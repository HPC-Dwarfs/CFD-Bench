/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The properties the embedded-boundary pressure operator has to have, checked
 * against a real solver on a real obstacle:
 *
 *   - symmetry, with and without a body, since that is what a conjugate
 *     gradient solver will need and it must not depend on the face weights;
 *   - a seven-point stencil, so no obstacle corner or edge introduces a
 *     diagonal coupling;
 *   - solid cells are identity rows whose pressure is exactly zero, including
 *     the deep interior cells the surface list deliberately leaves out;
 *   - the residual norm describes the fluid, so growing the solid fraction
 *     while leaving the fluid problem alone does not change it;
 *   - the surface list scales with the obstacle's surface rather than with its
 *     volume.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "discretization.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "surface-list.h"
#include "util.h"

#define IMAX 24
#define JMAX 24
#define KMAX 24

static Parameter Params;
static Discretization D;
static Solver S;

static void setup(CommType *base, const char *geometry, int boundary)
{
  initParameter(&Params);

  Params.name         = "check-obstacle";
  Params.imax         = IMAX;
  Params.jmax         = JMAX;
  Params.kmax         = KMAX;
  Params.xlength      = 4.0;
  Params.ylength      = 4.0;
  Params.zlength      = 4.0;
  /* Tight enough that the identity rows are unambiguous, loose enough that the
   * driver stays a check rather than a benchmark: this problem bottoms out near
   * 1e-16 in the mean square anyway. */
  Params.eps          = 1e-8;
  Params.omg          = 1.7;
  Params.itermax      = 20000;
  Params.levels       = 2;
  Params.presmooth    = 4;
  Params.postsmooth   = 4;
  Params.re           = 100.0;
  Params.tau          = 0.5;
  Params.gamma        = 0.9;
  Params.dt           = 0.02;
  Params.te           = 0.0;
  Params.gx = Params.gy = Params.gz = 0.0;
  Params.u_init = Params.v_init = Params.w_init = Params.p_init = 0.0;
  Params.geometryFile = (char *)geometry;

  Params.bcLeft = Params.bcRight = boundary;
  Params.bcBottom = Params.bcTop = boundary;
  Params.bcFront = Params.bcBack = boundary;

  D.comm = *base;
  commPartition(&D.comm, KMAX, JMAX, IMAX);
  initDiscretization(&D, &Params);
  initSolver(&S, &D, &Params);
  initProfiler(&D.comm);
}

/* Apply the operator once: out = A x, over the fluid unknowns. */
static void applyOperator(double *x, double *out)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;

  const double *Ax     = D.Ax;
  const double *Ay     = D.Ay;
  const double *Az     = D.Az;
  const double *Lambda = D.Lambda;

  double idx2          = 1.0 / (D.grid.dx * D.grid.dx);
  double idy2          = 1.0 / (D.grid.dy * D.grid.dy);
  double idz2          = 1.0 / (D.grid.dz * D.grid.dz);

  double *p            = x;

  commExchange(&D.comm, x);
  pressureBcApply(&S.bc, &D.comm, x, imaxLocal, jmaxLocal, kmaxLocal);

  size_t size = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);
  for (size_t i = 0; i < size; i++) {
    out[i] = 0.0;
  }

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }

        double pc = P(i, j, k);
        out[(size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
            (size_t)j * (imaxLocal + 2) + (size_t)i] =
            AX(i, j, k) * (P(i + 1, j, k) - pc) * idx2 +
            AX(i - 1, j, k) * (P(i - 1, j, k) - pc) * idx2 +
            AY(i, j, k) * (P(i, j + 1, k) - pc) * idy2 +
            AY(i, j - 1, k) * (P(i, j - 1, k) - pc) * idy2 +
            AZ(i, j, k) * (P(i, j, k + 1) - pc) * idz2 +
            AZ(i, j, k - 1) * (P(i, j, k - 1) - pc) * idz2;
      }
    }
  }
}

static double dotFluid(const double *a, const double *b)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;
  double sum           = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        size_t idx = (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                     (size_t)j * (imaxLocal + 2) + (size_t)i;
        sum += a[idx] * b[idx];
      }
    }
  }

  commReduceAll(&sum, SUM);
  return sum;
}

/* A count summed over the ranks. A surface list is per-rank, so a rank that
 * holds none of the body legitimately has an empty one. */
static double globalCount(int local)
{
  double v = (double)local;
  commReduceAll(&v, SUM);
  return v;
}

static void checkSymmetry(const char *label)
{
  int imaxLocal = D.comm.imaxLocal;
  int jmaxLocal = D.comm.jmaxLocal;
  int kmaxLocal = D.comm.kmaxLocal;
  size_t size   = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

  double *x     = calloc(size, sizeof(double));
  double *y     = calloc(size, sizeof(double));
  double *ax    = calloc(size, sizeof(double));
  double *ay    = calloc(size, sizeof(double));

  srand(12345);
  for (size_t i = 0; i < size; i++) {
    x[i] = (double)rand() / RAND_MAX - 0.5;
    y[i] = (double)rand() / RAND_MAX - 0.5;
  }

  applyOperator(x, ax);
  applyOperator(y, ay);

  double xAy   = dotFluid(x, ay);
  double yAx   = dotFluid(y, ax);
  double scale = fabs(xAy) + fabs(yAx) + 1e-300;

  CHECK_TRUE(fabs(xAy - yAx) / scale < 1e-12,
      "%s: operator is not symmetric, x.Ay = %.17g against y.Ax = %.17g",
      label,
      xAy,
      yAx);

  free(x);
  free(y);
  free(ax);
  free(ay);
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("obstacle");

  /* Symmetry without a body: the property must not come from the geometry. */
  setup(&comm, NULL, NOSLIP);
  checkSymmetry("obstacle-free");

  CHECK_TRUE(globalCount(S.surface.count) == 0.0,
      "an obstacle-free domain produced %.0f surface cells",
      globalCount(S.surface.count));

  /* Symmetry with one, and everything else the body implies. */
  setup(&comm, "sphere:2.0,2.0,2.0,1.0", NOSLIP);
  checkSymmetry("sphere");

  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  size_t size          = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);
  const double *Lambda = D.Lambda;

  CHECK_TRUE(globalCount(S.surface.count) > 0.0, "a sphere produced no surface cells");

  /* Seven-point stencil: a unit vector at one cell next to the body must
   * produce a result only at that cell and its six face neighbours. */
  {
    double *x  = calloc(size, sizeof(double));
    double *ax = calloc(size, sizeof(double));

    /* A fluid cell with at least one closed face, so it sits on the body. */
    int ti = -1, tj = -1, tk = -1;
    for (int k = 2; k < kmaxLocal && ti < 0; k++) {
      for (int j = 2; j < jmaxLocal && ti < 0; j++) {
        for (int i = 2; i < imaxLocal && ti < 0; i++) {
          if (LAM(i, j, k) > 0.0 && D.Ax[(size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                                         (size_t)j * (imaxLocal + 2) + (size_t)i] == 0.0) {
            ti = i;
            tj = j;
            tk = k;
          }
        }
      }
    }

    /* applyOperator and the reduction below are collective, so every rank has to
     * reach them whether or not it happens to hold a cell next to the body. A
     * rank that holds none simply contributes nothing. */
    int strideJ   = imaxLocal + 2;
    int strideK   = (imaxLocal + 2) * (jmaxLocal + 2);
    size_t centre = 0;
    int haveCell  = (ti > 0);

    if (haveCell) {
      centre    = (size_t)tk * strideK + (size_t)tj * strideJ + (size_t)ti;
      x[centre] = 1.0;
    }

    applyOperator(x, ax);

    double strays = 0.0;
    double found  = haveCell ? 1.0 : 0.0;

    if (haveCell) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          for (int i = 1; i < imaxLocal + 1; i++) {
            size_t idx = (size_t)k * strideK + (size_t)j * strideJ + (size_t)i;

            if (idx == centre || idx == centre + 1 || idx == centre - 1 ||
                idx == centre + strideJ || idx == centre - strideJ ||
                idx == centre + strideK || idx == centre - strideK) {
              continue;
            }

            if (ax[idx] != 0.0) {
              strays += 1.0;
            }
          }
        }
      }
    }

    commReduceAll(&strays, SUM);
    commReduceAll(&found, SUM);

    CHECK_TRUE(found > 0.0, "no rank holds a fluid cell next to the body");
    CHECK_TRUE(strays == 0.0,
        "the operator couples a cell next to the body to %.0f cells that are not its "
        "six face neighbours",
        strays);

    free(x);
    free(ax);
  }

  /* Solve, then check the identity rows. */
  {
    int offsets[NDIMS] = { 0, 0, 0 };
    commGetOffsets(&D.comm, offsets, KMAX, JMAX, IMAX);

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          double x = ((i - 1 + offsets[IDIM]) + 0.5) / (double)IMAX;
          double y = ((j - 1 + offsets[JDIM]) + 0.5) / (double)JMAX;
          double z = ((k - 1 + offsets[KDIM]) + 0.5) / (double)KMAX;

          double *rhs = D.rhs;
          RHS(i, j, k) =
              LAM(i, j, k) * cos(M_PI * x) * cos(M_PI * y) * cos(M_PI * z);
        }
      }
    }

    for (size_t i = 0; i < size; i++) {
      D.p[i] = 0.0;
    }

    PressureLevelType lv;
    pressureLevelFromSolver(&S, &lv);
    double before = pressureResidualNorm(&lv, D.p, D.rhs);

    double res    = solve(&S, D.p, D.rhs);

    /* An absolute tolerance would be a statement about double precision rather
     * than about the solver: this problem bottoms out around 1e-16 in the mean
     * square, which is where the operator stops being able to tell the residual
     * from rounding. The reduction is the thing worth asserting. */
    CHECK_TRUE(res < before * 1e-10,
        "the obstacle solve reduced the residual only from %.3e to %.3e",
        before,
        res);

    /* Every solid cell, not only the ones on the surface list, must be exactly
     * zero. The deep interior is the interesting half: nothing corrects it, and
     * it stays right only because relaxing a zero neighbourhood with a zero
     * right-hand side is exactly zero. */
    int nonzeroSolid = 0, deepSolid = 0, deepNonzero = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) > 0.0) {
            continue;
          }

          double *p = D.p;
          if (P(i, j, k) != 0.0) {
            ++nonzeroSolid;
          }

          int surrounded = LAM(i + 1, j, k) == 0.0 && LAM(i - 1, j, k) == 0.0 &&
                           LAM(i, j + 1, k) == 0.0 && LAM(i, j - 1, k) == 0.0 &&
                           LAM(i, j, k + 1) == 0.0 && LAM(i, j, k - 1) == 0.0;

          if (surrounded) {
            ++deepSolid;
            if (P(i, j, k) != 0.0) {
              ++deepNonzero;
            }
          }
        }
      }
    }

    double totals[3] = { (double)nonzeroSolid, (double)deepSolid, (double)deepNonzero };
    for (int i = 0; i < 3; i++) {
      commReduceAll(&totals[i], SUM);
    }

    CHECK_TRUE(totals[0] == 0.0,
        "%.0f solid cells hold a nonzero pressure after a converged solve",
        totals[0]);
    CHECK_TRUE(totals[1] > 0.0,
        "the sphere has no solid cell away from its surface, so the "
        "self-maintaining interior is untested");
    CHECK_TRUE(totals[2] == 0.0,
        "%.0f solid cells away from the surface drifted off zero",
        totals[2]);
  }

  /*
   * Perturbing the pressure inside the body must not move the fluid. This is
   * the scenario the aperture factor in the velocity correction exists for: a
   * face touching a solid cell is closed, so no pressure gradient crosses it.
   */
  {
    double *uBefore = calloc(size, sizeof(double));
    double *vBefore = calloc(size, sizeof(double));
    double *wBefore = calloc(size, sizeof(double));

    for (size_t i = 0; i < size; i++) {
      D.p[i]   = 0.1 * (double)(i % 7);
      D.f[i]   = 0.01 * (double)(i % 5);
      D.g[i]   = 0.02 * (double)(i % 3);
      D.h[i]   = 0.03 * (double)(i % 11);
    }

    adaptUV(&D);
    memcpy(uBefore, D.u, size * sizeof(double));
    memcpy(vBefore, D.v, size * sizeof(double));
    memcpy(wBefore, D.w, size * sizeof(double));

    /* Now put something arbitrary inside the body and correct again. */
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) == 0.0) {
            double *p = D.p;
            P(i, j, k) = 1234.5;
          }
        }
      }
    }

    adaptUV(&D);

    double moved = 0.0;
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          size_t idx = (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                       (size_t)j * (imaxLocal + 2) + (size_t)i;
          if (D.u[idx] != uBefore[idx] || D.v[idx] != vBefore[idx] ||
              D.w[idx] != wBefore[idx]) {
            moved += 1.0;
          }
        }
      }
    }
    commReduceAll(&moved, SUM);

    CHECK_TRUE(moved == 0.0,
        "perturbing the pressure inside the body moved %.0f velocity components",
        moved);

    free(uBefore);
    free(vBefore);
    free(wBefore);
  }

  /*
   * The residual describes the fluid problem. Whatever is left lying in the
   * solid cells -- pressure or right-hand side -- must not reach it, which is
   * what makes the reported norm independent of how much of the domain is
   * solid.
   */
  {
    PressureLevelType lv;
    pressureLevelFromSolver(&S, &lv);

    for (size_t i = 0; i < size; i++) {
      D.p[i]   = 0.05 * (double)(i % 13);
      D.rhs[i] = 0.07 * (double)(i % 9);
    }

    double clean = pressureResidualNorm(&lv, D.p, D.rhs);

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) > 0.0) {
            continue;
          }
          size_t idx = (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                       (size_t)j * (imaxLocal + 2) + (size_t)i;
          D.rhs[idx] = 999.0;
        }
      }
    }

    double dirtied = pressureResidualNorm(&lv, D.p, D.rhs);

    CHECK_NEAR(dirtied, clean, 0.0,
        "the right-hand side in the solid cells changed the reported residual");
  }

  /*
   * The surface list has to scale with the body's surface, not its volume.
   * That is the whole reason the interior sweep can stay geometry-free: a body
   * twice as wide has four times the surface and eight times the volume, so a
   * list that tracked volume would grow twice as fast.
   */
  {
    double small, large;

    setup(&comm, "sphere:2.0,2.0,2.0,0.5", NOSLIP);
    small = globalCount(S.surface.count);

    setup(&comm, "sphere:2.0,2.0,2.0,1.0", NOSLIP);
    large = globalCount(S.surface.count);

    double ratio = large / (small > 0.0 ? small : 1.0);

    CHECK_TRUE(small > 0.0 && large > 0.0,
        "one of the spheres produced no surface cells (%.0f and %.0f)",
        small,
        large);

    /* Doubling the radius should roughly quadruple the list and would octuple
     * it if it tracked volume. Allow a wide band -- the point is to tell an
     * area law from a volume law, not to measure the constant. */
    CHECK_TRUE(ratio > 2.5 && ratio < 6.0,
        "doubling the radius changed the surface list by a factor of %.2f; an area "
        "law gives about 4 and a volume law about 8",
        ratio);

    if (commIsMaster(&D.comm)) {
      printf("surface list: %.0f cells at radius 0.5, %.0f at radius 1.0, ratio %.2f\n",
          small,
          large,
          ratio);
    }
  }

  /*
   * Geometries a cell-type model cannot express. Both have to run; the plate in
   * particular has to stop the flow while leaving the cells on both sides
   * fluid, which is the case a destination-cell test would miss entirely.
   */
  {
    setup(&comm, "plate-x:2.0,1.0,1.0,3.0,3.0", NOSLIP);

    int closed = 0, solids = 0;
    int im = D.comm.imaxLocal, jm = D.comm.jmaxLocal, km = D.comm.kmaxLocal;

    for (int k = 1; k < km + 1; k++) {
      for (int j = 1; j < jm + 1; j++) {
        for (int i = 1; i < im + 1; i++) {
          size_t idx = (size_t)k * (im + 2) * (jm + 2) + (size_t)j * (im + 2) + (size_t)i;
          if (D.Ax[idx] == 0.0) {
            ++closed;
          }
          if (D.Lambda[idx] == 0.0) {
            ++solids;
          }
        }
      }
    }

    double t[2] = { (double)closed, (double)solids };
    commReduceAll(&t[0], SUM);
    commReduceAll(&t[1], SUM);

    CHECK_TRUE(t[0] > 0.0, "the plate closed no faces");
    CHECK_TRUE(t[1] == 0.0,
        "the plate produced %.0f solid cells; a zero-thickness body has none",
        t[1]);
    CHECK_TRUE(globalCount(S.surface.count) > 0.0, "the plate produced no surface cells");

    /* One isolated solid cell: six closed faces and no admissible mirror, the
     * configuration that has no representation at all as a cell type. */
    double dx = 4.0 / IMAX;
    char spec[128];
    snprintf(spec,
        sizeof(spec),
        "box:%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",
        10 * dx,
        10 * dx,
        10 * dx,
        11 * dx,
        11 * dx,
        11 * dx);
    setup(&comm, spec, NOSLIP);

    double solid = 0.0;
    for (int k = 1; k < D.comm.kmaxLocal + 1; k++) {
      for (int j = 1; j < D.comm.jmaxLocal + 1; j++) {
        for (int i = 1; i < D.comm.imaxLocal + 1; i++) {
          size_t idx = (size_t)k * (D.comm.imaxLocal + 2) * (D.comm.jmaxLocal + 2) +
                       (size_t)j * (D.comm.imaxLocal + 2) + (size_t)i;
          if (D.Lambda[idx] == 0.0) {
            solid += 1.0;
          }
        }
      }
    }
    commReduceAll(&solid, SUM);

    CHECK_TRUE(solid == 1.0,
        "the single-cell box produced %.0f solid cells, expected exactly 1",
        solid);
    CHECK_TRUE(globalCount(S.surface.count) > 0.0,
        "an isolated solid cell produced no surface cells");
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
