/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The matrix-free operator application, checked directly.
 *
 * tests/checks/obstacle.c checks the same properties of the operator as a
 * mathematical object, using its own geometry-reading apply. This driver checks
 * pressureApplyOperator -- the routine a Krylov solver actually calls, built
 * from a geometry-free bulk pass plus the surface list -- which is a different
 * claim: that the fast path computes what the slow one does.
 *
 *   - symmetry, x.Ay == y.Ax, with and without a body, because a conjugate
 *     gradient method needs it and the split must not break it;
 *   - positive definiteness on the fluid subspace, since the application
 *     negates the assembled operator precisely to obtain it;
 *   - a seven-point footprint from a unit vector next to an obstacle corner, so
 *     the surface pass introduces no coupling the stencil does not have;
 *   - agreement with pressureResidualNorm, which reads the geometry arrays, at
 *     cut cells -- the check that keeps pressureApplySurface from drifting away
 *     from pressureCorrectSurface.
 *
 * It lives entirely in solverbase.c, which every variant links, so it runs
 * unmodified under all of them and needs no #ifdef seam.
 */
#include <math.h>
#include <stdlib.h>

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
static PressureLevelType Lv;

static void setup(CommType *base, const char *geometry, int boundary)
{
  initParameter(&Params);

  Params.name         = "check-operator";
  Params.imax         = IMAX;
  Params.jmax         = JMAX;
  Params.kmax         = KMAX;
  Params.xlength      = 4.0;
  Params.ylength      = 4.0;
  Params.zlength      = 4.0;
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
  pressureLevelFromSolver(&S, &Lv);
}

static size_t fieldSize(void)
{
  return (size_t)(D.comm.imaxLocal + 2) * (D.comm.jmaxLocal + 2) *
         (D.comm.kmaxLocal + 2);
}

static size_t linear(int i, int j, int k)
{
  int imaxLocal = D.comm.imaxLocal;
  int jmaxLocal = D.comm.jmaxLocal;
  return (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) + (size_t)j * (imaxLocal + 2) +
         (size_t)i;
}

/* Over the fluid unknowns only, which is where the operator is definite and the
 * only place a Krylov method forms an inner product. */
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
        sum += a[linear(i, j, k)] * b[linear(i, j, k)];
      }
    }
  }

  commReduceAll(&sum, SUM);
  return sum;
}

/* A field that is zero in the solid, which is the invariant every vector a
 * Krylov iteration carries has to satisfy. */
static void randomFluidField(double *x, unsigned seed)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;
  size_t size          = fieldSize();

  for (size_t i = 0; i < size; i++) {
    x[i] = 0.0;
  }

  /* Seeded from the global position, not from rand(), so the field is the same
   * however the domain is divided and the checks mean the same thing at every
   * rank count. */
  int offsets[NDIMS] = { 0, 0, 0 };
  commGetOffsets(&D.comm, offsets, KMAX, JMAX, IMAX);

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }

        unsigned gi = (unsigned)(i - 1 + offsets[IDIM]);
        unsigned gj = (unsigned)(j - 1 + offsets[JDIM]);
        unsigned gk = (unsigned)(k - 1 + offsets[KDIM]);
        unsigned h  = seed + gi * 73856093u + gj * 19349663u + gk * 83492791u;

        h ^= h >> 13;
        h *= 1274126177u;
        h ^= h >> 16;

        x[linear(i, j, k)] = (double)h / (double)0xffffffffu - 0.5;
      }
    }
  }
}

static void checkSymmetry(const char *label)
{
  size_t size = fieldSize();
  double *x   = calloc(size, sizeof(double));
  double *y   = calloc(size, sizeof(double));
  double *ax  = calloc(size, sizeof(double));
  double *ay  = calloc(size, sizeof(double));

  randomFluidField(x, 12345u);
  randomFluidField(y, 67890u);

  pressureApplyOperator(&Lv, x, ax);
  pressureApplyOperator(&Lv, y, ay);

  double xAy   = dotFluid(x, ay);
  double yAx   = dotFluid(y, ax);
  double scale = fabs(xAy) + fabs(yAx) + 1e-300;

  CHECK_TRUE(fabs(xAy - yAx) / scale < 1e-12,
      "%s: the application is not symmetric, x.Ay = %.17g against y.Ax = %.17g",
      label,
      xAy,
      yAx);

  /* Positive definiteness on the fluid subspace. The application negates the
   * assembled operator for exactly this reason, so a negative quadratic form
   * here means a Krylov method would break down rather than converge slowly. */
  double xAx = dotFluid(x, ax);

  CHECK_TRUE(xAx > 0.0,
      "%s: the quadratic form is %.17g, so the negated operator is not positive "
      "definite on the fluid subspace",
      label,
      xAx);

  free(x);
  free(y);
  free(ax);
  free(ay);
}

/*
 * The residual of an arbitrary field, computed the way the solvers report it,
 * against the same residual assembled from the application.
 *
 * pressureResidualNorm reads Ax, Ay, Az and Lambda per cell; the application
 * reads no geometry in the bulk and the surface list at cut cells. They are two
 * expressions of one operator, and this is what keeps them from drifting apart.
 *
 *   residual   = Lambda*rhs - A p        (assembled sign)
 *   application = -A p                   (negated sign)
 *
 * so  sum over fluid of (Lambda*rhs + apply)^2 / fluidCells  must equal what
 * pressureResidualNorm returns.
 */
static void checkAgreesWithResidual(const char *label)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;

  size_t size          = fieldSize();
  double *p            = calloc(size, sizeof(double));
  double *ap           = calloc(size, sizeof(double));
  double *rhs          = D.rhs;

  randomFluidField(p, 24680u);

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        RHS(i, j, k) = 0.01 * (double)((i + 3 * j + 7 * k) % 13);
      }
    }
  }

  /* pressureResidualNorm exchanges and applies the boundary condition to p, and
   * so does the application; both see the same field either way. */
  double reported = pressureResidualNorm(&Lv, p, D.rhs);

  pressureApplyOperator(&Lv, p, ap);

  double sum      = 0.0;
  double worstCut = 0.0;
  int cutCells    = 0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          /* A solid cell is an identity row: the application must write exactly
           * zero there, and it contributes nothing to the residual. */
          CHECK_NEAR(ap[linear(i, j, k)],
              0.0,
              0.0,
              "%s: the application left a solid cell nonzero",
              label);
          continue;
        }

        double r = LAM(i, j, k) * RHS(i, j, k) + ap[linear(i, j, k)];
        sum += r * r;

        /* A cut cell: at least one face closed, so the bulk pass got it wrong
         * and the surface pass had to replace it. This is where a divergence
         * between the apply coefficients and the relaxation coefficients would
         * show up, and nowhere else. */
        int cut = LAM(i, j, k) < 1.0 || LAM(i + 1, j, k) == 0.0 ||
                  LAM(i - 1, j, k) == 0.0 || LAM(i, j + 1, k) == 0.0 ||
                  LAM(i, j - 1, k) == 0.0 || LAM(i, j, k + 1) == 0.0 ||
                  LAM(i, j, k - 1) == 0.0;

        if (cut) {
          ++cutCells;
          if (fabs(r) > worstCut) {
            worstCut = fabs(r);
          }
        }
      }
    }
  }

  commReduceAll(&sum, SUM);
  double assembled  = sum / S.fluidCells;

  double cutTotal   = (double)cutCells;
  commReduceAll(&cutTotal, SUM);

  double scale = fabs(reported) + fabs(assembled) + 1e-300;

  CHECK_TRUE(fabs(reported - assembled) / scale < 1e-12,
      "%s: the application implies a residual norm of %.17g where the residual "
      "routine reports %.17g",
      label,
      assembled,
      reported);

  /* Restating the comparison per cell rather than in the sum, so that a
   * cancelling pair of errors at two cut cells cannot hide. */
  {
    double worst = 0.0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) == 0.0) {
            continue;
          }

          double pc = p[linear(i, j, k)];
          double direct =
              -(D.Ax[linear(i, j, k)] * (p[linear(i + 1, j, k)] - pc) /
                      (D.grid.dx * D.grid.dx) +
                  D.Ax[linear(i - 1, j, k)] * (p[linear(i - 1, j, k)] - pc) /
                      (D.grid.dx * D.grid.dx) +
                  D.Ay[linear(i, j, k)] * (p[linear(i, j + 1, k)] - pc) /
                      (D.grid.dy * D.grid.dy) +
                  D.Ay[linear(i, j - 1, k)] * (p[linear(i, j - 1, k)] - pc) /
                      (D.grid.dy * D.grid.dy) +
                  D.Az[linear(i, j, k)] * (p[linear(i, j, k + 1)] - pc) /
                      (D.grid.dz * D.grid.dz) +
                  D.Az[linear(i, j, k - 1)] * (p[linear(i, j, k - 1)] - pc) /
                      (D.grid.dz * D.grid.dz));

          double diff = fabs(direct - ap[linear(i, j, k)]);
          if (diff > worst) {
            worst = diff;
          }
        }
      }
    }

    commReduceAll(&worst, MAX);

    CHECK_TRUE(worst < 1e-9,
        "%s: the application disagrees with the geometry-reading operator by %.3e at "
        "some cell",
        label,
        worst);
  }

  if (commIsMaster(&D.comm)) {
    printf("%s: residual %.3e from the application against %.3e reported, %.0f cut "
           "cells, worst cut-cell residual %.3e\n",
        label,
        assembled,
        reported,
        cutTotal,
        worstCut);
  }

  free(p);
  free(ap);
}

/* A unit vector at one cell may reach that cell and its six face neighbours and
 * nothing else, however the body cuts the cell. */
static void checkFootprint(void)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;

  size_t size          = fieldSize();
  double *x            = calloc(size, sizeof(double));
  double *ax           = calloc(size, sizeof(double));

  /* A corner of the body: a fluid cell with two closed faces in different
   * directions, which is where a nine-point coupling would appear if one were
   * ever introduced. */
  int ti = -1, tj = -1, tk = -1;

  for (int k = 2; k < kmaxLocal && ti < 0; k++) {
    for (int j = 2; j < jmaxLocal && ti < 0; j++) {
      for (int i = 2; i < imaxLocal && ti < 0; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }

        int closedI = LAM(i + 1, j, k) == 0.0 || LAM(i - 1, j, k) == 0.0;
        int closedJ = LAM(i, j + 1, k) == 0.0 || LAM(i, j - 1, k) == 0.0;
        int closedK = LAM(i, j, k + 1) == 0.0 || LAM(i, j, k - 1) == 0.0;

        if (closedI + closedJ + closedK >= 2) {
          ti = i;
          tj = j;
          tk = k;
        }
      }
    }
  }

  /* The application and the reduction are collective, so every rank reaches
   * them whether or not it happens to hold a corner. */
  int strideJ   = imaxLocal + 2;
  int strideK   = (imaxLocal + 2) * (jmaxLocal + 2);
  size_t centre = 0;
  int haveCell  = (ti > 0);

  if (haveCell) {
    centre    = linear(ti, tj, tk);
    x[centre] = 1.0;
  }

  pressureApplyOperator(&Lv, x, ax);

  double strays = 0.0;
  double found  = haveCell ? 1.0 : 0.0;

  if (haveCell) {
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          size_t idx = linear(i, j, k);

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

  CHECK_TRUE(found > 0.0, "no rank holds a fluid cell at a corner of the body");
  CHECK_TRUE(strays == 0.0,
      "the application couples a cell at a corner of the body to %.0f cells that are "
      "not its six face neighbours",
      strays);

  free(x);
  free(ax);
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  /* commFinalize frees the derived datatypes, which only commPartition creates,
   * so the communicator this driver finalizes has to have been partitioned. */
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("operator");

  /* Without a body: the split must be right before the surface list is even
   * populated, so that a failure below is attributable to the surface pass. */
  setup(&comm, NULL, NOSLIP);
  checkSymmetry("obstacle-free");
  checkAgreesWithResidual("obstacle-free");

  /* With one, and with an outflow boundary as well, so that the odd reflection
   * is exercised alongside the mirrored halo -- both have to leave the operator
   * symmetric or a Krylov method has no business running on it. */
  setup(&comm, "sphere:2.0,2.0,2.0,1.0", NOSLIP);
  checkSymmetry("sphere, walls");
  checkAgreesWithResidual("sphere, walls");
  checkFootprint();

  setup(&comm, "sphere:2.0,2.0,2.0,1.0", OUTFLOW);
  checkSymmetry("sphere, outflow");
  checkAgreesWithResidual("sphere, outflow");

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
