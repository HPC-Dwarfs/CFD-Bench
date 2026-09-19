/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Exercises the multigrid transfer operators and cycle one step at a time.
 * Only meaningful when the multigrid solver is the one linked in, so it skips
 * under the others.
 */
#include <math.h>
#include <stdlib.h>

#include "check.h"

#if !defined(SOLVER_mg)

#include <stdio.h>

int main(void)
{
  printf("== multigrid ==\nmultigrid: skipped, linked solver is not mg\n");
  return EXIT_SUCCESS;
}

#else

#include "discretization.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "util.h"

#define IMAX 32
#define JMAX 16
#define KMAX 16

#define AT(f, i, j, k, im, jm)                                                           \
  f[(size_t)(k) * ((im) + 2) * ((jm) + 2) + (size_t)(j) * ((im) + 2) + (size_t)(i)]

static Parameter Params;
static Discretization D;
static Solver S;

static void setup(CommType *base, int boundary)
{
  initParameter(&Params);

  Params.name       = "check-multigrid";
  Params.imax       = IMAX;
  Params.jmax       = JMAX;
  Params.kmax       = KMAX;
  Params.xlength    = 2.0;
  Params.ylength    = 1.0;
  Params.zlength    = 1.0;
  Params.eps        = 1e-11;
  Params.omg        = 1.7;
  Params.itermax    = 200;
  Params.levels     = 3;
  Params.presmooth  = 4;
  Params.postsmooth = 4;
  Params.re         = 100.0;
  Params.tau        = 0.5;
  Params.gamma      = 0.9;
  Params.dt         = 0.02;
  Params.te         = 0.0;
  Params.gx = Params.gy = Params.gz = 0.0;
  Params.u_init = Params.v_init = Params.w_init = Params.p_init = 0.0;

  Params.bcLeft = Params.bcRight = boundary;
  Params.bcBottom = Params.bcTop = boundary;
  Params.bcFront = Params.bcBack = boundary;

  D.comm = *base;
  commPartition(&D.comm, KMAX, JMAX, IMAX);
  initDiscretization(&D, &Params);
  initSolver(&S, &D, &Params);
  initProfiler(&D.comm);
}

static void setupWithGeometry(CommType *base, const char *geometry)
{
  Params.geometryFile = (char *)geometry;
  D.comm              = *base;
  commPartition(&D.comm, KMAX, JMAX, IMAX);
  initDiscretization(&D, &Params);
  initSolver(&S, &D, &Params);
}

static void fill(double *f, int im, int jm, int km, double value)
{
  size_t size = (size_t)(im + 2) * (jm + 2) * (km + 2);

  for (size_t i = 0; i < size; i++) {
    f[i] = value;
  }
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);

  CHECK_BEGIN("multigrid");

  setup(&comm, NOSLIP);

  int levels = mgTestLevels(&S);
  CHECK_TRUE(levels >= 2, "need at least two levels to test transfers, have %d", levels);

  if (levels >= 2) {
    int fi, fj, fk, ci, cj, ck;
    mgTestLevelExtents(&S, 0, &fi, &fj, &fk);
    mgTestLevelExtents(&S, 1, &ci, &cj, &ck);

    double fdx, fdy, fdz, cdx, cdy, cdz;
    mgTestLevelMesh(&S, 0, &fdx, &fdy, &fdz);
    mgTestLevelMesh(&S, 1, &cdx, &cdy, &cdz);

    /* Each level has its own mesh: the coarse spacing is twice the fine one.
     * Both used to be the finest spacing, so a coarse operator was wrong by a
     * factor of four per level. */
    CHECK_NEAR(cdx, 2.0 * fdx, 1e-15, "coarse dx is not twice the fine dx");
    CHECK_NEAR(cdy, 2.0 * fdy, 1e-15, "coarse dy is not twice the fine dy");
    CHECK_NEAR(cdz, 2.0 * fdz, 1e-15, "coarse dz is not twice the fine dz");

    double *fineR   = mgTestLevelR(&S, 0);
    double *coarseR = mgTestLevelR(&S, 1);

    /* Restriction of a constant is that constant. */
    fill(fineR, fi, fj, fk, 3.25);
    fill(coarseR, ci, cj, ck, -999.0);
    mgTestRestrict(&S, 0);

    double worst = 0.0;
    for (int k = 1; k < ck + 1; k++) {
      for (int j = 1; j < cj + 1; j++) {
        for (int i = 1; i < ci + 1; i++) {
          double d = fabs(AT(coarseR, i, j, k, ci, cj) - 3.25);
          if (d > worst) {
            worst = d;
          }
        }
      }
    }
    commReduceAll(&worst, MAX);
    CHECK_NEAR(worst, 0.0, 1e-15, "restriction of a constant is not that constant");

    /* Restriction lands in the cells it should: mark exactly the eight fine
     * children of one coarse cell and check only that coarse cell sees them. */
    fill(fineR, fi, fj, fk, 0.0);
    fill(coarseR, ci, cj, ck, -999.0);
    int ti = 2, tj = 2, tk = 2;
    for (int dk = 0; dk < 2; dk++) {
      for (int dj = 0; dj < 2; dj++) {
        for (int di = 0; di < 2; di++) {
          AT(fineR, 2 * ti - 1 + di, 2 * tj - 1 + dj, 2 * tk - 1 + dk, fi, fj) = 8.0;
        }
      }
    }
    mgTestRestrict(&S, 0);
    CHECK_NEAR(AT(coarseR, ti, tj, tk, ci, cj),
        8.0,
        1e-15,
        "the eight fine children did not land in their own coarse cell");

    int strays = 0;
    for (int k = 1; k < ck + 1; k++) {
      for (int j = 1; j < cj + 1; j++) {
        for (int i = 1; i < ci + 1; i++) {
          if (i == ti && j == tj && k == tk) {
            continue;
          }
          if (fabs(AT(coarseR, i, j, k, ci, cj)) > 1e-15) {
            ++strays;
          }
        }
      }
    }
    CHECK_TRUE(strays == 0, "%d coarse cells outside the target received a value", strays);

    /* Prolongation defines every fine cell: pre-fill with a sentinel and check
     * none of it survives. It used to be piecewise injection into the wrong
     * array, so half the fine cells kept whatever was there before. */
    double *fineE   = mgTestLevelE(&S, 0);
    double *coarseE = mgTestLevelE(&S, 1);

    fill(coarseE, ci, cj, ck, 2.5);
    fill(fineE, fi, fj, fk, -12345.0);
    mgTestProlongate(&S, 0);

    int sentinels = 0;
    for (int k = 1; k < fk + 1; k++) {
      for (int j = 1; j < fj + 1; j++) {
        for (int i = 1; i < fi + 1; i++) {
          if (AT(fineE, i, j, k, fi, fj) == -12345.0) {
            ++sentinels;
          }
        }
      }
    }
    CHECK_TRUE(sentinels == 0, "%d fine cells kept the sentinel after prolongation",
        sentinels);

    /* Restriction then prolongation of a constant returns that constant, away
     * from the boundary where the wall condition mirrors it anyway. */
    fill(fineR, fi, fj, fk, 1.75);
    mgTestRestrict(&S, 0);
    fill(coarseE, ci, cj, ck, 0.0);
    for (int k = 1; k < ck + 1; k++) {
      for (int j = 1; j < cj + 1; j++) {
        for (int i = 1; i < ci + 1; i++) {
          AT(coarseE, i, j, k, ci, cj) = AT(coarseR, i, j, k, ci, cj);
        }
      }
    }
    fill(fineE, fi, fj, fk, 0.0);
    mgTestProlongate(&S, 0);

    worst = 0.0;
    for (int k = 1; k < fk + 1; k++) {
      for (int j = 1; j < fj + 1; j++) {
        for (int i = 1; i < fi + 1; i++) {
          double d = fabs(AT(fineE, i, j, k, fi, fj) - 1.75);
          if (d > worst) {
            worst = d;
          }
        }
      }
    }
    commReduceAll(&worst, MAX);
    CHECK_NEAR(worst, 0.0, 1e-14,
        "restriction followed by prolongation did not return the constant");
  }

  /* A constant is a fixed point of a full V-cycle with a zero right-hand side
   * under the all-wall condition: the operator annihilates constants and the
   * boundary condition preserves them, at every level. This only holds when the
   * condition is applied on the coarse levels too. */
  {
    int fi, fj, fk;
    mgTestLevelExtents(&S, 0, &fi, &fj, &fk);
    size_t size = (size_t)(fi + 2) * (fj + 2) * (fk + 2);

    for (size_t i = 0; i < size; i++) {
      D.p[i]   = 4.0;
      D.rhs[i] = 0.0;
    }

    mgTestVcycle(&S, D.p, D.rhs);

    double worst = 0.0;
    for (int k = 1; k < fk + 1; k++) {
      for (int j = 1; j < fj + 1; j++) {
        for (int i = 1; i < fi + 1; i++) {
          double d = fabs(AT(D.p, i, j, k, fi, fj) - 4.0);
          if (d > worst) {
            worst = d;
          }
        }
      }
    }
    commReduceAll(&worst, MAX);
    CHECK_NEAR(worst, 0.0, 1e-12, "a constant is not a fixed point of a V-cycle");
  }

  /* A V-cycle must reduce the residual by more than pre- plus post-smoothing
   * alone would, which is the whole point of the coarse correction. Before, the
   * correction was prolongated into the residual array while correct() read the
   * error array, so the cycle degenerated to smoothing. */
  {
    int fi, fj, fk;
    mgTestLevelExtents(&S, 0, &fi, &fj, &fk);
    double dx, dy, dz;
    mgTestLevelMesh(&S, 0, &dx, &dy, &dz);

    size_t size  = (size_t)(fi + 2) * (fj + 2) * (fk + 2);
    double cells = (double)IMAX * JMAX * KMAX;

    double *rhs  = calloc(size, sizeof(double));
    double *p0   = calloc(size, sizeof(double));

    PressureLevelType desc;
    pressureLevelFromSolver(&S, &desc);

    int offsets[NDIMS] = { 0, 0, 0 };
    commGetOffsets(&D.comm, offsets, KMAX, JMAX, IMAX);

    /* A smooth, mean-free right-hand side, so the singular all-wall system is
     * compatible. */
    for (int k = 1; k < fk + 1; k++) {
      for (int j = 1; j < fj + 1; j++) {
        for (int i = 1; i < fi + 1; i++) {
          double x = ((i - 1 + offsets[IDIM]) + 0.5) / (double)IMAX;
          double y = ((j - 1 + offsets[JDIM]) + 0.5) / (double)JMAX;
          double z = ((k - 1 + offsets[KDIM]) + 0.5) / (double)KMAX;
          AT(rhs, i, j, k, fi, fj) =
              cos(M_PI * x) * cos(M_PI * y) * cos(M_PI * z);
        }
      }
    }

    double before = pressureResidualNorm(&desc, p0, rhs);

    /* Smoothing only, with as many sweeps as a whole cycle spends. */
    for (size_t i = 0; i < size; i++) {
      D.p[i] = 0.0;
    }
    mgTestSmooth(&S, 0, D.p, rhs, Params.presmooth + Params.postsmooth);
    double smoothed = pressureResidualNorm(&desc, D.p, rhs);

    /* A full cycle, same starting point. */
    for (size_t i = 0; i < size; i++) {
      D.p[i] = 0.0;
    }
    mgTestVcycle(&S, D.p, rhs);
    double cycled = pressureResidualNorm(&desc, D.p, rhs);

    CHECK_TRUE(cycled < before, "a V-cycle did not reduce the residual (%.3e -> %.3e)",
        before, cycled);
    CHECK_TRUE(cycled < smoothed,
        "a V-cycle did no better than smoothing alone (%.3e against %.3e), so the "
        "coarse correction is not reaching the solution",
        cycled,
        smoothed);

    if (commIsMaster(&D.comm)) {
      printf("residual: start %.3e, smoothing only %.3e, full cycle %.3e\n",
          before,
          smoothed,
          cycled);
    }

    free(rhs);
    free(p0);
  }

  /*
   * A body resolved on the finest grid has to survive coarsening, otherwise the
   * coarse correction solves a different problem from the one it is correcting.
   * The apertures and volume fractions are coarsened geometrically -- four fine
   * faces to a coarse face, eight fine cells to a coarse cell -- so a body
   * several cells across stays represented for several levels.
   */
  {
    setupWithGeometry(&comm, "sphere:1.0,0.5,0.5,0.4");

    int levelsG = mgTestLevels(&S);

    for (int l = 0; l < levelsG; l++) {
      double solid   = (double)mgTestLevelSolidCount(&S, l);
      double surface = (double)mgTestLevelSurfaceCount(&S, l);
      commReduceAll(&solid, SUM);
      commReduceAll(&surface, SUM);

      CHECK_TRUE(solid > 0.0,
          "level %d of %d has no solid cells, so the body is not represented there",
          l,
          levelsG);
      CHECK_TRUE(surface > 0.0,
          "level %d of %d has an empty surface list despite holding the body",
          l,
          levelsG);

      if (commIsMaster(&D.comm)) {
        printf("level %d: %.0f solid cells, %.0f surface cells\n", l, solid, surface);
      }
    }
  }

  /*
   * A body too small to survive coarsening is not an error: the cycle still has
   * to converge, with the coarse levels simply blind to it. initSolver says so
   * on the run header rather than leaving it to be discovered as a convergence
   * problem.
   */
  {
    setupWithGeometry(&comm, "sphere:1.0,0.5,0.5,0.12");

    int levelsG = mgTestLevels(&S);
    double coarsest = (double)mgTestLevelSolidCount(&S, levelsG - 1);
    commReduceAll(&coarsest, SUM);

    CHECK_TRUE(coarsest == 0.0,
        "the small sphere was still resolved at the coarsest level, so the "
        "unresolved case is untested (%.0f solid cells)",
        coarsest);

    int fi, fj, fk;
    mgTestLevelExtents(&S, 0, &fi, &fj, &fk);
    size_t sz = (size_t)(fi + 2) * (fj + 2) * (fk + 2);

    PressureLevelType desc;
    pressureLevelFromSolver(&S, &desc);

    int offsets[NDIMS] = { 0, 0, 0 };
    commGetOffsets(&D.comm, offsets, KMAX, JMAX, IMAX);

    for (int k = 1; k < fk + 1; k++) {
      for (int j = 1; j < fj + 1; j++) {
        for (int i = 1; i < fi + 1; i++) {
          double x = ((i - 1 + offsets[IDIM]) + 0.5) / (double)IMAX;
          double y = ((j - 1 + offsets[JDIM]) + 0.5) / (double)JMAX;
          double z = ((k - 1 + offsets[KDIM]) + 0.5) / (double)KMAX;
          AT(D.rhs, i, j, k, fi, fj) =
              cos(M_PI * x) * cos(M_PI * y) * cos(M_PI * z);
        }
      }
    }

    for (size_t i = 0; i < sz; i++) {
      D.p[i] = 0.0;
    }

    double before = pressureResidualNorm(&desc, D.p, D.rhs);
    for (int c = 0; c < 20; c++) {
      mgTestVcycle(&S, D.p, D.rhs);
    }
    double after = pressureResidualNorm(&desc, D.p, D.rhs);

    CHECK_TRUE(after < before,
        "multigrid did not converge with a body that vanishes on the coarse levels "
        "(%.3e to %.3e)",
        before,
        after);

    if (commIsMaster(&D.comm)) {
      printf("unresolved body: residual %.3e to %.3e over 20 cycles\n", before, after);
    }
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&D.comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

#endif
