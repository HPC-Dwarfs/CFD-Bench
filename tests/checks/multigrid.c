/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Exercises the multigrid transfer operators and cycle one step at a time.
 *
 * The hierarchy is built here, straight from the discretization, rather than
 * borrowed from whichever solver happens to be linked. The cycle lives in
 * multigrid.c, which every variant links, so this driver runs under all of
 * them and needs no #ifdef seam.
 */
#include <math.h>
#include <stdlib.h>

#include "check.h"

#include "discretization.h"
#include "multigrid.h"
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
static MultigridType Mg;

/* The hierarchy the checks below drive, described from the discretization the
 * setup just built. Rebuilt per setup, and the previous one released, so the
 * geometry cases do not inherit the level arrays of the case before them. */
static int HierarchyBuilt = 0;

static void buildHierarchy(void)
{
  if (HierarchyBuilt) {
    multigridFree(&Mg);
  }

  MultigridSpecType spec = { .comm = &D.comm,
    .bc                          = &S.bc,
    .grid                        = &D.grid,
    .Ax                          = D.Ax,
    .Ay                          = D.Ay,
    .Az                          = D.Az,
    .Lambda                      = D.Lambda,
    .iOffset                     = S.iOffset,
    .jOffset                     = S.jOffset,
    .kOffset                     = S.kOffset,
    .fluidCells                  = S.fluidCells,
    .levels                      = Params.levels,
    .presmooth                   = Params.presmooth,
    .postsmooth                  = Params.postsmooth,
    .smoothOmega                 = Params.smoothOmega };

  multigridBuild(&Mg, &spec);
  HierarchyBuilt = 1;
}

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
  buildHierarchy();
}

static void setupWithGeometry(CommType *base, const char *geometry)
{
  Params.geometryFile = (char *)geometry;
  D.comm              = *base;
  commPartition(&D.comm, KMAX, JMAX, IMAX);
  initDiscretization(&D, &Params);
  initSolver(&S, &D, &Params);
  buildHierarchy();
}



static void fill(double *f, int im, int jm, int km, double value)
{
  size_t size = (size_t)(im + 2) * (jm + 2) * (km + 2);

  for (size_t i = 0; i < size; i++) {
    f[i] = value;
  }
}


/*
 * One cycle from a zero guess, as a linear operator on the right-hand side.
 * That is exactly how a preconditioner uses it, so it is what has to be
 * symmetric.
 */
static void cycleApply(const double *rhs, double *out, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    D.p[i]  = 0.0;
    D.rhs[i] = rhs[i];
  }

  mgTestVcycle(&Mg, D.p, D.rhs);

  for (size_t i = 0; i < size; i++) {
    out[i] = D.p[i];
  }
}

/* Over the fluid unknowns, which is where the operator is defined. */
static double dotFluid(const double *a, const double *b, int im, int jm, int km)
{
  const double *Lambda = D.Lambda;
  int imaxLocal        = im;
  int jmaxLocal        = jm;
  double sum           = 0.0;

  for (int k = 1; k < km + 1; k++) {
    for (int j = 1; j < jm + 1; j++) {
      for (int i = 1; i < im + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        sum += AT(a, i, j, k, im, jm) * AT(b, i, j, k, im, jm);
      }
    }
  }

  commReduceAll(&sum, SUM);
  return sum;
}

/*
 * x . M y == y . M x, which is the property the whole symmetry work exists to
 * establish: a Krylov method preconditioned by an asymmetric operator is not
 * the method it claims to be. It holds only when all three of the smoother
 * pair, the coarsest solve and the transfer pair are symmetric, so this fails
 * if any one of them is reverted.
 */
static void checkCycleSymmetry(const char *label)
{
  int im, jm, km;
  mgTestLevelExtents(&Mg, 0, &im, &jm, &km);
  size_t size = (size_t)(im + 2) * (jm + 2) * (km + 2);

  double *x  = calloc(size, sizeof(double));
  double *y  = calloc(size, sizeof(double));
  double *mx = calloc(size, sizeof(double));
  double *my = calloc(size, sizeof(double));

  int offs[NDIMS] = { 0, 0, 0 };
  commGetOffsets(&D.comm, offs, KMAX, JMAX, IMAX);

  const double *Lambda = D.Lambda;
  int imaxLocal        = im;
  int jmaxLocal        = jm;

  for (int k = 1; k < km + 1; k++) {
    for (int j = 1; j < jm + 1; j++) {
      for (int i = 1; i < im + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }

        unsigned gi = (unsigned)(i - 1 + offs[IDIM]);
        unsigned gj = (unsigned)(j - 1 + offs[JDIM]);
        unsigned gk = (unsigned)(k - 1 + offs[KDIM]);

        unsigned h  = 101u + gi * 73856093u + gj * 19349663u + gk * 83492791u;
        h ^= h >> 13; h *= 1274126177u; h ^= h >> 16;
        AT(x, i, j, k, im, jm) = (double)h / (double)0xffffffffu - 0.5;

        unsigned g = 997u + gi * 73856093u + gj * 19349663u + gk * 83492791u;
        g ^= g >> 13; g *= 1274126177u; g ^= g >> 16;
        AT(y, i, j, k, im, jm) = (double)g / (double)0xffffffffu - 0.5;
      }
    }
  }

  cycleApply(x, mx, size);
  cycleApply(y, my, size);

  double xMy   = dotFluid(x, my, im, jm, km);
  double yMx   = dotFluid(y, mx, im, jm, km);
  double scale = fabs(xMy) + fabs(yMx) + 1e-300;

  CHECK_TRUE(fabs(xMy - yMx) / scale < 1e-10,
      "%s: the cycle is not symmetric, x.My = %.17g against y.Mx = %.17g",
      label,
      xMy,
      yMx);

  if (commIsMaster(&D.comm)) {
    printf("%s: relative asymmetry %.3e\n", label, fabs(xMy - yMx) / scale);
  }

  free(x); free(y); free(mx); free(my);
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);

  CHECK_BEGIN("multigrid");

  setup(&comm, NOSLIP);

  int levels = mgTestLevels(&Mg);
  CHECK_TRUE(levels >= 2, "need at least two levels to test transfers, have %d", levels);

  if (levels >= 2) {
    int fi, fj, fk, ci, cj, ck;
    mgTestLevelExtents(&Mg, 0, &fi, &fj, &fk);
    mgTestLevelExtents(&Mg, 1, &ci, &cj, &ck);

    double fdx, fdy, fdz, cdx, cdy, cdz;
    mgTestLevelMesh(&Mg, 0, &fdx, &fdy, &fdz);
    mgTestLevelMesh(&Mg, 1, &cdx, &cdy, &cdz);

    /* Each level has its own mesh: the coarse spacing is twice the fine one.
     * Both used to be the finest spacing, so a coarse operator was wrong by a
     * factor of four per level. */
    CHECK_NEAR(cdx, 2.0 * fdx, 1e-15, "coarse dx is not twice the fine dx");
    CHECK_NEAR(cdy, 2.0 * fdy, 1e-15, "coarse dy is not twice the fine dy");
    CHECK_NEAR(cdz, 2.0 * fdz, 1e-15, "coarse dz is not twice the fine dz");

    double *fineR   = mgTestLevelR(&Mg, 0);
    /* What restriction writes and what prolongation writes: separate arrays
     * from the level's own residual and solution, which they used to share. */
    double *coarseR = mgTestLevelB(&Mg, 1);
    double *fineE   = mgTestLevelCorr(&Mg, 0);
    double *coarseE = mgTestLevelE(&Mg, 1);

    /* Restriction of a constant is that constant. */
    fill(fineR, fi, fj, fk, 3.25);
    fill(coarseR, ci, cj, ck, -999.0);
    mgTestRestrict(&Mg, 0);

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

    /*
     * Restriction is a scalar multiple of the transpose of prolongation, which
     * is what makes the cycle symmetric and therefore usable as a
     * preconditioner. Stated as the identity it actually is:
     *
     *   <R f, c>_coarse  ==  (1/8) <f, P c>_fine
     *
     * for arbitrary f and c. This replaced a check that a unit fine impulse
     * restricted into exactly one coarse cell, which was true of the eight-cell
     * average and is false of the transpose by construction -- its footprint is
     * wider because prolongation's is. The identity is the stronger statement,
     * and unlike the old one it fails if either operator is changed alone.
     *
     * Both fields are seeded from global position rather than rand(), so the
     * identity means the same thing however the domain is divided.
     */
    fill(coarseR, ci, cj, ck, 0.0);
    fill(fineE, fi, fj, fk, 0.0);

    int offsetsT[NDIMS] = { 0, 0, 0 };
    commGetOffsets(&D.comm, offsetsT, KMAX, JMAX, IMAX);

    for (int k = 1; k < fk + 1; k++) {
      for (int j = 1; j < fj + 1; j++) {
        for (int i = 1; i < fi + 1; i++) {
          unsigned gi = (unsigned)(i - 1 + offsetsT[IDIM]);
          unsigned gj = (unsigned)(j - 1 + offsetsT[JDIM]);
          unsigned gk = (unsigned)(k - 1 + offsetsT[KDIM]);
          unsigned h  = 11u + gi * 73856093u + gj * 19349663u + gk * 83492791u;
          h ^= h >> 13;
          h *= 1274126177u;
          h ^= h >> 16;
          AT(fineR, i, j, k, fi, fj) = (double)h / (double)0xffffffffu - 0.5;
        }
      }
    }

    for (int k = 1; k < ck + 1; k++) {
      for (int j = 1; j < cj + 1; j++) {
        for (int i = 1; i < ci + 1; i++) {
          unsigned gi = (unsigned)(i - 1 + offsetsT[IDIM] / 2);
          unsigned gj = (unsigned)(j - 1 + offsetsT[JDIM] / 2);
          unsigned gk = (unsigned)(k - 1 + offsetsT[KDIM] / 2);
          unsigned h  = 29u + gi * 73856093u + gj * 19349663u + gk * 83492791u;
          h ^= h >> 13;
          h *= 1274126177u;
          h ^= h >> 16;
          AT(coarseE, i, j, k, ci, cj) = (double)h / (double)0xffffffffu - 0.5;
        }
      }
    }

    mgTestRestrict(&Mg, 0);
    mgTestProlongate(&Mg, 0);

    double lhs = 0.0, rhs = 0.0;

    for (int k = 1; k < ck + 1; k++) {
      for (int j = 1; j < cj + 1; j++) {
        for (int i = 1; i < ci + 1; i++) {
          lhs += AT(coarseR, i, j, k, ci, cj) * AT(coarseE, i, j, k, ci, cj);
        }
      }
    }

    for (int k = 1; k < fk + 1; k++) {
      for (int j = 1; j < fj + 1; j++) {
        for (int i = 1; i < fi + 1; i++) {
          rhs += AT(fineR, i, j, k, fi, fj) * AT(fineE, i, j, k, fi, fj);
        }
      }
    }

    rhs *= 0.125;

    commReduceAll(&lhs, SUM);
    commReduceAll(&rhs, SUM);

    double scaleT = fabs(lhs) + fabs(rhs) + 1e-300;

    CHECK_TRUE(fabs(lhs - rhs) / scaleT < 1e-12,
        "restriction is not the transpose of prolongation: <Rf, c> = %.17g against "
        "(1/8)<f, Pc> = %.17g",
        lhs,
        rhs);

    /* Prolongation defines every fine cell: pre-fill with a sentinel and check
     * none of it survives. It used to be piecewise injection into the wrong
     * array, so half the fine cells kept whatever was there before. */
    fill(coarseE, ci, cj, ck, 2.5);
    fill(fineE, fi, fj, fk, -12345.0);
    mgTestProlongate(&Mg, 0);

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

    /*
     * Prolongation is the 3/4 - 1/4 blend it is documented to be, which a
     * constant cannot show: any two weights summing to one preserve one.
     *
     * A coarse field equal to its own x index injects to a fine field whose
     * pair members differ by exactly 1, so the blend must leave them differing
     * by exactly 1/2 -- the pair (t, t+1) becomes (t + 1/4, t + 3/4). The
     * in-place form this replaced read a neighbour it had already written and
     * produced 0.5625 here instead.
     */
    fill(fineE, fi, fj, fk, 0.0);
    for (int k = 1; k < ck + 1; k++) {
      for (int j = 1; j < cj + 1; j++) {
        for (int i = 1; i < ci + 1; i++) {
          AT(coarseE, i, j, k, ci, cj) = (double)i;
        }
      }
    }
    mgTestProlongate(&Mg, 0);

    double worstPair = 0.0;
    int pairs        = 0;

    for (int k = 2; k < fk; k++) {
      for (int j = 2; j < fj; j++) {
        /* Interior pairs only: cell 1 and the last cell blend against a halo,
         * whose value the boundary condition rather than the ramp decides. */
        for (int i = 2; i + 1 < fi; i += 2) {
          double d = AT(fineE, i + 1, j, k, fi, fj) - AT(fineE, i, j, k, fi, fj);
          double e = fabs(d - 0.5);
          if (e > worstPair) {
            worstPair = e;
          }
          ++pairs;
        }
      }
    }

    commReduceAll(&worstPair, MAX);

    CHECK_TRUE(pairs > 0, "no interior pair was available to check the blend");
    CHECK_NEAR(worstPair, 0.0, 1e-12,
        "prolongation is not the 3/4 - 1/4 blend: a pair of fine cells whose "
        "injected values differ by 1 should differ by 1/2 afterwards");

    /* Restriction then prolongation of a constant returns that constant, away
     * from the boundary where the wall condition mirrors it anyway. */
    fill(fineR, fi, fj, fk, 1.75);
    mgTestRestrict(&Mg, 0);
    fill(coarseE, ci, cj, ck, 0.0);
    for (int k = 1; k < ck + 1; k++) {
      for (int j = 1; j < cj + 1; j++) {
        for (int i = 1; i < ci + 1; i++) {
          AT(coarseE, i, j, k, ci, cj) = AT(coarseR, i, j, k, ci, cj);
        }
      }
    }
    fill(fineE, fi, fj, fk, 0.0);
    mgTestProlongate(&Mg, 0);

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
    mgTestLevelExtents(&Mg, 0, &fi, &fj, &fk);
    size_t size = (size_t)(fi + 2) * (fj + 2) * (fk + 2);

    for (size_t i = 0; i < size; i++) {
      D.p[i]   = 4.0;
      D.rhs[i] = 0.0;
    }

    mgTestVcycle(&Mg, D.p, D.rhs);

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
    mgTestLevelExtents(&Mg, 0, &fi, &fj, &fk);
    double dx, dy, dz;
    mgTestLevelMesh(&Mg, 0, &dx, &dy, &dz);

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
    mgTestSmooth(&Mg, 0, D.p, rhs, Params.presmooth + Params.postsmooth);
    double smoothed = pressureResidualNorm(&desc, D.p, rhs);

    /* A full cycle, same starting point. */
    for (size_t i = 0; i < size; i++) {
      D.p[i] = 0.0;
    }
    mgTestVcycle(&Mg, D.p, rhs);
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

    int levelsG = mgTestLevels(&Mg);

    for (int l = 0; l < levelsG; l++) {
      double solid   = (double)mgTestLevelSolidCount(&Mg, l);
      double surface = (double)mgTestLevelSurfaceCount(&Mg, l);
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

    int levelsG = mgTestLevels(&Mg);
    double coarsest = (double)mgTestLevelSolidCount(&Mg, levelsG - 1);
    commReduceAll(&coarsest, SUM);

    CHECK_TRUE(coarsest == 0.0,
        "the small sphere was still resolved at the coarsest level, so the "
        "unresolved case is untested (%.0f solid cells)",
        coarsest);

    int fi, fj, fk;
    mgTestLevelExtents(&Mg, 0, &fi, &fj, &fk);
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

          /* Zero in the body, which is what computeRHS produces: a solid cell
           * is an identity row and its right-hand side is zero. Feeding it one
           * that is not asks the smoother to drive a solid cell away from the
           * zero it is required to hold, which is a property of the input and
           * not of the cycle. */
          AT(D.rhs, i, j, k, fi, fj) =
              (AT(D.Lambda, i, j, k, fi, fj) == 0.0)
                  ? 0.0
                  : cos(M_PI * x) * cos(M_PI * y) * cos(M_PI * z);
        }
      }
    }

    for (size_t i = 0; i < sz; i++) {
      D.p[i] = 0.0;
    }

    double before = pressureResidualNorm(&desc, D.p, D.rhs);
    for (int c = 0; c < 20; c++) {
      mgTestVcycle(&Mg, D.p, D.rhs);
    }
    double after = pressureResidualNorm(&desc, D.p, D.rhs);

    CHECK_TRUE(after < before,
        "multigrid did not converge with a body that vanishes on the coarse levels "
        "(%.3e to %.3e)",
        before,
        after);

    /*
     * And the body is still at zero afterwards.
     *
     * This is the sharpest case for it: the coarse levels do not represent this
     * sphere at all, so the error prolongated from them is an ordinary nonzero
     * field right where the body is. A correction applied without consulting the
     * geometry lands in the body and stays there -- the cut-cell pass repairs the
     * solid cells on the surface list, but a cell in the interior is deliberately
     * not on that list, so nothing takes it back off.
     */
    int nonzeroSolid = 0, deepSolid = 0, deepNonzero = 0;
    const double *Lambda = D.Lambda;
    int imaxLocal = fi, jmaxLocal = fj;

    for (int k = 1; k < fk + 1; k++) {
      for (int j = 1; j < fj + 1; j++) {
        for (int i = 1; i < fi + 1; i++) {
          if (LAM(i, j, k) > 0.0) {
            continue;
          }

          if (AT(D.p, i, j, k, fi, fj) != 0.0) {
            ++nonzeroSolid;
          }

          int surrounded = LAM(i + 1, j, k) == 0.0 && LAM(i - 1, j, k) == 0.0 &&
                           LAM(i, j + 1, k) == 0.0 && LAM(i, j - 1, k) == 0.0 &&
                           LAM(i, j, k + 1) == 0.0 && LAM(i, j, k - 1) == 0.0;

          if (surrounded) {
            ++deepSolid;
            if (AT(D.p, i, j, k, fi, fj) != 0.0) {
              ++deepNonzero;
            }
          }
        }
      }
    }

    double totals[3] = { (double)nonzeroSolid, (double)deepSolid, (double)deepNonzero };
    for (int t = 0; t < 3; t++) {
      commReduceAll(&totals[t], SUM);
    }

    CHECK_TRUE(totals[1] > 0.0,
        "the body has no solid cell away from its surface, so the cells the surface "
        "list does not cover are untested");
    CHECK_TRUE(totals[0] == 0.0,
        "%.0f solid cells hold a nonzero pressure after cycles whose coarse levels "
        "do not represent the body",
        totals[0]);
    CHECK_TRUE(totals[2] == 0.0,
        "%.0f solid cells away from the surface were moved off zero by the coarse "
        "correction",
        totals[2]);

    if (commIsMaster(&D.comm)) {
      printf("unresolved body: residual %.3e to %.3e over 20 cycles, %.0f interior "
             "solid cells, %.0f of them nonzero\n",
          before,
          after,
          totals[1],
          totals[2]);
    }
  }

  /* ---- The cycle is a symmetric operator ---- */
  setup(&comm, NOSLIP);
  checkCycleSymmetry("cycle symmetry, obstacle-free");

  setupWithGeometry(&comm, "sphere:1.0,0.5,0.5,0.3");
  checkCycleSymmetry("cycle symmetry, sphere");

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&D.comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

