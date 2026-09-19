/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The multigrid hierarchy and the cycle, see multigrid.h for why they live
 * here rather than in solver-mg.c.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "allocate.h"
#include "coloring.h"
#include "multigrid.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "surface-list.h"
#include "timing.h"
#include "util.h"

#define FINEST_LEVEL 0

/* util.h supplies P and RHS; the level's own error and residual need the same
 * treatment, and like those they read imaxLocal and jmaxLocal from scope. */
#define E(i, j, k)                                                                       \
  e[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define R(i, j, k)                                                                       \
  r[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

static void levelDescribe(MultigridType *mg, MgLevelType *lv, PressureLevelType *out)
{
  out->comm       = &lv->comm;
  out->bc         = mg->bc;
  out->Ax         = lv->Ax;
  out->Ay         = lv->Ay;
  out->Az         = lv->Az;
  out->Lambda     = lv->Lambda;
  out->list       = &lv->surface;
  out->imaxLocal  = lv->imaxLocal;
  out->jmaxLocal  = lv->jmaxLocal;
  out->kmaxLocal  = lv->kmaxLocal;
  out->dx         = lv->dx;
  out->dy         = lv->dy;
  out->dz         = lv->dz;
  out->fluidCells = lv->fluidCells;
}

static size_t levelSize(MgLevelType *lv)
{
  return (size_t)(lv->imaxLocal + 2) * (lv->jmaxLocal + 2) * (lv->kmaxLocal + 2);
}

static void zeroField(MgLevelType *lv, double *f)
{
  size_t size = levelSize(lv);

  for (size_t i = 0; i < size; i++) {
    f[i] = 0.0;
  }
}

static void prolongateBlend(double *f, int im, int jm, int km, int fi, int fj, int axis);

/*
 * Restriction: a scalar multiple of the transpose of prolongation.
 *
 * That is not a nicety. A V-cycle is symmetric only when the post-smoother is
 * the transpose of the pre-smoother, the coarsest solve is symmetric, and
 * R = c P^T; with prolongation trilinear, an eight-cell average is not that
 * multiple and the cycle cannot precondition a Krylov method.
 *
 * Prolongation is injection followed by three one-dimensional blends,
 *
 *   P    =  blend_z . blend_y . blend_x . inject
 *   P^T  =  inject^T . blend_x . blend_y . blend_z
 *
 * using that each blend is symmetric -- it updates a pair of cells from their
 * own two values with the same weight either way -- so each is its own
 * transpose, and only their order reverses. inject^T is the sum over a coarse
 * cell's eight children, which is what the averaging form already computed
 * before scaling.
 *
 * The scale is 1/8, chosen so that restriction still maps a constant to that
 * constant: the column sums of P are 8, so (1/8) P^T returns a constant
 * unchanged. That keeps the existing constant-preservation checks meaningful
 * rather than rewritten around a new normalization.
 *
 * Each blend reads only face neighbours, so restriction inherits prolongation's
 * property of needing no communication the halo exchange does not already
 * provide.
 */
static void restrictMG(MultigridType *mg, MgLevelType *fine, MgLevelType *coarse)
{
  double *fineField   = fine->r;
  double *work        = fine->scratch;
  double *coarseField = coarse->b;

  int fi              = fine->imaxLocal + 2;
  int fj              = fine->jmaxLocal + 2;
  int ci              = coarse->imaxLocal + 2;
  int cj              = coarse->jmaxLocal + 2;

  int im              = fine->imaxLocal;
  int jm              = fine->jmaxLocal;
  int km              = fine->kmaxLocal;

  size_t size         = levelSize(fine);

  /* The blends run on a copy: the residual they are given is still the level's
   * own and must survive. */
  for (size_t i = 0; i < size; i++) {
    work[i] = fineField[i];
  }

  /* The three blends, in the reverse of the order prolongation applies them.
   * Each touches one axis only, so in fact they commute and the order is free;
   * it is written reversed because that is what the transpose is, and a reader
   * checking the derivation should not have to notice the commutation first. */
  commExchange(&fine->comm, work);
  pressureBcApply(mg->bc, &fine->comm, work, im, jm, km);
  prolongateBlend(work, im, jm, km, fi, fj, 2);

  commExchange(&fine->comm, work);
  pressureBcApply(mg->bc, &fine->comm, work, im, jm, km);
  prolongateBlend(work, im, jm, km, fi, fj, 1);

  commExchange(&fine->comm, work);
  pressureBcApply(mg->bc, &fine->comm, work, im, jm, km);
  prolongateBlend(work, im, jm, km, fi, fj, 0);

  /* inject^T, scaled: the sum over each coarse cell's eight children. */
  commExchange(&fine->comm, work);

  for (int k = 1; k < coarse->kmaxLocal + 1; k++) {
    for (int j = 1; j < coarse->jmaxLocal + 1; j++) {
      for (int i = 1; i < coarse->imaxLocal + 1; i++) {
        double sum = 0.0;

        for (int dk = 0; dk < 2; dk++) {
          for (int dj = 0; dj < 2; dj++) {
            for (int di = 0; di < 2; di++) {
              size_t idx = (size_t)(2 * k - 1 + dk) * fi * fj +
                           (size_t)(2 * j - 1 + dj) * fi + (size_t)(2 * i - 1 + di);
              sum += work[idx];
            }
          }
        }

        coarseField[(size_t)k * ci * cj + (size_t)j * ci + (size_t)i] = sum * 0.125;
      }
    }
  }
}


/*
 * One of prolongation's three one-dimensional blends, along the axis whose
 * neighbour is `stride` away.
 *
 * `nb` is an involution: it pairs interior cells (2,3), (4,5) and so on, and
 * leaves cell 1 -- and the last cell when the extent is even -- blending against
 * a halo, which no pass writes. Updating a pair from its own two values is what
 * makes this the 3/4 - 1/4 stencil it is documented to be.
 *
 * Done that way rather than in place because the in-place form read a neighbour
 * it had already overwritten: for odd i the neighbour is i-1, written one
 * iteration earlier, which gave every other fine cell effective weights of
 * 0.8125 and 0.1875. A constant is invariant under either, which is why the
 * transfer checks did not see it. It also matters structurally -- the paired
 * form is symmetric, and restriction is its transpose.
 */
static void prolongateBlend(double *f, int im, int jm, int km, int fi, int fj, int axis)
{
  size_t stride = (axis == 0) ? 1 : (axis == 1) ? (size_t)fi : (size_t)fi * fj;

  for (int k = 1; k < km + 1; k++) {
    for (int j = 1; j < jm + 1; j++) {
      for (int i = 1; i < im + 1; i++) {
        /* Walk the axis in the layout's own index so the pairing below is the
         * same arithmetic the original nb computed. */
        int n         = (axis == 0) ? i : (axis == 1) ? j : k;
        size_t idx    = (size_t)k * fi * fj + (size_t)j * fi + (size_t)i;

        if (n & 1) {
          /* Lower member of a pair, or cell 1 against the halo. Cell 1 is the
           * only odd cell whose partner the loop does not also visit. */
          if (n == 1) {
            f[idx] = 0.75 * f[idx] + 0.25 * f[idx - stride];
          }
          continue;
        }

        /* Even cell: it and the odd cell above it are a pair, unless that one
         * is the halo. */
        int limit = (axis == 0) ? im : (axis == 1) ? jm : km;

        if (n == limit) {
          f[idx] = 0.75 * f[idx] + 0.25 * f[idx + stride];
          continue;
        }

        double t0     = f[idx];
        double t1     = f[idx + stride];
        f[idx]        = 0.75 * t0 + 0.25 * t1;
        f[idx + stride] = 0.75 * t1 + 0.25 * t0;
      }
    }
  }
}

/*
 * Trilinear prolongation for cell-centred data, writing every fine cell.
 *
 * A fine cell centre sits either a quarter or three quarters of the way between
 * two coarse centres in each direction, so the one-dimensional weights are
 * always 3/4 and 1/4; which coarse neighbour carries the 1/4 depends on which
 * of the two children the fine cell is.
 *
 * Done as injection followed by three one-dimensional passes rather than as one
 * tensor-product stencil. The two are the same operator, because it is
 * separable, but the tensor-product form reads coarse cells that are diagonal
 * neighbours, and the halo exchange transfers faces only -- so at a subdomain
 * edge that form would read whatever happened to be in memory, and the result
 * would depend on the decomposition. Each pass here reads only face
 * neighbours, which is exactly what the exchange provides.
 *
 * Which child a fine cell is follows from its local index alone: the fine
 * offset is twice the coarse one, so it is even, and local index 1 is always a
 * lower child.
 *
 * What was here before was piecewise injection, and it wrote into the residual
 * array while the correction read the error array, so the coarse correction
 * never reached the solution at all.
 */
static void prolongate(MultigridType *mg, MgLevelType *coarse, MgLevelType *fine)
{
  double *coarseField = coarse->e;
  double *fineField   = fine->corr;

  int ci              = coarse->imaxLocal + 2;
  int cj              = coarse->jmaxLocal + 2;
  int fi              = fine->imaxLocal + 2;
  int fj              = fine->jmaxLocal + 2;

  int im              = fine->imaxLocal;
  int jm              = fine->jmaxLocal;
  int km              = fine->kmaxLocal;

  /* Inject: every fine cell takes the value of the coarse cell containing it.
   * Reads no halo. */
  for (int k = 1; k < km + 1; k++) {
    for (int j = 1; j < jm + 1; j++) {
      for (int i = 1; i < im + 1; i++) {
        fineField[(size_t)k * fi * fj + (size_t)j * fi + (size_t)i] =
            coarseField[(size_t)((k + 1) / 2) * ci * cj + (size_t)((j + 1) / 2) * ci +
                        (size_t)((i + 1) / 2)];
      }
    }
  }

  /* Three one-dimensional corrections. After injection the neighbour across the
   * coarse cell boundary holds the neighbouring coarse value, so the 3/4 - 1/4
   * blend is a plain face-neighbour stencil. */
  commExchange(&fine->comm, fineField);
  pressureBcApply(mg->bc, &fine->comm, fineField, im, jm, km);

  prolongateBlend(fineField, im, jm, km, fi, fj, 0);

  commExchange(&fine->comm, fineField);
  pressureBcApply(mg->bc, &fine->comm, fineField, im, jm, km);

  prolongateBlend(fineField, im, jm, km, fi, fj, 1);

  commExchange(&fine->comm, fineField);
  pressureBcApply(mg->bc, &fine->comm, fineField, im, jm, km);

  prolongateBlend(fineField, im, jm, km, fi, fj, 2);
}

/*
 * p += e on this level, over the fluid only.
 *
 * A coarse level represents the body only approximately -- a feature thinner
 * than a coarse cell rounds away entirely -- so the error prolongated from it is
 * defined over cells the body occupies and is not zero there. Adding it would
 * move a solid cell off the zero its identity row requires.
 *
 * Nothing later in the cycle would take it back off. The post-smoother's
 * cut-cell pass repairs the solid cells on the surface list, which is why this
 * showed up only in the body's interior: those cells are deliberately absent
 * from that list, because relaxing a zero neighbourhood with a zero right-hand
 * side keeps them at zero without any help. That argument assumes nothing writes
 * into them from outside, and this is the one place that did.
 */
static void correct(MgLevelType *lv, double *p)
{
  int imaxLocal        = lv->imaxLocal;
  int jmaxLocal        = lv->jmaxLocal;
  int kmaxLocal        = lv->kmaxLocal;
  double *e            = lv->corr;
  const double *Lambda = lv->Lambda;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        P(i, j, k) += E(i, j, k);
      }
    }
  }
}

/*
 * Red-black SOR on one level, with that level's mesh and the global
 * checkerboard.
 *
 * reverse visits the two colours black-first instead of red-first. As an error
 * propagation operator, a red-then-black sweep has a black-then-red sweep as its
 * transpose, so running the post-smoother reversed is what makes the pair of
 * smoothing phases around the coarse correction a transpose pair -- and that,
 * with a symmetric coarsest solve and transfer operators that are a transpose
 * pair, is what makes the whole cycle symmetric.
 */
static void smooth(MultigridType *mg,
    MgLevelType *lv,
    double *p,
    const double *rhs,
    int sweeps,
    int reverse)
{
  int imaxLocal = lv->imaxLocal;
  int jmaxLocal = lv->jmaxLocal;
  int kmaxLocal = lv->kmaxLocal;

  double dx2    = lv->dx * lv->dx;
  double dy2    = lv->dy * lv->dy;
  double dz2    = lv->dz * lv->dz;
  double idx2   = 1.0 / dx2;
  double idy2   = 1.0 / dy2;
  double idz2   = 1.0 / dz2;
  double factor =
      mg->smoothOmega * 0.5 * (dx2 * dy2 * dz2) / (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);

  PressureLevelType desc;
  levelDescribe(mg, lv, &desc);
  double ignored = 0.0;

  for (int sweep = 0; sweep < sweeps; sweep++) {
    for (int c = 0; c < 2; c++) {
      int color = reverse ? 1 - c : c;

      commExchange(&lv->comm, p);

      /* Same split as the relaxation solvers: a geometry-free interior sweep,
       * then a correction over this level's own surface list. */
      pressureSaveSurface(&desc, p, color);

#ifdef PROFILING
      double bulkStart = getTimeStamp();
#endif

      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          int iStart = colorRowStart(color, j, k, lv->iOffset, lv->jOffset, lv->kOffset);

          for (int i = iStart; i < imaxLocal + 1; i += 2) {
            P(i, j, k) -=
                factor *
                (RHS(i, j, k) -
                    ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                        (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
                        (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2));
          }
        }
      }

#ifdef PROFILING
      T[SWEEP_BULK] += getTimeStamp() - bulkStart;
      C[SWEEP_BULK]++;
#endif

      pressureCorrectSurface(&desc, p, rhs, color, mg->smoothOmega, &ignored);
      pressureBcApply(mg->bc, &lv->comm, p, imaxLocal, jmaxLocal, kmaxLocal);
    }
  }
}

/* r = rhs - A p on this level, written into lv->r. */
static void residualField(MultigridType *mg, MgLevelType *lv, double *p, const double *rhs)
{
  int imaxLocal = lv->imaxLocal;
  int jmaxLocal = lv->jmaxLocal;
  int kmaxLocal = lv->kmaxLocal;

  double idx2   = 1.0 / (lv->dx * lv->dx);
  double idy2   = 1.0 / (lv->dy * lv->dy);
  double idz2   = 1.0 / (lv->dz * lv->dz);
  double *r     = lv->r;

  commExchange(&lv->comm, p);
  pressureBcApply(mg->bc, &lv->comm, p, imaxLocal, jmaxLocal, kmaxLocal);

  const double *Ax     = lv->Ax;
  const double *Ay     = lv->Ay;
  const double *Az     = lv->Az;
  const double *Lambda = lv->Lambda;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        /* A solid cell is an identity row already satisfied by p = 0, so it has
         * no residual to pass down the hierarchy. */
        if (LAM(i, j, k) == 0.0) {
          R(i, j, k) = 0.0;
          continue;
        }

        double pc  = P(i, j, k);
        R(i, j, k) = LAM(i, j, k) * RHS(i, j, k) -
                     (AX(i, j, k) * (P(i + 1, j, k) - pc) * idx2 +
                         AX(i - 1, j, k) * (P(i - 1, j, k) - pc) * idx2 +
                         AY(i, j, k) * (P(i, j + 1, k) - pc) * idy2 +
                         AY(i, j - 1, k) * (P(i, j - 1, k) - pc) * idy2 +
                         AZ(i, j, k) * (P(i, j, k + 1) - pc) * idz2 +
                         AZ(i, j, k - 1) * (P(i, j, k - 1) - pc) * idz2);
      }
    }
  }
}

static void vcycle(MultigridType *mg, int level, double *p, const double *rhs)
{
  MgLevelType *levels = mg->level;
  MgLevelType *lv     = &levels[level];

  if (level == mg->levels - 1) {
    /*
     * Coarsest level: relax hard enough to stand in for a solve, and
     * symmetrically, since the cycle is only symmetric if every part of it is.
     * Equal numbers of forward and backward sweeps make the error propagation
     * a product of the form A^T A, which is symmetric for any sweep count --
     * the same argument that makes the post-smoother the transpose of the
     * pre-smoother, applied to a level with nothing below it.
     */
    smooth(mg, lv, p, rhs, mg->presmooth, 0);
    smooth(mg, lv, p, rhs, mg->postsmooth, 1);
    return;
  }

  MgLevelType *coarse = &levels[level + 1];

  smooth(mg, lv, p, rhs, mg->presmooth, 0);
  residualField(mg, lv, p, rhs);
  restrictMG(mg, lv, coarse);

  /* The coarse problem is for the error, so it starts from zero. Leaving the
   * previous cycle's error in place made every cycle after the first solve a
   * different problem from the one it was given. */
  zeroField(coarse, coarse->e);
  vcycle(mg, level + 1, coarse->e, coarse->b);

  prolongate(mg, coarse, lv);
  correct(lv, p);
  pressureBcApply(mg->bc, &lv->comm, p, lv->imaxLocal, lv->jmaxLocal, lv->kmaxLocal);

  smooth(mg, lv, p, rhs, mg->postsmooth, 1);
}

/*
 * Coarsen the geometry by one level.
 *
 * A coarse face aperture is the mean of the four fine faces it covers; a coarse
 * volume fraction is the mean of the eight fine cells. Cheap, local, and
 * consistent with a flux-balance assembly, which is what makes it the right
 * first choice over a Galerkin coarse operator -- that would be 27-point in
 * three dimensions and would change the coarse kernels and their layout.
 *
 * The means are then rounded, because this phase carries binary apertures. A
 * feature thinner than a coarse cell rounds away, which is reported rather than
 * left to be discovered as a convergence problem.
 */
static void coarsenGeometry(MgLevelType *fine, MgLevelType *coarse)
{
  int fi = fine->imaxLocal + 2;
  int fj = fine->jmaxLocal + 2;
  int ci = coarse->imaxLocal + 2;
  int cj = coarse->jmaxLocal + 2;

#define FINE(f, i, j, k) (f)[(size_t)(k) * fi * fj + (size_t)(j) * fi + (size_t)(i)]
#define CRS(f, i, j, k) (f)[(size_t)(k) * ci * cj + (size_t)(j) * ci + (size_t)(i)]

  size_t size = (size_t)ci * cj * (coarse->kmaxLocal + 2);

  for (size_t i = 0; i < size; i++) {
    coarse->Ax[i] = coarse->Ay[i] = coarse->Az[i] = coarse->Lambda[i] = 1.0;
  }

  for (int k = 1; k < coarse->kmaxLocal + 1; k++) {
    for (int j = 1; j < coarse->jmaxLocal + 1; j++) {
      for (int i = 1; i < coarse->imaxLocal + 1; i++) {
        int fiI = 2 * i - 1, fjI = 2 * j - 1, fkI = 2 * k - 1;

        /* Eight fine cells for the volume fraction. */
        double vol = 0.0;
        for (int dk = 0; dk < 2; dk++) {
          for (int dj = 0; dj < 2; dj++) {
            for (int di = 0; di < 2; di++) {
              vol += FINE(fine->Lambda, fiI + di, fjI + dj, fkI + dk);
            }
          }
        }
        CRS(coarse->Lambda, i, j, k) = (vol >= 4.0) ? 1.0 : 0.0;

        /* Four fine faces for each aperture. The coarse x-face of cell i is the
         * fine x-face at 2i, which is the far side of the upper fine child. */
        double ax = 0.0, ay = 0.0, az = 0.0;
        for (int b = 0; b < 2; b++) {
          for (int a = 0; a < 2; a++) {
            ax += FINE(fine->Ax, 2 * i, fjI + a, fkI + b);
            ay += FINE(fine->Ay, fiI + a, 2 * j, fkI + b);
            az += FINE(fine->Az, fiI + a, fjI + b, 2 * k);
          }
        }
        CRS(coarse->Ax, i, j, k) = (ax >= 2.0) ? 1.0 : 0.0;
        CRS(coarse->Ay, i, j, k) = (ay >= 2.0) ? 1.0 : 0.0;
        CRS(coarse->Az, i, j, k) = (az >= 2.0) ? 1.0 : 0.0;
      }
    }
  }

  /* The same rule the producer enforces: a solid cell has no open face. */
  for (int k = 0; k < coarse->kmaxLocal + 2; k++) {
    for (int j = 0; j < coarse->jmaxLocal + 2; j++) {
      for (int i = 0; i < coarse->imaxLocal + 2; i++) {
        if (CRS(coarse->Lambda, i, j, k) > 0.0) {
          continue;
        }
        CRS(coarse->Ax, i, j, k) = 0.0;
        CRS(coarse->Ay, i, j, k) = 0.0;
        CRS(coarse->Az, i, j, k) = 0.0;
        if (i > 0) {
          CRS(coarse->Ax, i - 1, j, k) = 0.0;
        }
        if (j > 0) {
          CRS(coarse->Ay, i, j - 1, k) = 0.0;
        }
        if (k > 0) {
          CRS(coarse->Az, i, j, k - 1) = 0.0;
        }
      }
    }
  }

  commExchange(&coarse->comm, coarse->Ax);
  commExchange(&coarse->comm, coarse->Ay);
  commExchange(&coarse->comm, coarse->Az);
  commExchange(&coarse->comm, coarse->Lambda);

#undef FINE
#undef CRS
}

/* Global count of cells this level calls fluid. */
static double countFluid(MgLevelType *lv)
{
  int ci       = lv->imaxLocal + 2;
  int cj       = lv->jmaxLocal + 2;
  double fluid = 0.0;

  for (int k = 1; k < lv->kmaxLocal + 1; k++) {
    for (int j = 1; j < lv->jmaxLocal + 1; j++) {
      for (int i = 1; i < lv->imaxLocal + 1; i++) {
        if (lv->Lambda[(size_t)k * ci * cj + (size_t)j * ci + (size_t)i] > 0.0) {
          fluid += 1.0;
        }
      }
    }
  }

  commReduceAll(&fluid, SUM);
  return (fluid > 0.0) ? fluid : 1.0;
}

void multigridBuild(MultigridType *mg, const MultigridSpecType *spec)
{
  /*
   * The cycle is symmetric only when the two smoothing phases are a transpose
   * pair, and a product of four operators is not the transpose of a product of
   * five. Refused here rather than quietly reconciled: a cycle that is not the
   * one the setup asked for would still converge, and nothing downstream would
   * say so -- but a Krylov method preconditioned by it would no longer be
   * conjugate gradients.
   */
  if (spec->presmooth != spec->postsmooth) {
    if (commIsMaster(spec->comm)) {
      fprintf(stderr,
          "Multigrid: presmooth is %d and postsmooth is %d. They must be equal, "
          "because the post-smoother is the transpose of the pre-smoother and a "
          "transpose runs the same number of sweeps.\n",
          spec->presmooth,
          spec->postsmooth);
    }
    exit(EXIT_FAILURE);
  }

  mg->bc         = spec->bc;
  mg->smoothOmega = spec->smoothOmega;
  mg->levels     = spec->levels;
  mg->presmooth  = spec->presmooth;
  mg->postsmooth = spec->postsmooth;

  MgLevelType *levels = malloc((size_t)mg->levels * sizeof(MgLevelType));

  /* Finest level: the decomposition the solver already has. */
  levels[0].comm         = *spec->comm;
  levels[0].imaxLocal    = spec->comm->imaxLocal;
  levels[0].jmaxLocal    = spec->comm->jmaxLocal;
  levels[0].kmaxLocal    = spec->comm->kmaxLocal;
  levels[0].iOffset      = spec->iOffset;
  levels[0].jOffset      = spec->jOffset;
  levels[0].kOffset      = spec->kOffset;
  levels[0].dx           = spec->grid->dx;
  levels[0].dy           = spec->grid->dy;
  levels[0].dz           = spec->grid->dz;
  levels[0].cells        = (double)spec->grid->imax * spec->grid->jmax * spec->grid->kmax;
  levels[0].Ax           = (double *)spec->Ax;
  levels[0].Ay           = (double *)spec->Ay;
  levels[0].Az           = (double *)spec->Az;
  levels[0].Lambda       = (double *)spec->Lambda;
  levels[0].ownsGeometry = 0;
  levels[0].fluidCells   = spec->fluidCells;

  int built              = 1;

  for (int l = 1; l < mg->levels; l++) {
    MgLevelType *fine = &levels[l - 1];

    /* Coarsening halves each extent, and halving the global offset is only the
     * same partition when every local extent is even. Ranks must agree, so the
     * test is reduced across them. */
    double ok = (fine->imaxLocal % 2 == 0 && fine->jmaxLocal % 2 == 0 &&
                    fine->kmaxLocal % 2 == 0 && fine->imaxLocal / 2 >= 2 &&
                    fine->jmaxLocal / 2 >= 2 && fine->kmaxLocal / 2 >= 2)
                    ? 1.0
                    : 0.0;
    commReduceAll(&ok, MIN);

    if (ok < 0.5) {
      break;
    }

    MgLevelType *lv = &levels[l];
    commUpdateDatatypes(
        &fine->comm, &lv->comm, fine->imaxLocal, fine->jmaxLocal, fine->kmaxLocal);

    lv->imaxLocal = fine->imaxLocal / 2;
    lv->jmaxLocal = fine->jmaxLocal / 2;
    lv->kmaxLocal = fine->kmaxLocal / 2;
    lv->iOffset   = fine->iOffset / 2;
    lv->jOffset   = fine->jOffset / 2;
    lv->kOffset   = fine->kOffset / 2;
    lv->dx        = fine->dx * 2.0;
    lv->dy        = fine->dy * 2.0;
    lv->dz        = fine->dz * 2.0;
    lv->cells     = fine->cells / 8.0;
    ++built;
  }

  if (built < mg->levels) {
    if (commIsMaster(spec->comm)) {
      printf("Multigrid: requested %d levels, the decomposition supports %d\n",
          mg->levels,
          built);
    }
    mg->levels = built;
  }

  for (int l = 0; l < mg->levels; l++) {
    size_t size = levelSize(&levels[l]);
    levels[l].e       = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
    levels[l].r       = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
    levels[l].b       = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
    levels[l].corr    = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
    levels[l].scratch = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
    zeroField(&levels[l], levels[l].e);
    zeroField(&levels[l], levels[l].r);
    zeroField(&levels[l], levels[l].b);
    zeroField(&levels[l], levels[l].corr);
    zeroField(&levels[l], levels[l].scratch);

    if (l > 0) {
      levels[l].Ax           = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
      levels[l].Ay           = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
      levels[l].Az           = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
      levels[l].Lambda       = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
      levels[l].ownsGeometry = 1;
      coarsenGeometry(&levels[l - 1], &levels[l]);
      levels[l].fluidCells = countFluid(&levels[l]);
    }

    surfaceListBuild(&levels[l].surface,
        levels[l].imaxLocal,
        levels[l].jmaxLocal,
        levels[l].kmaxLocal,
        levels[l].iOffset,
        levels[l].jOffset,
        levels[l].kOffset,
        levels[l].Ax,
        levels[l].Ay,
        levels[l].Az,
        levels[l].Lambda,
        1.0 / (levels[l].dx * levels[l].dx),
        1.0 / (levels[l].dy * levels[l].dy),
        1.0 / (levels[l].dz * levels[l].dz));
  }

  /* A body that stops being represented on a coarse level is reported rather
   * than left to show up as a convergence problem. The cycle still runs; the
   * coarse correction is simply blind to the feature. */
  if (commIsMaster(spec->comm)) {
    double finestSolid = levels[0].cells - levels[0].fluidCells;

    for (int l = 1; l < mg->levels; l++) {
      double solid = levels[l].cells - levels[l].fluidCells;

      if (finestSolid > 0.0 && solid == 0.0) {
        printf("Multigrid: the obstacle is unresolved at level %d of %d and is not "
               "represented there\n",
            l,
            mg->levels);
        break;
      }
    }
  }


  mg->level = levels;
}

void multigridFree(MultigridType *mg)
{
  if (mg->level == NULL) {
    return;
  }

  for (int l = 0; l < mg->levels; l++) {
    surfaceListFree(&mg->level[l].surface);

    if (mg->level[l].ownsGeometry) {
      free(mg->level[l].Ax);
      free(mg->level[l].Ay);
      free(mg->level[l].Az);
      free(mg->level[l].Lambda);
    }

    free(mg->level[l].e);
    free(mg->level[l].r);
    free(mg->level[l].b);
    free(mg->level[l].corr);
    free(mg->level[l].scratch);
  }

  free(mg->level);
  mg->level  = NULL;
  mg->levels = 0;
}

void multigridCycle(MultigridType *mg, double *p, const double *rhs)
{
  vcycle(mg, FINEST_LEVEL, p, rhs);
}

SurfaceListType *multigridFinestSurface(MultigridType *mg)
{
  return &mg->level[FINEST_LEVEL].surface;
}

#if defined(TEST)
int mgTestLevels(MultigridType *mg) { return mg->levels; }

void mgTestLevelExtents(MultigridType *mg, int level, int *im, int *jm, int *km)
{
  MgLevelType *lv = &mg->level[level];
  *im             = lv->imaxLocal;
  *jm             = lv->jmaxLocal;
  *km             = lv->kmaxLocal;
}

void mgTestLevelMesh(MultigridType *mg, int level, double *dx, double *dy, double *dz)
{
  MgLevelType *lv = &mg->level[level];
  *dx             = lv->dx;
  *dy             = lv->dy;
  *dz             = lv->dz;
}

double *mgTestLevelE(MultigridType *mg, int level) { return mg->level[level].e; }

double *mgTestLevelR(MultigridType *mg, int level) { return mg->level[level].r; }

double *mgTestLevelB(MultigridType *mg, int level) { return mg->level[level].b; }

double *mgTestLevelCorr(MultigridType *mg, int level) { return mg->level[level].corr; }

void mgTestRestrict(MultigridType *mg, int level)
{
  restrictMG(mg, &mg->level[level], &mg->level[level + 1]);
}

void mgTestProlongate(MultigridType *mg, int level)
{
  prolongate(mg, &mg->level[level + 1], &mg->level[level]);
}

void mgTestResidualField(MultigridType *mg, int level, double *p, const double *rhs)
{
  residualField(mg, &mg->level[level], p, rhs);
}

void mgTestVcycle(MultigridType *mg, double *p, const double *rhs)
{
  vcycle(mg, FINEST_LEVEL, p, rhs);
}

void mgTestSmooth(MultigridType *mg, int level, double *p, const double *rhs, int sweeps)
{
  smooth(mg, &mg->level[level], p, rhs, sweeps, 0);
}

int mgTestLevelSolidCount(MultigridType *mg, int level)
{
  MgLevelType *lv      = &mg->level[level];
  int imaxLocal        = lv->imaxLocal;
  int jmaxLocal        = lv->jmaxLocal;
  const double *Lambda = lv->Lambda;
  int solid            = 0;

  for (int k = 1; k < lv->kmaxLocal + 1; k++) {
    for (int j = 1; j < lv->jmaxLocal + 1; j++) {
      for (int i = 1; i < lv->imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          ++solid;
        }
      }
    }
  }

  return solid;
}

int mgTestLevelSurfaceCount(MultigridType *mg, int level)
{
  return mg->level[level].surface.count;
}
#endif
