/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "allocate.h"
#include "coloring.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "timing.h"
#include "util.h"

#define FINEST_LEVEL 0

/*
 * One level of the hierarchy. Each carries its own extents, mesh size, global
 * offsets and communicator, because the operator, the colouring and the
 * boundary condition all depend on them and all three used to be taken from
 * the finest level regardless of which level was being relaxed.
 *
 * e and r are this level's error and residual. The finest level relaxes the
 * caller's pressure against the caller's right-hand side and uses only r; every
 * coarser level relaxes its own e against its own r.
 */
typedef struct {
  CommType comm;
  int imaxLocal, jmaxLocal, kmaxLocal;
  int iOffset, jOffset, kOffset;
  double dx, dy, dz;
  double cells; /* global interior cells at this level */
  double *e, *r;
} MgLevelType;

/* util.h supplies P and RHS; the level's own error and residual need the same
 * treatment, and like those they read imaxLocal and jmaxLocal from scope. */
#define E(i, j, k)                                                                       \
  e[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define R(i, j, k)                                                                       \
  r[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

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

/*
 * Restriction for cell-centred data: a coarse cell is the mean of the eight
 * fine cells it covers.
 *
 * What was here before applied a 27-point vertex-centred full-weighting stencil
 * to cell-centred data, and indexed two of its three loops with the bounds of
 * the wrong axis.
 */
static void restrictMG(MgLevelType *fine, MgLevelType *coarse)
{
  double *fineField   = fine->r;
  double *coarseField = coarse->r;

  int fi = fine->imaxLocal + 2;
  int fj = fine->jmaxLocal + 2;
  int ci = coarse->imaxLocal + 2;
  int cj = coarse->jmaxLocal + 2;

  commExchange(&fine->comm, fineField);

  for (int k = 1; k < coarse->kmaxLocal + 1; k++) {
    for (int j = 1; j < coarse->jmaxLocal + 1; j++) {
      for (int i = 1; i < coarse->imaxLocal + 1; i++) {
        double sum = 0.0;

        for (int dk = 0; dk < 2; dk++) {
          for (int dj = 0; dj < 2; dj++) {
            for (int di = 0; di < 2; di++) {
              size_t idx = (size_t)(2 * k - 1 + dk) * fi * fj +
                           (size_t)(2 * j - 1 + dj) * fi + (size_t)(2 * i - 1 + di);
              sum += fineField[idx];
            }
          }
        }

        coarseField[(size_t)k * ci * cj + (size_t)j * ci + (size_t)i] = sum * 0.125;
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
static void prolongate(Solver *s, MgLevelType *coarse, MgLevelType *fine)
{
  double *coarseField = coarse->e;
  double *fineField   = fine->e;

  int ci = coarse->imaxLocal + 2;
  int cj = coarse->jmaxLocal + 2;
  int fi = fine->imaxLocal + 2;
  int fj = fine->jmaxLocal + 2;

  int im = fine->imaxLocal;
  int jm = fine->jmaxLocal;
  int km = fine->kmaxLocal;

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
  pressureBcApply(&s->bc, &fine->comm, fineField, im, jm, km);

  for (int k = 1; k < km + 1; k++) {
    for (int j = 1; j < jm + 1; j++) {
      for (int i = 1; i < im + 1; i++) {
        size_t idx = (size_t)k * fi * fj + (size_t)j * fi + (size_t)i;
        size_t nb  = (i & 1) ? idx - 1 : idx + 1;
        fineField[idx] = 0.75 * fineField[idx] + 0.25 * fineField[nb];
      }
    }
  }

  commExchange(&fine->comm, fineField);
  pressureBcApply(&s->bc, &fine->comm, fineField, im, jm, km);

  for (int k = 1; k < km + 1; k++) {
    for (int j = 1; j < jm + 1; j++) {
      for (int i = 1; i < im + 1; i++) {
        size_t idx = (size_t)k * fi * fj + (size_t)j * fi + (size_t)i;
        size_t nb  = (j & 1) ? idx - fi : idx + fi;
        fineField[idx] = 0.75 * fineField[idx] + 0.25 * fineField[nb];
      }
    }
  }

  commExchange(&fine->comm, fineField);
  pressureBcApply(&s->bc, &fine->comm, fineField, im, jm, km);

  for (int k = 1; k < km + 1; k++) {
    for (int j = 1; j < jm + 1; j++) {
      for (int i = 1; i < im + 1; i++) {
        size_t idx     = (size_t)k * fi * fj + (size_t)j * fi + (size_t)i;
        size_t nb      = (k & 1) ? idx - (size_t)fi * fj : idx + (size_t)fi * fj;
        fineField[idx] = 0.75 * fineField[idx] + 0.25 * fineField[nb];
      }
    }
  }
}

/* p += e on this level. */
static void correct(MgLevelType *lv, double *p)
{
  int imaxLocal = lv->imaxLocal;
  int jmaxLocal = lv->jmaxLocal;
  int kmaxLocal = lv->kmaxLocal;
  double *e     = lv->e;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        P(i, j, k) += E(i, j, k);
      }
    }
  }
}

/* Red-black SOR on one level, with that level's mesh and the global
 * checkerboard. */
static void smooth(
    Solver *s, MgLevelType *lv, double *p, const double *rhs, int sweeps)
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
      s->omega * 0.5 * (dx2 * dy2 * dz2) / (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);

  for (int sweep = 0; sweep < sweeps; sweep++) {
    for (int color = 0; color < 2; color++) {
      commExchange(&lv->comm, p);

      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          int iStart =
              colorRowStart(color, j, k, lv->iOffset, lv->jOffset, lv->kOffset);

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

      pressureBcApply(&s->bc, &lv->comm, p, imaxLocal, jmaxLocal, kmaxLocal);
    }
  }
}

/* r = rhs - A p on this level, written into lv->r. */
static void residualField(Solver *s, MgLevelType *lv, double *p, const double *rhs)
{
  int imaxLocal = lv->imaxLocal;
  int jmaxLocal = lv->jmaxLocal;
  int kmaxLocal = lv->kmaxLocal;

  double idx2 = 1.0 / (lv->dx * lv->dx);
  double idy2 = 1.0 / (lv->dy * lv->dy);
  double idz2 = 1.0 / (lv->dz * lv->dz);
  double *r   = lv->r;

  commExchange(&lv->comm, p);
  pressureBcApply(&s->bc, &lv->comm, p, imaxLocal, jmaxLocal, kmaxLocal);

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        R(i, j, k) = RHS(i, j, k) -
                     ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                         (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
                         (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2);
      }
    }
  }
}

static void vcycle(Solver *s, int level, double *p, const double *rhs)
{
  MgLevelType *levels = (MgLevelType *)s->mgLevels;
  MgLevelType *lv     = &levels[level];

  if (level == s->levels - 1) {
    /* Coarsest level: relax hard enough to stand in for a solve. */
    smooth(s, lv, p, rhs, s->presmooth + s->postsmooth);
    return;
  }

  MgLevelType *coarse = &levels[level + 1];

  smooth(s, lv, p, rhs, s->presmooth);
  residualField(s, lv, p, rhs);
  restrictMG(lv, coarse);

  /* The coarse problem is for the error, so it starts from zero. Leaving the
   * previous cycle's error in place made every cycle after the first solve a
   * different problem from the one it was given. */
  zeroField(coarse, coarse->e);
  vcycle(s, level + 1, coarse->e, coarse->r);

  prolongate(s, coarse, lv);
  correct(lv, p);
  pressureBcApply(&s->bc, &lv->comm, p, lv->imaxLocal, lv->jmaxLocal, lv->kmaxLocal);

  smooth(s, lv, p, rhs, s->postsmooth);
}

void initSolver(Solver *s, Discretization *d, Parameter *p)
{
  solverBaseInit(s, d, p);

  s->levels     = p->levels;
  s->presmooth  = p->presmooth;
  s->postsmooth = p->postsmooth;

  MgLevelType *levels = malloc((size_t)s->levels * sizeof(MgLevelType));

  /* Finest level: the decomposition the solver already has. */
  levels[0].comm      = *s->comm;
  levels[0].imaxLocal = s->comm->imaxLocal;
  levels[0].jmaxLocal = s->comm->jmaxLocal;
  levels[0].kmaxLocal = s->comm->kmaxLocal;
  levels[0].iOffset   = s->iOffset;
  levels[0].jOffset   = s->jOffset;
  levels[0].kOffset   = s->kOffset;
  levels[0].dx        = s->grid->dx;
  levels[0].dy        = s->grid->dy;
  levels[0].dz        = s->grid->dz;
  levels[0].cells     = (double)s->grid->imax * s->grid->jmax * s->grid->kmax;

  int built           = 1;

  for (int l = 1; l < s->levels; l++) {
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

  if (built < s->levels) {
    if (commIsMaster(s->comm)) {
      printf("Multigrid: requested %d levels, the decomposition supports %d\n",
          s->levels,
          built);
    }
    s->levels = built;
  }

  for (int l = 0; l < s->levels; l++) {
    size_t size   = levelSize(&levels[l]);
    levels[l].e   = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
    levels[l].r   = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
    zeroField(&levels[l], levels[l].e);
    zeroField(&levels[l], levels[l].r);
  }

  s->mgLevels = levels;

  if (commIsMaster(s->comm)) {
    printf("Using Multigrid solver with %d levels\n", s->levels);
  }
}

double solve(Solver *s, double *p, const double *rhs)
{
  MgLevelType *levels = (MgLevelType *)s->mgLevels;
  MgLevelType *fine   = &levels[FINEST_LEVEL];

  double epssq        = s->eps * s->eps;
  int cycles          = 0;

  TIMESTART

  /* Multigrid used to run exactly one V-cycle per call, reading eps and itermax
   * into locals it never used, so there was no tolerance for it to converge to
   * and no way for it to agree with the relaxation solvers. */
  double res = pressureResidualNorm(&fine->comm,
      &s->bc,
      p,
      rhs,
      fine->imaxLocal,
      fine->jmaxLocal,
      fine->kmaxLocal,
      fine->dx,
      fine->dy,
      fine->dz,
      fine->cells);

  while ((res >= epssq) && (cycles < s->itermax)) {
    vcycle(s, FINEST_LEVEL, p, rhs);

    res = pressureResidualNorm(&fine->comm,
        &s->bc,
        p,
        rhs,
        fine->imaxLocal,
        fine->jmaxLocal,
        fine->kmaxLocal,
        fine->dx,
        fine->dy,
        fine->dz,
        fine->cells);
    cycles++;
  }
  TIMESTOP(SOLVER);

#ifdef VERBOSE
  if (commIsMaster(s->comm)) {
    printf("Multigrid took %d cycles to reach %e\n", cycles, sqrt(res));
  }

  printProfile(s->comm, cycles);
#endif

  return res;
}

#if defined(TEST) && defined(SOLVER_mg)
/* Test seam. See the declarations in solver.h for why these exist. */

int mgTestLevels(Solver *s) { return s->levels; }

void mgTestLevelExtents(Solver *s, int level, int *im, int *jm, int *km)
{
  MgLevelType *lv = &((MgLevelType *)s->mgLevels)[level];
  *im             = lv->imaxLocal;
  *jm             = lv->jmaxLocal;
  *km             = lv->kmaxLocal;
}

void mgTestLevelMesh(Solver *s, int level, double *dx, double *dy, double *dz)
{
  MgLevelType *lv = &((MgLevelType *)s->mgLevels)[level];
  *dx             = lv->dx;
  *dy             = lv->dy;
  *dz             = lv->dz;
}

double *mgTestLevelE(Solver *s, int level)
{
  return ((MgLevelType *)s->mgLevels)[level].e;
}

double *mgTestLevelR(Solver *s, int level)
{
  return ((MgLevelType *)s->mgLevels)[level].r;
}

void mgTestRestrict(Solver *s, int level)
{
  MgLevelType *levels = (MgLevelType *)s->mgLevels;
  restrictMG(&levels[level], &levels[level + 1]);
}

void mgTestProlongate(Solver *s, int level)
{
  MgLevelType *levels = (MgLevelType *)s->mgLevels;
  prolongate(s, &levels[level + 1], &levels[level]);
}

void mgTestResidualField(Solver *s, int level, double *p, const double *rhs)
{
  residualField(s, &((MgLevelType *)s->mgLevels)[level], p, rhs);
}

void mgTestVcycle(Solver *s, double *p, const double *rhs)
{
  vcycle(s, FINEST_LEVEL, p, rhs);
}

void mgTestSmooth(Solver *s, int level, double *p, const double *rhs, int sweeps)
{
  smooth(s, &((MgLevelType *)s->mgLevels)[level], p, rhs, sweeps);
}
#endif
