/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Not named solver-base.c on purpose: the Makefile treats solver-*.c as the
 * mutually exclusive solver variants and links exactly one of them.
 */
#include <stddef.h>

#include "profiler.h"
#include "timing.h"
#include "solver.h"
#include "util.h"

#ifdef PROFILING
void pressureSweepTimes(double *bulk, double *surface)
{
  *bulk    = T[SWEEP_BULK];
  *surface = T[SWEEP_SURFACE];
}

void pressureSweepTimesReset(void)
{
  T[SWEEP_BULK]    = 0.0;
  T[SWEEP_SURFACE] = 0.0;
  C[SWEEP_BULK]    = 0;
  C[SWEEP_SURFACE] = 0;
}
#endif

/*
 * The pressure operator, in one place.
 *
 *   sum over the six faces of  A_f (p_nb - p_c) / h^2  =  Lambda_c * rhs_c
 *
 * Each face contributes an equal and opposite coefficient to the two cells it
 * separates, so the operator is symmetric for any face weights -- which is why
 * fractional apertures will need no change here -- and its diagonal is
 * -sum_f A_f / h^2, per cell rather than domain-wide.
 *
 * A cell with zero volume fraction is an identity row with a zero right-hand
 * side. It stays in the vector so the arrays stay rectangular and the interior
 * sweep needs no index compression; its solution is zero.
 */

double pressureResidualNorm(const PressureLevelType *lv, double *p, const double *rhs)
{
  int imaxLocal        = lv->imaxLocal;
  int jmaxLocal        = lv->jmaxLocal;
  int kmaxLocal        = lv->kmaxLocal;

  const double *Ax     = lv->Ax;
  const double *Ay     = lv->Ay;
  const double *Az     = lv->Az;
  const double *Lambda = lv->Lambda;

  double idx2          = 1.0 / (lv->dx * lv->dx);
  double idy2          = 1.0 / (lv->dy * lv->dy);
  double idz2          = 1.0 / (lv->dz * lv->dz);
  double res           = 0.0;

  commExchange(lv->comm, p);
  pressureBcApply(lv->bc, lv->comm, p, imaxLocal, jmaxLocal, kmaxLocal);

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        /* Solid cells contribute nothing: the residual describes the fluid
         * problem, so growing the solid fraction must not change it. */
        if (LAM(i, j, k) == 0.0) {
          continue;
        }

        double pc = P(i, j, k);
        double r  = LAM(i, j, k) * RHS(i, j, k) -
                    (AX(i, j, k) * (P(i + 1, j, k) - pc) * idx2 +
                        AX(i - 1, j, k) * (P(i - 1, j, k) - pc) * idx2 +
                        AY(i, j, k) * (P(i, j + 1, k) - pc) * idy2 +
                        AY(i, j - 1, k) * (P(i, j - 1, k) - pc) * idy2 +
                        AZ(i, j, k) * (P(i, j, k + 1) - pc) * idz2 +
                        AZ(i, j, k - 1) * (P(i, j, k - 1) - pc) * idz2);

        res += r * r;
      }
    }
  }

  commReduceAll(&res, SUM);
  return res / lv->fluidCells;
}

/*
 * Remember what the listed cells of one colour hold, before the bulk sweep
 * overwrites them.
 *
 * The interior sweep is deliberately blind to geometry, so it computes the
 * wrong value for every cell on the list. Keeping the old value is what lets
 * the correction pass compute the right one from the right starting point,
 * instead of trying to undo an update it did not make.
 */
void pressureSaveSurface(const PressureLevelType *lv, const double *p, int color)
{
  SurfaceListType *list = lv->list;

  if (list == NULL || list->count == 0) {
    return;
  }

  int begin = (color == 0) ? 0 : list->colorCount[0];
  int end   = begin + list->colorCount[color];

  for (int e = begin; e < end; e++) {
    list->saved[e] = p[list->index[e]];
  }
}

/*
 * Relax the listed cells of one colour with the coefficients their geometry
 * actually implies, replacing whatever the bulk sweep left there.
 *
 * sweepRes is the residual the bulk sweep accumulated. For a listed cell that
 * contribution is wrong, so the bulk residual is recomputed here and swapped
 * for the correct one, leaving the total describing the fluid problem.
 */
void pressureCorrectSurface(const PressureLevelType *lv,
    double *p,
    const double *rhs,
    int color,
    double omega,
    double *sweepRes)
{
  SurfaceListType *list = lv->list;

  if (list == NULL || list->count == 0) {
    return;
  }

  int imaxLocal = lv->imaxLocal;
  int jmaxLocal = lv->jmaxLocal;

  double idx2   = 1.0 / (lv->dx * lv->dx);
  double idy2   = 1.0 / (lv->dy * lv->dy);
  double idz2   = 1.0 / (lv->dz * lv->dz);

  int strideJ   = imaxLocal + 2;
  int strideK   = (imaxLocal + 2) * (jmaxLocal + 2);

  int begin     = (color == 0) ? 0 : list->colorCount[0];
  int end       = begin + list->colorCount[color];

  double delta  = 0.0;

#ifdef PROFILING
  double surfaceStart = getTimeStamp();
#endif

  for (int e = begin; e < end; e++) {
    int idx     = list->index[e];
    double pOld = list->saved[e];

    double pE   = p[idx + 1];
    double pW   = p[idx - 1];
    double pN   = p[idx + strideJ];
    double pS   = p[idx - strideJ];
    double pT   = p[idx + strideK];
    double pB   = p[idx - strideK];

    /* Undo this cell's contribution to the bulk sweep's residual. */
    double rBulk =
        rhs[idx] - ((pE - 2.0 * pOld + pW) * idx2 + (pN - 2.0 * pOld + pS) * idy2 +
                       (pT - 2.0 * pOld + pB) * idz2);
    delta -= rBulk * rBulk;

    if (list->solid[e]) {
      /* An identity row with a zero right-hand side. Its residual is zero and
       * it contributes nothing to the fluid problem. */
      p[idx] = 0.0;
      continue;
    }

    double r = list->lambda[e] * rhs[idx] -
               (list->aE[e] * (pE - pOld) + list->aW[e] * (pW - pOld) +
                   list->aN[e] * (pN - pOld) + list->aS[e] * (pS - pOld) +
                   list->aT[e] * (pT - pOld) + list->aB[e] * (pB - pOld));

    p[idx]   = pOld - omega * r * list->invDiag[e];
    delta += r * r;
  }

  *sweepRes += delta;

#ifdef PROFILING
  T[SWEEP_SURFACE] += getTimeStamp() - surfaceStart;
  C[SWEEP_SURFACE]++;
#endif
}

/*
 * Apply the operator to the listed cells with the coefficients their geometry
 * actually implies, replacing whatever the bulk pass left there.
 *
 * Deliberately a sibling of pressureCorrectSurface rather than a reuse of it:
 * that one relaxes -- reads saved, writes a new p, and swaps a residual
 * contribution -- while this one applies, reading x and overwriting y. What
 * they share is the coefficients, which is what tests/checks/operator.c pins
 * by checking the application against pressureResidualNorm at cut cells.
 *
 * One pass over every entry, no colour split: an apply reads only x, so no
 * entry can see another's output and the traversal order cannot matter.
 */
void pressureApplySurface(const PressureLevelType *lv, const double *x, double *y)
{
  SurfaceListType *list = lv->list;

  if (list == NULL || list->count == 0) {
    return;
  }

  int imaxLocal = lv->imaxLocal;
  int jmaxLocal = lv->jmaxLocal;

  int strideJ   = imaxLocal + 2;
  int strideK   = (imaxLocal + 2) * (jmaxLocal + 2);

#ifdef PROFILING
  double surfaceStart = getTimeStamp();
#endif

  for (int e = 0; e < list->count; e++) {
    int idx = list->index[e];

    if (list->solid[e]) {
      /* An identity row applied to a vector whose solid entries are zero. */
      y[idx] = 0.0;
      continue;
    }

    double xc = x[idx];

    y[idx]    = -(list->aE[e] * (x[idx + 1] - xc) + list->aW[e] * (x[idx - 1] - xc) +
                list->aN[e] * (x[idx + strideJ] - xc) +
                list->aS[e] * (x[idx - strideJ] - xc) +
                list->aT[e] * (x[idx + strideK] - xc) +
                list->aB[e] * (x[idx - strideK] - xc));
  }

#ifdef PROFILING
  T[SWEEP_SURFACE] += getTimeStamp() - surfaceStart;
  C[SWEEP_SURFACE]++;
#endif
}

/* The operand and the result of an application, indexed like every other field.
 * util.h names a macro per field it knows about; these two are local because
 * only the application uses them. */
#define X(i, j, k)                                                                       \
  x[((k) * (imaxLocal + 2) * (jmaxLocal + 2)) + ((j) * (imaxLocal + 2)) + (i)]
#define Y(i, j, k)                                                                       \
  y[((k) * (imaxLocal + 2) * (jmaxLocal + 2)) + ((j) * (imaxLocal + 2)) + (i)]

/*
 * y = -A x, matrix-free, in the same two passes the relaxation sweep uses.
 *
 * The bulk pass reads no geometry at all, so an application costs the same with
 * a body and without one -- the property the whole embedded-boundary
 * representation exists to preserve, and the reason a Krylov solver, which
 * applies the operator every iteration rather than once per solve, does not
 * forfeit it. What the bulk pass gets wrong is confined to the surface list,
 * which is O(body surface).
 */
void pressureApplyOperator(const PressureLevelType *lv, double *x, double *y)
{
  int imaxLocal = lv->imaxLocal;
  int jmaxLocal = lv->jmaxLocal;
  int kmaxLocal = lv->kmaxLocal;

  double idx2   = 1.0 / (lv->dx * lv->dx);
  double idy2   = 1.0 / (lv->dy * lv->dy);
  double idz2   = 1.0 / (lv->dz * lv->dz);

  PROFILE(COMM, commExchange(lv->comm, x));
  pressureBcApply(lv->bc, lv->comm, x, imaxLocal, jmaxLocal, kmaxLocal);

#ifdef PROFILING
  double bulkStart = getTimeStamp();
#endif

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        double xc = X(i, j, k);

        Y(i, j, k) = -((X(i + 1, j, k) - 2.0 * xc + X(i - 1, j, k)) * idx2 +
                       (X(i, j + 1, k) - 2.0 * xc + X(i, j - 1, k)) * idy2 +
                       (X(i, j, k + 1) - 2.0 * xc + X(i, j, k - 1)) * idz2);
      }
    }
  }

#ifdef PROFILING
  T[SWEEP_BULK] += getTimeStamp() - bulkStart;
  C[SWEEP_BULK]++;
#endif

  pressureApplySurface(lv, x, y);
}

void solverBaseInit(Solver *s, Discretization *d, Parameter *p)
{
  s->eps      = p->eps;
  s->omega    = p->omg;
  s->itermax  = p->itermax;
  s->grid     = &d->grid;
  s->comm     = &d->comm;
  s->problem  = p->name;

  s->Ax       = d->Ax;
  s->Ay       = d->Ay;
  s->Az       = d->Az;
  s->Lambda   = d->Lambda;

  s->bcLeft   = p->bcLeft;
  s->bcRight  = p->bcRight;
  s->bcBottom = p->bcBottom;
  s->bcTop    = p->bcTop;
  s->bcFront  = p->bcFront;
  s->bcBack   = p->bcBack;

  pressureBcInit(
      &s->bc, p->bcLeft, p->bcRight, p->bcBottom, p->bcTop, p->bcFront, p->bcBack);

  int offsets[NDIMS] = { 0, 0, 0 };
  commGetOffsets(s->comm, offsets, p->kmax, p->jmax, p->imax);
  s->iOffset = offsets[IDIM];
  s->jOffset = offsets[JDIM];
  s->kOffset = offsets[KDIM];

  /* Residuals are normalized by the fluid cell count, not the cell count, so
   * that the reported norm does not shrink just because more of the domain is
   * solid. */
  int imaxLocal        = s->comm->imaxLocal;
  int jmaxLocal        = s->comm->jmaxLocal;
  int kmaxLocal        = s->comm->kmaxLocal;
  const double *Lambda = s->Lambda;
  double fluid         = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) > 0.0) {
          fluid += 1.0;
        }
      }
    }
  }

  commReduceAll(&fluid, SUM);
  s->fluidCells = (fluid > 0.0) ? fluid : 1.0;
}

void pressureLevelFromSolver(const Solver *s, PressureLevelType *lv)
{
  lv->comm       = s->comm;
  lv->bc         = &s->bc;
  lv->Ax         = s->Ax;
  lv->Ay         = s->Ay;
  lv->Az         = s->Az;
  lv->Lambda     = s->Lambda;
  lv->list       = (SurfaceListType *)&s->surface;
  lv->imaxLocal  = s->comm->imaxLocal;
  lv->jmaxLocal  = s->comm->jmaxLocal;
  lv->kmaxLocal  = s->comm->kmaxLocal;
  lv->dx         = s->grid->dx;
  lv->dy         = s->grid->dy;
  lv->dz         = s->grid->dz;
  lv->fluidCells = s->fluidCells;
}
