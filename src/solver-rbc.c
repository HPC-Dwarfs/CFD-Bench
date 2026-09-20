/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "coloring.h"
#include "comm.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "surface-list.h"
#include "timing.h"
#include "util.h"

/* Compressed-layout macros for color-split arrays.
 * Nc is the number of compressed columns (ceil of half the row width incl. ghosts).
 * ic = i/2 maps original column index to compressed index. */
#define Nc ((imaxLocal + 3) / 2)
#define PRED(ic, j, k) pRed[(k) * Nc * (jmaxLocal + 2) + (j) * Nc + (ic)]
#define PBLACK(ic, j, k) pBlack[(k) * Nc * (jmaxLocal + 2) + (j) * Nc + (ic)]
#define RHSRED(ic, j, k) rhsRed[(k) * Nc * (jmaxLocal + 2) + (j) * Nc + (ic)]
#define RHSBLACK(ic, j, k) rhsBlack[(k) * Nc * (jmaxLocal + 2) + (j) * Nc + (ic)]

void initSolver(Solver *s, Discretization *d, Parameter *p)
{
  solverBaseInit(s, d, p);

  double dx2 = s->grid->dx * s->grid->dx;
  double dy2 = s->grid->dy * s->grid->dy;
  double dz2 = s->grid->dz * s->grid->dz;

  surfaceListBuild(&s->surface,
      s->comm->imaxLocal,
      s->comm->jmaxLocal,
      s->comm->kmaxLocal,
      s->iOffset,
      s->jOffset,
      s->kOffset,
      s->Ax,
      s->Ay,
      s->Az,
      s->Lambda,
      1.0 / dx2,
      1.0 / dy2,
      1.0 / dz2);
}

/*
 * The cut-cell correction in the compressed colour-split layout.
 *
 * Same operator as pressureCorrectSurface, addressed differently: a cell of one
 * colour lives in its own array and all six of its face neighbours live in the
 * other. The x neighbours share or straddle the compressed column depending on
 * the parity of the original index, which is the same bookkeeping the bulk
 * passes do with ioff.
 *
 * This is why every list entry carries its compressed index as well as its
 * linear one: one list serves the natural layout and this one.
 */
static void correctSurfaceCompressed(const SurfaceListType *list,
    int color,
    double *pRed,
    double *pBlack,
    const double *rhsRed,
    const double *rhsBlack,
    int imaxLocal,
    int jmaxLocal,
    double omega,
    double idx2,
    double idy2,
    double idz2,
    double *sweepRes)
{
  if (list->count == 0) {
    return;
  }

  /* Colour 0 is the even parity of (i + j + k), which is what this layout calls
   * red. */
  double *own          = (color == 0) ? pRed : pBlack;
  double *other        = (color == 0) ? pBlack : pRed;
  const double *rhsOwn = (color == 0) ? rhsRed : rhsBlack;

  int begin            = (color == 0) ? 0 : list->colorCount[0];
  int end              = begin + list->colorCount[color];
  int strideJ          = imaxLocal + 2;
  double delta         = 0.0;

  for (int e = begin; e < end; e++) {
    int i        = list->index[e] % strideJ;
    int j        = list->jc[e];
    int k        = list->kc[e];
    int ic       = list->ic[e];

    int icE      = (i & 1) ? ic + 1 : ic;
    int icW      = (i & 1) ? ic : ic - 1;

    double pOld  = list->saved[e];

    double pE    = other[(k)*Nc * (jmaxLocal + 2) + (j)*Nc + (icE)];
    double pW    = other[(k)*Nc * (jmaxLocal + 2) + (j)*Nc + (icW)];
    double pN    = other[(k)*Nc * (jmaxLocal + 2) + (j + 1) * Nc + (ic)];
    double pS    = other[(k)*Nc * (jmaxLocal + 2) + (j - 1) * Nc + (ic)];
    double pT    = other[(k + 1) * Nc * (jmaxLocal + 2) + (j)*Nc + (ic)];
    double pB    = other[(k - 1) * Nc * (jmaxLocal + 2) + (j)*Nc + (ic)];

    double rhsc  = rhsOwn[(k)*Nc * (jmaxLocal + 2) + (j)*Nc + (ic)];

    double rBulk = rhsc - ((pE - 2.0 * pOld + pW) * idx2 + (pN - 2.0 * pOld + pS) * idy2 +
                              (pT - 2.0 * pOld + pB) * idz2);
    delta -= rBulk * rBulk;

    if (list->solid[e]) {
      own[(k)*Nc * (jmaxLocal + 2) + (j)*Nc + (ic)] = 0.0;
      continue;
    }

    double r = list->lambda[e] * rhsc -
               (list->aE[e] * (pE - pOld) + list->aW[e] * (pW - pOld) +
                   list->aN[e] * (pN - pOld) + list->aS[e] * (pS - pOld) +
                   list->aT[e] * (pT - pOld) + list->aB[e] * (pB - pOld));

    own[(k)*Nc * (jmaxLocal + 2) + (j)*Nc + (ic)] = pOld - omega * r * list->invDiag[e];
    delta += r * r;
  }

  *sweepRes += delta;
}

/* Keep the listed cells of one colour before the bulk pass overwrites them. */
static void saveSurfaceCompressed(const SurfaceListType *list,
    int color,
    const double *pRed,
    const double *pBlack,
    int imaxLocal,
    int jmaxLocal)
{
  if (list->count == 0) {
    return;
  }

  const double *own = (color == 0) ? pRed : pBlack;
  int begin         = (color == 0) ? 0 : list->colorCount[0];
  int end           = begin + list->colorCount[color];

  for (int e = begin; e < end; e++) {
    ((SurfaceListType *)list)->saved[e] =
        own[(list->kc[e]) * Nc * (jmaxLocal + 2) + (list->jc[e]) * Nc + (list->ic[e])];
  }
}

double solve(Solver *s, double *p, const double *rhs)
{
  int imaxLocal = s->comm->imaxLocal;
  int jmaxLocal = s->comm->jmaxLocal;
  int kmaxLocal = s->comm->kmaxLocal;

  int imax      = s->grid->imax;
  int jmax      = s->grid->jmax;
  int kmax      = s->grid->kmax;

  int iOffset   = s->iOffset;
  int jOffset   = s->jOffset;
  int kOffset   = s->kOffset;

  double eps    = s->eps;
  int itermax   = s->itermax;
  double dx2    = s->grid->dx * s->grid->dx;
  double dy2    = s->grid->dy * s->grid->dy;
  double dz2    = s->grid->dz * s->grid->dz;
  double idx2   = 1.0 / dx2;
  double idy2   = 1.0 / dy2;
  double idz2   = 1.0 / dz2;

  double factor =
      s->omega * 0.5 * (dx2 * dy2 * dz2) / (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);
  double epssq = eps * eps;
  double cells = (double)imax * jmax * kmax;
  int it       = 0;

  PressureLevelType lv;
  pressureLevelFromSolver(s, &lv);

  /* See solverbase.c: accumulated in-sweep for loop control only, reset every
   * iteration, with the reported value computed once at the end. */
  double sweepRes = DBL_MAX;

#ifdef _MPI

  TIMESTART
  while ((sweepRes >= epssq) && (it < itermax)) {
    sweepRes = 0.0;

    for (int color = 0; color < 2; color++) {
      PROFILE(COMM, commExchange(s->comm, p));

      pressureSaveSurface(&lv, p, color);

#ifdef PROFILING
      double bulkStart = getTimeStamp();
#endif

      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          int iStart = colorRowStart(color, j, k, iOffset, jOffset, kOffset);

          for (int i = iStart; i < imaxLocal + 1; i += 2) {

            double r = RHS(i, j, k) -
                       ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                           (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
                           (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2);

            P(i, j, k) -= (factor * r);
            sweepRes += (r * r);
          }
        }
      }

#ifdef PROFILING
      T[SWEEP_BULK] += getTimeStamp() - bulkStart;
      C[SWEEP_BULK]++;
#endif

      pressureCorrectSurface(&lv, p, rhs, color, s->omega, &sweepRes);
      pressureBcApply(&s->bc, s->comm, p, imaxLocal, jmaxLocal, kmaxLocal);
    }

    commReduceAll(&sweepRes, SUM);
    sweepRes = sweepRes / s->fluidCells;
#ifdef DEBUG
    if (commIsMaster(s->comm)) {
      printf("%d Residuum: %e\n", it, sweepRes);
    }
#endif

    PROFILE(COMM, commExchange(s->comm, p));
    it++;
  }
  TIMESTOP(SOLVER);

#else /* !_MPI: SIMD-friendly color-split layout with stride-1 inner loops */

  int colorSize    = Nc * (jmaxLocal + 2) * (kmaxLocal + 2);
  double *pRed     = (double *)calloc(colorSize, sizeof(double));
  double *pBlack   = (double *)calloc(colorSize, sizeof(double));
  double *rhsRed   = (double *)calloc(colorSize, sizeof(double));
  double *rhsBlack = (double *)calloc(colorSize, sizeof(double));

  /* Gather: split p and rhs into separate red/black arrays */
  for (int k = 0; k <= kmaxLocal + 1; k++) {
    for (int j = 0; j <= jmaxLocal + 1; j++) {
      for (int i = 0; i <= imaxLocal + 1; i++) {
        int ic = i / 2;
        if ((i + j + k) % 2 == 0) {
          PRED(ic, j, k)   = P(i, j, k);
          RHSRED(ic, j, k) = RHS(i, j, k);
        } else {
          PBLACK(ic, j, k)   = P(i, j, k);
          RHSBLACK(ic, j, k) = RHS(i, j, k);
        }
      }
    }
  }

  /* A wall mirrors the interior into the halo, an outflow reflects it oddly.
   * Every boundary of this layout is a physical one, because the compressed
   * path is the serial path. */
  double sgnLeft   = (s->bc.type[LEFT] == OUTFLOW) ? -1.0 : 1.0;
  double sgnRight  = (s->bc.type[RIGHT] == OUTFLOW) ? -1.0 : 1.0;
  double sgnBottom = (s->bc.type[BOTTOM] == OUTFLOW) ? -1.0 : 1.0;
  double sgnTop    = (s->bc.type[TOP] == OUTFLOW) ? -1.0 : 1.0;
  double sgnFront  = (s->bc.type[FRONT] == OUTFLOW) ? -1.0 : 1.0;
  double sgnBack   = (s->bc.type[BACK] == OUTFLOW) ? -1.0 : 1.0;

  TIMESTART
  while ((sweepRes >= epssq) && (it < itermax)) {
    sweepRes = 0.0;

    /* Pass 0: update red cells (neighbors live in pBlack) */
    saveSurfaceCompressed(&s->surface, 0, pRed, pBlack, imaxLocal, jmaxLocal);

    for (int k = 1; k <= kmaxLocal; k++) {
      for (int j = 1; j <= jmaxLocal; j++) {
        int pjk   = (j + k) % 2;
        int ioff  = pjk - 1;
        int icBeg = 1 - pjk;
        int icEnd = (imaxLocal - pjk) / 2;

        for (int ic = icBeg; ic <= icEnd; ic++) {
          double r =
              RHSRED(ic, j, k) -
              ((PBLACK(ic + ioff + 1, j, k) - 2.0 * PRED(ic, j, k) +
                   PBLACK(ic + ioff, j, k)) *
                      idx2 +
                  (PBLACK(ic, j + 1, k) - 2.0 * PRED(ic, j, k) + PBLACK(ic, j - 1, k)) *
                      idy2 +
                  (PBLACK(ic, j, k + 1) - 2.0 * PRED(ic, j, k) + PBLACK(ic, j, k - 1)) *
                      idz2);

          PRED(ic, j, k) -= factor * r;
          sweepRes += r * r;
        }
      }
    }

    correctSurfaceCompressed(&s->surface,
        0,
        pRed,
        pBlack,
        rhsRed,
        rhsBlack,
        imaxLocal,
        jmaxLocal,
        s->omega,
        idx2,
        idy2,
        idz2,
        &sweepRes);

    /* Pass 1: update black cells (neighbors live in pRed) */
    saveSurfaceCompressed(&s->surface, 1, pRed, pBlack, imaxLocal, jmaxLocal);

    for (int k = 1; k <= kmaxLocal; k++) {
      for (int j = 1; j <= jmaxLocal; j++) {
        int pjk   = (j + k) % 2;
        int ioff  = -pjk;
        int icBeg = pjk;
        int icEnd = (imaxLocal - 1 + pjk) / 2;

        for (int ic = icBeg; ic <= icEnd; ic++) {
          double r =
              RHSBLACK(ic, j, k) -
              ((PRED(ic + ioff + 1, j, k) - 2.0 * PBLACK(ic, j, k) +
                   PRED(ic + ioff, j, k)) *
                      idx2 +
                  (PRED(ic, j + 1, k) - 2.0 * PBLACK(ic, j, k) + PRED(ic, j - 1, k)) *
                      idy2 +
                  (PRED(ic, j, k + 1) - 2.0 * PBLACK(ic, j, k) + PRED(ic, j, k - 1)) *
                      idz2);

          PBLACK(ic, j, k) -= factor * r;
          sweepRes += r * r;
        }
      }
    }

    correctSurfaceCompressed(&s->surface,
        1,
        pRed,
        pBlack,
        rhsRed,
        rhsBlack,
        imaxLocal,
        jmaxLocal,
        s->omega,
        idx2,
        idy2,
        idz2,
        &sweepRes);

    /* The pressure boundary condition in the compressed layout. The sign is
     * what makes it the setup's condition rather than a wall everywhere: +1
     * mirrors the interior into the halo (zero normal gradient), -1 reflects it
     * oddly (zero pressure on the boundary face).
     *
     * j/k boundaries: halo and interior differ by 1 in j or k, so they have
     * opposite colours but the same compressed index ic. */

    /* k-direction */
    for (int j = 1; j <= jmaxLocal; j++) {
      for (int ic = 0; ic < Nc; ic++) {
        PRED(ic, j, 0)               = sgnFront * PBLACK(ic, j, 1);
        PBLACK(ic, j, 0)             = sgnFront * PRED(ic, j, 1);
        PRED(ic, j, kmaxLocal + 1)   = sgnBack * PBLACK(ic, j, kmaxLocal);
        PBLACK(ic, j, kmaxLocal + 1) = sgnBack * PRED(ic, j, kmaxLocal);
      }
    }

    /* j-direction */
    for (int k = 1; k <= kmaxLocal; k++) {
      for (int ic = 0; ic < Nc; ic++) {
        PRED(ic, 0, k)               = sgnBottom * PBLACK(ic, 1, k);
        PBLACK(ic, 0, k)             = sgnBottom * PRED(ic, 1, k);
        PRED(ic, jmaxLocal + 1, k)   = sgnTop * PBLACK(ic, jmaxLocal, k);
        PBLACK(ic, jmaxLocal + 1, k) = sgnTop * PRED(ic, jmaxLocal, k);
      }
    }

    /* i-direction: i=0 and i=1 have opposite colours, both at ic=0.
     * i=imax+1 and i=imax also have opposite colours. */
    for (int k = 1; k <= kmaxLocal; k++) {
      for (int j = 1; j <= jmaxLocal; j++) {
        int pjk = (j + k) % 2;

        /* Left halo from the first interior cell */
        if (pjk == 0) {
          PRED(0, j, k) = sgnLeft * PBLACK(0, j, k);
        } else {
          PBLACK(0, j, k) = sgnLeft * PRED(0, j, k);
        }

        /* Right halo from the last interior cell */
        int icg = (imaxLocal + 1) / 2;
        int ici = imaxLocal / 2;
        if ((imaxLocal + 1 + j + k) % 2 == 0) {
          PRED(icg, j, k) = sgnRight * PBLACK(ici, j, k);
        } else {
          PBLACK(icg, j, k) = sgnRight * PRED(ici, j, k);
        }
      }
    }

    sweepRes = sweepRes / s->fluidCells;
#ifdef DEBUG
    if (commIsMaster(s->comm)) {
      printf("%d Residuum: %e\n", it, sweepRes);
    }
#endif
    it++;
  }
  TIMESTOP(SOLVER);

  /* Scatter: merge pRed/pBlack back into p */
  for (int k = 0; k <= kmaxLocal + 1; k++) {
    for (int j = 0; j <= jmaxLocal + 1; j++) {
      for (int i = 0; i <= imaxLocal + 1; i++) {
        int ic = i / 2;
        if ((i + j + k) % 2 == 0) {
          P(i, j, k) = PRED(ic, j, k);
        } else {
          P(i, j, k) = PBLACK(ic, j, k);
        }
      }
    }
  }

  free(pRed);
  free(pBlack);
  free(rhsRed);
  free(rhsBlack);

#endif /* _MPI */

  /* Computed from the field that is about to be returned -- in the compressed
   * path that means after the scatter above, not from the split arrays. */
  double res = pressureResidualNorm(&lv, p, rhs);

#ifdef VERBOSE
  if (commIsMaster(s->comm)) {
    printf("Solver took %d iterations to reach %e\n", it, sqrt(res));
  }

  printProfile(s->comm, it);
#endif

  return res;
}
