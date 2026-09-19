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

void initSolver(Solver *s, Discretization *d, Parameter *p) { solverBaseInit(s, d, p); }

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
  double epssq    = eps * eps;
  double cells    = (double)imax * jmax * kmax;
  int it          = 0;

  /* See solverbase.c: accumulated in-sweep for loop control only, reset every
   * iteration, with the reported value computed once at the end. */
  double sweepRes = DBL_MAX;

#ifdef _MPI

  TIMESTART
  while ((sweepRes >= epssq) && (it < itermax)) {
    sweepRes = 0.0;

    for (int color = 0; color < 2; color++) {
      PROFILE(COMM, commExchange(s->comm, p));

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

      pressureBcApply(&s->bc, s->comm, p, imaxLocal, jmaxLocal, kmaxLocal);
    }

    commReduceAll(&sweepRes, SUM);
    sweepRes = sweepRes / cells;
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

    /* Pass 1: update black cells (neighbors live in pRed) */
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

    sweepRes = sweepRes / cells;
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
  double res = pressureResidualNorm(s->comm,
      &s->bc,
      p,
      rhs,
      imaxLocal,
      jmaxLocal,
      kmaxLocal,
      s->grid->dx,
      s->grid->dy,
      s->grid->dz,
      cells);

#ifdef VERBOSE
  if (commIsMaster(s->comm)) {
    printf("Solver took %d iterations to reach %e\n", it, sqrt(res));
  }

  printProfile(s->comm, it);
#endif

  return res;
}
