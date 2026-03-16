/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "comm.h"
#include "parameter.h"
#include "profiler.h"
#include "solver.h"
#include "timing.h"
#include "util.h"

/* Compressed-layout macros for color-split arrays.
 * Nc is the number of compressed columns (ceil of half the row width incl. ghosts).
 * ic = i/2 maps original column index to compressed index. */
#define Nc              ((imaxLocal + 3) / 2)
#define PRED(ic, j, k)  pRed[(k) * Nc * (jmaxLocal + 2) + (j) * Nc + (ic)]
#define PBLACK(ic, j, k) pBlack[(k) * Nc * (jmaxLocal + 2) + (j) * Nc + (ic)]
#define RHSRED(ic, j, k) rhsRed[(k) * Nc * (jmaxLocal + 2) + (j) * Nc + (ic)]
#define RHSBLACK(ic, j, k) rhsBlack[(k) * Nc * (jmaxLocal + 2) + (j) * Nc + (ic)]

void initSolver(Solver *s, Discretization *d, Parameter *p)
{
  s->eps     = p->eps;
  s->omega   = p->omg;
  s->itermax = p->itermax;
  s->grid    = &d->grid;
  s->comm    = &d->comm;
  s->problem = p->name;
}

double solve(Solver *s, double *p, const double *rhs)
{
  int imaxLocal = s->comm->imaxLocal;
  int jmaxLocal = s->comm->jmaxLocal;
  int kmaxLocal = s->comm->kmaxLocal;

  int imax      = s->grid->imax;
  int jmax      = s->grid->jmax;
  int kmax      = s->grid->kmax;

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
  int it       = 0;
  double res   = 1.0;

#ifdef _MPI
  int pass, ksw, jsw, isw;

  TIMESTART
  while ((res >= epssq) && (it < itermax)) {
    ksw = 1;

    for (pass = 0; pass < 2; pass++) {
      jsw = ksw;
      PROFILE(COMM, commExchange(s->comm, p));

      for (int k = 1; k < kmaxLocal + 1; k++) {
        isw = jsw;
        for (int j = 1; j < jmaxLocal + 1; j++) {
          for (int i = isw; i < imaxLocal + 1; i += 2) {

            double r = RHS(i, j, k) -
                       ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                           (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
                           (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2);

            P(i, j, k) -= (factor * r);
            res += (r * r);
          }
          isw = 3 - isw;
        }
        jsw = 3 - jsw;
      }
      ksw = 3 - ksw;
    }

    if (commIsBoundary(s->comm, FRONT)) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          P(i, j, 0) = P(i, j, 1);
        }
      }
    }

    if (commIsBoundary(s->comm, BACK)) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          P(i, j, kmaxLocal + 1) = P(i, j, kmaxLocal);
        }
      }
    }

    if (commIsBoundary(s->comm, BOTTOM)) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          P(i, 0, k) = P(i, 1, k);
        }
      }
    }

    if (commIsBoundary(s->comm, TOP)) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          P(i, jmaxLocal + 1, k) = P(i, jmaxLocal, k);
        }
      }
    }

    if (commIsBoundary(s->comm, LEFT)) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          P(0, j, k) = P(1, j, k);
        }
      }
    }

    if (commIsBoundary(s->comm, RIGHT)) {
      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          P(imaxLocal + 1, j, k) = P(imaxLocal, j, k);
        }
      }
    }

    commReduceAll(&res, SUM);
    res = res / (double)(imax * jmax * kmax);
#ifdef DEBUG
    if (commIsMaster(&s->comm)) {
      printf("%d Residuum: %e\n", it, res);
    }
#endif

    PROFILE(COMM, commExchange(s->comm, p));
    it++;
  }
  TIMESTOP(SOLVER);

#else /* !_MPI: SIMD-friendly color-split layout with stride-1 inner loops */

  int colorSize = Nc * (jmaxLocal + 2) * (kmaxLocal + 2);
  double *pRed      = (double *)calloc(colorSize, sizeof(double));
  double *pBlack    = (double *)calloc(colorSize, sizeof(double));
  double *rhsRed    = (double *)calloc(colorSize, sizeof(double));
  double *rhsBlack  = (double *)calloc(colorSize, sizeof(double));

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

  TIMESTART
  while ((res >= epssq) && (it < itermax)) {

    /* Pass 0: update red cells (neighbors live in pBlack) */
    for (int k = 1; k <= kmaxLocal; k++) {
      for (int j = 1; j <= jmaxLocal; j++) {
        int pjk    = (j + k) % 2;
        int ioff   = pjk - 1;
        int icBeg = 1 - pjk;
        int icEnd = (imaxLocal - pjk) / 2;

        for (int ic = icBeg; ic <= icEnd; ic++) {
          double r = RHSRED(ic, j, k) -
              ((PBLACK(ic + ioff + 1, j, k) - 2.0 * PRED(ic, j, k) +
                   PBLACK(ic + ioff, j, k)) * idx2 +
               (PBLACK(ic, j + 1, k) - 2.0 * PRED(ic, j, k) +
                   PBLACK(ic, j - 1, k)) * idy2 +
               (PBLACK(ic, j, k + 1) - 2.0 * PRED(ic, j, k) +
                   PBLACK(ic, j, k - 1)) * idz2);

          PRED(ic, j, k) -= factor * r;
          res += r * r;
        }
      }
    }

    /* Pass 1: update black cells (neighbors live in pRed) */
    for (int k = 1; k <= kmaxLocal; k++) {
      for (int j = 1; j <= jmaxLocal; j++) {
        int pjk    = (j + k) % 2;
        int ioff   = -pjk;
        int icBeg = pjk;
        int icEnd = (imaxLocal - 1 + pjk) / 2;

        for (int ic = icBeg; ic <= icEnd; ic++) {
          double r = RHSBLACK(ic, j, k) -
              ((PRED(ic + ioff + 1, j, k) - 2.0 * PBLACK(ic, j, k) +
                   PRED(ic + ioff, j, k)) * idx2 +
               (PRED(ic, j + 1, k) - 2.0 * PBLACK(ic, j, k) +
                   PRED(ic, j - 1, k)) * idy2 +
               (PRED(ic, j, k + 1) - 2.0 * PBLACK(ic, j, k) +
                   PRED(ic, j, k - 1)) * idz2);

          PBLACK(ic, j, k) -= factor * r;
          res += r * r;
        }
      }
    }

    /* Boundary conditions in compressed layout (Neumann: ghost = interior).
     * j/k BCs: ghost and interior differ by 1 in j or k, so they have
     * opposite colors but the same compressed index ic. */

    /* k-direction BCs */
    for (int j = 1; j <= jmaxLocal; j++) {
      for (int ic = 0; ic < Nc; ic++) {
        PRED(ic, j, 0)              = PBLACK(ic, j, 1);
        PBLACK(ic, j, 0)            = PRED(ic, j, 1);
        PRED(ic, j, kmaxLocal + 1)  = PBLACK(ic, j, kmaxLocal);
        PBLACK(ic, j, kmaxLocal + 1) = PRED(ic, j, kmaxLocal);
      }
    }

    /* j-direction BCs */
    for (int k = 1; k <= kmaxLocal; k++) {
      for (int ic = 0; ic < Nc; ic++) {
        PRED(ic, 0, k)              = PBLACK(ic, 1, k);
        PBLACK(ic, 0, k)            = PRED(ic, 1, k);
        PRED(ic, jmaxLocal + 1, k)  = PBLACK(ic, jmaxLocal, k);
        PBLACK(ic, jmaxLocal + 1, k) = PRED(ic, jmaxLocal, k);
      }
    }

    /* i-direction BCs: i=0 and i=1 have opposite colors, both at ic=0.
     * i=imax+1 and i=imax also have opposite colors. */
    for (int k = 1; k <= kmaxLocal; k++) {
      for (int j = 1; j <= jmaxLocal; j++) {
        int pjk = (j + k) % 2;

        /* Left: P(0,j,k) = P(1,j,k) */
        if (pjk == 0) {
          PRED(0, j, k) = PBLACK(0, j, k);
        } else {
          PBLACK(0, j, k) = PRED(0, j, k);
        }

        /* Right: P(imax+1,j,k) = P(imax,j,k) */
        int icg = (imaxLocal + 1) / 2;
        int ici = imaxLocal / 2;
        if ((imaxLocal + 1 + j + k) % 2 == 0) {
          PRED(icg, j, k) = PBLACK(ici, j, k);
        } else {
          PBLACK(icg, j, k) = PRED(ici, j, k);
        }
      }
    }

    res = res / (double)(imax * jmax * kmax);
#ifdef DEBUG
    if (commIsMaster(&s->comm)) {
      printf("%d Residuum: %e\n", it, res);
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

#ifdef VERBOSE
  if (commIsMaster(s->comm)) {
    printf("Solver took %d iterations to reach %f\n", it, sqrt(res));
  }

  printProfile(s->comm, it);
#endif

  return res;
}
