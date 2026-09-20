/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <stdlib.h>

#include "coloring.h"
#include "surface-list.h"

#define IDX(i, j, k)                                                                     \
  ((k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i))

/* True when the plain stencil would get this cell wrong. */
static int needsCorrection(int i,
    int j,
    int k,
    int imaxLocal,
    int jmaxLocal,
    const double *Ax,
    const double *Ay,
    const double *Az,
    const double *Lambda)
{
  double aE = Ax[IDX(i, j, k)];
  double aW = Ax[IDX(i - 1, j, k)];
  double aN = Ay[IDX(i, j, k)];
  double aS = Ay[IDX(i, j - 1, k)];
  double aT = Az[IDX(i, j, k)];
  double aB = Az[IDX(i, j, k - 1)];

  if (Lambda[IDX(i, j, k)] > 0.0) {
    /* A fluid cell matters once any of its faces is shut. */
    return (aE == 0.0 || aW == 0.0 || aN == 0.0 || aS == 0.0 || aT == 0.0 || aB == 0.0);
  }

  /* A solid cell matters only where it touches fluid. Its own faces are all
   * closed, so the neighbours have to be looked at directly. */
  return (Lambda[IDX(i + 1, j, k)] > 0.0 || Lambda[IDX(i - 1, j, k)] > 0.0 ||
          Lambda[IDX(i, j + 1, k)] > 0.0 || Lambda[IDX(i, j - 1, k)] > 0.0 ||
          Lambda[IDX(i, j, k + 1)] > 0.0 || Lambda[IDX(i, j, k - 1)] > 0.0);
}

void surfaceListBuild(SurfaceListType *list,
    int imaxLocal,
    int jmaxLocal,
    int kmaxLocal,
    int iOffset,
    int jOffset,
    int kOffset,
    const double *Ax,
    const double *Ay,
    const double *Az,
    const double *Lambda,
    double idx2,
    double idy2,
    double idz2)
{
  int count = 0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (needsCorrection(i, j, k, imaxLocal, jmaxLocal, Ax, Ay, Az, Lambda)) {
          ++count;
        }
      }
    }
  }

  list->count         = count;
  list->colorCount[0] = 0;
  list->colorCount[1] = 0;

  if (count == 0) {
    list->index = NULL;
    list->ic = list->jc = list->kc = NULL;
    list->color = list->solid = NULL;
    list->aE = list->aW = list->aN = list->aS = list->aT = list->aB = NULL;
    list->invDiag = list->lambda = list->saved = NULL;
    return;
  }

  list->index   = malloc((size_t)count * sizeof(int));
  list->ic      = malloc((size_t)count * sizeof(int));
  list->jc      = malloc((size_t)count * sizeof(int));
  list->kc      = malloc((size_t)count * sizeof(int));
  list->color   = malloc((size_t)count * sizeof(unsigned char));
  list->solid   = malloc((size_t)count * sizeof(unsigned char));
  list->aE      = malloc((size_t)count * sizeof(double));
  list->aW      = malloc((size_t)count * sizeof(double));
  list->aN      = malloc((size_t)count * sizeof(double));
  list->aS      = malloc((size_t)count * sizeof(double));
  list->aT      = malloc((size_t)count * sizeof(double));
  list->aB      = malloc((size_t)count * sizeof(double));
  list->invDiag = malloc((size_t)count * sizeof(double));
  list->lambda  = malloc((size_t)count * sizeof(double));
  list->saved   = malloc((size_t)count * sizeof(double));

  int at        = 0;

  /* Colour 0 first, then colour 1; each segment row-major, so the pass over it
   * walks memory forwards and is the same on every run. */
  for (int color = 0; color < 2; color++) {
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (colorOf(i, j, k, iOffset, jOffset, kOffset) != color) {
            continue;
          }
          if (!needsCorrection(i, j, k, imaxLocal, jmaxLocal, Ax, Ay, Az, Lambda)) {
            continue;
          }

          double aE        = Ax[IDX(i, j, k)] * idx2;
          double aW        = Ax[IDX(i - 1, j, k)] * idx2;
          double aN        = Ay[IDX(i, j, k)] * idy2;
          double aS        = Ay[IDX(i, j - 1, k)] * idy2;
          double aT        = Az[IDX(i, j, k)] * idz2;
          double aB        = Az[IDX(i, j, k - 1)] * idz2;
          double diag      = aE + aW + aN + aS + aT + aB;
          double lambda    = Lambda[IDX(i, j, k)];

          list->index[at]  = IDX(i, j, k);
          list->ic[at]     = i / 2;
          list->jc[at]     = j;
          list->kc[at]     = k;
          list->color[at]  = (unsigned char)color;
          list->solid[at]  = (lambda == 0.0);
          list->aE[at]     = aE;
          list->aW[at]     = aW;
          list->aN[at]     = aN;
          list->aS[at]     = aS;
          list->aT[at]     = aT;
          list->aB[at]     = aB;
          list->lambda[at] = lambda;
          /* A solid cell has no open face, so it has no diagonal to invert; it
           * is an identity row and is simply held at zero. A fluid cell with no
           * open face would be a sealed pocket, which the connectivity check
           * refuses before a solver ever sees it. */
          list->invDiag[at] = (diag > 0.0) ? 1.0 / diag : 0.0;
          list->saved[at]   = 0.0;

          ++at;
          ++list->colorCount[color];
        }
      }
    }
  }
}

void surfaceListFree(SurfaceListType *list)
{
  free(list->index);
  free(list->ic);
  free(list->jc);
  free(list->kc);
  free(list->color);
  free(list->solid);
  free(list->aE);
  free(list->aW);
  free(list->aN);
  free(list->aS);
  free(list->aT);
  free(list->aB);
  free(list->invDiag);
  free(list->lambda);
  free(list->saved);

  list->count         = 0;
  list->colorCount[0] = 0;
  list->colorCount[1] = 0;
}
