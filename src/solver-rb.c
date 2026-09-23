/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "coloring.h"
#include "comm.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "surface-list.h"
#include "timing.h"
#include "util.h"

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

double solve(Solver *s, double *p, const double *rhs)
{
  int imaxLocal = s->comm->imaxLocal;
  int jmaxLocal = s->comm->jmaxLocal;
  int kmaxLocal = s->comm->kmaxLocal;

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
  int it       = 0;

  PressureLevelType lv;
  pressureLevelFromSolver(s, &lv);

  /* Accumulated inside the sweep, so it describes the field part-way through
   * the iteration. That is what decides when to stop; the number reported at
   * the end comes from a dedicated pass. Reset every iteration -- it used to be
   * seeded at 1.0 and never cleared, so it only ever grew. */
  double sweepRes = DBL_MAX;

  TIMESTART
  while (solveContinues(sweepRes, epssq, it, itermax)) {
    sweepRes = 0.0;

    for (int color = 0; color < 2; color++) {
      PROFILE(COMM, commExchange(s->comm, p));

      /* The interior sweep reads no geometry at all, so its cost and its
       * instruction mix are the same with an obstacle and without one. What it
       * gets wrong is confined to the surface list, which is O(obstacle
       * surface) rather than O(domain). */
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

  double res = pressureResidualNorm(&lv, p, rhs);

  res = solveReport(s, "Solver", "iterations", it, sweepRes, res);

#ifdef VERBOSE
  printProfile(s->comm, it);
#endif

  return res;
}
