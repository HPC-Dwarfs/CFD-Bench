/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Multigrid as a solver: build a hierarchy, then run cycles on it until the
 * residual meets the tolerance.
 *
 * The hierarchy, the transfer operators and the cycle itself live in
 * multigrid.c, which every build links. This file is only the driver -- the
 * part that is specific to multigrid being the solver rather than a
 * preconditioner for one.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "multigrid.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "timing.h"
#include "util.h"

void initSolver(Solver *s, Discretization *d, Parameter *p)
{
  solverBaseInit(s, d, p);

  s->levels                 = p->levels;
  s->presmooth              = p->presmooth;
  s->postsmooth             = p->postsmooth;

  MultigridType *mg         = malloc(sizeof(MultigridType));

  MultigridSpecType spec    = { .comm = s->comm,
       .bc                          = &s->bc,
       .grid                        = s->grid,
       .Ax                          = s->Ax,
       .Ay                          = s->Ay,
       .Az                          = s->Az,
       .Lambda                      = s->Lambda,
       .iOffset                     = s->iOffset,
       .jOffset                     = s->jOffset,
       .kOffset                     = s->kOffset,
       .fluidCells                  = s->fluidCells,
       .levels                      = p->levels,
       .presmooth                   = p->presmooth,
       .postsmooth                  = p->postsmooth,
       .smoothOmega                 = p->smoothOmega };

  multigridBuild(mg, &spec);

  /* The hierarchy may support fewer levels than the setup asked for. */
  s->levels   = mg->levels;

  /* Solver.surface means the finest grid's surface list whatever the variant,
   * so anything outside the solver can ask about the body without knowing which
   * one is linked. Multigrid keeps a list per level; the finest of them is that
   * list, shared rather than copied. */
  s->surface  = *multigridFinestSurface(mg);

  s->mgLevels = mg;

  if (commIsMaster(s->comm)) {
    printf("Using Multigrid solver with %d levels\n", s->levels);
  }
}

double solve(Solver *s, double *p, const double *rhs)
{
  MultigridType *mg = (MultigridType *)s->mgLevels;

  double epssq      = s->eps * s->eps;
  int cycles        = 0;

  TIMESTART

  /* Multigrid used to run exactly one V-cycle per call, reading eps and itermax
   * into locals it never used, so there was no tolerance for it to converge to
   * and no way for it to agree with the relaxation solvers. */
  PressureLevelType desc;
  pressureLevelFromSolver(s, &desc);

  double res = pressureResidualNorm(&desc, p, rhs);

  while ((res >= epssq) && (cycles < s->itermax)) {
    multigridCycle(mg, p, rhs);

    res = pressureResidualNorm(&desc, p, rhs);
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
