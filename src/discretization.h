/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __DISCRETIZATION_H_
#define __DISCRETIZATION_H_

#include "comm.h"
#include "grid.h"
#include "parameter.h"
#include "pressure-bc.h"

enum BC { NOSLIP = 1, SLIP, OUTFLOW, PERIODIC };

typedef struct {
  /* geometry and grid information */
  Grid grid;
  /* arrays */
  double *p, *rhs;
  double *f, *g, *h;
  double *u, *v, *w;
  /* Obstacle geometry, as the rest of the solver sees it. Ax, Ay and Az are
   * face apertures co-located with u, v and w; Lambda is the cell volume
   * fraction, co-located with p. 0 is fully solid, 1 fully fluid. These four
   * arrays are the only channel between geometry and everything else -- nothing
   * downstream reads a file, a shape or a cell type. */
  double *Ax, *Ay, *Az, *Lambda;
  /* parameters */
  double eps, omega;
  double re, tau, gamma;
  double gx, gy, gz;
  /* time stepping */
  int itermax;
  double dt, te;
  double dtBound;
  char *problem;
  int bcLeft, bcRight, bcBottom, bcTop, bcFront, bcBack;
  /* the same boundary configuration the solvers use, so that the null-space
   * handling and the solvers agree on whether the operator is singular */
  PressureBcType pressureBc;
  CommType comm;
} Discretization;

extern void initDiscretization(Discretization *, Parameter *);
extern void computeRHS(Discretization *);
extern void normalizePressure(Discretization *);
extern void computeTimestep(Discretization *);
extern void setBoundaryConditions(Discretization *);
extern void setSpecialBoundaryCondition(Discretization *);
extern void computeFG(Discretization *);
extern void adaptUV(Discretization *);
#endif
