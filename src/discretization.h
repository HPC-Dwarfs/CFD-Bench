/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#ifndef __DISCRETIZATION_H_
#define __DISCRETIZATION_H_

#include "grid.h"
#include "parameter.h"
#include "comm.h"

enum BC { NOSLIP = 1, SLIP, OUTFLOW, PERIODIC };

typedef struct {
    /* geometry and grid information */
    Grid grid;
    /* arrays */
    double *p, *rhs;
    double *f, *g, *h;
    double *u, *v, *w;
    /* parameters */
    double eps, omega;
    double re, tau, gamma;
    double gx, gy, gz;
    /* time stepping */
    int itermax;
    double dt, te;
    double dtBound;
    char* problem;
    int bcLeft, bcRight, bcBottom, bcTop, bcFront, bcBack;
    Comm comm;
} Discretization;

extern void initDiscretization(Discretization*, Parameter*);
extern void computeRHS(Discretization*);
extern void normalizePressure(Discretization*);
extern void computeTimestep(Discretization*);
extern void setBoundaryConditions(Discretization*);
extern void setSpecialBoundaryCondition(Discretization*);
extern void computeFG(Discretization*);
extern void adaptUV(Discretization*);
#endif
