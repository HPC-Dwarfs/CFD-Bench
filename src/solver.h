/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __SOLVER_H_
#define __SOLVER_H_
#include "comm.h"
#include "discretization.h"
#include "grid.h"
#include "parameter.h"

typedef struct {
    /* geometry and grid information */
    Grid *grid;
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
    char *problem;
    int bcLeft, bcRight, bcBottom, bcTop, bcFront, bcBack;
    /* communication */
    double **r, **e;
    int levels, presmooth, postsmooth;
    Comm *comm;
} Solver;

extern double solve(Solver *, double *, double *);
extern void initSolver(Solver *, Discretization *, Parameter *);
#endif
