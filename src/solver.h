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
#include "pressure-bc.h"
#include "surface-list.h"

typedef struct {
  /* geometry and grid information */
  Grid *grid;
  /* arrays */
  double *p, *rhs;
  double *f, *g, *h;
  double *u, *v, *w;
  /* obstacle geometry, owned by the discretization */
  const double *Ax, *Ay, *Az, *Lambda;
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
  /* the pressure boundary condition the setup asks for, shared by every solver */
  PressureBcType bc;
  /* global index of this rank's first interior cell, minus one, so that the
   * red-black colouring can be taken from a cell's global position */
  int iOffset, jOffset, kOffset;
  /* the cells the geometry-free sweep gets wrong, on the finest grid */
  SurfaceListType surface;
  /* global fluid cell count, which residual norms are divided by */
  double fluidCells;
  /* communication */
  double **r, **e;
  int levels, presmooth, postsmooth;
  /* the multigrid level hierarchy: extents, mesh, communicator and work arrays
   * per level. Opaque outside solver-mg.c, which is the only file that builds
   * or reads it. */
  void *mgLevels;
  CommType *comm;
} Solver;

extern double solve(Solver *, double *, const double *);
extern void initSolver(Solver *, Discretization *, Parameter *);

/* Everything the three solvers set up identically. Each initSolver calls this
 * first and then adds whatever only it needs. */
extern void solverBaseInit(Solver *s, Discretization *d, Parameter *p);

#if defined(TEST) && defined(SOLVER_mg)
/*
 * Exposed for the check drivers only. The multigrid transfer operators and the
 * V-cycle are otherwise internal to solver-mg.c; a driver needs to drive them
 * one step at a time to check a constant restricts to a constant, that
 * prolongation leaves no fine cell untouched, and that a coarse level uses its
 * own mesh.
 */
extern int mgTestLevels(Solver *s);
extern void mgTestLevelExtents(Solver *s, int level, int *im, int *jm, int *km);
extern void mgTestLevelMesh(Solver *s, int level, double *dx, double *dy, double *dz);
extern double *mgTestLevelE(Solver *s, int level);
extern double *mgTestLevelR(Solver *s, int level);
extern void mgTestRestrict(Solver *s, int level);
extern void mgTestProlongate(Solver *s, int level);
extern void mgTestResidualField(Solver *s, int level, double *p, const double *rhs);
extern void mgTestVcycle(Solver *s, double *p, const double *rhs);
extern void mgTestSmooth(Solver *s, int level, double *p, const double *rhs, int sweeps);
/* Solid cells and surface-list length at one level, for checking that the
 * coarsened geometry still represents the body. */
extern int mgTestLevelSolidCount(Solver *s, int level);
extern int mgTestLevelSurfaceCount(Solver *s, int level);
#endif

/*
 * Everything the pressure operator needs on one grid. Multigrid builds one per
 * level; the relaxation solvers build a single one for the finest.
 */
typedef struct {
  CommType *comm;
  const PressureBcType *bc;
  const double *Ax, *Ay, *Az, *Lambda;
  SurfaceListType *list;
  int imaxLocal, jmaxLocal, kmaxLocal;
  double dx, dy, dz;
  double fluidCells;
} PressureLevelType;

/* Mean square residual over the fluid unknowns, for comparison with eps*eps. */
extern double pressureResidualNorm(
    const PressureLevelType *lv, double *p, const double *rhs);

/* Keep the listed cells of one colour before the bulk sweep overwrites them. */
extern void pressureSaveSurface(
    const PressureLevelType *lv, const double *p, int color);

/* Relax the listed cells of one colour with their real coefficients, and swap
 * their wrong contribution to sweepRes for the right one. */
extern void pressureCorrectSurface(const PressureLevelType *lv,
    double *p,
    const double *rhs,
    int color,
    double omega,
    double *sweepRes);

/* Fill a level description from the solver's finest grid. */
extern void pressureLevelFromSolver(const Solver *s, PressureLevelType *lv);
#endif
