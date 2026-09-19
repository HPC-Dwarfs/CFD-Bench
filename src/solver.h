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
  /* The cells the geometry-free sweep gets wrong, on the finest grid. Every
   * variant fills this, including the multigrid one, which also keeps a list
   * per coarse level. */
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
  /* The conjugate gradient solver's iteration vectors and its preconditioner.
   * Opaque outside solver-cg.c, the way mgLevels is outside solver-mg.c: four
   * more named field pointers would be carried by three variants that never
   * touch them. */
  void *cgState;
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

#if defined(TEST) && defined(SOLVER_cg)
/*
 * Exposed for the check drivers only. The preconditioner and the iteration
 * vectors are otherwise internal to solver-cg.c; a driver needs the
 * preconditioner as a bare linear operator to check its linearity and
 * symmetry, the iteration count to compare a preconditioned solve against an
 * unpreconditioned one, and a way to stop the recurrence part-way so that the
 * invariants the spec states about an intermediate iterate can be checked at
 * one.
 */
extern void cgTestPrecon(Solver *s, const double *r, double *z);
extern int cgTestIterations(Solver *s);
extern double *cgTestResidual(Solver *s);
extern double *cgTestDirection(Solver *s);
extern int cgTestIsSingular(Solver *s);
/* The same solve with the iteration limit lowered to steps, so an intermediate
 * iterate comes out of the real recurrence rather than a reimplementation. */
extern double cgTestSolveSteps(Solver *s, double *p, const double *rhs, int steps);
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

/*
 * A preconditioner, as a Krylov solver sees it.
 *
 * apply reads r and writes z outright, with no tolerance and no initial guess
 * to pass in, which is the linearity requirement made structural rather than
 * documented: there is nowhere for an inner convergence test or a warm start to
 * hide. It must also leave solid cells at exactly zero.
 *
 * The indirection exists so the multigrid preconditioner can be added without
 * touching solver-cg.c -- that change supplies a ctx holding a level hierarchy
 * and an apply that runs one symmetric V-cycle from zero.
 */
typedef struct {
  void (*apply)(void *ctx, const PressureLevelType *lv, const double *r, double *z);
  void *ctx;
} PreconType;

/* Mean square residual over the fluid unknowns, for comparison with eps*eps. */
extern double pressureResidualNorm(
    const PressureLevelType *lv, double *p, const double *rhs);

/*
 * Apply the operator to x, writing the result into y, matrix-free.
 *
 * The sign is negated against the operator as assembled above: the assembled
 * operator is negative semi-definite, and a Krylov method needs a positive
 * definite one, so this computes
 *
 *   y_c = -sum_f A_f (x_nb - x_c) / h^2   for a fluid cell,   0 for a solid one
 *
 * A solver using it therefore also negates its right-hand side, b_c =
 * -Lambda_c * rhs_c. Negating both sides leaves the solution and the magnitude
 * of the residual unchanged, so r.r is still exactly the numerator
 * pressureResidualNorm computes and the same eps means the same thing.
 *
 * x is exchanged and has the boundary condition applied to it, so it is not
 * const; y is written in full over the interior, halos untouched.
 *
 * Split the way the relaxation sweep is: a geometry-free 7-point bulk pass,
 * then a correction over the surface list alone.
 */
extern void pressureApplyOperator(const PressureLevelType *lv, double *x, double *y);

/* The surface half of the application: overwrite the listed cells with the
 * value their real coefficients imply, and solid cells with exactly zero.
 * Sibling of pressureCorrectSurface -- that one relaxes, this one applies -- so
 * the two carry the same coefficients and must stay in step. Uncoloured: an
 * apply reads only x, so traversal order cannot matter. */
extern void pressureApplySurface(
    const PressureLevelType *lv, const double *x, double *y);

/* Keep the listed cells of one colour before the bulk sweep overwrites them. */
extern void pressureSaveSurface(const PressureLevelType *lv, const double *p, int color);

/* Relax the listed cells of one colour with their real coefficients, and swap
 * their wrong contribution to sweepRes for the right one. */
extern void pressureCorrectSurface(const PressureLevelType *lv,
    double *p,
    const double *rhs,
    int color,
    double omega,
    double *sweepRes);

#ifdef PROFILING
/* Seconds spent in the interior sweep and in the cut-cell correction since the
 * last reset, for the timing check. */
extern void pressureSweepTimes(double *bulk, double *surface);
extern void pressureSweepTimesReset(void);
#endif

/* Fill a level description from the solver's finest grid. */
extern void pressureLevelFromSolver(const Solver *s, PressureLevelType *lv);
#endif
