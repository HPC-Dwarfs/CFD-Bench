/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __SOLVER_H_
#define __SOLVER_H_
#include <stdbool.h>

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
  /* the multigrid hierarchy, a MultigridType built by whichever solver wants
   * one. Opaque here: multigrid.h owns the type, and a build that never asks
   * for a hierarchy leaves this null. */
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

/*
 * Whether res is a finite number. Out of line and in its own object, compiled
 * without the fast-math assumption that nothing ever is anything else; an
 * inline isfinite() here would be compiled under that assumption and could be
 * folded to true. See finite.c.
 */
extern bool solveResidualIsFinite(double res);

/*
 * Whether an iteration goes on: the residual is finite, still at or above the
 * tolerance, and the budget is not spent.
 *
 * The finiteness test is not optional. "Continue while res >= epssq" alone
 * reads a NaN as the tolerance having been met, because a comparison against a
 * NaN is false whichever way it is written, so a diverged iteration used to
 * leave its loop and be reported as converged. An infinite residual that is not
 * a NaN would instead have kept iterating to the limit.
 *
 * The residual every solver tests is formed by a global reduction, so a NaN on
 * one rank reaches all of them and every rank leaves the loop together.
 */
static inline bool solveContinues(double res, double epssq, int it, int itermax)
{
  return solveResidualIsFinite(res) && (res >= epssq) && (it < itermax);
}

/*
 * Report how a solve stopped -- converged, stopped at the iteration limit, or
 * diverged -- and return the residual the caller should see.
 *
 * loopRes is the residual the loop tested, res the one computed from the field
 * being returned; unit names what the solver counts ("cycles", "iterations").
 * The converged message is the one every solver printed before, unchanged. A
 * divergence is reported whether or not the build is verbose, and the value
 * returned for it is never finite, so a caller can tell it from a converged
 * solve by the residual alone. Stopping at the limit returns res as it is: the
 * iterate is still the caller's to use.
 */
extern double solveReport(const Solver *s,
    const char *solver,
    const char *unit,
    int it,
    double loopRes,
    double res);


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
