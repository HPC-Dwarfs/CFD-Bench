/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The multigrid level hierarchy and the cycle that runs on it.
 *
 * Deliberately not named solver-*.c: the Makefile treats those as the mutually
 * exclusive solver variants and links exactly one of them, so a SOLVER=cg build
 * would otherwise contain no multigrid code at all and could not call a cycle
 * even to precondition with. This is the same reason solverbase.c carries the
 * operator.
 *
 * Nothing here knows what a Solver is. solver-mg.c is the driver that iterates
 * cycles to a tolerance; precon-mg.c applies a single cycle as a preconditioner.
 * Both build a hierarchy from the description below and neither reaches inside
 * it.
 */
#ifndef __MULTIGRID_H_
#define __MULTIGRID_H_

#include "comm.h"
#include "grid.h"
#include "pressure-bc.h"
#include "solver.h"
#include "surface-list.h"

/*
 * One level of the hierarchy. Each carries its own extents, mesh size, global
 * offsets and communicator, because the operator, the colouring and the
 * boundary condition all depend on them and all three used to be taken from
 * the finest level regardless of which level was being relaxed.
 *
 * e and r are this level's error and residual. The finest level relaxes the
 * caller's pressure against the caller's right-hand side and uses only r; every
 * coarser level relaxes its own e against its own r.
 */
typedef struct {
  CommType comm;
  int imaxLocal, jmaxLocal, kmaxLocal;
  int iOffset, jOffset, kOffset;
  double dx, dy, dz;
  double cells;      /* global interior cells at this level */
  double fluidCells; /* global fluid cells at this level */
  /* The geometry, coarsened from the level above: a coarse face aperture is the
   * mean of the four fine faces it covers, a coarse volume fraction the mean of
   * the eight fine cells. Level 0 borrows the caller's arrays. */
  double *Ax, *Ay, *Az, *Lambda;
  int ownsGeometry;
  SurfaceListType surface;
  double *e, *r;
  /* The right-hand side this level solves against, restricted from the level
   * above. Kept apart from r, which residualField overwrites with this level's
   * own residual partway through the cycle -- when the two were one array, the
   * post-smoother on every intermediate level relaxed against the residual
   * instead of against the right-hand side it was given. */
  double *b;
  /* The correction prolongated from the level below. Kept apart from e, which
   * is the solution this level is computing: at any level but the finest the
   * recursion passes e in as p, so when prolongation wrote e it overwrote the
   * very iterate it was about to be added to, and correct() doubled the iterate
   * instead of correcting it. */
  double *corr;
  /* Scratch for restriction, which blends on the fine grid before summing into
   * the coarse one and must not destroy the residual it is given. */
  double *scratch;
} MgLevelType;

/*
 * Everything a hierarchy needs to exist, and nothing about who wants one.
 *
 * The finest level borrows comm, the geometry arrays and the grid rather than
 * copying them, so they must outlive the hierarchy. Every coarser level owns
 * its own.
 */
typedef struct {
  CommType *comm;
  const PressureBcType *bc;
  const Grid *grid;
  const double *Ax, *Ay, *Az, *Lambda;
  int iOffset, jOffset, kOffset;
  double fluidCells;
  int levels;
  int presmooth, postsmooth;
  /* The smoother's relaxation factor. Deliberately not the SOR solvers' omg:
   * see the smoothing note in the change's design. */
  double smoothOmega;
} MultigridSpecType;

/*
 * A built hierarchy. levels may be fewer than the spec asked for, when the
 * decomposition cannot be coarsened that far; multigridBuild reports that and
 * records what it actually built.
 */
typedef struct {
  const PressureBcType *bc;
  double smoothOmega;
  int levels;
  int presmooth, postsmooth;
  MgLevelType *level;
} MultigridType;

/* Build the hierarchy described by spec. Reports on the master rank when it
 * supports fewer levels than asked, and when the body stops being represented
 * on a coarse level. */
extern void multigridBuild(MultigridType *mg, const MultigridSpecType *spec);

extern void multigridFree(MultigridType *mg);

/* One cycle on the finest level: p is relaxed against rhs and updated in
 * place. */
extern void multigridCycle(MultigridType *mg, double *p, const double *rhs);

/* The finest level's surface list, which a solver publishes as its own so that
 * anything outside can ask about the body without knowing which variant is
 * linked. Shared, not copied. */
extern SurfaceListType *multigridFinestSurface(MultigridType *mg);

/* The cycle as a preconditioner: one cycle from a zero guess. Lives in
 * precon-mg.c, which every build links. */
extern void preconMgInit(PreconType *precon, MultigridType *mg);

#if defined(TEST)
/*
 * Exposed for the check drivers only. The transfer operators and the cycle are
 * otherwise internal; a driver needs to drive them one step at a time to check
 * that a constant restricts to a constant, that prolongation leaves no fine
 * cell untouched, that restriction and prolongation are a transpose pair, and
 * that a coarse level uses its own mesh.
 */
extern int mgTestLevels(MultigridType *mg);
extern void mgTestLevelExtents(MultigridType *mg, int level, int *im, int *jm, int *km);
extern void mgTestLevelMesh(
    MultigridType *mg, int level, double *dx, double *dy, double *dz);
extern double *mgTestLevelE(MultigridType *mg, int level);
extern double *mgTestLevelR(MultigridType *mg, int level);
/* The right-hand side restriction writes, and the correction prolongation
 * writes. Separate arrays from e and r, which they used to share. */
extern double *mgTestLevelB(MultigridType *mg, int level);
extern double *mgTestLevelCorr(MultigridType *mg, int level);
extern void mgTestRestrict(MultigridType *mg, int level);
extern void mgTestProlongate(MultigridType *mg, int level);
extern void mgTestResidualField(
    MultigridType *mg, int level, double *p, const double *rhs);
extern void mgTestVcycle(MultigridType *mg, double *p, const double *rhs);
extern void mgTestSmooth(
    MultigridType *mg, int level, double *p, const double *rhs, int sweeps);
/* Solid cells and surface-list length at one level, for checking that the
 * coarsened geometry still represents the body. */
extern int mgTestLevelSolidCount(MultigridType *mg, int level);
extern int mgTestLevelSurfaceCount(MultigridType *mg, int level);
#endif

#endif // __MULTIGRID_H_
