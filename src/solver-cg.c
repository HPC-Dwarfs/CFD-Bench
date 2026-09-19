/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Preconditioned conjugate gradients on the embedded-boundary pressure system.
 *
 * The system as assembled is negative semi-definite, so this solver works with
 * both sides negated -- see pressureApplyOperator in solver.h. That leaves the
 * solution and the magnitude of the residual untouched, so r.r over the fluid
 * cells divided by the global fluid cell count is exactly what
 * pressureResidualNorm returns and the same eps means the same thing here as it
 * does for the relaxation solvers.
 *
 * Solid cells are identity rows with a zero right-hand side, so the system is
 * block diagonal and the two blocks have opposite definiteness. The iteration
 * is therefore confined to the fluid block: every vector it carries holds
 * exactly zero in the solid, every inner product skips those cells, and the
 * preconditioner writes zero there. That is a correctness statement rather than
 * an optimization -- a Krylov method cannot self-heal the way a colour sweep
 * can, so a solid cell that drifts off zero poisons every subsequent step
 * length.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "comm.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "surface-list.h"
#include "timing.h"
#include "util.h"

/*
 * The iteration's own storage, hung off Solver as an opaque pointer the way
 * the multigrid hierarchy is: four more named field pointers would be carried
 * by three variants that never look at them.
 *
 *   r   residual, b - A_neg x
 *   z   preconditioned residual, M r
 *   d   search direction
 *   q   A_neg d
 *   b   the negated right-hand side, -Lambda*rhs, projected once where the
 *       operator is singular. Kept rather than recomputed because the
 *       projection must be applied to it exactly once per solve.
 */
typedef struct {
  double *r, *z, *d, *q, *b;
  PreconType precon;
  /* The reciprocal of the diagonal an interior cell has, which is the same for
   * every one of them. A cut cell's own reciprocal is on the surface list. */
  double invDiagBulk;
  /* True where the boundary configuration leaves the operator singular, so the
   * constant over the fluid has to be projected out. Read once. */
  int singular;
  /* Iterations the last solve took, for the test seam. */
  int iterations;
} CgStateType;

#define CGSTATE(s) ((CgStateType *)(s)->cgState)

/* -------------------------------------------------------------------------
 * Preconditioners
 *
 * Both are fixed linear operators applied in the same split the operator
 * itself uses: a bulk pass that reads no geometry, then an override confined to
 * the surface list. Neither has a tolerance, an iteration count or a starting
 * guess, so z depends linearly on r by construction.
 *
 * Solid cells: a listed one is overridden with exactly zero. One in the body's
 * deep interior is not on the list and is not touched by the override, but r is
 * exactly zero there -- nothing in the iteration ever writes a solid cell -- so
 * a scalar multiple of it is exactly zero too.
 * ------------------------------------------------------------------------- */

/* z = r. The unpreconditioned case, spelled out rather than special-cased in
 * the recurrence, so that the two paths differ only in the operator applied. */
static void preconNone(void *ctx, const PressureLevelType *lv, const double *r, double *z)
{
  (void)ctx;

  SurfaceListType *list = lv->list;
  size_t size = (size_t)(lv->imaxLocal + 2) * (lv->jmaxLocal + 2) * (lv->kmaxLocal + 2);

#ifdef PROFILING
  double start = getTimeStamp();
#endif

  for (size_t i = 0; i < size; i++) {
    z[i] = r[i];
  }

  if (list != NULL) {
    for (int e = 0; e < list->count; e++) {
      if (list->solid[e]) {
        z[list->index[e]] = 0.0;
      }
    }
  }

#ifdef PROFILING
  T[PRECON] += getTimeStamp() - start;
  C[PRECON]++;
#endif
}

/*
 * z = D^-1 r, with D the operator's diagonal, sum_f A_f / h^2.
 *
 * For an interior cell that is the constant 2*(idx2 + idy2 + idz2); for a cut
 * cell the surface list already carries its reciprocal as invDiag, which is 0
 * for a solid cell. So the preconditioner needs no stored field at all -- a
 * scalar multiply over the bulk and an override over the list -- and reads no
 * geometry array, exactly like the operator application it accompanies.
 *
 * Diagonal, hence trivially symmetric, which is the point: a convergence
 * failure in this change is attributable to the conjugate gradient recurrence
 * rather than to the preconditioner.
 */
static void preconJacobi(
    void *ctx, const PressureLevelType *lv, const double *r, double *z)
{
  CgStateType *cg       = (CgStateType *)ctx;
  SurfaceListType *list = lv->list;

  int imaxLocal         = lv->imaxLocal;
  int jmaxLocal         = lv->jmaxLocal;
  int kmaxLocal         = lv->kmaxLocal;
  double scale          = cg->invDiagBulk;

#ifdef PROFILING
  double start = getTimeStamp();
#endif

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        size_t idx = (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                     (size_t)j * (imaxLocal + 2) + (size_t)i;
        z[idx] = scale * r[idx];
      }
    }
  }

  if (list != NULL) {
    for (int e = 0; e < list->count; e++) {
      int idx = list->index[e];
      z[idx]  = list->invDiag[e] * r[idx];
    }
  }

#ifdef PROFILING
  T[PRECON] += getTimeStamp() - start;
  C[PRECON]++;
#endif
}

/* -------------------------------------------------------------------------
 * Fluid-restricted reductions
 * ------------------------------------------------------------------------- */

/*
 * a.b over the cells with a nonzero volume fraction.
 *
 * Every vector the iteration carries is exactly zero in the solid, so skipping
 * those cells changes no sum. It is done anyway, and checked, because it is the
 * statement that the iteration lives on the fluid block: if a solid entry ever
 * became nonzero the dot would be the first thing to hide it.
 *
 * Local only -- the caller reduces, so that two dots taken in the same
 * iteration can be fused into one allreduce.
 */
static double dotLocal(const PressureLevelType *lv, const double *a, const double *b)
{
  int imaxLocal        = lv->imaxLocal;
  int jmaxLocal        = lv->jmaxLocal;
  int kmaxLocal        = lv->kmaxLocal;
  const double *Lambda = lv->Lambda;
  double sum           = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        size_t idx = (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                     (size_t)j * (imaxLocal + 2) + (size_t)i;
        sum += a[idx] * b[idx];
      }
    }
  }

  return sum;
}

/*
 * Remove the component along the operator's null vector.
 *
 * Where every boundary imposes a zero normal gradient the operator is singular
 * and its null space is spanned by the vector that is one over the fluid and
 * zero over the solid -- not by the constant over every cell, which is what
 * makes the averaging fluid-only rather than a matter of taste.
 *
 * Applied to the right-hand side and the initial guess once each, and to the
 * preconditioned residual every iteration, because a preconditioner is free to
 * inject a constant and the recurrence will then amplify it.
 *
 * This is the one place in the iteration that reads Lambda over the whole
 * domain. The sum does not need it -- solid entries are zero and contribute
 * nothing -- but the subtraction does, because a constant subtracted from a
 * solid cell would leave it holding that constant and nothing else would take
 * it back off. It runs only in singular setups.
 */
static void projectOutConstant(const PressureLevelType *lv, double *v)
{
  int imaxLocal        = lv->imaxLocal;
  int jmaxLocal        = lv->jmaxLocal;
  int kmaxLocal        = lv->kmaxLocal;
  const double *Lambda = lv->Lambda;
  double mean          = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        size_t idx = (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                     (size_t)j * (imaxLocal + 2) + (size_t)i;
        mean += v[idx];
      }
    }
  }

  commReduceAll(&mean, SUM);
  mean /= lv->fluidCells;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        size_t idx = (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                     (size_t)j * (imaxLocal + 2) + (size_t)i;
        v[idx] -= mean;
      }
    }
  }
}

/* -------------------------------------------------------------------------
 * Setup
 * ------------------------------------------------------------------------- */

void initSolver(Solver *s, Discretization *d, Parameter *p)
{
  solverBaseInit(s, d, p);

  double dx2       = s->grid->dx * s->grid->dx;
  double dy2       = s->grid->dy * s->grid->dy;
  double dz2       = s->grid->dz * s->grid->dz;

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

  CgStateType *cg  = malloc(sizeof(CgStateType));

  size_t size      = (size_t)(s->comm->imaxLocal + 2) * (s->comm->jmaxLocal + 2) *
                (s->comm->kmaxLocal + 2);

  cg->r = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
  cg->z = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
  cg->d = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
  cg->q = allocate(ARRAY_ALIGNMENT, size * sizeof(double));
  cg->b = allocate(ARRAY_ALIGNMENT, size * sizeof(double));

  for (size_t i = 0; i < size; i++) {
    cg->r[i] = 0.0;
    cg->z[i] = 0.0;
    cg->d[i] = 0.0;
    cg->q[i] = 0.0;
    cg->b[i] = 0.0;
  }

  cg->invDiagBulk = 1.0 / (2.0 * (1.0 / dx2 + 1.0 / dy2 + 1.0 / dz2));
  cg->singular    = pressureBcIsSingular(&s->bc) ? 1 : 0;
  cg->iterations  = 0;

  const char *precon = (p->precon != NULL) ? p->precon : "jacobi";

  if (strcmp(precon, "none") == 0) {
    cg->precon.apply = preconNone;
    cg->precon.ctx   = cg;
  } else if (strcmp(precon, "jacobi") == 0) {
    cg->precon.apply = preconJacobi;
    cg->precon.ctx   = cg;
  } else {
    /* Refused here rather than silently falling back, the way pressureBcInit
     * refuses a boundary code: a run that quietly used a different
     * preconditioner from the one asked for would be a benchmark result
     * attributed to the wrong method. "mg" is named because it is the one the
     * follow-up change adds, and a parameter file written for that change would
     * otherwise run unpreconditioned here. */
    if (commIsMaster(s->comm)) {
      fprintf(stderr,
          "Unsupported preconditioner \"%s\". Supported: none, jacobi.%s\n",
          precon,
          strcmp(precon, "mg") == 0
              ? " The multigrid preconditioner needs a symmetric V-cycle and is"
                " not implemented yet."
              : "");
    }
    exit(EXIT_FAILURE);
  }

  s->cgState = cg;

  if (commIsMaster(s->comm)) {
    printf("Using Conjugate Gradient solver with the %s preconditioner, %.1f MB of "
           "extra field storage\n",
        precon,
        (double)(5 * size * sizeof(double)) * 1.0e-6);
  }
}

/* -------------------------------------------------------------------------
 * The iteration
 * ------------------------------------------------------------------------- */

double solve(Solver *s, double *p, const double *rhs)
{
  CgStateType *cg = CGSTATE(s);

  PressureLevelType lv;
  pressureLevelFromSolver(s, &lv);

  int imaxLocal        = lv.imaxLocal;
  int jmaxLocal        = lv.jmaxLocal;
  int kmaxLocal        = lv.kmaxLocal;
  const double *Lambda = lv.Lambda;

  double *r            = cg->r;
  double *z            = cg->z;
  double *dir          = cg->d;
  double *q            = cg->q;
  double *b            = cg->b;

  size_t size   = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);
  double epssq  = s->eps * s->eps;
  int itermax   = s->itermax;
  int it        = 0;

  TIMESTART

  /* b = -Lambda*rhs over the fluid, zero in the solid: the negated system's
   * right-hand side, and an identity row's zero. */
  for (size_t i = 0; i < size; i++) {
    b[i] = 0.0;
  }

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        size_t idx = (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) +
                     (size_t)j * (imaxLocal + 2) + (size_t)i;
        b[idx] = -LAM(i, j, k) * rhs[idx];
      }
    }
  }

  /* The warm start. The caller's p carries the previous step's pressure, which
   * is a good guess; what it must not carry is anything in the solid, which
   * would enter the very first inner product. Listed solid cells are cleared
   * here; a deep interior one is already zero and stays so. */
  if (lv.list != NULL) {
    for (int e = 0; e < lv.list->count; e++) {
      if (lv.list->solid[e]) {
        p[lv.list->index[e]] = 0.0;
      }
    }
  }

  if (cg->singular) {
    /* Once each: the right-hand side, so the singular system is consistent and
     * a solution exists at all, and the initial guess, so the warm start
     * carries no constant of its own. */
    projectOutConstant(&lv, b);
    projectOutConstant(&lv, p);
  }

  /* r = b - A_neg p */
  pressureApplyOperator(&lv, p, q);

  for (size_t i = 0; i < size; i++) {
    r[i] = b[i] - q[i];
  }

  cg->precon.apply(cg->precon.ctx, &lv, r, z);

  if (cg->singular) {
    projectOutConstant(&lv, z);
  }

  /* One fused reduction for the stopping norm and the first rho. */
  double dots[2];
  dots[0] = dotLocal(&lv, r, r);
  dots[1] = dotLocal(&lv, r, z);
  { PROFILE(CG_DOT, commReduceAllN(dots, 2, SUM)); }

  double rr  = dots[0];
  double rho = dots[1];

  for (size_t i = 0; i < size; i++) {
    dir[i] = z[i];
  }

  while ((rr / lv.fluidCells >= epssq) && (it < itermax)) {
    pressureApplyOperator(&lv, dir, q);

    double dq = dotLocal(&lv, dir, q);
    { PROFILE(CG_DOT, commReduceAllN(&dq, 1, SUM)); }

    /* A breakdown rather than a slow run: the operator is positive definite on
     * the fluid subspace, so d.q is zero only if d is, which means the residual
     * already is. Stopping is the honest response; continuing would divide by
     * zero. */
    if (dq <= 0.0) {
      break;
    }

    double alpha = rho / dq;

#ifdef PROFILING
    double axpyStart = getTimeStamp();
#endif

    for (size_t i = 0; i < size; i++) {
      p[i] += alpha * dir[i];
      r[i] -= alpha * q[i];
    }

#ifdef PROFILING
    T[CG_AXPY] += getTimeStamp() - axpyStart;
    C[CG_AXPY]++;
#endif

    cg->precon.apply(cg->precon.ctx, &lv, r, z);

    if (cg->singular) {
      projectOutConstant(&lv, z);
    }

    /* The stopping norm and the next rho, fused: commReduceAllN takes a count,
     * so the norm costs no latency beyond the r.z the algorithm needs anyway.
     * Two allreduces an iteration is what a Krylov method costs and the
     * relaxation solvers do not -- it gets its own profiler region rather than
     * being lumped into the sweep. */
    dots[0] = dotLocal(&lv, r, r);
    dots[1] = dotLocal(&lv, r, z);
    { PROFILE(CG_DOT, commReduceAllN(dots, 2, SUM)); }

    double rhoNew = dots[1];
    double beta   = rhoNew / rho;

    rr            = dots[0];
    rho           = rhoNew;

#ifdef PROFILING
    double dirStart = getTimeStamp();
#endif

    for (size_t i = 0; i < size; i++) {
      dir[i] = z[i] + beta * dir[i];
    }

#ifdef PROFILING
    T[CG_AXPY] += getTimeStamp() - dirStart;
    C[CG_AXPY]++;
#endif

    it++;

#ifdef DEBUG
    if (commIsMaster(s->comm)) {
      printf("%d Residuum: %e\n", it, rr / lv.fluidCells);
    }
#endif
  }
  TIMESTOP(SOLVER);

  cg->iterations = it;

  /* Reported from a dedicated pass, the way every other solver reports it, so
   * the number in residual.dat is the residual of the field actually returned
   * rather than of the one the recurrence believed it had. */
  double res     = pressureResidualNorm(&lv, p, rhs);

#ifdef VERBOSE
  if (commIsMaster(s->comm)) {
    printf("Conjugate Gradient took %d iterations to reach %e\n", it, sqrt(res));
  }

  printProfile(s->comm, it);
#endif

  return res;
}

#if defined(TEST) && defined(SOLVER_cg)
/*
 * Exposed for the check drivers only. A driver needs to reach the
 * preconditioner as a bare linear operator, to read the iteration count the
 * last solve took, and to stop the recurrence part-way through so that the
 * invariants the spec states about an intermediate iterate can be checked at
 * one. Same pattern as the mgTest* entry points.
 */
void cgTestPrecon(Solver *s, const double *r, double *z)
{
  CgStateType *cg = CGSTATE(s);

  PressureLevelType lv;
  pressureLevelFromSolver(s, &lv);

  cg->precon.apply(cg->precon.ctx, &lv, r, z);
}

int cgTestIterations(Solver *s) { return CGSTATE(s)->iterations; }

double *cgTestResidual(Solver *s) { return CGSTATE(s)->r; }

double *cgTestDirection(Solver *s) { return CGSTATE(s)->d; }

int cgTestIsSingular(Solver *s) { return CGSTATE(s)->singular; }

/* Solve with the iteration limit temporarily lowered, which is how a driver
 * reaches an intermediate iterate: the recurrence is the same one, stopped
 * early rather than reimplemented. */
double cgTestSolveSteps(Solver *s, double *p, const double *rhs, int steps)
{
  int saved  = s->itermax;
  s->itermax = steps;

  double res = solve(s, p, rhs);

  s->itermax = saved;
  return res;
}
#endif
