/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The conjugate gradient solver, and the properties it needs that the
 * relaxation solvers never did.
 *
 *   - both preconditioners are linear and symmetric, which is what makes the
 *     preconditioned recurrence a conjugate gradient method at all rather than
 *     a plausible-looking iteration;
 *   - a fluid-only dot really is fluid-only: over a domain with a body, the
 *     unit vector dotted with itself is the fluid cell count;
 *   - convergence to the exact discrete solution on the Dirichlet and Neumann
 *     Poisson configurations, with a body and without one;
 *   - solid cells hold exactly zero at an intermediate iteration, not only at
 *     convergence, and adding solid volume away from the fluid changes neither
 *     the iteration count nor the fluid field;
 *   - the preconditioned solve takes no more iterations than the
 *     unpreconditioned one, so the preconditioner is demonstrably working
 *     rather than merely being applied.
 *
 * Reaches solver-cg.c's internals through the cgTest* seam, so it only builds
 * under SOLVER=cg and skips cleanly under the others.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "discretization.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "surface-list.h"
#include "util.h"

#if !defined(SOLVER_cg)

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("cg");
  printf("cg: skipped, this driver is linked against another solver\n");
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, CheckFailures);

  commFinalize(&comm);
  return EXIT_SUCCESS;
}

#else

#define IMAX 32
#define JMAX 24
#define KMAX 16

static Parameter Params;
static Discretization D;
static Solver S;

static void setup(CommType *base, int boundary, const char *geometry, const char *precon)
{
  initParameter(&Params);

  Params.name         = "check-cg";
  Params.imax         = IMAX;
  Params.jmax         = JMAX;
  Params.kmax         = KMAX;
  Params.xlength      = 2.0;
  Params.ylength      = 1.5;
  Params.zlength      = 1.0;

  Params.eps          = 1e-11;
  Params.omg          = 1.7;
  Params.itermax      = 20000;
  Params.levels       = 3;
  Params.presmooth    = 4;
  Params.postsmooth   = 4;

  Params.re           = 100.0;
  Params.tau          = 0.5;
  Params.gamma        = 0.9;
  Params.dt           = 0.02;
  Params.te           = 0.0;
  Params.gx = Params.gy = Params.gz = 0.0;
  Params.u_init = Params.v_init = Params.w_init = Params.p_init = 0.0;
  Params.geometryFile = (char *)geometry;
  Params.precon       = (char *)precon;

  Params.bcLeft = Params.bcRight = boundary;
  Params.bcBottom = Params.bcTop = boundary;
  Params.bcFront = Params.bcBack = boundary;

  D.comm = *base;
  commPartition(&D.comm, KMAX, JMAX, IMAX);
  initDiscretization(&D, &Params);
  initSolver(&S, &D, &Params);
  initProfiler(&D.comm);
}

static size_t fieldSize(void)
{
  return (size_t)(D.comm.imaxLocal + 2) * (D.comm.jmaxLocal + 2) *
         (D.comm.kmaxLocal + 2);
}

static size_t linear(int i, int j, int k)
{
  int imaxLocal = D.comm.imaxLocal;
  int jmaxLocal = D.comm.jmaxLocal;
  return (size_t)k * (imaxLocal + 2) * (jmaxLocal + 2) + (size_t)j * (imaxLocal + 2) +
         (size_t)i;
}

static double dotFluid(const double *a, const double *b)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;
  double sum           = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        sum += a[linear(i, j, k)] * b[linear(i, j, k)];
      }
    }
  }

  commReduceAll(&sum, SUM);
  return sum;
}

/* Zero in the solid, and the same field however the domain is divided: seeded
 * from the global position rather than from rand(), so an assertion means the
 * same thing at one rank and at four. */
static void randomFluidField(double *x, unsigned seed)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;
  size_t size          = fieldSize();

  for (size_t i = 0; i < size; i++) {
    x[i] = 0.0;
  }

  int offsets[NDIMS] = { 0, 0, 0 };
  commGetOffsets(&D.comm, offsets, KMAX, JMAX, IMAX);

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }

        unsigned gi = (unsigned)(i - 1 + offsets[IDIM]);
        unsigned gj = (unsigned)(j - 1 + offsets[JDIM]);
        unsigned gk = (unsigned)(k - 1 + offsets[KDIM]);
        unsigned h  = seed + gi * 73856093u + gj * 19349663u + gk * 83492791u;

        h ^= h >> 13;
        h *= 1274126177u;
        h ^= h >> 16;

        x[linear(i, j, k)] = (double)h / (double)0xffffffffu - 0.5;
      }
    }
  }
}

/* -------------------------------------------------------------------------
 * The preconditioner as a bare linear operator
 * ------------------------------------------------------------------------- */

static void checkPreconditioner(const char *label)
{
  size_t size = fieldSize();

  double *x   = calloc(size, sizeof(double));
  double *y   = calloc(size, sizeof(double));
  double *sum = calloc(size, sizeof(double));
  double *mx  = calloc(size, sizeof(double));
  double *my  = calloc(size, sizeof(double));
  double *ms  = calloc(size, sizeof(double));

  randomFluidField(x, 13579u);
  randomFluidField(y, 24681u);

  const double a = 2.5;
  const double b = -0.75;

  for (size_t i = 0; i < size; i++) {
    sum[i] = a * x[i] + b * y[i];
  }

  cgTestPrecon(&S, x, mx);
  cgTestPrecon(&S, y, my);
  cgTestPrecon(&S, sum, ms);

  /* Linearity: M(ax + by) == a Mx + b My. A preconditioner with an inner
   * tolerance or a warm start fails this, which is the point of checking it
   * rather than the seam's shape. */
  double worst = 0.0;
  double scale = 0.0;

  for (size_t i = 0; i < size; i++) {
    double want = a * mx[i] + b * my[i];
    double diff = fabs(ms[i] - want);

    if (diff > worst) {
      worst = diff;
    }
    if (fabs(want) > scale) {
      scale = fabs(want);
    }
  }

  commReduceAll(&worst, MAX);
  commReduceAll(&scale, MAX);

  CHECK_TRUE(worst <= 1e-12 * scale + 1e-300,
      "%s: the preconditioner is not linear, worst deviation %.3e against a largest "
      "value of %.3e",
      label,
      worst,
      scale);

  /* Symmetry over the fluid unknowns: x.My == y.Mx. */
  double xMy    = dotFluid(x, my);
  double yMx    = dotFluid(y, mx);
  double dscale = fabs(xMy) + fabs(yMx) + 1e-300;

  CHECK_TRUE(fabs(xMy - yMx) / dscale < 1e-12,
      "%s: the preconditioner is not symmetric, x.My = %.17g against y.Mx = %.17g",
      label,
      xMy,
      yMx);

  /* Positive definite, so the preconditioned system is one a conjugate gradient
   * method can be run on at all. */
  double xMx = dotFluid(x, mx);

  CHECK_TRUE(xMx > 0.0,
      "%s: the preconditioner's quadratic form is %.17g, so it is not positive "
      "definite",
      label,
      xMx);

  /* And it leaves the solid alone. Fed a residual that is already zero there,
   * which is the invariant the iteration maintains, every solid cell must come
   * back exactly zero -- not nearly zero. */
  {
    int imaxLocal        = D.comm.imaxLocal;
    int jmaxLocal        = D.comm.jmaxLocal;
    int kmaxLocal        = D.comm.kmaxLocal;
    const double *Lambda = D.Lambda;
    int nonzero          = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) > 0.0) {
            continue;
          }
          if (mx[linear(i, j, k)] != 0.0) {
            ++nonzero;
          }
        }
      }
    }

    double total = (double)nonzero;
    commReduceAll(&total, SUM);

    CHECK_TRUE(total == 0.0,
        "%s: the preconditioner left %.0f solid cells nonzero",
        label,
        total);
  }

  free(x);
  free(y);
  free(sum);
  free(mx);
  free(my);
  free(ms);
}

/* -------------------------------------------------------------------------
 * Convergence against a known discrete solution
 *
 * The right-hand side is built by applying the operator to a chosen field, so
 * that field is the exact discrete solution and any difference afterwards is
 * the solver's rather than the mesh's. Same construction tests/checks/poisson.c
 * uses, run here through the conjugate gradient solver and through a body.
 * ------------------------------------------------------------------------- */

static double meanFluid(const double *f)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;
  double sum           = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        sum += f[linear(i, j, k)];
      }
    }
  }

  commReduceAll(&sum, SUM);
  return sum / S.fluidCells;
}

/* Fill exact with a smooth field, zero in the solid, and D.rhs with the
 * right-hand side that makes it the discrete solution. */
static void buildExactProblem(double *exact)
{
  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;

  int offsets[NDIMS]   = { 0, 0, 0 };
  commGetOffsets(&D.comm, offsets, KMAX, JMAX, IMAX);

  size_t size = fieldSize();
  for (size_t i = 0; i < size; i++) {
    exact[i] = 0.0;
  }

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }

        double x    = ((i - 1 + offsets[IDIM]) + 0.5) * D.grid.dx / Params.xlength;
        double y    = ((j - 1 + offsets[JDIM]) + 0.5) * D.grid.dy / Params.ylength;
        double z    = ((k - 1 + offsets[KDIM]) + 0.5) * D.grid.dz / Params.zlength;

        /* The smooth part is the field tests/checks/poisson.c uses. On its own
         * it is a poor test of a Krylov method here: with every boundary an
         * outflow it is precisely the lowest eigenvector of the discrete
         * operator, so conjugate gradients terminate in one iteration and
         * neither the convergence check nor the comparison of preconditioners
         * says anything. The rough part spreads the right-hand side across the
         * spectrum so the solver has to work for its answer. It is seeded from
         * the global position, so the problem is the same however the domain is
         * divided. */
        unsigned gi = (unsigned)(i - 1 + offsets[IDIM]);
        unsigned gj = (unsigned)(j - 1 + offsets[JDIM]);
        unsigned gk = (unsigned)(k - 1 + offsets[KDIM]);
        unsigned h  = 987654321u + gi * 73856093u + gj * 19349663u + gk * 83492791u;

        h ^= h >> 13;
        h *= 1274126177u;
        h ^= h >> 16;

        double rough = (double)h / (double)0xffffffffu - 0.5;

        exact[linear(i, j, k)] = sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) +
                                 0.1 * rough;
      }
    }
  }

  /* rhs is what the solvers take: A p = Lambda*rhs. The application computes
   * -A p, so Lambda*rhs = -apply(exact), and a fluid cell's rhs is that divided
   * by its volume fraction. */
  PressureLevelType lv;
  pressureLevelFromSolver(&S, &lv);

  double *ap = calloc(size, sizeof(double));
  pressureApplyOperator(&lv, exact, ap);

  double *rhs = D.rhs;
  for (size_t i = 0; i < size; i++) {
    rhs[i] = 0.0;
  }

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          continue;
        }
        RHS(i, j, k) = -ap[linear(i, j, k)] / LAM(i, j, k);
      }
    }
  }

  free(ap);
}

/* Returns the iteration count the solve took. */
static int runExactCase(int boundary, const char *label)
{
  size_t size   = fieldSize();
  double *exact = calloc(size, sizeof(double));

  buildExactProblem(exact);

  for (size_t i = 0; i < size; i++) {
    D.p[i] = 0.0;
  }

  double res = solve(&S, D.p, D.rhs);
  int iters  = cgTestIterations(&S);

  CHECK_TRUE(res < Params.eps * Params.eps * 100.0,
      "%s: the solve did not converge, mean square residual %.3e against eps^2 %.3e",
      label,
      res,
      Params.eps * Params.eps);

  /* The reported residual has to be the residual of the field returned, which
   * is what makes CG's number comparable with every other solver's. */
  PressureLevelType lv;
  pressureLevelFromSolver(&S, &lv);
  double independent = pressureResidualNorm(&lv, D.p, D.rhs);

  CHECK_TRUE(fabs(independent - res) <= 1e-12 * fabs(res) + 1e-30,
      "%s: the reported residual %.17g is not the residual of the returned field "
      "%.17g",
      label,
      res,
      independent);

  /* The all-Neumann operator fixes the solution only up to a constant over the
   * fluid, so both fields have their fluid mean removed before comparison. */
  double shiftP = 0.0, shiftE = 0.0;

  if (boundary == NOSLIP) {
    shiftP = meanFluid(D.p);
    shiftE = meanFluid(exact);
  }

  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;
  const double *Lambda = D.Lambda;
  double maxDiff       = 0.0;
  int nonzeroSolid     = 0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          if (D.p[linear(i, j, k)] != 0.0) {
            ++nonzeroSolid;
          }
          continue;
        }

        double diff = fabs((D.p[linear(i, j, k)] - shiftP) -
                           (exact[linear(i, j, k)] - shiftE));
        if (diff > maxDiff) {
          maxDiff = diff;
        }
      }
    }
  }

  commReduceAll(&maxDiff, MAX);

  double solidTotal = (double)nonzeroSolid;
  commReduceAll(&solidTotal, SUM);

  CHECK_TRUE(maxDiff < 1e-4,
      "%s: the converged field differs from the exact discrete solution by %.3e",
      label,
      maxDiff);
  CHECK_TRUE(solidTotal == 0.0,
      "%s: %.0f solid cells hold a nonzero pressure after the solve",
      label,
      solidTotal);

  if (commIsMaster(&D.comm)) {
    printf("%s: %d iterations, residual %.3e, max deviation from the exact solution "
           "%.3e\n",
        label,
        iters,
        sqrt(res),
        maxDiff);
  }

  free(exact);
  return iters;
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  /* commFinalize frees the derived datatypes, which only commPartition creates,
   * so the communicator this driver finalizes has to have been partitioned. */
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("cg");

  const char *sphere = "sphere:1.0,0.75,0.5,0.3";

  /* ---- The preconditioners, as bare linear operators ---- */

  setup(&comm, OUTFLOW, NULL, "none");
  checkPreconditioner("none, obstacle-free");

  setup(&comm, OUTFLOW, NULL, "jacobi");
  checkPreconditioner("jacobi, obstacle-free");

  setup(&comm, OUTFLOW, sphere, "none");
  checkPreconditioner("none, sphere");

  setup(&comm, OUTFLOW, sphere, "jacobi");
  checkPreconditioner("jacobi, sphere");

  /* ---- A fluid-only dot really is fluid-only ---- */
  {
    size_t size = fieldSize();
    double *one = calloc(size, sizeof(double));

    /* The unit vector over the fluid: its dot with itself is the fluid cell
     * count and nothing else, which is the statement that the reduction skips
     * the body rather than merely relying on it holding zeros. */
    int imaxLocal        = D.comm.imaxLocal;
    int jmaxLocal        = D.comm.jmaxLocal;
    int kmaxLocal        = D.comm.kmaxLocal;
    const double *Lambda = D.Lambda;

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          one[linear(i, j, k)] = 1.0;
        }
      }
    }

    double solidCells = 0.0;
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) == 0.0) {
            solidCells += 1.0;
          }
        }
      }
    }
    commReduceAll(&solidCells, SUM);

    CHECK_TRUE(solidCells > 0.0, "the body produced no solid cells");
    CHECK_NEAR(dotFluid(one, one),
        S.fluidCells,
        0.0,
        "a unit-vector dot over a domain with %.0f solid cells is not the fluid cell "
        "count",
        solidCells);
    CHECK_TRUE(S.fluidCells + solidCells == (double)IMAX * JMAX * KMAX,
        "the fluid and solid counts %.0f and %.0f do not add up to the domain",
        S.fluidCells,
        solidCells);

    free(one);
  }

  /* ---- Convergence to the exact discrete solution ---- */

  setup(&comm, OUTFLOW, NULL, "jacobi");
  runExactCase(OUTFLOW, "dirichlet, obstacle-free");

  setup(&comm, NOSLIP, NULL, "jacobi");
  runExactCase(NOSLIP, "neumann, obstacle-free");

  setup(&comm, OUTFLOW, sphere, "jacobi");
  runExactCase(OUTFLOW, "dirichlet, sphere");

  setup(&comm, NOSLIP, sphere, "jacobi");
  runExactCase(NOSLIP, "neumann, sphere");

  /* ---- Solid cells are zero at an intermediate iteration ---- */
  {
    setup(&comm, OUTFLOW, sphere, "jacobi");

    size_t size   = fieldSize();
    double *exact = calloc(size, sizeof(double));

    buildExactProblem(exact);

    for (size_t i = 0; i < size; i++) {
      D.p[i] = 0.0;
    }

    /* Few enough that the solve is nowhere near converged: the claim is about
     * the iterate, not about the answer. */
    cgTestSolveSteps(&S, D.p, D.rhs, 5);

    CHECK_TRUE(cgTestIterations(&S) == 5,
        "the stopped solve took %d iterations rather than the 5 it was given",
        cgTestIterations(&S));

    int imaxLocal        = D.comm.imaxLocal;
    int jmaxLocal        = D.comm.jmaxLocal;
    int kmaxLocal        = D.comm.kmaxLocal;
    const double *Lambda = D.Lambda;
    const double *r      = cgTestResidual(&S);
    const double *dir    = cgTestDirection(&S);

    int nonzeroP = 0, nonzeroR = 0, nonzeroD = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) > 0.0) {
            continue;
          }
          size_t idx = linear(i, j, k);
          if (D.p[idx] != 0.0) {
            ++nonzeroP;
          }
          if (r[idx] != 0.0) {
            ++nonzeroR;
          }
          if (dir[idx] != 0.0) {
            ++nonzeroD;
          }
        }
      }
    }

    double totals[3] = { (double)nonzeroP, (double)nonzeroR, (double)nonzeroD };
    for (int i = 0; i < 3; i++) {
      commReduceAll(&totals[i], SUM);
    }

    CHECK_TRUE(totals[0] == 0.0,
        "%.0f solid cells hold a nonzero pressure five iterations in",
        totals[0]);
    CHECK_TRUE(totals[1] == 0.0,
        "%.0f solid cells hold a nonzero residual five iterations in",
        totals[1]);
    CHECK_TRUE(totals[2] == 0.0,
        "%.0f solid cells hold a nonzero search direction five iterations in",
        totals[2]);

    free(exact);
  }

  /* ---- Solid content does not reach the iteration ----
   *
   * The requirement: growing the solid part of the domain, with the fluid
   * region and its solution unchanged, must change neither the iteration count
   * nor the converged fluid field.
   *
   * Two bodies cannot both enclose the same fluid region and differ in volume
   * -- a hollow one would leave a sealed pocket, which the geometry validation
   * refuses -- so the comparison is made the other way round, holding the
   * geometry fixed and varying what the solid part of the domain contains. The
   * right-hand side is given large values inside the body between two otherwise
   * identical solves. A solver whose reductions, preconditioner or operator
   * reached into the solid would see a different problem; this one must not,
   * and the assertion is bit equality rather than a tolerance, because the two
   * runs are meant to execute the same arithmetic and not merely to agree.
   */
  {
    size_t size   = fieldSize();

    setup(&comm, OUTFLOW, sphere, "jacobi");

    double *exact = calloc(size, sizeof(double));
    double *clean = calloc(size, sizeof(double));

    buildExactProblem(exact);

    for (size_t i = 0; i < size; i++) {
      D.p[i] = 0.0;
    }
    solve(&S, D.p, D.rhs);

    int itersClean = cgTestIterations(&S);
    memcpy(clean, D.p, size * sizeof(double));

    int imaxLocal        = D.comm.imaxLocal;
    int jmaxLocal        = D.comm.jmaxLocal;
    int kmaxLocal        = D.comm.kmaxLocal;
    const double *Lambda = D.Lambda;
    double solidCells    = 0.0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) == 0.0) {
            D.rhs[linear(i, j, k)] = 1.0e6 * (double)((i + j + k) % 7 + 1);
            solidCells += 1.0;
          }
        }
      }
    }

    commReduceAll(&solidCells, SUM);

    for (size_t i = 0; i < size; i++) {
      D.p[i] = 0.0;
    }
    solve(&S, D.p, D.rhs);

    int itersLoaded  = cgTestIterations(&S);

    int differing    = 0;
    int nonzeroSolid = 0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          size_t idx = linear(i, j, k);

          if (LAM(i, j, k) == 0.0) {
            if (D.p[idx] != 0.0) {
              ++nonzeroSolid;
            }
            continue;
          }

          if (D.p[idx] != clean[idx]) {
            ++differing;
          }
        }
      }
    }

    double totals[2] = { (double)differing, (double)nonzeroSolid };
    for (int i = 0; i < 2; i++) {
      commReduceAll(&totals[i], SUM);
    }

    CHECK_TRUE(solidCells > 0.0, "the body produced no solid cells to load");
    CHECK_TRUE(itersLoaded == itersClean,
        "loading the body changed the iteration count from %d to %d",
        itersClean,
        itersLoaded);
    CHECK_TRUE(totals[0] == 0.0,
        "loading the body changed the converged fluid field at %.0f cells",
        totals[0]);
    CHECK_TRUE(totals[1] == 0.0,
        "%.0f solid cells hold a nonzero pressure after the loaded solve",
        totals[1]);

    if (commIsMaster(&D.comm)) {
      printf("solid content: %.0f solid cells loaded, %d iterations either way, "
             "%.0f fluid cells changed\n",
          solidCells,
          itersClean,
          totals[0]);
    }

    free(clean);
    free(exact);
  }

  /* ---- Preconditioning does not cost iterations ---- */
  {
    struct {
      int boundary;
      const char *geometry;
      const char *label;
    } cases[] = {
      { OUTFLOW, NULL, "dirichlet" },
      { NOSLIP, NULL, "neumann" },
      { OUTFLOW, "sphere:1.0,0.75,0.5,0.3", "dirichlet, sphere" },
    };

    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
      setup(&comm, cases[c].boundary, cases[c].geometry, "none");
      size_t size   = fieldSize();
      double *exact = calloc(size, sizeof(double));

      buildExactProblem(exact);
      for (size_t i = 0; i < size; i++) {
        D.p[i] = 0.0;
      }
      solve(&S, D.p, D.rhs);
      int plain = cgTestIterations(&S);

      setup(&comm, cases[c].boundary, cases[c].geometry, "jacobi");
      buildExactProblem(exact);
      for (size_t i = 0; i < size; i++) {
        D.p[i] = 0.0;
      }
      solve(&S, D.p, D.rhs);
      int preconditioned = cgTestIterations(&S);

      CHECK_TRUE(preconditioned <= plain,
          "%s: the jacobi solve took %d iterations against %d unpreconditioned, so "
          "the preconditioner costs rather than saves",
          cases[c].label,
          preconditioned,
          plain);

      if (commIsMaster(&comm)) {
        printf("%s: %d iterations unpreconditioned, %d with jacobi\n",
            cases[c].label,
            plain,
            preconditioned);
      }

      free(exact);
    }
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

#endif
