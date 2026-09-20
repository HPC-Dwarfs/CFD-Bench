/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The benchmark claim: the interior pressure sweep costs the same whether the
 * domain contains a body or not, and does not care how complicated that body
 * is. Obstacle handling lives in a separate pass whose work scales with the
 * body's surface rather than with the domain, so on a grid of this size it is a
 * few per cent of the sweep at most.
 *
 * Three configurations on the same grid:
 *
 *   none      no body at all
 *   sphere    one sphere, the simplest body with a curved surface
 *   lattice   125 small cubes of the same total volume and about six times the
 *             surface area
 *
 * The tolerance is set by the measurement rather than by the solver: a run is
 * repeated and the fastest taken, because the fastest is the one least
 * disturbed by everything else on the machine.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "check.h"
#include "discretization.h"
#include "parameter.h"
#include "profiler.h"
#include "solver.h"
#include "surface-list.h"
#include "timing.h"
#include "util.h"

#define IMAX 48
#define JMAX 48
#define KMAX 48

#define ITERATIONS 60
#define REPEATS 5

/* How far apart two timings may be before the difference is called real. The
 * surface pass is a few per cent of the sweep at this size, and a quiet machine
 * still varies by a few per cent between runs. */
#define TOLERANCE 0.20

static Parameter Params;
static Discretization D;
static Solver S;

typedef struct {
  const char *label;
  double nsPerUpdate;    /* interior sweep only */
  double nsPerUpdateAll; /* interior sweep plus the cut-cell correction */
  /* The same two numbers for pressureApplyOperator. A relaxation solver touches
   * the operator once per sweep; a Krylov solver applies it once per iteration,
   * so the claim has to hold for the application in its own right and not only
   * as a property of the sweep it was modelled on. */
  double nsPerApply;
  double nsPerApplyAll;
  double surfaceCells;
  double solidCells;
} ResultType;

static ResultType measure(CommType *base, const char *label, const char *geometry)
{
  initParameter(&Params);

  Params.name         = "sweeptime";
  Params.imax         = IMAX;
  Params.jmax         = JMAX;
  Params.kmax         = KMAX;
  Params.xlength      = 4.0;
  Params.ylength      = 4.0;
  Params.zlength      = 4.0;
  /* A tolerance of zero means the loop never exits early, so every
   * configuration performs exactly the same number of sweeps. */
  Params.eps          = 0.0;
  Params.omg          = 1.7;
  Params.itermax      = ITERATIONS;
  Params.levels       = 1;
  Params.presmooth    = 5;
  Params.postsmooth   = 5;
  Params.re           = 100.0;
  Params.tau          = 0.5;
  Params.gamma        = 0.9;
  Params.dt           = 0.02;
  Params.te           = 0.0;
  Params.gx = Params.gy = Params.gz = 0.0;
  Params.u_init = Params.v_init = Params.w_init = Params.p_init = 0.0;
  Params.geometryFile = (char *)geometry;

  Params.bcLeft = Params.bcRight = NOSLIP;
  Params.bcBottom = Params.bcTop = NOSLIP;
  Params.bcFront = Params.bcBack = NOSLIP;

  D.comm = *base;
  commPartition(&D.comm, KMAX, JMAX, IMAX);
  initDiscretization(&D, &Params);
  initSolver(&S, &D, &Params);
  initProfiler(&D.comm);

  int imaxLocal = D.comm.imaxLocal;
  int jmaxLocal = D.comm.jmaxLocal;
  int kmaxLocal = D.comm.kmaxLocal;
  size_t size   = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

  double best        = 0.0;
  double bestSurface = 0.0;

  for (int r = 0; r < REPEATS; r++) {
    for (size_t i = 0; i < size; i++) {
      D.p[i]   = 0.0;
      D.rhs[i] = 1e-3 * (double)(i % 31);
    }

    pressureSweepTimesReset();
    solve(&S, D.p, D.rhs);

    double bulk, surface;
    pressureSweepTimes(&bulk, &surface);

    if (r == 0 || bulk < best) {
      best        = bulk;
      bestSurface = surface;
    }
  }

  commReduceAll(&best, MAX);
  commReduceAll(&bestSurface, MAX);

  /* The operator application, timed the same way: the same number of passes
   * over the same grid, fastest repeat taken. */
  double bestApply        = 0.0;
  double bestApplySurface = 0.0;

  {
    PressureLevelType lv;
    pressureLevelFromSolver(&S, &lv);

    double *x = calloc(size, sizeof(double));
    double *y = calloc(size, sizeof(double));

    for (int r = 0; r < REPEATS; r++) {
      for (size_t i = 0; i < size; i++) {
        x[i] = 1e-3 * (double)(i % 31);
        y[i] = 0.0;
      }

      pressureSweepTimesReset();

      for (int it = 0; it < ITERATIONS; it++) {
        pressureApplyOperator(&lv, x, y);
      }

      double bulk, surface;
      pressureSweepTimes(&bulk, &surface);

      if (r == 0 || bulk < bestApply) {
        bestApply        = bulk;
        bestApplySurface = surface;
      }
    }

    free(x);
    free(y);
  }

  commReduceAll(&bestApply, MAX);
  commReduceAll(&bestApplySurface, MAX);

  double surface = (double)S.surface.count;
  double solid   = 0.0;

  {
    const double *Lambda = D.Lambda;
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) == 0.0) {
            solid += 1.0;
          }
        }
      }
    }
  }

  commReduceAll(&surface, SUM);
  commReduceAll(&solid, SUM);

  double updates = (double)IMAX * JMAX * KMAX * ITERATIONS;

  ResultType out;
  out.label          = label;
  out.nsPerUpdate    = best * 1e9 / updates;
  out.nsPerUpdateAll = (best + bestSurface) * 1e9 / updates;
  out.nsPerApply     = bestApply * 1e9 / updates;
  out.nsPerApplyAll  = (bestApply + bestApplySurface) * 1e9 / updates;
  out.surfaceCells = surface;
  out.solidCells   = solid;

  return out;
}

static void report(const ResultType *r)
{
  printf("%-8s sweep %7.3f ns/update (with correction %7.3f), apply %7.3f (%7.3f), "
         "%6.0f solid, %6.0f surface cells\n",
      r->label,
      r->nsPerUpdate,
      r->nsPerUpdateAll,
      r->nsPerApply,
      r->nsPerApplyAll,
      r->solidCells,
      r->surfaceCells);
}

/* The bulk time of one pass, picked out of a result so the two comparisons
 * below differ only in which pass they are about. */
typedef double (*SelectorType)(const ResultType *);

static double sweepBulk(const ResultType *r) { return r->nsPerUpdate; }
static double applyBulk(const ResultType *r) { return r->nsPerApply; }

static void compareWith(const ResultType *a,
    const ResultType *b,
    SelectorType pick,
    const char *pass,
    const char *what)
{
  double ta  = pick(a);
  double tb  = pick(b);
  double rel = fabs(ta - tb) / (0.5 * (ta + tb));

  CHECK_TRUE(rel < TOLERANCE,
      "%s (%s): %s took %.3f ns per update and %s took %.3f, a difference of %.1f%% "
      "against a %.0f%% tolerance",
      what,
      pass,
      a->label,
      ta,
      b->label,
      tb,
      100.0 * rel,
      100.0 * TOLERANCE);
}

static void compare(const ResultType *a, const ResultType *b, const char *what)
{
  compareWith(a, b, sweepBulk, "sweep", what);
  compareWith(a, b, applyBulk, "apply", what);
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("sweeptime");

  ResultType none    = measure(&comm, "none", NULL);
  ResultType sphere  = measure(&comm, "sphere", "sphere:2.0,2.0,2.0,1.0");
  ResultType lattice = measure(&comm, "lattice", "tests/geom/lattice.vox");

  if (commIsMaster(&comm)) {
    report(&none);
    report(&sphere);
    report(&lattice);
  }

  CHECK_TRUE(none.surfaceCells == 0.0,
      "the obstacle-free configuration has %.0f surface cells",
      none.surfaceCells);
  CHECK_TRUE(sphere.solidCells > 0.0 && lattice.solidCells > 0.0,
      "one of the bodies produced no solid cells");

  /* The two bodies have to be comparable in volume and very different in
   * surface, otherwise the comparison below says nothing. */
  double volRatio = lattice.solidCells / sphere.solidCells;
  CHECK_TRUE(volRatio > 0.6 && volRatio < 1.6,
      "the two bodies differ in solid volume by a factor of %.2f, so they are not "
      "comparable",
      volRatio);

  double surfRatio = lattice.surfaceCells / sphere.surfaceCells;
  CHECK_TRUE(surfRatio > 2.0,
      "the lattice has only %.2f times the surface of the sphere, so it is not the "
      "harder case it is meant to be",
      surfRatio);

  compare(&none, &sphere, "a body must not change the interior sweep");
  compare(&sphere, &lattice, "the shape of the body must not change it either");

  if (commIsMaster(&comm)) {
    printf("solid volume ratio %.2f, surface ratio %.2f\n", volRatio, surfRatio);
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
