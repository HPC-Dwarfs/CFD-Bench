/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Checks that the pressure null space is handled where there is one and left
 * alone where there is not:
 *
 *  - an enclosed setup, every boundary a wall, is singular, and advancing it
 *    over many steps must not let the mean pressure grow;
 *  - a setup with an outflow is not singular, and the handling must leave its
 *    pressure untouched;
 *  - the right-hand side of a singular system must be made compatible, which
 *    means its mean is removed before the solve.
 */
#include <math.h>
#include <stdlib.h>

#include "check.h"
#include "discretization.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "util.h"

#define IMAX 16
#define JMAX 16
#define KMAX 16
#define STEPS 40

static void buildParameter(
    Parameter *params, int leftRight, const char *name, const char *geometry)
{
  initParameter(params);

  params->name       = (char *)name;
  params->imax       = IMAX;
  params->jmax       = JMAX;
  params->kmax       = KMAX;
  params->xlength    = 1.0;
  params->ylength    = 1.0;
  params->zlength    = 1.0;
  params->eps        = 1e-6;
  params->omg        = 1.7;
  params->itermax    = 5000;
  params->levels     = 2;
  params->presmooth  = 4;
  params->postsmooth = 4;
  params->re         = 100.0;
  params->tau        = 0.5;
  params->gamma      = 0.9;
  params->dt         = 0.005;
  params->te         = 0.0;
  params->gx = params->gy = params->gz = 0.0;
  params->u_init = params->v_init = params->w_init = params->p_init = 0.0;

  params->geometryFile = (char *)geometry;

  params->bcLeft = params->bcRight = leftRight;
  params->bcBottom = params->bcTop = NOSLIP;
  params->bcFront = params->bcBack = NOSLIP;
}

static double meanPressure(Discretization *d)
{
  int imaxLocal = d->comm.imaxLocal;
  int jmaxLocal = d->comm.jmaxLocal;
  int kmaxLocal = d->comm.kmaxLocal;
  double *p     = d->p;
  double sum    = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        sum += P(i, j, k);
      }
    }
  }

  commReduceAll(&sum, SUM);
  return sum / ((double)IMAX * JMAX * KMAX);
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  /* commFinalize frees the derived datatypes, which only commPartition creates,
   * so the communicator this driver finalizes has to have been partitioned. */
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("nullspace");

  /* The enclosed case: singular, and it must stay bounded. */
  {
    Parameter params;
    Discretization d;
    Solver s;

    buildParameter(&params, NOSLIP, "check-nullspace-closed", NULL);
    d.comm = comm;
    commPartition(&d.comm, KMAX, JMAX, IMAX);
    initDiscretization(&d, &params);
    initSolver(&s, &d, &params);
    initProfiler(&d.comm);

    CHECK_TRUE(pressureBcIsSingular(&s.bc),
        "an all-wall setup was not detected as singular");

    int imaxLocal      = d.comm.imaxLocal;
    int jmaxLocal      = d.comm.jmaxLocal;
    int kmaxLocal      = d.comm.kmaxLocal;
    int offsets[NDIMS] = { 0, 0, 0 };
    commGetOffsets(&d.comm, offsets, KMAX, JMAX, IMAX);

    /* Run the same sequence twice: once with the null-space handling and once
     * without it. The absolute size of the mean after a solve is set by the
     * solve tolerance, so what matters is not a fixed bound but that the
     * handling stops the mean from walking away step after step. */
    double worstMean = 0.0;
    double worstMeanUnhandled = 0.0;

    for (int pass = 0; pass < 2; pass++) {
      int handled = (pass == 0);
      double worst = 0.0;

      for (size_t i = 0; i < (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);
           i++) {
        d.p[i] = 0.0;
      }

      for (int step = 0; step < STEPS; step++) {
        /* A right-hand side with a deliberate constant component, which is
         * exactly what the range projection has to remove. Without it the
         * singular solve has no solution and the iterate walks away. */
        double *rhs = d.rhs;

        for (int k = 1; k < kmaxLocal + 1; k++) {
          for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
              double x = ((i - 1 + offsets[IDIM]) + 0.5) / (double)IMAX;
              double y = ((j - 1 + offsets[JDIM]) + 0.5) / (double)JMAX;
              double z = ((k - 1 + offsets[KDIM]) + 0.5) / (double)KMAX;

              RHS(i, j, k) = cos(M_PI * x) * cos(M_PI * y) * cos(M_PI * z) + 1.0;
            }
          }
        }

        if (handled) {
          normalizePressure(&d);
        }
        solve(&s, d.p, d.rhs);

        double mean = fabs(meanPressure(&d));
        if (mean > worst) {
          worst = mean;
        }
      }

      if (handled) {
        worstMean = worst;
      } else {
        worstMeanUnhandled = worst;
      }
    }

    CHECK_TRUE(worstMean < 1e-3,
        "mean pressure of an enclosed setup reached %.3e over %d steps",
        worstMean,
        STEPS);

    /* The comparative form only means something when there is a drift to
     * remove. A relaxation solver left to itself walks away over these forty
     * steps, so the ratio is the assertion that the handling is what stops it.
     * A Krylov solver projects the constant out at every iteration of its own
     * accord and so never accumulates one, leaving the ratio a comparison of
     * two numbers that are both rounding -- which is a pass, not a defect, and
     * the absolute bound above is what carries the requirement there. */
    if (worstMeanUnhandled > 1e-6) {
      CHECK_TRUE(worstMean < 0.01 * worstMeanUnhandled,
          "null-space handling barely changed the drift: %.3e with it against %.3e "
          "without",
          worstMean,
          worstMeanUnhandled);
    } else {
      CHECK_TRUE(worstMean < 1e-6,
          "this solver does not accumulate a constant of its own (%.3e without the "
          "handling), but the handled run still reached %.3e",
          worstMeanUnhandled,
          worstMean);
    }

    if (commIsMaster(&d.comm)) {
      printf("enclosed over %d steps: worst mean pressure %.3e with the handling, "
             "%.3e without\n",
          STEPS,
          worstMean,
          worstMeanUnhandled);
    }
  }

  /*
   * The enclosed case with a body in it. The null vector of the operator is the
   * constant over the fluid and zero over the solid, not the constant over
   * every cell, so the projection has to average over fluid cells alone and
   * leave the solid ones at zero.
   *
   * A relaxation sweep repairs a solid cell the projection disturbed on its
   * next cut-cell pass, which is why this went unnoticed; a Krylov iteration
   * has no such pass and carries the error into every inner product it forms.
   */
  {
    Parameter params;
    Discretization d;
    Solver s;

    buildParameter(&params,
        NOSLIP,
        "check-nullspace-body",
        "sphere:0.5,0.5,0.5,0.2");
    d.comm = comm;
    commPartition(&d.comm, KMAX, JMAX, IMAX);
    initDiscretization(&d, &params);
    initSolver(&s, &d, &params);
    initProfiler(&d.comm);

    CHECK_TRUE(pressureBcIsSingular(&s.bc),
        "an all-wall setup with a body was not detected as singular");

    int imaxLocal        = d.comm.imaxLocal;
    int jmaxLocal        = d.comm.jmaxLocal;
    int kmaxLocal        = d.comm.kmaxLocal;
    const double *Lambda = d.Lambda;
    size_t size = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

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

    CHECK_TRUE(solidCells > 0.0,
        "the body produced no solid cells, so the projection is untested here");
    CHECK_TRUE(d.fluidCells < (double)IMAX * JMAX * KMAX,
        "the fluid cell count %.0f is the whole domain despite a body",
        d.fluidCells);

    /* A field whose fluid mean is a known constant, with the solid cells where
     * the solver leaves them. */
    {
      double *p            = d.p;
      double *rhs          = d.rhs;
      double expectedFluid = 0.0;

      for (size_t i = 0; i < size; i++) {
        d.p[i]   = 0.0;
        d.rhs[i] = 0.0;
      }

      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          for (int i = 1; i < imaxLocal + 1; i++) {
            if (LAM(i, j, k) == 0.0) {
              continue;
            }
            P(i, j, k) = 3.0 + 0.5 * (double)((i + j + k) % 4);
            expectedFluid += P(i, j, k);
          }
        }
      }

      commReduceAll(&expectedFluid, SUM);
      expectedFluid /= d.fluidCells;

      /* What the projection has to remove: the fluid mean, not the mean over
       * every cell, which the body's zeros would drag down. */
      normalizePressure(&d);

      double afterMean  = 0.0;
      int nonzeroSolid  = 0;

      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          for (int i = 1; i < imaxLocal + 1; i++) {
            if (LAM(i, j, k) == 0.0) {
              if (P(i, j, k) != 0.0) {
                ++nonzeroSolid;
              }
              continue;
            }
            afterMean += P(i, j, k);
          }
        }
      }

      commReduceAll(&afterMean, SUM);
      afterMean /= d.fluidCells;

      double solidTotal = (double)nonzeroSolid;
      commReduceAll(&solidTotal, SUM);

      CHECK_TRUE(solidTotal == 0.0,
          "%.0f solid cells hold a nonzero pressure after the projection",
          solidTotal);
      CHECK_NEAR(afterMean,
          0.0,
          1e-12,
          "the fluid mean after the projection is not zero, so the constant "
          "removed was not the fluid mean (it was %.17g)",
          expectedFluid);
    }

    /* And the solve itself: this solver converges, ends with the body still at
     * exactly zero, and agrees with whatever the others produce. The
     * cross-solver comparison is the shell gate's job -- what a single driver
     * can assert is that the field is the unique one with zero fluid mean,
     * which every solver must therefore reach. */
    {
      int offsets[NDIMS] = { 0, 0, 0 };
      commGetOffsets(&d.comm, offsets, KMAX, JMAX, IMAX);

      double *p   = d.p;
      double *rhs = d.rhs;

      for (size_t i = 0; i < size; i++) {
        d.p[i] = 0.0;
      }

      for (int step = 0; step < 10; step++) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
          for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
              double x     = ((i - 1 + offsets[IDIM]) + 0.5) / (double)IMAX;
              double y     = ((j - 1 + offsets[JDIM]) + 0.5) / (double)JMAX;
              double z     = ((k - 1 + offsets[KDIM]) + 0.5) / (double)KMAX;

              RHS(i, j, k) = LAM(i, j, k) *
                             (cos(M_PI * x) * cos(M_PI * y) * cos(M_PI * z) + 1.0);
            }
          }
        }

        normalizePressure(&d);
        solve(&s, d.p, d.rhs);
      }

      double res           = 0.0;
      int nonzeroSolid     = 0;
      double worstFluidAbs = 0.0;

      PressureLevelType lv;
      pressureLevelFromSolver(&s, &lv);
      res = pressureResidualNorm(&lv, d.p, d.rhs);

      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          for (int i = 1; i < imaxLocal + 1; i++) {
            if (LAM(i, j, k) == 0.0) {
              if (P(i, j, k) != 0.0) {
                ++nonzeroSolid;
              }
              continue;
            }
            if (fabs(P(i, j, k)) > worstFluidAbs) {
              worstFluidAbs = fabs(P(i, j, k));
            }
          }
        }
      }

      double solidTotal = (double)nonzeroSolid;
      commReduceAll(&solidTotal, SUM);
      commReduceAll(&worstFluidAbs, MAX);

      CHECK_TRUE(res < params.eps * params.eps,
          "a singular setup with a body did not converge: %.3e against eps^2 %.3e",
          res,
          params.eps * params.eps);
      CHECK_TRUE(worstFluidAbs < 1e3,
          "the fluid pressure reached %.3e over ten steps, so it is drifting",
          worstFluidAbs);

      /* Reported rather than asserted. Whether a solver leaves the body at zero
       * after a solve is that solver's property and tests/checks/obstacle.c
       * already asserts it; repeating it here would only make one defect fail
       * two drivers. What this driver is about is the projection, which is
       * checked above and holds under every solver. */
      if (commIsMaster(&d.comm)) {
        printf("enclosed with a body: %.0f solid cells, residual %.3e, largest fluid "
               "pressure %.3e, %.0f solid cells left nonzero by the solve\n",
            solidCells,
            sqrt(res),
            worstFluidAbs,
            solidTotal);
      }
    }
  }

  /* The outflow case: not singular, and the handling must change nothing. */
  {
    Parameter params;
    Discretization d;
    Solver s;

    buildParameter(&params, OUTFLOW, "check-nullspace-open", NULL);
    d.comm = comm;
    commPartition(&d.comm, KMAX, JMAX, IMAX);
    initDiscretization(&d, &params);
    initSolver(&s, &d, &params);

    CHECK_TRUE(!pressureBcIsSingular(&s.bc),
        "a setup with an outflow was reported as singular");

    int imaxLocal = d.comm.imaxLocal;
    int jmaxLocal = d.comm.jmaxLocal;
    int kmaxLocal = d.comm.kmaxLocal;
    size_t size = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

    double *before = malloc(size * sizeof(double));
    double *rhsBefore = malloc(size * sizeof(double));

    for (size_t i = 0; i < size; i++) {
      d.p[i]   = 0.5 + 0.001 * (double)(i % 17);
      d.rhs[i] = 1.0 + 0.002 * (double)(i % 13);
      before[i]    = d.p[i];
      rhsBefore[i] = d.rhs[i];
    }

    normalizePressure(&d);

    int changed = 0;
    for (size_t i = 0; i < size; i++) {
      if (d.p[i] != before[i] || d.rhs[i] != rhsBefore[i]) {
        ++changed;
      }
    }

    CHECK_TRUE(changed == 0,
        "null-space handling changed %d values in a non-singular setup",
        changed);

    free(before);
    free(rhsBefore);
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
