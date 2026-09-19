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

static void buildParameter(Parameter *params, int leftRight, const char *name)
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

    buildParameter(&params, NOSLIP, "check-nullspace-closed");
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
    CHECK_TRUE(worstMean < 0.01 * worstMeanUnhandled,
        "null-space handling barely changed the drift: %.3e with it against %.3e "
        "without",
        worstMean,
        worstMeanUnhandled);

    if (commIsMaster(&d.comm)) {
      printf("enclosed over %d steps: worst mean pressure %.3e with the handling, "
             "%.3e without\n",
          STEPS,
          worstMean,
          worstMeanUnhandled);
    }
  }

  /* The outflow case: not singular, and the handling must change nothing. */
  {
    Parameter params;
    Discretization d;
    Solver s;

    buildParameter(&params, OUTFLOW, "check-nullspace-open");
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
