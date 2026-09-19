/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Solves a Poisson problem whose answer is known exactly, so that the solver is
 * measured against the discrete operator rather than against a discretization
 * of a continuous problem.
 *
 * The right-hand side is built by applying the solver's own operator to a
 * chosen field, so the exact discrete solution is that field. Any difference
 * after a converged solve is the solver's, not the mesh's.
 *
 * Two configurations:
 *
 *   Dirichlet   every boundary an outflow, so the operator is non-singular and
 *               the solution is unique.
 *   Neumann     every boundary a wall, so the operator is singular and the
 *               solution is unique only up to a constant; both fields are
 *               compared with their means removed.
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

#define IMAX 32
#define JMAX 24
#define KMAX 16

static void buildParameter(Parameter *params, int boundary)
{
  initParameter(params);

  params->name       = "check-poisson";
  params->imax       = IMAX;
  params->jmax       = JMAX;
  params->kmax       = KMAX;
  params->xlength    = 2.0;
  params->ylength    = 1.5;
  params->zlength    = 1.0;

  params->eps        = 1e-11;
  params->omg        = 1.7;
  params->itermax    = 20000;
  params->levels     = 3;
  params->presmooth  = 4;
  params->postsmooth = 4;

  params->re         = 100.0;
  params->tau        = 0.5;
  params->gamma      = 0.9;
  params->dt         = 0.02;
  params->te         = 0.0;
  params->gx = params->gy = params->gz = 0.0;
  params->u_init = params->v_init = params->w_init = params->p_init = 0.0;

  params->bcLeft = params->bcRight = boundary;
  params->bcBottom = params->bcTop = boundary;
  params->bcFront = params->bcBack = boundary;
}

/* The same seven-point operator the solvers relax, applied once. */
static void applyOperator(Solver *s,
    double *p,
    double *out,
    int imaxLocal,
    int jmaxLocal,
    int kmaxLocal,
    double dx,
    double dy,
    double dz)
{
  double idx2 = 1.0 / (dx * dx);
  double idy2 = 1.0 / (dy * dy);
  double idz2 = 1.0 / (dz * dz);

  commExchange(s->comm, p);
  pressureBcApply(&s->bc, s->comm, p, imaxLocal, jmaxLocal, kmaxLocal);

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        out[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)] =
            (P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
            (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
            (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2;
      }
    }
  }
}

static double meanOf(
    const double *f, int imaxLocal, int jmaxLocal, int kmaxLocal, double cells)
{
  double sum = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        sum += f[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)];
      }
    }
  }

  commReduceAll(&sum, SUM);
  return sum / cells;
}

static int runCase(CommType *base, int boundary, const char *label)
{
  Parameter params;
  Discretization d;
  Solver s;

  buildParameter(&params, boundary);

  /* MPI is initialized once by main; each case gets its own partition of the
   * same world. */
  d.comm = *base;
  commPartition(&d.comm, KMAX, JMAX, IMAX);
  initDiscretization(&d, &params);
  initSolver(&s, &d, &params);
  initProfiler(&d.comm);

  int imaxLocal = d.comm.imaxLocal;
  int jmaxLocal = d.comm.jmaxLocal;
  int kmaxLocal = d.comm.kmaxLocal;

  int offsets[NDIMS] = { 0, 0, 0 };
  commGetOffsets(&d.comm, offsets, KMAX, JMAX, IMAX);

  size_t size    = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);
  double *exact  = calloc(size, sizeof(double));
  double cells   = (double)IMAX * JMAX * KMAX;

  double dx      = d.grid.dx;
  double dy      = d.grid.dy;
  double dz      = d.grid.dz;

  /* A smooth field that vanishes on every boundary face, so it satisfies the
   * Dirichlet configuration exactly and is a perfectly ordinary field for the
   * Neumann one. */
  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        double x = ((i - 1 + offsets[IDIM]) + 0.5) * dx / params.xlength;
        double y = ((j - 1 + offsets[JDIM]) + 0.5) * dy / params.ylength;
        double z = ((k - 1 + offsets[KDIM]) + 0.5) * dz / params.zlength;

        exact[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)] =
            sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
      }
    }
  }

  /* rhs = A * exact, so exact is the discrete solution by construction. */
  applyOperator(
      &s, exact, d.rhs, imaxLocal, jmaxLocal, kmaxLocal, dx, dy, dz);

  for (size_t i = 0; i < size; i++) {
    d.p[i] = 0.0;
  }

  double res = solve(&s, d.p, d.rhs);

  CHECK_TRUE(res < params.eps * params.eps * 100.0,
      "%s: solve did not converge, reported mean square residual %.3e against eps^2 "
      "%.3e",
      label,
      res,
      params.eps * params.eps);

  /* An independently computed residual of the returned field must agree with
   * what the solver reported. */
  double independent = pressureResidualNorm(&d.comm,
      &s.bc,
      d.p,
      d.rhs,
      imaxLocal,
      jmaxLocal,
      kmaxLocal,
      dx,
      dy,
      dz,
      cells);

  CHECK_TRUE(fabs(independent - res) <= 1e-12 * fabs(res) + 1e-30,
      "%s: reported residual %.17g is not the residual of the returned field %.17g",
      label,
      res,
      independent);

  /* Compare the fields. The all-Neumann operator fixes the solution only up to
   * a constant, so remove the mean from both before comparing. */
  double shiftP = 0.0, shiftE = 0.0;

  if (boundary == NOSLIP) {
    shiftP = meanOf(d.p, imaxLocal, jmaxLocal, kmaxLocal, cells);
    shiftE = meanOf(exact, imaxLocal, jmaxLocal, kmaxLocal, cells);
  }

  double maxDiff = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        size_t idx =
            (size_t)(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i);
        double diff = fabs((d.p[idx] - shiftP) - (exact[idx] - shiftE));
        if (diff > maxDiff) {
          maxDiff = diff;
        }
      }
    }
  }

  commReduceAll(&maxDiff, MAX);

  CHECK_TRUE(maxDiff < 1e-4,
      "%s: converged field differs from the exact discrete solution by %.3e",
      label,
      maxDiff);

  if (commIsMaster(&d.comm)) {
    printf("%s: residual %.3e, max deviation from exact solution %.3e\n",
        label,
        sqrt(res),
        maxDiff);
  }

  free(exact);
  return 0;
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  /* commFinalize frees the derived datatypes, which only commPartition creates,
   * so the communicator this driver finalizes has to have been partitioned. */
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("poisson");

  runCase(&comm, OUTFLOW, "dirichlet");
  runCase(&comm, NOSLIP, "neumann");

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
