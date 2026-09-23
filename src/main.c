/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file. */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "allocate.h"
#include "discretization.h"
#include "forces.h"
#ifdef TEST
#include "fielddump.h"
#endif
#include "parameter.h"
#include "particletracing.h"
#include "profiler.h"
#include "progress.h"
#include "solver.h"
#include "timing.h"
#include "vtkWriter.h"

int main(int argc, char **argv)
{
  double timeStart;
  double timeStop;
  Parameter p;
  Solver s;
  Discretization d;

  commInit(&d.comm, argc, argv);
  initParameter(&p);
  FILE *fp = NULL;
  if (commIsMaster(&d.comm)) {
    fp = initResidualWriter();
  }

  if (argc != 2) {
    printf("Usage: %s <configFile>\n", argv[0]);
    exit(EXIT_SUCCESS);
  }

  readParameter(&p, argv[1]);
  commPartition(&d.comm, p.kmax, p.jmax, p.imax);

  if (commIsMaster(&d.comm)) {
    printParameter(&p);
  }

  initDiscretization(&d, &p);
  initSolver(&s, &d, &p);
  initProfiler(&d.comm);

  /* A body in the domain means the force on it is worth recording: it is what
   * the reference benchmarks are defined in terms of, and it costs a pass over
   * the body's surface. */
  FILE *forceFile = NULL;
  int haveBody    = forcesHaveBody(&d);

  if (haveBody && commIsMaster(&d.comm)) {
    forceFile = fopen("forces.dat", "w");
    if (forceFile != NULL) {
      fprintf(forceFile, "# time fx fy fz\n");
    }
  }

  ParticleTracerType tracer;
  particleTracerInit(&tracer, &d, &p);
#ifndef VERBOSE
  initProgress(&d.comm, d.te);
#endif

  double tau = d.tau;
  double te  = d.te;
  double t   = 0.0;
  int nt     = 0;
  double res = 0.0;

  timeStart  = getTimeStamp();
  while (t <= te) {
    if (tau > 0.0) {
      computeTimestep(&d);
    }
    setBoundaryConditions(&d);
    setSpecialBoundaryCondition(&d);
    computeFG(&d);
    computeRHS(&d);
    /* Every step, not every hundredth: where the operator is singular the
     * right-hand side has to be made compatible before each solve, not
     * occasionally. It is a no-op where a boundary pins the pressure. */
    normalizePressure(&d);
    res = solve(&s, d.p, d.rhs);

    /* A diverged solve returns a residual that is not a finite number, the same
     * one on every rank since it is a global reduction, so every rank stops
     * here together. Nothing after this point would be a result: the field is
     * NaNs, and writing it out is how a diverged run used to pass for a good
     * one. Stopping at the iteration limit is not caught here -- its iterate is
     * finite and a bounded budget is a legitimate setup. */
    if (!solveResidualIsFinite(res)) {
      if (commIsMaster(&d.comm)) {
        fprintf(stderr,
            "The pressure solve diverged at time step %d (t = %f). Stopping "
            "without writing any output.\n",
            nt + 1,
            t);
      }
      exit(EXIT_FAILURE);
    }

    adaptUV(&d);

    if (commIsMaster(&d.comm)) {
      writeResidual(fp, t, res);
    }

    if (haveBody) {
      double fx, fy, fz;
      forcesCompute(&d, &fx, &fy, &fz);

      if (forceFile != NULL) {
        fprintf(forceFile, "%.10e %.10e %.10e %.10e\n", t, fx, fy, fz);
      }
    }

    particleTracerStep(&tracer, &d, t);

    t += d.dt;
    nt++;

#ifdef VERBOSE
    if (commIsMaster(s.comm)) {
      printf("TIME %f , TIMESTEP %f\n", t, d.dt);
    }
#else
    printProgress(t);
#endif
  }
  timeStop = getTimeStamp();
#ifndef VERBOSE
  stopProgress();
#endif
  if (commIsMaster(s.comm)) {
    printf("Solution took %.2fs\n", timeStop - timeStart);
  }

#ifdef TEST
  const char *dumpPath = fieldDumpPath();
  if (dumpPath != NULL) {
    fieldDumpWrite(&d.comm, s.grid, dumpPath, d.p, d.u, d.v, d.w);
  }
#endif

  timeStart = getTimeStamp();
#ifdef _VTK_WRITER_MPI
  VtkOptions opts = { .grid = s.grid, .comm = s.comm };
  vtkOpen(&opts, s.problem);
  vtkScalar(&opts, "pressure", d.p);
  vtkVector(&opts, "velocity", (VtkVector) { d.u, d.v, d.w });
  vtkClose(&opts);
#else
  if (fp != NULL)
    fclose(fp);

  double *pg;
  double *ug;
  double *vg;
  double *wg;

  if (commIsMaster(s.comm)) {
    size_t bytesize = s.grid->imax * s.grid->jmax * s.grid->kmax * sizeof(double);

    pg              = allocate(ARRAY_ALIGNMENT, bytesize);
    ug              = allocate(ARRAY_ALIGNMENT, bytesize);
    vg              = allocate(ARRAY_ALIGNMENT, bytesize);
    wg              = allocate(ARRAY_ALIGNMENT, bytesize);
  }

  commCollectResult(s.comm,
      ug,
      vg,
      wg,
      pg,
      d.u,
      d.v,
      d.w,
      d.p,
      s.grid->kmax,
      s.grid->jmax,
      s.grid->imax);

  if (commIsMaster(s.comm)) {
    VtkOptions opts = { .grid = s.grid };
    vtkOpen(&opts, s.problem);
    vtkScalar(&opts, "pressure", pg);
    vtkVector(&opts, "velocity", (VtkVector) { ug, vg, wg });
    vtkClose(&opts);
  }

#endif

  timeStop = getTimeStamp();

  if (commIsMaster(s.comm)) {
    printf("Result output took %.2fs\n", timeStop - timeStart);
  }

  particleTracerFinalize(&tracer, &d);

  if (forceFile != NULL) {
    fclose(forceFile);
  }

  finalizeProfiler();
  commFinalize(s.comm);
  return EXIT_SUCCESS;
}
