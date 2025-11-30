/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file. */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "allocate.h"
#include "discretization.h"
#include "parameter.h"
#include "progress.h"
#include "solver.h"
#include "timing.h"
#include "vtkWriter.h"

int main(int argc, char **argv)
{
  double timeStart, timeStop;
  Parameter p;
  Solver s;
  Discretization d;

  commInit(&d.comm, argc, argv);
  initParameter(&p);
  FILE *fp;
  if (commIsMaster(&d.comm))
    fp = initResidualWriter();

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

#ifndef VERBOSE
  initProgress(d.te);
#endif

  double tau = d.tau;
  double te  = d.te;
  double t   = 0.0;
  int nt     = 0;
  double res = 0.0;

  timeStart  = getTimeStamp();
  while (t <= te) {
    if (tau > 0.0)
      computeTimestep(&d);
    setBoundaryConditions(&d);
    setSpecialBoundaryCondition(&d);
    computeFG(&d);
    computeRHS(&d);
    if (nt % 100 == 0)
      normalizePressure(&d);
    res = solve(&s, d.p, d.rhs);
    adaptUV(&d);

    if (commIsMaster(&d.comm))
      writeResidual(fp, t, res);

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

  timeStart = getTimeStamp();
#ifdef _VTK_WRITER_MPI
  VtkOptions opts = { .grid = s.grid, .comm = s.comm };
  vtkOpen(&opts, s.problem);
  vtkScalar(&opts, "pressure", d.p);
  vtkVector(&opts, "velocity", (VtkVector) { d.u, d.v, d.w });
  vtkClose(&opts);
#else
  if (commIsMaster(&d.comm))
    fclose(fp);

  double *pg, *ug, *vg, *wg;

  if (commIsMaster(s.comm)) {
    size_t bytesize = s.grid->imax * s.grid->jmax * s.grid->kmax * sizeof(double);

    pg              = allocate(64, bytesize);
    ug              = allocate(64, bytesize);
    vg              = allocate(64, bytesize);
    wg              = allocate(64, bytesize);
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

  commFinalize(s.comm);
  return EXIT_SUCCESS;
}
