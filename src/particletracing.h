/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __PARTICLETRACING_H_
#define __PARTICLETRACING_H_

#include "discretization.h"

/*
 * Massless tracer particles, so that flow around an embedded body can be seen
 * as streaklines rather than only as a field.
 *
 * A particle is obstructed by the apertures of the faces its path crosses, not
 * by the volume fraction of the cell it lands in. The difference is the whole
 * point: a body one face thick has fluid on both sides, so a destination-cell
 * test sees nothing wrong and the particle passes straight through it.
 *
 * Tracing is off unless a setup asks for it, and a setup that does not ask
 * produces exactly the fields and files it would with no tracer present.
 */

typedef struct {
  double x, y, z;
} ParticleType;

typedef struct {
  int enabled;

  /* configuration */
  int perBatch;
  double startTime, injectPeriod, writePeriod;
  double x1, y1, z1, x2, y2, z2;

  /* state */
  double lastInject, lastWrite;
  long batch;
  int writeIndex;

  ParticleType *pool;
  int count;
  int capacity;

  /* running totals, summed across ranks when reported */
  long injected;
  long removedBoundary;
  long removedBody;

  /* this rank's slice of the domain, in physical coordinates */
  double xLo, xHi, yLo, yHi, zLo, zHi;
} ParticleTracerType;

extern void particleTracerInit(ParticleTracerType *t, Discretization *d, Parameter *p);

/* Inject if a batch is due, advance every particle by the step the flow solver
 * just took, hand over the ones that left this rank, and write if due. */
extern void particleTracerStep(ParticleTracerType *t, Discretization *d, double time);

/* Report the totals and release the pool. */
extern void particleTracerFinalize(ParticleTracerType *t, Discretization *d);

/* Hand every particle that left this rank to the one that now owns it. Part of
 * a step; separate so a check driver can advance and migrate on its own. */
extern void particleTracerMigrate(ParticleTracerType *t, Discretization *d);

#ifdef TEST
/* Exposed for the check drivers. */
extern void particleTracerSeed(
    const ParticleTracerType *t, long batch, int index, double *x, double *y, double *z);
extern int particleTracerLiveCount(const ParticleTracerType *t);
extern void particleTracerInject(ParticleTracerType *t, Discretization *d);
extern void particleTracerAdvance(ParticleTracerType *t, Discretization *d, double dt);
extern void particleTracerTotals(const ParticleTracerType *t,
    Discretization *d,
    double *injected,
    double *alive,
    double *removedBoundary,
    double *removedBody);
#endif

#endif // __PARTICLETRACING_H_
