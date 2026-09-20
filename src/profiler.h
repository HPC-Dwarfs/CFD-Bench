/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <stddef.h>
#include <stdio.h>

#include "comm.h"

#define MAX_STR_LENGTH 30
/* BULK is the geometry-free interior sweep and SURFACE the cut-cell correction
 * that follows it. They are timed apart because the benchmark's claim is about
 * the interior sweep specifically: the correction is a separate pass whose work
 * scales with the body's surface, so lumping the two together would measure
 * something the claim does not make. */
/* CG_DOT, CG_AXPY and PRECON are the Krylov solver's own work: the global
 * inner products, the vector updates, and the preconditioner application. The
 * first is the interesting one -- two allreduces an iteration is a cost a
 * Krylov method has and the relaxation solvers do not, so it is measured rather
 * than inferred. They stay at zero in a build that links another solver. */
typedef enum {
  SOLVER = 0,
  COMM,
  SWEEP_BULK,
  SWEEP_SURFACE,
  CG_DOT,
  CG_AXPY,
  PRECON,
  NUMREGIONS
} RegionsType;

#ifdef PROFILING
#define PROFILE(tag, call)                                                               \
  const double ts = getTimeStamp();                                                      \
  call;                                                                                  \
  T[tag] += (getTimeStamp() - ts);                                                       \
  C[tag]++;

#define TIMESTART double timeStart = getTimeStamp();

#define TIMESTOP(tag)                                                                    \
  T[tag] += (getTimeStamp() - timeStart);                                                \
  C[tag]++;

extern double T[NUMREGIONS];
extern size_t C[NUMREGIONS];
#else
#define PROFILE(tag, call) call;
#define TIMESTART
#define TIMESTOP(tag)

#endif

extern void printProfile(CommType *c, int it);
extern void initProfiler(CommType *c);
extern void finalizeProfiler(void);
