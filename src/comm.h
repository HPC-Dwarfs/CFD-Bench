/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __COMM_H_
#define __COMM_H_

#include "stdbool.h"

#if defined(_MPI)
#include <mpi.h>
#endif
/*
 * Spatial directions:
 * ICORD (0) from 0 (LEFT) to imax (RIGHT)
 * JCORD (1) from 0 (BOTTOM) to jmax (TOP)
 * KCORD (2) from 0 (FRONT) to kmax (BACK)
 * All derived Subarray types are in C ordering
 * with indices KDIM (0), JDIM(1) and IDIM(2)
 * */
typedef enum direction { LEFT = 0, RIGHT, BOTTOM, TOP, FRONT, BACK, NDIRS } DirectionType;
typedef enum coordinates { ICORD = 0, JCORD, KCORD, NCORDS } CoordinatesType;
typedef enum dimension { KDIM = 0, JDIM, IDIM, NDIMS } DimensionType;
enum layer { HALO = 0, BULK };
enum op { MIN = 0, MAX, SUM };

typedef struct {
  int rank;
  int size;
#if defined(_MPI)
  MPI_Comm comm;
  MPI_Datatype sbufferTypes[NDIRS];
  MPI_Datatype rbufferTypes[NDIRS];
#endif
  int neighbours[NDIRS];
  int coords[NDIMS], dims[NDIMS];
  int imaxLocal, jmaxLocal, kmaxLocal;
} CommType;

extern void commInit(CommType *c, int argc, char **argv);
extern void commPartition(CommType *c, int kmax, int jmax, int imax);
extern void commFinalize(CommType *comm);
extern void commPrintConfig(CommType *);
extern void commExchange(CommType *, double *);
extern void commShift(CommType *c, double *f, double *g, double *h);
extern void commReduceAll(double *v, int op);
/* The same all-reduce over several values at once, in place. A Krylov solver
 * forms two inner products per iteration and needs both on every rank, so it
 * fuses them into one call and pays one latency instead of two. Not
 * commReduce: that one reduces to the master and exists for the profiler's
 * reporting, which is the one place a result nobody else sees is enough. */
extern void commReduceAllN(double *v, int count, int op);
extern void commReduce(double *v, double *o, int count, int op);
extern int commIsBoundary(CommType *c, DirectionType direction);
extern void commGetOffsets(CommType *c, int offsets[], int kmax, int jmax, int imax);
/* Which rank holds the global cell (gi, gj, gk). Particle migration needs it:
 * a particle can land anywhere in one step, so its new owner has to be found
 * from its position rather than by asking the neighbours. */
extern int commRankOfCell(
    CommType *c, int gi, int gj, int gk, int imax, int jmax, int kmax);
extern void commFreeCommunicator(CommType *comm);
extern void commUpdateDatatypes(
    CommType *oldcomm, CommType *newcomm, int imaxLocal, int jmaxLocal, int kmaxLocal);
extern void commCollectResult(CommType *c,
    double *ug,
    double *vg,
    double *wg,
    double *pg,
    double *u,
    double *v,
    double *w,
    double *p,
    int kmax,
    int jmax,
    int imax);

static inline bool commIsMaster(CommType *c)
{
  return c->rank == 0;
}
#endif // __COMM_H_
