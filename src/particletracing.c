/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "particletracing.h"
#include "util.h"

#define PARTICLE_DIR "vis_files"

/* ---------------------------------------------------------------------- */

/*
 * Seed positions come from a counter-based hash of the global particle index
 * and the batch number, not from a random number generator.
 *
 * A generator would have to be seeded, drawn from in the same order on every
 * rank, and would still differ between runs. A hash of the index is the same
 * everywhere by construction: every rank can compute the position of every
 * particle without talking to anyone, and two runs at different rank counts
 * inject exactly the same particles.
 */
static uint64_t mix(uint64_t z)
{
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static double unitInterval(uint64_t h)
{
  /* 53 bits is all a double can hold, and taking the high ones avoids the
   * weaker low bits of the mix. */
  return (double)(h >> 11) * (1.0 / 9007199254740992.0);
}

void particleTracerSeed(
    const ParticleTracerType *t, long batch, int index, double *x, double *y, double *z)
{
  uint64_t base = mix((uint64_t)batch * 0x100000001b3ULL + (uint64_t)index);

  *x            = t->x1 + (t->x2 - t->x1) * unitInterval(mix(base + 1));
  *y            = t->y1 + (t->y2 - t->y1) * unitInterval(mix(base + 2));
  *z            = t->z1 + (t->z2 - t->z1) * unitInterval(mix(base + 3));
}

/* ---------------------------------------------------------------------- */

static void reserve(ParticleTracerType *t, int want)
{
  if (want <= t->capacity) {
    return;
  }

  int capacity = (t->capacity > 0) ? t->capacity : 1;

  while (capacity < want) {
    capacity *= 2;
  }

  t->pool     = realloc(t->pool, (size_t)capacity * sizeof(ParticleType));
  t->capacity = capacity;

  if (t->pool == NULL) {
    fprintf(stderr, "particle tracing: cannot grow the pool to %d particles\n", capacity);
    exit(EXIT_FAILURE);
  }
}

/*
 * Removing a particle moves the last one into its place. That reclaims the slot
 * at once rather than leaving a hole for a later compaction pass to find, and
 * keeps the array dense so the advance loop stays a straight walk over live
 * particles.
 */
static void removeAt(ParticleTracerType *t, int i)
{
  t->pool[i] = t->pool[t->count - 1];
  --t->count;
}

int particleTracerLiveCount(const ParticleTracerType *t)
{
  return t->count;
}

/* ---------------------------------------------------------------------- */

static int ownsPosition(const ParticleTracerType *t, double x, double y, double z)
{
  /* Half-open, so a position on a subdomain boundary belongs to exactly one
   * rank however the domain was divided. */
  return x >= t->xLo && x < t->xHi && y >= t->yLo && y < t->yHi && z >= t->zLo &&
         z < t->zHi;
}

static int insideDomain(Discretization *d, double x, double y, double z)
{
  return x >= 0.0 && x < d->grid.xlength && y >= 0.0 && y < d->grid.ylength && z >= 0.0 &&
         z < d->grid.zlength;
}

/* The local cell holding a position, in this rank's 1-based interior indices. */
static void localCell(ParticleTracerType *t,
    Discretization *d,
    double x,
    double y,
    double z,
    int *ci,
    int *cj,
    int *ck)
{
  (void)t;
  *ci = (int)floor(x / d->grid.dx) - d->iOffset + 1;
  *cj = (int)floor(y / d->grid.dy) - d->jOffset + 1;
  *ck = (int)floor(z / d->grid.dz) - d->kOffset + 1;
}

static int cellIsFluid(Discretization *d, int i, int j, int k)
{
  int imaxLocal        = d->comm.imaxLocal;
  int jmaxLocal        = d->comm.jmaxLocal;
  int kmaxLocal        = d->comm.kmaxLocal;
  const double *Lambda = d->Lambda;

  if (i < 0 || i > imaxLocal + 1 || j < 0 || j > jmaxLocal + 1 || k < 0 ||
      k > kmaxLocal + 1) {
    return 0;
  }

  return LAM(i, j, k) > 0.0;
}

/* ---------------------------------------------------------------------- */

/*
 * Trilinear interpolation of one velocity component at its own face location.
 *
 * Each component sits on a different face of the cell, so one stencil cannot
 * serve all three: u is offset by half a cell in y and z relative to the cell
 * centre and sits exactly on the x-faces, and likewise for v and w. The offsets
 * below are what that difference amounts to in index space.
 */
static double interpolate(Discretization *d,
    const double *field,
    double x,
    double y,
    double z,
    double offX,
    double offY,
    double offZ)
{
  int imaxLocal = d->comm.imaxLocal;
  int jmaxLocal = d->comm.jmaxLocal;
  int kmaxLocal = d->comm.kmaxLocal;

  double a      = x / d->grid.dx - d->iOffset + offX;
  double b      = y / d->grid.dy - d->jOffset + offY;
  double c      = z / d->grid.dz - d->kOffset + offZ;

  int i0        = (int)floor(a);
  int j0        = (int)floor(b);
  int k0        = (int)floor(c);

  double fx     = a - i0;
  double fy     = b - j0;
  double fz     = c - k0;

  /* Clamp into the array including its halo. A particle this rank owns is never
   * more than one cell from its interior, so this only ever trims the very edge
   * of the stencil. */
  if (i0 < 0) {
    i0 = 0;
    fx = 0.0;
  }
  if (j0 < 0) {
    j0 = 0;
    fy = 0.0;
  }
  if (k0 < 0) {
    k0 = 0;
    fz = 0.0;
  }
  if (i0 > imaxLocal) {
    i0 = imaxLocal;
    fx = 0.0;
  }
  if (j0 > jmaxLocal) {
    j0 = jmaxLocal;
    fy = 0.0;
  }
  if (k0 > kmaxLocal) {
    k0 = kmaxLocal;
    fz = 0.0;
  }

  int stride = imaxLocal + 2;
  int slab   = (imaxLocal + 2) * (jmaxLocal + 2);

  double sum = 0.0;

  for (int dk = 0; dk < 2; dk++) {
    double wk = dk ? fz : 1.0 - fz;
    for (int dj = 0; dj < 2; dj++) {
      double wj = dj ? fy : 1.0 - fy;
      for (int di = 0; di < 2; di++) {
        double wi = di ? fx : 1.0 - fx;
        size_t idx =
            (size_t)(k0 + dk) * slab + (size_t)(j0 + dj) * stride + (size_t)(i0 + di);
        sum += wi * wj * wk * field[idx];
      }
    }
  }

  return sum;
}

/*
 * Walk the cell faces the segment from (x0,y0,z0) to (x1,y1,z1) crosses, in
 * order, and stop at the first one that is shut.
 *
 * Testing the cell the particle lands in instead would be cheaper and wrong: a
 * body one face thick has fluid on both sides, so the destination cell says
 * nothing was in the way. It would also depend on the step being short enough
 * not to jump a thin body entirely, which an adaptive time step does not
 * promise.
 *
 * Returns 1 when the path is blocked.
 */
static int pathBlocked(ParticleTracerType *t,
    Discretization *d,
    double x0,
    double y0,
    double z0,
    double x1,
    double y1,
    double z1)
{
  int imaxLocal    = d->comm.imaxLocal;
  int jmaxLocal    = d->comm.jmaxLocal;
  const double *Ax = d->Ax;
  const double *Ay = d->Ay;
  const double *Az = d->Az;

  double dx        = d->grid.dx;
  double dy        = d->grid.dy;
  double dz        = d->grid.dz;

  int ci, cj, ck;
  localCell(t, d, x0, y0, z0, &ci, &cj, &ck);

  int ti, tj, tk;
  localCell(t, d, x1, y1, z1, &ti, &tj, &tk);

  double sx = x1 - x0;
  double sy = y1 - y0;
  double sz = z1 - z0;

  int stepI = (sx > 0.0) ? 1 : (sx < 0.0 ? -1 : 0);
  int stepJ = (sy > 0.0) ? 1 : (sy < 0.0 ? -1 : 0);
  int stepK = (sz > 0.0) ? 1 : (sz < 0.0 ? -1 : 0);

  /* Parameter along the segment at which the next face in each direction is
   * reached, and how far apart successive ones are. */
  double tMaxX = INFINITY, tMaxY = INFINITY, tMaxZ = INFINITY;
  double tDeltaX = INFINITY, tDeltaY = INFINITY, tDeltaZ = INFINITY;

  int gi = ci - 1 + d->iOffset;
  int gj = cj - 1 + d->jOffset;
  int gk = ck - 1 + d->kOffset;

  if (stepI != 0) {
    double next = (stepI > 0) ? (gi + 1) * dx : gi * dx;
    tMaxX       = (next - x0) / sx;
    tDeltaX     = dx / fabs(sx);
  }
  if (stepJ != 0) {
    double next = (stepJ > 0) ? (gj + 1) * dy : gj * dy;
    tMaxY       = (next - y0) / sy;
    tDeltaY     = dy / fabs(sy);
  }
  if (stepK != 0) {
    double next = (stepK > 0) ? (gk + 1) * dz : gk * dz;
    tMaxZ       = (next - z0) / sz;
    tDeltaZ     = dz / fabs(sz);
  }

  int budget = abs(ti - ci) + abs(tj - cj) + abs(tk - ck) + 3;

  for (int n = 0; n < budget; n++) {
    if (tMaxX > 1.0 && tMaxY > 1.0 && tMaxZ > 1.0) {
      return 0; /* the segment ends before the next face */
    }

    if (tMaxX <= tMaxY && tMaxX <= tMaxZ) {
      /* Crossing the x-face of this cell, or of the one below it. */
      int face = (stepI > 0) ? ci : ci - 1;
      if (Ax[(size_t)ck * (imaxLocal + 2) * (jmaxLocal + 2) +
              (size_t)cj * (imaxLocal + 2) + (size_t)face] == 0.0) {
        return 1;
      }
      ci += stepI;
      tMaxX += tDeltaX;
    } else if (tMaxY <= tMaxZ) {
      int face = (stepJ > 0) ? cj : cj - 1;
      if (Ay[(size_t)ck * (imaxLocal + 2) * (jmaxLocal + 2) +
              (size_t)face * (imaxLocal + 2) + (size_t)ci] == 0.0) {
        return 1;
      }
      cj += stepJ;
      tMaxY += tDeltaY;
    } else {
      int face = (stepK > 0) ? ck : ck - 1;
      if (Az[(size_t)face * (imaxLocal + 2) * (jmaxLocal + 2) +
              (size_t)cj * (imaxLocal + 2) + (size_t)ci] == 0.0) {
        return 1;
      }
      ck += stepK;
      tMaxZ += tDeltaZ;
    }

    /* The walk should never leave this rank's halo: the time step is limited by
     * the same CFL condition the flow is, so a particle moves well under a cell
     * per step and the segment spans at most one face in each direction. This
     * bounds the loop rather than describing a path that is expected to be
     * taken -- if it ever triggers, the step is longer than the grid and the
     * obstruction test is not the only thing that has stopped being true. */
    if (ci < 0 || ci > imaxLocal + 1 || cj < 0 || cj > jmaxLocal + 1 || ck < 0 ||
        ck > d->comm.kmaxLocal + 1) {
      return 0;
    }
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

void particleTracerInject(ParticleTracerType *t, Discretization *d)
{
  for (int index = 0; index < t->perBatch; index++) {
    double x, y, z;
    particleTracerSeed(t, t->batch, index, &x, &y, &z);

    /* Every rank walks the whole batch and keeps what is its own, which needs
     * no communication and gives the same answer whatever the decomposition. */
    if (!ownsPosition(t, x, y, z)) {
      continue;
    }

    int ci, cj, ck;
    localCell(t, d, x, y, z, &ci, &cj, &ck);

    /* A particle inside the body would have nowhere to go. */
    if (!cellIsFluid(d, ci, cj, ck)) {
      continue;
    }

    reserve(t, t->count + 1);
    t->pool[t->count].x = x;
    t->pool[t->count].y = y;
    t->pool[t->count].z = z;
    ++t->count;
    ++t->injected;
  }

  ++t->batch;
}

void particleTracerAdvance(ParticleTracerType *t, Discretization *d, double dt)
{
  /* adaptUV leaves the velocity halo as the last exchange left it, and the
   * interpolation reads one cell beyond the interior. */
  commExchange(&d->comm, d->u);
  commExchange(&d->comm, d->v);
  commExchange(&d->comm, d->w);

  int i = 0;

  while (i < t->count) {
    double x = t->pool[i].x;
    double y = t->pool[i].y;
    double z = t->pool[i].z;

    /* u sits on the x-faces and at cell centres in y and z; v and w are the
     * same statement rotated. */
    double u  = interpolate(d, d->u, x, y, z, 0.0, 0.5, 0.5);
    double v  = interpolate(d, d->v, x, y, z, 0.5, 0.0, 0.5);
    double w  = interpolate(d, d->w, x, y, z, 0.5, 0.5, 0.0);

    double nx = x + dt * u;
    double ny = y + dt * v;
    double nz = z + dt * w;

    if (pathBlocked(t, d, x, y, z, nx, ny, nz)) {
      ++t->removedBody;
      removeAt(t, i);
      continue;
    }

    if (!insideDomain(d, nx, ny, nz)) {
      ++t->removedBoundary;
      removeAt(t, i);
      continue;
    }

    t->pool[i].x = nx;
    t->pool[i].y = ny;
    t->pool[i].z = nz;
    ++i;
  }
}

/* ---------------------------------------------------------------------- */

#if defined(_MPI)
static MPI_Datatype ParticleDatatype;
static int ParticleDatatypeReady = 0;

static MPI_Datatype particleDatatype(void)
{
  if (!ParticleDatatypeReady) {
    MPI_Type_contiguous(3, MPI_DOUBLE, &ParticleDatatype);
    MPI_Type_commit(&ParticleDatatype);
    ParticleDatatypeReady = 1;
  }

  return ParticleDatatype;
}
#endif

/*
 * Hand every particle that left this rank to the one that owns where it is now.
 *
 * A step can carry a particle anywhere, so the destination is worked out from
 * its position through the Cartesian topology rather than by asking the six
 * neighbours. One collective with variable counts then moves them all; the
 * pairs with nothing to say cost nothing.
 */
void particleTracerMigrate(ParticleTracerType *t, Discretization *d)
{
#if defined(_MPI)
  int size = d->comm.size;

  if (size == 1) {
    return;
  }

  int *sendCount = calloc((size_t)size, sizeof(int));
  int *recvCount = calloc((size_t)size, sizeof(int));
  int *sendDispl = calloc((size_t)size, sizeof(int));
  int *recvDispl = calloc((size_t)size, sizeof(int));

  /* Which rank each departing particle belongs to. */
  int *owner  = malloc((size_t)(t->count > 0 ? t->count : 1) * sizeof(int));
  int leaving = 0;

  for (int i = 0; i < t->count; i++) {
    double x = t->pool[i].x;
    double y = t->pool[i].y;
    double z = t->pool[i].z;

    if (ownsPosition(t, x, y, z)) {
      owner[i] = -1;
      continue;
    }

    int gi = (int)floor(x / d->grid.dx);
    int gj = (int)floor(y / d->grid.dy);
    int gk = (int)floor(z / d->grid.dz);

    owner[i] =
        commRankOfCell(&d->comm, gi, gj, gk, d->grid.imax, d->grid.jmax, d->grid.kmax);
    ++sendCount[owner[i]];
    ++leaving;
  }

  MPI_Alltoall(sendCount, 1, MPI_INT, recvCount, 1, MPI_INT, d->comm.comm);

  int totalSend = 0, totalRecv = 0;
  for (int r = 0; r < size; r++) {
    sendDispl[r] = totalSend;
    recvDispl[r] = totalRecv;
    totalSend += sendCount[r];
    totalRecv += recvCount[r];
  }

  ParticleType *outbox =
      malloc((size_t)(totalSend > 0 ? totalSend : 1) * sizeof(ParticleType));
  ParticleType *inbox =
      malloc((size_t)(totalRecv > 0 ? totalRecv : 1) * sizeof(ParticleType));

  int *at = calloc((size_t)size, sizeof(int));

  for (int i = 0; i < t->count; i++) {
    if (owner[i] < 0) {
      continue;
    }
    outbox[sendDispl[owner[i]] + at[owner[i]]] = t->pool[i];
    ++at[owner[i]];
  }

  /* Drop the departed from this rank's pool, back to front so the swap with the
   * last element cannot move one that has not been looked at yet. */
  for (int i = t->count - 1; i >= 0; i--) {
    if (owner[i] >= 0) {
      removeAt(t, i);
    }
  }

  MPI_Alltoallv(outbox,
      sendCount,
      sendDispl,
      particleDatatype(),
      inbox,
      recvCount,
      recvDispl,
      particleDatatype(),
      d->comm.comm);

  reserve(t, t->count + totalRecv);
  for (int i = 0; i < totalRecv; i++) {
    t->pool[t->count++] = inbox[i];
  }

  (void)leaving;

  free(sendCount);
  free(recvCount);
  free(sendDispl);
  free(recvDispl);
  free(owner);
  free(outbox);
  free(inbox);
  free(at);
#else
  (void)t;
  (void)d;
#endif
}

/* ---------------------------------------------------------------------- */

/* A ParaView file-series index next to the particle files. The legacy reader
 * would otherwise number the files 0, 1, 2, ...; opening this file instead
 * animates them in simulation time, so they line up with anything else keyed
 * to it. Rewritten whole after every file, so it is complete if a run stops. */
static void writeSeries(ParticleTracerType *t, double time)
{
  if (t->writeIndex >= t->seriesCapacity) {
    int capacity = t->seriesCapacity > 0 ? 2 * t->seriesCapacity : 64;
    double *grown = realloc(t->seriesTimes, (size_t)capacity * sizeof(double));

    if (grown == NULL) {
      fprintf(stderr, "particle tracing: cannot record the time of file %d\n",
          t->writeIndex);
      return;
    }
    t->seriesTimes    = grown;
    t->seriesCapacity = capacity;
  }
  t->seriesTimes[t->writeIndex] = time;

  char path[256];
  snprintf(path, sizeof(path), "%s/particles.vtk.series", PARTICLE_DIR);

  FILE *fp = fopen(path, "w");

  if (fp == NULL) {
    fprintf(stderr, "particle tracing: cannot write %s\n", path);
    return;
  }

  fprintf(fp, "{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n");
  for (int i = 0; i <= t->writeIndex; i++) {
    fprintf(fp,
        "    { \"name\" : \"particles_%05d.vtk\", \"time\" : %.10g }%s\n",
        i,
        t->seriesTimes[i],
        i < t->writeIndex ? "," : "");
  }
  fprintf(fp, "  ]\n}\n");

  fclose(fp);
}

static void writeParticles(ParticleTracerType *t, Discretization *d, double time)
{
  int rank          = d->comm.rank;
  int size          = d->comm.size;

  int local         = t->count;
  int *counts       = NULL;
  int *displs       = NULL;
  ParticleType *all = NULL;
  int total         = local;

#if defined(_MPI)
  if (rank == 0) {
    counts = malloc((size_t)size * sizeof(int));
    displs = malloc((size_t)size * sizeof(int));
  }

  MPI_Gather(&local, 1, MPI_INT, counts, 1, MPI_INT, 0, d->comm.comm);

  if (rank == 0) {
    total = 0;
    for (int r = 0; r < size; r++) {
      displs[r] = total;
      total += counts[r];
    }
    all = malloc((size_t)(total > 0 ? total : 1) * sizeof(ParticleType));
  }

  MPI_Gatherv(t->pool,
      local,
      particleDatatype(),
      all,
      counts,
      displs,
      particleDatatype(),
      0,
      d->comm.comm);
#else
  (void)size;
  all = t->pool;
#endif

  if (rank == 0) {
    mkdir(PARTICLE_DIR, 0755);

    char path[256];
    snprintf(path, sizeof(path), "%s/particles_%05d.vtk", PARTICLE_DIR, t->writeIndex);

    FILE *fp = fopen(path, "w");

    if (fp == NULL) {
      fprintf(stderr, "particle tracing: cannot write %s\n", path);
    } else {
      /* Legacy VTK polydata, which is what the field output already is, so a
       * particle file and the fields for the same time open together. */
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "particles\n");
      fprintf(fp, "ASCII\n");
      fprintf(fp, "DATASET POLYDATA\n");
      fprintf(fp, "POINTS %d double\n", total);

      for (int i = 0; i < total; i++) {
        fprintf(fp, "%.10e %.10e %.10e\n", all[i].x, all[i].y, all[i].z);
      }

      fprintf(fp, "VERTICES %d %d\n", total, 2 * total);
      for (int i = 0; i < total; i++) {
        fprintf(fp, "1 %d\n", i);
      }

      fclose(fp);
    }

    writeSeries(t, time);
  }

#if defined(_MPI)
  if (rank == 0) {
    free(counts);
    free(displs);
    free(all);
  }
#endif

  ++t->writeIndex;
}

/* ---------------------------------------------------------------------- */

void particleTracerInit(ParticleTracerType *t, Discretization *d, Parameter *p)
{
  memset(t, 0, sizeof(*t));

  t->enabled = (p->numberOfParticles > 0);

  if (!t->enabled) {
    return;
  }

  t->perBatch     = p->numberOfParticles;
  t->startTime    = p->startTime;
  t->injectPeriod = p->injectTimePeriod;
  t->writePeriod  = p->writeTimePeriod;
  t->x1           = p->x1;
  t->y1           = p->y1;
  t->z1           = p->z1;
  t->x2           = p->x2;
  t->y2           = p->y2;
  t->z2           = p->z2;

  t->lastInject   = p->startTime - p->injectTimePeriod;
  t->lastWrite    = p->startTime - p->writeTimePeriod;

  t->xLo          = d->iOffset * d->grid.dx;
  t->xHi          = (d->iOffset + d->comm.imaxLocal) * d->grid.dx;
  t->yLo          = d->jOffset * d->grid.dy;
  t->yHi          = (d->jOffset + d->comm.jmaxLocal) * d->grid.dy;
  t->zLo          = d->kOffset * d->grid.dz;
  t->zHi          = (d->kOffset + d->comm.kmaxLocal) * d->grid.dz;

  /* One batch to start with; the pool grows when it has to. Sizing it from the
   * grid instead -- which is what the two-dimensional tracer does -- would
   * reserve millions of slots in three dimensions for a pool that holds a few
   * thousand. */
  reserve(t, t->perBatch);

  if (commIsMaster(&d->comm)) {
    printf("Particle tracing:\n");
    printf("\t%d particles per batch from t = %g every %g\n",
        t->perBatch,
        t->startTime,
        t->injectPeriod);
    printf("\tseed region (%g, %g, %g) to (%g, %g, %g)\n",
        t->x1,
        t->y1,
        t->z1,
        t->x2,
        t->y2,
        t->z2);
    printf("\twriting to %s every %g\n", PARTICLE_DIR, t->writePeriod);
  }
}

void particleTracerStep(ParticleTracerType *t, Discretization *d, double time)
{
  if (!t->enabled || time < t->startTime) {
    return;
  }

  if (t->injectPeriod > 0.0 && (time - t->lastInject) >= t->injectPeriod) {
    particleTracerInject(t, d);
    t->lastInject = time;
  }

  particleTracerAdvance(t, d, d->dt);
  particleTracerMigrate(t, d);

  if (t->writePeriod > 0.0 && (time - t->lastWrite) >= t->writePeriod) {
    writeParticles(t, d, time);
    t->lastWrite = time;
  }
}

void particleTracerTotals(const ParticleTracerType *t,
    Discretization *d,
    double *injected,
    double *alive,
    double *removedBoundary,
    double *removedBody)
{
  double a = (double)t->injected;
  double b = (double)t->count;
  double c = (double)t->removedBoundary;
  double e = (double)t->removedBody;

  (void)d;
  commReduceAll(&a, SUM);
  commReduceAll(&b, SUM);
  commReduceAll(&c, SUM);
  commReduceAll(&e, SUM);

  *injected        = a;
  *alive           = b;
  *removedBoundary = c;
  *removedBody     = e;
}

void particleTracerFinalize(ParticleTracerType *t, Discretization *d)
{
  if (!t->enabled) {
    return;
  }

  double injected, alive, atBoundary, atBody;
  particleTracerTotals(t, d, &injected, &alive, &atBoundary, &atBody);

  if (commIsMaster(&d->comm)) {
    printf("Particles: %.0f injected, %.0f still in the domain, %.0f left through a "
           "boundary, %.0f stopped at the body\n",
        injected,
        alive,
        atBoundary,
        atBody);
  }

  free(t->pool);
  free(t->seriesTimes);
  t->pool     = NULL;
  t->capacity = 0;
  t->count    = 0;
}
