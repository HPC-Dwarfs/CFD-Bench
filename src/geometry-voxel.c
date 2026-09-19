/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry-voxel.h"

/* Minimum voxels per grid cell per direction. Below this the represented
 * geometry changes under grid refinement and a resolution sweep stops
 * measuring convergence. */
#define MIN_VOXELS_PER_CELL 4

#define SOLID_THRESHOLD 128

/* Aspect ratios are compared with this relative tolerance before warning. */
#define ASPECT_TOLERANCE 0.01

typedef struct {
  int nx, ny, nz;                 /* voxels in the whole volume */
  double xlength, ylength, zlength;
  /* this rank's slab, in voxel indices: [ox, ox+sx) etc. */
  int ox, oy, oz;
  int sx, sy, sz;
  unsigned char *data;
  size_t bytes;
  unsigned long long checksum;
} VoxelVolumeType;

static VoxelVolumeType Volume;

static void fail(const char *path, const char *reason)
{
  fprintf(stderr, "geometryFile '%s': %s\n", path, reason);
  exit(EXIT_FAILURE);
}

/* Read the next header token, skipping whitespace and comment lines. */
static int nextToken(FILE *fp, char *buf, size_t n)
{
  int c;

  for (;;) {
    c = fgetc(fp);

    if (c == EOF) {
      return 0;
    }

    if (c == '#') {
      while (c != EOF && c != '\n') {
        c = fgetc(fp);
      }
      continue;
    }

    if (!isspace(c)) {
      break;
    }
  }

  size_t i = 0;
  while (c != EOF && !isspace(c) && c != '#' && i + 1 < n) {
    buf[i++] = (char)c;
    c        = fgetc(fp);
  }
  buf[i] = '\0';

  if (c == '#') {
    ungetc(c, fp);
  }

  return i > 0;
}

/* FNV-1a over the whole file, streamed so that no rank has to hold it. */
static unsigned long long fileChecksum(FILE *fp)
{
  unsigned long long hash = 1469598103934665603ULL;
  unsigned char buf[65536];
  size_t got;

  rewind(fp);

  while ((got = fread(buf, 1, sizeof(buf), fp)) > 0) {
    for (size_t i = 0; i < got; i++) {
      hash ^= buf[i];
      hash *= 1099511628211ULL;
    }
  }

  return hash;
}

void geometryVoxelLoad(const char *path, const GeometryDomainType *domain)
{
  FILE *fp = fopen(path, "rb");

  if (fp == NULL) {
    fail(path, "cannot be opened");
  }

  char tok[64];

  if (!nextToken(fp, tok, sizeof(tok)) || strcmp(tok, "P5V") != 0) {
    fail(path, "does not start with the P5V magic");
  }

  int dims[3];
  for (int d = 0; d < 3; d++) {
    if (!nextToken(fp, tok, sizeof(tok))) {
      fail(path, "header ends before the three voxel counts");
    }
    dims[d] = atoi(tok);
    if (dims[d] <= 0) {
      fail(path, "has a non-positive voxel count");
    }
  }

  if (!nextToken(fp, tok, sizeof(tok))) {
    fail(path, "header ends before the maximum value");
  }

  if (atoi(tok) != 255) {
    fail(path, "has a maximum value other than 255");
  }

  /* nextToken stops on the whitespace that ends a token and has already taken
   * it from the stream, so the data starts exactly here. */
  long dataOffset = ftell(fp);

  Volume.nx       = dims[0];
  Volume.ny       = dims[1];
  Volume.nz       = dims[2];
  Volume.xlength  = domain->xlength;
  Volume.ylength  = domain->ylength;
  Volume.zlength  = domain->zlength;

  /* The payload has to be all there before anything is believed about it. */
  if (fseek(fp, 0, SEEK_END) != 0) {
    fail(path, "cannot be seeked");
  }

  long fileSize  = ftell(fp);
  long wantBytes = (long)Volume.nx * Volume.ny * Volume.nz;

  if (fileSize - dataOffset < wantBytes) {
    char msg[160];
    snprintf(msg,
        sizeof(msg),
        "is truncated: the header declares %d x %d x %d voxels, which needs %ld bytes, "
        "but only %ld follow the header",
        Volume.nx,
        Volume.ny,
        Volume.nz,
        wantBytes,
        fileSize - dataOffset);
    fail(path, msg);
  }

  /* Resolution: the represented geometry must not change as the grid is
   * refined, which needs the volume to out-resolve the grid. */
  double perCellX = (double)Volume.nx / domain->imax;
  double perCellY = (double)Volume.ny / domain->jmax;
  double perCellZ = (double)Volume.nz / domain->kmax;

  if (perCellX < MIN_VOXELS_PER_CELL || perCellY < MIN_VOXELS_PER_CELL ||
      perCellZ < MIN_VOXELS_PER_CELL) {
    fprintf(stderr,
        "geometryFile '%s': too coarse for the grid. The volume is %d x %d x %d "
        "voxels and the grid is %d x %d x %d cells, giving %.2f x %.2f x %.2f voxels "
        "per cell; at least %d per cell per direction are required.\n",
        path,
        Volume.nx,
        Volume.ny,
        Volume.nz,
        domain->imax,
        domain->jmax,
        domain->kmax,
        perCellX,
        perCellY,
        perCellZ,
        MIN_VOXELS_PER_CELL);
    exit(EXIT_FAILURE);
  }

  /* Aspect: sampling an anisotropic volume is allowed but is almost always a
   * mistake, so say so and carry on. */
  double volXY = (double)Volume.nx / Volume.ny;
  double domXY = domain->xlength / domain->ylength;
  double volXZ = (double)Volume.nx / Volume.nz;
  double domXZ = domain->xlength / domain->zlength;

  if (fabs(volXY - domXY) > ASPECT_TOLERANCE * domXY ||
      fabs(volXZ - domXZ) > ASPECT_TOLERANCE * domXZ) {
    fprintf(stderr,
        "geometryFile '%s': aspect ratio mismatch. The volume is %.4f : %.4f (x:y, "
        "x:z) and the domain is %.4f : %.4f; the volume will be sampled "
        "anisotropically.\n",
        path,
        volXY,
        volXZ,
        domXY,
        domXZ);
  }

  Volume.checksum = fileChecksum(fp);

  /*
   * This rank's slab: the voxels covering its interior and halo cells, plus one
   * voxel of margin so that a sub-cell sample landing exactly on a boundary
   * cannot fall outside.
   */
  double vx = (double)Volume.nx / domain->imax;
  double vy = (double)Volume.ny / domain->jmax;
  double vz = (double)Volume.nz / domain->kmax;

  int lo[3], hi[3];
  lo[0]    = (int)((domain->iOffset - 1) * vx) - 1;
  hi[0]    = (int)((domain->iOffset + domain->imaxLocal + 1) * vx) + 1;
  lo[1]    = (int)((domain->jOffset - 1) * vy) - 1;
  hi[1]    = (int)((domain->jOffset + domain->jmaxLocal + 1) * vy) + 1;
  lo[2]    = (int)((domain->kOffset - 1) * vz) - 1;
  hi[2]    = (int)((domain->kOffset + domain->kmaxLocal + 1) * vz) + 1;

  int n[3] = { Volume.nx, Volume.ny, Volume.nz };

  for (int d = 0; d < 3; d++) {
    if (lo[d] < 0) {
      lo[d] = 0;
    }
    if (hi[d] > n[d] - 1) {
      hi[d] = n[d] - 1;
    }
  }

  Volume.ox    = lo[0];
  Volume.oy    = lo[1];
  Volume.oz    = lo[2];
  Volume.sx    = hi[0] - lo[0] + 1;
  Volume.sy    = hi[1] - lo[1] + 1;
  Volume.sz    = hi[2] - lo[2] + 1;

  Volume.bytes = (size_t)Volume.sx * Volume.sy * Volume.sz;
  Volume.data  = malloc(Volume.bytes);

  if (Volume.data == NULL) {
    fail(path, "slab does not fit in memory");
  }

  /* One positioned read per contiguous voxel row. No rank ever holds more than
   * its own slab, and the rows a given cell is sampled from are the same
   * whatever the decomposition, so the apertures come out bit-identical. */
  for (int z = 0; z < Volume.sz; z++) {
    for (int y = 0; y < Volume.sy; y++) {
      long off = dataOffset +
                 ((long)(z + Volume.oz) * Volume.ny + (y + Volume.oy)) * Volume.nx +
                 Volume.ox;

      if (fseek(fp, off, SEEK_SET) != 0) {
        fail(path, "cannot be seeked to a voxel row");
      }

      unsigned char *row = Volume.data + ((size_t)z * Volume.sy + y) * Volume.sx;

      if (fread(row, 1, (size_t)Volume.sx, fp) != (size_t)Volume.sx) {
        fail(path, "ends in the middle of a voxel row");
      }
    }
  }

  fclose(fp);
}

int geometryVoxelIsSolid(double x, double y, double z)
{
  if (Volume.data == NULL) {
    return 0;
  }

  int vx = (int)(x / Volume.xlength * Volume.nx);
  int vy = (int)(y / Volume.ylength * Volume.ny);
  int vz = (int)(z / Volume.zlength * Volume.nz);

  /* Outside the volume is outside the domain, which is the boundary
   * condition's business rather than the obstacle's. */
  if (vx < 0 || vx >= Volume.nx || vy < 0 || vy >= Volume.ny || vz < 0 ||
      vz >= Volume.nz) {
    return 0;
  }

  int lx = vx - Volume.ox;
  int ly = vy - Volume.oy;
  int lz = vz - Volume.oz;

  if (lx < 0 || lx >= Volume.sx || ly < 0 || ly >= Volume.sy || lz < 0 ||
      lz >= Volume.sz) {
    return 0;
  }

  return Volume.data[((size_t)lz * Volume.sy + ly) * Volume.sx + lx] < SOLID_THRESHOLD;
}

size_t geometryVoxelSlabBytes(void)
{
  return Volume.data == NULL ? 0 : Volume.bytes;
}

void geometryVoxelFree(void)
{
  free(Volume.data);
  Volume.data  = NULL;
  Volume.bytes = 0;
}

unsigned long long geometryChecksum(void)
{
  return Volume.checksum;
}

void geometryVolumeSize(int *nx, int *ny, int *nz)
{
  *nx = Volume.nx;
  *ny = Volume.ny;
  *nz = Volume.nz;
}

#ifdef TEST
int geometryVoxelAt(int vx, int vy, int vz)
{
  int lx = vx - Volume.ox;
  int ly = vy - Volume.oy;
  int lz = vz - Volume.oz;

  if (Volume.data == NULL || lx < 0 || lx >= Volume.sx || ly < 0 || ly >= Volume.sy ||
      lz < 0 || lz >= Volume.sz) {
    return -1;
  }

  return Volume.data[((size_t)lz * Volume.sy + ly) * Volume.sx + lx];
}

void geometryVoxelSlabOrigin(int *ox, int *oy, int *oz)
{
  *ox = Volume.ox;
  *oy = Volume.oy;
  *oz = Volume.oz;
}
#endif
