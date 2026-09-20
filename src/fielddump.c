/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "fielddump.h"

#define NFIELDS 4

#define L(v, i, j, k)                                                                    \
  v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

const char *fieldDumpPath(void)
{
  return getenv("NUSIF_FIELD_DUMP");
}

/* Copy this rank's interior of one field into a contiguous buffer, in the same
 * k-slowest order the global array uses. */
static void packInterior(
    const double *src, double *dst, int imaxLocal, int jmaxLocal, int kmaxLocal)
{
  size_t idx = 0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        dst[idx++] = L(src, i, j, k);
      }
    }
  }
}

/* Place a contiguous subdomain block at its global position. */
static void placeBlock(double *global,
    const double *block,
    int imax,
    int jmax,
    int iOffset,
    int jOffset,
    int kOffset,
    int ni,
    int nj,
    int nk)
{
  size_t idx = 0;

  for (int k = 0; k < nk; k++) {
    for (int j = 0; j < nj; j++) {
      double *row = global + (size_t)(k + kOffset) * imax * jmax +
                    (size_t)(j + jOffset) * imax + iOffset;
      memcpy(row, block + idx, (size_t)ni * sizeof(double));
      idx += ni;
    }
  }
}

void fieldDumpWrite(CommType *comm,
    const Grid *grid,
    const char *filename,
    const double *p,
    const double *u,
    const double *v,
    const double *w)
{
  int imaxLocal              = comm->imaxLocal;
  int jmaxLocal              = comm->jmaxLocal;
  int kmaxLocal              = comm->kmaxLocal;

  int imax                   = grid->imax;
  int jmax                   = grid->jmax;
  int kmax                   = grid->kmax;

  const double *src[NFIELDS] = { p, u, v, w };

  size_t localCount          = (size_t)imaxLocal * jmaxLocal * kmaxLocal;
  size_t globalCount         = (size_t)imax * jmax * kmax;

  double *block              = allocate(ARRAY_ALIGNMENT, localCount * sizeof(double));
  double *global             = NULL;

  if (commIsMaster(comm)) {
    global = allocate(ARRAY_ALIGNMENT, NFIELDS * globalCount * sizeof(double));
  }

  int offsets[NDIMS] = { 0, 0, 0 };
  commGetOffsets(comm, offsets, kmax, jmax, imax);

  for (int f = 0; f < NFIELDS; f++) {
    packInterior(src[f], block, imaxLocal, jmaxLocal, kmaxLocal);

    if (commIsMaster(comm)) {
      placeBlock(global + (size_t)f * globalCount,
          block,
          imax,
          jmax,
          offsets[IDIM],
          offsets[JDIM],
          offsets[KDIM],
          imaxLocal,
          jmaxLocal,
          kmaxLocal);
    }

#if defined(_MPI)
    if (!commIsMaster(comm)) {
      int meta[6] = {
        offsets[IDIM], offsets[JDIM], offsets[KDIM], imaxLocal, jmaxLocal, kmaxLocal
      };
      MPI_Send(meta, 6, MPI_INT, 0, 0, comm->comm);
      MPI_Send(block, (int)localCount, MPI_DOUBLE, 0, 1, comm->comm);
    } else {
      for (int r = 1; r < comm->size; r++) {
        int meta[6];
        MPI_Recv(meta, 6, MPI_INT, r, 0, comm->comm, MPI_STATUS_IGNORE);

        size_t count = (size_t)meta[3] * meta[4] * meta[5];
        double *recv = allocate(ARRAY_ALIGNMENT, count * sizeof(double));
        MPI_Recv(recv, (int)count, MPI_DOUBLE, r, 1, comm->comm, MPI_STATUS_IGNORE);

        placeBlock(global + (size_t)f * globalCount,
            recv,
            imax,
            jmax,
            meta[0],
            meta[1],
            meta[2],
            meta[3],
            meta[4],
            meta[5]);
        free(recv);
      }
    }
#endif
  }

  if (commIsMaster(comm)) {
    FILE *fp = fopen(filename, "wb");

    if (fp == NULL) {
      fprintf(stderr, "fieldDumpWrite: cannot open %s for writing\n", filename);
      exit(EXIT_FAILURE);
    }

    int32_t dims[3] = { imax, jmax, kmax };
    fwrite("NUSIFD01", 1, 8, fp);
    fwrite(dims, sizeof(int32_t), 3, fp);
    fwrite(global, sizeof(double), NFIELDS * globalCount, fp);
    fclose(fp);

    printf("Wrote field dump %s (%d x %d x %d)\n", filename, imax, jmax, kmax);
    free(global);
  }

  free(block);
}
