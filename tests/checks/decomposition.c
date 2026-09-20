/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Checks that the domain decomposition and the global offsets derived from it
 * describe the same partition:
 *
 *  - every rank's offsets are in range and its block fits inside the domain;
 *  - the blocks tile the domain exactly, with no gap and no overlap;
 *  - rank 0 sits at the origin, which is the only thing a serial build can
 *    report and therefore the thing both builds have to agree on.
 *
 * The tiling check is what would catch a coords/dims pairing that named one
 * axis and indexed another.
 */
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "coloring.h"

#define IMAX 37
#define JMAX 19
#define KMAX 23

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  commPartition(&comm, KMAX, JMAX, IMAX);

  CHECK_BEGIN("decomposition");

  int offsets[NDIMS];
  memset(offsets, 0x7f, sizeof(offsets));
  commGetOffsets(&comm, offsets, KMAX, JMAX, IMAX);

  int iOffset = offsets[IDIM];
  int jOffset = offsets[JDIM];
  int kOffset = offsets[KDIM];

  CHECK_TRUE(iOffset >= 0 && iOffset + comm.imaxLocal <= IMAX,
      "x block [%d,%d) does not fit in [0,%d)",
      iOffset,
      iOffset + comm.imaxLocal,
      IMAX);
  CHECK_TRUE(jOffset >= 0 && jOffset + comm.jmaxLocal <= JMAX,
      "y block [%d,%d) does not fit in [0,%d)",
      jOffset,
      jOffset + comm.jmaxLocal,
      JMAX);
  CHECK_TRUE(kOffset >= 0 && kOffset + comm.kmaxLocal <= KMAX,
      "z block [%d,%d) does not fit in [0,%d)",
      kOffset,
      kOffset + comm.kmaxLocal,
      KMAX);

  if (commIsMaster(&comm)) {
    CHECK_TRUE(iOffset == 0 && jOffset == 0 && kOffset == 0,
        "rank 0 is not at the origin, offsets are (%d,%d,%d)",
        iOffset,
        jOffset,
        kOffset);
  }

  /* Mark every global cell this rank owns, then sum the marks across ranks.
   * A partition that tiles the domain marks every cell exactly once. */
  size_t cells  = (size_t)IMAX * JMAX * KMAX;
  double *marks = malloc(cells * sizeof(double));

  for (size_t i = 0; i < cells; i++) {
    marks[i] = 0.0;
  }

  for (int k = 0; k < comm.kmaxLocal; k++) {
    for (int j = 0; j < comm.jmaxLocal; j++) {
      for (int i = 0; i < comm.imaxLocal; i++) {
        size_t idx = (size_t)(k + kOffset) * IMAX * JMAX + (size_t)(j + jOffset) * IMAX +
                     (size_t)(i + iOffset);
        marks[idx] += 1.0;
      }
    }
  }

#if defined(_MPI)
  double *total = malloc(cells * sizeof(double));
  MPI_Allreduce(marks, total, (int)cells, MPI_DOUBLE, MPI_SUM, comm.comm);
  memcpy(marks, total, cells * sizeof(double));
  free(total);
#endif

  int unclaimed = 0, duplicated = 0;

  for (size_t i = 0; i < cells; i++) {
    if (marks[i] < 0.5) {
      ++unclaimed;
    } else if (marks[i] > 1.5) {
      ++duplicated;
    }
  }

  CHECK_TRUE(unclaimed == 0, "%d global cells are owned by no rank", unclaimed);
  CHECK_TRUE(duplicated == 0, "%d global cells are owned by more than one rank", duplicated);

  /*
   * The red-black colour of a global cell must not depend on which rank holds
   * it. Record the colour this rank assigns to each cell it owns, encoded so
   * that a disagreement between decompositions is visible in a single number,
   * and compare the total against what the global position alone predicts.
   *
   * The colouring used to start from the local index, so a cell's colour
   * flipped whenever a subdomain boundary landed on an odd offset.
   */
  int mismatches = 0;

  for (int k = 1; k < comm.kmaxLocal + 1; k++) {
    for (int j = 1; j < comm.jmaxLocal + 1; j++) {
      for (int i = 1; i < comm.imaxLocal + 1; i++) {
        int local  = colorOf(i, j, k, iOffset, jOffset, kOffset);
        /* The same cell, addressed as though one rank owned everything. */
        int global = colorOf(i + iOffset, j + jOffset, k + kOffset, 0, 0, 0);

        if (local != global) {
          ++mismatches;
        }
      }
    }
  }

  CHECK_TRUE(mismatches == 0,
      "%d cells are coloured differently than their global position implies",
      mismatches);

  free(marks);

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures (rank %d of %d, local %dx%dx%d at %d,%d,%d)\n",
      CheckName,
      CheckCount,
      failures,
      comm.rank,
      comm.size,
      comm.imaxLocal,
      comm.jmaxLocal,
      comm.kmaxLocal,
      iOffset,
      jOffset,
      kOffset);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
