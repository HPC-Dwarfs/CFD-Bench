/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Produces geometry under the real decomposition and reports quantities that
 * must not depend on how the domain was divided:
 *
 *   - a checksum of each aperture field, keyed by global cell index so that it
 *     is independent of the order cells are visited in and of which rank holds
 *     them. Summed as exact integers, so the reduction order cannot change it.
 *   - the number of connected fluid regions.
 *   - the bytes each rank allocated for its voxel slab, which has to scale with
 *     the subdomain rather than with the whole volume.
 *
 * tests/check-geometry-ranks.sh runs this at several rank counts and compares.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "geometry.h"
#include "geometry-voxel.h"

#define IMAX 24
#define JMAX 24
#define KMAX 24

#define VOLUME "tests/geom/sphere.vox"

#define G(v, i, j, k)                                                                    \
  v[(size_t)(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (size_t)(j) * (imaxLocal + 2) +    \
      (size_t)(i)]

/*
 * A checksum that depends on which global cells are solid and on nothing else.
 * Each contribution is a small integer, so the sum is exact in double precision
 * and the reduction order cannot matter.
 */
static double fieldChecksum(const double *field, int imaxLocal, int jmaxLocal,
    int kmaxLocal, int iOffset, int jOffset, int kOffset)
{
  double sum = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (G(field, i, j, k) != 0.0) {
          continue;
        }

        long gi  = i - 1 + iOffset;
        long gj  = j - 1 + jOffset;
        long gk  = k - 1 + kOffset;
        long idx = (gk * JMAX + gj) * IMAX + gi;

        sum += (double)(1 + (idx * 2654435761L) % 1000003L);
      }
    }
  }

  commReduceAll(&sum, SUM);
  return sum;
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  commPartition(&comm, KMAX, JMAX, IMAX);

  CHECK_BEGIN("geometry-mpi");

  int imaxLocal      = comm.imaxLocal;
  int jmaxLocal      = comm.jmaxLocal;
  int kmaxLocal      = comm.kmaxLocal;

  int offsets[NDIMS] = { 0, 0, 0 };
  commGetOffsets(&comm, offsets, KMAX, JMAX, IMAX);

  GeometryDomainType domain;
  memset(&domain, 0, sizeof(domain));
  domain.imaxLocal = imaxLocal;
  domain.jmaxLocal = jmaxLocal;
  domain.kmaxLocal = kmaxLocal;
  domain.iOffset   = offsets[IDIM];
  domain.jOffset   = offsets[JDIM];
  domain.kOffset   = offsets[KDIM];
  domain.imax      = IMAX;
  domain.jmax      = JMAX;
  domain.kmax      = KMAX;
  domain.xlength = domain.ylength = domain.zlength = 4.0;
  domain.dx = domain.dy = domain.dz = 4.0 / IMAX;

  size_t size    = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);
  double *Ax     = calloc(size, sizeof(double));
  double *Ay     = calloc(size, sizeof(double));
  double *Az     = calloc(size, sizeof(double));
  double *Lambda = calloc(size, sizeof(double));

  GeometrySpecType spec;
  geometryParseSpec(&spec, VOLUME);
  geometryVoxelLoad(spec.file, &domain);
  geometryProduce(&spec, &domain, Ax, Ay, Az, Lambda);

  double sumAx = fieldChecksum(Ax, imaxLocal, jmaxLocal, kmaxLocal, domain.iOffset,
      domain.jOffset, domain.kOffset);
  double sumAy = fieldChecksum(Ay, imaxLocal, jmaxLocal, kmaxLocal, domain.iOffset,
      domain.jOffset, domain.kOffset);
  double sumAz = fieldChecksum(Az, imaxLocal, jmaxLocal, kmaxLocal, domain.iOffset,
      domain.jOffset, domain.kOffset);
  double sumL  = fieldChecksum(Lambda, imaxLocal, jmaxLocal, kmaxLocal, domain.iOffset,
      domain.jOffset, domain.kOffset);

  int regions  = geometryValidateConnectivity(&comm, &domain, Ax, Ay, Az, Lambda, 0);

  CHECK_TRUE(regions == 1, "the sphere leaves %d fluid regions, expected 1", regions);

  /* The slab has to be a fraction of the whole volume once the domain is
   * divided; on one rank it is naturally the whole thing. */
  size_t slab = geometryVoxelSlabBytes();
  int nx, ny, nz;
  geometryVolumeSize(&nx, &ny, &nz);
  size_t whole = (size_t)nx * ny * nz;

  if (comm.size > 1) {
    CHECK_TRUE(slab < whole,
        "a rank of %d allocated %zu bytes for its slab, the whole volume is %zu",
        comm.size,
        slab,
        whole);
  }

  if (commIsMaster(&comm)) {
    printf("APERTURE-CHECKSUM Ax %.0f Ay %.0f Az %.0f Lambda %.0f\n",
        sumAx,
        sumAy,
        sumAz,
        sumL);
    printf("REGIONS %d\n", regions);
  }

  printf("SLAB rank %d of %d: %zu bytes of %zu in the whole volume\n",
      comm.rank,
      comm.size,
      slab,
      whole);

  geometryVoxelFree();
  free(Ax);
  free(Ay);
  free(Az);
  free(Lambda);

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
