/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __GEOMETRY_VOXEL_H_
#define __GEOMETRY_VOXEL_H_

#include <stddef.h>

#include "geometry.h"

/*
 * The voxel volume format: the three-dimensional counterpart of the binary PGM
 * the two-dimensional solver reads.
 *
 *   P5V
 *   # optional comment lines, anywhere in the header
 *   <nx> <ny> <nz>
 *   <maxval>
 *   <nx*ny*nz raw bytes, x fastest, then y, then z>
 *
 * A voxel below 128 is solid, 128 or above is fluid. Exactly one whitespace
 * character separates the maximum value from the data.
 *
 * Voxel (vx, vy, vz) covers the physical box
 *
 *   x in [vx, vx+1) * xlength / nx,  and likewise for y and z,
 *
 * so the index order matches the axis order and the origin is the domain
 * origin. No axis is flipped; a generator writes x fastest and counts every
 * axis upward, the same direction the grid is indexed.
 *
 * A rank reads only the voxels covering its own subdomain plus the sampling
 * margin, by computing byte offsets from the header and seeking. Memory is
 * O(local voxels), which matters because a volume at the required sampling
 * density is 64 bytes per grid cell -- a benchmark-sized grid would be
 * gigabytes if every rank held all of it.
 */

/* Open the volume, validate it against the grid, and read this rank's slab.
 * Aborts on any problem, naming the file and the reason. */
extern void geometryVoxelLoad(const char *path, const GeometryDomainType *domain);

/* Is the point inside a solid voxel? Valid only after geometryVoxelLoad, and
 * only for points inside this rank's slab. */
extern int geometryVoxelIsSolid(double x, double y, double z);

/* Bytes this rank allocated for its slab, for the per-rank memory check. */
extern size_t geometryVoxelSlabBytes(void);

extern void geometryVoxelFree(void);

#ifdef TEST
/* Exposed for the check drivers: read one voxel straight from the slab, and
 * report the header without sampling. */
extern int geometryVoxelAt(int vx, int vy, int vz);
extern void geometryVoxelSlabOrigin(int *ox, int *oy, int *oz);
#endif

#endif // __GEOMETRY_VOXEL_H_
