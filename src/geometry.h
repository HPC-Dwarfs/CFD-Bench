/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __GEOMETRY_H_
#define __GEOMETRY_H_

#include "comm.h"

/*
 * Obstacle geometry, as the rest of the solver sees it.
 *
 * Geometry is carried by four fields co-located with the staggered unknowns:
 *
 *   Lambda(i,j,k)  volume fraction of the cell, co-located with p
 *   Ax(i,j,k)      aperture of the face between (i,j,k) and (i+1,j,k), with u
 *   Ay(i,j,k)      aperture of the face between (i,j,k) and (i,j+1,k), with v
 *   Az(i,j,k)      aperture of the face between (i,j,k) and (i,j,k+1), with w
 *
 * 0 means fully solid and 1 fully fluid. Values in between are reserved for a
 * later phase; a producer here emits only 0 or 1.
 *
 * These four fields are the only channel between geometry and the rest of the
 * solver. Nothing downstream reads a file, a shape or a cell type.
 */

typedef enum {
  GEOMETRY_NONE = 0, /* obstacle-free */
  GEOMETRY_VOXEL,    /* a voxel volume file, thresholded */
  GEOMETRY_SPHERE,   /* analytic sphere */
  GEOMETRY_CYLINDER, /* analytic cylinder, axis-aligned */
  GEOMETRY_BOX,      /* analytic axis-aligned box */
  GEOMETRY_PLATE     /* analytic zero-thickness plate, axis-aligned */
} GeometryKindType;

/* Which axis an analytic body is aligned with. */
typedef enum { AXIS_X = 0, AXIS_Y, AXIS_Z } GeometryAxisType;

typedef struct {
  GeometryKindType kind;
  const char *file;                     /* GEOMETRY_VOXEL */
  double xCenter, yCenter, zCenter;     /* SPHERE, CYLINDER */
  double radius;                        /* SPHERE, CYLINDER */
  GeometryAxisType axis;                /* CYLINDER, PLATE */
  double x0, y0, z0, x1, y1, z1;        /* BOX, PLATE */
} GeometrySpecType;

/*
 * What a producer needs to know about the grid it is filling. Deliberately not
 * Discretization or CommType: a producer depends on neither, which keeps it
 * callable from a check driver with a grid that no solver owns.
 */
typedef struct {
  int imaxLocal, jmaxLocal, kmaxLocal;  /* interior cells held by this rank */
  int iOffset, jOffset, kOffset;        /* global index of this rank's first
                                         * interior cell */
  int imax, jmax, kmax;                 /* interior cells in the whole domain */
  double dx, dy, dz;
  double xlength, ylength, zlength;
} GeometryDomainType;

/*
 * Fill Ax, Ay, Az and Lambda over the local subdomain including the halo layer,
 * so that no exchange is needed to make them consistent. Each array is
 * (imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2) doubles.
 *
 * Aborts the run if the spec names geometry that cannot be honoured.
 */
extern void geometryProduce(const GeometrySpecType *spec,
    const GeometryDomainType *domain,
    double *Ax,
    double *Ay,
    double *Az,
    double *Lambda);

/*
 * Verify that the fluid part of the domain forms exactly one region under face
 * connectivity -- six neighbours per cell. Each sealed fluid pocket adds an
 * independent constant to the pressure null space, which the solvers here do
 * not carry.
 *
 * Unlike the producers this needs the communicator: connectivity is a global
 * property and a region can span any number of ranks.
 *
 * Returns the number of fluid regions found; aborts before returning when that
 * is not one and abortOnFailure is set.
 */
extern int geometryValidateConnectivity(CommType *comm,
    const GeometryDomainType *domain,
    const double *Ax,
    const double *Ay,
    const double *Az,
    const double *Lambda,
    int abortOnFailure);

/* Build a spec from the `geometryFile` parameter value. A name of the form
 * "sphere:xc,yc,zc,r", "cylinder-z:xc,yc,r", "box:x0,y0,z0,x1,y1,z1" or
 * "plate-x:x,y0,z0,y1,z1" selects an analytic body; anything else is a path to
 * a voxel volume; NULL or empty means obstacle-free. */
extern void geometryParseSpec(GeometrySpecType *spec, const char *name);

/* One line describing the geometry in use, for the run header. */
extern void geometryPrintHeader(
    const GeometrySpecType *spec, const GeometryDomainType *domain);

/* Checksum and dimensions of the volume last loaded, 0 when none was. */
extern unsigned long long geometryChecksum(void);
extern void geometryVolumeSize(int *nx, int *ny, int *nz);

#endif // __GEOMETRY_H_
