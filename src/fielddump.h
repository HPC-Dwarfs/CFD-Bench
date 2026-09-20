/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __FIELDDUMP_H_
#define __FIELDDUMP_H_

#include "comm.h"
#include "grid.h"

/*
 * A raw binary dump of the interior of p, u, v, w over the whole domain, for
 * comparing two runs against each other. The solver's only other output is VTK,
 * which averages the staggered velocities to cell centres; that loses exactly
 * the face values the embedded-boundary checks are about, so this writes the
 * unmodified face values instead.
 *
 * Layout, little-endian host order throughout:
 *
 *   char   magic[8]   "NUSIFD01"
 *   int32  imax, jmax, kmax
 *   double p[imax*jmax*kmax]      cell centres
 *   double u[imax*jmax*kmax]      x-faces, u(i,j,k) is the face at i
 *   double v[imax*jmax*kmax]      y-faces
 *   double w[imax*jmax*kmax]      z-faces
 *
 * Gathered to rank 0 and written there. Collective: every rank must call it.
 */
extern void fieldDumpWrite(CommType *comm,
    const Grid *grid,
    const char *filename,
    const double *p,
    const double *u,
    const double *v,
    const double *w);

/* The dump path a run was asked for, or NULL when none was. Taken from the
 * environment rather than the parameter file so that a setup file stays the
 * description of a physical problem. */
extern const char *fieldDumpPath(void);

#endif // __FIELDDUMP_H_
