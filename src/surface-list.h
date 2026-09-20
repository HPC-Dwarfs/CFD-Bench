/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __SURFACE_LIST_H_
#define __SURFACE_LIST_H_

/*
 * The cells where a geometry-free relaxation sweep gets the wrong answer.
 *
 * Two classes, and in three dimensions both scale with the obstacle's surface
 * area rather than with its volume:
 *
 *  - fluid cells with at least one closed face, where the plain stencil couples
 *    across a wall it should not;
 *  - solid cells with a fluid cell next to them, which the plain stencil would
 *    drag off the zero they are required to hold.
 *
 * A solid cell whose six neighbours are all solid needs nothing: relaxing a
 * zero neighbourhood with a zero right-hand side leaves it at exactly zero, so
 * deep interiors maintain themselves and cost nothing at runtime. That is what
 * makes the list O(surface) against O(volume) of bulk work, and it is why the
 * interior sweep can stay branch-free and geometry-free.
 *
 * Stored as arrays rather than as an array of structs, and split into the two
 * colours of the global checkerboard: entries [0, colorCount[0]) first, then
 * the rest, each segment row-major.
 *
 * The colouring matters. Relaxing the list in one pass would be a Gauss-Seidel
 * over it, so a cut cell next to another cut cell would see its neighbour's new
 * value or its old one depending on the traversal order -- and across a rank
 * boundary never the new one, which would make the iteration depend on how the
 * domain was divided. Two cut cells of the same colour are never face
 * neighbours, so within a colour the order cannot matter.
 *
 * Every entry carries its index in both forms the solvers need: the natural
 * linear index, and the compressed (ic, j, k) of the colour-split layout the
 * vectorized variant relaxes in.
 */
typedef struct {
  int count;
  int colorCount[2]; /* entries of each colour; colour 0 occupies [0, c0) */

  int *index;        /* linear index into the field arrays */
  int *ic, *jc, *kc; /* the same cell in the compressed colour-split layout */
  unsigned char *color;
  unsigned char *solid; /* 1 where the cell is an identity row */

  /* Face coefficients, aperture over spacing squared, in the order
   * east, west, north, south, top, bottom. */
  double *aE, *aW, *aN, *aS, *aT, *aB;
  double *invDiag; /* 1 / sum of the face coefficients; 0 for a solid cell */
  double *lambda;  /* volume fraction, which weights the right-hand side */

  /* Scratch: the value each listed cell held before the bulk sweep overwrote
   * it, so the correction can be computed from the right starting point. */
  double *saved;
} SurfaceListType;

/*
 * Build the list for one grid. idx2, idy2 and idz2 are 1/dx^2, 1/dy^2 and
 * 1/dz^2 of the level being built; the offsets place this rank on the global
 * checkerboard.
 */
extern void surfaceListBuild(SurfaceListType *list,
    int imaxLocal,
    int jmaxLocal,
    int kmaxLocal,
    int iOffset,
    int jOffset,
    int kOffset,
    const double *Ax,
    const double *Ay,
    const double *Az,
    const double *Lambda,
    double idx2,
    double idy2,
    double idz2);

extern void surfaceListFree(SurfaceListType *list);

#endif // __SURFACE_LIST_H_
