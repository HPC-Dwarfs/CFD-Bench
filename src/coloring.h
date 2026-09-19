/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __COLORING_H_
#define __COLORING_H_

/*
 * The red-black checkerboard, defined on the global grid.
 *
 * Every sweep used to start its colour at the local index 1, so which colour a
 * cell belonged to depended on which rank held it and therefore on the number
 * of ranks. Deriving the colour from the cell's global position instead makes
 * the iteration -- and the answer -- the same under any decomposition.
 *
 * Local indices are 1-based over the interior and iOffset is the global index
 * of local cell 1 minus one, so i + iOffset identifies the cell globally up to
 * a constant shift. A constant shift only renames the two colours, and the
 * shift chosen here is the one that makes a serial run agree with the
 * compressed solver's existing (i + j + k) parity, so that layout's index
 * arithmetic keeps working unchanged.
 *
 * Everything that colours cells uses these helpers so that the bulk sweeps and
 * the obstacle surface pass cannot drift apart.
 */

/* The colour, 0 or 1, of the interior cell at local (i, j, k). */
static inline int colorOf(int i, int j, int k, int iOffset, int jOffset, int kOffset)
{
  return ((i + iOffset) + (j + jOffset) + (k + kOffset)) & 1;
}

/*
 * The first local i in a row of given colour, for the row at local (j, k).
 * Returns 1 or 2; the caller steps by 2 from there.
 */
static inline int colorRowStart(
    int color, int j, int k, int iOffset, int jOffset, int kOffset)
{
  /* colorOf(1, j, k, ...) is the colour of the first interior cell in the row.
   * If it already matches, start there, otherwise start one cell further in. */
  return (colorOf(1, j, k, iOffset, jOffset, kOffset) == color) ? 1 : 2;
}

#endif // __COLORING_H_
