/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file. */
#ifndef __UTIL_H_
#define __UTIL_H_
#define HLINE "--------------------------------------------------------------\n"

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef ABS
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#endif

#define P(i, j, k)                                                                       \
  p[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define F(i, j, k)                                                                       \
  f[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define G(i, j, k)                                                                       \
  g[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define H(i, j, k)                                                                       \
  h[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define U(i, j, k)                                                                       \
  u[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define V(i, j, k)                                                                       \
  v[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define W(i, j, k)                                                                       \
  w[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]
#define RHS(i, j, k)                                                                     \
  rhs[(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (j) * (imaxLocal + 2) + (i)]

#endif // __UTIL_H_
