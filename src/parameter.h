/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __PARAMETER_H_
#define __PARAMETER_H_

typedef struct {
  int imax, jmax, kmax;
  double xlength, ylength, zlength;
  int itermax;
  double eps, omg;
  double re, tau, gamma;
  double te, dt;
  double gx, gy, gz;
  char *name;
  int bcLeft, bcRight, bcBottom, bcTop, bcFront, bcBack;
  double u_init, v_init, w_init, p_init;
  int levels, presmooth, postsmooth;
} Parameter;

void initParameter(Parameter *);
void readParameter(Parameter *, const char *);
void printParameter(Parameter *);
#endif
