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
  /* Obstacle geometry: a path to a voxel volume, or an analytic body such as
   * "sphere:xc,yc,zc,r". Absent means an obstacle-free domain. */
  char *geometryFile;
  /* Particle tracing. Absent, or a count of zero, means no tracing at all. The
   * seed region is the axis-aligned box (x1,y1,z1) to (x2,y2,z2). */
  int numberOfParticles;
  double startTime, injectTimePeriod, writeTimePeriod;
  double x1, y1, z1, x2, y2, z2;
} Parameter;

void initParameter(Parameter *);
void readParameter(Parameter *, const char *);
void printParameter(Parameter *);
#endif
