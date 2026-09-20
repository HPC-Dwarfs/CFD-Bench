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
  /* Relaxation factor the multigrid smoother uses. Not omg: that is the
   * optimum for SOR as a solver, and a factor above 1 amplifies exactly the
   * high-frequency modes a smoother exists to damp. The two were only ever
   * confused because a one-directional cycle hid the difference.
   *
   * 1.3 is the fastest value measured that both converges and leaves the
   * converged field inside the cross-solver agreement gate: 1.6 is stable and
   * quicker but stops marginally outside it, and 1.8 -- which the setups carry
   * as omg -- diverges. Lower is safer and slower. */
  double smoothOmega;
  /* Which preconditioner the conjugate gradient solver uses: "none" or
   * "jacobi". Read by every build and acted on only by SOLVER=cg, the way
   * levels and the smoothing counts are read by every build and acted on only
   * by SOLVER=mg. An unsupported value is refused at initialization. */
  char *precon;
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
