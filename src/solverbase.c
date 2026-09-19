/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Not named solver-base.c on purpose: the Makefile treats solver-*.c as the
 * mutually exclusive solver variants and links exactly one of them.
 */
#include "solver.h"
#include "util.h"

/*
 * The residual norm of a field, as a dedicated pass.
 *
 * The relaxation sweeps accumulate a residual as they go, which costs nothing
 * extra but describes the field as it was part-way through the sweep, not the
 * field that comes back. That is good enough to decide when to stop iterating
 * and wrong to report. So the sweeps keep their cheap accumulation for loop
 * control, and this runs once per solve to produce the number that is reported
 * and returned -- one extra pass per solve rather than per iteration, which
 * keeps it off the benchmark's headline kernel.
 *
 * Returns the mean square residual, to be compared against eps * eps, which is
 * the convention the solvers already used.
 */
double pressureResidualNorm(CommType *comm,
    const PressureBcType *bc,
    double *p,
    const double *rhs,
    int imaxLocal,
    int jmaxLocal,
    int kmaxLocal,
    double dx,
    double dy,
    double dz,
    double globalCells)
{
  double idx2 = 1.0 / (dx * dx);
  double idy2 = 1.0 / (dy * dy);
  double idz2 = 1.0 / (dz * dz);
  double res  = 0.0;

  commExchange(comm, p);
  pressureBcApply(bc, comm, p, imaxLocal, jmaxLocal, kmaxLocal);

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        double r = RHS(i, j, k) -
                   ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                       (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) * idy2 +
                       (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) * idz2);
        res += r * r;
      }
    }
  }

  commReduceAll(&res, SUM);
  return res / globalCells;
}

void solverBaseInit(Solver *s, Discretization *d, Parameter *p)
{
  s->eps      = p->eps;
  s->omega    = p->omg;
  s->itermax  = p->itermax;
  s->grid     = &d->grid;
  s->comm     = &d->comm;
  s->problem  = p->name;

  s->Ax       = d->Ax;
  s->Ay       = d->Ay;
  s->Az       = d->Az;
  s->Lambda   = d->Lambda;

  s->bcLeft   = p->bcLeft;
  s->bcRight  = p->bcRight;
  s->bcBottom = p->bcBottom;
  s->bcTop    = p->bcTop;
  s->bcFront  = p->bcFront;
  s->bcBack   = p->bcBack;

  pressureBcInit(&s->bc,
      p->bcLeft,
      p->bcRight,
      p->bcBottom,
      p->bcTop,
      p->bcFront,
      p->bcBack);

  int offsets[NDIMS] = { 0, 0, 0 };
  commGetOffsets(s->comm, offsets, p->kmax, p->jmax, p->imax);
  s->iOffset = offsets[IDIM];
  s->jOffset = offsets[JDIM];
  s->kOffset = offsets[KDIM];
}
