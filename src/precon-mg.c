/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The multilevel cycle as a preconditioner: one cycle, from a zero initial
 * guess, every time.
 *
 * That is the whole implementation, and it is deliberately that small. The
 * requirement a preconditioner is held to is that it be a fixed linear
 * operator -- the same work on every application, no inner convergence test,
 * and a result that depends linearly on the residual it is given. One cycle
 * from zero satisfies all three by construction rather than by discipline:
 * there is no tolerance to pass in and nowhere for a warm start to hide.
 *
 * Not named solver-*.c, so every build links it; see multigrid.h.
 */
#include <stddef.h>

#include "multigrid.h"
#include "solver.h"

/*
 * z = M r, with M one V-cycle of the hierarchy in ctx, negated.
 *
 * The negation is the sign convention, not an adjustment. The cycle relaxes the
 * operator as assembled, which is negative semi-definite, so one cycle applied
 * to r approximates A^-1 r. The conjugate gradient solver works with the
 * negated system -- see pressureApplyOperator -- and needs an approximate
 * inverse of -A, which is what a Krylov method requires to be positive
 * definite. So M = -cycle, and without the sign the preconditioner is symmetric
 * and linear but negative definite, which CG cannot use.
 *
 * Solid cells come back at exactly zero without being masked here. The cycle
 * starts from a zero guess, the smoother keeps a solid cell at zero when its
 * right-hand side is zero there -- which it is, since r carries no residual in
 * the body -- and the coarse correction is confined to the fluid. That is the
 * invariant a Krylov solver needs at every iteration, not only at convergence.
 */
static void preconMgApply(
    void *ctx, const PressureLevelType *lv, const double *r, double *z)
{
  MultigridType *mg = (MultigridType *)ctx;

  size_t size = (size_t)(lv->imaxLocal + 2) * (lv->jmaxLocal + 2) * (lv->kmaxLocal + 2);

  for (size_t i = 0; i < size; i++) {
    z[i] = 0.0;
  }

  multigridCycle(mg, z, r);

  for (size_t i = 0; i < size; i++) {
    z[i] = -z[i];
  }
}

void preconMgInit(PreconType *precon, MultigridType *mg)
{
  precon->apply = preconMgApply;
  precon->ctx   = mg;
}
