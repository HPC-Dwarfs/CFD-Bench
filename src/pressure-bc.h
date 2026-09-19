/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __PRESSURE_BC_H_
#define __PRESSURE_BC_H_

#include "comm.h"
#include <stdbool.h>

/*
 * The pressure boundary condition, in one place.
 *
 * Every solver used to impose a zero normal gradient on all six physical
 * boundaries whatever the setup said, so the velocity boundary condition never
 * reached the pressure operator and the operator was singular in every setup --
 * including the ones with an outflow, which the null-space handling would then
 * decline to treat.
 *
 * The condition now follows the discretization:
 *
 *   NOSLIP, SLIP   a wall. No flux crosses it, so the normal pressure gradient
 *                  is zero and the halo mirrors the interior.
 *   OUTFLOW        the pressure is pinned to zero on the boundary face. For a
 *                  cell-centred field the face value is the mean of the halo
 *                  and the interior cell, so this is an odd reflection.
 *   PERIODIC       rejected. setBoundaryConditions implements it as an empty
 *                  case for velocity and the Cartesian communicator is built
 *                  non-periodic, so treating it as a wall here would be a
 *                  silent disagreement with a velocity condition that does not
 *                  exist.
 *
 * The homogeneous form of either condition is the same, which is why the same
 * routine serves a multigrid level carrying an error as well as the finest
 * level carrying the pressure.
 */

/* The six boundary condition codes of a setup, indexed by DirectionType. */
typedef struct {
  int type[NDIRS];
} PressureBcType;

/* Collect the setup's six boundary codes. Aborts naming the boundary and the
 * code when one is unsupported, before any time step runs. */
extern void pressureBcInit(PressureBcType *bc,
    int bcLeft,
    int bcRight,
    int bcBottom,
    int bcTop,
    int bcFront,
    int bcBack);

/* Apply the condition to the halo of p on this rank's physical boundaries.
 * Ranks that own no physical boundary in a direction leave that halo to the
 * exchange. Valid for any level, given that level's local extents. */
extern void pressureBcApply(const PressureBcType *bc,
    CommType *comm,
    double *p,
    int imaxLocal,
    int jmaxLocal,
    int kmaxLocal);

/* True when every boundary imposes a zero normal gradient, which leaves the
 * operator singular with the constants in its null space. Read from the
 * configuration rather than from the matrix. */
extern bool pressureBcIsSingular(const PressureBcType *bc);

#endif // __PRESSURE_BC_H_
