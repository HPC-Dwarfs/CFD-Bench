/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#ifndef __FORCES_H_
#define __FORCES_H_

#include "discretization.h"

/*
 * The force the flow exerts on the embedded body.
 *
 * Integrated over the closed faces, which are the body's surface: a face with
 * zero aperture separating a fluid cell from a solid one is a piece of wall, and
 * its area and orientation are known exactly. Two contributions:
 *
 *   pressure   the fluid pressure next to the face, acting along the face
 *              normal. For a bluff body at these Reynolds numbers this is most
 *              of the drag.
 *   viscous    the wall shear, taken from the tangential velocity one half cell
 *              from the wall. On a staircase surface this is a first-order
 *              estimate, which is the accuracy binary apertures allow.
 *
 * The staircase is why the result is worth reporting rather than trusting: a
 * body whose surface is made of axis-aligned faces has a larger area than the
 * smooth body it approximates, and no amount of care in this integration fixes
 * that. Fractional apertures are what would.
 */
extern void forcesCompute(Discretization *d, double *fx, double *fy, double *fz);

/* True when the domain actually contains a body, so a caller can skip the
 * whole business for an obstacle-free setup. */
extern int forcesHaveBody(Discretization *d);

#endif // __FORCES_H_
