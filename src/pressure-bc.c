/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <stdio.h>
#include <stdlib.h>

#include "discretization.h"
#include "pressure-bc.h"
#include "util.h"

static const char *const DIRECTION_NAME[NDIRS] = {
  "bcLeft", "bcRight", "bcBottom", "bcTop", "bcFront", "bcBack"
};

static const char *typeName(int type)
{
  switch (type) {
  case NOSLIP:
    return "NOSLIP";
  case SLIP:
    return "SLIP";
  case OUTFLOW:
    return "OUTFLOW";
  case PERIODIC:
    return "PERIODIC";
  default:
    return "unknown";
  }
}

void pressureBcInit(PressureBcType *bc,
    int bcLeft,
    int bcRight,
    int bcBottom,
    int bcTop,
    int bcFront,
    int bcBack)
{
  bc->type[LEFT]   = bcLeft;
  bc->type[RIGHT]  = bcRight;
  bc->type[BOTTOM] = bcBottom;
  bc->type[TOP]    = bcTop;
  bc->type[FRONT]  = bcFront;
  bc->type[BACK]   = bcBack;

  for (int d = 0; d < NDIRS; d++) {
    switch (bc->type[d]) {
    case NOSLIP:
    case SLIP:
    case OUTFLOW:
      break;
    case PERIODIC:
      fprintf(stderr,
          "Boundary %s is %s, which the pressure solver does not implement. The "
          "velocity boundary condition does not implement it either, so the run "
          "would silently solve a different problem.\n",
          DIRECTION_NAME[d],
          typeName(bc->type[d]));
      exit(EXIT_FAILURE);
    default:
      fprintf(stderr,
          "Boundary %s has unknown boundary condition %d. Expected 1 (NOSLIP), 2 "
          "(SLIP) or 3 (OUTFLOW).\n",
          DIRECTION_NAME[d],
          bc->type[d]);
      exit(EXIT_FAILURE);
    }
  }
}

bool pressureBcIsSingular(const PressureBcType *bc)
{
  for (int d = 0; d < NDIRS; d++) {
    if (bc->type[d] == OUTFLOW) {
      return false;
    }
  }

  return true;
}

/* A wall mirrors the interior into the halo, giving a zero normal gradient. An
 * outflow reflects it oddly, putting zero on the boundary face. */
static double haloValue(int type, double interior)
{
  return (type == OUTFLOW) ? -interior : interior;
}

void pressureBcApply(const PressureBcType *bc,
    CommType *comm,
    double *p,
    int imaxLocal,
    int jmaxLocal,
    int kmaxLocal)
{
  if (commIsBoundary(comm, FRONT)) {
    int type = bc->type[FRONT];
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        P(i, j, 0) = haloValue(type, P(i, j, 1));
      }
    }
  }

  if (commIsBoundary(comm, BACK)) {
    int type = bc->type[BACK];
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        P(i, j, kmaxLocal + 1) = haloValue(type, P(i, j, kmaxLocal));
      }
    }
  }

  if (commIsBoundary(comm, BOTTOM)) {
    int type = bc->type[BOTTOM];
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        P(i, 0, k) = haloValue(type, P(i, 1, k));
      }
    }
  }

  if (commIsBoundary(comm, TOP)) {
    int type = bc->type[TOP];
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        P(i, jmaxLocal + 1, k) = haloValue(type, P(i, jmaxLocal, k));
      }
    }
  }

  if (commIsBoundary(comm, LEFT)) {
    int type = bc->type[LEFT];
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        P(0, j, k) = haloValue(type, P(1, j, k));
      }
    }
  }

  if (commIsBoundary(comm, RIGHT)) {
    int type = bc->type[RIGHT];
    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        P(imaxLocal + 1, j, k) = haloValue(type, P(imaxLocal, j, k));
      }
    }
  }
}
