/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include "forces.h"
#include "util.h"

int forcesHaveBody(Discretization *d)
{
  int imaxLocal        = d->comm.imaxLocal;
  int jmaxLocal        = d->comm.jmaxLocal;
  int kmaxLocal        = d->comm.kmaxLocal;
  const double *Lambda = d->Lambda;
  double solid         = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (LAM(i, j, k) == 0.0) {
          solid += 1.0;
        }
      }
    }
  }

  commReduceAll(&solid, SUM);
  return solid > 0.0;
}

void forcesCompute(Discretization *d, double *fx, double *fy, double *fz)
{
  int imaxLocal        = d->comm.imaxLocal;
  int jmaxLocal        = d->comm.jmaxLocal;
  int kmaxLocal        = d->comm.kmaxLocal;

  const double *Ax     = d->Ax;
  const double *Ay     = d->Ay;
  const double *Az     = d->Az;
  const double *Lambda = d->Lambda;

  double *p            = d->p;
  double *u            = d->u;
  double *v            = d->v;
  double *w            = d->w;

  double dx            = d->grid.dx;
  double dy            = d->grid.dy;
  double dz            = d->grid.dz;

  double areaX         = dy * dz;
  double areaY         = dx * dz;
  double areaZ         = dx * dy;

  double nu            = 1.0 / d->re;

  /* The shear terms read a velocity one cell outside this rank's interior, and
   * adaptUV leaves the velocity halo as the last exchange left it. On one rank
   * that halo is the physical boundary and is already right; once the domain is
   * divided it is a neighbour's stale value, so it has to be refreshed before
   * the surface of a body that straddles a rank boundary can be integrated. */
  commExchange(&d->comm, u);
  commExchange(&d->comm, v);
  commExchange(&d->comm, w);

  double sx = 0.0, sy = 0.0, sz = 0.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {

        /* ---- faces normal to x ---- */
        if (AX(i, j, k) == 0.0) {
          /* A wall only where fluid meets solid; a closed face between two
           * solid cells is inside the body and carries no load. */
          if (LAM(i, j, k) > 0.0 && LAM(i + 1, j, k) == 0.0) {
            /* Body lies at higher x, so its outward normal is -x and the
             * pressure pushes it towards +x. */
            sx += P(i, j, k) * areaX;

            /* Shear from the tangential velocity half a cell away.
             *
             * Unlike the pressure term, this does not change sign with the side
             * the body is on. Flow along a wall drags the wall the way it is
             * going whether the fluid is above or below it, so both branches
             * add: the normal reverses and so does the velocity gradient, and
             * the two reversals cancel. */
            double vt = 0.5 * (V(i, j, k) + V(i, j - 1, k));
            double wt = 0.5 * (W(i, j, k) + W(i, j, k - 1));
            sy += nu * vt / (0.5 * dx) * areaX;
            sz += nu * wt / (0.5 * dx) * areaX;
          } else if (LAM(i, j, k) == 0.0 && LAM(i + 1, j, k) > 0.0) {
            sx -= P(i + 1, j, k) * areaX;

            double vt = 0.5 * (V(i + 1, j, k) + V(i + 1, j - 1, k));
            double wt = 0.5 * (W(i + 1, j, k) + W(i + 1, j, k - 1));
            sy += nu * vt / (0.5 * dx) * areaX;
            sz += nu * wt / (0.5 * dx) * areaX;
          }
        }

        /* ---- faces normal to y ---- */
        if (AY(i, j, k) == 0.0) {
          if (LAM(i, j, k) > 0.0 && LAM(i, j + 1, k) == 0.0) {
            sy += P(i, j, k) * areaY;

            double ut = 0.5 * (U(i, j, k) + U(i - 1, j, k));
            double wt = 0.5 * (W(i, j, k) + W(i, j, k - 1));
            sx += nu * ut / (0.5 * dy) * areaY;
            sz += nu * wt / (0.5 * dy) * areaY;
          } else if (LAM(i, j, k) == 0.0 && LAM(i, j + 1, k) > 0.0) {
            sy -= P(i, j + 1, k) * areaY;

            double ut = 0.5 * (U(i, j + 1, k) + U(i - 1, j + 1, k));
            double wt = 0.5 * (W(i, j + 1, k) + W(i, j + 1, k - 1));
            sx += nu * ut / (0.5 * dy) * areaY;
            sz += nu * wt / (0.5 * dy) * areaY;
          }
        }

        /* ---- faces normal to z ---- */
        if (AZ(i, j, k) == 0.0) {
          if (LAM(i, j, k) > 0.0 && LAM(i, j, k + 1) == 0.0) {
            sz += P(i, j, k) * areaZ;

            double ut = 0.5 * (U(i, j, k) + U(i - 1, j, k));
            double vt = 0.5 * (V(i, j, k) + V(i, j - 1, k));
            sx += nu * ut / (0.5 * dz) * areaZ;
            sy += nu * vt / (0.5 * dz) * areaZ;
          } else if (LAM(i, j, k) == 0.0 && LAM(i, j, k + 1) > 0.0) {
            sz -= P(i, j, k + 1) * areaZ;

            double ut = 0.5 * (U(i, j, k + 1) + U(i - 1, j, k + 1));
            double vt = 0.5 * (V(i, j, k + 1) + V(i, j - 1, k + 1));
            sx += nu * ut / (0.5 * dz) * areaZ;
            sy += nu * vt / (0.5 * dz) * areaZ;
          }
        }
      }
    }
  }

  commReduceAll(&sx, SUM);
  commReduceAll(&sy, SUM);
  commReduceAll(&sz, SUM);

  *fx = sx;
  *fy = sy;
  *fz = sz;
}
