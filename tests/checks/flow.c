/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * What a completed time step has to be true of once the domain contains a body:
 *
 *   - the aperture-weighted divergence vanishes in every fluid cell, and is no
 *     worse next to the body than anywhere else;
 *   - the velocity on every closed face is zero to machine precision;
 *   - the volume flux through any cross-section of the channel is the same,
 *     which follows from the divergence being zero and no flux crossing either
 *     the walls or the body;
 *   - the operator couples a cell only to its six face neighbours, at an
 *     obstacle edge where two faces are shut and at a corner where three are.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "discretization.h"
#include "parameter.h"
#include "pressure-bc.h"
#include "profiler.h"
#include "solver.h"
#include "util.h"

#define IMAX 48
#define JMAX 24
#define KMAX 24
#define STEPS 12

static Parameter Params;
static Discretization D;
static Solver S;

static void setup(CommType *base, const char *geometry, int inletBc)
{
  initParameter(&Params);

  Params.name         = "canal";
  Params.imax         = IMAX;
  Params.jmax         = JMAX;
  Params.kmax         = KMAX;
  Params.xlength      = 8.0;
  Params.ylength      = 4.0;
  Params.zlength      = 4.0;
  Params.eps          = 1e-9;
  Params.omg          = 1.8;
  Params.itermax      = 50000;
  Params.levels       = 3;
  Params.presmooth    = 5;
  Params.postsmooth   = 5;
  Params.re           = 100.0;
  Params.tau          = 0.5;
  Params.gamma        = 0.9;
  Params.dt           = 0.02;
  Params.te           = 0.0;
  Params.gx = Params.gy = Params.gz = 0.0;
  Params.u_init       = 1.0;
  Params.v_init = Params.w_init = Params.p_init = 0.0;
  Params.geometryFile = (char *)geometry;

  Params.bcLeft  = inletBc;
  Params.bcRight = OUTFLOW;
  Params.bcBottom = Params.bcTop = NOSLIP;
  Params.bcFront = Params.bcBack = NOSLIP;

  D.comm = *base;
  commPartition(&D.comm, KMAX, JMAX, IMAX);
  initDiscretization(&D, &Params);
  initSolver(&S, &D, &Params);
  initProfiler(&D.comm);
}

/* One time step, exactly as main.c runs it. */
static void step(void)
{
  computeTimestep(&D);
  setBoundaryConditions(&D);
  setSpecialBoundaryCondition(&D);
  computeFG(&D);
  computeRHS(&D);
  normalizePressure(&D);
  solve(&S, D.p, D.rhs);
  adaptUV(&D);
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("flow");

  /* A sphere in a channel, off-centre so the wake is not symmetric. */
  setup(&comm, "sphere:2.0,2.0,2.0,0.6", NOSLIP);

  for (int n = 0; n < STEPS; n++) {
    step();
  }

  /*
   * adaptUV updates the interior faces only; the velocity halo stays whatever
   * the last exchange left there, because nothing in the time loop reads it
   * again until computeFG exchanges it afresh. The divergence below does read
   * it, at the first interior cell of every rank, so it has to be brought up to
   * date here. On one rank that halo is the physical boundary and is already
   * right, which is exactly why this only shows up once the domain is divided.
   */
  commExchange(&D.comm, D.u);
  commExchange(&D.comm, D.v);
  commExchange(&D.comm, D.w);

  int imaxLocal        = D.comm.imaxLocal;
  int jmaxLocal        = D.comm.jmaxLocal;
  int kmaxLocal        = D.comm.kmaxLocal;

  const double *Ax     = D.Ax;
  const double *Ay     = D.Ay;
  const double *Az     = D.Az;
  const double *Lambda = D.Lambda;
  double *u            = D.u;
  double *v            = D.v;
  double *w            = D.w;

  double idx           = 1.0 / D.grid.dx;
  double idy           = 1.0 / D.grid.dy;
  double idz           = 1.0 / D.grid.dz;

  /*
   * Divergence, split into the cells that touch the body and the rest. A body
   * that the pressure solve could not see would show up here as a divergence
   * next to it far larger than anywhere else.
   */
  {
    double nearBody = 0.0, awayFromBody = 0.0, nearCount = 0.0;

    for (int k = 1; k < kmaxLocal + 1; k++) {
      for (int j = 1; j < jmaxLocal + 1; j++) {
        for (int i = 1; i < imaxLocal + 1; i++) {
          if (LAM(i, j, k) == 0.0) {
            continue;
          }

          double div =
              (AX(i, j, k) * U(i, j, k) - AX(i - 1, j, k) * U(i - 1, j, k)) * idx +
              (AY(i, j, k) * V(i, j, k) - AY(i, j - 1, k) * V(i, j - 1, k)) * idy +
              (AZ(i, j, k) * W(i, j, k) - AZ(i, j, k - 1) * W(i, j, k - 1)) * idz;

          div        = fabs(div);

          int cut    = (AX(i, j, k) == 0.0 || AX(i - 1, j, k) == 0.0 ||
                     AY(i, j, k) == 0.0 || AY(i, j - 1, k) == 0.0 ||
                     AZ(i, j, k) == 0.0 || AZ(i, j, k - 1) == 0.0);

          if (cut) {
            nearCount += 1.0;
            if (div > nearBody) {
              nearBody = div;
            }
          } else if (div > awayFromBody) {
            awayFromBody = div;
          }
        }
      }
    }

    commReduceAll(&nearBody, MAX);
    commReduceAll(&awayFromBody, MAX);
    commReduceAll(&nearCount, SUM);

    CHECK_TRUE(nearCount > 0.0, "no fluid cell touches the body");

    /* The pressure solve is converged to eps, and the divergence it leaves is
     * of that order divided by the time step. Both halves of the domain are
     * held to the same bound, which is the point: the body must not be worse. */
    double bound = 1e-4;

    CHECK_TRUE(awayFromBody < bound,
        "divergence away from the body is %.3e",
        awayFromBody);
    CHECK_TRUE(nearBody < bound,
        "divergence next to the body is %.3e against %.3e away from it",
        nearBody,
        awayFromBody);

    if (commIsMaster(&D.comm)) {
      printf("divergence: %.3e next to the body, %.3e away from it, over %.0f cut "
             "cells\n",
          nearBody,
          awayFromBody,
          nearCount);
    }
  }

  /* No velocity on a closed face, exactly. */
  {
    double worst = 0.0, faces = 0.0;

    for (int k = 0; k < kmaxLocal + 2; k++) {
      for (int j = 0; j < jmaxLocal + 2; j++) {
        for (int i = 0; i < imaxLocal + 2; i++) {
          if (AX(i, j, k) == 0.0) {
            faces += 1.0;
            worst = fmax(worst, fabs(U(i, j, k)));
          }
          if (AY(i, j, k) == 0.0) {
            faces += 1.0;
            worst = fmax(worst, fabs(V(i, j, k)));
          }
          if (AZ(i, j, k) == 0.0) {
            faces += 1.0;
            worst = fmax(worst, fabs(W(i, j, k)));
          }
        }
      }
    }

    commReduceAll(&worst, MAX);
    commReduceAll(&faces, SUM);

    CHECK_TRUE(faces > 0.0, "the sphere closed no faces");
    CHECK_NEAR(worst, 0.0, 0.0,
        "a velocity on a closed face is nonzero over %.0f closed faces",
        faces);
  }

  /*
   * Volume flux through cross-sections. With a divergence-free field, no flux
   * through the walls and none through the body, the flux through every plane
   * of constant x has to be the same -- at any time, not only at a steady
   * state, because it follows from the divergence rather than from the flow
   * having settled.
   */
  {
    int offsets[NDIMS] = { 0, 0, 0 };
    commGetOffsets(&D.comm, offsets, KMAX, JMAX, IMAX);

    double area    = D.grid.dy * D.grid.dz;
    double lo = 0.0, hi = 0.0;
    int haveAny    = 0;

    for (int gi = IMAX / 4; gi < 3 * IMAX / 4; gi += IMAX / 8) {
      double flux = 0.0;
      int i       = gi - offsets[IDIM] + 1;

      if (i >= 1 && i < imaxLocal + 1) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
          for (int j = 1; j < jmaxLocal + 1; j++) {
            flux += AX(i, j, k) * U(i, j, k) * area;
          }
        }
      }

      commReduceAll(&flux, SUM);

      if (!haveAny) {
        lo = hi = flux;
        haveAny = 1;
      } else {
        lo = fmin(lo, flux);
        hi = fmax(hi, flux);
      }
    }

    CHECK_TRUE(haveAny, "no cross-section was sampled");
    CHECK_TRUE(fabs(hi - lo) < 1e-6 * (fabs(hi) + 1e-30),
        "volume flux varies between cross-sections, from %.17g to %.17g",
        lo,
        hi);

    if (commIsMaster(&D.comm)) {
      printf("cross-section flux: %.6f to %.6f, spread %.3e\n", lo, hi, hi - lo);
    }
  }

  /*
   * The stencil stays compact at an edge and at a corner of the body, where two
   * and three of a cell's faces are shut. A body with a staircase surface has
   * plenty of both.
   */
  {
    size_t size   = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);
    double *x     = calloc(size, sizeof(double));
    double *ax    = calloc(size, sizeof(double));
    double *p     = x;

    int strideJ   = imaxLocal + 2;
    int strideK   = (imaxLocal + 2) * (jmaxLocal + 2);

    double idx2   = 1.0 / (D.grid.dx * D.grid.dx);
    double idy2   = 1.0 / (D.grid.dy * D.grid.dy);
    double idz2   = 1.0 / (D.grid.dz * D.grid.dz);

    for (int want = 2; want <= 3; want++) {
      for (size_t i = 0; i < size; i++) {
        x[i] = 0.0;
      }

      size_t centre = 0;
      int have      = 0;

      for (int k = 2; k < kmaxLocal && !have; k++) {
        for (int j = 2; j < jmaxLocal && !have; j++) {
          for (int i = 2; i < imaxLocal && !have; i++) {
            if (LAM(i, j, k) == 0.0) {
              continue;
            }

            int shut = (AX(i, j, k) == 0.0) + (AX(i - 1, j, k) == 0.0) +
                       (AY(i, j, k) == 0.0) + (AY(i, j - 1, k) == 0.0) +
                       (AZ(i, j, k) == 0.0) + (AZ(i, j, k - 1) == 0.0);

            if (shut == want) {
              centre = (size_t)k * strideK + (size_t)j * strideJ + (size_t)i;
              have   = 1;
            }
          }
        }
      }

      if (have) {
        x[centre] = 1.0;
      }

      commExchange(&D.comm, x);
      pressureBcApply(&S.bc, &D.comm, x, imaxLocal, jmaxLocal, kmaxLocal);

      for (size_t i = 0; i < size; i++) {
        ax[i] = 0.0;
      }

      for (int k = 1; k < kmaxLocal + 1; k++) {
        for (int j = 1; j < jmaxLocal + 1; j++) {
          for (int i = 1; i < imaxLocal + 1; i++) {
            if (LAM(i, j, k) == 0.0) {
              continue;
            }
            double pc = P(i, j, k);
            ax[(size_t)k * strideK + (size_t)j * strideJ + (size_t)i] =
                AX(i, j, k) * (P(i + 1, j, k) - pc) * idx2 +
                AX(i - 1, j, k) * (P(i - 1, j, k) - pc) * idx2 +
                AY(i, j, k) * (P(i, j + 1, k) - pc) * idy2 +
                AY(i, j - 1, k) * (P(i, j - 1, k) - pc) * idy2 +
                AZ(i, j, k) * (P(i, j, k + 1) - pc) * idz2 +
                AZ(i, j, k - 1) * (P(i, j, k - 1) - pc) * idz2;
          }
        }
      }

      double strays = 0.0;
      double found  = have ? 1.0 : 0.0;

      if (have) {
        for (int k = 1; k < kmaxLocal + 1; k++) {
          for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imaxLocal + 1; i++) {
              size_t idx = (size_t)k * strideK + (size_t)j * strideJ + (size_t)i;

              if (idx == centre || idx == centre + 1 || idx == centre - 1 ||
                  idx == centre + strideJ || idx == centre - strideJ ||
                  idx == centre + strideK || idx == centre - strideK) {
                continue;
              }

              if (ax[idx] != 0.0) {
                strays += 1.0;
              }
            }
          }
        }
      }

      commReduceAll(&strays, SUM);
      commReduceAll(&found, SUM);

      CHECK_TRUE(found > 0.0,
          "no cell has exactly %d closed faces, so that case is untested",
          want);
      CHECK_TRUE(strays == 0.0,
          "a cell with %d closed faces couples to %.0f cells beyond its six face "
          "neighbours",
          want,
          strays);
    }

    free(x);
    free(ax);
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
