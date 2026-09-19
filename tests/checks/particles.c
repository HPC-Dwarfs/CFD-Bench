/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The particle tracer:
 *
 *   - seed positions depend on the setup and nothing else, so two runs and two
 *     rank counts inject the same particles;
 *   - a particle is injected only inside the seed region and only into fluid;
 *   - it follows a uniform flow exactly, and follows the time step the solver
 *     actually took rather than the one the parameter file asked for;
 *   - it cannot cross a closed face, including a body one face thick, which is
 *     the case a destination-cell test misses entirely;
 *   - it is conserved: everything injected is either still in the domain, gone
 *     through a boundary, or stopped at the body;
 *   - storage tracks the number of particles rather than the size of the grid.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "discretization.h"
#include "parameter.h"
#include "particletracing.h"
#include "profiler.h"
#include "solver.h"
#include "util.h"

#define IMAX 32
#define JMAX 16
#define KMAX 16

static Parameter Params;
static Discretization D;
static Solver S;
static ParticleTracerType Tracer;

static void setup(CommType *base, const char *geometry, int particles, int grid)
{
  initParameter(&Params);

  Params.name              = "check-particles";
  Params.imax              = grid * 2;
  Params.jmax              = grid;
  Params.kmax              = grid;
  Params.xlength           = 2.0;
  Params.ylength           = 1.0;
  Params.zlength           = 1.0;
  Params.eps               = 1e-7;
  Params.omg               = 1.7;
  Params.itermax           = 5000;
  Params.levels            = 1;
  Params.presmooth         = 4;
  Params.postsmooth        = 4;
  Params.re                = 100.0;
  Params.tau               = 0.5;
  Params.gamma             = 0.9;
  Params.dt                = 0.01;
  Params.te                = 0.0;
  Params.gx = Params.gy = Params.gz = 0.0;
  Params.u_init = Params.v_init = Params.w_init = Params.p_init = 0.0;
  Params.geometryFile      = (char *)geometry;

  Params.bcLeft = Params.bcRight = NOSLIP;
  Params.bcBottom = Params.bcTop = NOSLIP;
  Params.bcFront = Params.bcBack = NOSLIP;

  Params.numberOfParticles = particles;
  Params.startTime         = 0.0;
  Params.injectTimePeriod  = 0.1;
  Params.writeTimePeriod   = 0.0; /* no files from a check */
  Params.x1                = 0.1;
  Params.y1                = 0.1;
  Params.z1                = 0.1;
  Params.x2                = 0.5;
  Params.y2                = 0.9;
  Params.z2                = 0.9;

  D.comm = *base;
  commPartition(&D.comm, Params.kmax, Params.jmax, Params.imax);
  initDiscretization(&D, &Params);
  initSolver(&S, &D, &Params);
  initProfiler(&D.comm);
  particleTracerInit(&Tracer, &D, &Params);
}

/* Set the whole velocity field to a constant. */
static void uniformFlow(double u, double v, double w)
{
  int imaxLocal = D.comm.imaxLocal;
  int jmaxLocal = D.comm.jmaxLocal;
  int kmaxLocal = D.comm.kmaxLocal;
  size_t size   = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

  for (size_t i = 0; i < size; i++) {
    D.u[i] = u;
    D.v[i] = v;
    D.w[i] = w;
  }
}

static double globalCount(int local)
{
  double v = (double)local;
  commReduceAll(&v, SUM);
  return v;
}

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);
  commPartition(&comm, 4, 4, 4);

  CHECK_BEGIN("particles");

  /* Seed positions are a function of the configuration alone. */
  {
    setup(&comm, NULL, 64, KMAX);

    int inRegion = 0;
    for (int i = 0; i < 64; i++) {
      double x, y, z;
      particleTracerSeed(&Tracer, 0, i, &x, &y, &z);

      if (x >= Params.x1 && x <= Params.x2 && y >= Params.y1 && y <= Params.y2 &&
          z >= Params.z1 && z <= Params.z2) {
        ++inRegion;
      }

      /* Asking twice has to give the same answer, which a generator with
       * hidden state would not. */
      double x2, y2, z2;
      particleTracerSeed(&Tracer, 0, i, &x2, &y2, &z2);
      CHECK_TRUE(x == x2 && y == y2 && z == z2,
          "seed %d is not reproducible within a run",
          i);
    }

    CHECK_TRUE(inRegion == 64, "%d of 64 seeds fell outside the seed region", inRegion);

    /* Different batches must not repeat the same positions. */
    double ax, ay, az, bx, by, bz;
    particleTracerSeed(&Tracer, 0, 0, &ax, &ay, &az);
    particleTracerSeed(&Tracer, 1, 0, &bx, &by, &bz);
    CHECK_TRUE(ax != bx || ay != by || az != bz,
        "two batches produced the same seed position");
  }

  /* A particle follows a uniform flow exactly. */
  {
    setup(&comm, NULL, 200, KMAX);
    uniformFlow(1.0, 0.0, 0.0);

    particleTracerInject(&Tracer, &D);

    double before = globalCount(particleTracerLiveCount(&Tracer));
    CHECK_TRUE(before > 0.0, "no particle was injected into an empty channel");

    /* Advance, then check every survivor against the seed it came from. The
     * seeds are reproducible, so a particle can be matched to its own by the
     * two coordinates a flow along x does not change. */
    double dt = 0.05;
    int steps = 4;

    for (int n = 0; n < steps; n++) {
      particleTracerAdvance(&Tracer, &D, dt);
      particleTracerMigrate(&Tracer, &D);
    }

    double worstX = 0.0;
    int unmatched = 0;

    for (int i = 0; i < particleTracerLiveCount(&Tracer); i++) {
      double px = Tracer.pool[i].x;
      double py = Tracer.pool[i].y;
      double pz = Tracer.pool[i].z;

      int found = 0;

      for (int sIdx = 0; sIdx < 200; sIdx++) {
        double sx, sy, sz;
        particleTracerSeed(&Tracer, 0, sIdx, &sx, &sy, &sz);

        if (sy == py && sz == pz) {
          double want = sx + steps * dt * 1.0;
          double err  = fabs(px - want);
          if (err > worstX) {
            worstX = err;
          }
          found = 1;
          break;
        }
      }

      if (!found) {
        ++unmatched;
      }
    }

    commReduceAll(&worstX, MAX);

    CHECK_TRUE(globalCount(unmatched) == 0.0,
        "%.0f particles could not be matched to the seed they came from",
        globalCount(unmatched));
    CHECK_NEAR(worstX, 0.0, 1e-12,
        "a particle in a uniform flow of 1.0 is off its exact position by %.3e "
        "after %d steps of %g",
        worstX,
        steps,
        dt);

    double after = globalCount(particleTracerLiveCount(&Tracer));
    CHECK_TRUE(after == before,
        "a uniform flow with no obstacle lost %.0f of %.0f particles",
        before - after,
        before);
  }

  /* The displacement follows the step it is given, not the one configured. */
  {
    setup(&comm, NULL, 1, KMAX);
    uniformFlow(1.0, 0.0, 0.0);

    double sx, sy, sz;
    particleTracerSeed(&Tracer, 0, 0, &sx, &sy, &sz);

    particleTracerInject(&Tracer, &D);

    int holder  = particleTracerLiveCount(&Tracer) > 0;
    double moved = 0.0;
    double dt    = 0.037; /* deliberately not Params.dt */

    particleTracerAdvance(&Tracer, &D, dt);
    particleTracerMigrate(&Tracer, &D);

    if (holder && particleTracerLiveCount(&Tracer) > 0) {
      moved = Tracer.pool[0].x - sx;
    }

    commReduceAll(&moved, SUM);

    CHECK_NEAR(moved, dt * 1.0, 1e-12,
        "the particle moved by %.12f where the step it was given implies %.12f",
        moved,
        dt * 1.0);
  }

  /* A body one face thick has to stop every particle, and none may end up
   * inside a solid cell. */
  {
    /* A plate at x = 1.0 covering the middle of the cross-section. It cannot
     * span the whole of it: that would cut the fluid into two regions, which
     * the connectivity check refuses before a solver ever sees it. The rim left
     * open is also what makes the check sharp -- particles behind the plate
     * must stop, and the ones beside it must not. */
    setup(&comm, "plate-x:1.0,0.2,0.2,0.8,0.8", 400, KMAX);
    uniformFlow(1.0, 0.0, 0.0);

    particleTracerInject(&Tracer, &D);
    double injected = globalCount(particleTracerLiveCount(&Tracer));
    CHECK_TRUE(injected > 0.0, "no particle was injected upstream of the plate");

    /* Long enough that an unobstructed particle would be well past the plate. */
    for (int n = 0; n < 60; n++) {
      particleTracerAdvance(&Tracer, &D, 0.03);
      particleTracerMigrate(&Tracer, &D);
    }

    int beyond = 0, past = 0;
    for (int i = 0; i < particleTracerLiveCount(&Tracer); i++) {
      double py = Tracer.pool[i].y;
      double pz = Tracer.pool[i].z;
      int behind = (py >= 0.2 && py <= 0.8 && pz >= 0.2 && pz <= 0.8);

      if (Tracer.pool[i].x > 1.0) {
        if (behind) {
          ++beyond;
        } else {
          ++past;
        }
      }
    }

    CHECK_TRUE(globalCount(beyond) == 0.0,
        "%.0f particles crossed a plate one face thick",
        globalCount(beyond));
    CHECK_TRUE(globalCount(past) > 0.0,
        "no particle got past the plate around its edge, so the check cannot tell "
        "obstruction from the particles simply not arriving");

    double injectedT, alive, atBoundary, atBody;
    particleTracerTotals(&Tracer, &D, &injectedT, &alive, &atBoundary, &atBody);

    CHECK_TRUE(atBody > 0.0, "the plate stopped no particles at all");
    CHECK_NEAR(injectedT,
        alive + atBoundary + atBody,
        0.0,
        "particles are not conserved: %.0f injected against %.0f alive, %.0f gone "
        "through a boundary and %.0f stopped at the body",
        injectedT,
        alive,
        atBoundary,
        atBody);
  }

  /* No particle may sit inside the body. */
  {
    setup(&comm, "sphere:1.0,0.5,0.5,0.25", 400, KMAX);
    uniformFlow(1.0, 0.0, 0.0);

    for (int n = 0; n < 40; n++) {
      particleTracerInject(&Tracer, &D);
      particleTracerAdvance(&Tracer, &D, 0.02);
      particleTracerMigrate(&Tracer, &D);
    }

    int imaxLocal        = D.comm.imaxLocal;
    int jmaxLocal        = D.comm.jmaxLocal;
    const double *Lambda = D.Lambda;
    int inside           = 0;

    for (int i = 0; i < particleTracerLiveCount(&Tracer); i++) {
      int ci = (int)floor(Tracer.pool[i].x / D.grid.dx) - D.iOffset + 1;
      int cj = (int)floor(Tracer.pool[i].y / D.grid.dy) - D.jOffset + 1;
      int ck = (int)floor(Tracer.pool[i].z / D.grid.dz) - D.kOffset + 1;

      if (ci < 0 || ci > imaxLocal + 1 || cj < 0 || cj > jmaxLocal + 1 || ck < 0 ||
          ck > D.comm.kmaxLocal + 1) {
        continue;
      }

      if (LAM(ci, cj, ck) == 0.0) {
        ++inside;
      }
    }

    CHECK_TRUE(globalCount(inside) == 0.0,
        "%.0f particles ended up inside the body",
        globalCount(inside));

    double injectedT, alive, atBoundary, atBody;
    particleTracerTotals(&Tracer, &D, &injectedT, &alive, &atBoundary, &atBody);
    CHECK_NEAR(injectedT,
        alive + atBoundary + atBody,
        0.0,
        "particles are not conserved around a sphere");
  }

  /* Storage follows the particle count, not the grid. */
  {
    setup(&comm, NULL, 500, 16);
    particleTracerInject(&Tracer, &D);
    int coarse = Tracer.capacity;

    setup(&comm, NULL, 500, 32);
    particleTracerInject(&Tracer, &D);
    int fine   = Tracer.capacity;

    CHECK_TRUE(coarse == fine,
        "the pool grew from %d to %d slots when the grid was refined, so it is "
        "sized by the grid rather than by the particles",
        coarse,
        fine);
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
