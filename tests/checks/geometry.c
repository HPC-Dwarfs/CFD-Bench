/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Exercises the geometry producers directly, without a solver: the producer
 * interface deliberately takes a plain description of the grid rather than the
 * solver structures, so a driver can fill aperture arrays for a grid nobody
 * owns and look at the result.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "check.h"
#include "geometry.h"
#include "geometry-voxel.h"

#define IMAX 24
#define JMAX 24
#define KMAX 24

#define AT(f, i, j, k, im, jm)                                                           \
  f[(size_t)(k) * ((im) + 2) * ((jm) + 2) + (size_t)(j) * ((im) + 2) + (size_t)(i)]

static GeometryDomainType makeDomain(int im, int jm, int km, double lx, double ly,
    double lz)
{
  GeometryDomainType d;
  memset(&d, 0, sizeof(d));

  d.imaxLocal = im;
  d.jmaxLocal = jm;
  d.kmaxLocal = km;
  d.imax      = im;
  d.jmax      = jm;
  d.kmax      = km;
  d.xlength   = lx;
  d.ylength   = ly;
  d.zlength   = lz;
  d.dx        = lx / im;
  d.dy        = ly / jm;
  d.dz        = lz / km;

  return d;
}

typedef struct {
  double *Ax, *Ay, *Az, *Lambda;
  size_t size;
} FieldsType;

static FieldsType makeFields(const GeometryDomainType *d)
{
  FieldsType f;
  f.size = (size_t)(d->imaxLocal + 2) * (d->jmaxLocal + 2) * (d->kmaxLocal + 2);
  f.Ax     = calloc(f.size, sizeof(double));
  f.Ay     = calloc(f.size, sizeof(double));
  f.Az     = calloc(f.size, sizeof(double));
  f.Lambda = calloc(f.size, sizeof(double));
  return f;
}

static void freeFields(FieldsType *f)
{
  free(f->Ax);
  free(f->Ay);
  free(f->Az);
  free(f->Lambda);
}

/* Solid volume the apertures represent, in physical units. */
static double solidVolume(const GeometryDomainType *d, const FieldsType *f)
{
  double cell = d->dx * d->dy * d->dz;
  double sum  = 0.0;

  for (int k = 1; k < d->kmaxLocal + 1; k++) {
    for (int j = 1; j < d->jmaxLocal + 1; j++) {
      for (int i = 1; i < d->imaxLocal + 1; i++) {
        sum += (1.0 - AT(f->Lambda, i, j, k, d->imaxLocal, d->jmaxLocal)) * cell;
      }
    }
  }

  return sum;
}

/* No cell may be solid and still have an open face. */
static int openFacesOnSolidCells(const GeometryDomainType *d, const FieldsType *f)
{
  int im = d->imaxLocal, jm = d->jmaxLocal;
  int bad = 0;

  for (int k = 0; k < d->kmaxLocal + 2; k++) {
    for (int j = 0; j < jm + 2; j++) {
      for (int i = 0; i < im + 2; i++) {
        if (AT(f->Lambda, i, j, k, im, jm) > 0.0) {
          continue;
        }

        if (AT(f->Ax, i, j, k, im, jm) > 0.0 || AT(f->Ay, i, j, k, im, jm) > 0.0 ||
            AT(f->Az, i, j, k, im, jm) > 0.0) {
          ++bad;
        }
        if (i > 0 && AT(f->Ax, i - 1, j, k, im, jm) > 0.0) {
          ++bad;
        }
        if (j > 0 && AT(f->Ay, i, j - 1, k, im, jm) > 0.0) {
          ++bad;
        }
        if (k > 0 && AT(f->Az, i, j, k - 1, im, jm) > 0.0) {
          ++bad;
        }
      }
    }
  }

  return bad;
}

int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;

  CHECK_BEGIN("geometry");

  /* An obstacle-free spec leaves everything open. */
  {
    GeometryDomainType d = makeDomain(IMAX, JMAX, KMAX, 1.0, 1.0, 1.0);
    FieldsType f         = makeFields(&d);

    GeometrySpecType spec;
    geometryParseSpec(&spec, NULL);
    CHECK_TRUE(spec.kind == GEOMETRY_NONE, "an absent geometryFile is not obstacle-free");

    geometryProduce(&spec, &d, f.Ax, f.Ay, f.Az, f.Lambda);

    int closed = 0;
    for (size_t i = 0; i < f.size; i++) {
      if (f.Ax[i] != 1.0 || f.Ay[i] != 1.0 || f.Az[i] != 1.0 || f.Lambda[i] != 1.0) {
        ++closed;
      }
    }
    CHECK_TRUE(closed == 0, "a null producer closed %d entries", closed);

    freeFields(&f);
  }

  /* Analytic sphere: the represented solid volume must match the analytic one
   * to within roughly one cell layer over the surface. */
  {
    double lx = 4.0, radius = 1.0;
    GeometryDomainType d = makeDomain(IMAX, JMAX, KMAX, lx, lx, lx);
    FieldsType f         = makeFields(&d);

    GeometrySpecType spec;
    geometryParseSpec(&spec, "sphere:2.0,2.0,2.0,1.0");
    CHECK_TRUE(spec.kind == GEOMETRY_SPHERE, "sphere spec was not parsed");
    CHECK_NEAR(spec.radius, radius, 0.0, "sphere radius was not parsed");

    geometryProduce(&spec, &d, f.Ax, f.Ay, f.Az, f.Lambda);

    double got    = solidVolume(&d, &f);
    double want   = 4.0 / 3.0 * M_PI * radius * radius * radius;
    /* One cell layer over the sphere's surface, which is the accuracy a binary
     * aperture model can have. */
    double layer  = 4.0 * M_PI * radius * radius * d.dx;

    CHECK_TRUE(fabs(got - want) < layer,
        "sphere volume %.6f differs from the analytic %.6f by more than one cell "
        "layer %.6f",
        got,
        want,
        layer);

    CHECK_TRUE(openFacesOnSolidCells(&d, &f) == 0,
        "a solid cell kept an open face for the analytic sphere");

    /* Fully fluid and fully solid cells, checked at named places. */
    int im = d.imaxLocal, jm = d.jmaxLocal;
    int cc = IMAX / 2;
    CHECK_NEAR(AT(f.Lambda, cc, cc, cc, im, jm), 0.0, 0.0,
        "the cell at the sphere centre is not fully solid");
    CHECK_NEAR(AT(f.Lambda, 1, 1, 1, im, jm), 1.0, 0.0,
        "the corner cell, far from the sphere, is not fully fluid");
    CHECK_NEAR(AT(f.Ax, 1, 1, 1, im, jm), 1.0, 0.0,
        "a face far from the sphere is not fully open");

    freeFields(&f);
  }

  /* A zero-thickness plate closes a plane of faces and leaves both sides
   * fluid, which is the body a cell-type model cannot express at all. */
  {
    GeometryDomainType d = makeDomain(IMAX, JMAX, KMAX, 1.0, 1.0, 1.0);
    FieldsType f         = makeFields(&d);

    /* Normal to x, sitting on the face between cells 11 and 12. */
    double at = 12.0 * d.dx;
    char spec_s[128];
    snprintf(spec_s, sizeof(spec_s), "plate-x:%.17g,0.2,0.2,0.8,0.8", at);

    GeometrySpecType spec;
    geometryParseSpec(&spec, spec_s);
    CHECK_TRUE(spec.kind == GEOMETRY_PLATE, "plate spec was not parsed");

    geometryProduce(&spec, &d, f.Ax, f.Ay, f.Az, f.Lambda);

    int im = d.imaxLocal, jm = d.jmaxLocal;
    int mid = JMAX / 2;

    CHECK_NEAR(AT(f.Ax, 12, mid, mid, im, jm), 0.0, 0.0,
        "the face on the plate is not closed");
    CHECK_NEAR(AT(f.Lambda, 12, mid, mid, im, jm), 1.0, 0.0,
        "the cell on one side of the plate is not fluid");
    CHECK_NEAR(AT(f.Lambda, 13, mid, mid, im, jm), 1.0, 0.0,
        "the cell on the other side of the plate is not fluid");

    int solids = 0;
    for (size_t i = 0; i < f.size; i++) {
      if (f.Lambda[i] == 0.0) {
        ++solids;
      }
    }
    CHECK_TRUE(solids == 0, "a zero-thickness plate produced %d solid cells", solids);

    freeFields(&f);
  }

  /* A voxel volume: orientation, sampling, and agreement with the analytic body
   * it was rasterized from. */
  {
    double lx            = 4.0;
    GeometryDomainType d = makeDomain(IMAX, JMAX, KMAX, lx, lx, lx);
    FieldsType f         = makeFields(&d);

    GeometrySpecType spec;
    geometryParseSpec(&spec, "tests/geom/sphere.vox");
    CHECK_TRUE(spec.kind == GEOMETRY_VOXEL, "a path was not read as a voxel volume");

    geometryVoxelLoad(spec.file, &d);

    int nx, ny, nz;
    geometryVolumeSize(&nx, &ny, &nz);
    CHECK_TRUE(nx == 128 && ny == 128 && nz == 128,
        "voxel dimensions read as %d x %d x %d, expected 128 x 128 x 128 (regenerate "
        "with tests/make-geom.sh)",
        nx,
        ny,
        nz);
    CHECK_TRUE(geometryChecksum() != 0, "no checksum was computed");

    geometryProduce(&spec, &d, f.Ax, f.Ay, f.Az, f.Lambda);

    double got   = solidVolume(&d, &f);
    double want  = 4.0 / 3.0 * M_PI;
    double layer = 4.0 * M_PI * d.dx;

    CHECK_TRUE(fabs(got - want) < layer,
        "rasterized sphere volume %.6f differs from the analytic %.6f by more than "
        "one cell layer %.6f",
        got,
        want,
        layer);

    CHECK_TRUE(openFacesOnSolidCells(&d, &f) == 0,
        "a solid cell kept an open face for the rasterized sphere");

    geometryVoxelFree();
    freeFields(&f);
  }

  /* Orientation: a box in one named octant must come back in that octant, which
   * is what catches a flipped or transposed axis. */
  {
    double lx            = 4.0;
    GeometryDomainType d = makeDomain(IMAX, JMAX, KMAX, lx, lx, lx);
    FieldsType f         = makeFields(&d);

    /* Low x, low y, high z. */
    GeometrySpecType spec;
    geometryParseSpec(&spec, "box:0.1,0.1,3.0,1.0,1.0,3.9");
    geometryProduce(&spec, &d, f.Ax, f.Ay, f.Az, f.Lambda);

    int im = d.imaxLocal, jm = d.jmaxLocal;
    int lowLowHigh = 0, elsewhere = 0;

    for (int k = 1; k < KMAX + 1; k++) {
      for (int j = 1; j < JMAX + 1; j++) {
        for (int i = 1; i < IMAX + 1; i++) {
          if (AT(f.Lambda, i, j, k, im, jm) != 0.0) {
            continue;
          }

          int lowX  = (i - 1) < IMAX / 2;
          int lowY  = (j - 1) < JMAX / 2;
          int highZ = (k - 1) >= KMAX / 2;

          if (lowX && lowY && highZ) {
            ++lowLowHigh;
          } else {
            ++elsewhere;
          }
        }
      }
    }

    CHECK_TRUE(lowLowHigh > 0, "the box produced no solid cells in its own octant");
    CHECK_TRUE(elsewhere == 0,
        "%d solid cells landed outside the octant the box was placed in, so an axis "
        "is flipped or transposed",
        elsewhere);

    freeFields(&f);
  }

  /* The same volume on two grids differing by a factor of two must represent
   * the same body, which is what the resolution rule exists to guarantee. */
  {
    double lx = 4.0;
    double volumes[2];

    for (int pass = 0; pass < 2; pass++) {
      int n                = pass == 0 ? 8 : 16;
      GeometryDomainType d = makeDomain(n, n, n, lx, lx, lx);
      FieldsType f         = makeFields(&d);

      GeometrySpecType spec;
      geometryParseSpec(&spec, "tests/geom/sphere.vox");
      geometryVoxelLoad(spec.file, &d);
      geometryProduce(&spec, &d, f.Ax, f.Ay, f.Az, f.Lambda);

      volumes[pass] = solidVolume(&d, &f);

      geometryVoxelFree();
      freeFields(&f);
    }

    /* The coarser grid's sampling tolerance is a cell layer on that grid. */
    double coarseLayer = 4.0 * M_PI * (lx / 8.0);

    CHECK_TRUE(fabs(volumes[0] - volumes[1]) < coarseLayer,
        "the same volume represents %.6f on a coarse grid and %.6f on a fine one, "
        "more than the coarse sampling tolerance %.6f apart",
        volumes[0],
        volumes[1],
        coarseLayer);
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
