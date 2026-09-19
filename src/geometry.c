/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.h"
#include "geometry-voxel.h"

/*
 * Sub-cell sampling resolution. A cell's volume fraction is counted over
 * SUBSAMPLES^3 points and a face aperture over SUBSAMPLES^2, then rounded.
 * The counting is what will produce fractional values in a later phase; only
 * the rounding at the end is specific to binary apertures.
 */
#define SUBSAMPLES 4

#define G(v, i, j, k)                                                                    \
  v[(size_t)(k) * (imaxLocal + 2) * (jmaxLocal + 2) + (size_t)(j) * (imaxLocal + 2) +    \
      (size_t)(i)]

/* Is the point solid under this spec? The one place geometry is defined. */
static int isSolid(const GeometrySpecType *spec, double x, double y, double z)
{
  switch (spec->kind) {
  case GEOMETRY_NONE:
    return 0;

  case GEOMETRY_VOXEL:
    return geometryVoxelIsSolid(x, y, z);

  case GEOMETRY_SPHERE: {
    double ddx = x - spec->xCenter;
    double ddy = y - spec->yCenter;
    double ddz = z - spec->zCenter;
    return (ddx * ddx + ddy * ddy + ddz * ddz) <= spec->radius * spec->radius;
  }

  case GEOMETRY_CYLINDER: {
    /* Infinite in its own axis; the domain boundaries cap it. */
    double a, b;
    switch (spec->axis) {
    case AXIS_X:
      a = y - spec->yCenter;
      b = z - spec->zCenter;
      break;
    case AXIS_Y:
      a = x - spec->xCenter;
      b = z - spec->zCenter;
      break;
    default:
      a = x - spec->xCenter;
      b = y - spec->yCenter;
      break;
    }
    return (a * a + b * b) <= spec->radius * spec->radius;
  }

  case GEOMETRY_BOX:
    return (x >= spec->x0 && x <= spec->x1 && y >= spec->y0 && y <= spec->y1 &&
            z >= spec->z0 && z <= spec->z1);

  case GEOMETRY_PLATE:
    /* A zero-thickness plate is never "solid" as a volume: it closes faces and
     * leaves the cells on both sides fluid. Handled in the face sampling below,
     * not here. */
    return 0;
  }

  return 0;
}

/*
 * A plate has no volume, so it cannot be found by sampling points. It is
 * defined instead as a rectangle lying in a grid-aligned plane, and it closes
 * exactly the faces that plane passes through. Returns 1 when the face whose
 * centre is (x, y, z) and whose normal is `axis` lies on the plate.
 */
static int faceOnPlate(const GeometrySpecType *spec, GeometryAxisType axis, double x,
    double y, double z, double dx, double dy, double dz)
{
  if (spec->kind != GEOMETRY_PLATE || spec->axis != axis) {
    return 0;
  }

  /* Within half a cell of the plane counts as on it, so that the plate lands on
   * the nearest face rather than falling between two. */
  switch (axis) {
  case AXIS_X:
    return fabs(x - spec->x0) <= 0.5 * dx && y >= spec->y0 && y <= spec->y1 &&
           z >= spec->z0 && z <= spec->z1;
  case AXIS_Y:
    return fabs(y - spec->y0) <= 0.5 * dy && x >= spec->x0 && x <= spec->x1 &&
           z >= spec->z0 && z <= spec->z1;
  default:
    return fabs(z - spec->z0) <= 0.5 * dz && x >= spec->x0 && x <= spec->x1 &&
           y >= spec->y0 && y <= spec->y1;
  }
}

/* Fraction of a cell box that is solid, by sub-sampling. */
static double solidFractionCell(const GeometrySpecType *spec, double x0, double y0,
    double z0, double dx, double dy, double dz)
{
  int solid = 0;

  for (int kk = 0; kk < SUBSAMPLES; kk++) {
    for (int jj = 0; jj < SUBSAMPLES; jj++) {
      for (int ii = 0; ii < SUBSAMPLES; ii++) {
        double x = x0 + (ii + 0.5) * dx / SUBSAMPLES;
        double y = y0 + (jj + 0.5) * dy / SUBSAMPLES;
        double z = z0 + (kk + 0.5) * dz / SUBSAMPLES;
        solid += isSolid(spec, x, y, z);
      }
    }
  }

  return (double)solid / (SUBSAMPLES * SUBSAMPLES * SUBSAMPLES);
}

/* Fraction of a face rectangle that is solid, by sub-sampling. */
static double solidFractionFace(const GeometrySpecType *spec, GeometryAxisType axis,
    double x, double y, double z, double dx, double dy, double dz)
{
  int solid = 0;

  for (int bb = 0; bb < SUBSAMPLES; bb++) {
    for (int aa = 0; aa < SUBSAMPLES; aa++) {
      double px = x, py = y, pz = z;

      switch (axis) {
      case AXIS_X:
        py = y + ((aa + 0.5) / SUBSAMPLES - 0.5) * dy;
        pz = z + ((bb + 0.5) / SUBSAMPLES - 0.5) * dz;
        break;
      case AXIS_Y:
        px = x + ((aa + 0.5) / SUBSAMPLES - 0.5) * dx;
        pz = z + ((bb + 0.5) / SUBSAMPLES - 0.5) * dz;
        break;
      default:
        px = x + ((aa + 0.5) / SUBSAMPLES - 0.5) * dx;
        py = y + ((bb + 0.5) / SUBSAMPLES - 0.5) * dy;
        break;
      }

      solid += isSolid(spec, px, py, pz);
    }
  }

  return (double)solid / (SUBSAMPLES * SUBSAMPLES);
}

void geometryProduce(const GeometrySpecType *spec,
    const GeometryDomainType *domain,
    double *Ax,
    double *Ay,
    double *Az,
    double *Lambda)
{
  int imaxLocal = domain->imaxLocal;
  int jmaxLocal = domain->jmaxLocal;
  int kmaxLocal = domain->kmaxLocal;

  double dx     = domain->dx;
  double dy     = domain->dy;
  double dz     = domain->dz;

  if (spec->kind == GEOMETRY_NONE) {
    size_t size = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);

    for (size_t i = 0; i < size; i++) {
      Ax[i] = Ay[i] = Az[i] = Lambda[i] = 1.0;
    }
    return;
  }

  /* Halo included, so that the geometry needs no exchange to be consistent.
   * A halo cell outside the domain is fluid: the domain boundary is the
   * boundary condition's business, not the obstacle's. */
  for (int k = 0; k < kmaxLocal + 2; k++) {
    for (int j = 0; j < jmaxLocal + 2; j++) {
      for (int i = 0; i < imaxLocal + 2; i++) {
        int gi     = i - 1 + domain->iOffset;
        int gj     = j - 1 + domain->jOffset;
        int gk     = k - 1 + domain->kOffset;

        /* Lower corner of the cell in physical coordinates. */
        double x0  = gi * dx;
        double y0  = gj * dy;
        double z0  = gk * dz;

        int inside = (gi >= 0 && gi < domain->imax && gj >= 0 && gj < domain->jmax &&
                      gk >= 0 && gk < domain->kmax);

        G(Lambda, i, j, k) =
            (inside && solidFractionCell(spec, x0, y0, z0, dx, dy, dz) >= 0.5) ? 0.0
                                                                              : 1.0;

        /* Face centres: the x-face of this cell is its upper x boundary. */
        double fxx = x0 + dx, fxy = y0 + 0.5 * dy, fxz = z0 + 0.5 * dz;
        double fyx = x0 + 0.5 * dx, fyy = y0 + dy, fyz = z0 + 0.5 * dz;
        double fzx = x0 + 0.5 * dx, fzy = y0 + 0.5 * dy, fzz = z0 + dz;

        G(Ax, i, j, k) =
            (solidFractionFace(spec, AXIS_X, fxx, fxy, fxz, dx, dy, dz) >= 0.5 ||
                faceOnPlate(spec, AXIS_X, fxx, fxy, fxz, dx, dy, dz))
                ? 0.0
                : 1.0;
        G(Ay, i, j, k) =
            (solidFractionFace(spec, AXIS_Y, fyx, fyy, fyz, dx, dy, dz) >= 0.5 ||
                faceOnPlate(spec, AXIS_Y, fyx, fyy, fyz, dx, dy, dz))
                ? 0.0
                : 1.0;
        G(Az, i, j, k) =
            (solidFractionFace(spec, AXIS_Z, fzx, fzy, fzz, dx, dy, dz) >= 0.5 ||
                faceOnPlate(spec, AXIS_Z, fzx, fzy, fzz, dx, dy, dz))
                ? 0.0
                : 1.0;
      }
    }
  }

  /*
   * A face is open only if both cells it separates are fluid.
   *
   * The volume fraction and the face apertures are sampled independently, so
   * rounding each to binary separately can disagree: a cell that is mostly
   * solid rounds to Lambda = 0 while a face along its edge, sampled on a plane
   * that falls in the fluid, rounds to open. That would leave a solid cell with
   * a face the flow can cross, and a pressure inside the body that reaches the
   * fluid -- exactly what the identity rows the solvers give solid cells assume
   * cannot happen.
   *
   * The converse is deliberately not imposed: a closed face between two fluid
   * cells is a zero-thickness plate, which is a body this model exists to
   * represent.
   */
  for (int k = 0; k < kmaxLocal + 2; k++) {
    for (int j = 0; j < jmaxLocal + 2; j++) {
      for (int i = 0; i < imaxLocal + 2; i++) {
        if (G(Lambda, i, j, k) > 0.0) {
          continue;
        }

        G(Ax, i, j, k) = 0.0;
        G(Ay, i, j, k) = 0.0;
        G(Az, i, j, k) = 0.0;

        if (i > 0) {
          G(Ax, i - 1, j, k) = 0.0;
        }
        if (j > 0) {
          G(Ay, i, j - 1, k) = 0.0;
        }
        if (k > 0) {
          G(Az, i, j, k - 1) = 0.0;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

static int parseAxisSuffix(const char *name, const char *prefix, GeometryAxisType *axis)
{
  size_t n = strlen(prefix);

  if (strncmp(name, prefix, n) != 0) {
    return 0;
  }

  if (name[n] == '-' && name[n + 2] == ':') {
    switch (name[n + 1]) {
    case 'x':
      *axis = AXIS_X;
      return 1;
    case 'y':
      *axis = AXIS_Y;
      return 1;
    case 'z':
      *axis = AXIS_Z;
      return 1;
    default:
      return 0;
    }
  }

  return 0;
}

void geometryParseSpec(GeometrySpecType *spec, const char *name)
{
  memset(spec, 0, sizeof(*spec));
  spec->kind = GEOMETRY_NONE;

  if (name == NULL || name[0] == '\0') {
    return;
  }

  if (strncmp(name, "sphere:", 7) == 0) {
    if (sscanf(name + 7,
            "%lf,%lf,%lf,%lf",
            &spec->xCenter,
            &spec->yCenter,
            &spec->zCenter,
            &spec->radius) != 4) {
      fprintf(stderr,
          "geometryFile: '%s' is not sphere:xc,yc,zc,r\n",
          name);
      exit(EXIT_FAILURE);
    }
    spec->kind = GEOMETRY_SPHERE;
    return;
  }

  if (parseAxisSuffix(name, "cylinder", &spec->axis)) {
    const char *args = strchr(name, ':') + 1;
    double a, b, r;

    if (sscanf(args, "%lf,%lf,%lf", &a, &b, &r) != 3) {
      fprintf(stderr,
          "geometryFile: '%s' is not cylinder-<axis>:c1,c2,r\n",
          name);
      exit(EXIT_FAILURE);
    }

    switch (spec->axis) {
    case AXIS_X:
      spec->yCenter = a;
      spec->zCenter = b;
      break;
    case AXIS_Y:
      spec->xCenter = a;
      spec->zCenter = b;
      break;
    default:
      spec->xCenter = a;
      spec->yCenter = b;
      break;
    }

    spec->radius = r;
    spec->kind   = GEOMETRY_CYLINDER;
    return;
  }

  if (strncmp(name, "box:", 4) == 0) {
    if (sscanf(name + 4,
            "%lf,%lf,%lf,%lf,%lf,%lf",
            &spec->x0,
            &spec->y0,
            &spec->z0,
            &spec->x1,
            &spec->y1,
            &spec->z1) != 6) {
      fprintf(stderr, "geometryFile: '%s' is not box:x0,y0,z0,x1,y1,z1\n", name);
      exit(EXIT_FAILURE);
    }
    spec->kind = GEOMETRY_BOX;
    return;
  }

  if (parseAxisSuffix(name, "plate", &spec->axis)) {
    const char *args = strchr(name, ':') + 1;
    double at, a0, b0, a1, b1;

    if (sscanf(args, "%lf,%lf,%lf,%lf,%lf", &at, &a0, &b0, &a1, &b1) != 5) {
      fprintf(stderr,
          "geometryFile: '%s' is not plate-<axis>:position,a0,b0,a1,b1\n",
          name);
      exit(EXIT_FAILURE);
    }

    switch (spec->axis) {
    case AXIS_X:
      spec->x0 = at;
      spec->y0 = a0;
      spec->z0 = b0;
      spec->y1 = a1;
      spec->z1 = b1;
      break;
    case AXIS_Y:
      spec->y0 = at;
      spec->x0 = a0;
      spec->z0 = b0;
      spec->x1 = a1;
      spec->z1 = b1;
      break;
    default:
      spec->z0 = at;
      spec->x0 = a0;
      spec->y0 = b0;
      spec->x1 = a1;
      spec->y1 = b1;
      break;
    }

    spec->kind = GEOMETRY_PLATE;
    return;
  }

  spec->kind = GEOMETRY_VOXEL;
  spec->file = name;
}

static const char *const AXIS_NAME[] = { "x", "y", "z" };

void geometryPrintHeader(
    const GeometrySpecType *spec, const GeometryDomainType *domain)
{
  printf("Obstacle geometry:\n");

  switch (spec->kind) {
  case GEOMETRY_NONE:
    printf("\tnone, the domain is obstacle-free\n");
    return;

  case GEOMETRY_VOXEL: {
    int nx, ny, nz;
    geometryVolumeSize(&nx, &ny, &nz);
    printf("\tvoxel volume %s\n", spec->file);
    printf("\tvoxels: %d x %d x %d\n", nx, ny, nz);
    printf("\tchecksum: %016llx\n", geometryChecksum());
    printf("\tvoxels per cell: %.2f x %.2f x %.2f\n",
        (double)nx / domain->imax,
        (double)ny / domain->jmax,
        (double)nz / domain->kmax);
    return;
  }

  case GEOMETRY_SPHERE:
    printf("\tanalytic sphere at (%g, %g, %g), radius %g\n",
        spec->xCenter,
        spec->yCenter,
        spec->zCenter,
        spec->radius);
    return;

  case GEOMETRY_CYLINDER:
    printf("\tanalytic cylinder along %s at (%g, %g, %g), radius %g\n",
        AXIS_NAME[spec->axis],
        spec->xCenter,
        spec->yCenter,
        spec->zCenter,
        spec->radius);
    return;

  case GEOMETRY_BOX:
    printf("\tanalytic box (%g, %g, %g) to (%g, %g, %g)\n",
        spec->x0,
        spec->y0,
        spec->z0,
        spec->x1,
        spec->y1,
        spec->z1);
    return;

  case GEOMETRY_PLATE:
    printf("\tanalytic plate normal to %s\n", AXIS_NAME[spec->axis]);
    return;
  }
}

/* ---------------------------------------------------------------------- */

/*
 * Connected components of the fluid region, by label propagation.
 *
 * Every fluid cell starts labelled with its own global index. A label spreads
 * to a face neighbour only through an open face, so the obstacle geometry is
 * what defines connectivity. Sweeping forward and backward in turn propagates a
 * label a long way per pass rather than one cell per pass, and the halo
 * exchange between passes carries labels across rank boundaries, so the answer
 * does not depend on the decomposition.
 *
 * At convergence each component's lowest-indexed cell is the only one still
 * carrying its own index, so counting those counts the components.
 */
int geometryValidateConnectivity(CommType *comm,
    const GeometryDomainType *domain,
    const double *Ax,
    const double *Ay,
    const double *Az,
    const double *Lambda,
    int abortOnFailure)
{
  int imaxLocal = domain->imaxLocal;
  int jmaxLocal = domain->jmaxLocal;
  int kmaxLocal = domain->kmaxLocal;

  size_t size   = (size_t)(imaxLocal + 2) * (jmaxLocal + 2) * (kmaxLocal + 2);
  double *label = malloc(size * sizeof(double));

  if (label == NULL) {
    fprintf(stderr, "geometry: cannot allocate the connectivity workspace\n");
    exit(EXIT_FAILURE);
  }

  /* SOLID is larger than any real label, so a min-propagation never adopts it. */
  const double SOLID = 1.0e18;

  for (size_t i = 0; i < size; i++) {
    label[i] = SOLID;
  }

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (G(Lambda, i, j, k) > 0.0) {
          double gi     = i - 1 + domain->iOffset;
          double gj     = j - 1 + domain->jOffset;
          double gk     = k - 1 + domain->kOffset;
          G(label, i, j, k) =
              (gk * domain->jmax + gj) * domain->imax + gi;
        }
      }
    }
  }

  for (int pass = 0; pass < 10000; pass++) {
    double changed = 0.0;

    commExchange(comm, label);

    for (int dir = 0; dir < 2; dir++) {
      int kBeg = dir ? kmaxLocal : 1, kEnd = dir ? 0 : kmaxLocal + 1;
      int step = dir ? -1 : 1;

      for (int k = kBeg; k != kEnd; k += step) {
        int jBeg = dir ? jmaxLocal : 1, jEnd = dir ? 0 : jmaxLocal + 1;

        for (int j = jBeg; j != jEnd; j += step) {
          int iBeg = dir ? imaxLocal : 1, iEnd = dir ? 0 : imaxLocal + 1;

          for (int i = iBeg; i != iEnd; i += step) {
            if (G(Lambda, i, j, k) == 0.0) {
              continue;
            }

            double best = G(label, i, j, k);

            /* A neighbour is reachable only through an open face. The x-face of
             * cell i-1 is the one between i-1 and i. */
            if (G(Ax, i - 1, j, k) > 0.0 && G(label, i - 1, j, k) < best) {
              best = G(label, i - 1, j, k);
            }
            if (G(Ax, i, j, k) > 0.0 && G(label, i + 1, j, k) < best) {
              best = G(label, i + 1, j, k);
            }
            if (G(Ay, i, j - 1, k) > 0.0 && G(label, i, j - 1, k) < best) {
              best = G(label, i, j - 1, k);
            }
            if (G(Ay, i, j, k) > 0.0 && G(label, i, j + 1, k) < best) {
              best = G(label, i, j + 1, k);
            }
            if (G(Az, i, j, k - 1) > 0.0 && G(label, i, j, k - 1) < best) {
              best = G(label, i, j, k - 1);
            }
            if (G(Az, i, j, k) > 0.0 && G(label, i, j, k + 1) < best) {
              best = G(label, i, j, k + 1);
            }

            if (best < G(label, i, j, k)) {
              G(label, i, j, k) = best;
              changed           = 1.0;
            }
          }
        }
      }
    }

    commReduceAll(&changed, MAX);

    if (changed == 0.0) {
      break;
    }
  }

  /* Each component's lowest-indexed cell still carries its own index. */
  double regions  = 0.0;
  double firstCell = -1.0;
  double secondCell = -1.0;

  for (int k = 1; k < kmaxLocal + 1; k++) {
    for (int j = 1; j < jmaxLocal + 1; j++) {
      for (int i = 1; i < imaxLocal + 1; i++) {
        if (G(Lambda, i, j, k) == 0.0) {
          continue;
        }

        double gi  = i - 1 + domain->iOffset;
        double gj  = j - 1 + domain->jOffset;
        double gk  = k - 1 + domain->kOffset;
        double own = (gk * domain->jmax + gj) * domain->imax + gi;

        if (G(label, i, j, k) == own) {
          regions += 1.0;

          if (firstCell < 0.0 || own < firstCell) {
            secondCell = firstCell;
            firstCell  = own;
          } else if (secondCell < 0.0 || own < secondCell) {
            secondCell = own;
          }
        }
      }
    }
  }

  commReduceAll(&regions, SUM);
  free(label);

  int count = (int)(regions + 0.5);

  if (count != 1 && abortOnFailure) {
    if (commIsMaster(comm)) {
      fprintf(stderr,
          "geometry: the fluid region is not connected. Found %d separate regions "
          "under face connectivity; the pressure solvers carry one constant in the "
          "null space, not %d. Each sealed pocket has to be removed or opened.\n",
          count,
          count);
    }

    /* Name a cell in each region, so the pockets can be found. */
    double myFirst = firstCell, mySecond = secondCell;
    for (int r = 0; r < comm->size; r++) {
      if (comm->rank == r && myFirst >= 0.0) {
        long idx = (long)myFirst;
        fprintf(stderr,
            "  region root at global cell (%ld, %ld, %ld)\n",
            idx % domain->imax,
            (idx / domain->imax) % domain->jmax,
            idx / ((long)domain->imax * domain->jmax));

        if (mySecond >= 0.0) {
          long idx2 = (long)mySecond;
          fprintf(stderr,
              "  region root at global cell (%ld, %ld, %ld)\n",
              idx2 % domain->imax,
              (idx2 / domain->imax) % domain->jmax,
              idx2 / ((long)domain->imax * domain->jmax));
        }
      }
    }

    exit(EXIT_FAILURE);
  }

  return count;
}
