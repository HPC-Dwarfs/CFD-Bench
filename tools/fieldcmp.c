/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Compare two field dumps written by src/fielddump.c.
 *
 *   tools/fieldcmp <a.dump> <b.dump> [tolerance] [--l2]
 *
 * Reports the max absolute and L2 differences per field. Exits 0 when every
 * field is within the tolerance, 1 when one is not, and 2 when the files
 * cannot be compared at all. With no tolerance given the comparison is exact,
 * so a file compared against itself reports zero and exits 0.
 *
 * --l2 judges by the L2 difference rather than the largest one. Two runs of a
 * multi-step simulation that differ only in the order their global sums were
 * formed agree closely almost everywhere and can differ by more at a single
 * point, so the L2 is the honest measure of whether they are the same run.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NFIELDS 4

static const char *const FIELD_NAME[NFIELDS] = { "p", "u", "v", "w" };

static double *readDump(const char *path, int32_t dims[3], size_t *count)
{
  FILE *fp = fopen(path, "rb");

  if (fp == NULL) {
    fprintf(stderr, "fieldcmp: cannot open %s\n", path);
    exit(2);
  }

  char magic[8];
  if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, "NUSIFD01", 8) != 0) {
    fprintf(stderr, "fieldcmp: %s is not a NUSIFD01 dump\n", path);
    exit(2);
  }

  if (fread(dims, sizeof(int32_t), 3, fp) != 3) {
    fprintf(stderr, "fieldcmp: %s has a truncated header\n", path);
    exit(2);
  }

  *count       = (size_t)dims[0] * dims[1] * dims[2];
  size_t total = NFIELDS * (*count);
  double *data = malloc(total * sizeof(double));

  if (data == NULL) {
    fprintf(stderr, "fieldcmp: out of memory reading %s\n", path);
    exit(2);
  }

  if (fread(data, sizeof(double), total, fp) != total) {
    fprintf(stderr, "fieldcmp: %s is truncated\n", path);
    exit(2);
  }

  fclose(fp);
  return data;
}

int main(int argc, char **argv)
{
  if (argc < 3 || argc > 5) {
    fprintf(stderr, "Usage: %s <a.dump> <b.dump> [tolerance] [--l2]\n", argv[0]);
    return 2;
  }

  double tol = 0.0;
  int useL2  = 0;

  for (int i = 3; i < argc; i++) {
    if (strcmp(argv[i], "--l2") == 0) {
      useL2 = 1;
    } else {
      tol = atof(argv[i]);
    }
  }

  int32_t dimsA[3], dimsB[3];
  size_t countA, countB;

  double *a = readDump(argv[1], dimsA, &countA);
  double *b = readDump(argv[2], dimsB, &countB);

  if (dimsA[0] != dimsB[0] || dimsA[1] != dimsB[1] || dimsA[2] != dimsB[2]) {
    fprintf(stderr,
        "fieldcmp: grid mismatch, %dx%dx%d against %dx%dx%d\n",
        dimsA[0],
        dimsA[1],
        dimsA[2],
        dimsB[0],
        dimsB[1],
        dimsB[2]);
    return 2;
  }

  int failed = 0;

  for (int f = 0; f < NFIELDS; f++) {
    const double *pa = a + (size_t)f * countA;
    const double *pb = b + (size_t)f * countA;

    double maxDiff = 0.0, sumSq = 0.0;
    size_t argMax = 0;

    for (size_t i = 0; i < countA; i++) {
      double d = fabs(pa[i] - pb[i]);
      if (d > maxDiff) {
        maxDiff = d;
        argMax  = i;
      }
      sumSq += d * d;
    }

    double l2 = sqrt(sumSq / (double)countA);
    int bad   = useL2 ? (l2 > tol) : (maxDiff > tol);

    size_t k  = argMax / ((size_t)dimsA[0] * dimsA[1]);
    size_t j  = (argMax / dimsA[0]) % dimsA[1];
    size_t i  = argMax % dimsA[0];

    printf("%s: max %.17g at (%zu,%zu,%zu), L2 %.17g%s\n",
        FIELD_NAME[f],
        maxDiff,
        i,
        j,
        k,
        l2,
        bad ? "   FAIL" : "");

    if (bad) {
      failed = 1;
    }
  }

  free(a);
  free(b);
  return failed;
}
