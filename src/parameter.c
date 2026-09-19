/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parameter.h"
#include "util.h"
#define MAXLINE 4096

void initParameter(Parameter *param)
{
  param->xlength           = 1.0;
  param->ylength           = 1.0;
  param->zlength           = 1.0;
  param->imax              = 100;
  param->jmax              = 100;
  param->kmax              = 100;
  param->itermax           = 1000;
  param->eps               = 0.0001;
  param->omg               = 1.7;
  param->re                = 100.0;
  param->gamma             = 0.9;
  param->tau               = 0.5;
  param->levels            = 5;
  param->presmooth         = 5;
  param->postsmooth        = 5;
  param->geometryFile      = NULL;
  param->numberOfParticles = 0;
  param->startTime         = 0.0;
  param->injectTimePeriod  = 0.0;
  param->writeTimePeriod   = 0.0;
  param->x1 = param->y1 = param->z1 = 0.0;
  param->x2 = param->y2 = param->z2 = 0.0;
}

void readParameter(Parameter *param, const char *filename)
{
  FILE *fp = fopen(filename, "r");
  char line[MAXLINE];
  int i;

  if (!fp) {
    fprintf(stderr, "Could not open parameter file: %s\n", filename);
    exit(EXIT_FAILURE);
  }

  while (!feof(fp)) {
    line[0] = '\0';
    if (fgets(line, MAXLINE, fp) == NULL && ferror(fp) != 0) {
      fprintf(stderr, "Error in fgets function\n");
    }

    for (i = 0; line[i] != '\0' && line[i] != '#'; i++)
      ;
    line[i] = '\0';

    /* Split on any whitespace, not on spaces alone, and trim what is left.
     * A value used to keep whatever trailing whitespace the line carried, so a
     * line such as "name dcavity" with no trailing comment produced the string
     * "dcavity\n". Every setup-specific boundary condition is selected by
     * strcmp on that name, so it silently did nothing -- the shipped setups
     * only worked because each happened to have a trailing comment for the
     * '#' strip above to remove. */
    char *tok = strtok(line, " \t\n\r\f\v");
    char *val = strtok(NULL, " \t\n\r\f\v");

#define PARSE_PARAM(p, f)                                                                \
  if (strncmp(tok, #p, sizeof(#p) / sizeof(#p[0]) - 1) == 0) {                           \
    param->p = f(val);                                                                   \
  }
#define PARSE_STRING(p) PARSE_PARAM(p, strdup)
#define PARSE_INT(p) PARSE_PARAM(p, atoi)
#define PARSE_REAL(p) PARSE_PARAM(p, atof)

    if (tok != NULL && val != NULL) {
      PARSE_REAL(xlength);
      PARSE_REAL(ylength);
      PARSE_REAL(zlength);
      PARSE_INT(imax);
      PARSE_INT(jmax);
      PARSE_INT(kmax);
      PARSE_INT(itermax);
      PARSE_INT(levels);
      PARSE_INT(presmooth);
      PARSE_INT(postsmooth);
      PARSE_REAL(eps);
      PARSE_REAL(omg);
      PARSE_REAL(re);
      PARSE_REAL(tau);
      PARSE_REAL(gamma);
      PARSE_REAL(dt);
      PARSE_REAL(te);
      PARSE_REAL(gx);
      PARSE_REAL(gy);
      PARSE_REAL(gz);
      PARSE_STRING(name);
      PARSE_STRING(geometryFile);
      PARSE_INT(bcLeft);
      PARSE_INT(bcRight);
      PARSE_INT(bcBottom);
      PARSE_INT(bcTop);
      PARSE_INT(bcFront);
      PARSE_INT(bcBack);
      PARSE_REAL(u_init);
      PARSE_REAL(v_init);
      PARSE_REAL(w_init);
      PARSE_REAL(p_init);
      PARSE_INT(numberOfParticles);
      PARSE_REAL(startTime);
      PARSE_REAL(injectTimePeriod);
      PARSE_REAL(writeTimePeriod);
      PARSE_REAL(x1);
      PARSE_REAL(y1);
      PARSE_REAL(z1);
      PARSE_REAL(x2);
      PARSE_REAL(y2);
      PARSE_REAL(z2);
    }
  }

  fclose(fp);
}

void printParameter(Parameter *param)
{
  printf("Parameters for %s\n", param->name);
  printf("Boundary conditions Left:%d Right:%d Bottom:%d Top:%d Front:%d "
         "Back:%d\n",
      param->bcLeft,
      param->bcRight,
      param->bcBottom,
      param->bcTop,
      param->bcFront,
      param->bcBack);
  printf("\tReynolds number: %.2f\n", param->re);
  printf("\tInit arrays: U:%.2f V:%.2f W:%.2f P:%.2f\n",
      param->u_init,
      param->v_init,
      param->w_init,
      param->p_init);
  printf("Geometry data:\n");
  printf("\tDomain box size (x, y, z): %.2f, %.2f, %.2f\n",
      param->xlength,
      param->ylength,
      param->zlength);
  printf("\tCells (x, y, z): %d, %d, %d\n", param->imax, param->jmax, param->kmax);
  printf("Timestep parameters:\n");
  printf("\tDefault stepsize: %.2f, Final time %.2f\n", param->dt, param->te);
  printf("\tTau factor: %.2f\n", param->tau);
  printf("Iterative solver parameters:\n");
  printf("\tMax iterations: %d\n", param->itermax);
  printf("\tepsilon (stopping tolerance) : %f\n", param->eps);
  printf("\tgamma (stopping tolerance) : %f\n", param->gamma);
  printf("\tomega (SOR relaxation): %f\n", param->omg);
}
