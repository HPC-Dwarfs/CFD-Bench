/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file. */

#include <stdio.h>

#include "comm.h"

#ifndef __PROGRESS_H_
#define __PROGRESS_H_
/* The bar is drawn by the master rank alone: every rank steps the same global
 * time, so one rank's view is the whole run's progress, and letting all of them
 * redraw would just interleave carriage returns on the shared stdout. */
extern void initProgress(CommType *, double);
extern void printProgress(double);
extern void stopProgress(void);
extern FILE *initResidualWriter(void);
extern void writeResidual(FILE *, double, double);
#endif
