/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file. */

#include <stdio.h>

#ifndef __PROGRESS_H_
#define __PROGRESS_H_
extern void initProgress(double);
extern void printProgress(double);
extern void stopProgress(void);
extern FILE *initResidualWriter(void);
extern void writeResidual(FILE *, double, double);
#endif
