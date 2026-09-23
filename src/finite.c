/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * The one test in the code that has to see a NaN, in a file of its own.
 *
 * Everything else is compiled with -ffast-math, which tells the compiler that
 * no value is ever a NaN or an infinity. Under that assumption it may fold
 * isfinite(x) to true, and it may equally fold a test on the bits of x, since a
 * finite double's exponent is never all ones. GCC does the first. So the test
 * is compiled here with -fno-finite-math-only (the Makefile says so for this
 * object alone) and called out of line, where the caller's assumptions cannot
 * reach it.
 */
#include <math.h>
#include <stdbool.h>

#include "solver.h"

bool solveResidualIsFinite(double res) { return isfinite(res); }
