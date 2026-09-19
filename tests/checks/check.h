/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Minimal assertion helpers for the check drivers. A driver is an ordinary
 * program linked against the solver objects; it opens with CHECK_BEGIN, makes
 * assertions, and ends with CHECK_END, which returns 0 when every assertion
 * held and 1 otherwise.
 *
 * Assertions do not abort. A driver runs all of its checks so that one run
 * reports everything that is wrong, not only the first thing.
 */
#ifndef __CHECK_H_
#define __CHECK_H_

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "comm.h"

static int CheckFailures = 0;
static int CheckCount    = 0;
static const char *CheckName = "";

#define CHECK_BEGIN(name)                                                                \
  do {                                                                                   \
    CheckName = (name);                                                                  \
    printf("== %s ==\n", CheckName);                                                     \
  } while (0)

#define CHECK_TRUE(cond, ...)                                                            \
  do {                                                                                   \
    ++CheckCount;                                                                        \
    if (!(cond)) {                                                                       \
      ++CheckFailures;                                                                   \
      printf("FAIL %s:%d: ", __FILE__, __LINE__);                                         \
      printf(__VA_ARGS__);                                                               \
      printf("\n");                                                                      \
    }                                                                                    \
  } while (0)

/* Absolute tolerance. Use 0.0 to demand bit equality. */
#define CHECK_NEAR(got, want, tol, ...)                                                  \
  do {                                                                                   \
    ++CheckCount;                                                                        \
    double checkGot_  = (got);                                                           \
    double checkWant_ = (want);                                                          \
    if (!(fabs(checkGot_ - checkWant_) <= (tol))) {                                      \
      ++CheckFailures;                                                                   \
      printf("FAIL %s:%d: ", __FILE__, __LINE__);                                         \
      printf(__VA_ARGS__);                                                               \
      printf(" (got %.17g, want %.17g, tol %.17g)\n", checkGot_, checkWant_, (tol));     \
    }                                                                                    \
  } while (0)

#define CHECK_END()                                                                      \
  do {                                                                                   \
    printf("%s: %d checks, %d failures\n", CheckName, CheckCount, CheckFailures);        \
    return CheckFailures == 0 ? 0 : 1;                                                   \
  } while (0)

#endif // __CHECK_H_
