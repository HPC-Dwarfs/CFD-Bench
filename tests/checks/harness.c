/* Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of NuSiF solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 *
 * Checks the check harness itself: that a driver links against the solver
 * objects, that MPI starts and stops cleanly at whatever rank count it was
 * launched with, and that a failing assertion really does produce a nonzero
 * exit status. Pass --fail to make it fail on purpose; tests/run-checks.sh uses
 * that to confirm it does not report success when a driver fails.
 */
#include <stdlib.h>

#include "check.h"

int main(int argc, char **argv)
{
  CommType comm;
  commInit(&comm, argc, argv);

  /* commFinalize frees the derived datatypes, which only commPartition creates,
   * so every driver has to partition before it finalizes even when it does not
   * care about the decomposition. */
  commPartition(&comm, 8, 8, 8);

  CHECK_BEGIN("harness");

  CHECK_TRUE(comm.size >= 1, "communicator size is %d", comm.size);
  CHECK_TRUE(comm.imaxLocal > 0 && comm.jmaxLocal > 0 && comm.kmaxLocal > 0,
      "local subdomain is %dx%dx%d",
      comm.imaxLocal,
      comm.jmaxLocal,
      comm.kmaxLocal);
  CHECK_TRUE(comm.rank >= 0 && comm.rank < comm.size,
      "rank %d out of range for size %d",
      comm.rank,
      comm.size);

  /* commIsMaster must agree with the rank on every rank, including serial
   * builds where it is the only rank. */
  CHECK_TRUE(commIsMaster(&comm) == (comm.rank == 0),
      "commIsMaster disagrees with rank %d",
      comm.rank);

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--fail") == 0) {
      CHECK_TRUE(0, "deliberate failure requested with --fail");
    }
  }

  int failures = CheckFailures;
  printf("%s: %d checks, %d failures\n", CheckName, CheckCount, failures);

  commFinalize(&comm);
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
