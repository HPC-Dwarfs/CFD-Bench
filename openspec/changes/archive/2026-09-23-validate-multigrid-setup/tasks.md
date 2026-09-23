# Tasks

## 1. Refuse a setup that cannot produce a correct solve

- [x] 1.1 Record the current state: run `tests/run-all.sh` and note every step passes, so the change has a known-good starting point; verify the recorded result is green end to end
- [x] 1.2 Refuse `levels < 1` in `multigridBuild`, beside the existing equal-count refusal, naming the parameter and the value; verify `levels 0` and `levels -1` now exit nonzero with a message instead of the segfault they produce today, under both `SOLVER=mg` and `SOLVER=cg` with `precon mg`
- [x] 1.3 Leave the existing clamp of a depth larger than the decomposition supports exactly as it is, and note in the code why the two directions differ: a request too deep can be honoured approximately and reported, a request below one cannot be honoured at all; verify a setup asking for more levels than its grid allows still runs and still reports the number built
- [x] 1.4 Refuse `presmooth < 1` or `postsmooth < 1` in `multigridBuild`; verify `presmooth 0` and `presmooth -1` are refused rather than running to `itermax` with a residual four orders of magnitude above `eps`
- [x] 1.5 Refuse `smoothOmega` outside the open interval `(0, 2)` in `multigridBuild`, with a message that names `smoothOmega` specifically so it is not confused with `omg`; verify `0`, `-1`, `2.0` and `3.0` are refused and that `1.3`, the default, and `1.9` are accepted
- [x] 1.6 Refuse `omg` outside `(0, 2)` in `solverBaseInit`, which is the one place every solver variant passes through, with a message naming `omg`; verify `omg 3.0` is refused under `rb`, `rbc`, `mg` and `cg`, and that the shipped values 1.7 and 1.8 are accepted
- [x] 1.7 Add the refused values to `tests/check-rejects.sh`, following the existing `expect_reject` pattern of a nonzero exit, a message naming the problem and no `TIME` line; verify the script passes and that each new case fails if its guard is reverted
- [x] 1.8 Confirm no shipped setup is caught by any of the new refusals: verify `tests/check-setups.sh` passes and that each setup in `testcases/` still starts and executes time steps

## 2. Report how a solve stopped

- [x] 2.1 Establish the converged path is unchanged before changing it: record `SOLVER=mg`, `rb`, `rbc` and `cg` iteration counts, residuals and field dumps on `sphere-baseline`; verify these are the reference the next tasks compare against
- [x] 2.2 Change the multigrid loop at `src/solver-mg.c:87` to stop on a residual that is not finite as well as on the tolerance and the limit, and to report the three outcomes distinguishably; verify a converged solve reports and returns exactly what 2.1 recorded, and that `smoothOmega 3.0` now reports divergence instead of `took 13 cycles to reach nan`
- [x] 2.3 Make the same change in `src/solver-rb.c:78`; verify `omg 3.0` reports divergence rather than `took 552 iterations to reach nan`, and the converged path matches 2.1
- [x] 2.4 Make the same change in both `rbc` loops, `src/solver-rbc.c:197` and `:283` — the compressed and the natural layout — since missing one leaves the serial path silently unfixed; verify both by running the divergence case under an MPI build and under the serial build `check-solver-obstacle.sh` creates
- [x] 2.5 Make the same change in `src/solver-cg.c:468`; verify the converged path matches 2.1 and that `tests/checks/cg.c` still passes at 1 and 4 ranks
- [x] 2.6 Confirm the residual test is consistent across ranks: verify a diverging setup is reported as diverged by every rank at 4 ranks, since the residual is formed by a global reduction and a NaN on one rank must propagate to all
- [x] 2.7 Confirm an exhausted iteration budget is still non-fatal and still returns its iterate: verify `tests/checks/nullspace.c` passes, since it runs a deliberately incompatible right-hand side to `itermax` and compares the result against the handled case

## 3. Act on a diverged solve

- [x] 3.1 Have `main.c` treat a non-finite residual from `solve()` as fatal, before `adaptUV` and before any output is written, reporting the time step at which it diverged; verify a diverging setup exits nonzero and that no VTK file is produced from it
- [x] 3.2 Confirm the three outcomes are distinguishable to a reader and to a script: verify a converged run, a run that stops at its iteration limit and a diverged run produce different messages and that only the diverged one exits nonzero
- [x] 3.3 Add a diverging setup to `tests/check-rejects.sh`, or to a sibling check if a divergence that is detected mid-run rather than at initialization does not fit that script's shape; verify it fails the run and that the check passes

## 4. Confirm nothing else moved

- [x] 4.1 Verify `tests/record-baseline.sh -v` is green, so no converging solve changed its field
- [x] 4.2 Verify `tests/check-solver-obstacle.sh` passes with all four solvers, so the converged fields and the cross-solver agreement are untouched
- [x] 4.3 Run `tests/run-all.sh` end to end and compare against the state recorded in task 1.1; verify every step passes and that any difference is attributable to a new refusal or to the changed reporting of a solve that did not converge
