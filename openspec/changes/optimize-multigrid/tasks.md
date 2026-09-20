# Tasks

## 1. Make the gate judge converged fields

- [ ] 1.1 Record the current state: run `tests/run-all.sh` and note every step passes, so this change has a known-good starting point; verify the recorded result is green end to end
- [ ] 1.2 Change `tests/check-solver-obstacle.sh` to solve each setup at a tolerance an order of magnitude below the setup's `eps` for the cross-solver comparison, and to keep comparing the resulting fields at `eps`; verify all four solvers still agree and the script exits zero
- [ ] 1.3 Confirm the change distinguishes what it is meant to: verify that `mg` at one level deeper than the setup asks for, which lands 5.1e-04 from `rb` under the old gate, now agrees, and that a deliberately wrong field still fails
- [ ] 1.4 Leave the rank-count L2 comparison as it is, and note in the script why it needs no tightening: it compares a solver against itself and carries no cross-solver stopping artifact; verify the rank-count section still passes for all four solvers

## 2. Solve the coarsest level

- [ ] 2.1 Replace the coarsest level's fixed `presmooth + postsmooth` relaxation with a much larger fixed sweep count, split evenly forward and backward so the symmetric shape keeps its `A^T A` form; verify `tests/checks/multigrid.c` still reports the cycle symmetric to machine precision with and without a body
- [ ] 2.2 Choose the sweep count from measurement rather than by guess: record the coarse-level residual reduction and the resulting cycle count on `sphere-baseline` for a few values, and pick the knee; verify the chosen value is recorded in the code next to the constant
- [ ] 2.3 Add a case to `tests/checks/multigrid.c` asserting the coarsest level's residual falls by orders of magnitude within one cycle, not by the factor a handful of sweeps would give; verify it passes at 1 and 4 ranks and fails with 2.1 reverted
- [ ] 2.4 Confirm the coarse solve performs identical work regardless of the residual it is given, so the preconditioner stays a fixed linear operator; verify `tests/checks/cg.c` still reports the `mg` preconditioner linear and symmetric
- [ ] 2.5 Measure the effect on its own, before any depth change: record `SOLVER=mg` cycles and wall clock on `sphere-baseline` at the setup's current `levels`; verify the cross-solver gate still passes

## 3. Use the depth the grids support

- [ ] 3.1 Raise `levels` to what each grid supports in `testcases/regression/sphere-baseline.par` (3 to 4), `testcases/regression/dcavity-baseline.par` (3 to 5), `testcases/flow/canal.par` (3 to 4), `testcases/flow/dcavity.par` (4 to 7) and `testcases/flow/schaefer-turek.par` (3 to 5), changing no grid dimension; verify each setup reports the level count it actually built and executes time steps
- [ ] 3.2 Leave `testcases/flow/karman.par`, `testcases/flow/backstep.par` and `testcases/regression/canal-baseline.par` alone, and record why: the first is limited by a grid dimension that halves to an odd number, the other two are already at their limit; verify `karman.par` still reports that the decomposition supports fewer levels than it requests
- [ ] 3.3 Add a case to `tests/checks/multigrid.c` asserting that a hierarchy built on a grid supporting more levels than requested honours the request, and that one requesting more than the decomposition allows is clamped and reports the number built; verify both at 1 and 4 ranks
- [ ] 3.4 Re-record the `mg` baselines the deeper hierarchies move, listing which setups changed and by how much; verify `tests/record-baseline.sh -v` is green afterwards and that the `rb`, `rbc` and `cg` baselines are untouched
- [ ] 3.5 Confirm the deeper hierarchies still solve the same system: verify `tests/check-solver-obstacle.sh` passes with all four solvers and that `mg` agrees with `rb`
- [ ] 3.6 Record the gain: `SOLVER=mg` cycles and wall clock on `sphere-baseline` before and after, and the same for `SOLVER=cg precon mg`; verify the converged fields agree with `rb`

## 4. A cycle shape the solver can choose

- [ ] 4.1 Add an `MgShapeType` field to `MultigridSpecType` and `MultigridType` selecting the symmetric or the fast shape, with every caller passing the symmetric one so nothing moves yet; verify the tree builds under all four solvers and `tests/checks/multigrid.c` reports no change
- [ ] 4.2 Branch the post-smoother direction, the coarsest solve's split and the restriction operator on the shape, keeping the eight-cell averaging restriction alongside the transposed one; verify the symmetric shape's results are bit-identical to before the branch was introduced
- [ ] 4.3 Move the equal-`presmooth`/`postsmooth` refusal behind the symmetric shape, since only that shape needs it; verify `tests/check-rejects.sh` still refuses unequal counts for a symmetric build and accepts them for a fast one
- [ ] 4.4 Have `solver-mg.c` request the fast shape and `precon-mg.c` the symmetric one; verify `tests/checks/cg.c` still reports the `mg` preconditioner symmetric and `tests/checks/multigrid.c` reports the solver's shape converging
- [ ] 4.5 Extend `tests/checks/multigrid.c` so the symmetry assertions apply to the symmetric shape and the fast shape gets its own assertion that it converges to the same field to the solve tolerance; verify both at 1 and 4 ranks

## 5. Decide whether the split is worth keeping

- [ ] 5.1 Measure the two shapes against each other at full depth and with the strengthened coarse solve: `SOLVER=mg` cycles and wall clock on `sphere-baseline` for the symmetric shape and for the fast shape at its best stable relaxation factor; verify both converge and agree with `rb` through the gate
- [ ] 5.2 Decide on that evidence, and record the decision in `design.md` with the numbers behind it: keep the split if the fast shape is meaningfully ahead, and remove it if it is not, reverting section 4 rather than shipping a second shape for its own sake; verify whichever tree results builds under all four solvers
- [ ] 5.3 Document the outcome in the README's Pressure Solvers section: the hierarchy depth setups should use, the coarse solve, and — if the split is kept — which shape each caller gets and why a solver does not need symmetry; verify a fresh `make SOLVER=mg` and `make SOLVER=cg` from a clean tree succeed
- [ ] 5.4 Run `tests/run-all.sh` end to end and compare against the state recorded in task 1.1; verify every step passes and that any difference is attributable to the deeper hierarchies, the coarse solve or the cycle shape
