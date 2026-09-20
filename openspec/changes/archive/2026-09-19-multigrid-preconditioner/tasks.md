# Tasks

## 1. Confine the coarse correction to the fluid

- [x] 1.1 Record the current state: run `tests/run-all.sh` and note that `check drivers, mg` fails on `tests/checks/obstacle.c` at 1 and 4 ranks with 576 solid cells nonzero, and that every other step passes; verify the failure reproduces at `HEAD` so the fix below is attributable
- [x] 1.2 Skip cells with zero volume fraction in `correct()` in `src/solver-mg.c`, reading the level's own `Lambda`; verify `tests/checks/obstacle.c` passes under `SOLVER=mg` at 1 and 4 ranks, where it failed before
- [x] 1.3 Add a case to `tests/checks/multigrid.c` that runs one cycle on a domain whose obstacle is unresolved at the coarse level and asserts every fully solid cell is exactly zero afterwards, including cells absent from the surface list; verify it fails with 1.2 reverted and passes with it applied
- [x] 1.4 Confirm the fix moves nothing else: verify `tests/record-baseline.sh -v` is still green and `tests/run-all.sh` now passes every step

## 2. Extract the hierarchy and the cycle

- [x] 2.1 Create `src/multigrid.h` declaring a context carrying the boundary condition, `omega`, the level count, the smoothing counts and the level array, plus entry points to build and free a hierarchy and to run one cycle; verify the tree still builds under all four solvers with nothing yet using it
- [x] 2.2 Move `MgLevelType`, `coarsenGeometry`, `countFluid`, `levelDescribe`, `zeroField`, `restrictMG`, `prolongate`, `correct`, `smooth`, `residualField` and `vcycle` from `src/solver-mg.c` to `src/multigrid.c` unchanged, taking the new context instead of `Solver *`; verify `make SOLVER=mg` links and the file is picked up with no Makefile change
- [x] 2.3 Reduce `src/solver-mg.c` to `initSolver` building the hierarchy and `solve` iterating cycles to the tolerance, keeping the hierarchy behind the existing opaque `mgLevels` pointer; verify `tests/checks/multigrid.c` passes unchanged at 1 and 4 ranks
- [x] 2.4 Move the `mgTest*` seam to `src/multigrid.c` and widen its guard from `defined(TEST) && defined(SOLVER_mg)` to any build that owns a hierarchy; verify `make SOLVER=mg tests` builds and the driver still passes
- [x] 2.5 Prove the extraction is a pure move: verify the `mg` baselines are bit-identical across it with `tests/record-baseline.sh -v`, and that `tests/checks/multigrid.c` and `tests/checks/obstacle.c` report the same counts as before

## 3. Make the cycle symmetric

- [x] 3.1 Give `smooth` a direction and have the post-smoother visit the two colours in the reverse order of the pre-smoother; verify a cycle still reduces the residual on the Poisson problem in `tests/checks/multigrid.c` by more than smoothing alone
- [x] 3.2 Add a `smoothOmega` parameter defaulting to 1.0 that the cycle uses in place of `omg`, which stays the SOR solvers' relaxation factor; verify a parameter file setting it is read back, that `sphere-baseline` no longer reaches NaN with the post-smoother reversed, and that `tests/checks/multigrid.c` passes
- [x] 3.3 Reject unequal `presmooth` and `postsmooth` at initialization, naming both values, the way `pressureBcInit` rejects a boundary code; verify `tests/check-rejects.sh` covers it and that equal counts are accepted
- [x] 3.4 Set `presmooth` and `postsmooth` to 5 in `canal.par` and `dcavity.par`, matching every other shipped setup and their own baseline variants; verify `tests/check-setups.sh` passes and both setups still execute time steps
- [x] 3.5 Replace the coarsest level's one-directional relaxation with `k` forward sweeps followed by `k` backward sweeps; verify the coarsest-level assertions in `tests/checks/multigrid.c` still hold
- [x] 3.6 Make each of prolongation's three one-dimensional passes a true 3/4-1/4 blend by updating the two cells of a pair from their own values, removing the aliasing that gives every other fine cell weights of 0.8125 and 0.1875; verify a coarse ramp prolongates to the exact trilinear result and that the existing constant-preservation checks still hold
- [x] 3.7 Rebuild restriction as the transpose of prolongation — the three one-dimensional passes reversed, then the eight-child sum, scaled by 1/8 — reading only face neighbours at a subdomain boundary; verify restriction of a constant is still that constant and that restriction followed by prolongation of a constant returns it
- [x] 3.8 Replace the single-target restriction assertion in `tests/checks/multigrid.c`, which the wider transpose footprint invalidates, with the transpose-pair identity: for arbitrary fine `f` and coarse `c`, `<R f, c>` equals `<f, P c>` up to the scale factor; verify it passes at 1 and 4 ranks and fails if either operator is changed alone
- [x] 3.9 Add a symmetry check to `tests/checks/multigrid.c` asserting `x·My == y·Mx` for one cycle from a zero guess, with and without an obstacle, judged relative to the magnitude of the two inner products; verify it passes at 1 and 4 ranks and fails with 3.1, 3.5 or 3.7 reverted
- [x] 3.10 Re-record the `mg` baselines the symmetric cycle moves, listing which setups changed and by how much; verify `tests/record-baseline.sh -v` is green afterwards and that the relaxation and `cg` baselines are untouched
- [x] 3.11 Confirm the new cycle still solves the same system: verify `tests/check-solver-obstacle.sh` passes with all four solvers and that `mg` still agrees with `rb` at the setup's tolerance

## 4. The multigrid preconditioner

- [x] 4.1 Create `src/precon-mg.c` exposing a `PreconType` whose `apply` zeroes `z` and runs one cycle from that zero guess, with `ctx` holding the hierarchy; verify it compiles into every build with no Makefile change
- [x] 4.2 Build the hierarchy in `src/solver-cg.c` when `precon mg` is selected, store it in the preconditioner's `ctx`, and accept `mg` where it was previously refused; verify `make SOLVER=cg` links and the banner reports the preconditioner and the level count
- [x] 4.3 Remove the `precon mg` rejection case from `tests/check-rejects.sh` and replace it with an acceptance case, keeping the unknown-value rejection; verify the script passes
- [x] 4.4 Extend `tests/checks/cg.c` with an `mg` case asserting the preconditioner is linear and symmetric, exactly as the `none` and `jacobi` cases do; verify it passes at 1 and 4 ranks
- [x] 4.5 Extend `tests/checks/cg.c` to assert the `mg` preconditioner leaves solid cells exactly zero, and that an `mg`-preconditioned solve converges to the exact discrete solution on the Dirichlet and Neumann Poisson configurations with and without an obstacle; verify both at 1 and 4 ranks
- [x] 4.6 Extend `tests/checks/cg.c` to assert the `mg` solve takes no more iterations than the `jacobi` solve of the same problem to the same tolerance; verify it holds on both Poisson configurations

## 5. Measure and integrate

- [x] 5.1 Add a refinement case to `tests/checks/cg.c` solving the same problem on a grid and on one refined in every direction, with `mg` and with `jacobi`, asserting the `mg` iteration count grows by a substantially smaller factor; verify it passes at 1 and 4 ranks and record both growth factors
- [x] 5.2 Measure and record, on `sphere-baseline`, the iteration counts and wall-clock of `SOLVER=mg`, `SOLVER=cg precon jacobi` and `SOLVER=cg precon mg`, so the cost of the symmetric cycle and of the two allreduces is visible rather than inferred; verify the three converged fields agree at the setup's tolerance
- [x] 5.3 Document the multigrid preconditioner, the equal-smoothing-count constraint and the measured refinement behaviour in the README's Pressure Solvers section; verify the `precon` table lists three values and that a fresh `make SOLVER=cg` from a clean tree succeeds
- [x] 5.4 Run `tests/run-all.sh` end to end and compare against the state recorded in task 1.1; verify every step passes and that every difference is attributable to the symmetric cycle or to the `correct()` fix
