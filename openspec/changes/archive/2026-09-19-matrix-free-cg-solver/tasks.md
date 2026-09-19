# Tasks

## 1. Repair the acceptance gates

- [x] 1.1 Replace `NusifSolver-$TOOLCHAIN` with `CFD-Solver-$TOOLCHAIN` in `tests/check-setups.sh`, `tests/record-baseline.sh`, `tests/check-schaefer-turek.sh`, `tests/check-solver-obstacle.sh`, `tests/check-rejects.sh` and `tests/check-particles.sh`, including the `-test` and `NusifSolver-SERIAL-*` forms; verify `grep -rn NusifSolver tests/` returns nothing
- [x] 1.2 Replace the stale `NusifSolver-*` line in `.gitignore` with `CFD-Solver-*` and drop the narrower `/CFD-Solver-CLANG` it was masking, so the test binary and every other toolchain's binaries stop showing up as untracked; verify `git status` reports no build artefact after `make tests`
- [x] 1.3 Record the selected solver in a `$(BUILD_DIR)/solver.sel` stamp and make it a prerequisite of the target, the test target and every check binary, so selecting a variant whose object file is already current still relinks; verify that building `rb`, then `cg`, then `rb` again leaves a binary that reports `rb`, and that two consecutive `tests/run-all.sh` runs agree
- [x] 1.4 Run `tests/run-all.sh` on the current tree and record which steps pass, so this change has a known-good starting point to compare against; verify every step reports OK or a pre-existing failure that is written down

## 2. Matrix-free operator application

- [x] 2.1 Declare `pressureApplyOperator` and `pressureApplySurface` in `src/solver.h` next to the existing `PressureLevelType` block, documenting the negated sign convention from design.md decision 2; verify the tree still builds under all three existing solvers
- [x] 2.2 Implement `pressureApplySurface` in `src/solverbase.c` as a sibling of `pressureCorrectSurface`: one uncoloured pass over all `count` entries reading `aE..aB` and `lambda`, writing `y = 0` for solid entries; verify it compiles and that no existing check driver regresses
- [x] 2.3 Implement `pressureApplyOperator` in `src/solverbase.c` as halo exchange, boundary condition, a geometry-free 7-point bulk loop, then `pressureApplySurface`; verify `tests/run-checks.sh` still passes under `rb`, `rbc` and `mg`
- [x] 2.4 Add `tests/checks/operator.c` covering operator symmetry (`x·Ay == y·Ax`) with and without an obstacle, positive definiteness on the fluid subspace, the seven-point footprint from a unit vector next to an obstacle corner, and agreement with `pressureResidualNorm` at cut cells; verify it passes under all three existing solvers at 1 and 4 ranks
- [x] 2.5 Add profiler regions for the application's bulk and surface passes reusing `SWEEP_BULK` and `SWEEP_SURFACE`, and extend `tests/checks/sweeptime.c` to measure `pressureApplyOperator` across the none/sphere/lattice configurations; verify the bulk time agrees across the three within the driver's existing tolerance

## 3. Fluid-only null space

- [x] 3.1 Change `removeMean` in `src/discretization.c` to sum over cells with nonzero volume fraction, divide by the global fluid cell count, and leave solid cells at zero; verify a singular obstacle setup ends each step with every solid cell exactly zero
- [x] 3.2 Add a null-space case to `tests/checks/nullspace.c` for a fully enclosed setup containing an obstacle: solid cells stay zero after projection, the removed constant is the fluid mean, and all three existing solvers converge and agree; verify the driver passes at 1 and 4 ranks
- [x] 3.3 Re-record the baselines the fix moves with `tests/record-baseline.sh`, listing which setups changed and by how much; verify `tests/record-baseline.sh -v` is green afterwards

## 4. The CG solver variant

- [x] 4.1 Add the `PreconType` function-pointer seam to `src/solver.h` and a `void *cgState` field to `Solver` alongside `mgLevels`; verify the tree builds under the three existing solvers with the new field unused
- [x] 4.2 Add a `precon` string parameter to `src/parameter.h`/`.c` with `PARSE_STRING`, defaulting to `jacobi`, printed by `printParameter`; verify a parameter file setting `precon none` is read back correctly and that `presmooth` still parses
- [x] 4.3 Create `src/solver-cg.c` with `initSolver` calling `solverBaseInit`, building the surface list the way `solver-rb.c` does, and allocating `r`, `z`, `d`, `q` into `cgState`; verify `make SOLVER=cg` links and the banner reports the solver and its extra memory
- [x] 4.4 Add `commReduceAllN(v, count, op)` to `src/comm.c`/`.h` — the in-place all-reduce `commReduceAll` already is, over several values at once — because the existing `commReduce` reduces to the master only and CG needs its step lengths on every rank; implement fluid-only dot products on it, with count 2 for the fused `{r·r, r·z}` and count 1 for `d·q`; verify a unit-vector dot over a domain with solid cells returns the fluid count
- [x] 4.5 Implement the PCG recurrence with the negated system, warm start from the caller's `p`, stopping on `r·r / fluidCells < eps*eps` or `itermax`, printing `took N iterations to reach E` in the form `tests/check-solver-obstacle.sh` greps for; verify an obstacle-free `dcavity` run converges and reports a count
- [x] 4.6 Apply the fluid-indicator projection to `b` and the initial guess once, and to `z` each iteration, gated on `pressureBcIsSingular`; verify a fully enclosed run over many steps keeps the fluid mean pressure bounded

## 5. Preconditioners

- [x] 5.1 Implement the `none` preconditioner as a copy of `r` into `z` with solid cells zeroed; verify unpreconditioned CG converges on the Poisson problem with a known discrete solution
- [x] 5.2 Implement the `jacobi` preconditioner as a scalar bulk multiply by `1 / (2*(idx2+idy2+idz2))` followed by a surface-list override using `invDiag`, reading no geometry array; verify solid cells receive exactly zero and that `tests/checks/sweeptime.c` sees no geometry-dependent bulk cost
- [x] 5.3 Reject an unsupported `precon` value at initialization, naming the value, the way `pressureBcInit` rejects a boundary type; verify `tests/check-rejects.sh` covers it and that `precon mg` is rejected until the follow-up change lands

## 6. CG test driver

- [x] 6.1 Add the `#if defined(TEST) && defined(SOLVER_cg)` seam in `src/solver.h` and `src/solver-cg.c` exposing the preconditioner application, the iteration state and a single-iteration step, following the `mgTest*` pattern; verify `make SOLVER=cg tests` builds
- [x] 6.2 Add `tests/checks/cg.c` asserting preconditioner linearity and symmetry (`x·My == y·Mx`) for both `none` and `jacobi`; verify the driver passes at 1 and 4 ranks
- [x] 6.3 Extend `tests/checks/cg.c` with convergence to the exact discrete solution on the Dirichlet and Neumann Poisson configurations from `tests/checks/poisson.c`, with and without an obstacle; verify the converged field matches to the driver's tolerance
- [x] 6.4 Extend `tests/checks/cg.c` with the invariant checks the spec names: solid cells exactly zero at an intermediate iteration, and an unchanged iteration count and fluid field when solid volume is added away from the fluid region; verify both pass at 1 and 4 ranks
- [x] 6.5 Extend `tests/checks/cg.c` to assert that the `jacobi` solve takes no more iterations than the `none` solve of the same problem to the same tolerance; verify the assertion holds on both Poisson configurations

## 7. Integration

- [x] 7.1 Add `cg` to the solver loops in `tests/run-all.sh` and `tests/check-solver-obstacle.sh`, and give the rank-count comparison a within-one tolerance for Krylov solvers rather than strict equality; verify `tests/check-solver-obstacle.sh` passes with all four solvers
- [x] 7.2 Confirm CG's converged field agrees with `rb`, `rbc` and `mg` on the `sphere-baseline` setup via `tools/fieldcmp` at the setup's tolerance; verify the comparison passes in both directions against each solver
- [x] 7.3 Add `cg` to the `SOLVER` comment in `config.mk` and `mk/config-default.mk`, and document the solver, the `precon` parameter and the two-allreduce cost in the README's Pressure Solvers section; verify the README lists four solvers and that a fresh `make SOLVER=cg` from a clean tree succeeds
- [x] 7.4 Run `tests/run-all.sh` end to end and compare against the baseline recorded in task 1.4; verify every step passes and that any difference is attributable to the `removeMean` fix
