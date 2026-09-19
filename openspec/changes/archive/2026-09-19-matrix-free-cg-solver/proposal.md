# Proposal

## Why

A preconditioned CG solver is the stated endpoint of the embedded-boundary work:
that change built face apertures and volume fractions specifically so all
solvers share one symmetric positive definite operator, and its proposal names
CG as the reason. The operator is now in place and unused by any Krylov method —
the project ships three relaxation solvers and nothing that exploits symmetry.

CG also needs guarantees the relaxation solvers never did. It applies the
operator every iteration rather than once per solve, it forms global inner
products, and it cannot self-heal a violated invariant the way a colour sweep
does. Two things in the current code are fine for SOR and wrong for CG:
`removeMean` averages over solid cells and leaves them holding a nonzero value,
breaking the `p = 0` identity rows; and the only operator application that
exists, multigrid's `residualField`, reads four geometry streams per cell,
discarding the geometry-free-bulk property the benchmark's central claim rests
on.

## What Changes

- **New `SOLVER=cg` build variant** (`src/solver-cg.c`), a first-class solver
  alongside `rb`, `rbc` and `mg`, selected the same way and subject to the same
  cross-solver agreement checks.
- **A matrix-free operator application in `solverbase.c`**, shared by any solver
  that wants one: a geometry-blind 7-point bulk pass followed by an
  `O(obstacle surface)` correction over the existing surface list. Same split,
  same instruction mix and same memory streams as the SOR bulk sweep, so CG's
  SpMV is directly comparable to it.
- **A diagonal (Jacobi) preconditioner**, built from each cell's open faces.
  Trivially symmetric, so a convergence failure in this change is attributable
  to CG itself rather than to the preconditioner.
- **Fluid-only inner products and null-space handling.** Dot products skip
  zero-volume cells; where the boundary configuration makes the operator
  singular, the right-hand side, the initial guess and the preconditioned
  residual are projected against the fluid indicator vector.
- **BREAKING (numerically): `removeMean` becomes fluid-only.** It currently sums
  over every interior cell, divides by the full cell count rather than the fluid
  cell count, and writes the shifted value into solid cells. Fixing it changes
  the pressure field of singular setups that contain an obstacle. Non-singular
  setups are unaffected, since the projection is already skipped there.
- **Stopping criterion identical to the other solvers.** CG's `r·r` over fluid
  cells divided by the global fluid cell count is exactly what
  `pressureResidualNorm` returns, so `residual.dat`, `tools/fieldcmp` and the
  shell gates keep working. `r·r` and `r·z` are reduced together in one
  all-reduce.
- **Profiler regions for the Krylov work** (`CG_DOT`, `CG_AXPY`, `PRECON`)
  alongside `SWEEP_BULK` and `SWEEP_SURFACE`, so the allreduce cost CG exposes
  and the relaxation solvers do not is measured rather than lumped in.
- **The cross-rank iteration-count requirement is relaxed for Krylov solvers**
  to equality within one iteration. CG's step lengths come from global sums,
  whose floating-point reduction order depends on the decomposition, so a run
  near the tolerance can legitimately take one more or one fewer iteration. The
  converged field must still agree to the solve tolerance, unchanged.
- **The check scripts are repaired and extended.** Six of them still invoke
  `./NusifSolver-$TOOLCHAIN`, which the binary rename replaced with
  `CFD-Solver-$(TOOLCHAIN)`; they fail today. They are fixed and `cg` is added to
  their solver loops, along with new check drivers for operator symmetry and for
  CG itself. `.gitignore` carries the same stale name, so the test binary shows
  up as untracked; it is fixed with them.

- **The build learns which solver a build directory holds.** `make` links one
  `solver-$(SOLVER).o` and compiles everything with `-DSOLVER_$(SOLVER)`, but
  none of that is a timestamp, so switching back to a variant whose object file
  is already current left the previously linked binaries in place.
  `tests/run-all.sh` loops over every solver, so its second run reported one
  variant's results under another's name — a defect this change would otherwise
  inherit while adding a fourth solver to that loop, and which its own
  acceptance run cannot be trusted without. A stamp file recording the selection
  makes the choice an ordinary prerequisite.

### Non-goals

- **The multigrid preconditioner.** It needs a symmetric V-cycle, which the
  current cycle is not: the smoother always sweeps colours in one order, the
  coarsest solve is one-directional, and restriction is not a scalar multiple of
  the transpose of prolongation. That is a separate change, which this one is
  shaped to accept — the preconditioner is applied through a seam, not inlined.
- Pipelined or communication-reducing CG variants.
- Any change to the operator itself, the geometry representation, or the
  relaxation solvers' sweeps.

## Capabilities

### New Capabilities

None. CG is another solver for a system this project already specifies.

### Modified Capabilities

- `pressure-solvers`: adds requirements for a Krylov solver — that the operator
  is applied matrix-free with the same geometry-free bulk split the relaxation
  sweeps use, that inner products and the null-space projection are restricted
  to fluid unknowns, and that solid cells hold exactly zero throughout the
  iteration. Amends the existing cross-decomposition requirement so a Krylov
  solver's iteration count is equal within one rather than exactly equal.

Note: `openspec/specs/` is empty because `embedded-boundary-geometry` is
complete but not yet archived. `pressure-solvers` is that change's delta path
and is reused here rather than duplicated under a new name.

## Impact

- **New**: `src/solver-cg.c`, `tests/checks/cg.c`, `tests/checks/operator.c`.
- **Modified**: `src/solverbase.c` and `src/solver.h` (operator application and
  its surface-list sibling, the fluid-only reductions, the diagonal
  preconditioner); `src/comm.c`/`.h` (a counted all-reduce); `src/discretization.c`
  (`removeMean`); `src/profiler.h` and `src/profiler.c` (regions);
  `src/parameter.c`/`.h` (a preconditioner selection); `config.mk`,
  `mk/config-default.mk` and `README.md` (the `cg` option); `Makefile` and
  `.gitignore` (the solver-selection stamp and the stale binary name);
  `tests/run-all.sh`, `tests/check-solver-obstacle.sh` and the five other
  scripts carrying the stale binary name.
- **Unmodified**: the operator's coefficients, `surface-list.c`,
  `pressure-bc.c`, `geometry*.c`, and all three existing solvers' sweeps. CG
  must reproduce their converged fields, which is the primary acceptance gate.
- **No new dependencies**, but one new communication entry point. `commReduce`
  accepts a count and looked like the fused reduction's home, but it is a
  reduction to the master — it exists for the profiler's reporting, which is the
  one place a result nobody else sees is enough. CG's `alpha` and `beta` are
  needed on every rank, so `comm.c` gains `commReduceAllN(v, count, op)`: the
  in-place all-reduce `commReduceAll` already is, over several values at once.
  No new library and no new MPI concept.
