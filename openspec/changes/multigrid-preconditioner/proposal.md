# Proposal

## Why

`matrix-free-cg-solver` shipped a Krylov framework and said plainly what it was
missing: Jacobi preconditioning does not fix the Poisson condition number, so
CG's iteration count still grows with the grid and the solver is not yet the
fastest thing the benchmark has. A multigrid preconditioner is the step that
makes it one, and that change was shaped to accept it — the preconditioner
arrives through a function-pointer seam, not inlined.

It was deferred for a concrete reason: the V-cycle is not symmetric, and a
Krylov method preconditioned by an asymmetric operator is not conjugate
gradients at all. Three things make it asymmetric. The smoother sweeps its two
colours in the same order before and after the coarse correction; the coarsest
level is relaxed one-directionally; and restriction is an eight-cell average
while prolongation is trilinear, which are not a transpose pair.

Blocking all of that is a defect the same change recorded but did not fix.
`tests/checks/obstacle.c` fails under `SOLVER=mg`, at one rank and at four, with
576 solid cells holding a nonzero pressure after a converged solve. The cause is
`correct()`: it adds the prolongated error to every interior cell without
consulting the geometry. The post-smoother's surface-list pass repairs the solid
cells it knows about, but a cell in the body's deep interior is deliberately not
on that list, so nothing takes the correction back off. That already violates a
requirement the project has specified since `embedded-boundary-geometry` — solid
cells hold exactly zero after a converged solve — and it makes `tests/run-all.sh`
exit non-zero on a clean tree.

The two are not merely adjacent. A Krylov solver requires solid cells to hold
exactly zero at *every* iteration, not only at convergence. A V-cycle that writes
into the body would break that invariant the first time it was applied, so the
defect has to be fixed before the preconditioner can exist at all.

## What Changes

- **The coarse-grid correction respects the geometry.** `correct()` leaves cells
  with zero volume fraction alone, so a correction prolongated over the body
  never reaches it. This fixes `tests/checks/obstacle.c` under `SOLVER=mg` and
  makes `tests/run-all.sh` green on a clean tree for the first time since the
  obstacle work landed.
- **The multigrid hierarchy and cycle move out of `solver-mg.c`** into a unit
  every build links, the way `solverbase.c` already is. The Makefile links
  exactly one `solver-*.c`, so today a `SOLVER=cg` binary contains no multigrid
  code whatsoever and could not call a V-cycle if it wanted to. `solver-mg.c`
  becomes the thin driver that iterates cycles to a tolerance.
- **The V-cycle becomes symmetric**, in all three places it is not:
  - the post-smoother sweeps its colours in the reverse order of the
    pre-smoother, making the smoother pair a transpose pair;
  - the coarsest level is relaxed symmetrically rather than one-directionally;
  - restriction becomes a scalar multiple of prolongation's transpose.
- **BREAKING (numerically): `SOLVER=mg` adopts the symmetric cycle.** There is
  one cycle, not two. Multigrid's cycle counts and its converged fields move, and
  its baselines are re-recorded deliberately. The alternative — a second cycle
  used only for preconditioning — was rejected: two cycles is two things to keep
  correct, and the one that is not exercised by `SOLVER=mg` is the one that would
  rot.
- **Restriction is rebuilt as the transpose of trilinear prolongation**: the same
  three one-dimensional passes, reversed, followed by the eight-child sum.
  Deliberately not a 27-point tensor-product form — the existing prolongation is
  written as three separable passes precisely so that each reads only face
  neighbours, which is what the halo exchange provides, and its transpose keeps
  that property. No new communication.
- **`precon mg` becomes a supported value**, applying one V-cycle from a zero
  initial guess with fixed sweep counts. That satisfies the fixed-linear-operator
  requirement the spec already states: no inner tolerance, no warm start, the
  same work every application.
- **New check coverage** for the cycle as an operator: that it is symmetric, that
  it leaves solid cells at zero, and that restriction and prolongation are a
  transpose pair.

### Non-goals

- **Making multigrid a Galerkin method.** Forming the coarse operator as `R A P`
  would make symmetry structural rather than argued, but it becomes a 27-point
  stencil in three dimensions and changes the coarse kernels and their layout.
  The existing flux-balance coarsening stays.
- **Changing the smoother itself.** It stays red-black SOR; only the order in
  which the post-smoother visits the colours changes.
- **W-cycles, F-cycles, or full multigrid.**
- **Any change to the operator, the geometry representation, or the relaxation
  solvers' sweeps.**

## Capabilities

### New Capabilities

None. This change makes an existing solver correct and reuses it.

### Modified Capabilities

- `pressure-solvers`: adds requirements that the multilevel cycle is symmetric
  and usable as a preconditioner, and that a multigrid preconditioner is shipped
  and accelerates convergence. Amends the transfer-operator requirement so that
  restriction and prolongation must be a transpose pair, and the coarse-correction
  requirement so that the correction is confined to the fluid unknowns.
- `embedded-boundary-discretization`: amends the identity-row requirement so that
  it binds a multilevel correction explicitly, not only a relaxation sweep. The
  existing scenarios are about a converged solve and about a deep solid cell under
  the geometry-free sweep; neither reaches a correction prolongated from a coarse
  level, which is the gap this defect fell through.

## Impact

- **New**: `src/multigrid.c` and `src/multigrid.h` (the hierarchy and the cycle,
  extracted); `src/precon-mg.c` (the seam implementation); a symmetric-cycle
  section in `tests/checks/multigrid.c` and an `mg` case in `tests/checks/cg.c`.
- **Modified**: `src/solver-mg.c` (reduced to a driver); `src/solver-cg.c` (the
  `mg` preconditioner becomes selectable rather than refused); `src/solver.h`;
  `tests/checks/multigrid.c` (the restriction footprint assertion changes with
  the transfer pair); `tests/check-rejects.sh` (`precon mg` is no longer a
  refusal); `README.md`.
- **Baselines**: `tests/baseline/*-mg.dump` are re-recorded. The relaxation
  solvers and CG with `precon jacobi` are unaffected.
- **Binary size**: every build gains the multigrid object, including `rb` and
  `rbc`, which never call it. This is the price of `SOLVER=cg` being able to
  reach a V-cycle, and it is what `solverbase.c` already does for the operator.
- **No new dependencies.**
