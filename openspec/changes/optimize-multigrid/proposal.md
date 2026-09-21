# Proposal

## Why

`multigrid-preconditioner` made the V-cycle correct and symmetric, and bought a
preconditioner worth 24x on CG. It also cost `SOLVER=mg` roughly twice the
cycles it used to take, and left three things unexamined.

Measuring them says the symmetry is not what is holding multigrid back.

- **The shipped setups run shallow hierarchies.** `dcavity.par` asks for 4
  levels on a grid that supports 7; `schaefer-turek.par` asks for 3 of 5;
  `sphere-baseline` asks for 3 of 4. On `sphere-baseline` the fourth level takes
  multigrid from 29 cycles to 8, and the wall clock from 0.36 s to 0.17 s. That
  is a larger win than anything symmetry costs, it needs no symmetry given up,
  and the preconditioner gets it too.
- **The coarsest level is not solved.** At three levels on `sphere-baseline` the
  coarsest grid is 12x6x6 and is "solved" by ten relaxation sweeps. That is the
  weakness the extra level is compensating for, and it gets worse the shallower
  the hierarchy is.
- **The acceptance gate rejects faster solvers that are correct.**
  `check-solver-obstacle.sh` compares the largest pointwise pressure difference
  against `eps`, but `eps` bounds a *mean square residual*. Those are different
  quantities and nothing relates them. Two configurations measured here converge
  demonstrably to the right answer and fail the gate anyway: the deeper
  hierarchy lands 5.1e-04 from the reference at `eps = 1e-04` and 6.5e-06 at
  `eps = 1e-05`, and a larger smoothing factor behaves the same way. Any
  multigrid speedup trips this, so it has to be fixed before the rest can land.

Separately, the symmetric cycle is a constraint the standalone solver does not
need. Only a preconditioner requires it. Two of the constraints cost nothing
(reversing the post-smoother is the same work; the symmetric coarsest solve is
the same sweep count), but two are real: restriction as the transpose of
prolongation costs roughly a fifth more work per cycle, and the reversed
smoother goes unstable above a relaxation factor of about 1.6 where the forward
one runs happily at 1.8. At three levels that was worth 14 cycles against 29.

## What Changes

- **The cross-solver gate converges before it compares.** Both solvers are run
  to a tolerance well below the setup's `eps`, and the comparison is then made
  at `eps`. The agreement being tested is between converged fields, not between
  two solvers' stopping points, which is what the requirement always meant.
- **The shipped setups get the hierarchies their grids support.** A `levels`
  value the decomposition cannot build is already clamped and reported, so this
  is safe by construction: the solver reports what it actually built.
- **The coarsest level is solved rather than smoothed.** Enough work that the
  coarse problem stops being the cycle's limiting factor, and cheap because the
  grid is tiny.
- **Two cycles, selected by who is asking.** `SOLVER=mg` runs a cycle free to be
  asymmetric: forward smoothing on both sides, the cheap eight-cell averaging
  restriction, and no equal-sweep-count constraint. The preconditioner keeps the
  symmetric cycle unchanged. They share every part that does not differ.
- **The asymmetric cycle earns its place or is dropped.** Its measured advantage
  was taken at three levels, before the coarsest solve was strengthened. The
  same comparison is run again at full depth first; if the gap has closed, the
  split is not worth two cycles to maintain and this change says so and removes
  it rather than shipping a second cycle for its own sake.

  **Outcome: dropped.** At full depth with the coarsest level solved, the two
  shapes take the same number of cycles, and the asymmetric one wins about 5% of
  wall clock -- the per-cycle cost of the transposed restriction, not
  convergence. It was built, measured and removed; the numbers are in
  `design.md` decision 5 and the implementation is in the history.

### Non-goals

- **Changing any setup's grid.** `karman.par` is 200x50x50 and supports only two
  levels because 50 halves to 25; the fix is a grid dimension, which would change
  its resolution and its recorded behaviour. `schaefer-turek.par` is a validation
  benchmark whose published coefficients depend on its grid. Both are left alone.
- **Agglomerating coarse levels onto fewer ranks.** The hierarchy stops when a
  *local* extent cannot halve, so depth is decomposition-dependent. Redistributing
  coarse levels would lift that, and is a much larger change.
- **W-cycles, F-cycles and full multigrid.** Cycle *shape* is a separate axis from
  cycle *depth*, and depth is the measured win.
- **A Galerkin coarse operator.** Already rejected in `multigrid-preconditioner`
  for the same reason: 27-point in three dimensions.
- **Any change to the preconditioner's numerics.** It keeps the symmetric cycle
  it has.

## Capabilities

### New Capabilities

None. This makes an existing solver faster and an existing gate honest.

### Modified Capabilities

- `pressure-solvers`: amends the symmetry requirement so that it binds the cycle
  used as a preconditioner rather than every multilevel cycle the project ships,
  since a solver has no need of it. Adds requirements that the coarsest level is
  solved rather than approximated, and that a solver uses the hierarchy depth its
  grid and decomposition allow. Amends the cross-solver agreement requirement so
  that agreement is judged between converged fields rather than at whatever point
  each solver happened to stop.

## Impact

- **Modified**: `src/multigrid.c` and `src/multigrid.h` (a cycle-shape selection,
  the coarsest solve); `src/solver-mg.c` and `src/solver-cg.c` (which shape each
  asks for); `tests/check-solver-obstacle.sh` (converge before comparing);
  `tests/checks/multigrid.c` (the symmetry check becomes specific to the
  preconditioner's shape, and the asymmetric shape gets its own convergence
  coverage); `testcases/flow/` and `testcases/regression/` (`levels` in five
  setups only, never a grid).
- **Baselines**: `tests/baseline/*-mg.dump` are re-recorded. `rb`, `rbc` and `cg`
  are unaffected.
- **Sequencing**: this lands after `move-testcases`, which relocates every setup
  file. The paths above are the post-move ones; against the old layout they
  would not be found.
- **No new dependencies.**
