# Design

## Context

See proposal.md — Why. What matters for the approach:

- `solverbase.c` already owns everything the three solvers share: the
  `PressureLevelType` description of the operator on one grid,
  `pressureResidualNorm`, and the save/correct pair that implements the
  geometry-free-bulk plus surface-list split. It is deliberately not named
  `solver-*.c` so the Makefile always links it.
- The Makefile links exactly one `solver-$(SOLVER).o`. A new `src/solver-cg.c`
  needs **no Makefile change** to be built — the existing pattern rule picks it
  up. The Makefile is touched for an unrelated reason: nothing about the
  selection is visible to `make` as a timestamp, so switching back to an
  already-built variant did not relink. See decision 10.
- The operator is symmetric with `-A` positive definite on the fluid subspace,
  and both boundary conditions (`wall` = mirrored halo, `outflow` = odd
  reflection) are homogeneous linear halo maps, so neither breaks symmetry.
  `embedded-boundary-discretization` already specifies this.
- Solid cells are identity rows with no open faces, so the system is block
  diagonal: `[I_solid | A_fluid]`. There is no coupling to eliminate — but the
  two blocks have opposite definiteness, so CG has to be confined to the fluid
  block rather than run on the whole vector.
- The only operator application in the tree is `residualField`, `static` inside
  `solver-mg.c`, and it reads four geometry streams per cell.

## Goals / Non-Goals

**Goals:**

- One operator application, in `solverbase.c`, level-generic through
  `PressureLevelType`, usable by CG now and by multigrid later.
- A preconditioner seam that the multigrid V-cycle can be dropped into without
  restructuring `solver-cg.c`.
- CG's stopping test producing the byte-identical quantity the other solvers
  report, so the existing shell gates and `tools/fieldcmp` need no new concepts.

**Non-Goals:**

- Extracting multigrid from `solver-mg.c`. That happens in the follow-up change
  that makes the cycle symmetric; this change only leaves room for it.
- Replacing `residualField` with the new application. Doing that changes
  multigrid's numerics and belongs with the multigrid work, not here.
- Reducing the reduction count below two per iteration.

## Decisions

### 1. A separate `SOLVER=cg` variant, not a flag on `mg`

**Chosen:** a new `src/solver-cg.c` exporting `initSolver`/`solve`, exactly like
the other three.

The alternative — a `krylov` parameter that wraps `solver-mg.c`'s cycle — avoids
a new file but breaks two established conventions: check drivers discriminate on
`-DSOLVER_$(SOLVER)`, and the shell gates loop `for solver in rb rbc mg`. A real
variant slots into both. It also means the `cg` binary does not carry multigrid
code, which matters for a benchmark whose binaries are profiled.

### 2. Solve the negated system

The operator as assembled is negative semi-definite: `A p = sum_f A_f (p_nb -
p_c)/h^2` and the equation is `A p = Lambda * rhs`. CG requires positive
definiteness, so the solver works with

```
  b_c = -Lambda_c * rhs_c      (fluid),   0  (solid)
  (A_neg x)_c = -sum_f A_f (x_nb - x_c)/h^2   (fluid),   0  (solid)
```

Negating both sides leaves the solution and the residual `r = b - A_neg x`
unchanged in magnitude, so `r·r` is still exactly the numerator
`pressureResidualNorm` computes. Flipping the sign in the apply rather than
post-negating a residual field avoids a whole extra pass.

Alternative considered: keep the native sign and use CG on `-A` implicitly by
negating `alpha`. Rejected — it makes every sign in the recurrence a thing to
re-derive when reading the code.

### 3. `pressureApplyOperator` and `pressureApplySurface` in `solverbase.c`

```
  void pressureApplyOperator(const PressureLevelType *lv,
                             double *x, double *y);
```

Structure mirrors the relaxation sweep exactly:

```
  commExchange(lv->comm, x); pressureBcApply(lv->bc, ..., x);

  # bulk: no geometry streams, no branches, 7-point
  for all interior cells:
      y = -((xE - 2x + xW)*idx2 + (xN - 2x + xS)*idy2 + (xT - 2x + xB)*idz2)

  # surface: O(obstacle area), real coefficients from the existing list
  pressureApplySurface(lv, x, y);
```

`pressureApplySurface` is a **sibling** of `pressureCorrectSurface`, not a reuse
of it: that one relaxes (reads `saved`, writes a new `p`, accumulates a residual
delta), this one applies (reads `x`, overwrites `y`). It reuses the same
`SurfaceListType` — `aE..aB`, `lambda`, `solid` — and unlike the relaxation pass
it needs no colour split, because an apply has no ordering constraint. It walks
all `count` entries in one loop and writes `y = 0` for solid entries.

This is what keeps the `Interior sweep cost is independent of geometry`
requirement true for CG, which the existing spec already binds "in any of the
solvers the project ships".

Alternative considered: expose `residualField` from `solver-mg.c`. Rejected —
it reads `Ax/Ay/Az/Lambda` for every cell, which is four extra streams on CG's
hot path and forfeits the benchmark's central claim.

### 4. Jacobi preconditioner with no diagonal array

The diagonal is `sum_f A_f / h^2`. For an interior cell that is the constant
`2*(idx2 + idy2 + idz2)`; for a cut cell the surface list already carries
`invDiag`, and it is `0` for solid cells.

So the preconditioner needs **no stored field at all**: a scalar multiply over
the bulk, then an override pass over the surface list. Same split, same
zero-geometry-stream property as the apply. Solid cells get `z = 0` for free
from `invDiag == 0`, which is exactly the invariant the spec demands.

### 5. A preconditioner seam sized for the multigrid follow-up

```c
typedef struct {
  void (*apply)(void *ctx, const PressureLevelType *lv,
                const double *r, double *z);
  void *ctx;
} PreconType;
```

`solver-cg.c` calls `precon->apply(...)` and knows nothing else. The follow-up
change adds a `precon-mg.c` whose `ctx` is the level hierarchy and whose `apply`
is one V-cycle from a zero start; `solver-cg.c` is untouched. Selected by a
`precon` parameter (`none` | `jacobi`, with `mg` rejected until it exists),
parsed with the existing `PARSE_STRING` the way `geometryFile` is.

The seam also encodes the linearity requirement structurally: `apply` receives
`r` and writes `z` outright, so there is nowhere for an inner tolerance or a
warm start to hide.

### 6. Two fused reductions per iteration

```
  q      = A_neg d
  d·q                                      <- allreduce 1 (count 1)
  alpha  = rho / (d·q)
  x     += alpha d ;  r -= alpha q
  z      = M r  [+ null projection]
  {r·r, r·z}                               <- allreduce 2 (count 2)
```

The stopping norm costs nothing beyond the `r·z` the algorithm needs anyway,
provided the two travel in one call. `commReduce(v, o, count, op)` takes a count
and looks like the place for that, but it is `MPI_Reduce` to the master: it
exists for `printProfile`, where only rank 0 reads the answer, and on every
other rank the result is untouched. CG's `alpha` and `beta` are needed
everywhere, so the count belongs on the all-reduce instead and `comm.c` gains

```c
  void commReduceAllN(double *v, int count, int op);
```

which is what `commReduceAll` already does, over several values at once.
Extending `commReduceAll` itself was the alternative and was rejected: it is
called from a dozen places that all pass one value, and widening its signature
would churn every one of them to describe a case only CG has. Two allreduces per
iteration is standard PCG; the relaxation solvers have one, and only for their
stopping test. That difference is real and is the point — it gets its own
profiler region rather than being hidden.

All three dots skip `Lambda == 0` cells. Since solid `x`, `r`, `z`, `d`, `q` are
all held at exactly zero, skipping is a correctness statement rather than an
optimization, and the check for it is cheap to assert in the test driver.

### 7. Null space: project the fluid indicator vector out of `b`, `x0` and `z`

Where `pressureBcIsSingular` holds, the null vector is the constant over fluid
cells and zero over solid cells. Three projection points:

- `b`, once, so the system is consistent;
- `x0`, once, so the warm start from the previous time step carries no constant;
- `z`, every iteration, because a preconditioner is free to inject a constant
  and CG will then amplify it.

`removeMean` in `discretization.c` is fixed in the same breath: it currently
sums over every interior cell, divides by `imax*jmax*kmax` instead of
`fluidCells`, and writes the shifted value into solid cells — leaving them at
`-mean`. SOR self-heals that on the next cut-cell pass; CG does not. This is the
one change with an externally visible numerical effect, and it is confined to
singular setups that contain an obstacle.

### 8. CG state hangs off `Solver` as an opaque pointer

`Solver` already carries `void *mgLevels` with the comment that it is opaque
outside `solver-mg.c`. CG follows that precedent with `void *cgState` holding
`r`, `z`, `d`, `q` and the `PreconType`, rather than adding four more named
field pointers that three of the four solvers would never touch.

Memory cost: four extra finest-level fields, `4 * (imaxLocal+2)(jmaxLocal+2)
(kmaxLocal+2)` doubles. Stated in the build banner alongside the existing
"Using Multigrid solver with N levels" line.

### 9. Test seams, and one driver that needs none

- `tests/checks/operator.c` tests `pressureApplyOperator` — symmetry,
  definiteness, the seven-point footprint, agreement with
  `pressureResidualNorm`. It lives in `solverbase.c`, which every variant links,
  so this driver runs unmodified under all four and needs **no `#ifdef` seam**.
  It is worth having on its own: it is the first direct test of the operator's
  symmetry, which until now was an argued property rather than a measured one.
- `tests/checks/cg.c` needs `#if defined(TEST) && defined(SOLVER_cg)` to reach
  the preconditioner and the iteration state, matching the existing `mgTest*`
  pattern in `solver.h`.

### 10. Repair the acceptance gates before extending them

`check-setups.sh`, `record-baseline.sh`, `check-schaefer-turek.sh`,
`check-solver-obstacle.sh`, `check-rejects.sh` and `check-particles.sh` all
invoke `./NusifSolver-$TOOLCHAIN`; commit `ef77c8e` renamed the binary to
`CFD-Solver-$(TOOLCHAIN)`. They fail today. `check-solver-obstacle.sh` is this
change's primary acceptance gate, so the rename is fixed here rather than left
for someone else — in `.gitignore` too, which carries the same stale name and so
leaves the test binary showing up as untracked. The gate is also taught that
`cg` is a Krylov solver and gets the within-one iteration-count tolerance the
spec now grants it.

A second defect in the same path is not a script but the Makefile. Exactly one
`solver-$(SOLVER).o` is linked and everything is compiled with
`-DSOLVER_$(SOLVER)`, but neither is visible to `make` as a file timestamp, so
selecting a variant whose object file is already current relinks nothing and the
previously built binaries are run instead. `tests/run-all.sh` loops over every
solver and ends on the last one, so its **second** run reports one variant's
results under another's name. This change adds a fourth solver to that loop and
takes `run-all.sh` as its own before-and-after evidence, so the defect has to go
first or task 7.4 means nothing.

**Chosen:** write the selection to `$(BUILD_DIR)/solver.sel` whenever it
changes, and make that file a prerequisite of the target, the test target and
every check binary. The selection then behaves like any other dependency.

Alternatives considered: a `clean` between variants, rejected as it turns each
gate into a full rebuild of everything; and one build directory per solver,
rejected because `BUILD_DIR` is already overridden by two of the gates for their
own purposes and nesting the two conventions is worse than a stamp file.

## Risks / Trade-offs

- **CG stalls or breaks down where SOR merely converges slowly.** → The
  acceptance path is staged so the cause is never ambiguous:
  `tests/checks/operator.c` proves symmetry and definiteness before any CG runs,
  and the Jacobi preconditioner is trivially symmetric, so a failure in
  `tests/checks/cg.c` is attributable to CG's own implementation and nothing
  else. This is the main reason the multigrid preconditioner is a separate
  change.

- **The `removeMean` fix changes published baselines.** → It only affects
  singular setups containing an obstacle. `tests/record-baseline.sh` will flag
  whichever shipped setups move; those baselines are re-recorded deliberately
  and the move is explained in the commit, rather than the fix being silently
  absorbed.

- **Jacobi-CG may be slower in wall-clock than SOR on these grids.** Diagonal
  preconditioning does not fix the Poisson condition number, so the iteration
  count still grows with the grid. → This is expected and is not a regression:
  the change's value is a correct, verified Krylov framework with the operator
  and reductions in place. The speed argument arrives with the multigrid
  preconditioner. Saying so up front stops the first benchmark run from reading
  as a failure.

- **Two allreduces per iteration make CG latency-sensitive at high rank
  counts.** → Not mitigated, and not a defect — it is a genuine property of the
  method and one of the reasons a CFD benchmark wants a Krylov solver at all.
  Measured via its own profiler region so it is visible rather than inferred.

- **`pressureApplySurface` is a near-copy of `pressureCorrectSurface`, inviting
  drift between the two.** → They are placed adjacent in `solverbase.c` under
  the shared operator comment, and `tests/checks/operator.c` checks the apply
  against `pressureResidualNorm` — which uses the geometry-reading path — at
  cells adjacent to the obstacle. A divergence between the relaxation
  coefficients and the apply coefficients fails that check.

- **`PARSE_PARAM` matches parameter names by prefix.** A new `precon` key could
  in principle shadow another. → `presmooth` differs at the fourth character, so
  there is no clash today; the rejects gate gains a case that an unknown
  `precon` value aborts at initialization rather than falling back to `none`.

## Migration Plan

No data or format migration. Rollout is a build-time selection: `SOLVER=cg` is
additive and the default in `config.mk` stays `mg`. Rollback is reverting the
`solverbase.c`, `discretization.c` and script changes; `src/solver-cg.c` is
inert unless selected.

The one non-additive step is the `removeMean` fix, which affects all four
solvers. It ships with re-recorded baselines in the same commit so that
`tests/record-baseline.sh` is green at every point in history.
