# Design

## Context

See proposal.md — Why. What shapes the approach:

- The Makefile links every `src/*.c` except `vtkWriter-*.c` and `solver-*.c`, of
  which it links exactly one each. `solverbase.c` is deliberately not named
  `solver-base.c` for that reason. A file that every variant must link therefore
  needs only a name outside that pattern — no Makefile change.
- `solver-mg.c` today holds four things that are not the same thing: the level
  hierarchy and its geometry coarsening, the transfer operators, the smoother
  and cycle, and the driver that iterates cycles to a tolerance. Only the last
  is specific to `SOLVER=mg`.
- `pressureApplyOperator` and `pressureApplySurface` already live in
  `solverbase.c` and are shared. `residualField` inside `solver-mg.c` is a
  second, geometry-reading operator application; `matrix-free-cg-solver`
  explicitly left it alone, and this change does too.
- The `PreconType` seam takes `(ctx, level, r, z)` and nothing else. A V-cycle
  fits it as long as it runs a fixed number of cycles from a zero guess.
- `presmooth` and `postsmooth` are independent parameters. Two shipped setups,
  `canal.par` and `dcavity.par`, set them to 50 and 5.

## Goals / Non-Goals

**Goals:**

- One cycle, symmetric, used by `SOLVER=mg` and as CG's preconditioner.
- The extraction leaves `solver-mg.c` a driver and introduces no behaviour of
  its own, so that the numerical change is attributable to the symmetry work and
  not to the move.
- Symmetry that is measured rather than argued, before any preconditioned solve
  runs.

**Non-Goals:**

- Replacing `residualField` with the matrix-free application. It would change
  multigrid's numerics in the same change that changes them for other reasons,
  making a regression impossible to attribute.
- Making the hierarchy usable from a build that did not ask for it. `rb` and
  `rbc` link the object and never construct a hierarchy.

## Decisions

### 1. Fix the correction before anything else

`correct()` adds `e` to `p` at every interior cell. The fix is one geometry test
in that loop: a cell with zero volume fraction keeps its zero.

**Chosen:** mask in `correct()`, reading the level's own `Lambda`.

Alternative considered: zero `e` inside `prolongate()` instead. Rejected, though
it is nearly equivalent — `prolongate` is a transfer operator and is about to
become half of a transpose pair, so keeping it a pure interpolation with no
geometry in it is what makes the transpose argument checkable. The geometry
belongs where the correction is applied.

Alternative considered: extend the surface list to cover deep solid cells.
Rejected outright — the list is `O(body surface)` by construction, and that is
the property the whole benchmark rests on.

This lands first, as its own step, because it is a bug fix against a requirement
the project already states and it can be verified on its own:
`tests/checks/obstacle.c` under `SOLVER=mg` goes from failing to passing with no
other change in the tree.

### 2. Extract into `src/multigrid.c`, leaving a driver behind

```
  multigrid.c   hierarchy, geometry coarsening, transfers, smoother, cycle
  solver-mg.c   initSolver builds the hierarchy; solve iterates cycles
  precon-mg.c   the PreconType whose apply is one cycle from zero
```

The cycle stops taking `Solver *`. It takes a context carrying only what it
uses — the boundary condition, `omega`, the level count, the smoothing counts
and the level array — so that a preconditioner can build one without pretending
to be a solver. `Solver` keeps an opaque pointer, as it does today.

The extraction is a pure move: same code, same numerics, verified by
`tests/checks/multigrid.c` and the mg baselines being untouched across it. Doing
it in the same step as the symmetry work would make any regression
unattributable.

Trade-off accepted: `rb` and `rbc` binaries gain the multigrid object without
calling it. This is what `solverbase.c` already does, and the alternative — a
conditional compile — buys a smaller binary at the cost of a build matrix.

### 3. Symmetry, in the three places it is missing

A V-cycle is symmetric when the post-smoother is the transpose of the
pre-smoother, the coarsest solve is symmetric, and `R = c P^T`.

**Smoother.** Red-black SOR sweeping red then black has as its transpose the
same sweep black then red. So `smooth` gains a direction, and the post-smoother
runs with the colours reversed. Nothing about the sweep itself changes.

**Equal counts.** `S_post = S_pre^T` also requires the same number of sweeps on
each side: a product of two operators is not the transpose of a product of
three. `presmooth` and `postsmooth` must therefore be equal, which two shipped
setups violate.

**Chosen:** refuse unequal counts at initialization, naming both values, the way
an unsupported `precon` and an unsupported boundary code are refused. `canal.par`
and `dcavity.par` move to 5 and 5, matching every other shipped setup and both of
their own baseline variants.

Alternative considered: silently use `min(presmooth, postsmooth)`. Rejected —
it produces a cycle that is symmetric but is not the one the setup asked for,
and the setup would never say so.

Alternative considered: let the counts differ and accept asymmetry for
`SOLVER=mg`, requiring equality only when preconditioning. Rejected, because
"one cycle" is the whole point of the choice made in this change; a cycle whose
symmetry depends on which caller invoked it is two cycles wearing one name.

**The relaxation factor.** Found during implementation, and not anticipated
above: the factor the smoother uses is the setup's `omg`, which is 1.7 — the
optimum for SOR *as a solver*. As a smoother that is a poor choice, because with
a factor above 1 the sweep amplifies some high-frequency modes instead of
damping them, and which modes survive depends on the colour order. The forward
sweep happens to survive it; its transpose does not. Reversing the post-smoother
with `omg = 1.7` drives `sphere-baseline` to NaN, and takes the Poisson check
from 1.25e-01 to 3.08e-01 in a single cycle. Forcing the smoothing factor to 1.0
makes the reversed cycle pass every existing check.

**Chosen:** a `smoothOmega` parameter, defaulting to 1.3, which the cycle reads
in place of `omg`. `omg` stays what it has always been, the relaxation factor of
the SOR solvers. The two were never the same quantity; they only looked like one
because a single cycle direction hid the difference.

Alternatives considered: hardcoding 1.0, rejected because the smoothing factor
is worth tuning per setup and a constant cannot be; and clamping `omg` to 1.0,
rejected because a setup asking for 1.7 would silently get something else.

Measured on `sphere-baseline`, cycles for the last solve: 41 at 1.0, 29 at 1.3,
21 at 1.6, and divergence at 1.8. Symmetry itself does not depend on the value —
it was verified at 1.0 and at 1.6 alike — so the factor trades convergence
against stability and nothing else. Stability is not the only bound. Compared against `rb` on the same setup, the
converged pressure differs by 9.6e-05 at 1.0, 6.0e-05 at 1.3 and 1.12e-04 at
1.6, against a gate of 1e-04 — all three converge to the same solution, as a
tighter `eps` confirms to 4.4e-06, but 1.6 stops marginally outside the
agreement the project requires. 1.3 is chosen as the fastest value that both
converges and passes that gate. It is a parameter rather than a constant because
the usable range is narrow and setup-dependent.

This moves multigrid's numerics for the pre-smoother as well as the post-,
which is a larger baseline move than the symmetry work alone implied.

**Coarsest level.** Today it is `presmooth + postsmooth` forward sweeps.
Symmetric form: `k` forward sweeps then `k` backward sweeps, which as an error
propagation operator is `A^T A` and therefore symmetric for any `k`.

### 4. Prolongation has to become trilinear before it can be transposed

Found during implementation. `prolongate` is documented as trilinear and is not:
each of its three one-dimensional passes updates in place and reads a face
neighbour that, for every other cell, the same loop has already overwritten.

```
  nb = (i & 1) ? idx - 1 : idx + 1;
  fineField[idx] = 0.75 * fineField[idx] + 0.25 * fineField[nb];
```

For odd `i` the neighbour is `idx - 1`, written one iteration earlier, so those
cells blend against a partly interpolated value. On a coarse ramp of 10, 20, 30
the pass produces 18.125 where trilinear interpolation gives 17.5 — effective
weights of 0.8125 and 0.1875. It predates this change and is invisible to the
existing checks, both of which restrict and prolongate a **constant**, which any
pair of weights summing to one preserves.

It blocks the decision below: an operator whose matrix form depends on traversal
order has no transpose of the shape that decision assumes.

**Chosen:** make each pass a true blend by updating the two cells of a pair
together from their own values, and only then take the transpose. The pairing is
already there in the index arithmetic — `nb` is an involution, pairing (2,3),
(4,5) and so on, with cell 1 and the last cell blending against a halo no pass
writes. Computing a pair from two temporaries removes the aliasing and, as a
side effect, makes each pass symmetric, which is what makes the transpose in
decision 5 tractable.

Alternative considered: transpose the sweep as it actually is. Its transpose
exists — it is the reversed traversal — but it would enshrine the wrong weights
and leave both operators defined by a loop order rather than by a stencil.

This moves multigrid's numerics a third time. Better interpolation should cost
fewer cycles, partly offsetting the smoothing-factor change, but the direction
is not assumed here; it is measured.

### 5. Restriction becomes the transpose of prolongation

Prolongation is injection followed by three one-dimensional 3/4–1/4 passes,
written that way so each pass reads only face neighbours — the tensor-product
form would read diagonal neighbours the halo exchange does not carry, making the
result depend on the decomposition.

Its transpose inherits that structure exactly, reversed:

```
  P  =  pass_z . pass_y . pass_x . inject
  P^T =  inject^T . pass_x^T . pass_y^T . pass_z^T
```

`inject^T` is the sum over a coarse cell's eight children, which is what the
current restriction already computes before scaling. Each `pass^T` is again a
one-dimensional face-neighbour stencil. So restriction stays separable, stays
face-only, and needs no new communication — it gains three passes and their
exchanges over the current single averaging pass.

**Scale factor.** `c = 1/8`, chosen so that restriction still maps a constant to
that constant: the column sums of `P` are 8, so `(1/8) P^T` applied to a
constant returns it. That keeps the existing `restriction of a constant is that
constant` check meaningful rather than needing to be rewritten around a new
normalization.

**What this breaks in the tests.** The current driver asserts that a unit fine
impulse restricts into exactly one coarse cell. The transpose has a wider
footprint by construction, so that assertion is replaced with the transpose-pair
identity itself, which is the stronger statement and the one the spec now names.

Alternative considered: drop prolongation to injection and keep the cheap
average, an exact transpose pair for free. Rejected — injection is a poor
interpolation and would cost convergence for `SOLVER=mg`, which is not a trade
this change is entitled to make on mg's behalf.

Alternative considered: a Galerkin coarse operator, `R A P`, making symmetry
structural. Rejected as out of scope: 27-point in three dimensions, changing the
coarse kernels and their storage.

### 6. The preconditioner is one cycle from zero

```c
  static void preconMg(void *ctx, const PressureLevelType *lv,
                       const double *r, double *z)
  {
      zero(z);                 /* a zero initial guess, every time */
      multigridCycle(ctx, z, r);
  }
```

Fixed sweep counts, one cycle, no residual test inside. That is what the
existing `A preconditioner is a fixed linear operator` requirement demands, and
the seam's shape already prevents an inner tolerance from being smuggled in.

Solid cells come out at zero because decision 1 confines the correction and the
smoother maintains the body at zero from a zero start.

`solver-cg.c` builds the hierarchy when `precon mg` is selected, holds it in the
preconditioner's `ctx`, and is otherwise untouched — which is what the seam was
built for. `precon mg` stops being refused.

**Levels.** The hierarchy uses the same `levels` parameter multigrid already
reads, and reports when the decomposition supports fewer, as it does today.

### 7. The preconditioner negates the cycle

Found during implementation. The cycle relaxes the operator as assembled, which
is negative semi-definite, so one cycle applied to `r` approximates `A^-1 r`.
The conjugate gradient solver works with the negated system and needs an
approximate inverse of `-A`. Without the sign the preconditioner is linear and
symmetric -- it passes both of those checks -- and negative definite, which is
exactly the case a Krylov method cannot use and which the definiteness check in
`tests/checks/cg.c` caught.

**Chosen:** negate the cycle's result in `precon-mg.c`. It is the sign
convention made explicit in the one place that bridges the two, rather than a
correction applied somewhere along the way.

### 8. Default stays `jacobi`

`precon` keeps its current default. Making `mg` the default is a benchmarking
decision that should follow the measurements this change produces, not precede
them.

## Risks / Trade-offs

- **Each cycle costs more.** The reversed post-smoother is the same work, and the
  symmetric coarsest solve is the same sweep count, but restriction gains three
  passes and three halo exchanges. → Cycle counts should fall enough to pay for
  it; if they do not, that is a finding worth having and the measurement is part
  of the change. Stated up front so the first timing run does not read as a
  regression.

- **`SOLVER=mg`'s numbers move and its baselines are re-recorded.** → Confined to
  mg. The relaxation solvers and `precon jacobi` are untouched, and the
  cross-solver agreement gate is what says the new cycle still solves the same
  system. The move ships with the re-recorded baselines in the same commit.

- **Two shipped setups change their smoothing counts.** → `canal.par` and
  `dcavity.par` go from 50/5 to 5/5. Both have baseline variants already at 5/5,
  so the shipped setup and its baseline stop disagreeing — which is arguably a
  fix in its own right.

- **Symmetry is exact only in exact arithmetic.** The check compares two inner
  products that are formed by different routes through the same operator. →
  Judged relative to their magnitude, the way the operator and preconditioner
  symmetry checks already are, not against an absolute epsilon.

- **mg-preconditioned CG may not beat `SOLVER=mg` in wall-clock**, since it pays
  a cycle plus two allreduces per iteration where mg pays a cycle. → Not a
  failure: the value is a Krylov method whose convergence no longer degrades with
  refinement, which is what the spec now requires and what the relaxation and
  multilevel solvers cannot offer. Both numbers get measured and reported.

- **The extraction is a large diff that touches a working solver.** → It is a
  separate step with its own gate: the mg baselines and
  `tests/checks/multigrid.c` must be bit-identical across it, before any
  symmetry change lands.

## Migration Plan

Four steps, each independently verifiable, in order:

1. The `correct()` fix — `tests/checks/obstacle.c` under `mg` passes, nothing
   else moves.
2. The extraction — a pure move, mg baselines bit-identical.
3. The symmetric cycle — mg baselines re-recorded, symmetry measured.
4. `precon mg` — additive; `precon jacobi` and every other solver untouched.

Rollback is per step. Step 1 is worth keeping regardless of what happens to the
rest, since it fixes a requirement violation that predates this change.

## Open Questions

- Whether the coarsest-level sweep count should stay `presmooth + postsmooth`
  split evenly, or become its own parameter. Deferrable: it changes a constant,
  not the specs, the approach or the task breakdown, and the measurements in step
  3 are what should decide it.
