# Design

## Context

See proposal.md — Why. What shapes the approach:

- `multigrid.c` already parameterises the parts that would differ. `smooth` takes
  a direction; `multigridBuild` takes a `MultigridSpecType`. A second cycle shape
  is a field on that spec and a branch in two places, not a second file.
- The hierarchy stops when a **local** extent cannot halve
  (`imaxLocal % 2 == 0 && imaxLocal / 2 >= 2`, reduced with `MIN` across ranks),
  so the depth available depends on the decomposition, not only on the grid. The
  solver already clamps to what it can build and prints
  `Multigrid: requested N levels, the decomposition supports M`.
- That clamp is what makes raising `levels` in a setup safe: a value the
  decomposition cannot honour is reported and reduced, never failed.
- The coarsest level today is `presmooth` forward sweeps then `postsmooth`
  backward ones. At three levels on `sphere-baseline` that is ten sweeps on a
  12x6x6 grid — a relaxation, not a solve.
- `check-solver-obstacle.sh` dumps each solver's field after a run at the setup's
  own `eps` and compares with `tools/fieldcmp` at that same `eps`.
- `tests/checks/multigrid.c` builds its own hierarchy from a
  `MultigridSpecType`, so it can drive either shape without touching a solver.

Measured before this change, on `sphere-baseline` at `eps = 1e-4`:

```
  mg, levels 3 (shipped)        29 cycles   0.36 s   6.0e-05 from rb   ok
  mg, levels 4                   8 cycles   0.17 s   5.1e-04 from rb   FAILS GATE
  mg, levels 4, eps 1e-5          -            -     6.5e-06 from rb   ok
  cg + precon mg, levels 3        5 iter    0.12 s
  cg + precon mg, levels 4        4 iter    0.12 s
```

## Goals / Non-Goals

**Goals:**

- A gate that distinguishes a wrong answer from an early stop, so that any
  later optimization can be judged at all.
- The depth the grids already support, without touching a grid.
- A coarse solve strong enough that depth is not compensating for it.
- Evidence on whether an asymmetric solver cycle is still worth having once the
  first three land.

**Non-Goals:**

- Changing the preconditioner's numerics. It keeps the symmetric shape.
- Reworking `tools/fieldcmp`. The gate decides what to compare and at what
  tolerance; the tool already takes both.

## Decisions

### 1. The gate converges past the tolerance it compares at

**Chosen:** the gate runs each solver at a tolerance an order of magnitude below
the setup's `eps`, and compares the resulting fields at `eps`.

A solver's `eps` bounds a mean square residual over the fluid cells.
`tools/fieldcmp` reports the largest pointwise difference. Nothing relates the
two: a solver that takes larger steps crosses the residual threshold from
further out and shows a larger pointwise difference while being no less
converged. Two configurations measured here do exactly that, and both agree to
6.5e-06 once converged.

This is not a loosening. The comparison tolerance is unchanged; what changes is
that both sides are converged before it is applied, which is what the
requirement always said ("each converged to a tight tolerance") and what the
script did not do.

Alternatives considered: comparing the L2 difference, which is the norm `eps`
actually bounds — rejected because it weakens the check against a single bad
point, which is the failure mode the gate most wants to catch. And giving the
comparison its own constant — rejected as an arbitrary number to justify per
setup, where converging further is a statement about the solvers.

The whole-run L2 comparison already in the gate for the rank-count check stays
as it is; it compares a solver against itself and has no cross-solver stopping
artifact.

### 2. Cycle shape is a field on the spec, not a second implementation

```c
typedef enum { MG_SHAPE_SYMMETRIC, MG_SHAPE_FAST } MgShapeType;
```

on `MultigridSpecType`, carried into `MultigridType`. It selects, in three
places:

| | symmetric | fast |
|---|---|---|
| post-smoother | reversed colours | forward |
| coarsest solve | k forward + k backward | 2k forward |
| restriction | `(1/8) P^T`, three blends | eight-cell average |
| equal sweep counts | required | not required |

Everything else — the hierarchy, the geometry coarsening, prolongation, the
residual, the correction masking, the three aliasing fixes — is shared and
unconditional. The two shapes differ in four branches, not in two code paths.

`solver-mg.c` asks for `MG_SHAPE_FAST`, `precon-mg.c` for `MG_SHAPE_SYMMETRIC`.
The equal-count refusal moves behind the symmetric branch, since only that shape
needs it.

Alternative considered: two cycle functions. Rejected — the shared part is the
overwhelming majority and is where every defect this project has found in
multigrid actually lived; duplicating it is how the copies drift.

Both shapes keep the eight-cell average path, which the symmetric shape no
longer uses. It is a handful of lines and it is what `fast` restores, so it is
kept rather than deleted and reintroduced.

### 3. The coarsest level gets a bounded, much stronger solve

**Chosen:** many more relaxation sweeps at the coarsest level, a fixed count,
tuned so the coarse residual falls by orders of magnitude.

The grid there is tiny — 6x3x3 at four levels on `sphere-baseline` — so a large
sweep count is cheap in absolute terms and is dwarfed by one pass over the
finest grid.

A fixed count keeps the preconditioner a fixed linear operator, which an inner
convergence test would break. For the symmetric shape the sweeps are split
evenly forward and backward, preserving the `A^T A` form that makes it
symmetric.

Alternatives considered: a direct solve at the coarsest level. It is the
textbook answer and would be exact, but it needs a factorization, a dense or
banded representation of the coarse operator, and a decision about which rank
owns it when the coarse level is still distributed. That is a larger change than
this one, and the measurements below will show whether it is needed — if a large
fixed sweep count already removes the coarse level as the limiting factor, a
direct solve buys nothing.

Alternative considered: a Krylov solve at the coarsest level. Rejected for the
preconditioner outright, since its iteration count would vary with the residual.

### 4. Setups get the depth their grids support, and no grid moves

`levels` is raised where the grid supports more, and nowhere else:

```
  dcavity.par        128x128x128   4 -> 7
  schaefer-turek.par 256x48x48     3 -> 5
  canal.par          200x40x40     3 -> 4
  dcavity-baseline   32x32x32      3 -> 5
  sphere-baseline    48x24x24      3 -> 4
  canal-baseline     50x10x10      2      already at its limit
  backstep.par       140x30x30     2      already at its limit
  karman.par         200x50x50     3 -> 2 already clamped at runtime
```

`karman.par` is the interesting one and is left alone: 50 halves to 25, so it
supports two levels whatever it asks for, and the solver already says so. Fixing
it means changing the grid to 200x48x48, which changes its resolution and the
flow it records. That is a setup decision, not a solver one.

**Revised while implementing: `schaefer-turek.par` keeps its three levels.**

The table above reasons from what each grid coarsens to, and that is only half
the question for a setup with a body. Measured on the shipped setup, the
cylinder is represented at every level at three, and stops being represented at
**level 3 of 4** — so every level this setup could gain is one whose coarse
correction is blind to the cylinder. `tests/check-setups.sh` already refuses a
shipped setup that reports `unresolved at level`, and it is right to: this is
the validation benchmark, its drag and lift are what it exists to produce, and
they are what a correction computed without the body would move.

The setup file had argued exactly this before the change and was right. Decision
3 is what makes leaving it cheap: the coarse solve is most of what the extra
depth would have bought, and it no longer needs depth to compensate for it.

So the depth requirement is read as "the depth its grid **and its body** support".
`karman.par` is left alone for its grid, `schaefer-turek.par` for its body, and
both say so in the setup file. The other five setups are unaffected — four have
no body at all, and `sphere-baseline` already reported an unresolved body at
three levels, so raising it to four changes nothing about that.

### 5. The asymmetric shape has to earn its place

Its advantage was measured at three levels, with a weak coarse solve, before
this change. 14 cycles against 29 is a 2.1x ratio in exactly the regime where
the smoother is doing most of the work. Deepening the hierarchy and strengthening
the coarse solve both reduce how much the smoother matters, so the gap should
narrow, and it may close.

**The comparison is run again at full depth before the split is kept.** If the
fast shape is not meaningfully ahead there, this change removes it and says so
rather than shipping a second shape for its own sake. The spec permits both but
requires neither, so dropping it costs nothing already promised.

This is a decision the measurements make, not one taken here. The task list
sequences it that way: depth and coarse solve first, the comparison next, the
split kept or dropped on the result.

## Risks / Trade-offs

- **Deeper hierarchies move the mg baselines again**, a third time in three
  changes. → Confined to mg, re-recorded with the change, and the cross-solver
  gate is what says the field is still right. Decision 1 makes that gate
  trustworthy first, which is why it lands first.

- **A deeper hierarchy on a coarse decomposition may build fewer levels than the
  setup now asks for**, so a run at high rank counts silently gets a shallower
  cycle than the same setup at one rank. → Already reported at initialization,
  and the depth requirement makes that reporting a requirement rather than a
  convenience. The gate's rank-count comparison is what catches a real
  divergence.

  **Realised, and the gate needed changing for it.** `sphere-baseline` at four
  levels builds 4 on 1 and 2 ranks and 3 on 4 and 8, reported every time, and
  takes 3 cycles against 5. The gate asserted that a non-Krylov solver's
  iteration count does not change with the rank count, which was true while
  every setup asked for a depth any decomposition could build and is not true
  now — so it failed a run whose fields agreed.

  Equality is the wrong requirement there: two different hierarchies are two
  different iterations, and a shallower one legitimately takes more cycles. The
  gate now compares counts only when both runs built the same depth, and when
  they did not, requires that the solver reported the clamp. The field
  comparison is untouched and is what says the two runs still solved the same
  system. Nothing is relaxed for `rb`, `rbc` or `cg`.

- **Two shapes is two things to keep correct**, which is the argument that was
  made against them when `multigrid-preconditioner` chose a single cycle. → The
  shared part is now large and the differing part is four branches, which is a
  much better ratio than it was then; and decision 5 keeps the split only if it
  is worth that cost.

- **A large fixed coarse sweep count is a tuned constant.** → Chosen from
  measurement and recorded, and bounded so the preconditioner stays a fixed
  linear operator. If it proves fragile, decision 3 names the direct solve as
  the next step.

- **The gate becomes slower**, since every solver in it runs to a tighter
  tolerance. → It is a correctness gate, not a benchmark, and the relaxation
  solvers dominate its runtime already.

## Migration Plan

Ordered so each step is judged by a gate that already works:

1. The gate converges before comparing. Nothing else changes; every solver still
   passes.
2. The coarsest solve is strengthened. Measured on its own, at the current depth.
3. Setups get their available depth. Baselines re-recorded.
4. The fast shape is added and compared at full depth, then kept or dropped.

Rollback is per step. Steps 1 to 3 are independent of the cycle split and are
worth keeping whatever step 4 concludes.

## Open Questions

- Whether a direct coarsest solve is needed, or whether a large fixed sweep
  count is enough. Deferrable: decision 3 ships the sweep count, and the
  measurements in step 2 are what would justify the larger change.
- Whether `karman.par` should move to a grid that can be coarsened. It is a
  setup and resolution decision rather than a solver one, and it changes what
  that setup records.
