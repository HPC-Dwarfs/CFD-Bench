# Design

## Context

See proposal.md — Why. What shapes the approach:

- The hierarchy stops when a **local** extent fails
  `imaxLocal % 2 == 0 && imaxLocal / 2 >= 2`, reduced with `MIN` across ranks.
  For an extent `n = 2^k * m` with `m` odd, split across `2^p` ranks in that
  direction, the depth is `k - p + 1`, bottoming out at `m`.
- `commPartition` calls `MPI_Dims_create(size, 3, dims)`, which factors the rank
  count as cubically as it can and knows nothing about the domain's aspect
  ratio. It returns the factors in non-increasing order and the largest is
  applied to `x`, which happens to suit elongated domains.
- `main.c` calls `computeTimestep` only when `tau > 0`. With `tau 0` the step is
  whatever `dt` says, so a fixed step count is already expressible.
- `cylinder-z:xc,yc,r` and `box:` analytic primitives exist and go through the
  same `isSolid` predicate the voxel path does, sampled at `SUBSAMPLES = 4` per
  direction and thresholded at 0.5.
- `testcases/bench/` exists and is empty but for a README describing what
  belongs there.

Measured depth against rank count, which is what the sizes below come from:

```
  tier                    Mcells   P=1   P=27  P=216  P=1728
  dcavity  48^3              0.1    5L     4L     3L      2L
  dcavity  96^3              0.9    6L     5L     4L      3L
  dcavity 192^3              7.1    7L     6L     5L      4L

  a 128-cube, for contrast:  2.1    7L     2L     1L      -
```

## Goals / Non-Goals

**Goals:**

- Sizes derived from a stated rule rather than picked, so a fourth tier or a
  fourth case can be added by applying it.
- Tiers of one case differing only in cells, so their run times can be compared.
- The sizing rule enforced by a check, since its violation is silent.

**Non-Goals:**

- Choosing what to do with the timings. No results format, no runner, no
  regression thresholds on performance.
- Making the existing setups benchmarkable. They stay as they are.

## Decisions

### 1. Size from the smallest extent and the target rank count

The rule, in one line: **the smallest extent must still be at least 16 cells
after the decomposition divides it**, because four levels need `2^3 * 2`.

```
  smallest global extent  >=  16 * (ranks per direction)

    P = 27   ->  d = 3   ->  n_min >= 48
    P = 216  ->  d = 6   ->  n_min >= 96
    P = 1728 ->  d = 12  ->  n_min >= 192
```

So the tiers are `n_min` of 48, 96 and 192, and the other extents follow the
domain's integer ratio. Four levels is the threshold because below it the
coarsest grid is large enough that the cycle is doing little coarse work — which
is precisely the weakness `optimize-multigrid` addresses on the validation
setups.

**Chosen sizes:**

```
  case      ratio    small           medium          large
  dcavity   1:1:1    48^3            96^3            192^3
  karman    4:1:1    192x48x48       384x96x96       768x192x192
  canal     8:1:1    384x48x48       768x96x96       1536x192x192

  cells     dcavity  0.1 M           0.9 M           7.1 M
            karman   0.4 M           3.5 M           28.3 M
            canal    0.9 M           7.1 M           56.6 M
```

Note that 192 is what 1728 ranks requires, not 384. An earlier sketch of this
change said 384³ for the large cavity; that reaches five levels at 1728 ranks
rather than four, at eight times the memory — 56.6 M cells against 7.1 M. Four
levels is the stated threshold, so 192 is the size the rule gives.

### 2. Extents of the form `3 * 2^k`, not `2^k`

The counter-intuitive part, and the reason to write the rule down.

`MPI_Dims_create` produces factors like 3, 6 and 12 for rank counts that are not
powers of two. A pure power of two divided by 6 is odd, and an odd local extent
has no hierarchy at all:

```
  128 / 6 = 21  ->  1 level        192 / 6 = 32  ->  6 levels
  128 / 12 = 10 ->  2 levels       192 / 12 = 16 ->  4 levels
```

48, 96 and 192 are `3 * 2^4`, `3 * 2^5` and `3 * 2^6`. They divide by 2, 3, 4,
6, 8 and 12, which covers the factorisations that actually occur.

Alternative considered: powers of two, and require power-of-two rank counts.
Rejected — it makes the benchmark silently wrong rather than slow when someone
runs it at 27 or 216 ranks, and a benchmark whose validity depends on an
unstated precondition is a trap.

### 3. Domains rounded to exact integer ratios

The tier grids must match the domain's aspect ratio or the cells stretch, and a
point-smoothed multigrid cycle degrades on stretched cells — which the benchmark
would then be measuring by accident.

`karman` is 30x8x8, a ratio of 3.75:1:1, which no `3 * 2^k` grid matches. So the
benchmark channel is **32x8x8**, exactly 4:1:1, and the cells are cubic at every
tier. Likewise `canal` becomes **32x4x4** rather than 30x4x4, which also removes
the 1.5x anisotropy `canal.par` carries today. `dcavity` is already cubic.

This is why these are separate setups rather than tiered copies. A benchmark
case is free to adjust its domain for a clean discretization; a validation case
is not, and none of them are touched.

The cylinder keeps its position and radius — `cylinder-z:5.0,4.0,1.0`, the same
body the voxel volume describes — so the flow is the one `karman` produces, in a
channel two units longer.

### 4. A fixed step count, via `tau 0`

`tau 0` disables the adaptive controller, so `dt` is taken literally and the run
executes `te / dt` steps regardless of the grid or the flow.

**Chosen:** the same `dt` and `te` across the three tiers of a case, so the tiers
differ only in cells. The count is set to amortise startup — geometry sampling,
hierarchy construction and the first solve are not representative — while
keeping the small tier quick enough to run often.

The exact count comes from measurement rather than from this document; the task
list picks it from the small tier's startup fraction. What is decided here is
that it is fixed and identical across tiers.

Alternative considered: keeping `tau 0.5` and a physical final time, as the flow
setups do. Rejected outright — the number of steps would then depend on the
velocity field and the cell size, so the large tier would run a different number
of steps from the small one and the two would not be comparable. That is the
specific defect this decision exists to avoid.

### 5. Analytic geometry

Both geometry paths evaluate the same `isSolid` predicate at the same `4^3`
subsample points and threshold at the same 0.5. For a shape with a closed form
the analytic path is therefore at least as faithful — strictly more so, since
the voxel path quantises to its own grid first — and it is invariant under
refinement, which a fixed voxel volume is not.

The voxel path's real advantage is geometry with no closed form, which a
cylinder is not. Its cost here is concrete: `MIN_VOXELS_PER_CELL` is 4, and the
large `karman` tier would need a volume of about 226 million voxels or the run
is refused before the first step.

`karman.par` and `backstep.par` keep their voxel volumes, so the loader and its
rejection paths stay covered by `check-rejects.sh` and `check-setups.sh`.

### 6. The depth check computes, then validates its own arithmetic

The sizing rule has to be enforced, but the rank counts the large tier targets
cannot be run in a check.

**Chosen:** compute the expected depth for each tier at each rank count it
targets, from the grid and `MPI_Dims_create`'s factorisation, and assert it meets
the tier's claim. Then, at the rank counts that *can* be run — 1, 8 and 27 —
compare the computed figure against the level count the solver reports, so the
arithmetic is checked against the implementation rather than merely restating it.

That second half is what stops the check from being a second, drifting copy of
the coarsening rule. If `multigrid.c`'s limit ever changes, the small-scale
comparison fails and the computed figures are known to be stale.

Alternative considered: running each tier at each target rank count. Not
possible — 1728 ranks is not available in a check, and the large tiers are too
slow regardless.

## Risks / Trade-offs

- **The depth arithmetic duplicates the solver's coarsening rule**, and could
  drift from it. → Decision 6 is the mitigation: the computed values are checked
  against what the solver actually reports at the rank counts that can be run,
  so drift is caught rather than assumed away.

- **`MPI_Dims_create`'s factorisation is implementation-defined.** A different
  MPI could factor 216 as something other than 6x6x6 and change the depths. →
  The check computes from the same balanced factorisation the tiers were sized
  against and compares against reality at small rank counts, which is where a
  differing implementation would first show.

- **The large tiers are big.** 56.6 M cells for `canal` is a cluster job, not
  something to run casually. → That is what the tiers are for; the small tier
  exists to be run often, and nothing runs the large one by default.

- **Three cases is narrow coverage.** No backward-facing step, no validation
  benchmark. → Deliberate, and the shapes that matter for a solver benchmark are
  covered: a cubic domain with no geometry, an elongated one, and one with an
  obstacle in the flow.

- **The benchmark cases diverge slightly from the flow cases they are named
  after**, in domain length. → They are separate files in a separate directory
  with a README saying what they are for; the risk is someone reading a
  benchmark number as a `karman` result, which naming and documentation address.

## Migration Plan

Purely additive — nothing existing changes, so there is no rollback beyond
deleting the new files.

1. The three small tiers, and the depth check that certifies them.
2. The step count, measured on the small tiers and applied to all three sizes.
3. The medium and large tiers, certified by the same check.
4. Documentation.

## Open Questions

- Whether the depth check belongs in `tests/run-all.sh`. It is fast and it
  guards a silent failure, which argues for including it; it also checks
  configuration rather than behaviour, which is unlike everything else there.
  Deferrable: the check is worth having either way, and where it is invoked
  changes no artifact.
