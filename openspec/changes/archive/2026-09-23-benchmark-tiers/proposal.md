# Proposal

## Why

The project has no setup that can be run at scale, and the ones it has are the
wrong shape for it.

The shipped flow cases are validation-shaped: a physical final time, an adaptive
time step, and in two cases a recorded baseline. Running one at 512 ranks tells
you little, because the amount of work depends on how many steps the adaptive
controller chooses, which depends on the grid. Two runs at different sizes are
not comparable.

Worse, their grids are hostile to multigrid once the domain is divided. The
hierarchy stops when a **local** extent can no longer be halved, so depth is a
property of the grid *and* the decomposition. `karman.par` is 200x50x50 and
gets two levels before any decomposition at all, because 50 halves to 25.
`MPI_Dims_create` then factors the rank count without knowing the aspect ratio,
so at 216 ranks it produces 6x6x6 and a grid of 128 cells per direction becomes
21 — odd, and the hierarchy collapses to a single level. A benchmark that does
that is not measuring multigrid; it is measuring a smoother with extra steps,
and nothing in the output says so.

That is the trap worth naming: **pure powers of two are the worst choice** for a
grid that will meet an arbitrary rank count. Extents of the form `3 * 2^k` — 24,
48, 96, 192 — divide cleanly by 2, 3, 4, 6, 8 and 12, so they survive the
factorisations `MPI_Dims_create` actually produces. Measured against rank counts
from 1 to 1728, a 192-cube keeps 4 levels where a 128-cube drops to 1.

## What Changes

- **Nine benchmark setups in `testcases/bench/`**, three tiers each for three
  cases, sized so the multigrid hierarchy survives to the rank count each tier
  targets:

  | case | domain | small | medium | large |
  |---|---|---|---|---|
  | `dcavity` | 1:1:1 | 48³ | 96³ | 192³ |
  | `karman` | 4:1:1 | 192x48x48 | 384x96x96 | 768x192x192 |
  | `canal` | 8:1:1 | 384x48x48 | 768x96x96 | 1536x192x192 |

  Smallest extent 48, 96 and 192, which holds four levels or more at roughly 27,
  216 and 1728 ranks respectively.

- **A fixed amount of work per run.** Each tier sets `tau 0`, which disables the
  adaptive controller, with a fixed `dt` and a final time that is an exact
  multiple of it. Every tier of a case runs the same number of steps, so the
  only thing that varies with the tier is the cells.

- **Analytic geometry.** The `karman` tiers use `cylinder-z:` rather than a voxel
  volume. Both paths funnel through the same point predicate and the same 4³
  subsampling, so a voxel volume adds a quantisation step and no fidelity; and
  at `MIN_VOXELS_PER_CELL = 4` a large tier would need a volume of roughly 226
  million voxels, regenerated per tier, or the run is refused at startup.

- **Domains rounded to exact integer ratios** so cells are cubic: the `karman`
  tiers use a 32x8x8 channel rather than karman's 30x8x8, and the `canal` tiers
  32x4x4 rather than 30x4x4. Anisotropic cells degrade a point-smoothed
  multigrid cycle, and a benchmark should not measure that by accident. This
  also removes the 1.5x anisotropy `canal.par` carries today.

- **A check that the tiers are actually benchmarkable**, asserting the hierarchy
  depth each one reaches at the rank counts it targets. Without it the sizing
  rule is a comment, and the failure it guards against is silent.

### Non-goals

- **No baselines.** These exist to be timed, not compared against a stored
  field. Numerical agreement is the regression setups' job and stays there.
- **No change to any existing setup.** The flow and regression setups keep their
  grids, their tolerances and their recorded fields. This change only adds.
- **`backstep` and `schaefer-turek` are not tiered.** The first has a 4.67:1:1
  domain that rounds to nothing convenient; the second is 250:41 with 41 prime,
  so it cannot be made both isotropic and rank-friendly, and its published
  coefficients depend on its discretization. Three cases — a cubic cavity, an
  elongated channel, and one with an obstacle — cover the shapes that matter.
- **No coarse-level agglomeration.** Redistributing coarse levels onto fewer
  ranks would lift the depth limit entirely and is the right long-term answer;
  it is also a much larger change, and the tiers are what make its absence
  measurable.
- **No benchmark runner or results format.** Invoking a setup is already one
  command. What to do with the timings is a separate question.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

- `pressure-solvers`: adds requirements that a setup intended for scaling
  measurements is sized so the multilevel hierarchy survives the decomposition
  at the scale it targets, and that such a setup performs a fixed amount of work
  independent of its grid. Both are properties of how the solver behaves under
  decomposition rather than of the files themselves, and both fail silently
  today — a collapsed hierarchy still converges, just slowly, and reports
  nothing unusual.

## Impact

- **New**: nine setups under `testcases/bench/`, and a check that certifies their
  hierarchy depth at the rank counts they target.
- **Modified**: `testcases/bench/README.md` (it currently says the tiers are yet
  to arrive); `tests/run-all.sh` if the depth check joins the suite; `README.md`
  to describe the tiers and how to run one.
- **Unmodified**: every solver source file, every existing setup, every recorded
  baseline. This change adds configuration and a check, and touches no code.
- **No new dependencies.**

## Sequencing

Independent of `optimize-multigrid`, and of no consequence to it. That change
raises `levels` in the validation setups; these are new files it does not touch.
Either may land first.

It does depend on `move-testcases`, which created `testcases/bench/` — already
archived.
