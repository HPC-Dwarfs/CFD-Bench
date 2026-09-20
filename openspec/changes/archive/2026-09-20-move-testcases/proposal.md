# Proposal

## Why

Setup files are scattered and their location no longer says anything about what
they are for. Five shipped setups sit in the repository root next to the
Makefile and the sources; three shortened regression variants sit in
`tests/setups/`. Nothing distinguishes a setup meant to be run from one that
exists only to be compared against a recorded field.

Two queued changes make that worse rather than better. `optimize-multigrid`
edits five of the existing files, and the benchmark tiers discussed alongside it
add around nine more. Both are cheaper against a layout that is already final,
and both would otherwise be written against paths that are about to move.

The groups also behave differently in ways the current layout hides:

- the shipped setups are physical cases, run at their own resolution and final
  time, with no recorded field;
- the regression variants are the same cases shortened and shrunk, and are the
  only ones with baseline dumps;
- the benchmark tiers, when they arrive, will have no baselines at all and exist
  to be run at scale.

A flat directory of eight files becomes seventeen once the tiers land, mixing
three purposes with no signal but a filename suffix.

## What Changes

- **All setup files move under `testcases/`**, grouped by what they are for:

  ```
  testcases/
    flow/        dcavity, canal, karman, backstep, schaefer-turek
    regression/  dcavity-baseline, canal-baseline, sphere-baseline
    bench/       (empty; the tiers land here in a later change)
  ```

- **The scripts follow.** Every shipped setup is reached through five sites, so
  this is a small edit: one path in `check-setups.sh`'s helper covers all five
  flow cases, two bespoke variants in the same file name `backstep` and `karman`
  directly, `check-schaefer-turek.sh` names its own setup twice, and
  `record-baseline.sh` and `check-solver-obstacle.sh` each name the regression
  directory once.
- **The README's setup table and its two run examples are updated**, and the two
  `../../canal.par`-style comments inside the regression files are corrected.
- **No setup content changes.** Not a grid, not a tolerance, not a final time,
  not a level count. The files move and nothing inside them is touched.

### Non-goals

- **Changing any setup's parameters.** `optimize-multigrid` raises `levels` in
  five of these files; that is its work and stays there. This change must be
  provably inert, and editing content would destroy the only verification a move
  has.
- **Adding the benchmark tiers.** `testcases/bench/` is created empty so the
  later change has somewhere to land, and nothing more.
- **Renaming any setup.** `dcavity.par` stays `dcavity.par`; only its directory
  changes, so git records renames and history follows.
- **Touching `$WORK`-generated parameter files.** The twenty-four `.par`
  references in `check-rejects.sh` are all temporaries written per run and are
  unaffected.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

None. The specs describe "the setup parameter file" as a concept and never name
a path or a directory, so no requirement changes. The change declares
`skip_specs: true` accordingly: it is a pure relocation with no observable
behaviour change, and inventing a requirement to satisfy validation would be
worse than declaring the truth.

## Impact

- **Moved**: `dcavity.par`, `canal.par`, `karman.par`, `backstep.par` and
  `schaefer-turek.par` from the repository root to `testcases/flow/`;
  `dcavity-baseline.par`, `canal-baseline.par` and `sphere-baseline.par` from
  `tests/setups/` to `testcases/regression/`.
- **Modified**: `tests/check-setups.sh`, `tests/check-schaefer-turek.sh`,
  `tests/record-baseline.sh`, `tests/check-solver-obstacle.sh` (paths only);
  `README.md` (the setup table and two run examples).
- **Unmodified**: every solver source file, the Makefile, the baselines
  themselves, and the contents of every setup that moves.
- **Verification**: `tests/run-all.sh` is green before and after, and the
  recorded baselines are bit-identical across the move. That is the whole
  acceptance criterion — a relocation that changes a result has not relocated,
  it has broken something.
- **No new dependencies.**

## Sequencing

This should land **before** `optimize-multigrid` and before the benchmark tiers.

Both of those touch setup files, and both are written against paths this change
replaces. More importantly, the repository is green now and its baselines were
re-recorded at `2698747`, which is the state a pure move most wants to be
verified against. `optimize-multigrid` moves the multigrid baselines again for
its deeper hierarchies, so doing the move afterwards would mean proving "nothing
changed" against a reference that had just changed for other reasons.

The cost of going first is one update to `optimize-multigrid`, whose tasks and
design name concrete setup paths. That is mechanical and smaller than the
alternative.
