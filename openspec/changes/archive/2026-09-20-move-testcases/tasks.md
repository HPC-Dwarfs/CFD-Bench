# Tasks

## 1. Establish the reference the move is judged against

- [x] 1.1 Run `tests/run-all.sh` and confirm it is green end to end, so the move has a known-good state to be compared against; verify every step reports OK
- [x] 1.2 Take a copy of the recorded baselines from `tests/baseline/` before anything moves, so the comparison in task 5.2 is against bytes rather than against a re-recording; verify the copy holds one dump per setup and solver currently recorded

## 2. Move the files

- [x] 2.1 Create `testcases/flow/`, `testcases/regression/` and `testcases/bench/`, with `bench/` left empty for the later tier change and carrying a short note saying so; verify the three directories exist and `bench/` is tracked despite being empty
- [x] 2.2 Move `dcavity.par`, `canal.par`, `karman.par`, `backstep.par` and `schaefer-turek.par` from the repository root to `testcases/flow/` with `git mv`, changing no file content; verify `git status` reports five renames and no content modifications
- [x] 2.3 Move `dcavity-baseline.par`, `canal-baseline.par` and `sphere-baseline.par` from `tests/setups/` to `testcases/regression/` with `git mv`, changing no file content, and remove `tests/setups/` once empty; verify `git status` reports three renames and that `tests/setups/` is gone
- [x] 2.4 Confirm nothing was altered in transit: verify every moved file is byte-identical to its pre-move content

## 3. Point the scripts at the new layout

- [x] 3.1 Update the setup path in `check_setup()` in `tests/check-setups.sh`, which reaches all five flow cases through one expression; verify the five shipped setups are still found and still execute time steps
- [x] 3.2 Update the two bespoke variants in `tests/check-setups.sh` that name `backstep.par` and `karman.par` directly for their flow assertions; verify the backstep blocked-flow and karman transverse-velocity checks both still run and pass
- [x] 3.3 Update the two references to `schaefer-turek.par` in `tests/check-schaefer-turek.sh`; verify the shortened Schaefer-Turek run still completes and reports its coefficients
- [x] 3.4 Update the regression directory in `tests/record-baseline.sh` and the `sphere-baseline.par` path in `tests/check-solver-obstacle.sh`; verify both scripts locate their setups and run
- [x] 3.5 Confirm nothing was missed: verify a search for the old locations finds no reference to a root-level `.par` or to `tests/setups` outside this change's own artifacts and the archived changes

## 4. Update the documentation

- [x] 4.1 Update the shipped-setup table and the two run examples in `README.md` to the new paths; verify the table lists five setups and every path in it exists
- [x] 4.2 Correct the `../../canal.par` and `../../dcavity.par` comments inside the moved regression files to point at their new locations; verify the referenced paths exist
- [x] 4.3 Replace the two stale `./NusifSolver-CLANG` invocations in `README.md` with `./CFD-Solver-CLANG`, which the binary rename in `ef77c8e` missed; verify no occurrence of the old binary name remains anywhere in the repository

## 5. Prove the move changed nothing

- [x] 5.1 Run `tests/run-all.sh` end to end and compare against task 1.1; verify every step that passed before passes now, with no step newly skipped for a setup it can no longer find
- [x] 5.2 Compare the baselines against the copy taken in task 1.2; verify every dump is bit-identical, since a relocation that moves a number has not relocated but broken something
- [x] 5.3 Confirm the history is followable: verify `git log --follow` reaches the pre-move history for one moved file from each of the two groups
