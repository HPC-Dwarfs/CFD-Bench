# Tasks

## 1. The small tiers and the check that certifies them

- [ ] 1.1 Add `testcases/bench/dcavity-small.par` (48³ on a 1x1x1 domain), `testcases/bench/karman-small.par` (192x48x48 on 32x8x8 with `cylinder-z:5.0,4.0,1.0`) and `testcases/bench/canal-small.par` (384x48x48 on 32x4x4), each carrying the physics of the flow case it is named after and a header saying it is a benchmark tier, not a validation setup; verify each starts, reports the level count it built, and executes time steps
- [ ] 1.2 Confirm the cells are cubic in all three: verify `xlength/imax`, `ylength/jmax` and `zlength/kmax` agree to within a per cent for each small tier, so no tier measures the cost of stretched cells
- [ ] 1.3 Add `tests/check-bench-tiers.sh`, which reads every setup in `testcases/bench/`, computes the hierarchy depth it would reach at each rank count it targets using the same balanced factorisation `MPI_Dims_create` produces, and fails when a tier falls below four levels at its target; verify it passes for the three small tiers and fails if a tier's grid is edited to an extent that halves to an odd number
- [ ] 1.4 Extend that script to validate its own arithmetic: at 1, 8 and 27 ranks, run each small tier and compare the level count the solver reports against the computed figure; verify the two agree, so the check is tied to the implementation rather than restating its coarsening rule
- [ ] 1.5 Confirm the check would have caught the problem that motivated this change: verify it rejects a 128³ grid targeting 216 ranks, which `MPI_Dims_create` divides into 21 per direction and which therefore has no hierarchy

## 2. A fixed, measured step count

- [ ] 2.1 Set `tau 0` and an explicit `dt` in the three small tiers so the adaptive controller is bypassed; verify each run executes exactly `te / dt` steps and that the step count does not change when the Reynolds number or the initial velocity is altered
- [ ] 2.2 Measure what fraction of a small-tier run is startup — geometry sampling, hierarchy construction and the first solve — and choose a step count that amortises it while keeping the small tier quick enough to run often; verify the chosen count is recorded in each tier's header alongside the measurement that produced it
- [ ] 2.3 Apply the same `dt` and `te` to every tier of a case, so its three sizes differ only in cells; verify all three tiers of a case report the same number of steps

## 3. The medium and large tiers

- [ ] 3.1 Add the medium tiers — `dcavity` 96³, `karman` 384x96x96, `canal` 768x96x96 — with the physics, step count and domains of their small counterparts; verify each starts and that `tests/check-bench-tiers.sh` certifies four levels or more at 216 ranks
- [ ] 3.2 Add the large tiers — `dcavity` 192³, `karman` 768x192x192, `canal` 1536x192x192; verify each starts and that the check certifies four levels or more at 1728 ranks
- [ ] 3.3 Confirm the tiers form a ladder rather than three unrelated sizes: verify each case's tiers have extents in the ratio 1:2:4 and cell counts in the ratio 1:8:64
- [ ] 3.4 Run the small and medium tiers of each case at 1 and 8 ranks and record cycles and wall clock per solver, so the tiers have a first set of numbers and an obvious regression would be visible later; verify every run converges to its tolerance

## 4. Documentation and integration

- [ ] 4.1 Replace `testcases/bench/README.md`, which currently says the tiers are yet to arrive, with what they are: the sizing rule, the rank count each tier targets, and why extents are `3 * 2^k` rather than powers of two; verify the file names every tier that exists
- [ ] 4.2 Describe the tiers in `README.md` — what they are for, how to run one, and that they carry no baselines because they exist to be timed rather than compared; verify the documented command runs a small tier successfully from a clean tree
- [ ] 4.3 Decide whether `tests/check-bench-tiers.sh` joins `tests/run-all.sh`, and record the reasoning either way in the script's header; verify the suite passes with whichever choice is made
- [ ] 4.4 Confirm the change is additive: verify no existing setup, recorded baseline or solver source file was modified, and that `tests/run-all.sh` passes every step exactly as it did before
