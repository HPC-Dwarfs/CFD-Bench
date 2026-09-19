# Tasks

## 1. Verification harness

- [x] 1.1 Add a `TEST` build path to the Makefile that builds a named check driver instead of the normal `main`, and verify it builds and runs without disturbing the normal build of any `SOLVER` variant
- [x] 1.2 Add a raw binary dump of `p`, `u`, `v`, `w` over the global domain, written under the test build path only, and verify two consecutive runs of `dcavity.par` produce byte-identical dumps
- [x] 1.3 Add `tools/fieldcmp.c` reporting the max and L2 difference between two field dumps, and verify it reports exactly zero for a dump compared against itself
- [x] 1.4 Add `tests/run-checks.sh` that builds and runs every check driver at a given rank count, and verify it reports a nonzero exit status when a driver fails
- [x] 1.5 Add `tests/record-baseline.sh` and record reference dumps for `dcavity.par` and `canal.par` under each of `rb`, `rbc` and `mg`, and verify each recorded dump is reproduced by a fresh run
- [x] 1.6 Add `tests/run-all.sh` aggregating every check with a pass/fail summary, and verify it runs end to end against the current unmodified solver

## 2. Geometry-independent correctness fixes

- [x] 2.1 Reset the residual accumulator at the start of each iteration in all three solvers, and verify in a check driver that the reported residual of a converged solve matches an independently computed residual norm
- [x] 2.2 Give `commGetOffsets` a serial path that writes zero offsets, and verify a check driver reads the same offsets from a one-rank MPI build and a serial build
- [x] 2.3 Correct the `coords`/`dims` enum pairing in `commPartition` to name the axis it actually indexes, and verify local subdomain sizes are unchanged on 1, 2, 4 and 8 ranks
- [x] 2.4 Derive the red-black colour from each cell's global position in all three solvers, and verify a check driver assigns a given global cell the same colour on 1, 2, 4 and 8 ranks
- [x] 2.5 Extract the pressure boundary condition into one routine honouring `bcLeft`…`bcBack` — zero normal gradient at a wall or slip boundary, `p = 0` at an outflow — used by all three solvers, and verify the four duplicated copies are gone and `canal.par` converges
- [x] 2.6 Abort at initialization when a pressure boundary is configured `PERIODIC`, naming the boundary and the type, and verify no time steps execute for such a setup
- [x] 2.7 Replace multigrid restriction with the eight-child average for cell-centred data, with its loop bounds addressing the axes they name, and verify a constant field restricts to the same constant to machine precision
- [x] 2.8 Replace multigrid prolongation with trilinear interpolation reading `e[level+1]` and writing `e[level]`, defining every fine cell, and verify no fine cell retains a pre-filled sentinel value after prolongation
- [x] 2.9 Verify in a check driver that restriction followed by prolongation of a constant returns the same constant to machine precision
- [x] 2.10 Zero `e[level+1]` before each coarse solve, and verify a V-cycle on a Poisson problem with a known solution reduces the residual by more than pre- plus post-smoothing alone
- [x] 2.11 Index every multigrid level with its own stride instead of the finest one, and verify the restriction of a known fine field lands in the expected coarse cells
- [x] 2.12 Use each level's own mesh size in the smoother and in the residual, and verify the coarse-level operator applied to a known quadratic reproduces its Laplacian on that level's mesh
- [x] 2.13 Apply the pressure boundary condition at every multigrid level rather than only the finest, and verify a constant field is a fixed point of a full V-cycle with zero right-hand side
- [x] 2.14 Iterate multigrid V-cycles until `eps` is reached, bounded by `itermax`, reporting the cycle count, and verify `dcavity.par` and `canal.par` converge to the stated tolerance
- [x] 2.15 Drive the pressure mean removal from the boundary-condition configuration instead of a fixed step interval, and extend it to project the right-hand side onto the range of the operator, and verify an enclosed setup keeps its mean pressure bounded over many steps while an outflow setup is unchanged by it
- [x] 2.16 Verify all three solvers converge to the same field on `dcavity.par` and `canal.par`, and re-record the baselines from task 1.5 as the reference for the geometry work
- [x] 2.17 Copy the Cartesian coordinates, dimensions and neighbours onto each coarse multigrid communicator, and verify the boundary condition reaches every level
- [x] 2.18 Split parameter lines on any whitespace and trim the value, and verify a setup whose `name` line carries no trailing comment still matches its setup-specific boundary condition

Tasks 2.17 and 2.18 were added during implementation.

`commUpdateDatatypes` copied neither `coords` nor `dims` onto the coarse
communicator, and `commIsBoundary` reads both, so every coarse level answered
"am I on a physical boundary?" from uninitialized memory. Task 2.13 cannot be
done without this.

The parameter reader split on spaces alone and kept whatever trailing
whitespace a line carried, so `name dcavity` with no trailing comment produced
the string `"dcavity\n"`. Every setup-specific boundary condition is selected by
`strcmp` on that name, so it silently did nothing; the shipped setups only
worked because each happened to carry a trailing comment for the `#` strip to
remove. It surfaced when a reduced `dcavity` baseline written without comments
ran with a stationary lid and a pressure right-hand side of exactly zero.

Two decisions were taken during implementation and are recorded here because
they change the shipped setups rather than only the code.

The residual accumulator was seeded at `1.0`, never reset, and divided by the
cell count each iteration, so it fell below `eps^2` after about two iterations
whatever the real residual was. Every shipped setup was therefore running two
pressure iterations per step and reporting convergence. With the residual
computed honestly, red-black SOR needs roughly 1300 iterations for `dcavity` and
5900 for `canal` at the baseline resolutions, against budgets of 1000 and 500.
`itermax` was raised in both shipped setups so that the pressure solve actually
reaches its stated tolerance.

`canal.par` also shipped `levels 1`, which makes the multigrid solver a plain
relaxation. It is now `3`, and `dcavity.par` is `4`, so that the multigrid
variant is actually multigrid and the three solvers are comparable. Multigrid
reaches the same tolerance in 3 to 19 cycles where red-black SOR needs hundreds
to thousands of iterations.

## 3. Aperture data model

- [x] 3.1 Add `Ax`, `Ay`, `Az` and `Lambda` as `double` arrays to the discretization, allocated alongside the existing fields, and verify allocation sizes and alignment match the existing arrays
- [x] 3.2 Initialize all four to 1 everywhere as the obstacle-free default, and verify an obstacle-free run is bit-identical to the task 2.16 baseline for every solver
- [x] 3.3 Define the geometry-producer interface that fills the four arrays including the halo layer, taking grid and offset information rather than the solver structures, and verify a null producer leaves the obstacle-free default intact
- [x] 3.4 Give the solvers access to the aperture arrays, and verify the build has no remaining solver that reads geometry through a back channel

## 4. Geometry producers

- [x] 4.1 Add the `geometryFile` parameter key and its absence-means-obstacle-free default, and verify `dcavity.par` and `canal.par` still run unchanged and the key does not collide with any existing key under the parser's prefix matching
- [x] 4.2 Implement the analytic producers for a sphere, an axis-aligned cylinder, a box and an axis-aligned plate, and verify the represented solid volume matches the analytic volume within one cell layer
- [x] 4.3 Implement the voxel volume header parser, and verify it accepts a header with comment lines and arbitrary whitespace and reports the declared dimensions for a known small file
- [x] 4.4 Reject unreadable, malformed or truncated geometry files at initialization with a message naming the file and the reason, and verify no time steps execute for each rejection case
- [x] 4.5 Map the volume to the domain with the documented axis order and orientation, and verify with an asymmetric test volume that the obstacle appears in the expected octant of the written output
- [x] 4.6 Derive apertures and volume fractions by sub-cell sampling, and verify a fully fluid cell yields 1, a fully solid cell yields 0, and a one-voxel solid plane aligned with a grid face yields zero apertures on that plane with volume fraction 1 on both sides
- [x] 4.7 Close every face of a solid cell as a final producer pass, and verify no cell in any shipped geometry has zero volume fraction together with a nonzero face aperture
- [x] 4.8 Have each rank read only the voxel slab covering its subdomain plus the sampling margin, by computing byte offsets from the header, and verify the gathered aperture fields are bit-identical on 1 and 8 ranks including for an obstacle straddling a subdomain boundary
- [x] 4.9 Verify the per-rank geometry allocation scales with the subdomain rather than the whole volume, by running a volume far larger than one rank's share and reporting the allocated bytes per rank
- [x] 4.10 Print the geometry source, voxel dimensions and a content checksum in the run header, and verify two volumes differing in one voxel produce different checksums

## 5. Geometry validation

- [x] 5.1 Enforce a minimum of 4 voxels per cell per direction, aborting with the volume resolution, the grid resolution and the required minimum, and verify the abort triggers for an under-resolved volume and not for a conforming one
- [x] 5.2 Warn when any of the volume's three side ratios differs from the corresponding domain ratio beyond tolerance, naming both ratios, and verify the run continues
- [x] 5.3 Implement the connected-component check over fluid cells under 6-connectivity, parallel-consistent across ranks, and verify a single connected region passes and reports the same count on 1 and 8 ranks
- [x] 5.4 Abort when more than one fluid region is found, reporting the count and one cell location per region, and verify with a volume containing a sealed pocket
- [x] 5.5 Verify the represented solid volume for one file agrees across two grids differing by a factor of two in resolution, within the coarser grid's sampling tolerance
- [x] 5.6 Verify the solid volume produced by the analytic path and by a high-resolution rasterization of the same body agree within the volume's sampling tolerance

## 6. Discretization against apertures

- [x] 6.1 Rewrite the momentum predictor to derive no-slip from face apertures, and verify an obstacle-free run stays bit-identical to the task 3.2 baseline
- [x] 6.2 Rewrite the velocity correction so it applies no pressure gradient across a closed face and does not update velocities on closed faces or inside solid regions, and verify that perturbing the pressure inside a solid region before the correction leaves the fluid velocity unchanged
- [x] 6.3 Weight the right-hand side by the cell volume fraction, and verify an obstacle-free run is unchanged
- [x] 6.4 Assemble the pressure operator as the aperture-weighted six-face flux balance with a per-cell diagonal, and verify symmetry to machine precision for random vectors, both with and without an obstacle
- [x] 6.5 Give cells with zero volume fraction identity rows with zero right-hand side, initialize their pressure to zero, and verify their pressure remains exactly zero after a converged solve
- [x] 6.6 Compute residual norms over fluid unknowns only, normalized by the fluid cell count, and verify the reported norm is unchanged when the solid fraction grows while the fluid region and its solution stay the same
- [x] 6.7 Verify a geometry with an isolated solid cell, and one with a zero-thickness plate, both run to completion and produce a solution with no flow across the plate

## 7. Kernel split and solvers

- [x] 7.1 Build the surface cell list at initialization — fluid cells with at least one closed face, plus solid cells with at least one fluid face neighbour — in two global-checkerboard colour segments, each sorted row-major, with precomputed face coefficients and inverse diagonal, and verify its length scales with obstacle surface area rather than solid volume across two obstacle sizes
- [x] 7.2 Store both the natural linear index and the compressed `(ic, j, k)` triple per list entry, and verify a check driver reads the same cell through either form
- [x] 7.3 Split the red-black sweep into a geometry-free bulk sweep over all cells plus a correction pass over the surface list, with a halo exchange between the two colour segments, and verify the converged result matches the pre-split implementation to the solve tolerance
- [x] 7.4 Verify that solid cells not on the surface list stay at exactly zero through a full solve, asserting it in a check driver
- [x] 7.5 Order the cut-cell pass after the colour sweeps, and verify red-black SOR converges on an obstacle setup with iteration count within the agreed factor of the obstacle-free case
- [x] 7.6 Apply the cut-cell pass in the compressed solver's `PRED`/`PBLACK` layout, and verify its converged field on an obstacle setup matches the red-black solver's to the solve tolerance
- [x] 7.7 Apply the same bulk-plus-correction structure inside the multigrid smoother, and verify a multigrid solve matches the red-black solve on the same obstacle setup to the solve tolerance
- [x] 7.8 Coarsen apertures and volume fractions through the multigrid hierarchy — four-face mean for apertures, eight-cell mean for volume fractions — and verify the obstacle is represented at every level for a resolved body
- [x] 7.9 Report the level at which an obstacle feature becomes unresolved, and verify the solver still converges for an obstacle that vanishes on the coarsest level
- [x] 7.10 Verify asymptotic multigrid residual reduction per cycle on an obstacle setup is within the agreed factor of the obstacle-free case
- [x] 7.11 Verify all three solvers produce the same converged pressure field and equal iteration counts on 1 and 8 ranks for the same obstacle setup

Three notes from implementing sections 6 and 7.

Guarding the velocity correction with a branch on the aperture, rather than
multiplying by it, changes what the compiler contracts under `-ffast-math` and
moved the last bit of every obstacle-free result. The correction multiplies by
the aperture instead: an aperture is 0 or 1, so the open case is bit-for-bit the
update the solver always did, and the closed case gives the zero the velocity on
a closed face is required to have.

The multigrid residual is now the aperture-weighted operator rather than the
plain seven-point stencil, which regroups the arithmetic. That moves
obstacle-free multigrid results by 2e-16 to 5e-15 -- rounding, not behaviour,
but it does mean the recorded multigrid baselines were re-recorded at this
point. The relaxation solvers are unaffected, because their sweep is unchanged
and the aperture-weighted residual only ever reports.

The bulk sweep writes a wrong value into every listed cell before the
correction runs, so the correction cannot simply adjust what it finds: it needs
the value the cell held beforehand. Each colour therefore saves its listed cells
before the sweep and recomputes from that. The saved array is O(surface), so
this costs nothing that scales, and it is what keeps the interior sweep free of
any branch or geometry read.

## 8. Setups, assets and documentation

- [ ] 8.1 Add `tools/genvox.py` rasterizing a body to the voxel format, and verify a volume it writes round-trips through the solver's parser to the expected voxel values
- [ ] 8.2 Add `tools/make-geometry.sh` regenerating the shipped volumes at resolutions that clear the 4-voxels-per-cell minimum on each setup's own grid and on one refinement of it, and verify the produced volumes pass the resolution and connectivity checks
- [ ] 8.3 Add `karman.par` with a cylinder in a channel, its volume and a particle block, and verify the run starts without parameter warnings and produces a wake
- [ ] 8.4 Add `backstep.par` with a box step and its volume, and verify the run starts and the step blocks flow through its volume
- [ ] 8.5 Add `schaefer-turek.par` using the analytic cylinder producer, and verify the run starts at the geometry the published benchmark specifies
- [ ] 8.6 Update `README.md` with the voxel format, the resolution rule, the analytic producers, the particle parameters and the test harness, and verify the documented commands run as written

## 9. Verification and benchmarking

- [ ] 9.1 Verify the aperture-weighted divergence over fluid cells adjacent to the obstacle is within the solve tolerance, matching the rest of the domain
- [ ] 9.2 Verify the velocity component normal to every zero-aperture face is zero to machine precision after a completed time step
- [ ] 9.3 Verify net volume flux through channel cross-sections is constant at steady state for flow around an obstacle, and that flux through the body surface is zero
- [ ] 9.4 Verify the operator couples each cell only to its six face neighbours by applying it to a unit vector at a cell adjacent to an obstacle corner and at one adjacent to an obstacle edge
- [ ] 9.5 Measure interior-sweep time with no obstacle and with a small obstacle on the same grid for each of the three solvers, and verify they agree within measurement noise
- [ ] 9.6 Measure interior-sweep time for a simple and a geometrically complex obstacle of comparable solid volume, and verify they agree within measurement noise
- [ ] 9.7 Run the Schäfer–Turek 3D cylinder case and verify drag, lift and Strouhal number are reported against their published reference ranges
## 10. Particle tracing

- [ ] 10.1 Add the particle parameters — count, start time, injection period, write period and the six seed-region bounds — and verify each is read back correctly and that none collides with `xlength`, `ylength` or `zlength` under the parser's prefix matching
- [ ] 10.2 Add the particle pool with heap storage sized from the injection rate and grown on demand, plus the compaction pass that reclaims removed particles, and verify storage per rank is substantially unchanged when the grid is refined by a factor of two in each direction
- [ ] 10.3 Implement deterministic seed-position generation from the global particle index and batch number, and verify two runs and two rank counts inject identical seed positions
- [ ] 10.4 Implement injection into the seed region from the start time at the injection period, skipping positions inside solid regions, and verify no particle is created before the start time, all lie in the seed region, and none lies in a cell of zero volume fraction
- [ ] 10.5 Implement trilinear interpolation of each velocity component at its own staggered face location, and verify a particle in a uniform flow of known velocity displaces by velocity times elapsed time
- [ ] 10.6 Advance particles with the flow solver's current time step, and verify displacement over a step where the adaptive step changed matches the step actually taken rather than the parameter-file value
- [ ] 10.7 Implement obstruction by walking the faces the particle path crosses and stopping at the first zero aperture, removing and counting the particle, and verify no particle appears on the far side of a zero-thickness plate
- [ ] 10.8 Verify no written particle position lies inside a cell of zero volume fraction over a long run around a closed body
- [ ] 10.9 Implement position-to-rank lookup through the Cartesian topology and migration by one collective exchange of variable counts, and verify a particle crossing a subdomain boundary is held afterwards by exactly one rank, and that one is the rank owning its new position
- [ ] 10.10 Remove and count particles advected past a domain boundary, and verify injected equals written plus removed-at-boundary plus removed-at-body at the end of a run
- [ ] 10.11 Write particle positions as VTK polydata at the write period, gathered to rank 0, and verify a particle file and the field output for the same time load together and the particles lie within the domain
- [ ] 10.12 Verify a setup that configures no particles produces fields, output files and run time bit-identical to the same setup before tracing existed
- [ ] 10.13 Verify the set of particle positions written at a given time agrees on 1 and 8 ranks to within the reproducibility limits of the flow field

## 11. Aggregate

- [ ] 11.1 Add every check from sections 4 through 10 to `tests/run-all.sh`, and verify the aggregate run reports a single pass/fail summary at 1 and 4 ranks
