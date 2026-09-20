# Proposal

## Why

This solver has no way to place a body in the domain. Every loop in
`computeFG`, `computeRHS`, `adaptUV` and all three pressure solvers is an
unconditional sweep over the full index range, and no parameter, array or file
format describes obstacle geometry. The benchmark can therefore run only the
lid-driven cavity and an empty canal, and has nothing to validate against the
published 3D cylinder references.

Adding obstacles the traditional way — a per-cell type field with mirror
conditions, the encoding the 2D sibling of this code used to carry — would buy a
representation that cannot express a thin plate or an isolated solid cell, and
that the pressure solvers would have to be taught about one at a time. This code
is an HPC benchmark with red-black SOR, a compressed-layout red-black variant
and multigrid all as first-class solvers, and a preconditioned CG solver as the
intended endpoint. All four need one operator, and CG needs it symmetric
positive definite.

Independently, the three pressure solvers carry defects that make geometry work
impossible to validate: the reported residual is not the residual, the red-black
colouring is derived from local indices so the answer depends on the rank count,
the pressure boundary condition ignores the setup's own configuration and leaves
the operator singular even where an outflow should pin it, and the multigrid
cycle never applies its coarse correction at all.

## What Changes

**Geometry representation**

- Add face apertures `Ax` (co-located with `u`), `Ay` (with `v`), `Az` (with
  `w`), and cell volume fractions `Lambda` (with `p`) — four `double` arrays of
  the same shape and allocation as the existing fields. 0 means fully solid, 1
  fully fluid; this change emits only 0 and 1.
- Add a `geometryFile` parameter key naming either a voxel volume file or an
  analytic body. Absent key means an obstacle-free domain, so both existing
  setups keep running unchanged.
- Define a raw voxel volume format — a minimal ASCII header followed by one byte
  per voxel — as the 3D counterpart of the binary PGM the 2D code reads. A voxel
  below 128 is solid. No third-party image or volume library.
- Each rank reads only the voxel slab covering its own subdomain plus the sample
  margin, by computing byte offsets from the header and seeking. Memory per rank
  scales with the rank's subdomain, not with the whole volume, which at the
  required sampling density reaches gigabytes.
- Keep a small analytic producer (sphere, axis-aligned cylinder, box,
  axis-aligned plate) behind the same interface, for order-of-accuracy work that
  a rasterized body would cap at its voxel size.
- Validate at initialization: reject a volume too coarse for the grid, and
  reject geometry whose fluid part is not a single connected region under
  6-connectivity, since each sealed pocket adds a null-space vector no solver
  here can handle.

**Discretization**

- Enforce no-slip through face apertures: a closed face carries no velocity
  unknown and no pressure flux. `computeFG` and `adaptUV` honour apertures;
  `adaptUV` no longer updates velocities on closed faces or inside solids.
- Formulate the pressure Poisson operator as a face-flux balance,
  `sum_f A_f (p_nb - p_c) / h^2 = Lambda * rhs_c`, over the six faces of a cell.
  Each face contributes a symmetric pair, so the operator is symmetric for any
  face weight — binary now, fractional later — and stays a compact 7-point
  stencil.
- Give fully solid cells identity rows (`aC = 1`, off-diagonals zero, `rhs = 0`)
  so the system stays SPD and the arrays stay rectangular, with no index
  compression in the hot kernel.
- Split the Poisson kernel into a branch-free bulk sweep over all cells reading
  no geometry, plus a correction pass over a compact list of cut cells. In 3D
  that list is O(obstacle surface area) against O(domain volume) of bulk work,
  so geometric complexity and solid volume are both free at runtime and the
  benchmark's headline kernel stays byte-for-byte comparable with and without
  obstacles.

**Solvers**

- All three solvers gain a per-cell diagonal in place of the hardcoded constant
  relaxation factor, and a third pass over the cut-cell list.
- The compressed red-black solver needs that pass in its own colour-split
  `(ic, j, k)` layout, so the cut-cell list carries both index forms and is
  built once at initialization rather than per `solve()` call.
- Multigrid applies the same correction inside its smoother and coarsens
  geometry: a coarse face aperture is the mean of the four fine faces it covers,
  a coarse volume fraction the mean of eight fine cells.
- Consolidate the pressure boundary condition, currently duplicated in four
  places, into one routine that follows the setup's own configuration — zero
  normal gradient at a wall or slip boundary, `p = 0` at an outflow — applied at
  every multigrid level.
- Drive `normalizePressure` from the boundary-condition configuration instead of
  a fixed step interval, and project the right-hand side onto the range of the
  operator, not just the iterate.
- Fix the pre-existing defects that block all of the above, each of which is
  reachable from the sources: the residual accumulator seeded at `1.0` and never
  reset; red-black colouring taken from local rather than global indices;
  the pressure boundary condition forced to a zero normal gradient on all six
  faces whatever the setup says; multigrid restriction applying a vertex-centred
  stencil to cell-centred data with two of its loop bounds swapped; prolongation
  writing into the residual array while the correction reads the error array, so
  the coarse correction never reaches the pressure; every level indexed with the
  finest stride and relaxed with the finest mesh size; the boundary condition
  applied only at the finest level; the coarse error never zeroed before a
  coarse solve; and one V-cycle per `solve()` with `eps` and `itermax` read into
  unused locals.

**Particle tracing**

- Add massless tracer particles, so that flow around an embedded body can be
  seen as streaklines rather than only as a pressure and velocity field. The 2D
  sibling of this code has this feature; this solver does not.
- Particles are injected from a configured seed region at a configured period
  from a configured start time, advected by the velocity interpolated from the
  staggered components, migrated between ranks as they move, and written
  periodically as VTK polydata that opens in the same tool as the field output.
- Obstruction is decided from the apertures of the faces a particle's path
  crosses, not from the volume fraction of the cell it lands in. Testing the
  destination cell — which is what the 2D code does — misses a body one face
  thick entirely, because such a body has volume fraction 1 on both sides, so a
  particle passes straight through it. This is the reason particle tracing
  belongs with the aperture work rather than after it.
- Four defects of the 2D implementation are not carried over, each of which is
  survivable in 2D and is not in 3D:
  - Migration buffers are declared as stack arrays sized by rank count times
    the subdomain cell count, and the particle pool is sized by the subdomain
    cell count. In 3D that is gigabytes of stack for a pool holding at most a
    few thousand particles. Storage is sized from the injection rate and grown
    on demand instead.
  - The obstruction test uses the destination cell's volume fraction, as above.
  - The tracer's time step is captured once at initialization from the
    parameter file and never updated, while the flow solver takes an adaptive
    step whenever `tau > 0`, so particles advance with the wrong step on most
    setups. The current step is used instead.
  - Migration sends to and receives from every rank every step, finding a
    particle's owner by scanning an all-gathered table of rank extents, and
    injection uses unseeded `rand()`, so seeds differ between runs and between
    rank counts. The Cartesian topology gives the owning rank directly, the
    exchange becomes a single collective, and injection becomes a deterministic
    sequence that does not depend on the decomposition.

**Out of scope (deliberate, follow-on work)**

- Fractional apertures in `(0,1)` and the second-order wall treatment and
  small-cell strategy they require. This change ships binary apertures, but
  ships the bulk/surface kernel split and the aperture data model that make
  fractional apertures purely additive.
- The preconditioned CG solver. This change establishes the SPD operator and the
  null-space handling it needs.
- Moving obstacles.
- Implementing `PERIODIC`, which is an empty case in every branch of
  `setBoundaryConditions` and is paired with `periods = {0,0,0}` in
  `commPartition`. The shared pressure boundary routine rejects it rather than
  silently treating it as a wall.
- OpenMP. `README.md` lists it and the build system has the flag, but no source
  file contains a single `omp` pragma.

## Capabilities

### New Capabilities

- `obstacle-geometry`: How obstacle geometry is specified (voxel volume or
  analytic producer), how it is turned into face apertures and cell volume
  fractions, how each rank obtains its own part without communication, and what
  geometry inputs are rejected.
- `embedded-boundary-discretization`: How the momentum predictor, the velocity
  correction and the pressure Poisson operator honour apertures; the symmetry,
  identity-row and stencil-compactness properties the operator must guarantee;
  the bulk/surface kernel split.
- `pressure-solvers`: What red-black SOR, compressed red-black SOR and multigrid
  must do to solve the embedded-boundary operator correctly, including the
  cut-cell pass, per-cell diagonals, aperture coarsening, transfer-operator
  consistency, null-space handling and decomposition independence.
- `particle-tracing`: How tracer particles are configured, seeded, advected,
  obstructed by embedded boundaries, exchanged between ranks and written out,
  and what determinism and storage bounds the implementation must hold to.

### Modified Capabilities

None. The project has no existing specs (`openspec list --specs` reports none).

## Impact

**Source**

- `src/discretization.h`, `src/discretization.c` — `Ax`, `Ay`, `Az`, `Lambda`
  added to `Discretization` and allocated alongside the existing fields
  (`:81-105`); `computeFG` (`:455`) and `adaptUV` (`:654`) rewritten against
  apertures; `computeRHS` (`:334`) weighted by `Lambda`; `normalizePressure`
  (`:403`) made conditional and extended to the right-hand side.
- `src/solver.h` — solvers need the aperture arrays and the cut-cell list.
- `src/solver-rb.c` — bulk/surface split, per-cell diagonal, global-checkerboard
  colouring, residual reset, pressure boundary condition moved out.
- `src/solver-rbc.c` — the same, plus the cut-cell pass in the compressed
  `PRED`/`PBLACK` layout; note its `_MPI` path is currently a copy of
  `solver-rb.c` and stays one.
- `src/solver-mg.c` — restriction, prolongation, per-level strides and mesh
  sizes, coarse-error zeroing, V-cycle loop, aperture coarsening per level.
- New `src/geometry.h`, `src/geometry.c` — voxel-volume reader, analytic
  producer, aperture derivation, resolution and connectivity validation,
  cut-cell list construction, run-header line and checksum.
- New `src/pressure-bc.h`, `src/pressure-bc.c` — the single pressure boundary
  routine replacing the four duplicated copies in the solvers.
- New `src/particletracing.h`, `src/particletracing.c` — particle pool,
  injection, advection, obstruction, rank migration and output, adapted from the
  2D implementation with the four defects above corrected.
- `src/parameter.h`, `src/parameter.c` — `geometryFile` added, plus the particle
  keys `numberOfParticles`, `startTime`, `injectTimePeriod`, `writeTimePeriod`
  and the six seed-region bounds. The parser matches keys with `strncmp` over
  the literal's length (`:56-58`), so a new key must not be a prefix of, nor
  prefixed by, an existing one; `geometryFile` satisfies that against the
  current key set, and the two-character seed bounds have to be checked against
  `xlength`, `ylength` and `zlength` rather than assumed safe.
- `src/comm.h`, `src/comm.c` — geometry needs each rank's global cell offsets.
  `commGetOffsets` (`:300-322`) already computes them but its whole body is
  inside `#if defined(_MPI)`, so a serial build leaves the array untouched; it
  gains a serial path writing zeros. `commPartition` (`:539-541`) pairs
  `coords[KDIM]` with `dims[ICORD]` and so on; this is numerically correct,
  because `KDIM == ICORD == 0`, but it is a trap for code that has to reason
  about which axis is which, so the enum use is corrected without behaviour
  change. Particle migration additionally needs a position-to-rank lookup over
  the Cartesian topology and one collective exchange of variable counts; neither
  exists today.
- `src/main.c` — geometry initialization before the time loop; the
  `normalizePressure` call site (`:68-70`); the particle tracer's
  initialization, its per-step call after `adaptUV` (`:72`), and its teardown.
- `Makefile` — a `TEST` build path for check drivers. New translation units are
  picked up by the existing `wildcard`, but any new *variant* family would need
  the same filter treatment as `solver-%.o` (`:30-32`).

**Inputs and assets**

- New `karman.par` (cylinder in a channel), `backstep.par` (3D backward-facing
  step) and `schaefer-turek.par` (the published 3D cylinder benchmark).
- New `geometry/` holding a versioned voxel volume per rasterized setup, with a
  content checksum printed in the run header so a silently edited volume cannot
  be mistaken for the reference case.
- New `tools/` — the volume generator, the generator driver script, and a field
  comparison utility.
- New `tests/` — check drivers and shell checks. The solver currently writes
  only VTK, so a raw field dump readable by the comparison utility is added
  under the test build path.
- New `vis_files/` output directory for the per-time particle files. `karman.par`
  gains a particle block; the other setups do not, so they stay untraced.
- `dcavity.par` and `canal.par` are unchanged and must stay bit-identical
  through the geometry work, which is what makes them the baseline.

**Behaviour**

- Obstacle-free results change where the defect fixes correct them — most
  visibly for multigrid, which currently applies no coarse correction, and for
  `canal`, whose outflow never reached the pressure operator. Red-black SOR
  results change where the colouring was rank-dependent.
- Reported residuals change in all three solvers, since the accumulator was
  never reset.

**Dependencies**

- None added. The voxel format is parsed directly.
