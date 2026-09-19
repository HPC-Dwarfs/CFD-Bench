# Design

## Context

See `proposal.md` — Why. The constraints that shape this design:

- The code is an HPC benchmark. The pressure sweep is the headline kernel and is
  memory-bound. In 3D each cell already streams a 7-point stencil against ~12
  flops, so any per-cell geometry data added to that loop is paid directly in
  bandwidth. Seven coefficient arrays is not an acceptable way to carry
  geometry.
- Three solvers are first-class and a preconditioned CG solver is intended. All
  four must share one operator, and CG requires it symmetric positive definite.
- One of the three, `src/solver-rbc.c`, exists specifically to get stride-1
  inner loops by splitting the unknowns into separate red and black arrays with
  a compressed x index. Anything that forces indirect addressing into its inner
  loop destroys the reason that file exists.
- Multigrid needs the operator to coarsen, which rules out approaches based on a
  large penalization coefficient.
- Geometry flexibility is wanted but is the least important axis, so the
  geometry producer should be small and dependency-free.

Existing state the design has to work around, all verified in the sources:

- There is no obstacle machinery at all — no flag array, no cell types, no
  parameter key. This change is additive; nothing is being replaced.
- `res` is seeded to `1.0` and never reset between iterations in any solver
  (`src/solver-rb.c:50`, `src/solver-rbc.c:60`, `src/solver-mg.c:197`, `:252`),
  so the reported residual is not the residual.
- The red-black colour starts from the local index 1 on every rank
  (`src/solver-rb.c:55`, `src/solver-mg.c:199`, `:254`) and `src/solver-rbc.c:167`
  splits on local `(i+j+k)%2`, so the colouring — and the answer — depend on the
  decomposition.
- All three solvers force a zero normal pressure gradient on all six boundaries
  regardless of `bcLeft`…`bcBack` (`src/solver-rb.c:81-148`,
  `src/solver-rbc.c:93-139` and `:236-277`, `src/solver-mg.c:102-175`), so the
  velocity boundary condition never reaches the pressure operator and the
  operator is singular in every setup — including `canal`, which has outflows on
  two faces.
- Multigrid has several defects that compound: `restrictMG` applies a 27-point
  vertex-centred full-weighting stencil to cell-centred data and swaps two of
  its loop bounds (`src/solver-mg.c:34-36`); `prolongate` writes piecewise
  injection into `s->r[level]` while `correct` reads `s->e[level]`
  (`:67-100`), so the coarse correction never reaches the pressure and the
  finest level degenerates to pre- plus post-smoothing; every level is allocated
  at the global finest size and indexed with the finest local stride (`:378`,
  `:14-21`); `smooth` and `calculateResidual` use the finest `dx, dy, dz` on
  every level (`:182-190`, `:236-246`), so a coarse operator is wrong by
  `4^level`; the boundary condition is applied only at `FINEST_LEVEL`
  (`:319-320`, `:343`, `:349`); `e[level+1]` is never zeroed before a coarse
  solve (`:334`); and `solve` runs exactly one V-cycle, reading `eps` and
  `itermax` into locals it never uses (`:391-402`).
- `normalizePressure` (`src/discretization.c:403`) exists and works, but runs
  every 100 steps unconditionally (`src/main.c:68-70`) whether or not the system
  is singular, and touches only the iterate, not the right-hand side.
- `commGetOffsets` (`src/comm.c:300-322`) computes each rank's global cell
  offsets but its entire body is inside `#if defined(_MPI)`, so a serial build
  leaves the caller's array untouched.
- `commPartition` (`src/comm.c:539-541`) pairs `coords[KDIM]` with
  `dims[ICORD]`, `coords[IDIM]` with `dims[KCORD]`. This is numerically correct
  — `KDIM`, `ICORD` are both 0 and `IDIM`, `KCORD` are both 2 — but it is a trap
  for code that has to reason about which axis is which, which geometry does.
- The solver writes VTK only. There is no field format two runs can be diffed
  in, and no test harness of any kind.
- There is no particle tracing. The 2D sibling has an implementation
  (`2D-mpi/src/particletracing.c`, 636 lines) that is the starting point, but
  four of its properties do not survive the move to 3D and are treated as
  defects to correct rather than behaviour to reproduce; each is called out at
  the decision that replaces it.

## Goals / Non-Goals

**Goals**

- One geometry representation that serves the momentum predictor, the velocity
  correction and the pressure operator, and that all three solvers can consume
  without any of them learning what an obstacle is.
- An interior pressure sweep whose cost and instruction mix are unaffected by
  geometry, in all three solvers, so benchmark numbers remain comparable across
  setups.
- A data model and kernel structure in which introducing fractional apertures
  later is a change of values, not of types, loops or solver structure.
- Geometry whose memory cost per rank scales with the rank's subdomain, since a
  3D volume at the required sampling density is two orders of magnitude larger
  than the 2D image this design is adapted from.
- A particle tracer that is inert unless a setup asks for it, so it cannot
  perturb a benchmark measurement, and whose obstruction test is a statement
  about apertures rather than about cells.

**Non-Goals**

- Second-order wall treatment, sub-cell reconstruction and the small-cell
  problem. Binary apertures only here.
- Changing the projection method, the time-stepping scheme, the momentum
  discretization away from the obstacle, or the MPI decomposition strategy.
- Performance tuning of the surface pass. It is O(surface area) and will not be
  on the critical path.
- Rewriting `src/solver-rbc.c`'s `_MPI` path, which is currently a copy of
  `src/solver-rb.c`. It stays a copy; only the shared changes reach it.
- Higher-order particle integration, particle-laden flow, or any particle that
  influences the flow. The tracers here are massless and passive.

## Decisions

### Geometry lives on faces, not cells

Store x-face apertures `Ax` (co-located with `u`), y-face apertures `Ay` (with
`v`), z-face apertures `Az` (with `w`) and cell volume fractions `Lambda` (with
`p`). Four `double` arrays with the same shape, allocation and alignment as the
existing fields.

The obvious alternative — a per-cell type field with mirror boundary conditions,
which is what the 2D sibling of this code carried before its own
embedded-boundary change — fails because a face carries one velocity degree of
freedom while a cell-type mirror condition needs exactly one fluid side to
mirror from. A thin wall demands two contradicting mirrors, and in 3D a solid
cell can have fluid on any of six sides, so the number of inadmissible
configurations is larger still. Putting geometry on faces removes the
contradiction rather than working around it: a zero-thickness plate is a plane
of closed faces with no solid cell at all, and an isolated solid cell is six
closed faces. The category "inadmissible cell" never comes into existence.

It also makes the two roles of geometry the same object: `Ax(i,j,k) == 0` means
simultaneously "`u(i,j,k)` is not an unknown" and "no pressure flux crosses this
face". Momentum and pressure see identical geometry by construction.

_Alternatives considered._ Cell types plus an admissibility repair pass —
rejected, it preserves the representational limit and forbids thin bodies.
Volume penalization — rejected, the stiff coefficient jump breaks multigrid
coarse-grid correction and moves the difficulty into the momentum solve.
Non-symmetric cut-cell flux reconstruction — rejected, it forecloses CG.

### The operator is a face-flux balance, which makes symmetry structural

Assemble the pressure operator as `sum_f A_f (p_nb - p_c) / h^2 = Lambda *
rhs_c` over the six faces. Each face contributes `+A_f` to one row and `-A_f` to
the other, so symmetry holds for _any_ face weight. This is why binary apertures
now and fractional apertures later need no change in the solvers: the SPD
property, the 7-point structure and the CG-readiness are properties of the
assembly, not of the values.

The diagonal becomes per-cell, `aC = -sum_f A_f`, replacing the domain-wide
constant currently computed once per solve in each of the three solvers
(`src/solver-rb.c:46`, `src/solver-rbc.c:56`, `src/solver-mg.c:191`).

### A face is open only if both cells it separates are fluid

The volume fraction and the face apertures are sampled independently, and
rounding each to binary separately can disagree: a cell that is mostly solid
rounds to `Lambda = 0` while a face along its edge, sampled on a plane that
falls in the fluid, rounds to open. That leaves a solid cell with a face the
flow can cross, and a pressure inside the body that reaches the fluid —
precisely what the identity rows below assume cannot happen.

The producer therefore closes every face of a solid cell as a final pass. The
2D change discovered this during implementation; it is adopted here up front.

The converse is deliberately not imposed: a closed face between two fluid cells
is a zero-thickness plate, which is a body this model exists to represent. One
consequence follows and is worth stating, because it is the configuration a
cell-type encoding could not express: a solid cell open on opposite sides can no
longer be produced. Its representable form in this model is a plane of closed
faces with no solid cell at all.

### Solid cells stay in the vector as identity rows

Cells with `Lambda == 0` get `aC = 1`, zero off-diagonals, `rhs = 0`, and stay
part of the rectangular arrays.

Compacting the unknown vector to fluid cells only would force indirect
addressing into the hot kernel. In `src/solver-rb.c` that costs the vectorized
sweep; in `src/solver-rbc.c` it costs the stride-1 inner loop over the
compressed index, which is the only reason that variant exists. Identity rows
keep SPD, keep the arrays rectangular, and cost only the flops already being
spent on those cells.

### Bulk sweep plus a surface list, with solid interiors self-maintaining

The pressure kernel splits into:

1. a bulk sweep over all cells applying the plain constant-coefficient 7-point
   stencil, reading no geometry at all;
2. a correction pass over a compact list of cells where the bulk result is
   wrong.

The load-bearing observation is that the list is O(surface area), not O(solid
volume). If pressure is initialized to zero in solid cells and their right-hand
side is zero, a bulk relaxation of a solid cell whose six neighbours are all
solid computes a zero residual and leaves the value at exactly zero. Deep solid
interiors therefore stay correct with no correction at all. Only two classes of
cell need the list:

- fluid cells with at least one closed face;
- solid cells with at least one fluid face neighbour.

In 3D this is O(N^2/3) against O(N) bulk work, so the margin is wider than in
2D, not narrower — geometric complexity and solid volume both become free at
runtime.

Store the list as arrays rather than an array of structs, sorted row-major
within each colour segment, with per-entry precomputed face coefficients and
inverse diagonal. Sorting keeps the pass cache-friendly and makes it
deterministic across runs.

_Alternative considered._ A masked bulk kernel reading one byte of packed face
flags per cell — simpler, no third relaxation pass. Rejected because it must be
rewritten into the split form when fractional apertures arrive, and the rewrite
would land on the benchmark's headline kernel. It is also worse in 3D than in
2D: six face flags per cell rather than four, against a stencil that is already
bandwidth-bound.

### The surface list carries two index forms, and is built once

`src/solver-rb.c` and `src/solver-mg.c` address the unknowns by the natural
linear index. `src/solver-rbc.c`'s serial path addresses them as `PRED(ic,j,k)`
and `PBLACK(ic,j,k)` with `ic = i/2`. One list serves all three by storing, per
entry, both the natural linear index and the compressed `(ic, j, k)` triple.

The list is built at initialization, not inside `solve()`. `src/solver-rbc.c`
already allocates and gathers its four colour arrays once per `solve()` call,
i.e. once per time step; adding a per-call geometry build on top of that would
put an O(surface area) allocation on the time-step path for no reason.

### The surface list is coloured on the global checkerboard

Relaxing the list in one pass makes it a Gauss-Seidel over the list, so a cut
cell next to another cut cell sees its neighbour's new value or its old one
depending on traversal order — and across a rank boundary the neighbour is in
the halo and is never the new one. The iteration would then depend on how the
domain was divided.

A Jacobi block would be order-independent but diverges at the relaxation factor
these solvers use. So the list is built in two colour segments of the global
checkerboard, with a halo exchange between them: two cut cells of the same
colour are never face neighbours, so within a colour the order cannot matter,
and the second colour sees the first identically on every rank.

The colour is `(i + iOffset + j + jOffset + k + kOffset) % 2`, and the bulk
sweep is recoloured the same way. That is also the fix for the pre-existing
rank-dependence, so the two are one change rather than two. It requires each
rank's global offsets, which `commGetOffsets` already computes; it gains a
serial path writing zeros, since its body is currently entirely inside `#if
defined(_MPI)`.

### The voxel volume is read by seeking, not by slurping

The format is the 3D counterpart of the binary PGM the 2D code reads: an ASCII
header — magic, optional comment lines, `nx ny nz`, maxval — followed by
`nx*ny*nz` raw bytes with x varying fastest, then y, then z. A voxel below 128
is solid. The parser is a few dozen lines and needs no library.

```
P5V
# karman cylinder, 30 x 8 x 8 domain
1600 800 800
255
<nx*ny*nz raw bytes>
```

The 2D design has every rank read the whole image and keep its own slice. That
does not carry over. At the minimum sampling density of 4 voxels per cell per
direction a volume is 64 bytes per cell, so a benchmark grid of 400 x 200 x 200
needs 1600 x 800 x 800 voxels — one gigabyte, per rank.

Instead a rank computes the byte offset of each voxel row directly from the
header and reads only the rows covering its own subdomain plus the sampling
margin. Memory per rank is O(local voxels). Apertures stay bit-identical across
rank counts because the sampling of a given cell reads exactly the same voxels
either way — the slab bounds change what is in memory, not what is read for any
particular cell. No communication is involved, so the rank-count independence is
a property of the construction rather than something to be checked for.

_Trade-off._ Many small positioned reads instead of one large sequential one.
The reads are row-contiguous and the whole thing happens once at initialization,
so this is not on any hot path. MPI-IO can replace it later without changing the
interface.

_Alternative considered._ A stack of 2D PGM slices, one per z-plane, which would
reuse the 2D parser verbatim and stay human-inspectable. Rejected: it turns one
file into hundreds or thousands, and turns a slab read into a per-file open.

### Sampling produces apertures and volume fractions, and nothing else reads the file

A cell's volume fraction is the fraction of covered voxels that are fluid,
counted over the sub-cell box; a face aperture is the same count over the
sub-face rectangle. With binary rounding this is a threshold on that fraction.
The same counting produces fractional values unchanged in a later phase, which
is why a raster input makes fractional apertures easier rather than harder.

The four arrays are the sole interface between geometry and the rest of the
solver, including the halo layer — filled directly from the volume rather than
exchanged, so no communication is needed to make them consistent.

### Apertures are stored as `double` from the start

Even though this change's values are only 0 or 1. They are never read in the
bulk sweep, so their footprint does not enter the streaming working set, and
keeping the type means a later phase changes values only.

### Analytic producers sit behind the same interface

Sphere, axis-aligned cylinder, box and axis-aligned plate, selected through the
same `geometryFile` value, because a rasterized body is only as accurate as its
voxel size, which would cap any order-of-accuracy study. They are also what
`schaefer-turek.par` uses, since the published reference ranges are for an exact
cylinder.

### Multigrid coarsens apertures geometrically

Coarse face aperture = mean of the four fine faces it covers; coarse volume
fraction = mean of the eight fine cells. Cheap, local, and consistent with the
flux-balance assembly.

Galerkin (RAP) coarsening would be more robust for under-resolved features but
produces a 27-point coarse stencil in 3D, changing the coarse kernels and their
memory layout. Given that multigrid-preconditioned CG is the intended endpoint,
and CG absorbs the deficiency of an imperfect coarse-grid correction, geometric
coarsening is the right first choice. Where a feature becomes unresolved at some
level, the solver reports the level rather than silently degrading.

### The multigrid cycle is repaired as part of this change

The defects listed in Context are in scope on one ground: the geometry work
cannot be validated on top of them. Specifically — restriction becomes the
eight-child average for cell-centred data with its loop bounds corrected;
prolongation becomes trilinear interpolation from `e[level+1]` into `e[level]`,
writing every fine cell, with `e[level+1]` zeroed before each coarse solve; each
level gets its own stride and its own mesh size; the boundary condition is
applied at every level; and `solve` gains a V-cycle loop bounded by `eps` and
`itermax`, reporting the cycle count. Without the last one there is no tolerance
for multigrid to converge to and therefore no way for it to agree with the
relaxation solvers, which is a requirement in
`specs/pressure-solvers/spec.md`.

The residual accumulator is reset at the start of each iteration in all three
solvers, independently of multigrid.

### An outflow pins the pressure; only then is the operator non-singular

All three solvers impose a zero normal pressure gradient on all six physical
boundaries whatever the setup says, so the velocity boundary condition never
reaches the pressure operator and the operator is singular in every setup —
including `canal`, which has outflows, and which the null-space handling would
therefore decline to treat.

The pressure boundary condition now follows the discretization: a wall or slip
boundary sets the normal gradient to zero and the halo mirrors the interior; an
outflow pins `p = 0` on the boundary face, which for a cell-centred field is an
odd reflection. `PERIODIC` is rejected at initialization, because
`setBoundaryConditions` implements it as an empty case for velocity and
`commPartition` creates the Cartesian communicator with `periods = {0,0,0}` —
treating it as a wall for pressure would be a silent disagreement with a
velocity condition that does not exist.

One routine shared by all three solvers applies it, at every multigrid level,
since the homogeneous form of either condition is the same and coarse levels
carry the error rather than the pressure. This replaces four duplicated copies,
which is also why it is worth extracting rather than fixing in place.

This is what makes "an outflow setup is non-singular" true of the code and not
only of the specification, and it is what lets the singular case be detected
from the boundary-condition configuration rather than from the matrix.

### Null space handled by range projection plus mean removal

Where every boundary imposes a zero normal pressure gradient the system is
singular. Project the right-hand side onto the range of the operator (subtract
its fluid-cell mean, `Lambda`-weighted) and remove the constant from the
iterate. Detect the singular case from the boundary-condition configuration
rather than testing the matrix.

`normalizePressure` already removes a mean, but over all cells rather than fluid
cells, unconditionally rather than when singular, every hundredth step rather
than every solve, and only from the iterate. All four of those change. The
relaxation solvers tolerate the drift today; multigrid coarse solves and CG do
not.

### More than one fluid region is rejected

Each sealed fluid pocket contributes an independent constant to the null space.
Pinning one degree of freedom per component would generalize this, but it puts
per-component bookkeeping into every solver and every multigrid level. For a
benchmark with curated inputs, rejecting at initialization with a clear
diagnostic is the better trade. Connectivity is checked under 6-connectivity and
is a global property, so unlike the producers it needs the communicator.

### Particles are obstructed by the faces they cross, not by where they land

The 2D implementation advances a particle, then deletes it if the cell it landed
in has `Lambda == 0` (`2D-mpi/src/particletracing.c:358-364`). That test cannot
see a body one face thick: a zero-thickness plate is a plane of closed faces
between cells whose volume fractions are both 1, so a particle crosses it and
the destination-cell test reports nothing wrong. It also depends on the step
being short enough that a particle cannot step over a thin solid region
entirely, which adaptive time stepping does not guarantee.

So obstruction is decided from the path: walk the cell faces the segment from
the old to the new position crosses, in order, and stop at the first whose
aperture is zero. The particle is then removed and counted.

This is the reason particle tracing belongs in this change rather than after it.
The property the 2D spec claims — "no particle is advected through a closed
face" — is a statement about apertures, and it is only checkable once apertures
exist. Bolting tracing on afterwards would mean writing the destination-cell
test first and then replacing it.

_Alternatives considered._ Stopping the particle at the face, which shows where
flow impinges but piles tracers up on the surface permanently; and reflecting
it, which is not physical for a tracer in a no-slip flow and needs a sub-step
search. Deletion matches what a streakline visualization wants.

### Velocity is interpolated at each component's own face location

Each of `u`, `v`, `w` lives on a different face of the cell, so a single
interpolation stencil cannot serve all three. Interpolate each component
trilinearly from the eight surrounding values of that component, using the
offset appropriate to its own staggering — the 3D form of what the 2D code does
with two bilinear interpolations at two different offsets
(`2D-mpi/src/particletracing.c:298-330`).

Integration stays forward Euler, matching the momentum equation's own time
discretization; a higher-order particle integrator would be more accurate than
the flow field it samples.

### The tracer uses the flow solver's current time step

The 2D tracer captures `dt` once at initialization from the parameter file
(`2D-mpi/src/particletracing.c:92`) and never updates it, while `computeTimestep`
recomputes the step every iteration whenever `tau > 0` — which is every shipped
setup. Particles therefore advance with a step the flow is not taking. The
tracer reads the current step instead. This is a correctness fix, not a design
choice, but it is recorded here because the tracer's structure invites it: the
tracer holds a copy of everything it needs, and a copy of an adaptive quantity
is wrong by construction.

### Particle storage is sized by particles, not by cells

The 2D implementation sizes its pool at `imaxLocal * jmaxLocal` and, worse,
declares its migration buffers as stack arrays of
`size * estimatedNumParticles` particles (`2D-mpi/src/particletracing.c:285-288`).
In 2D at 400 x 200 per rank that is already tens of megabytes of stack; the 3D
equivalent multiplies it by `kmaxLocal` and by the larger particle record, which
does not run.

Storage is instead sized from what the configuration actually produces — the
injected count per batch times the number of batches in flight — on the heap,
grown when it fills, and compacted when removals accumulate. The 2D `compress`
pass is kept; the sizing and the buffers are not.

### Migration goes through the Cartesian topology and one collective

The 2D code finds a migrating particle's destination by scanning an all-gathered
table of every rank's extents (`2D-mpi/src/particletracing.c:341-354`), then
sends to and receives from every other rank every step
(`2D-mpi/src/particletracing.c:368-395`) — `O(size^2)` messages per step, almost
all of them empty.

A particle's owner follows directly from its position: divide by the cell size,
find which rank's index range contains the result, and ask the Cartesian
communicator for that rank. That is `O(1)` per particle and needs no gathered
table. The exchange itself becomes one collective with variable counts, so empty
pairs cost nothing.

Unlike the geometry, which is derived locally and needs no exchange at all, a
particle can in principle land in any rank's subdomain within one step, so this
is a genuine all-to-all in shape even though it is nearly empty in practice. A
neighbour-only exchange would be cheaper but would be wrong for a step long
enough to cross a whole subdomain, and the time step is adaptive.

### Injection is deterministic and decomposition-independent

The 2D code seeds with unseeded `rand()` on every rank
(`2D-mpi/src/particletracing.c:620-635`), so the particles differ between runs,
and each rank draws its own sequence, so they differ between rank counts too.
Every rank then filters the same-indexed batch against its own extents, which
only works because all ranks happen to draw the same sequence from the same
default seed — an accident of the C library, not a property.

Seed positions are instead derived from the global particle index and the batch
number by a small deterministic function, so a given particle has the same seed
position on every rank and in every run. Ranks still filter by their own
extents; the difference is that the sequence is now defined rather than
inherited from `rand()`. This is what makes "same particles on different rank
counts" checkable.

### Particles are written as VTK polydata

The project's entire visualization story is legacy VTK opened in ParaView
(`src/vtkWriter-seq.c`, `src/vtkWriter-mpi.c`). Particles are written the same
way, as `DATASET POLYDATA` with one vertex per particle, so a particle file and
the field output for the same time load together without a custom reader. The 2D
code writes a bare list of coordinates to `.dat`, which needs one.

Output is gathered to rank 0 and written there, matching how the serial VTK path
already collects fields (`src/comm.c:327-477`). Particle counts are small enough
that a collective write buys nothing.

### The verification harness has to be built before anything else

There is nothing to compare against today: no test directory, no check driver,
no reference output, and VTK as the only output format. Every acceptance
criterion in the three specs needs something that can run it.

So the change starts by adding a `TEST` build path for check drivers, a raw
field dump of `p, u, v, w` that two runs can be diffed in, a field comparison
utility, and recorded baselines for `dcavity.par` and `canal.par` under each
solver. The baselines are recorded _before_ the defect fixes, so that the fixes'
effect on the obstacle-free setups is visible rather than assumed, and
re-recorded after step 2 as the reference the geometry work must not disturb.

## Risks / Trade-offs

- **A 3D geometry volume is 64 bytes per cell at the minimum sampling density,
  so generating and versioning benchmark geometry is a real offline cost** →
  Ship a generator and a driver script rather than only the files, so a volume
  can be regenerated rather than stored at every resolution; keep the analytic
  producer for the cases where a raster buys nothing.
- **Slab reading is more code than slurping, and an off-by-one in the offset
  arithmetic would produce a subtly shifted body rather than an error** → The
  rank-count bit-identity check is what catches it, and it is a specified
  acceptance criterion rather than an incidental test.
- **Geometric coarsening degrades multigrid convergence for under-resolved
  features** → Report the level at which a feature is lost; make the retained
  convergence rate an explicit acceptance criterion; multigrid-preconditioned CG
  later covers the residual weakness.
- **The third relaxation pass changes the SOR splitting and may shift its
  convergence rate** → Measure iteration counts against the obstacle-free case
  as an acceptance criterion rather than assuming they match.
- **The compressed solver's colour-split layout is a second place the cut-cell
  pass can be wrong, and wrong differently** → It is held to producing the same
  converged field as the other two solvers, which is a specified requirement, not
  an implementation detail left to inspection.
- **The self-maintaining solid interior depends on solid pressure being exactly
  zero** → It holds because relaxing a zero neighbourhood with a zero right-hand
  side is exactly zero in floating point, but it is a property worth asserting in
  a check driver rather than trusting; any code path that writes a nonzero value
  into a solid cell breaks it silently.
- **Volume orientation and aspect ratio are easy to get wrong and fail quietly**
  → The axis mapping is specified explicitly and gets a dedicated check with an
  asymmetric volume; an aspect mismatch produces a warning naming both ratios.
- **The defect fixes change obstacle-free results, so there is no "unchanged"
  baseline to hold the geometry work against until step 2 completes** → Baselines
  are recorded twice, before and after step 2, and the geometry work is held
  against the second.
- **Four extra `double` arrays per grid** → They are not read in the bulk sweep,
  so they do not enter the streaming working set.
- **Benchmark inputs become binary assets** → The volume is versioned with the
  code and its checksum is printed in the run header, so a silently edited file
  cannot be mistaken for the reference case.
- **Particle tracing is a second feature in one change, and it is the part with
  no benchmark riding on it** → It is confined to its own translation unit and is
  inert unless a setup configures it, so an obstacle-free or untraced run is
  bit-identical either way; that is a specified acceptance criterion rather than
  an assumption.
- **Path-based obstruction is more work per particle than a destination-cell
  test, and a face-walk has more edge cases — a particle starting exactly on a
  face, or a step crossing three faces at once** → The count is small relative to
  the flow solve, so the cost is not the concern; the edge cases are covered by
  the thin-plate scenario, which is precisely the case the cheaper test gets
  wrong.
- **Deterministic injection removes the only source of variety between runs, so
  a tracing bug that depends on particle placement will reproduce rather than
  show up intermittently** → That is the intended trade: reproducibility is worth
  more here than coverage-by-accident, and the seed sequence is a function of the
  configuration, so varying it is a parameter change.

## Migration Plan

Order matters, because the later steps cannot be validated on top of the earlier
defects.

1. **Verification harness.** `TEST` build path, check drivers, raw field dump,
   comparison utility, recorded baselines for `dcavity.par` and `canal.par`
   under all three solvers. Nothing about the solver changes.
2. **Correctness fixes that are independent of geometry.** Residual reset;
   global-checkerboard colouring; the shared pressure boundary routine; the
   multigrid transfer, stride, mesh-size, boundary-level, coarse-error and
   V-cycle-loop fixes; conditional null-space handling. Validate against step 1
   and re-record the baselines.
3. **Aperture data model and geometry producers.** Introduce `Ax`, `Ay`, `Az`,
   `Lambda`, the analytic producers and the voxel reader, plus the validation
   rules. Verify against an obstacle-free volume, which must reproduce step 2
   exactly.
4. **Discretization against apertures.** Momentum predictor, velocity
   correction, operator assembly, identity rows, fluid-only residual norms.
5. **Kernel split and solver updates.** Bulk sweep, surface list in both index
   forms, per-cell diagonal, cut-cell pass in all three solvers, multigrid
   aperture coarsening.
6. **Particle tracing.** The tracer, its parameters and its output, with
   obstruction against the apertures from step 3. It depends on step 5 only for
   having a converged flow field to trace, so it can be built alongside it.
7. **Setup migration.** Ship the generator, the volumes, and `karman.par`,
   `backstep.par` and `schaefer-turek.par`; document the format, the resolution
   rule and the particle parameters.

Rollback is per-step: steps 1, 3 and 6 are additive and independently
revertible — step 6 in particular touches nothing the flow solve reads. Step 2
changes obstacle-free results by design and is the one that has to be reverted
as a unit.

## Open Questions

- Which checksum to print for the geometry volume. Any cheap non-cryptographic
  checksum satisfies the requirement; the choice does not affect specs, approach
  or tasks.
- Whether the analytic producer is selected through a dedicated parameter key or
  through a reserved `geometryFile` value. Cosmetic, decided during
  implementation.
- Whether apertures should later be stored as `float` to halve their footprint.
  Only worth revisiting if the surface pass or the coarse levels turn out to be
  memory-limited, which the O(surface area) cost makes unlikely.
