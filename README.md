# NuSiF CFD Solver

The NuSiF-Solver implements a 3D structured incompressible Navier-Stokes solver.
This solver uses finite difference discretization on a staggered grid as
described in [this book](https://epubs.siam.org/doi/10.1137/1.9780898719703) by
Michael Griebel. The program supports sequential and MPI IO output of the
results in VTK file format, which can be visualized using the ParaView
application.

## Overview

The solver uses Chorin's projection method, which is explicit in the velocities
and implicit in the pressure. Physical quantities are stored on a staggered
grid: pressure at cell centers and velocities at cell faces. The time derivative
is discretized with Euler's method. Spatial derivatives use central differences;
the Donor cell differencing scheme is used for convective terms.

### Supported Tool Chains

- GCC
- CLANG
- ICX

### Parallelism

- Sequential
- MPI (derived datatypes, Cartesian topology, neighborhood collectives, MPI IO)
- OpenMP

### Pressure Solvers

- Red-Black SOR
- Compressed Red-Black SOR
- Geometric Multigrid
- Preconditioned Conjugate Gradient

All four solve the same system, including when the domain contains an
obstacle, and agree to the solve tolerance. Select one at build time with
`SOLVER=` (see [Configure](#1-configure)); the binary carries exactly one.

The conjugate gradient variant applies the pressure operator matrix-free once
per iteration, in the same two passes the relaxation sweeps use: a bulk pass
that reads no geometry at all, then a correction confined to the cells the
obstacle cuts. Its inner products, its preconditioner and its null-space
projection are all restricted to the fluid unknowns, and solid cells hold
exactly zero at every iteration rather than only at convergence.

It is the one solver that costs global communication per iteration: two
all-reductions, one for the direction's energy and one fusing the residual norm
with the preconditioned inner product. The relaxation solvers have a single
reduction, and only for their stopping test. That difference is real and is part
of what a Krylov solver is; it is measured in its own profiler regions
(`CG_DOT`, `CG_AXPY`, `PRECON`) rather than lumped into the sweep.

The preconditioner is chosen with a `precon` line in the parameter file:

| `precon` | meaning                                                          |
|----------|------------------------------------------------------------------|
| `jacobi` | the operator's diagonal, the default. Built from each cell's open faces, so it needs no stored field and reads no geometry array. |
| `none`   | no preconditioning, for comparison.                              |

An unrecognised value is refused at initialization rather than falling back to
a default. `mg` is named in that message because a multigrid preconditioner is
the intended next step and needs a symmetric V-cycle, which the current cycle is
not; a parameter file asking for it is refused until that exists.

Diagonal preconditioning does not fix the Poisson condition number, so on these
grids the iteration count still grows with resolution and the conjugate
gradient solver is not necessarily the fastest in wall-clock terms. What it
provides is a Krylov framework with the operator, the reductions and the
null-space handling in place and checked.

### Obstacles

The domain may contain an embedded body. Geometry is carried by face apertures
co-located with the velocities and a volume fraction co-located with the
pressure, so the momentum predictor, the velocity correction and the pressure
operator all see the same body. See [Obstacle geometry](#obstacle-geometry).

## Build

### 1. Configure

Copy or edit `config.mk` to select the tool chain and options:

```make
# Supported: GCC, CLANG, ICX
TOOLCHAIN ?= GCC
# Supported: true, false
ENABLE_MPI ?= true
ENABLE_OPENMP ?= false
# Supported: rb, rbc, mg, cg
SOLVER ?= rb
# Supported: seq, mpi
VTK_OUTPUT_FMT ?= seq

OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE
#OPTIONS +=  -DVERBOSE_AFFINITY
#OPTIONS +=  -DVERBOSE_DATASIZE
#OPTIONS +=  -DVERBOSE_TIMER
```

The verbosity options enable detailed output about the solver, affinity
settings, allocation sizes, and timer resolution. For debugging, enable
`-DVERBOSE`.

### 2. Build

```sh
make
```

Multiple tool chains can coexist in the same directory. Intermediate build
results are stored in `./build/<TOOLCHAIN>/`. The executable is named
`NusifSolver-<TOOLCHAIN>`.

To see all executed commands:

```sh
make Q=
```

### 3. Clean

Remove intermediate build results for the active tool chain:

```sh
make clean
```

Remove all build results for all tool chains, including data and visualization
output:

```sh
make distclean
```

### 4. (Optional) Generate Assembler

```sh
make asm
```

Assembler files are placed in `./build/<TOOLCHAIN>/`.

## CLANG Tooling Support

The Makefile generates a `.clangd` configuration with the correct compiler
flags for the clang language server. This requires GNU Make 4.0 or newer.
Note: the default Make version on macOS is 3.81 — install a newer version via
[Homebrew](https://brew.sh). Without GNU Make 4.0+, LSP support will be
restricted to previously opened buffers.

Alternatively, use [Bear](https://github.com/rizsotto/Bear) to generate a
`compile_commands.json` compilation database, which also enables jump-to-definition
without a previously opened buffer:

```sh
bear -- make
```

The repository includes `.clang-format` and `.clang-tidy` files to enforce
consistent formatting and naming conventions. To reformat all source files:

```sh
make format
```

This requires `clang-format` in your `PATH`.

## Usage

Provide a parameter file describing the problem to solve:

```sh
./NusifSolver-CLANG dcavity.par
```

Five example setups are included:

| Setup | Flow | Geometry |
|---|---|---|
| `dcavity.par` | lid-driven cavity | none |
| `canal.par` | empty canal flow | none |
| `karman.par` | cylinder in a channel | voxel volume |
| `backstep.par` | backward-facing step | voxel volume |
| `schaefer-turek.par` | the published 3D cylinder benchmark, case 3D-2Z | analytic cylinder |

The two setups that reference a voxel volume need it generated first:

```sh
tools/make-geometry.sh
```

To plot the pressure solver residual as a function of iteration:

```sh
make plot
```

## Obstacle geometry

A setup places a body in the domain with a `geometryFile` entry. Without one the
domain is obstacle-free and the solver behaves exactly as it did before
obstacles existed.

### Analytic bodies

The value may name a body directly, which is not limited by any sampling
resolution and is what the reference benchmarks use:

```
geometryFile  sphere:xc,yc,zc,r
geometryFile  cylinder-z:xc,yc,r          # also cylinder-x and cylinder-y
geometryFile  box:x0,y0,z0,x1,y1,z1
geometryFile  plate-z:z,x0,y0,x1,y1       # also plate-x and plate-y
```

A plate has no volume: it closes a plane of faces and leaves the cells on both
sides fluid. That is a body a per-cell obstacle type cannot express at all.

### Voxel volumes

Otherwise the value is a path to a voxel volume, the three-dimensional
counterpart of the binary PGM the two-dimensional solver reads:

```
P5V
# optional comment lines, anywhere in the header
<nx> <ny> <nz>
255
<nx*ny*nz raw bytes, x fastest, then y, then z>
```

A voxel below 128 is solid, 128 or above is fluid. Voxel `(vx, vy, vz)` covers
the box `[vx, vx+1) * xlength / nx` and likewise in y and z, so the index order
matches the axis order and the origin is the domain origin. No axis is flipped.

Write one with the generator:

```sh
tools/genvox.py sphere --out geometry/sphere.vox \
    --size 256 256 256 --domain 4 4 4 --center 2 2 2 --radius 0.5

tools/genvox.py cylinder --out geometry/karman.vox \
    --size 810 216 216 --domain 30 8 8 --axis z --center 5 4 --radius 1

tools/genvox.py box --out geometry/backstep.vox \
    --size 560 120 120 --domain 7 1.5 1.5 --corner 0 0 0 --extent 1 0.5 1.5
```

The volumes the shipped setups use are generated rather than committed: at the
required sampling density a volume is 64 bytes per grid cell, so the karman one
is 32 MB. `tools/make-geometry.sh` rebuilds them in a few seconds.

### Rules the solver enforces

A volume must provide **at least 4 voxels per grid cell in every direction**.
Below that the represented body changes as the grid is refined, which would make
a resolution sweep measure something other than convergence. An under-resolved
volume is refused rather than sampled.

The **fluid region must be connected** under face connectivity. Each sealed
pocket contributes an independent constant to the pressure null space, which the
solvers do not carry, so geometry that encloses one is refused.

A volume whose aspect ratio disagrees with the domain produces a warning and is
then sampled anisotropically. A volume that cannot be opened, is not a P5V file,
or is shorter than its header declares is refused. Every refusal happens at
initialization, before any time step runs, and names the file and the reason.

The run header records which geometry is in use, the voxel dimensions and a
checksum of the file, so a silently edited volume cannot be mistaken for the
reference case.

## Particle tracing

A setup can release massless tracer particles, so that flow around a body shows
up as streaklines. It is off unless asked for: a setup with no particle block
produces exactly the fields and files it would with no tracer present.

```
numberOfParticles   400     # released per batch; 0 or absent means no tracing
startTime           20.0    # when the first batch is released
injectTimePeriod    1.0     # between batches
writeTimePeriod     0.5     # between output files

x1  0.2                     # the seed box, one corner
y1  0.5
z1  0.5
x2  0.4                     # and the other
y2  7.5
z2  7.5
```

Positions are written to `vis_files/particles_NNNNN.vtk` as VTK polydata, which
opens in ParaView alongside the field output.

Seed positions come from a hash of the particle's index and batch number rather
than from a random number generator, so the same setup releases the same
particles in every run and at every rank count, with no communication needed to
arrange it.

A particle is stopped by the **apertures of the faces its path crosses**, not by
the volume fraction of the cell it lands in. That distinction is the point: a
body one face thick has fluid on both sides, so a destination-cell test sees
nothing in the way and the particle passes straight through it. Particles that
meet a body, or leave through a domain boundary, are removed and counted; the
run reports how many of each at the end.

## Tests

```sh
tests/run-all.sh
```

This builds each solver variant in turn and runs everything: the check drivers
at one and several ranks, the inputs the solver has to refuse, geometry against
rank count, the solvers against an obstacle, and the recorded field baselines.

The pieces can also be run on their own:

| Command | What it covers |
|---|---|
| `tests/run-checks.sh [-n RANKS]` | the check drivers in `tests/checks/` |
| `tests/check-setups.sh` | every shipped setup starts and its body reaches the flow |
| `tests/check-rejects.sh` | inputs that must be refused |
| `tests/check-geometry-ranks.sh` | apertures are identical whatever the rank count |
| `tests/check-solver-obstacle.sh` | the three solvers agree with a body present |
| `tests/check-schaefer-turek.sh` | drag, lift and Strouhal against their published ranges |
| `tests/check-particles.sh` | particle output, and that it does not depend on the rank count |
| `tests/record-baseline.sh [-v]` | record or verify the field baselines |

Baselines and test geometry are generated, not committed; `record-baseline.sh`
and `tests/make-geom.sh` rebuild them.

The check drivers link against the solver objects, so they need the test build:

```sh
make tests
```

which also produces `NusifSolver-<TOOLCHAIN>-test`, a solver that writes a raw
dump of `p`, `u`, `v` and `w` when `NUSIF_FIELD_DUMP` names a path. `tools/fieldcmp`
compares two dumps -- `--l2` judges by the L2 difference rather than the largest
one, which is the honest measure when two runs differ only in the order their
global sums were formed -- and `tools/fieldprobe.py` reports statistics over a
box of cells.

### Forces on a body

A run whose domain contains an obstacle writes `forces.dat`, a time series of
the force the flow exerts on it, integrated over the closed faces that make up
its surface. `tools/stcoeffs.py` turns that into the drag and lift coefficients
and the Strouhal number of the Schaefer-Turek benchmark and reports each against
its published range:

```sh
./NusifSolver-CLANG schaefer-turek.par
tools/stcoeffs.py forces.dat
```

That check reports rather than asserts. Binary apertures give a body a staircase
surface and a first-order wall treatment, so the drag is expected to come out
high and a wake resolved by ten cells across the cylinder may not shed at all.
`STRICT=1 tests/check-schaefer-turek.sh` turns the report into a pass-or-fail
check, which is the measurement a later fractional-aperture change has to
satisfy.
