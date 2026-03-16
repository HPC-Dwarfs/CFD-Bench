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

## Build

### 1. Configure

Copy or edit `config.mk` to select the tool chain and options:

```make
# Supported: GCC, CLANG, ICX
TOOLCHAIN ?= GCC
# Supported: true, false
ENABLE_MPI ?= true
ENABLE_OPENMP ?= false
# Supported: rb, rbc, mg
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

Two example test cases are included:

- `dcavity.par` — lid-driven cavity
- `canal.par` — empty canal flow

To plot the pressure solver residual as a function of iteration:

```sh
make plot
```
