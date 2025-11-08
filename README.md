# NuSiF CFD Solver

The NuSiF-Solver implements a 3D structured incompressible Navier Stokes solver.
This solver uses finite difference discretization using a staggered grid as
described in [this book](https://epubs.siam.org/doi/10.1137/1.9780898719703) by
Michael Griebel. The program supports sequential and MPI IO output of the
results in VTK file format, that can be visualized using the Paraview
application.

## Infos

### Preconfigured tool chains

- GCC
- CLANG
- ICX

### Parallelism

- MPI
- OpenMP

### Solvers

- Red-Black SOR
- Geometric Multigrid

## Build

1. Configure the tool chain and additional options in `config.mk`:

```make
# Supported: GCC, CLANG, ICX
TOOLCHAIN ?= GCC
# Supported: true, false
ENABLE_MPI ?= true
ENABLE_OPENMP ?= false
# Supported: rb, mg
SOLVER ?= rb
# Supported: seq, mpi
VTK_OUTPUT_FMT ?= seq

ENABLE_OPENMP ?= false

OPTIONS +=  -DARRAY_ALIGNMENT=64
#OPTIONS +=  -DVERBOSE
#OPTIONS +=  -DVERBOSE_AFFINITY
#OPTIONS +=  -DVERBOSE_DATASIZE
#OPTIONS +=  -DVERBOSE_TIMER
```

The verbosity options enable detailed output about solver, affinity settings, allocation sizes and timer resolution.
For debugging you may want to set the VERBOSE option:

```make
# Supported: GCC, CLANG, ICX
TAG ?= GCC
# Supported: true, false
ENABLE_MPI ?= true
ENABLE_OPENMP ?= false
# Supported: rb, mg
SOLVER ?= rb
# Supported: seq, mpi
VTK_OUTPUT_FMT ?= seq

ENABLE_OPENMP ?= false

OPTIONS +=  -DARRAY_ALIGNMENT=64
OPTIONS +=  -DVERBOSE
#OPTIONS +=  -DVERBOSE_AFFINITY
#OPTIONS +=  -DVERBOSE_DATASIZE
#OPTIONS +=  -DVERBOSE_TIMER
```

1. Build with:

```sh
make
```

You can build multiple tool chains in the same directory, but notice that the
Makefile is only acting on the one currently set. Intermediate build results are
located in the `./build/<TOOLCHAIN>` directory. The executable is named
`NusifSolver-<TOOLCHAIN>`.

To output all executed commands use:

```sh
make Q=
```

1. Clean up intermediate build results of active tool chain use:

```sh
make clean
```

To clean all build results of all tool chains including all data and
visualization output use:

```sh
make distclean
```

1. (Optional) Generate assembler for all source files:

```sh
make asm
```

The assembler files will also be located in the `./build/<TOOLCHAIN>` directory.

## CLANG tooling support

The Makefile will generate a .clangd configuration to correctly set all options
for the clang language server. This is only important if you use an editor with
LSP support and want to edit or explore the source code. It is required to use
GNU Make 4.0 or newer. While building with older make versions will work, the
generation of the .clangd configuration for the clang language server will not
work. The default Make version included in MacOS is 3.81! Newer make versions
can be easily installed on MacOS using the Homebrew package manager. This will
restrict LSP support for all loaded compilation units.

An alternative is to use [Bear](https://github.com/rizsotto/Bear), a tool that
generates a compilation database for clang tooling. This method also will enable
to jump to any definition without a previously opened buffer. You have to build
one time with Bear as a wrapper:

```sh
bear -- make
```

The repository includes `.clang-format` and `.clang-tidy` files to enforce
consistent formatting and variable naming.

To reformat all source files use:

```sh
make format
```

This required `clang-format` in your `PATH`.

## Usage

You have to provide a parameter file describing the problem you want to solve:

```

./NusifSolver-CLANG dcavity.par

```

Example test cases are given in the `dcavity.par` (a lid driven cavity test
case) and `canal.par` (simulating a empty canal) files.

You can plot the resulting residual as a function of iterations using:

```sh
make plot
```
