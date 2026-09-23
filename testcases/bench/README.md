# Benchmark tiers

Setups sized for scaling runs rather than for physics: three cases, three
sizes each. They exist to be timed. None carries a recorded baseline and none
of their flows is a result -- numerical agreement is the regression setups'
job, in `../regression`, and the physical cases are in `../flow`.

| case | domain | small | medium | large |
|---|---|---|---|---|
| `dcavity` | 1 x 1 x 1 | 48³ | 96³ | 192³ |
| `karman` | 32 x 8 x 8 | 192 x 48 x 48 | 384 x 96 x 96 | 768 x 192 x 192 |
| `canal` | 32 x 4 x 4 | 384 x 48 x 48 | 768 x 96 x 96 | 1536 x 192 x 192 |
| | **cells** | 0.1 / 0.4 / 0.9 M | 0.9 / 3.5 / 7.1 M | 7.1 / 28.3 / 56.6 M |
| | **targets** | up to 27 ranks | up to 216 ranks | up to 1728 ranks |

Run one like any other setup:

```sh
mpirun -n 8 ./CFD-Bench-CLANG testcases/bench/dcavity-small.par
```

## Sizing rule

Multigrid coarsens until a *local* extent can no longer be halved, so the depth
of the hierarchy depends on the grid and the decomposition together. A grid that
has seven levels on one rank can have none at all on 216. Nothing fails when
that happens -- the solve converges more slowly and reports nothing unusual --
so a benchmark on such a grid would time a smoother and call it multigrid.

Four levels is the threshold, which needs every local extent to be at least 16.
`MPI_Dims_create` divides the ranks as cubically as it can, so at the rank
count a tier targets the smallest extent is divided by its cube root:

```
  smallest extent  >=  16 * (ranks per direction)

      27 ranks  ->  3 per direction   ->  48
     216 ranks  ->  6 per direction   ->  96
    1728 ranks  ->  12 per direction  ->  192
```

The other extents follow the domain's ratio, and the domains are rounded to
exact integer ratios -- `karman` is 32 x 8 x 8 rather than `../flow/karman.par`'s
30 x 8 x 8, `canal` 32 x 4 x 4 rather than 30 x 4 x 4 -- so every cell is cubic.
Stretched cells would degrade the point smoother, and the benchmark would be
measuring that by accident.

### Why `3 * 2^k` and not a power of two

The counter-intuitive part. 48, 96 and 192 are `3 * 2^4`, `3 * 2^5` and
`3 * 2^6`, and divide evenly by 2, 3, 4, 6, 8 and 12 -- the factors
`MPI_Dims_create` actually produces. A power of two divided by 3 or 6 is odd,
and an odd local extent has no hierarchy at all:

```
  128 / 6  = 21  ->  1 level         192 / 6  = 32  ->  6 levels
  128 / 12 = 10  ->  2 levels        192 / 12 = 16  ->  4 levels
```

So a 128³ grid, which looks like the natural choice, collapses to a single
level at 216 ranks.

### Where a body binds first

A body limits depth too: coarsening past the level at which it stops being
represented computes a correction to a problem without it. The `karman`
cylinder is twelve cells across in the small tier and is lost below level 2, so
`karman-small` is built with three levels, not four, and says so. The medium
tier's body allows exactly four; at the large tier the grid binds first.

## A fixed amount of work

Every tier sets `tau 0`, which bypasses the adaptive time step, and a fixed
`dt`. The three sizes of a case share the same `dt` and `te` and run the same
200 steps, so they differ only in cells. An adaptive step would let the grid
and the flow choose the step count, and two sizes would then not be comparable.

`dt` is the one the large tier's stability bound allows with margin, and is a
power of two so that the time accumulates exactly; `te` is 199 steps of it,
because the loop runs while `t <= te`. How 200 was chosen is recorded in each
file: it amortises the first two solves, which are far more expensive than the
rest because the initial flow is not divergence-free, and stops short of the
regime where every solve is a single cycle and a longer run times the solver
less rather than more.

## What is checked

Each file states the rank counts it is meant for, the depth it keeps at every
one of them, and whether the grid or the body sets that depth:

```
# bench-ranks       1 8 27
# bench-min-levels  4
# bench-bound       grid
```

`tests/check-bench-tiers.sh` computes the depth every tier reaches at every
rank count it names and fails any that falls short. It also runs the small
tiers at 1, 8 and 27 ranks and compares the solver's reported level count with
the computed one, so the computation is checked against the implementation
rather than restating it. Adding a tier means writing those three lines and
running the check.

## First measurements

Multigrid (`SOLVER=mg`), small and medium tiers, 200 steps each, on an Apple M2
Pro (6 performance cores) with Apple clang 21, `-O3 -ffast-math` and Open MPI
5.1. Wall clock is the time-stepping loop only; cycles are summed over all 200
pressure solves. Every run converged at every step.

| tier | ranks | levels | cycles | wall clock |
|---|---|---|---|---|
| `dcavity-small` | 1 | 5 | 402 | 1.4 s |
| | 8 | 4 | 402 | 1.1 s |
| `karman-small` | 1 | 3 | 519 | 8.9 s |
| | 8 | 3 | 519 | 3.8 s |
| `canal-small` | 1 | 5 | 1206 | 26.9 s |
| | 8 | 4 | 3969 | 41.7 s |
| `dcavity-medium` | 1 | 6 | 417 | 11.9 s |
| | 8 | 5 | 417 | 5.6 s |
| `karman-medium` | 1 | 4 | 627 | 65.2 s |
| | 8 | 4 | 627 | 27.9 s |
| `canal-medium` | 1 | 6 | 1766 | 316.5 s |
| | 8 | 5 | 5872 | 423.8 s |

The one to look at is `canal`. Divided over 8 ranks it needs 3.3 times the
cycles it does on one, at both sizes, and is slower in wall clock than the
single-rank run. `dcavity` also loses a level at 8 ranks and its cycle count
does not move, so the lost level alone does not explain it. Not yet
investigated; it is the kind of thing these tiers exist to show.

The relaxation solvers and conjugate gradients are not in the table. On
`canal-medium` they take hours to days at these rank counts, which is a job for
a machine that can afford it rather than for the first set of numbers.
