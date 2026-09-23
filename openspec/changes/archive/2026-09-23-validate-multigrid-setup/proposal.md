# Proposal

## Why

Setup files are the interface users edit. A value that cannot produce a correct
solve is accepted today, and the run either crashes or — worse — finishes,
reports success and writes a field that is wrong.

Measured on a 16x16x16 enclosed cavity, `eps 1e-3`, `itermax 1000`:

```
  levels 0         SEGFAULT (exit 139)     also under SOLVER=cg with precon mg
  levels -1        SEGFAULT (exit 139)
  presmooth 0      exit 0   1000 cycles, residual 6.81e+00
  smoothOmega 0    exit 0   1000 cycles, residual 7.01e+00
  smoothOmega 2.0  exit 0   1000 cycles, residual 6.26e-02
  smoothOmega -1   exit 0      7 cycles, residual nan
  smoothOmega 3.0  exit 0     13 cycles, residual nan
  omg 3.0 (rb)     exit 0    552 iterations, residual nan
```

Two separate defects produce that table.

**A setup that cannot work is accepted.** `levels 0` leaves the hierarchy with no
levels allocated and the cycle indexes past the end of it. The existing clamp is
one-sided: it catches a depth too large and reports it, and does nothing about a
depth too small. `presmooth 0` builds a cycle with no smoother, which cannot
converge no matter how many cycles it is given. `smoothOmega` outside `(0, 2)`
is outside the range where the relaxation it scales converges at all, and `omg`
has the same bound and the same gap.

**A diverged solve is reported as a converged one.** Every solver stops on
`while (res >= eps*eps)`, and `NaN >= x` is false for any `x`, so the moment the
residual becomes NaN the loop exits, the solver prints `took 552 iterations to
reach nan`, returns success, and the NaN field is written to VTK. The exit status
is 0. Nothing downstream can tell that run apart from a good one.

The second is why the first matters. A bad setup value would be a nuisance if it
reliably crashed; it is a correctness problem because the run completes and
claims to have worked.

## What Changes

- **A setup that cannot produce a correct solve is refused at initialization**,
  naming the parameter and the value, before any time step runs. This follows
  the pattern the project already uses for an unsupported boundary code and an
  unknown preconditioner. It covers `levels`, `presmooth`/`postsmooth`,
  `smoothOmega` and `omg`.
- **A solve reports how it stopped** — converged, stopped at the iteration
  limit, or diverged — in three distinguishable forms, in every solver. A
  divergence stops the run with a failing status and produces no output. Hitting
  the iteration limit is reported as such but stays non-fatal: a bounded budget
  is a legitimate configuration, and `tests/checks/nullspace.c` deliberately runs
  a non-converging solve and needs its iterate back.
- **The multigrid tolerance requirement says what reaching the tolerance means**
  — a finite residual below it — so that the NaN case is excluded by the
  requirement and not only by the implementation.

### Scope note

`omg` is a relaxation-solver parameter rather than a multigrid one, so it sits
slightly outside this change's name. It is included deliberately: `omg` and
`smoothOmega` are the same quantity under the same `(0, 2)` bound, and refusing
one while accepting `omg 3.0` — which is the very value that produced the NaN
above — would be an arbitrary line.

### Non-goals

- **Validating the rest of the setup.** `eps`, `itermax`, `dt`, `tau`, the grid
  extents and the geometry parameters are not swept here. The two requirements
  added are about solves that cannot converge; a general parameter-validation
  pass is a larger and different piece of work.
- **Changing what any valid setup computes.** Every refusal is for a value that
  cannot produce a correct solve today. No converging setup changes its result,
  its iteration count or its baselines.
- **Detecting slow convergence or stagnation.** Exhausting `itermax` is reported
  honestly, but deciding that a converging-but-slow solve should be abandoned
  early is a policy question this does not take.
- **Recovering from divergence.** Falling back to a smaller relaxation factor, or
  restarting, is out of scope. The solve fails and says why.

## Capabilities

### New Capabilities

None. This makes an existing solver honest about input it cannot use and output
it did not achieve.

### Modified Capabilities

- `pressure-solvers`: adds a requirement that a setup which cannot produce a
  correct solve is refused at initialization rather than run, and a requirement
  that a solve reports how it stopped, distinguishing convergence from a
  divergence and from an exhausted iteration budget — which the current
  `while (res >= eps*eps)` test cannot do, since it reads a non-finite residual
  as the tolerance having been met. Amends the multigrid tolerance requirement so
  that "reaches that tolerance" means a finite residual below it.

## Impact

- **Modified**: `src/multigrid.c` (the guards on `levels`, the smoothing counts
  and `smoothOmega`, beside the existing equal-count refusal); `src/solverbase.c`
  or `src/parameter.c` (the `omg` bound, wherever the shared solver setup
  validates); `src/solver-mg.c:87`, `src/solver-rb.c:78`,
  `src/solver-rbc.c:197` and `:283`, `src/solver-cg.c:468` (the convergence
  test and what each reports).
- **Tests**: `tests/check-rejects.sh` gains the refused values; a diverging setup
  that must now fail rather than report success needs a home there too.
- **Baselines**: unaffected. Every shipped setup is inside the accepted range,
  and no converging solve changes.
- **Sequencing**: independent of `optimize-multigrid`, which is complete. It
  touches the same file and is easier to land after it.
- **No new dependencies.**
