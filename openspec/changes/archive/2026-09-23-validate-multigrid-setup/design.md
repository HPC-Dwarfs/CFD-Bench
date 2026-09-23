# Design

## Context

See proposal.md — Why. What shapes the approach:

- The project already has a refusal idiom: `pressureBcInit` refuses an
  unsupported boundary code, `solver-cg.c` refuses an unknown preconditioner,
  and `multigridBuild` refuses unequal smoothing counts. All three print to
  `stderr` on the master rank and `exit(EXIT_FAILURE)` before any time step.
  `tests/check-rejects.sh` tests them by requiring a nonzero exit, a message
  naming the problem, and no `TIME` line in the output.
- The guards fall in two places. `levels`, `presmooth`/`postsmooth` and
  `smoothOmega` reach `multigridBuild` through `MultigridSpecType`, next to the
  refusal already there. `omg` is a relaxation-solver parameter and reaches
  `solverBaseInit` as `s->omega`, which is the only place every solver variant
  passes through.
- All four solver loops have the same shape:
  `while ((res >= epssq) && (it < itermax))` — `src/solver-mg.c:87`,
  `src/solver-rb.c:78`, `src/solver-rbc.c:197` and `:283`,
  `src/solver-cg.c:468`. The `rbc` variant has two because it carries a
  compressed and a natural layout.
- `main.c:92` calls `solve()` once per step and uses the returned residual only
  for `writeResidual`. Nothing inspects it for success. Output is written after
  the loop.
- **`tests/checks/nullspace.c:154` calls `solve()` directly on a deliberately
  incompatible right-hand side, runs it to `itermax`, and compares the result
  against the handled case.** It needs a non-converging solve to return its
  iterate. This is the constraint that separates the two non-converged outcomes:
  divergence may abort, exhausting the budget may not.

## Goals / Non-Goals

**Goals:**

- No setup value that cannot produce a correct solve reaches a time step.
- No run that diverged can be mistaken for one that converged, by a person
  reading the output or by a script reading the exit status.
- One place that decides what to do about a failed solve, rather than four.

**Non-Goals:**

- Changing any converging solve's result, iteration count, message or baseline.
- A general parameter-validation pass over the whole setup file.
- Deciding when a slow-but-converging solve should be given up on.

## Decisions

### 1. Refusals go where the parameter is already consumed, not in a new validator

**Chosen:** `levels`, the smoothing counts and `smoothOmega` are checked in
`multigridBuild` beside the existing equal-count refusal; `omg` is checked in
`solverBaseInit`.

A single `validateParameters()` over the whole `Parameter` struct was the
alternative. Rejected for now: most of the struct has no validation and writing
one function that checks four fields invites it to grow into the general pass
this change lists as a non-goal, without the requirement that would say what
belongs in it. Putting each guard where the value is consumed also keeps the
message next to the code that explains why the bound exists.

The cost is that `multigridBuild` is reached only by a build that constructs a
hierarchy, so `SOLVER=rb` never sees the `levels` guard. That is the same
property the existing equal-count refusal has, and `check-rejects.sh` already
builds a dedicated multigrid binary to test it.

### 2. `levels < 1` is refused rather than clamped up to 1

The depth requirement says a request larger than the decomposition supports is
clamped and reported, so clamping `0` up to `1` would be symmetric. Rejected:
the two are not the same case. Asking for more depth than the grid allows is a
reasonable thing to write in a setup — the author cannot always know the
decomposition — and the solver can honour it approximately and say so. Asking
for zero or fewer levels is not a request the solver can approximate; it is a
value with no meaning, and silently substituting a different one hides a typo in
a file the user will read back later.

### 3. Divergence is detected in the loop and acted on by the caller

**Chosen:** each solver's loop breaks when the residual stops being finite and
reports it; `main.c` treats a non-finite return from `solve()` as fatal, before
`adaptUV` and before any output.

Detection has to be in the loop because that is where the residual is, and
because the loop must stop — continuing to iterate on NaNs wastes the whole
iteration budget producing nothing. The policy — stop the run, write nothing —
belongs in `main.c` because it is the same policy for every solver, and because
the check drivers are a second caller that wants a different one.

That second caller is the reason not to `exit()` from inside `solve()`, which
would be a smaller change: `tests/checks/nullspace.c` runs a non-converging
solve on purpose. It is not a *diverging* one, so an `exit()` on NaN alone would
not break it today — but a solver that aborts the process is a solver the check
drivers cannot use to study failure, and they are the natural home for tests of
exactly this behaviour.

Alternative considered: keeping the loop condition and letting NaN fall out of
it as now, then testing `isfinite` once afterwards. It gives the same answer for
the NaN case and is a one-line change per solver. Rejected because it leaves the
loop condition itself wrong — a reader still has to know that `NaN >= x` is
false to understand why the loop terminates — and because an infinite residual
that is not NaN would keep iterating to the limit.

### 4. Exhausting the iteration limit stays non-fatal

Required by the constraint in Context: `nullspace.c` needs the iterate from a
solve that ran out of budget. It is also a legitimate configuration — a bounded
budget is a normal thing to ask for.

What changes is only the reporting. Today every solve prints the same
`took N iterations to reach X` whether it converged or not, so the two are
distinguished only by comparing `X` against the setup's `eps` by eye. Each
solver reports the three outcomes distinguishably instead.

Whether `main.c` should *also* stop on an exhausted budget is a policy question
this change does not take, and the non-goals say so. Divergence is the case
where continuing is meaningless; a budget-limited solve produces an iterate that
a time step can legitimately use.

### 5. The three outcomes are named in the message, not encoded in the return

`solve()` keeps returning the residual. The outcome is distinguishable from it —
non-finite means diverged — and `main.c` needs no new type to act on that.

Returning a status enum alongside was considered and rejected as a larger
interface change across four variants and their check drivers for one bit that
the residual already carries.

## Risks / Trade-offs

- **A setup that currently runs starts being refused.** → Every refusal is for a
  value that cannot produce a correct solve, so anything refused was already
  producing a wrong answer or a crash. `check-setups.sh` runs every shipped
  setup and is what says none of them is affected; a task exists for it.

- **A run that currently exits 0 starts exiting nonzero.** This is the point,
  but it could surface in a harness that does not expect it. → Only a diverged
  run changes status, and a diverged run's output was NaN. Anything relying on
  that exit code was relying on a wrong answer.

- **`isfinite` on a residual computed by a global reduction.** A NaN on one rank
  propagates through `MPI_Allreduce` to every rank, so the test is consistent
  across ranks without extra communication. Worth a check at more than one rank
  rather than assuming.

- **The `omg` guard sits in `solverBaseInit`, which multigrid also calls**, so a
  multigrid build validates both factors. That is correct but means the
  `smoothOmega` message and the `omg` message must name their parameter clearly
  or a user will fix the wrong line in the setup file.

- **Four loops to change, and `rbc` has two of them.** → The rbc pair is the
  known trap; a task names both line numbers, and the check drivers run under
  every variant.

## Migration Plan

Ordered so each step is verifiable on its own:

1. The refusals, with `check-rejects.sh` cases and a `check-setups.sh` run to
   show no shipped setup is caught.
2. The convergence test and reporting, per solver, with the converged path
   asserted unchanged first.
3. `main.c` acts on a diverged solve.

Rollback is per step. Step 1 is independent of steps 2 and 3.

## Open Questions

None that affect the specs, the approach or the tasks. One deferred:
`tests/checks/nullspace.c` is currently the only caller that depends on a
non-converging solve returning, and it does so incidentally rather than as a
stated requirement of that driver. If a later change makes budget exhaustion
fatal, that dependency is where it will surface.
