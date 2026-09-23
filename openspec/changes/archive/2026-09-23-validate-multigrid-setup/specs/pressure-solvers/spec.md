# Spec Delta

## ADDED Requirements

### Requirement: A setup that cannot produce a correct solve is refused

A setup value that makes a correct solve impossible SHALL be refused at
initialization, naming the parameter and the value it was given, before any time
step is executed. This is the treatment an unsupported boundary code and an
unknown preconditioner already receive, and for the same reason: a run that
cannot produce a right answer is worse than no run, because its output is
indistinguishable from a good one.

This applies to a value that is outside the range in which the method works at
all, not to one that merely makes it slow. In particular:

- a multilevel hierarchy SHALL be at least one level deep;
- a cycle SHALL apply at least one smoothing sweep on each side of its coarse
  correction, since a cycle with no smoother cannot converge however many cycles
  it is given;
- a relaxation factor SHALL lie strictly between zero and two, which is the
  range in which the relaxation it scales converges. This binds every relaxation
  factor a setup carries, whether it drives a solver or a smoother.

A value that is merely larger than the problem can use SHALL NOT be refused
where the solver can adapt to it and report what it did — a requested hierarchy
depth greater than the decomposition supports is clamped and reported, and
remains so.

#### Scenario: A hierarchy depth below one is refused

- **WHEN** a setup requests fewer than one multigrid level
- **THEN** the solver aborts at initialization naming the parameter and the value, rather than building a hierarchy it cannot cycle on

#### Scenario: A cycle with no smoother is refused

- **WHEN** a setup requests zero or fewer smoothing sweeps on either side of the coarse correction
- **THEN** the solver aborts at initialization, rather than running a cycle that cannot reduce the residual

#### Scenario: A relaxation factor outside the convergent range is refused

- **WHEN** a setup gives a relaxation factor that is zero or negative, or two or greater, for either a solver or a smoother
- **THEN** the solver aborts at initialization naming the parameter and the value

#### Scenario: A depth the decomposition cannot provide is still accepted

- **WHEN** a setup requests more levels than the local extents can be halved to provide
- **THEN** the solver builds as many as it can and reports the number built, as before, because it can honour the request approximately and say so

#### Scenario: Every shipped setup is accepted

- **WHEN** each setup the project ships is run
- **THEN** none of them is refused, so the refusals describe values no correct setup uses

### Requirement: A solve reports how it stopped

A solve SHALL end in one of three distinguishable outcomes — it converged, it
stopped at the iteration limit, or it diverged — and SHALL report which. The
three are not interchangeable and SHALL NOT be reported in the same form.

Divergence SHALL be judged on the residual being a finite number. This is a
requirement about the convergence test and not only about the message: a test of
the form "continue while the residual exceeds the tolerance" treats a residual
that is not a number as having met the tolerance, because a comparison against a
NaN is false whichever way it is written. A solver whose iteration has diverged
therefore leaves that test, reports the iteration count it happened to reach, and
returns a field of NaNs as though it had converged — the one outcome a solver
must never produce, since nothing downstream can tell it from a correct run.

A diverged solve SHALL stop the run with a failing status, and no field output
SHALL be produced from it.

Stopping at the iteration limit SHALL NOT stop the run. It is a legitimate
configuration — a deliberately bounded budget, or a problem being studied
precisely because it does not converge — and the iterate SHALL still be returned
to the caller. What changes is only that it is reported as having stopped at the
limit rather than as having converged.

#### Scenario: A diverged solve fails

- **WHEN** a solve's residual becomes infinite or not a number
- **THEN** the solver reports divergence and the iteration at which it occurred, and the run exits with a failing status

#### Scenario: A field of NaNs is never written as output

- **WHEN** an iteration diverges
- **THEN** no field output is produced from that run, so a diverged run cannot be mistaken for a converged one by anything reading the results

#### Scenario: Exhausting the iteration limit is reported but not fatal

- **WHEN** a solve reaches its iteration limit with the residual still above the tolerance
- **THEN** the solver reports that it stopped at the limit and the residual it reached, distinctly from having converged, and returns the iterate so that a caller studying a non-converging problem still receives it

#### Scenario: Every solver is held to it

- **WHEN** any shipped solver diverges
- **THEN** it reports the divergence, so the guarantee does not depend on which solver variant the binary was built with

#### Scenario: A converging solve is unaffected

- **WHEN** a solve reaches its tolerance
- **THEN** it reports the same iteration count and residual it reported before, and its converged field is unchanged

## MODIFIED Requirements

### Requirement: Multigrid iterates to the requested tolerance

Multigrid SHALL repeat cycles until the configured tolerance is reached or the configured iteration limit is exhausted, and SHALL report the number of cycles used.

Reaching the tolerance SHALL mean a finite residual below it. A residual that is
not a finite number does not reach any tolerance, and SHALL be reported as a
divergence rather than as a solve that finished early.

#### Scenario: Multigrid honours the tolerance

- **WHEN** a setup is solved by multigrid with a given tolerance and iteration limit
- **THEN** the solve either reaches that tolerance or stops at the limit, and the cycle count is reported

#### Scenario: A non-finite residual is not a tolerance being met

- **WHEN** a multigrid solve's residual becomes not a number
- **THEN** the solve is reported as diverged rather than as having reached the tolerance in the cycles taken so far
