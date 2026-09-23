# Spec Delta

## ADDED Requirements

### Requirement: The coarsest level is solved

The coarsest level of a multilevel hierarchy SHALL be solved, not merely
relaxed. The correction a cycle carries upward is only as good as the coarse
problem it came from, so a coarsest level that is still far from its own
solution limits the whole cycle regardless of how many levels sit above it.

The work SHALL be bounded and independent of the residual, so that a cycle used
as a preconditioner remains a fixed linear operator.

#### Scenario: The coarsest level reaches its own solution

- **WHEN** the coarsest level's problem is solved as part of a cycle
- **THEN** its residual is reduced by orders of magnitude, not by the factor a fixed handful of relaxation sweeps would achieve

#### Scenario: The coarse solve does not vary with the residual

- **WHEN** the coarsest level is solved for right-hand sides of widely differing magnitude
- **THEN** each solve performs the same operations, with no inner convergence test

#### Scenario: A deeper hierarchy does not need a stronger coarse solve to converge

- **WHEN** the same problem is solved with the hierarchy the grid supports and with one level fewer
- **THEN** both converge, and the deeper one takes no more cycles

### Requirement: A multilevel solver uses the depth available to it

A multilevel solver SHALL build the hierarchy its grid and its decomposition
allow, up to the depth the setup requests. Where the requested depth cannot be
built, the solver SHALL report the depth it actually built rather than failing,
since coarsening is limited by the local extents and therefore by how the domain
was divided.

A setup SHALL NOT be configured with a depth materially shallower than its grid
and its body support, because a shallow hierarchy leaves the coarsest problem
large and the cycle weak.

Where coarsening would stop representing an embedded body, the level at which
that happens bounds the depth worth building, since a coarse correction computed
on a domain that no longer contains the body is not a correction to the problem
being solved. A setup SHALL be configured to the shallower of the two limits,
and SHALL record which limit binds it.

#### Scenario: Requested depth exceeds what the decomposition allows

- **WHEN** a setup requests more levels than the local extents can be halved to provide
- **THEN** the solver builds as many as it can, reports the number it built, and solves to the requested tolerance

#### Scenario: Depth is reported, not assumed

- **WHEN** a multilevel solver initializes
- **THEN** the number of levels it built is visible in its output, so a shallow hierarchy is apparent rather than inferred from poor convergence

#### Scenario: A body bounds the depth before the grid does

- **WHEN** a setup contains a body that coarsening stops representing at a shallower level than the grid stops halving at
- **THEN** the setup is configured to the depth at which the body is still represented, rather than to the depth the grid alone would allow

#### Scenario: Depth differing with the decomposition is reported rather than compared

- **WHEN** the same setup is solved on two rank counts whose local extents allow different depths
- **THEN** each run reports the depth it built, their converged fields agree, and their iteration counts are not required to match, because two hierarchies of different depth are two different iterations

## MODIFIED Requirements

### Requirement: All pressure solvers solve the same system

Every pressure solver the project ships SHALL solve the same embedded-boundary pressure system, with the same treatment of closed faces and solid cells. Given the same setup and a sufficiently tight tolerance, the solvers SHALL produce the same pressure field.

Agreement between solvers SHALL be judged between converged fields. Each solver SHALL be converged to a tolerance tighter than the one their agreement is judged at, because a solver's tolerance bounds a mean square residual over the domain and does not bound the largest pointwise difference between two fields. Comparing two solvers at the tolerance they each stopped at measures how they stop as much as what they converge to, and a faster iteration that takes larger steps can cross that threshold from further away while being no less correct.

#### Scenario: Solvers agree

- **WHEN** the same setup with an obstacle is solved by each shipped solver, each converged well past the tolerance the comparison is made at
- **THEN** the resulting pressure fields agree with one another to within the comparison tolerance

#### Scenario: A faster solver is not penalised for stopping sooner

- **WHEN** two solvers converge to the same problem, one taking markedly fewer or larger steps than the other
- **THEN** their converged fields agree, and neither is reported as disagreeing on account of where its iteration crossed the stopping threshold

#### Scenario: Obstacle is visible to every solver

- **WHEN** an obstacle is introduced into a setup
- **THEN** every solver produces a pressure field whose gradient across zero-aperture faces is zero to the solve tolerance

#### Scenario: A solver using a different internal data layout gives the same answer

- **WHEN** a solver stores the unknowns in a layout other than the domain's natural one in order to vectorize its sweep
- **THEN** its converged field with an obstacle present matches the other solvers' to the solve tolerance, so the layout is not observable in the result

### Requirement: A Krylov pressure solver is shipped

The project SHALL ship a conjugate gradient pressure solver as a first-class
solver variant, selected the same way as the relaxation solvers, and subject to
the same requirement that every shipped solver solves the same embedded-boundary
system and agrees with the others once converged.

#### Scenario: Krylov solver agrees with the relaxation solvers

- **WHEN** a setup containing an obstacle is solved by the Krylov solver and by each relaxation solver, each converged well past the tolerance the comparison is made at
- **THEN** all of the resulting pressure fields agree with one another to within the comparison tolerance

#### Scenario: Krylov solver honours the tolerance

- **WHEN** a setup is solved by the Krylov solver with a given tolerance and iteration limit
- **THEN** the solve either reaches that tolerance or stops at the limit, and the iteration count is reported in the same form as the other solvers report theirs

#### Scenario: The reported residual is comparable across solvers

- **WHEN** the Krylov solver reports a residual norm for a converged solve
- **THEN** that norm is the same quantity the other solvers report — the mean square residual over the fluid unknowns — so the same tolerance means the same thing for every solver

### Requirement: The multilevel cycle is a symmetric operator

A multilevel cycle used to precondition a Krylov solver SHALL be a symmetric
linear operator on the fluid unknowns, so that the method preconditioned by it
remains the method it claims to be. Every part of that cycle SHALL be symmetric:
the smoothing applied after the coarse correction SHALL be the transpose of the
smoothing applied before it, the coarsest level SHALL be solved by a symmetric
process, and the transfer operators SHALL be a transpose pair.

A cycle used directly as a solver SHALL NOT be required to be symmetric. Nothing
about a stationary iteration needs it, and the constraints that buy it — a
restriction that is the transpose of prolongation rather than the cheaper
average, equal smoothing counts, and a reversed smoother that is unstable at
relaxation factors the forward one tolerates — cost that solver convergence and
work per cycle. Where the project ships both, each SHALL be identified by which
of the two it is, and the symmetric one SHALL remain available unchanged to the
preconditioner.

#### Scenario: The cycle is symmetric

- **WHEN** one cycle of the preconditioning shape is applied from a zero initial guess to arbitrary right-hand sides `x` and `y`, over the fluid unknowns
- **THEN** `x` applied to the cycle of `y` equals `y` applied to the cycle of `x`, to machine precision

#### Scenario: The cycle is symmetric with a body present

- **WHEN** the same comparison is made on a domain containing an obstacle
- **THEN** the two inner products agree to machine precision, so the embedded boundary does not break symmetry

#### Scenario: Smoothing is reversed around the coarse correction

- **WHEN** the preconditioning cycle smooths before and after the coarse correction
- **THEN** the second pass visits the unknowns in the reverse of the order the first pass used

#### Scenario: The solver's cycle converges whether or not it is symmetric

- **WHEN** a multilevel solver runs a cycle shape that is not constrained to be symmetric
- **THEN** it converges to the same field the symmetric shape converges to, to the solve tolerance, and reports the cycles it took
