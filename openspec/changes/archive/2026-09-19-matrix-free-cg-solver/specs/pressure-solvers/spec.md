# Spec Delta

## ADDED Requirements

### Requirement: A Krylov pressure solver is shipped

The project SHALL ship a conjugate gradient pressure solver as a first-class
solver variant, selected the same way as the relaxation solvers, and subject to
the same requirement that every shipped solver solves the same embedded-boundary
system and agrees with the others to the solve tolerance.

#### Scenario: Krylov solver agrees with the relaxation solvers

- **WHEN** a setup containing an obstacle is solved by the Krylov solver and by each relaxation solver, each converged to a tight tolerance
- **THEN** all of the resulting pressure fields agree with one another to within that tolerance

#### Scenario: Krylov solver honours the tolerance

- **WHEN** a setup is solved by the Krylov solver with a given tolerance and iteration limit
- **THEN** the solve either reaches that tolerance or stops at the limit, and the iteration count is reported in the same form as the other solvers report theirs

#### Scenario: The reported residual is comparable across solvers

- **WHEN** the Krylov solver reports a residual norm for a converged solve
- **THEN** that norm is the same quantity the other solvers report — the mean square residual over the fluid unknowns — so the same tolerance means the same thing for every solver

### Requirement: Operator application obeys the geometry-free bulk split

Applying the pressure operator to a vector SHALL use the same split the
relaxation sweeps use: a bulk pass that reads no geometry, followed by a
correction pass confined to the cells adjacent to the obstacle surface. A
solver that applies the operator once per iteration SHALL NOT read per-cell
geometry over the whole domain to do so.

#### Scenario: Operator application cost is independent of geometry

- **WHEN** the operator is applied on the same grid with no obstacle and with an obstacle occupying a small fraction of the domain
- **THEN** the measured time of the bulk pass agrees between the two runs within measurement noise

#### Scenario: Application and relaxation agree on the operator

- **WHEN** the operator is applied to a field, and the residual of that same field is computed by the independent residual routine the solvers report from
- **THEN** the two agree to machine precision, including at cells adjacent to the obstacle and at solid cells

### Requirement: Krylov inner products are taken over the fluid subspace

All inner products, norms and preconditioner applications in a Krylov solver
SHALL be restricted to cells with nonzero volume fraction. Solid cells SHALL
hold exactly zero at every iteration, not only at convergence, and SHALL
contribute nothing to any scalar the iteration computes.

#### Scenario: Solid volume does not affect the iteration

- **WHEN** the solid fraction of a domain is increased while the fluid region and its solution are unchanged
- **THEN** the Krylov solver's iteration count and its converged fluid field are unchanged

#### Scenario: Solid cells are zero throughout

- **WHEN** the Krylov iteration is stopped at any intermediate iteration
- **THEN** the pressure in every fully solid cell is exactly zero at that point

### Requirement: A preconditioner is a fixed linear operator

A preconditioner used by a Krylov solver SHALL be a fixed linear operator: it
SHALL perform the same amount of work on every application, SHALL NOT use an
inner convergence test, and SHALL start from a zero initial guess so that its
result depends linearly on the residual it is given.

#### Scenario: Preconditioner is linear

- **WHEN** the preconditioner is applied to two vectors and to a linear combination of them
- **THEN** the result for the combination equals the same combination of the individual results, to machine precision

#### Scenario: Preconditioner is symmetric

- **WHEN** the preconditioner is applied to arbitrary vectors `x` and `y` over the fluid unknowns
- **THEN** `x` applied to `M y` equals `y` applied to `M x` to machine precision

#### Scenario: Preconditioner work does not vary with the residual

- **WHEN** the same preconditioner is applied to residuals of widely differing magnitude
- **THEN** each application performs the same operations, so the Krylov recurrence sees one fixed operator

### Requirement: A preconditioner accelerates convergence

A preconditioned Krylov solve SHALL take no more iterations than the
unpreconditioned solve of the same problem to the same tolerance, so that the
preconditioner is demonstrably doing work rather than merely being applied.

#### Scenario: Preconditioning reduces the iteration count

- **WHEN** the same setup is solved to the same tolerance with and without preconditioning
- **THEN** the preconditioned solve takes no more iterations than the unpreconditioned one

## MODIFIED Requirements

### Requirement: The pressure null space is handled explicitly

Where all boundaries impose a zero normal pressure gradient, the system is singular with a one-dimensional null space. The solver SHALL detect that case from the boundary-condition configuration, project the right-hand side onto the range of the operator, and remove the constant component from the pressure field, so that the iterate does not drift.

The null vector SHALL be the constant over the fluid cells and zero over the solid cells, which is what the operator's null space actually contains. Every projection SHALL therefore average over fluid cells only and SHALL leave solid cells at zero. A solver that regenerates the constant during its iteration SHALL project it out at each iteration rather than only before the solve.

#### Scenario: Pressure does not drift

- **WHEN** a fully enclosed setup is advanced over many time steps
- **THEN** the mean pressure over fluid cells remains bounded and does not grow with the number of steps

#### Scenario: Right-hand side is made compatible

- **WHEN** discretization produces a right-hand side whose fluid-cell sum is nonzero through rounding
- **THEN** the solver removes that component before iterating, and the solve converges to the stated tolerance

#### Scenario: Null-space handling does not disturb a non-singular setup

- **WHEN** the setup includes an outflow boundary, making the system non-singular
- **THEN** the pressure field is unchanged by the null-space handling

#### Scenario: Projection does not disturb solid cells

- **WHEN** the null-space projection runs on a singular setup containing an obstacle
- **THEN** solid cells still hold exactly zero afterwards, and the constant removed is the mean over the fluid cells rather than over all cells

#### Scenario: A singular setup with an obstacle converges

- **WHEN** a fully enclosed setup containing an obstacle is solved by each shipped solver
- **THEN** each converges to the stated tolerance and the converged fields agree with one another to within it

### Requirement: Solver results are independent of the MPI decomposition

For a given setup and tolerance, each solver SHALL produce the same converged solution regardless of the number of MPI ranks, up to the reproducibility limits of floating-point reduction order.

A solver whose iterates do not depend on a global reduction SHALL produce the same iteration count on every rank count. A Krylov solver, whose step lengths are formed from global inner products and therefore depend on reduction order, SHALL produce iteration counts that differ by at most one across rank counts. The converged field requirement is the same for both.

#### Scenario: Same answer on different rank counts

- **WHEN** the same obstacle setup is solved on 1 rank and on 8 ranks
- **THEN** the converged pressure fields agree to within the solve tolerance

#### Scenario: Relaxation and multilevel iteration counts are equal

- **WHEN** the same obstacle setup is solved on 1 rank and on 8 ranks by a solver whose iterates do not depend on a global reduction
- **THEN** the iteration counts are equal

#### Scenario: Krylov iteration counts agree within one

- **WHEN** the same obstacle setup is solved on 1 rank and on 8 ranks by a Krylov solver
- **THEN** the iteration counts differ by at most one
