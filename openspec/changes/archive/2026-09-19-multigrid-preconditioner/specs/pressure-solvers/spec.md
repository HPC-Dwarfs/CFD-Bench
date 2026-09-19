# Spec Delta

## ADDED Requirements

### Requirement: The multilevel cycle is a symmetric operator

The multilevel cycle SHALL be a symmetric linear operator on the fluid
unknowns, so that a Krylov method preconditioned by it remains the method it
claims to be. Every part of the cycle SHALL be symmetric: the smoothing applied
after the coarse correction SHALL be the transpose of the smoothing applied
before it, the coarsest level SHALL be solved by a symmetric process, and the
transfer operators SHALL be a transpose pair.

#### Scenario: The cycle is symmetric

- **WHEN** one cycle is applied from a zero initial guess to arbitrary right-hand sides `x` and `y`, over the fluid unknowns
- **THEN** `x` applied to the cycle of `y` equals `y` applied to the cycle of `x`, to machine precision

#### Scenario: The cycle is symmetric with a body present

- **WHEN** the same comparison is made on a domain containing an obstacle
- **THEN** the two inner products agree to machine precision, so the embedded boundary does not break symmetry

#### Scenario: Smoothing is reversed around the coarse correction

- **WHEN** the cycle smooths before and after the coarse correction
- **THEN** the second pass visits the unknowns in the reverse of the order the first pass used

### Requirement: A multigrid preconditioner is shipped

The project SHALL offer the multilevel cycle as a preconditioner for its Krylov
solver, selectable in the same way as the other preconditioners. It SHALL
perform one cycle from a zero initial guess with a fixed amount of work, so that
it satisfies the requirements every preconditioner is held to.

#### Scenario: Multigrid preconditioning is selectable

- **WHEN** a setup selects the multigrid preconditioner
- **THEN** the Krylov solver uses it, reports that it is doing so, and converges to the configured tolerance

#### Scenario: The multigrid preconditioner is a fixed linear operator

- **WHEN** the multigrid preconditioner is applied to residuals of widely differing magnitude
- **THEN** each application performs the same number of cycles and the same amount of work, with no inner convergence test

### Requirement: Multigrid preconditioning is insensitive to grid refinement

A Krylov solve preconditioned by the multilevel cycle SHALL NOT lose
effectiveness as the grid is refined the way a diagonal preconditioner does.
Refining the grid SHALL increase its iteration count by substantially less than
it increases the diagonally preconditioned count for the same problem and
tolerance.

#### Scenario: Iteration count grows slowly under refinement

- **WHEN** the same problem is solved to the same tolerance on a grid and on a grid refined in every direction, with multigrid preconditioning and with diagonal preconditioning
- **THEN** the multigrid-preconditioned iteration count grows by a substantially smaller factor than the diagonally preconditioned one

#### Scenario: Multigrid preconditioning costs fewer iterations

- **WHEN** a setup containing an obstacle is solved to the same tolerance with multigrid and with diagonal preconditioning
- **THEN** the multigrid-preconditioned solve takes fewer iterations

## MODIFIED Requirements

### Requirement: Multigrid transfer operators match the cell-centred layout

Restriction and prolongation SHALL be consistent with the cell-centred placement of pressure in three dimensions, and prolongation SHALL define a value at every fine cell it is responsible for.

Restriction SHALL be a scalar multiple of the transpose of prolongation. This is what makes the cycle symmetric, and it is a property of the pair rather than of either operator alone, so neither may be changed without the other.

Both operators SHALL read only face neighbours at a subdomain boundary, since that is what the halo exchange provides; an operator that reads a diagonal neighbour would make the result depend on the decomposition.

#### Scenario: Prolongation leaves no stale values

- **WHEN** a correction is prolongated from a coarse level to a fine level
- **THEN** every fine cell in the prolongated region receives a value derived from the coarse correction, with none retaining a value from a previous cycle

#### Scenario: Transfer operators are consistent

- **WHEN** a constant field is restricted and then prolongated
- **THEN** the result is the same constant to machine precision

#### Scenario: Restriction is the transpose of prolongation

- **WHEN** an arbitrary fine-level field is restricted, and an arbitrary coarse-level field is prolongated
- **THEN** the coarse field applied to the restricted fine field equals, up to the fixed scale factor, the fine field applied to the prolongated coarse field, to machine precision

#### Scenario: Transfers agree across decompositions

- **WHEN** the same field is restricted and prolongated on one rank and on several
- **THEN** the results agree to machine precision

### Requirement: The coarse correction reaches the solution

Multigrid SHALL apply the correction computed on each coarse level to the solution on the level below it, starting each coarse solve from a zero error guess.

The correction SHALL be confined to cells with nonzero volume fraction. A coarse level represents the body only approximately, so a correction prolongated from it is defined over cells the body occupies; applying it there would move a solid cell off the zero its identity row requires, and nothing later in the cycle would take it back off.

#### Scenario: A cycle reduces the residual

- **WHEN** a multigrid cycle is run on a Poisson problem with a known solution
- **THEN** the residual after the cycle is smaller than before it, by more than pre- and post-smoothing alone would achieve

#### Scenario: Each level uses its own mesh

- **WHEN** the coarse-level operator is applied to a field whose exact Laplacian is known
- **THEN** the result matches that Laplacian on that level's mesh, not on the finest mesh

#### Scenario: The correction does not disturb the body

- **WHEN** a cycle is run on a domain containing an obstacle whose coarse representation differs from its fine one
- **THEN** every fully solid cell holds exactly zero after the correction, including cells in the body's interior that no surface correction pass visits
