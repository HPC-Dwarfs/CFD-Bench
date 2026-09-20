# Spec Delta

## MODIFIED Requirements

### Requirement: Solid cells are carried as identity rows

Cells with zero volume fraction SHALL remain part of the unknown vector and SHALL be assigned identity rows with zero right-hand side, so that the operand arrays stay rectangular and require no index compression. Their pressure SHALL be zero and SHALL not influence any fluid cell.

This SHALL hold for every update a solver makes, not only for the relaxation sweep and not only at convergence. A solver that adds a correction computed elsewhere — on a coarser grid, or by a preconditioner — SHALL confine that correction to cells with nonzero volume fraction. The self-maintaining property of the geometry-free sweep covers a cell it relaxes from a zero neighbourhood; it does not cover a value written into the body from outside, and a cell in the body's interior is deliberately absent from the surface list that would otherwise repair it.

#### Scenario: Solid pressure stays zero

- **WHEN** the pressure solve converges on a domain containing solid cells
- **THEN** the pressure in every fully solid cell is exactly zero

#### Scenario: Deep solid interiors maintain themselves

- **WHEN** a solid cell all of whose six neighbours are also solid is relaxed by the geometry-free sweep
- **THEN** its value remains exactly zero, so it needs no correction

#### Scenario: Corrections from another grid do not reach the body

- **WHEN** a solver adds a correction computed on a coarser grid, whose representation of the body differs from the fine one
- **THEN** every fully solid cell still holds exactly zero afterwards, including cells in the body's interior

#### Scenario: Residual measures fluid only

- **WHEN** the solver reports a residual norm
- **THEN** solid cells contribute nothing to it, and the reported norm is normalized by the number of fluid cells
