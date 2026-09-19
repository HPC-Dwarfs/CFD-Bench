# Spec Delta

## Purpose

Defines how the projection method on the 3D staggered grid honours embedded obstacle boundaries, and the numerical properties the resulting pressure operator must guarantee so that every pressure solver — current and planned — can share it.

## ADDED Requirements

### Requirement: Geometry enters the discretization only through apertures

The momentum predictor, the velocity correction and the pressure operator SHALL derive all obstacle information from face apertures and cell volume fractions. The discretization SHALL NOT classify cells into boundary orientation types, and SHALL NOT impose any restriction on which combinations of open and closed faces a cell may have.

#### Scenario: Isolated solid cell

- **WHEN** the geometry produces a single solid cell surrounded by fluid on all six sides
- **THEN** the run proceeds and produces a solution, with no classification error and no special case

#### Scenario: Zero-thickness plate

- **WHEN** the geometry produces a plane of closed faces with fluid cells on both sides and no solid cell at all
- **THEN** the run proceeds and produces a solution in which no flow crosses the plane

#### Scenario: Obstacle-free domain is unaffected

- **WHEN** all apertures and volume fractions are 1
- **THEN** the computed fields are bit-identical to those of the solver before apertures were introduced

### Requirement: No-slip is enforced at closed faces

The normal velocity at a face with zero aperture SHALL be zero at every stage of the time step, including immediately after the velocity correction. Tangential velocity adjacent to a closed face SHALL satisfy the no-slip condition to the accuracy of the scheme.

#### Scenario: Normal velocity vanishes on the obstacle surface

- **WHEN** a time step completes with an obstacle present
- **THEN** the velocity component normal to every zero-aperture face is zero to machine precision

#### Scenario: No flow through a closed body

- **WHEN** a channel flow around a closed body is advanced to a steady state
- **THEN** the volume flux through the body's surface is zero to machine precision

#### Scenario: Downstream consumers see a consistent field

- **WHEN** particle tracing samples the velocity field after the velocity correction
- **THEN** it observes zero normal velocity at obstacle faces, and no particle is advected through a closed face

### Requirement: Velocity correction respects closed faces

The velocity correction SHALL NOT apply a pressure gradient across a closed face, and SHALL NOT update velocity components that lie on closed faces or inside solid regions.

#### Scenario: Pressure inside a body does not drive the flow

- **WHEN** the pressure field inside a solid region is perturbed arbitrarily before the velocity correction
- **THEN** the corrected velocity field in the fluid is unchanged

### Requirement: Fluid velocity is discretely divergence-free

After the velocity correction, the aperture-weighted divergence of the velocity field SHALL vanish in every fluid cell, to the tolerance of the pressure solve, including in cells adjacent to an obstacle.

#### Scenario: Divergence near the body

- **WHEN** a time step completes with the pressure solve converged to its tolerance
- **THEN** the maximum aperture-weighted divergence over fluid cells adjacent to the obstacle is within the same tolerance as over the rest of the domain

#### Scenario: Mass is conserved around a closed body

- **WHEN** a steady channel flow around an obstacle is advanced to a steady state
- **THEN** the net volume flux through any cross-section of the channel is constant to the tolerance of the pressure solve

### Requirement: Pressure operator is symmetric and positive definite

The pressure operator SHALL be assembled as a balance of face fluxes over the six faces of each cell, so that each face contributes an equal and opposite coefficient to the two cells it separates. The resulting operator SHALL be symmetric, and positive definite once its constant null vector is removed. This property SHALL hold for any face weights, so that later introduction of fractional apertures does not change it.

#### Scenario: Operator symmetry

- **WHEN** the pressure operator is applied to arbitrary vectors `x` and `y` over the fluid unknowns
- **THEN** `x` applied to `A y` equals `y` applied to `A x` to machine precision

#### Scenario: Symmetry with an obstacle present

- **WHEN** the same test is run on a domain containing an obstacle with closed faces
- **THEN** symmetry still holds to machine precision

### Requirement: Solid cells are carried as identity rows

Cells with zero volume fraction SHALL remain part of the unknown vector and SHALL be assigned identity rows with zero right-hand side, so that the operand arrays stay rectangular and require no index compression. Their pressure SHALL be zero and SHALL not influence any fluid cell.

#### Scenario: Solid pressure stays zero

- **WHEN** the pressure solve converges on a domain containing solid cells
- **THEN** the pressure in every fully solid cell is exactly zero

#### Scenario: Deep solid interiors maintain themselves

- **WHEN** a solid cell all of whose six neighbours are also solid is relaxed by the geometry-free sweep
- **THEN** its value remains exactly zero, so it needs no correction

#### Scenario: Residual measures fluid only

- **WHEN** the solver reports a residual norm
- **THEN** solid cells contribute nothing to it, and the reported norm is normalized by the number of fluid cells

### Requirement: Interior sweep cost is independent of geometry

The cost of the interior pressure sweep SHALL NOT depend on the presence or complexity of obstacle geometry, in any of the solvers the project ships. Obstacle handling SHALL be confined to a separate pass whose work scales with the number of cells adjacent to the obstacle surface, not with the number of cells in the domain nor with the solid volume.

#### Scenario: Benchmark comparability

- **WHEN** the same grid is run with no obstacle and with an obstacle occupying a small fraction of the domain
- **THEN** the measured time of the interior sweep agrees between the two runs within measurement noise, for each shipped solver

#### Scenario: Geometry complexity does not affect the interior sweep

- **WHEN** a simple obstacle is replaced by a geometrically complex one of comparable solid volume
- **THEN** the measured time of the interior sweep is unchanged within measurement noise

#### Scenario: Surface work scales with surface area

- **WHEN** the same body is represented at two sizes on the same grid
- **THEN** the number of cells in the correction pass scales with the body's surface area rather than with its solid volume

### Requirement: Pressure operator remains a compact seven-point stencil

On the finest grid the pressure operator SHALL couple each cell only to its six face neighbours, for any obstacle geometry.

#### Scenario: No diagonal or edge coupling

- **WHEN** the operator is applied to a unit vector at a single cell adjacent to an obstacle corner or edge
- **THEN** the result is nonzero only at that cell and its six face neighbours
