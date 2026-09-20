# Spec Delta

## Purpose

Defines what the pressure solvers must deliver when the domain contains embedded obstacle boundaries, so that red-black SOR, compressed red-black SOR and multigrid remain interchangeable first-class solvers and a preconditioned CG solver can be added without changing the operator.

## ADDED Requirements

### Requirement: All pressure solvers solve the same system

Every pressure solver the project ships SHALL solve the same embedded-boundary pressure system, with the same treatment of closed faces and solid cells. Given the same setup and a sufficiently tight tolerance, the solvers SHALL produce the same pressure field.

#### Scenario: Solvers agree

- **WHEN** the same setup with an obstacle is solved by each shipped solver, each converged to a tight tolerance
- **THEN** the resulting pressure fields agree with one another to within that tolerance

#### Scenario: Obstacle is visible to every solver

- **WHEN** an obstacle is introduced into a setup
- **THEN** every solver produces a pressure field whose gradient across zero-aperture faces is zero to the solve tolerance

#### Scenario: A solver using a different internal data layout gives the same answer

- **WHEN** a solver stores the unknowns in a layout other than the domain's natural one in order to vectorize its sweep
- **THEN** its converged field with an obstacle present matches the other solvers' to the solve tolerance, so the layout is not observable in the result

### Requirement: Relaxation uses the correct per-cell diagonal

Relaxation-based solvers SHALL use the diagonal implied by each cell's open faces rather than a single domain-wide constant, so that cells adjacent to an obstacle are relaxed correctly.

#### Scenario: Relaxation converges next to an obstacle

- **WHEN** a relaxation solver is run on a setup with an obstacle
- **THEN** the residual decreases monotonically in the asymptotic regime, with no divergence or stagnation localized at the obstacle surface

### Requirement: Cut cells are updated consistently by relaxation

Cells adjacent to the obstacle SHALL be relaxed in a dedicated pass, ordered consistently with the colour sweeps, so that the overall iteration remains a convergent splitting of the operator and does not depend on the traversal order within a pass.

#### Scenario: Cut-cell pass preserves convergence

- **WHEN** a relaxation solver with the cut-cell pass is run to convergence on a setup with an obstacle
- **THEN** it converges, and its iteration count is within a small factor of the same setup without an obstacle

#### Scenario: Neighbouring cut cells do not see each other within a pass

- **WHEN** two cut cells that are face neighbours are relaxed
- **THEN** they are relaxed in different passes, so neither result depends on which was visited first

### Requirement: Colouring is derived from global position

The red-black colouring of every sweep, and the colour segments of the cut-cell pass, SHALL be derived from each cell's position in the global domain rather than from its index within a rank's subdomain.

#### Scenario: Colouring does not shift with the decomposition

- **WHEN** the same cell is relaxed under two different domain decompositions
- **THEN** it is assigned the same colour in both, and is therefore relaxed in the same pass

### Requirement: Multigrid coarsens the geometry consistently

Multigrid SHALL construct a geometry hierarchy by coarsening apertures and volume fractions, so that every level solves a consistent embedded-boundary problem.

#### Scenario: Coarse levels see the obstacle

- **WHEN** a multigrid cycle runs on a setup with an obstacle resolved by several cells
- **THEN** the obstacle is represented at every level of the hierarchy

#### Scenario: Convergence rate is retained

- **WHEN** multigrid is run on a setup with a resolved obstacle
- **THEN** the asymptotic residual reduction per cycle is within a small factor of the same setup without an obstacle

#### Scenario: Obstacle too small for the coarsest level

- **WHEN** an obstacle feature becomes unresolved at a coarse level
- **THEN** the solver still converges, and reports the level at which the feature was lost

### Requirement: Multigrid transfer operators match the cell-centred layout

Restriction and prolongation SHALL be consistent with the cell-centred placement of pressure in three dimensions, and prolongation SHALL define a value at every fine cell it is responsible for.

#### Scenario: Prolongation leaves no stale values

- **WHEN** a correction is prolongated from a coarse level to a fine level
- **THEN** every fine cell in the prolongated region receives a value derived from the coarse correction, with none retaining a value from a previous cycle

#### Scenario: Transfer operators are consistent

- **WHEN** a constant field is restricted and then prolongated
- **THEN** the result is the same constant to machine precision

### Requirement: The coarse correction reaches the solution

Multigrid SHALL apply the correction computed on each coarse level to the solution on the level below it, starting each coarse solve from a zero error guess.

#### Scenario: A cycle reduces the residual

- **WHEN** a multigrid cycle is run on a Poisson problem with a known solution
- **THEN** the residual after the cycle is smaller than before it, by more than pre- and post-smoothing alone would achieve

#### Scenario: Each level uses its own mesh

- **WHEN** the coarse-level operator is applied to a field whose exact Laplacian is known
- **THEN** the result matches that Laplacian on that level's mesh, not on the finest mesh

### Requirement: Multigrid iterates to the requested tolerance

Multigrid SHALL repeat cycles until the configured tolerance is reached or the configured iteration limit is exhausted, and SHALL report the number of cycles used.

#### Scenario: Multigrid honours the tolerance

- **WHEN** a setup is solved by multigrid with a given tolerance and iteration limit
- **THEN** the solve either reaches that tolerance or stops at the limit, and the cycle count is reported

### Requirement: The pressure boundary condition follows the setup

The pressure boundary condition imposed by every solver SHALL follow the setup's own velocity boundary configuration: a wall or slip boundary imposes a zero normal pressure gradient, and an outflow boundary pins the pressure on that boundary face. A boundary type the solver does not support SHALL be rejected at initialization rather than silently treated as a wall. The condition SHALL be applied at every level of a multilevel solver.

#### Scenario: An outflow setup is non-singular

- **WHEN** a setup with an outflow boundary is solved
- **THEN** the outflow pins the pressure, the system is non-singular, and every solver converges to the stated tolerance

#### Scenario: An unsupported boundary type is rejected

- **WHEN** a setup requests a pressure boundary type the solver does not implement
- **THEN** the solver aborts at initialization naming the boundary and the type, rather than applying a different condition

### Requirement: The pressure null space is handled explicitly

Where all boundaries impose a zero normal pressure gradient, the system is singular with a one-dimensional null space. The solver SHALL detect that case from the boundary-condition configuration, project the right-hand side onto the range of the operator, and remove the constant component from the pressure field, so that the iterate does not drift.

#### Scenario: Pressure does not drift

- **WHEN** a fully enclosed setup is advanced over many time steps
- **THEN** the mean pressure over fluid cells remains bounded and does not grow with the number of steps

#### Scenario: Right-hand side is made compatible

- **WHEN** discretization produces a right-hand side whose fluid-cell sum is nonzero through rounding
- **THEN** the solver removes that component before iterating, and the solve converges to the stated tolerance

#### Scenario: Null-space handling does not disturb a non-singular setup

- **WHEN** the setup includes an outflow boundary, making the system non-singular
- **THEN** the pressure field is unchanged by the null-space handling

### Requirement: Solver reporting reflects the fluid domain

Residuals, iteration counts and convergence messages SHALL be computed over fluid unknowns only, SHALL be independent of the number of solid cells in the domain, and SHALL reflect the current iterate rather than accumulating across iterations.

#### Scenario: Solid cells do not flatter the residual

- **WHEN** the solid fraction of a domain is increased while the fluid region and its solution are unchanged
- **THEN** the reported residual norm is unchanged

#### Scenario: The reported residual is the residual

- **WHEN** a solve converges and reports a residual norm
- **THEN** that norm matches an independently computed norm of the residual of the returned field

### Requirement: Solver results are independent of the MPI decomposition

For a given setup and tolerance, each solver SHALL produce the same converged solution and the same iteration count regardless of the number of MPI ranks, up to the reproducibility limits of floating-point reduction order.

#### Scenario: Same answer on different rank counts

- **WHEN** the same obstacle setup is solved on 1 rank and on 8 ranks
- **THEN** the converged pressure fields agree to within the solve tolerance, and the iteration counts are equal
