# pressure-solvers Spec

## Purpose

Defines what the pressure solvers must deliver when the domain contains embedded obstacle boundaries, so that red-black SOR, compressed red-black SOR and multigrid remain interchangeable first-class solvers and a preconditioned CG solver can be added without changing the operator.

## Requirements

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

### Requirement: Multigrid iterates to the requested tolerance

Multigrid SHALL repeat cycles until the configured tolerance is reached or the configured iteration limit is exhausted, and SHALL report the number of cycles used.

#### Scenario: Multigrid honours the tolerance

- **WHEN** a setup is solved by multigrid with a given tolerance and iteration limit
- **THEN** the solve either reaches that tolerance or stops at the limit, and the cycle count is reported

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

### Requirement: Solver reporting reflects the fluid domain

Residuals, iteration counts and convergence messages SHALL be computed over fluid unknowns only, SHALL be independent of the number of solid cells in the domain, and SHALL reflect the current iterate rather than accumulating across iterations.

#### Scenario: Solid cells do not flatter the residual

- **WHEN** the solid fraction of a domain is increased while the fluid region and its solution are unchanged
- **THEN** the reported residual norm is unchanged

#### Scenario: The reported residual is the residual

- **WHEN** a solve converges and reports a residual norm
- **THEN** that norm matches an independently computed norm of the residual of the returned field

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
