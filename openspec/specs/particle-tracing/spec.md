# particle-tracing Spec

## Purpose

Defines how massless tracer particles are seeded, advected, exchanged between ranks and written out, so that flow around an embedded body can be visualized as streaklines, and so that a particle never passes through a body the pressure solve treats as solid.

## Requirements

### Requirement: Particle tracing is configured per setup and is off by default

Particle tracing SHALL be enabled by the presence of its parameters in the setup file and SHALL be absent from a run that does not configure it. A setup that does not configure tracing SHALL produce the same fields, the same output files and the same run time it would produce with no tracing code present.

#### Scenario: Setup does not configure tracing

- **WHEN** a setup parameter file configures no particles
- **THEN** no particle is created, no particle output file is written, and the computed fields are bit-identical to a run of the same setup before tracing existed

#### Scenario: Setup configures tracing

- **WHEN** a setup parameter file gives a particle count, a start time, an injection period, a write period and a seed region
- **THEN** the solver traces particles from the start time onwards and reports the configuration in the run header

### Requirement: Particles are injected from a seed region into fluid only

Particles SHALL be injected within a configured axis-aligned seed region of the domain, repeatedly at a configured period, beginning at a configured start time. A particle SHALL NOT be created at a position inside a solid region.

#### Scenario: Injection begins at the start time

- **WHEN** the simulated time reaches the configured start time
- **THEN** the first batch of particles is injected, and none existed before that time

#### Scenario: Seed positions lie in the seed region

- **WHEN** a batch of particles is injected
- **THEN** every particle's position lies within the configured seed region

#### Scenario: No particle is seeded inside a body

- **WHEN** the seed region overlaps a solid region of the geometry
- **THEN** the particles that would fall inside the body are not created, and the remaining particles are injected normally

### Requirement: Injection is deterministic and independent of the decomposition

The sequence of seed positions SHALL depend only on the setup, not on the number of MPI ranks, the shape of the decomposition, or the order in which ranks execute. Repeating a run SHALL reproduce the same seed positions.

#### Scenario: Same particles on different rank counts

- **WHEN** the same setup is run on 1 rank and on 8 ranks
- **THEN** the set of particle positions injected at a given time is the same, regardless of which rank owns each one

#### Scenario: Repeated runs agree

- **WHEN** the same setup is run twice
- **THEN** the two runs inject identical seed positions

### Requirement: Particles are advected by the interpolated velocity field

Each particle SHALL be advanced using the velocity interpolated to its position from the staggered velocity components, sampling each component at its own face location. The advance SHALL use the time step the flow solver is currently taking, including when that time step is chosen adaptively.

#### Scenario: A particle follows a uniform flow

- **WHEN** particles are traced in a uniform flow of known velocity over a known interval
- **THEN** each particle's displacement matches velocity times elapsed time to the accuracy of the integration scheme

#### Scenario: Advection follows the adaptive time step

- **WHEN** a setup with adaptive time stepping changes its time step during a run
- **THEN** particle displacement over that step reflects the time step the flow solver actually took

### Requirement: Particles do not cross closed faces

A particle SHALL NOT be advected across a face whose aperture is zero, nor come to rest inside a cell whose volume fraction is zero. Whether a particle is obstructed SHALL be determined from the faces its path crosses, not from the cell it lands in, so that a body one cell face thick obstructs it.

#### Scenario: A particle meeting a body is removed

- **WHEN** a particle's path over one step would cross a face of zero aperture
- **THEN** the particle is removed from the simulation and counted as removed

#### Scenario: A zero-thickness plate obstructs particles

- **WHEN** particles are traced towards a plane of closed faces with fluid cells on both sides
- **THEN** no particle appears on the far side of the plane

#### Scenario: No particle ends inside a body

- **WHEN** a long run around a closed body completes
- **THEN** no written particle position lies inside a cell of zero volume fraction

### Requirement: Particles are exchanged between ranks and conserved

A particle whose new position lies in another rank's subdomain SHALL be transferred to that rank and continue to be traced there. A particle SHALL NOT be duplicated, and SHALL NOT be lost except by leaving the domain or being obstructed.

#### Scenario: A particle crosses a subdomain boundary

- **WHEN** a particle's new position lies outside the owning rank's subdomain but inside the domain
- **THEN** exactly one rank holds it afterwards, and it is the rank whose subdomain contains the new position

#### Scenario: Particle count is conserved

- **WHEN** a run completes
- **THEN** the total particles written, plus those removed at the domain boundary, plus those removed at a body, equals the total injected

#### Scenario: Tracing agrees across rank counts

- **WHEN** the same setup is traced on 1 rank and on 8 ranks
- **THEN** the set of particle positions written at a given time agrees to within the reproducibility limits of the flow field itself

### Requirement: Particles leaving the domain are removed

A particle advected past a domain boundary SHALL be removed and counted, and SHALL NOT be transferred to any rank.

#### Scenario: A particle leaves through an outflow

- **WHEN** a particle is advected past an outflow boundary
- **THEN** it is removed, the removal is counted, and no rank continues to trace it

### Requirement: Particle storage is bounded and does not depend on grid size

The memory a rank uses for particles SHALL be governed by the number of particles it holds, not by the number of cells in its subdomain, and SHALL grow as needed rather than being fixed at a worst-case size. Storage for removed particles SHALL be reclaimed.

#### Scenario: Storage tracks particle count, not grid size

- **WHEN** the same tracing configuration is run on two grids differing by a factor of two in each direction
- **THEN** the reported particle storage per rank is substantially the same

#### Scenario: Removed particles are reclaimed

- **WHEN** a run removes many particles over time
- **THEN** particle storage does not grow without bound, and the traced particle count matches the number of live particles

### Requirement: Particles are written in a format the project's visualization reads

Particle positions SHALL be written at a configured period, gathered across ranks, into a file per output time that the same visualization tool as the field output can open.

#### Scenario: Output is written at the configured period

- **WHEN** the configured write period elapses
- **THEN** a particle file for that time is written containing every live particle in the domain exactly once

#### Scenario: Output opens alongside the field output

- **WHEN** a particle file and the corresponding field output are loaded together
- **THEN** both open in the project's visualization tool and the particles lie within the domain the fields cover
