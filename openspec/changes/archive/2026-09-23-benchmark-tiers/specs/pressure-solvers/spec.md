# Spec Delta

## ADDED Requirements

### Requirement: A setup intended for scaling measurements keeps its hierarchy

A setup shipped for measuring performance at scale SHALL be sized so that the
multilevel hierarchy retains a useful depth at every rank count it is intended
for. Coarsening is limited by the local extents, so depth is a property of the
grid and the decomposition together, and a grid chosen without regard to the
rank counts it will meet can lose its hierarchy entirely.

This SHALL be checked rather than assumed. A hierarchy that collapses does not
fail: it converges more slowly and reports nothing unusual, so a benchmark run
on a collapsed hierarchy measures a smoother and presents it as multigrid.

#### Scenario: Depth survives the decomposition

- **WHEN** a setup intended for scaling is divided across any rank count it targets
- **THEN** the hierarchy built on each rank still has at least the depth the setup claims, and the solver reports the depth it built

#### Scenario: A grid that cannot be coarsened is rejected as a benchmark

- **WHEN** a setup intended for scaling has an extent that stops being divisible by two before the depth it targets
- **THEN** that is reported as a defect of the setup, rather than being discovered later as poor convergence

#### Scenario: The claimed depth is verified against the solver

- **WHEN** the depth a setup is expected to reach is computed from its grid and rank count
- **THEN** it agrees with the depth the solver reports at the rank counts where that can be run directly, so the expectation is checked against the implementation rather than restating it

### Requirement: A scaling measurement performs a fixed amount of work

A setup shipped for measuring performance SHALL perform an amount of work that
depends on its grid and on nothing else. The number of time steps SHALL be fixed
by configuration rather than chosen by an adaptive controller, so that two sizes
of the same case differ only in the cells they cover.

#### Scenario: Two sizes of the same case are comparable

- **WHEN** the same case is run at two grid sizes intended for comparison
- **THEN** both execute the same number of time steps, so the difference in run time is attributable to the problem size

#### Scenario: The step count does not depend on the flow

- **WHEN** a setup intended for scaling is run
- **THEN** the time step is the one its configuration states, not one derived from the velocity field, and the run executes the number of steps that configuration implies
