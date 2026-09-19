# obstacle-geometry Spec

## Purpose

Defines how obstacle geometry is supplied to the 3D solver, how it is converted into the face apertures and cell volume fractions the discretization consumes, how each rank obtains its own part of it, and which geometry inputs are accepted or rejected.

## Requirements

### Requirement: Geometry is supplied as a voxel volume

The solver SHALL accept obstacle geometry as a binary voxel volume named by a `geometryFile` entry in the setup parameter file. The format SHALL consist of an ASCII header giving the magic identifier, the three voxel counts and the maximum value, followed by one byte per voxel with x varying fastest, then y, then z. A voxel value below 128 SHALL denote solid; a value of 128 or above SHALL denote fluid. The solver SHALL NOT require any third-party image or volume library.

#### Scenario: Setup names a geometry volume

- **WHEN** a setup parameter file sets `geometryFile` to a readable voxel volume
- **THEN** the solver uses that volume as the obstacle geometry for the run

#### Scenario: Setup omits the geometry file

- **WHEN** a setup parameter file has no `geometryFile` entry
- **THEN** the solver runs with an obstacle-free domain, producing the same result as before obstacles existed

#### Scenario: Geometry file is missing or malformed

- **WHEN** `geometryFile` names a file that cannot be opened, or whose header is not valid, or whose payload is shorter than its header declares
- **THEN** the solver aborts at initialization with a message naming the file and the reason
- **AND** no time steps are executed

### Requirement: Volume maps to the physical domain with fixed orientation

The volume SHALL cover the full physical domain `[0, xlength] x [0, ylength] x [0, zlength]`. The mapping from voxel index to physical position SHALL be fixed and documented in each axis, so that the same file produces the same body on every grid and every rank count.

#### Scenario: Obstacle position matches the volume

- **WHEN** a volume places a solid region in one identified octant of its index space
- **THEN** the resulting obstacle occupies the corresponding octant of the physical domain in the computed fields and in the written output

#### Scenario: Volume aspect ratio disagrees with the domain

- **WHEN** any of the volume's three side ratios differs from the corresponding domain side ratio by more than a small tolerance
- **THEN** the solver reports a warning naming both ratios at initialization
- **AND** continues, sampling the volume anisotropically

### Requirement: Geometry is resolved into apertures and volume fractions

For every cell the solver SHALL derive a volume fraction, and for every x-face, y-face and z-face an aperture, by sampling the volume at sub-cell resolution. A value of 0 SHALL mean fully solid and 1 fully fluid. These four quantities SHALL be the sole interface between geometry and the rest of the solver; no other part of the solver SHALL consult the volume file.

#### Scenario: Fully fluid cell

- **WHEN** every voxel covering a cell and its faces is fluid
- **THEN** that cell's volume fraction and all six of its face apertures are 1

#### Scenario: Fully solid cell

- **WHEN** every voxel covering a cell and its faces is solid
- **THEN** that cell's volume fraction and all six of its face apertures are 0

#### Scenario: A solid cell never has an open face

- **WHEN** independent sampling of a cell and of one of its faces would round the cell to solid and the face to open
- **THEN** the face is closed, so that no pressure inside a body can reach the fluid

#### Scenario: Zero-thickness plate

- **WHEN** the volume contains a solid plane one voxel thick with fluid on both sides, aligned with a grid face
- **THEN** the apertures of the faces on that plane are 0 and the volume fractions of the cells on both sides remain 1
- **AND** the run proceeds without error

### Requirement: Volume resolution must exceed the grid resolution

The solver SHALL require at least a configured minimum number of voxels per grid cell in each direction, so that the represented geometry does not change as the grid is refined across a benchmark resolution sweep. The default minimum SHALL be 4 voxels per cell per direction.

#### Scenario: Volume is too coarse for the grid

- **WHEN** the volume provides fewer than the minimum voxels per cell in any direction
- **THEN** the solver aborts at initialization, reporting the volume resolution, the grid resolution and the required minimum

#### Scenario: Geometry is grid-independent across a refinement sweep

- **WHEN** the same volume is used for two grids whose resolutions differ by a factor of two, both satisfying the minimum
- **THEN** the total solid volume represented on the two grids agrees to within the sampling tolerance of the coarser grid

### Requirement: Apertures are independent of MPI decomposition

Derived apertures and volume fractions SHALL be bit-identical for a given volume and grid regardless of the number of MPI ranks or the shape of the domain decomposition. Deriving them SHALL require no communication between ranks.

#### Scenario: Same geometry under different rank counts

- **WHEN** the same setup is run on 1 rank and on 8 ranks
- **THEN** the gathered aperture and volume-fraction fields are bit-identical

#### Scenario: Obstacle straddles a subdomain boundary

- **WHEN** an obstacle spans the boundary between two subdomains
- **THEN** its apertures are identical to those produced when the same obstacle lies in the interior of a single subdomain

### Requirement: Per-rank geometry reading is bounded by the subdomain

The memory a rank uses to read geometry SHALL scale with the voxels covering that rank's own subdomain and its sampling margin, not with the size of the whole volume. A rank SHALL NOT be required to hold the complete volume in memory.

#### Scenario: Large volume on many ranks

- **WHEN** a volume far larger than one rank's share is used for a run on many ranks
- **THEN** each rank's geometry reading allocates only its own slab, and the run completes

#### Scenario: Slab reading does not change the result

- **WHEN** the same geometry is produced for a decomposition that splits the volume differently
- **THEN** the apertures and volume fractions are unchanged, because each cell is sampled from exactly the same voxels either way

### Requirement: Fluid region must be simply connected

The solver SHALL verify at initialization that the fluid part of the domain forms exactly one connected region under face connectivity — six neighbours per cell — and SHALL abort when it does not. A sealed fluid pocket adds an independent constant to the pressure null space, which the pressure solvers do not support.

#### Scenario: Volume encloses a fluid pocket

- **WHEN** the volume contains a fluid region fully surrounded by solid
- **THEN** the solver aborts at initialization, reporting the number of fluid regions found and the location of one cell in each

#### Scenario: Single fluid region

- **WHEN** all fluid cells are connected through open faces
- **THEN** initialization completes and the run proceeds

#### Scenario: Region count does not depend on the decomposition

- **WHEN** the same geometry is checked on 1 rank and on 8 ranks
- **THEN** the reported number of fluid regions is the same

### Requirement: Analytic geometry producer for verification

The solver SHALL additionally offer analytically defined obstacles — at minimum a sphere and an axis-aligned plate — producing apertures and volume fractions through the same interface as the volume-file path. This exists so that order-of-accuracy and reference benchmark studies are not limited by the voxel resolution of a rasterized body.

#### Scenario: Analytic sphere selected

- **WHEN** a setup selects the analytic sphere producer with a centre and radius
- **THEN** the solver produces apertures and volume fractions for that sphere, and the rest of the solver behaves identically to the volume-file path

#### Scenario: Analytic and rasterized geometry agree

- **WHEN** a volume is generated by rasterizing the same body at high resolution
- **THEN** the solid volume represented by the file path matches the analytic path to within the sampling tolerance of the volume

### Requirement: Geometry input is recorded for reproducibility

The solver SHALL print, in the run header, the geometry source in use and — for the volume-file path — the voxel dimensions and a checksum of the file contents.

#### Scenario: Run header identifies the geometry

- **WHEN** a run using a volume file starts
- **THEN** the header contains the filename, the voxel dimensions and a checksum sufficient to distinguish it from a modified file

#### Scenario: A modified volume is distinguishable

- **WHEN** two volumes differing in a single voxel are each used for a run
- **THEN** the checksums printed in the two run headers differ
