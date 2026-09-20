#!/bin/sh
#
# Regenerate the voxel volumes the shipped setups reference.
#
# Usage: tools/make-geometry.sh
#
# The volumes are generated rather than committed. In two dimensions the
# equivalent image is a few megabytes and lives in the repository; a volume of
# the same relative resolution is 64 bytes per grid cell, so the karman volume
# below is 32 MB and a refinement of it would be a quarter of a gigabyte. That
# does not belong in version control, and it costs seconds to rebuild.
#
# Each resolution is chosen so that the volume clears the solver's minimum of 4
# voxels per cell per direction on the setup's own grid, and so that the volume
# aspect ratio matches the domain exactly and no anisotropic-sampling warning is
# produced. A refinement study needs a volume regenerated at a matching
# resolution; the minimum is enforced, so an under-resolved one is refused
# rather than silently sampled.
#
# schaefer-turek.par is not listed here: it uses the analytic cylinder producer,
# because a rasterized body is only as accurate as its voxel size and that would
# cap the very quantity the benchmark measures.

set -eu

ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT="$ROOT/geometry"
GEN="$ROOT/tools/genvox.py"

mkdir -p "$OUT"

# karman: a 30 x 8 x 8 channel with a cylinder of radius 1 along z at (5, 4),
#   benchmark grid 200 x 50 x 50.
#   The domain is 15:4:4, so the volume has to be a multiple of that to avoid an
#   anisotropic-sampling warning. 54 x (15, 4, 4) is the smallest multiple that
#   also clears 4 voxels per cell on this grid, giving 4.05 x 4.32 x 4.32.
"$GEN" cylinder --out "$OUT/karman.vox" --size 810 216 216 \
    --domain 30 8 8 --axis z --center 5 4 --radius 1

# backstep: a 7 x 1.5 x 1.5 channel with a step filling x in [0,1], y in [0,0.5]
#   and the full z extent, benchmark grid 140 x 30 x 30.
#   560 x 120 x 120 is 14:3:3 like the domain, and gives 4 voxels per cell.
"$GEN" box --out "$OUT/backstep.vox" --size 560 120 120 \
    --domain 7 1.5 1.5 --corner 0 0 0 --extent 1 0.5 1.5

ls -l "$OUT"
