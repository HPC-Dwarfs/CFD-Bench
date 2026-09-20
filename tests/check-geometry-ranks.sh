#!/bin/sh
#
# Geometry must not depend on how the domain was divided.
#
# Usage: tests/check-geometry-ranks.sh [ranks ...]     (default 1 2 4 8)
#
# Runs the geometry driver at each rank count and requires that the aperture
# checksums and the number of connected fluid regions are identical, and that a
# rank's voxel slab shrinks as the domain is divided.
#
# Environment:
#   TOOLCHAIN  selects the build directory (default CLANG)
#   MPIRUN     mpirun command (default mpirun)

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${TOOLCHAIN:-CLANG}
MPIRUN=${MPIRUN:-mpirun}
BIN="$ROOT/build/$TOOLCHAIN/check-geometry-mpi"

if [ $# -gt 0 ]; then
    ranks="$*"
else
    ranks="1 2 4 8"
fi

if ! make -C "$ROOT" ${SOLVER:+SOLVER="$SOLVER"} tests >/dev/null; then
    echo "check-geometry-ranks.sh: build failed" >&2
    exit 2
fi

if [ ! -f "$ROOT/tests/geom/sphere.vox" ]; then
    "$ROOT/tests/make-geom.sh" >/dev/null || exit 2
fi

status=0
reference=""
refSlab=""

for n in $ranks; do
    out=$( cd "$ROOT" && "$MPIRUN" -n "$n" "$BIN" 2>/dev/null )

    if [ $? -ne 0 ]; then
        printf '%s ranks: FAILED, the driver exited nonzero\n' "$n"
        status=1
        continue
    fi

    line=$(printf '%s' "$out" | grep -E '^(APERTURE-CHECKSUM|REGIONS)' | tr '\n' ' ')
    slab=$(printf '%s' "$out" | grep '^SLAB rank 0 ' | awk '{print $6}')

    if [ -z "$reference" ]; then
        reference=$line
        refSlab=$slab
        printf '%s ranks: %s(reference, slab %s bytes)\n' "$n" "$line" "$slab"
        continue
    fi

    if [ "$line" = "$reference" ]; then
        printf '%s ranks: identical (slab %s bytes)\n' "$n" "$slab"
    else
        printf '%s ranks: DIFFERS\n  got %s\n  want %s\n' "$n" "$line" "$reference"
        status=1
    fi

    if [ "$slab" -ge "$refSlab" ]; then
        printf '%s ranks: FAILED, slab did not shrink (%s bytes against %s on one rank)\n' \
            "$n" "$slab" "$refSlab"
        status=1
    fi
done

printf '\n==========================================\n'
if [ $status -eq 0 ]; then
    echo "geometry is independent of the decomposition"
else
    echo "GEOMETRY DEPENDS ON THE DECOMPOSITION"
fi
exit $status
