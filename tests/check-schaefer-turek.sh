#!/bin/sh
#
# Run the Schaefer-Turek 3D cylinder case and report drag, lift and Strouhal
# number against their published ranges.
#
# Usage: tests/check-schaefer-turek.sh
#
# This reports rather than asserts by default. The quantities depend on accuracy
# this phase of the work deliberately does not deliver: binary apertures give
# the cylinder a staircase surface and a first-order wall treatment, so the
# drag is expected to come out high and the wake may not shed at all. Set
# STRICT=1 to make a quantity outside its published range a failure, which is
# the measurement a later fractional-aperture change has to satisfy.
#
# The shipped setup at its own resolution run to its own final time is hours of
# work. By default this runs a shortened, coarser version, which is enough to
# show the force integration works and to give a drag figure; set FULL=1 for the
# real thing.
#
# Environment:
#   TOOLCHAIN  selects the binary suffix (default CLANG)
#   FULL       1 to run the shipped setup unmodified
#   STRICT     1 to fail when a quantity is outside its published range
#   TE         override the final time of the shortened run

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${TOOLCHAIN:-CLANG}
FULL=${FULL:-0}
STRICT=${STRICT:-0}
TE=${TE:-0.5}

BIN="$ROOT/CFD-Solver-$TOOLCHAIN"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

if ! make -C "$ROOT" ${SOLVER:+SOLVER="$SOLVER"} >/dev/null 2>&1; then
    echo "check-schaefer-turek.sh: build failed" >&2
    exit 2
fi

par="$WORK/st.par"

if [ "$FULL" = "1" ]; then
    cp "$ROOT/testcases/flow/schaefer-turek.par" "$par"
    echo "Running the shipped setup unmodified. This takes hours."
else
    # Half the resolution in every direction and a short final time. The
    # geometry, the Reynolds number and the inlet profile are untouched, so the
    # drag is still the benchmark's drag, measured on a coarser mesh.
    sed -e "s/^te  *[0-9.]*/te       $TE/" \
        -e 's/^imax  *[0-9]*/imax          128/' \
        -e 's/^jmax  *[0-9]*/jmax          24/' \
        -e 's/^kmax  *[0-9]*/kmax          24/' \
        -e 's/^levels  *[0-9]*/levels        2/' \
        "$ROOT/testcases/flow/schaefer-turek.par" > "$par"
    echo "Running a shortened, half-resolution version; set FULL=1 for the shipped setup."
fi

( cd "$WORK" && "$BIN" "$par" > "$WORK/run.log" 2>&1 )
rc=$?

if [ $rc -ne 0 ]; then
    echo "check-schaefer-turek.sh: the run failed" >&2
    tail -5 "$WORK/run.log" >&2
    exit 1
fi

if [ ! -f "$WORK/forces.dat" ]; then
    echo "check-schaefer-turek.sh: no force series was written" >&2
    exit 1
fi

printf '\n'

if [ "$STRICT" = "1" ]; then
    "$ROOT/tools/stcoeffs.py" "$WORK/forces.dat" --strict
    status=$?
else
    "$ROOT/tools/stcoeffs.py" "$WORK/forces.dat"
    status=$?
fi

printf '\n==========================================\n'
if [ $status -eq 0 ]; then
    echo "Schaefer-Turek measured"
else
    echo "SCHAEFER-TUREK OUTSIDE ITS PUBLISHED RANGES"
fi
exit $status
