#!/bin/sh
#
# Every shipped setup starts, reads the geometry it names, and says nothing
# alarming while doing it.
#
# Usage: tests/check-setups.sh
#
# Each setup is run for a couple of time steps with its final time overridden,
# because the point is that it starts correctly rather than that it reaches a
# steady state. What is checked:
#
#   - the run reaches a time step, so nothing was refused at initialization;
#   - the geometry header names what the setup asked for;
#   - a setup with a voxel volume reports at least the required voxels per cell;
#   - no aspect-ratio warning, since the shipped volumes are generated to match
#     their domains exactly;
#   - the obstacle survives to the coarsest multigrid level the setup asks for.
#
# Environment:
#   TOOLCHAIN  selects the binary suffix (default CLANG)

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${TOOLCHAIN:-CLANG}
BIN="$ROOT/CFD-Solver-$TOOLCHAIN"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

status=0

if ! make -C "$ROOT" ${SOLVER:+SOLVER="$SOLVER"} >/dev/null 2>&1; then
    echo "check-setups.sh: build failed" >&2
    exit 2
fi

# The obstacle setups need their volumes; they are generated, not committed.
if [ ! -f "$ROOT/geometry/karman.vox" ] || [ ! -f "$ROOT/geometry/backstep.vox" ]; then
    "$ROOT/tools/make-geometry.sh" >/dev/null || exit 2
fi

# $1 setup, $2 substring the geometry header has to contain
check_setup() {
    setup=$1
    want=$2

    par="$WORK/$setup.par"
    # Two steps is enough to show it runs; the shipped final times are minutes.
    sed 's/^te  *[0-9.]*/te       0.0001/' "$ROOT/testcases/flow/$setup.par" > "$par"

    printf -- '-- %s\n' "$setup"

    out=$( cd "$ROOT" && "$BIN" "$par" 2>&1 )
    rc=$?

    if [ $rc -ne 0 ]; then
        printf '%s: FAILED, exited %d\n' "$setup" "$rc"
        printf '%s\n' "$out" | tail -3
        status=1
        return
    fi

    if ! printf '%s' "$out" | grep -q "^TIME "; then
        printf '%s: FAILED, no time step was executed\n' "$setup"
        status=1
        return
    fi

    if ! printf '%s' "$out" | grep -q "$want"; then
        printf '%s: FAILED, geometry header does not mention "%s"\n' "$setup" "$want"
        printf '%s\n' "$out" | grep -A5 "Obstacle geometry"
        status=1
        return
    fi

    if printf '%s' "$out" | grep -q "aspect ratio mismatch"; then
        printf '%s: FAILED, the shipped volume does not match its domain\n' "$setup"
        status=1
        return
    fi

    if printf '%s' "$out" | grep -q "unresolved at level"; then
        printf '%s: FAILED, the obstacle does not survive the multigrid hierarchy\n' \
            "$setup"
        printf '%s\n' "$out" | grep "unresolved at level"
        status=1
        return
    fi

    perCell=$(printf '%s' "$out" | sed -n 's/.*voxels per cell: \([0-9.]*\).*/\1/p')
    if [ -n "$perCell" ]; then
        printf '%s: OK (%s voxels per cell)\n' "$setup" "$perCell"
    else
        printf '%s: OK (analytic geometry)\n' "$setup"
    fi
}

# A body has to do something to the flow. These are cheap, bounded runs: the
# point is that the geometry reaches the physics, not that the flow is
# converged.
check_physics() {
    if ! make -C "$ROOT" ${SOLVER:+SOLVER="$SOLVER"} tests >/dev/null 2>&1; then
        echo "check-setups.sh: test build failed" >&2
        status=1
        return
    fi

    printf -- '-- backstep blocks the flow through the step\n'

    par="$WORK/backstep-flow.par"
    sed 's/^te  *[0-9.]*/te       0.2/' "$ROOT/testcases/flow/backstep.par" > "$par"

    ( cd "$ROOT" && NUSIF_FIELD_DUMP="$WORK/backstep.dump" \
        "$BIN-test" "$par" >/dev/null 2>&1 )

    # The step fills x in [0,1] of 7 and y in [0,0.5] of 1.5, on a 140 x 30 x 30
    # grid: cells 0..19 in x and 0..9 in y. Probe well inside it.
    probe=$("$ROOT/tools/fieldprobe.py" "$WORK/backstep.dump" --field u \
        --box 2 17 2 7 2 27)
    absmax=$(printf '%s' "$probe" | awk '{print $3}')

    if awk -v v="$absmax" 'BEGIN { exit !(v < 1e-30) }'; then
        printf 'backstep: OK, the flow inside the step is exactly zero\n'
    else
        printf 'backstep: FAILED, |u| inside the step reaches %s\n' "$absmax"
        status=1
    fi

    printf -- '-- karman disturbs the flow behind the cylinder\n'

    par="$WORK/karman-flow.par"
    sed 's/^te  *[0-9.]*/te       2.0/' "$ROOT/testcases/flow/karman.par" > "$par"

    ( cd "$ROOT" && NUSIF_FIELD_DUMP="$WORK/karman.dump" \
        "$BIN-test" "$par" >/dev/null 2>&1 )

    # The cylinder sits at x = 5 of 30 on a 200-cell grid, so around cell 33.
    # Downstream of it the transverse velocity can only be nonzero because the
    # body is there: the inflow is uniform and the channel is straight.
    probe=$("$ROOT/tools/fieldprobe.py" "$WORK/karman.dump" --field v \
        --box 40 90 10 40 10 40)
    absmax=$(printf '%s' "$probe" | awk '{print $3}')

    if awk -v v="$absmax" 'BEGIN { exit !(v > 1e-3) }'; then
        printf 'karman: OK, transverse velocity behind the cylinder reaches %s\n' \
            "$absmax"
    else
        printf 'karman: FAILED, the flow behind the cylinder is undisturbed (|v| max %s)\n' \
            "$absmax"
        status=1
    fi
}

check_setup dcavity        "obstacle-free"
check_setup canal          "obstacle-free"
check_setup karman         "geometry/karman.vox"
check_setup backstep       "geometry/backstep.vox"
check_setup schaefer-turek "analytic cylinder"

check_physics

printf '\n==========================================\n'
if [ $status -eq 0 ]; then
    echo "every shipped setup starts"
else
    echo "SOME SHIPPED SETUPS DO NOT START"
fi
exit $status
