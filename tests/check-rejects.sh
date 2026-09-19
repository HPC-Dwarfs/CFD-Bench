#!/bin/sh
#
# Inputs the solver has to refuse, and refuse before doing any work.
#
# Usage: tests/check-rejects.sh
#
# Each case builds a setup file, runs it, and requires a nonzero exit status, a
# message naming the problem, and no time step having been executed.
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

if [ ! -x "$BIN" ]; then
    if ! make -C "$ROOT" >/dev/null; then
        echo "check-rejects.sh: build failed" >&2
        exit 2
    fi
fi

# $1 label, $2 setup file, $3 expected substring in the message
expect_reject() {
    label=$1
    par=$2
    want=$3

    out=$( cd "$WORK" && "$BIN" "$par" 2>&1 )
    rc=$?

    printf -- '-- %s\n' "$label"

    if [ $rc -eq 0 ]; then
        printf '%s: FAILED, the run was accepted\n' "$label"
        status=1
        return
    fi

    if ! printf '%s' "$out" | grep -q "$want"; then
        printf '%s: FAILED, message did not mention "%s"\n' "$label" "$want"
        printf '   got: %s\n' "$(printf '%s' "$out" | tail -3)"
        status=1
        return
    fi

    if printf '%s' "$out" | grep -q "^TIME "; then
        printf '%s: FAILED, time steps were executed before the refusal\n' "$label"
        status=1
        return
    fi

    printf '%s: OK\n' "$label"
}

base() {
    cat <<EOF
name check-reject
bcLeft    $1
bcRight   $2
bcBottom  1
bcTop     1
bcFront   1
bcBack    1
gx 0.0
gy 0.0
gz 0.0
re 100.0
u_init 0.0
v_init 0.0
w_init 0.0
p_init 0.0
xlength 1.0
ylength 1.0
zlength 1.0
imax 8
jmax 8
kmax 8
te 0.01
dt 0.01
tau 0.5
levels 1
presmooth 2
postsmooth 2
itermax 100
eps 0.001
omg 1.7
gamma 0.9
EOF
}

base 4 1 > "$WORK/periodic.par"
expect_reject "periodic boundary is refused" "$WORK/periodic.par" "PERIODIC"

base 7 1 > "$WORK/unknown.par"
expect_reject "unknown boundary code is refused" "$WORK/unknown.par" "unknown boundary condition"

# Geometry files the solver has to refuse. Each needs its own setup because the
# geometry name is a parameter.
with_geometry() {
    base 1 1
    printf 'geometryFile %s\n' "$1"
}

"$ROOT/tests/make-geom.sh" >/dev/null 2>&1 || true

with_geometry "$WORK/does-not-exist.vox" > "$WORK/missing.par"
expect_reject "a missing geometry file is refused" "$WORK/missing.par" "cannot be opened"

printf 'not a voxel volume at all\n' > "$WORK/garbage.vox"
with_geometry "$WORK/garbage.vox" > "$WORK/garbage.par"
expect_reject "a file that is not a volume is refused" "$WORK/garbage.par" "P5V"

# A header that promises more data than the file holds.
printf 'P5V\n64 64 64\n255\n' > "$WORK/short.vox"
head -c 1000 /dev/zero >> "$WORK/short.vox"
with_geometry "$WORK/short.vox" > "$WORK/short.par"
expect_reject "a truncated volume is refused" "$WORK/short.par" "truncated"

# Fewer than the required voxels per cell: an 8^3 volume on an 8^3 grid.
"$ROOT/tools/genvox.py" blank --out "$WORK/coarse.vox" --size 8 8 8 --domain 1 1 1 \
    >/dev/null 2>&1
with_geometry "$WORK/coarse.vox" > "$WORK/coarse.par"
expect_reject "an under-resolved volume is refused" "$WORK/coarse.par" "too coarse"

# A sealed fluid pocket adds a null-space vector the solvers do not carry.
if [ -f "$ROOT/tests/geom/pocket.vox" ]; then
    with_geometry "$ROOT/tests/geom/pocket.vox" > "$WORK/pocket.par"
    sed -i.bak -e 's/^imax 8$/imax 16/' -e 's/^jmax 8$/jmax 16/' -e 's/^kmax 8$/kmax 16/' \
        -e 's/^xlength 1.0$/xlength 4.0/' -e 's/^ylength 1.0$/ylength 4.0/' \
        -e 's/^zlength 1.0$/zlength 4.0/' "$WORK/pocket.par"
    expect_reject "a sealed fluid pocket is refused" "$WORK/pocket.par" "not connected"
fi

# A warning, not a refusal: an anisotropic volume is sampled anyway.
if [ -f "$ROOT/tests/geom/skewed.vox" ]; then
    with_geometry "$ROOT/tests/geom/skewed.vox" > "$WORK/skewed.par"
    sed -i.bak -e 's/^xlength 1.0$/xlength 4.0/' -e 's/^ylength 1.0$/ylength 4.0/' \
        -e 's/^zlength 1.0$/zlength 4.0/' -e 's/^imax 8$/imax 16/' \
        -e 's/^jmax 8$/jmax 16/' -e 's/^kmax 8$/kmax 16/' "$WORK/skewed.par"

    printf -- '-- aspect ratio mismatch warns and continues\n'
    out=$( cd "$WORK" && "$BIN" "$WORK/skewed.par" 2>&1 )
    rc=$?

    if [ $rc -ne 0 ]; then
        echo "aspect ratio mismatch: FAILED, the run was refused rather than warned"
        status=1
    elif ! printf '%s' "$out" | grep -q "aspect ratio mismatch"; then
        echo "aspect ratio mismatch: FAILED, no warning was printed"
        status=1
    elif ! printf '%s' "$out" | grep -q "^TIME "; then
        echo "aspect ratio mismatch: FAILED, the run did not proceed"
        status=1
    else
        echo "aspect ratio mismatch: OK"
    fi
fi

# Unequal smoothing counts leave the multigrid cycle asymmetric, which is a
# cycle the setup did not ask for and which a Krylov method cannot be
# preconditioned by. Refused rather than reconciled.
uneven() {
    base 1 1 | sed -e 's/^presmooth 2$/presmooth 4/' -e 's/^postsmooth 2$/postsmooth 2/'
    printf 'levels 2\n'
}

# Needs a build that actually constructs a hierarchy: the refusal lives where
# the cycle is built, so a relaxation solver -- or CG with a cheap
# preconditioner -- never reaches it and is right not to.
MGBIN="$ROOT/CFD-Solver-SOLVERMG"

if make -C "$ROOT" SOLVER=mg BUILD_DIR=./build/SOLVERMG \
    TARGET="CFD-Solver-SOLVERMG" >/dev/null 2>&1; then

    uneven > "$WORK/uneven-smoothing.par"

    saved_bin=$BIN
    BIN=$MGBIN
    expect_reject "unequal smoothing counts are refused" "$WORK/uneven-smoothing.par" \
        "must be equal"
    BIN=$saved_bin

    rm -rf "$ROOT/build/SOLVERMG" "$MGBIN"
else
    echo "FAILED: could not build the multigrid variant"
    status=1
fi

# The preconditioner selection, which only the conjugate gradient variant acts
# on, so it needs its own binary. An unsupported value has to abort at
# initialization rather than fall back to a default: a benchmark number
# attributed to a preconditioner the run never used is worse than no number.
CGBIN="$ROOT/CFD-Solver-SOLVERCG"

if make -C "$ROOT" SOLVER=cg BUILD_DIR=./build/SOLVERCG \
    TARGET="CFD-Solver-SOLVERCG" >/dev/null 2>&1; then

    with_precon() {
        base 1 1
        printf 'precon %s\n' "$1"
    }

    saved_bin=$BIN
    BIN=$CGBIN

    with_precon "bogus" > "$WORK/precon-bogus.par"
    expect_reject "an unknown preconditioner is refused" "$WORK/precon-bogus.par" \
        "Unsupported preconditioner"

    BIN=$saved_bin

    # And the supported values are accepted, so the check above is about the
    # value and not about the parameter being rejected outright. mg is among
    # them now; it was refused by name until the cycle became symmetric.
    for good in none jacobi mg; do
        with_precon "$good" > "$WORK/precon-$good.par"
        printf -- '-- precon %s is accepted\n' "$good"

        if ( cd "$WORK" && "$CGBIN" "$WORK/precon-$good.par" >/dev/null 2>&1 ); then
            printf 'precon %s is accepted: OK\n' "$good"
        else
            printf 'precon %s is accepted: FAILED, the run was refused\n' "$good"
            status=1
        fi
    done

    rm -rf "$ROOT/build/SOLVERCG" "$CGBIN"
else
    echo "FAILED: could not build the conjugate gradient variant"
    status=1
fi

printf '\n==========================================\n'
if [ $status -eq 0 ]; then
    echo "all rejection checks passed"
else
    echo "SOME REJECTION CHECKS FAILED"
fi
exit $status
