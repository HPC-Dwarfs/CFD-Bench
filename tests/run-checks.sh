#!/bin/sh
#
# Build and run every check driver in tests/checks/.
#
# Usage: tests/run-checks.sh [-n RANKS] [driver ...]
#
#   -n RANKS   run each driver under mpirun with this many ranks (default 1,
#              which runs the binary directly so a serial build works too)
#   driver     run only the named drivers, e.g. "harness"; default is all
#
# Exits 0 only when every driver exits 0.
#
# Environment:
#   TOOLCHAIN  selects the build directory (default CLANG)
#   MPIRUN     mpirun command (default mpirun)

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${TOOLCHAIN:-CLANG}
MPIRUN=${MPIRUN:-mpirun}
BUILD_DIR="$ROOT/build/$TOOLCHAIN"

ranks=1
while [ $# -gt 0 ]; do
    case $1 in
        -n) ranks=$2; shift 2 ;;
        -n*) ranks=${1#-n}; shift ;;
        --) shift; break ;;
        -*) echo "run-checks.sh: unknown option $1" >&2; exit 2 ;;
        *) break ;;
    esac
done

# Honour a solver selection from the environment, so a caller can run the
# drivers against each variant in turn without the build here resetting it.
if ! make -C "$ROOT" ${SOLVER:+SOLVER="$SOLVER"} tests >/dev/null; then
    echo "run-checks.sh: build failed" >&2
    exit 2
fi

if [ $# -gt 0 ]; then
    drivers=""
    for name in "$@"; do
        bin="$BUILD_DIR/check-$name"
        if [ ! -x "$bin" ]; then
            echo "run-checks.sh: no such check driver: $name" >&2
            exit 2
        fi
        drivers="$drivers $bin"
    done
else
    drivers=$(ls "$BUILD_DIR"/check-* 2>/dev/null)
fi

if [ -z "$drivers" ]; then
    echo "run-checks.sh: no check drivers found in $BUILD_DIR" >&2
    exit 2
fi

status=0
passed=0
failed=0

for bin in $drivers; do
    name=$(basename "$bin" | sed 's/^check-//')
    printf '\n---------- %s (%s rank(s)) ----------\n' "$name" "$ranks"

    if [ "$ranks" -gt 1 ]; then
        "$MPIRUN" -n "$ranks" "$bin"
    else
        "$bin"
    fi

    if [ $? -eq 0 ]; then
        printf '%s: OK\n' "$name"
        passed=$((passed + 1))
    else
        printf '%s: FAILED\n' "$name"
        failed=$((failed + 1))
        status=1
    fi
done

printf '\n==========================================\n'
printf 'check drivers: %d passed, %d failed\n' "$passed" "$failed"
exit $status
