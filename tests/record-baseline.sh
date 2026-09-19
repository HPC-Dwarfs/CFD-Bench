#!/bin/sh
#
# Record reference field dumps for the obstacle-free setups under every solver.
#
# Usage: tests/record-baseline.sh [-v] [-o DIR]
#
#   -v       verify instead of record: compare a fresh run against the stored
#            dumps and fail when any differs
#   -o DIR   where the dumps live (default tests/baseline)
#
# The setups in tests/setups/ are the shipped ones at their own resolution with
# a shortened final time; see the comment at the top of each. Every solver is
# recorded because the three do not agree today and the whole point of the
# baseline is to show what the correctness fixes change.
#
# Environment:
#   TOOLCHAIN  selects the build directory (default CLANG)
#   SOLVERS    solvers to record (default "rb rbc mg")
#   SETUPS     setup names to record (default "dcavity canal")

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${TOOLCHAIN:-CLANG}
SOLVERS=${SOLVERS:-"rb rbc mg"}
SETUPS=${SETUPS:-"dcavity canal"}

OUT="$ROOT/tests/baseline"
verify=0

while [ $# -gt 0 ]; do
    case $1 in
        -v) verify=1; shift ;;
        -o) OUT=$2; shift 2 ;;
        *) echo "record-baseline.sh: unknown option $1" >&2; exit 2 ;;
    esac
done

mkdir -p "$OUT"
status=0

for solver in $SOLVERS; do
    printf '\n===== solver %s =====\n' "$solver"

    if ! make -C "$ROOT" SOLVER="$solver" tests >/dev/null; then
        echo "record-baseline.sh: build failed for SOLVER=$solver" >&2
        exit 2
    fi

    for setup in $SETUPS; do
        par="$ROOT/tests/setups/$setup-baseline.par"
        ref="$OUT/$setup-$solver.dump"

        if [ ! -f "$par" ]; then
            echo "record-baseline.sh: missing setup $par" >&2
            exit 2
        fi

        if [ "$verify" -eq 1 ]; then
            if [ ! -f "$ref" ]; then
                echo "record-baseline.sh: no recorded baseline $ref" >&2
                status=1
                continue
            fi
            got="$OUT/$setup-$solver.new"
        else
            got="$ref"
        fi

        ( cd "$ROOT" && NUSIF_FIELD_DUMP="$got" \
            "./CFD-Solver-$TOOLCHAIN-test" "$par" >/dev/null ) || {
            echo "record-baseline.sh: run failed for $setup under $solver" >&2
            status=1
            continue
        }

        if [ "$verify" -eq 1 ]; then
            printf -- '-- %s / %s --\n' "$setup" "$solver"
            if "$ROOT/tools/fieldcmp" "$ref" "$got" 0.0; then
                printf '%s/%s: unchanged\n' "$setup" "$solver"
            else
                printf '%s/%s: CHANGED\n' "$setup" "$solver"
                status=1
            fi
            rm -f "$got"
        else
            printf 'recorded %s\n' "$ref"
        fi
    done
done

if [ "$verify" -eq 1 ]; then
    printf '\n==========================================\n'
    if [ $status -eq 0 ]; then
        echo "baselines unchanged"
    else
        echo "BASELINES CHANGED"
    fi
fi

exit $status
