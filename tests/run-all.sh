#!/bin/sh
#
# Every check in one place. Run from anywhere.
#
# Usage: tests/run-all.sh
#
# Environment:
#   TOOLCHAIN  selects the build directory (default CLANG)
#   RANKS      rank count for the multi-rank pass (default 4)

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
RANKS=${RANKS:-4}
status=0

step() {
    label=$1
    shift

    if [ ! -x "$1" ] && [ ! -f "$1" ]; then
        printf '\n========== %s ==========\n' "$label"
        printf '%s: SKIPPED (not present yet)\n' "$label"
        return
    fi

    printf '\n========== %s ==========\n' "$label"

    if "$@"; then
        printf '%s: OK\n' "$label"
    else
        printf '%s: FAILED\n' "$label"
        status=1
    fi
}

# The check drivers link whichever solver variant is selected, and some of them
# exercise one variant's internals, so run them under each in turn.
for solver in rb rbc mg cg; do
    export SOLVER=$solver
    step "check drivers, $solver, 1 rank"       "$ROOT/tests/run-checks.sh"
    step "check drivers, $solver, $RANKS ranks" "$ROOT/tests/run-checks.sh" -n "$RANKS"
done
unset SOLVER

step "rejected inputs"          "$ROOT/tests/check-rejects.sh"
step "geometry vs rank count"   "$ROOT/tests/check-geometry-ranks.sh"
step "solvers vs obstacle"      "$ROOT/tests/check-solver-obstacle.sh"
step "shipped setups start"     "$ROOT/tests/check-setups.sh"
step "particle tracing"         "$ROOT/tests/check-particles.sh"
step "Schaefer-Turek cylinder"  "$ROOT/tests/check-schaefer-turek.sh"
step "baselines unchanged"      "$ROOT/tests/record-baseline.sh" -v

printf '\n==================================================\n'
if [ $status -eq 0 ]; then
    echo "all checks passed"
else
    echo "SOME CHECKS FAILED"
fi
exit $status
