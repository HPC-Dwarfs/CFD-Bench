#!/bin/sh
#
# What every solver has to deliver once the domain contains a body.
#
# Usage: tests/check-solver-obstacle.sh
#
#   - all four solvers converge to the same field;
#   - the cut-cell pass does not cost a disproportionate number of iterations
#     against the same setup without the body;
#   - the answer and the iteration count do not depend on the rank count.
#
# Environment:
#   TOOLCHAIN  selects the binary suffix (default CLANG)
#   MPIRUN     mpirun command (default mpirun)
#   RANKS      rank count for the decomposition comparison (default 8)
#   FACTOR     iteration-count factor allowed against the obstacle-free case
#              (default 3)

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${TOOLCHAIN:-CLANG}
MPIRUN=${MPIRUN:-mpirun}
RANKS=${RANKS:-8}
FACTOR=${FACTOR:-3}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

PAR="$ROOT/testcases/regression/sphere-baseline.par"
FREE="$WORK/sphere-free.par"

# The same setup without the body, for the iteration-count comparison.
grep -v '^geometryFile' "$PAR" > "$FREE"

status=0

# Iterations or cycles the last solve of a run took.
iterations() {
    awk '/took .* (iterations|cycles) to reach/ { for (i = 1; i < NF; i++) if ($(i+1) == "iterations" || $(i+1) == "cycles") n = $i } END { print n }' "$1"
}

# Krylov solvers form their step lengths from global inner products, whose
# floating-point reduction order depends on the decomposition, so a run near the
# tolerance can legitimately take one more or one fewer iteration on a different
# rank count. The converged field requirement is unchanged.
is_krylov() {
    case $1 in
        cg) return 0 ;;
        *) return 1 ;;
    esac
}

printf '========== solvers agree with a body ==========\n'

for solver in rb rbc mg cg; do
    if ! make -C "$ROOT" SOLVER="$solver" >/dev/null 2>&1 ||
       ! make -C "$ROOT" SOLVER="$solver" tests >/dev/null 2>&1; then
        echo "check-solver-obstacle.sh: build failed for $solver" >&2
        exit 2
    fi

    ( cd "$ROOT" && NUSIF_FIELD_DUMP="$WORK/$solver.dump" \
        "./CFD-Bench-$TOOLCHAIN-test" "$PAR" > "$WORK/$solver.log" 2>&1 )

    ( cd "$ROOT" && "./CFD-Bench-$TOOLCHAIN" "$PAR" > "$WORK/$solver-plain.log" 2>&1 )
    ( cd "$ROOT" && "./CFD-Bench-$TOOLCHAIN" "$FREE" > "$WORK/$solver-free.log" 2>&1 )

    withBody=$(iterations "$WORK/$solver-plain.log")
    without=$(iterations "$WORK/$solver-free.log")

    printf '%s: %s with the body, %s without\n' "$solver" "$withBody" "$without"

    if [ -z "$withBody" ] || [ -z "$without" ]; then
        printf '%s: FAILED, could not read an iteration count\n' "$solver"
        status=1
        continue
    fi

    limit=$((without * FACTOR + 10))
    if [ "$withBody" -gt "$limit" ]; then
        printf '%s: FAILED, the body cost %s iterations against %s without it, more than %sx\n' \
            "$solver" "$withBody" "$without" "$FACTOR"
        status=1
    fi
done

# The solve tolerance the setup asks for; fields are compared against it.
TOL=$(sed -n 's/^eps  *\([0-9.eE+-]*\).*/\1/p' "$PAR" | head -1)

for solver in rbc mg cg; do
    printf -- '-- rb against %s\n' "$solver"
    if "$ROOT/tools/fieldcmp" "$WORK/rb.dump" "$WORK/$solver.dump" "$TOL"; then
        printf 'rb and %s agree to the solve tolerance\n' "$solver"
    else
        printf 'rb and %s DIFFER by more than the solve tolerance\n' "$solver"
        status=1
    fi
done

printf '\n========== the answer does not depend on the rank count ==========\n'

for solver in rb rbc mg cg; do
    make -C "$ROOT" SOLVER="$solver" tests >/dev/null 2>&1

    ( cd "$ROOT" && NUSIF_FIELD_DUMP="$WORK/$solver-r1.dump" \
        "./CFD-Bench-$TOOLCHAIN-test" "$PAR" > "$WORK/$solver-r1.log" 2>&1 )
    ( cd "$ROOT" && NUSIF_FIELD_DUMP="$WORK/$solver-rn.dump" \
        "$MPIRUN" -n "$RANKS" "./CFD-Bench-$TOOLCHAIN-test" "$PAR" \
        > "$WORK/$solver-rn.log" 2>&1 )

    one=$(iterations "$WORK/$solver-r1.log")
    many=$(iterations "$WORK/$solver-rn.log")

    printf -- '-- %s: %s on 1 rank, %s on %s ranks\n' "$solver" "$one" "$many" "$RANKS"

    if is_krylov "$solver"; then
        drift=$((one - many))
        [ "$drift" -lt 0 ] && drift=$((-drift))

        if [ "$drift" -gt 1 ]; then
            printf '%s: FAILED, iteration count moved by %s with the rank count, more than the one a Krylov solver is allowed\n' \
                "$solver" "$drift"
            status=1
        fi
    elif [ "$one" != "$many" ]; then
        printf '%s: FAILED, iteration count changed with the rank count\n' "$solver"
        status=1
    fi

    # Judged by the L2 difference, not the largest one. A whole multi-step run
    # cannot be bit-reproducible across rank counts and is not meant to be: the
    # residual norm, the mean the null-space handling removes and the velocity
    # maximum the time step comes from are all global sums, formed in an order
    # the decomposition decides. Those differences start at machine precision
    # and the flow amplifies them, so two runs agree closely almost everywhere
    # and can differ by more at a single point. The equal iteration counts above
    # are what says the solver itself is decomposition-independent.
    if "$ROOT/tools/fieldcmp" "$WORK/$solver-r1.dump" "$WORK/$solver-rn.dump" \
        "$TOL" --l2; then
        printf '%s: same field on 1 and %s ranks (L2 within the solve tolerance)\n' \
            "$solver" "$RANKS"
    else
        printf '%s: FAILED, the field changed with the rank count\n' "$solver"
        status=1
    fi

    # And no single point may be wildly off, which would mean something other
    # than reduction order.
    if ! "$ROOT/tools/fieldcmp" "$WORK/$solver-r1.dump" "$WORK/$solver-rn.dump" \
        "$(awk -v t="$TOL" 'BEGIN { print t * 20 }')" >/dev/null; then
        printf '%s: FAILED, a single point differs by more than twenty times the solve tolerance\n' \
            "$solver"
        status=1
    fi
done

printf '\n========== the compressed layout agrees with the natural one ==========\n'

# The compressed colour-split layout is the serial path of the rbc variant: an
# MPI build of rbc relaxes in the natural layout and so never exercises it. Its
# cut-cell pass addresses PRED and PBLACK through the second index form every
# list entry carries, which is the thing being checked here.
serialBuild() {
    make -C "$ROOT" ENABLE_MPI=false BUILD_DIR=./build/SERIAL SOLVER="$1" \
        TARGET="CFD-Bench-SERIAL-$1" "CFD-Bench-SERIAL-$1-test" >/dev/null 2>&1
}

if serialBuild rb && serialBuild rbc; then
    for solver in rb rbc; do
        ( cd "$ROOT" && NUSIF_FIELD_DUMP="$WORK/serial-$solver.dump" \
            "./CFD-Bench-SERIAL-$solver-test" "$PAR" > "$WORK/serial-$solver.log" 2>&1 )
    done

    a=$(iterations "$WORK/serial-rb.log")
    b=$(iterations "$WORK/serial-rbc.log")
    printf -- '-- serial: rb %s iterations, rbc %s iterations\n' "$a" "$b"

    if "$ROOT/tools/fieldcmp" "$WORK/serial-rb.dump" "$WORK/serial-rbc.dump" "$TOL"; then
        echo "the compressed layout agrees with the natural one"
    else
        echo "FAILED: the compressed cut-cell pass disagrees with the natural layout"
        status=1
    fi

    rm -f "$ROOT/CFD-Bench-SERIAL-rb" "$ROOT/CFD-Bench-SERIAL-rbc" \
        "$ROOT/CFD-Bench-SERIAL-rb-test" "$ROOT/CFD-Bench-SERIAL-rbc-test"
else
    echo "FAILED: could not build the serial variants"
    status=1
fi

printf '\n==========================================\n'
if [ $status -eq 0 ]; then
    echo "every solver handles the body the same way"
else
    echo "SOLVERS DISAGREE ON THE BODY"
fi
exit $status
