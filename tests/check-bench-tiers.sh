#!/bin/sh
#
# The benchmark tiers keep their multigrid hierarchy at the rank counts they are
# meant for.
#
# Usage: tests/check-bench-tiers.sh
#
# A hierarchy that collapses does not fail. It converges more slowly and says
# nothing unusual, so a benchmark run on it measures a smoother and presents it
# as multigrid. Coarsening stops when a *local* extent can no longer be halved,
# so depth is a property of the grid and the decomposition together, and a grid
# that is fine on one rank can have no hierarchy at all on 216.
#
# Each setup in testcases/bench/ states, in comment lines, what it claims:
#
#   # bench-ranks       1 8 27     the rank counts it is meant for
#   # bench-min-levels  4          the depth it keeps at every one of them
#   # bench-bound       grid       what sets that depth: grid, or body
#
# Four levels is the threshold. A tier may claim fewer only by saying its body
# binds it -- coarsening past the level where a body stops being represented
# computes a correction to a different problem, so there the shallower hierarchy
# is the right one and the tier has to say so.
#
# Three parts:
#
#   1. For every tier and every rank count it names, compute the depth from its
#      grid, its levels request and the factorisation MPI_Dims_create gives,
#      and require it to meet the claim. The large tiers target 1728 ranks, so
#      this has to be computed rather than run.
#   2. Check the computation against the solver: run each small tier at 1, 8
#      and 27 ranks and compare the level count it reports with the computed
#      one. This is what keeps part 1 from being a second copy of the
#      coarsening rule that drifts from the first; if multigrid.c's limit
#      changes, this fails.
#   3. Check the check: it has to reject the grids that motivated it.
#
# Deliberately not part of tests/run-all.sh. Everything there checks what the
# code does; this checks configuration, and part 2 runs 27 ranks, which is more
# than a quick pass should ask of a laptop. Run it when a tier or the coarsening
# rule changes.
#
# Environment:
#   TOOLCHAIN     selects the build directory (default CLANG)
#   MPIRUN        mpirun command (default mpirun)
#   OVERSUBSCRIBE flag letting mpirun start more ranks than cores
#                 (default --oversubscribe, Open MPI's spelling; set it empty
#                 for an MPI that does not need or know it)

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${TOOLCHAIN:-CLANG}
MPIRUN=${MPIRUN:-mpirun}
OVERSUBSCRIBE=${OVERSUBSCRIBE---oversubscribe}
BENCH="$ROOT/testcases/bench"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

status=0

# depth <setup> <ranks>: the number of levels the solver builds for this setup
# at this rank count, or a claim check with --check.
#
# The factorisation is the balanced one MPI_Dims_create produces -- prime
# factors, largest first, each given to the dimension with the smallest product
# so far, then sorted non-increasing -- and the largest factor goes to x, which
# is how commPartition applies it. A rank's extent is N / d, plus one on the
# first N % d ranks. A level exists while every local extent on every rank is
# even and halves to at least 2.
depth() {
    python3 - "$@" <<'EOF'
import sys

def factorise(n):
    primes, p = [], 2
    while n > 1:
        while n % p == 0:
            primes.append(p)
            n //= p
        p += 1
    return sorted(primes, reverse=True)

def dims_create(ranks):
    dims = [1, 1, 1]
    for p in factorise(ranks):
        dims[dims.index(min(dims))] *= p
    return sorted(dims, reverse=True)

def read(path):
    params, claims = {}, {}
    for line in open(path):
        words = line.split()
        if len(words) >= 2 and words[0] == "#" and words[1].startswith("bench-"):
            claims[words[1]] = words[2:]
        elif words and not words[0].startswith("#"):
            params[words[0]] = words[1] if len(words) > 1 else ""
    return params, claims

def levels(params, ranks):
    extents = [int(params["imax"]), int(params["jmax"]), int(params["kmax"])]
    local = [
        [n // d + (1 if n % d > r else 0) for r in range(d)]
        for n, d in zip(extents, dims_create(ranks))
    ]
    built, requested = 1, int(params.get("levels", 5))
    while built < requested:
        sizes = [s for direction in local for s in direction]
        if not all(s % 2 == 0 and s // 2 >= 2 for s in sizes):
            break
        local = [[s // 2 for s in direction] for direction in local]
        built += 1
    return built

path = sys.argv[1]
params, claims = read(path)

if sys.argv[2] != "--check":
    print(levels(params, int(sys.argv[2])))
    sys.exit(0)

missing = [k for k in ("bench-ranks", "bench-min-levels", "bench-bound") if k not in claims]
if missing:
    print("missing " + ", ".join(missing))
    sys.exit(1)

claim = int(claims["bench-min-levels"][0])
bound = claims["bench-bound"][0]
failed = False

if bound not in ("grid", "body"):
    print(f"bench-bound is {bound}, expected grid or body")
    sys.exit(1)

if claim < 4 and bound != "body":
    print(f"claims {claim} levels; fewer than 4 needs the body to be what binds it")
    failed = True

for ranks in map(int, claims["bench-ranks"]):
    got = levels(params, ranks)
    dims = "x".join(map(str, dims_create(ranks)))
    mark = "ok" if got >= claim else "BELOW CLAIM"
    print(f"   {ranks:5d} ranks ({dims:>8}): {got} levels  {mark}")
    failed |= got < claim

sys.exit(1 if failed else 0)
EOF
}

printf -- '-- 1. every tier keeps the depth it claims at the rank counts it names\n'

tiers=$(ls "$BENCH"/*.par 2>/dev/null)

if [ -z "$tiers" ]; then
    echo "check-bench-tiers.sh: no setups in $BENCH" >&2
    exit 2
fi

for par in $tiers; do
    name=$(basename "$par" .par)
    out=$(depth "$par" --check)
    rc=$?

    if [ $rc -eq 0 ]; then
        printf '%s: OK\n%s\n' "$name" "$out"
    else
        printf '%s: FAILED\n%s\n' "$name" "$out"
        status=1
    fi
done

printf -- '\n-- 2. the computed depth is the depth the solver builds\n'

# Only multigrid reports its level count, so this needs a multigrid build.
MGBIN="$ROOT/CFD-Bench-BENCHMG"

if make -C "$ROOT" SOLVER=mg BUILD_DIR=./build/BENCHMG TARGET="CFD-Bench-BENCHMG" \
    >/dev/null 2>&1; then

    for par in "$BENCH"/*-small.par; do
        name=$(basename "$par" .par)

        # One step and one cycle: the hierarchy is built at initialization, so
        # nothing after the first solve changes the count.
        sed -e 's/^te .*/te 0.0/' -e 's/^itermax .*/itermax 1/' "$par" \
            > "$WORK/$name.par"

        for ranks in 1 8 27; do
            want=$(depth "$par" "$ranks")
            out=$( cd "$WORK" && $MPIRUN $OVERSUBSCRIBE -n "$ranks" "$MGBIN" \
                "$WORK/$name.par" 2>&1 )
            got=$(printf '%s' "$out" | sed -n 's/^Using Multigrid solver with \([0-9]*\) levels$/\1/p')

            if [ -z "$got" ]; then
                printf '%s at %s ranks: FAILED, the solver reported no level count\n' \
                    "$name" "$ranks"
                printf '%s\n' "$out" | tail -3
                status=1
            elif [ "$got" != "$want" ]; then
                printf '%s at %s ranks: FAILED, computed %s levels, the solver built %s\n' \
                    "$name" "$ranks" "$want" "$got"
                status=1
            else
                printf '%s at %s ranks: OK, %s levels computed and built\n' \
                    "$name" "$ranks" "$got"
            fi

            # A body-bound tier is only right if the body really is still there
            # at the depth it was bounded to.
            if printf '%s' "$out" | grep -q "unresolved at level"; then
                printf '%s at %s ranks: FAILED, the body is not represented on every level\n' \
                    "$name" "$ranks"
                status=1
            fi
        done
    done

    rm -rf "$ROOT/build/BENCHMG" "$MGBIN"
else
    echo "FAILED: could not build the multigrid variant"
    status=1
fi

printf -- '\n-- 3. the check rejects the grids it exists to catch\n'

# $1 label, $2 grid "imax jmax kmax", $3 ranks, $4 claimed levels, $5 bound
expect_rejected() {
    set -- "$1" $2 "$3" "$4" "$5"
    cat > "$WORK/bad.par" <<EOF
# bench-ranks       $5
# bench-min-levels  $6
# bench-bound       $7
imax $2
jmax $3
kmax $4
levels 8
EOF

    if depth "$WORK/bad.par" --check >/dev/null; then
        printf '%s: FAILED, the check accepted it\n' "$1"
        status=1
    else
        printf '%s: OK, rejected\n' "$1"
    fi
}

# What motivated the check: MPI_Dims_create divides 216 ranks 6x6x6, and 128 / 6
# is 21 -- odd, so there is no hierarchy at all.
expect_rejected "a 128^3 grid at 216 ranks" "128 128 128" 216 4 grid

# karman.par's own grid: 50 halves to 25, which is odd, so it has two levels on
# one rank whatever it asks for.
expect_rejected "an extent that halves to an odd number" "200 50 50" 1 4 grid

# Claiming fewer than four levels without saying the body binds is not allowed,
# even where the grid meets the smaller claim.
expect_rejected "a shallow claim with no body to justify it" "48 48 48" 27 3 grid

printf '\n==========================================\n'
if [ $status -eq 0 ]; then
    echo "every benchmark tier keeps its hierarchy"
else
    echo "SOME BENCHMARK TIERS DO NOT KEEP THEIR HIERARCHY"
fi
exit $status
