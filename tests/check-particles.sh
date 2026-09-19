#!/bin/sh
#
# Particle tracing end to end: the files it writes, and that what it writes does
# not depend on how the domain was divided.
#
# Usage: tests/check-particles.sh
#
#   - a setup that traces writes VTK polydata that names the right number of
#     points and stays inside the domain;
#   - the same setup on 1 and several ranks writes the same set of positions;
#   - a setup that does not configure tracing writes no particle files at all.
#
# Environment:
#   TOOLCHAIN  selects the binary suffix (default CLANG)
#   MPIRUN     mpirun command (default mpirun)
#   RANKS      rank count for the comparison (default 4)

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${TOOLCHAIN:-CLANG}
MPIRUN=${MPIRUN:-mpirun}
RANKS=${RANKS:-4}
BIN="$ROOT/CFD-Solver-$TOOLCHAIN"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

status=0

if ! make -C "$ROOT" ${SOLVER:+SOLVER="$SOLVER"} >/dev/null 2>&1; then
    echo "check-particles.sh: build failed" >&2
    exit 2
fi

# A small traced setup: a channel with a sphere, seeded upstream of it.
cat > "$WORK/traced.par" <<'EOF'
name canal
bcLeft    1
bcRight   3
bcBottom  1
bcTop     1
bcFront   1
bcBack    1
gx 0.0
gy 0.0
gz 0.0
re 100.0
u_init 1.0
v_init 0.0
w_init 0.0
p_init 0.0
xlength 8.0
ylength 4.0
zlength 4.0
imax 48
jmax 24
kmax 24
te 0.4
dt 0.02
tau 0.5
levels 2
presmooth 5
postsmooth 5
itermax 50000
eps 0.0001
omg 1.8
gamma 0.9
geometryFile sphere:2.0,2.0,2.0,0.5
numberOfParticles 300
startTime 0.0
injectTimePeriod 0.1
writeTimePeriod 0.1
x1 0.2
y1 0.5
z1 0.5
x2 0.6
y2 3.5
z2 3.5
EOF

# The same setup with no particle block at all.
grep -vE '^(numberOfParticles|startTime|injectTimePeriod|writeTimePeriod|x1|y1|z1|x2|y2|z2) ' \
    "$WORK/traced.par" > "$WORK/untraced.par"

run() {
    outdir=$1
    ranks=$2
    par=$3

    mkdir -p "$outdir"
    if [ "$ranks" -gt 1 ]; then
        ( cd "$outdir" && "$MPIRUN" -n "$ranks" "$BIN" "$par" > run.log 2>&1 )
    else
        ( cd "$outdir" && "$BIN" "$par" > run.log 2>&1 )
    fi
}

printf -- '-- one rank writes particle files\n'
run "$WORK/r1" 1 "$WORK/traced.par"

files=$(ls "$WORK/r1/vis_files"/particles_*.vtk 2>/dev/null | wc -l | tr -d ' ')

if [ "$files" -lt 2 ]; then
    printf 'FAILED: %s particle files were written, expected several\n' "$files"
    status=1
else
    printf 'wrote %s particle files\n' "$files"
fi

last=$(ls "$WORK/r1/vis_files"/particles_*.vtk 2>/dev/null | tail -1)

if [ -n "$last" ]; then
    # The header has to agree with the body, and every point has to be inside
    # the domain: a stray coordinate would mean a particle escaped or a file was
    # assembled wrongly.
    if python3 - "$last" 8.0 4.0 4.0 <<'PYEOF'
import sys

path, lx, ly, lz = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
lines = open(path).read().splitlines()

if lines[3] != "DATASET POLYDATA":
    sys.exit("not polydata")

n = int(lines[4].split()[1])
pts = [tuple(float(v) for v in lines[5 + i].split()) for i in range(n)]

if len(pts) != n:
    sys.exit("point count does not match the header")

verts = lines[5 + n].split()
if verts[0] != "VERTICES" or int(verts[1]) != n:
    sys.exit("vertex count does not match the header")

for x, y, z in pts:
    if not (0.0 <= x <= lx and 0.0 <= y <= ly and 0.0 <= z <= lz):
        sys.exit(f"point ({x}, {y}, {z}) is outside the domain")

print(f"{n} points, header consistent, all inside the domain")
PYEOF
    then
        :
    else
        printf 'FAILED: %s is not well formed\n' "$last"
        status=1
    fi
fi

printf -- '\n-- the same positions on %s ranks\n' "$RANKS"
run "$WORK/rn" "$RANKS" "$WORK/traced.par"

# Compare the final file as a set of positions. Order is not meaningful: the
# particles are gathered rank by rank, so the same set comes out permuted.
a=$(ls "$WORK/r1/vis_files"/particles_*.vtk 2>/dev/null | tail -1)
b=$(ls "$WORK/rn/vis_files"/particles_*.vtk 2>/dev/null | tail -1)

if [ -z "$a" ] || [ -z "$b" ]; then
    printf 'FAILED: one of the runs wrote no particle file\n'
    status=1
elif python3 - "$a" "$b" <<'PYEOF'
import sys

def points(path):
    lines = open(path).read().splitlines()
    n = int(lines[4].split()[1])
    return [tuple(float(v) for v in lines[5 + i].split()) for i in range(n)]

a, b = points(sys.argv[1]), points(sys.argv[2])

if len(a) != len(b):
    sys.exit(f"{len(a)} particles on one rank against {len(b)} on several")

a.sort()
b.sort()

worst = max((max(abs(p - q) for p, q in zip(pa, pb)) for pa, pb in zip(a, b)),
            default=0.0)

# The flow field itself is only reproducible across rank counts to the solve
# tolerance, and a particle integrates that field, so the positions inherit it.
if worst > 1e-3:
    sys.exit(f"positions differ by {worst:.3e}")

print(f"{len(a)} particles, positions agree to {worst:.3e}")
PYEOF
then
    :
else
    printf 'FAILED: the particles differ between rank counts\n'
    status=1
fi

printf -- '\n-- a setup without a particle block writes none\n'
run "$WORK/none" 1 "$WORK/untraced.par"

if [ -d "$WORK/none/vis_files" ] &&
   [ "$(ls "$WORK/none/vis_files" 2>/dev/null | wc -l | tr -d ' ')" -gt 0 ]; then
    printf 'FAILED: a setup with no particle block wrote particle files\n'
    status=1
else
    printf 'no particle files, as asked\n'
fi

printf '\n==========================================\n'
if [ $status -eq 0 ]; then
    echo "particle tracing behaves"
else
    echo "PARTICLE TRACING DOES NOT BEHAVE"
fi
exit $status
