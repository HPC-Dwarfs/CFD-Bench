#!/usr/bin/env python3
"""
Report statistics of one field over a box of cells in a dump written by
src/fielddump.c.

    tools/fieldprobe.py <dump> --field v --box I0 I1 J0 J1 K0 K1

Indices are global and zero-based, inclusive at both ends; a negative index
counts back from the far side, so -1 is the last cell. Omitting --box covers
the whole domain.

Prints one line of "min max absmax mean", which is enough for a shell check to
assert that a wake exists or that a body stops the flow.
"""
import argparse
import struct
import sys

FIELDS = ("p", "u", "v", "w")


def read_dump(path):
    with open(path, "rb") as fp:
        magic = fp.read(8)
        if magic != b"NUSIFD01":
            sys.exit(f"fieldprobe: {path} is not a NUSIFD01 dump")

        nx, ny, nz = struct.unpack("<3i", fp.read(12))
        count = nx * ny * nz
        data = {}

        for name in FIELDS:
            raw = fp.read(count * 8)
            if len(raw) != count * 8:
                sys.exit(f"fieldprobe: {path} is truncated in field {name}")
            data[name] = struct.unpack(f"<{count}d", raw)

    return nx, ny, nz, data


def resolve(value, n):
    return value + n if value < 0 else value


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dump")
    ap.add_argument("--field", choices=FIELDS, default="u")
    ap.add_argument("--box", nargs=6, type=int,
                    metavar=("I0", "I1", "J0", "J1", "K0", "K1"))
    args = ap.parse_args(argv)

    nx, ny, nz, data = read_dump(args.dump)

    if args.box:
        i0, i1, j0, j1, k0, k1 = args.box
        i0, i1 = resolve(i0, nx), resolve(i1, nx)
        j0, j1 = resolve(j0, ny), resolve(j1, ny)
        k0, k1 = resolve(k0, nz), resolve(k1, nz)
    else:
        i0, i1, j0, j1, k0, k1 = 0, nx - 1, 0, ny - 1, 0, nz - 1

    for lo, hi, n, name in ((i0, i1, nx, "i"), (j0, j1, ny, "j"), (k0, k1, nz, "k")):
        if not (0 <= lo <= hi < n):
            sys.exit(f"fieldprobe: {name} range {lo}..{hi} is outside 0..{n - 1}")

    field = data[args.field]
    lo = hi = None
    absmax = 0.0
    total = 0.0
    count = 0

    for k in range(k0, k1 + 1):
        for j in range(j0, j1 + 1):
            base = (k * ny + j) * nx
            for i in range(i0, i1 + 1):
                v = field[base + i]
                lo = v if lo is None or v < lo else lo
                hi = v if hi is None or v > hi else hi
                absmax = max(absmax, abs(v))
                total += v
                count += 1

    print(f"{lo!r} {hi!r} {absmax!r} {total / count!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
