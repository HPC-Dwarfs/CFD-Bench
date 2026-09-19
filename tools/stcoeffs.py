#!/usr/bin/env python3
"""
Turn the force time series into the coefficients the Schaefer-Turek benchmark
is defined in terms of, and report each against its published range.

    tools/stcoeffs.py forces.dat [--umean 1.0] [--diameter 0.1] [--span 0.41]

The coefficients are

    cD = 2 Fx / (rho Umean^2 D H)
    cL = 2 Fy / (rho Umean^2 D H)

with rho = 1, D the cylinder diameter and H the spanwise extent. The Strouhal
number is taken from the lift signal, St = f D / Umean, with f estimated from
the mean interval between upward zero crossings of the lift about its own mean
over the second half of the record. That is deliberately crude: a record that
has not begun to shed periodically will not produce a frequency, and the script
says so rather than inventing one.

Exits 1 when any quantity is outside its published range and --strict is given.
Without --strict it reports and exits 0, because this is a measurement of what
the solver currently does rather than a promise about it.
"""
import argparse
import sys

# Published ranges for case 3D-2Z.
REFERENCE = {
    "cD_max": (3.2000, 3.3000),
    "cL_max": (0.0020, 0.0040),
    "St": (0.2600, 0.3400),
}


def read_series(path):
    times, fx, fy = [], [], []

    with open(path) as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            times.append(float(parts[0]))
            fx.append(float(parts[1]))
            fy.append(float(parts[2]))

    if len(times) < 4:
        sys.exit(f"stcoeffs: {path} has too few samples to say anything")

    return times, fx, fy


def strouhal(times, cl, diameter, umean):
    """Upward zero crossings of the lift about its mean, over the tail."""
    half = len(times) // 2
    t = times[half:]
    y = cl[half:]

    if len(t) < 8:
        return None, "the record is too short"

    mean = sum(y) / len(y)
    dev = [v - mean for v in y]
    span = max(dev) - min(dev)
    scale = max(abs(v) for v in y) + 1e-300

    # A signal that never moves is not oscillating, it is steady.
    if span < 1e-6 * scale or span == 0.0:
        return None, "the lift does not oscillate, so there is nothing to measure"

    crossings = []
    for i in range(1, len(dev)):
        if dev[i - 1] <= 0.0 < dev[i]:
            # Linear interpolation to the crossing time.
            frac = -dev[i - 1] / (dev[i] - dev[i - 1])
            crossings.append(t[i - 1] + frac * (t[i] - t[i - 1]))

    if len(crossings) < 3:
        return None, f"only {len(crossings)} lift cycles in the record"

    periods = [b - a for a, b in zip(crossings, crossings[1:])]
    period = sum(periods) / len(periods)

    if period <= 0.0:
        return None, "the lift period came out non-positive"

    return (1.0 / period) * diameter / umean, None


def report(name, value, note=None):
    lo, hi = REFERENCE[name]

    if value is None:
        print(f"  {name:8s} not measured   reference {lo:.4f} to {hi:.4f}   ({note})")
        return False

    inside = lo <= value <= hi
    mark = "within" if inside else "OUTSIDE"
    if inside:
        detail = ""
    else:
        off = (value - lo) / lo * 100.0 if value < lo else (value - hi) / hi * 100.0
        detail = f", {off:+.1f}% past the {'low' if value < lo else 'high'} end"

    print(f"  {name:8s} {value:10.4f}   reference {lo:.4f} to {hi:.4f}   "
          f"{mark}{detail}")
    return inside


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("forces")
    ap.add_argument("--umean", type=float, default=1.0)
    ap.add_argument("--diameter", type=float, default=0.1)
    ap.add_argument("--span", type=float, default=0.41)
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args(argv)

    times, fx, fy = read_series(args.forces)

    scale = 2.0 / (args.umean ** 2 * args.diameter * args.span)
    cd = [f * scale for f in fx]
    cl = [f * scale for f in fy]

    # The start-up transient is not what the benchmark describes.
    half = len(times) // 2

    print(f"Schaefer-Turek 3D-2Z, {len(times)} samples over t = {times[0]:.3f} "
          f"to {times[-1]:.3f}")
    print(f"  (coefficients taken over the second half, from t = {times[half]:.3f})")

    st, note = strouhal(times, cl, args.diameter, args.umean)

    ok = True
    ok &= report("cD_max", max(cd[half:]))
    ok &= report("cL_max", max(cl[half:]))
    ok &= report("St", st, note)

    print(f"  cD mean over the tail: {sum(cd[half:]) / len(cd[half:]):.4f}")

    if not ok:
        print("\nAt least one quantity is outside its published range. With binary "
              "apertures the body has a staircase surface and a first-order wall "
              "treatment, which is the known limit of this phase; fractional "
              "apertures are what would close the gap.")

    if args.strict and not ok:
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
