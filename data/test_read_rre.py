"""
Tests for read_rre.py -- that the binary is the same measurement as the export.

The .rre now supplies the READINGS for 277 of the 402 fittable curves, so a
misread block would not announce itself: it would quietly change every rate in
the package. These tests pin the one thing that licenses the substitution --
that where both sources exist they agree to within the export's own rounding
step -- and the resolution gain that motivated it.

    python data/test_read_rre.py
"""
import sys

import numpy as np

from curve_metrics import ABSORBANCE_QUANTUM, QUANTISATION_SIGMA
from fit_dataset import PRIMARY_SCOPE, build_curves
from read_rre import RRE_SIGMA, covered, experiment_number, read_all

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def _unfloored_noise(values):
    """curve_noise without its floor -- the floor is what is under test."""
    values = np.asarray(values, dtype=float)
    curvature = values[2:] - 2 * values[1:-1] + values[:-2]
    return 1.4826 * float(np.median(np.abs(curvature))) / np.sqrt(6)


def test_agreement():
    """
    Every in-scope curve must reproduce from the .rre to the export's rounding.

    This is the check the substitution rests on. If it fails, do NOT widen the
    tolerance: the block offset or the sample mapping has moved, and a curve
    silently assigned to the wrong cuvette is worse than no upgrade at all.
    """
    print("\nthe binary agrees with the export")
    curves = [c for c in build_curves()[0] if c.experiment in PRIMARY_SCOPE]
    check("119 in-scope curves", len(curves) == 119, f"got {len(curves)}")

    instrument = read_all()
    worst, compared = 0.0, 0
    for curve in curves:
        block = instrument.get(curve.experiment, {}).get(curve.sample)
        if block is None or len(block) != len(curve.absorbance):
            continue
        drift = np.abs((block - block[0])
                       - (curve.absorbance - curve.absorbance[0])).max()
        worst = max(worst, float(drift))
        compared += 1
    check("all 119 have a matching .rre block", compared == 119, f"got {compared}")
    check("and none drifts by more than the export's 0.001 AU step",
          worst < ABSORBANCE_QUANTUM, f"worst {worst:.5f} AU")


def test_source_selection():
    """The .rre is preferred where it exists, and only where it exists."""
    print("\nwhich source each curve came from")
    curves, _ = build_curves()
    used = {c.source for c in curves}
    check("sources are only 'rre' and 'txt'", used <= {"rre", "txt"}, str(used))
    from_rre = sum(1 for c in curves if c.source == "rre")
    check("277 of 402 fittable curves read from the instrument file",
          from_rre == 277, f"got {from_rre}")
    check("every in-scope curve reads from the instrument file",
          all(c.source == "rre" for c in curves if c.experiment in PRIMARY_SCOPE))
    # A run with no .rre must keep the export rather than losing its curves.
    without = [c for c in curves if c.experiment not in covered()]
    check("runs with no .rre keep their export",
          without and all(c.source == "txt" for c in without),
          f"{len(without)} curves")


def test_resolution():
    """The gain that motivated all of this, stated as numbers."""
    print("\nresolution")
    check("the .rre floor is far below the export's",
          RRE_SIGMA < QUANTISATION_SIGMA / 100,
          f"{RRE_SIGMA:.2e} vs {QUANTISATION_SIGMA:.2e}")

    # The export's rounding, not the instrument, is what zeroes the scatter on
    # most in-scope curves. Pinned because it is the whole reason for the swap.
    scoped = [c for c in build_curves()[0] if c.experiment in PRIMARY_SCOPE]
    real = np.array([_unfloored_noise(c.absorbance) for c in scoped])
    check("no in-scope curve reads zero noise any more",
          (real > 0).all(), f"{int((real == 0).sum())} still do")
    check("and their real noise sits below the export's floor",
          np.median(real) < QUANTISATION_SIGMA,
          f"median {np.median(real):.2e} vs floor {QUANTISATION_SIGMA:.2e}")

    instrument = read_all()
    check("read_all covers the whole scope",
          set(PRIMARY_SCOPE) <= set(instrument),
          f"missing {sorted(set(PRIMARY_SCOPE) - set(instrument))}")
    check("experiment_number parses rate<n>.rre",
          experiment_number("data/Mads/rate139.rre") == 139
          and experiment_number("rate033.rre") == 33
          and experiment_number("5h mads_t.rme") is None)


def test_no_reference_channel():
    """
    The files carry seven sample channels and no baseline trace.

    Recorded as a test because it is the answer to a question that will be
    asked again: the absolute absorbance of the cuvette contents is NOT
    recoverable from these files, so nothing downstream may assume it is.
    """
    print("\nwhat the files do not contain")
    instrument = read_all()
    counts = {len(samples) for samples in instrument.values()}
    check("no run exposes more than eight channels", max(counts) <= 8, str(sorted(counts)))
    raw = open("data/Mads/rate139.rre", "rb").read().lower()
    check("no reference, blank or baseline block",
          not any(tag in raw for tag in (b"\x00ref", b"blank", b"baseline")))
    first = [samples[1][0] for samples in instrument.values() if 1 in samples]
    check("curves start at ~0 AU, i.e. already referenced",
          np.median(np.abs(first)) < 0.05, f"median |A0| {np.median(np.abs(first)):.4f}")


if __name__ == "__main__":
    test_agreement()
    test_source_selection()
    test_resolution()
    test_no_reference_channel()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
