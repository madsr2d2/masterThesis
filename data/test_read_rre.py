"""
Tests for read_rre.py -- that the binary is the same measurement as the export.

The .rre now supplies the READINGS for all 402 fittable curves, so a
misread block would not announce itself: it would quietly change every rate in
the package. These tests pin the one thing that licenses the substitution --
that where both sources exist they agree to within the export's own rounding
step -- and the resolution gain that motivated it.

    python data/test_read_rre.py
"""
import glob
import os
import re
import sys

import numpy as np

from curve_metrics import ABSORBANCE_QUANTUM, QUANTISATION_SIGMA
from fit_dataset import DROP_FIRST_READING, PRIMARY_SCOPE, build_curves
from read_rre import (ARCHIVE_DIR, RRE_SIGMA, covered,
                      experiment_number, read_all)

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
        if block is not None and DROP_FIRST_READING:
            # build_curves discards the first reading of every run, so the
            # curve is one shorter than the block it came from. Comparing the
            # raw lengths would silently compare nothing -- as it did, scoring
            # 0 of 119 -- so the block is trimmed the same way.
            block = block[1:]
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
    check("all 402 fittable curves read from the instrument file",
          from_rre == len(curves), f"got {from_rre} of {len(curves)}")
    check("every in-scope curve reads from the instrument file",
          all(c.source == "rre" for c in curves if c.experiment in PRIMARY_SCOPE))
    # THIS BLOCK USED TO ASSERT THE BUG. Until 2026-09-02 it read "28 cuvettes
    # have no block in their run's binary" and explained that exp 6 holds
    # Sample001, 002 and 004 and no 003. The block was there all along: the
    # instrument wrote sample 3's label lowercase, `sample003`, in 31 files
    # across exps 1-32, and `read_rre`'s pattern was case-sensitive. A missing
    # cuvette was observed, rationalised as absent data, and frozen into a
    # check -- which is worse than not checking, because it made the defect
    # look adjudicated. The assertion is now that NOTHING falls back.
    available = covered()
    without = [c for c in curves if c.experiment not in available]
    check("every fittable experiment now has a readable .rre",
          not without, f"{len(without)} curves in uncovered runs")
    fallbacks = [c for c in curves if c.source == "txt"]
    check("no cuvette falls back to the export",
          not fallbacks,
          f"{len(fallbacks)} on .txt: "
          f"{sorted({c.experiment for c in fallbacks})}")
    # And the label the pattern has to keep matching, in both cases. If a
    # future edit tightens it back to `Sample00\\d`, 28 curves go quietly back
    # onto the rounded export and no other check would notice.
    lowercase = [path for path in glob.glob(os.path.join(ARCHIVE_DIR, "*.rre"))
                 if re.search(rb"sample00\d",
                              open(path, "rb").read())]
    check("the lowercase sample label is still matched",
          len(lowercase) >= 25 and not fallbacks,
          f"{len(lowercase)} files carry a lowercase label")


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


def test_both_naming_conventions_are_read():
    """
    Both `rate<n>.rre` and `mads_t<n>.rre` must map to their experiment.

    Only the first was matched until 2026-09-01, so 32 files covering exps
    2-32 were never opened and 97 curves kept the .txt export's floor -- 1096x
    coarser than the instrument's own -- for no reason but a regex.

    The `mads_t` pairing is not asserted from the filename. It is licensed by
    `test_agreement`'s rule, applied per sample: `mads_t003.rre` reproduces all
    seven of exp 3's exported cuvettes to within the export's rounding, at 227
    points each. This test pins the mapping AND that agreement, so a future
    edit cannot quietly re-break either half.
    """
    print("\nboth .rre naming conventions are read")
    check("rate<n>.rre maps to n", experiment_number("rate148.rre") == 148)
    check("mads_t<n>.rre maps to n", experiment_number("mads_t003.rre") == 3,
          f"got {experiment_number('mads_t003.rre')}")
    check("leading zeros are stripped",
          experiment_number("mads_t032.rre") == 32)
    check("an unrelated file maps to nothing",
          experiment_number("5h mads_t.rme") is None)

    available = covered()
    check("exps 2-32 now have a readable .rre",
          {3, 6, 23, 26, 32} <= available,
          f"missing {sorted({3, 6, 23, 26, 32} - available)}")

    # Exp 3 is the archive's widest enzyme-free buffer titration and was
    # entirely on the export's floor. All seven cuvettes must come from the
    # instrument now, which only happens if every one passed the agreement test.
    curves, _ = build_curves()
    exp3 = [c for c in curves if c.experiment == 3]
    check("all seven of exp 3's cuvettes now read from the instrument",
          len(exp3) == 7 and all(c.source == "rre" for c in exp3),
          f"{len(exp3)} curves, sources {sorted({c.source for c in exp3})}")
    # The gain is not a smaller number, it is a MEASURED one. Exp 3's real
    # scatter is about 3.0e-4, just above the export's 2.887e-4 floor, so four
    # of its seven cuvettes previously reported the floor exactly -- a constant
    # standing in for a measurement. None does now.
    check("no exp 3 cuvette reports the export's floor as its noise",
          not any(abs(c.noise - QUANTISATION_SIGMA) < 1e-12 for c in exp3),
          f"noises {[f'{c.noise:.2e}' for c in exp3]}")

    # The .rre may only ever ADD resolution. A curve that had one must not
    # lose it, and the in-scope block must be untouched by this change.
    in_scope = [c for c in curves if c.experiment in PRIMARY_SCOPE]
    check("the in-scope block is still entirely .rre",
          len(in_scope) == 119 and all(c.source == "rre" for c in in_scope),
          f"{len(in_scope)} curves")


if __name__ == "__main__":
    test_agreement()
    test_source_selection()
    test_resolution()
    test_no_reference_channel()
    test_both_naming_conventions_are_read()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
