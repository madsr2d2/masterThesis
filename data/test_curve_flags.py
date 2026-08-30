"""
Tests for the progress-curve findings in build_dossier.py.

Three kinds of test live here. The first builds synthetic curves whose defect
is known by construction, so the rules have a validation gate rather than only
agreeing with the archive. The second is a regression test for each way the
previous rule was wrong -- it flagged 33 healthy curves and missed 25 damaged
ones, and each test here fails if that behaviour comes back. The third checks
the rules against the real archive, where the counts are pinned so a silent
change in what the dossier reports cannot pass unnoticed.

    python data/test_curve_flags.py      # from the repository root
"""
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

from build_manifest import KNOWN_EXCLUSIONS
from build_dossier import (BACKTRACK_ABSORBANCE, FLAT_ABSORBANCE,
                           QUANTISATION_SIGMA, curve_backtrack,
                           curve_defects, curve_diagnostics, curve_findings,
                           curve_noise)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def curve(values, points=None, duration_min=60.0):
    """Wraps a value array in the row shape curve_findings expects."""
    values = np.asarray(values, dtype=float)
    return {"sample": 1, "points": points if points is not None else len(values),
            "duration_min": duration_min, "start": values[0], "end": values[-1],
            "max": values.max(), "net": values[-1] - values[0],
            "noise": curve_noise(values), "backtrack": curve_backtrack(values),
            "values": values}


def kinds(row):
    return {("backwards" if "backwards" in m else
             "thrashing" if "backtracks" in m else
             "points" if "points" in m else "flat"): level
            for level, m in curve_findings(row)}


def rising(net=0.3, n=120, offset=0.0):
    return offset + net * (1 - np.exp(-3 * np.linspace(0, 1, n)))


def _spiked(values, width, amplitude=0.15):
    """Adds excursions of a given width, to probe the median filter's edge."""
    values = np.asarray(values, dtype=float).copy()
    for start in range(20, len(values) - width, 40):
        values[start:start + width] += amplitude
    return values


# --- noise ----------------------------------------------------------------

def test_noise():
    print("\nnoise estimate")
    check("a perfectly smooth curve floors at the quantisation sigma",
          curve_noise(rising()) == QUANTISATION_SIGMA)
    check("a constant curve floors too, rather than reporting zero",
          curve_noise(np.zeros(200)) == QUANTISATION_SIGMA)

    rng = np.random.default_rng(0)
    for sigma in (0.002, 0.01):
        estimate = curve_noise(rising() + rng.normal(0, sigma, 120))
        check(f"recovers an injected sigma of {sigma}",
              0.7 * sigma < estimate < 1.4 * sigma, f"got {estimate:.5f}")

    steep = curve_noise(np.linspace(0, 2.0, 200))
    check("a steep but smooth ramp is not mistaken for noise",
          steep == QUANTISATION_SIGMA, f"got {steep:.5f}")


# --- backtracking ---------------------------------------------------------

def test_backtrack():
    print("\nbacktracking")
    check("a monotone curve backtracks zero",
          curve_backtrack(rising()) < 1e-9)
    check("a steep monotone ramp backtracks zero",
          curve_backtrack(np.linspace(0, 2.0, 200)) < 1e-9)

    # One excursion of known amplitude: up 0.1, back down, then on as before.
    values = rising(n=200).copy()
    values[100:110] += 0.1
    amplitude = curve_backtrack(values)
    check("a single 0.1 AU excursion is measured at about 0.1",
          0.08 < amplitude < 0.12, f"got {amplitude:.4f}")

    rng = np.random.default_rng(1)
    noisy = curve_backtrack(rising() + rng.normal(0, 0.001, 120))
    check("ordinary point noise does not accumulate into an excursion",
          noisy < BACKTRACK_ABSORBANCE, f"got {noisy:.4f}")


# --- the rules ------------------------------------------------------------

def test_synthetic():
    print("\nrules, on curves whose defect is known by construction")
    check("a healthy rise has no finding at all",
          kinds(curve(rising())) == {})
    check("a healthy rise from a NEGATIVE baseline is still clean",
          kinds(curve(rising(offset=-0.05))) == {},
          "this is the case the old max<0 rule got wrong 33 times")
    check("a curve that falls is a defect",
          kinds(curve(-rising())).get("backwards") == "defect")
    check("a curve going nowhere is a note, not a defect",
          kinds(curve(rising(net=0.002))).get("flat") == "note",
          "the slowest rung of a titration must survive the dossier")
    check("a flat curve at a negative offset is also only a note",
          kinds(curve(rising(net=0.002, offset=-0.05))).get("flat") == "note")

    # Excursions must be wider than the 5-point median window, or the filter
    # removes them -- which is exactly what it is there to do to point noise.
    thrashing = rising(net=0.1, n=200).copy()
    for start in range(20, 180, 40):
        thrashing[start:start + 12] += 0.15
    check("a curve that thrashes is a defect even though it ends up right",
          kinds(curve(thrashing)).get("thrashing") == "defect",
          "the old endpoint-only rule could not see this")
    check("a spike narrower than the median window is filtered, not flagged",
          "thrashing" not in kinds(curve(_spiked(rising(net=0.1, n=200), 2))),
          "single-point glitches are noise, not a disturbed measurement")
    check("a short curve is a defect",
          kinds(curve(rising(n=10))).get("points") == "defect")

    check("a fall smaller than the flat threshold is not called backwards",
          "backwards" not in kinds(curve(rising(net=-0.004))),
          "a -0.002 drift is not a reaction-direction failure")
    check("a large fall clears both the sigma and the absorbance bar",
          abs(curve(-rising(net=0.3))["net"]) > FLAT_ABSORBANCE)


# --- regressions ----------------------------------------------------------

def test_regressions():
    print("\nregressions against the previous rule")

    def old_rule(row):
        return (row["net"] < -1e-9) or (row["max"] < 0) or (row["points"] < 20)

    healthy_but_offset = curve(rising(net=0.037, offset=-0.044))
    check("exp 67 s2's shape: rises 0.037 from a -0.044 baseline",
          old_rule(healthy_but_offset) and not curve_defects(healthy_but_offset),
          "the old rule called this dead; it is a 60-sigma rise")

    drift = curve(rising(net=0.3) - np.linspace(0, 0.002, 120))
    check("a 0.002 downward drift on a 0.3 rise is not a defect",
          not curve_defects(drift))

    check("the sigma bar alone cannot condemn a curve: a tiny but significant "
          "fall stays a note",
          kinds(curve(np.linspace(0, -0.005, 300))).get("flat") == "note")


# --- the archive ----------------------------------------------------------

EXPECTED_DEFECTS = {
    26: 4,                                           # 10 points over 7 minutes
    13: 1, 55: 1, 80: 1, 128: 2, 130: 2, 131: 2,     # thrashing
    135: 4, 138: 3, 140: 2, 141: 1, 142: 2, 146: 1,
    149: 1, 150: 1, 151: 1,
}

# Taken from the manifest rather than restated, so an exclusion added there
# cannot leave this test asserting against a stale set.
REMOVED = set(KNOWN_EXCLUSIONS)


def test_archive():
    print("\nthe archive")
    import pandas as pd
    manifest = pd.read_csv("data/manifest.csv")

    defects, flat = {}, 0
    for number in sorted(manifest["experiment"].unique()):
        diagnostics, _ = curve_diagnostics(int(number))
        if diagnostics is None:
            continue
        for _, row in diagnostics.iterrows():
            if int(number) in REMOVED:
                continue
            if curve_defects(row):
                defects[int(number)] = defects.get(int(number), 0) + 1
            if any(level == "note" for level, _ in curve_findings(row)):
                flat += 1

    check("the live experiments carrying a curve defect are exactly the "
          "surveyed ones", defects == EXPECTED_DEFECTS,
          f"unexpected {set(defects) - set(EXPECTED_DEFECTS)}, "
          f"missing {set(EXPECTED_DEFECTS) - set(defects)}")

    for number in (50, 85, 58, 77, 78):
        diagnostics, _ = curve_diagnostics(number)
        every = all(curve_defects(row) for _, row in diagnostics.iterrows())
        check(f"exp {number}, excluded for running backwards, still has every "
              f"curve defective", every)

    print(f"        {sum(defects.values())} defective curves across "
          f"{len(defects)} live experiments; {flat} noted flat")


if __name__ == "__main__":
    test_noise()
    test_backtrack()
    test_synthetic()
    test_regressions()
    test_archive()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
