"""
Tests for curve_metrics.py, and the guard that keeps it the only definition.

test_no_duplicate_definitions is the one that earns its place. Six functions
were defined in two modules each and four had diverged; the lag statistic's two
copies disagreed on 96 of 402 curves and would have put 21% in a thesis where
the evidence says 34%. That class of bug is invisible in review and expensive
in print, so it is now a test failure.

    python data/test_curve_metrics.py
"""
import ast
import glob
import os
import sys

import numpy as np

from curve_metrics import (ACCELERATION_SIGMA, INITIAL_WINDOW, LAG_THRESHOLD,
                           QUANTISATION_SIGMA, acceleration, curve_noise,
                           initial_rate, line_fit, line_slope, peak_position,
                           window_size)
from fit_dataset import source_floor
from read_rre import RRE_SIGMA

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


# Names every module is allowed to define for itself: per-module entry points
# and test harness helpers, not shared measurements.
PERMITTED_DUPLICATES = {
    "main", "build", "analyse", "report", "check", "close",
    "test_regressions",
}


def _defined_names(path):
    """Top-level functions and classes a module defines, by name."""
    tree = ast.parse(open(path).read())
    return {node.name: node.lineno for node in tree.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                                 ast.ClassDef))}


def test_no_duplicate_definitions():
    """
    No name may be defined at top level in two modules.

    A shared measurement belongs in curve_metrics and is imported. If this
    fails, do not rename one of the copies -- delete one and import the other,
    or the two will drift apart exactly as they did before.
    """
    print("\nno duplicate definitions")
    seen = {}
    duplicates = {}
    for path in sorted(glob.glob(os.path.join(os.path.dirname(__file__) or ".",
                                              "*.py"))):
        module = os.path.basename(path)[:-3]
        for name, line in _defined_names(path).items():
            if name in PERMITTED_DUPLICATES:
                continue
            if name in seen:
                duplicates.setdefault(name, [seen[name]]).append(f"{module}:{line}")
            else:
                seen[name] = f"{module}:{line}"
    check("no shared name is defined in two modules",
          not duplicates,
          "; ".join(f"{n} in {', '.join(w)}" for n, w in sorted(duplicates.items())))


def test_lag_statistic():
    """
    The published 151/402 has to come out of the canonical implementation.

    It was 136/402 (34%) until 2026-08-31, when the readings moved from the
    .txt exports to the instrument's own .rre files. Nothing about the
    statistic changed; the export rounds to 0.001 AU and that rounding was
    flattening 15 curves' lags below the threshold. See read_rre.py.
    """
    print("\nthe lag statistic")
    from fit_dataset import build_curves
    curves, _ = build_curves()
    positions = np.array([peak_position(c.absorbance, c.times) for c in curves],
                         dtype=float)
    lagging = int(np.nansum(positions > LAG_THRESHOLD))
    check("402 fittable curves", len(curves) == 402, f"got {len(curves)}")
    check("151 of them lag, as MECHANISM.md and FITTING.md report",
          lagging == 151, f"got {lagging}")

    # The shape the statistic is supposed to detect, and the shape it is not.
    times = np.linspace(0, 3000, 300)
    straight = 1e-5 * times
    sigmoid = 0.1 / (1 + np.exp(-(times - 2000) / 200))
    check("a straight line does not lag", peak_position(straight, times) == 0.0)
    check("a sigmoid whose inflection is late does lag",
          peak_position(sigmoid, times) > LAG_THRESHOLD)


def test_noise_and_rate():
    print("\nnoise and initial rate")
    rng = np.random.default_rng(0)
    times = np.linspace(0, 3000, 400)
    check("noise floors at the quantisation sigma on a smooth curve",
          curve_noise(1e-5 * times) == QUANTISATION_SIGMA)
    check("noise floors rather than dividing by nothing on a short curve",
          curve_noise([0.1, 0.2, 0.3]) == QUANTISATION_SIGMA)
    sigma = 0.004
    estimate = curve_noise(1e-5 * times + rng.normal(0, sigma, len(times)))
    check("noise recovers a known sigma through a strong trend",
          abs(estimate - sigma) < 0.25 * sigma, f"{estimate:.5f} vs {sigma}")

    slope, stderr, rms = initial_rate(times, 2e-5 * times)
    check("initial_rate recovers a known slope",
          abs(slope - 2e-5) < 1e-9, f"got {slope}")
    check("initial_rate returns (slope, stderr, rms)",
          np.isfinite([slope, stderr, rms]).all())
    check("initial_rate reads only the leading window",
          window_size(len(times), INITIAL_WINDOW) == int(len(times) * INITIAL_WINDOW))
    check("a floor keeps stderr finite on a perfectly straight window",
          stderr >= 0 and np.isfinite(stderr))
    check("line_slope is line_fit without the intercept",
          line_slope(times, 2e-5 * times) == line_fit(times, 2e-5 * times)[1:])


def test_acceleration():
    """
    The autocatalysis statistic, and the case peak_position gets wrong.

    A lag phase starts flat, and a flat start makes the first point-wise
    gradient a coin flip about zero -- which trips peak_position's
    `slope[0] <= 0` guard and scores the curve as having no lag at all. On the
    in-scope block that guard silences 31 of 96 live curves, including all six
    live curves of exp 142, whose lag is visible by eye. `acceleration` fits
    slopes over blocks instead, so it survives a flat start.
    """
    print("\nthe acceleration statistic")
    times = np.linspace(0, 1000, 80)

    z, _ = acceleration(times, 1e-4 * times)
    check("a straight line does not accelerate", abs(z) < ACCELERATION_SIGMA,
          f"z={z:.2f}")
    z, _ = acceleration(times, 0.05 * (1 - np.exp(-times / 200)))
    check("a saturating curve does not accelerate", z < 0, f"z={z:.2f}")
    z, where = acceleration(times, 0.05 / (1 + np.exp(-(times - 600) / 80)))
    check("a sigmoid accelerates", z > ACCELERATION_SIGMA, f"z={z:.2f}")
    check("and its steepest block sits near the inflection",
          0.4 < where < 0.8, f"where={where:.2f}")

    # A lag phase read at three decimals: exactly flat, then a ramp.
    flat_then_ramp = np.round(np.concatenate(
        [np.zeros(40), 4e-5 * (times[40:] - times[40])]), 3)
    z, _ = acceleration(times, flat_then_ramp)
    check("a flat start does not defeat it", z > ACCELERATION_SIGMA,
          f"z={z:.2f}")
    check("where peak_position's first-point guard does defeat it",
          peak_position(flat_then_ramp, times) == 0.0,
          "guard no longer fires -- peak_position changed, revisit this test")

    check("too few points for two blocks returns nan",
          not np.isfinite(acceleration(times[:6], times[:6] * 1e-4)[0]))


def test_floor_belongs_to_the_source():
    """
    The variance floor is an argument, and it changes the verdict.

    `curve_noise` has taken its floor as an argument since 2026-08-31, but
    `line_fit` hardcoded QUANTISATION_SIGMA until 2026-09-01 -- so every
    standard error in the package, and the acceleration z-score that divides
    by two of them, was floored at the .txt export's 0.001 AU rounding even on
    .rre curves read a thousand times finer. It bound on 52 of the 110 live
    in-scope curves and cost 3 of them their acceleration verdict (48/110 read
    where the instrument says 51/110).

    These checks fail if the floor is ever hardcoded again.
    """
    print("\nthe floor is a property of the source")
    times = np.linspace(0, 1000, 80)
    fine = RRE_SIGMA

    # A signal below the export's quantisation but far above the .rre's: real
    # to the instrument, invisible to the export. This is the regime the whole
    # 2026-08-31 .rre swap was about.
    straight = 1e-7 * times

    _, _, coarse_stderr, _ = line_fit(times, straight)
    _, _, fine_stderr, _ = line_fit(times, straight, fine)
    check("line_fit takes a floor and a smaller one gives a smaller stderr",
          fine_stderr < coarse_stderr,
          f"{fine_stderr:.3e} vs {coarse_stderr:.3e}")
    check("and the ratio is the ratio of the floors, since both are floored",
          abs(coarse_stderr / fine_stderr
              - QUANTISATION_SIGMA / RRE_SIGMA) < 1e-6 * (QUANTISATION_SIGMA / RRE_SIGMA),
          f"ratio {coarse_stderr / fine_stderr:.1f}")

    check("line_fit's default is still the export's floor",
          line_fit(times, straight) == line_fit(times, straight,
                                                QUANTISATION_SIGMA))

    # The verdict itself moves: a rise of a few .rre quanta accelerates when
    # judged against the instrument's own floor and does not when judged
    # against the export's.
    # Total rise 1e-4 AU: a tenth of one export quantum, ~380 .rre quanta.
    ramp = np.concatenate([np.zeros(40), 2e-7 * (times[40:] - times[40])])
    check("a sub-quantum acceleration is invisible at the export's floor",
          acceleration(times, ramp)[0] < ACCELERATION_SIGMA,
          f"z={acceleration(times, ramp)[0]:.2f}")
    check("and visible at the instrument's",
          acceleration(times, ramp, floor=fine)[0] > ACCELERATION_SIGMA,
          f"z={acceleration(times, ramp, floor=fine)[0]:.2f}")

    check("source_floor maps the sources to the two constants",
          source_floor("rre") == RRE_SIGMA
          and source_floor("txt") == QUANTISATION_SIGMA)


if __name__ == "__main__":
    test_no_duplicate_definitions()
    test_lag_statistic()
    test_noise_and_rate()
    test_acceleration()
    test_floor_belongs_to_the_source()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
