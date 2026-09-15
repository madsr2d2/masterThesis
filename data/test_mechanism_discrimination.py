"""
The gate for `data/mechanism_discrimination.py`: the summaries and their anchors.

Stage C's plan (`PLAN_MECHANISM_DISCRIMINATION.md`) fixes the summary counts and
scatter before any candidate exists, because those are the numbers a candidate
is scored against. Task 1's anchors A1 are the curves and runs per substrate and
buffer, the admitted counts and quantiles of L, E and D, their within-run and
replicate spread. If one differs, stop -- the plan's A1.

    python data/test_mechanism_discrimination.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import mechanism_discrimination as md

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


_SUBSTRATES = ("4OMe-BnOH", "BnOH")
_SUMMARY_TABLES = {}
_SUMMARY_SCATTER_CACHE = {}


def _summary_tables():
    """`summary_table` on both substrates, built once per process."""
    if not _SUMMARY_TABLES:
        for substrate in _SUBSTRATES:
            _SUMMARY_TABLES[substrate] = md.summary_table(substrate)
    return _SUMMARY_TABLES


def _scatter_of(substrate):
    """`summary_scatter` on one substrate, built once per process."""
    if substrate not in _SUMMARY_SCATTER_CACHE:
        _SUMMARY_SCATTER_CACHE[substrate] = md.summary_scatter(substrate)
    return _SUMMARY_SCATTER_CACHE[substrate]


def test_window_slopes_on_a_line():
    print("\nwindow slopes on a line")
    t = np.arange(0, 3600.0 + 1, 60.0)
    y = 2e-5 * t + 0.01
    noise = 1e-3
    slopes = md.window_slopes(t, y, noise)
    for name, _, _ in md.SUMMARY_WINDOWS:
        slope, _ = slopes[name]
        check(f"{name}'s slope is 2e-5",
              abs(slope - 2e-5) <= 1e-10 * 2e-5, f"{slope!r}")
    chosen = t[(t >= 1800.0) & (t <= 2700.0)]
    tc = chosen - chosen.mean()
    expected = 1e-3 / np.sqrt(np.sum(tc * tc))
    _, se_q3 = slopes["Q3"]
    check("se_Q3 is the interval's own width",
          abs(se_q3 - expected) <= 1e-12 * expected, f"{se_q3!r}")


def test_summaries_are_nan_when_not_admitted():
    print("\nsummaries are NaN when not admitted")
    t = np.arange(0, 3600.0 + 1, 60.0)
    y = np.full(t.shape, 0.01)
    got = md.curve_summaries(t, y, 1e-3)
    for name in ("L", "E", "D"):
        check(f"{name} is NaN on a flat line", np.isnan(got[name]),
              f"{got[name]!r}")


def test_summary_anchors():
    print("\nthe summary anchors")
    curves = {
        "4OMe-BnOH": {"curves": 147, "runs": 38,
                      "Boric": (40, 10), "Phosphate": (92, 23),
                      "Pyrophosphate": (15, 5)},
        "BnOH": {"curves": 164, "runs": 31,
                 "Boric": (27, 7), "Phosphate": (19, 5),
                 "Pyrophosphate": (118, 19)},
    }
    tables = _summary_tables()
    for substrate, want in curves.items():
        rows = tables[substrate]["rows"]
        check(f"{substrate} has {want['curves']} curves",
              len(rows) == want["curves"], str(len(rows)))
        check(f"{substrate} has {want['runs']} runs",
              rows.experiment.nunique() == want["runs"],
              str(rows.experiment.nunique()))
        for buffer in ("Boric", "Phosphate", "Pyrophosphate"):
            part = rows[rows.buffer == buffer]
            got = (len(part), part.experiment.nunique())
            check(f"{substrate} {buffer} is {want[buffer]}",
                  got == want[buffer], f"{got}")

    anchors = {
        ("4OMe-BnOH", "L"): (147, -13.227, -10.576, -8.727, 0.008, 0.739, 109),
        ("4OMe-BnOH", "E"): (145, -0.623, -0.025, 0.539, 0.004, 0.279, 107),
        ("4OMe-BnOH", "D"): (146, -0.308, -0.024, 0.204, 0.011, 0.129, 108),
        ("BnOH", "L"): (158, -14.627, -11.978, -10.212, 0.024, 0.919, 127),
        ("BnOH", "E"): (158, -0.867, 0.224, 0.781, 0.012, 0.810, 127),
        ("BnOH", "D"): (156, -0.416, -0.095, 0.435, 0.034, 0.447, 125),
    }
    for (substrate, name), want in anchors.items():
        got = _scatter_of(substrate).loc[name]
        check(f"{substrate} {name} admitted is {want[0]}",
              int(got.admitted) == want[0], str(int(got.admitted)))
        check(f"{substrate} {name} quantiles are {want[1:4]}",
              max(abs(got.q10 - want[1]), abs(got.q50 - want[2]),
                  abs(got.q90 - want[3])) <= 0.001,
              f"{got.q10:.3f} / {got.q50:.3f} / {got.q90:.3f}")
        check(f"{substrate} {name} median_se is {want[4]}",
              abs(got.median_se - want[4]) <= 0.001, f"{got.median_se:.3f}")
        check(f"{substrate} {name} within_run_sd is {want[5]} (df {want[6]})",
              abs(got.within_run_sd - want[5]) <= 0.001
              and int(got.within_run_df) == want[6],
              f"{got.within_run_sd:.3f} ({int(got.within_run_df)})")

    replicate = {"L": (0.236, 12, 16), "E": (0.705, 12, 16),
                 "D": (0.132, 12, 16)}
    scatter = _scatter_of("4OMe-BnOH")
    for name, (sd, df, n) in replicate.items():
        got = scatter.loc[name]
        check(f"4OMe-BnOH {name} replicate SD is {sd} (df {df}, n {n})",
              abs(got.replicate_sd - sd) <= 0.001
              and int(got.replicate_df) == df and int(got.replicate_n) == n,
              f"{got.replicate_sd:.3f} ({int(got.replicate_df)}, "
              f"{int(got.replicate_n)})")
    for name in ("L", "E", "D"):
        got = _scatter_of("BnOH").loc[name]
        check(f"BnOH {name} has no replicate rows", np.isnan(got.replicate_sd)
              and int(got.replicate_df) == 0 and int(got.replicate_n) == 0,
              f"{got.replicate_sd!r} ({int(got.replicate_df)}, "
              f"{int(got.replicate_n)})")


if __name__ == "__main__":
    test_window_slopes_on_a_line()
    test_summaries_are_nan_when_not_admitted()
    test_summary_anchors()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
