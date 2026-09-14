"""
The gate for `data/rate_laws.py`: the curve-parameter tables and their anchors.

Anchors are the counts the plan fixes before any fitting code trusts the
tables: 311 live catalysed curves, 94 the activation-sink form cannot hold, 217
that remain, and per substrate and element the curves and runs each table
carries. If one differs, stop -- the plan's A1.

    python data/test_rate_laws.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import rate_laws
import scope

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


_SUBSTRATES = ("4OMe-BnOH", "BnOH")
_TABLES = {}


def _tables():
    """`curve_parameters` on both substrates, built once per process."""
    if not _TABLES:
        for substrate in _SUBSTRATES:
            _TABLES[substrate] = rate_laws.curve_parameters(substrate)
    return _TABLES


def test_the_archive_anchors():
    print("\nthe archive anchors")
    frame = scope.frame(scope.archive())
    catalysed = frame[frame.live & (frame.e0 > 0)]
    check("311 live catalysed curves in the frame", len(catalysed) == 311,
          str(len(catalysed)))
    tables = _tables()
    removed = sum(
        int(t["excluded"].query(
            "element == 'v_act' and reason == 'form cannot hold the curve'"
        )["count"].iloc[0]) for t in tables.values())
    check("94 removed as form cannot hold the curve", removed == 94,
          str(removed))
    retained = sum(len(t["rows"]) for t in tables.values())
    check("217 remain", retained == 217, str(retained))


def test_the_substrate_anchors():
    print("\nthe substrate anchors (before the positivity exclusions)")
    tables = _tables()
    expected = {
        "4OMe-BnOH": {"curves": 105, "runs": 34, "v_act": 73, "tau": 81,
                      "k": 44, "bubble_load": 1},
        "BnOH": {"curves": 112, "runs": 31, "v_act": 65, "tau": 74,
                 "k": 27, "bubble_load": 8},
    }
    for substrate, want in expected.items():
        rows = tables[substrate]["rows"]
        got = {
            "curves": len(rows),
            "runs": rows.experiment.nunique(),
            "v_act": int(rows.v_act_resolved.sum()),
            "tau": int(rows.tau_resolved.sum()),
            "k": int((rows.k_state == "resolved").sum()),
            "bubble_load": int((rows.bubble_load > 1).sum()),
        }
        for key, value in want.items():
            check(f"{substrate} {key} is {value}", got[key] == value,
                  str(got[key]))


def test_the_element_table_anchors():
    print("\nthe element tables after every exclusion (curves / runs)")
    tables = _tables()
    expected = {
        "v_act": {"4OMe-BnOH": (73, 32), "BnOH": (65, 29)},
        "k_act_lag": {"4OMe-BnOH": (70, 29), "BnOH": (28, 18)},
        "k_act_burst": {"4OMe-BnOH": (11, 5), "BnOH": (46, 24)},
        "k_sink": {"4OMe-BnOH": (44, 21), "BnOH": (27, 18)},
    }
    for element, per_substrate in expected.items():
        for substrate, (curves, runs) in per_substrate.items():
            table = tables[substrate]["elements"][element]
            got = (len(table), table.experiment.nunique())
            check(f"{substrate} {element} is {curves} / {runs}",
                  got == (curves, runs), str(got))
    for substrate, count in (("4OMe-BnOH", 4), ("BnOH", 10)):
        rows = tables[substrate]["rows"]
        check(f"{substrate} has {count} curves with v_act <= 0",
              int((rows.v_act <= 0).sum()) == count,
              str(int((rows.v_act <= 0).sum())))


def test_lag_and_burst_share_runs():
    print("\nruns carrying both families")
    tables = _tables()
    expected = {
        "4OMe-BnOH": {16, 21, 42, 46, 47, 48, 130},
        "BnOH": {51, 66, 71, 75, 76, 83},
    }
    for substrate, want in expected.items():
        rows = tables[substrate]["rows"]
        both = {int(experiment) for experiment, group in
                rows.groupby("experiment")
                if set(group.family) == {"lag", "burst"}}
        check(f"{substrate} runs include {sorted(want)}", want <= both,
              str(sorted(both)))


def test_se_is_the_interval_width():
    print("\nthe log-scale standard error")
    interval = (np.exp(-1.0), np.exp(1.0))
    check("on (e^-1, e^1) se is 2 / (2 * 1.96)",
          np.isclose(rate_laws._log_se(interval), 2.0 / (2 * 1.96)),
          str(rate_laws._log_se(interval)))


if __name__ == "__main__":
    test_the_archive_anchors()
    test_the_substrate_anchors()
    test_the_element_table_anchors()
    test_lag_and_burst_share_runs()
    test_se_is_the_interval_width()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
