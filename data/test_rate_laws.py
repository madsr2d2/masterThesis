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
import pandas as pd

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


def _synthetic_table(n, buffer=None, experiment=None, **columns):
    """A table with every column `fit_model` and its health checks read."""
    table = pd.DataFrame({
        "buffer": ["Phosphate"] * n,
        "experiment": list(range(1, n + 1)),
        "sample": list(range(1, n + 1)),
        "s0": np.linspace(1.0, 10.0, n),
        "h2o2": np.linspace(10.0, 100.0, n),
        "hoo": np.linspace(0.01, 1.0, n),
        "buf": np.linspace(10.0, 100.0, n),
        "e0": np.linspace(0.01, 0.1, n),
        "kelvin": np.linspace(288.0, 313.0, n),
    })
    if buffer is not None:
        table["buffer"] = buffer
    if experiment is not None:
        table["experiment"] = experiment
    for key, value in columns.items():
        table[key] = value
    table["se"] = 0.1
    table["weight"] = 1.0 / (table.se ** 2 + rate_laws.SE_FLOOR ** 2)
    table["y"] = 0.0
    return table


def _flags(table, model, starts=(-2.0, 0.0, 2.0)):
    return rate_laws.fit_model(table, model, starts=starts)["health"]["flags"]


def test_a_strong_planted_model_is_recovered():
    """A3: the machinery test. The true model must come back, off the truth."""
    print("\nthe strong planted model (A3)")
    table = _tables()["4OMe-BnOH"]["elements"]["v_act"].copy()
    s0 = table.s0.to_numpy(dtype=float)
    hoo = table.hoo.to_numpy(dtype=float)
    buf = table.buf.to_numpy(dtype=float)
    intercept = {"Boric": 0.0, "Phosphate": -0.5, "Pyrophosphate": 0.5}
    planted = (np.array([intercept[b] for b in table.buffer])
               + 0.8 * np.log(s0) + 0.5 * np.log(hoo)
               + np.log(0.02 * buf / (1.0 + 0.02 * buf)))
    generator = np.random.default_rng(0)
    table["y"] = planted + generator.normal(0.0, 0.01, len(table))
    result = rate_laws.fit_model(
        table, {"S": "power", "HOO": "power", "BUF": "bind"},
        starts=(-2.0, 0.0, 2.0))
    coefficients = result["coefficients"]
    check("a_S within 0.03 of 0.8",
          abs(coefficients["a_S"] - 0.8) <= 0.03,
          f"{coefficients['a_S']:.4f}")
    check("a_HOO within 0.03 of 0.5",
          abs(coefficients["a_HOO"] - 0.5) <= 0.03,
          f"{coefficients['a_HOO']:.4f}")
    check("log10 K_B within 0.1 of log10 0.02",
          abs(coefficients["log10(K_B)"] - np.log10(0.02)) <= 0.1,
          f"{coefficients['log10(K_B)']:.4f}")
    flags = result["health"]["flags"]
    check("no at-bound flag",
          not any(flag.startswith("at bound") for flag in flags), str(flags))
    check("no rank-deficient flag", "rank deficient" not in flags, str(flags))


def test_health_flags():
    """One hand-built table per flag, each producing exactly that flag."""
    print("\nthe health flags, planted one at a time")
    at_bound = _synthetic_table(20)
    at_bound["y"] = np.log(at_bound.buf.to_numpy(dtype=float))
    check("at bound",
          _flags(at_bound, {"BUF": "bind"}) == ["at bound: K_B"],
          str(_flags(at_bound, {"BUF": "bind"})))

    collinear = _synthetic_table(20)
    first = np.linspace(-1.0, 1.0, 20)
    second = first + 0.02 * np.sin(np.arange(20))
    collinear["s0"] = np.exp(first)
    collinear["h2o2"] = np.exp(second)
    collinear["y"] = first + 0.01 * np.sin(3.0 * np.arange(20))
    check("collinear",
          _flags(collinear, {"S": "power", "H": "power"})
          == ["collinear: a_S/a_H r=-1.00"],
          str(_flags(collinear, {"S": "power", "H": "power"})))

    deficient = _synthetic_table(20)
    x1 = np.linspace(-1.0, 1.0, 20)
    x2 = (-1.0) ** np.arange(20) * 0.9
    deficient["s0"] = np.exp(x1)
    deficient["h2o2"] = np.exp(x2)
    deficient["hoo"] = np.exp(x1 + x2)
    check("rank deficient",
          _flags(deficient, {"S": "power", "H": "power", "HOO": "power"})
          == ["rank deficient"],
          str(_flags(deficient, {"S": "power", "H": "power",
                                 "HOO": "power"})))

    few_runs = _synthetic_table(20)
    hoo = np.full(20, 1.0)
    hoo[0], hoo[1] = np.exp(0.1), np.exp(-0.1)
    few_runs["hoo"] = hoo
    check("carried by fewer than 3 runs",
          _flags(few_runs, {"HOO": "power"})
          == ["carried by fewer than 3 runs: a_HOO"],
          str(_flags(few_runs, {"HOO": "power"})))

    short = _synthetic_table(10)
    check("too few curves",
          _flags(short, {}) == ["too few curves"],
          str(_flags(short, {})))


def test_identifiability_rules():
    print("\nthe identifiability rules")
    table = _synthetic_table(20)
    table["s0"] = 5.0
    ok, reason = rate_laws.term_identifiable(table, "S", "power")
    check("S:power at constant [S] is not identifiable", not ok)
    check("the reason names the distinct values",
          reason == "fewer than 2 distinct values", reason)


def test_within_run_check_catches_a_between_run_confound():
    print("\nthe within-run check, planted")
    experiments, s0, y = [], [], []
    for run in range(1, 13):
        for step in (0.0, 0.5, 1.0):
            experiments.append(run)
            s0.append(run + step)
            y.append(2.0 * run)
    table = _synthetic_table(len(y), experiment=experiments, s0=s0)
    generator = np.random.default_rng(0)
    table["y"] = np.array(y) + generator.normal(0.0, 0.01, len(y))
    result = rate_laws.within_run_check(table, {"S": "power"})
    row = result["table"].loc["a_S"]
    check("the between-run a_S is non-zero", abs(row.between_runs) > 0.5,
          f"{row.between_runs:.4f}")
    check("the within-run a_S is near zero", abs(row.within_runs) < 0.1,
          f"{row.within_runs:.4f}")
    check("the two disagree", bool(row.disagree))


if __name__ == "__main__":
    test_the_archive_anchors()
    test_the_substrate_anchors()
    test_the_element_table_anchors()
    test_lag_and_burst_share_runs()
    test_se_is_the_interval_width()
    test_a_strong_planted_model_is_recovered()
    test_health_flags()
    test_identifiability_rules()
    test_within_run_check_catches_a_between_run_confound()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
