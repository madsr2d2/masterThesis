"""
The gate for `data/rate_laws.py`: the curve-parameter tables and their anchors.

Anchors are the counts the plan fixes before any fitting code trusts the
tables: 311 live catalysed curves, 94 the activation-sink form cannot hold, 217
that remain, and per substrate and element the curves and runs each table
carries. If one differs, stop -- the plan's A1.

    python data/test_rate_laws.py
"""
import json
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


def test_run_clustered_errors_cover_a_between_run_null():
    """200 planted nulls: the naive error over-rejects, the clustered does not."""
    print("\nthe run-clustered errors, 200 planted nulls")
    naive, clustered = 0, 0
    for seed in range(200):
        generator = np.random.default_rng(seed)
        hoo_run = np.exp(generator.normal(0.0, 1.0, 12))
        offset = generator.normal(0.0, 0.3, 12)
        experiments = np.repeat(np.arange(1, 13), 4)
        noise = generator.normal(0.0, 0.02, 48)
        table = _synthetic_table(48, experiment=experiments)
        table["hoo"] = hoo_run[experiments - 1]
        table["y"] = offset[experiments - 1] + noise
        table["weight"] = 1.0
        result = rate_laws.fit_model(table, {"HOO": "power"})
        coefficient = result["coefficients"]["a_HOO"]
        if abs(coefficient / result["stderr"]["a_HOO"]) > 2.0:
            naive += 1
        if abs(coefficient / result["stderr_clustered"]["a_HOO"]) > 2.0:
            clustered += 1
    check("the naive errors reject more than 22% of the nulls",
          naive / 200 > 0.22, f"{naive}/200")
    check("the clustered errors reject fewer than 17% of the nulls",
          clustered / 200 < 0.17, f"{clustered}/200")


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


def test_temperature_needs_four_temperatures():
    print("\nthe temperature floor")
    table = _synthetic_table(20)
    table["kelvin"] = 273.15 + np.array([25.0, 30.0, 35.0] * 6 + [25.0, 30.0])
    ok, reason = rate_laws.term_identifiable(table, "T", "arrhenius")
    check("three temperatures leave T unidentifiable", not ok)
    check("the reason names the temperature count",
          reason == "fewer than 4 temperatures", reason)


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


def test_within_run_check_never_keeps_between_run_families():
    print("\nthe within-run families")
    experiments = list(np.repeat(np.arange(1, 7), 4))
    hoo = np.tile(np.array([1.0, 1.07, 1.13, 1.20]), 6)
    table = _synthetic_table(24, experiment=experiments, hoo=hoo)
    table["y"] = np.log(table.hoo.to_numpy(dtype=float))
    result = rate_laws.within_run_check(table, {"HOO": "power"})
    check("the within-run design has no HOO coefficient",
          "a_HOO" not in result["table"].index, str(result["table"].index))
    check("HOO is dropped as a between-run family",
          dict(result["dropped"]).get("a_HOO") == "between-run family",
          str(result["dropped"]))


def test_every_run_is_held_out_once_and_never_trained_on():
    print("\nthe leave-one-run-out folds")
    experiments, s0 = [], []
    for run in range(1, 7):
        for step in (0.0, 1.0, 2.0):
            experiments.append(run)
            s0.append(run + step)
    table = _synthetic_table(len(s0), experiment=experiments, s0=s0)
    table["y"] = np.log(table.s0.to_numpy(dtype=float))
    result = rate_laws.cross_validate(table, {"S": "power"})
    scored = set(result["scores"].index)
    check("every run is held out exactly once",
          scored == set(range(1, 7)) and result["skipped"] == 0, str(scored))
    trained_on_self = [run for run in scored
                       if run in result["trained"][run]]
    check("no fold trained on the run it scores", not trained_on_self,
          str(trained_on_self))
    missing = [run for run in scored
               if set(result["trained"][run]) != set(range(1, 7)) - {run}]
    check("every fold trained on every other run", not missing, str(missing))


def test_the_tie_rule():
    print("\nthe tie rule")
    folds = range(10)
    alternate = np.array([(-1.0) ** i for i in folds])
    scores = pd.DataFrame(
        [np.ones(10),
         1.0 + 0.2 + 0.5 * alternate,
         1.0 + 1.0 + 0.5 * alternate],
        index=["best", "tied", "not"], columns=folds)
    tied = rate_laws.tie_set(scores)
    check("the best model and the within-error model are tied",
          set(tied) == {"best", "tied"}, str(tied))
    check("the model outside the error is not tied", "not" not in tied,
          str(tied))


def test_the_tie_rule_needs_both_tests():
    print("\nthe tie rule needs both tests")
    folds = np.arange(30)
    base = 10.0 ** (folds / 5.0)
    scores = pd.DataFrame(
        [base, 1.10 * base,
         base + np.where(folds < 6, 40.0 * base.max(), 0.0),
         base * (1.0 + 0.02 * (-1.0) ** folds)],
        index=["best", "consistent", "concentrated", "equal"],
        columns=folds)
    statistics = rate_laws.tie_statistics(scores)
    tied = rate_laws.tie_set(scores)
    check("a constant proportional excess is not tied",
          "consistent" not in tied, str(tied))
    check("a six-fold absolute excess is not tied",
          "concentrated" not in tied, str(tied))
    check("the raw test is the one that catches it",
          not bool(statistics.loc["concentrated", "raw_kept"]),
          str(statistics.loc["concentrated"].to_dict()))
    check("a within-error wobble is tied", "equal" in tied, str(tied))
    check("the best model is tied", "best" in tied, str(tied))


def test_the_tie_rule_cannot_see_three_fold_failures():
    print("\nthe tie rule's blind spot")
    folds = np.arange(30)
    base = 10.0 ** (folds / 5.0)
    scores = pd.DataFrame(
        [base, base + np.where(folds < 3, 40.0 * base.max(), 0.0)],
        index=["best", "hidden"], columns=folds)
    statistics = rate_laws.tie_statistics(scores)
    tied = rate_laws.tie_set(scores)
    check("a three-fold absolute excess is tied", "hidden" in tied, str(tied))
    check("its worst fold is more than 10 log units out",
          float(statistics.loc["hidden", "worst_fold_log_ratio"]) > 10.0,
          str(statistics.loc["hidden", "worst_fold_log_ratio"]))


def test_term_verdicts():
    print("\nthe verdict strings")
    dropped = pd.DataFrame(
        [{"family": "T", "option": "arrhenius",
          "reason": "fewer than 2 distinct values"}])
    none = {"S": None, "H": None, "HOO": None, "BUF": None, "E": None,
            "T": None}
    cases = (
        ("S: power in every tied model",
         [dict(none, S="power")], dropped),
        ("S: a dependence in every tied model, option undecided (mm, power)",
         [dict(none, S="power"), dict(none, S="mm")], dropped),
        ("S: absent from every tied model", [dict(none)], dropped),
        ("S: undecided", [dict(none, S="power"), dict(none)], dropped),
        ("T: not identifiable on this table (fewer than 2 distinct values)",
         [dict(none)], dropped),
    )
    for expected, tie, table in cases:
        family = expected.split(":")[0]
        verdicts = rate_laws.term_verdicts(tie, table, "v_act")
        check(expected, verdicts[family] == expected, str(verdicts[family]))


def _planted_identifiability(workers=1):
    """
    What the real 4OMe-BnOH design can identify when the truth is known.

    The realistic planted model of the plan: `S:mm` with Km = 2 mM,
    `HOO:power` 0.5, `BUF:bind` with K_B = 0.02 /mM and the Task 3
    intercepts, with noise at each row's own sqrt(se^2 + SE_FLOOR^2), seeds
    0, 1 and 2. The search is run on each; the verdict tables are saved to
    `data/fits/rate_laws/planted_identifiability.json` and read back on later
    runs of the gate (the plan's skip-existing rule for long jobs).
    """
    path = os.path.join(rate_laws.RATE_LAW_DIR, "planted_identifiability.json")
    if os.path.exists(path):
        with open(path) as handle:
            return json.load(handle)
    os.makedirs(rate_laws.RATE_LAW_DIR, exist_ok=True)
    table = _tables()["4OMe-BnOH"]["elements"]["v_act"].copy()
    s0 = table.s0.to_numpy(dtype=float)
    hoo = table.hoo.to_numpy(dtype=float)
    buf = table.buf.to_numpy(dtype=float)
    intercept = {"Boric": 0.0, "Phosphate": -0.5, "Pyrophosphate": 0.5}
    signal = (np.array([intercept[b] for b in table.buffer])
              + np.log(s0 / (2.0 + s0)) + 0.5 * np.log(hoo)
              + np.log(0.02 * buf / (1.0 + 0.02 * buf)))
    spread = np.sqrt(table.se.to_numpy(dtype=float) ** 2
                     + rate_laws.SE_FLOOR ** 2)
    report = {"substrate": "4OMe-BnOH", "element": "v_act", "seeds": {}}
    for seed in (0, 1, 2):
        planted = table.copy()
        planted["y"] = signal + np.random.default_rng(seed).normal(
            0.0, spread)
        found = rate_laws._search_table(planted, "v_act", workers=workers)
        report["seeds"][str(seed)] = {
            "verdicts": found["verdicts"],
            "tie": found["tie"],
            "tie_size": len(found["tie"]),
            "models": found["models"],
        }
    with open(path, "w") as handle:
        json.dump(report, handle, default=float, indent=1, sort_keys=True)
    return report


def test_the_search_finds_a_realistic_planted_model():
    """Recorded, not asserted: the planted identifiability of the design."""
    print("\nthe realistic planted model (recorded, not asserted)")
    report = _planted_identifiability()
    check("the planted record covers the three seeds",
          set(report["seeds"]) == {"0", "1", "2"}, str(set(report["seeds"])))
    for seed, found in sorted(report["seeds"].items()):
        summary = ", ".join(sorted(found["verdicts"].values()))
        print(f"  seed {seed}: tie {found['tie_size']} of {found['models']} "
              f"models")
        print(f"    {summary}")


def _planted_pre_equilibrium(k_rate, k_clock, seed=0):
    """The real 4OMe-BnOH v_act and k_act_lag tables with the pre-equilibrium
    planted: v_act ~ K buf/(1 + K buf), k_act ~ (1 + K buf)."""
    tables = _tables()
    rate = tables["4OMe-BnOH"]["elements"]["v_act"].copy()
    clock = tables["4OMe-BnOH"]["elements"]["k_act_lag"].copy()
    intercept = {"Boric": 0.0, "Phosphate": -0.5, "Pyrophosphate": 0.5}
    generator = np.random.default_rng(seed)
    for table, constant, relax in ((rate, k_rate, False),
                                   (clock, k_clock, True)):
        buf = table.buf.to_numpy(dtype=float)
        base = np.array([intercept[b] for b in table.buffer])
        base = base + (np.log1p(constant * buf) if relax
                       else np.log(constant * buf / (1.0 + constant * buf)))
        table["y"] = base + generator.normal(0.0, 0.01, len(table))
    return rate, clock


def test_shared_binding_on_a_planted_pre_equilibrium():
    print("\nthe shared binding constant, planted")
    rate, clock = _planted_pre_equilibrium(0.03, 0.03)
    found = rate_laws.shared_binding_law("4OMe-BnOH", "BUF",
                                         rate_table=rate, clock_table=clock)
    check("one K ties with separate K",
          found["verdict"] == "one K ties with separate K",
          found["verdict"])
    check("the shared K is within 0.1 decade of 0.03",
          abs(np.log10(found["shared_k"]) - np.log10(0.03)) <= 0.1,
          f"{found['shared_k']:.4f}")


def test_separate_binding_is_detected():
    print("\nseparate binding constants, planted")
    rate, clock = _planted_pre_equilibrium(0.003, 0.3)
    found = rate_laws.shared_binding_law("4OMe-BnOH", "BUF",
                                         rate_table=rate, clock_table=clock)
    check("separate K predicts better",
          found["verdict"] == "separate K predicts better",
          found["verdict"])


def test_shared_binding_at_bound_is_not_identified():
    print("\nthe shared binding constant at a bound")
    tables = _tables()
    rate = tables["4OMe-BnOH"]["elements"]["v_act"].copy()
    clock = tables["4OMe-BnOH"]["elements"]["k_act_lag"].copy()
    intercept = {"Boric": 0.0, "Phosphate": -0.5, "Pyrophosphate": 0.5}
    generator = np.random.default_rng(0)
    rate["y"] = (np.array([intercept[b] for b in rate.buffer])
                 + generator.normal(0.0, 0.01, len(rate)))
    buf = clock.buf.to_numpy(dtype=float)
    clock["y"] = (np.array([intercept[b] for b in clock.buffer])
                  + np.log(buf) + generator.normal(0.0, 0.01, len(clock)))
    found = rate_laws.shared_binding_law("4OMe-BnOH", "BUF",
                                         rate_table=rate, clock_table=clock)
    check("not identified (K at bound)",
          found["verdict"] == "not identified (K at bound)",
          found["verdict"])


if __name__ == "__main__":
    test_the_archive_anchors()
    test_the_substrate_anchors()
    test_the_element_table_anchors()
    test_lag_and_burst_share_runs()
    test_se_is_the_interval_width()
    test_a_strong_planted_model_is_recovered()
    test_run_clustered_errors_cover_a_between_run_null()
    test_health_flags()
    test_identifiability_rules()
    test_temperature_needs_four_temperatures()
    test_within_run_check_catches_a_between_run_confound()
    test_within_run_check_never_keeps_between_run_families()
    test_every_run_is_held_out_once_and_never_trained_on()
    test_the_tie_rule()
    test_the_tie_rule_needs_both_tests()
    test_the_tie_rule_cannot_see_three_fold_failures()
    test_term_verdicts()
    test_the_search_finds_a_realistic_planted_model()
    test_shared_binding_on_a_planted_pre_equilibrium()
    test_separate_binding_is_detected()
    test_shared_binding_at_bound_is_not_identified()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
