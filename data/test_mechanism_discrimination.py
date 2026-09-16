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
import pandas as pd

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


def test_candidate_curves_match_the_ode():
    print("\nthe candidate curve matches the ODE")
    from scipy.integrate import solve_ivp
    t = np.linspace(0.0, 6000.0, 101)
    rows = pd.DataFrame([{
        "buffer": "Phosphate", "temperature": 25.0, "pH": 4.0, "s0": 1.0,
        "h2o2": 1.0, "hoo": 1.0, "buf": 0.1, "e0": 1.0}])

    def parameters(candidate, V, tau, theta0, theta_ss, k, X):
        """Parameters that make the rate expressions return these targets."""
        q = theta_ss / tau
        k_r = (1.0 - theta_ss) / tau
        base = 1.0 / (1.0 + 10.0 ** (4.0 - 14.0))
        Y = 1.0 / (1.0 + 1e-8)
        lk_cat = np.log10(V / (1.0 / 1.01 * Y))
        out = {"lk_cat[Phosphate]": lk_cat, "lK_S": -2.0, "lK_O": -8.0,
               "lk_f": np.log10(q / (X * base)), "lk_r": np.log10(k_r),
               "pKa": 14.0, "theta0": theta0}
        out["lk_d" if candidate == "C4" else "lk_s"] = np.log10(k)
        return out

    cases = (
        ("C0", 2e-5, 1500.0, 0.1, 0.6, 2e-4, 1.0),
        ("C0", 2e-5, 900.0, 0.9, 0.2, 1e-9, 1.0),
        ("C4", 2e-5, 1500.0, 0.8, 0.3, 2e-4, 0.1),
    )
    for candidate, V, tau, theta0, theta_ss, k, X in cases:
        curves = md.candidate_curves(
            candidate, parameters(candidate, V, tau, theta0, theta_ss, k, X),
            rows, [t])
        for label, got, want in (
                ("V", curves["V"][0], V), ("tau", curves["tau"][0], tau),
                ("theta_ss", curves["theta_ss"][0], theta_ss),
                ("k", curves["k"][0], k)):
            check(f"{candidate} {label} is {want:g}",
                  abs(got - want) <= 1e-9 * want, f"{got!r}")
        if candidate == "C4":
            def rhs(time, state):
                theta, _ = state
                return [(theta_ss - theta) / tau, V * theta * np.exp(-k * time)]
        else:
            def rhs(time, state):
                theta, A = state
                return [(theta_ss - theta) / tau, V * theta - k * A]
        solution = solve_ivp(rhs, (t[0], t[-1]), [theta0, 0.0], t_eval=t,
                             rtol=1e-11, atol=1e-16)
        got = curves["A"][0][1:]
        want = solution.y[1][1:]
        worst = float(np.max(np.abs(got - want)
                             / np.maximum(np.abs(want), 1e-30)))
        check(f"{candidate} A(t) matches the ODE to 1e-8", worst <= 1e-8,
              f"{worst:.2e}")


def test_c5_matches_the_ode():
    print("\nthe C5 curve matches the ODE")
    from scipy.integrate import solve_ivp
    t = np.linspace(0.0, 6000.0, 601)
    V_base, tau, theta0, theta_ss = 2e-5, 1500.0, 0.1, 0.6
    k_ox, lK_O, hoo = 2e-4, -3.0, 1e-2
    rows = pd.DataFrame([{
        "buffer": "Phosphate", "temperature": 25.0, "pH": 4.0, "s0": 1.0,
        "h2o2": 1.0, "hoo": hoo, "buf": 0.1, "e0": 1.0}])
    q = theta_ss / tau
    k_r = (1.0 - theta_ss) / tau
    base = 1.0 / (1.0 + 10.0 ** (4.0 - 14.0))
    parameters = {
        "lk_cat[Phosphate]": np.log10(V_base / (1.0 / 1.01)),
        "lK_S": -2.0, "lK_O": lK_O, "lk_f": np.log10(q / (0.1 * base)),
        "lk_r": np.log10(k_r), "pKa": 14.0, "theta0": theta0,
        "lk_ox": np.log10(k_ox), "Ea_cat": 0.0, "Ea_act": 0.0, "Ea_ox": 0.0,
    }
    curves = md.candidate_curves("C5", parameters, rows, [t])

    def rhs(time, state):
        decayed = hoo * np.exp(-k_ox * time)
        Y = decayed / (10.0 ** lK_O + decayed)
        theta = theta_ss + (theta0 - theta_ss) * np.exp(-time / tau)
        return [V_base * Y * theta]

    solution = solve_ivp(rhs, (t[0], t[-1]), [0.0], t_eval=t,
                         rtol=1e-11, atol=1e-16)
    got = curves["A"][0][1:]
    want = solution.y[0][1:]
    worst = float(np.max(np.abs(got - want)
                         / np.maximum(np.abs(want), 1e-30)))
    check("C5 A(t) matches the ODE to 1e-4", worst <= 1e-4, f"{worst:.2e}")


def test_run_likelihood_matches_the_dense_normal():
    print("\nthe run likelihood matches the dense normal")
    from scipy.stats import multivariate_normal
    rng = np.random.default_rng(0)
    e = rng.normal(size=5)
    se = rng.uniform(0.1, 0.3, 5)
    runs = np.array([1, 1, 1, 2, 2])
    covariance = (np.diag(se ** 2 + 0.4 ** 2)
                  + 0.7 ** 2 * (runs[:, None] == runs[None, :]))
    dense = -multivariate_normal(
        mean=np.zeros(5), cov=covariance).logpdf(e)
    series = float(md.run_likelihood(e, se, 0.4, 0.7, runs).sum())
    check("the Series' sum equals the dense normal's negative logpdf",
          abs(series - dense) <= 1e-10, f"{series!r} vs {dense!r}")
    check("and both are 3.9404970364222",
          abs(series - 3.9404970364222) <= 1e-10, f"{series!r}")


def test_midpoint_likelihood_anchors():
    print("\nthe midpoint likelihood anchors")
    anchors = {
        "4OMe-BnOH": {"C0": 4873.435372, "C1": 4730.711800,
                      "C2": 6297.876950, "C3": 3082.772937,
                      "C4": 4730.544204, "C5": 4829.060949},
        "BnOH": {"C0": 15622.845649, "C1": 15005.891228,
                 "C2": 15975.468150, "C3": 14998.498241,
                 "C4": 15005.932644, "C5": 15031.666318},
    }
    tables = _summary_tables()
    for substrate, per_candidate in anchors.items():
        table = tables[substrate]
        buffers = sorted(set(table["rows"].buffer))
        for candidate, want in per_candidate.items():
            lower, upper = md.candidate_bounds(candidate, substrate, buffers)
            got = md.candidate_nll(candidate, 0.5 * (lower + upper), table)
            check(f"{substrate} {candidate} midpoint NLL is {want}",
                  abs(got - want) <= 1e-6, f"{got:.6f}")


def test_a6_discrimination_anchors():
    print("\nthe A6 discrimination anchors, five original candidates")
    five = ("C0", "C1", "C2", "C3", "C4")
    status = {
        "4OMe-BnOH": {
            "C0": ("best", "tied", "excluded", "excluded", "excluded"),
            "C1": ("tied", "best", "excluded", "excluded", "excluded"),
            "C2": ("tied", "tied", "best", "tied", "tied"),
            "C3": ("tied", "tied", "tied", "best", "tied"),
            "C4": ("tied", "excluded", "tied", "excluded", "best"),
        },
        "BnOH": {
            "C0": ("excluded", "tied", "tied", "tied", "best"),
            "C1": ("tied", "tied", "best", "tied", "tied"),
            "C2": ("tied", "tied", "best", "tied", "tied"),
            "C3": ("tied", "tied", "best", "tied", "tied"),
            "C4": ("tied", "tied", "best", "tied", "tied"),
        },
    }
    recovered = {
        "4OMe-BnOH": {candidate: "truth recovered" for candidate in five},
        "BnOH": {"C0": "truth not recovered",
                 "C1": "truth recovered", "C2": "truth recovered",
                 "C3": "truth recovered", "C4": "truth recovered"},
    }
    four_ome_pairs = {
        "C0 vs C1": "not distinguishable",
        "C0 vs C2": "one-way (C0 as truth excludes C2)",
        "C0 vs C3": "one-way (C0 as truth excludes C3)",
        "C0 vs C4": "one-way (C0 as truth excludes C4)",
        "C1 vs C2": "one-way (C1 as truth excludes C2)",
        "C1 vs C3": "one-way (C1 as truth excludes C3)",
        "C1 vs C4": "distinguishable",
        "C2 vs C3": "not distinguishable",
        "C2 vs C4": "not distinguishable",
        "C3 vs C4": "one-way (C4 as truth excludes C3)",
    }
    bnoh_pairs = {f"{first} vs {second}": "not distinguishable"
                  for index, first in enumerate(five)
                  for second in five[index + 1:]}
    for substrate in _SUBSTRATES:
        table = md.discrimination_table(substrate, candidates=five)
        check(f"{substrate} status table matches A6",
              all(table["status"].loc[truth, candidate]
                  == status[substrate][truth][index]
                  for truth in five for index, candidate in enumerate(five)),
              "\n" + table["status"].to_string())
        check(f"{substrate} recovered matches A6",
              all(table["recovered"].loc[truth] == recovered[substrate][truth]
                  for truth in five),
              table["recovered"].to_string())
        want = four_ome_pairs if substrate == "4OMe-BnOH" else bnoh_pairs
        check(f"{substrate} pair table matches A6",
              table["pairs"].to_dict() == want, str(table["pairs"].to_dict()))


def test_candidate_tie_rule():
    print("\nthe candidate tie rule")
    folds = list(range(10))
    A = np.arange(1.0, 11.0)
    B = A + 1.0
    C = A + np.tile([0.6, -0.5], 5)
    scores = pd.DataFrame([A, B, C], index=("A", "B", "C"), columns=folds)
    result = md.candidate_tie(scores)
    check("A is best", result.loc["A", "status"] == "best",
          result.loc["A", "status"])
    check("B is excluded", result.loc["B", "status"] == "excluded",
          result.loc["B", "status"])
    check("C is tied", result.loc["C", "status"] == "tied",
          result.loc["C", "status"])


def test_degeneracy_flags():
    print("\nthe degeneracy flags")
    table = _summary_tables()["4OMe-BnOH"]
    buffers = sorted(set(table["rows"].buffer))
    names = md.candidate_names("C1", "4OMe-BnOH", buffers)
    lower, upper = md.candidate_bounds("C1", "4OMe-BnOH", buffers)
    x = 0.5 * (lower + upper)
    x[names.index("Ea_act")] = 199.0
    verdict, _ = md.candidate_degeneracy(
        "C1", {"success": True, "x": x}, table, "4OMe-BnOH")
    check("an off-midpoint parameter is at bound", "at bound: Ea_act" in verdict,
          verdict)
    x[names.index("lk_f")] = -12.0
    x[names.index("lk_r")] = -10.0
    verdict, _ = md.candidate_degeneracy(
        "C1", {"success": True, "x": x}, table, "4OMe-BnOH")
    check("the clocks are outside the window on every curve",
          "clocks outside the run window on 147 of 147 curves" in verdict,
          verdict)


if __name__ == "__main__":
    test_window_slopes_on_a_line()
    test_summaries_are_nan_when_not_admitted()
    test_summary_anchors()
    test_candidate_curves_match_the_ode()
    test_c5_matches_the_ode()
    test_run_likelihood_matches_the_dense_normal()
    test_midpoint_likelihood_anchors()
    test_a6_discrimination_anchors()
    test_candidate_tie_rule()
    test_degeneracy_flags()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
