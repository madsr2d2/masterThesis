"""
Tests for the M0-M4 ladder's bookkeeping in fit_kinetics.py: which model a
term is read against, and which comparisons are not comparisons at all.

The ladder was first read against whichever rung was LISTED before each model
(DATA_VERIFICATION.md 2026-09-14), so M2's K4 was differenced against M1b, a
model it does not contain, and two stage-2 costs fitted on different frozen
backgrounds were differenced as if nested. And a larger model that fitted
WORSE than the one it contains -- which no optimum can do -- was read as a term
that "does not earn". These pin the reading on planted saves, then on the real
ones.

Fast: it reads saved fits and never runs the optimiser.

    python data/test_fit_ladder.py
"""
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from fit_dataset import build_curves
from fit_kinetics import (MODEL_PARENTS, MODEL_STAGES, STAGE_ONE, STAGE_TWO,
                          ladder_checks, ladder_f_test)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


PLANTED_BLOCK = ("BnOH", 25.0, "Phosphate")
PLANTED_POINTS = {1: 1000, 2: 1000}
BACKGROUND = {"k_can": 1.0, "k3": 1.0, "k0": 1e-9, "k5": 0.0, "k6": 0.0, "r": 0.1}


def _stage(free, cost, constants=BACKGROUND, converged=True):
    return {"free_parameters": list(free), "cost": float(cost),
            "rms_absorbance": 0.01, "converged": converged,
            "constants": dict(constants)}


def _planted_model_ladder(directory, m4_converged=True):
    """
    A ladder whose costs are chosen so that the listed-order reading and the
    parent reading disagree on every verdict the review found.
    """
    sink = BACKGROUND | {"k_sink": 1e-3}
    saturated = BACKGROUND | {"km_s_background": 0.1}
    stages = {
        "M0": (_stage(STAGE_ONE, 1000), _stage(STAGE_TWO, 500)),
        "M1": (_stage(STAGE_ONE, 1000), _stage(STAGE_TWO + ("km_s",), 400)),
        "M1b": (_stage(STAGE_ONE + ("km_s_background",), 600, saturated),
                _stage(STAGE_TWO + ("km_s",), 300, saturated)),
        # fits WORSE than M1, which it contains
        "M2": (_stage(STAGE_ONE, 1000),
               _stage(STAGE_TWO + ("km_s", "K4"), 450)),
        # a different background, so its stage 2 is not M2's
        "M3": (_stage(STAGE_ONE + ("k_sink",), 990, sink),
               _stage(STAGE_TWO + ("km_s", "K4"), 200, sink)),
        "M4": (_stage(STAGE_ONE + ("k_sink",), 990, sink),
               _stage(STAGE_TWO + ("km_s", "K4", "k_act_r", "K_act"), 100,
                      sink, converged=m4_converged)),
        "M4b": (_stage(STAGE_ONE + ("km_s_background",), 600, saturated),
                _stage(STAGE_TWO + ("km_s", "k_act_r", "K_act"), 299,
                       saturated)),
    }
    for model, (one, two) in stages.items():
        path = os.path.join(directory, f"BnOH_25C_Phosphate_{model}.json")
        with open(path, "w") as handle:
            json.dump({"stage_1": one, "stage_2": two}, handle)


def _row(table, model, stage):
    return table[(table.model == model) & (table.stage == stage)].iloc[0]


def _planted_table(**options):
    with tempfile.TemporaryDirectory() as directory:
        _planted_model_ladder(directory, **options)
        return ladder_f_test(PLANTED_BLOCK, directory=directory,
                             points=PLANTED_POINTS)


def test_each_term_is_read_against_the_model_it_extends():
    print("\nthe parent, not the previous rung")
    check("every rung but M0 names its parent",
          set(MODEL_PARENTS) == set(MODEL_STAGES) - {"M0"},
          f"{sorted(MODEL_PARENTS)}")
    table = _planted_table()
    for model, stage, parent in (("M2", 2, "M1"), ("M3", 1, "M2"),
                                 ("M4b", 2, "M1b"), ("M1b", 1, "M1")):
        row = _row(table, model, stage)
        check(f"{model} stage {stage} is read against {parent}",
              row.vs == parent, f"read against {row.vs}")


def test_a_larger_model_that_fits_worse_is_an_optimiser_failure():
    print("\na larger model cannot fit worse at its optimum")
    row = _row(_planted_table(), "M2", 2)
    check("M2 above M1's cost is an optimiser failure, not a term that fails",
          row.verdict == "optimiser failure", f"{row.verdict}, F {row.f:.1f}")


def test_stage_two_is_not_compared_across_backgrounds():
    print("\nstage 2 on different frozen backgrounds")
    table = _planted_table()
    for model in ("M3", "M1b"):
        row = _row(table, model, 2)
        check(f"{model} stage 2 is not differenced against its parent's",
              row.verdict == "not nested: stage 1 differs" and row.f != row.f,
              f"{row.verdict}, F {row.f}")


def test_the_bar_is_applied_only_to_comparable_rows():
    print("\nthe F bar")
    table = _planted_table()
    expected = {("M1b", 1): "earns", ("M3", 1): "does not earn",
                ("M4", 2): "earns", ("M4b", 2): "does not earn"}
    for (model, stage), verdict in expected.items():
        row = _row(table, model, stage)
        check(f"{model} stage {stage} {verdict}", row.verdict == verdict,
              f"{row.verdict}, F {row.f:.2f}")
    unconverged = _row(_planted_table(m4_converged=False), "M4", 2)
    check("an unconverged fit gets no verdict on the bar",
          unconverged.verdict == "not converged", unconverged.verdict)


def test_the_saved_ladders():
    """The real saves, read the way the 2026-09-14 review read them."""
    print("\nthe saved M0-M4b ladders")
    curves, _ = build_curves()
    bnoh = ("BnOH", 25.0, "Phosphate")
    anisyl = ("4OMe-BnOH", 40.0, "Phosphate")
    tables = {}
    for block in (bnoh, anisyl):
        scoped = [c for c in curves if c.group == block]
        points = {1: sum(len(c) for c in scoped if c.conditions.e0 == 0),
                  2: sum(len(c) for c in scoped if c.conditions.e0 > 0)}
        tables[block] = ladder_f_test(block, points=points)
    expected = (
        (anisyl, "M3", 1, "earns"),                     # the sink, stage 1
        (anisyl, "M4b", 2, "earns"),                    # activation on M1b
        (anisyl, "M1", 2, "optimiser failure"),
        (anisyl, "M2", 2, "optimiser failure"),
        (bnoh, "M3", 1, "earns"),                       # the sink, stage 1
        (bnoh, "M3", 2, "not nested: stage 1 differs"),
        (bnoh, "M4", 2, "not converged"),
        (bnoh, "M4b", 2, "does not earn"),              # activation on M1b
    )
    for block, model, stage, verdict in expected:
        row = _row(tables[block], model, stage)
        check(f"{block[0]} {model} stage {stage}: {verdict}",
              row.verdict == verdict, f"{row.verdict}, F {row.f:.1f}")

    checks = {block: ladder_checks(block, curves=curves)
              for block in (bnoh, anisyl)}
    background = _row(checks[bnoh], "M1b", 1)
    check("BnOH M1b's saturated background leaves the model's substrate order"
          " below the data's",
          background.order_model < background.order_data
          - background.stderr_data,
          f"model {background.order_model:+.3f}, data "
          f"{background.order_data:+.3f} +/- {background.stderr_data:.3f}")
    confound = _row(checks[anisyl], "M1b", 1)
    check("the 4OMe enzyme-free substrate ladders move [buf] with [S]",
          abs(confound.s0_buf_median_r) > 0.9,
          f"median r {confound.s0_buf_median_r:+.3f} over "
          f"{confound.s0_buf_runs} runs")
    lagging = _row(checks[anisyl], "M4b", 2)
    check("4OMe M4b's activation lags more catalysed curves than the data do",
          lagging.lag_model > lagging.lag_data,
          f"model {lagging.lag_model}, data {lagging.lag_data}")


if __name__ == "__main__":
    test_each_term_is_read_against_the_model_it_extends()
    test_a_larger_model_that_fits_worse_is_an_optimiser_failure()
    test_stage_two_is_not_compared_across_backgrounds()
    test_the_bar_is_applied_only_to_comparable_rows()
    test_the_saved_ladders()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
