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

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from fit_dataset import REFERENCE_OMITS_UNRULED, TWO_AXIS_BLOCK, build_curves
from fit_kinetics import (MODEL_PARENTS, MODEL_STAGES, STAGE_ONE, STAGE_TWO,
                          _assert_observation_design, ladder_checks,
                          ladder_f_test)
from verify_enzyme import reference_design

FAILURES = []

MANIFEST_PATH = "data/manifest.csv"
SHEET_DIRECTORY = "data/data"


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


def test_every_fitted_curve_has_the_design_its_observable_assumes():
    """
    The observation switch is only right where the design agrees with [enz].

    Over the two blocks the sequential fit supports, every catalysed curve's
    reference omits the enzyme (so it is an increment) and every enzyme-free
    curve's reference omits the H2O2 (so it is the raw background). The two
    exceptions are exps 3 and 6, whose layouts `verify_enzyme.reference_design`
    cannot classify; they stay None until R0.6 and fall back to the absolute
    observable, which is what stage 1 always used.
    """
    print("\nthe observation matches each curve's reference design")
    blocks = (("BnOH", 25.0, "Phosphate"), ("4OMe-BnOH", 40.0, "Phosphate"))
    curves, _ = build_curves()
    for block in blocks:
        scoped = [c for c in curves if c.group == block]
        catalysed = [c for c in scoped if c.conditions.e0 > 0]
        background = [c for c in scoped if c.conditions.e0 == 0]
        wrong = [(c.experiment, c.reference_omits) for c in catalysed
                 if c.reference_omits != "enzyme"]
        check(f"{block[0]}: every catalysed curve's reference omits the enzyme",
              not wrong, str(wrong))
        wrong = [(c.experiment, c.reference_omits) for c in background
                 if c.reference_omits not in ("h2o2", None)]
        check(f"{block[0]}: every enzyme-free curve is a raw background",
              not wrong, str(wrong))
        unruled = {c.experiment for c in scoped if c.reference_omits is None}
        check(f"{block[0]}: the only unruled designs are exps 3 and 6",
              unruled <= REFERENCE_OMITS_UNRULED, str(unruled))


def test_the_design_guard_rejects_a_mislabelled_curve():
    """
    `sequential_fit` must not silently fit a curve whose design contradicts
    its [enz]. The guard is tested directly so this stays fast -- it is the
    one rule that keeps the observation switch honest.
    """
    print("\nthe observation guard")
    curves, _ = build_curves()
    catalysed = next(c for c in curves
                     if c.conditions.e0 > 0 and c.reference_omits == "enzyme")
    background = next(c for c in curves
                      if c.conditions.e0 == 0 and c.reference_omits == "h2o2")

    saved = catalysed.reference_omits
    catalysed.reference_omits = None
    try:
        _assert_observation_design([catalysed])
        raised = False
    except ValueError:
        raised = True
    finally:
        catalysed.reference_omits = saved
    check("a catalysed curve with an unknown design is rejected", raised)

    saved = background.reference_omits
    background.reference_omits = "enzyme"
    try:
        _assert_observation_design([background])
        raised = False
    except ValueError:
        raised = True
    finally:
        background.reference_omits = saved
    check("an enzyme-free curve labelled catalysed is rejected", raised)

    check("the matching designs pass",
          _assert_observation_design([catalysed, background]) is None)


def _sheet(experiment, manifest):
    return pd.read_excel(os.path.join(SHEET_DIRECTORY,
                                      manifest.loc[experiment, "xls_file"]),
                         sheet_name="Sheet1", header=None)


def test_the_two_axis_sheets_read_as_enzyme_references():
    """
    The Sum rows under the two-axis cuvette tables are not cuvettes.

    Exps 135-151 state each sum with a real total volume, so the old reader
    took `Sum:` and `Sum*9:` for two more cuvettes, split the table evenly and
    returned "other". Their `Ref.` rows carry `Enz` 0.000 like every catalysed
    sheet, so stopping at the sum row leaves seven measured against seven
    references and the reference omits only the enzyme.
    """
    print("\nthe two-axis sheets are enzyme references")
    manifest = pd.read_csv(MANIFEST_PATH).set_index("experiment")
    for number in (135, 140, 151):
        design = reference_design(_sheet(number, manifest))
        check(f"exp {number} classifies as enzyme", design == "enzyme",
              str(design))


def test_the_sum_rows_change_only_the_two_axis_designs():
    """A2: no experiment outside the two-axis block changes design."""
    print("\nthe Sum rows change only the two-axis designs")
    manifest = pd.read_csv(MANIFEST_PATH).set_index("experiment")
    changed = []
    for number in sorted(manifest.index):
        if number in TWO_AXIS_BLOCK:
            continue
        sheet = _sheet(number, manifest)
        before = reference_design(sheet, stop_at_sum=False)
        after = reference_design(sheet)
        if before != after:
            changed.append((int(number), before, after))
    check("no experiment outside the two-axis block changes design",
          not changed, str(changed))


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
    test_every_fitted_curve_has_the_design_its_observable_assumes()
    test_the_design_guard_rejects_a_mislabelled_curve()
    test_the_two_axis_sheets_read_as_enzyme_references()
    test_the_sum_rows_change_only_the_two_axis_designs()
    test_the_saved_ladders()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
