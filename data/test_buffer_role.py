"""
Tests for buffer_role.py.

The load-bearing one is `test_the_species_test_can_read_a_planted_answer`:
this module's headline is a NEGATIVE -- that the archive cannot tell general
acid from general base -- and a negative is worthless unless the machinery
could have read a positive. So the test plants both and checks that each is
recovered and that the other is then excluded.

    python data/test_buffer_role.py
"""
import sys
import os

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import buffer_role
import scope
import solution_chemistry
from buffer_role import (catalytic_coefficient, identity_overlap,
                         overlap_width, separable, species_prediction,
                         titration_table)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def test_the_species_prediction_is_henderson_hasselbalch():
    """The two predicted ratios are the base and acid fractions, and nothing else."""
    print("\nwhat each species predicts")
    got = species_prediction(7.00, 7.53)
    for pH, fraction in ((7.00, got["low_base_fraction"]),
                         (7.53, got["high_base_fraction"])):
        expected = 1.0 / (1.0 + 10.0 ** (got["pka"] - pH))
        check(f"the base fraction at pH {pH:.2f} is Henderson-Hasselbalch",
              abs(fraction - expected) < 1e-9,
              f"{fraction:.4f} against {expected:.4f}")
    # To four decimals, not to machine precision: `general_base` is the ratio
    # of the base CONCENTRATIONS and the fractions divide by (acid + base),
    # which is not the same total at two pHs when a third species exists.
    check("general base is the ratio of the base fractions",
          abs(got["general_base"]
              - got["high_base_fraction"] / got["low_base_fraction"]) < 1e-4,
          f"{got['general_base']:.6f} against "
          f"{got['high_base_fraction'] / got['low_base_fraction']:.6f}")
    # Not "1 - the base fraction": phosphate has more than two species, so
    # acid + base is not the total and the complement is not the acid. What
    # must hold exactly is the ratio of the two predictions, which is
    # 10^(dpH) by Henderson-Hasselbalch whatever else is in solution.
    check("base over acid is 10^(delta pH), as it has to be",
          abs(got["general_base"] / got["general_acid"]
              - 10.0 ** (7.53 - 7.00)) < 1e-3,
          f"{got['general_base'] / got['general_acid']:.4f} against "
          f"{10.0 ** 0.53:.4f}")
    check("and they straddle 1, which is what makes the test a test",
          got["general_acid"] < 1.0 < got["general_base"],
          f"{got['general_acid']:.2f} and {got['general_base']:.2f}")
    check("a titration at ONE pH cannot separate them: the ratio is 1",
          abs(species_prediction(7.0, 7.0)["general_base"] - 1.0) < 1e-9
          and abs(species_prediction(7.0, 7.0)["general_acid"] - 1.0) < 1e-9)


def _planted_ratio(ratio, levels=(50.0, 100.0, 150.0, 200.0), noise=0.0, seed=4):
    """
    Two runs at pH 7.00 and two at pH 7.53, built from a known coefficient
    ratio and the substrate concentrations the real runs sit at.
    """
    generator = np.random.default_rng(seed)
    rows = []
    design = ((32, 7.00, 8.251, 1.0), (34, 7.00, 8.251, 1.0),
              (36, 7.53, 57.900, ratio), (37, 7.53, 12.228, ratio))
    for experiment, pH, s0, factor in design:
        base = 3.0e-5 if pH < 7.2 else 2.6e-4
        for index, buffer_mM in enumerate(levels):
            scale = s0 ** buffer_role.SUBSTRATE_ORDER
            rate = (base + 1.9e-7 * factor * buffer_mM) * scale
            if noise:
                rate *= float(np.exp(generator.normal(0.0, noise)))
            rows.append({"experiment": experiment, "sample": index + 1,
                         "pH": pH, "s0": s0, "buf": buffer_mM,
                         "v_peak": rate, "live": True,
                         "buffer": "Phosphate", "h2o2": 82.5})
    return pd.DataFrame(rows)


def test_the_species_test_can_read_a_planted_answer():
    """
    The negative is only worth something if a positive would have been read.
    """
    print("\nthe species test, on planted data")
    prediction = species_prediction(7.00, 7.53)
    for name, planted in (("general_base", prediction["general_base"]),
                          ("general_acid", prediction["general_acid"]),
                          ("spectator", 1.0)):
        exact = catalytic_coefficient(frame=_planted_ratio(planted))
        check(f"a planted {name.replace('_', ' ')} is recovered exactly "
              f"with no noise",
              abs(exact["ratio"] - planted) < 1e-6,
              f"{exact['ratio']:.6f} against {planted:.6f}")
        got = catalytic_coefficient(frame=_planted_ratio(planted, noise=0.02))
        verdict = separable(got, prediction)
        # Within the fit's OWN error, not within an absolute tolerance: the
        # high-pH runs sit ten times higher, so an unweighted fit gives their
        # residuals ten times the leverage and the ratio is loose by design.
        check(f"...and within its own error with noise",
              abs(got["ratio"] - planted) < 3 * got["ratio_stderr"],
              f"{got['ratio']:.3f} +- {got['ratio_stderr']:.3f} against "
              f"{planted:.3f}")
        check(f"...and {name.replace('_', ' ')} survives its own data",
              name in verdict["survivors"], f"{verdict['survivors']}")
        others = [other for other in ("general_base", "general_acid")
                  if other != name]
        check(f"...while {others[0].replace('_', ' ')} is excluded",
              others[0] not in verdict["survivors"]
              or name == "spectator",
              f"{verdict['survivors']}")


def test_the_archive_excludes_nothing():
    """The published negative, and the reason for it."""
    print("\nthe archive's own answer")
    prediction = species_prediction(7.00, 7.53)
    got = catalytic_coefficient(drop=(35,))
    verdict = separable(got, prediction)
    check("all three hypotheses survive",
          len(verdict["survivors"]) == 3, f"{verdict['survivors']}")
    check("because the measurement sits between the two predictions",
          prediction["general_acid"] < got["ratio"]
          < prediction["general_base"],
          f"{prediction['general_acid']:.2f} < {got['ratio']:.2f} < "
          f"{prediction['general_base']:.2f}")
    check("and two of its own errors reach past both",
          2 * got["ratio_stderr"] > max(
              abs(got["ratio"] - prediction["general_base"]),
              abs(got["ratio"] - prediction["general_acid"])),
          f"2 x {got['ratio_stderr']:.2f}")
    check("the high-pH coefficient is the poorly determined one",
          got["high_stderr"] / abs(got["high"])
          > 0.3, f"{got['high_stderr'] / abs(got['high']):.2f}")


def test_the_titration_table_and_the_maps():
    """The five runs, and the two confound maps."""
    print("\nthe titrations and the maps")
    table = titration_table()
    check("five titrations, twenty curves",
          len(table) == 5 and int(table.curves.sum()) == 20,
          f"{len(table)} runs, {int(table.curves.sum())} curves")
    check("every one holds [S] fixed inside itself",
          all(scope.frame((int(e),)).s0.nunique() == 1
              for e in table.experiment))
    check("they sit at three pH values", table.pH.nunique() == 3,
          f"{sorted(table.pH.unique())}")
    low = table[table.pH < buffer_role.TITRATION_PH_SPLIT]
    check("and the order saturates within the low-pH pair",
          float(low[low.buffer_low < 10].order.iloc[0])
          > float(low[low.buffer_low >= 10].order.iloc[0]),
          f"{low.order.round(3).tolist()}")

    identity = identity_overlap()
    widest = overlap_width(identity)
    check("the widest two-buffer pH overlap is about two units",
          abs(widest["width"] - 2.01) < 0.02, f"{widest['width']:.2f}")
    check("and those two cells share no peroxide at all",
          not widest["shares_peroxide"], f"{widest['peroxide']}")
    check("the peroxide check is on VALUES, not on ranges",
          min(widest["peroxide"][1]) < min(widest["peroxide"][0])
          < max(widest["peroxide"][1]),
          "a range test would have called this a match")


def test_regressions():
    """The numbers buffer/ANALYSIS.md quotes."""
    print("\nthe published numbers")
    table = titration_table().set_index("experiment")
    for experiment, order, error in ((32, 0.371, 0.018), (34, 0.792, 0.185),
                                     (37, 0.122, 0.039)):
        row = table.loc[experiment]
        check(f"exp {experiment}: order {order:+.3f} +- {error:.3f}",
              abs(row.order - order) < 0.001 and abs(row.order_stderr - error)
              < 0.001,
              f"{row.order:+.3f} +- {row.order_stderr:.3f}")
    # NOT "no order is resolved": exp 37's +0.122 +- 0.039 is three standard
    # errors from zero. What collapses is the SIZE.
    check("above the pKa the largest order is a third of the low-pH one",
          max(table.loc[e].order for e in (35, 36, 37))
          < 0.4 * table.loc[32].order,
          f"{max(table.loc[e].order for e in (35, 36, 37)):.3f} against "
          f"{table.loc[32].order:.3f}")
    check("the buffer-independent term is what grew, roughly eightfold",
          6.0 < table.loc[[36, 37]].intercept.mean() / table.loc[32].intercept
          < 14.0,
          f"{table.loc[[36, 37]].intercept.mean() / table.loc[32].intercept:.1f}x")
    check("while the coefficient did not",
          0.5 < table.loc[[36, 37]].coefficient.mean()
          / table.loc[32].coefficient < 3.0,
          f"{table.loc[[36, 37]].coefficient.mean() / table.loc[32].coefficient:.2f}x")
    got = catalytic_coefficient(drop=(35,))
    check("the measured ratio is +1.06 +- 0.77",
          abs(got["ratio"] - 1.06) < 0.02
          and abs(got["ratio_stderr"] - 0.77) < 0.02,
          f"{got['ratio']:+.3f} +- {got['ratio_stderr']:.3f}")
    check("exp 35 is dropped on its own internal consistency",
          float(table.loc[35].order_r2) < 0.2,
          f"R2 {table.loc[35].order_r2:.3f}")


if __name__ == "__main__":
    test_the_species_prediction_is_henderson_hasselbalch()
    test_the_species_test_can_read_a_planted_answer()
    test_the_archive_excludes_nothing()
    test_the_titration_table_and_the_maps()
    test_regressions()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
