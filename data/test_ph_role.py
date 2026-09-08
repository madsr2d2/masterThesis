"""
Tests for ph_role.py.

The load-bearing one is `test_a_planted_pair_of_orders_is_recovered`: the
folder's headline is that three of the four ladders agree on an order in
[HOO-] and the fourth does not, and that split is only worth reporting if the
fitting machinery can read back a planted order at all, on data shaped like a
pH ladder rather than like the two-axis block `scope.orders` was built for.

    python data/test_ph_role.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import induction
import ph_role
import scope
from ph_role import (boric_turnover, pooled_rate_order, rate_ladder,
                     rate_ladder_table, rate_ladders)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def _planted_ladder(order_s0, order_hoo, runs=(1, 2, 3, 4, 5),
                    pH_values=(6.0, 6.75, 7.5, 8.25, 9.0),
                    s0_values=(1.85, 3.70, 5.55, 7.40), noise=0.0, seed=7):
    """
    A synthetic pH ladder: `runs` experiments, each a four-rung substrate
    ladder at its own pH, [HOO-] standing in for pH the way it does for the
    two-axis block (`PAX_AXIS = "hoo"`) so the same fit that reads that
    block's cuvette-matched ladder can be pointed at this one.
    """
    generator = np.random.default_rng(seed)
    rows = []
    for run, pH in zip(runs, pH_values):
        # A toy [HOO-]: rises 10x per pH unit, as it must while [H2O2] and
        # the peroxide's own pKa are held fixed -- the qualitative fact the
        # fit is allowed to lean on, not a claim about the acid's real pKa.
        hoo = 1.0e-4 * 10.0 ** (pH - 7.0)
        for sample, s0 in enumerate(s0_values, start=1):
            rate = 1.0e-5 * (s0 ** order_s0) * (hoo / 1.0e-4) ** order_hoo
            if noise:
                rate *= float(np.exp(generator.normal(0.0, noise)))
            rows.append({"experiment": run, "sample": sample, "pH": pH,
                        "s0": s0, "hoo": hoo, "h2o2": 82.5,
                        "vmax": rate, "live": True})
    return pd.DataFrame(rows)


def test_a_planted_pair_of_orders_is_recovered():
    """Both orders, exactly, with no noise; within error once noise is added."""
    print("\na planted (order_s0, order_hoo) pair, on data shaped like a ladder")
    for order_s0, order_hoo in ((0.5, 0.6), (0.0, -0.03), (-0.2, 0.9)):
        exact = rate_ladder((901, 902, 903, 904, 905), parameter="vmax",
                            frame=_planted_ladder(order_s0, order_hoo))
        check(f"order_s0={order_s0:+.2f}: recovered exactly with no noise",
              abs(exact["s0"]["slope"] - order_s0) < 1e-6,
              f"{exact['s0']['slope']:.6f}")
        check(f"order_hoo={order_hoo:+.2f}: recovered exactly with no noise",
              abs(exact["hoo"]["slope"] - order_hoo) < 1e-6,
              f"{exact['hoo']['slope']:.6f}")
        noisy = rate_ladder((901, 902, 903, 904, 905), parameter="vmax",
                            frame=_planted_ladder(order_s0, order_hoo,
                                                  noise=0.05))
        check(f"order_s0={order_s0:+.2f}: within 3 stderr with 5% noise",
              abs(noisy["s0"]["slope"] - order_s0)
              < 3 * noisy["s0"]["stderr"],
              f"{noisy['s0']['slope']:+.3f} +/- {noisy['s0']['stderr']:.3f}")
        check(f"order_hoo={order_hoo:+.2f}: within 3 stderr with 5% noise",
              abs(noisy["hoo"]["slope"] - order_hoo)
              < 3 * noisy["hoo"]["stderr"],
              f"{noisy['hoo']['slope']:+.3f} +/- {noisy['hoo']['stderr']:.3f}")
        check("a planted frame gets no schedule control",
              np.isnan(exact["schedule_collinearity"]))


def test_a_flat_ladder_reads_as_flat_and_not_as_noise():
    """The boric result is a real zero, not an artefact of too few points."""
    print("\na planted zero order, at the boric ladder's own size and pH span")
    flat = rate_ladder((910, 911, 912, 913, 914, 915, 916, 917, 918),
                       parameter="vmax",
                       frame=_planted_ladder(
                           0.6, 0.0,
                           runs=range(910, 919),
                           pH_values=(8.46, 8.98, 9.23, 9.40, 9.50,
                                     9.51, 9.70, 10.07, 10.34),
                           noise=0.05, seed=11))
    check("a planted zero is not mistaken for a strong order",
          abs(flat["hoo"]["slope"]) < 3 * flat["hoo"]["stderr"],
          f"{flat['hoo']['slope']:+.3f} +/- {flat['hoo']['stderr']:.3f}")


def test_the_two_axis_ladders_are_filtered_to_strong_runs():
    """`_ladder_scope` drops the weak runs on the two-axis ladders, and only there."""
    print("\nwhich runs each ladder is actually read over")
    strong = set(scope.strong_runs())
    for label, ladder in (("low", scope.PH_LADDER_TWO_AXIS_LOW),
                          ("high", scope.PH_LADDER_TWO_AXIS_HIGH)):
        kept = ph_role._ladder_scope(ladder)
        check(f"two-axis {label}: every kept run is a strong run",
              set(kept) <= strong, f"{kept}")
        check(f"two-axis {label}: at least one run was actually dropped",
              len(kept) < len(ladder), f"{len(kept)} of {len(ladder)}")
    for label, ladder in (("phosphate", scope.PH_LADDER_PHOSPHATE),
                          ("boric", scope.PH_LADDER_BORIC)):
        kept = ph_role._ladder_scope(ladder)
        check(f"{label}: not filtered -- concentration_agreement has "
              f"nothing to say off the two-axis block's own L design",
              kept == tuple(ladder), f"{kept}")


def test_the_boric_ladder_turns_over_rather_than_saturating():
    """A single log-log slope over a rise-then-fall ladder hides the rise."""
    print("\nthe boric ladder split at its own peak")
    turnover = boric_turnover()
    check("the peak sits at exp 43, pH 9.50",
          turnover["peak_experiment"] == 43
          and abs(turnover["peak_pH"] - 9.50) < 0.01,
          f"exp {turnover['peak_experiment']}, pH {turnover['peak_pH']}")
    check("below the peak the order is positive, not flat",
          turnover["below"]["hoo"]["slope"]
          > 3 * turnover["below"]["hoo"]["stderr"] > 0.0,
          f"{turnover['below']['hoo']['slope']:+.4f} +/- "
          f"{turnover['below']['hoo']['stderr']:.4f}")
    check("...and weaker than phosphate or pyrophosphate's",
          turnover["below"]["hoo"]["slope"] < 0.4)
    above = turnover["above_by_experiment"]
    check("both runs past pH 10 sit below the ladder's own peak",
          all(row["median_vmax"] < 1.42e-4 for row in above.values()),
          f"{above}")
    check("vmax_corrected does not rescue the decline -- not the O2 side "
          "reaction",
          all(abs(row["median_vmax_corrected"] - row["median_vmax"])
              < 0.25 * row["median_vmax"] for row in above.values()),
          f"{above}")


def test_regressions():
    """The published numbers, locked down. `python data/ph_role.py` prints them."""
    print("\nthe archive's own four ladders")
    table = rate_ladder_table()
    check("four ladders", len(table) == 4, f"{len(table)}")
    expected = {
        "phosphate 4OMe": (0.5959, 0.0321),
        "boric 4OMe": (-0.0335, 0.0400),
        "pyrophosphate BnOH 136-142": (0.6158, 0.0655),
        "pyrophosphate BnOH 143-151": (0.5704, 0.0565),
    }
    for ladder, (order, stderr) in expected.items():
        row = table.loc[ladder]
        check(f"{ladder}: order_hoo {order:+.4f} +/- {stderr:.4f}",
              abs(row.order_hoo - order) < 0.001
              and abs(row.stderr_hoo - stderr) < 0.001,
              f"{row.order_hoo:+.4f} +/- {row.stderr_hoo:.4f}")
    check("boric alone sits at zero; the other three do not",
          abs(table.loc["boric 4OMe"].order_hoo)
          < 3 * table.loc["boric 4OMe"].stderr_hoo
          and all(abs(table.loc[l].order_hoo) > 5 * table.loc[l].stderr_hoo
                  for l in expected if l != "boric 4OMe"))

    print("\npooling: boric is the reason four ladders do not agree")
    all_four = pooled_rate_order()
    three = pooled_rate_order(drop=("boric 4OMe",))
    check("all four: chi2 is large (boric does not belong)",
          all_four["chi2"] > 50.0, f"{all_four['chi2']:.1f} on {all_four['dof']}")
    check("without boric: chi2 is small (the other three agree)",
          three["chi2"] < 2.0, f"{three['chi2']:.2f} on {three['dof']}")
    check("without boric: pooled order +0.594 +/- 0.026",
          abs(three["pooled"] - 0.594) < 0.002
          and abs(three["stderr"] - 0.026) < 0.002,
          f"{three['pooled']:+.4f} +/- {three['stderr']:.4f}")

    print("\nthe two-axis block's own cuvette-matched reading, for comparison")
    published = scope.ph_order(parameter="vmax", scope=scope.strong_runs())
    check("the published pooled two-axis pH order is +0.554 +/- 0.040",
          abs(published.loc["pooled", "order"] - 0.5536) < 0.001
          and abs(published.loc["pooled", "stderr"] - 0.0397) < 0.001,
          f"{published.loc['pooled', 'order']:+.4f} +/- "
          f"{published.loc['pooled', 'stderr']:.4f}")
    check("this module's cruder pooled fit lands within its own error of it",
          abs(three["pooled"] - published.loc["pooled", "order"])
          < 2 * (three["stderr"] + published.loc["pooled", "stderr"]))


if __name__ == "__main__":
    test_a_planted_pair_of_orders_is_recovered()
    test_a_flat_ladder_reads_as_flat_and_not_as_noise()
    test_the_two_axis_ladders_are_filtered_to_strong_runs()
    test_the_boric_ladder_turns_over_rather_than_saturating()
    test_regressions()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
