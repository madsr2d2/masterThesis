"""
Tests for scope.py's selection and order machinery.

Separate from `test_fit_kinetics`, which owns the same module's block
definitions, because that suite runs the optimiser and takes minutes. What is
here is fast enough to run every time, and one of the two things it guards
against had already happened twice.

    python data/test_scope.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import induction
import scope

FAILURES = []


def check(name, condition, detail=""):
    print(f"  {'pass' if condition else 'FAIL'}  {name}"
          f"{': ' + detail if detail else ''}")
    if not condition:
        FAILURES.append(f"{name} {detail}")


def test_an_axis_the_offsets_absorb_is_not_reported():
    """
    The identification guard, planted at the level where it failed.

    An axis that varies over the data but is CONSTANT INSIDE EVERY RUN is
    collinear with the per-experiment indicators `within=True` adds, so `lstsq`
    splits its coefficient through the pseudo-inverse and returns a number with
    a small standard error attached. Testing the whole column instead of the
    within-run spread let that through: `scope.arm_orders` on the block's
    peroxide arm returned a substrate order of -6.094 +/- 0.078 on 63 curves,
    on an arm where every run holds the substrate fixed.
    """
    print("\nan axis the experiment offsets absorb")
    rows = []
    for experiment, fixed in enumerate((1.0, 4.0, 16.0), start=1):
        for step in (1.0, 2.0, 4.0, 8.0):
            rows.append({"experiment": experiment, "s0": fixed, "h2o2": step,
                         "live": True, "v0": 1e-5 * step ** 0.5,
                         "vmax": 1e-5 * step ** 0.5})
    planted = pd.DataFrame(rows)

    within = scope.orders("vmax", frame=planted, within=True)
    check("the axis that moves inside runs is measured",
          abs(within["order_h2o2"] - 0.5) < 1e-6,
          f"{within['order_h2o2']:+.4f}")
    check("and the axis that only moves between them is not",
          np.isnan(within["order_s0"]), f"{within['order_s0']}")

    # Without offsets the same axis IS identified, and the guard must not
    # refuse it -- `orders(within=False)` is how the pooled row of
    # `order_table` is built.
    pooled = scope.orders("vmax", frame=planted, within=False)
    check("but with no offsets to hide behind it is identified again",
          np.isfinite(pooled["order_s0"]), f"{pooled['order_s0']}")


def test_the_block_still_measures_both_axes():
    """The guard must not have made the block's own orders disappear."""
    print("\nthe block itself is unaffected")
    table = scope.order_table()
    finite = table[["order_s0", "order_h2o2"]].notna().all().all()
    check("every order in the block's table is a number", finite,
          f"{table[['order_s0', 'order_h2o2']].isna().sum().to_dict()}")
    row = table.loc[("vmax", "within-experiment")]
    check("and v_max's substrate order is the published one",
          abs(row.order_s0 - 0.0906) < 5e-4, f"{row.order_s0:+.4f}")


def test_each_arm_agrees_with_the_joint_fit():
    """
    The L has no interior, so the joint fit assumes additivity. Measure it.

    Not a tolerance plucked from the air: each arm holds the other axis fixed
    within a run and so needs no such assumption, and if the two ever parted
    company the joint numbers would be the ones to stop quoting.
    """
    print("\nthe arms against the joint fit")
    table = scope.arm_orders()
    check("each arm reports only its own axis",
          len(table) == 8, f"{len(table)} rows")
    check("and no arm disagrees with the joint fit past 2 sigma",
          table.sigma.max() < 2.0,
          f"worst {table.sigma.max():.2f} at {table.sigma.idxmax()}")


def test_a_ph_ladder_is_derived_and_not_listed():
    """
    Every run of a ladder shares a composition AND an enzyme loading.

    Both halves matter. Dropping the composition would compare cuvettes that
    are not the same cuvette; dropping the loading would pool exps 141-142 at
    0.014 mM with exps 136-140 at 0.034 on the same composition, and since
    those two sit at the TOP of that group's pH range the pooled slope would be
    biased by whatever the rate does with enzyme.
    """
    print("\nthe pH ladders")
    ladders = scope.ph_ladders()
    check("two ladders over the whole block", len(ladders) == 2,
          f"{sorted(ladders)}")
    for label, group in ladders.items():
        check(f"{label}: one enzyme loading", group.e0.nunique() == 1,
              f"{sorted(group.e0.unique())}")
        designs = {tuple(sorted(zip(block.s0.round(6), block.h2o2.round(6))))
                   for _, block in group.groupby("experiment")}
        check(f"{label}: one composition", len(designs) == 1,
              f"{len(designs)} designs")
        check(f"{label}: at least three pH values",
              group.pH.nunique() >= scope.PH_LADDER_MINIMUM,
              f"{group.pH.nunique()}")

    check("exps 141 and 142 are not in any ladder",
          not any(set(group.experiment) & {141, 142}
                  for group in ladders.values()),
          "they share exps 136-140's composition at a different loading")


def test_the_weak_runs_flatten_the_ph_ladder():
    """
    The agreement filter is not cosmetic here, and the direction is the point.

    Exps 149-151 sit at the BOTTOM of the pH ladder and measure the cell's own
    wander rather than a rate, so keeping them stops the ladder falling where
    the chemistry does. Filtered, the two ladders agree; unfiltered they do not.
    """
    print("\nthe agreement filter, on the pH order")
    strong = scope.ph_order("vmax", scope=scope.strong_runs())
    every = scope.ph_order("vmax")
    check("the filtered slope is the steeper one",
          strong.loc["pooled", "order"] > every.loc["pooled", "order"],
          f"{strong.loc['pooled', 'order']:+.3f} against "
          f"{every.loc['pooled', 'order']:+.3f}")

    def gap(table):
        rows = table.drop(index="pooled")
        spread = float(np.hypot(*rows.stderr.to_numpy()))
        return abs(float(np.diff(rows.order.to_numpy())[0])) / spread

    check("and the two ladders agree once it is applied",
          gap(strong) < 2.5, f"{gap(strong):.1f} sigma")
    check("where unfiltered they do not", gap(every) > 4.0,
          f"{gap(every):.1f} sigma")


def test_the_schedule_control_is_a_real_control():
    """
    It only works because the two ladders were run in opposite directions.

    If they had both climbed in pH, a stock ageing over the twelve days would
    be indistinguishable from the pH order and this function would be
    decoration. Assert the premise, not just the conclusion.
    """
    print("\nthe schedule control")
    table, verdict = scope.ph_schedule_control("vmax")
    check("the two ladders are opposed against the schedule",
          verdict["opposed_schedules"],
          f"{table.pH_vs_schedule.to_dict()}")
    check("each correlation is a strong one, not a shrug",
          bool((table.pH_vs_schedule.abs() > 0.5).all()),
          f"{table.pH_vs_schedule.abs().min():.2f}")
    check("and both pH orders come out with the same sign",
          verdict["orders_agree_in_sign"], f"{table.order.to_dict()}")

    dates = scope.run_dates()
    check("every run in the block has a collection date",
          dates.date.notna().all(), f"{dates.date.isna().sum()} missing")
    check("and experiment number runs with the calendar",
          list(dates.sort_index().order) == sorted(dates.order),
          f"{list(dates.order)}")


def test_the_two_levers_are_measured_on_different_contrasts():
    """
    `hoo_consistency` is only a test if the two sides are independent.

    They are, structurally: one is measured within runs against per-experiment
    offsets, the other between runs against per-cuvette offsets. What can still
    go wrong is that one side quietly stops being measured, and then the gap is
    a comparison with a constant.
    """
    print("\nthe two levers on [HOO-]")
    result = scope.hoo_consistency("vmax")
    check("both sides carry curves",
          result["within_curves"] > 20 and result["across_curves"] > 20,
          f"{result['within_curves']}, {result['across_curves']}")
    check("both sides carry an error bar",
          result["within_stderr"] > 0 and result["across_stderr"] > 0,
          f"{result['within_stderr']:.3f}, {result['across_stderr']:.3f}")
    check("and they part company past 3 sigma",
          result["sigma"] > 3.0, f"{result['sigma']:.2f}")

    arm = induction.ladder_arms(scope.frame(scope.strong_runs()))["peroxide arm"]
    check("the within-run side is the peroxide arm and holds [S] fixed per run",
          all(block.s0.nunique() == 1
              for _, block in arm.groupby("experiment")),
          "so its substrate order is the one the guard refuses")


def test_the_acceleration_band_uses_live_curves_only():
    """A dead curve's `accelerates` is its quantisation staircase stepping late."""
    print("\nthe acceleration bands")
    bands = scope.acceleration_by_ph()
    frame = scope.frame()
    check("the bands cover every live curve and no dead one",
          int(bands.curves.sum()) == int(frame.live.sum()),
          f"{int(bands.curves.sum())} against {int(frame.live.sum())}")
    check("acceleration is commoner above pH 9",
          bands.loc["pH >= 9", "share"] > bands.loc["pH < 9", "share"],
          f"{bands.loc['pH >= 9', 'share']:.2f} against "
          f"{bands.loc['pH < 9', 'share']:.2f}")


if __name__ == "__main__":
    test_an_axis_the_offsets_absorb_is_not_reported()
    test_the_block_still_measures_both_axes()
    test_each_arm_agrees_with_the_joint_fit()
    test_a_ph_ladder_is_derived_and_not_listed()
    test_the_weak_runs_flatten_the_ph_ladder()
    test_the_schedule_control_is_a_real_control()
    test_the_two_levers_are_measured_on_different_contrasts()
    test_the_acceleration_band_uses_live_curves_only()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
