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

import curve_metrics
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


def test_the_burst_amplitude_is_read_off_predictions_not_parameters():
    """
    Planted so the parameters are degenerate and the prediction is not.

    Two exponentials with nearly the same time constant let the linear solve
    trade enormous opposite amplitudes between them without moving the fitted
    curve -- exp 135 sample 3 comes back with B_fast = -241 against B_slow =
    +303 on a curve that moves 0.06 AU. Any burst read off `B_fast`, or off
    their sum, inherits that. Reading the curve does not, and this test is the
    statement of why.
    """
    print("\nthe burst is a prediction, not a parameter")
    times = np.linspace(0.0, 1000.0, 400)
    rise = 0.02 * (1 - np.exp(-times / 100.0)) - 3e-6 * times
    amplitude, when, bounded = curve_metrics.burst_amplitude(times, rise)
    check("a curve that turns over inside the run is bounded", bounded,
          f"peak at {when:.0f} s of {times[-1]:.0f}")
    check("and its amplitude is the height it reached",
          abs(amplitude - (rise.max() - rise[0])) < 1e-12,
          f"{amplitude:.6f}")

    climbing = 0.02 * (1 - np.exp(-times / 100.0))
    _, _, still_rising = curve_metrics.burst_amplitude(times, climbing)
    check("a curve still rising at the last reading is NOT bounded",
          not still_rising, "so it cannot be compared between runs")

    # The degeneracy itself: two shapes with identical predictions and wildly
    # different parameter splits must give the same burst.
    slow = 0.02 * (1 - np.exp(-times / 100.0)) - 3e-6 * times
    traded = (slow + 5.0 * (1 - np.exp(-times / 1e9))
              - 5.0 * (1 - np.exp(-times / 1e9)))
    check("and a parameter trade that leaves the curve alone changes nothing",
          abs(curve_metrics.burst_amplitude(times, traded)[0] - amplitude)
          < 1e-9, "predictions are what is read")


def test_the_burst_bound_bites_between_runs_and_not_within_them():
    """
    Which is why `burst_drivers` may use every live curve and `enzyme_pair` may not.

    Every cuvette of a run shares its length, so a truncated rise is truncated
    identically across the run and a per-experiment offset absorbs it. Between
    runs it is not, and the block spans 3000 to 28740 s.
    """
    print("\nwhere the bound bites")
    frame = scope.frame()
    live = frame[frame.live]
    # NOT exactly constant, and the reason is the instrument rather than the
    # design: it reads the seven cuvettes in sequence, so six runs end with
    # some cuvettes one 60 s interval short of the others. That is 0.7-1.6%,
    # against 9.6x between runs, so the argument survives -- but it is 1.6%
    # and not zero, and asserting zero here failed the first time it was run.
    inside = live.groupby("experiment").duration_s.agg(lambda d: d.max() / d.min())
    check("run length is constant inside a run to within one reading",
          bool((inside <= 1.02).all()),
          f"worst {inside.max():.4f}x in exp {inside.idxmax()}")
    check("and varies a long way between them",
          live.duration_s.max() / live.duration_s.min() > 5,
          f"{live.duration_s.max() / live.duration_s.min():.1f}x")
    check("so the within-run spread is two orders below the between-run one",
          (inside.max() - 1) * 100 < live.duration_s.max() / live.duration_s.min(),
          f"{(inside.max() - 1) * 100:.1f}% against "
          f"{live.duration_s.max() / live.duration_s.min():.0f}00%")

    drivers = scope.burst_drivers()
    check("so the burst's concentration orders use every live curve",
          drivers["curves"] == int(live.burst.gt(0).sum()),
          f"{drivers['curves']} of {int(live.live.sum())}")
    check("most of which never finish their rise",
          drivers["bounded"] < drivers["curves"] / 2,
          f"{drivers['bounded']} bounded of {drivers['curves']}")
    check("and the catalyst is not identified there",
          not drivers["enzyme_identified"],
          "[enz] never moves inside a run")


def test_the_enzyme_pair_is_derived_and_survives_its_window():
    """
    The pair, the correction it needs, and the sweep that says it is not a window.

    The pair is chosen by the smallest pH gap among runs that share a
    composition and differ in loading -- not listed -- so a change in the
    archive moves it rather than silently invalidating it.
    """
    print("\nthe enzyme lever")
    table, verdict = scope.enzyme_pair()
    check("the pair shares a composition and steps the catalyst",
          verdict["expected"] > 2.0, f"{verdict['expected']:.2f}x in [enz]")
    check("with the smallest pH gap the block offers",
          verdict["pH_gap"] < 0.15, f"{verdict['pH_gap']:.2f} pH units")
    check("and all seven cuvettes pair up", verdict["cuvettes"] == 7,
          f"{verdict['cuvettes']}")
    check("the pH gap is corrected, and is small beside the enzyme step",
          1.0 < verdict["pH_correction"] < 1.2,
          f"{verdict['pH_correction']:.3f} against "
          f"{verdict['expected']:.2f}")

    check("the rise scales with the catalyst",
          verdict["sigma_no_dependence"] > 3.0,
          f"no dependence excluded at {verdict['sigma_no_dependence']:.1f} sigma")
    check("and is not distinguishable from first order in it",
          verdict["sigma_first_order"] < 2.0,
          f"{verdict['sigma_first_order']:.1f} sigma from {verdict['expected']:.2f}x")

    sweep = scope.enzyme_pair_sensitivity()
    check("every window in the sweep says the same thing",
          bool((sweep.order > 0.5).all() and (sweep.order < 1.6).all()),
          f"{sweep.order.round(2).to_dict()}")
    check("and the published window is the largest one, not a chosen one",
          verdict["window_s"] == sweep.window_s.max(),
          f"{verdict['window_s']:.0f} s")


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
    test_the_burst_amplitude_is_read_off_predictions_not_parameters()
    test_the_burst_bound_bites_between_runs_and_not_within_them()
    test_the_enzyme_pair_is_derived_and_survives_its_window()
    test_the_acceleration_band_uses_live_curves_only()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
