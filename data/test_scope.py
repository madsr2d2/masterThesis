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


def test_the_gas_is_gas_and_not_the_instrument():
    """
    The three controls that make the chop O2 rather than anything else.

    A ladder in [H2O2] alone would also fit a peroxide-dependent CHEMISTRY, so
    the ladder is not enough on its own; the two controls beside it are what
    close it. Turnover: two runs sit at the block's top peroxide and carry no
    detachment, and they are its two weakest. Synchrony: a lamp, a shutter or
    the carousel would take all seven cuvettes of a run at once.
    """
    print("\nthe chop is gas")
    ladder = scope.bubble_ladder()
    shares = (ladder.with_drops / ladder.curves).to_numpy()
    check("no curve at the bottom of the peroxide range detaches",
          shares[0] == 0.0, f"{shares[0]:.2f}")
    check("every curve at the top of it does",
          shares[-1] == 1.0, f"{shares[-1]:.2f}")
    check("and the drop count rises monotonically between",
          bool((np.diff(ladder.mean_drops.to_numpy()) > 0).all()),
          ", ".join(f"{v:.2f}" for v in ladder.mean_drops))

    runs = scope.bubble_turnover_control()
    top = runs[np.isclose(runs.top_h2o2, 73.424)]
    quiet = top[top.drops == 0]
    check("peroxide alone does not do it -- some top-peroxide runs are quiet",
          len(quiet) > 0, f"{sorted(quiet.index)}")
    check("and the quiet ones are the weakest runs at that peroxide",
          float(quiet.agreement.max()) < float(
              top[top.drops > 0].agreement.min()),
          f"{quiet.agreement.max():.2f} against "
          f"{top[top.drops > 0].agreement.min():.2f}")

    together = scope.bubble_synchrony()
    check("detachments are not synchronised between cuvettes of a run",
          abs(together["observed"] - together["expected"])
          < 3 * np.sqrt(max(together["expected"], 1.0)),
          f"{together['observed']} observed against "
          f"{together['expected']:.1f} expected over {together['pairs']} pairs")


def test_stitching_is_refused_on_the_real_curves():
    """
    Mass balance, on the archive rather than on a synthetic curve.

    Joining the pieces keeps the artefact's rise and drops its fall, so it
    inflates by the sum of the drops. On this block that carries at least one
    curve past the most absorbance its own substrate could ever make -- which
    the subtraction, being conserving, cannot do to any curve.
    """
    print("\nstitching, against the substrate")
    balance = scope.bubble_mass_balance()
    check("no curve as read exceeds its substrate ceiling",
          int((balance.raw > 1.0).sum()) == 0)
    check("stitching pushes at least one curve past it",
          int((balance.stitched > 1.0).sum()) > 0,
          f"worst {balance.stitched.max():.2f} on exp "
          f"{int(balance.loc[balance.stitched.idxmax(), 'experiment'])} "
          f"sample {int(balance.loc[balance.stitched.idxmax(), 'sample'])}")
    check("the subtraction leaves every curve inside it",
          int((balance.corrected > 1.0).sum()) == 0)
    # It cannot do otherwise, and saying so here is the point: the gas starts
    # at nothing and ends at nothing or more, so the repair can only ever
    # LOWER the net rise. Stitching is the one direction the mass balance
    # forbids, and it is the direction stitching always goes.
    check("because the repair can only lower a curve, never raise it",
          bool((balance.corrected <= balance.raw + 1e-12).all()),
          f"worst {float((balance.corrected - balance.raw).max()):+.2e}")
    check("and stitching only ever raises one",
          bool((balance.stitched >= balance.raw - 1e-12).all()))
    # The curves stitching ruins are the ones that were still speeding up, so
    # the conversion it claims cannot be believed on its own terms either.
    worst = balance.nlargest(5, "stitched")
    check("and the curves it ruins were accelerating, not spending substrate",
          float(worst.ramp_gain.min()) > 1.0,
          f"gains {', '.join(f'{v:.1f}' for v in worst.ramp_gain)}")


def test_the_correction_recovers_a_planted_rate():
    """
    The recovery table, on real curves with a known truth.

    Donors are the block's own clean curves; the artefact is planted, so the
    rate before it is the truth. It is planted TWO WAYS, because the choice
    between them is the whole difference between the segment ramp and the
    model that replaced it: with each bubble emptying at its detachment, which
    is what the ramp assumed, and with each emptying only partly, which is what
    `bubble_record(141, 3)` says actually happens. A repair that only works
    under its own assumption is not a repair.

    And scored two ways within each, because the reconstruction subtracts only
    gas it watched leave. `ends_holding=False` stops production at the last
    release, so all the planted gas went and the repair is answerable for all
    of it; `ends_holding=True` leaves the run still making gas at the last
    reading, so it ends holding a bubble nothing ever saw go and the repair
    keeps it. The FIRST is the recovery. The second is the stated systematic,
    and the gap between them is its size.
    """
    print("\nrecovery of a planted rate")
    for emptying in (True, False):
        how = "emptying" if emptying else "partly emptying"
        table = scope.bubble_recovery(emptying=emptying)
        shed = scope.bubble_recovery(emptying=emptying, ends_holding=False)
        check(f"{how}: there are clean donor curves to plant into",
              int(table.n.iloc[0]) >= 40, f"{int(table.n.iloc[0])} plantings")
        for severity, row in table.iterrows():
            check(f"{how}, {severity:g}x: stitching never beats leaving it "
                  f"alone", row.stitched >= row.raw - 1e-9,
                  f"{row.stitched:.4f} against {row.raw:.4f}")
            if severity >= 0.5:
                # At 0.25x with the bubbles emptying only partly the two are
                # equal to four decimals: barely a drop clears the detector,
                # so there is nothing for stitching to add back. From 0.5x on
                # it is strictly worse.
                check(f"{how}, {severity:g}x: and past the smallest artefact "
                      f"it is strictly worse", row.stitched > row.raw,
                      f"{row.stitched:.4f} against {row.raw:.4f}")
            check(f"{how}, {severity:g}x: the reconstruction beats both",
                  abs(row.rebuilt - 1) < abs(row.raw - 1)
                  and abs(row.rebuilt - 1) < abs(row.stitched - 1),
                  f"{row.rebuilt:.2f} against {row.raw:.2f} and "
                  f"{row.stitched:.2f}")
            # THE CLAIM THE SEGMENT RAMP COULD NOT MAKE. It held to 1.01 while
            # the artefact was small and reached 1.60 by 2x; this stays inside
            # a tenth at every severity, under both plantings -- on the gas
            # that was seen to leave.
            check(f"{how}, {severity:g}x: and on gas that left it is within a "
                  f"tenth of the truth",
                  abs(shed.rebuilt[severity] - 1.0) < 0.10,
                  f"{shed.rebuilt[severity]:.2f}")
            # The systematic runs one way: keeping a bubble means removing too
            # little, never too much.
            check(f"{how}, {severity:g}x: and the bubble it keeps only ever "
                  f"raises the answer",
                  row.rebuilt >= shed.rebuilt[severity] - 0.01,
                  f"{row.rebuilt:.2f} against {shed.rebuilt[severity]:.2f}")
        check(f"{how}: and the systematic is real at the largest artefact",
              table.rebuilt.iloc[-1] - shed.rebuilt.iloc[-1] > 0.2,
              f"{table.rebuilt.iloc[-1]:.2f} against "
              f"{shed.rebuilt.iloc[-1]:.2f}")

    # A repair that moves a curve it was not needed on is a repair that has to
    # be defended on every curve. This one is the identity there.
    untouched = []
    for curve in scope.curves():
        values = np.asarray(curve.absorbance, dtype=float)
        if len(curve_metrics.bubble_drops(values, curve.noise)):
            continue
        untouched.append(np.array_equal(curve_metrics.debubble(
            np.asarray(curve.times, dtype=float), values, curve.noise)[0],
            values))
    check("and on a curve with no detachment it changes nothing at all",
          all(untouched) and len(untouched) > 40,
          f"{sum(untouched)} of {len(untouched)}")


def test_every_detachment_is_corrected_and_nothing_else_is():
    """
    The repair removes gas, and ONLY gas.

    Two halves. Every detachment has to be gone from the reconstruction --
    that is what the segment ramp failed, leaving five consecutive steps past
    8 sigma in exp 141 cuvette 3 and the whole of a -0.0165 AU fall at 60
    sigma in exp 144 cuvette 2. And a fall the model decided was NOT gas has
    to be left exactly where it is, because laundering an instrument
    excursion through a gas model is how a real early rise gets removed.

    So the guarantee is `worst_at_event` and not `rebuilt_worst`. The
    reconstructions still fall by up to -61.1 sigma somewhere, and the curves
    doing that are behaving correctly: those are the 34 rejected excursions,
    which `isolated_outliers` nominates and nothing removes automatically.
    """
    print("\nevery detachment corrected, and nothing else touched")
    table = scope.rebuild_smoothness()
    repaired = table[~table.clean]
    frame = scope.frame().set_index(["experiment", "sample"])
    untouched = [(int(r.experiment), int(r["sample"]))
                 for _, r in repaired.iterrows()
                 if not np.isfinite(frame.loc[(int(r.experiment),
                                               int(r["sample"])), "gas_rate"])]
    check("there are curves carrying gas to repair", len(repaired) > 20,
          f"{len(repaired)} of {len(table)} live")
    check("the readings carry falls no noise can explain",
          float(repaired.raw_worst.min()) < -100,
          f"{repaired.raw_worst.min():.1f} sigma")

    fixable = repaired[[(int(r.experiment), int(r["sample"])) not in untouched
                        for _, r in repaired.iterrows()]]
    check("every detachment is corrected in full",
          float(fixable.worst_at_event.min()) >= -1e-9,
          f"worst {fixable.worst_at_event.min():+.3f} sigma over "
          f"{len(fixable)} curves")
    check("and the only curve left uncorrected is the first-interval one",
          untouched == [(135, 6)], f"{untouched}")

    # THE OTHER HALF. A fall that is not gas must survive the repair, so that
    # it stays visible as the instrument problem it is.
    check("falls rejected as excursions are left in place",
          int(table.excursions.sum()) > 0,
          f"{int(table.excursions.sum())} rejected")
    check("and they are why the rebuilt curves still fall somewhere",
          float(repaired.rebuilt_worst.min()) < -20,
          f"{repaired.rebuilt_worst.min():.1f} sigma")
    for experiment, sample in ((149, 5), (149, 1)):
        row = table[(table.experiment == experiment)
                    & (table["sample"] == sample)].iloc[0]
        check(f"exp {experiment} cuvette {sample} is left entirely alone",
              int(row.bubble_events) == 0 and int(row.excursions) > 0,
              f"{int(row.bubble_events)} detachments, "
              f"{int(row.excursions)} excursions")
        check(f"  and its reconstruction is the readings",
              abs(row.raw_worst - row.rebuilt_worst) < 1e-9)


def test_the_excursion_test_on_the_curve_that_forced_it():
    """
    Exp 149 cuvette 5, the real curve the excursion test was written for.

    Its two "detachments" are 9.3 and 8.2 sigma and neither is gas. The first
    falls 0.00206 AU and the NEXT READING climbs 0.00222 straight back; the
    second falls off a reading that is an isolated spike. Between them they
    set a production rate of 6.2e-6 AU/s, and the repair then removed 0.0097
    AU from a curve that rose 0.0262 -- turning a real early rise into a flat
    line, while staying perfectly monotone and passing every test there was.
    """
    print("\na fall that comes straight back is not gas")
    curve = {c.sample: c for c in scope.curves_of(149)}[5]
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    falls = curve_metrics.bubble_drops(values, curve.noise)
    check("the detector still sees both falls", len(falls) == 2, f"{len(falls)}")
    check("and they are past the threshold on their own",
          float(np.diff(values)[falls].min() / curve.noise)
          < -curve_metrics.BUBBLE_DROP_SIGMA,
          f"{np.diff(values)[falls].min() / curve.noise:.1f} sigma")
    check("the first is undone by the very next reading",
          float(values[10] - values[9]) > float(values[8] - values[9]),
          f"+{values[10] - values[9]:.5f} against a fall of "
          f"{values[8] - values[9]:.5f}")
    check("neither survives as a detachment",
          curve_metrics.detachments(values, curve.noise) == [], "")
    rebuilt, _ = curve_metrics.debubble(times, values, curve.noise)
    check("so the curve is returned exactly as it was read",
          np.array_equal(rebuilt, values))

    # WITHOUT the test, the repair takes a third of the curve.
    unfiltered = curve_metrics.detachments(values, curve.noise,
                                           recovery=np.inf)
    rate = curve_metrics.bubble_rate(times, values, unfiltered)
    held = curve_metrics.bubble_profile(times, values, unfiltered, rate)
    check("and that is not a small correction it is refusing",
          held.max() / float(values[-1] - values[0]) > 0.3,
          f"{held.max():.4f} AU off a rise of {values[-1] - values[0]:.4f}")

    # The one-sided shape of the test matters: local_outlier_z cannot do this,
    # because its window spans the fall and a genuine step flags itself.
    worst = max(scope.curves(), key=lambda c: float(
        np.max(np.maximum.accumulate(np.asarray(c.absorbance, dtype=float))
               - np.asarray(c.absorbance, dtype=float))))
    z = curve_metrics.local_outlier_z(
        np.asarray(worst.times, dtype=float),
        np.asarray(worst.absorbance, dtype=float), worst.noise)
    events = curve_metrics.detachments(
        np.asarray(worst.absorbance, dtype=float), worst.noise)
    check("a genuine step change flags ITSELF under local_outlier_z",
          any(z[start] > curve_metrics.OUTLIER_SIGMA for start, _ in events),
          "which is why the excursion test never looks across the fall")


def test_the_gas_may_not_outlast_the_evidence():
    """
    A MONOTONE RECONSTRUCTION CAN STILL BE THE WRONG ONE.

    The smoothness test above passes on a curve that has been dragged far too
    low, because pulling a curve down by a smooth ramp leaves it smooth. That
    is exactly what happened: the production rate is fixed by the detachments,
    and where those are early and the run is long it was extrapolated across
    hours in which nothing detached. Exp 149 cuvette 4 sheds 0.0031 AU once,
    1920 s into an 8 h run, and the uncapped profile grew to 0.0273 over the
    remaining 7.5 h -- 8.8x the largest bubble that curve ever shed -- taking
    the reconstruction to -0.0209 against a raw rise of +0.0064. Twelve of 49
    curves held more than twice their own largest bubble and three finished
    below zero, and every one of them passed the smoothness test.

    Capping the tail at the most the beam had carried bounded that without
    curing it: the reconstruction still sat a flat 0.0022 AU under exp 149
    cuvette 3's own readings for 82% of the run, and under exp 150 cuvette 1's
    by 99% of everything it rose. `unreleased_gas` is the rule that cures it --
    THE BEAM MAY HOLD AT MOST THE DETACHMENTS STILL TO COME -- so a reading
    after the last fall carries no gas at all and every reconstruction lands
    back ON the readings. This is the test that would have caught it, and the
    curves are named so it cannot come back quietly.
    """
    print("\nthe gas may not outlast the evidence for it")
    table = scope.rebuild_smoothness()
    repaired = table[~table.clean]
    check("no reconstruction ends holding gas",
          float(repaired.gas_at_end.abs().max()) == 0.0,
          f"worst {repaired.gas_at_end.abs().max():.2e}")
    ratio = repaired.gas_held / repaired.biggest_bubble
    check("no curve holds many times its own largest bubble",
          float(ratio.max()) < 6.0, f"worst {ratio.max():.1f}x")
    check("and the typical curve holds about one",
          float(ratio.median()) < 2.0, f"median {ratio.median():.1f}x")
    for experiment, sample in ((149, 2), (149, 3), (149, 4), (150, 1)):
        row = table[(table.experiment == experiment)
                    & (table["sample"] == sample)]
        if not len(row):
            continue
        check(f"exp {experiment} cuvette {sample} ends on its readings",
              float(row.gas_at_end.iloc[0]) == 0.0,
              f"{row.gas_at_end.iloc[0]:+.2e} AU below them")

    # THE PRICE, STATED. Subtracting only gas that was watched to leave means
    # a run still making gas when the recording stopped keeps the bubble it
    # never shed. Against a planting that ends at a release the repair is
    # exact; against one that ends mid-bubble it is not, and the gap between
    # the two is the systematic rather than a failure.
    shed = scope.bubble_recovery(ends_holding=False).rebuilt
    holding = scope.bubble_recovery(ends_holding=True).rebuilt
    check("gas that was seen to leave is removed in full",
          float(np.abs(shed - 1.0).max()) < 0.05,
          ", ".join(f"{x:.2f}" for x in shed))
    check("and gas that never left is left in, which is the stated cost",
          float(holding.max()) > 1.2,
          ", ".join(f"{x:.2f}" for x in holding))
    check("the curves at risk of it are the ones with a short quiet tail",
          int((repaired.quiet_tail > 1).sum()) >= 10,
          f"{int((repaired.quiet_tail > 1).sum())} of {len(repaired)} run "
          f"more than a full shedding interval past their last detachment")

    # A reconstruction that ends below where it started has removed more than
    # the curve ever rose. One curve is allowed to graze it and it has to be
    # the one the load already flags.
    negative = table[table.rebuilt_net < 0]
    check("at most one curve rebuilds to a falling curve",
          len(negative) <= 1, f"{len(negative)}")
    for _, row in negative.iterrows():
        check(f"and exp {int(row.experiment)} cuvette {int(row['sample'])} "
              f"is flagged by its load",
              row.bubble_load > scope.BUBBLE_LOAD_CEILING,
              f"load {row.bubble_load:.2f}")
        check("and it is zero within the curve's own noise",
              abs(row.rebuilt_net) < 0.001, f"{row.rebuilt_net:+.4f}")

    check("every live curve still carries a rate after the repair",
          int((scope.frame().query("live").vmax_corrected <= 0).sum()) == 0)


def test_the_gas_rate_belongs_to_the_peroxide():
    """
    `bubble_rate` never sees a concentration, so what it correlates with is a
    prediction rather than a fit. It is first order in peroxide.
    """
    print("\nwhat the fitted gas rate depends on")
    drivers = scope.gas_rate_drivers()
    check("the gas rate is first order in peroxide",
          abs(drivers["pooled_h2o2"] - 1.0) < 2 * drivers["pooled_stderr_h2o2"],
          f"{drivers['pooled_h2o2']:+.3f} +/- "
          f"{drivers['pooled_stderr_h2o2']:.3f}")
    check("and that is not nothing -- it is excluded from zero",
          drivers["pooled_h2o2"] > 3 * drivers["pooled_stderr_h2o2"],
          f"{drivers['pooled_h2o2'] / drivers['pooled_stderr_h2o2']:.1f} sigma")
    check("the substrate carries far less of it than the peroxide does",
          abs(drivers["order_s0"]) < 0.5 * abs(drivers["pooled_h2o2"]),
          f"{drivers['order_s0']:+.3f} against {drivers['pooled_h2o2']:+.3f}")


def test_no_published_order_rests_on_the_gas():
    """
    The reason the block's conclusions survive a defect this large.

    A substrate-blind, peroxide-driven additive artefact is exactly what would
    manufacture the flat substrate order of section 2, so the orders have to
    be read under every repair before that order can stand.

    THE SUBSTRATE ORDER IS THE ONE AT RISK AND IT DOES NOT MOVE. The peroxide
    order does: the reconstruction takes it from +0.794 to +0.666, which is
    1.2 sigma of its own error and the largest shift any repair here produces.
    That is not a failure of the repair, it is the repair working -- the gas
    is made from peroxide, so an uncorrected artefact must inflate the
    peroxide order, and taking the gas out must bring it down. It is reported
    rather than asserted away, and it is not significant.
    """
    print("\nthe orders under each repair")
    table = scope.bubble_sensitivity()
    published = table.loc[("vmax", "all live")]

    def shift(row, axis):
        gap = abs(row[f"order_{axis}"] - published[f"order_{axis}"])
        spread = float(np.hypot(row[f"stderr_{axis}"],
                                published[f"stderr_{axis}"]))
        return gap, spread, gap / spread

    for (treatment, subset), row in table.iterrows():
        if treatment == "vmax" and subset == "all live":
            continue
        gap, spread, sigma = shift(row, "s0")
        check(f"{treatment} on {subset}: the substrate order does not move",
              gap < spread, f"{gap:.3f} against {spread:.3f}")
        gap, spread, sigma = shift(row, "h2o2")
        check(f"{treatment} on {subset}: the peroxide order moves less than "
              f"1.5 sigma", sigma < 1.5, f"{sigma:.1f} sigma")

    corrected = table.loc[("vmax_corrected", "all live")]
    check("the peroxide order comes DOWN under the reconstruction, as an "
          "artefact made from peroxide requires",
          corrected["order_h2o2"] < published["order_h2o2"],
          f"{corrected['order_h2o2']:+.3f} against "
          f"{published['order_h2o2']:+.3f}")
    check("and it is the largest shift of the three repairs",
          shift(corrected, "h2o2")[2]
          > max(shift(table.loc[(t, "all live")], "h2o2")[2]
                for t in ("vmax_monotone",)),
          f"{shift(corrected, 'h2o2')[2]:.1f} sigma")
    check("the substrate order does not move toward zero under the bound",
          abs(table.loc[("vmax_monotone", "all live"), "order_s0"])
          >= abs(published["order_s0"]))

    # ONE curve of 110 loses its rate to the repair, and it has to be named
    # rather than counted: a rate that goes negative when the gas comes out
    # had no rate to begin with.
    lost = int(published["n"]) - int(corrected["n"])
    check("at most one live curve loses a positive rate to the repair",
          lost <= 1, f"{lost}")
    data = scope.frame()
    data = data[data.live & (data.vmax_corrected <= 0)]
    check("and it is a weak run, outside every quoted number in the block",
          set(data.experiment) <= set(scope.TWO_AXIS_BLOCK)
          - set(scope.strong_runs()),
          f"{[(int(r.experiment), int(r['sample'])) for _, r in data.iterrows()]}")


def test_the_load_ceiling_flags_and_does_not_exclude():
    print("\nthe load ceiling")
    table = scope.bubble_table()
    frame = scope.frame()
    check("the table holds every live curve",
          len(table) == int(frame.live.sum()),
          f"{len(table)} against {int(frame.live.sum())}")
    check("some curves are beyond repair", int((~table.repairable).sum()) > 0,
          f"{int((~table.repairable).sum())}")
    check("and they are the high-peroxide ones",
          float(table[~table.repairable].h2o2.median())
          > float(table[table.repairable].h2o2.median()),
          f"{table[~table.repairable].h2o2.median():.1f} against "
          f"{table[table.repairable].h2o2.median():.1f}")
    check("the ceiling is a flag, not an exclusion -- the frame keeps them all",
          int(frame.live.sum()) == len(table))
    check("a curve with no detachment carries no load",
          float(table[table.bubble_drops == 0].bubble_load.max()) == 0.0)


def test_the_large_steps_say_which_beam_the_gas_is_in():
    """
    The sign of a large step separates the two cuvettes.

    A bubble in the SAMPLE beam scatters light out of the aperture, so the
    difference climbs while it grows and drops when it goes: slow up, sudden
    down. A bubble in the REFERENCE beam does the same to the reference, and
    the subtraction inverts it. So the two populations are separable by sign,
    and if they were balanced the reading would be in doubt.
    """
    print("\nwhich beam the gas is in")
    asymmetry = scope.bubble_step_asymmetry()
    check("large steps fall far more often than they rise",
          asymmetry["falls"] > 3 * asymmetry["rises"],
          f"{asymmetry['falls']} falls against {asymmetry['rises']} rises, "
          f"ratio {asymmetry['ratio']:.1f}")
    check("and the largest fall dwarfs the largest rise",
          abs(asymmetry["largest_fall"]) > 3 * asymmetry["largest_rise"],
          f"{asymmetry['largest_fall']:.0f} against "
          f"+{asymmetry['largest_rise']:.0f}")
    check("the tail threshold is well above the detachment threshold",
          scope.BUBBLE_STEP_SIGMA > 2 * curve_metrics.BUBBLE_DROP_SIGMA,
          f"{scope.BUBBLE_STEP_SIGMA:g} against "
          f"{curve_metrics.BUBBLE_DROP_SIGMA:g}")
    counted = sum(max(len(np.asarray(c.absorbance, float)) - 1, 0)
                  for c in scope.curves())
    check("every step in the block is counted exactly once",
          asymmetry["steps"] == counted,
          f"{asymmetry['steps']} against {counted}")
    # And the statistic has to be able to see the other sign, or it is not a
    # test. A curve mirrored top to bottom must swap the two counts.
    values = np.asarray(scope.curves_of(141)[2].absorbance, dtype=float)
    noise = scope.curves_of(141)[2].noise
    forward = np.diff(values) / noise
    mirrored = np.diff(-values) / noise
    check("mirroring a curve swaps its rises and falls",
          int((mirrored > scope.BUBBLE_STEP_SIGMA).sum())
          == int((forward < -scope.BUBBLE_STEP_SIGMA).sum()))


def test_the_worked_curve_carries_more_than_one_bubble():
    """
    `bubble_record` on the curve section 5 walks through.

    The repair charges each drop to growth since the previous drop. That is
    right on average and this curve is the counter-example: the largest drop
    follows the SHORTEST growth window, and the level after the last release
    sits below every earlier one -- which a monotone chemistry beneath a single
    bubble cannot produce. The test asserts the counter-example, because it is
    the honest limit of the correction and a document quoting it must not be
    able to drift away from it.
    """
    print("\nthe worked curve carries more than one bubble")
    experiment, sample = scope.BUBBLE_WORKED_EXAMPLE
    record = scope.bubble_record(experiment, sample)
    check("the record has a row per detachment",
          len(record) == int(scope.frame()[
              (scope.frame().experiment == experiment)
              & (scope.frame()["sample"] == sample)].bubble_drops.iloc[0]),
          f"{len(record)} rows")
    curve = {c.sample: c for c in scope.curves_of(experiment)}[sample]
    values = np.asarray(curve.absorbance, dtype=float)
    drops = curve_metrics.bubble_drops(values, curve.noise)
    check("and its drops are the readings' own",
          np.allclose(record["lost"].to_numpy(), -np.diff(values)[drops]))
    check("every drop clears the detachment threshold",
          bool((record.sigma <= -curve_metrics.BUBBLE_DROP_SIGMA).all()),
          f"worst {record.sigma.max():.1f}")

    check("the largest drop follows the shortest growth window",
          record.loc[record["lost"].idxmax(), "grew_s"] == record.grew_s.min(),
          f"{record.loc[record['lost'].idxmax(), 'grew_s']:.0f} s of growth shed "
          f"{record["lost"].max():.4f} AU")
    check("so the drops do not scale with the time available to grow",
          float(np.corrcoef(record.grew_s, record["lost"])[0, 1]) < 0.5,
          f"r = {float(np.corrcoef(record.grew_s, record["lost"])[0, 1]):+.2f}")
    check("and the last level sits below every earlier one",
          record.after.iloc[-1] < record.after.iloc[:-1].min(),
          f"{record.after.iloc[-1]:+.4f} against "
          f"{record.after.iloc[:-1].min():+.4f}")
    check("which a monotone curve under one bubble could not do",
          not record.after.is_monotonic_increasing)
    check("the curve sheds more than it nets", record["lost"].sum() >
          float(values[-1] - values[0]),
          f"{record["lost"].sum():.4f} against {values[-1] - values[0]:.4f}")


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
    test_the_gas_is_gas_and_not_the_instrument()
    test_stitching_is_refused_on_the_real_curves()
    test_the_correction_recovers_a_planted_rate()
    test_no_published_order_rests_on_the_gas()
    test_the_load_ceiling_flags_and_does_not_exclude()
    test_the_large_steps_say_which_beam_the_gas_is_in()
    test_the_worked_curve_carries_more_than_one_bubble()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
