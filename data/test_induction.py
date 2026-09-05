"""
Tests for induction.py.

Three of these would have caught a real mistake. The landmark has to return a
number for a curve with no induction rather than nothing, because dropping
those rows is the censoring the statistic exists to avoid. The driver
regression has to read back BOTH a planted clock and a planted product
threshold, because a regression that always says "clock" would agree with every
conclusion in this module for the wrong reason. And `scope.orders` has to
refuse an axis the block does not move: exps 127-131 hold [S] fixed and
returned an order of +2.14 +- 0.18 in it before the guard existed.

    python data/test_induction.py
"""
import sys
import os

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import induction
import scope
from curve_metrics import lag_time
# `FakeCurve` is test_slowdown's, imported rather than written again: it
# already carries the four attributes `induction_point` reads, and two copies
# of a fixture drift the way two copies of a measurement do.
from test_slowdown import FakeCurve
import pandas as pd
from induction import (DEPTH_FLOOR, INDUCTION_CLOCK_SLOPE, INDUCTION_FLOOR,
                       INDUCTION_PRODUCT_SLOPE, PERHYDRATE_ORDER_GAP,
                       common_window, induction_blocks,
                       induction_drivers, induction_point,
                       joint_buffer_order, joint_order,
                       joint_peroxide_order,
                       PRODUCT_CLOCK_SLOPE, PRODUCT_THRESHOLD_SLOPE,
                       induction_table, product_at_landmark,
                       product_recovery,
                       lag_arrhenius_sweep, lag_ladder, lag_orders,
                       lag_ph_ladders, lag_signal_control, lag_window_frame,
                       order_ratio, pooled_ladder,
                       peroxide_geometric_mean, peroxide_saturation,
                       composition_collinearity, ladder_arms,
                       lag_identifiability, replicate_floor, WHOLE_ARCHIVE,
                       saturation_fraction, sign_drivers,
                       signal_control, substrate_lever, trap_constant)
from summary_kinetics import fit_progress

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def lag_curve(v=1e-5, tau=800.0, span=8000.0, points=400, noise=0.0, seed=0):
    """A = v(t - tau(1 - e^(-t/tau))): rate 0 at t=0, rising to v."""
    times = np.linspace(0.0, span, points)
    values = v * (times - tau * (1.0 - np.exp(-times / tau)))
    if noise:
        values = values + np.random.default_rng(seed).normal(0.0, noise, points)
    return FakeCurve(times, values)


def test_the_landmark_measures_what_it_claims():
    """Half-rise of a known relaxation, and zero where there is no rise."""
    print("\nthe landmark")
    tau = 800.0
    found = induction_point(lag_curve(tau=tau))
    # The rate is v(1 - e^(-t/tau)), so it reaches half of v at tau.ln2.
    # The rolling window averages that over a tenth of the run, which pushes
    # the crossing later; a factor of two is the tolerance this earns.
    expected = tau * np.log(2.0)
    check("a pure lag's half-rise is within a factor of two of tau.ln2",
          0.5 * expected <= found.t_ind <= 2.5 * expected,
          f"{found.t_ind:.0f} s against {expected:.0f} s")
    # depth is smeared by the window: a relaxation as long as the window reads
    # 1 - (tau/w)(1 - e^(-w/tau)), which for tau = w is 0.632, not 1.
    window = 0.1 * 8000.0
    for relaxation in (tau, 4000.0):
        smeared = (relaxation / window) * (1.0 - np.exp(-window / relaxation))
        seen = induction_point(lag_curve(tau=relaxation)).depth
        check(f"depth reads the window-smeared amplitude at tau = {relaxation:.0f}",
              abs(seen - smeared) < 0.02, f"{seen:.3f} against {smeared:.3f}")
    check("and a long relaxation recovers most of it",
          induction_point(lag_curve(tau=4000.0)).depth > 0.85,
          f"{induction_point(lag_curve(tau=4000.0)).depth:.3f}")

    times = np.linspace(0.0, 8000.0, 400)
    straight = FakeCurve(times, 1e-5 * times)
    found = induction_point(straight)
    check("a straight line has no induction, and returns 0 rather than nan",
          found.t_ind == 0.0 and found.depth == 0.0,
          f"t_ind {found.t_ind}, depth {found.depth}")
    shallow = lag_curve(tau=800.0, span=8000.0)
    shallow.absorbance = shallow.absorbance + 1e-5 * 500.0 * shallow.times
    check("and so does a rise below the depth resolution",
          induction_point(shallow).t_ind == 0.0,
          f"{induction_point(shallow).t_ind}, "
          f"depth {induction_point(shallow).depth:.5f}")

    # A decelerating curve: fastest in its first window, like every
    # enzyme-free 4OMe run in the archive.
    falling = FakeCurve(times, 5e-2 * (1.0 - np.exp(-times / 3000.0)))
    found = induction_point(falling)
    check("a decelerating curve has no induction either",
          found.t_ind == 0.0, f"{found.t_ind}")

    check("longer tau gives a longer induction",
          induction_point(lag_curve(tau=2000.0)).t_ind
          > induction_point(lag_curve(tau=400.0)).t_ind)
    check("the induction does not move when the rate does",
          abs(induction_point(lag_curve(v=4e-5, tau=tau)).t_ind
              - induction_point(lag_curve(v=1e-5, tau=tau)).t_ind) < 1.0,
          "the shape is v-independent by construction and must read that way")


def test_it_agrees_with_the_gated_statistic():
    """
    Where `curve_metrics.lag_time` is defined, the two must be the same number.

    They are the same landmark; `induction_point` only removes the acceleration
    gate. If they ever disagree, one of them has drifted.
    """
    print("\nagreement with curve_metrics.lag_time")
    from fit_dataset import source_floor
    gaps = []
    for curve in scope.curves(scope.TEMPERATURE_SERIES):
        gated = lag_time(curve.times, curve.absorbance,
                         floor=source_floor(curve.source))
        if not np.isfinite(gated):
            continue
        gaps.append(abs(gated - induction_point(curve).t_ind))
    check("both are defined on at least half the temperature series",
          len(gaps) >= 12, f"{len(gaps)} curves")
    check("and they agree exactly wherever both are defined",
          max(gaps) < 1e-9, f"largest gap {max(gaps):.3g} s")


def _planted(rule, rates, spans=(6000.0,), seed=11):
    """
    A frame whose induction time follows a known rule.

    `rule` is "clock" (a fixed 900 s whatever the rate) or "product" (a fixed
    0.02 AU of product, so the time is 0.02/rate). Each experiment holds one
    span and carries a four-rung rate ladder, which is the archive's design.
    """
    generator = np.random.default_rng(seed)
    rows = []
    for index, span in enumerate(spans):
        for rung, rate in enumerate(rates):
            noise = float(np.exp(generator.normal(0.0, 0.12)))
            time = 900.0 * noise if rule == "clock" else 0.02 / rate * noise
            rows.append({"experiment": 100 + index, "sample": rung + 1,
                         "live": True, "t_ind": min(time, span),
                         "peak_rate": rate, "v_peak": rate, "vmax": rate,
                         "s0": 1.85 * (rung + 1), "h2o2": 82.5,
                         "duration_s": span, "depth": 0.5,
                         "net": 0.05, "noise": 1e-5})
    return pd.DataFrame(rows)


def test_induction_drivers_tell_a_clock_from_a_product():
    """Both answers have to be readable, or neither answer means anything."""
    print("\nthe driver regression")
    rates = np.array([2e-6, 4e-6, 8e-6, 1.6e-5])
    spans = (4000.0, 8000.0, 16000.0)
    clock = induction_drivers(_planted("clock", rates, spans))
    product = induction_drivers(_planted("product", rates, spans))
    check("a planted clock reads as a clock",
          abs(clock["slope"] - INDUCTION_CLOCK_SLOPE) < 3 * clock["stderr"],
          f"{clock['slope']:+.3f} +- {clock['stderr']:.3f}")
    check("a planted clock excludes product control",
          abs(clock["slope"] - INDUCTION_PRODUCT_SLOPE) > 3 * clock["stderr"],
          f"{clock['slope']:+.3f} +- {clock['stderr']:.3f}")
    check("a planted product threshold reads as one",
          abs(product["slope"] - INDUCTION_PRODUCT_SLOPE)
          < 3 * product["stderr"],
          f"{product['slope']:+.3f} +- {product['stderr']:.3f}")
    check("a planted product threshold excludes the clock",
          abs(product["slope"] - INDUCTION_CLOCK_SLOPE) > 3 * product["stderr"],
          f"{product['slope']:+.3f} +- {product['stderr']:.3f}")

    # The same pair through the order ratio, which is the route with no
    # errors-in-variables. `v_peak` here is proportional to s0^1 by
    # construction, so the ratio is the induction's order over 1.
    ratio = order_ratio(_planted("product", rates, spans))
    check("the order ratio agrees with the driver regression on a product rule",
          abs(ratio["ratio"] - INDUCTION_PRODUCT_SLOPE)
          < 3 * ratio["ratio_stderr"],
          f"{ratio['ratio']:+.2f} +- {ratio['ratio_stderr']:.2f}")


def _peroxide_frame(constant, rule, levels=(2.5, 10.0, 40.0, 160.0),
                    runs=4, seed=3, baseline=4000.0, substrate=(8.0,),
                    substrate_order=0.0, unresolved=0):
    """
    Curves built from a known `K + H2O2 <=> KP`, one peroxide ladder per run.

    `rule` says what the induction is. "adduct": the relaxation is the approach
    to that equilibrium, 1/tau = k_f h + k_r, and the rate follows [KP]. "trap":
    the activation is unimolecular from FREE catalyst, 1/tau = k_act/(1 + K h),
    while the rate still follows [KP] -- which is the reading section 4a says
    the archive's sign would imply.

    `substrate` and `substrate_order` give the rate a factor in [S] that the
    CLOCK does not carry, which is what makes the substrate axis a control:
    the joint slope there is the substrate order and not the adduct's +1.
    `unresolved` plants that many rows per run whose fitted clock is nonsense
    and whose `tau_slow_resolved` is False, so a gate can be shown to bite.
    """
    generator = np.random.default_rng(seed)
    rows = []
    for run in range(runs):
        level = float(np.exp(generator.normal(0.0, 0.3)))
        sample = 0
        for peroxide in levels:
            for concentration in substrate:
                sample += 1
                saturation = (constant * peroxide
                              / (1.0 + constant * peroxide))
                rate = (3e-6 * level * saturation
                        * (concentration / substrate[0]) ** substrate_order
                        * float(np.exp(generator.normal(0.0, 0.05))))
                if rule == "adduct":
                    # k_r = 1/baseline, k_f = K.k_r, so 1/tau = (1 + K h)/
                    # baseline. `baseline` is large enough that no time lands
                    # on INDUCTION_FLOOR, which would break the identity below.
                    time = baseline / (1.0 + constant * peroxide)
                else:
                    time = baseline * (1.0 + constant * peroxide)
                noise = float(np.exp(generator.normal(0.0, 0.05)))
                rows.append({"experiment": 200 + run, "sample": sample,
                             "live": True, "t_ind": time * noise,
                             "peak_rate": rate, "v_peak": rate, "vmax": rate,
                             "tau_slow": time * noise,
                             "tau_slow_resolved": True,
                             "s0": concentration, "h2o2": peroxide,
                             "duration_s": 9000.0,
                             "depth": 0.5, "net": 0.05, "noise": 1e-5})
        for spare in range(unresolved):
            sample += 1
            # A clock the fit could not pin: the value is whatever the
            # optimiser stopped at, and the flag is the only thing that knows.
            rows.append({"experiment": 200 + run, "sample": sample,
                         "live": True, "t_ind": baseline,
                         "peak_rate": 3e-6, "v_peak": 3e-6, "vmax": 3e-6,
                         "tau_slow": baseline * 40.0 * (spare + 1),
                         "tau_slow_resolved": False,
                         "s0": substrate[0], "h2o2": levels[spare % len(levels)],
                         "duration_s": 9000.0,
                         "depth": 0.5, "net": 0.05, "noise": 1e-5})
    return pd.DataFrame(rows)


def test_the_joint_order_reads_back_the_scheme_it_is_built_from():
    """
    The whole point of the joint test: +1 for an adduct, and not +1 otherwise.

    A test that only checked the adduct case would pass on a regression that
    always returned +1, which is exactly the failure mode worth guarding --
    the archive's answer is that the constraint is VIOLATED, so the machinery
    has to be able to see it satisfied too.
    """
    print("\nthe joint peroxide order")
    for constant in (0.005, 0.03, 0.2):
        adduct = joint_peroxide_order(_peroxide_frame(constant, "adduct"))
        check(f"an adduct with K = {constant} /mM gives the required "
              f"{PERHYDRATE_ORDER_GAP:+.0f}",
              abs(adduct["slope"] - PERHYDRATE_ORDER_GAP)
              < 3 * adduct["stderr"] + 0.02,
              f"{adduct['slope']:+.3f} +- {adduct['stderr']:.3f}")
    # The identity is EXACT and carries no curvature: log(v/tau) is
    # log K + log h - log(k_r) whatever K is, so a global log-log slope over a
    # wide ladder has to return +1 and not an average of local slopes.
    exact = joint_peroxide_order(_peroxide_frame(0.2, "adduct", runs=12))
    check("and the identity carries no curvature, so a wide ladder is fine",
          abs(exact["slope"] - PERHYDRATE_ORDER_GAP) < 0.03,
          f"{exact['slope']:+.3f} over a 64-fold ladder")
    # But the floor breaks it, and in the direction of this module's own
    # conclusion, so the bias has to be visible rather than assumed away.
    floored = joint_peroxide_order(
        _peroxide_frame(0.2, "adduct", runs=12, baseline=1200.0))
    check("an induction driven below the floor biases the joint order DOWN",
          floored["slope"] < exact["slope"] - 0.05,
          f"{floored['slope']:+.3f} against {exact['slope']:+.3f}")

    trapped = joint_peroxide_order(_peroxide_frame(0.03, "trap"))
    check("a trap does not, and falls short rather than over",
          trapped["slope"] < PERHYDRATE_ORDER_GAP - 3 * trapped["stderr"],
          f"{trapped['slope']:+.3f} +- {trapped['stderr']:.3f}")
    check("and the shortfall is reported in standard errors",
          trapped["sigma"] > 3, f"{trapped['sigma']:.1f}")


def _planted_buffer_ladder(constant, rule="pre-equilibrium",
                           levels=((3.125, 6.25, 12.5, 25.0),
                                   (50.0, 100.0, 150.0, 200.0)),
                           runs=(34, 32), day=(1.0, 0.55), baseline=4000.0):
    """
    Two runs on one buffer ladder, built from the scheme the joint test asserts.

    `pre-equilibrium` puts the buffer in the activating equilibrium, so
    `v ~ K b/(1 + K b)` and `1/tau ~ 1 + K b` and the joint order is +1 by
    construction. `spectator` leaves the induction alone and keeps the same
    rate, so the joint order collapses to the rate's own order -- the negative
    control, without which a regression stuck at +1 would pass.

    The two runs get DIFFERENT levels, a 1.8-fold day step in the rate, because
    that is what the archive's pair does and what the free level per run is for.

    `baseline` is `tau` at zero buffer, and it has to be large enough that no
    planted time lands on `INDUCTION_FLOOR` -- the clip breaks the identity and
    biases the answer DOWN. At 400 s it did: K = 0.1 /mM over a ladder reaching
    200 mM puts tau at 19 s and the joint order at +0.743 rather than +1. That
    is a property of the floor and not of the regression, and the test below
    now asserts it deliberately instead of tripping over it.
    """
    rows = []
    for experiment, rungs, scale in zip(runs, levels, day):
        for index, buffer_mM in enumerate(rungs):
            bound = constant * buffer_mM / (1.0 + constant * buffer_mM)
            rate = scale * 1e-4 * bound
            if rule == "pre-equilibrium":
                tau = baseline / (1.0 + constant * buffer_mM)
            else:
                tau = baseline
            rows.append({"experiment": experiment, "sample": index + 1,
                         "buf": buffer_mM, "t_ind": tau, "v_peak": rate,
                         "live": True, "net": 0.5, "noise": 1e-4})
    return pd.DataFrame(rows)


def test_the_joint_order_takes_a_fitted_clock_and_a_gate():
    """
    The +1 does not belong to the landmark, and the block that can't use one.

    `t_ind` is a ROLLING WINDOW, so it is not comparable between runs of
    different length, and on the two-axis block `signal_control` fails at
    +0.619 +- 0.228 -- the landmark there is partly measuring the
    spectrophotometer. A time constant from the PROGRESS FIT carries no window
    and escapes both objections, and the identity holds for it just the same,
    so `joint_order` has to take one.

    With it comes the gate. An unresolved time constant is whatever the
    optimiser stopped at, and there is no floor that makes it honest -- so the
    flag has to drop the row, and this checks that it does by planting rows
    that would wreck the answer if it did not.
    """
    print("\nthe joint order on a fitted clock")
    frame = _peroxide_frame(0.2, "adduct", runs=12)
    fitted = joint_order(frame, axis="h2o2", rate="vmax",
                         timescale="tau_slow", covariates=("s0",),
                         floor=None, gate="tau_slow_resolved")
    check("a fitted clock reads back the planted adduct",
          abs(fitted["slope"] - PERHYDRATE_ORDER_GAP) < 0.03,
          f"{fitted['slope']:+.3f} +- {fitted['stderr']:.3f}")
    landmark = joint_peroxide_order(frame)
    check("and it agrees with the landmark where both are available",
          abs(fitted["slope"] - landmark["slope"]) < 0.03,
          f"{fitted['slope']:+.3f} against {landmark['slope']:+.3f}")

    spoilt = _peroxide_frame(0.2, "adduct", runs=12, unresolved=2)
    gated = joint_order(spoilt, axis="h2o2", rate="vmax",
                        timescale="tau_slow", covariates=("s0",),
                        floor=None, gate="tau_slow_resolved")
    ungated = joint_order(spoilt, axis="h2o2", rate="vmax",
                          timescale="tau_slow", covariates=("s0",),
                          floor=None)
    check("the gate drops the unresolved rows",
          gated["points"] == fitted["points"]
          and ungated["points"] > gated["points"],
          f"{gated['points']} gated against {ungated['points']} ungated")
    check("and it has to, because they would carry the answer",
          abs(fitted["slope"] - gated["slope"]) < 1e-9
          and abs(ungated["slope"] - PERHYDRATE_ORDER_GAP) > 0.1,
          f"gated {gated['slope']:+.3f}, ungated {ungated['slope']:+.3f}")

    # A floor is for a landmark, not for a fit: it would put the unresolved
    # constants ON the floor and call them the fastest curves in the block.
    floored = joint_order(spoilt, axis="h2o2", rate="vmax",
                          timescale="tau_slow", covariates=("s0",))
    check("a floor is no substitute for the gate",
          abs(floored["slope"] - PERHYDRATE_ORDER_GAP)
          > abs(gated["slope"] - PERHYDRATE_ORDER_GAP),
          f"floored {floored['slope']:+.3f} against gated {gated['slope']:+.3f}")


def test_the_joint_order_control_axis_reads_the_substrate_order():
    """
    THE AXIS IS THE CONTROL, and a test that only ever confirms is not a test.

    The +1 belongs to the species that draws the catalyst into its active form.
    Ask it of an axis whose species does not -- the alcohol -- and the answer
    is that axis's own order on the rate, because the clock does not carry it.
    So the control is not "it should be zero"; it is "it should be the
    substrate order", and missing +1 is the consequence.
    """
    print("\nthe joint order's control axis")
    for planted in (0.0, 0.6, 1.4):
        frame = _peroxide_frame(0.2, "adduct", runs=10,
                                substrate=(2.0, 4.0, 8.0, 16.0),
                                substrate_order=planted)
        control = joint_order(frame, axis="s0", rate="vmax",
                              timescale="tau_slow", covariates=("h2o2",),
                              floor=None, gate="tau_slow_resolved")
        check(f"the substrate axis returns its planted order {planted:+.1f}",
              abs(control["slope"] - planted) < 3 * control["stderr"] + 0.02,
              f"{control['slope']:+.3f} +- {control['stderr']:.3f}")
        peroxide = joint_order(frame, axis="h2o2", rate="vmax",
                               timescale="tau_slow", covariates=("s0",),
                               floor=None, gate="tau_slow_resolved")
        check("while the peroxide axis still returns the adduct's +1",
              abs(peroxide["slope"] - PERHYDRATE_ORDER_GAP) < 0.03,
              f"{peroxide['slope']:+.3f}")
    # And with a substrate order well away from 1 the control MISSES the +1,
    # which is the whole use of it.
    frame = _peroxide_frame(0.2, "adduct", runs=10,
                            substrate=(2.0, 4.0, 8.0, 16.0),
                            substrate_order=0.0)
    control = joint_order(frame, axis="s0", rate="vmax", timescale="tau_slow",
                          covariates=("h2o2",), floor=None,
                          gate="tau_slow_resolved")
    check("so a control axis rejects the constraint it was never subject to",
          control["sigma"] > 3, f"{control['sigma']:.1f} sigma from +1")


def test_the_two_named_orders_are_the_general_one():
    """
    `joint_peroxide_order` and `joint_buffer_order` were two copies of one
    regression, differing in their filters and their minimum. They are now one
    function with the species filled in, and this is what stops a third copy:
    each wrapper has to equal the general call it claims to be.
    """
    print("\nthe named orders are the general one")
    frame = _peroxide_frame(0.03, "adduct", runs=6)
    named = joint_peroxide_order(frame)
    general = joint_order(frame, axis="h2o2", rate="v_peak",
                          timescale="t_ind", covariates=("s0",),
                          floor=INDUCTION_FLOOR, minimum=8)
    check("the peroxide order IS the general order on h2o2",
          all(abs(named[k] - general[k]) < 1e-12
              for k in ("slope", "stderr", "sigma")),
          f"{named['slope']:+.6f} against {general['slope']:+.6f}")

    ladder = _planted_buffer_ladder(0.02)
    named = joint_buffer_order(ladder)
    general = joint_order(ladder, axis="buf", rate="v_peak",
                          timescale="t_ind", floor=INDUCTION_FLOOR,
                          minimum=4, live_only=False)
    check("and the buffer order IS the general order on buf",
          all(abs(named[k] - general[k]) < 1e-12
              for k in ("slope", "stderr", "sigma")),
          f"{named['slope']:+.6f} against {general['slope']:+.6f}")


def test_the_joint_buffer_order_reads_back_its_own_scheme():
    """
    +1 when the buffer is in the activating equilibrium, and not otherwise.

    The archive's answer here is that the constraint is SATISFIED -- which is
    the opposite of the peroxide axis and so the more dangerous way round: a
    regression that returned +1 on anything would confirm the conclusion for
    the wrong reason. Hence the spectator control.
    """
    print("\nthe joint buffer order")
    for constant in (0.004, 0.02, 0.1):
        got = joint_buffer_order(_planted_buffer_ladder(constant))
        check(f"a buffer pre-equilibrium at K = {constant} /mM gives "
              f"{PERHYDRATE_ORDER_GAP:+.0f}",
              abs(got["slope"] - PERHYDRATE_ORDER_GAP) < 0.02,
              f"{got['slope']:+.3f}")
    spectator = joint_buffer_order(
        _planted_buffer_ladder(0.02, rule="spectator"))
    check("a buffer that does not touch the induction falls short of +1",
          spectator["slope"] < PERHYDRATE_ORDER_GAP - 0.1,
          f"{spectator['slope']:+.3f}")
    check("and the shortfall is what the rate order alone would give",
          spectator["slope"] > 0.0, f"{spectator['slope']:+.3f}")
    # The day step is what the free level per run absorbs. Double it and the
    # answer must not move; pool the runs instead and it does.
    stepped = joint_buffer_order(
        _planted_buffer_ladder(0.02, day=(1.0, 0.2)))
    check("a larger between-run step does not move it, because of the levels",
          abs(stepped["slope"] - PERHYDRATE_ORDER_GAP) < 0.02,
          f"{stepped['slope']:+.3f}")

    # THE FLOOR BREAKS THE IDENTITY, and downward -- towards the shortfall the
    # peroxide axis reports -- so it is asserted here rather than avoided.
    # This test did not run until 2026-09-03 and tripped over it at K = 0.1.
    clipped = joint_buffer_order(_planted_buffer_ladder(0.1, baseline=400.0))
    check("an induction driven below the floor biases the buffer order DOWN",
          clipped["slope"] < PERHYDRATE_ORDER_GAP - 0.1,
          f"{clipped['slope']:+.3f} with tau clipped at "
          f"{INDUCTION_FLOOR:.0f} s")
    clear = joint_buffer_order(_planted_buffer_ladder(0.1))
    check("and clear of the floor the same ladder reads back +1",
          abs(clear["slope"] - PERHYDRATE_ORDER_GAP) < 0.02,
          f"{clear['slope']:+.3f}")
    # And the archive's own ladder is clear of it: its shortest induction is
    # 120 s and the answer does not move for any floor up to that.
    check("the archive's buffer ladder is not floor-limited",
          all(abs(joint_buffer_order(_planted_buffer_ladder(0.02),
                                     floor=floor)["slope"]
                  - PERHYDRATE_ORDER_GAP) < 0.02
              for floor in (1.0, 30.0, 60.0, 120.0)))


def test_the_saturation_fit_recovers_a_planted_constant():
    """K profiled on the scheme's own form, with the exponent held at 1."""
    print("\nthe peroxide saturation fit")
    for constant in (0.01, 0.05):
        table = _peroxide_frame(constant, "adduct", runs=8)
        got = peroxide_saturation(table, parameter="v_peak")
        check(f"K = {constant} /mM is inside the profile interval",
              got["constant_low"] <= constant <= got["constant_high"],
              f"{got['constant_low']:.4f}-{got['constant_high']:.4f}, "
              f"best {got['constant']:.4f}")
    # A rate that really is first order must not be called saturating.
    straight = _peroxide_frame(1e-6, "adduct", runs=8)
    got = peroxide_saturation(straight, parameter="v_peak")
    check("a genuinely first-order rate reads back as first order",
          abs(got["order"] - 1.0) < 0.1, f"a = {got['order']:.3f}")
    check("and is not rejected as one",
          got["first_order_f"] < 4.0, f"F = {got['first_order_f']:.2f}")


def test_the_trap_constant_inverts_its_own_order():
    """trap_constant is the algebraic inverse of the order it is given."""
    print("\nthe trap constant")
    for constant, peroxide in ((0.01, 40.0), (0.05, 25.0), (0.2, 10.0)):
        order = constant * peroxide / (1.0 + constant * peroxide)
        got = trap_constant(order, 0.05, peroxide)
        check(f"K = {constant} /mM recovered from its own order "
              f"{order:.3f} at {peroxide:.0f} mM",
              abs(got["constant"] - constant) < 1e-9,
              f"{got['constant']:.6f}")
    check("an order at or above 1 has no finite K",
          not np.isfinite(trap_constant(1.0, 0.1, 40.0)["constant"]))
    check("and neither does one at or below 0",
          not np.isfinite(trap_constant(-0.2, 0.1, 40.0)["constant"]))
    # 29 /M is 0.029 /mM, and -RT ln(29) is -8.3 kJ/mol at 298 K.
    got = trap_constant(0.029 * 28.3 / (1.0 + 0.029 * 28.3), 0.05, 28.3)
    check("the free energy is quoted from a MOLAR constant",
          abs(got["free_energy_kJ"] + 8.3) < 0.1,
          f"{got['free_energy_kJ']:+.2f} kJ/mol from {got['molar']:.0f} /M")


def test_the_floor_does_not_carry_the_answer():
    """Zeros go in at a floor, and the conclusion must survive moving it."""
    print("\nthe floor")
    table = _induction_table()
    block = induction.induction_blocks(table)["4OMe catalysed"]
    slopes = [induction_drivers(block, floor=floor)["slope"]
              for floor in induction.FLOOR_SWEEP]
    check("every floor from 1 s to 300 s keeps the slope far from -1",
          all(slope - INDUCTION_PRODUCT_SLOPE > 0.6 for slope in slopes),
          " ".join(f"{s:+.3f}" for s in slopes))
    check("and none of them is more than 0.4 from the clock",
          all(abs(slope - INDUCTION_CLOCK_SLOPE) < 0.4 for slope in slopes),
          " ".join(f"{s:+.3f}" for s in slopes))
    check("the floor sits at the shortest induction the readings can resolve",
          INDUCTION_FLOOR >= float(block[block.t_ind > 0].t_ind.min()) - 1e-9,
          f"floor {INDUCTION_FLOOR} s, smallest measured "
          f"{float(block[block.t_ind > 0].t_ind.min()):.0f} s")


def test_orders_refuse_an_axis_the_block_does_not_move():
    """The guard in scope.orders, which this module is the reason for."""
    print("\nunidentified axes")
    table = _induction_table()
    fixed = table[table.experiment.isin(induction.PEROXIDE_LEVER)]
    check("exps 127-131 hold [S] at one value", fixed.s0.nunique() == 1,
          f"{sorted(fixed.s0.unique())}")
    got = scope.orders("v_peak", frame=fixed)
    check("so no substrate order is reported for them",
          not np.isfinite(got["order_s0"]), f"{got['order_s0']}")
    check("and the peroxide order still is",
          np.isfinite(got["order_h2o2"]), f"{got['order_h2o2']}")
    whole = scope.orders("vmax")
    check("the two-axis orders are untouched by the guard",
          abs(whole["order_h2o2"] - 0.7936) < 5e-4
          and abs(whole["order_s0"] - 0.0906) < 5e-4,
          f"{whole['order_s0']:.4f}, {whole['order_h2o2']:.4f}")


_CACHE = {}


def _induction_table():
    if "table" not in _CACHE:
        _CACHE["table"] = induction.induction_table(induction.WHOLE_ARCHIVE)
    return _CACHE["table"]


def test_the_sign_comes_off_the_fit_and_not_off_the_depth():
    """
    `progress_kind` has to distinguish shapes `depth` cannot.

    depth = 1 - start/peak is floored at zero, so a curve that is fastest in
    its first window and one that starts just below its peak read the same. The
    sign has to separate them, and it has to come from the fitted form rather
    than from the readings' extremes.
    """
    print("\nwhich way the curve points")
    from summary_kinetics import fit_progress
    times = np.linspace(0.0, 8000.0, 400)
    # With noise, because a NOISELESS single exponential is not a realistic
    # input and does not behave like one: its residual falls to numerical dust,
    # the F test always buys the second phase, and a pure lag comes back as
    # "two lags" with the one relaxation split across both terms.
    grain = np.random.default_rng(5).normal(0.0, 2e-4, len(times))
    lag = 1e-5 * (times - 900.0 * (1.0 - np.exp(-times / 900.0))) + grain
    burst = 1e-5 * (times + 900.0 * (1.0 - np.exp(-times / 900.0))) + grain
    check("a lag curve is called one of the lag-first kinds",
          fit_progress(times, lag).kind in induction.LAG_FIRST_KINDS,
          fit_progress(times, lag).kind)
    check("a burst curve is not",
          fit_progress(times, burst).kind not in induction.LAG_FIRST_KINDS,
          fit_progress(times, burst).kind)
    check("the lag's fast amplitude is positive and the burst's negative",
          fit_progress(times, lag).amplitudes[0] > 0
          > fit_progress(times, burst).amplitudes[0],
          f"{fit_progress(times, lag).amplitudes[0]:.4g} and "
          f"{fit_progress(times, burst).amplitudes[0]:.4g}")
    clean = 1e-5 * (times + 900.0 * (1.0 - np.exp(-times / 900.0)))
    check("and induction_point cannot tell a burst from a straight line",
          induction_point(FakeCurve(times, clean)).depth
          == induction_point(FakeCurve(times, 1e-5 * times)).depth == 0.0)

    table = _induction_table()
    check("every live curve in the archive carries a sign",
          bool(table[table.live].progress_kind.notna().all()))
    check("and lag_first is exactly the kinds that begin below their rate",
          set(table[table.lag_first].progress_kind.unique())
          == set(induction.LAG_FIRST_KINDS)
          - {"two lags"} | ({"two lags"} if
                            (table.progress_kind == "two lags").any() else set()),
          f"{sorted(table[table.lag_first].progress_kind.unique())}")


def test_the_L_has_to_be_split_before_it_is_read():
    """
    Pooling the two-axis L reports a peroxide effect the peroxide arm denies.

    This is the failure the split exists to prevent, so it is checked rather
    than described: the pooled fit and the arm have to disagree, and the arm
    has to be the one that holds [S] fixed.
    """
    print("\nthe two arms of the L")
    table = _induction_table()
    scoped = induction.induction_blocks(table)["BnOH two-axis (135-151)"]
    arms = ladder_arms(scoped)
    check("the substrate arm holds [H2O2] fixed inside every run",
          bool((arms["substrate arm"].groupby("experiment").h2o2.nunique()
                == 1).all()))
    check("the peroxide arm holds [S] fixed inside every run",
          bool((arms["peroxide arm"].groupby("experiment").s0.nunique()
                == 1).all()))
    pooled = sign_drivers(scoped, axis="h2o2", control=True)
    arm = sign_drivers(arms["peroxide arm"], axis="h2o2", control=True)
    check("the pooled peroxide effect is larger than the arm's",
          abs(pooled["h2o2"]) > abs(arm["h2o2"]),
          f"pooled {pooled['h2o2']:+.3f} against arm {arm['h2o2']:+.3f}")
    check("and inside the arm it is consistent with nothing",
          abs(arm["h2o2"]) < 2 * arm["h2o2_stderr"],
          f"{arm['h2o2']:+.3f} +- {arm['h2o2_stderr']:.3f}")


def test_the_blocks_differ_in_their_buffer_collinearity():
    """The caveat section 3 and section 6 both rest on, as a number."""
    print("\nthe substrate/buffer collinearity")
    table = _induction_table()
    blocks = induction.induction_blocks(table)
    four = composition_collinearity(blocks["4OMe catalysed"])
    scoped = composition_collinearity(blocks["BnOH two-axis (135-151)"])
    check("[S] and [buf] move together in every 4OMe run",
          four["median"] < -0.9 and four["constant_buffer"] == 0,
          f"median {four['median']:+.2f}, {four['constant_buffer']} constant "
          f"of {four['runs']}")
    check("and in no two-axis run",
          scoped["median"] == 0.0
          and scoped["constant_buffer"] == scoped["runs"],
          f"median {scoped['median']:+.2f}, {scoped['constant_buffer']} "
          f"constant of {scoped['runs']}")


def test_the_window_free_lag_recovers_a_planted_relaxation():
    """`lag_profile` has to read back the tau and depth it was given."""
    print("\nthe window-free lag, against a planting")
    for tau in (200.0, 800.0, 3000.0):
        curve = lag_curve(tau=tau, span=20000.0, points=400)
        fit = fit_progress(curve.times, curve.absorbance)
        depth, half, peak, start = fit.lag_profile(float(curve.times[-1]))
        check(f"tau {tau:.0f} s reads back as tau.ln2 = {tau * np.log(2):.0f} s",
              abs(half - tau * np.log(2)) < 0.06 * tau * np.log(2),
              f"{half:.1f}")
        check(f"tau {tau:.0f} s: the depth of a curve starting at rest is 1",
              abs(depth - 1.0) < 0.02, f"{depth:.4f}")

    # A curve that begins at its fastest has NO lag, and the statistic has to
    # say zero rather than something small and positive.
    times = np.linspace(0.0, 8000.0, 400)
    burst = 1e-5 * times + 5e-3 * (1.0 - np.exp(-times / 500.0))
    fit = fit_progress(times, burst)
    depth, half, _, _ = fit.lag_profile(float(times[-1]))
    check("a burst-first curve returns depth 0 and clock 0",
          depth == 0.0 and half == 0.0, f"depth {depth}, half {half}")

    # And it is read INSIDE the window. A run truncated before the rise is over
    # must under-read, not extrapolate: that under-reading is the whole reason
    # `lag_window_frame` exists, so it has to be real. A PARTIAL lag is needed
    # to see it -- a curve that starts at rest has depth 1 through any window,
    # because the depth is measured against the initial rate and that is zero.
    times = np.linspace(0.0, 30000.0, 600)
    rate, tau = 1e-5, 3000.0
    partial = rate * times - 0.5 * rate * tau * (1.0 - np.exp(-times / tau))
    fit = fit_progress(times, partial)
    full, _, _, _ = fit.lag_profile(30000.0)
    short, _, _, _ = fit.lag_profile(2000.0)
    check("the planted half-depth reads back over the whole run",
          abs(full - 0.5) < 0.02, f"{full:.3f}")
    check("and a window shorter than the rise under-reads it",
          short < full - 0.10, f"{short:.3f} against {full:.3f}")


def test_the_window_free_lag_is_not_windowed():
    """
    The property the whole section rests on: no run-fraction anywhere.

    `induction_point` reads through a window a tenth of the run, so the SAME
    curve recorded for twice as long gives a different landmark. `lag_profile`
    must not: read the same relaxation over two spans and the clock is the
    same number, because it comes off the form and not off a window.
    """
    print("\nthe window-free lag carries no window")
    clocks = []
    for span, points in ((24000.0, 800), (12000.0, 400)):
        curve = lag_curve(tau=600.0, span=span, points=points)
        fit = fit_progress(curve.times, curve.absorbance)
        clocks.append(fit.lag_profile(float(curve.times[-1]))[1])
    check("the fitted clock is the same on a run twice as long",
          abs(clocks[0] - clocks[1]) < 0.05 * clocks[0],
          f"{clocks[0]:.1f} against {clocks[1]:.1f}")

    # WHERE IT IS NOT INVARIANT, AND WHY THAT IS THE DESIGN. The profile is
    # read inside the window, so a run that ends before the rise is over
    # reports a shorter clock -- honestly, since it has not seen the rest. It
    # is the SCHEDULE the statistic still carries, not a window, and
    # `lag_window_frame` is the answer to it.
    truncated, complete = [], []
    for span, points in ((8000.0, 300), (40000.0, 1500)):
        curve = lag_curve(tau=3000.0, span=span, points=points)
        fit = fit_progress(curve.times, curve.absorbance)
        (truncated if span < 4 * 3000.0 else complete).append(
            fit.lag_profile(float(curve.times[-1]))[1])
    check("a run only 2.7 time constants long under-reads its own clock",
          truncated[0] < complete[0] - 0.05 * complete[0],
          f"{truncated[0]:.1f} against {complete[0]:.1f}")

    # And the landmark's window dependence is not synthetic: on a noiseless
    # planting the two agree, and it is the ARCHIVE where the landmark tracks
    # run length at +0.437 +/- 0.181. `test_regressions` holds that number.
    landmark = induction_point(lag_curve(tau=600.0, span=24000.0, points=800))
    check("the landmark agrees with the fitted clock on a clean planting",
          abs(landmark.t_ind - clocks[0]) < 0.10 * clocks[0],
          f"{landmark.t_ind:.1f} against {clocks[0]:.1f}")


def test_a_between_run_ladder_refuses_offsets_that_absorb_it():
    """
    pH is one value per run, so a per-experiment offset IS the pH axis.

    This is `scope._moves` one level out, and it was live for an afternoon:
    `lag_ph_ladders` first ran with one offset per experiment and reported
    +0.549 +/- 0.014 on the boric ladder -- forty sigma of pure collinearity.
    """
    print("\na between-run ladder and the offsets that would absorb it")
    raised = False
    try:
        lag_ladder(scope.PH_LADDER_BORIC, axis="pH", offsets="experiment")
    except ValueError as problem:
        raised = "absorbed" in str(problem)
    check("one offset per experiment on a pH ladder raises", raised)

    got = lag_ladder(scope.PH_LADDER_BORIC, axis="pH")
    check("with a single level it returns a coefficient",
          np.isfinite(got["lag_half_s"]["slope"]), f"{got.get('lag_half_s')}")
    check("and it reports both collinearities the reader needs",
          abs(got["schedule_collinearity"]) > 0.5
          and np.isfinite(got["signal_collinearity"]),
          f"{got['schedule_collinearity']:+.2f}, "
          f"{got['signal_collinearity']:+.2f}")

    # The two-axis ladders repeat seven compositions across runs, so a CUVETTE
    # offset is matched and must be allowed.
    matched = lag_ladder(scope.PH_LADDER_TWO_AXIS_HIGH, axis="pH",
                         offsets="sample")
    check("a cuvette offset is not collinear with pH and is allowed",
          np.isfinite(matched["lag_half_s"]["slope"]),
          f"{matched.get('lag_half_s')}")


def test_the_common_window_removes_the_schedule():
    """
    Every run in a `lag_window_frame` is fitted over the same span.

    That is the property, and it is checkable directly: the boric ladder's runs
    span 1260 to 17940 s and their fitted windows must all be the shortest.
    """
    print("\nthe common window")
    table = lag_window_frame(scope.PH_LADDER_BORIC)
    spans = table.window_s.to_numpy(dtype=float)
    durations = table.duration_s.to_numpy(dtype=float)
    check("the runs really do differ in length",
          durations.max() / durations.min() > 10,
          f"{durations.min():.0f} to {durations.max():.0f} s")
    check("and every curve is fitted over the same window",
          float(spans.max() - spans.min()) < 120.0,
          f"{spans.min():.0f} to {spans.max():.0f} s")
    check("which is the shortest run",
          abs(float(spans.max()) - common_window(scope.PH_LADDER_BORIC)) < 120.0,
          f"{spans.max():.0f} against "
          f"{common_window(scope.PH_LADDER_BORIC):.0f}")

    # And it has to change the answer, or there was nothing to remove.
    whole = lag_ladder(scope.PH_LADDER_BORIC, axis="pH",
                       window=common_window(scope.PH_LADDER_BORIC) * 20)
    shared = lag_ladder(scope.PH_LADDER_BORIC, axis="pH")
    check("the schedule was carrying the boric ladder's pH coefficient",
          abs(whole["lag_half_s"]["slope"] - shared["lag_half_s"]["slope"]) > 0.5,
          f"{whole['lag_half_s']['slope']:+.3f} whole runs against "
          f"{shared['lag_half_s']['slope']:+.3f} shared")


def test_the_barrier_survives_the_schedule_and_the_naive_control_does_not():
    """
    The one between-run result this section rescues rather than withdraws.

    On whole runs the induction's barrier is 83.7 +/- 8.9 kJ/mol; putting
    log(run length) into that regression drops it to about 55 on a six-point
    collinearity of +0.66. At a window all six runs share it is 73-84, which is
    the number to quote -- and both bounds have to be checked, or the sweep is
    decoration.
    """
    print("\nthe induction's barrier, against the schedule")
    sweep = lag_arrhenius_sweep()
    full = [row for row in sweep if row["fraction"] == 1.0][0]
    check("at the common window all six temperatures contribute",
          full["temperatures"] == 6, f"{full['temperatures']}")
    check("and the barrier is covalent-sized: 60-100 kJ/mol",
          60.0 < full["activation_kj"] < 100.0,
          f"{full['activation_kj']:.1f} +/- {full['stderr_kj']:.1f}")
    wide = [row["activation_kj"] for row in sweep if row["fraction"] >= 0.75]
    check("it is stable from three quarters of the window upwards",
          max(wide) - min(wide) < 15.0, f"{[round(v, 1) for v in wide]}")
    narrow = [row for row in sweep if row["fraction"] <= 0.5]
    check("and it collapses below it, because the cold clocks are censored",
          all(row["activation_kj"] < 30.0 for row in narrow),
          f"{[round(row['activation_kj'], 1) for row in narrow]}")
    check("which the temperature count shows rather than hides",
          all(row["temperatures"] < 6 for row in narrow),
          f"{[row['temperatures'] for row in narrow]}")


def test_the_replicate_floor_is_the_bar_between_runs():
    """Four runs that differ in nothing still differ; by how much."""
    print("\nthe replicate floor")
    floor = replicate_floor()
    check("it is four runs of one composition", floor["runs"] == 4,
          f"{floor['runs']}")
    check("the clock reproduces between runs to better than 1.5x",
          1.0 < floor["clock_ratio"] < 1.5, f"{floor['clock_ratio']:.2f}x")
    check("the DEPTH does not -- its replicate spread is over half a unit",
          floor["depth_spread"] > 0.5, f"{floor['depth_spread']:.3f}")


def test_lag_orders_agree_with_the_package_order_machinery():
    """
    `lag_orders` must BE `scope.orders`, not a second copy of it.

    The lag statistic already had two definitions once and they disagreed on
    96 of 402 curves. This checks the wrapper adds the signal pairing and
    changes nothing else.
    """
    print("\nlag_orders is scope.orders")
    block = scope.frame(scope.TWO_AXIS_BLOCK)
    got = lag_orders(block)
    direct = scope.orders("lag_half_s", frame=block, floor=INDUCTION_FLOOR)
    check("the bare substrate order is scope.orders' own",
          abs(got["lag_half_s"]["s0"]["order"] - direct["order_s0"]) < 1e-12,
          f"{got['lag_half_s']['s0']['order']} against {direct['order_s0']}")
    check("and so is the peroxide order",
          abs(got["lag_half_s"]["h2o2"]["order"] - direct["order_h2o2"]) < 1e-12)

    # The joint fit is not the single-axis fit, and the difference is the point.
    single = scope.orders("lag_half_s", frame=block, floor=INDUCTION_FLOOR,
                          terms=("s0",))
    check("fitting one axis of an L is a different number",
          abs(single["order_s0"] - direct["order_s0"]) > 0.1,
          f"{single['order_s0']:+.3f} against {direct['order_s0']:+.3f}")

    # And the covariate has to be identified the way an order is.
    held = scope.orders("lag_half_s", frame=block, floor=INDUCTION_FLOOR,
                        covariates=("e0",))
    check("a covariate constant inside every run is refused, not fitted",
          "held_e0" not in held, f"{sorted(held)}")


def test_the_signal_collinearity_says_which_axis_can_be_asked():
    """
    The two-axis block fails its signal control and can still be asked about
    substrate, because in that block substrate buys no signal. That is the
    structural argument the whole section rests on, so it is measured.
    """
    print("\nwhich axis a block that fails the signal control can carry")
    block = scope.frame(scope.TWO_AXIS_BLOCK)
    got = lag_orders(block)
    check("peroxide IS the signal here",
          got["signal_collinearity_h2o2"] > 0.4,
          f"{got['signal_collinearity_h2o2']:+.3f}")
    check("substrate is not",
          abs(got["signal_collinearity_s0"]) < 0.2,
          f"{got['signal_collinearity_s0']:+.3f}")
    peroxide = got["lag_half_s"]["h2o2"]
    substrate = got["lag_half_s"]["s0"]
    check("so holding the signal moves the peroxide order a long way",
          abs(peroxide["order"] - peroxide["controlled"]) > 0.2,
          f"{peroxide['order']:+.3f} -> {peroxide['controlled']:+.3f}")
    check("and the substrate order stays negative either way",
          substrate["order"] < 0 and substrate["controlled"] < 0,
          f"{substrate['order']:+.3f} -> {substrate['controlled']:+.3f}")

    control = lag_signal_control(block)
    check("the block fails the signal control outright",
          control["lag_half_s"]["slope"]
          > 3 * control["lag_half_s"]["stderr"],
          f"{control['lag_half_s']['slope']:+.3f} +/- "
          f"{control['lag_half_s']['stderr']:.3f}")

    catalysed = induction_blocks(scope.frame(scope.archive()))["4OMe catalysed"]
    passes = lag_signal_control(catalysed)
    check("the catalysed 4OMe block passes it",
          abs(passes["lag_half_s"]["slope"])
          < 2 * passes["lag_half_s"]["stderr"],
          f"{passes['lag_half_s']['slope']:+.3f} +/- "
          f"{passes['lag_half_s']['stderr']:.3f}")


def test_the_pooled_ladder_notices_a_disagreement():
    """A pooled coefficient that hides a disagreement is worse than none."""
    print("\npooling four ladders")
    agreeing = [{"x": {"controlled": 0.30, "controlled_stderr": 0.10}},
                {"x": {"controlled": 0.34, "controlled_stderr": 0.10}}]
    disagreeing = [{"x": {"controlled": -0.60, "controlled_stderr": 0.10}},
                   {"x": {"controlled": +0.60, "controlled_stderr": 0.10}}]
    good = pooled_ladder(agreeing, "x")
    bad = pooled_ladder(disagreeing, "x")
    check("two ladders that agree give a small chi2", good["chi2"] < 1.0,
          f"{good['chi2']:.2f}")
    check("two that do not give a large one", bad["chi2"] > 50.0,
          f"{bad['chi2']:.2f}")
    check("and the pooled value is the inverse-variance mean",
          abs(good["pooled"] - 0.32) < 1e-9, f"{good['pooled']}")

    real = pooled_ladder(lag_ph_ladders(), "lag_half_s")
    check("the four real pH ladders agree with each other",
          real["chi2"] < 7.81, f"chi2 {real['chi2']:.2f} on {real['dof']}")


def test_the_saturation_fraction_bounds_both_schemes():
    """A coefficient outside (0,1) falsifies the trap; inside, it sizes it."""
    print("\nthe bounded schemes")
    trap = saturation_fraction(1.0, 0.2)
    check("a pH coefficient is divided by ln 10 first",
          abs(trap["fraction"] - 1.0 / np.log(10.0)) < 1e-12,
          f"{trap['fraction']:.4f}")
    check("and +1.0 per pH unit is a half-saturated trap",
          trap["inside_trap"] and not trap["inside_activating"],
          f"{trap}")
    outside = saturation_fraction(4.0, 0.2)
    check("a coefficient past ln 10 per pH unit falsifies both schemes",
          not outside["inside_trap"] and not outside["inside_activating"],
          f"{outside['fraction']:.3f}")
    activating = saturation_fraction(-0.5, 0.1, per_ph=False)
    check("a negative order is the activating branch",
          activating["inside_activating"] and not activating["inside_trap"],
          f"{activating['fraction']:.3f}")


def test_the_identifiability_table_names_the_between_run_axes():
    """Five of the seven variables are one value per run. Say so."""
    print("\nwhat the archive moves, and where")
    table = lag_identifiability()
    between = sorted(table.index[table.runs_moving_it == 0])
    check("temperature, pH, the buffer salt and the substrate move only "
          "between runs",
          between == ["buffer", "pH", "substrate", "temperature"],
          f"{between}")
    check("[S] moves inside most runs",
          table.loc["s0", "runs_moving_it"] > 70,
          f"{table.loc['s0', 'runs_moving_it']}")
    check("[buf] inside about half",
          40 < table.loc["buf", "runs_moving_it"] < 70,
          f"{table.loc['buf', 'runs_moving_it']}")
    check("[H2O2] inside a fifth",
          10 < table.loc["h2o2", "runs_moving_it"] < 30,
          f"{table.loc['h2o2', 'runs_moving_it']}")
    check("a pH span is reported in units and not as a ratio",
          table.loc["pH", "unit"] == "units"
          and table.loc["pH", "widest_within_run"] == 0.0,
          f"{table.loc['pH', 'unit']}, {table.loc['pH', 'widest_within_run']}")


def test_the_product_at_the_landmark_tells_a_clock_from_a_threshold():
    """
    Plant both hypotheses and check the statistic reads each one back.

    A CLOCK: every curve relaxes with the same tau and only the rate differs,
    so the landmark falls at the same time and the product made by it is
    proportional to the rate -- slope +1.

    A THRESHOLD: tau is set to C/rate, so the landmark moves inversely with the
    rate and the product made by it is the same on every curve -- slope 0.

    Both plantings use the same curve form and the same landmark, so a
    statistic that could not separate them would return the same number twice.
    """
    print("\nthe product at the landmark, against both plantings")
    # Twelve rates over eight-fold, which is what `product_at_landmark`
    # needs to fit at all -- it refuses under ten curves, because a slope
    # from four points with a run offset is one degree of freedom.
    rates = tuple(2e-6 * 8.0 ** (index / 11.0) for index in range(12))

    def planted(tau_of):
        rows = []
        for index, rate in enumerate(rates, start=1):
            tau = tau_of(rate)
            times = np.linspace(0.0, 40000.0, 800)
            values = rate * (times - tau * (1.0 - np.exp(-times / tau)))
            curve = FakeCurve(times, values, sample=index)
            found = induction_point(curve)
            rows.append({"experiment": 1, "sample": index, "live": True,
                         "t_ind": found.t_ind, "made": found.made,
                         "v_peak": found.peak_rate, "epsilon": curve.epsilon,
                         "s0": 5.0})
        return pd.DataFrame(rows)

    clock = product_at_landmark(planted(lambda rate: 3000.0),
                                minimum_per_run=4, calibrate=False)
    threshold = product_at_landmark(planted(lambda rate: 6e-3 / rate),
                                    minimum_per_run=4, calibrate=False)
    check("a planted clock reads back as +1, exactly",
          abs(clock["slope"] - PRODUCT_CLOCK_SLOPE) < 0.02,
          f"{clock['slope']:+.4f}")
    check("a planted product threshold does NOT read back as 0",
          threshold["slope"] > PRODUCT_THRESHOLD_SLOPE + 0.15,
          f"{threshold['slope']:+.4f}")
    check("and the two are still far apart",
          clock["slope"] - threshold["slope"] > 0.5,
          f"{clock['slope']:+.3f} against {threshold['slope']:+.3f}")

    # The spread comparison points the same way on the same plantings and uses
    # no regression at all -- and carries the same bias, which is why it is
    # quoted as a direction and not as a discriminator on its own.
    check("the clock planting holds the TIME constant",
          clock["spread_time"] < 0.05 * clock["spread_product"],
          f"{clock['spread_time']:.4f} against {clock['spread_product']:.4f}")
    check("the threshold planting holds the PRODUCT more nearly constant",
          threshold["spread_product"] < 0.6 * threshold["spread_time"],
          f"{threshold['spread_product']:.4f} against "
          f"{threshold['spread_time']:.4f}")


def test_the_product_test_is_calibrated_against_its_own_bias():
    """
    A planted threshold reads +0.4 through this landmark, and the sigma has
    to be measured from there. Quoting the nominal 0.0 overstates the
    exclusion by more than a factor of two.
    """
    print("\nthe product test, calibrated")
    recovery = product_recovery()
    low, high = recovery["threshold_range"]
    check("a planted threshold reads back well above zero, at every geometry",
          low > 0.25 and high < 0.55, f"{low:.3f} to {high:.3f}")
    check("and a planted clock reads back at +1 at every geometry",
          abs(recovery["clock_range"][0] - 1.0) < 0.01
          and abs(recovery["clock_range"][1] - 1.0) < 0.01,
          f"{recovery['clock_range']}")

    table = induction_table(WHOLE_ARCHIVE)
    got = product_at_landmark(induction_blocks(table)["4OMe catalysed"])
    check("the archive's answer is quoted against the recovered value",
          abs(got["threshold_reads"] - high) < 1e-9,
          f"{got['threshold_reads']:.3f} against {high:.3f}")
    naive = abs(got["slope"] - PRODUCT_THRESHOLD_SLOPE) / got["stderr"]
    check("which nearly halves the exclusion the nominal value would give",
          naive > 1.8 * got["threshold_sigma"],
          f"{naive:.1f} sigma nominal against "
          f"{got['threshold_sigma']:.1f} calibrated")
    check("and the exclusion survives the correction",
          got["threshold_sigma"] > 3.0, f"{got['threshold_sigma']:.1f}")


def test_the_product_test_agrees_with_the_driver_regression():
    """
    It is `induction_drivers` in other units, and that is checkable.

    Over the induction the product is about the rate times the time, so the
    two slopes must differ by 1. If they ever stop doing so, one of them has
    started measuring something else.
    """
    print("\nthe product test is the driver regression restated")
    table = induction_table(WHOLE_ARCHIVE)
    block = induction_blocks(table)["4OMe catalysed"]
    product = product_at_landmark(block)
    driver = induction_drivers(block, rate="v_peak")
    check("the two slopes differ by 1, as the algebra requires",
          abs((product["slope"] - driver["slope"]) - 1.0) < 0.10,
          f"{product['slope']:+.3f} against {driver['slope']:+.3f}")
    check("so it is quoted as a restatement and not as a second witness",
          "not independent" in product_at_landmark.__doc__.lower()
          or "NOT INDEPENDENT" in product_at_landmark.__doc__)

    check("it drops the curves with no landmark, and says how many",
          product["dropped"] > 0
          and product["curves"] + product["dropped"] == int(block.live.sum()),
          f"{product['curves']} kept, {product['dropped']} dropped")
    check("the archive's own answer is the clock, not the threshold",
          product["clock_sigma"] < 1.0 and product["threshold_sigma"] > 3.0,
          f"{product['clock_sigma']:.1f} from a clock, "
          f"{product['threshold_sigma']:.1f} from a threshold")
    check("and the time is the more nearly conserved quantity",
          product["spread_time"] < product["spread_product"],
          f"sd(log t) {product['spread_time']:.3f} against "
          f"sd(log P) {product['spread_product']:.3f}")


def test_regressions():
    """The numbers induction/ANALYSIS.md quotes."""
    print("\nthe published numbers")
    table = _induction_table()
    named = induction.induction_blocks(table)
    summary = induction.channel_summary(table)

    check("49 enzyme-free 4OMe curves, and not one accelerates",
          summary["enzyme_free"]["curves"] == 49
          and summary["enzyme_free"]["accelerates"] == 0,
          f"{summary['enzyme_free']}")
    check("none of them is deep either",
          summary["enzyme_free"]["deep"] == 0)
    check("the largest acceleration z in the enzyme-free block is 2.25",
          abs(summary["enzyme_free"]["max_accel_z"] - 2.25) < 0.01,
          f"{summary['enzyme_free']['max_accel_z']:.2f}")
    check("151 catalysed 4OMe curves, 83 of them deep",
          summary["catalysed"]["curves"] == 151
          and summary["catalysed"]["deep"] == 83,
          f"{summary['catalysed']}")
    check("the enzyme-free block includes a run of 17934 s",
          summary["enzyme_free"]["longest_s"] == 17934.0,
          f"{summary['enzyme_free']['longest_s']}")

    contrast = induction.channel_contrast(table)
    free40 = contrast[(contrast.temperature == 40.0)
                      & (contrast.channel == "enzyme-free")].iloc[0]
    check("37 matched enzyme-free curves at 40 C, none accelerating",
          free40.curves == 37 and free40.accelerates == 0,
          f"{free40.curves}, {free40.accelerates}")

    control = induction.schedule_control(table)
    check("exps 19 and 14 share a run length",
          control["same_span"] and control[19]["duration_s"] == 17934.0)
    check("and their induction time constants differ 6.9-fold",
          abs(control["tau_ratio"] - 6.9) < 0.1, f"{control['tau_ratio']:.2f}")

    drivers = induction_drivers(named["4OMe catalysed"], rate="v_peak")
    check("the clock/product coefficient is -0.025 +- 0.109 on 147 curves",
          abs(drivers["slope"] + 0.025) < 0.005
          and abs(drivers["stderr"] - 0.109) < 0.005
          and drivers["points"] == 147,
          f"{drivers['slope']:+.3f} +- {drivers['stderr']:.3f}, "
          f"n={drivers['points']}")
    check("product control is excluded by more than five standard errors",
          (drivers["slope"] - INDUCTION_PRODUCT_SLOPE) / drivers["stderr"] > 5,
          f"{(drivers['slope'] + 1) / drivers['stderr']:.1f} sigma")

    ratio = order_ratio(named["4OMe catalysed"])
    check("the order route agrees: -0.26 +- 0.32",
          abs(ratio["ratio"] + 0.26) < 0.02
          and abs(ratio["ratio_stderr"] - 0.32) < 0.02,
          f"{ratio['ratio']:+.3f} +- {ratio['ratio_stderr']:.3f}")
    check("the substrate order of the rate is +0.471 +- 0.084",
          abs(ratio["rate_order"] - 0.471) < 0.005
          and abs(ratio["rate_stderr"] - 0.084) < 0.005,
          f"{ratio['rate_order']:+.3f} +- {ratio['rate_stderr']:.3f}")

    depth = scope.orders("depth", frame=named["4OMe catalysed"],
                         floor=DEPTH_FLOOR)
    check("the amplitude has no substrate order: -0.114 +- 0.169",
          abs(depth["order_s0"] + 0.114) < 0.005
          and abs(depth["stderr_s0"] - 0.169) < 0.005,
          f"{depth['order_s0']:+.3f} +- {depth['stderr_s0']:.3f}")

    signal = signal_control(named["4OMe catalysed"])
    check("and it does not track signal-to-noise: +0.003 +- 0.149",
          abs(signal["signal_slope"] - 0.003) < 0.005
          and abs(signal["signal_stderr"] - 0.149) < 0.005,
          f"{signal['signal_slope']:+.3f} +- {signal['signal_stderr']:.3f}")
    for name, expected in (("BnOH two-axis (135-151)", 0.619),):
        got = signal_control(named[name])
        check(f"{name} fails the same control at {expected:+.3f}",
              abs(got["signal_slope"] - expected) < 0.005,
              f"{got['signal_slope']:+.3f}")

    ladder = induction.buffer_lever(table)
    joint = joint_buffer_order(ladder)
    check("the buffer axis MEETS the joint constraint: +1.094 +- 0.150",
          abs(joint["slope"] - 1.094) < 0.005
          and abs(joint["stderr"] - 0.150) < 0.005,
          f"{joint['slope']:+.3f} +- {joint['stderr']:.3f}")
    check("which is where the peroxide axis fails, so the two differ",
          joint["sigma"] < 1.0, f"{joint['sigma']:.1f} sigma")
    buffer_signal = signal_control(ladder)
    check("and the ladder passes its own signal control at -0.275 +- 0.268",
          abs(buffer_signal["signal_slope"] + 0.275) < 0.005
          and abs(buffer_signal["signal_stderr"] - 0.268) < 0.005,
          f"{buffer_signal['signal_slope']:+.3f} +- "
          f"{buffer_signal['signal_stderr']:.3f}")
    wide = signal_control(induction.buffer_lever(table, width=900.0))
    check("a 900 s window does not: the control fails and the joint order "
          "overshoots",
          wide["signal_slope"] < -0.5
          and joint_buffer_order(
              induction.buffer_lever(table, width=900.0))["slope"] > 1.2,
          f"{wide['signal_slope']:+.3f}")

    lever = induction.peroxide_lever(table)
    check("the peroxide block's induction order is the wrong sign for an adduct",
          lever["induction_order"] > 0,
          f"{lever['induction_order']:+.3f}")
    check("but its landmark tracks signal-to-noise, so it settles nothing",
          lever["signal_slope"] > 2 * lever["signal_stderr"],
          f"{lever['signal_slope']:+.3f} +- {lever['signal_stderr']:.3f}")

    gap = induction.activation_contrast()
    check("the induction is 12.0 kJ/mol below the turnover in free energy",
          abs(gap["gibbs_gap_kJ"] + 11.99) < 0.05,
          f"{gap['gibbs_gap_kJ']:+.2f}")
    check("which is 126x faster at 298 K",
          abs(gap["rate_ratio"] - 126) < 1, f"{gap['rate_ratio']:.0f}")
    check("the free-energy gap is many standard errors wide",
          abs(gap["gibbs_gap_kJ"]) > 10 * gap["gibbs_gap_stderr"],
          f"{gap['gibbs_gap_kJ']:+.2f} +- {gap['gibbs_gap_stderr']:.2f}")
    check("the entropy gap is NOT: it is about one standard error",
          abs(gap["entropy_gap_J"]) < 1.5 * gap["entropy_gap_stderr"],
          f"{gap['entropy_gap_J']:+.1f} +- {gap['entropy_gap_stderr']:.1f}")
    check("and neither is the enthalpy gap",
          abs(gap["enthalpy_gap_kJ"]) < 1.5 * gap["enthalpy_gap_stderr"],
          f"{gap['enthalpy_gap_kJ']:+.1f} +- {gap['enthalpy_gap_stderr']:.1f}")
    check("the induction's barrier is far too large to be physical",
          gap["induction"]["activation_kJ"] > 60.0,
          f"{gap['induction']['activation_kJ']:.1f} kJ/mol")

    windows = induction.landmark_window()
    narrow = windows[windows.window == "300 s"].iloc[0]
    tenth = windows[windows.window.str.startswith("0.10")].iloc[0]
    check("a 300 s window destroys the cold end: 529 s against 4289 s",
          abs(narrow.cold_s - 529) < 2 and abs(tenth.cold_s - 4289) < 2,
          f"{narrow.cold_s:.0f} against {tenth.cold_s:.0f}")
    check("and takes the activation energy with it, 24 against 86",
          narrow.activation_kJ < 40 and tenth.activation_kJ > 80,
          f"{narrow.activation_kJ:.1f} against {tenth.activation_kJ:.1f}")

    drift = induction.schedule_dependence(table)
    check("between runs the induction time tracks the schedule, +0.44",
          abs(drift["span"] - 0.437) < 0.01,
          f"{drift['span']:+.3f} +- {drift['span_stderr']:.3f}")
    check("and not pH",
          abs(drift["pH"]) < drift["pH_stderr"],
          f"{drift['pH']:+.3f} +- {drift['pH_stderr']:.3f}")

    for name, expected in (("4OMe catalysed", (147, 98)),
                           ("temperature series", (24, 22)),
                           ("BnOH two-axis (135-151)", (110, 46)),
                           ("4OMe enzyme-free", (49, 10))):
        block = named[name]
        block = block[block.live]
        check(f"{name}: {expected[1]} of {expected[0]} curves begin below "
              f"their eventual rate",
              (len(block), int(block.lag_first.sum())) == expected,
              f"{len(block)}, {int(block.lag_first.sum())}")
    arms = ladder_arms(named["BnOH two-axis (135-151)"])
    substrate = sign_drivers(arms["substrate arm"], axis="s0", control=True)
    check("two-axis: the sign tracks substrate at -0.112 +- 0.052",
          abs(substrate["s0"] + 0.112) < 0.005
          and abs(substrate["s0_stderr"] - 0.052) < 0.005,
          f"{substrate['s0']:+.3f} +- {substrate['s0_stderr']:.3f}")
    check("against a signal-to-noise lean of the opposite sign",
          substrate["signal"] > 0,
          f"{substrate['signal']:+.3f} +- {substrate['signal_stderr']:.3f}")
    four = sign_drivers(named["4OMe catalysed"], axis="s0", control=True)
    check("and the 4OMe block gives the opposite sign, +0.182 +- 0.073",
          abs(four["s0"] - 0.182) < 0.005
          and abs(four["s0_stderr"] - 0.073) < 0.005,
          f"{four['s0']:+.3f} +- {four['s0_stderr']:.3f}")
    ladder = induction.buffer_lever(table)
    check("the two buffer runs earn different model forms",
          set(ladder[ladder.experiment == 34].phases) == {2}
          and set(ladder[ladder.experiment == 32].phases) == {1},
          f"{ladder.groupby('experiment').phases.unique().to_dict()}")
    check("and the break is at the run boundary because the schedule is",
          ladder[ladder.experiment == 34].span_s.min()
          > 2.5 * ladder[ladder.experiment == 32].span_s.max(),
          f"{ladder.groupby('experiment').span_s.first().to_dict()}")
    step = induction.buffer_join_step(ladder)
    check("the rate falls 1.80x across the join, so the runs need levels",
          abs(step["step"] - 1.80) < 0.01, f"{step['step']:.3f}")
    order = induction.buffer_order(ladder)
    check("read on a common footing the two runs agree in sign",
          order["slope_32"] < 0 and order["slope_34"] < 0,
          f"{order['slope_34']:+.3f} and {order['slope_32']:+.3f}")
    check("and the pooled slope is -0.433 +- 0.201",
          abs(order["slope"] + 0.433) < 0.005
          and abs(order["stderr"] - 0.201) < 0.005,
          f"{order['slope']:+.3f} +- {order['stderr']:.3f}")
    swept = [induction.buffer_order(
        induction.buffer_lever(table, width=width))["slope"]
        for width in induction.BUFFER_WINDOW_SWEEP]
    check("no window from 300 s to 1200 s changes its sign",
          all(value < 0 for value in swept),
          " ".join(f"{value:+.3f}" for value in swept))

    fixed = induction.substrate_order_corrected(
        named["4OMe catalysed"], order["slope"], order["stderr"])
    check("correcting route two for the buffer moves it towards the threshold",
          abs(fixed["corrected"] - fixed["threshold"])
          < abs(fixed["measured"] - fixed["threshold"]),
          f"{fixed['measured']:+.3f} -> {fixed['corrected']:+.3f} against "
          f"{fixed['threshold']:+.3f}")
    check("and it stops excluding it, which the document says",
          abs(fixed["corrected"] - fixed["threshold"])
          < 2 * fixed["corrected_stderr"],
          f"{abs(fixed['corrected'] - fixed['threshold']) / fixed['corrected_stderr']:.2f} sigma")
    drivers_again = induction_drivers(named["4OMe catalysed"], rate="v_peak")
    check("while route one, whose regressor is the rate, is untouched",
          (drivers_again["slope"] - INDUCTION_PRODUCT_SLOPE)
          / drivers_again["stderr"] > 5)

    lever = substrate_lever(named["4OMe catalysed"])
    check("the in-run substrate lever is 4.0x on 28 of 38 experiments",
          abs(lever["median_lever"] - 4.0) < 0.01
          and lever["laddered"] == 28 and lever["experiments"] == 38,
          f"{lever}")


if __name__ == "__main__":
    test_the_landmark_measures_what_it_claims()
    test_it_agrees_with_the_gated_statistic()
    test_induction_drivers_tell_a_clock_from_a_product()
    test_the_joint_order_reads_back_the_scheme_it_is_built_from()
    test_the_joint_buffer_order_reads_back_its_own_scheme()
    test_the_joint_order_takes_a_fitted_clock_and_a_gate()
    test_the_joint_order_control_axis_reads_the_substrate_order()
    test_the_two_named_orders_are_the_general_one()
    test_the_saturation_fit_recovers_a_planted_constant()
    test_the_trap_constant_inverts_its_own_order()
    test_the_floor_does_not_carry_the_answer()
    test_orders_refuse_an_axis_the_block_does_not_move()
    test_the_sign_comes_off_the_fit_and_not_off_the_depth()
    test_the_L_has_to_be_split_before_it_is_read()
    test_the_blocks_differ_in_their_buffer_collinearity()
    test_the_window_free_lag_recovers_a_planted_relaxation()
    test_the_window_free_lag_is_not_windowed()
    test_a_between_run_ladder_refuses_offsets_that_absorb_it()
    test_the_common_window_removes_the_schedule()
    test_the_barrier_survives_the_schedule_and_the_naive_control_does_not()
    test_the_replicate_floor_is_the_bar_between_runs()
    test_lag_orders_agree_with_the_package_order_machinery()
    test_the_signal_collinearity_says_which_axis_can_be_asked()
    test_the_pooled_ladder_notices_a_disagreement()
    test_the_saturation_fraction_bounds_both_schemes()
    test_the_identifiability_table_names_the_between_run_axes()
    test_the_product_at_the_landmark_tells_a_clock_from_a_threshold()
    test_the_product_test_is_calibrated_against_its_own_bias()
    test_the_product_test_agrees_with_the_driver_regression()
    test_regressions()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
