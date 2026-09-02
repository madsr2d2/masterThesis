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
from induction import (DEPTH_FLOOR, INDUCTION_CLOCK_SLOPE, INDUCTION_FLOOR,
                       INDUCTION_PRODUCT_SLOPE, induction_drivers,
                       induction_point, order_ratio, signal_control,
                       substrate_lever)

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
    check("the in-scope orders are untouched by the guard",
          abs(whole["order_h2o2"] - 0.7936) < 5e-4
          and abs(whole["order_s0"] - 0.0906) < 5e-4,
          f"{whole['order_s0']:.4f}, {whole['order_h2o2']:.4f}")


_CACHE = {}


def _induction_table():
    if "table" not in _CACHE:
        _CACHE["table"] = induction.induction_table(induction.WHOLE_ARCHIVE)
    return _CACHE["table"]


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
    for name, expected in (("BnOH in scope (135-151)", 0.619),):
        got = signal_control(named[name])
        check(f"{name} fails the same control at {expected:+.3f}",
              abs(got["signal_slope"] - expected) < 0.005,
              f"{got['signal_slope']:+.3f}")

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

    lever = substrate_lever(named["4OMe catalysed"])
    check("the in-run substrate lever is 4.0x on 28 of 38 experiments",
          abs(lever["median_lever"] - 4.0) < 0.01
          and lever["laddered"] == 28 and lever["experiments"] == 38,
          f"{lever}")


if __name__ == "__main__":
    test_the_landmark_measures_what_it_claims()
    test_it_agrees_with_the_gated_statistic()
    test_induction_drivers_tell_a_clock_from_a_product()
    test_the_floor_does_not_carry_the_answer()
    test_orders_refuse_an_axis_the_block_does_not_move()
    test_regressions()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
