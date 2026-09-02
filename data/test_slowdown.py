"""
Tests for slowdown.py.

The claims this module makes are all of the form "these data are better
described by X than by Y", so the tests are of the same form: build a curve
from a known X, and check that the machinery says X and not Y. Two of them
would have caught real mistakes made while writing it -- the landmark that
reported a fall on curves that accelerate throughout, and a driver regression
that cannot tell a clock from a product when the design does not separate them.

    python data/test_slowdown.py
"""
import sys
import os

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import slowdown
from slowdown import (SinkFit, _clock_shape, _inhibition_shape, _sink_shape,
                      deceleration_drivers, fit_slowdown, plateau_scaling,
                      selectivity, sink_fit, verdict)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


class FakeCurve:
    """The five attributes slowdown reads off a Curve, and nothing else."""

    def __init__(self, times, absorbance, noise=1e-5, source="rre",
                 experiment=1, sample=1, s0=5.0, e0=0.25):
        self.times = np.asarray(times, dtype=float)
        self.absorbance = np.asarray(absorbance, dtype=float)
        self.noise = noise
        self.source = source
        self.experiment = experiment
        self.sample = sample
        self.substrate = "4OMe-BnOH"
        self.buffer = "Phosphate"
        self.pH = 7.0
        self.temperature = 25.0
        self.epsilon = 7.53
        self.buf = 65.0
        self.conditions = type("C", (), {"s0": s0, "e0": e0, "h2o2": 82.5})()


def sink_curve(v=5e-5, k=1e-4, tau=300.0, span=6000.0, points=360, seed=0):
    """A' = v(1 - e^(-t/tau)) - k A, with noise."""
    times = np.linspace(0.0, span, points)
    values = v * _sink_shape(times, tau, k)
    rng = np.random.default_rng(seed)
    return times, values + rng.normal(0.0, 1e-5, points)


def test_shapes_are_the_solutions_they_claim():
    """Each closed form must satisfy its own differential equation."""
    print("\nclosed forms")
    times = np.linspace(1e-3, 4000.0, 40001)
    tau, k, v = 400.0, 3e-4, 4e-5
    sink = v * _sink_shape(times, tau, k)
    residual = (np.gradient(sink, times)
                - (v * (1 - np.exp(-times / tau)) - k * sink))
    check("sink solves A' = v(1-e^(-t/tau)) - kA",
          float(np.max(np.abs(residual[10:]))) < 1e-9 * v / 1e-5,
          f"max {np.max(np.abs(residual[10:])):.3e}")
    clock = v * _clock_shape(times, tau, k)
    residual = (np.gradient(clock, times)
                - v * (1 - np.exp(-times / tau)) * np.exp(-k * times))
    check("clock solves A' = v(1-e^(-t/tau))e^(-kt)",
          float(np.max(np.abs(residual[10:]))) < 1e-9 * v / 1e-5,
          f"max {np.max(np.abs(residual[10:])):.3e}")
    ki, depth = 0.02, v / 0.02
    inhibited = ki * _inhibition_shape(times, tau, depth)
    residual = (np.gradient(inhibited, times)
                - v * (1 - np.exp(-times / tau)) / (1 + inhibited / ki))
    check("inhibition solves A' = v(1-e^(-t/tau))/(1+A/Ki)",
          float(np.max(np.abs(residual[10:]))) < 1e-9 * v / 1e-5,
          f"max {np.max(np.abs(residual[10:])):.3e}")


def test_the_contest_separates_what_it_can_and_admits_what_it_cannot():
    """
    A plateauing curve must not be called inhibition, and the reverse.

    `clock` and `sink` are NOT separated, and this test asserts that they are
    not: both approach a horizontal asymptote through two exponentials, and on
    one curve's shape there is nothing to choose between them. A clock curve
    here is fitted slightly better by the sink form than by the form it was
    built from. The module says so in its docstring and separates them by
    where else the decay appears -- in the enzyme-free channel -- rather than
    by residuals. What the residuals CAN decide is whether the curve reaches a
    plateau at all, which is the only thing they are used for.
    """
    print("\nmodel recovery")
    times = np.linspace(0.0, 6000.0, 360)
    rng = np.random.default_rng(1)
    plateaus = {
        "sink": 5e-5 * _sink_shape(times, 300.0, 2e-4),
        "clock": 5e-5 * _clock_shape(times, 300.0, 2e-4),
    }
    for name, clean in plateaus.items():
        values = clean + rng.normal(0.0, 2e-5, len(times))
        sse = {n: fit_slowdown(times, values, n).sse
               for n in ("clock", "sink", "inhibition")}
        best = min(sse, key=sse.get)
        check(f"a {name} curve is not called inhibition",
              best != "inhibition" and sse["inhibition"] > 2 * sse[best],
              ", ".join(f"{k}={v:.3e}" for k, v in sse.items()))
    inhibited = 0.02 * _inhibition_shape(times, 300.0, 5e-5 / 0.02)
    values = inhibited + rng.normal(0.0, 2e-5, len(times))
    sse = {n: fit_slowdown(times, values, n).sse
           for n in ("clock", "sink", "inhibition")}
    won, _ = verdict({k: np.sqrt(v) for k, v in sse.items()})
    check("an inhibition curve is won by inhibition", won == "inhibition",
          f"got {won}: " + ", ".join(f"{k}={v:.3e}" for k, v in sse.items()))
    check("clock and sink are not separable on one curve",
          abs(sse["clock"] - sse["sink"]) < max(sse["clock"], sse["sink"]),
          f"clock {sse['clock']:.3e}, sink {sse['sink']:.3e}")


def test_sink_fit_recovers_its_constants():
    """The linearisation has to return the k and plateau it was built from."""
    print("\nsink linearisation")
    v, k = 5e-5, 1.5e-4
    times, values = sink_curve(v=v, k=k, span=20000.0, points=400)
    fitted = sink_fit(FakeCurve(times, values))
    check("k within 15%", abs(fitted.k - k) / k < 0.15,
          f"{fitted.k:.3e} against {k:.3e}")
    check("plateau within 15%", abs(fitted.plateau - v / k) / (v / k) < 0.15,
          f"{fitted.plateau:.4f} against {v / k:.4f}")
    check("prefers the sink", fitted.prefers == "sink",
          f"{fitted.prefers}: rate {fitted.rate_r2:.4f} "
          f"reciprocal {fitted.reciprocal_r2:.4f}")

    ki, rate = 0.02, 5e-5
    inhibited = ki * _inhibition_shape(times, 300.0, rate / ki)
    rng = np.random.default_rng(3)
    fitted = sink_fit(FakeCurve(times, inhibited + rng.normal(0, 1e-5,
                                                              len(times))))
    check("an inhibition curve is not called a sink",
          fitted.prefers in ("inhibition", "tied"),
          f"{fitted.prefers}: rate {fitted.rate_r2:.4f} "
          f"reciprocal {fitted.reciprocal_r2:.4f}")


def test_a_rising_curve_has_no_fall():
    """
    A lag curve must not report a decay landmark.

    `decay_point` did, before the rate was read through a window: the fitted
    form is free to put a spike at t = 0 and the landmark was measured against
    it, so curves that accelerate from start to finish were reported as
    falling to three quarters of their peak inside the first minute.
    """
    print("\nno fall on a rising curve")
    times = np.linspace(0.0, 6000.0, 360)
    values = 5e-5 * (times - 900.0 * (1 - np.exp(-times / 900.0)))
    rng = np.random.default_rng(5)
    curve = FakeCurve(times, values + rng.normal(0.0, 1e-5, len(times)))
    t_fall, a_fall, peak, _ = slowdown.decay_point(curve, fraction=0.75)
    check("no landmark on a pure lag", not np.isfinite(t_fall),
          f"t_fall {t_fall}, a_fall {a_fall}")
    check("its peak rate is still reported", np.isfinite(peak) and peak > 0)
    check("sink_fit declines to read it", sink_fit(curve).points == 0)


def _driver_frame(rule, spans, rates, seed=7):
    """A frame of synthetic curves obeying one law, for the regression."""
    rng = np.random.default_rng(seed)
    rows = []
    for index, (span, rate) in enumerate(zip(spans, rates)):
        product = rate * span
        if rule == "clock":
            ratio = np.exp(-span / 8000.0)
        else:
            ratio = 1.0 / (1.0 + product / 0.05)
        ratio *= np.exp(rng.normal(0.0, 0.05))
        rows.append({"experiment": 100 + index // 4, "live": True,
                     "late_over_early": ratio, "net": product,
                     "duration_s": span, "epsilon": 7.53,
                     "substrate": "4OMe-BnOH", "buffer": "Phosphate",
                     "e0": 0.25, "temperature": 25.0, "pH": 7.0,
                     "differential": True})
    return pd.DataFrame(rows)


def test_the_driver_regression_tells_a_clock_from_a_product():
    """
    On a design that separates them, each law must be read back.

    This is the test that matters: everything the module concludes rests on
    this regression being able to distinguish "fell with time" from "fell with
    product", and it can only do that if run length and product are varied
    independently -- which the archive does by accident and this test does on
    purpose.
    """
    print("\ndriver regression")
    rng = np.random.default_rng(11)
    spans = np.exp(rng.uniform(np.log(1500), np.log(18000), 80))
    rates = np.exp(rng.uniform(np.log(3e-6), np.log(3e-4), 80))
    for rule, expected in (("clock", "span"), ("product", "product")):
        frame = _driver_frame(rule, spans, rates)
        row = deceleration_drivers(frame)
        other = "product" if expected == "span" else "span"
        check(f"a {rule} law loads on {expected}",
              abs(row[expected]) > 4 * row[expected + "_stderr"]
              and row[expected] < 0,
              f"{row[expected]:+.3f} +/- {row[expected + '_stderr']:.3f}")
        check(f"a {rule} law does not load on {other}",
              abs(row[other]) < 4 * row[other + "_stderr"]
              or row[other] > 0,
              f"{row[other]:+.3f} +/- {row[other + '_stderr']:.3f}")


def test_selectivity_and_plateau_scaling():
    """Both are arithmetic on quantities defined elsewhere; pin them."""
    print("\nderived quantities")
    check("k_A/k_S is [S] over the plateau in mM",
          abs(selectivity(0.753, 5.0, 7.53) - 50.0) < 1e-9,
          f"{selectivity(0.753, 5.0, 7.53)}")
    check("a nonpositive plateau has no selectivity",
          not np.isfinite(selectivity(-0.1, 5.0, 7.53)))
    s0 = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    table = pd.DataFrame({"s0": s0, "plateau": 0.3 * s0 ** 0.6,
                          "points": 10, "epsilon": 7.53})
    order = plateau_scaling(table)
    check("plateau_scaling recovers a planted order",
          abs(order["order"] - 0.6) < 1e-6, f"{order['order']}")


def test_pinning_k_recovers_the_production_rate():
    """
    With k supplied, the sink model must return the v it was built from.

    This is the whole basis of `production_frame`: the published rate is net
    of the sink, and `v` is recovered by fitting the sink model with its rate
    constant pinned. If the pinned fit could not recover a planted v there
    would be nothing to correct with.
    """
    print("\npinned sink fit")
    v, k = 5e-5, 1.5e-4
    times, values = sink_curve(v=v, k=k, span=12000.0, points=400, seed=9)
    pinned = fit_slowdown(times, values, "sink", decay=k)
    check("v within 10%", abs(pinned.extra["rate"] - v) / v < 0.10,
          f"{pinned.extra['rate']:.3e} against {v:.3e}")
    check("k is the one supplied", abs(pinned.decay - k) < 1e-12,
          f"{pinned.decay:.3e}")
    # And the free fit is the thing pinning exists to avoid: on a run that has
    # not turned over, (v, k) trade off freely along v/k = plateau.
    short, values = sink_curve(v=v, k=k, span=1200.0, points=200, seed=9)
    free = fit_slowdown(short, values, "sink")
    check("the free fit on a short run is not to be trusted",
          not (0.5 < free.extra["rate"] / v < 2.0)
          or not (0.5 < free.decay / k < 2.0),
          f"free v {free.extra['rate']:.3e}, k {free.decay:.3e} -- if this "
          f"passes, pinning may no longer be necessary")


def test_the_tail_has_to_be_longer_than_its_own_window():
    """
    SINK_MINIMUM_POINTS counts overlapping window positions and cannot catch a
    tail shorter than one window. This is the guard that can.

    Planted: a curve whose rate turns over only at the very end, so its tail is
    a fraction of a wide window. Without the guard it yields a `k` drawn
    through less than one independent measurement of the rate.
    """
    print("\nthe tail's own length")
    times, values = sink_curve(v=5e-5, k=8e-5, tau=300.0, span=6000.0,
                               points=360, seed=3)
    curve = FakeCurve(times, values)
    wide = sink_fit(curve, fraction=0.45)
    narrow = sink_fit(curve, fraction=0.15)
    check("a window wide enough to swallow the tail is refused",
          not np.isfinite(wide.k), f"k = {wide.k}")
    check("and the same curve at a workable window is not",
          np.isfinite(narrow.k) and narrow.windows >= slowdown.SINK_MINIMUM_WINDOWS,
          f"k = {narrow.k:.3g}, {narrow.windows:.2f} windows")
    check("the threshold is stated, not implied",
          slowdown.SINK_MINIMUM_WINDOWS >= 1.0,
          f"{slowdown.SINK_MINIMUM_WINDOWS}")


def test_the_sink_window_does_not_carry_the_null():
    """
    The published null must hold at every window, and the systematic must be
    reported rather than hidden -- it is three times the statistical error.
    """
    print("\nthe sink constant's window")
    sweep = slowdown.sink_window_sensitivity()
    check("the null holds at every width", sweep.attrs["null_holds"],
          f"{sweep.shift_kJ.round(2).tolist()}")
    check("the activation energy does not: it spreads by 30 kJ/mol",
          abs(sweep.attrs["activation_spread"] - 30.1) < 0.5,
          f"{sweep.attrs['activation_spread']:.1f}")
    check("which is larger than the error quoted at the published window",
          sweep.attrs["activation_spread"]
          > 2 * float(sweep[np.isclose(sweep.window,
                                       slowdown.SINK_WINDOW)].stderr_kJ.iloc[0]))
    check("and the published width still admits its 8 curves",
          int(sweep[np.isclose(sweep.window, slowdown.SINK_WINDOW)]
              .curves.iloc[0]) == 8, f"{sweep.curves.tolist()}")


def test_regressions():
    """The published numbers, so a refactor cannot move them quietly."""
    print("\npublished numbers")
    import scope
    frame = scope.frame(tuple(range(1, 152)))
    named = slowdown.substrate_blocks(frame)
    series = deceleration_drivers(named["temperature series"])
    check("temperature series loads on product, not on span",
          series["product"] < -0.3 and abs(series["span"]) < 0.3,
          f"product {series['product']:+.3f}, span {series['span']:+.3f}")
    free = deceleration_drivers(named["4OMe enzyme-free"])
    check("the enzyme-free block loads on span, not on product",
          free["span"] < -0.2 and free["product"] > -0.1,
          f"product {free['product']:+.3f}, span {free['span']:+.3f}")
    bnoh = deceleration_drivers(named["BnOH catalysed, all buffers"], fixed=True)
    four = deceleration_drivers(named["4OMe catalysed, phosphate"], fixed=True)
    gap = bnoh["product"] - four["product"]
    error = np.hypot(bnoh["product_stderr"], four["product_stderr"])
    check("the two substrates differ by more than 4 sigma", gap > 4 * error,
          f"{gap:.3f} +/- {error:.3f}")

    effect = slowdown.sink_effect_on_activation()
    check("naming the fall does not move the activation energy",
          abs(effect["activation_shift"]) < effect["activation_shift_stderr"],
          f"{effect['activation_shift']:+.2f} +/- "
          f"{effect['activation_shift_stderr']:.2f} kJ/mol")
    check("the corrected route is the noisier one, so it is not the headline",
          effect["corrected"]["rms"] > effect["published"]["rms"],
          f"{effect['corrected']['rms']:.3f} against "
          f"{effect['published']['rms']:.3f}")
    check("the production rate is above the published one",
          effect["lift"] > 1.0, f"{effect['lift']:.3f}")


if __name__ == "__main__":
    test_shapes_are_the_solutions_they_claim()
    test_the_contest_separates_what_it_can_and_admits_what_it_cannot()
    test_sink_fit_recovers_its_constants()
    test_a_rising_curve_has_no_fall()
    test_the_driver_regression_tells_a_clock_from_a_product()
    test_selectivity_and_plateau_scaling()
    test_pinning_k_recovers_the_production_rate()
    test_the_tail_has_to_be_longer_than_its_own_window()
    test_the_sink_window_does_not_carry_the_null()
    test_regressions()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
