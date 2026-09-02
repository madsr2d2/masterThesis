"""
Tests for summary_kinetics.py.

Two kinds of check. The first is recovery: synthetic curves built from a known
power law, with the real between-experiment offsets and the real noise put back
in, and the regression has to find the exponents. A regression that has never
recovered a known answer is evidence about linear algebra, not about chemistry.

The second kind guards the diagnostics, which is where this module earns its
keep. The regression will always return numbers; the diagnostics are what say
whether the numbers mean anything. Two of them are pinned against designs
constructed to be degenerate on purpose -- an axis that never moves inside an
experiment, and a cuvette that died -- because both occur in the real data and
both silently corrupt an offset model.

    python data/test_summary_kinetics.py
"""
import sys

import numpy as np
import pandas as pd

from fit_dataset import build_curves, group_curves
from summary_kinetics import (ABSORPTION_FLOOR, BURST_TAU_CAP, DEAD_CURVE_SNR,
                              DESIGN_AXES, INITIAL_WINDOW, fit_burst,
                              fit_progress, fit_two_phase, line_fit,
                              _two_phase_design,
                              OUTLIER_FACTOR, REPLICATE_RSD, absorbed_axes,
                              buffer_concentrations, experiment_outliers,
                              initial_rate, line_slope, profile_km, regress,
                              replicate_scatter, slope_ratio, summarise,
                              to_frame, window_sensitivity,
                              within_experiment_variation)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def synthetic_frame(order_s=0.30, order_hoo=1.10, order_buf=0.0,
                    offset_spread=0.7, sigma=0.045, seed=0):
    """
    Curve summaries from a known law, with a per-experiment offset put back in.

    log10 v0 = c_experiment + order_s.log[S] + order_hoo.log[HOO] + order_buf.log[buf]

    The offsets are drawn at the spread the real data shows (about 500x between
    the extreme experiments), so a regression that ignores them faces the same
    problem the real one does.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for experiment in range(1, 16):
        offset = rng.normal(-4.0, offset_spread)
        for sample in range(1, 5):
            s0 = float(10 ** rng.uniform(-1.0, 1.5))
            hoo = float(10 ** rng.uniform(-2.8, -2.0))
            buf = float(10 ** rng.uniform(0.5, 2.3))
            log_v = (offset + order_s * np.log10(s0) + order_hoo * np.log10(hoo)
                     + order_buf * np.log10(buf) + rng.normal(0.0, sigma))
            rows.append(dict(experiment=experiment, sample=sample, s0=s0,
                             h2o2=82.5, hoo=hoo, buf=buf, e0=0.0,
                             v0=float(10 ** log_v), v0_stderr=0.0,
                             window_rms=0.0, noise=6e-4, slope_ratio=1.0,
                             points=200, duration=1000.0, amplitude=0.05,
                             conversion=0.005))
    return pd.DataFrame(rows)


def test_line_slope():
    print("line fitting")
    times = np.linspace(0.0, 1000.0, 200)
    slope, stderr, rms = line_slope(times, 3.0 + 2.5e-5 * times)
    check("recovers an exact slope", abs(slope - 2.5e-5) < 1e-12,
          f"got {slope:.6e}")
    check("an exact line has zero residual", rms < 1e-14, f"got {rms:.3e}")

    rng = np.random.default_rng(1)
    noisy = 3.0 + 2.5e-5 * times + rng.normal(0.0, 6e-4, len(times))
    slope, stderr, rms = line_slope(times, noisy)
    check("recovers a noisy slope inside 3 sigma",
          abs(slope - 2.5e-5) < 3 * stderr, f"got {slope:.3e} +/- {stderr:.3e}")
    check("residual matches the noise put in", 4e-4 < rms < 8e-4, f"got {rms:.3e}")
    check("standard error is positive", stderr > 0)

    slope, stderr, rms = line_slope(np.array([0.0, 0.0, 0.0]), np.array([1.0, 2.0, 3.0]))
    check("zero time span returns nan rather than dividing by zero", np.isnan(slope))


def test_window():
    print("windowing")
    times = np.linspace(0.0, 1000.0, 100)
    # rate halves halfway through: early slope 2e-5, late slope 1e-5
    values = np.where(times <= 500, 2e-5 * times, 1e-2 + 1e-5 * (times - 500))
    early, *_ = initial_rate(times, values, 0.2)
    check("initial rate uses the early window only",
          abs(early - 2e-5) < 1e-9, f"got {early:.3e}")
    ratio = slope_ratio(times, values, 0.2)
    check("slope ratio sees the deceleration", abs(ratio - 0.5) < 1e-6,
          f"got {ratio:.4f}")
    check("slope ratio is 1 on a straight line",
          abs(slope_ratio(times, 3.0 + 2e-5 * times, 0.2) - 1.0) < 1e-9)
    check("the default window is the one the docstring justifies",
          INITIAL_WINDOW == 0.20)

    short = np.linspace(0.0, 10.0, 6)
    slope, *_ = initial_rate(short, 1e-4 * short, 0.01)
    check("a tiny fraction still leaves enough points for a slope",
          np.isfinite(slope))


def test_line_intercept_and_floor():
    print("line intercept and the quantisation floor")
    times = np.linspace(0.0, 1000.0, 40)
    intercept, slope, stderr, rms = line_fit(times, -0.002 + 2.5e-5 * times)
    check("line_fit returns the intercept the fit chose",
          abs(intercept - (-0.002)) < 1e-12, f"got {intercept:.6f}")
    check("the intercept is what a drawn fit line needs",
          abs((intercept + slope * times[-1]) - (-0.002 + 2.5e-5 * times[-1])) < 1e-12)

    # Points sitting exactly on a line: residual is zero, but the instrument
    # only reports three decimals, so the standard error must not be zero too.
    exact = np.arange(5) * 0.004
    _, _, floored, residual = line_fit(np.arange(5) * 49.0, exact)
    check("an exactly collinear window still gets a positive standard error",
          floored > 1e-9, f"got {floored:.3e}")
    check("...and its residual is genuinely zero", residual < 1e-15,
          f"got {residual:.3e}")
    check("the floored error is of a believable size",
          1e-7 < floored < 1e-4, f"got {floored:.3e}")


def test_burst_fit():
    print("the burst / lag form")
    times = np.linspace(0.0, 3000.0, 150)

    def build(c, v_ss, B, tau, sigma=0.0, seed=0):
        clean = c + v_ss * times - B * (1 - np.exp(-times / tau))
        if sigma:
            clean = clean + np.random.default_rng(seed).normal(0.0, sigma, len(times))
        return clean

    truth = dict(c=0.01, v_ss=4e-5, B=0.02, tau=400.0)
    fit = fit_burst(times, build(**truth))
    for name, value in truth.items():
        got = getattr(fit, name)
        check(f"recovers {name}", abs(got - value) <= 0.05 * abs(value) + 1e-9,
              f"got {got:.5g}, want {value:.5g}")
    check("tau is resolved when the transient sits inside the window", fit.resolved,
          f"interval {fit.tau_interval}")
    check("B > 0 is reported as a lag", fit.kind == "lag", f"got {fit.kind}")
    check("v0 obeys v0 = v_ss - B/tau",
          abs(fit.v0 - (fit.v_ss - fit.B / fit.tau)) < 1e-12)
    check("v0 is below v_ss for a lag", fit.v0 < fit.v_ss)
    check("lag time is B/v_ss",
          abs(fit.lag_time - fit.B / fit.v_ss) < 1e-9)
    # tau is profiled on a discrete grid, so it lands near 400 rather than on
    # it. The scale that matters is the instrument's 0.001 AU quantum, not
    # machine precision: a reconstruction closer than that is exact as far as
    # this data can tell.
    deviation = float(np.max(np.abs(fit.predict(times) - build(**truth))))
    check("predict reproduces the curve well inside one quantisation step",
          deviation < 1e-4, f"max deviation {deviation:.2e} AU")

    burst = fit_burst(times, build(0.01, 4e-5, -0.02, 400.0))
    check("B < 0 is reported as a burst", burst.kind == "burst", f"got {burst.kind}")
    check("v0 is above v_ss for a burst", burst.v0 > burst.v_ss)

    # tau far beyond the run: the transient is collinear with the linear term
    slow = fit_burst(times, build(0.01, 4e-5, 0.4, 60000.0, sigma=6e-4, seed=2))
    check("tau beyond the run length is reported UNRESOLVED", not slow.resolved,
          f"tau={slow.tau:.0f} interval {slow.tau_interval}")
    check("an unresolved fit says so through .kind",
          slow.kind == "unresolved", f"got {slow.kind}")

    straight = fit_burst(times, 0.01 + 4e-5 * times + np.random.default_rng(4)
                         .normal(0.0, 6e-4, len(times)))
    check("a straight line leaves tau unresolved", not straight.resolved)
    check("a straight line's v_ss is still its slope",
          abs(straight.v_ss - 4e-5) < 6e-6, f"got {straight.v_ss:.3e}")

    check("too few points returns a blank fit rather than raising",
          not np.isfinite(fit_burst(times[:4], build(**truth)[:4]).tau))
    check("the grid cap is above 1, so tau may exceed the run length",
          BURST_TAU_CAP > 1.0)


def test_bounded_v0_does_not_mean_a_resolved_shape():
    """
    `bounded` asks about v0 only. tau can still be pinned at a grid end.

    BoundedBurstFit was written without tau's profile, replacing
    `fit_burst.resolved` with a v0-only criterion, and that lost real
    information: on the enzyme-free BnOH set tau is pinned on 15 of 27 curves
    and 11 of those still report `bounded`. Once tau collapses to the floor or
    runs to the cap the model degenerates -- to a step, to a straight line --
    and v0 -> v_ss becomes exactly determined, so v0 is trustworthy while the
    BURST is not. Restored as `tau_resolved` / `shape_is_meaningful`.
    """
    print("\nbounded v0 does not imply a resolved shape")
    from fit_dataset import build_curves, source_floor
    from summary_kinetics import fit_burst_bounded

    curves = [c for c in build_curves()[0]
              if c.experiment in (3, 6, 65, 67, 69, 70)]
    fits = [fit_burst_bounded(c.times, c.absorbance,
                              noise_floor=source_floor(c.source))
            for c in curves]

    bounded = sum(f.bounded for f in fits)
    resolved = sum(f.shape_is_meaningful for f in fits)
    both = sum(1 for f in fits if f.bounded and f.shape_is_meaningful)
    check("more curves bound v0 than resolve tau",
          bounded > resolved, f"bounded {bounded}, resolved {resolved}")
    check("and the two are not the same set",
          both < bounded, f"both {both}, bounded {bounded}")
    check("some curves report bounded v0 with an unresolved tau",
          any(f.bounded and not f.shape_is_meaningful for f in fits),
          "none found")

    # A pinned tau must be reported pinned, not silently quoted.
    for fit in fits:
        if not fit.shape_is_meaningful:
            continue
        check_once = (fit.tau_low <= fit.tau <= fit.tau_high)
        if not check_once:
            check("a resolved tau lies inside its own profile interval", False,
                  f"tau {fit.tau:.4g} outside [{fit.tau_low:.4g}, "
                  f"{fit.tau_high:.4g}]")
            break
    else:
        check("every resolved tau lies inside its own profile interval", True)

    # The two curves with a negative settled rate are exactly the trap: both
    # report a bounded v0, and neither has a meaningful shape.
    backwards = [f for f in fits if f.settles_backwards]
    check("the v_ss < 0 fits all report bounded v0",
          backwards and all(f.bounded for f in backwards),
          f"{len(backwards)} such fits")
    check("...and none of them has a resolved tau",
          all(not f.shape_is_meaningful for f in backwards),
          f"{sum(f.shape_is_meaningful for f in backwards)} resolved")


def test_lag_branch_is_gated_per_curve():
    """
    The B <= 0 constraint applies only where a curve does not accelerate.

    A blanket `B <= 0` was added on 2026-09-01 to stop the burst form returning
    negative initial rates, justified by "0 of 16 pass the acceleration test".
    That was measured on the constant-buffer runs and generalised to all 27
    without checking. (The 16 now score 1, exp 67 sample 3 having crossed the
    gate when the first reading of every run was dropped -- which is a second
    reason not to hard-code a count as a constraint.) Exps 3 and 6 hold four curves that DO accelerate, two at
    z = +8.4 and +11.8, and the blanket rule bound on two of them -- forcing a
    decelerating shape onto curves whose own z-score says they rise.

    "auto" asks each curve instead. This test pins that it asks, that it
    forbids where it should, and that it still admits no negative v0.
    """
    print("\nthe lag branch is gated per curve")
    from curve_metrics import ACCELERATION_SIGMA, acceleration
    from fit_dataset import build_curves, source_floor
    from summary_kinetics import fit_burst_bounded

    curves = [c for c in build_curves()[0] if c.experiment in (3, 6, 65, 67, 69, 70)]
    check("the enzyme-free BnOH set is 27 cuvettes", len(curves) == 27,
          f"got {len(curves)}")

    accelerating = [c for c in curves
                    if acceleration(c.times, c.absorbance,
                                    floor=source_floor(c.source))[0]
                    > ACCELERATION_SIGMA]
    check("some of them genuinely accelerate, so a blanket B<=0 is wrong",
          len(accelerating) >= 3, f"{len(accelerating)} accelerate")
    # They were all in the titration runs until 2026-09-01, when dropping the
    # first reading moved exp 67 sample 3 from z = +2.9 to +3.5 and across the
    # gate. That is the 3-sigma cut being a hard edge on a continuous
    # statistic, not a change in the chemistry.
    check("most of them are in the titration runs",
          sum(1 for c in accelerating if c.experiment in (3, 6))
          >= len(accelerating) - 1,
          f"experiments {sorted(c.experiment for c in accelerating)}")

    def run(mode):
        return [fit_burst_bounded(c.times, c.absorbance, constrain=mode,
                                  noise_floor=source_floor(c.source))
                for c in curves]

    auto, always, never = run("auto"), run(True), run(False)
    check("unconstrained, some v0 come back negative",
          sum(f.v0 < 0 for f in never) == 4,
          f"{sum(f.v0 < 0 for f in never)} negative")
    # "auto" admits a negative v0 only where the curve's own z says a lag is
    # real -- exp 67 sample 3, at z = +3.5. It is not quotable: `bounded` is
    # False and tau is unresolved, so the safety net that matters still holds.
    negative_auto = [f for f in auto if f.v0 < 0]
    check("auto admits far fewer negative v0 than leaving the branch open",
          len(negative_auto) < sum(f.v0 < 0 for f in never),
          f"auto {len(negative_auto)}, unconstrained "
          f"{sum(f.v0 < 0 for f in never)}")
    check("and any it does admit is flagged unbounded",
          all(not f.bounded for f in negative_auto),
          f"{sum(f.bounded for f in negative_auto)} bounded")
    check("auto bounds more curves than leaving the branch open",
          sum(f.bounded for f in auto) > sum(f.bounded for f in never),
          f"{sum(f.bounded for f in auto)} vs {sum(f.bounded for f in never)}")
    check("and fewer than shutting it everywhere, which is the point",
          sum(f.bounded for f in auto) < sum(f.bounded for f in always),
          f"{sum(f.bounded for f in auto)} vs {sum(f.bounded for f in always)}")

    # On an accelerating curve, "auto" must leave the branch open, so it has to
    # agree with the unconstrained fit and may differ from the blanket one.
    index = {id(c): i for i, c in enumerate(curves)}
    for curve in accelerating:
        i = index[id(curve)]
        check(f"exp {curve.experiment} sample {curve.sample} keeps its lag branch",
              auto[i].v0 == never[i].v0,
              f"auto {auto[i].v0:.3e} vs unconstrained {never[i].v0:.3e}")


def test_burst_on_real_curves():
    print("the burst / lag form on real curves")
    curves, _ = build_curves()
    # Experiment 26 is the only true replicate set in the enzyme-free data:
    # four cuvettes at identical conditions whose initial rates agree to 9%.
    # Every one of them is sound, which is what makes it the right test -- any
    # misbehaviour of the burst form here cannot be blamed on the data.
    block = [c for c in curves if c.experiment == 26]
    check("experiment 26 has its four replicate cuvettes", len(block) == 4,
          f"got {len(block)}")
    fits = [fit_burst(c.times, c.absorbance) for c in block]
    check("the fits describe these curves inside their own noise",
          all(f.rms < c.noise for f, c in zip(fits, block)),
          f"rms/noise {[round(f.rms / c.noise, 2) for f, c in zip(fits, block)]}")
    check("tau is unresolved on most of even these clean replicates",
          sum(f.resolved for f in fits) <= 1,
          f"resolved: {[f.resolved for f in fits]}")
    # The point of the whole diagnostic: fitting well is not the same as being
    # identified.
    #
    # This was pinned on the SIGN of v0 until 2026-09-01 -- one of these four
    # sound, mutually agreeing cuvettes returned a negative initial rate from a
    # collapsed tau. On the instrument's own readings none of them does, so
    # that much was partly an artefact of the export's 0.001 AU rounding. The
    # underlying failure is unchanged and shows up in the PROFILE: every one of
    # these four fits inside its own noise, and every one still fails to bound
    # v0, three of them with an interval that spans zero.
    from summary_kinetics import fit_burst_bounded
    profiles = [fit_burst_bounded(c.times, c.absorbance, constrain=False)
                for c in block]
    check("no v0 is bounded on even these clean replicates",
          not any(p.bounded for p in profiles),
          f"half-widths {[round(p.half_width, 2) for p in profiles]}")
    check("and most profiles reach below zero",
          sum(p.v0_low < 0 for p in profiles) >= 3,
          f"lows {[f'{p.v0_low:.1e}' for p in profiles]}")
    check("...on curves whose straight-line rate is firmly positive",
          all(initial_rate(c.times, c.absorbance, INITIAL_WINDOW)[0] > 0
              for c in block))

    survivors = [c for c in curves if c.experiment == 25]
    check("experiment 25 keeps only its two sound cuvettes after exclusion",
          sorted(c.sample for c in survivors) == [1, 3],
          f"got {sorted(c.sample for c in survivors)}")


def test_recovery():
    print("exponent recovery from synthetic curves")
    frame = synthetic_frame(order_s=0.30, order_hoo=1.10, order_buf=0.0)
    fitted = regress(frame, "v0", None, per_experiment=True)
    for name, truth in (("log[S]", 0.30), ("log[hoo]", 1.10), ("log[buf]", 0.0)):
        value, stderr = fitted.coefficient(name)
        check(f"recovers {name} = {truth:+.2f}", abs(value - truth) < 3 * stderr,
              f"got {value:+.3f} +/- {stderr:.3f}")

    pooled = regress(frame, "v0", None, per_experiment=False)
    check("ignoring the offsets inflates the residual",
          pooled.residual_scatter > fitted.residual_scatter,
          f"pooled x{pooled.residual_scatter:.2f} vs "
          f"offsets x{fitted.residual_scatter:.2f}")
    check("the offset model is preferred by AIC", fitted.aic < pooled.aic,
          f"{fitted.aic:.1f} vs {pooled.aic:.1f}")


def test_saturation_recovery():
    print("saturation recovery")
    rng = np.random.default_rng(3)
    km = 0.5
    rows = []
    for experiment in range(1, 13):
        offset = rng.normal(-4.0, 0.5)
        for s0 in np.logspace(-1.2, 1.4, 6):
            v = 10 ** (offset + rng.normal(0.0, 0.045)) * (s0 / (km + s0))
            rows.append(dict(experiment=experiment, sample=len(rows), s0=float(s0),
                             h2o2=82.5, hoo=3e-3, buf=75.0, e0=0.0, v0=float(v),
                             v0_stderr=0.0, window_rms=0.0, noise=6e-4,
                             slope_ratio=1.0, points=200, duration=1000.0,
                             amplitude=0.05, conversion=0.005))
    frame = pd.DataFrame(rows)
    saturating = profile_km(frame, "v0", per_experiment=True)
    low, high = saturating.km_interval
    # The saturation term carries a free exponent, so Km and that exponent
    # trade off and the point estimate comes back biased low. A factor of two
    # is the honest precision, and matches what the real block delivers
    # (interval 0.001-7.2 mM).
    check("recovers Km to within a factor of two",
          0.5 * km <= saturating.km <= 2.0 * km, f"got {saturating.km:.3f} vs {km}")
    check("the profile interval is resolved, not the whole grid",
          np.isfinite(low) and high < 100.0 and high / low < 20.0,
          f"[{low:.3f}, {high:.3f}]")
    check("the profile interval contains the estimate", low <= saturating.km <= high)
    power = regress(frame, "v0", None, per_experiment=True)
    check("saturation beats a power law when the truth saturates",
          saturating.aic < power.aic, f"{saturating.aic:.1f} vs {power.aic:.1f}")

    # ... and does NOT win when the truth is a clean power law
    straight = synthetic_frame(order_s=0.30, sigma=0.02, seed=5)
    check("saturation does not win on data that has none",
          profile_km(straight, "v0", True).aic >= regress(straight, "v0", None, True).aic - 2.0)


def test_absorption_diagnostic():
    print("collinearity with the per-experiment offsets")
    # [HOO-] constant inside each experiment but different between them: an
    # offset model cannot identify it, and must say so.
    rows = []
    for experiment in range(1, 11):
        hoo = 10 ** (-2.0 - 0.08 * experiment)
        for sample, s0 in enumerate(np.logspace(-1, 1, 4), start=1):
            rows.append(dict(experiment=experiment, sample=sample, s0=float(s0),
                             h2o2=82.5, hoo=float(hoo), buf=75.0, e0=0.0,
                             v0=1e-4 * float(s0) ** 0.3, v0_stderr=0.0,
                             window_rms=0.0, noise=6e-4, slope_ratio=1.0,
                             points=200, duration=1000.0, amplitude=0.05,
                             conversion=0.005))
    frame = pd.DataFrame(rows)
    surviving = within_experiment_variation(frame)
    check("an axis constant within experiments keeps no within-variance",
          surviving["hoo"] < 1e-9, f"got {surviving['hoo']:.3g}")
    check("an axis that varies within experiments keeps its variance",
          surviving["s0"] > 0.99, f"got {surviving['s0']:.3g}")
    check("the absorbed axis is flagged", "hoo" in absorbed_axes(frame))
    check("the identifiable axis is not flagged", "s0" not in absorbed_axes(frame))
    check("regress reports the absorption it is suffering",
          "hoo" in regress(frame, "v0", None, per_experiment=True).absorbed)
    check("the pooled fit reports no absorption, having no offsets",
          regress(frame, "v0", None, per_experiment=False).absorbed == ())
    check("the floor is a fraction, not a count", 0.0 < ABSORPTION_FLOOR < 1.0)


def test_outliers():
    print("dead cuvettes")
    frame = synthetic_frame(seed=11)
    clean = regress(frame, "v0", None, per_experiment=True)
    killed = frame.copy()
    # A real dead cuvette is slow AND barely clears the noise; the rule requires
    # both, so the synthetic one has to have both.
    killed.loc[killed.index[0], "v0"] = float(killed.v0.iloc[0]) / 200.0
    killed.loc[killed.index[0], "amplitude"] = 0.003
    dead = (int(killed["experiment"].iloc[0]), int(killed["sample"].iloc[0]))
    check("a dead cuvette is flagged against its own experiment",
          dead in {(e, s) for e, s, _, _ in experiment_outliers(killed)})
    check("healthy curves are not flagged",
          len(experiment_outliers(frame)) == 0,
          f"flagged {len(experiment_outliers(frame))}")
    spoiled = regress(killed, "v0", None, per_experiment=True)
    check("one dead cuvette measurably wrecks the fit",
          spoiled.sse > 3 * clean.sse,
          f"{spoiled.sse:.4f} vs {clean.sse:.4f}")
    check("the outlier factor is above 1", OUTLIER_FACTOR > 1.0)
    check("a slow but clean curve is not called dead",
          not experiment_outliers(
              frame.assign(v0=np.where(frame.index == frame.index[0],
                                       frame.v0 / 200.0, frame.v0),
                           amplitude=1.0, noise=1e-4)),
          "a high-SNR curve must not be flagged however slow it is")


def test_real_block():
    print("the real enzyme-free block")
    curves, _ = build_curves()
    buffers = buffer_concentrations()
    block = [c for c in curves if c.conditions.e0 == 0
             and c.group == ("4OMe-BnOH", 40.0, "Phosphate")]
    # Sourced from fit_dataset rather than typed, so the two cannot drift apart:
    # the block lost 20 curves on 2026-08-31 when exps 32 and 34-37 were ruled
    # catalysed, and a hardcoded count would have hidden that instead of failing.
    expected = sum(len(g) for k, g in group_curves(curves, enzyme_free=True).items()
                   if k == ("4OMe-BnOH", 40.0, "Phosphate"))
    check("the block is the curves fit_dataset reports",
          len(block) == expected, f"got {len(block)}, fit_dataset says {expected}")
    frame = to_frame(summarise(block, buffers))
    quality = (frame.window_rms / frame.noise).replace([np.inf, -np.inf], np.nan)
    # A ratio near or below 1 means a straight line describes the window to
    # within the readings' own scatter. The lower bound is loose because the
    # point-noise estimate is conservative on the short, coarsely sampled runs
    # that dominate this block; the claim being tested is that the window is
    # linear, and only the upper bound can falsify it.
    check("a straight line over 20% fits to about the point noise",
          0.5 <= np.nanmedian(quality) <= 1.3, f"got {np.nanmedian(quality):.2f}")
    full = to_frame(summarise(block, buffers, 1.00))
    full_quality = (full.window_rms / full.noise).replace([np.inf, -np.inf], np.nan)
    check("a straight line over the whole run does not",
          np.nanmedian(full_quality) > 4.0, f"got {np.nanmedian(full_quality):.2f}")

    sets = replicate_scatter(frame)
    check("experiment 26 is found as a replicate set",
          any(s["experiment"] == 26 and s["n"] == 4 for s in sets))
    check("the replicate floor is near the documented default",
          any(abs(s["rsd"] - REPLICATE_RSD) < 0.05 for s in sets),
          f"got {[round(s['rsd'], 3) for s in sets]}")
    check("h2o2 is constant, so it is not a design axis",
          frame.h2o2.nunique() == 1 and "h2o2" not in DESIGN_AXES)
    check("most curves decelerate", (frame.slope_ratio < 0.9).sum() > len(frame) / 2,
          f"{(frame.slope_ratio < 0.9).sum()}/{len(frame)}")

    flagged = {(e, s) for e, s, _, _ in experiment_outliers(frame)}
    kept = frame[~frame.set_index(["experiment", "sample"]).index.isin(flagged)]
    offsets = regress(kept, "v0", None, per_experiment=True)
    pooled = regress(kept, "v0", None, per_experiment=False)
    check("per-experiment offsets cut the residual sum of squares several-fold",
          pooled.sse > 4 * offsets.sse, f"{pooled.sse:.3f} vs {offsets.sse:.3f}")
    # Until 2026-08-31 this asserted the opposite -- that [S] was absorbed by the
    # per-experiment offsets, at 6.4% surviving variance. The 20 curves removed
    # that day were buffer titrations at a FIXED substrate, so they contributed
    # offsets and no substrate contrast; dropping them roughly doubled what is
    # left to fit with, and the substrate order became interpretable.
    survives = within_experiment_variation(kept)["s0"]
    check("the substrate order now rests on real within-experiment contrast",
          survives > ABSORPTION_FLOOR, f"got {survives:.3f}, floor {ABSORPTION_FLOOR}")

    # Also inverted on 2026-08-31, and for the same reason. With the catalysed
    # curves out, the order stops depending on where the window is drawn: the
    # spread across 10-100% is now far smaller than the fit's own standard error,
    # so the window is no longer a systematic worth carrying.
    rows = window_sensitivity(block, buffers, drop_outliers=True)
    orders = [r["order"] for r in rows]
    spread = max(orders) - min(orders)
    check("the substrate order no longer depends on the window",
          spread < min(r["stderr"] for r in rows),
          f"{min(orders):+.2f} to {max(orders):+.2f}, spread {spread:.2f}")
    check("wider windows fit the curve worse",
          rows[-1]["rms_over_noise"] > rows[0]["rms_over_noise"])


def test_two_phase_and_selection():
    """
    The two-phase form recovers a rise-then-fall, and is refused when it is not one.

    The forms are exactly nested -- B2 = 0 is the one-phase form -- so the
    selection is an F test and the thing that must hold is BOTH directions: the
    second phase is taken when it is there and refused when it is not. A
    selector that only ever accepts is not a selector.
    """
    print("\ntwo-phase form and nested selection")
    times = np.arange(0, 130) * 50.0
    v_ss, B1, tau1, B2, tau2 = 2e-5, 6e-3, 300.0, -1.2e-2, 4000.0
    generator = np.random.default_rng(0)

    truth = (0.001 + v_ss * times - B1 * (1 - np.exp(-times / tau1))
             - B2 * (1 - np.exp(-times / tau2)))
    fitted = fit_progress(times, truth + generator.normal(0, 2e-4, len(times)))
    check("a rise-then-fall curve selects two phases", fitted.phases == 2,
          f"phases {fitted.phases}, F {fitted.f_statistic:.1f}")
    check("and recovers the steady rate", abs(fitted.two.v_ss - v_ss) < 0.1 * v_ss,
          f"{fitted.two.v_ss:.3e} against {v_ss:.3e}")
    check("and the fast time constant", abs(fitted.two.tau1 - tau1) < 0.3 * tau1,
          f"{fitted.two.tau1:.0f} against {tau1:.0f}")
    check("and names the shape", fitted.two.kind == "lag then fall",
          fitted.two.kind)
    peak, when = fitted.peak_rate
    check("the peak is interior and above the steady rate",
          peak > v_ss and 0 < when < times[-1], f"{peak:.3e} at {when:.0f} s")

    lag = 0.001 + v_ss * times - B1 * (1 - np.exp(-times / tau1))
    refused = fit_progress(times, lag + generator.normal(0, 2e-4, len(times)))
    check("a pure lag is refused the second phase", refused.phases == 1,
          f"phases {refused.phases}, F {refused.f_statistic:.1f}")
    check("and its peak rate is the one-phase asymptote",
          abs(refused.peak_rate[0] - refused.one.v_ss) < 1e-12)

    straight = 0.001 + v_ss * times
    line = fit_progress(times, straight + generator.normal(0, 2e-4, len(times)))
    check("a straight line is refused too", line.phases == 1,
          f"phases {line.phases}, F {line.f_statistic:.1f}")

    # The nesting itself: at B2 = 0 the two-phase form IS the one-phase form,
    # so it can never fit worse. A negative F would mean the grid search had
    # missed its own special case.
    for values in (truth, lag, straight):
        noisy = values + generator.normal(0, 2e-4, len(times))
        both = fit_progress(times, noisy)
        check("the two-phase fit never does worse than the one-phase",
              both.two.sse <= both.one.sse * 1.0001,
              f"{both.two.sse:.3e} against {both.one.sse:.3e}")

    # The fast route -- precomputed normal equations rather than an n-row
    # least squares per node -- must give the same answer. Checked on a NOISY
    # curve: normal equations square the condition number, so on a noiseless
    # one the residual is pure rounding (about 1e-8 here) and the comparison
    # measures floating point rather than the algebra.
    noisy = truth + generator.normal(0, 2e-4, len(times))
    exact = fit_two_phase(times, noisy)
    design = _two_phase_design(times - times[0], exact.tau1, exact.tau2)
    beta, *_ = np.linalg.lstsq(design, noisy, rcond=None)
    residual = noisy - design @ beta
    check("the precomputed normal equations match a direct least squares",
          abs(exact.sse - float(residual @ residual))
          <= 1e-8 * abs(exact.sse),
          f"{exact.sse:.9e} against {float(residual @ residual):.9e}")


if __name__ == "__main__":
    test_line_slope()
    test_line_intercept_and_floor()
    test_window()
    test_burst_fit()
    test_burst_on_real_curves()
    test_lag_branch_is_gated_per_curve()
    test_bounded_v0_does_not_mean_a_resolved_shape()
    test_recovery()
    test_saturation_recovery()
    test_absorption_diagnostic()
    test_outliers()
    test_real_block()
    test_two_phase_and_selection()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
