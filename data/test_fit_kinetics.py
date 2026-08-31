"""
Tests for fit_dataset.py and fit_kinetics.py.

The test that matters here is parameter recovery: synthetic curves are
generated from known rate constants at the real experimental conditions, noise
is added at the level the instrument actually delivers, and the fitter has to
find the constants back. A fitter that has never been shown to recover a known
answer is not evidence about the mechanism, only about the optimiser.

The rest pin the row selection against the counts the notebook's
clean_experiment_dataframe produces, so the two selections cannot drift apart.

    python data/test_fit_kinetics.py
"""
import sys

import numpy as np
import pandas as pd

from fit_dataset import (BASELINE_POINTS, PRIMARY_SCOPE, PRIMARY_SCOPE_BLOCK,
                         QUANTISATION_SIGMA, Curve, build_curves, curve_noise,
                         group_curves, in_scope, select_fittable)
from fit_kinetics import (BOUNDS, FAILURE_RESIDUAL, INITIAL, STAGE_ONE,
                          STAGE_TWO, fit_group, residuals, sequential_fit)
from kinetic_model import Conditions, RateConstants, observable, simulate

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


# What the notebook's clean_experiment_dataframe produces, and what README.md
# reports. Pinned so a change to either selection shows up as a failure here
# rather than as a quietly different fit.
EXPECTED_ROWS = 402
EXPECTED_EXPERIMENTS = 88

# The fitting scope, pinned. See fit_dataset.PRIMARY_SCOPE for why these runs
# and not others.
EXPECTED_SCOPE_CURVES = 119
EXPECTED_SCOPE_EXPERIMENTS = 17

# What "well designed" means here, in numbers rather than in a folder name: a
# run whose own cuvettes vary BOTH the substrate and the peroxide, so both
# orders are measurable inside the run rather than across per-experiment
# offsets. A two-fold spread is the bar, and the separation it draws is not
# marginal -- inside the scope the smallest ladders are 40x in substrate and
# 6.9x in peroxide, while outside it every single run holds one of the two
# axes exactly constant. The bar could sit anywhere between 1x and 6.9x
# without changing which runs qualify.
LADDER_MINIMUM = np.log(2.0)
SCOPE_SUBSTRATE_LADDER = np.log(39.9)   # the weakest substrate ladder in scope is 40.0x
SCOPE_PEROXIDE_LADDER = np.log(6.8)     # the weakest peroxide ladder in scope


def test_selection():
    print("\nrow selection matches clean_experiment_dataframe")
    data = pd.read_csv("data/experiment_data.csv")
    selected, report = select_fittable(data)
    check(f"{EXPECTED_ROWS} rows survive", len(selected) == EXPECTED_ROWS,
          f"got {len(selected)}")
    check(f"{EXPECTED_EXPERIMENTS} experiments survive",
          selected.experiment.nunique() == EXPECTED_EXPERIMENTS,
          f"got {selected.experiment.nunique()}")
    check("every rule removed something (none is dead)",
          all(report[k] > 0 for k in
              ("excluded_experiments", "excluded_samples", "excluded_buffers")),
          str(report))
    check("no carbonate survives", not (selected.buffer == "Carbonate").any())
    check("selection is idempotent",
          len(select_fittable(selected)[0]) == len(selected))


def test_curves():
    print("\ncurve assembly")
    curves, report = build_curves()
    check("every fittable row produced a curve",
          report["curves"] == EXPECTED_ROWS,
          f"{report['curves']} curves from {report['rows_out']} rows, "
          f"dropped {report['dropped']}")
    check("nothing was dropped silently",
          not any(report["dropped"].values()), str(report["dropped"]))

    check("time starts at zero for every curve",
          all(curve.times[0] == 0.0 for curve in curves))
    check("time is ascending for every curve",
          all(np.all(np.diff(curve.times) > 0) for curve in curves))
    check("the baseline puts every curve near zero at t = 0",
          max(abs(curve.absorbance[0]) for curve in curves) < 0.05,
          f"worst {max(abs(c.absorbance[0]) for c in curves):.4f} AU")
    check("noise never falls below the quantisation floor",
          all(curve.noise >= QUANTISATION_SIGMA - 1e-15 for curve in curves))
    check("every curve carries a positive extinction coefficient",
          all(curve.epsilon > 0 for curve in curves))
    check("[HOO-] is positive everywhere",
          all(curve.conditions.hoo > 0 for curve in curves))

    data = pd.read_csv("data/experiment_data.csv")
    lookup = {(int(r["experiment"]), int(r["sample"])): r
              for r in data.to_dict("records")}
    mismatched = [c for c in curves
                  if not np.isclose(c.conditions.s0, lookup[(c.experiment, c.sample)]["[sub]"])
                  or not np.isclose(c.conditions.e0, lookup[(c.experiment, c.sample)]["[enz]"])]
    check("each curve's conditions come from its own dataset row",
          not mismatched, f"{len(mismatched)} mismatched")

    grouped = group_curves(curves)
    check("no group mixes substrates, temperatures or buffers",
          all(len({(c.substrate, c.temperature, c.buffer) for c in block}) == 1
              for block in grouped.values()))
    free = group_curves(curves, enzyme_free=True)
    catalysed = group_curves(curves, enzyme_free=False)
    check("the two stages partition the curves",
          sum(map(len, free.values())) + sum(map(len, catalysed.values())) == len(curves))


def test_noise_estimator():
    print("\nnoise estimator")
    generator = np.random.default_rng(0)
    trend = np.linspace(0.0, 0.5, 200)
    known = 0.004
    noisy = trend + generator.normal(0.0, known, 200)
    estimate = curve_noise(noisy)
    check("recovers a known sigma through a strong linear trend",
          abs(estimate - known) / known < 0.15, f"{estimate:.5f} vs {known}")
    check("a perfectly smooth curve floors at the quantisation sigma",
          curve_noise(trend) == QUANTISATION_SIGMA)
    check("a short curve floors rather than dividing by nothing",
          curve_noise([0.1, 0.2, 0.3]) == QUANTISATION_SIGMA)


def _synthetic(constants, conditions_list, sigma, seed=0, points=60, duration=3000.0):
    """Curves generated from known constants, with instrument-level noise."""
    generator = np.random.default_rng(seed)
    times = np.linspace(0.0, duration, points)
    curves = []
    for index, conditions in enumerate(conditions_list):
        signal = observable(constants, conditions, times)
        epsilon = 1.23
        measured = epsilon * signal + generator.normal(0.0, sigma, points)
        measured -= np.median(measured[:BASELINE_POINTS])
        curves.append(Curve(
            experiment=900 + index, sample=1, substrate="BnOH", buffer="Phosphate",
            pH=8.0, temperature=25.0, epsilon=epsilon, times=times,
            absorbance=measured, baseline=0.0, noise=max(sigma, QUANTISATION_SIGMA),
            conditions=conditions,
        ))
    return curves


def test_parameter_recovery():
    """
    The fitter's own validation gate. Note the truth uses r = 2.0: at r <= 1 the
    model cannot make a lag (see kinetic_model.test_acceleration_requires_r_above_one),
    so a recovery test at r <= 1 would only ever exercise the flat corner of the
    parameter space.
    """
    print("\nparameter recovery from synthetic curves")
    truth = RateConstants(k_can=5.0, k3=2e-2, k0=2e-8, k5=0.0, k6=0.0, r=2.0)
    conditions = [Conditions(s0=s, h2o2=100.0, e0=0.0, hoo=2e-3)
                  for s in (0.5, 1.0, 2.5, 5.0, 8.0, 12.0)]

    noiseless = _synthetic(truth, conditions, sigma=0.0)
    result = fit_group(noiseless, STAGE_ONE, restarts=1, seed=1)
    for name in STAGE_ONE:
        found, expected = getattr(result.constants, name), getattr(truth, name)
        if name == "r":
            ok = abs(found - expected) < 0.15
            detail = f"{found:.3f} vs {expected:.3f}"
        else:
            ok = abs(np.log10(found) - np.log10(expected)) < 0.20
            detail = f"{found:.3g} vs {expected:.3g} ({np.log10(found/expected):+.2f} decades)"
        check(f"noiseless: {name} recovered", ok, detail)

    # The noise level and signal size are matched to the real BnOH/25/phosphate
    # block, which runs 0.005-0.088 AU net against a median per-curve noise of
    # 0.0006 AU. An earlier version of this test used 0.002 AU noise on curves
    # topping out at 0.014 AU -- a signal-to-noise of about 1, some thirty times
    # worse than the instrument delivers -- and its failure said nothing about
    # the fitter.
    noisy = _synthetic(truth, conditions, sigma=0.0006, seed=7)
    result = fit_group(noisy, STAGE_ONE, restarts=1, seed=1)
    # k_can, k0 and r are recovered; k3 is not, and is excluded here rather than
    # given a loose tolerance, because test_k3_is_a_lower_bound_only explains
    # why. Across noise seeds 7, 11 and 23 the recovered values move by
    # -0.02/-0.10/+0.16 decades for k_can, 0.00/+0.01/+0.01 for k0 and
    # +0.03/+0.01/-0.13 for r, against +0.10/+1.06/-0.97 decades for k3.
    for name in ("k_can", "k0", "r"):
        found, expected = getattr(result.constants, name), getattr(truth, name)
        if name == "r":
            ok, detail = abs(found - expected) < 0.4, f"{found:.3f} vs {expected:.3f}"
        else:
            ok = abs(np.log10(found) - np.log10(expected)) < 0.4
            detail = f"{found:.3g} vs {expected:.3g} ({np.log10(found/expected):+.2f} decades)"
        check(f"at instrument noise: {name} recovered", ok, detail)
    check("the fit reports residuals near the noise it was given",
          0.3 < result.rms_sigma < 4.0, f"{result.rms_sigma:.2f} sigma")


def test_k3_is_a_lower_bound_only():
    """
    Why k3 is not in the recovery test above.

    Step 3 consumes the peracid that steps 1-2 make. Once k3 is large enough
    that PBA reaches quasi-steady state, v3 -> v_can, and since benzoate is
    produced at v3 the observable is then set by v_can alone -- k3 has dropped
    out of it. Raising k3 further only lowers the standing [PBA] in exact
    proportion, which nothing measures.

    So progress curves bound k3 from below and say nothing above that bound.
    The same is true of k5 in stage 2, for the same kind of reason, and it is
    why the real BnOH/25/phosphate block reports k_can and k3 anticorrelated at
    -0.979: near the threshold the two trade off.
    """
    print("\nk3 drops out of the observable above the quasi-steady-state threshold")
    conditions = Conditions(s0=8.0, h2o2=100.0, e0=0.0, hoo=2e-3)
    times = np.linspace(0.0, 3000.0, 60)
    reference = observable(RateConstants(k_can=5.0, k3=2e-2, k0=2e-8, r=2.0),
                           conditions, times)

    def deviation(k3):
        signal = observable(RateConstants(k_can=5.0, k3=k3, k0=2e-8, r=2.0),
                            conditions, times)
        return float(np.max(np.abs(signal - reference)) / reference[-1])

    high = deviation(2e1)      # 1000x the reference
    low = deviation(2e-4)      # 1/100th of it
    check("raising k3 a thousandfold barely moves the signal",
          high < 0.01, f"{high * 100:.2f}% change")
    check("lowering k3 a hundredfold moves it a great deal",
          low > 0.1, f"{low * 100:.1f}% change")
    check("so the curves bound k3 from below, not from above", low > 20 * high,
          f"low {low:.3g} vs high {high:.3g}")

    peracid = [simulate(RateConstants(k_can=5.0, k3=k3, k0=2e-8, r=2.0),
                        conditions, times)["PBA"][-1] for k3 in (2e-1, 2e0, 2e1)]
    check("above the threshold [PBA] falls in inverse proportion to k3",
          np.allclose(np.array(peracid[:-1]) / np.array(peracid[1:]), 10.0, rtol=0.02),
          f"ratios {np.array(peracid[:-1]) / np.array(peracid[1:])}")


def test_stage_two_recovery():
    print("\nstage 2 recovery with stage 1 frozen")
    truth = RateConstants(k_can=5.0, k3=2e-2, k0=2e-8, k5=4.0, k6=30.0, r=2.0)
    background = [Conditions(s0=s, h2o2=100.0, e0=0.0, hoo=2e-3)
                  for s in (0.5, 2.0, 5.0, 10.0)]
    catalysed = [Conditions(s0=s, h2o2=100.0, e0=e, hoo=2e-3)
                 for s in (1.0, 4.0, 8.0) for e in (0.05, 0.15, 0.28)]
    curves = _synthetic(truth, background, sigma=0.0) + \
             _synthetic(truth, catalysed, sigma=0.0, seed=3)

    first, second = sequential_fit(curves, restarts=3, seed=2)
    check("stage 1 reports k5 and k6 as zero rather than as their seed values",
          first.constants.k5 == 0.0 and first.constants.k6 == 0.0,
          f"k5 = {first.constants.k5}, k6 = {first.constants.k6}")
    check("stage 2 carried stage 1's constants forward unchanged",
          second.constants.k_can == first.constants.k_can
          and second.constants.r == first.constants.r)

    found = second.constants.k6
    check("stage 2: k6 recovered",
          abs(np.log10(found) - np.log10(truth.k6)) < 0.35,
          f"{found:.3g} vs {truth.k6:.3g} ({np.log10(found / truth.k6):+.2f} decades)")

    # k5 is a different story, and the reason is chemistry rather than
    # optimisation: once the seed is fast enough that starting the loop is not
    # what limits the observable, making it faster still changes almost nothing.
    # The cost surface is therefore steep below the true k5 and nearly flat
    # above it, so progress curves give a LOWER BOUND on the seed constant, not
    # a value. This is asserted on the cost surface directly rather than on
    # where the optimiser lands, because the flat direction means it can stop
    # anywhere along it.
    base = truth.replace(k5=truth.k5, k6=truth.k6)
    catalysed_curves = [c for c in curves if c.conditions.e0 > 0]

    def cost_at(factor):
        vector = np.array([np.log10(truth.k5 * factor), np.log10(truth.k6)])
        return float(0.5 * np.sum(residuals(vector, STAGE_TWO, base,
                                            catalysed_curves) ** 2))

    below, at, above = cost_at(0.1), cost_at(1.0), cost_at(100.0)
    check("stage 2: reducing k5 tenfold is strongly penalised",
          below > 100.0, f"cost {below:.3g}")
    check("stage 2: raising k5 a hundredfold is barely penalised",
          above < below / 100.0, f"cost {above:.3g} vs {below:.3g} below")
    check("stage 2: so k5 is identified only as a lower bound",
          second.constants.k5 > truth.k5 / 3.0,
          f"fit gave {second.constants.k5:.3g}, truth {truth.k5:.3g}")
    print(f"        k5 cost surface: x0.1 -> {below:.3g}, x1 -> {at:.3g}, "
          f"x100 -> {above:.3g}")


def test_residual_machinery():
    print("\nresidual machinery")
    truth = RateConstants(k_can=5.0, k3=2e-2, k0=3e-9, r=2.0)
    conditions = [Conditions(s0=s, h2o2=100.0, e0=0.0, hoo=2e-3) for s in (1.0, 5.0)]
    curves = _synthetic(truth, conditions, sigma=0.0)

    at_truth = residuals(np.array([np.log10(truth.k_can), np.log10(truth.k3),
                                   np.log10(truth.k0), truth.r]),
                         STAGE_ONE, RateConstants(), curves)
    check("residuals vanish at the generating parameters",
          float(np.max(np.abs(at_truth))) < 1e-6,
          f"max {np.max(np.abs(at_truth)):.2e}")

    absurd = residuals(np.array([8.0, 6.0, -2.0, 5.0]), STAGE_ONE,
                       RateConstants(), curves)
    check("an unintegrable parameter set gives finite penalty residuals",
          np.all(np.isfinite(absurd)))
    check("the penalty is large enough to be rejected",
          float(np.max(np.abs(absurd))) >= 1.0)
    check("residual length always matches the data",
          len(absurd) == sum(len(c) for c in curves))

    short, long_ = curves[0], curves[1]
    long_.times = np.linspace(0, 3000, 60)
    per_curve = residuals(np.array([0.0, -2.0, -9.0, 0.3]), STAGE_ONE,
                          RateConstants(), [short])
    per_point = residuals(np.array([0.0, -2.0, -9.0, 0.3]), STAGE_ONE,
                          RateConstants(), [short], weighting="point")
    check("per-curve weighting is per-point weighting divided by sqrt(n)",
          np.allclose(per_curve * np.sqrt(len(short)), per_point))



def _ladder_spread(curves, attribute):
    """Log-range of `attribute` across one run's own cuvettes."""
    values = np.log([max(getattr(c.conditions, attribute), 1e-12) for c in curves])
    return values.max() - values.min()


def _two_axis_runs(curves):
    """Experiments carrying both an internal substrate and an internal peroxide ladder."""
    by_experiment = {}
    for curve in curves:
        by_experiment.setdefault(curve.experiment, []).append(curve)
    return {experiment for experiment, group in by_experiment.items()
            if _ladder_spread(group, "s0") >= LADDER_MINIMUM
            and _ladder_spread(group, "h2o2") >= LADDER_MINIMUM}


def test_scope():
    """
    Pins the fitting scope, and re-derives it from the designs rather than
    trusting the list.

    The last check is the one that matters: if a run outside the scope ever
    turns out to carry the same two-axis design -- because an exclusion was
    lifted, a sheet was re-read, or an unexported run was recovered -- this
    fails and forces a decision, instead of leaving it silently out of the fit.
    """
    print("\nfitting scope")
    curves, _ = build_curves()
    scoped = in_scope(curves)

    check(f"{EXPECTED_SCOPE_CURVES} curves in scope",
          len(scoped) == EXPECTED_SCOPE_CURVES, f"got {len(scoped)}")
    check(f"{EXPECTED_SCOPE_EXPERIMENTS} experiments in scope",
          len({c.experiment for c in scoped}) == EXPECTED_SCOPE_EXPERIMENTS,
          f"got {len({c.experiment for c in scoped})}")
    check("every scoped experiment survived the exclusions",
          {c.experiment for c in scoped} == set(PRIMARY_SCOPE),
          f"missing {sorted(set(PRIMARY_SCOPE) - {c.experiment for c in scoped})}")
    check("the scope is one block, so rate constants may be pooled across it",
          {c.group for c in scoped} == {PRIMARY_SCOPE_BLOCK})
    check("every scoped run carries the full seven cuvettes",
          all(len([c for c in scoped if c.experiment == e]) == 7 for e in PRIMARY_SCOPE))

    # The reason for the scope, restated as measurements.
    check("every scoped run measures its substrate order internally, over >= 39.9x",
          all(_ladder_spread([c for c in scoped if c.experiment == e], "s0")
              >= SCOPE_SUBSTRATE_LADDER for e in PRIMARY_SCOPE))
    check("every scoped run measures its peroxide order internally, over >= 6.8x",
          all(_ladder_spread([c for c in scoped if c.experiment == e], "h2o2")
              >= SCOPE_PEROXIDE_LADDER for e in PRIMARY_SCOPE))
    pH_values = {c.pH for c in scoped}
    check("the scope spans at least three pH units in one block",
          max(pH_values) - min(pH_values) >= 3.0,
          f"{min(pH_values)}-{max(pH_values)}")
    hoo = [c.conditions.hoo for c in scoped]
    check("[HOO-] spans at least three decades across the scope",
          np.log10(max(hoo) / max(min(hoo), 1e-12)) >= 3.0)

    # The guard: nothing outside the scope has the design the scope was chosen for.
    outside = _two_axis_runs(curves) - set(PRIMARY_SCOPE)
    check("no run outside the scope varies both axes at all",
          not outside, f"unscoped two-axis runs: {sorted(outside)}")


def test_configuration():
    print("\nconfiguration")
    check("every fitted parameter has bounds",
          set(STAGE_ONE + STAGE_TWO) <= set(BOUNDS))
    check("every fitted parameter has a starting value",
          set(STAGE_ONE + STAGE_TWO) <= set(INITIAL))
    check("every starting value lies inside its bounds",
          all(BOUNDS[n][0] <= INITIAL[n] <= BOUNDS[n][1] for n in INITIAL))
    check("r is allowed above 1, so the fit can reach the only regime that lags",
          BOUNDS["r"][1] > 1.0)
    check("r is not allowed negative",
          BOUNDS["r"][0] >= 0.0)
    check("the two stages are disjoint",
          not set(STAGE_ONE) & set(STAGE_TWO))
    check("the failure penalty is finite",
          np.isfinite(FAILURE_RESIDUAL))


if __name__ == "__main__":
    test_selection()
    test_scope()
    test_curves()
    test_noise_estimator()
    test_parameter_recovery()
    test_k3_is_a_lower_bound_only()
    test_stage_two_recovery()
    test_residual_machinery()
    test_configuration()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
