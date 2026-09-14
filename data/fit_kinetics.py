"""
Fits the reduced mechanism of MECHANISM.md to the progress curves.

The strategy is the sequential one the reduction was done to enable. The model
is exactly linear in E0, so the enzyme-free controls determine the background
constants on their own, and freezing those leaves only the two catalysed
constants for the catalysed runs:

    stage 1   E0 = 0   ->  k_can, k3, k0, r      (2 ODEs' worth of chemistry)
    stage 2   E0 > 0   ->  k5, k6                (stage 1 frozen)

Both stages pool curves only within one (substrate, temperature, buffer) cell.
Temperature moves every rate constant through Arrhenius, the two substrates are
different molecules, and MECHANISM.md's buffer section argues at length that the
four buffers are chemically different reagents rather than four ways of setting
pH -- so pooling across any of those three would be fitting one constant to two
different quantities.

Residuals are taken in ABSORBANCE, not in mM. The conversion to concentration
divides by an extinction coefficient that differs sixfold between the two
substrates, which would rescale the instrument's roughly-constant noise along
with the signal and quietly weight one substrate six times the other.

    python data/fit_kinetics.py                       # the default block
    python data/fit_kinetics.py --substrate 4OMe-BnOH --temperature 40
    python data/fit_kinetics.py --list                # what blocks exist
    python data/fit_kinetics.py --save results.json
"""
import json
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from curve_metrics import peak_position
from fit_dataset import (BASELINE_POINTS, TWO_AXIS_BLOCK,
                         TWO_AXIS_GROUP, build_curves, group_curves,
                         in_block)
from kinetic_model import (LOG_PARAMETERS, PARAMETER_NAMES, Conditions,
                           RateConstants, increment, observable, pack, simulate,
                           unpack)
from summary_kinetics import TWO_PHASE_F

# Bounds in the optimiser's own coordinates: log10 for rate constants, linear
# for r. Wide enough not to shape the answer, tight enough that the integrator
# is not asked to do something absurd. A fitted value sitting ON a bound is
# reported as such -- for r especially, that is the result, not a detail.
#
# The extensions are bounded so that the OFF value is the FAR bound, and a
# fitted value AT that bound is "off" (see report and to_dict): k_sink, K4 and
# K_act are off at 0 (lower), km_s, k_act_r and km_s_background at infinity
# (upper). k_act_r and km_s are physical at small positive values, so their
# lower bound is just wide, not off.
BOUNDS = {
    "k_can": (-8.0, 8.0),
    "k3": (-10.0, 6.0),
    "k0": (-16.0, -2.0),
    "k5": (-8.0, 10.0),
    "k6": (-8.0, 10.0),
    "r": (0.0, 5.0),
    "k_sink": (-10.0, 2.0),
    "K4": (-8.0, 6.0),
    "km_s": (-2.0, 4.0),
    "k_act_r": (-8.0, 2.0),
    "K_act": (-8.0, 4.0),
    "km_s_background": (-2.0, 4.0),
}

# Starting points, in the same coordinates. The scaling behind k_can: the seed
# makes aldehyde at v0 = k0[H2O2][S], so over a run of length T the pool reaches
# ~v0*T, and steps 1-2 only matter once k_can*[HOO-]*(v0*T)^2 is comparable to
# v0 -- which for this dataset's numbers puts k_can near 1 mM^-2 s^-1.
INITIAL = {"k_can": 0.0, "k3": -2.0, "k0": -9.0, "k5": 0.0, "k6": 0.0, "r": 0.3,
           # extensions. k_sink from the product sink (~1e-4 s^-1, slowdown),
           # km_s from the per-run Michaelis fits (2-7 mM, saturation), k_act_r
           # from the activation clock (tau ~ hundreds of s), K4 and K_act
           # weakly determined.
           "k_sink": -4.0, "K4": -1.0, "km_s": 0.5, "k_act_r": -2.0,
           "K_act": -1.0, "km_s_background": 0.5}

# The extended stages' starting table. The original k5 = 1 is six decades above
# the real seed scale, and the extensions add k6, K4, km_s, k_act_r and K_act,
# whose priors come from the same fits: k6 from F5 (4.4e-3), K4 from Step 1
# (K ~ 0.04-0.05 /mM), km_s from `saturation.michaelis_by_element` (2-7 mM),
# k_act_r from the activation clock (tau ~ hundreds of s) and K_act from Step 2
# (K_shared 0.018 /mM). The old path never sees this table.
EXTENDED_INITIAL = {"k_can": 0.0, "k3": -2.0, "k0": -9.0, "k5": -6.0,
                    "k6": -2.0, "r": 0.3, "k_sink": -4.0, "K4": -1.3,
                    "km_s": 0.5, "k_act_r": -2.5, "K_act": -1.5,
                    "km_s_background": 0.5}

# The original six, which is what the base RateConstants is built from. The
# extensions must stay at their OFF defaults in `base` unless a stage actually
# frees them, or the "untouched" path would silently run the extended model.
BASE_PARAMETERS = ("k_can", "k3", "k0", "k5", "k6", "r")

STAGE_ONE = ("k_can", "k3", "k0", "r")
STAGE_TWO = ("k5", "k6")

# The extended splits (Stage 3.2). Stage 1 is all E0 = 0, so only the sink and
# a saturable BACKGROUND seed can be asked of it; the catalysed binding,
# saturation and activation are stage 2's.
STAGE_ONE_EXTENDED = STAGE_ONE + ("k_sink", "km_s_background")
STAGE_TWO_EXTENDED = STAGE_TWO + ("K4", "km_s", "k_act_r", "K_act")

# The nested ladder of Stage 3.3, one tuple pair per model. Each model's stage
# 1 is the background model and stage 2 adds the catalysed term; a term is
# "earned" by the weighted F test and is not carried into the next model if it
# is not (that is why M2 builds on M1 and not M1b). `M4` is what `--extended`
# selects.
MODEL_STAGES = {
    "M0": (STAGE_ONE, STAGE_TWO),
    "M1": (STAGE_ONE, STAGE_TWO + ("km_s",)),
    "M1b": (STAGE_ONE + ("km_s_background",), STAGE_TWO + ("km_s",)),
    "M2": (STAGE_ONE, STAGE_TWO + ("km_s", "K4")),
    "M3": (STAGE_ONE + ("k_sink",), STAGE_TWO + ("km_s", "K4")),
    "M4": (STAGE_ONE + ("k_sink",),
           STAGE_TWO + ("km_s", "K4", "k_act_r", "K_act")),
    # The EARNED-term model: M4 built on the terms that cleared the F bar and
    # dropping the ones that did not (K4 adds nothing on either block; k_sink
    # only moves stage 1's local minimum). It is what Stage 3.3's "do not carry
    # them into the next model" asks for, and it is the model the activation
    # question is actually asked of.
    "M4b": (STAGE_ONE + ("km_s_background",),
            STAGE_TWO + ("km_s", "k_act_r", "K_act")),
}

# What a failed integration costs. Large enough that the optimiser walks away,
# finite so it never poisons the Jacobian with a NaN.
FAILURE_RESIDUAL = 1e3


def baseline_like_data(model):
    """
    Puts a modelled curve on the same footing as the measurement.

    Public because anything that compares model to data -- the fit, and
    plot_fit.py -- must apply exactly this transformation. A plot that skipped it
    would show a disagreement the fit was never asked to remove.

    fit_dataset subtracts a baseline from each measured curve: the median of its
    first few readings. That is not the same as the curve's value at t = 0 --
    the reaction has already moved during those readings -- so a model that is
    exactly zero at t = 0 sits systematically above data that has had its own
    early points subtracted. Applying the identical operation to the model
    removes the bias. It is small, but it is a bias rather than noise: without
    it the residual at the true parameters does not vanish even on noiseless
    synthetic curves, which is how it was found.
    """
    if len(model) < 2:
        return model
    return model - np.median(model[:max(1, min(BASELINE_POINTS, len(model) // 10))])

# Integration tolerances used during fitting. Looser than kinetic_model's own
# defaults because a fit runs tens of thousands of integrations: at these
# settings a curve costs 17 ms instead of 33 ms, and the trajectory differs from
# a rtol=1e-10 reference by 5e-8 relative -- five orders of magnitude below the
# 0.001 AU the instrument records, so the loosening cannot move a fitted
# constant. Reporting and per-curve statistics use the tight defaults.
FIT_RTOL = 1e-6
FIT_ATOL = 1e-12

# And a much tighter evaluation cap than kinetic_model's own. A healthy curve
# integrates in a few thousand rhs calls; an optimiser probing log-space will
# propose sets that need millions, and at ~1 us per call the default 200,000
# turns each of those into a third of a second -- which is what turns a fit from
# minutes into hours. Abandoning them after 20,000 costs nothing: a parameter set
# that cannot be integrated to this accuracy is not a candidate answer.
FIT_MAX_EVALUATIONS = 20_000


@dataclass
class FitResult:
    constants: RateConstants
    free_names: tuple
    curves: list
    cost: float
    rms_absorbance: float
    rms_sigma: float               # in units of each curve's own noise
    at_bound: dict                 # name -> 'lower' | 'upper'
    standard_errors: dict
    correlation: np.ndarray
    condition_number: float
    per_curve: list = field(default_factory=list)
    success: bool = True
    message: str = ""
    observation: str = "design"


def _weights(curve, weighting):
    """
    Per-point residual denominator.

    'curve' (default) divides by sqrt(n) as well as by the curve's noise, so
    every curve contributes the same total weight. Points within a progress
    curve are strongly autocorrelated -- systematic model error dwarfs the
    reading noise -- so a 368-point 5-hour run does not carry 37 times the
    independent information of a 10-point one, and letting it carry 37 times the
    weight would let a handful of long runs decide the answer.

    'point' weights each reading equally, which is the right choice only if the
    residuals really are independent. Both are offered so the fit can be
    reported either way.
    """
    if weighting == "point":
        return curve.noise
    return curve.noise * np.sqrt(len(curve))


def observed_signal(curve, constants, observation="design", **kwargs):
    """
    The modelled signal the fit is compared against, for one curve.

    `observation="design"` matches the observable to the run's own design:
    a curve whose reference cuvette omitted the enzyme is a catalytic
    INCREMENT (`kinetic_model.increment`), everything else is the absolute
    signal. `observation="absolute"` reproduces the fitter's behaviour before
    R0.0 exactly, which is kept so the old M0-M4b saves can be reproduced and
    the two compared.

    Public because `plot_fit.py` must draw the same curve the fit was scored
    against -- a plot that used the other choice would show a disagreement the
    fit was never asked to remove.
    """
    if observation == "design" and curve.reference_omits == "enzyme":
        return increment(constants, curve.conditions, curve.times, **kwargs)
    return observable(constants, curve.conditions, curve.times, **kwargs)


def residuals(vector, free_names, base, curves, weighting="curve",
              observation="design"):
    """Stacked weighted residuals, model minus measurement, in absorbance."""
    constants = unpack(vector, free_names, base)
    stacked = []
    for curve in curves:
        denominator = _weights(curve, weighting)
        signal = observed_signal(curve, constants, observation,
                                 rtol=FIT_RTOL, atol=FIT_ATOL,
                                 max_evaluations=FIT_MAX_EVALUATIONS)
        if signal is None or not np.all(np.isfinite(signal)):
            stacked.append(np.full(len(curve), FAILURE_RESIDUAL))
            continue
        stacked.append((baseline_like_data(curve.epsilon * signal) - curve.absorbance)
                       / denominator)
    return np.concatenate(stacked) if stacked else np.zeros(0)


def _statistics(result, free_names, curves, weighting):
    """
    Standard errors, the correlation matrix and the Jacobian's condition
    number, all from J at the solution.

    The correlation matrix is the point of this function. With E0 varying
    across barely any of the dataset, k5 and k6 can be strongly correlated with
    each other, and a fit that reports only best-fit values would hide that.
    """
    jacobian = result.jac
    count = len(free_names)
    residual_count = jacobian.shape[0]
    degrees = max(residual_count - count, 1)
    variance = 2.0 * result.cost / degrees

    hessian = jacobian.T @ jacobian
    singular = np.linalg.svd(hessian, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else np.inf
    try:
        covariance = np.linalg.inv(hessian) * variance
        errors = np.sqrt(np.abs(np.diag(covariance)))
        scale = np.outer(errors, errors)
        correlation = np.divide(covariance, scale, out=np.zeros_like(covariance),
                                where=scale > 0)
    except np.linalg.LinAlgError:
        errors = np.full(count, np.nan)
        correlation = np.full((count, count), np.nan)
    return dict(zip(free_names, errors)), correlation, condition


def _per_curve(constants, curves, observation="design"):
    """
    One row per curve: how well it is fitted, and whether the model reproduces
    its SHAPE. `peak_data` and `peak_model` are where the steepest point sits as
    a fraction of the run -- the statistic MECHANISM.md used to falsify the
    aldehyde reading, so the fit is checked against it rather than only against
    a residual norm.
    """
    rows = []
    for curve in curves:
        signal = observed_signal(curve, constants, observation)
        if signal is None:
            rows.append({"experiment": curve.experiment, "sample": curve.sample,
                         "rms": np.nan, "sigma": np.nan,
                         "peak_data": np.nan, "peak_model": np.nan})
            continue
        model = baseline_like_data(curve.epsilon * signal)
        error = model - curve.absorbance
        rows.append({
            "experiment": curve.experiment,
            "sample": curve.sample,
            "rms": float(np.sqrt(np.mean(error ** 2))),
            "sigma": float(np.sqrt(np.mean(error ** 2)) / curve.noise),
            "net_data": float(curve.absorbance[-1]),
            "net_model": float(model[-1]),
            "peak_data": peak_position(curve.absorbance, curve.times),
            "peak_model": peak_position(model, curve.times),
        })
    return rows


# Starting points are SCREENED rather than guessed. The cost surface has long
# narrow valleys -- k_can/r come out ~0.999 anticorrelated on synthetic data and
# k_can/k3 ~0.979 on the real curves -- and a handful of blind starts slides
# along one and stops. That is not hypothetical: with blind starts, stage 2 on
# BnOH/25/phosphate returned a "converged" fit at cost 4.1e7 when cost 7.7e3 was
# available at k5, k6 -> 0, a factor of 5300. Sampling the whole box, ranking by
# one cheap residual evaluation each, and only then running the optimiser from
# the best few costs about 25 s and removes that failure mode.
SCREEN_SAMPLES = 64

# r is additionally swept on a fixed ladder because it is the decisive
# parameter, and because r = 1 is a boundary in the model's behaviour rather
# than just a number: below it the observable cannot produce a lag at all
# (see kinetic_model). Both sides deserve a start of their own.
R_SWEEP = (0.1, 0.5, 1.0, 2.0, 3.5)


def _guaranteed_points(free_names, start, lower, upper):
    """
    Starts that always run, whatever the screen thinks of them.

    The nominal start and the r ladder are structural: the nominal one is where
    the recovery tests converge exactly, and r = 1 is a boundary in the model's
    behaviour rather than just a number, so both sides of it deserve a start.
    Screening must never drop these -- an earlier version let it, and noiseless
    recovery went from exact to four decades out, because a low cost AT a point
    is a poor predictor of where an optimiser starting there ends up.
    """
    points = [np.asarray(start, dtype=float)]
    if "r" in free_names:
        index = free_names.index("r")
        for value in R_SWEEP:
            guess = np.asarray(start, dtype=float).copy()
            guess[index] = float(np.clip(value, lower[index], upper[index]))
            points.append(guess)
    return points


def _hypercube(free_names, lower, upper, seed):
    """
    A Latin-hypercube sample of the whole bounded box: one stratified draw per
    parameter, independently shuffled, so each parameter's range is covered
    evenly instead of clumping the way uniform sampling does in this many
    dimensions. This is what reaches corners a jittered start never does --
    small k5 and k6, in the case that motivated it.
    """
    generator = np.random.default_rng(seed)
    count = len(free_names)
    strata = (np.arange(SCREEN_SAMPLES)[:, None]
              + generator.random((SCREEN_SAMPLES, count))) / SCREEN_SAMPLES
    for column in range(count):
        generator.shuffle(strata[:, column])
    return [lower + row * (upper - lower) for row in strata]


def _screen(points, free_names, base, curves, weighting, keep, observation="design"):
    """
    Ranks candidate starting points by cost and returns the best `keep`.

    One residual evaluation each. Points the integrator cannot handle come back
    at the failure penalty and sort themselves to the bottom, so no special
    casing is needed.
    """
    scored = []
    for point in points:
        residual = residuals(point, free_names, base, curves, weighting, observation)
        scored.append((float(0.5 * np.sum(residual ** 2)), point))
    scored.sort(key=lambda pair: pair[0])
    return [point for _, point in scored[:keep]]


def fit_group(curves, free_names, base=None, weighting="curve", restarts=4,
              seed=0, initial=None, observation="design"):
    """
    Fits `free_names` to `curves`, holding everything else at `base`.

    Two kinds of start are used, and the best result over all of them is
    returned: the guaranteed ones (nominal, plus the r ladder when r is free),
    which always run, and the `restarts` best of a screened Latin-hypercube
    sample. Neither alone is enough -- see _guaranteed_points and
    SCREEN_SAMPLES, each of which records the failure that put it here.

    `initial` overrides the starting table (default `INITIAL`). The extended
    fit needs it because the original nominal k5 = 1 is six decades above the
    real seed scale (~1e-6): from there the catalysed fit cannot find the
    planted or the real basin, which is a starting problem and not a chemistry
    one.

    `observation` is "design" (the default) or "absolute"; see
    `observed_signal`. It is carried on the result and written by `to_dict`.
    """
    base = base or RateConstants(**{name: 10.0 ** INITIAL[name] if name in LOG_PARAMETERS
                                    else INITIAL[name] for name in BASE_PARAMETERS})
    table = INITIAL if initial is None else initial
    lower = np.array([BOUNDS[name][0] for name in free_names])
    upper = np.array([BOUNDS[name][1] for name in free_names])
    start = np.array([table[name] for name in free_names])

    starts = _guaranteed_points(free_names, start, lower, upper)
    starts += _screen(_hypercube(free_names, lower, upper, seed),
                      free_names, base, curves, weighting, restarts, observation)

    best = None
    for guess in starts:
        try:
            trial = least_squares(
                residuals, guess, bounds=(lower, upper), method="trf",
                args=(free_names, base, curves, weighting, observation),
                x_scale="jac", max_nfev=300,
            )
        except Exception as error:  # a solver blow-up must not kill the run
            print(f"    restart: {type(error).__name__}: {error}")
            continue
        if best is None or trial.cost < best.cost:
            best = trial
    if best is None:
        return FitResult(base, tuple(free_names), curves, np.inf, np.nan, np.nan,
                         {}, {}, np.zeros((0, 0)), np.inf,
                         success=False, message="every restart failed",
                         observation=observation)

    constants = unpack(best.x, free_names, base)
    at_bound = {}
    for name, value in zip(free_names, best.x):
        low, high = BOUNDS[name]
        span = high - low
        if value <= low + 1e-6 * span:
            at_bound[name] = "lower"
        elif value >= high - 1e-6 * span:
            at_bound[name] = "upper"

    errors, correlation, condition = _statistics(best, free_names, curves, weighting)
    per_curve = _per_curve(constants, curves, observation)
    finite = [row["rms"] for row in per_curve if np.isfinite(row["rms"])]
    sigmas = [row["sigma"] for row in per_curve if np.isfinite(row.get("sigma", np.nan))]

    return FitResult(
        constants=constants,
        free_names=tuple(free_names),
        curves=curves,
        cost=float(best.cost),
        rms_absorbance=float(np.sqrt(np.mean(np.square(finite)))) if finite else np.nan,
        rms_sigma=float(np.mean(sigmas)) if sigmas else np.nan,
        at_bound=at_bound,
        standard_errors=errors,
        correlation=correlation,
        condition_number=condition,
        per_curve=per_curve,
        success=bool(best.success),
        message=str(best.message),
        observation=observation,
    )


def _assert_observation_design(curves):
    """
    The observation switch is only meaningful where the design matches it.

    A catalysed curve (e0 > 0) whose reference did NOT omit the enzyme is not a
    catalytic increment, so `increment` would be the wrong observable for it --
    and a curve with e0 = 0 whose reference DID omit the enzyme is a catalysed
    curve mislabelled as a background. Neither is fitted as either: the design
    is read off the cuvette table, independently of [enz], so a disagreement is
    a defect to rule on, not a curve to guess at.
    """
    for curve in curves:
        if curve.conditions.e0 > 0 and curve.reference_omits != "enzyme":
            raise ValueError(
                f"exp {curve.experiment} sample {curve.sample} has e0 = "
                f"{curve.conditions.e0} but its reference omits "
                f"{curve.reference_omits!r}, not the enzyme: it is not a "
                f"catalytic increment and is not fitted as one")
        if curve.conditions.e0 == 0 and curve.reference_omits == "enzyme":
            raise ValueError(
                f"exp {curve.experiment} sample {curve.sample} has e0 = 0 but "
                f"its reference omits the enzyme: it is a catalysed curve "
                f"mislabelled as a background")


def sequential_fit(curves, weighting="curve", restarts=4, seed=0, on_stage=None,
                   extended=False, stages=None, observation="design"):
    """
    Stage 1 on the enzyme-free curves, then stage 2 on the catalysed ones with
    stage 1 frozen. Returns (stage_one, stage_two); stage_two is None when the
    block has no catalysed curves.

    `on_stage(name, result)` is called as each stage finishes. A full fit takes
    tens of minutes, and stage 1's answer is worth having on screen before
    stage 2 starts rather than after it ends.

    `extended=True` swaps in `STAGE_ONE_EXTENDED` / `STAGE_TWO_EXTENDED`; a
    `stages=(one, two)` pair names a specific rung of the M0-M4 ladder
    (`MODEL_STAGES`). The old path is untouched: stage 1 zeroes k5 and k6 as
    before and the base constants keep the extensions OFF.

    `observation` is passed through to `fit_group` and defaults to "design"
    (the increment for catalysed curves). The design guard runs first: a curve
    whose reference design contradicts its [enz] is a defect and stops the fit.
    """
    _assert_observation_design(curves)
    if stages is not None:
        stage_one, stage_two = stages
    else:
        stage_one = STAGE_ONE_EXTENDED if extended else STAGE_ONE
        stage_two = STAGE_TWO_EXTENDED if extended else STAGE_TWO
    enzyme_free = [c for c in curves if c.conditions.e0 == 0]
    catalysed = [c for c in curves if c.conditions.e0 > 0]
    if not enzyme_free:
        raise ValueError("no enzyme-free curves in this block: "
                         "stage 1 has nothing to determine the background from")

    # The extended starting table only for a model that actually has
    # extensions free; M0 through `stages` is still the old path.
    old_path = stages is None and not extended
    start_table = None if old_path or stages == MODEL_STAGES["M0"] else EXTENDED_INITIAL
    first = fit_group(enzyme_free, stage_one, weighting=weighting,
                      restarts=restarts, seed=seed, initial=start_table,
                      observation=observation)
    # Enzyme-free curves carry no information about k5 or k6 -- both are
    # multiplied by E0 = 0 -- so stage 1 must not report whatever value they
    # happened to be seeded with. Zeroing them keeps the reported constants
    # honest; stage 2's own starting values come from INITIAL, not from here.
    first = FitResult(**{**first.__dict__,
                         "constants": first.constants.replace(k5=0.0, k6=0.0)})
    if on_stage:
        on_stage("stage_1", first)
    if not catalysed:
        return first, None
    second = fit_group(catalysed, stage_two, base=first.constants,
                       weighting=weighting, restarts=restarts, seed=seed,
                       initial=start_table, observation=observation)
    if on_stage:
        on_stage("stage_2", second)
    return first, second


def profile_r(curves, weighting="curve", restarts=3, seed=0, values=None):
    """
    Profile likelihood over r: fix r, fit the rest, record the cost.

    This is how a parameter correlated at -0.999 with another has to be
    reported. The fit determines the product k_can*r far better than it
    determines either factor, so a single best-fit r with a standard error would
    overstate what the data says. The profile shows the whole range of r the
    curves tolerate, and -- since r <= 1 makes a lag impossible -- whether they
    tolerate any value the spectroscopy would accept.

    Returns a list of {r, cost, k_can, k3, k0} ordered by r.
    """
    values = values or (0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0)
    free = tuple(name for name in STAGE_ONE if name != "r")
    profile = []
    for value in values:
        base = RateConstants(**{name: 10.0 ** INITIAL[name] for name in
                                ("k_can", "k3", "k0", "k5", "k6")}, r=value)
        result = fit_group(curves, free, base=base, weighting=weighting,
                           restarts=restarts, seed=seed)
        profile.append({
            "r": value, "cost": result.cost,
            **{name: getattr(result.constants, name) for name in free},
            "k_can_times_r": result.constants.k_can * value,
        })
    return profile


def report_profile(profile):
    best = min(profile, key=lambda row: row["cost"])
    print("\n  profile over r (r fixed, the rest refitted)")
    print(f"  {'r':>6s} {'cost':>12s} {'cost/min':>9s} {'k_can':>12s} "
          f"{'k_can*r':>12s} {'k3':>11s} {'k0':>11s}")
    for row in profile:
        ratio = row["cost"] / best["cost"] if best["cost"] > 0 else np.inf
        mark = "  <- best" if row is best else ""
        print(f"  {row['r']:6.2f} {row['cost']:12.5g} {ratio:9.2f} "
              f"{row['k_can']:12.4g} {row['k_can_times_r']:12.4g} "
              f"{row['k3']:11.4g} {row['k0']:11.4g}{mark}")
    products = [row["k_can_times_r"] for row in profile if row["r"] > 0]
    if products:
        print(f"  k_can*r varies by {max(products) / min(products):.2g}x across the "
              f"profile while k_can alone varies by "
              f"{max(r['k_can'] for r in profile) / min(r['k_can'] for r in profile):.2g}x")
    return best


# --- reporting -------------------------------------------------------------

UNITS = {"k_can": "mM^-2 s^-1", "k3": "mM^-1 s^-1", "k0": "mM^-1 s^-1",
         "k5": "mM^-2 s^-1", "k6": "mM^-1 s^-1", "r": "",
         "k_sink": "s^-1", "K4": "mM^-1", "km_s": "mM",
         "k_act_r": "s^-1", "K_act": "mM^-1", "km_s_background": "mM"}


def report(result, title):
    print(f"\n{title}")
    print(f"  {len(result.curves)} curves / "
          f"{len({c.experiment for c in result.curves})} experiments, "
          f"{sum(len(c) for c in result.curves):,} points")
    if not result.success:
        print(f"  optimiser did not converge: {result.message}")

    print(f"\n  {'parameter':10s} {'value':>12s} {'+- (log10)':>11s} {'units':12s}")
    for name in result.free_names:
        value = getattr(result.constants, name)
        error = result.standard_errors.get(name, np.nan)
        flag = f"  <- AT {result.at_bound[name].upper()} BOUND" if name in result.at_bound else ""
        shown = f"{value:12.4g}" if name != "r" else f"{value:12.4f}"
        print(f"  {name:10s} {shown} {error:11.3f} {UNITS[name]:12s}{flag}")

    print(f"\n  rms residual   {result.rms_absorbance:.5f} AU "
          f"({result.rms_sigma:.1f}x the curves' own noise on average)")
    print(f"  Jacobian condition number  {result.condition_number:.3g}"
          + ("   <- badly determined" if result.condition_number > 1e8 else ""))

    if result.correlation.size and np.all(np.isfinite(result.correlation)):
        print("\n  parameter correlations")
        print("            " + " ".join(f"{n:>8s}" for n in result.free_names))
        for i, name in enumerate(result.free_names):
            row = " ".join(f"{result.correlation[i, j]:8.3f}"
                           for j in range(len(result.free_names)))
            print(f"  {name:8s}  {row}")
        off = [(abs(result.correlation[i, j]), result.free_names[i], result.free_names[j])
               for i in range(len(result.free_names))
               for j in range(i + 1, len(result.free_names))]
        worst = max(off) if off else None
        if worst and worst[0] > 0.95:
            print(f"  {worst[1]} and {worst[2]} are {worst[0]:.3f} correlated: "
                  f"this block constrains their combination, not each separately")

    shapes = [(row["peak_data"], row["peak_model"]) for row in result.per_curve
              if np.isfinite(row.get("peak_data", np.nan))]
    if shapes:
        late_data = sum(1 for d, _ in shapes if d > LAG_PEAK_FRACTION)
        late_model = sum(1 for _, m in shapes if m > LAG_PEAK_FRACTION)
        print(f"\n  shape: {late_data}/{len(shapes)} measured curves reach peak slope "
              f"past 15% into the run")
        print(f"         {late_model}/{len(shapes)} modelled curves do")
        if late_data and not late_model:
            print("         the model reproduces no lag at all -- see the r > 1 bound "
                  "in kinetic_model.py")

    worst_curves = sorted((row for row in result.per_curve if np.isfinite(row["rms"])),
                          key=lambda row: -row["sigma"])[:5]
    if worst_curves:
        print("\n  worst-fitted curves")
        for row in worst_curves:
            print(f"    exp {row['experiment']:3d} sample {row['sample']}: "
                  f"rms {row['rms']:.4f} AU = {row['sigma']:6.1f} sigma  "
                  f"(net: data {row['net_data']:+.4f}, model {row['net_model']:+.4f})")


def constants_from_record(record):
    """
    A saved stage's constant table as a `RateConstants`.

    `to_dict` writes every name in PARAMETER_NAMES and writes None for a
    non-finite OFF value (km_s, k_act_r, km_s_background); an older save
    carries the original six only. Both are filled from the dataclass defaults,
    so any save this module has written stays readable.
    """
    defaults = RateConstants()
    values = {}
    for name in PARAMETER_NAMES:
        value = record["constants"].get(name)
        values[name] = getattr(defaults, name) if value is None else float(value)
    return RateConstants(**values)


def to_dict(result, title):
    """A JSON-safe record of one stage, for --save."""
    return {
        "block": title,
        "free_parameters": list(result.free_names),
        "constants": {name: (None if not np.isfinite(getattr(result.constants, name))
                             else float(getattr(result.constants, name)))
                      for name in PARAMETER_NAMES},
        "standard_errors_log10": {k: (None if not np.isfinite(v) else float(v))
                                  for k, v in result.standard_errors.items()},
        "at_bound": result.at_bound,
        "cost": result.cost,
        "rms_absorbance": result.rms_absorbance,
        "rms_sigma": result.rms_sigma,
        "condition_number": (None if not np.isfinite(result.condition_number)
                             else result.condition_number),
        "correlation": (result.correlation.tolist()
                        if np.all(np.isfinite(result.correlation)) else None),
        "curves": len(result.curves),
        "experiments": sorted({c.experiment for c in result.curves}),
        "per_curve": result.per_curve,
        "converged": result.success,
        "observation": result.observation,
    }


# Which model each rung of the ladder EXTENDS. A term is read against the model
# it adds to, never against whichever rung was listed before it. Until
# 2026-09-14 `ladder_f_test` did the latter, so M2's K4 was differenced against
# M1b -- a model M2 does not contain -- and came back at F = -197.
MODEL_PARENTS = {"M1": "M0", "M1b": "M1", "M2": "M1", "M3": "M2", "M4": "M3",
                 "M4b": "M1b"}

# How far, relatively, a larger model's cost may sit above its parent's before
# it is called an optimiser failure. At its optimum a model that CONTAINS
# another can never fit worse, so anything past round-off means the optimiser
# missed the larger model's optimum, and no F built on that pair means anything.
NESTING_TOLERANCE = 1e-6

# A curve "lags" when its steepest slope comes later than this share of the
# run -- the count `report` prints and FITTING.md's lag fraction reads.
LAG_PEAK_FRACTION = 0.15


def _ladder_saves(block, models, directory):
    """The `--model` saves for one block, keyed by model; absent rungs skipped."""
    import os

    substrate, temperature, buffer_name = block
    saves = {}
    for model in models:
        path = os.path.join(directory, f"{substrate}_{temperature:.0f}C_"
                            f"{buffer_name}_{model}.json")
        if os.path.exists(path):
            with open(path) as handle:
                saves[model] = json.load(handle)
    return saves


def _same_background(one, two):
    """
    Whether two saved stage-1 constant tables are the same background.

    A name a save does not carry is at its OFF default (the M0 saves predate
    the extensions), and None is how `to_dict` writes an infinite OFF value.
    """
    defaults = RateConstants()
    for name in PARAMETER_NAMES:
        values = []
        for table in (one, two):
            value = table.get(name, getattr(defaults, name))
            values.append(np.inf if value is None else float(value))
        first, second = values
        if np.isinf(first) or np.isinf(second):
            if first != second:
                return False
        elif not np.isclose(first, second, rtol=1e-9, atol=0.0):
            return False
    return True


def ladder_f_test(block, models=tuple(MODEL_STAGES), directory="data/fits",
                  points=None):
    """
    The weighted F of each rung of the M0-M4 ladder against the model it
    EXTENDS (`MODEL_PARENTS`), per stage, with a verdict on whether that F is
    a comparison at all.

    `F = ((cost_small - cost_big)/delta_p) / (cost_big/(N - p_big))` on the
    weighted least-squares cost the fitter minimises, N the residual points
    and p_big the larger model's free parameters, against the package's
    `TWO_PHASE_F`. The verdicts, in the order they are tested:

      - "not nested": the parent frees a parameter the child does not.
      - "not nested: stage 1 differs": a STAGE-2 cost is only comparable when
        both models froze the SAME background. M3 adds `k_sink` to stage 1,
        so its stage 2 is fitted on a different base from M2's, and
        differencing the two costs prices the background, not a stage-2 term.
      - "nothing added": the child frees no new parameter in this stage.
      - "optimiser failure": the larger model fits WORSE than the one it
        contains (beyond `NESTING_TOLERANCE`), which cannot happen at its
        optimum.
      - "parent not converged" / "not converged": the optimiser's own flag.
      - "earns" / "does not earn": F against the bar.

    READ "EARNS" WITH THE MISFIT BESIDE IT. These fits sit tens to hundreds of
    times the noise over thousands of serially correlated readings, and there
    any term that absorbs any systematic misfit clears F = 12 -- one cleared
    it at 259,200. "Earns" says the term helps this model; it does not say the
    term is the missing chemistry.

    `points` is {1: N_stage1, 2: N_stage2}; left None it is counted from the
    `build_curves` the fits read, so the denominator is the fit's own.
    """
    saves = _ladder_saves(block, models, directory)
    if points is None:
        curves, _ = build_curves()
        scoped = [c for c in curves if c.group == block]
        points = {1: sum(len(c) for c in scoped if c.conditions.e0 == 0),
                  2: sum(len(c) for c in scoped if c.conditions.e0 > 0)}
    rows = []
    for model in models:
        if model not in saves:
            continue
        for stage, key in ((1, "stage_1"), (2, "stage_2")):
            entry = saves[model].get(key)
            if entry is None:
                continue
            cost = float(entry["cost"])
            free = list(entry["free_parameters"])
            row = {"stage": stage, "model": model, "free": len(free),
                   "cost": cost, "rms_absorbance": entry["rms_absorbance"],
                   "converged": bool(entry.get("converged", True)),
                   "vs": None, "added": "", "f": np.nan, "verdict": "baseline"}
            parent = MODEL_PARENTS.get(model)
            small = saves.get(parent, {}).get(key) if parent else None
            if small is None:
                rows.append(row)
                continue
            row["vs"] = parent
            added = [name for name in free if name not in small["free_parameters"]]
            row["added"] = ",".join(added)
            if not set(small["free_parameters"]) <= set(free):
                row["verdict"] = "not nested"
            elif stage == 2 and not _same_background(
                    saves[parent]["stage_1"]["constants"],
                    saves[model]["stage_1"]["constants"]):
                row["verdict"] = "not nested: stage 1 differs"
            elif not added:
                row["verdict"] = "nothing added"
            else:
                degrees = points[stage] - len(free)
                small_cost = float(small["cost"])
                row["f"] = ((small_cost - cost) / len(added) / (cost / degrees)
                            if degrees > 0 and cost > 0 else np.nan)
                if cost > small_cost * (1.0 + NESTING_TOLERANCE):
                    row["verdict"] = "optimiser failure"
                elif not small.get("converged", True):
                    row["verdict"] = "parent not converged"
                elif not row["converged"]:
                    row["verdict"] = "not converged"
                else:
                    row["verdict"] = ("earns" if row["f"] > TWO_PHASE_F
                                      else "does not earn")
            rows.append(row)
    return pd.DataFrame(rows)


def ladder_checks(block, models=tuple(MODEL_STAGES), directory="data/fits",
                  curves=None, saves=None):
    """
    What the F test cannot see, per model and stage, off the `--model` saves.

    - `converged`, whether a correlation matrix was saved (`to_dict` drops a
      non-finite one), the largest |correlation| and its pair, and `at_bound`.
      A pair at 1.000 is one identified combination, not two constants.
    - The SUBSTRATE ORDER the data and the model give, through
      `scope.orders` with one offset per run, on each curve's NET RISE. That
      is a proxy for the initial-rate order FITTING.md F1 quotes, and good for
      the comparison it is used for here: whether a fitted saturation leaves
      the model's order where the data's is.
    - The lag counts, `peak_* > LAG_PEAK_FRACTION`, data against model.
    - `induction.composition_collinearity` on the same curves: where [buf]
      moves with [S] inside the runs, an order or a saturation "in [S]" is one
      in the PAIR, and a background that saturates in substrate is equally
      one that falls with buffer.
    """
    import induction
    import scope

    if saves is None:
        saves = _ladder_saves(block, models, directory)
    if curves is None:
        curves, _ = build_curves()
    conditions = pd.DataFrame([
        {"experiment": c.experiment, "sample": c.sample, "s0": c.conditions.s0,
         "buf": c.conditions.buf, "live": True}
        for c in curves if c.group == block])
    rows = []
    for model in models:
        if model not in saves:
            continue
        for stage, key in ((1, "stage_1"), (2, "stage_2")):
            entry = saves[model].get(key)
            if entry is None:
                continue
            table = pd.DataFrame(entry["per_curve"]).merge(
                conditions, on=["experiment", "sample"])
            data = scope.orders("net_data", frame=table, terms=("s0",),
                                live_only=False)
            fitted = scope.orders("net_model", frame=table, terms=("s0",),
                                  live_only=False)
            collinear = induction.composition_collinearity(table)
            names = entry["free_parameters"]
            raw = entry.get("correlation")
            largest, pair = np.nan, None
            if raw is not None and len(names) > 1:
                magnitude = np.abs(np.array(raw, dtype=float))
                np.fill_diagonal(magnitude, 0.0)
                i, j = np.unravel_index(int(np.nanargmax(magnitude)),
                                        magnitude.shape)
                largest, pair = float(magnitude[i, j]), f"{names[i]}/{names[j]}"
            rows.append({
                "stage": stage, "model": model, "curves": int(len(table)),
                "converged": bool(entry.get("converged", True)),
                "correlation_saved": raw is not None,
                "largest_correlation": largest, "pair": pair,
                "at_bound": ",".join(f"{name}:{side}" for name, side
                                     in (entry.get("at_bound") or {}).items()),
                "order_data": data["order_s0"], "stderr_data": data["stderr_s0"],
                "order_model": fitted["order_s0"],
                "stderr_model": fitted["stderr_s0"],
                "lag_data": int((table.peak_data > LAG_PEAK_FRACTION).sum()),
                "lag_model": int((table.peak_model > LAG_PEAK_FRACTION).sum()),
                "s0_buf_runs": collinear.get("runs", 0),
                "s0_buf_median_r": collinear.get("median", np.nan),
                "s0_buf_constant_runs": collinear.get("constant_buffer", 0)})
    return pd.DataFrame(rows)


def increment_comparison(block, directory="data/fits", curves=None):
    """
    R0.0's report: the M0 refit on the catalytic INCREMENT beside the old
    absolute M0 (`<block>_R_M0.json` against `<block>_M0.json`).

    Per stage it prints the rms in AU and the cost for both saves, and for
    stage 2 the fitted k5, k6 and their correlation. It also prints the peak
    [PBA] the refit's constants reach on three catalysed curves -- F6's hidden
    intermediate, the test of whether the autocatalytic loop does anything --
    and `ladder_checks`' substrate order and lag counts for the refit.

    Stage 1 must be unchanged by R0.0: the observation switch only touches
    curves whose reference omitted the enzyme, and every stage-1 curve is
    enzyme-free.
    """
    import os

    substrate, temperature, buffer_name = block
    stem = f"{substrate}_{temperature:.0f}C_{buffer_name}"
    with open(os.path.join(directory, f"{stem}_M0.json")) as handle:
        old = json.load(handle)
    with open(os.path.join(directory, f"{stem}_R_M0.json")) as handle:
        refit = json.load(handle)

    print(f"\n{stem}: M0 on the increment (R_M0) beside M0 absolute (M0)")
    for stage in ("stage_1", "stage_2"):
        if stage not in old or stage not in refit:
            continue
        for label, save in (("absolute", old), ("increment", refit)):
            entry = save[stage]
            print(f"  {stage} {label}: rms {entry['rms_absorbance']:.5f} AU "
                  f"({entry['rms_sigma']:.1f}x noise), cost {entry['cost']:.5g}")
        if stage == "stage_2":
            for name in ("k5", "k6"):
                print(f"    {name}: absolute {old[stage]['constants'].get(name):.4g}"
                      f"  ->  increment {refit[stage]['constants'].get(name):.4g}")
            correlation = refit[stage].get("correlation")
            names = refit[stage]["free_parameters"]
            if correlation is not None and "k5" in names and "k6" in names:
                i, j = names.index("k5"), names.index("k6")
                print(f"    correlation k5/k6 (increment): "
                      f"{correlation[i][j]:+.3f}")

    if curves is None:
        curves, _ = build_curves()
    scoped = [c for c in curves if c.group == block and c.conditions.e0 > 0]
    constants = constants_from_record(refit["stage_2"])
    experiments = sorted({c.experiment for c in scoped})[:3]
    print("  peak [PBA] (mM) at the increment constants (F6):")
    for experiment in experiments:
        samples = sorted([c for c in scoped if c.experiment == experiment],
                         key=lambda c: c.sample)
        peaks = []
        for curve in samples:
            trajectory = simulate(constants, curve.conditions, curve.times)
            peaks.append(np.nan if trajectory is None
                         else float(trajectory["PBA"].max()))
        print(f"    exp {experiment}: " + ", ".join(f"{p:.3g}" for p in peaks))

    checks = ladder_checks(block, models=("M0R",), curves=curves,
                           saves={"M0R": refit})
    columns = ["stage", "curves", "order_data", "order_model",
               "lag_data", "lag_model", "s0_buf_runs", "s0_buf_median_r"]
    print("  ladder_checks on the increment refit:")
    print(checks[columns].to_string(index=False))
    return checks


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--substrate", default="BnOH")
    parser.add_argument("--temperature", type=float, default=25.0)
    parser.add_argument("--buffer", default="Phosphate")
    parser.add_argument("--weighting", choices=("curve", "point"), default="curve")
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--observation", choices=("design", "absolute"),
                        default="design",
                        help="'design' (default) compares each catalysed curve "
                             "with the catalytic increment its reference cuvette "
                             "leaves; 'absolute' reproduces the pre-R0.0 fitter "
                             "exactly and is kept only for comparison")
    parser.add_argument("--profile-r", action="store_true",
                        help="profile the cost over r instead of fitting it, "
                             "which is how a -0.999-correlated parameter should "
                             "be reported")
    parser.add_argument("--list", action="store_true",
                        help="list the blocks that have both stages, then exit")
    parser.add_argument("--extended", action="store_true",
                        help="fit the extended model (sink, saturation, "
                             "activation) instead of the reduced one")
    parser.add_argument("--model", choices=tuple(MODEL_STAGES), default=None,
                        help="a rung of the nested M0-M4 ladder (Stage 3.3); "
                             "overrides --extended")
    parser.add_argument("--save", default=None)
    parser.add_argument("--scope", choices=("two-axis", "all"),
                        default="two-axis",
                        help="'two-axis' restricts to "
                             "fit_dataset.TWO_AXIS_BLOCK (exps 135-151, the "
                             "archive's only two-axis designs); 'all' uses "
                             "every fittable curve. See FITTING.md.")
    arguments = parser.parse_args()

    curves, _ = build_curves()
    if arguments.scope == "two-axis":
        # Restricting here rather than at the block key matters: exps 75 and 76
        # share the block's (BnOH, 25 C, Pyrophosphate) key but carry the
        # unresolved hexametaphosphate speciation question, so selecting the
        # block alone would silently pull them in.
        scoped = in_block(curves)
        if scoped:
            curves = scoped

    if arguments.list:
        free = group_curves(curves, enzyme_free=True)
        catalysed = group_curves(curves, enzyme_free=False)
        print(f"  {'substrate':11s} {'T':>4s} {'buffer':14s} {'E0=0':>6s} {'E0>0':>6s}  "
              f"sequential fit")
        for key in sorted(set(free) | set(catalysed), key=str):
            substrate, temperature, buffer_name = key
            first, second = len(free.get(key, [])), len(catalysed.get(key, []))
            if not first:
                verdict = "no background data"
            elif not second:
                verdict = "background only"
            elif first < 10:
                verdict = f"weak background ({first} curves)"
            else:
                verdict = "yes"
            print(f"  {substrate:11s} {temperature:4.0f} {buffer_name:14s} "
                  f"{first:6d} {second:6d}  {verdict}")
        return 0

    key = (arguments.substrate, arguments.temperature, arguments.buffer)
    block = [c for c in curves if c.group == key]
    if not block:
        print(f"no curves in block {key}" +
              (f" within --scope two-axis (exps {min(TWO_AXIS_BLOCK)}-"
               f"{max(TWO_AXIS_BLOCK)}, block {TWO_AXIS_GROUP})"
               if arguments.scope == "two-axis" else "") + "; try --list")
        return 1

    title = f"{arguments.substrate}, {arguments.temperature:.0f} C, {arguments.buffer}"
    print(f"=== {title} ===")
    print(f"weighting: per-{arguments.weighting}, {arguments.restarts} restarts")

    if arguments.profile_r:
        enzyme_free = [c for c in block if c.conditions.e0 == 0]
        print(f"\nSTAGE 1 profile: {len(enzyme_free)} enzyme-free curves")
        report_profile(profile_r(enzyme_free, weighting=arguments.weighting,
                                 restarts=arguments.restarts))
        return 0

    if arguments.model:
        stages = MODEL_STAGES[arguments.model]
    elif arguments.extended:
        stages = (STAGE_ONE_EXTENDED, STAGE_TWO_EXTENDED)
    else:
        stages = None
    stage_one, stage_two = stages if stages else (STAGE_ONE, STAGE_TWO)
    titles = {
        "stage_1": f"STAGE 1  enzyme-free  ->  {', '.join(stage_one)}",
        "stage_2": f"STAGE 2  catalysed  ->  {', '.join(stage_two)}   (stage 1 frozen)",
    }
    first, second = sequential_fit(block, weighting=arguments.weighting,
                                   restarts=arguments.restarts,
                                   extended=arguments.extended, stages=stages,
                                   observation=arguments.observation,
                                   on_stage=lambda name, result: report(result, titles[name]))
    if second is None:
        print("\nSTAGE 2 skipped: no catalysed curves in this block")

    if arguments.save:
        payload = {"block": title, "weighting": arguments.weighting,
                   "observation": arguments.observation,
                   "stage_1": to_dict(first, title)}
        if second is not None:
            payload["stage_2"] = to_dict(second, title)
        with open(arguments.save, "w") as handle:
            json.dump(payload, handle, indent=2)
        print(f"\nsaved to {arguments.save}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
