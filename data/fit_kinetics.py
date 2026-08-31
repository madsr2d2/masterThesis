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
from scipy.optimize import least_squares

from fit_dataset import BASELINE_POINTS, build_curves, group_curves
from kinetic_model import (LOG_PARAMETERS, Conditions, RateConstants,
                           observable, pack, unpack)

# Bounds in the optimiser's own coordinates: log10 for rate constants, linear
# for r. Wide enough not to shape the answer, tight enough that the integrator
# is not asked to do something absurd. A fitted value sitting ON a bound is
# reported as such -- for r especially, that is the result, not a detail.
BOUNDS = {
    "k_can": (-8.0, 8.0),
    "k3": (-10.0, 6.0),
    "k0": (-16.0, -2.0),
    "k5": (-8.0, 10.0),
    "k6": (-8.0, 10.0),
    "r": (0.0, 5.0),
}

# Starting points, in the same coordinates. The scaling behind k_can: the seed
# makes aldehyde at v0 = k0[H2O2][S], so over a run of length T the pool reaches
# ~v0*T, and steps 1-2 only matter once k_can*[HOO-]*(v0*T)^2 is comparable to
# v0 -- which for this dataset's numbers puts k_can near 1 mM^-2 s^-1.
INITIAL = {"k_can": 0.0, "k3": -2.0, "k0": -9.0, "k5": 0.0, "k6": 0.0, "r": 0.3}

STAGE_ONE = ("k_can", "k3", "k0", "r")
STAGE_TWO = ("k5", "k6")

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


def residuals(vector, free_names, base, curves, weighting="curve"):
    """Stacked weighted residuals, model minus measurement, in absorbance."""
    constants = unpack(vector, free_names, base)
    stacked = []
    for curve in curves:
        denominator = _weights(curve, weighting)
        signal = observable(constants, curve.conditions, curve.times,
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


def _per_curve(constants, curves):
    """
    One row per curve: how well it is fitted, and whether the model reproduces
    its SHAPE. `peak_data` and `peak_model` are where the steepest point sits as
    a fraction of the run -- the statistic MECHANISM.md used to falsify the
    aldehyde reading, so the fit is checked against it rather than only against
    a residual norm.
    """
    rows = []
    for curve in curves:
        signal = observable(constants, curve.conditions, curve.times)
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
            "peak_data": _peak_position(curve.absorbance, curve.times),
            "peak_model": _peak_position(model, curve.times),
        })
    return rows


def _peak_position(values, times):
    """
    Where the steepest point sits, as a fraction of the run. Returns 0.0 for a
    curve whose slope never rises above its initial value by more than 5% --
    without that guard a straight line's gradient is a flat array of ties and
    argmax picks an arbitrary index, which reads as a late peak that is not
    there.
    """
    if len(values) < 5:
        return np.nan
    slope = np.gradient(np.asarray(values, dtype=float), times)
    if slope[0] <= 0 or slope.max() <= 1.05 * slope[0]:
        return 0.0
    return float(times[np.argmax(slope)] / times[-1])


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


def _screen(points, free_names, base, curves, weighting, keep):
    """
    Ranks candidate starting points by cost and returns the best `keep`.

    One residual evaluation each. Points the integrator cannot handle come back
    at the failure penalty and sort themselves to the bottom, so no special
    casing is needed.
    """
    scored = []
    for point in points:
        residual = residuals(point, free_names, base, curves, weighting)
        scored.append((float(0.5 * np.sum(residual ** 2)), point))
    scored.sort(key=lambda pair: pair[0])
    return [point for _, point in scored[:keep]]


def fit_group(curves, free_names, base=None, weighting="curve", restarts=4, seed=0):
    """
    Fits `free_names` to `curves`, holding everything else at `base`.

    Two kinds of start are used, and the best result over all of them is
    returned: the guaranteed ones (nominal, plus the r ladder when r is free),
    which always run, and the `restarts` best of a screened Latin-hypercube
    sample. Neither alone is enough -- see _guaranteed_points and
    SCREEN_SAMPLES, each of which records the failure that put it here.
    """
    base = base or RateConstants(**{name: 10.0 ** INITIAL[name] if name in LOG_PARAMETERS
                                    else INITIAL[name] for name in INITIAL})
    lower = np.array([BOUNDS[name][0] for name in free_names])
    upper = np.array([BOUNDS[name][1] for name in free_names])
    start = np.array([INITIAL[name] for name in free_names])

    starts = _guaranteed_points(free_names, start, lower, upper)
    starts += _screen(_hypercube(free_names, lower, upper, seed),
                      free_names, base, curves, weighting, restarts)

    best = None
    for guess in starts:
        try:
            trial = least_squares(
                residuals, guess, bounds=(lower, upper), method="trf",
                args=(free_names, base, curves, weighting),
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
                         success=False, message="every restart failed")

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
    per_curve = _per_curve(constants, curves)
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
    )


def sequential_fit(curves, weighting="curve", restarts=4, seed=0, on_stage=None):
    """
    Stage 1 on the enzyme-free curves, then stage 2 on the catalysed ones with
    stage 1 frozen. Returns (stage_one, stage_two); stage_two is None when the
    block has no catalysed curves.

    `on_stage(name, result)` is called as each stage finishes. A full fit takes
    tens of minutes, and stage 1's answer is worth having on screen before
    stage 2 starts rather than after it ends.
    """
    enzyme_free = [c for c in curves if c.conditions.e0 == 0]
    catalysed = [c for c in curves if c.conditions.e0 > 0]
    if not enzyme_free:
        raise ValueError("no enzyme-free curves in this block: "
                         "stage 1 has nothing to determine the background from")

    first = fit_group(enzyme_free, STAGE_ONE, weighting=weighting,
                      restarts=restarts, seed=seed)
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
    second = fit_group(catalysed, STAGE_TWO, base=first.constants,
                       weighting=weighting, restarts=restarts, seed=seed)
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
         "k5": "mM^-2 s^-1", "k6": "mM^-1 s^-1", "r": ""}


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
        late_data = sum(1 for d, _ in shapes if d > 0.15)
        late_model = sum(1 for _, m in shapes if m > 0.15)
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


def to_dict(result, title):
    """A JSON-safe record of one stage, for --save."""
    return {
        "block": title,
        "free_parameters": list(result.free_names),
        "constants": {name: getattr(result.constants, name)
                      for name in ("k_can", "k3", "k0", "k5", "k6", "r")},
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
    }


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--substrate", default="BnOH")
    parser.add_argument("--temperature", type=float, default=25.0)
    parser.add_argument("--buffer", default="Phosphate")
    parser.add_argument("--weighting", choices=("curve", "point"), default="curve")
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--profile-r", action="store_true",
                        help="profile the cost over r instead of fitting it, "
                             "which is how a -0.999-correlated parameter should "
                             "be reported")
    parser.add_argument("--list", action="store_true",
                        help="list the blocks that have both stages, then exit")
    parser.add_argument("--save", default=None)
    arguments = parser.parse_args()

    curves, _ = build_curves()

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
        print(f"no curves in block {key}; try --list")
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

    titles = {
        "stage_1": f"STAGE 1  enzyme-free  ->  {', '.join(STAGE_ONE)}",
        "stage_2": f"STAGE 2  catalysed  ->  {', '.join(STAGE_TWO)}   (stage 1 frozen)",
    }
    first, second = sequential_fit(block, weighting=arguments.weighting,
                                   restarts=arguments.restarts,
                                   on_stage=lambda name, result: report(result, titles[name]))
    if second is None:
        print("\nSTAGE 2 skipped: no catalysed curves in this block")

    if arguments.save:
        payload = {"block": title, "weighting": arguments.weighting,
                   "stage_1": to_dict(first, title)}
        if second is not None:
            payload["stage_2"] = to_dict(second, title)
        with open(arguments.save, "w") as handle:
            json.dump(payload, handle, indent=2)
        print(f"\nsaved to {arguments.save}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
