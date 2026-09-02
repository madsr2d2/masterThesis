"""
Why the 4OMe-BnOH progress curves slow down.

The temperature series (exps 14-19) and every enzyme-free 4OMe run in the
archive share one feature: the rate rises to a maximum and then falls, at
conversions of 0.3-3%. `temperature_series/ANALYSIS.md` 3a describes the fall
with a second exponential and says explicitly that the description must stay a
description. This module is the attempt to name the process.

Four candidates survive `ANALYSIS.md` 3a's eliminations:

    depletion    the substrate runs out          EXCLUDED, 0.3-1.1% conversion
    death        the catalyst dies               EXCLUDED for BnOH (Selwyn)
    clock        a shared reagent is consumed    dA/dt = v(t) . exp(-k t)
    sink         the chromophore reacts onward   dA/dt = v(t) - k . A
    inhibition   the product blocks turnover     dA/dt = v(t) / (1 + A/Ki)

`clock`, `sink` and `inhibition` are written here as PROGRESS-CURVE MODELS with
the same parameter count -- an offset, an amplitude, an induction time constant
and one decay parameter -- so "which of these is the fall" is a comparison of
residuals at equal cost rather than an argument about plausibility. Every one
of them carries the same induction term `v(t) = v (1 - exp(-t/tau))`, because
these curves have an induction period whatever else they have, and a model that
cannot represent it will spend its decay parameter on the rise.

The three differ in one thing only, and it is the thing in question:

    clock       the CAPACITY to react decays with time; the product level is
                irrelevant. Plateaus at v/k - v/(k + 1/tau).
    sink        the SIGNAL is consumed in proportion to how much has built up.
                Plateaus at v/k. Shape is independent of v.
    inhibition  the RATE is divided by how much product has built up.
                Never plateaus: A grows as sqrt(2 Ki v t) forever.

So the discriminating question is whether the curve approaches a horizontal
asymptote, and the discriminating design is one that changes v while holding
time fixed -- which is what the substrate ladder inside each run does (2.2x)
and what the six temperatures do (20x).

    plateau reached, shape independent of v   ->  clock or sink
    no plateau, shape depends on v            ->  inhibition

Distinguishing `clock` from `sink` is NOT possible from one curve's shape --
both are two exponentials approaching a plateau -- and this module does not
pretend otherwise. They are separated by where else the decay appears:
`background_decay` asks whether the enzyme-free cuvette, which has no catalyst
to inhibit and no catalytic cycle to poison, decays the same way.

    python data/slowdown.py

Names here are unique to this module; test_curve_metrics.test_no_duplicate_
definitions enforces that.
"""
import sys
import os
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from curve_metrics import model_residual
from fit_dataset import source_floor
import scope


# Grid geometry, shared by all three models so their residuals are comparable.
# Both axes are profiled on a log grid and the remaining two parameters are
# linear at every node, which is the same device fit_burst and fit_two_phase
# use: no optimiser, so no local minimum to fall into and no starting guess to
# defend. 40 x 40 is coarser than fit_burst's 240 because it is two dimensional
# and because nothing downstream quotes a time constant from here -- the
# comparison is between RESIDUALS, and a residual surface is flat near its
# minimum where a parameter is not.
SLOWDOWN_GRID = 40
INDUCTION_FLOOR = 1 / 300.0     # grid start, as a fraction of the run
INDUCTION_CAP = 2.0             # grid end, as a multiple of the run
# The inhibition model's second axis is dimensionless: u = v/Ki has units of
# 1/time, so u * (run length) says how far into inhibition the run reaches.
# Below 1e-3 the model is a straight line and Ki is unidentifiable; above 1e3
# the first reading is already inhibited. Both ends are meant to be reachable
# and neither is meant to be chosen.
INHIBITION_DEPTH = (1e-3, 1e3)

# A model has to beat the others by more than rounding to be reported as
# beating them. 1% in the residual sum of squares on ~370 readings is well
# outside what the grid's own coarseness can produce and well inside what a
# real difference in shape produces -- the falling curves separate the models
# by 20-80%.
SLOWDOWN_MARGIN = 0.01

# The models cost four parameters each; the descriptive two-phase form costs
# six. Residuals are reported per curve in units of the curve's own noise,
# through curve_metrics.model_residual, which divides by (n - parameters).
SLOWDOWN_PARAMETERS = 4


@dataclass(frozen=True)
class SlowdownFit:
    """
    One mechanistic progress-curve model fitted to one curve.

    `offset` and `amplitude` are the two linear parameters, `induction` the
    rise time constant and `decay` the one parameter that differs between the
    models -- a rate constant for `clock` and `sink`, and v/Ki for
    `inhibition`. `extra` carries whatever else the model defines: Ki, or the
    plateau the model is heading for.
    """
    name: str
    offset: float
    amplitude: float
    induction: float
    decay: float
    sse: float
    points: int
    extra: dict

    @property
    def rms(self):
        return float(np.sqrt(self.sse / self.points)) if self.points else np.nan


def _induction_integral(times, induction):
    """The production integral of a first-order rise: t - tau(1 - e^(-t/tau))."""
    return times - induction * (1.0 - np.exp(-times / induction))


def _clock_shape(times, induction, decay):
    """
    Integral of v(1 - e^(-t/tau)) e^(-k t), per unit v.

    The capacity to react is what decays, so the decay multiplies the RATE and
    the accumulated signal stops growing because nothing is being made any
    more -- the shape a consumed oxidant or a dying catalyst produces.
    """
    combined = decay + 1.0 / induction
    return ((1.0 - np.exp(-decay * times)) / decay
            - (1.0 - np.exp(-combined * times)) / combined)


def _sink_shape(times, induction, decay):
    """
    Solution of A' + kA = v(1 - e^(-t/tau)) with A(0) = 0, per unit v.

    Production never slows; what is made is consumed in proportion to how much
    of it there is. The shape is INDEPENDENT of v, which is the property that
    separates it from inhibition on a substrate ladder.
    """
    gap = decay - 1.0 / induction
    return ((1.0 - np.exp(-decay * times)) / decay
            - (np.exp(-times / induction) - np.exp(-decay * times)) / gap)


def _inhibition_shape(times, induction, depth):
    """
    Solution of A' = v(1 - e^(-t/tau))/(1 + A/Ki), per unit Ki.

    Separating gives A + A^2/(2 Ki) = v W(t), so A = Ki(sqrt(1 + 2 u W) - 1)
    with u = v/Ki. Written this way the model is LINEAR in Ki at fixed
    (tau, u), which is what lets it be profiled on the same grid as the other
    two. It has no asymptote: A keeps growing as sqrt(2 Ki v t).
    """
    return np.sqrt(1.0 + 2.0 * depth * _induction_integral(times, induction)) - 1.0


_SHAPES = {"clock": _clock_shape, "sink": _sink_shape,
           "inhibition": _inhibition_shape}


def _profile(times, values, shape, first_grid, second_grid):
    """
    Least squares over `offset + amplitude * shape(t; a, b)` on a 2-D grid.

    Both linear parameters are solved in closed form at every node from the
    node's own sufficient statistics, and the whole grid is one array
    operation. The alternative -- an optimiser over four parameters -- was not
    used anywhere else in this package for the reason given in fit_burst, and
    is not used here.
    """
    count = len(times)
    nodes = [(a, b) for a in first_grid for b in second_grid]
    design = np.empty((len(nodes), count))
    for index, (a, b) in enumerate(nodes):
        with np.errstate(all="ignore"):
            design[index] = shape(times, a, b)
    design[~np.isfinite(design)] = 0.0
    sum_y = float(values.sum())
    total = float(values @ values)
    sum_g = design.sum(axis=1)
    sum_gg = np.einsum("ij,ij->i", design, design)
    sum_gy = design @ values
    determinant = count * sum_gg - sum_g ** 2
    with np.errstate(all="ignore"):
        offset = (sum_gg * sum_y - sum_g * sum_gy) / determinant
        amplitude = (count * sum_gy - sum_g * sum_y) / determinant
        cost = total - offset * sum_y - amplitude * sum_gy
    cost = np.where(np.isfinite(cost) & (determinant > 0), cost, np.inf)
    best = int(np.argmin(cost))
    return (nodes[best][0], nodes[best][1], float(offset[best]),
            float(amplitude[best]), float(max(cost[best], 0.0)))


def fit_slowdown(times, values, name, points=SLOWDOWN_GRID):
    """
    Fit one of the three mechanistic forms to one curve.

    `name` is "clock", "sink" or "inhibition". Returns a SlowdownFit whose
    `sse` is directly comparable across the three: same data, same parameter
    count, same grid resolution.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    times = times - times[0]
    span = float(times[-1])
    blank = SlowdownFit(name, *([np.nan] * 5), len(times), {})
    if span <= 0 or len(times) < SLOWDOWN_PARAMETERS + 4:
        return blank

    induction = np.logspace(np.log10(span * INDUCTION_FLOOR),
                            np.log10(span * INDUCTION_CAP), points)
    if name == "inhibition":
        second = np.logspace(np.log10(INHIBITION_DEPTH[0] / span),
                             np.log10(INHIBITION_DEPTH[1] / span), points)
    else:
        # The decay constant is gridded as its own time constant on the same
        # range as the induction, so neither process is given a head start.
        # The 1.03 offset keeps k off 1/tau exactly, where _sink_shape divides
        # by zero; the grid is logarithmic so the shift is invisible.
        second = 1.0 / (np.logspace(np.log10(span * INDUCTION_FLOOR),
                                    np.log10(span * INDUCTION_CAP),
                                    points) * 1.03)
    first, decay, offset, amplitude, cost = _profile(
        times, values, _SHAPES[name], induction, second)

    extra = {}
    if name == "inhibition":
        # amplitude IS Ki; the rate constant follows from the profiled u.
        extra["ki"] = amplitude
        extra["rate"] = amplitude * decay
        extra["plateau"] = np.inf
    else:
        extra["rate"] = amplitude
        extra["plateau"] = float(
            amplitude * _SHAPES[name](np.array([1e9 * first]), first, decay)[0])
    return SlowdownFit(name, offset, amplitude, float(first), float(decay),
                       cost, len(times), extra)


def contest(curve, points=SLOWDOWN_GRID):
    """
    All three models against one curve, with residuals in units of its noise.

    Returns a dict: model name -> (SlowdownFit, residual). The residual is
    curve_metrics.model_residual, so it is the same quantity `progress_resid`
    reports for the descriptive form and the two can be compared directly.
    """
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    noise = curve.noise
    out = {}
    for name in _SHAPES:
        fitted = fit_slowdown(times, values, name, points=points)
        if not np.isfinite(fitted.sse):
            out[name] = (fitted, np.nan)
            continue
        shifted = times - times[0]
        with np.errstate(all="ignore"):
            predicted = fitted.offset + fitted.amplitude * _SHAPES[name](
                shifted, fitted.induction, fitted.decay)
        out[name] = (fitted, model_residual(values, predicted,
                                            SLOWDOWN_PARAMETERS, noise))
    return out


def verdict(residuals, margin=SLOWDOWN_MARGIN):
    """
    Which model won, or "tied" when the best two are inside `margin`.

    Reported on the residual, not on a p-value: the three forms are not nested,
    so there is no F test between them, and a likelihood ratio on serially
    correlated progress-curve residuals would be an overstatement dressed as a
    number. See TWO_PHASE_F for the same caution stated at greater length.
    """
    ranked = sorted((r, n) for n, r in residuals.items() if np.isfinite(r))
    if not ranked:
        return "unresolved", np.nan
    if len(ranked) == 1:
        return ranked[0][1], np.nan
    best, runner = ranked[0], ranked[1]
    lead = (runner[0] - best[0]) / best[0]
    return (best[1] if lead > margin else "tied"), float(lead)


# The rate has to fall by a stated amount before "when did it fall" means
# anything, and the amount has to be reachable inside a run. Half is the
# natural landmark for `inhibition` -- the rate halves exactly when A = Ki --
# and it is reached inside the run by most of the enzyme-free curves and by
# none of the 15-30 C catalysed ones, which is itself part of the answer.
DECAY_FRACTIONS = (0.9, 0.75, 0.5)

# The rate is read as the slope over this fraction of the run, which is the
# resolution at which these curves have a shape at all.
DECAY_WINDOW = 0.05

# Where the two hypotheses put the slope of log(time to fall) on log(peak
# rate). Product control fixes the PRODUCT at which the rate falls, so the
# time to reach it is inversely proportional to the rate; a clock fixes the
# TIME and lets the product land where it may.
PRODUCT_SLOPE = -1.0
CLOCK_SLOPE = 0.0


def decay_point(curve, fraction=0.5, samples=4000):
    """
    Where the fitted rate falls to `fraction` of its maximum: (t, A, v, t_max).

    Read off `summary_kinetics.fit_progress`'s CHOSEN form, not off the
    readings: the readings' own slope is too noisy to locate a landmark on,
    and the chosen form is already the one the F test says the curve earned.
    Evaluated only inside the observed span -- a landmark extrapolated past the
    last reading is a property of the model, not of the experiment -- so a
    curve whose rate never falls that far returns NaN and is counted as
    censored rather than silently given the run's end.
    """
    from summary_kinetics import fit_progress          # local: heavy import
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    shifted = times - times[0]
    fitted = fit_progress(times, values)
    grid = np.linspace(0.0, float(shifted[-1]), samples)
    predicted = fitted.predict(grid + times[0])
    if not np.all(np.isfinite(predicted)):
        return (np.nan,) * 4
    # A WINDOWED slope, not the analytic derivative. The fitted form is free
    # to put a spike at t = 0 -- a two-phase fit with a short tau1 and B1 < 0
    # does, and its instantaneous rate there can be ten times anything the
    # readings show -- and a landmark measured against that spike reports a
    # "fall" on curves that accelerate throughout. Averaging the fit over the
    # same window the block-slope tables use removes it, and costs nothing on
    # a curve whose fall is real and slow.
    step = max(1, int(round(samples * DECAY_WINDOW)))
    rate = (predicted[step:] - predicted[:-step]) / (grid[step:] - grid[:-step])
    grid = grid[:-step]
    predicted = predicted[:-step]
    top = int(np.argmax(rate))
    peak = float(rate[top])
    if not peak > 0:
        return (np.nan,) * 4
    # The fall has to still be there at the last reading. Without this a lag
    # curve whose fitted rate dips for a moment early on -- exps 9, 10 and 11
    # do, at 40-200 s and 0.0004 AU -- reports a "fall" that is a transient in
    # the fit, and those rows land in the regression below at a thousandth of
    # the product and a hundredth of the time of a real one.
    sustained = rate[-1] <= fraction * peak
    below = np.nonzero(rate[top:] <= fraction * peak)[0]
    if not len(below) or not sustained:
        return np.nan, np.nan, peak, float(grid[top])
    index = top + int(below[0])
    return (float(grid[index]), float(predicted[index] - predicted[0]),
            peak, float(grid[top]))


def decay_table(experiments, fraction=0.5):
    """
    One row per curve: the fall landmark, with everything needed to test it.

    Columns: experiment, sample, substrate, buffer, pH, temperature, s0, e0,
    span_s, v_peak, t_peak_s, t_fall_s, a_fall, censored.
    """
    import pandas as pd
    rows = []
    for curve in scope.curves(tuple(experiments)):
        t_fall, a_fall, peak, t_peak = decay_point(curve, fraction=fraction)
        rows.append({
            "experiment": curve.experiment, "sample": curve.sample,
            "substrate": curve.substrate, "buffer": curve.buffer,
            "pH": curve.pH, "temperature": curve.temperature,
            "s0": curve.conditions.s0, "e0": curve.conditions.e0,
            "h2o2": curve.conditions.h2o2, "buf": curve.buf,
            "epsilon": curve.epsilon, "noise": curve.noise,
            "span_s": float(curve.times[-1] - curve.times[0]),
            "net": float(curve.absorbance[-1] - curve.absorbance[0]),
            "v_peak": peak, "t_peak_s": t_peak,
            "t_fall_s": t_fall, "a_fall": a_fall,
            "censored": not np.isfinite(t_fall),
        })
    return pd.DataFrame(rows)


def decay_scaling(table):
    """
    Regress log(time to fall) on log(peak rate), and say which law it matches.

    This is the whole discrimination in one number. Both hypotheses predict a
    straight line through these points and they predict DIFFERENT SLOPES:

        product control   the rate falls at a fixed A, reached at t = A/v,
                          so the slope is -1
        a clock           the rate falls at a fixed t whatever v is,
                          so the slope is 0

    and nothing else in the design has to be held constant for that to be the
    test, which is why it is worth more than any single matched pair. The
    lever is the range of `v_peak` in the table; report it, because a slope
    fitted over a factor of two says very little and one fitted over sixty says
    a great deal.

    Censored rows -- the rate never fell that far inside the run -- are
    dropped, and their count is returned. They are not missing at random: a
    slow curve runs out of time before a fast one does, which biases the slope
    TOWARDS zero. So a slope near -1 survives the censoring; a slope near 0
    does not establish a clock on its own.
    """
    live = table[~table.censored & (table.v_peak > 0) & (table.t_fall_s > 0)]
    if len(live) < 3:
        return {"points": len(live), "censored": int(table.censored.sum()),
                "slope": np.nan, "stderr": np.nan, "lever": np.nan,
                "intercept": np.nan, "r2": np.nan}
    x = np.log(live.v_peak.to_numpy())
    y = np.log(live.t_fall_s.to_numpy())
    design = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    variance = float(resid @ resid) / max(1, len(x) - 2)
    covariance = variance * np.linalg.pinv(design.T @ design)
    spread = float(np.var(y))
    return {
        "points": int(len(live)), "censored": int(table.censored.sum()),
        "slope": float(beta[1]), "stderr": float(np.sqrt(covariance[1, 1])),
        "intercept": float(beta[0]),
        "lever": float(np.exp(x.max() - x.min())),
        "r2": float(1 - (variance * max(1, len(x) - 2)) / (spread * len(x)))
        if spread > 0 else np.nan,
    }


# The rolling window the sink linearisation reads its rate from. Wider than
# curve_metrics.LAG_WINDOW because the question here is the SLOPE of the rate
# against product over the whole tail, not where an induction ends, and a
# tenth of the run leaves the tail of a short curve with too few independent
# windows to regress.
SINK_WINDOW = 0.15
SINK_MINIMUM_POINTS = 8
# A curve has to be described by the straight line before which straight line
# it is can be argued about. 0.95 keeps the curves whose decline is deep
# enough to have a shape and drops the ones whose "decline" is a slope through
# the noise of a nearly flat tail.
SINK_CLEAN_R2 = 0.95
# The tail has to have fallen this far below the maximum before "the rate
# declines with product" is a measurement rather than a slope through noise.
SINK_DECLINE = 0.85


@dataclass(frozen=True)
class SinkFit:
    """
    The rate read against the product already made, after the rate maximum.

    Two hypotheses put a STRAIGHT LINE through different transforms of the
    same two columns, which is the whole reason to do it this way rather than
    by fitting whole-curve forms:

        sink        A' = v0 - k A        so the RATE is linear in A
        inhibition  A' = v0/(1 + A/Ki)   so 1/RATE is linear in A

    Both are 2-parameter lines on the same points, so `rate_r2` against
    `reciprocal_r2` is a fair comparison, and it is decided by the curves that
    decline furthest -- over a 30% decline the two transforms are both nearly
    straight and neither wins.

    `plateau` is v0/k, the absorbance the sink would settle at, and it is a
    prediction the run can be checked against. `ki` is the inhibition
    constant the reciprocal line implies. A NEGATIVE ki means the reciprocal
    line has the wrong sign to be an inhibition constant at all.
    """
    v0: float
    k: float
    plateau: float
    ki: float
    rate_r2: float
    reciprocal_r2: float
    decline: float
    points: int

    @property
    def prefers(self):
        if not np.isfinite(self.rate_r2) or not np.isfinite(self.reciprocal_r2):
            return "unresolved"
        if abs(self.rate_r2 - self.reciprocal_r2) < SLOWDOWN_MARGIN:
            return "tied"
        return "sink" if self.rate_r2 > self.reciprocal_r2 else "inhibition"


def _line_r2(x, y):
    """Slope, intercept and R^2 of a straight line, or nans."""
    if len(x) < 3 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        return np.nan, np.nan, np.nan
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    spread = float(((y - y.mean()) ** 2).sum())
    return (float(slope), float(intercept),
            float(1 - float(resid @ resid) / spread) if spread > 0 else np.nan)


def sink_fit(curve, fraction=SINK_WINDOW):
    """
    Regress the rolling rate on the accumulated product, after the maximum.

    The rate comes from `curve_metrics.rolling_slope` at the curve's own noise
    floor, and the product from the readings themselves at the same window
    centres -- no model, no extrapolation, and nothing that has to be believed
    before the answer means something.
    """
    from curve_metrics import rolling_slope                # local: cycle-free
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    floor = source_floor(curve.source)
    centres, slopes = rolling_slope(times, values, fraction, floor)
    blank = SinkFit(*([np.nan] * 7), 0)
    if len(slopes) < SINK_MINIMUM_POINTS:
        return blank
    product = np.interp(centres, times, values) - values[0]
    top = int(np.argmax(slopes))
    rate, made = slopes[top:], product[top:]
    if len(rate) < SINK_MINIMUM_POINTS or not np.all(rate > 0):
        return blank
    decline = float(rate[-1] / rate[0])
    if decline > SINK_DECLINE:
        return blank
    slope, intercept, rate_r2 = _line_r2(made, rate)
    inverse_slope, inverse_intercept, reciprocal_r2 = _line_r2(made, 1.0 / rate)
    k = -slope
    return SinkFit(
        v0=intercept, k=float(k),
        plateau=float(intercept / k) if k > 0 else np.nan,
        ki=(float(inverse_intercept / inverse_slope)
            if inverse_slope != 0 else np.nan),
        rate_r2=rate_r2, reciprocal_r2=reciprocal_r2,
        decline=decline, points=int(len(rate)))


def sink_table(experiments, fraction=SINK_WINDOW):
    """One row per curve that declines far enough to be read this way."""
    import pandas as pd
    rows = []
    for curve in scope.curves(tuple(experiments)):
        fitted = sink_fit(curve, fraction=fraction)
        rows.append({
            "experiment": curve.experiment, "sample": curve.sample,
            "substrate": curve.substrate, "buffer": curve.buffer,
            "pH": curve.pH, "temperature": curve.temperature,
            "s0": curve.conditions.s0, "e0": curve.conditions.e0,
            "h2o2": curve.conditions.h2o2, "buf": curve.buf,
            "epsilon": curve.epsilon,
            "span_s": float(curve.times[-1] - curve.times[0]),
            "net": float(curve.absorbance[-1] - curve.absorbance[0]),
            "v0": fitted.v0, "k_sink": fitted.k, "plateau": fitted.plateau,
            "ki": fitted.ki, "rate_r2": fitted.rate_r2,
            "reciprocal_r2": fitted.reciprocal_r2, "decline": fitted.decline,
            "prefers": fitted.prefers, "points": fitted.points,
        })
    return pd.DataFrame(rows)


def selectivity(plateau, s0, epsilon):
    """
    How much faster the oxidant attacks the product than the substrate.

    If the measured species is made from the substrate and destroyed by the
    same oxidant, the steady state is where the two are equal:

        k_S [S] [Ox]  =  k_A [A] [Ox]      ->     k_A/k_S = [S] / [A]_inf

    and [A]_inf is the plateau in mM. So one number falls out of the plateau
    that a calculation can be checked against: the barrier difference between
    oxidising the alcohol and oxidising its own aldehyde, in the same medium,
    with the concentration term divided out.

    It is an ESTIMATE, not a measurement, for one reason that cannot be fixed
    from these runs: the plateau is read off the catalytic INCREMENT, and the
    aldehyde in the cuvette is the increment plus whatever the enzyme-free
    background made. So this is an upper bound on k_A/k_S, and the correction
    goes the same way for every curve.
    """
    if not (plateau > 0 and s0 > 0 and epsilon > 0):
        return np.nan
    return float(s0 / (plateau / epsilon))


# A curve has to have made enough signal for "how much product" to mean
# something. 0.005 AU is about 5x the largest baseline drift in the archive
# and 0.0007 mM of 4OMe product; below it the ratio below is a ratio of noise.
DRIVER_FLOOR = 0.005
# Where each hypothesis puts the two coefficients of the driver regression.
DRIVER_CLOCK = ("span", "product")


def deceleration_drivers(frame, fixed=False, extra=()):
    """
    Regress how much a curve slowed on how LONG it ran and how MUCH it made.

        log(late rate / early rate) = a . log(span) + b . log(product) + c

    This is the discrimination the shapes cannot make. A single progress curve
    cannot separate "the rate fell with time" from "the rate fell with
    product", because within one curve the product only ever grows with time.
    Across curves the two come apart -- the archive holds 1470 s runs that
    reached 0.27 AU and 17934 s runs that reached 0.045 -- and in the catalysed
    4OMe phosphate set the two regressors correlate only -0.39.

        a < 0, b = 0     a clock: something decays on its own schedule
        a = 0, b < 0     product control: the rate is set by what has been made
        a = b = 0        neither, over the range the design covers

    `late_over_early` is `scope`'s own statistic (last fifth against first
    fifth) and is not recomputed here. Product is `net/epsilon` in mM, so the
    two substrates are compared as CONCENTRATIONS and not as absorbances --
    they differ by a factor of 6.1 in epsilon, and on absorbance the BnOH set
    would look like it reached a third of what it reaches.

    `fixed=True` adds one dummy per experiment, which absorbs temperature, pH,
    buffer, enzyme, run length and cell, and leaves only the substrate ladder
    inside each run to carry the product term. The span term is dropped there:
    every cuvette of an experiment shares a run length, so it is collinear with
    the dummies and its coefficient in that model means nothing.

    Beware the direction of the one bias that cannot be removed: `net` is the
    integral of the rate, so a curve that decelerates less makes MORE product
    at the same starting rate. That pushes b UPWARDS, towards zero. A negative
    b is therefore a floor on the effect, not a ceiling.
    """
    live = frame[frame.live & (frame.late_over_early > 0)
                 & (frame.net > DRIVER_FLOOR)].copy()
    if len(live) < 6:
        return {"points": len(live)}
    columns = {"product": np.log(live.net.to_numpy() / live.epsilon.to_numpy())}
    if not fixed:
        columns["span"] = np.log(live.duration_s.to_numpy())
    for name in extra:
        columns[name] = live[name].to_numpy(dtype=float)
    design = [np.ones(len(live))] + [columns[n] for n in columns]
    names = list(columns)
    if fixed:
        for experiment in sorted(live.experiment.unique())[1:]:
            design.append((live.experiment == experiment).to_numpy(float))
            names.append(f"experiment_{experiment}")
    design = np.column_stack(design)
    response = np.log(live.late_over_early.to_numpy())
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    resid = response - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(live) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    spread = float(((response - response.mean()) ** 2).sum())
    out = {"points": int(len(live)),
           "r2": float(1 - float(resid @ resid) / spread) if spread else np.nan,
           "collinearity": float(np.corrcoef(columns["product"],
                                             columns["span"])[0, 1])
           if "span" in columns else np.nan}
    for index, name in enumerate(names):
        if name.startswith("experiment_"):
            continue
        out[name] = float(beta[index + 1])
        out[name + "_stderr"] = float(np.sqrt(max(covariance[index + 1,
                                                             index + 1], 0.0)))
    return out


def plateau_scaling(table):
    """
    Regress log(plateau) on log([S]), and compare with the substrate order.

    If the signal is a species the oxidant makes from the substrate and then
    destroys, the stationary level is `A_inf = v(S)/k`, so the plateau has to
    carry the SAME substrate order the rate does -- about +0.58 here, not +1,
    because the production is partly saturated. That is a prediction with
    nothing fitted to it: the substrate order was measured on the rates, in a
    different block, before any of this.
    """
    live = table[(table.points > 0) & (table.plateau > 0) & (table.s0 > 0)]
    if len(live) < 4:
        return {"points": len(live), "order": np.nan, "stderr": np.nan}
    x = np.log(live.s0.to_numpy())
    y = np.log(live.plateau.to_numpy())
    design = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    variance = float(resid @ resid) / max(1, len(x) - 2)
    covariance = variance * np.linalg.pinv(design.T @ design)
    return {"points": int(len(live)), "order": float(beta[1]),
            "stderr": float(np.sqrt(covariance[1, 1])),
            "lever": float(np.exp(x.max() - x.min()))}


def product_budget(frame):
    """
    How much product exists, against how much catalyst there is to inhibit.

    Inhibition of the catalyst has an arithmetic to answer before it has a
    mechanism. Blocking a fraction f of the catalyst takes f.[enz] of
    inhibitor BOUND, whatever the binding constant -- tightening the binding
    moves the equilibrium, it does not create inhibitor. So a run that loses
    half its rate to product inhibition needs product of order half the
    catalyst concentration, and if it has a fiftieth, the explanation is not
    inhibition of the catalyst.

    Returns the product made (mM), as a fraction of [enz], per experiment.
    """
    import pandas as pd
    live = frame[frame.live & (frame.e0 > 0)].copy()
    live["product_mM"] = live.net / live.epsilon
    live["of_enzyme"] = live.product_mM / live.e0
    return (live.groupby(["substrate", "experiment"])
            .agg(temperature=("temperature", "first"), pH=("pH", "first"),
                 e0=("e0", "first"), product_mM=("product_mM", "max"),
                 of_enzyme=("of_enzyme", "max"),
                 late_over_early=("late_over_early", "median")))


def substrate_blocks(frame=None):
    """
    The archive cut the way this question needs it, from the frame's own
    columns rather than from a list of experiment numbers.

    `differential` is the structural classification -- what the reference
    channel omitted -- not the filename, which is why exps 32 and 34-37 sit in
    the catalysed blocks despite being called `with_NO_E`.
    """
    if frame is None:
        frame = scope.frame(tuple(range(1, 152)))
    four = frame.substrate == "4OMe-BnOH"
    bnoh = frame.substrate == "BnOH"
    phosphate = frame.buffer == "Phosphate"
    return {
        "4OMe catalysed, phosphate": frame[four & frame.differential & phosphate],
        "4OMe catalysed, all buffers": frame[four & frame.differential],
        "4OMe enzyme-free": frame[four & ~frame.differential],
        "BnOH catalysed, all buffers": frame[bnoh & frame.differential],
        "BnOH in scope (135-151)":
            frame[frame.experiment.isin(scope.PRIMARY_SCOPE)],
        "BnOH enzyme-free": frame[bnoh & ~frame.differential],
        "temperature series":
            frame[frame.experiment.isin(scope.TEMPERATURE_SERIES)],
    }


def report(frame=None):
    """Print the whole argument, in the order it has to be made."""
    if frame is None:
        frame = scope.frame(tuple(range(1, 152)))
    named = substrate_blocks(frame)

    print("\n1. WHAT THE SLOWDOWN TRACKS")
    print("   log(late rate / early rate) on log(run length) and "
          "log(product, mM)")
    print(f"   {'block':30s} {'n':>4s} {'R2':>6s} {'rho':>6s} "
          f"{'span':>16s} {'product':>16s}")
    for name, block in named.items():
        row = deceleration_drivers(block)
        if row.get("points", 0) < 6:
            continue
        print(f"   {name:30s} {row['points']:4d} {row['r2']:6.3f} "
              f"{row['collinearity']:+6.2f} "
              f"{row['span']:+8.3f} +/- {row['span_stderr']:.3f} "
              f"{row['product']:+8.3f} +/- {row['product_stderr']:.3f}")
    print("\n   with every condition absorbed by one dummy per experiment, so"
          "\n   only the substrate ladder inside each run carries the product:")
    for name in ("4OMe catalysed, phosphate", "4OMe enzyme-free",
                 "BnOH catalysed, all buffers", "BnOH in scope (135-151)"):
        row = deceleration_drivers(named[name], fixed=True)
        print(f"   {name:30s} {row['points']:4d} {row['r2']:6.3f} "
              f"{'':6s} {'':16s} {row['product']:+8.3f} "
              f"+/- {row['product_stderr']:.3f}")

    print("\n2. WHAT SHAPE THE FALL HAS")
    catalysed = named["4OMe catalysed, phosphate"]
    table = sink_table(sorted(catalysed.experiment.unique()))
    good = table[(table.points > 0) & (table.rate_r2 > SINK_CLEAN_R2)]
    counts = good.prefers.value_counts().to_dict()
    print(f"   {len(good)} catalysed 4OMe phosphate curves are fitted well "
          f"enough to choose on (R2 > {SINK_CLEAN_R2})")
    print(f"   rate linear in product (a sink):    {counts.get('sink', 0)}")
    print(f"   1/rate linear in product (inhibition): "
          f"{counts.get('inhibition', 0)}")
    print(f"   tied within {SLOWDOWN_MARGIN:.0%}: {counts.get('tied', 0)}")
    print(f"   median R2   rate {good.rate_r2.median():.4f}   "
          f"reciprocal {good.reciprocal_r2.median():.4f}")

    print("\n3. WHAT THE SINK LAW PREDICTS, AND WHETHER IT HOLDS")
    order = plateau_scaling(good)
    print(f"   A_inf = v(S)/k, so the plateau must carry the rate's own "
          f"substrate order")
    print(f"   plateau order  {order['order']:+.3f} +/- {order['stderr']:.3f} "
          f"over {order['lever']:.0f}x in [S], {order['points']} curves")
    from arrhenius import substrate_order                # local: heavy import
    measured = float(substrate_order().corrected.mean())
    print(f"   rate's substrate order, measured on the rates of the "
          f"temperature series and not fitted here: {measured:+.3f}")
    selectivities = np.array([selectivity(p, s, e) for p, s, e
                              in zip(good.plateau, good.s0, good.epsilon)])
    selectivities = selectivities[np.isfinite(selectivities)]
    print(f"   k_A/k_S = [S]/[A]_inf   median {np.median(selectivities):.0f}"
          f"  IQR {np.percentile(selectivities, 25):.0f}-"
          f"{np.percentile(selectivities, 75):.0f}   (an UPPER bound)")

    print("\n4. WHETHER THERE IS ENOUGH PRODUCT TO INHIBIT THE CATALYST")
    budget = product_budget(frame)
    series = budget.loc["4OMe-BnOH"].reindex(scope.TEMPERATURE_SERIES)
    for experiment, row in series.iterrows():
        print(f"   exp {experiment:3d}  {row.temperature:4.0f} C   product "
              f"{row.product_mM * 1000:6.1f} uM   [enz] {row.e0 * 1000:5.0f} uM"
              f"   {row.of_enzyme:6.1%} of it   late/early "
              f"{row.late_over_early:.2f}")


def main():
    report()


if __name__ == "__main__":
    main()
