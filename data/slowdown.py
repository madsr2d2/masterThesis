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

A ROUTE THAT WAS TRIED AND WITHDRAWN. The obvious statistic is a landmark --
when, and at what product, does the rate fall to three quarters of its peak --
regressed against the peak rate: product control makes that time go as 1/v and
a clock makes it flat. It does not work, and the reason is censoring. A curve
whose rate never falls that far inside its run has no landmark, slow curves run
out of time before fast ones do, and dropping them biases the slope towards the
clock's prediction by an amount that cannot be estimated from the survivors.
`decay_point` remains as a diagnostic; nothing in this module's conclusions
rests on it. What replaced it is `deceleration_drivers`, which is defined for
every curve and so cannot be censored at all.

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


def fit_slowdown(times, values, name, points=SLOWDOWN_GRID, decay=None):
    """
    Fit one of the three mechanistic forms to one curve.

    `name` is "clock", "sink" or "inhibition". Returns a SlowdownFit whose
    `sse` is directly comparable across the three: same data, same parameter
    count, same grid resolution.

    `decay` PINS the second parameter instead of profiling it, which turns a
    four-parameter fit into a three-parameter one and removes the degeneracy
    that makes the free fit useless on short runs: with k free, a curve that
    has not yet turned over is fitted equally well by any (v, k) pair holding
    v/k at the plateau, and three of these 24 curves come back with k a
    hundredfold out. Supplying k from `sink_activation` costs the fit nothing
    it could have determined for itself.
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
    if decay is not None:
        second = np.array([float(decay)])
    elif name == "inhibition":
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


# The rate is read as the slope over this fraction of the run, which is the
# resolution at which these curves have a shape at all.
DECAY_WINDOW = 0.05

# Where the two hypotheses put the slope of log(time to fall) on log(peak
# rate). Product control fixes the PRODUCT at which the rate falls, so the
# time to reach it is inversely proportional to the rate; a clock fixes the
# TIME and lets the product land where it may.
PRODUCT_SLOPE = -1.0
CLOCK_SLOPE = 0.0


def _windowed_rate(curve, samples=4000):
    """
    The chosen fitted form's shape, sampled: (t, A(t), rate(t)).

    Read off `summary_kinetics.fit_progress`'s CHOSEN form, not off the
    readings: the readings' own slope is too noisy to locate a landmark on,
    and the chosen form is already the one the F test says the curve earned.
    Evaluated only inside the observed span -- anything past the last reading
    is a property of the model, not of the experiment.

    A WINDOWED slope, not the analytic derivative. The fitted form is free to
    put a spike at t = 0 -- a two-phase fit with a short tau1 and B1 < 0 does,
    and its instantaneous rate there can be ten times anything the readings
    show -- and a landmark measured against that spike reports a "fall" on
    curves that accelerate throughout. Averaging the fit over the same window
    the block-slope tables use removes it, and costs nothing on a curve whose
    fall is real and slow.
    """
    from summary_kinetics import fit_progress          # local: heavy import
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    shifted = times - times[0]
    fitted = fit_progress(times, values)
    grid = np.linspace(0.0, float(shifted[-1]), samples)
    predicted = fitted.predict(grid + times[0])
    if not np.all(np.isfinite(predicted)):
        return None
    step = max(1, int(round(samples * DECAY_WINDOW)))
    rate = (predicted[step:] - predicted[:-step]) / (grid[step:] - grid[:-step])
    return grid[:-step], predicted[:-step], rate


def decay_point(curve, fraction=0.5, samples=4000):
    """
    Where the fitted rate falls to `fraction` of its maximum: (t, A, v, t_max).

    A curve whose rate never falls that far returns NaN and is counted as
    censored rather than silently given the run's end.
    """
    sampled = _windowed_rate(curve, samples=samples)
    if sampled is None:
        return (np.nan,) * 4
    grid, predicted, rate = sampled
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


# The rolling window the sink linearisation reads its rate from. Wider than
# curve_metrics.LAG_WINDOW because the question here is the SLOPE of the rate
# against product over the whole tail, not where an induction ends, and a
# tenth of the run leaves the tail of a short curve with too few independent
# windows to regress.
SINK_WINDOW = 0.15
SINK_MINIMUM_POINTS = 8
# HOW MUCH INDEPENDENT INFORMATION THE TAIL HAS TO CONTAIN, as a multiple of
# the window's own width. `SINK_MINIMUM_POINTS` counts window POSITIONS, and
# consecutive positions overlap by all but one reading, so it can be satisfied
# by a tail shorter than a single window -- a straight line drawn through less
# than one independent measurement of the rate.
#
# That is not hypothetical. At a window of 0.30 of the run the 15 and 20 C
# cuvettes are admitted with tails of 0.32, 0.47 and 0.54 windows, they return
# a k at 15 C larger than the k at 25 C, and the sink's activation energy
# collapses from 72.3 +- 10.0 to 7.9 +- 33.0 -- which then breaks
# `sink_effect_on_activation`'s null, the one this package's section 5 rests
# on. 1.5 keeps every curve the published window admits (the tightest is 2.08)
# and rejects all three of those. `sink_window_sensitivity` prints the sweep.
SINK_MINIMUM_WINDOWS = 1.5
# A curve has to be described by the straight line before which straight line
# it is can be argued about. 0.95 keeps the curves whose decline is deep
# enough to have a shape and drops the ones whose "decline" is a slope through
# the noise of a nearly flat tail.
SINK_CLEAN_R2 = 0.95
# The tail has to have fallen this far below the maximum before "the rate
# declines with product" is a measurement rather than a slope through noise.
SINK_DECLINE = 0.85
# ...but only when the question is WHICH FORM. To read k off the line, a
# shallow decline is enough as long as the line itself is good, and the R2
# gate already asks that. The two thresholds answer different questions and
# the 25 C run sits between them: it declines to 0.87-0.97, too little to
# choose a functional form on and quite enough to fit a slope to.
SINK_SLOPE_DECLINE = 0.98


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
    windows: float
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


def sink_fit(curve, fraction=SINK_WINDOW, decline=SINK_DECLINE):
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
    blank = SinkFit(*([np.nan] * 8), 0)
    if len(slopes) < SINK_MINIMUM_POINTS:
        return blank
    product = np.interp(centres, times, values) - values[0]
    top = int(np.argmax(slopes))
    rate, made = slopes[top:], product[top:]
    if len(rate) < SINK_MINIMUM_POINTS or not np.all(rate > 0):
        return blank
    # The tail in units of the window it is read through. See
    # SINK_MINIMUM_WINDOWS: the point count above cannot catch this, because
    # window positions overlap.
    width = fraction * float(times[-1] - times[0])
    independent = (float(centres[-1] - centres[top]) / width
                   if width > 0 else np.nan)
    if not independent >= SINK_MINIMUM_WINDOWS:
        return blank
    fallen = float(rate[-1] / rate[0])
    if fallen > decline:
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
        decline=fallen, windows=float(independent), points=int(len(rate)))


def sink_table(experiments, fraction=SINK_WINDOW, decline=SINK_DECLINE):
    """One row per curve that declines far enough to be read this way."""
    import pandas as pd
    rows = []
    for curve in scope.curves(tuple(experiments)):
        fitted = sink_fit(curve, fraction=fraction, decline=decline)
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
            "windows": fitted.windows,
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

    print("\n5. WHAT IT DOES TO THE ACTIVATION PARAMETERS")
    effect = sink_effect_on_activation()
    sweep = sink_window_sensitivity()
    print("\n   what the sink constant owes to its window:")
    print(sweep.to_string(index=False))
    print(f"   activation spread across windows "
          f"{sweep.attrs['activation_spread']:.1f} kJ/mol; "
          f"the null holds at every width: {sweep.attrs['null_holds']}")
    print(f"   the sink's own Ea  {effect['sink_activation_kJ']:.1f} +/- "
          f"{effect['sink_stderr_kJ']:.1f} kJ/mol, on "
          + ", ".join(f"{t - 273.15:.0f}" for t in effect["sink_temperatures"])
          + " C")
    print(f"   {'estimator':22s} {'Ea kJ/mol':>16s} {'dS J/mol/K':>16s} "
          f"{'rms':>7s}")
    for label, key in (("v_peak, published", "published"),
                       ("v_prod, k pinned", "corrected")):
        row = effect[key]
        print(f"   {label:22s} {row['activation_kJ']:8.2f} +/- "
              f"{row['activation_stderr']:.2f} {row['entropy_J']:9.1f} +/- "
              f"{row['entropy_stderr']:.1f} {row['rms']:7.3f}")
    print(f"   the production rate sits {effect['lift'] - 1:+.1%} above "
          f"v_peak on the median curve,")
    print(f"   between {effect['lift_low'] - 1:+.1%} and "
          f"{effect['lift_high'] - 1:+.1%} across the six temperatures, with "
          f"no order to it,")
    print(f"   so Ea moves by {effect['activation_shift']:+.2f} +/- "
          f"{effect['activation_shift_stderr']:.2f} kJ/mol -- nothing.")


# The sink's rate constant is read off the line in `sink_fit`, so a curve
# qualifies when that LINE is good, not when the decline is deep. 0.85 is
# looser than SINK_CLEAN_R2 because this needs the coldest runs it can get and
# they decline least; anything looser lets a slope through noise in.
SINK_ARRHENIUS_R2 = 0.85
SINK_ARRHENIUS_TEMPERATURES = 3


def sink_constants(experiments, floor=SINK_ARRHENIUS_R2,
                   decline=SINK_SLOPE_DECLINE):
    """
    Per-curve sink rate constant k, for curves where the line is trustworthy.

    Returns a DataFrame with temperature, kelvin, k and the line's R². Curves
    whose rate never turns over -- every 15 and 20 °C cuvette here -- have no
    post-maximum tail to regress and are absent, which is the honest outcome
    and the reason `sink_activation` has to extrapolate to reach them.
    """
    table = sink_table(experiments, decline=decline)
    live = table[(table.points > 0) & (table.k_sink > 0)
                 & (table.rate_r2 > floor)].copy()
    live["kelvin"] = live.temperature + 273.15
    return live


def sink_activation(experiments, floor=SINK_ARRHENIUS_R2):
    """
    An Arrhenius fit to the sink's own rate constant, and a predictor for k(T).

    Fitted on the per-temperature MEDIAN of k rather than on every curve,
    because the four cuvettes of a run share a cell, a day and a peroxide
    aliquot and are not four independent draws -- the same argument
    `pooled_arrhenius` makes about the four rungs.

    Returns the `arrhenius.arrhenius_fit` dict with an extra `predict`
    callable, or None if fewer than three temperatures survive. Three is the
    minimum a line can be fitted to and it is exactly what this block gives:
    treat the activation energy as an order of magnitude and the predictor as
    an interpolation, which is all `production_frame` asks of it.
    """
    from arrhenius import arrhenius_fit                # local: heavy import
    live = sink_constants(experiments, floor=floor)
    if live.empty:
        return None
    grouped = live.groupby("kelvin").k_sink.median()
    if len(grouped) < SINK_ARRHENIUS_TEMPERATURES:
        return None
    fit = arrhenius_fit(grouped.index.to_numpy(dtype=float),
                        grouped.to_numpy(dtype=float))
    if fit is None:
        return None
    slope, intercept = fit["slope"], fit["intercept"]
    fit["predict"] = lambda kelvin: np.exp(intercept + slope / np.asarray(
        kelvin, dtype=float))
    fit["temperatures"] = [float(k) for k in grouped.index]
    fit["curves"] = int(len(live))
    return fit


# The windows `sink_window_sensitivity` sweeps. 0.15 is what everything here is
# quoted at; the others are there to put a systematic on it.
SINK_WINDOW_SWEEP = (0.15, 0.20, 0.25, 0.30)


def sink_window_sensitivity(experiments=None, widths=SINK_WINDOW_SWEEP,
                            floor=SINK_ARRHENIUS_R2):
    """
    What the sink's rate constant, and everything pinned to it, owe to the window.

    `k` is read as the slope of the rolling rate against the product already
    made, and that rolling rate comes through a window that is a FRACTION of
    the run -- so a comparison of k across runs of 1470 and 17934 s is a
    comparison through windows that differ twelvefold. This sweeps the width
    and reports what moves.

    WHAT IT FOUND, and the reason it exists. Two things, and only one of them
    was a defect:

      the sink's own activation energy is WINDOW-DEPENDENT, rising 72 -> 88 ->
      95 -> 102 kJ/mol across 0.15 to 0.30 because a wider window smooths the
      slow cold curves more than the fast warm ones and flattens their
      rate-against-product slope. That systematic is larger than the +- 10.0
      statistical error and belongs beside it.

      the NULL that section 5 rests on -- that naming the fall does not move
      the activation parameters -- holds at every width, and holds more
      strongly as the window widens. That is the load-bearing result and it is
      not affected.

    Before `SINK_MINIMUM_WINDOWS` existed a width of 0.30 also admitted three
    cold cuvettes whose tail was shorter than a single window, returned a k at
    15 C larger than the k at 25 C, and broke the null. That is a guard now
    rather than a caveat.
    """
    import pandas as pd
    experiments = (scope.TEMPERATURE_SERIES if experiments is None
                   else experiments)
    original = sink_table.__defaults__
    rows = []
    try:
        for width in widths:
            sink_table.__defaults__ = (width, original[1])
            live = sink_constants(experiments, floor=floor)
            row = {"window": float(width), "curves": int(len(live)),
                   "temperatures": int(live.kelvin.nunique())}
            try:
                fit = sink_activation(experiments, floor=floor)
                row["activation_kJ"] = fit["activation_kJ"]
                row["stderr_kJ"] = fit["stderr_kJ"]
            except ValueError:
                row["activation_kJ"] = np.nan
                row["stderr_kJ"] = np.nan
            try:
                effect = sink_effect_on_activation(experiments, floor=floor)
                row["shift_kJ"] = effect["activation_shift"]
                row["shift_stderr"] = effect["activation_shift_stderr"]
            except ValueError:
                row["shift_kJ"] = np.nan
                row["shift_stderr"] = np.nan
            rows.append(row)
    finally:
        sink_table.__defaults__ = original
    table = pd.DataFrame(rows)
    good = table[np.isfinite(table.activation_kJ)]
    table.attrs["activation_spread"] = (float(good.activation_kJ.max()
                                              - good.activation_kJ.min())
                                        if len(good) else np.nan)
    table.attrs["null_holds"] = bool(
        (table.shift_kJ.abs() < table.shift_stderr).all())
    return table


def production_frame(experiments=None, floor=SINK_ARRHENIUS_R2):
    """
    The temperature series with the sink taken back out of the rate.

    WHY THE PUBLISHED RATE IS TOO LOW. If the signal obeys `A' = v - kA` then
    every rate the instrument reports is already net of the loss, and the
    production rate `v` -- the quantity an activation energy is wanted for --
    is larger than anything read off the curve.

    HOW IT IS RECOVERED, and why not the obvious way. The obvious way is to
    add `k.A` back at the moment the rate peaks. It does not work: on the four
    coldest runs the fitted rate peaks at the END of the run rather than
    inside it, so that recipe silently swaps `v_peak` -- which is an
    asymptote, and truncation-free by construction (see ANALYSIS.md 3a) -- for
    a value read at the last reading, which reintroduces exactly the cold-end
    truncation `v_peak` exists to avoid. It moved the activation energy by
    +2.4 kJ/mol of pure bookkeeping.

    What is done instead is to fit the sink model itself,

        A' = v (1 - exp(-t/tau)) - k A

    with **k pinned** at `sink_activation`'s value for that temperature. `v`
    is then a parameter of the model rather than a reconstruction, and pinning
    k removes the degeneracy that makes the free four-parameter fit useless
    here: a curve that has not turned over is fitted equally well by any
    (v, k) holding v/k at the plateau, and three of these 24 come back with k
    a hundredfold out.

    Adds `sink_k`, `sink_tau`, `sink_plateau` and `v_prod`.
    """
    from arrhenius import series_frame, TEMPERATURE_SERIES   # local
    experiments = TEMPERATURE_SERIES if experiments is None else experiments
    frame = series_frame(experiments)
    fit = sink_activation(experiments, floor=floor)
    if fit is None:
        raise ValueError("no sink Arrhenius: fewer than "
                         f"{SINK_ARRHENIUS_TEMPERATURES} temperatures give a "
                         "readable rate-against-product line")
    frame = frame.copy()
    frame["sink_k"] = fit["predict"](frame.kelvin.to_numpy(dtype=float))
    rates, taus, plateaus = [], [], []
    for curve, (_, row) in zip(scope.curves(tuple(experiments)),
                               frame.iterrows()):
        fitted = fit_slowdown(np.asarray(curve.times, dtype=float),
                              np.asarray(curve.absorbance, dtype=float),
                              "sink", decay=float(row.sink_k))
        rates.append(fitted.extra.get("rate", np.nan))
        taus.append(fitted.induction)
        plateaus.append(fitted.extra.get("plateau", np.nan))
    frame["v_prod"] = rates
    frame["sink_tau"] = taus
    frame["sink_plateau"] = plateaus
    frame.attrs["sink_activation_kJ"] = fit["activation_kJ"]
    frame.attrs["sink_stderr_kJ"] = fit["stderr_kJ"]
    frame.attrs["sink_temperatures"] = fit["temperatures"]
    frame.attrs["sink_curves"] = fit["curves"]
    return frame


def sink_effect_on_activation(experiments=None, floor=SINK_ARRHENIUS_R2):
    """
    Does naming the fall change the activation parameters? Answer: no.

    Returns a dict comparing `v_peak` -- the published estimator, the maximum
    of a DESCRIPTIVE fit and therefore net of the sink -- with `v_prod`, the
    production rate of the sink model itself. The two differ by a factor that
    is systematically greater than one and only weakly ordered in temperature,
    and **a factor that does not order in temperature cancels out of a slope**.
    So the correction lands on the LEVEL: it raises every rate a few per cent,
    which moves the Eyring entropy and the absolute free energy a little and
    the enthalpy not at all.

    Reported rather than applied. `v_prod` is the more nearly correct
    quantity and the noisier one -- the sink model imposes a shape the four
    coldest runs cannot test, and its Arrhenius scatter is twice `v_peak`'s --
    so quoting it as the headline would trade a bias smaller than the error
    for a variance larger than it.
    """
    from arrhenius import activation_parameters, GAS_CONSTANT   # local
    frame = production_frame(experiments, floor=floor)
    published = activation_parameters("v_peak", frame=frame)
    corrected = activation_parameters("v_prod", frame=frame)
    ratio = (frame.v_prod / frame.v_peak).to_numpy(dtype=float)
    ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
    by_temperature = (frame.v_prod / frame.v_peak).groupby(
        frame.temperature).median()
    lift = float(np.median(ratio))
    return {
        "published": published, "corrected": corrected,
        "lift": lift,
        "lift_low": float(np.min(by_temperature)),
        "lift_high": float(np.max(by_temperature)),
        "by_temperature": by_temperature,
        "activation_shift": float(corrected["activation_kJ"]
                                  - published["activation_kJ"]),
        "activation_shift_stderr": float(np.hypot(
            corrected["activation_stderr"], published["activation_stderr"])),
        "entropy_from_lift": float(GAS_CONSTANT * np.log(lift)),
        "sink_activation_kJ": frame.attrs["sink_activation_kJ"],
        "sink_stderr_kJ": frame.attrs["sink_stderr_kJ"],
        "sink_temperatures": frame.attrs["sink_temperatures"],
    }


def main():
    report()


if __name__ == "__main__":
    main()


