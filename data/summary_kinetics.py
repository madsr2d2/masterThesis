"""
Model-free summaries of the progress curves, and the scaling laws they obey.

This is the layer beneath `fit_kinetics.py`. It fits no mechanism. It reduces
each curve to two numbers and then asks how those numbers scale with what was
in the cuvette, which is a question the data can answer without committing to
an ODE.

Three facts about this dataset make that the right first move:

  * Conversion is under 1% on most curves, so a progress curve is essentially
    its own initial-rate regime. The information an ODE fit normally extracts
    from a progress curve lives in the curvature that appears as species
    deplete, and there is almost none of that here.
  * A straight line over the first 20% of a run fits to 1.10x the curve's own
    point noise. Over the whole run it fits to 5.85x. So the early window is
    described down to the noise floor, and everything an ODE could tell us
    about that window is already in its slope.
  * Pooling experiments with one shared amplitude -- which is what a fit over a
    block does -- leaves 80% of the variance in log v0 unexplained. Giving each
    experiment its own free amplitude drops the residual sum of squares 5.1x.
    That between-experiment offset is not kinetics, and until it is modelled it
    will dominate any residual, however right the chemistry is.

WHY A STRAIGHT LINE, AND WHY 20%

Chosen by measurement, not taste. Polynomials of degree 1-3 over windows of
10-100% were judged on whether they fit to the noise floor, how precise the
slope was, and how well it reproduced across the four identical cuvettes of
experiment 26:

    deg  window  rms/noise  v0 rel.SE  exp26 RSD
      1     10%       0.93       2.4%       9.1%
      1     20%       1.10       2.0%       9.1%   <-- used
      1     50%       1.89       1.2%       9.1%
      1    100%       5.86       1.3%       3.8%
      2     20%       0.67       5.3%      21.5%
      3     20%       0.53       8.7%      44.9%

Degree 1 at 20% sits at rms/noise = 1: the line is an adequate description of
that window and no more. Higher degrees fit BELOW the noise floor on short
windows and multiply the replicate scatter by 2.4 and 4.9, because the
curvature term takes its variance out of the slope. A full-run line is
reproducible but is not measuring an initial rate at all.

The replicate column is flat from 10% to 50% for a mundane reason: experiment
26's curves are only ten points long, so `window_size`'s five-point floor binds
at every one of those fractions and the same five points are used throughout.
That set therefore pins the reproducibility of the METHOD, not of the window
choice, and the window is justified by the rms and standard-error columns
alone. Reproduce this table with `python data/summary_kinetics.py --group
4OMe-BnOH,40`, whose window-sensitivity section prints the same rms column.

The obvious alternative -- the textbook induction form
`A = c + v_ss t - B(1 - exp(-t/tau))` -- was tried and rejected. It fits to 1.1x
noise, but it is not identifiable on these runs: for tau near or beyond the run
length the exponential IS the linear term, and the unconstrained fit exploits
that (median burst amplitude 2.1x the entire observed signal). Capping tau at a
third of the run pins 69% of curves at the cap, and the "steady rate" it returns
gives a substrate order of -2.64 inside a single experiment. There is no plateau
at 1% conversion to anchor it against.

WINDOW SENSITIVITY IS A REAL SYSTEMATIC

The substrate order measured this way moves from +0.13 at a 10% window to +0.43
at 100% -- a 3.3x swing from an analysis choice, against a +/-0.07 statistical
error on any one of them. The curves are curved and the curvature is itself
[S]-dependent, so an order is only meaningful with its window attached.
`window_sensitivity()` measures this, and `report_windows()` prints it; treat
its spread as the systematic error on any order quoted from here.

    python data/summary_kinetics.py                        # enzyme-free, per group
    python data/summary_kinetics.py --group 4OMe-BnOH,40   # one block, in full
    python data/summary_kinetics.py --catalysed
    python data/summary_kinetics.py --save summaries.csv
"""
import json
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from curve_metrics import (ACCELERATION_SIGMA, INITIAL_WINDOW, LAG_THRESHOLD,
                           MINIMUM_WINDOW_POINTS, QUANTISATION_SIGMA,
                           acceleration, initial_rate, line_fit,
                           line_slope, peak_position, window_size)
from fit_dataset import (ABSORBANCE_QUANTUM, DATASET_PATH,
                         QUANTISATION_SIGMA, build_curves)
from solution_chemistry import add_solution_columns

# Fraction of each run the rate is measured over. See the module docstring for
# how this was chosen; it is not a free knob, and moving it moves every order.
# INITIAL_WINDOW, MINIMUM_WINDOW_POINTS: imported from curve_metrics above.

# A window shorter than this cannot support a slope worth reporting.


# Curve-to-curve scatter in v0 at IDENTICAL conditions, from the four cuvettes
# of experiment 26 -- the only true replicate set in the enzyme-free data. The
# slope's own standard error is 2.0%, so reproducibility is the binding
# constraint by a factor of four and this, not the fit error, is what a
# regression residual should be judged against.
#
# Weak evidence, and it should be treated as such: one experiment, four
# cuvettes, ten points each, at a single [S]. It is a floor on the achievable
# residual, not a calibrated error bar. `replicate_scatter()` recomputes it on
# whatever subset is in hand and `report_replicates()` prints both.
REPLICATE_RSD = 0.091

# The axes a rate is regressed against, all in log10. [H2O2] is deliberately
# absent: it is constant at 82.5 mM across all 59 enzyme-free curves of the
# largest block and varies 2x over three levels in the only block that moves it
# at all, so its order is not measurable and asking for it only adds a column
# of near-zero variance.
DESIGN_AXES = ("s0", "hoo", "buf")

# chi-square, 1 degree of freedom, 95% -- the profile-likelihood threshold for
# the interval on Km.
CHI2_95 = 3.84


@dataclass(frozen=True)
class CurveSummary:
    """One cuvette reduced to what the design can actually resolve."""
    experiment: int
    sample: int
    substrate: str
    temperature: float
    buffer_name: str
    pH: float
    s0: float            # mM substrate
    h2o2: float          # mM peroxide
    hoo: float           # mM hydroperoxide anion
    buf: float           # mM buffer
    e0: float            # mM catalyst
    v0: float            # AU/s, slope over the first INITIAL_WINDOW of the run
    v0_stderr: float     # AU/s, from the line fit alone -- NOT reproducibility
    window_rms: float    # AU, residual of that line
    window_quanta: float # how far that window rose, in 0.001 AU steps
    noise: float         # AU, the curve's own point noise
    slope_ratio: float   # v(last window) / v(first window); <1 burst, >1 lag
    points: int
    duration: float      # s
    amplitude: float     # AU
    conversion: float    # fraction of substrate turned over, if all signal is A

    @property
    def group(self):
        return (self.substrate, self.temperature, self.buffer_name)

    @property
    def window_quality(self):
        """Residual of the straight line in units of the curve's point noise."""
        return self.window_rms / self.noise if self.noise > 0 else np.nan


@dataclass(frozen=True)
class Regression:
    """One weighted least-squares fit of log10(response) against the axes."""
    response: str
    names: tuple
    coefficients: np.ndarray
    stderrs: np.ndarray
    n: int
    parameters: int
    sse: float
    aic: float
    residual_scatter: float   # multiplicative, 10**sigma
    condition_number: float
    per_experiment: bool
    km: float = np.nan
    km_interval: tuple = (np.nan, np.nan)
    absorbed: tuple = ()      # axes the per-experiment offsets swallow

    def coefficient(self, name):
        if name not in self.names:
            return np.nan, np.nan
        index = self.names.index(name)
        return float(self.coefficients[index]), float(self.stderrs[index])


def window_quanta(times, values, fraction=INITIAL_WINDOW):
    """
    How far the fitted window rises, in instrument quanta.

    The measure of whether a rate is measurable at all. Absorbance is reported
    to three decimals, so a window that climbs less than a few 0.001 AU steps
    cannot constrain a slope however clean the cuvette was -- experiment 25's
    dead sample 2 rises 0.8 quanta, and so does the legitimate bottom rung of a
    titration. That is the point: this says the RATE carries no information, not
    that the cuvette failed, and the two need different consequences. See
    `curve_screen.eligibility`.
    """
    count = window_size(len(times), fraction)
    _, slope, _, _ = line_fit(times[:count], values[:count])
    span = times[count - 1] - times[0]
    if not np.isfinite(slope) or span <= 0:
        return np.nan
    return float(abs(slope) * span / ABSORBANCE_QUANTUM)


def slope_ratio(times, values, fraction=INITIAL_WINDOW):
    """
    Late slope over early slope: the curve's shape as one dimensionless number.

    Below 1 the curve decelerates (a burst), above 1 it accelerates (a lag).
    Being a ratio of two rates from the same cuvette, it divides out the
    per-experiment amplitude offset that dominates v0 -- which is exactly why
    the shape is reported this way rather than as a fitted time constant. The
    identifiable alternative would be tau from an exponential, and tau is not
    identifiable here.
    """
    count = window_size(len(times), fraction)
    early, *_ = line_slope(times[:count], values[:count])
    late, *_ = line_slope(times[-count:], values[-count:])
    if not np.isfinite(early) or early == 0:
        return np.nan
    return float(late / early)


# --- the burst / lag form -------------------------------------------------
#
#     A(t) = c + v_ss.t - B (1 - exp(-t/tau))
#
#     dA/dt = v_ss - (B/tau) exp(-t/tau)
#     v0    = v_ss - B/tau          rate at t = 0
#     v_ss  = dA/dt as t -> inf     rate the curve settles to
#     lag   = B / v_ss              where the asymptote crosses A = c
#
#     B > 0  ->  v0 < v_ss  ->  the rate rises: a LAG
#     B < 0  ->  v0 > v_ss  ->  the rate falls: a BURST
#     B = 0  ->  a straight line
#
# For fixed tau the model is LINEAR in (c, v_ss, B), so fitting is a
# one-dimensional profile over tau wrapped around a linear solve: no starting
# values, no local minima in three of the four parameters.
#
# It is degenerate at both ends of the tau range, which is why this is a
# diagnostic here and not the source of any rate:
#
#   tau -> inf   1 - exp(-t/tau) ~ t/tau - t^2/(2 tau^2), whose leading term is
#                exactly collinear with v_ss.t. Only the t^2 term separates
#                them and it is suppressed by 1/tau^2, so v_ss and B/tau slide
#                along a flat valley once tau approaches the run length.
#   tau -> 0     1 - exp(-t/tau) -> 1, a step: B merges into c while
#                v0 = v_ss - B/tau diverges.
#
# So tau is profiled, not just optimised, and a fit whose interval runs to
# either end of the grid is reported UNRESOLVED rather than quoted.
BURST_TAU_FLOOR = 1 / 300.0    # grid start, as a fraction of the run length
BURST_TAU_CAP = 2.0            # grid end, as a multiple of the run length
BURST_GRID_POINTS = 240
BURST_MINIMUM_POINTS = 6       # four parameters need more than four points


@dataclass(frozen=True)
class BurstFit:
    """One curve described as a transient decaying onto a steady rate."""
    c: float
    v_ss: float
    B: float
    tau: float
    v0: float
    lag_time: float
    sse: float
    rms: float
    tau_interval: tuple
    resolved: bool
    points: int

    @property
    def kind(self):
        if not np.isfinite(self.B) or not self.resolved:
            return "unresolved"
        if self.B > 0:
            return "lag"
        if self.B < 0:
            return "burst"
        return "linear"

    def predict(self, times):
        times = np.asarray(times, dtype=float)
        if not np.isfinite(self.tau) or self.tau <= 0:
            return np.full(len(times), np.nan)
        return self.c + self.v_ss * times - self.B * (1 - np.exp(-times / self.tau))


def _burst_design(times, tau):
    return np.column_stack([np.ones(len(times)), times,
                            -(1 - np.exp(-times / tau))])


def fit_burst(times, values, cap=BURST_TAU_CAP, floor=BURST_TAU_FLOOR,
              points=BURST_GRID_POINTS):
    """
    Fit A = c + v_ss.t - B(1 - exp(-t/tau)), profiling tau on a log grid.

    Returns a BurstFit. `resolved` is False when the 95% profile interval on
    tau reaches either end of the grid, which means the data do not locate the
    transient and none of v0, B or tau from that fit should be quoted.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    times = times - times[0]
    span = times[-1]
    blank = BurstFit(*( [np.nan] * 8 + [(np.nan, np.nan), False, len(times)] ))
    if span <= 0 or len(times) < BURST_MINIMUM_POINTS:
        return blank

    grid = np.logspace(np.log10(span * floor), np.log10(span * cap), points)
    costs = np.empty(len(grid))
    solutions = []
    for index, tau in enumerate(grid):
        design = _burst_design(times, tau)
        beta, *_ = np.linalg.lstsq(design, values, rcond=None)
        residual = values - design @ beta
        costs[index] = residual @ residual
        solutions.append(beta)

    best = int(np.argmin(costs))
    tau = float(grid[best])
    c, v_ss, B = (float(v) for v in solutions[best])
    degrees = max(1, len(times) - 4)
    inside = grid[costs <= costs[best] * (1 + CHI2_95 / degrees)]
    low, high = (float(inside.min()), float(inside.max())) if len(inside) else (np.nan, np.nan)
    resolved = bool(np.isfinite(low) and low > grid[0] * 1.05 and high < grid[-1] * 0.95)
    return BurstFit(
        c=c, v_ss=v_ss, B=B, tau=tau,
        v0=v_ss - B / tau,
        lag_time=B / v_ss if v_ss != 0 else np.nan,
        sse=float(costs[best]),
        rms=float(np.sqrt(costs[best] / len(times))),
        tau_interval=(low, high),
        resolved=resolved,
        points=len(times),
    )


# --- the two-phase form, and choosing between it and the one above ---------
#
# WHY A SECOND TERM. `fit_burst` has ONE exponential, so its rate
# v_ss - (B/tau) exp(-t/tau) is MONOTONE: it approaches v_ss from below (lag,
# B > 0) or from above (burst, B < 0) and never turns. Fourteen of the
# temperature series' 24 curves have a rate that rises to a maximum and then
# falls -- exp 16's upper rungs peak at 1666-2303 s and end about 30% down --
# and no value of B represents that. Forced through the one-phase form those
# curves sit at up to 3.8x their noise.
#
# WHAT THE SECOND TERM IS NOT. Its time constant must NOT be read as a decay
# rate constant or given an activation energy. The process behind the fall is
# unidentified: substrate depletion is excluded (0.5-1.1% conversion), catalyst
# inactivation is excluded by the Selwyn test on exps 59/60
# (`scope.selwyn_test` -- inactivation needs the ratio below 1 and it is 1.10
# to 2.20), and a single product threshold is unsupported. See MECHANISM.md.
# The term is DESCRIPTIVE: it keeps the shape from corrupting v_ss, and that is
# all it is for.
TWO_PHASE_GRID_POINTS = 48     # per axis; the pair is profiled on tau1 < tau2
# The second phase costs two parameters, so it has to earn them. This is the
# F threshold on (RSS1 - RSS2)/2 over RSS2/(n - 6). F(2, ~125) at alpha = 0.001
# is about 7.3; 12.0 is deliberately well above it because progress-curve
# residuals are serially correlated -- consecutive readings share whatever the
# model is missing -- and a nominal F over-rejects the simpler model when they
# are. Treat it as "clearly better", not as a p-value.
TWO_PHASE_F = 12.0


@dataclass
class TwoPhaseFit:
    """
    One curve as a fast transient onto a steady rate, then a slow drift off it.

        A = c + v_ss.t - B1(1 - exp(-t/tau1)) - B2(1 - exp(-t/tau2))

    Sign convention follows BurstFit: B > 0 is a LAG (the rate starts low and
    rises), B < 0 a burst. The shape this exists for is B1 > 0 > B2 with
    tau1 < tau2 -- switch on, then fall away -- and `kind` names what was
    actually found rather than assuming it.

    `resolved` says whether TAU2 IS QUOTABLE -- whether its 95% profile
    interval stays inside the grid. It does NOT say whether the second phase is
    real, and the two must not be confused: at 40 C and 5.549 mM the slow phase
    is unmistakable (F = 791, residual 3.81x noise down to 1.05x) while tau2 is
    entirely unlocated, because the run is 6517 s and any tau2 from 2545 s to
    the grid cap fits equally well. Using resolution to decide model FORM is
    the mistake `model_residual` was written to prevent -- `bounded` asks
    whether a parameter is pinned, the residual asks whether the form fits, and
    exp 65 is where they came apart.

    `v_ss_stderr` is the companion warning. As tau2 grows past the run the
    second term flattens into a straight line and becomes degenerate with
    v_ss.t, so v_ss stops being an asymptotic rate and starts absorbing the
    drift. A large stderr here says exactly that.
    """
    c: float
    v_ss: float
    B1: float
    tau1: float
    B2: float
    tau2: float
    sse: float
    rms: float
    v_ss_stderr: float
    tau2_interval: tuple
    resolved: bool
    points: int

    @property
    def kind(self):
        if not (np.isfinite(self.B1) and np.isfinite(self.B2)):
            return "unresolved"
        if self.B1 > 0 and self.B2 < 0:
            return "lag then fall"
        if self.B1 < 0 and self.B2 < 0:
            return "burst then fall"
        if self.B1 > 0 and self.B2 > 0:
            return "two lags"
        return "mixed"

    @property
    def peak_rate(self):
        """The maximum of the fitted rate, and when it occurs: (rate, time)."""
        if not np.isfinite(self.tau1) or not np.isfinite(self.tau2):
            return np.nan, np.nan
        grid = np.linspace(0.0, self.tau2 * 6.0, 2000)
        rate = (self.v_ss - (self.B1 / self.tau1) * np.exp(-grid / self.tau1)
                - (self.B2 / self.tau2) * np.exp(-grid / self.tau2))
        best = int(np.argmax(rate))
        return float(rate[best]), float(grid[best])

    def predict(self, times):
        times = np.asarray(times, dtype=float)
        if not (np.isfinite(self.tau1) and np.isfinite(self.tau2)):
            return np.full(len(times), np.nan)
        return (self.c + self.v_ss * times
                - self.B1 * (1 - np.exp(-times / self.tau1))
                - self.B2 * (1 - np.exp(-times / self.tau2)))


def _two_phase_design(times, tau1, tau2):
    return np.column_stack([np.ones(len(times)), times,
                            -(1 - np.exp(-times / tau1)),
                            -(1 - np.exp(-times / tau2))])


def fit_two_phase(times, values, cap=BURST_TAU_CAP, floor=BURST_TAU_FLOOR,
                  points=TWO_PHASE_GRID_POINTS):
    """
    Fit the two-phase form, profiling (tau1, tau2) on a log grid, tau1 < tau2.

    Given both time constants the model is LINEAR in c, v_ss, B1 and B2, so
    every node is one least-squares solve and the search has no local minima to
    fall into -- the same reason `fit_burst` profiles tau rather than handing
    all four parameters to an optimiser. The grid is coarser per axis than the
    one-phase one (48 against 240) because it is two-dimensional; that costs
    resolution on the time constants and none on v_ss, which is what the form
    exists to protect.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    times = times - times[0]
    span = times[-1]
    blank = TwoPhaseFit(*([np.nan] * 9 + [(np.nan, np.nan), False, len(times)]))
    # Six parameters need more than six points, and the profile needs slack.
    if span <= 0 or len(times) < BURST_MINIMUM_POINTS + 2:
        return blank

    grid = np.logspace(np.log10(span * floor), np.log10(span * cap), points)
    # The design's columns are [1, t, -E_i, -E_j] and only the last two depend
    # on the grid, so every inner product the normal equations need can be
    # computed ONCE and indexed. That turns each node from an n-row least
    # squares into a 4x4 solve: about 400x faster, which is what makes this
    # affordable on all 402 curves rather than on one block.
    basis = -(1.0 - np.exp(-np.outer(1.0 / grid, times)))     # (points, n)
    ones = np.ones(len(times))
    fixed = np.column_stack([ones, times])                    # (n, 2)
    ff = fixed.T @ fixed                                      # (2, 2)
    fb = fixed.T @ basis.T                                    # (2, points)
    bb = basis @ basis.T                                      # (points, points)
    fy = fixed.T @ values                                     # (2,)
    by = basis @ values                                       # (points,)
    total = float(values @ values)

    best_cost, best = np.inf, None
    profile = np.full(len(grid), np.inf)
    for second in range(1, len(grid)):
        # All `first < second` at once: a stack of 4x4 systems, solved in one
        # call. The loop over `second` stays because the stack is what makes
        # each step cheap, not the Python.
        first = np.arange(second)
        normal = np.empty((second, 4, 4))
        normal[:, :2, :2] = ff
        normal[:, :2, 2] = fb[:, first].T
        normal[:, :2, 3] = fb[:, second]
        normal[:, 2, :2] = fb[:, first].T
        normal[:, 3, :2] = fb[:, second]
        normal[:, 2, 2] = bb[first, first]
        normal[:, 3, 3] = bb[second, second]
        normal[:, 2, 3] = normal[:, 3, 2] = bb[first, second]
        target = np.empty((second, 4))
        target[:, :2] = fy
        target[:, 2] = by[first]
        target[:, 3] = by[second]
        with np.errstate(all="ignore"):
            try:
                # target[..., None] because numpy reads a (M, 4) right-hand
                # side against a (M, 4, 4) stack as a single matrix, not as M
                # vectors.
                beta = np.linalg.solve(normal, target[..., None])[..., 0]
            except np.linalg.LinAlgError:
                # One singular node in the stack fails the whole batch, so fall
                # back to the least-squares route for this `second` only.
                beta = np.stack([np.linalg.lstsq(n, t, rcond=None)[0]
                                 for n, t in zip(normal, target)])
        # SSE = y.y - 2 b.X'y + b'X'X b, from the same precomputed pieces.
        cost = (total - 2 * np.einsum("ij,ij->i", beta, target)
                + np.einsum("ij,ijk,ik->i", beta, normal, beta))
        cost = np.where(np.isfinite(cost), cost, np.inf)
        pick = int(np.argmin(cost))
        if cost[pick] < profile[second]:
            profile[second] = float(cost[pick])
        if cost[pick] < best_cost:
            best_cost = float(cost[pick])
            best = (int(first[pick]), second, beta[pick].copy())
    if best is None or not np.isfinite(best_cost):
        return blank
    # Floating point can drive a residual sum of squares slightly negative when
    # the fit is near-perfect; it is a cost, so clamp rather than propagate.
    best_cost = max(best_cost, 0.0)

    first, second, beta = best
    degrees = max(1, len(times) - 6)
    inside = grid[1:][profile[1:] <= best_cost * (1 + CHI2_95 / degrees)]
    low, high = ((float(inside.min()), float(inside.max())) if len(inside)
                 else (np.nan, np.nan))
    resolved = bool(np.isfinite(low) and low > grid[1] * 1.05
                    and high < grid[-1] * 0.95)
    c, v_ss, B1, B2 = (float(v) for v in beta)
    # v_ss's own standard error, from the linear solve at the chosen pair. It
    # is conditional on (tau1, tau2) and so understates the true uncertainty,
    # but it catches the failure that matters: a tau2 long enough that its term
    # is a straight line, which makes v_ss and B2/tau2 trade off freely.
    design = _two_phase_design(times, grid[first], grid[second])
    variance = best_cost / max(1, len(times) - 4)
    covariance = variance * np.linalg.pinv(design.T @ design)
    return TwoPhaseFit(
        c=c, v_ss=v_ss, B1=B1, tau1=float(grid[first]),
        B2=B2, tau2=float(grid[second]),
        sse=best_cost, rms=float(np.sqrt(best_cost / len(times))),
        v_ss_stderr=float(np.sqrt(max(covariance[1, 1], 0.0))),
        tau2_interval=(low, high), resolved=resolved, points=len(times))


@dataclass
class ProgressFit:
    """Whichever form the curve earned, and the evidence for the choice."""
    phases: int
    one: BurstFit
    two: TwoPhaseFit
    f_statistic: float
    reason: str

    @property
    def chosen(self):
        return self.two if self.phases == 2 else self.one

    @property
    def v_ss(self):
        """
        The t -> infinity asymptote of whichever form was chosen.

        DO NOT PUT THIS ON AN ARRHENIUS PLOT when two phases were selected. It
        is an extrapolation far outside the data, and with a decay in the fit
        it is not constrained by anything: on the temperature series it comes
        out NEGATIVE on two of the 35 C curves and gives an Arrhenius scatter
        of 0.19 to 1.18 in ln units, against 0.08 for the peak rate. Use
        `peak_rate`.
        """
        return self.chosen.v_ss

    @property
    def kind(self):
        """
        WHICH WAY THE CURVE POINTS, from whichever form it earned.

        The one-phase form has two shapes (`lag`, `burst`) and the two-phase
        form four (`lag then fall`, `burst then fall`, `two lags`, `mixed`),
        and the sign convention is the same in both: B > 0 means the rate
        starts BELOW its eventual value and rises.

        This is computed on every curve in the archive and was, until
        2026-09-02, thrown away -- `scope.frame` carried the time constants and
        the peak rate and no sign at all, so nothing could ask whether a curve
        began slow or began fast. It now carries `progress_kind`, `B_fast` and
        `B_slow`, and `induction.sign_table` reads them.
        """
        if self.phases == 2:
            return self.two.kind
        one = self.one
        if not np.isfinite(one.B):
            return "unresolved"
        return "lag" if one.B > 0 else "burst"

    @property
    def amplitudes(self):
        """The chosen form's (B_fast, B_slow); B_slow is nan on one phase."""
        if self.phases == 2:
            return float(self.two.B1), float(self.two.B2)
        return float(self.one.B), float("nan")

    @property
    def peak_rate(self):
        """
        The largest rate the fitted model reaches, and when: (rate, time).

        THE OBSERVABLE TO USE, and it is defined the same way on both forms:
        the maximum of the fitted dA/dt. What that means differs by shape and
        the difference is the point --

          one phase, lag    the rate rises monotonically to v_ss, so the
                            supremum IS v_ss, reached as t -> infinity
          one phase, burst  the rate falls from v0, so the maximum is v0 at 0
          two phases        an interior maximum, which is what the second term
                            exists to locate

        This is what `curve_metrics.peak_rate` measures off the raw readings,
        with the truncation removed: the block statistic can only find a
        maximum inside the window, and at 15 and 20 C the window ends before
        the rate levels off. The fit knows the shape and does not have to stop
        where the readings do.
        """
        if self.phases == 2:
            return self.two.peak_rate
        one = self.one
        if not np.isfinite(one.tau) or not np.isfinite(one.B):
            return np.nan, np.nan
        if one.B > 0:                      # lag: rises to v_ss
            return float(one.v_ss), float("inf")
        return float(one.v_ss - one.B / one.tau), 0.0   # burst: highest at t=0

    def predict(self, times):
        return self.chosen.predict(times)


def fit_progress(times, values, threshold=TWO_PHASE_F, **kwargs):
    """
    Fit both forms and return the one the curve earns, with the reason.

    NESTED SELECTION, and the nesting is exact: B2 = 0 turns the two-phase form
    into the one-phase form, so the comparison is a plain F test on two extra
    parameters, and F > `threshold` is the whole rule.

    IT IS DELIBERATELY NOT GATED ON TAU2 BEING RESOLVED. The first version of
    this function required both, and it rejected the clearest two-phase curves
    in the archive -- 40 C at 5.549 mM, F = 791, residual falling from 3.81x
    noise to 1.05x -- because tau2 is 7562 s in a 6517 s run and anything above
    2545 s fits as well. Whether a phase EXISTS and whether its time constant
    is PINNED are different questions, and answering the first with the second
    is the error `model_residual` exists to prevent. `TwoPhaseFit.resolved`
    still reports the second question, and tau2 should not be quoted where it
    is False.

    Returns a ProgressFit. Read `phases`, then `chosen`.
    """
    one = fit_burst(times, values, **kwargs)
    two = fit_two_phase(times, values, **kwargs)
    count = len(np.asarray(times))
    if not np.isfinite(one.sse) or not np.isfinite(two.sse):
        return ProgressFit(1, one, two, np.nan, "a form did not fit")
    degrees = max(1, count - 6)
    gain = one.sse - two.sse
    statistic = (float((gain / 2.0) / (two.sse / degrees))
                 if two.sse > 0 else np.nan)
    if not np.isfinite(statistic) or statistic <= threshold:
        return ProgressFit(1, one, two, statistic,
                           "the second phase does not pay for its parameters")
    return ProgressFit(2, one, two, statistic, "the second phase is earned")


# How wide the v0 profile interval may be, as a fraction of v0, before the
# initial rate counts as unbounded. 0.30 is where the enzyme-free BnOH curves
# separate: the curves inside it agree with their own line fit to 1.02-2.14x,
# the ones outside disagree by up to 28x.
BURST_V0_HALFWIDTH = 0.30


@dataclass(frozen=True)
class BoundedBurstFit:
    """
    The burst/lag form fitted with B <= 0, and v0 profiled rather than quoted.

    TWO CHANGES FROM `fit_burst`, both aimed at the initial rate.

    1. `B <= 0`, BUT ONLY WHERE THE CURVE DOES NOT ACCELERATE. The
       unconstrained fit is free to choose a LAG -- a rate that rises -- and on
       a near-straight curve it sometimes does, with a collapsed tau, returning
       a NEGATIVE v0 where the line gives a firmly positive one (4 of the 27
       enzyme-free BnOH curves; exp 67 sample 3 gives -2.07e-4 against a line's
       +3.35e-6, with a profile interval that never reaches zero).

       A blanket `B <= 0` was the first fix and was wrong. It rested on "0 of
       16 pass the `acceleration` test", which is true of the constant-buffer
       runs and was generalised to all 27 without checking: exps 3 and 6 hold
       four curves that DO accelerate, two of them at z = +8.4 and +11.8, and
       the blanket rule bound on two of them -- forcing a decelerating shape
       onto curves the data says are rising.

       `constrain="auto"` therefore asks each curve. Where `acceleration`
       clears ACCELERATION_SIGMA the lag branch stays open; everywhere else it
       is shut. On the enzyme-free BnOH set that bounds 19 of 27 against 13
       unconstrained, with no negative v0 and no curve made to fit a shape its
       own z-score contradicts. `constrain=True`/`False` force the old
       behaviours.

       `noise_floor` is the SOURCE's floor and reaches `acceleration`, which
       divides by two standard errors that `line_fit` floors. Passing the .txt
       default on .rre data suppresses the z this decision is made on -- see
       `fit_dataset.source_floor`.

    2. `v0_low`/`v0_high` come from profiling v0 across every tau whose cost is
       within the 95% band, NOT from a standard error at the optimum. This is
       the diagnostic that matters and `fit_burst.resolved` is not a substitute
       for it: `resolved` asks whether TAU is located, and those are different
       questions. Exps 6's four cuvettes have tau unresolved yet v0 pinned to a
       0.00 half-width -- as tau -> inf the curve is a straight line, which
       kills tau and B but leaves v0 -> v_ss perfectly determined.

    WHAT THE CONSTRAINT DOES NOT FIX. It raises the number of bounded curves
    from 13 to 21 of 27 and removes every negative v0, but it trades the
    tau -> inf degeneracy for the tau -> 0 one: exp 65 samples 1 and 2 come
    back at 10x and 28x their line rate, with intervals to match. v0 is an
    extrapolation to the boundary, so its uncertainty is set by the assumed
    shape near t = 0, and no fitting refinement recovers information the curve
    does not carry. Read `bounded` before quoting `v0`.
    """
    c: float
    v_ss: float
    B: float
    tau: float
    v0: float
    v0_low: float
    v0_high: float
    tau_low: float
    tau_high: float
    tau_resolved: bool
    sse: float
    rms: float
    points: int
    bounded: bool

    @property
    def half_width(self):
        """v0's profile half-width as a fraction of v0."""
        if not np.isfinite(self.v0) or self.v0 == 0:
            return np.inf
        return float((self.v0_high - self.v0_low) / (2 * abs(self.v0)))

    @property
    def shape_is_meaningful(self):
        """
        True when tau is located, so `kind`, `B` and `v_ss` mean something.

        `bounded` ASKS ONLY ABOUT v0 AND IS NOT A SUBSTITUTE. On the
        enzyme-free BnOH set tau is pinned at a grid end on 15 of 27 curves and
        11 of those still report `bounded` -- because once tau collapses to the
        floor or runs to the cap the model degenerates (to a step, to a
        straight line) and v0 -> v_ss becomes exactly determined. v0 is then
        trustworthy and the BURST is not: exps 69 sample 3 and 70 sample 4 are
        bounded, with tau at the cap and a NEGATIVE v_ss.

        This restores what `fit_burst.resolved` reports and what this class
        dropped when it was written. See DATA_VERIFICATION.md, 2026-09-01.
        """
        return bool(self.tau_resolved)

    @property
    def kind(self):
        """
        What shape was fitted: "burst", "lag", "clamped" or "unresolved".

        WITHOUT THIS, v0 IS AMBIGUOUS. On a burst (B < 0) the rate falls, so v0
        is the MAXIMUM and v_ss the settled rate. On a lag (B > 0) the rate
        rises, so v0 is the INDUCTION rate -- the reaction before it gets going
        -- and v_ss is the developed one. Reporting v0 from both as "the
        initial rate" puts two different quantities in one column, which is the
        same trap `curve_metrics.peak_rate` warns about for vmax against v0.

        "clamped" is B == 0 reached because the lag branch was SHUT, not
        because the curve is straight: the constrained optimum is a line
        whenever the free optimum wanted B > 0. It is reported separately so a
        clamped fit is never read as evidence of linearity.
        """
        if not np.isfinite(self.B):
            return "unresolved"
        if self.B > 0:
            return "lag"
        if self.B < 0:
            return "burst"
        return "clamped"

    @property
    def settles_backwards(self):
        """
        True when v_ss is negative: a late rate that runs the reaction backwards.

        Physically impossible here, and it happens -- exps 69 sample 3 and 70
        sample 4 -- on fits whose v0 profile is otherwise tight. tau has
        collapsed far enough that the exponential absorbs the whole curve and
        v_ss is extrapolating past it. `bounded` does not catch this because it
        asks only about v0, so check both before quoting either.
        """
        return bool(np.isfinite(self.v_ss) and self.v_ss < 0)

    def predict(self, times):
        times = np.asarray(times, dtype=float)
        if not np.isfinite(self.tau) or self.tau <= 0:
            return np.full(len(times), np.nan)
        return self.c + self.v_ss * times - self.B * (1 - np.exp(-times / self.tau))


def _burst_solve(times, values, tau, constrain):
    """(sse, beta) for one tau, clamping B to 0 if a positive B is forbidden."""
    design = _burst_design(times, tau)
    beta, *_ = np.linalg.lstsq(design, values, rcond=None)
    if constrain and beta[2] > 0:
        # B > 0 is a lag. Forbidden: refit as a straight line, which is the
        # constrained optimum whenever the unconstrained one wants B > 0.
        line, *_ = np.linalg.lstsq(design[:, :2], values, rcond=None)
        beta = np.array([line[0], line[1], 0.0])
    residual = values - design @ beta
    return float(residual @ residual), beta


def fit_burst_bounded(times, values, cap=BURST_TAU_CAP, floor=BURST_TAU_FLOOR,
                      points=BURST_GRID_POINTS, constrain="auto",
                      half_width=BURST_V0_HALFWIDTH,
                      noise_floor=QUANTISATION_SIGMA):
    """
    Fit the burst/lag form and profile v0. Returns a BoundedBurstFit.

    `constrain` is "auto" (shut the lag branch only where this curve does not
    accelerate), True (always shut it) or False (never). See BoundedBurstFit.

    `bounded` is True only when the whole 95% profile interval on v0 is
    positive AND narrower than `half_width` either side. Quote `v0` only then.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    times = times - times[0]
    span = times[-1] if len(times) else 0.0
    blank = BoundedBurstFit(*([np.nan] * 9 + [False, np.nan, np.nan,
                                              len(times), False]))
    if span <= 0 or len(times) < BURST_MINIMUM_POINTS:
        return blank

    if constrain == "auto":
        # The curve's own verdict, not a rule about the block it sits in.
        z, _ = acceleration(times, values, floor=noise_floor)
        forbid_lag = not (np.isfinite(z) and z > ACCELERATION_SIGMA)
    else:
        forbid_lag = bool(constrain)

    grid = np.logspace(np.log10(span * floor), np.log10(span * cap), points)
    costs = np.empty(len(grid))
    rates = np.empty(len(grid))
    solutions = []
    for index, tau in enumerate(grid):
        cost, beta = _burst_solve(times, values, tau, forbid_lag)
        costs[index] = cost
        rates[index] = beta[1] - beta[2] / tau
        solutions.append(beta)

    best = int(np.argmin(costs))
    degrees = max(1, len(times) - 4)
    inside = costs <= costs[best] * (1 + CHI2_95 / degrees)
    low, high = float(rates[inside].min()), float(rates[inside].max())
    # tau's own profile interval, and whether it clears both ends of the grid.
    # Running to either end means the transient is not located and the model
    # has degenerated -- to a straight line at the cap, to a step at the floor.
    # The same test `fit_burst.resolved` applies, restored here.
    taus = grid[inside]
    tau_low, tau_high = float(taus.min()), float(taus.max())
    resolved = bool(tau_low > grid[0] * 1.05 and tau_high < grid[-1] * 0.95)
    c, v_ss, B = (float(v) for v in solutions[best])
    v0 = float(rates[best])
    relative = (high - low) / (2 * abs(v0)) if v0 != 0 else np.inf
    return BoundedBurstFit(
        c=c, v_ss=v_ss, B=B, tau=float(grid[best]), v0=v0,
        v0_low=low, v0_high=high,
        tau_low=tau_low, tau_high=tau_high, tau_resolved=resolved,
        sse=float(costs[best]),
        rms=float(np.sqrt(costs[best] / len(times))),
        points=len(times),
        bounded=bool(low > 0 and relative < half_width),
    )


def buffer_concentrations(dataset_path=DATASET_PATH):
    """{(experiment, sample): [buf] in mM} -- the one axis `Curve` does not carry."""
    data = add_solution_columns(pd.read_csv(dataset_path))
    return {(int(row["experiment"]), int(row["sample"])): float(row["[buf]"])
            for row in data.to_dict("records")}


def summarise(curves, buffers=None, fraction=INITIAL_WINDOW):
    """Reduce every curve to a CurveSummary."""
    if buffers is None:
        buffers = buffer_concentrations()
    summaries = []
    for curve in curves:
        v0, stderr, rms = initial_rate(curve.times, curve.absorbance, fraction)
        amplitude = float(curve.absorbance.max() - curve.absorbance.min())
        denominator = curve.epsilon * curve.conditions.s0
        summaries.append(CurveSummary(
            experiment=curve.experiment,
            sample=curve.sample,
            substrate=curve.substrate,
            temperature=curve.temperature,
            buffer_name=curve.buffer,
            pH=curve.pH,
            s0=curve.conditions.s0,
            h2o2=curve.conditions.h2o2,
            hoo=curve.conditions.hoo,
            buf=buffers.get((curve.experiment, curve.sample), np.nan),
            e0=curve.conditions.e0,
            v0=v0,
            v0_stderr=stderr,
            window_rms=rms,
            window_quanta=window_quanta(curve.times, curve.absorbance, fraction),
            noise=curve.noise,
            slope_ratio=slope_ratio(curve.times, curve.absorbance, fraction),
            points=len(curve),
            duration=float(curve.times[-1]),
            amplitude=amplitude,
            conversion=float(amplitude / denominator) if denominator > 0 else np.nan,
        ))
    return summaries


def to_frame(summaries):
    """CurveSummary list -> DataFrame, with `group` as a printable string."""
    frame = pd.DataFrame([asdict(s) for s in summaries])
    if len(frame):
        frame["group"] = [str(s.group) for s in summaries]
    return frame


def replicate_scatter(frame, minimum=3):
    """
    Relative scatter of v0 across cuvettes run at IDENTICAL conditions.

    This is the noise floor a regression should be judged against: no model,
    right or wrong, can explain variation the experiment does not reproduce.
    Grouped inside an experiment, so a day-to-day offset cannot inflate it.
    """
    keys = ["experiment", "s0", "h2o2", "hoo", "buf", "e0"]
    sets = []
    for key, block in frame.groupby(keys, dropna=False):
        rates = block.v0.to_numpy()
        if len(block) < minimum or not np.all(np.isfinite(rates)) or np.mean(rates) == 0:
            continue
        sets.append(dict(experiment=int(key[0]), n=len(block), s0=float(key[1]),
                         buf=float(key[4]), mean=float(np.mean(rates)),
                         rsd=float(np.std(rates, ddof=1) / abs(np.mean(rates))),
                         noise_ratio=float(np.std(rates, ddof=1) / abs(np.mean(rates))
                                           / max(np.mean(block.v0_stderr / block.v0.abs()),
                                                 1e-12))))
    return sets


# A cuvette this far BELOW its own experiment's median rate, on a curve that
# also barely cleared the noise, is not a slow reaction -- it is a failed one,
# the dead runs already documented in DATA_VERIFICATION.md. Two such curves
# contribute more to the residual sum of squares of the offset model than the
# other fifty-six combined, so they are always reported, and `--drop-outliers`
# excludes them.
#
# Both conditions are required, and neither alone would do. The median ratio on
# its own assumes conditions barely move inside an experiment, which is true of
# the 4OMe/40 block (6% of the substrate variance is within-experiment) and
# would misfire badly on a design with wide ladders. The signal-to-noise floor
# on its own would condemn the legitimately slow low-substrate curves --
# experiment 26's four healthy replicates sit at SNR 56-63, barely above the
# dead curves. Together they are specific.
OUTLIER_FACTOR = 5.0
DEAD_CURVE_SNR = 30.0


def experiment_outliers(frame, factor=OUTLIER_FACTOR, response="v0",
                        snr_floor=DEAD_CURVE_SNR):
    """
    Curves that look like failed cuvettes rather than slow ones.

    A curve is flagged when its response is non-positive, or when it is more
    than `factor` BELOW its own experiment's median AND the whole curve rose
    less than `snr_floor` times its own point noise. Judged inside the
    experiment, so the between-experiment offset that dominates this dataset
    cannot make a healthy cuvette look dead. Only low outliers are flagged: a
    dead cuvette reads low, never high.

    This is a screen, not a rule. It is always reported and never applied
    unless the caller asks (`--drop-outliers`).
    """
    flagged = []
    for _, block in frame.groupby("experiment"):
        rates = block[response].to_numpy(dtype=float)
        good = rates[np.isfinite(rates) & (rates > 0)]
        if len(good) < 3:
            continue
        centre = float(np.median(good))
        for row in block.to_dict("records"):
            value = float(row[response])
            noise = float(row.get("noise", 0.0) or 0.0)
            snr = (float(row.get("amplitude", np.inf)) / noise
                   if noise > 0 else np.inf)
            dead = (not np.isfinite(value) or value <= 0
                    or (value * factor < centre and snr < snr_floor))
            if dead:
                flagged.append((int(row["experiment"]), int(row["sample"]),
                                value, centre))
    return flagged


def _dummies(frame):
    """One indicator column per experiment, and their names."""
    experiments = sorted(frame.experiment.unique())
    columns = [(frame.experiment.to_numpy() == e).astype(float) for e in experiments]
    return columns, [f"exp{e}" for e in experiments]


# Fraction of an axis's log-variance that must survive per-experiment centring
# for its coefficient to mean anything. An axis that barely moves inside an
# experiment is nearly collinear with that experiment's indicator: lstsq still
# returns a number, and the number is a ratio of two small quantities. The
# symptom is a coefficient far larger than any chemistry allows carrying a
# standard error larger than itself.
ABSORPTION_FLOOR = 0.10


def within_experiment_variation(frame, axes=DESIGN_AXES):
    """
    Fraction of each axis's log-variance surviving per-experiment centring.

    Per-experiment offsets can only be paid for out of WITHIN-experiment
    contrast. This is how much of each axis is left to pay with, and it is the
    number that decides whether an order from the offset model is real.
    """
    surviving = {}
    for axis in axes:
        values = np.log10(frame[axis].to_numpy(dtype=float))
        if not np.all(np.isfinite(values)) or np.var(values) <= 0:
            continue
        means = frame.groupby("experiment")[axis].transform(
            lambda v: np.log10(v.astype(float)).mean()).to_numpy()
        surviving[axis] = float(np.var(values - means) / np.var(values))
    return surviving


def absorbed_axes(frame, axes=DESIGN_AXES, floor=ABSORPTION_FLOOR):
    """
    Axes the per-experiment offsets swallow, wholly or nearly.

    An axis held constant WITHIN every experiment is perfectly collinear with
    that experiment's indicator and its coefficient is not identifiable at all;
    one that varies only slightly is identifiable in principle and worthless in
    practice. Both are reported rather than silently dropped.
    """
    return tuple(axis for axis, fraction in within_experiment_variation(frame, axes).items()
                 if fraction < floor)


def design_matrix(frame, km=None, per_experiment=True, axes=DESIGN_AXES):
    """
    Columns of the log-log regression.

    The substrate column is log10([S]) for a power law, or log10([S]/(Km+[S]))
    for a saturating law -- the only nonlinearity in the whole model, and it is
    one-dimensional, so it is profiled on a grid rather than handed to an
    optimiser.
    """
    columns, names = [], []
    substrate = frame.s0.to_numpy()
    # A column with no variance is exactly collinear with the intercept, or
    # with the per-experiment offsets. lstsq via pinv will still split the
    # coefficient between them and report a small standard error on a number
    # that means nothing -- a constant-[S] block reported an order of -3.906
    # +/- 0.080 before this check existed. Drop it instead.
    if np.std(np.log10(substrate)) > 1e-9:
        if km is None:
            columns.append(np.log10(substrate))
            names.append("log[S]")
        else:
            columns.append(np.log10(substrate / (km + substrate)))
            names.append("saturation([S])")
    for axis in axes:
        if axis == "s0":
            continue
        values = np.log10(frame[axis].to_numpy())
        if np.std(values) > 1e-9:
            columns.append(values)
            names.append(f"log[{axis}]")
    if per_experiment:
        extra, extra_names = _dummies(frame)
        columns.extend(extra)
        names.extend(extra_names)
    else:
        columns.append(np.ones(len(frame)))
        names.append("const")
    return np.column_stack(columns), names


def _solve(design, response):
    """lstsq plus the covariance a rank-deficient design still permits."""
    coefficients, *_ = np.linalg.lstsq(design, response, rcond=None)
    residual = response - design @ coefficients
    sse = float(residual @ residual)
    degrees = max(1, len(response) - np.linalg.matrix_rank(design))
    covariance = (sse / degrees) * np.linalg.pinv(design.T @ design)
    return coefficients, np.sqrt(np.abs(np.diag(covariance))), sse, residual


def _condition(design, experiments=None):
    """
    Condition number of the part of the design the fit actually identifies.

    Indicator columns are dropped, and when per-experiment offsets are present
    the remaining columns are centred within each experiment first -- the
    offsets absorb every between-experiment contrast, so conditioning measured
    on the uncentred columns flatters a design that has nothing left to fit
    with. This is why an axis can look well conditioned and still come back
    with a standard error larger than its coefficient.
    """
    keep = [i for i in range(design.shape[1])
            if design[:, i].std() > 1e-12 and not set(np.unique(design[:, i])) <= {0.0, 1.0}]
    if len(keep) < 2:
        return 1.0
    block = design[:, keep].copy()
    if experiments is not None:
        for value in np.unique(experiments):
            mask = experiments == value
            block[mask] -= block[mask].mean(axis=0)
    spread = block.std(axis=0)
    if np.any(spread <= 1e-12):
        return np.inf
    return float(np.linalg.cond((block - block.mean(0)) / spread))


def regress(frame, response="v0", km=None, per_experiment=True, axes=DESIGN_AXES):
    """
    Fit log10(response) against the log design. Returns a Regression.

    Rows whose response is not positive are dropped -- a log-log law cannot
    represent them, and on this data a non-positive v0 means a dead cuvette
    rather than a negative rate.
    """
    usable = frame[np.isfinite(frame[response]) & (frame[response] > 0)]
    for axis in axes:
        usable = usable[np.isfinite(usable[axis]) & (usable[axis] > 0)]
    if len(usable) < 4:
        raise ValueError(f"only {len(usable)} usable rows for {response!r}")

    values = np.log10(usable[response].to_numpy())
    design, names = design_matrix(usable, km, per_experiment, axes)
    coefficients, stderrs, sse, residual = _solve(design, values)
    parameters = int(np.linalg.matrix_rank(design)) + (1 if km is not None else 0)
    labels = usable.experiment.to_numpy() if per_experiment else None
    return Regression(
        response=response,
        names=tuple(names),
        coefficients=coefficients,
        stderrs=stderrs,
        n=len(usable),
        parameters=parameters,
        sse=sse,
        aic=len(values) * np.log(sse / len(values)) + 2 * parameters,
        residual_scatter=float(10 ** np.std(residual)),
        condition_number=_condition(design, labels),
        per_experiment=per_experiment,
        absorbed=absorbed_axes(usable, axes) if per_experiment else (),
    )


def profile_km(frame, response="v0", per_experiment=True, axes=DESIGN_AXES,
               low=1e-3, high=3e2, points=700):
    """
    Profile the fit over Km, returning the best Regression with its interval.

    Km is the one parameter the model is nonlinear in, so it gets a grid and a
    profile-likelihood interval rather than a standard error off a Jacobian.
    An interval that runs to both ends of the grid means the data do not
    resolve saturation at all, which is a result and is reported as one.
    """
    if np.std(np.log10(frame.s0.to_numpy())) <= 1e-9:
        # No substrate contrast at all: there is no saturation curve to profile.
        return regress(frame, response, None, per_experiment, axes)
    grid = np.logspace(np.log10(low), np.log10(high), points)
    costs = np.array([regress(frame, response, km, per_experiment, axes).sse
                      for km in grid])
    best_index = int(np.argmin(costs))
    km = float(grid[best_index])
    fitted = regress(frame, response, km, per_experiment, axes)
    degrees = max(1, fitted.n - fitted.parameters)
    inside = grid[costs <= costs[best_index] * (1 + CHI2_95 / degrees)]
    interval = ((float(inside.min()), float(inside.max())) if len(inside)
                else (np.nan, np.nan))
    return Regression(**{**asdict(fitted), "km": km, "km_interval": interval,
                         "coefficients": fitted.coefficients,
                         "stderrs": fitted.stderrs})


def window_sensitivity(curves, buffers=None, response="v0",
                       fractions=(0.10, 0.20, 0.30, 0.50, 1.00),
                       per_experiment=True, drop_outliers=False):
    """
    The substrate order as a function of the window it was measured over.

    The orders this module reports are not window-independent -- they move by a
    factor of three across the range below, which is far larger than their
    statistical error. Any order quoted from here has to carry its window, and
    the spread across windows is the honest systematic on it.
    """
    if buffers is None:
        buffers = buffer_concentrations()
    rows = []
    for fraction in fractions:
        frame = to_frame(summarise(curves, buffers, fraction))
        if drop_outliers:
            flagged = {(e, s) for e, s, _, _ in experiment_outliers(frame)}
            if flagged:
                frame = frame[~frame.set_index(["experiment", "sample"]).index.isin(flagged)]
        try:
            fitted = regress(frame, response, None, per_experiment)
        except ValueError:
            continue
        order, stderr = fitted.coefficient("log[S]")
        quality = frame.window_rms / frame.noise
        rows.append(dict(fraction=fraction, order=order, stderr=stderr,
                         rms_over_noise=float(np.nanmedian(quality)),
                         n=fitted.n, scatter=fitted.residual_scatter))
    return rows


def report_design(frame, title):
    """What the design can and cannot separate, before any fit is quoted."""
    print(f"\n{title}")
    print(f"  {len(frame)} curves, {frame.experiment.nunique()} experiments")
    for axis in DESIGN_AXES + ("h2o2",):
        values = frame[axis].dropna()
        values = values[values > 0]
        if not len(values):
            print(f"    {axis:5s}  --")
            continue
        span = values.max() / values.min()
        flag = "   CONSTANT" if span < 1.01 else ""
        print(f"    {axis:5s}  {values.min():10.4g} - {values.max():<10.4g} "
              f"span {span:8.1f}x  {values.nunique():3d} levels{flag}")
    logs = {a: np.log10(frame[a].to_numpy()) for a in DESIGN_AXES
            if np.all(np.isfinite(np.log10(frame[a].to_numpy().astype(float))))
            and np.std(np.log10(frame[a].to_numpy())) > 1e-9}
    if len(logs) > 1:
        names = list(logs)
        matrix = np.corrcoef(np.array([logs[n] for n in names]))
        print("    log-design correlations:")
        print("         " + "".join(f"{n:>9s}" for n in names))
        for i, name in enumerate(names):
            print(f"      {name:5s}" + "".join(f"{matrix[i, j]:9.2f}"
                                               for j in range(len(names))))
    surviving = within_experiment_variation(frame)
    if surviving and frame.experiment.nunique() > 1:
        print("    log-variance surviving per-experiment centring "
              "(what an offset model has left to fit with):")
        for axis, fraction in surviving.items():
            verdict = "  <-- absorbed by the offsets" if fraction < ABSORPTION_FLOOR else ""
            print(f"      {axis:5s} {fraction:6.1%}{verdict}")


def report_curves(frame, title):
    """Curve-level quality: does the straight line describe its window?"""
    print(f"\n{title}")
    quality = (frame.window_rms / frame.noise).replace([np.inf, -np.inf], np.nan)
    print(f"  straight line over the first {INITIAL_WINDOW:.0%}: "
          f"median rms = {np.nanmedian(quality):.2f}x the curve's point noise")
    relative = (frame.v0_stderr / frame.v0.abs()).replace([np.inf, -np.inf], np.nan)
    print(f"  v0 relative standard error: median {np.nanmedian(relative):.1%} "
          f"(reproducibility floor is {REPLICATE_RSD:.1%})")
    print(f"  conversion: median {frame.conversion.median():.2%}, "
          f"max {frame.conversion.max():.2%}")
    flagged = experiment_outliers(frame)
    if flagged:
        print(f"  {len(flagged)} curve(s) look like failed cuvettes -- non-positive, or "
              f"{OUTLIER_FACTOR:.0f}x below their experiment median with whole-curve "
              f"SNR < {DEAD_CURVE_SNR:.0f} (--drop-outliers excludes them):")
        for experiment, sample, value, centre in flagged:
            print(f"    exp {experiment} sample {sample}  v0={value:.3e}  "
                  f"vs experiment median {centre:.3e}")
    shape = frame.slope_ratio.replace([np.inf, -np.inf], np.nan).dropna()
    if len(shape):
        burst = int((shape < 0.9).sum())
        lag = int((shape > 1.1).sum())
        print(f"  shape (late slope / early slope): median {shape.median():.2f}; "
              f"{burst} decelerating, {lag} accelerating, "
              f"{len(shape) - burst - lag} linear within 10%")


def report_replicates(frame):
    """The reproducibility floor, measured rather than assumed."""
    sets = replicate_scatter(frame)
    print(f"\nREPLICATE SETS (identical conditions inside one experiment)")
    if not sets:
        print("  none -- the reproducibility floor cannot be measured on this subset")
        return
    for entry in sets:
        print(f"  exp {entry['experiment']:3d}  n={entry['n']}  "
              f"[S]={entry['s0']:.2f} mM  [buf]={entry['buf']:.0f} mM  "
              f"v0 = {entry['mean']:.3e} +/- {100 * entry['rsd']:.1f}%")
    print(f"  -> reproducibility floor {100 * np.mean([e['rsd'] for e in sets]):.1f}% "
          f"across {len(sets)} set(s); module default is {100 * REPLICATE_RSD:.1f}%")


def report_regression(fitted, title):
    """One fit, with the axes that matter first and the offsets suppressed."""
    print(f"\n{title}")
    axis_of = {"log[S]": "s0", "saturation([S])": "s0",
               "log[hoo]": "hoo", "log[buf]": "buf"}
    interesting = [n for n in fitted.names if not n.startswith("exp")]
    for name in interesting:
        value, stderr = fitted.coefficient(name)
        if axis_of.get(name) in fitted.absorbed:
            # Nearly collinear with the offsets. The standard error here is
            # about a contrast the design barely contains, so a small one is
            # not evidence; say so rather than printing "resolved".
            verdict = "ABSORBED by the offsets -- not interpretable"
        elif abs(value) > 2 * stderr:
            verdict = "resolved"
        else:
            verdict = "NOT resolved"
        print(f"    {name:18s} {value:+7.3f} +/- {stderr:.3f}   {verdict}")
    if np.isfinite(fitted.km):
        low, high = fitted.km_interval
        unresolved = "   (interval spans the grid: saturation not resolved)" \
            if not np.isfinite(low) or (low <= 1.1e-3 and high >= 2.7e2) else ""
        print(f"    {'Km':18s} {fitted.km:7.3f} mM  "
              f"[95% {low:.3f} - {high:.2f}]{unresolved}")
    offsets = [n for n in fitted.names if n.startswith("exp")]
    if offsets:
        values = np.array([fitted.coefficient(n)[0] for n in offsets])
        print(f"    per-experiment offsets: {len(offsets)}, "
              f"spread {10 ** (values.max() - values.min()):.1f}x")
    print(f"    n={fitted.n}  parameters={fitted.parameters}  "
          f"SSE={fitted.sse:.4f}  AIC={fitted.aic:.1f}  "
          f"residual x{fitted.residual_scatter:.2f}  "
          f"design condition {fitted.condition_number:.1f}")
    if fitted.absorbed:
        print(f"    WARNING: {', '.join(fitted.absorbed)} absorbed by the offsets")


def report_windows(rows):
    """The systematic that any quoted order has to carry."""
    print(f"\nWINDOW SENSITIVITY OF THE SUBSTRATE ORDER")
    print(f"  {'window':>7s} {'order':>8s} {'+/-':>7s} {'rms/noise':>10s} {'scatter':>8s}")
    for row in rows:
        print(f"  {row['fraction']:7.0%} {row['order']:+8.3f} {row['stderr']:7.3f} "
              f"{row['rms_over_noise']:10.2f} {row['scatter']:7.2f}x")
    if len(rows) > 1:
        orders = [r["order"] for r in rows]
        print(f"  -> order moves {min(orders):+.2f} to {max(orders):+.2f} across windows; "
              f"treat +/-{(max(orders) - min(orders)) / 2:.2f} as the systematic")


def compare(frame, response="v0"):
    """Power law vs saturating, pooled vs per-experiment: the four-way table."""
    results = {}
    for per_experiment in (False, True):
        label = "per-experiment" if per_experiment else "pooled"
        try:
            results[(label, "power law")] = regress(frame, response, None, per_experiment)
            results[(label, "saturating")] = profile_km(frame, response, per_experiment)
        except ValueError:
            continue
    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--group", default=None,
                        help="comma-separated substring filter on "
                             "(substrate, temperature, buffer), e.g. '4OMe-BnOH,40'")
    parser.add_argument("--catalysed", action="store_true",
                        help="summarise E0 > 0 curves instead of the enzyme-free ones")
    parser.add_argument("--all", action="store_true", help="both, pooled")
    parser.add_argument("--response", default="v0", choices=("v0", "slope_ratio"))
    parser.add_argument("--window", type=float, default=INITIAL_WINDOW)
    parser.add_argument("--drop-outliers", action="store_true",
                        help="exclude curves that look like failed cuvettes; see "
                         "experiment_outliers() for the rule")
    parser.add_argument("--save", default=None, help="write per-curve summaries to CSV")
    parser.add_argument("--save-json", default=None, help="write regression results to JSON")
    arguments = parser.parse_args()

    curves, _ = build_curves()
    if not arguments.all:
        curves = [c for c in curves
                  if (c.conditions.e0 > 0) == bool(arguments.catalysed)]
    if arguments.group:
        wanted = [t.strip().lower() for t in arguments.group.split(",")]
        curves = [c for c in curves
                  if all(t in str(c.group).lower() for t in wanted)]
    if not curves:
        print("no curves match that filter")
        return 1

    buffers = buffer_concentrations()
    frame = to_frame(summarise(curves, buffers, arguments.window))
    kind = "ALL" if arguments.all else ("CATALYSED" if arguments.catalysed
                                        else "ENZYME-FREE")
    print(f"{kind}: {len(frame)} curves, {frame.experiment.nunique()} experiments, "
          f"{frame.group.nunique()} group(s), window {arguments.window:.0%}")

    if arguments.save:
        frame.to_csv(arguments.save, index=False)
        print(f"per-curve summaries -> {arguments.save}")

    saved = {}
    for name, block in sorted(frame.groupby("group")):
        report_design(block, f"=== {name} ===")
        report_curves(block, "CURVE QUALITY")
        report_replicates(block)
        if arguments.drop_outliers:
            flagged = {(e, s) for e, s, _, _ in experiment_outliers(block)}
            if flagged:
                keep = ~block.set_index(["experiment", "sample"]).index.isin(flagged)
                block = block[keep]
                print(f"\n  dropped {len(flagged)} outlier curve(s); "
                      f"{len(block)} remain")
        if block.experiment.nunique() < 2 or len(block) < 8:
            print("\n  too few curves or experiments here to regress")
            continue
        for (scope, form), fitted in sorted(compare(block, arguments.response).items()):
            report_regression(fitted, f"REGRESSION [{scope}, {form}]")
            saved[f"{name}|{scope}|{form}"] = {
                "response": fitted.response, "names": list(fitted.names),
                "coefficients": [float(v) for v in fitted.coefficients],
                "stderrs": [float(v) for v in fitted.stderrs],
                "n": fitted.n, "sse": fitted.sse, "aic": fitted.aic,
                "km": fitted.km, "km_interval": list(fitted.km_interval),
                "residual_scatter": fitted.residual_scatter,
                "absorbed": list(fitted.absorbed)}
        block_curves = [c for c in curves if str(c.group) == name]
        report_windows(window_sensitivity(block_curves, buffers, arguments.response,
                                          drop_outliers=arguments.drop_outliers))

    if arguments.save_json:
        with open(arguments.save_json, "w") as handle:
            json.dump(saved, handle, indent=2)
        print(f"\nregression results -> {arguments.save_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
