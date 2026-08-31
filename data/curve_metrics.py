"""
Curve-level measurements, defined once.

Every number that describes the *shape* of a progress curve -- its noise, its
initial rate, where its slope peaks -- is defined here and nowhere else. The
modules that plot, fit, screen and summarise curves all import from this one.

This module exists because they did not. Six functions were defined twice
across the package and four had silently diverged; the worst was the lag
statistic, whose two copies disagreed on the verdict for 96 of 402 curves and
would have reported 21% instead of the published 34%. A duplicate here is a
wrong number in a thesis, not a style problem, and `test_no_duplicate_defs`
in test_curve_metrics.py fails if one reappears.

Pure numpy. No I/O, no dataset, no pandas -- so anything may import it.
"""
import numpy as np

# Absorbance is recorded to three decimals, so a reading's quantisation alone
# contributes this much standard deviation. Used as a floor on every noise and
# residual estimate: a curve sitting on three or four distinct levels otherwise
# reports zero noise, and a zero denominator would give it infinite weight.
ABSORBANCE_QUANTUM = 0.001
QUANTISATION_SIGMA = ABSORBANCE_QUANTUM / np.sqrt(12)

# The leading fraction of a run an initial rate is measured over, and the
# fewest points that fraction may be floored to.
INITIAL_WINDOW = 0.20
MINIMUM_WINDOW_POINTS = 5

# A curve "lags" when its slope peaks later than this fraction of the run.
# 0.15 is the threshold behind the 34% (136/402) figure in MECHANISM.md.
LAG_THRESHOLD = 0.15


def curve_noise(values, floor=QUANTISATION_SIGMA):
    """
    Point-to-point noise from the median absolute second difference.

    `floor` is a property of the SOURCE, not of the chemistry: the .txt
    exports are rounded to 0.001 AU and need QUANTISATION_SIGMA, while a
    curve read from a .rre carries ~1000x finer readings and must be given
    rre_io.RRE_SIGMA instead. Leaving the default in place on .rre data
    would report 2.89e-4 AU for curves whose real noise is 1.2e-4.

    The second difference annihilates any linear trend, so this measures a
    progress curve's scatter without being inflated by the progress itself.
    1.4826 converts a median absolute deviation to a standard deviation and the
    sqrt(6) undoes the variance the second difference introduces.
    """
    values = np.asarray(values, dtype=float)
    if len(values) < 5:
        return floor
    curvature = values[2:] - 2 * values[1:-1] + values[:-2]
    estimate = 1.4826 * np.median(np.abs(curvature)) / np.sqrt(6)
    return max(float(estimate), floor)

def peak_position(values, times, smooth_fraction=0.05):
    """
    Where the steepest point sits, as a fraction of the run, smoothed over ~5%
    of the run first.

    THIS IS THE CANONICAL LAG STATISTIC. The 34% (136/402) quoted in
    MECHANISM.md and FITTING.md is exactly `peak_position(...) > LAG_THRESHOLD`
    counted over the fittable curves. An unsmoothed variant lived in
    fit_kinetics until 2026-08-31 and reported 84 rather than 136 at the same
    threshold -- it disagreed with this one on 96 of 402 curves. Do not write
    another one; call this.
    """
    values = np.asarray(values, dtype=float)
    window = max(3, int(len(values) * smooth_fraction) | 1)
    if len(values) < window + 4:
        return np.nan
    smoothed = np.convolve(values, np.ones(window) / window, mode="valid")
    centres = times[window // 2:window // 2 + len(smoothed)]
    slope = np.gradient(smoothed, centres)
    if slope[0] <= 0 or slope.max() <= 1.05 * slope[0]:
        return 0.0
    return float((centres[np.argmax(slope)] - centres[0]) / (centres[-1] - centres[0]))

def line_fit(times, values):
    """
    Ordinary least squares line: (intercept, slope, slope stderr, residual rms).

    Time is rescaled to the window before solving. The design matrix is
    [1, t] and t runs to a few thousand seconds, so on the raw scale the
    normal equations are conditioned around 1e7 for no reason at all.

    The residual variance is floored at the quantisation variance. Readings
    come out of the instrument at three decimals, so a window whose points
    happen to sit exactly on a line has a residual of zero and would otherwise
    report an impossible standard error -- experiment 25 sample 3 rises by
    exactly 0.004 AU per reading and produced 1.5e-20 AU/s before this floor
    existed. A weighted fit would have given that curve essentially infinite
    weight on the strength of the instrument's rounding.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    scale = times[-1] - times[0]
    if scale <= 0 or len(times) < 3:
        return np.nan, np.nan, np.nan, np.nan
    design = np.column_stack([np.ones(len(times)), (times - times[0]) / scale])
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ coefficients
    degrees = max(1, len(times) - 2)
    variance = max(float(residual @ residual) / degrees, QUANTISATION_SIGMA ** 2)
    covariance = variance * np.linalg.pinv(design.T @ design)
    return (float(coefficients[0]),
            float(coefficients[1] / scale),
            float(np.sqrt(covariance[1, 1]) / scale),
            float(np.sqrt(residual @ residual / len(times))))

def line_slope(times, values):
    """(slope, stderr, rms) -- `line_fit` without the intercept."""
    _, slope, stderr, rms = line_fit(times, values)
    return slope, stderr, rms

def window_size(count, fraction):
    """How many leading points `fraction` of a curve is, floored so a slope exists."""
    return max(MINIMUM_WINDOW_POINTS, min(count, int(count * fraction)))

def initial_rate(times, values, fraction=INITIAL_WINDOW):
    """Slope over the first `fraction` of a run, in absorbance per second."""
    count = window_size(len(times), fraction)
    return line_slope(times[:count], values[:count])

# A curve accelerates when a later block of it is steeper than its first block
# by this many combined standard errors. 3.0 leaves room for the fact that the
# steepest block is a maximum over four candidates, not a single pre-chosen one.
ACCELERATION_SIGMA = 3.0

def acceleration(times, values, fraction=INITIAL_WINDOW):
    """
    How much steeper a curve gets after its start, in standard errors.

    Returns `(z, where)`. `z` is (steepest block's slope - first block's slope)
    divided by the two slopes' combined standard error; `where` is the centre
    of the steepest block as a fraction of the run. `z > ACCELERATION_SIGMA` is
    the autocatalysis signature: the curve is measurably steeper later than it
    was at the start.

    THIS IS THE STATISTIC TO USE ON SMALL-AMPLITUDE CURVES, and most of the
    in-scope block is small-amplitude -- exp 139 changes by 0.005-0.065 AU
    against 0.0006 AU of noise. `peak_position` asks where a point-wise
    gradient is largest, and on curves like those the gradient's own noise is
    a median 28% of the largest gradient in the curve, so its argmax is partly
    a draw. Slopes over blocks average that noise down and carry a standard
    error, so a verdict here can be given a threshold rather than a hope.

    The blocks are consecutive and non-overlapping, `fraction` of the run each.
    Overlapping windows would let the maximum roam over many correlated
    candidates and inflate every z.
    """
    slopes, stderrs, centres = block_slopes(times, values, fraction)
    if len(slopes) < 2:
        return np.nan, np.nan
    if not np.isfinite(slopes[0]) or not np.isfinite(stderrs[0]):
        return np.nan, np.nan
    best = 1 + int(np.argmax(slopes[1:]))
    combined = float(np.hypot(stderrs[0], stderrs[best]))
    if combined <= 0:
        return np.nan, np.nan
    return float((slopes[best] - slopes[0]) / combined), float(centres[best])


def block_slopes(times, values, fraction=INITIAL_WINDOW):
    """
    Straight-line slopes over consecutive, non-overlapping blocks of a curve.

    Returns `(slopes, stderrs, centres)`, each of length `1/fraction` or so,
    with `centres` given as fractions of the run. This is the shared spine of
    `acceleration` and `peak_rate`: both need the same blocks, and computing
    them twice from two places is how the measurements in this package drifted
    apart the first time.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    count = window_size(len(times), fraction)
    blocks = len(times) // count
    span = times[-1] - times[0]
    if blocks < 1 or span <= 0:
        return np.array([]), np.array([]), np.array([])
    slopes, stderrs, centres = [], [], []
    for index in range(blocks):
        a, b = index * count, (index + 1) * count
        slope, stderr, _ = line_slope(times[a:b], values[a:b])
        slopes.append(slope)
        stderrs.append(stderr)
        centres.append((0.5 * (times[a] + times[b - 1]) - times[0]) / span)
    return np.array(slopes), np.array(stderrs), np.array(centres)


def peak_rate(times, values, fraction=INITIAL_WINDOW):
    """
    The steepest block's slope, its standard error, and where it sits.

    On an autocatalytic curve this is the *developed* rate, and `initial_rate`
    is the rate before the catalyst has built up. The two carry different
    reaction orders and neither substitutes for the other: 25 of the 27
    substrate-ladder cuvettes in exps 139-145 are still accelerating when the
    initial window closes, so their `initial_rate` describes the induction
    period rather than the reaction.

    Returns `(slope, stderr, where)`, `where` as a fraction of the run.
    """
    slopes, stderrs, centres = block_slopes(times, values, fraction)
    if not len(slopes):
        return np.nan, np.nan, np.nan
    best = int(np.argmax(slopes))
    return float(slopes[best]), float(stderrs[best]), float(centres[best])


# The rolling window a lag time is read through, as a fraction of the run, and
# the fraction of the rise that defines "the reaction has started".
LAG_WINDOW = 0.10
LAG_LEVEL = 0.5

def rolling_slope(times, values, fraction=LAG_WINDOW):
    """
    Least-squares slope through a window of `fraction` of the run, at every
    point it fits. Returns `(centres, slopes)`.

    A window, not `np.gradient`: on this dataset a point-wise gradient's own
    noise is a median 28% of the largest gradient in the curve, so its shape is
    partly a draw. A window of a tenth of the run averages that down while
    still resolving an induction period, which 20% blocks cannot.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    count = window_size(len(times), fraction)
    if len(times) < count + 2:
        return np.array([]), np.array([])
    centres, slopes = [], []
    for start in range(0, len(times) - count + 1):
        stop = start + count
        slope, _, _ = line_slope(times[start:stop], values[start:stop])
        centres.append(0.5 * (times[start] + times[stop - 1]))
        slopes.append(slope)
    return np.array(centres), np.array(slopes)


def lag_time(times, values, fraction=LAG_WINDOW, level=LAG_LEVEL):
    """
    When the reaction reaches half its eventual speed, in SECONDS.

    Defined as the first time the rolling slope crosses
    `initial + level * (maximum - initial)`, and only for curves that pass
    the `acceleration` test. Returns nan otherwise: a curve whose slope does
    not measurably rise has no induction period to time, and on a straight
    line the crossing is set by floating-point dust.

    Reported in seconds, never as a fraction of the run. `vmax_where` is a
    fraction, and fractions are not comparable across runs of 51 and 480
    minutes: among the accelerating in-scope curves the fraction correlates
    +0.84 with run length, so it measures the schedule as much as the
    chemistry.
    """
    z, _ = acceleration(times, values)
    if not np.isfinite(z) or z <= ACCELERATION_SIGMA:
        return np.nan
    centres, slopes = rolling_slope(times, values, fraction)
    if len(slopes) < 3:
        return np.nan
    start, peak = slopes[0], slopes.max()
    if not np.isfinite(start) or peak <= start:
        return np.nan
    crossings = np.flatnonzero(slopes >= start + level * (peak - start))
    if not len(crossings):
        return np.nan
    return float(centres[crossings[0]] - centres[0])
