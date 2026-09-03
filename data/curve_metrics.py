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

# The .txt EXPORT records absorbance to three decimals, so a reading's
# quantisation alone contributes this much standard deviation. It is the floor
# on a noise or residual estimate for a curve read from an export: one sitting
# on three or four distinct levels otherwise reports zero noise, and a zero
# denominator would give it infinite weight.
#
# IT IS NOT THE FLOOR FOR EVERY CURVE. Since 2026-08-31 most readings come from
# the instrument's own .rre at ~1e-6 AU, and read_rre.RRE_SIGMA -- 1096x smaller
# -- is the floor those need. Every function here that floors takes the value as
# an argument for that reason; fit_dataset.source_floor maps a Curve.source to
# the right one. Passing the default on .rre data overstates a curve's noise by
# up to 2.4x and its slope error by a median 1.4x.
ABSORBANCE_QUANTUM = 0.001
QUANTISATION_SIGMA = ABSORBANCE_QUANTUM / np.sqrt(12)

# Bisection steps for `bubble_rate`. 80 halvings take the bracket below any
# absorbance the instrument can resolve, and the search is over a scalar, so
# the cost is 80 forward passes of a curve rather than anything to tune.
BISECTION_ROUNDS = 80

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
    `read_rre.RRE_SIGMA` instead -- `fit_dataset.source_floor` maps a
    `Curve.source` to the right one. Leaving the default in place on .rre
    data would report 2.89e-4 AU for curves whose real noise is 1.2e-4.

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

def line_fit(times, values, floor=QUANTISATION_SIGMA):
    """
    Ordinary least squares line: (intercept, slope, slope stderr, residual rms).

    Time is rescaled to the window before solving. The design matrix is
    [1, t] and t runs to a few thousand seconds, so on the raw scale the
    normal equations are conditioned around 1e7 for no reason at all.

    The residual variance is floored at `floor ** 2`, so that a window whose
    points happen to sit exactly on a line cannot report an impossible standard
    error -- experiment 25 sample 3 rises by exactly 0.004 AU per reading and
    produced 1.5e-20 AU/s before this floor existed. A weighted fit would have
    given that curve essentially infinite weight on the strength of the
    instrument's rounding.

    `floor` IS A PROPERTY OF THE SOURCE, exactly as it is in `curve_noise`, and
    the default here is the .txt export's. A .rre curve must be given
    `read_rre.RRE_SIGMA` instead -- `fit_dataset.source_floor` maps one to the
    other. The floor's only job is to stop a degenerate zero, so the right
    value is the source's own digitisation limit, which then almost never
    binds and lets the measured scatter speak.

    Until 2026-09-01 this floor was hardcoded at QUANTISATION_SIGMA for every
    curve, which is 1096x RRE_SIGMA. It bound on 52 of the 110 live two-axis
    curves and inflated their slope errors by a median 1.4x, and since
    `acceleration` divides by exactly these errors it was suppressing the
    z-scores it is measured by: the two-axis acceleration count read 48/110
    where the instrument's own readings say 51/110.
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
    variance = max(float(residual @ residual) / degrees, float(floor) ** 2)
    covariance = variance * np.linalg.pinv(design.T @ design)
    return (float(coefficients[0]),
            float(coefficients[1] / scale),
            float(np.sqrt(covariance[1, 1]) / scale),
            float(np.sqrt(residual @ residual / len(times))))

def line_slope(times, values, floor=QUANTISATION_SIGMA):
    """(slope, stderr, rms) -- `line_fit` without the intercept."""
    _, slope, stderr, rms = line_fit(times, values, floor)
    return slope, stderr, rms

def window_size(count, fraction):
    """How many leading points `fraction` of a curve is, floored so a slope exists."""
    return max(MINIMUM_WINDOW_POINTS, min(count, int(count * fraction)))

def initial_rate(times, values, fraction=INITIAL_WINDOW,
                 floor=QUANTISATION_SIGMA):
    """Slope over the first `fraction` of a run, in absorbance per second."""
    count = window_size(len(times), fraction)
    return line_slope(times[:count], values[:count], floor)

# A curve accelerates when a later block of it is steeper than its first block
# by this many combined standard errors. 3.0 leaves room for the fact that the
# steepest block is a maximum over four candidates, not a single pre-chosen one.
ACCELERATION_SIGMA = 3.0

def acceleration(times, values, fraction=INITIAL_WINDOW,
                 floor=QUANTISATION_SIGMA):
    """
    How much steeper a curve gets after its start, in standard errors.

    Returns `(z, where)`. `z` is (steepest block's slope - first block's slope)
    divided by the two slopes' combined standard error; `where` is the centre
    of the steepest block as a fraction of the run. `z > ACCELERATION_SIGMA` is
    the autocatalysis signature: the curve is measurably steeper later than it
    was at the start.

    THIS IS THE STATISTIC TO USE ON SMALL-AMPLITUDE CURVES, and most of the
    two-axis block is small-amplitude -- exp 139 changes by 0.005-0.065 AU
    against 0.0006 AU of noise. `peak_position` asks where a point-wise
    gradient is largest, and on curves like those the gradient's own noise is
    a median 28% of the largest gradient in the curve, so its argmax is partly
    a draw. Slopes over blocks average that noise down and carry a standard
    error, so a verdict here can be given a threshold rather than a hope.

    The blocks are consecutive and non-overlapping, `fraction` of the run each.
    Overlapping windows would let the maximum roam over many correlated
    candidates and inflate every z.

    `floor` reaches the two standard errors this z divides by, so it changes
    the verdict directly -- pass the curve's SOURCE floor, not the default.
    See `line_fit`.
    """
    slopes, stderrs, centres = block_slopes(times, values, fraction, floor)
    if len(slopes) < 2:
        return np.nan, np.nan
    if not np.isfinite(slopes[0]) or not np.isfinite(stderrs[0]):
        return np.nan, np.nan
    best = 1 + int(np.argmax(slopes[1:]))
    combined = float(np.hypot(stderrs[0], stderrs[best]))
    if combined <= 0:
        return np.nan, np.nan
    return float((slopes[best] - slopes[0]) / combined), float(centres[best])


# A two-segment fit needs enough points each side for either slope to mean
# anything, and the same MINIMUM_WINDOW_POINTS the block statistics use is the
# right number: below it a "slope" is two readings and their noise.
SEGMENT_MINIMUM_POINTS = MINIMUM_WINDOW_POINTS
# Above this the curve is steeper AFTER its breakpoint than before it. A real
# reaction run to partial conversion decelerates, so the whole enzyme-free
# BnOH set sits at 0.12-1.23 -- 1.5 is well clear of that and of the noise on
# it, without being so high that only a step change trips it.
SEGMENT_RATIO_STEEP = 1.5


# A two-break fit costs three more parameters than a one-break fit -- one
# breakpoint and one line -- so the same discipline applies as to the fitting
# forms: it has to earn them on an F test. Set above the nominal F(3, ~120) at
# alpha = 0.001 (about 5.8) for the same reason `summary_kinetics.TWO_PHASE_F`
# is: the residuals of a progress curve are serially correlated, and a nominal
# F over-rejects the simpler description when they are.
SEGMENT_F = 10.0


def _segment_errors(times, values, minimum):
    """
    SSE of a straight line through every contiguous stretch, as a matrix.

    error[i, j] is the residual sum of squares of the readings from i to j
    exclusive, or inf where the stretch is too short. Built from prefix sums so
    each entry is O(1): a two-break search visits O(n^2) stretch pairs, and at
    368 readings the naive route refits a line 68000 times.
    """
    count = len(times)
    zero = np.zeros(1)
    sx = np.concatenate([zero, np.cumsum(times)])
    sy = np.concatenate([zero, np.cumsum(values)])
    sxx = np.concatenate([zero, np.cumsum(times * times)])
    syy = np.concatenate([zero, np.cumsum(values * values)])
    sxy = np.concatenate([zero, np.cumsum(times * values)])
    start = np.arange(count + 1)[:, None]
    stop = np.arange(count + 1)[None, :]
    length = stop - start
    with np.errstate(divide="ignore", invalid="ignore"):
        n = np.where(length > 0, length, np.nan)
        gx = sx[stop] - sx[start]
        gy = sy[stop] - sy[start]
        cxx = (sxx[stop] - sxx[start]) - gx * gx / n
        cyy = (syy[stop] - syy[start]) - gy * gy / n
        cxy = (sxy[stop] - sxy[start]) - gx * gy / n
        error = cyy - np.where(cxx > 0, cxy * cxy / cxx, 0.0)
    error = np.where(length >= minimum, error, np.inf)
    # Rounding can push a residual sum of squares slightly negative on a
    # stretch that lies exactly on a line; it is a cost, so clamp it.
    return np.where(np.isfinite(error), np.maximum(error, 0.0), np.inf)


def _segment_slope(times, values, start, stop):
    if stop - start < 2:
        return np.nan
    x = times[start:stop]
    y = values[start:stop]
    centred = x - x.mean()
    denominator = float(centred @ centred)
    return float(centred @ (y - y.mean()) / denominator) if denominator > 0 \
        else np.nan


def segment_breaks(times, values, breaks=1, minimum=SEGMENT_MINIMUM_POINTS):
    """
    The best split of a curve into `breaks` + 1 straight stretches.

    Returns (break_times, slopes, sse). `breaks` may be 1 or 2; two is what a
    curve whose rate rises and then falls needs, and one cannot describe it --
    a single breakpoint on such a curve lands on whichever change is stronger
    and says nothing about the other. That is how the early break on the
    temperature series' 40 C curves went unreported: `segmented_fit` looked for
    one and found the late one.

    Exhaustive over breakpoints, on the O(1) stretch errors from
    `_segment_errors`, so there is no optimiser and no local minimum.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    count = len(times)
    if breaks not in (1, 2) or count < minimum * (breaks + 1):
        return (), (), np.nan
    error = _segment_errors(times, values, minimum)

    if breaks == 1:
        candidates = np.arange(minimum, count - minimum + 1)
        if not len(candidates):
            return (), (), np.nan
        total = error[0, candidates] + error[candidates, count]
        best = int(candidates[np.argmin(total)])
        return ((float(times[best]),),
                (_segment_slope(times, values, 0, best),
                 _segment_slope(times, values, best, count)),
                float(total.min()))

    first = np.arange(minimum, count - 2 * minimum + 1)
    second = np.arange(2 * minimum, count - minimum + 1)
    if not len(first) or not len(second):
        return (), (), np.nan
    grid = (error[0, first][:, None] + error[np.ix_(first, second)]
            + error[second, count][None, :])
    grid = np.where(first[:, None] + minimum <= second[None, :], grid, np.inf)
    if not np.isfinite(grid).any():
        return (), (), np.nan
    flat = int(np.argmin(grid))
    row, column = divmod(flat, grid.shape[1])
    a, b = int(first[row]), int(second[column])
    return ((float(times[a]), float(times[b])),
            (_segment_slope(times, values, 0, a),
             _segment_slope(times, values, a, b),
             _segment_slope(times, values, b, count)),
            float(grid[row, column]))


def segment_selection(times, values, threshold=SEGMENT_F,
                      minimum=SEGMENT_MINIMUM_POINTS):
    """
    One breakpoint or two? Nested, chosen on an F test.

    A k-break fit has 3k + 2 free parameters -- two per stretch plus the
    breakpoints -- so the second costs three, and F is
    (SSE1 - SSE2)/3 over SSE2/(n - 8).

    Returns a dict: `breaks`, `times`, `slopes`, `f_statistic`, `pattern`.

    READ `pattern`, NOT `breaks`. A progress curve is smooth, and three
    straight lines fit a smooth bend better than two whether or not anything
    happened -- so the F test accepts a second breakpoint on 20 of the
    temperature series' 24 curves, including the 15 C ones whose slopes simply
    rise 0.08 -> 0.15 -> 0.17. The COUNT is a statement about piecewise-linear
    approximation. The SEQUENCE OF SLOPES is the statement about the curve:
    "rise then fall" needs a maximum in the middle and cannot be produced by
    bending one way.

    Use it to say WHERE a curve changes, not to decide what it is:
    `summary_kinetics.fit_progress` is the model.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    one_times, one_slopes, one_sse = segment_breaks(times, values, 1, minimum)
    two_times, two_slopes, two_sse = segment_breaks(times, values, 2, minimum)
    result = {"breaks": 1, "times": one_times, "slopes": one_slopes,
              "f_statistic": np.nan, "pattern": _slope_pattern(one_slopes)}
    if not np.isfinite(one_sse):
        return result
    if not np.isfinite(two_sse) or two_sse <= 0:
        return result
    degrees = max(1, len(times) - 8)
    statistic = float(((one_sse - two_sse) / 3.0) / (two_sse / degrees))
    result["f_statistic"] = statistic
    if statistic > threshold:
        result.update(breaks=2, times=two_times, slopes=two_slopes)
    result["pattern"] = _slope_pattern(result["slopes"])
    return result


def _slope_pattern(slopes):
    """Name the sequence of slopes: what the curve does, not how many lines."""
    slopes = [v for v in slopes if np.isfinite(v)]
    if len(slopes) < 2:
        return "unresolved"
    rises = [b > a for a, b in zip(slopes, slopes[1:])]
    if all(rises):
        return "rising"
    if not any(rises):
        return "falling"
    if rises[0] and not rises[-1]:
        return "rise then fall"
    return "fall then rise"


def segmented_fit(times, values, floor=QUANTISATION_SIGMA,
                  minimum=SEGMENT_MINIMUM_POINTS):
    """
    The best split of a curve into two straight lines: where, and how the slope
    changes across it.

    Returns `(break_time, slope_before, slope_after, ratio)`. The breakpoint is
    chosen by exhaustive search over interior readings -- there are at most a
    few hundred candidates and the objective is not convex, so a search is both
    cheaper and safer than an optimiser.

    WHY THIS EXISTS, AND WHAT IT CATCHES THAT NOTHING ELSE HERE DOES. Every
    other shape statistic in this module compares a curve's START to its END:
    `acceleration` takes the first block against the steepest, `peak_position`
    the point-wise argmax, `late_over_early` the last fifth against the first.
    A curve that runs flat, breaks upward in the MIDDLE, and then plateaus
    defeats all three, because its first and last stretches look ordinary and
    the event is between them. Exp 65 is exactly that curve and it ranked
    213-340 of 386 on `late_over_early` -- mid-pack -- while all four of its
    cuvettes turn out to steepen by 1.8-15.9x across a break the four share to
    within 56 s. It was found by eye, off a plot, which is the argument for
    having the number.

    `ratio` is the statistic to read: slope_after / slope_before. A reaction
    followed to partial conversion DECELERATES, so healthy curves sit below 1
    -- across the enzyme-free BnOH set, 0.12 to 1.23. A ratio above
    SEGMENT_RATIO_STEEP means the curve got materially faster part way through,
    which for a background reaction at a few percent conversion needs an
    explanation outside the reaction.

    The break TIME matters as much as the ratio, and only across a run. Four
    cuvettes of one run differ only in substrate, so a break they SHARE is
    driven by something that is not the substrate -- the buffer, the peroxide,
    the cell or the instrument. See `scope.synchronised_break`.

    `floor` is passed through to nothing here: this function reports slopes and
    their ratio, not standard errors, so it takes the argument only to keep the
    signature uniform with the rest of the module. It is deliberately unused.
    """
    # One break, via the shared search. Kept as its own name because the
    # ratio it returns is what `scope.frame` and `synchronised_break` read.
    breaks, slopes, _ = segment_breaks(times, values, 1, minimum)
    if not breaks:
        return np.nan, np.nan, np.nan, np.nan
    before, after = slopes
    ratio = after / before if before else np.nan
    return breaks[0], before, after, float(ratio)


def block_slopes(times, values, fraction=INITIAL_WINDOW,
                 floor=QUANTISATION_SIGMA):
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
        slope, stderr, _ = line_slope(times[a:b], values[a:b], floor)
        slopes.append(slope)
        stderrs.append(stderr)
        centres.append((0.5 * (times[a] + times[b - 1]) - times[0]) / span)
    return np.array(slopes), np.array(stderrs), np.array(centres)


def peak_rate(times, values, fraction=INITIAL_WINDOW,
              floor=QUANTISATION_SIGMA):
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
    slopes, stderrs, centres = block_slopes(times, values, fraction, floor)
    if not len(slopes):
        return np.nan, np.nan, np.nan
    best = int(np.argmax(slopes))
    return float(slopes[best]), float(stderrs[best]), float(centres[best])


# A reading is suspect when it sits this many of the curve's own noise away
# from a local fit through its neighbours. 5 is where the isolated flags stop
# growing quickly: 0.66% of all readings are isolated single spikes at 5 sigma.
OUTLIER_SIGMA = 5.0
# Neighbours each side, and the local polynomial degree. DEGREE 2 IS LOAD
# BEARING. A local LINE is wrong wherever the curve bends, so it reads real
# kinetics as outliers -- on exp 65 sample 3 it flags the flat-to-rise
# transition at -6.5 sigma. The quadratic scores that same point -2.7 and
# leaves it alone while still catching the 11 sigma bump four points earlier.
OUTLIER_NEIGHBOURS = 4
OUTLIER_DEGREE = 2


def local_outlier_z(times, values, noise, half=OUTLIER_NEIGHBOURS,
                    degree=OUTLIER_DEGREE):
    """
    Each reading's leave-one-out residual against its neighbours, in `noise`.

    The point being scored is EXCLUDED from the fit that predicts it, so a
    genuine spike cannot drag the curve toward itself and hide. The window is
    one-sided at the ends, which is the whole reason this exists: the leading
    reading carries more leverage than any interior point, because v0 is an
    extrapolation to t = 0.

    The instrument's own first reading was the worst-behaved point in the
    archive -- 15.9% of curves beyond 5 sigma against 7.5% for the last on the
    identical test -- and since 2026-09-01 `fit_dataset.DROP_FIRST_READING`
    removes it before anything here sees a curve. On the curves this function
    is now given the gap is much smaller and has changed sign on this test:
    14.7% of leading readings are flagged at 5 sigma against 16.2% of last
    ones. Point 0 is no longer a special case statistically; it is still where
    an error hurts most.

    Returns an array of z-scores, nan where a window could not be formed.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    count = len(times)
    out = np.full(count, np.nan)
    if count < degree + 3 or noise <= 0:
        return out
    for index in range(count):
        low, high = max(0, index - half), min(count, index + half + 1)
        window = [j for j in range(low, high) if j != index]
        if len(window) < degree + 2:
            low, high = max(0, index - 2 * half), min(count, index + 2 * half + 1)
            window = [j for j in range(low, high) if j != index]
        if len(window) < degree + 2:
            continue
        coefficients = np.polyfit(times[window], values[window], degree)
        out[index] = (values[index]
                      - np.polyval(coefficients, times[index])) / noise
    return out


def isolated_outliers(times, values, noise, sigma=OUTLIER_SIGMA, **kwargs):
    """
    Indices of suspect readings that stand ALONE, and of those that do not.

    Returns `(isolated, in_runs)`. Only the first is evidence of an artefact,
    and the distinction is a timescale argument rather than a statistical one.
    At 30-60 s sampling nothing chemical can move in one interval and revert in
    the next, so a single reading out of line with both neighbours is not
    chemistry. Two or more consecutive ones are not separated from chemistry at
    all, and this dataset's striking shapes are live hypotheses -- see
    curve_screen.py, "CURVE SHAPE IS NEVER A DEFECT". Across the archive the
    split is 445 isolated against 1429 in runs, the longest run being 16
    consecutive readings. (It was 463 against 1470 before the first reading of
    every run was discarded. The DATA_VERIFICATION.md entry quoting those is
    the pre-drop measurement and stands as the evidence for the drop.)

    NOTHING HERE EXCLUDES ANYTHING. This nominates; convictions go into
    build_manifest.KNOWN_SAMPLE_EXCLUSIONS by hand, with their evidence.

    THREE KNOWN LIMITATIONS, none fatal for an advisory flag:

      masking     two adjacent spikes each sit in the other's fitting window
                  and pull it toward themselves, so the second often falls
                  under the threshold. Injected at +9 sigma each they score
                  +5.4 and +4.9.
      endpoint    a bad first reading drags its neighbour past the threshold,
                  and the pair then reads as a run rather than as one isolated
                  spike. That happened on 21 of the 86 real curves whose
                  leading reading was flagged before the first-reading drop and
                  on 8 of the 56 after it, which is why `first_point_flagged`
                  in scope.frame is taken from z[0] directly rather than from
                  membership of `isolated`.
      sharp kink  an INSTANTANEOUS change of slope is flagged, scoring -6.3 on
                  a synthetic one-reading kink. Real transitions here are
                  gradual -- exp 65 sample 3's flat-to-rise spans several
                  readings and scores -2.7 -- but a genuinely abrupt feature
                  would be nominated, which is one more reason nothing is
                  removed automatically.
    """
    z = local_outlier_z(times, values, noise, **kwargs)
    flagged = np.flatnonzero(np.isfinite(z) & (np.abs(z) > sigma))
    if not len(flagged):
        return np.array([], dtype=int), np.array([], dtype=int)
    groups = np.split(flagged, np.flatnonzero(np.diff(flagged) != 1) + 1)
    isolated = np.array([g[0] for g in groups if len(g) == 1], dtype=int)
    in_runs = np.array([i for g in groups if len(g) > 1 for i in g], dtype=int)
    return isolated, in_runs


def whole_slope(times, values, floor=QUANTISATION_SIGMA):
    """
    Least-squares slope through the ENTIRE curve. Returns `(slope, stderr)`.

    No window to choose, and it uses every point -- but it is the average rate
    over the run, not the initial one, so on a curve that decelerates it reads
    low by however much the curve bends. Quote it beside `initial_rate` rather
    than instead of it: the gap between the two IS the curvature.

    NOT `line_fit(...)[:2]`. `line_fit` returns (intercept, slope, stderr, rms),
    so slicing its first two hands back the INTERCEPT as the rate -- an
    absorbance, not a rate at all. That bug shipped in c41f459 on 2026-09-01,
    and the nonsense buffer order it produced was mistakenly written up as
    evidence that a whole-curve line is biased by deceleration. See
    DATA_VERIFICATION.md.
    """
    slope, stderr, _ = line_slope(times, values, floor)
    return slope, stderr


def model_residual(values, predicted, parameters, noise):
    """
    A fitted model's residual RMS in units of the curve's own noise.

    One definition, used for every form, so "which model fits this curve"
    is a comparison rather than an argument -- the quadratic and the burst
    form differ in parameter count and an unadjusted RMS would favour the
    larger one for free.

    This asks a DIFFERENT question from `BoundedBurstFit.bounded`, and the
    two come apart on real curves. `bounded` asks whether the data pins the
    parameter; this asks whether the model describes the data at all. Exp
    65's four cuvettes report bounded v0 on fits sitting 7-8x above their
    noise: the parameter is perfectly determined, of a form that is wrong.
    Near 1 the model is at the noise; much above it, the extrapolation to
    t = 0 is an extrapolation of the wrong shape.
    """
    values = np.asarray(values, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    degrees = max(1, len(values) - parameters)
    if not noise > 0:
        return np.nan
    return float(np.sqrt(float((values - predicted) @ (values - predicted))
                         / degrees) / noise)


def quadratic_rate(times, values, floor=QUANTISATION_SIGMA):
    """
    Initial rate from a whole-curve quadratic: `(v0, stderr, curvature_t)`.

    Fits A = c + v0 t + a t^2 over every point, so v0 is the slope at t = 0
    with no window chosen anywhere -- the answer to the objection that
    `INITIAL_WINDOW` is arbitrary. `curvature_t` is a/se(a); |t| > 3 means the
    curve measurably bends, and its SIGN says which way (negative decelerates).

    WHAT IT IS NOT. A quadratic is a two-term Taylor series, so it is only the
    right shape while the bend is gentle. On the enzyme-free BnOH curves it
    misses by 2-7x the noise wherever the deceleration is strong, and a form
    that does not fit within noise cannot have its extrapolation to t = 0
    trusted far. Read `curvature_t` and the residual before quoting v0.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    times = times - times[0]
    if len(times) < 4:
        return np.nan, np.nan, np.nan
    design = np.column_stack([np.ones(len(times)), times, times ** 2])
    beta, *_ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ beta
    degrees = max(1, len(times) - 3)
    # Floored exactly as line_fit does, and for the same reason: below the
    # source's digitisation the residual variance is quantisation, not scatter.
    variance = max(float(residual @ residual) / degrees, floor ** 2)
    covariance = variance * np.linalg.pinv(design.T @ design)
    stderr = np.sqrt(np.diag(covariance))
    curvature_t = float(beta[2] / stderr[2]) if stderr[2] > 0 else np.nan
    return float(beta[1]), float(stderr[1]), curvature_t


# The rolling window a lag time is read through, as a fraction of the run, and
# the fraction of the rise that defines "the reaction has started".
LAG_WINDOW = 0.10
LAG_LEVEL = 0.5

def rolling_slope(times, values, fraction=LAG_WINDOW,
                  floor=QUANTISATION_SIGMA):
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
        slope, _, _ = line_slope(times[start:stop], values[start:stop], floor)
        centres.append(0.5 * (times[start] + times[stop - 1]))
        slopes.append(slope)
    return np.array(centres), np.array(slopes)


def lag_time(times, values, fraction=LAG_WINDOW, level=LAG_LEVEL,
             floor=QUANTISATION_SIGMA):
    """
    When the reaction reaches half its eventual speed, in SECONDS.

    Defined as the first time the rolling slope crosses
    `initial + level * (maximum - initial)`, and only for curves that pass
    the `acceleration` test. Returns nan otherwise: a curve whose slope does
    not measurably rise has no induction period to time, and on a straight
    line the crossing is set by floating-point dust.

    Reported in seconds, never as a fraction of the run. `vmax_where` is a
    fraction, and fractions are not comparable across runs of 51 and 480
    minutes: among the accelerating two-axis curves the fraction correlates
    +0.84 with run length, so it measures the schedule as much as the
    chemistry.
    """
    z, _ = acceleration(times, values, floor=floor)
    if not np.isfinite(z) or z <= ACCELERATION_SIGMA:
        return np.nan
    centres, slopes = rolling_slope(times, values, fraction, floor)
    if len(slopes) < 3:
        return np.nan
    start, peak = slopes[0], slopes.max()
    if not np.isfinite(start) or peak <= start:
        return np.nan
    crossings = np.flatnonzero(slopes >= start + level * (peak - start))
    if not len(crossings):
        return np.nan
    return float(centres[crossings[0]] - centres[0])


# A peak this close to either end of the run is AT THE BOUNDARY, not inside it.
# Two of the block's 480-reading runs place their maximum on the last reading;
# calling that a turnover would make the amplitude a measure of how long the
# instrument was left running.
BURST_EDGE = 0.02


def burst_amplitude(times, fitted, edge=BURST_EDGE):
    """
    How far a fitted curve climbs above its own start before it turns over.

    Returns `(amplitude, time, bounded)`. `bounded` is False when the maximum
    sits within `edge` of either end of the run, and then the amplitude is a
    LOWER BOUND -- the curve had not stopped rising when the reading stopped.
    Nothing may be compared across runs on an unbounded value, because the
    two-axis block spans 3000 to 28740 s and the comparison would be of run
    lengths.

    READ OFF THE FIT, NOT THE READINGS, and for two reasons.

    The small one is noise: the maximum of 480 readings is biased upward by
    about one noise excursion, which on this block is 1.5-2.7e-4 AU against
    features of 5e-3 to 1e-1 -- up to 5% of the quantity, all of it in one
    direction.

    The large one is that the fit's PREDICTIONS are stable where its
    PARAMETERS are not. When the two exponentials are nearly degenerate the
    linear solve trades enormous opposite amplitudes between them -- exp 135
    sample 3 returns B_fast = -241 against B_slow = +303 on a curve that moves
    0.06 AU -- so `-B_fast` is not the burst and neither is `-(B_fast +
    B_slow)`: the trade leaves the curve alone and moves only the split. What
    the trade cannot touch is the value of the fitted curve at a time inside
    the window, which is what this reads.
    """
    times = np.asarray(times, dtype=float)
    fitted = np.asarray(fitted, dtype=float)
    if len(times) < 3 or len(fitted) != len(times) or not np.all(
            np.isfinite(fitted)):
        return np.nan, np.nan, False
    index = int(np.argmax(fitted))
    span = len(times) - 1
    bounded = bool(edge * span < index < (1 - edge) * span)
    return (float(fitted[index] - fitted[0]), float(times[index]), bounded)


# A downward step this many of the curve's own noise is a DETACHMENT, not a
# reading. 8 is chosen from the block's own step distribution rather than from
# a table: across the two-axis block's 28827 steps, 23 rise by more than +20
# sigma and 122 fall by more than -20 sigma, and the largest fall (-260 sigma)
# is 4.6x the largest rise. That asymmetry is the artefact's signature -- slow
# growth, sudden release -- and 8 sits well above the noise while staying below
# the smallest step the asymmetry is visible in.
BUBBLE_DROP_SIGMA = 8.0


def bubble_drops(values, noise, sigma=BUBBLE_DROP_SIGMA):
    """
    Indices `i` where the reading falls from `i` to `i + 1` by more than
    `sigma` of the curve's own noise: the detachments.

    ABSORBANCE THAT GOES AWAY WAS NEVER PRODUCT. Benzaldehyde does not
    un-form, so a fall of this size is not the reaction running backwards; it
    is something leaving the light path. In this archive that something is O2:
    the drops scale with [H2O2] (0 of 8 curves below 5 mM carry one, 5 of 5
    above 80 mM do), they need TURNOVER as well as peroxide (exps 136 and 137
    sit at 73.4 mM and carry none, being the block's two weakest runs), they
    are NOT synchronised between the cuvettes of a run, and their onset is
    delayed to a median 40% of the run while the solution supersaturates.
    See DATA_VERIFICATION.md 2026-09-03 and MECHANISM.md refs 34-35, which are
    the ketone-catalysed decomposition of peroxide to O2.

    Pass `noise` floored by the curve's SOURCE -- `fit_dataset.source_floor` --
    or a .rre curve is compared against the export's quantisation and every
    detachment in it is 1096x under-counted.
    """
    values = np.asarray(values, dtype=float)
    if len(values) < 2 or not np.isfinite(noise) or noise <= 0:
        return np.array([], dtype=int)
    return np.flatnonzero(np.diff(values) < -sigma * noise)


def monotone_bound(values):
    """
    The greatest non-decreasing function lying under the readings.

    THE ASSUMPTION-FREE BOUND. A bubble only ever ADDS absorbance and product
    only ever accumulates, so with `A_obs = A_chem + b` and `b >= 0`,
    `A_chem(t) <= min(A_obs(s) : s >= t)` -- which is this. It needs no drop
    detector, no bubble model and no free parameter, and it is an UPPER bound
    on the chemistry, so a rate read off it cannot be inflated by gas.

    It is the systematic to quote beside `debubble`, not a replacement for it:
    it also removes real curvature, so it under-reads a clean curve. Over the
    two-axis block it costs a median 0.6% of `vmax` and up to 83% on the worst
    curve, which is the point -- the spread IS the contamination.
    """
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return values.copy()
    return np.minimum.accumulate(values[::-1])[::-1]


# How much of a fall a single ADJACENT reading may undo before the fall is read
# as an instrument excursion rather than as gas leaving. Half.
#
# It is the same timescale argument `isolated_outliers` rests on, applied to
# the other artefact: at 60 s sampling a bubble cannot grow half its size in
# one interval, and absorbance that goes away because gas left the beam does
# not come back at all, let alone in one reading. The block separates cleanly
# on it -- its ten largest detachments have adjacent readings moving -0.001 to
# 0.110 of the fall, while the excursions sit at 0.8 to 1.66.
BUBBLE_RECOVERY_FRACTION = 0.5


def detachments(values, noise, sigma=BUBBLE_DROP_SIGMA,
                recovery=BUBBLE_RECOVERY_FRACTION):
    """
    The falls of `bubble_drops`, grouped into events: `(start, stop)` index
    pairs, where the reading falls from `start` to `stop`.

    ONE BUBBLE CAN TAKE MORE THAN ONE READING TO LEAVE. The instrument reads
    every 60 s and a bubble sliding out of the beam is not obliged to finish
    inside one interval, so a run of consecutive falling steps is one
    detachment, not several. Treating them separately is not a rounding
    error: `bubble_ramp` did, and because it dated each bubble's growth from
    the previous fall it gave the second of a pair a growth window of ZERO
    seconds, skipped it, and left the whole of it in the corrected curve --
    -0.0165 AU at 60 sigma on exp 144 cuvette 2, and -0.0200 at 29 sigma on
    exp 140 cuvette 4.
    """
    drops = bubble_drops(values, noise, sigma=sigma)
    if not len(drops):
        return []
    values = np.asarray(values, dtype=float)
    events = []
    start = previous = int(drops[0])
    for index in drops[1:]:
        index = int(index)
        if index == previous + 1:
            previous = index
        else:
            events.append((start, previous + 1))
            start = previous = index
    events.append((start, previous + 1))
    return [event for event in events
            if not _is_excursion(values, event, recovery)]


def _is_excursion(values, event, recovery=BUBBLE_RECOVERY_FRACTION):
    """
    Is this fall an instrument excursion rather than gas leaving?

    GAS THAT LEAVES DOES NOT COME BACK, and a bubble does not grow half its
    size in one 60 s reading. So a fall flanked by a single reading that
    climbs a comparable amount is a spike -- either the fall departs from an
    anomalously high reading, in which case it is the return off one, or it
    lands on an anomalously low one, in which case the level is back next
    reading.

    `local_outlier_z` CANNOT BE USED FOR THIS, though it is the obvious tool:
    its window spans the fall, so a genuine step change flags itself. That is
    the "sharp kink" limitation its own docstring records, and it is not
    marginal here -- exp 135 cuvette 2's 0.1196 AU detachment scores +130.
    This looks only at the two readings immediately either side, which no step
    change can make anomalous.

    Exp 149 cuvette 5 is the curve that forced it. Its two "detachments" are
    9.3 and 8.2 sigma, both instrument excursions: the first falls 0.00206 and
    the next reading climbs 0.00222 straight back, the second falls off a
    reading that is an isolated spike. Between them they set a production rate
    of 6.2e-6 AU/s, and the repair then removed 0.0097 AU from a curve that
    rose 0.0262 -- flattening a real early rise into a straight line.
    """
    start, stop = event
    drop = float(values[start] - values[stop])
    if drop <= 0:
        return True
    into = (float(values[start] - values[start - 1]) if start >= 1 else 0.0)
    out = (float(values[stop + 1] - values[stop])
           if stop + 1 < len(values) else 0.0)
    return max(into, out) > recovery * drop


def unreleased_gas(values, events):
    """
    At every reading, the gas that is still going to be seen leaving.

    ONLY GAS THE RUN WAS WATCHED TO SHED IS SUBTRACTED. A detachment is the
    whole of the evidence that the beam held gas, so at any reading the beam
    may hold at most the total of the detachments still to come -- and after
    the last one, nothing at all. The cap never blocks a detachment: at the
    reading before the jth fall it is `d_j` plus everything later, which is
    already more than `d_j` has to pay.

    THE READINGS THEMSELVES REFUTE THE ALTERNATIVE. Carrying the rate across a
    quiet tail asserts gas that nothing ever saw leave, and on 18 of the 45
    detaching curves the rate would have made more gas over that tail than the
    trace rose in total -- exp 149 cuvette 4 sheds 0.0031 AU once, 1920 s into
    an 8.0 h run, and the rate carried over the remaining 7.5 h makes 0.0736
    against a tail that rises 0.0041. Gas that accumulates detaches, and 11 of
    the 45 ran on for more than a full detachment interval without one. So
    production stopped, and the beam is empty because the last thing it did
    was shed.

    This is what keeps the reconstruction ON the readings wherever a run has
    stopped shedding. The ceiling it replaced held the tail at the most the
    beam ever carried, which put exp 149 cuvette 3 a flat 0.0022 AU below its
    own readings for 82% of the run, and exp 150 cuvette 1 below them by 99%
    of everything it rose.

    The price is stated rather than assumed: a run that ENDS holding a bubble
    it never shed keeps that bubble, so its `vmax` is the readings' and not the
    chemistry's. `scope.bubble_recovery(ends_holding=...)` is both halves.
    """
    values = np.asarray(values, dtype=float)
    owed = np.zeros(len(values))
    for start, stop in events:
        owed[:start + 1] += values[start] - values[stop]
    return owed


def bubble_profile(times, values, events, rate):
    """
    The gas held in the beam at each reading, `b(t) >= 0`.

    Gas is made at a steady `rate` -- the peroxide is in enough excess that
    its decomposition does not slow over a run -- and leaves in the whole of
    each detachment. Between detachments `b` climbs, and THREE clauses bound
    the climb.

      it may not outrun    `b` grows by at most what the reading itself gained,
      the curve            so `f = A_obs - b` can never fall across an ordinary
                           step. This is the whole of `f' >= 0`: the chemistry
                           does not run backwards, so the gas cannot grow
                           faster than the trace it rides on.
      it may not go        a detachment takes `b` down by the size of the fall
      negative             and no further.
      it may not exceed    `unreleased_gas`: the beam holds at most the total
      what leaves          of the detachments still to come, so a reading after
                           the last one carries no gas at all. A rate fixed by
                           the drops says nothing about a stretch that has
                           none, and the trace over such a stretch says the
                           production stopped.

    THE FIRST CLAUSE IS WHAT `bubble_ramp` LACKED. That model dated each
    bubble from the previous detachment and made it reach the full size of the
    drop, so a large drop after a short window implied a bubble growing faster
    than the curve: over exp 141 cuvette 3's last six readings it subtracted
    +0.0061 per reading from a trace rising +0.0018, and the "corrected" curve
    fell at every one of them -- five steps past 8 sigma, to -20. Here that
    curve's worst remaining step is -2.1 sigma.

    Carrying gas over is the other half. A detachment empties ONE bubble and
    the window may hold several -- exp 141 cuvette 3 sheds its largest drop
    after its shortest growth window (r = -0.91), which one bubble cannot do
    -- so `b` is not reset to zero at a detachment, only reduced. That is why
    the rate is set by the CUMULATIVE demand of `bubble_rate` and not drop by
    drop.
    """
    values = np.asarray(values, dtype=float)
    times = np.asarray(times, dtype=float)
    held = np.zeros(len(values))
    owed = unreleased_gas(values, events)
    # The growth increment does not depend on how much is already held, so a
    # whole stretch between detachments is a cumulative sum rather than a
    # loop. Only the detachments themselves are stepped through, and there are
    # at most a few dozen -- which matters, because `bubble_rate` bisects and
    # so calls this some tens of times per curve.
    room = np.maximum(np.diff(values), 0.0)
    growth = np.minimum(rate * np.diff(times), room)
    # `np.minimum(start + cumsum(growth), owed)` IS the saturating recursion,
    # not an approximation of it: every increment is non-negative and `owed`
    # never rises, so once the running sum meets the cap it stays at it, and
    # clipping the sum gives the same answer as clipping at every step.
    position = 0
    for start, stop in events:
        if start > position:
            held[position + 1:start + 1] = np.minimum(
                held[position] + np.cumsum(growth[position:start]),
                owed[position + 1:start + 1])
        held[start + 1:stop + 1] = np.maximum(
            held[start] - (values[start] - values[start + 1:stop + 1]), 0.0)
        position = stop
    if position < len(values) - 1:
        held[position + 1:] = np.minimum(
            held[position] + np.cumsum(growth[position:]), owed[position + 1:])
    return held


def quiet_tail(times, events):
    """
    The stretch after the last detachment, in the run's own shedding intervals.

    WHERE THE SUBTRACTED GAS ENDS AND THE STATED SYSTEMATIC BEGINS. A bubble
    that keeps growing detaches, so a run that sheds every 1900 s and then goes
    27000 s without a fall has stopped making gas, and `unreleased_gas` is
    exactly right to leave its tail alone. A run that stops shedding one
    interval before the last reading may simply have ended mid-bubble, and
    there the reconstruction keeps gas it should have removed.

    Over the two-axis block 11 of the 45 detaching curves exceed 1, and the
    four the eye picks out -- exps 149 cuvette 4 at 14.4, 150 cuvette 1 at 7.7,
    149 cuvette 2 at 5.7 and 149 cuvette 3 at 4.6 -- are the top of the list.
    The first interval is dated from the start of the run, because the first
    bubble grew from nothing.
    """
    if not events:
        return np.nan
    times = np.asarray(times, dtype=float)
    marks = [times[start] for start, _ in events]
    cadence = (marks[-1] - times[0]) / len(marks)
    if cadence <= 0:
        return np.inf
    return float((times[-1] - times[events[-1][1]]) / cadence)


def bubble_shortfall(times, values, events, rate):
    """
    The largest detachment this `rate` cannot pay for, in absorbance.

    A bubble cannot shed gas that was never made. Zero or less means every
    detachment is affordable and `A_obs - bubble_profile` is non-decreasing
    across every one of them.
    """
    if not events:
        return 0.0
    held = bubble_profile(times, values, events, rate)
    return max(float((values[start] - values[stop]) - held[start])
               for start, stop in events)


def bubble_rate(times, values, events, rounds=BISECTION_ROUNDS):
    """
    The least steady production rate that pays for every detachment, AU/s.

    THE ONE FREE PARAMETER, and it is pinned rather than fitted. Gas that
    leaves the beam was made before it left, so the rate is bounded below by
    the cumulative demand of the detachments so far; and the least such rate
    is the least gas that explains them, which makes the reconstruction an
    UPPER bound on the chemistry -- the same direction as `monotone_bound`,
    and the safe one.

    It is not a nuisance parameter. Over the two-axis block's 50 detaching
    curves it rises with the CATALYST, +1.97 +/- 0.61 per decade of `[enz]`,
    and is flat in substrate, -0.09 +/- 0.21 -- which is the peroxide
    decomposition the gas was argued to be, recovered by a fit that never saw
    either concentration.

    Returns `inf` when no rate suffices. That is one curve in the block, exp
    135 cuvette 6, whose fall is in the FIRST interval: a bubble that grew
    before the run began leaves no rise in the data to date it from, and
    `debubble` returns such a curve untouched.
    """
    if not events:
        return 0.0
    steps = np.diff(np.asarray(values, dtype=float))
    intervals = np.diff(np.asarray(times, dtype=float))
    if not len(steps) or not np.any(intervals > 0):
        return np.inf
    high = 4.0 * max(float(np.max(steps) / np.min(intervals[intervals > 0])),
                     1e-12)
    if bubble_shortfall(times, values, events, high) > 0:
        return np.inf
    low = 0.0
    for _ in range(rounds):
        middle = 0.5 * (low + high)
        if bubble_shortfall(times, values, events, middle) > 0:
            low = middle
        else:
            high = middle
    return high


def debubble(times, values, noise, sigma=BUBBLE_DROP_SIGMA):
    """
    The readings with the gas taken out. Returns `(reconstructed, events)`.

    `A_obs = f + b`: a non-decreasing chemistry and a non-negative gas that
    climbs steadily and leaves in jumps. `bubble_rate` fixes the one
    parameter, `bubble_profile` builds `b`, and this returns `A_obs - b`.

    WHAT IT GUARANTEES. `f` is non-decreasing across every detachment and
    every ordinary step, so every fall the model calls gas is gone from the
    result: over the block's 180 detachments `scope.rebuild_smoothness`'s
    `worst_at_event` is zero or above on all 44 repairable curves, the one
    exception being exp 135 cuvette 6, the first-interval case `bubble_rate`
    returns `inf` for and this leaves alone. `f` also never rises further than
    the readings do, because `b >= 0` at both ends -- which is the mass
    balance stitching breaks.

    WHAT IT DELIBERATELY DOES NOT TOUCH is a fall `detachments` rejected as an
    instrument excursion, so a repaired curve can still carry a large single
    fall. That is not a failure: absorbance that comes straight back was never
    gas leaving, and running it through a gas model is how exp 149 cuvette 5
    lost a third of a real early rise.

    AGAINST A PLANTED TRUTH it recovers `vmax` at 1.00, 0.99, 0.97 and 0.97 as
    the artefact grows from 0.25 to 2x the chemistry when each bubble empties,
    and 1.02, 1.02, 1.03, 1.08 when they empty only partly. Stitching gives
    1.15, 1.32, 1.64, 2.34 and the old segment ramp 1.01, 1.01, 1.12, 1.60.
    On a curve with no detachment it returns the readings UNCHANGED, so it
    cannot move a clean curve: the null over 19 of them is 1.000000 exactly.

    STITCHING IS THIS MODEL AT `rate = 0` -- every bubble springing into being
    full-sized at the instant it leaves. That is why it fails, and it fails by
    the amount of the artefact's whole upward half: stitched, exp 135 cuvette
    4 ends at 126% of the most absorbance its 0.219 mM of substrate could ever
    make. Read `bubble_load` before quoting a rate, and `monotone_bound` for
    the assumption-free bracket on the other side.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    events = detachments(values, noise, sigma=sigma)
    rate = bubble_rate(times, values, events)
    if not np.isfinite(rate):
        return values.copy(), events
    return values - bubble_profile(times, values, events, rate), events


def bubble_load(values, drops):
    """
    Absorbance lost to detachments, divided by the curve's net rise.

    The severity axis, and the one that decides what a curve may be used for.
    Over the two-axis block's 110 live curves: 60 carry no detachment at all,
    31 sit below 0.5 where `debubble` is unbiased, 6 between 0.5 and 1 where it
    leaves about a tenth, and 13 above 1 -- all four substrate rungs of exp 135
    at 4.0-8.7, plus inner rungs of 138, 140, 141, 142 and 150 -- where the
    rate is not measurable by any means here.

    Returns 0.0 for a curve with no detachment and nan where the net rise is
    not positive, since the ratio has no meaning on a curve that went nowhere.

    Takes the raw falling steps of `bubble_drops`, not the grouped events of
    `detachments`: the total shed is the same either way, and this is the
    severity axis rather than a count.
    """
    values = np.asarray(values, dtype=float)
    drops = np.asarray(drops, dtype=int)
    if len(values) < 2:
        return np.nan
    net = float(values[-1] - values[0])
    if net <= 0:
        return np.nan
    if not len(drops):
        return 0.0
    return float(-np.diff(values)[drops].sum() / net)
