"""
The one way to get at the in-scope experiments.

Fitting is scoped to exps 135-151 (see fit_dataset.PRIMARY_SCOPE and
FITTING.md). This module is the front door to them: it loads the curves, hangs
every derived per-curve quantity off them in one frame, and answers the
questions that keep getting asked about the block's design.

WHY THIS EXISTS. Every analysis in this repository's history has re-derived
its own initial rate, its own noise estimate, its own lag position -- and the
copies drifted, most damagingly in the lag statistic, whose two versions
disagreed on 96 of 402 curves. The measurements now live in curve_metrics and
the selection lives here. If an analysis needs a quantity this module does not
expose, ADD IT HERE rather than computing it in a script; that is the whole
point of the module.

    from scope import curves, frame, design

    cs = curves()                     # 119 Curve objects, exps 135-151
    df = frame()                      # one row per curve, derived columns filled
    design()                          # the block's design, as a table

    python data/scope.py              # print the summary
    python data/scope.py --design     # the per-experiment design table
"""
import numpy as np
import pandas as pd

from curve_metrics import (ACCELERATION_SIGMA, INITIAL_WINDOW, LAG_THRESHOLD,
                           OUTLIER_SIGMA, acceleration, initial_rate,
                           isolated_outliers, lag_time, local_outlier_z,
                           model_residual, peak_position, peak_rate,
                           quadratic_rate, segmented_fit, SEGMENT_RATIO_STEEP,
                           whole_slope)
from fit_dataset import (PRIMARY_SCOPE, PRIMARY_SCOPE_BLOCK, build_curves,
                         in_scope, source_floor)
from solution_chemistry import dominant_buffer_pair
from summary_kinetics import fit_burst_bounded

# A run's own cuvettes have to move an axis by at least this much before that
# axis counts as measured inside the run rather than across experiments.
LADDER_MINIMUM = 2.0

# Net change below this multiple of a curve's own noise is not a measurement of
# anything. Exps 150 and 151 are mostly flat by this rule. That does NOT make
# them a background: every run in PRIMARY_SCOPE carries enzyme, and their
# cuvettes carry no concentration information either (concentration_agreement
# 0.61 and 0.005). The block has no enzyme-free control -- see
# DATA_VERIFICATION.md, 2026-08-31.
LIVE_SIGNAL_NOISE_MULTIPLE = 20.0


def curves(scope=PRIMARY_SCOPE):
    """The in-scope Curve objects, in (experiment, sample) order."""
    all_curves, _ = build_curves()
    return sorted(in_scope(all_curves, scope),
                  key=lambda c: (c.experiment, c.sample))


def frame(scope=PRIMARY_SCOPE):
    """
    One row per in-scope curve, with every derived quantity already attached.

    Columns: experiment, source, sample, substrate, buffer, temperature, buf,
    pH, s0, h2o2, e0, hoo, duration_s,
    points, outliers, outliers_in_runs, first_point_z, first_point_flagged,
    noise, net, live, v0, v0_stderr, v0_rms, vmax, vmax_stderr, vmax_where,
    gain, vmax_time_s, lag_time_s, conversion, peak, lags, accel_z, accel_where, accelerates,
    late_over_early.

    `v0` is the rate before the catalyst has built up and `vmax` the rate
    after; on this block they are different measurements with different
    orders, so pick deliberately rather than reaching for v0 by habit.

    `accelerates` is the autocatalysis verdict; prefer it to `lags` on this
    block, whose curves are small-amplitude enough that the point-wise
    gradient behind `lags` is part noise. See curve_metrics.acceleration.

    Reach for this before writing a loop over `curves()`. The columns are named
    for what they are, so a question like "does the substrate order hold at low
    peroxide" is a groupby on this frame rather than a new script.
    """
    rows = []
    for curve in curves(scope):
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        # curve.noise, not a fresh curve_noise call: build_curves floors it
        # by the curve's SOURCE, and a .rre curve floored at the .txt
        # export's quantisation reports 2.4x its real noise.
        noise = curve.noise
        # The same floor has to reach the RATES, not just the noise. Every one
        # of these divides by a standard error that line_fit floors, and until
        # 2026-09-01 that floor was hardcoded at the export's quantisation for
        # every curve -- suppressing the acceleration z on the .rre data this
        # scope is entirely made of. See fit_dataset.source_floor.
        floor = source_floor(curve.source)
        net = float(values[-1] - values[0])
        v0, v0_stderr, v0_rms = initial_rate(times, values, floor=floor)
        peak = peak_position(values, times)
        accel_z, accel_where = acceleration(times, values, floor=floor)
        vmax, vmax_stderr, vmax_where = peak_rate(times, values, floor=floor)
        # Three more rate estimators, so that "does this conclusion depend on
        # how the rate was measured" is a groupby rather than an argument.
        # v0 uses the first 20% of the run, v0_whole every point with no bend
        # allowed, v0_quad every point with one bend allowed.
        v0_whole, v0_whole_stderr = whole_slope(times, values, floor=floor)
        v0_quad, v0_quad_stderr, curvature_t = quadratic_rate(
            times, values, floor=floor)
        # And the burst/lag form, the only estimator here that is a SHAPE
        # rather than a polynomial: A = c + v_ss t - B(1 - exp(-t/tau)), with
        # v0 = v_ss - B/tau. It is carried with three guards, because its v0
        # is the one number in this frame that can look confident while
        # meaning nothing:
        #
        #   v0_burst_bounded    the profile-likelihood interval over tau is
        #                       narrow enough to quote at all
        #   v0_burst_kind       burst / lag / clamped / unresolved. On a LAG
        #                       v0 is the INDUCTION rate, not the maximum, so
        #                       pooling the two into one column mixes two
        #                       different quantities.
        #   v0_burst_resid      residual RMS in units of the curve's own
        #                       noise. `bounded` asks whether the PARAMETER is
        #                       determined, not whether the MODEL fits, and
        #                       those come apart: exp 65's four cuvettes all
        #                       report bounded=True on fits that collapse to
        #                       B = 0 with tau at the floor of its grid --
        #                       a straight line wearing a four-parameter form
        #                       -- and sit at 7-8x noise.
        burst = fit_burst_bounded(times, values, noise_floor=floor)
        burst_pred = (burst.c + burst.v_ss * times
                      - burst.B * (1.0 - np.exp(-times / burst.tau)))
        burst_resid = model_residual(values, burst_pred, 4, noise)
        # The quadratic scored the same way, so the two are comparable. Its
        # own docstring has warned since it was written that it misses by
        # 2-7x the noise where the deceleration is strong; this puts a number
        # on it per curve instead of leaving it as a caveat.
        quad_design = np.column_stack([np.ones(len(times)), times, times ** 2])
        quad_beta, *_ = np.linalg.lstsq(quad_design, values, rcond=None)
        v0_quad_resid = model_residual(values, quad_design @ quad_beta, 3, noise)
        # Suspect readings. `outliers` counts the ISOLATED ones -- a single
        # reading out of line with both neighbours, which nothing chemical can
        # produce at this sampling rate. `outliers_in_runs` counts flagged
        # readings with a flagged neighbour, which may be real structure and
        # are never treated as artefacts. Neither excludes anything.
        break_time, slope_before, slope_after, break_ratio = segmented_fit(
            times, values, floor=floor)
        isolated, in_runs = isolated_outliers(times, values, noise)
        # The leading reading gets its own flag, taken from z[0] rather than
        # from membership of `isolated`: a bad leading point drags its
        # neighbour past the threshold and the pair then reads as a run. That
        # happened on 21 of the 86 flagged curves before the first-reading drop
        # and on 8 of the 56 after it.
        #
        # The flag is kept although its statistical case has weakened. The
        # instrument's own first reading was flagged on 21.4% of curves against
        # 14.7% for the last; since that reading is discarded the leading one
        # is flagged on 13.9% against 15.2%, so it is no longer the outlier
        # class it was. What survives is structural: t = 0 is where v0 is
        # extrapolated to, so a bad point there costs more than anywhere else,
        # and the run-masking above still hides it from `isolated`.
        outlier_z = local_outlier_z(times, values, noise)
        first_z = float(outlier_z[0]) if len(outlier_z) else np.nan
        buf_acid, buf_base, buf_pka = dominant_buffer_pair(
            curve.buffer, curve.pH, curve.buf)
        rows.append({
            "experiment": curve.experiment,
            "source": curve.source,
            "sample": curve.sample,
            # The block a curve belongs to, carried per row rather than
            # assumed: once `scope` is a free parameter a frame may span
            # several (substrate, temperature, buffer) cells, and a caller
            # that pools across one without meaning to has no way to notice.
            "substrate": curve.substrate,
            "buffer": curve.buffer,
            "temperature": curve.temperature,
            # Buffer CONCENTRATION, mM -- distinct from `buffer`, the salt.
            # In every enzyme-free titration this falls as `s0` rises, because
            # substrate volume displaced buffer volume; see BUFFER_CONFOUNDED.
            "buf": curve.buf,
            # The conjugate pair straddling this pH, in mM. General acid/base
            # catalysis and a peroxo-adduct route are both first order in a
            # SPECIES, not in the total -- see
            # solution_chemistry.dominant_buffer_pair, which also says why
            # this dataset cannot separate them.
            "buf_acid": buf_acid,
            "buf_base": buf_base,
            "buf_pka": buf_pka,
            "pH": curve.pH,
            "s0": curve.conditions.s0,
            "h2o2": curve.conditions.h2o2,
            "e0": curve.conditions.e0,
            "hoo": curve.conditions.hoo,
            "duration_s": float(times[-1] - times[0]),
            "points": len(times),
            "noise": noise,
            "net": net,
            "live": net > LIVE_SIGNAL_NOISE_MULTIPLE * noise,
            "v0": v0,
            "v0_stderr": v0_stderr,
            "v0_rms": v0_rms,
            "v0_whole": v0_whole,
            "v0_whole_stderr": v0_whole_stderr,
            "v0_quad": v0_quad,
            "v0_quad_stderr": v0_quad_stderr,
            "curvature_t": curvature_t,
            "v0_burst": burst.v0,
            "v0_burst_low": burst.v0_low,
            "v0_burst_high": burst.v0_high,
            "v0_burst_bounded": bool(burst.bounded),
            "v0_burst_kind": burst.kind,
            "v0_burst_resid": burst_resid,
            "v0_quad_resid": v0_quad_resid,
            "tau": burst.tau,
            "tau_resolved": bool(burst.tau_resolved),
            "outliers": len(isolated),
            "outliers_in_runs": len(in_runs),
            "first_point_z": first_z,
            "first_point_flagged": bool(np.isfinite(first_z)
                                        and abs(first_z) > OUTLIER_SIGMA),
            "vmax": vmax,
            "vmax_stderr": vmax_stderr,
            "vmax_where": vmax_where,
            # vmax/v0: how many times over the reaction sped up. accel_z
            # says whether the speed-up is real, this says how big it is --
            # z also carries the curve's noise and length, so z alone is not
            # comparable between cuvettes.
            "gain": vmax / v0 if np.isfinite(v0) and v0 > 0 else np.nan,
            # Fraction of the cuvette's substrate that has turned over by the
            # end of the run. This is the axis an autocatalysis driven by
            # PRODUCT lives on: a loop that feeds on product closes at a
            # given product/substrate ratio, not at a given absolute rate.
            # When the curve is steepest, in seconds. The FRACTION is what
            # `vmax_where` reports, and fractions are not comparable between
            # runs of 51 and 480 minutes; a mechanism predicts a time.
            "vmax_time_s": vmax_where * float(times[-1] - times[0]),
            "lag_time_s": lag_time(times, values, floor=floor),
            "conversion": (net / (curve.epsilon * curve.conditions.s0)
                           if curve.epsilon > 0 and curve.conditions.s0 > 0
                           else np.nan),
            "peak": peak,
            "lags": bool(peak > LAG_THRESHOLD) if np.isfinite(peak) else False,
            "accel_z": accel_z,
            "accel_where": accel_where,
            "accelerates": bool(accel_z > ACCELERATION_SIGMA)
            if np.isfinite(accel_z) else False,
            "late_over_early": _late_over_early(times, values),
            # WHETHER THIS CURVE CAN SHOW A BACKGROUND FEATURE AT ALL.
            # Every run is double-beam and what the reference channel omits
            # decides what the curve means (kinetics_io, DATA_VERIFICATION.md
            # 2026-08-31): an enzyme run's reference omits the ENZYME, so the
            # background appears in both beams and CANCELS, and the curve is a
            # catalytic increment. A background run's reference omits the
            # H2O2, so the curve is the raw reaction with nothing subtracted.
            # A background shape can therefore only be looked for on
            # `differential == False` rows, and comparing the two populations
            # for one is a category error -- which is how the boric probe's
            # first control set was chosen wrongly on 2026-09-01.
            "differential": bool(curve.conditions.e0 > 0),
            # The two-line split. Every other shape column here compares the
            # curve's start to its end, and a curve that breaks upward in the
            # MIDDLE and then plateaus defeats all of them -- exp 65 sat
            # mid-pack on `late_over_early` while steepening 1.8-15.9x across
            # a break its four cuvettes share to 56 s. `break_ratio` is the
            # one to read; see curve_metrics.segmented_fit and
            # `synchronised_break`.
            "break_time": break_time,
            "slope_before": slope_before,
            "slope_after": slope_after,
            "break_ratio": break_ratio,
        })
    return pd.DataFrame(rows)


def _late_over_early(times, values, fraction=0.2):
    """
    Late-window slope divided by early-window slope.

    Above 1 the curve is accelerating, which is the autocatalysis signature.
    Kept separate from `peak_position` because it answers a different question:
    peak_position asks *where* the curve is steepest, this asks whether it is
    still getting steeper at the end.
    """
    count = max(4, int(len(times) * fraction))
    early = np.polyfit(times[:count], values[:count], 1)[0]
    late = np.polyfit(times[-count:], values[-count:], 1)[0]
    return float(late / early) if early > 0 else np.nan


def ladder(axis, scope=PRIMARY_SCOPE):
    """
    How far each run moves `axis` within its own cuvettes, as a factor.

    `axis` is a column of `frame()` -- "s0", "h2o2", "hoo", "e0".
    """
    data = frame(scope)
    return data.groupby("experiment")[axis].agg(
        lambda values: float(values.max() / max(values.min(), 1e-12)))


def within_experiment_share(axis, scope=PRIMARY_SCOPE):
    """
    The fraction of an axis's log-variance that lives inside experiments.

    This is the number the scope was chosen on. Near 1 means the order in that
    axis is measured within runs and cannot be absorbed by a per-experiment
    offset; near 0 means it rests entirely on comparing one run to another.
    """
    data = frame(scope)
    logged = np.log(np.maximum(data[axis].to_numpy(dtype=float), 1e-12))
    total = logged.var()
    if total <= 0:
        return 0.0
    groups = data.experiment.to_numpy()
    within = np.average(
        [logged[groups == e].var() for e in np.unique(groups)],
        weights=[int((groups == e).sum()) for e in np.unique(groups)])
    return float(within / total)


def design(scope=PRIMARY_SCOPE):
    """One row per experiment: what it varies, how long it ran, what it saw."""
    data = frame(scope)
    rows = []
    for experiment, group in data.groupby("experiment"):
        rows.append({
            "experiment": experiment,
            "pH": group.pH.iloc[0],
            "hoo_mM": group.hoo.iloc[0],
            "cuvettes": len(group),
            "s0_ladder": group.s0.max() / max(group.s0.min(), 1e-12),
            "h2o2_ladder": group.h2o2.max() / max(group.h2o2.min(), 1e-12),
            "duration_min": group.duration_s.max() / 60.0,
            "live": int(group.live.sum()),
            "lagging": int(group.lags.sum()),
            # Live curves only. A dead curve's "acceleration" is its
            # quantisation staircase happening to step late: exp 151 has
            # one live curve and would otherwise report two accelerating.
            "accelerating": int(group.loc[group.live, "accelerates"].sum()),
            "median_late_over_early": float(group.loc[group.live,
                                                      "late_over_early"].median()),
        })
    return pd.DataFrame(rows).set_index("experiment")


def blocks(scope=PRIMARY_SCOPE):
    """
    The (substrate, temperature, buffer) cells a scope spans, with counts.

    Rate constants may be pooled only within one cell (FITTING.md F7), so a
    scope that returns more than one row here is a scope no fit may be run on
    as a unit. PRIMARY_SCOPE returns exactly one; the enzyme-free BnOH set
    returns two, phosphate and boric, which is why it is a background
    characterisation and not a fit.
    """
    data = frame(scope)
    return (data.groupby(["substrate", "temperature", "buffer"])
            .agg(curves=("experiment", "size"),
                 experiments=("experiment", "nunique"))
            .sort_values("curves", ascending=False))


def summary(scope=PRIMARY_SCOPE):
    """The scope in one paragraph of numbers, for printing."""
    data = frame(scope)
    cells = blocks(scope)
    return {
        # Read off the curves, not assumed: with `scope` a free parameter this
        # was reporting PRIMARY_SCOPE_BLOCK for every scope it was handed.
        "block": cells.index[0] if len(cells) == 1 else tuple(cells.index),
        "experiments": int(data.experiment.nunique()),
        "curves": len(data),
        "live_curves": int(data.live.sum()),
        "pH_range": (float(data.pH.min()), float(data.pH.max())),
        "hoo_decades": float(np.log10(data.hoo.max() / max(data.hoo.min(), 1e-12))),
        "within_experiment_s0": within_experiment_share("s0", scope),
        "within_experiment_h2o2": within_experiment_share("h2o2", scope),
        "lagging": int(data.lags.sum()),
        "accelerating": int(data.accelerates.sum()),
        "accelerating_live": int(data.loc[data.live, "accelerates"].sum()),
    }


# The other axis of a two-axis run is held constant along each ladder, so a
# ladder is the smallest set of cuvettes that isolates one concentration.
AXIS_PARTNER = {"s0": "h2o2", "h2o2": "s0"}


def ladder_groups(axis, scope=PRIMARY_SCOPE, live_only=True, minimum=3):
    """
    The runs of cuvettes that vary `axis` alone, as (label, sub-frame) pairs.

    Within one experiment, cuvettes sharing a value of the partner axis differ
    only in `axis` -- and in nothing else at all, since pH, [HOO-], enzyme,
    cell and day are properties of the run. A slope read along one of these is
    therefore an order, not a correlation. This is the same argument the scope
    rests on, applied to a single ladder instead of the whole block.

    `minimum` is the fewest cuvettes a ladder must have to be worth a slope.
    """
    partner = AXIS_PARTNER[axis]
    data = frame(scope)
    if live_only:
        data = data[data.live]
    groups = []
    for (experiment, held), group in data.groupby(["experiment", partner]):
        group = group.sort_values(axis)
        if len(group) < minimum or group[axis].nunique() < minimum:
            continue
        groups.append((f"exp {experiment}, {partner}={held:.4g}", group))
    return groups


def ladder_trend(parameter, axis, scope=PRIMARY_SCOPE):
    """
    A parameter's median value at the bottom and top rung of the ladders.

    For quantities that cannot be logged -- `accel_z` is a signed statistic, not
    a rate -- an order is meaningless, but the question "does this track the
    axis" still has an answer. Taking the bottom and top rung of each ladder
    keeps the comparison inside a run, so it carries the same protection from
    between-run confounding that `orders` gets from its offsets.

    Returns (low, high, n_ladders).
    """
    low, high = [], []
    for _, group in ladder_groups(axis, scope):
        ordered = group.sort_values(axis)
        low.append(float(ordered[parameter].iloc[0]))
        high.append(float(ordered[parameter].iloc[-1]))
    if not low:
        return np.nan, np.nan, 0
    # nanmedian, not median: `gain` is nan wherever v0 was not positive, and
    # one such rung would otherwise erase the whole trend.
    return float(np.nanmedian(low)), float(np.nanmedian(high)), len(low)


def local_orders(parameter, axis, scope=PRIMARY_SCOPE):
    """
    The log-log slope between each ADJACENT pair of rungs, ladder by ladder.

    `orders` fits one slope to the whole block and so cannot tell saturation
    from inhibition: both flatten the average. This returns the order as a
    function of concentration instead. Saturation approaches zero and stays
    there; inhibition goes negative at the top of the range. That distinction
    is the whole question when a curve stops responding to its substrate.

    Returns a frame with columns: experiment, axis value (geometric mean of the
    pair), order, and the pH of the run.
    """
    rows = []
    for _, group in ladder_groups(axis, scope):
        group = group.sort_values(axis)
        x = np.log(group[axis].to_numpy(dtype=float))
        y = group[parameter].to_numpy(dtype=float)
        for i in range(len(group) - 1):
            if not (y[i] > 0 and y[i + 1] > 0) or x[i + 1] <= x[i]:
                continue
            rows.append({
                "experiment": int(group.experiment.iloc[i]),
                axis: float(np.exp(0.5 * (x[i] + x[i + 1]))),
                "order": float((np.log(y[i + 1]) - np.log(y[i]))
                               / (x[i + 1] - x[i])),
                "pH": float(group.pH.iloc[i]),
            })
    return pd.DataFrame(rows)


def concentration_agreement(scope=PRIMARY_SCOPE, parameter="vmax"):
    """
    Per experiment: does the rate follow the concentrations across its cuvettes?

    A run's cuvettes share pH, [HOO-], enzyme, cell and day, so the only thing
    that may move `parameter` between them is [S] and [H2O2]. This correlates
    each run's observed log rate against the log rate the block's own fitted
    orders predict for its cuvettes. Near +1 the run obeys the rate law; near 0
    it does not, and whatever is moving its cuvettes is not the reaction.

    THIS IS THE SCREEN FOR DRIFT-DOMINATED RUNS. A cuvette measured at the level
    of the cell's own wander still produces a rate and a standard error, and
    nothing else in this package would call it out.

    Mildly circular: the orders are fitted over the same curves, so a block of
    drift-dominated runs would flatten the orders they are then judged against.
    Read it as a ranking between runs, not as an absolute score.

    Returns a frame indexed by experiment: hoo, median rate, agreement, n.
    """
    data = frame(scope)
    data = data[data.live & (data[parameter] > 0)]
    fitted = orders(parameter, scope)
    rows = []
    for experiment, group in data.groupby("experiment"):
        if len(group) < 4:
            continue
        predicted = (fitted["order_s0"] * np.log(group.s0)
                     + fitted["order_h2o2"] * np.log(group.h2o2))
        observed = np.log(group[parameter].to_numpy(dtype=float))
        if predicted.std() <= 0 or observed.std() <= 0:
            continue
        rows.append({"experiment": experiment,
                     "hoo": float(group.hoo.iloc[0]),
                     "median_rate": float(group[parameter].median()),
                     "agreement": float(np.corrcoef(observed, predicted)[0, 1]),
                     "cuvettes": len(group)})
    return pd.DataFrame(rows).set_index("experiment")


# The per-curve quantities it is meaningful to ask an order of. `v0` and `vmax`
# are rates; `net` is an extent, whose "order" is a shape statement about how
# far a cuvette gets, not a rate law; `gain` is the speed-up factor, whose
# order says how the autocatalysis itself depends on each concentration.
ORDER_PARAMETERS = ("v0", "vmax", "net", "gain")

# A run has to predict its own cuvettes this well before its orders are worth
# reading. concentration_agreement correlates each run's observed log vmax
# against the log rate its own cuvette concentrations imply; a run that scores
# low is telling you its cuvettes differ by something other than what was put
# in them -- drift, a bad blank, or a signal too small to carry the ladder.
#
# The archive separates cleanly here, which is why the threshold is not fine-
# tuned: the eleven runs above it score 0.724 to 0.974, and the five below
# score 0.609 down to 0.005. Exp 151 does not appear at all -- its cuvettes
# scatter 234-fold with two negative rates.
AGREEMENT_FLOOR = 0.70


def strong_runs(scope=PRIMARY_SCOPE, floor=AGREEMENT_FLOOR):
    """
    The experiments whose own cuvettes predict their own rates.

    Returns a sorted tuple of experiment numbers. Orders quoted in MECHANISM.md
    are measured over these; quoting them over all 17 runs moves the substrate
    order of vmax from +0.01 to +0.11 and drops the fit's R2 from 0.88 to 0.81,
    because the runs that fail this test contribute scatter and no signal.
    """
    table = concentration_agreement(scope)
    return tuple(sorted(int(e) for e in
                        table.index[table.agreement >= floor]))


def orders(parameter="v0", scope=PRIMARY_SCOPE, within=True, live_only=True):
    """
    Apparent reaction orders in [S] and [H2O2], from a log-log regression.

    Fits log(parameter) = intercept + a*log[S] + b*log[H2O2], with a free
    offset per experiment when `within` is True.

    THE OFFSETS ARE THE POINT. A run holds pH, [HOO-], enzyme batch, cell and
    day constant across its own cuvettes, so a per-experiment offset absorbs
    all of them and leaves `a` and `b` measured only from contrast between
    cuvettes of the same run. That is exactly what this scope was chosen to
    provide -- 100.0% of its log[S] variance and 94.1% of its log[H2O2]
    variance is within-experiment. Set `within=False` to see what the same
    numbers look like when between-run differences are allowed to carry the
    fit; the gap between the two is the confounding the scope exists to avoid.

    Returns a dict: order_s0, stderr_s0, order_h2o2, stderr_h2o2, n, r2.
    """
    data = frame(scope)
    if live_only:
        data = data[data.live]
    data = data[(data[parameter] > 0) & np.isfinite(data[parameter])]
    if len(data) < 4:
        return {"order_s0": np.nan, "stderr_s0": np.nan, "order_h2o2": np.nan,
                "stderr_h2o2": np.nan, "n": len(data), "r2": np.nan}

    y = np.log(data[parameter].to_numpy(dtype=float))
    columns = [np.log(data.s0.to_numpy(dtype=float)),
               np.log(data.h2o2.to_numpy(dtype=float))]
    if within:
        # One indicator per experiment, and no separate intercept -- the
        # indicators already span it. Adding both would make X rank-deficient.
        labels = np.unique(data.experiment.to_numpy())
        columns += [(data.experiment.to_numpy() == e).astype(float)
                    for e in labels]
    else:
        columns.append(np.ones(len(data)))
    design_matrix = np.column_stack(columns)

    coefficients, *_ = np.linalg.lstsq(design_matrix, y, rcond=None)
    residual = y - design_matrix @ coefficients
    degrees = max(1, len(y) - np.linalg.matrix_rank(design_matrix))
    variance = float(residual @ residual) / degrees
    covariance = variance * np.linalg.pinv(design_matrix.T @ design_matrix)
    total = float(((y - y.mean()) ** 2).sum())
    return {
        "order_s0": float(coefficients[0]),
        "stderr_s0": float(np.sqrt(covariance[0, 0])),
        "order_h2o2": float(coefficients[1]),
        "stderr_h2o2": float(np.sqrt(covariance[1, 1])),
        "n": int(len(y)),
        "r2": float(1 - (residual @ residual) / total) if total > 0 else np.nan,
    }


def order_table(scope=PRIMARY_SCOPE, parameters=ORDER_PARAMETERS):
    """Orders for each parameter, within-experiment and pooled, side by side."""
    rows = []
    for parameter in parameters:
        for within in (True, False):
            result = orders(parameter, scope, within=within)
            rows.append({"parameter": parameter,
                         "fit": "within-experiment" if within else "pooled",
                         **result})
    return pd.DataFrame(rows).set_index(["parameter", "fit"])



# ---------------------------------------------------------------------------
# The +/- chemzyme controls.
#
# Exps 65-71 are consecutive runs from June 2010 in which the same substrate
# ladder, [H2O2], buffer, pH and temperature were run twice, once without the
# chemzyme and once with it at 0.028 mM. They are the only paired controls in
# the archive on BnOH: every other enzyme-free BnOH run (exps 3 and 6) sits at
# a pH and [H2O2] no catalysed run shares.
#
# They are NOT in PRIMARY_SCOPE and cannot be pooled with it -- phosphate and
# boric buffer, and only one value of [H2O2] per run, so they carry no
# peroxide order. What they carry is the one comparison the primary scope
# cannot make at all: the same chemistry with and without the catalyst.
#
# (enzyme-free experiments, catalysed partner, label)
PAIRED_CONTROLS = (
    ((65,), 66, "boric pH 8.51"),
    ((67,), 68, "phosphate pH 8.01"),
    ((69, 70), 71, "phosphate pH 8.01, low rungs"),
)

# Enzyme-free BnOH runs whose [buf] is CONSTANT along the substrate ladder,
# and so the only ones from which a substrate order may be read.
FREE_BNOH = (65, 67, 69, 70)

# The archive holds two more enzyme-free BnOH runs, and they are a trap. Exps 3
# and 6 are buffer titrations: [buf] falls 85 -> 25 mM as [sub] rises 1.28 ->
# 8.98 mM, r = -0.91 and -0.98 against log[sub]. Their rate falls with
# substrate, which looks like the turnover the catalysed block shows and is a
# buffer effect wearing substrate's clothes -- FITTING.md F1 has said so since
# 2026-08-29. Pooling them into FREE_BNOH turns 2 clean rung-pairs above 3 mM
# into 7 and swings the median order from -0.245 to -0.431. Do not.
FREE_BNOH_BUFFER_TITRATIONS = (3, 6)

# Rungs count as the same rung if their [BnOH] agrees to this relative
# tolerance. The ladders were made from one dilution series, so matching is
# exact in practice; the tolerance only absorbs rounding in the sheets.
RUNG_TOLERANCE = 0.02


def paired_controls(controls=PAIRED_CONTROLS):
    """
    One row per substrate rung that exists both with and without the chemzyme.

    Columns: pair, s0, v0_free, v0_enz, vmax_free, vmax_enz, ratio (the
    catalysed vmax over the enzyme-free one), accel_free, accel_enz, live --
    where `live` is True only if BOTH sides carry a signal. Read `ratio` only
    on live rows: a ratio against a dead curve is a ratio against noise, and
    three of the twelve catalysed cuvettes are dead where their enzyme-free
    partner is alive.

    For scale when reading `ratio`: exps 69 and 70 are the same experiment run
    twice, and their vmax disagrees by up to 1.55x rung for rung. Nothing
    inside that factor is a measurement of anything.
    """
    scope = tuple(sorted({e for free, cat, _ in controls
                          for e in (*free, cat)}))
    data = frame(scope)
    rows = []
    for free, catalysed, label in controls:
        left = data[data.experiment.isin(free)]
        right = data[data.experiment == catalysed]
        for s0 in sorted(right.s0.unique()):
            a = left[np.isclose(left.s0, s0, rtol=RUNG_TOLERANCE)]
            b = right[np.isclose(right.s0, s0, rtol=RUNG_TOLERANCE)]
            if not len(a) or not len(b):
                continue
            b = b.iloc[0]
            both_live = bool(a.live.all() and b.live)
            rows.append({
                "pair": label, "s0": float(s0),
                "v0_free": float(a.v0.median()), "v0_enz": float(b.v0),
                "vmax_free": float(a.vmax.median()), "vmax_enz": float(b.vmax),
                "ratio": float(b.vmax / a.vmax.median()) if both_live else np.nan,
                "accel_free": float(a.accel_z.max()),
                "accel_enz": float(b.accel_z),
                "live": both_live,
            })
    return pd.DataFrame(rows)


# The same experiment run twice: exps 69 and 70 share every declared
# condition. Their disagreement is the archive's only direct measure of
# run-to-run reproducibility on BnOH, and so the yardstick any ratio has to
# clear before it means anything.
REPLICATE_PAIR = (69, 70)


def catalytic_effect(controls=PAIRED_CONTROLS, replicate=REPLICATE_PAIR):
    """
    What adding 0.028 mM chemzyme does to BnOH oxidation, on the paired runs.

    Returns the median vmax ratio over the live matched rungs, its range, the
    count, and -- as the yardstick that decides whether the ratio means
    anything -- the largest vmax disagreement between exps 69 and 70, which
    are the same experiment run twice.
    """
    table = paired_controls(controls)
    live = table[table.live]
    first, second = replicate
    data = frame(replicate)
    repeats = []
    for s0 in sorted(data[data.experiment == first].s0.unique()):
        a = float(data[(data.experiment == first) & (data.s0 == s0)].vmax.iloc[0])
        b = float(data[(data.experiment == second) & (data.s0 == s0)].vmax.iloc[0])
        repeats.append(max(a, b) / min(a, b))
    return {
        "rungs": int(len(live)),
        "median_ratio": float(live.ratio.median()),
        "ratio_range": (float(live.ratio.min()), float(live.ratio.max())),
        "replicate_scatter": float(max(repeats)),
        "dead_with_enzyme": int((~table.live).sum()),
    }


# ---------------------------------------------------------------------------
# The literature's own kinetics for this reaction, for scale.
#
# From MECHANISM.md reference 4 (ChemCatChem 2025), the only Bols-group paper
# retrieved in full. Its catalyst (diketone 8) is not necessarily the
# "a-diesterketon" exp 66's sheet names, so this is an order-of-magnitude
# comparison and nothing finer.
LITERATURE = {
    "source": "ChemCatChem 2025, diketone 8 (MECHANISM.md ref 4)",
    "kcat_per_s": 44e-5,
    "km_mM": 1.25,
    "kcat_over_kuncat": 28000,
    "catalyst_mM": 0.4,
    "h2o2_mM": 72.0,
    "pH": 7.0,
}

# The sheets' extinction coefficient for BnOH at 285 nm, used only to turn our
# AU/s into mM/s. What the absorbance actually measures is MECHANISM.md's
# leading open question, so treat any rate in mM/s as uncertain by whatever
# that answer turns out to be.
BNOH_EPSILON = 1.23

# The enzyme-free BnOH runs at NEAR-NEUTRAL pH -- the only ones comparable to
# the literature's pH 7 without extrapolating across a decade of [HOO-].
#
# These are exps 3 and 6, the same buffer titrations FREE_BNOH excludes. The
# collinearity between [buf] and [sub] destroys any SUBSTRATE ORDER read from
# them; it does not destroy the RATE of an individual cuvette, which is what
# this comparison needs. Do not use them for an order. See FITTING.md F1.
FREE_BNOH_NEUTRAL = (3, 6)


# Every enzyme-free BnOH run the manifest keeps: the four with a constant-[buf]
# substrate ladder and the two buffer titrations. Exp 64 is NOT here -- it is
# excluded as an aborted run (7 minutes at dt = 28 s), which is why the boric
# pair 64/65 cannot supply a within-pair [H2O2] contrast.
FREE_BNOH_ALL = FREE_BNOH + FREE_BNOH_NEUTRAL

# The boric-buffer runs, isolated so that "what does borate carry" is a scope
# rather than an argument. In the enzyme-free BnOH data this is exp 65 alone.
#
# MECHANISM.md says to treat boric points as suspect and gives three separate
# reasons, all of which bite hardest exactly where exp 65 sits (pH 8.51):
# borate forms PEROXOBORATE with H2O2 (K = 2.0e-8, significant above pH ~7.7),
# a much faster oxidant than H2O2 itself; it generates DIOXABORIRANE, a
# competing electrophilic oxidant with a 2.8 kcal/mol barrier; and boric acid
# CATALYSES PEROXYACID HYDROLYSIS ~12-fold with a maximum at pH 8.4-9, which
# would be actively destroying the intermediate the mechanism runs through.
#
# The data agreed before that was consulted. Exp 65 is the only run in the set
# that NEITHER the quadratic NOR the burst/lag form fits -- 4.8-5.7x noise and
# 6.9-8.3x respectively -- and its samples 3 and 4 return a negative v0_quad,
# so they silently leave any log-log fit. See DATA_VERIFICATION.md 2026-09-01.
#
# EXCLUDING IT IS NOT FREE. Exp 65 is the only run at pH 8.51 and carries the
# top of the [HOO-] range on its own (0.089 mM against 0.041 and 0.0012), so
# without it the pH axis has two levels and a [HOO-] order taken across them
# is a two-point line that cannot be checked for curvature. That is why this
# is a named alternative scope and not a deletion.
BORIC_BUFFER = (65,)
FREE_BNOH_PHOSPHATE = tuple(e for e in FREE_BNOH_ALL if e not in BORIC_BUFFER)

# EXP 65 HAS NO USABLE RATE, and every default in this module that fits one now
# excludes it. Ruled 2026-09-01. The reason is not that boric is a different
# buffer -- that would be an argument for a sensitivity, which is what
# `boric_sensitivity` was -- it is that the run's curves cannot be reduced to a
# rate at all:
#
#   * all four cuvettes break upward mid-run, at 504-560 s, and steepen by
#     1.82-15.94x across the break (`synchronised_break`). So `v0` measures the
#     pre-break stretch, `vmax` the post-break one, and they are not estimates
#     of one quantity that disagree -- they are two different quantities.
#   * `v0_quad` returns a NEGATIVE rate on two of the four, which silently
#     drops them from any log-log fit, so the "all six runs" law rested on two
#     of exp 65's four cuvettes without saying so.
#   * neither rate form fits it: 6.9-8.3x noise for the burst, 4.8-5.7x for the
#     quadratic, against a median 1.08-1.18x over the rest of the block.
#   * all four burst fits collapse to B = 0 with tau at its grid floor while
#     still reporting bounded = True.
#
# There is no defensible choice of estimator, so there is no defensible number,
# and a sensitivity that reports the fit both ways implies one of them is right.
#
# THIS IS NOT A KNOWN_EXCLUSIONS ENTRY, deliberately. That would drop exp 65
# from the dataset, and its SHAPE is the most informative thing in the boric
# block: it is the only background run in the archive that breaks, which is the
# only direct evidence bearing on whether the buffer makes an oxidant. The
# shape is evidence; the rate is not. `FREE_BNOH_ALL` therefore still exists
# and still contains exp 65, for shape work and for `boric_sensitivity`, which
# is now a record of what the exclusion bought rather than a live alternative.
BORIC_RATE_UNUSABLE = BORIC_BUFFER


def background_orders(scope=FREE_BNOH_PHOSPHATE, terms=("s0", "h2o2", "hoo"),
                      within=False, parameter="vmax", drop_accelerating=False):
    """
    How the enzyme-free rate depends on substrate, peroxide and pH.

    Fitted across all six enzyme-free BnOH runs, which span pH 6.71 to 8.51.
    Returns a dict of orders with standard errors.

    The [HOO-] order is the one that matters here: it comes out near +0.84,
    i.e. the uncatalysed reaction is close to first order in the peroxide
    anion, so its rate climbs about tenfold per pH unit and a background
    measured at pH 8 says almost nothing about a background at pH 7. Ignoring
    that is how this module first reported an 876-fold discrepancy against the
    literature where the honest figure is nearer 34.

    [H2O2] and [HOO-] are collinear (hoo = h2o2 * f(pH)), so read their sum as
    the peroxide dependence and the [HOO-] term as the pH part.

    `terms` is the list of axes to fit. IT DOES NOT DEFAULT TO EVERYTHING, and
    the omission that matters is `buf`. Exps 3 and 6 are buffer titrations --
    [buf] falls 85 -> 25 mM as [sub] rises -- so a fit without a `buf` term
    does not drop the buffer effect, it RELABELS it as a substrate order. See
    BUFFER_CONFOUNDED and `buffer_dependence`. Every returned order carries a
    `vif_` alongside it for the same reason: an order whose variance is
    inflated tenfold is a number the design cannot support, and it looks
    exactly like one that can.
    """
    data = frame(scope)
    data = data[data.live]
    # A log-log order is undefined for a non-positive rate, and some
    # estimators produce them: `v0_quad` extrapolates to t = 0 and returns a
    # negative rate on the two exp 65 cuvettes whose curvature the quadratic
    # cannot hold. Dropping them here rather than propagating a nan means the
    # `n` this function reports is the count it actually fitted.
    data = data[(data[parameter] > 0) & np.isfinite(data[parameter])]
    if drop_accelerating:
        # An "initial rate" read off a curve whose rate is still RISING is the
        # induction rate, not the reaction at the stated concentrations -- the
        # same distinction curve_metrics.peak_rate draws between v0 and vmax.
        # Four enzyme-free BnOH curves accelerate and all four are in the
        # titrations, so this is a sensitivity worth reporting rather than a
        # default: it removes 2 of the 10 live curves the buffer order's
        # titration arm rests on.
        data = data[~data.accelerates]
    y = np.log(data[parameter].to_numpy(dtype=float))
    axes = [np.log(data[c].to_numpy(dtype=float)) for c in terms]
    if within:
        # One indicator per experiment instead of a single intercept, as in
        # `orders`. Anything constant across a run -- pH, [HOO-], [H2O2], the
        # cell, the day -- is absorbed, so the remaining orders are measured
        # only from contrast BETWEEN CUVETTES OF THE SAME RUN. A term that is
        # itself constant within every run is then unidentifiable and must not
        # be in `terms`; it would be collinear with the indicators.
        labels = np.unique(data.experiment.to_numpy())
        design_matrix = np.column_stack(
            axes + [(data.experiment.to_numpy() == e).astype(float)
                    for e in labels])
    else:
        design_matrix = np.column_stack([np.ones(len(data))] + axes)
    coefficients, *_ = np.linalg.lstsq(design_matrix, y, rcond=None)
    residual = y - design_matrix @ coefficients
    dof = max(1, len(data) - np.linalg.matrix_rank(design_matrix))
    variance = float(residual @ residual) / dof
    stderr = np.sqrt(np.diag(variance
                             * np.linalg.pinv(design_matrix.T @ design_matrix)))
    result = {"n": int(len(data)), "terms": tuple(terms), "within": bool(within),
              "dof": int(dof),
              "r2": float(1 - (residual @ residual)
                          / ((y - y.mean()) ** 2).sum())}
    # The orders are the LEADING coefficients in both layouts: `within` puts
    # the indicators after them, the pooled fit puts the intercept first.
    offset = 0 if within else 1
    inflation = variance_inflation(data, terms, within=within)
    for name, value, error in zip(terms, coefficients[offset:], stderr[offset:]):
        result[f"order_{name}"] = float(value)
        result[f"stderr_{name}"] = float(error)
        result[f"vif_{name}"] = inflation[name]
    return result


def variance_inflation(data, terms, within=False):
    """
    Each term's variance inflation factor against the others, on log axes.

    VIF = 1/(1 - R2) where R2 is from regressing one term on the rest. It is
    the factor by which collinearity widens that coefficient's standard error,
    so it says whether an order is measured or merely reported. Above about 10
    the coefficient is arithmetic, not evidence.
    """
    logged = {name: np.log(data[name].to_numpy(dtype=float)) for name in terms}
    result = {}
    for name in terms:
        others = [logged[o] for o in terms if o != name]
        if within:
            labels = np.unique(data.experiment.to_numpy())
            others = others + [(data.experiment.to_numpy() == e).astype(float)
                               for e in labels]
        elif not others:
            result[name] = 1.0
            continue
        matrix = np.column_stack([np.ones(len(data))] + others)
        target = logged[name]
        coefficients, *_ = np.linalg.lstsq(matrix, target, rcond=None)
        residual = target - matrix @ coefficients
        total = float(((target - target.mean()) ** 2).sum())
        r2 = 1 - float(residual @ residual) / total if total > 0 else 0.0
        result[name] = float(1.0 / (1.0 - r2)) if r2 < 1 else float("inf")
    return result


# The two enzyme-free BnOH designs, named by what they can and cannot measure.
#
# BUFFER_CONFOUNDED (exps 3, 6) move [buf] 85 -> 25 mM DOWN as [sub] moves
# 1.28 -> 8.98 mM UP, because substrate volume displaced buffer volume in the
# cuvette. BUFFER_FIXED (exps 67, 69, 70; exp 65 left it on 2026-09-01,
# see BORIC_RATE_UNUSABLE) holds [buf] at 85 mM along
# their substrate ladder. Neither design varies [buf] at constant [sub] --
# no enzyme-free run in the archive does, in any block -- so the buffer order
# is not directly measurable and has to be recovered from the disagreement
# between these two, which is what `buffer_dependence` does.
BUFFER_CONFOUNDED = FREE_BNOH_NEUTRAL
BUFFER_FIXED = tuple(e for e in FREE_BNOH if e not in BORIC_RATE_UNUSABLE)


def buffer_dependence(anchor=BUFFER_FIXED, titration=BUFFER_CONFOUNDED,
                      parameter="vmax", drop_accelerating=False):
    """
    The order in buffer concentration, from the substrate order it corrupts.

    THE DESIGN PROBLEM. No enzyme-free run varies [buf] at constant [sub], so
    no single run measures a buffer order. Two things are measurable instead:

      a   the substrate order where [buf] is CONSTANT (`anchor`), and
      a'  the substrate order where [buf] falls as [sub] rises (`titration`),

    both read within experiments, so both are free of pH, [H2O2], cell and day.
    If the rate goes as [sub]^a [buf]^d, and within the titration runs
    log[buf] = g log[sub] + constant, then fitting the titrations without a
    buffer term returns a' = a + d g, and so

        d = (a' - a) / g.

    WHY NOT JUST FIT BOTH TERMS. Because the titrations cannot carry it: on
    those eight live curves [sub] and [buf] run at VIF 11-14, and the joint fit
    returns d = +0.73 +/- 0.62, consistent with anything from zero to two.
    Pooling all six runs instead does return a tight d, but there [buf] is 85+
    in every pH 8.0-8.5 run and sweeps only in the pH 6.71 ones, so [buf] is
    partly a label for pH -- and the [HOO-] order drops from +0.84 to +0.74
    when the buffer term is added, which is that theft made visible. This
    route uses only within-run contrast on both sides and never asks one
    regression to separate the two.

    THE ASSUMPTION IT RESTS ON, stated because it is not testable here: that
    the substrate order is the same at the titrations' pH 6.71 as at the
    anchor's pH 8.01-8.51. The archive has no run that could check it.

    Returns a dict: order_s0_fixed, order_s0_titration, coupling (g),
    order_buf and its standard error, and the counts behind each.
    """
    clean = background_orders(anchor, terms=("s0",), within=True,
                              parameter=parameter,
                              drop_accelerating=drop_accelerating)
    dirty = background_orders(titration, terms=("s0",), within=True,
                              parameter=parameter,
                              drop_accelerating=drop_accelerating)

    # g: how [buf] tracks [sub] inside the titration runs, on log axes, with a
    # free offset per run -- the same within-experiment contrast the orders use.
    data = frame(titration)
    data = data[data.live]
    labels = np.unique(data.experiment.to_numpy())
    design_matrix = np.column_stack(
        [np.log(data.s0.to_numpy(dtype=float))]
        + [(data.experiment.to_numpy() == e).astype(float) for e in labels])
    coefficients, *_ = np.linalg.lstsq(
        design_matrix, np.log(data.buf.to_numpy(dtype=float)), rcond=None)
    coupling = float(coefficients[0])

    gap = dirty["order_s0"] - clean["order_s0"]
    # g is measured from concentrations that were pipetted, not from a noisy
    # observable, so its own error is negligible beside the two rate orders'.
    error = float(np.hypot(clean["stderr_s0"], dirty["stderr_s0"]) / abs(coupling))
    return {
        "order_s0_fixed": clean["order_s0"],
        "stderr_s0_fixed": clean["stderr_s0"],
        "n_fixed": clean["n"],
        "order_s0_titration": dirty["order_s0"],
        "stderr_s0_titration": dirty["stderr_s0"],
        "n_titration": dirty["n"],
        "coupling": coupling,
        "order_buf": float(gap / coupling),
        "stderr_buf": error,
    }


# The 4OMe-BnOH / 40 C enzyme-free runs. Exp 31 is deliberately NOT here: it
# is the same design at 35 C, and temperature moves every rate constant through
# Arrhenius, so including it would pool two cells (FITTING.md F7).
FREE_4OME_40C = (23, 24, 25, 26, 27, 28, 29, 30, 38, 39)


def boric_sensitivity(estimators=("v0_quad", "v0_burst", "vmax", "v0",
                                  "v0_whole"),
                      terms=("s0", "h2o2", "hoo")):
    """
    Every order in the enzyme-free BnOH set, with and without the boric run.

    Returns a DataFrame indexed by (estimator, term) with columns `all`,
    `all_stderr`, `phosphate`, `phosphate_stderr`, plus a `buf` pair from
    `buffer_dependence` -- the buffer order is computed differently, from the
    contrast between the fixed-buffer and titration designs, so it cannot come
    out of the same regression.

    THE POINT IS THE SPREAD ACROSS ESTIMATORS, not any single row. Five ways of
    measuring a rate should give one order; where they do not, the disagreement
    is the measurement's, not the chemistry's. On the substrate order the five
    span 0.24 in log units with exp 65 in and 0.06 without it -- so the boric
    run was the disagreement. The buffer order, by contrast, barely notices it,
    which is the robustness the headline actually rests on.

    See BORIC_BUFFER for why excluding it is defensible and what it costs.
    """
    rows = []
    for estimator in estimators:
        both = {}
        for label, scope_ in (("all", FREE_BNOH_ALL),
                              ("phosphate", FREE_BNOH_PHOSPHATE)):
            both[label] = background_orders(scope_, terms=terms,
                                            parameter=estimator)
        for term in terms:
            rows.append({
                "estimator": estimator, "term": term,
                "all": both["all"].get(f"order_{term}", np.nan),
                "all_stderr": both["all"].get(f"stderr_{term}", np.nan),
                "phosphate": both["phosphate"].get(f"order_{term}", np.nan),
                "phosphate_stderr": both["phosphate"].get(f"stderr_{term}",
                                                          np.nan),
            })
        # The buffer order comes from the two-design contrast, so its two
        # versions are the anchor WITH and WITHOUT exp 65. Both anchors are
        # named explicitly: BUFFER_FIXED stopped containing exp 65 on
        # 2026-09-01, and reading the "all" column off the default silently
        # made both columns the same number.
        full = buffer_dependence(anchor=FREE_BNOH, parameter=estimator)
        cut = buffer_dependence(anchor=BUFFER_FIXED, parameter=estimator)
        rows.append({
            "estimator": estimator, "term": "buf",
            "all": full["order_buf"], "all_stderr": full["stderr_buf"],
            "phosphate": cut["order_buf"], "phosphate_stderr": cut["stderr_buf"],
        })
    table = pd.DataFrame(rows).set_index(["estimator", "term"])
    return table


def boric_spread(table=None):
    """
    How far the five estimators disagree on each order, with boric and without.

    Returns a DataFrame indexed by term with columns `all`, `phosphate` --
    the max-minus-min across estimators. This is the summary that decides
    whether excluding boric helped: a spread that shrinks means the estimators
    were disagreeing about exp 65 rather than about the chemistry.
    """
    if table is None:
        table = boric_sensitivity()
    spread = table.groupby("term").agg(
        all=("all", lambda v: float(v.max() - v.min())),
        phosphate=("phosphate", lambda v: float(v.max() - v.min())))
    return spread


def buffer_cross_check(scope=FREE_4OME_40C):
    """
    The buffer order again, on the 4OMe-BnOH / 40 C block, independently.

    That block is the only other place with buffer contrast: three substrate
    values recur at two buffer levels each across its experiments (0.38 mM at
    75 and 90, 2.06 at 80 and 85, 4.12 at 70 and 75). The contrast is BETWEEN
    experiments, so this fit takes no per-experiment offsets -- offsets would
    absorb the very thing being measured. pH (6.97-7.00) and [H2O2] (82.5 mM)
    are constant across the block, so there is no pH for [buf] to proxy.

    Different substrate and temperature, so this checks the SIGN and rough
    SIZE of `buffer_dependence`, not its value.
    """
    return background_orders(scope, terms=("s0", "buf"), within=False)

# The one matched pair in the enzyme-free archive that changes the buffer SALT
# and almost nothing else. Exps 65 and 67 ran the same substrate ladder
# (7.310 / 3.655 / 1.827 / 0.365 mM) at the same [H2O2] (122.426 mM), the same
# temperature, the same instrument, and the same .rre source, at 87.5 against
# 85.0 mM buffer. Only the salt -- boric against phosphate -- and the pH
# (8.51 against 8.01) differ.
PEROXO_PAIR = (65, 67)
# The SAME pair with enzyme, and the reason it is needed: exp 65's curves have
# a synchronised mid-run break that no other run in the block has, so its rate
# numbers do not describe one process (`synchronised_break`,
# DATA_VERIFICATION.md 2026-09-01). Exps 66 and 68 match on enzyme (0.028 mM),
# buffer (85.0 mM), peroxide (122.426 mM), substrate (2.741 mM BnOH) and
# temperature, differing only in salt and pH, and both run smooth.
#
# WITHDRAWN AS A PEROXO PROBE, same day, one hour later. Both runs are
# catalysed, so their reference channels omit the ENZYME and the background is
# subtracted out of both. A buffer-made oxidant acts on the background, which
# means it is present in BOTH beams of each run and CANCELS. This pair cannot
# detect the thing it was added to detect -- not confounded, blind. It is kept
# only as the catalysed boric-vs-phosphate comparison it actually is, which is
# a statement about the CATALYSED reaction in the two buffers.
CATALYSED_PEROXO_PAIR = (66, 68)


def synchronised_break(scope=FREE_BNOH_ALL):
    """
    Per run: do its cuvettes break at the same TIME, and do they steepen?

    A run's four cuvettes differ only in substrate -- same buffer, same
    peroxide, same cell, same day. So a breakpoint they SHARE cannot be driven
    by the substrate, and a break that is also a STEEPENING cannot be the
    reaction decelerating toward conversion. The two together are the
    signature this function reports.

    Returns a DataFrame indexed by experiment: `span` (max - min break time,
    in seconds), `median_break`, `max_ratio`, `steep` (how many cuvettes
    exceed SEGMENT_RATIO_STEEP), `n`, and the sorted break times and ratios.

    ONLY ASK IT OF BACKGROUND RUNS. A shape in the background is invisible in a
    catalysed run BY CONSTRUCTION: an enzyme run's reference channel omits the
    enzyme, so the background is in both beams and cancels, and the curve is a
    catalytic increment (`frame`'s `differential` column, kinetics_io,
    DATA_VERIFICATION.md 2026-08-31). This function raises if `scope` mixes the
    two, because the first control set chosen for the boric probe was four
    catalysed runs and the conclusion drawn from their smoothness -- that exp
    65's break was not borate chemistry -- did not follow.

    WHAT IT FOUND, on the 20 background experiments -- the whole un-subtracted
    population -- 18 of which have live curves. Seventeen are phosphate and
    none of them breaks: ratios 0.22-1.23, with one isolated cuvette of exp
    3's six at 2.44. Exactly one boric run has live curves -- exp 65 -- and
    all four of its cuvettes steepen, by 1.82,
    2.04, 5.59 and 15.94, across breaks spanning 56 s, two of its 28 s
    sampling intervals. The steepening is LARGEST at the lowest substrate
    (15.94 at 0.365 mM against 1.82 at 7.310 mM), so the clock is not the
    substrate.

    That is 1 boric run of 1 showing it against 17 of 17 not, which is
    consistent with borate chemistry and rests on a single run. The other
    boric background run, exp 64, was aborted at 448 s -- BEFORE exp 65's
    break -- and is dead besides, so it cannot test it. The missing experiment
    is a repeat of exp 65. See DATA_VERIFICATION.md 2026-09-01.
    """
    data = frame(scope)
    data = data[data.live & np.isfinite(data.break_time)]
    if data.differential.nunique() > 1:
        mixed = sorted(data[data.differential].experiment.unique())
        raise ValueError(
            "synchronised_break was given both catalysed and background runs "
            f"({mixed} are catalysed). A catalysed curve is an increment whose "
            "reference channel already subtracted the background, so it cannot "
            "show a background shape and its smoothness is not evidence about "
            "one. Pass background runs only.")
    rows = []
    for experiment, group in data.groupby("experiment"):
        order = group.sort_values("break_time")
        rows.append({
            "experiment": int(experiment),
            "n": int(len(group)),
            "span": float(group.break_time.max() - group.break_time.min()),
            "median_break": float(group.break_time.median()),
            "max_ratio": float(group.break_ratio.max()),
            "steep": int((group.break_ratio > SEGMENT_RATIO_STEEP).sum()),
            "breaks": [int(v) for v in order.break_time],
            "ratios": [round(float(v), 2) for v in order.break_ratio],
        })
    return pd.DataFrame(rows).set_index("experiment")


def peroxo_buffer_test(pair=PEROXO_PAIR, orders_scope=FREE_BNOH_PHOSPHATE,
                       estimators=("v0", "vmax", "v0_whole", "v0_quad")):
    """
    Does a buffer that DOES form a peroxo species run faster than the law?

    THE QUESTION. The enzyme-free rate is first order in buffer
    (`buffer_dependence`), and two mechanisms give that: general acid/base
    catalysis, or the buffer making an oxidant -- phosphate + H2O2 ->
    peroxomonophosphate. Within the phosphate runs those are indistinguishable
    (see background_reaction/ANALYSIS.md section 6b: log[buf], log[H2PO4-] and
    log[HPO4^2-] are the same variable, correlation 1.000000).

    THE WAY ROUND IT. Borate is the buffer where the peroxo route is not a
    hypothesis. MECHANISM.md item 39 has B(OH)3 + H2O2 -> peroxoborate with
    K = 2.0e-8, "significant above pH ~ 7.7", and the anionic peroxoborates
    are much faster oxidants than H2O2 itself. Exp 65 is boric buffer at
    pH 8.51 -- above that threshold, at 122 mM H2O2 -- so a substantial part
    of its boron is peroxoborate. If a buffer-derived peroxo oxidant is what
    carries a first-order buffer term, exp 65 must run far above a rate law
    fitted without one. It does not.

    READ THIS BEFORE QUOTING THAT. On the default pair the test is WEAK, and
    it is weak for a reason found after it was written: exp 65's four cuvettes
    share a mid-run breakpoint at 504-560 s across which every one of them
    STEEPENS, by 1.82 to 15.94x, most at the LOWEST substrate. No other run in
    the block does this -- every other one decelerates -- so a single rate
    number does not describe exp 65 at all: `vmax` reads the post-break
    stretch and `v0` the pre-break one, and they are not the same process. The
    excesses below are therefore a comparison between two different things,
    and the right reading of them is "borate is nowhere fast", not "borate
    matches the law".

    AND THERE IS NO SECOND PROBE. CATALYSED_PEROXO_PAIR was added as one and
    withdrawn the same day: both its runs are catalysed, so a buffer-made
    oxidant sits in both beams and cancels. Among the 21 enzyme-free
    background experiments -- the entire population in which a background
    feature is even visible -- exactly one is boric, and it is exp 65. So this test rests on
    one run and that run is the one with the break. See `synchronised_break`
    and DATA_VERIFICATION.md 2026-09-01. The missing experiment is a repeat of
    exp 65: enzyme-free, boric, run long.

    HOW THE PREDICTION IS MADE. The rate law is fitted on `orders_scope`,
    which excludes the boric run, so exp 65 is out of sample. Because the two
    runs share [S] and [H2O2] exactly, the predicted ratio depends only on the
    [buf] and [HOO-] orders; the substrate and peroxide orders, the two worst
    determined here, drop out. Matching is per cuvette on `s0`, not on run
    medians.

    WHAT IT DOES NOT SETTLE. pH is not matched, so the prediction leans on the
    [HOO-] order holding from 8.01 up to 8.51 -- an extrapolation, and exp 65
    is the only run there. It is one run of four cuvettes against one run of
    four, on a day and a cell that are not controlled. Exp 65 is also the run
    neither rate form fits (section 3a) and its noise runs 1.5-2.8x exp 67's.
    A null result here is evidence against the peroxo route, not proof, and it
    says nothing directly about phosphate: it says the mechanism does not show
    up where it certainly operates.

    Returns a DataFrame indexed by estimator: `predicted` (boric/phosphate
    from the phosphate-only law), `observed` (median over matched cuvettes),
    `excess` = observed / predicted, `n` matched cuvettes, and the per-cuvette
    ratios in `ratios`.
    """
    boric_experiment, phosphate_experiment = pair
    data = frame(tuple(sorted(set(pair))))
    data = data[data.live]
    left = data[data.experiment == boric_experiment]
    right = data[data.experiment == phosphate_experiment]

    rows = []
    for estimator in estimators:
        if orders_scope is None:
            # UNCORRECTED, for the catalysed pair. The enzyme-free rate law
            # does not describe a catalysed comparison and there is no
            # catalysed law in this buffer to put in its place, so quote the
            # raw ratio against the [HOO-] ratio and let the reader correct:
            # boric carries 2.19x the hydroperoxide, so anything at or below
            # 1.00 is boric running SLOWER than pH alone would give it.
            law = {"order_buf": 0.0, "order_hoo": 0.0}
        else:
            law = background_orders(orders_scope, parameter=estimator,
                                    terms=("s0", "h2o2", "hoo", "buf"))
        # [S] and [H2O2] are identical between the two runs, so only these two
        # terms survive the ratio. Asserted rather than assumed: a pair that
        # did not match on them would need the other two orders as well.
        predicted = ((float(left.buf.median()) / float(right.buf.median()))
                     ** law["order_buf"]
                     * (float(left.hoo.median()) / float(right.hoo.median()))
                     ** law["order_hoo"])
        ratios = {}
        for s0 in sorted(set(left.s0) & set(right.s0), reverse=True):
            # A log-log order is undefined for a non-positive rate, and
            # `v0_quad` returns one on two of exp 65's cuvettes: skip the
            # cuvette rather than propagate a nan into the median.
            a = left[left.s0 == s0][estimator].to_numpy(dtype=float)
            b = right[right.s0 == s0][estimator].to_numpy(dtype=float)
            if len(a) != 1 or len(b) != 1 or not (a[0] > 0 and b[0] > 0):
                continue
            ratios[float(s0)] = float(a[0] / b[0])
        observed = float(np.median(list(ratios.values()))) if ratios else np.nan
        rows.append({"estimator": estimator,
                     "order_buf": law["order_buf"],
                     "order_hoo": law["order_hoo"],
                     "predicted": float(predicted),
                     "observed": observed,
                     "excess": observed / float(predicted),
                     "n": len(ratios),
                     "ratios": ratios})
    return pd.DataFrame(rows).set_index("estimator")


def literature_comparison(scope=FREE_BNOH_NEUTRAL,
                          orders_scope=FREE_BNOH_PHOSPHATE,
                          orders_terms=("s0", "h2o2", "hoo")):
    """
    Our enzyme-free background against the literature's uncatalysed rate.

    Both are put at the literature's conditions -- pH 7.0, 72 mM H2O2 -- by
    scaling our rates with the orders `background_orders` measures. Using the
    pH 8.0-8.5 runs instead and correcting only [H2O2] inflates the answer
    about 25-fold, which is a statement about [HOO-] and not about our cuvettes.

    Returns a frame with one row per live cuvette plus a `summary` attribute
    holding the median excess and the enhancement the literature's kcat would
    produce at OUR catalyst loading -- which is the number that explains why no
    enhancement is visible anywhere in this archive.
    """
    orders = background_orders(orders_scope, terms=orders_terms)
    h2o2_order, hoo_order = orders["order_h2o2"], orders["order_hoo"]
    kuncat = LITERATURE["kcat_per_s"] / LITERATURE["kcat_over_kuncat"]

    data = frame(scope)
    rows = []
    for _, row in data[data.live].iterrows():
        # to the literature's [H2O2], then across the pH gap through [HOO-]
        factor = ((LITERATURE["h2o2_mM"] / row.h2o2) ** h2o2_order
                  * (10 ** (LITERATURE["pH"] - row.pH)) ** hoo_order)
        ours = float(row.vmax) / BNOH_EPSILON * factor
        theirs = kuncat * float(row.s0)
        rows.append({"experiment": int(row.experiment), "pH": float(row.pH),
                     "s0": float(row.s0), "ours_mM_s": ours,
                     "literature_mM_s": theirs, "excess": ours / theirs})
    table = pd.DataFrame(rows)
    return table


def background_model(scope=FREE_BNOH_NEUTRAL,
                     orders_scope=FREE_BNOH_PHOSPHATE,
                     orders_terms=("s0", "h2o2", "hoo")):
    """
    An amplitude and three orders that predict the enzyme-free rate, in mM/s.

    The orders come from all six enzyme-free runs (`background_orders`); the
    amplitude is anchored on the near-neutral ones, so the model is pinned
    where it is compared to the literature and extrapolated -- not fitted --
    across the pH gap to the catalysed runs.
    """
    orders = background_orders(orders_scope, terms=orders_terms)
    data = frame(scope)
    data = data[data.live]
    exponents = (orders["order_s0"], orders["order_h2o2"], orders["order_hoo"])
    predicted = (data.s0 ** exponents[0] * data.h2o2 ** exponents[1]
                 * data.hoo ** exponents[2])
    amplitude = float(np.median(data.vmax / BNOH_EPSILON / predicted))
    return amplitude, exponents


# The catalysed runs that have an enzyme-free counterpart to be judged against:
# the three paired-control partners plus exps 73 and 83.
CATALYSED_WITH_BACKGROUND = (66, 68, 71, 73, 83)


def predicted_enhancement(scope=CATALYSED_WITH_BACKGROUND,
                          background_scope=FREE_BNOH_NEUTRAL,
                          orders_scope=FREE_BNOH_PHOSPHATE):
    """
    What the literature's kcat would show at THIS archive's catalyst loading.

    For every live cuvette of the catalysed runs that have an enzyme-free
    counterpart -- exps 66, 68, 71 (the paired controls) and 73, 83 -- this
    predicts the background at that cuvette's own conditions, adds the
    catalytic contribution kcat*E0*S/(Km+S), and reports the ratio to the
    background alone. That ratio is what the experiment would have had to
    resolve.

    It comes out at a median 1.3x, range about 1.15-1.9x. Exps 69 and 70 are
    the SAME experiment run twice and their vmax disagrees by up to 1.55x. So
    the enhancement these runs were capable of detecting is smaller than their
    own reproducibility, and the observed 0.63x is not evidence about the
    catalyst. At the literature's own 0.4 mM the predicted ratio is above 40x,
    which nothing could miss; no BnOH run in this archive exceeds 0.069 mM.
    """
    amplitude, (a, b, c) = background_model(background_scope, orders_scope)
    rows = []
    for experiment in scope:
        data = frame((experiment,))
        for _, row in data[data.live].iterrows():
            background = amplitude * row.s0 ** a * row.h2o2 ** b * row.hoo ** c
            catalysed = (LITERATURE["kcat_per_s"] * row.e0 * row.s0
                         / (LITERATURE["km_mM"] + row.s0))
            rows.append({
                "experiment": experiment, "pH": float(row.pH),
                "s0": float(row.s0), "e0": float(row.e0),
                "observed_mM_s": float(row.vmax) / BNOH_EPSILON,
                "background_mM_s": float(background),
                "catalysed_mM_s": float(catalysed),
                "expected_ratio": float((background + catalysed) / background),
                "at_literature_loading": float(
                    (background + catalysed * LITERATURE["catalyst_mM"] / row.e0)
                    / background),
            })
    return pd.DataFrame(rows)


# The scopes worth a name. Anything else is spelled out on the command line as
# experiment numbers, so a one-off question does not need a constant.
NAMED_SCOPES = {
    "primary": PRIMARY_SCOPE,
    "free-bnoh": FREE_BNOH,
    "free-bnoh-all": FREE_BNOH_ALL,
    "free-bnoh-neutral": FREE_BNOH_NEUTRAL,
    "free-bnoh-phosphate": FREE_BNOH_PHOSPHATE,
    "boric": BORIC_BUFFER,
    "paired": tuple(sorted({e for free, cat, _ in PAIRED_CONTROLS
                            for e in (*free, cat)})),
}


def parse_scope(text):
    """
    A scope from a name, or from experiment numbers: "3,6" or "135-151".

    Returns a frozenset. Raises ValueError on anything it cannot read, rather
    than quietly returning an empty scope -- an empty scope produces an empty
    frame, and an empty frame is a table of zeroes that looks like a result.
    """
    if text in NAMED_SCOPES:
        return frozenset(NAMED_SCOPES[text])
    experiments = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part.lstrip("-"):
            low, high = part.split("-", 1)
            experiments.update(range(int(low), int(high) + 1))
        else:
            experiments.add(int(part))
    if not experiments:
        raise ValueError(f"empty scope: {text!r}")
    return frozenset(experiments)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--scope", default="primary",
                        help="a named scope (%s) or experiment numbers "
                             "(\"3,6\", \"135-151\"); default primary"
                             % ", ".join(NAMED_SCOPES))
    parser.add_argument("--design", action="store_true",
                        help="print the per-experiment design table")
    parser.add_argument("--orders", action="store_true",
                        help="print the apparent reaction orders")
    parser.add_argument("--controls", action="store_true",
                        help="print the +/- chemzyme paired controls")
    parser.add_argument("--literature", action="store_true",
                        help="compare the paired controls against the literature")
    parser.add_argument("--buffer", action="store_true",
                        help="the enzyme-free rate's dependence on [buf]")
    arguments = parser.parse_args()
    chosen = parse_scope(arguments.scope)

    if arguments.design:
        table = design(chosen)
        with pd.option_context("display.width", 200, "display.max_columns", 20):
            print(table.to_string(float_format=lambda v: f"{v:.3g}"))
        return 0

    if arguments.orders:
        with pd.option_context("display.width", 200):
            print(order_table(chosen).to_string(
                float_format=lambda v: f"{v:.3f}"))
        return 0

    if arguments.controls:
        with pd.option_context("display.width", 200):
            print(paired_controls().to_string(
                index=False, float_format=lambda v: f"{v:.4g}"))
        effect = catalytic_effect()
        print(f"\nvmax(+chemzyme)/vmax(-chemzyme) over "
              f"{effect['rungs']} live matched rungs: "
              f"median {effect['median_ratio']:.2f}x, range "
              f"{effect['ratio_range'][0]:.2f}-{effect['ratio_range'][1]:.2f}x")
        print(f"same experiment run twice (exps 69 vs 70) disagrees by up to "
              f"{effect['replicate_scatter']:.2f}x -- the ratio is not "
              f"resolved from no effect")
        return 0

    if arguments.buffer:
        print("enzyme-free BnOH: the substrate order depends on what the "
              "BUFFER was doing\n")
        result = buffer_dependence()
        print(f"  [buf] held constant   exps {BUFFER_FIXED}, "
              f"n={result['n_fixed']:2d}:  order in [sub] "
              f"{result['order_s0_fixed']:+.3f} +/- {result['stderr_s0_fixed']:.3f}")
        print(f"  [buf] falling         exps {BUFFER_CONFOUNDED}, "
              f"n={result['n_titration']:2d}:  order in [sub] "
              f"{result['order_s0_titration']:+.3f} +/- "
              f"{result['stderr_s0_titration']:.3f}")
        print(f"\n  The same reaction reads a POSITIVE substrate order where "
              f"[buf] is held and a\n  NEGATIVE one where [buf] falls as "
              f"[sub] rises. The difference is the buffer.\n")
        print(f"  coupling  dlog[buf]/dlog[sub] within the titrations: "
              f"{result['coupling']:+.3f}")
        print(f"  => order in [buf] = {result['order_buf']:+.2f} +/- "
              f"{result['stderr_buf']:.2f}   (approximately FIRST order)\n")
        check = buffer_cross_check()
        print(f"  independent cross-check, 4OMe-BnOH / 40 C / phosphate, "
              f"n={check['n']} (between-run\n  contrast at fixed pH and "
              f"[H2O2], different substrate and temperature):")
        print(f"      order in [buf]  {check['order_buf']:+.2f} +/- "
              f"{check['stderr_buf']:.2f}   (VIF {check['vif_buf']:.1f})")
        print(f"      order in [sub]  {check['order_s0']:+.2f} +/- "
              f"{check['stderr_s0']:.2f}")
        with_boric = buffer_dependence(anchor=FREE_BNOH)
        print(f"\n  The anchor is exps {BUFFER_FIXED} -- PHOSPHATE ONLY. Exp 65 "
              f"left it on 2026-09-01:\n  its curves break mid-run and have no "
              f"usable rate (BORIC_RATE_UNUSABLE,\n  `synchronised_break`). "
              f"With it in, this read {with_boric['order_buf']:+.2f} +/- "
              f"{with_boric['stderr_buf']:.2f}.")

        print(f"\n  WHAT MAKES IT FIRST ORDER is NOT SETTLED. Catalysis by a "
              f"buffer species, or the\n  buffer making an oxidant? Both are "
              f"first order in a buffer species and the\n  phosphate design "
              f"cannot separate them. Borate is where the second is not a\n  "
              f"hypothesis, and exp {PEROXO_PAIR[0]} is the only boric "
              f"BACKGROUND run in the archive -- so the\n  test rests on the "
              f"one run whose curves are unusable. Predicting it from a law\n  "
              f"fitted on phosphate alone ({PEROXO_PAIR[0]} against "
              f"{PEROXO_PAIR[1]}, cuvette for cuvette):\n")
        peroxo = peroxo_buffer_test()
        for estimator, row in peroxo.iterrows():
            print(f"      {estimator:9s} predicted {row['predicted']:.2f}x   "
                  f"observed {row['observed']:.2f}x   "
                  f"excess {row['excess']:.2f}x   (n={int(row['n'])})")
        print(f"\n  Read as 'borate is nowhere fast', not as 'borate matches "
              f"the law': exp 65's\n  vmax reads the post-break stretch and "
              f"its v0 the pre-break one. The catalysed\n  runs cannot stand "
              f"in -- their reference channel already subtracted the\n  "
              f"background, so a buffer-made oxidant cancels. See "
              f"background_reaction/\n  ANALYSIS.md section 6b; the decisive "
              f"tests are 31P NMR and a repeat of exp 65.")

        print(f"\n  The in-scope block is UNAFFECTED: [buf] = 75.013 mM in "
              f"all 119 of its curves,\n  so no buffer variation can reach "
              f"its substrate order.")
        return 0

    if arguments.literature:
        orders = background_orders()
        print(f"enzyme-free BnOH background, {orders['n']} live curves "
              f"(R2 {orders['r2']:.2f}):")
        for name in ("s0", "h2o2", "hoo"):
            print(f"  order in {name:5s} {orders['order_' + name]:+.2f} "
                  f"+/- {orders['stderr_' + name]:.2f}")
        print(f"\n  the [HOO-] order is why pH matters: the background climbs "
              f"about tenfold per pH unit,\n  so a background measured at pH 8 "
              f"says little about one at pH 7.\n")

        table = literature_comparison()
        print(f"at the literature's pH {LITERATURE['pH']:.1f} / "
              f"{LITERATURE['h2o2_mM']:.0f} mM H2O2, our near-neutral "
              f"enzyme-free runs {FREE_BNOH_NEUTRAL} against its uncatalysed rate:")
        with pd.option_context("display.width", 200):
            print(table.to_string(index=False,
                                  float_format=lambda v: f"{v:.4g}"))
        print(f"\n  median excess: {table.excess.median():.0f}x  "
              f"(range {table.excess.min():.0f}-{table.excess.max():.0f}x)")

        enhancement = predicted_enhancement()
        print(f"\nwhat the literature's kcat would show at this archive's "
              f"loading:")
        print(f"  predicted enhancement over background: median "
              f"{enhancement.expected_ratio.median():.2f}x "
              f"(range {enhancement.expected_ratio.min():.2f}-"
              f"{enhancement.expected_ratio.max():.2f}x)")
        print(f"  observed:                              median "
              f"{(enhancement.observed_mM_s / enhancement.background_mM_s).median():.2f}x")
        print(f"  at the literature's {LITERATURE['catalyst_mM']} mM loading:   "
              f"    median "
              f"{enhancement.at_literature_loading.median():.0f}x")
        print(f"\n  exps 69 and 70 are the same experiment run twice and "
              f"disagree by up to "
              f"{catalytic_effect()['replicate_scatter']:.2f}x, so the "
              f"predicted\n  enhancement is smaller than the reproducibility. "
              f"These runs could not have seen it.")
        print(f"\n  source: {LITERATURE['source']}")
        return 0

    facts = summary(chosen)
    cells = blocks(chosen)
    span = (f"exps {min(chosen)}-{max(chosen)}"
            if sorted(chosen) == list(range(min(chosen), max(chosen) + 1))
            else "exps " + ",".join(str(e) for e in sorted(chosen)))
    print(f"scope        {span}")
    for (substrate, temperature, buffer_name), row in cells.iterrows():
        print(f"block        {substrate} / {temperature:.0f} C / {buffer_name}"
              f"  ({row.curves} curves, {row.experiments} experiments)")
    print(f"curves       {facts['curves']} over {facts['experiments']} experiments, "
          f"{facts['live_curves']} with a live signal")
    print(f"pH           {facts['pH_range'][0]:.2f} to {facts['pH_range'][1]:.2f}, "
          f"{facts['hoo_decades']:.1f} decades of [HOO-]")
    print(f"contrast     log[S] {100 * facts['within_experiment_s0']:.1f}% "
          f"within-experiment, log[H2O2] "
          f"{100 * facts['within_experiment_h2o2']:.1f}%")
    print(f"lag          {facts['lagging']} curves peak after "
          f"{LAG_THRESHOLD:.0%} of the run")
    print(f"acceleration {facts['accelerating_live']} of "
          f"{facts['live_curves']} live curves are steeper later than at "
          f"the start, by >{ACCELERATION_SIGMA:.0f} sigma")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
