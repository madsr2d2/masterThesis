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
                           acceleration, initial_rate, lag_time,
                           peak_position, peak_rate)
from fit_dataset import (PRIMARY_SCOPE, PRIMARY_SCOPE_BLOCK, build_curves,
                         in_scope)

# A run's own cuvettes have to move an axis by at least this much before that
# axis counts as measured inside the run rather than across experiments.
LADDER_MINIMUM = 2.0

# Net change below this multiple of a curve's own noise is not a measurement of
# anything. Exps 150 and 151 are mostly flat by this rule, which is the point:
# they are the block's in-cell background.
LIVE_SIGNAL_NOISE_MULTIPLE = 20.0


def curves(scope=PRIMARY_SCOPE):
    """The in-scope Curve objects, in (experiment, sample) order."""
    all_curves, _ = build_curves()
    return sorted(in_scope(all_curves, scope),
                  key=lambda c: (c.experiment, c.sample))


def frame(scope=PRIMARY_SCOPE):
    """
    One row per in-scope curve, with every derived quantity already attached.

    Columns: experiment, source, sample, pH, s0, h2o2, e0, hoo, duration_s,
    points,
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
        net = float(values[-1] - values[0])
        v0, v0_stderr, v0_rms = initial_rate(times, values)
        peak = peak_position(values, times)
        accel_z, accel_where = acceleration(times, values)
        vmax, vmax_stderr, vmax_where = peak_rate(times, values)
        rows.append({
            "experiment": curve.experiment,
            "source": curve.source,
            "sample": curve.sample,
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
            "lag_time_s": lag_time(times, values),
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


def summary(scope=PRIMARY_SCOPE):
    """The scope in one paragraph of numbers, for printing."""
    data = frame(scope)
    return {
        "block": PRIMARY_SCOPE_BLOCK,
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


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--design", action="store_true",
                        help="print the per-experiment design table")
    parser.add_argument("--orders", action="store_true",
                        help="print the apparent reaction orders")
    arguments = parser.parse_args()

    if arguments.design:
        table = design()
        with pd.option_context("display.width", 200, "display.max_columns", 20):
            print(table.to_string(float_format=lambda v: f"{v:.3g}"))
        return 0

    if arguments.orders:
        with pd.option_context("display.width", 200):
            print(order_table().to_string(float_format=lambda v: f"{v:.3f}"))
        return 0

    facts = summary()
    print(f"scope        exps {min(PRIMARY_SCOPE)}-{max(PRIMARY_SCOPE)}, "
          f"block {facts['block'][0]} / {facts['block'][1]:.0f} C / {facts['block'][2]}")
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
