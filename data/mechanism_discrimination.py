"""
Which mechanisms survive the progress curves -- Stage C of the mechanism work.

The plan is `PLAN_MECHANISM_DISCRIMINATION.md`. This module holds it: the three
SUMMARIES read off every live catalysed curve (L, E, D), the five candidate
mechanisms, the random-effect likelihood, the fitter and the cross-validation.

The summaries (plan section 3). A summary is a shape number read off a curve
without a fit, so it says the same thing about the readings and about a
candidate's predicted curve. `L` is the level, `ln` of the slope over the Q3
window (0.50-0.75 of the run); `E` is the early shape, `ln(slope_H1/slope_H2)`
over the first and second halves; `D` is the late shape,
`ln(slope_Q4/slope_Q3)`. Every slope is an ordinary least-squares slope of the
values on time inside its window, with error `noise/sqrt(sum(tc*tc))`, and a
summary is admitted only where the slopes it takes the log of clear twice
their own error -- otherwise it is NaN and is never scored. A candidate's
predicted summaries use the same windows and the same slope formula.

The candidates (plan section 4). C0-C4 are five hand-written mechanisms for
one catalysed curve, sharing an activation clock, one turnover constant per
buffer, saturation in [S], a rate proportional to [enz] and one initial active
fraction per substrate. They differ in what draws the catalyst active (`X`:
nothing but base, the buffer, or HOO-), in what the rate saturates in (`z`),
and in whether the product is lost (C0-C3) or the catalyst is (C4). Fitting
scores the SUMMARIES through a normal random effect per run (plan section
4.6), never a fitted curve parameter -- those scatter twice as much between
repeat runs as the summaries do.

A number that reaches a document is returned by a function here.
"""
import os

import numpy as np
import pandas as pd

import scope


SUMMARY_WINDOWS = (
    ("Q3", 0.50, 0.75),
    ("H1", 0.00, 0.50),
    ("H2", 0.50, 1.00),
    ("Q4", 0.75, 1.00),
)


def window_slopes(times, values, noise):
    """
    `{name: (slope, se)}` for the four windows (plan section 3).

    A window `(a, b)` holds the readings with
    `t[0] + a*T <= t <= t[0] + b*T`, both ends inclusive, where
    `T = t[-1] - t[0]`. The slope is the ordinary least-squares slope of
    `values` on `t` inside it (`slope = sum(tc*(y - mean y)) / sum(tc*tc)`),
    and the error is `noise/sqrt(sum(tc*tc))`. An empty or single-reading
    window returns `(nan, inf)`.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    span = float(times[-1] - times[0])
    out = {}
    for name, a, b in SUMMARY_WINDOWS:
        low = times[0] + a * span
        high = times[0] + b * span
        chosen = (times >= low) & (times <= high)
        t = times[chosen]
        y = values[chosen]
        if t.size < 2:
            out[name] = (np.nan, np.inf)
            continue
        tc = t - t.mean()
        sxx = float(np.sum(tc * tc))
        if sxx <= 0.0:
            out[name] = (np.nan, np.inf)
            continue
        out[name] = (float(np.sum(tc * (y - y.mean())) / sxx),
                     float(noise) / np.sqrt(sxx))
    return out


def curve_summaries(times, values, noise):
    """
    `{"L", "E", "D", "L_se", "E_se", "D_se"}` for one curve (plan section 3).

    A summary that is not admitted is NaN and is never scored; its `se` is
    computed for every curve, admitted or not.
    """
    slopes = window_slopes(times, values, noise)
    q3, q3_se = slopes["Q3"]
    h1, h1_se = slopes["H1"]
    h2, h2_se = slopes["H2"]
    q4, q4_se = slopes["Q4"]
    out = {"L": np.nan, "E": np.nan, "D": np.nan,
           "L_se": np.nan, "E_se": np.nan, "D_se": np.nan}
    if q3 > 2.0 * q3_se:
        out["L"] = float(np.log(q3))
        out["L_se"] = q3_se / abs(q3)
    if h1 > 2.0 * h1_se and h2 > 2.0 * h2_se:
        out["E"] = float(np.log(h1 / h2))
        out["E_se"] = float(np.hypot(h1_se / abs(h1), h2_se / abs(h2)))
    if q4 > 2.0 * q4_se and q3 > 2.0 * q3_se:
        out["D"] = float(np.log(q4 / q3))
        out["D_se"] = float(np.hypot(q4_se / abs(q4), q3_se / abs(q3)))
    return out


def summary_table(substrate):
    """
    Section 3's curves for one substrate, with the six summary columns.

    One row per curve of `scope.frame(scope.archive())` with `live` true,
    `e0 > 0` and this substrate, in frame order with the index reset and no
    other exclusion -- curves the activation-sink form cannot hold ARE
    included. `times` holds each row's own time array. The summaries are read
    from `CurveFit.corrected` and `CurveFit.noise`.
    """
    frame = scope.frame(scope.archive())
    rows = frame[frame.live & (frame.e0 > 0)
                 & (frame.substrate == substrate)].reset_index(drop=True)
    fits = scope.fits(scope.archive())
    times = []
    summaries = []
    for record in rows.itertuples(index=False):
        fit = fits[(int(record.experiment), int(record.sample))]
        t = np.asarray(fit.times, dtype=float)
        times.append(t)
        summaries.append(
            curve_summaries(t, np.asarray(fit.corrected, dtype=float),
                            float(fit.noise)))
    keep = ["experiment", "sample", "buffer", "temperature", "pH", "s0",
            "h2o2", "hoo", "buf", "e0"]
    out_rows = rows[keep].copy()
    for name in ("L", "E", "D", "L_se", "E_se", "D_se"):
        out_rows[name] = [summary[name] for summary in summaries]
    return {"rows": out_rows, "times": times}


def _pooled_sd(values, admitted, groups, keep=None):
    """
    Pooled within-group SD over the admitted values.

    `groups` labels each row (a run, or a cuvette for the replicate pooling);
    a group contributes only with two or more admitted values, and `keep`
    optionally restricts to a subset of rows. Returns `(sd, df, n)`, with
    `(nan, 0, 0)` when nothing pools. The denominator is `df`.
    """
    if keep is None:
        keep = np.ones(values.shape, dtype=bool)
    use = admitted & keep
    total = 0.0
    df = 0
    n = 0
    for group in pd.unique(groups[use]):
        vals = values[use & (groups == group)]
        if vals.size < 2:
            continue
        total += float(np.sum((vals - vals.mean()) ** 2))
        df += vals.size - 1
        n += vals.size
    if df == 0:
        return float("nan"), 0, 0
    return float(np.sqrt(total / df)), df, n


def summary_scatter(substrate):
    """
    The scatter of each summary, within runs and across replicates.

    Indexed by `L`, `E`, `D` (plan section 7). `within_run_sd` pools runs with
    two or more admitted curves, `sqrt(sum((value - run mean)^2) /
    sum(n_run - 1))`, with `within_run_df` the denominator; it still holds
    every design effect, so it is a RAW spread rather than a noise estimate.
    `replicate_*` pools the same way over `scope.REPLICATE_RUNS`, grouped by
    cuvette, and is NaN where there are none (BnOH).
    """
    table = summary_table(substrate)
    rows = table["rows"]
    experiments = rows.experiment.to_numpy()
    cuvettes = rows["sample"].to_numpy()
    replicate_keep = np.isin(experiments, scope.REPLICATE_RUNS)
    out = []
    for name in ("L", "E", "D"):
        values = rows[name].to_numpy(dtype=float)
        admitted = np.isfinite(values)
        se = rows[f"{name}_se"].to_numpy(dtype=float)
        within_sd, within_df, _ = _pooled_sd(values, admitted, experiments)
        rep_sd, rep_df, rep_n = _pooled_sd(values, admitted, cuvettes,
                                           keep=replicate_keep)
        out.append({
            "admitted": int(admitted.sum()),
            "curves": int(values.size),
            "q10": float(np.nanquantile(values, 0.10)) if admitted.any()
            else np.nan,
            "q50": float(np.nanquantile(values, 0.50)) if admitted.any()
            else np.nan,
            "q90": float(np.nanquantile(values, 0.90)) if admitted.any()
            else np.nan,
            "median_se": float(np.median(se[admitted])) if admitted.any()
            else np.nan,
            "within_run_sd": within_sd,
            "within_run_df": int(within_df),
            "replicate_sd": rep_sd,
            "replicate_df": int(rep_df),
            "replicate_n": int(rep_n),
        })
    return pd.DataFrame(out, index=pd.Index(("L", "E", "D"), name="summary"))
