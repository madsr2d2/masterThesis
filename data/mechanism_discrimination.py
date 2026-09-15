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
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import qmc

import scope
import summary_kinetics


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
    return {"rows": out_rows, "times": times, "substrate": substrate}


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


CANDIDATES = ("C0", "C1", "C2", "C3", "C4")

# (X, z, loss) per candidate (plan section 4.3). `X` names what draws the
# catalyst active, `z` what the turnover rate saturates in, `loss` whether the
# product (C0-C3) or the catalyst (C4) is lost.
_LAW_PARTS = {
    "C0": ("one", "hoo", "product"),
    "C1": ("buf", "hoo", "product"),
    "C2": ("buf", "buf_hoo", "product"),
    "C3": ("hoo", "hoo", "product"),
    "C4": ("buf", "hoo", "catalyst"),
}

_GAS_CONSTANT = 8.314462618e-3      # kJ / (mol K)
_REFERENCE_TEMPERATURE = 298.15     # K


def candidate_names(candidate, substrate, buffers):
    """
    Section 4.5's parameter vector, in exactly that order.

    Buffers are sorted alphabetically, so every fit and every fold share one
    layout. 4OMe-BnOH carries the three Arrhenius parameters; BnOH carries
    none.
    """
    names = [f"lk_cat[{buffer}]" for buffer in sorted(buffers)]
    names += ["lK_S", "lK_O", "lk_f", "lk_r", "pKa", "theta0"]
    names.append("lk_d" if candidate == "C4" else "lk_s")
    if substrate == "4OMe-BnOH":
        names += ["Ea_cat", "Ea_act",
                  "Ea_d" if candidate == "C4" else "Ea_s"]
    names += ["ls_w_L", "ls_w_E", "ls_w_D", "ls_b_L", "ls_b_E", "ls_b_D"]
    return names


def candidate_bounds(candidate, substrate, buffers):
    """
    Section 4.5's L-BFGS-B box bounds, in `candidate_names` order.
    """
    lower, upper = [], []

    def add(lo, hi):
        lower.append(lo)
        upper.append(hi)

    for _ in sorted(buffers):
        add(-8.0, 2.0)                       # lk_cat
    add(-2.0, 4.0)                           # lK_S
    add(-8.0, 4.0 if candidate == "C2" else 2.0)   # lK_O
    add(-12.0, 6.0)                          # lk_f
    add(-10.0, -1.0)                         # lk_r
    add(4.0, 14.0)                           # pKa
    add(0.0, 1.0)                            # theta0
    add(-9.0, -2.0)                          # lk_s / lk_d
    if substrate == "4OMe-BnOH":
        add(0.0, 200.0)                      # Ea_cat
        add(0.0, 200.0)                      # Ea_act
        add(0.0, 200.0)                      # Ea_s / Ea_d
    for _ in range(3):
        add(-3.0, 1.0)                       # ls_w_L/E/D
    for _ in range(3):
        add(-3.0, 1.0)                       # ls_b_L/E/D
    return np.array(lower, dtype=float), np.array(upper, dtype=float)


def candidate_curves(candidate, parameters, rows, times):
    """
    Every curve's `A(t)` and the mechanism quantities behind it (4.2, 4.4).

    `parameters` is a `{name: value}` dict; `rows` and `times` come from a
    `summary_table`. Returns `{"A": [array per curve], "tau", "theta_ss", "V",
    "k"}`, each of the last four an array over curves.
    """
    buffer = rows["buffer"].to_numpy()
    temperature = rows["temperature"].to_numpy(dtype=float)
    pH = rows["pH"].to_numpy(dtype=float)
    s0 = rows["s0"].to_numpy(dtype=float)
    hoo = rows["hoo"].to_numpy(dtype=float)
    buf = rows["buf"].to_numpy(dtype=float)
    e0 = rows["e0"].to_numpy(dtype=float)
    invT = 1.0 / (temperature + 273.15) - 1.0 / _REFERENCE_TEMPERATURE
    base = 1.0 / (1.0 + 10.0 ** (pH - parameters["pKa"]))
    with np.errstate(all="ignore"):
        if "Ea_act" in parameters:
            activation = np.exp(-parameters["Ea_act"] * invT / _GAS_CONSTANT)
            catalyst = np.exp(-parameters["Ea_cat"] * invT / _GAS_CONSTANT)
            loss = np.exp(
                -(parameters["Ea_d"] if candidate == "C4"
                  else parameters["Ea_s"]) * invT / _GAS_CONSTANT)
        else:
            activation = 1.0
            catalyst = 1.0
            loss = 1.0
        draw, saturate, _ = _LAW_PARTS[candidate]
        X = {"one": 1.0, "buf": buf, "hoo": hoo}[draw]
        z = {"buf_hoo": buf * hoo,
             "hoo": hoo}[saturate]
        q = 10.0 ** parameters["lk_f"] * X * base * activation
        k_r = 10.0 ** parameters["lk_r"] * activation
        tau = 1.0 / (q + k_r)
        theta_ss = q / (q + k_r)
        Y = z / (10.0 ** parameters["lK_O"] + z)
        lk_cat = np.array([parameters[f"lk_cat[{b}]"] for b in buffer])
        V = (10.0 ** lk_cat * e0 * s0 / (10.0 ** parameters["lK_S"] + s0)
             * Y * catalyst)
        k = np.broadcast_to(
            np.asarray(10.0 ** parameters[
                "lk_d" if candidate == "C4" else "lk_s"] * loss, dtype=float),
            (len(rows),))
    A = []
    for index in range(len(rows)):
        t = np.asarray(times[index], dtype=float)
        if candidate == "C4":
            r = k[index] + 1.0 / tau[index]
            curve = V[index] * (
                theta_ss[index] * (-np.expm1(-k[index] * t)) / k[index]
                + (parameters["theta0"] - theta_ss[index])
                * (-np.expm1(-r * t)) / r)
        else:
            h, g = summary_kinetics._activation_sink_columns(
                t, tau[index], k[index])
            curve = (V[index] * parameters["theta0"] * h[0]
                     + V[index] * theta_ss[index] * g[0])
        A.append(curve)
    return {"A": A, "tau": tau, "theta_ss": theta_ss, "V": V, "k": k}


def candidate_predict(candidate, parameters, rows, times):
    """
    The predicted `{"L", "E", "D"}` arrays (plan section 3).

    The same windows and slope formula as the observed summaries, on the
    candidate's `A(t)` at each curve's own times. A slope a summary takes the
    log of that is not positive gives `-50.0`.
    """
    curves = candidate_curves(candidate, parameters, rows, times)
    out = {"L": [], "E": [], "D": []}
    for t, curve in zip(times, curves["A"]):
        slopes = window_slopes(t, curve, 1.0)
        q3 = slopes["Q3"][0]
        h1 = slopes["H1"][0]
        h2 = slopes["H2"][0]
        q4 = slopes["Q4"][0]
        out["L"].append(np.log(q3) if q3 > 0 else -50.0)
        out["E"].append(np.log(h1 / h2) if h1 > 0 and h2 > 0 else -50.0)
        out["D"].append(np.log(q4 / q3) if q4 > 0 and q3 > 0 else -50.0)
    return {name: np.array(values, dtype=float)
            for name, values in out.items()}


def run_likelihood(errors, se, sigma_w, sigma_b, runs):
    """
    `NLL_g` per run, a `pandas.Series` indexed by run (plan section 4.6).

    The exact negative log-density of a normal with a shared random offset per
    run, over the curves being scored.
    """
    errors = np.asarray(errors, dtype=float)
    se = np.asarray(se, dtype=float)
    runs = np.asarray(runs)
    out = {}
    for run in pd.unique(runs):
        mask = runs == run
        e = errors[mask]
        d = se[mask] ** 2 + sigma_w ** 2
        A = float(np.sum(1.0 / d))
        B = float(np.sum(e / d))
        C = float(np.sum(e ** 2 / d))
        out[run] = 0.5 * (
            C - sigma_b ** 2 * B ** 2 / (1.0 + sigma_b ** 2 * A)
            + float(np.sum(np.log(d))) + float(np.log(1.0 + sigma_b ** 2 * A))
            + e.size * float(np.log(2.0 * np.pi)))
    return pd.Series(out)


def candidate_nll(candidate, x, table, keep=None):
    """
    The total negative log-likelihood of a candidate's summaries (section 4.6).

    Scores the rows admitted AND kept, for each of L, E, D, with
    `sigma_w = 10**ls_w_q` and `sigma_b = 10**ls_b_q`. `keep` is a boolean mask
    over the table's rows (None = all). Returns `1e12` if the total is not
    finite. The buffers are always the full table's sorted buffers, so the
    vector layout never changes between folds.
    """
    rows = table["rows"]
    times = table["times"]
    substrate = table["substrate"]
    buffers = sorted(pd.unique(rows["buffer"]))
    names = candidate_names(candidate, substrate, buffers)
    parameters = dict(zip(names, np.asarray(x, dtype=float)))
    predicted = candidate_predict(candidate, parameters, rows, times)
    experiments = rows["experiment"].to_numpy()
    if keep is None:
        keep = np.ones(len(rows), dtype=bool)
    total = 0.0
    for name in ("L", "E", "D"):
        observed = rows[name].to_numpy(dtype=float)
        se = rows[f"{name}_se"].to_numpy(dtype=float)
        scored = np.isfinite(observed) & keep
        if not scored.any():
            continue
        sigma_w = 10.0 ** float(parameters[f"ls_w_{name}"])
        sigma_b = 10.0 ** float(parameters[f"ls_b_{name}"])
        terms = run_likelihood(observed[scored] - predicted[name][scored],
                               se[scored], sigma_w, sigma_b,
                               experiments[scored])
        total += float(terms.sum())
    return total if np.isfinite(total) else 1e12


def candidate_degeneracy(candidate, fit, table, substrate):
    """
    The degeneracy flags, in the plan's order (section 8.9).

    Returns `(verdict, flags)`: `"not degenerate"` or
    `"degenerate: " + "; ".join(flags)`, plus the flag list. "clocks outside
    the run window" counts a curve whose `tau` is above ten run lengths or
    below a three-hundredth of one.
    """
    rows = table["rows"]
    times = table["times"]
    buffers = sorted(pd.unique(rows["buffer"]))
    names = candidate_names(candidate, substrate, buffers)
    lower, upper = candidate_bounds(candidate, substrate, buffers)
    x = np.asarray(fit["x"], dtype=float)
    flags = []
    if not fit["success"]:
        flags.append("not converged")
    span = upper - lower
    for name, value, lo, hi in zip(names, x, lower, upper):
        if min(value - lo, hi - value) <= 0.01 * (hi - lo):
            flags.append(f"at bound: {name}")
    parameters = dict(zip(names, x))
    tau = candidate_curves(candidate, parameters, rows, times)["tau"]
    outside = 0
    for index, t in enumerate(times):
        duration = float(np.asarray(t, dtype=float)[-1] - np.asarray(t, dtype=float)[0])
        if tau[index] > 10.0 * duration or tau[index] < duration / 300.0:
            outside += 1
    if outside > len(rows) / 2:
        flags.append(
            f"clocks outside the run window on {outside} of {len(rows)} curves")
    if flags:
        return "degenerate: " + "; ".join(flags), flags
    return "not degenerate", flags


def fit_candidate(candidate, table, substrate, keep=None, start=None,
                  starts=16, seed=0, maxiter=1000):
    """
    Fit one candidate (plan section 8.8), from one start or many.

    With `start`, one L-BFGS-B run from it; otherwise one from each row of a
    Latin-hypercube design. Returns the lowest-NLL start's `names`, `x`, `nll`,
    `success`, `nit`, the per-start `(nll, success)` list, `refit`, and the
    degeneracy verdict in `flags` (with the flag list in `flag_list`). If the
    best start did not succeed it is refit ONCE from its own `x`.
    """
    rows = table["rows"]
    buffers = sorted(pd.unique(rows["buffer"]))
    names = candidate_names(candidate, substrate, buffers)
    lower, upper = candidate_bounds(candidate, substrate, buffers)

    def objective(vector):
        return candidate_nll(candidate, vector, table, keep)

    if start is not None:
        starts_x = [np.asarray(start, dtype=float)]
    else:
        unit = qmc.LatinHypercube(d=len(names), seed=seed).random(starts)
        starts_x = [lower + (upper - lower) * row for row in unit]
    results = []
    for x0 in starts_x:
        result = minimize(objective, x0, method="L-BFGS-B",
                          bounds=list(zip(lower, upper)),
                          options={"maxiter": maxiter})
        results.append(result)
    best = min(results, key=lambda result: float(result.fun))
    refit = False
    if not bool(best.success):
        retry = minimize(objective, best.x, method="L-BFGS-B",
                         bounds=list(zip(lower, upper)),
                         options={"maxiter": 3000})
        best = retry
        refit = True
    verdict, flags = candidate_degeneracy(
        candidate, {"success": bool(best.success), "x": best.x}, table,
        substrate)
    return {
        "names": names,
        "x": best.x,
        "nll": float(best.fun),
        "success": bool(best.success),
        "nit": int(best.nit),
        "starts": [(float(result.fun), bool(result.success))
                   for result in results],
        "refit": refit,
        "flags": verdict,
        "flag_list": flags,
    }


def cross_validate_candidate(candidate, table, substrate, full):
    """
    Score each held-out run under a fit to the others (plan section 8.10).

    A run whose buffer appears in no other run is skipped and counted, because
    the fold would have no `lk_cat` for it. Otherwise the fold is fitted with
    `keep = experiment != run`, starting at the full fit's `x`, and the
    held-out run is scored -- never refitted.
    """
    rows = table["rows"]
    experiments = rows["experiment"].to_numpy()
    buffers = rows["buffer"].to_numpy()
    folds = {}
    fold_success = {}
    skipped = 0
    for run in sorted(pd.unique(experiments)):
        held = set(buffers[experiments == run])
        others = set(buffers[experiments != run])
        if held - others:
            skipped += 1
            continue
        fold = fit_candidate(candidate, table, substrate,
                             keep=experiments != run, start=full["x"],
                             maxiter=1000)
        folds[int(run)] = float(candidate_nll(
            candidate, fold["x"], table, keep=experiments == run))
        fold_success[int(run)] = bool(fold["success"])
    return {"folds": folds, "skipped": skipped,
            "fold_success": fold_success}


def candidate_tie(scores):
    """
    The tie rule (plan section 8.11), one row per candidate.

    `scores` is candidates x folds; a fold with any NaN is dropped. `best` is
    the lowest total; a candidate is `"tied"` if its mean difference from the
    best is within `2*std(ddof=1)/sqrt(folds)`, else `"excluded"`.
    `worst_fold_difference` shows a candidate that fails on a few folds only.
    """
    scores = scores.dropna(axis=1)
    totals = scores.sum(axis=1)
    best = totals.idxmin()
    difference = scores.subtract(scores.loc[best], axis=1)
    folds = scores.shape[1]
    out = []
    for candidate in scores.index:
        mean = float(difference.loc[candidate].mean())
        bar = (2.0 * float(difference.loc[candidate].std(ddof=1))
               / np.sqrt(folds))
        if candidate == best:
            status = "best"
        elif mean <= bar:
            status = "tied"
        else:
            status = "excluded"
        out.append({"total": float(totals.loc[candidate]),
                    "mean_difference": mean, "bar": bar, "status": status,
                    "worst_fold_difference":
                        float(difference.loc[candidate].max())})
    return pd.DataFrame(out, index=scores.index)
