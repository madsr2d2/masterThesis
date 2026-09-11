"""
Which rate should an order be measured on? Every headline result, on every
candidate, side by side.

    python data/rate_choice.py

WHY THIS EXISTS. `vmax` -- `curve_metrics.peak_rate`, the steepest 20% block
slope of the READINGS -- is the rate almost every headline order in this
package is measured on, and it is not a parameter of the function the curves
pages now draw. That function,

    A(t) = c + v_ss.t - sum B_i (1 - e^(-t/tau_i)),

defines three rates of its own, and they sit at different stages of the
mechanism: `v0_fit_corrected` is the rate at mixing (on a lag, before the
catalyst has activated; on a burst, the fast first turnover),
`v_peak_corrected` is the fitted maximum (full activation, before the
extent-dependent decline), and `v_ss_fit_corrected` is the rate the fit
settles onto (the activated steady rate on a one-phase lag; the residual
after the decline on a two-phase curve). `lag_peak` is the fitted maximum
read INSIDE the run, which is `v_peak_corrected` wherever that maximum is
reached before the last reading. All four come off the GAS-CORRECTED
`fit_progress`, because the gas is made from peroxide and would inflate the
peroxide order of any rate read with it still in.

Replacing `vmax` is a scientific choice and it is not made here. This module
is the evidence the choice is made from: it routes each candidate through the
SAME functions that produced the published numbers -- nothing is re-fitted or
re-implemented -- so the `vmax` column reproduces what the documents say and
every other column is the same analysis on a different rate.

TWO WAYS TO READ IT, and both are needed. `order_comparison()` lets each
candidate keep whatever curves it can -- which is what a replacement would
actually do -- but `scope.orders` drops non-positive rates, and the fitted
v(0) is negative on the early-trough curves while a two-phase `v_ss` can be
too, so each column is then fitted on a slightly different set of curves.
`order_comparison(common=True)` pairs each candidate with the reference on
the curves where THAT candidate is finite and positive, so a difference
between the two is a difference between the STATISTICS and not the samples.

IT IS PAIRWISE ON PURPOSE. A single mask over all six at once was tried first
and is not neutral: between them the fitted v(0) and v_ss go non-positive on
every curve of exps 44 and 49 -- the boric ladder's two runs above pH 10,
which ARE its turnover -- so a joint mask quietly truncated that ladder below
its peak and moved even plain `vmax`'s boric order from -0.034 to +0.186.
Neither candidate does that alone: v(0) fails on two cuvettes of exp 44 and
v_ss on the other two. Pairing confines each loss to the candidate that
causes it, and `candidate_coverage`'s `runs_lost` names the whole runs a
candidate cannot be read on, so a reader knows which rows sit on a shorter
ladder.

`shape_divergence` asks the remaining question: where the candidates part
company, is it the curve's SHAPE that decides -- and since shape tracks pH
(`scope.acceleration_by_ph`), would a shape-dependent choice of rate quietly
become a pH confound?

Nothing here excludes, replaces or reweights anything. No number from it is
in any document yet.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import arrhenius
import induction
import ph_role
import scope

# The candidates, and what each is. `vmax` is the published baseline and
# `vmax_corrected` the same statistic on the rebuilt series; the four fitted
# ones are all off `scope.CurveFit.progress_corrected`, the fit every curves
# page draws.
RATE_CANDIDATES = {
    "vmax": "steepest 20% block slope of the readings (published baseline)",
    "vmax_corrected": "the same, on the gas-corrected series",
    "v0_fit_corrected": "fitted rate at t = 0: v_ss - sum B_i/tau_i",
    "v_peak_corrected": "fitted rate's maximum, over all t",
    "lag_peak": "fitted rate's maximum inside the run",
    "v_ss_fit_corrected": "fitted rate as t -> infinity",
    # The activation-sink form's (`summary_kinetics.fit_activation_sink`),
    # off `scope.CurveFit.activation_sink`: the activated rate with no product
    # made, on every curve, and then only where its interval is narrow.
    "v_act_corrected": "activation-sink form: activated rate, no product",
    "v_act_where_resolved_corrected": "the same, blank where unresolved",
}

# The rate the shape comparison measures every other candidate against: the
# baseline, on the same gas-corrected series the fitted candidates come from.
DIVERGENCE_REFERENCE = "vmax_corrected"


def _common(data, columns):
    """
    `data` with `columns` blanked on any curve where one of them fails.

    A curve survives only if every one of `columns` is finite and positive on
    it. Blanking, rather than dropping, keeps the row -- so `live`, the design
    columns and anything an analysis reads besides the rate are untouched, and
    each analysis's own filter then removes the same curves for every column
    in the set. Only the named columns are touched; see the module docstring
    for why the set is a PAIR and not every candidate at once.
    """
    data = data.copy()
    # Deduplicated: the reference paired with itself is a pair of one.
    columns = list(dict.fromkeys(c for c in columns if c in data))
    values = data[columns].to_numpy(dtype=float)
    keep = np.all(np.isfinite(values) & (values > 0), axis=1)
    data.loc[~keep, columns] = np.nan
    return data


def _frames(common=None, shape=None):
    """
    The tables each analysis reads, standard or restricted.

    `common` is None for the standard tables, or a tuple of rate columns to
    hold to one set of curves -- `order_comparison` passes (candidate,
    reference). `shape` names a `SHAPE_CLASSES` entry: every rate column is
    blanked on curves outside that class, so each analysis sees only the class
    and nothing else about the frame changes.
    """
    strong = scope.strong_runs()
    ladder_runs = tuple(sorted(set().union(*scope.PH_LADDERS.values())))
    frames = {
        "block": scope.frame(),
        "strong": scope.frame(strong),
        "series": arrhenius.series_frame(),
        "ladders": None,
        "ph_ladders": scope.ph_ladders(strong),
    }
    if not common and shape is None:
        return frames

    def restrict(data):
        if shape is not None:
            data = data.copy()
            columns = [c for c in SHAPE_RATE_COLUMNS if c in data]
            data.loc[~shape_members(data, shape).to_numpy(), columns] = np.nan
        return _common(data, common) if common else data
    frames["block"] = restrict(frames["block"])
    frames["strong"] = restrict(frames["strong"])
    frames["series"] = restrict(frames["series"])
    frames["ladders"] = restrict(scope.frame(ladder_runs))
    frames["ph_ladders"] = {label: restrict(group) for label, group
                            in frames["ph_ladders"].items()}
    return frames


# Each analysis takes (rate column, frames) and returns
# (estimate, stderr, curves, extra) where `extra` is the one further number
# that analysis is read by -- named in RATE_ANALYSES beside it.
def _two_axis(term, frame_key):
    def analysis(rate, frames):
        fit = scope.orders(rate, frame=frames[frame_key], within=True)
        return (fit[f"order_{term}"], fit[f"stderr_{term}"], fit["n"],
                np.nan)
    return analysis


def _ph_order(rate, frames):
    table = scope.ph_order(rate, ladders=frames["ph_ladders"])
    if "pooled" not in table.index:
        return np.nan, np.nan, 0, np.nan
    row = table.loc["pooled"]
    return row.order, row.stderr, int(row.curves), np.nan


def _three_ladders(rate, frames):
    rows = ph_role.rate_ladders(rate, frame=frames["ladders"])
    pooled = ph_role.pooled_rate_order(rows, drop=("boric 4OMe",))
    curves = sum(int(r["n"]) for r in rows if r["ladder"] != "boric 4OMe")
    return (pooled.get("pooled", np.nan), pooled.get("stderr", np.nan),
            curves, pooled.get("chi2", np.nan))


def _boric(rate, frames):
    rows = ph_role.rate_ladders(rate, frame=frames["ladders"])
    row = next(r for r in rows if r["ladder"] == "boric 4OMe")
    return row["hoo"]["slope"], row["hoo"]["stderr"], int(row["n"]), np.nan


def _arrhenius(rate, frames):
    # A RATIO OF TWO RATES IS NOT DIVIDED BY THE ENZYME. `pooled_arrhenius`
    # turns a rate into a turnover by dividing by [enz]; a ratio is already
    # dimensionless, and dividing it adds a stray -log[enz] to the response,
    # which moves the slope wherever [enz] tracks temperature. It put the lag
    # class's ratio at 52.6 kJ/mol where early minus late is 57.9, and
    # `test_rate_choice` caught it on the identity the ratio exists for.
    fit = arrhenius.pooled_arrhenius(rate, frame=frames["series"],
                                     per_enzyme=rate not in RATIO_COLUMNS)
    return fit["activation_kJ"], fit["stderr_kJ"], fit["n"], np.nan


def _plus_one(rate, frames):
    table = induction.joint_clocks(frames["block"], rate=rate)
    row = table.loc[("tau_slow_corrected", "axis")]
    return row.order, row.stderr, int(row.curves), row.sigma


def _saturation(rate, frames):
    fit = induction.peroxide_saturation(frames["block"], parameter=rate)
    width = (fit.get("order_high", np.nan) - fit.get("order_low", np.nan)) / 2
    return (fit.get("order", np.nan), width, int(fit.get("points", 0)),
            fit.get("first_order_f", np.nan))


# (name, what `extra` is, the analysis). The names say the block and the cut,
# because "the peroxide order" is four different numbers in this package.
RATE_ANALYSES = (
    ("two-axis [S] order, all live", "", _two_axis("s0", "block")),
    ("two-axis [H2O2] order, all live", "", _two_axis("h2o2", "block")),
    ("two-axis [H2O2] order, strong runs", "", _two_axis("h2o2", "strong")),
    ("two-axis pH order, strong runs", "", _ph_order),
    ("[HOO-] order, three ladders pooled", "chi2", _three_ladders),
    ("[HOO-] order, boric ladder", "", _boric),
    ("activation energy, kJ/mol", "", _arrhenius),
    ("+1 rule through tau_slow", "sigma from +1", _plus_one),
    # Its error is HALF THE PROFILE INTERVAL, not a standard error: the
    # saturation fit profiles the order on a grid and reports the range the
    # F test admits.
    ("peroxide saturation order (+/- half profile)", "F against a = 1",
     _saturation),
)


def order_comparison(common=False, rates=tuple(RATE_CANDIDATES),
                     reference=DIVERGENCE_REFERENCE, shape=None):
    """
    Every analysis in `RATE_ANALYSES`, on every candidate rate.

    Returns a long table: one row per (analysis, rate) with `estimate`,
    `stderr`, `curves` and `extra` (what `extra` means is in `extra_is`).
    `order_comparison_wide` pivots it for reading.

    `common=False` lets each candidate keep the curves it can -- the table a
    replacement would actually produce. `common=True` reads each candidate AND
    `reference` on the curves where that candidate is usable, and adds
    `reference_estimate` / `reference_stderr`, so each row is a like-for-like
    comparison of two statistics. On a candidate positive on every live curve
    the mask is a no-op and the row equals the standard one.

    THE STRONG RUNS ARE HELD FIXED at `scope.strong_runs()` in both, so the
    comparison isolates the rate. Which runs would be strong under each
    candidate is a separate question, and `strong_runs_by_rate` asks it.

    `shape` restricts every analysis to one `SHAPE_CLASSES` class; see
    `shape_split`, which is the reading of it that means something.
    """
    shared = _frames(shape=shape)
    rows = []
    for rate in rates:
        frames = _frames((rate, reference), shape) if common else shared
        for name, extra_is, analysis in RATE_ANALYSES:
            estimate, stderr, curves, extra = analysis(rate, frames)
            row = {"analysis": name, "rate": rate,
                   "estimate": float(estimate), "stderr": float(stderr),
                   "curves": int(curves), "extra": float(extra),
                   "extra_is": extra_is}
            if common:
                ref, ref_err, _, _ = analysis(reference, frames)
                row.update(reference_estimate=float(ref),
                           reference_stderr=float(ref_err))
            rows.append(row)
    return pd.DataFrame(rows)


def order_comparison_wide(table=None):
    """`order_comparison` as `estimate +/- stderr (n)`, analyses by rates."""
    table = order_comparison() if table is None else table
    paired = "reference_estimate" in table

    def cell(r):
        text = (f"{r.estimate:+.3f} +/- {r.stderr:.3f} ({r.curves})"
                if np.isfinite(r.estimate) else f"-- ({r.curves})")
        if np.isfinite(r.extra):
            text += f" [{r.extra:.2f}]"
        if paired and np.isfinite(r.reference_estimate):
            text += f"  ref {r.reference_estimate:+.3f}"
        return text
    cells = table.assign(cell=[cell(r) for r in table.itertuples()])
    order = [name for name, _, _ in RATE_ANALYSES]
    wide = cells.pivot(index="analysis", columns="rate", values="cell")
    return wide.loc[order, [r for r in RATE_CANDIDATES if r in wide.columns]]


def candidate_coverage(block=None):
    """
    How many live curves each candidate can be read on, and its one caveat.

    Columns: `live`, `usable` (finite and positive), `non_positive`, and
    `caveat` / `caveat_curves` -- the specific way that candidate fails to be
    a measurement inside the run:

      v0_fit_corrected     the fitted form is not resolved, so its slope at
                           t = 0 is not pinned either
      v_peak_corrected     the fitted maximum lies at t = infinity or past the
                           last reading, so it is an extrapolation
      lag_peak             the same curves: there the in-run maximum is the
                           rate at the LAST reading, which grows with run
                           length and so carries the schedule
      v_ss_fit_corrected   two-phase curves, where the asymptote lies past the
                           decline rather than after the activation
      v_act_corrected      the activation-sink form's v_act is not pinned --
                           mostly curves still accelerating at their last
                           reading, where only v_act/tau is

    `block` is a scope; the default is the whole archive.
    """
    data = scope.frame(scope.archive() if block is None else block)
    live = data[data.live]
    peak_time = live.v_peak_corrected_time.to_numpy(dtype=float)
    beyond = ~np.isfinite(peak_time) | (peak_time > live.duration_s.to_numpy())
    caveats = {
        "vmax": ("", np.zeros(len(live), dtype=bool)),
        "vmax_corrected": ("", np.zeros(len(live), dtype=bool)),
        "v0_fit_corrected": ("fitted form unresolved",
                             ~live.v0_fit_resolved_corrected.to_numpy(bool)),
        "v_peak_corrected": ("maximum outside the run", beyond),
        "lag_peak": ("maximum at the last reading", beyond),
        "v_ss_fit_corrected": ("two-phase: past the decline",
                               live.phases_corrected.to_numpy() == 2),
        "v_act_corrected": ("v_act unresolved",
                            ~live.v_act_resolved_corrected.to_numpy(bool)),
        "v_act_where_resolved_corrected": ("", np.zeros(len(live), dtype=bool)),
    }
    rows = []
    for rate, description in RATE_CANDIDATES.items():
        values = live[rate].to_numpy(dtype=float)
        finite = np.isfinite(values)
        usable = finite & (values > 0)
        caveat, flagged = caveats[rate]
        # WHOLE RUNS, because a between-run axis loses a rung when every
        # cuvette of one run drops -- which is how exps 44 and 49, the top of
        # the boric ladder, vanish under the fitted v(0) and v_ss.
        by_run = pd.Series(usable, index=live.index).groupby(
            live.experiment.to_numpy()).any()
        rows.append({"rate": rate, "is": description, "live": len(live),
                     "usable": int(usable.sum()),
                     "non_positive": int((finite & (values <= 0)).sum()),
                     "caveat": caveat, "caveat_curves": int(flagged.sum()),
                     "runs_lost": [int(e) for e in by_run.index[~by_run]]})
    return pd.DataFrame(rows).set_index("rate")


def shape_divergence(block=None, reference=DIVERGENCE_REFERENCE):
    """
    Where each candidate parts from `reference`, grouped by the curve's shape.

    One row per (rate, shape): `curves`, the median and interquartile range of
    log(rate / reference), and the median pH of those curves. Shape is
    `progress_kind_corrected` -- the corrected fit's own kind, since that is
    the fit the candidates come from.

    READ THE pH COLUMN BESIDE THE RATIO. If one shape carries a large ratio
    and also sits at high pH, then choosing a rate that behaves differently on
    that shape moves a pH order, and the shape and the pH cannot be told apart
    in the result. `divergence_against_ph` puts a number on that.
    """
    data = scope.frame(scope.archive() if block is None else block)
    live = data[data.live]
    rows = []
    for rate in RATE_CANDIDATES:
        if rate in (reference, "vmax"):
            continue
        both = live[np.isfinite(live[rate]) & (live[rate] > 0)
                    & np.isfinite(live[reference]) & (live[reference] > 0)]
        ratio = np.log(both[rate] / both[reference])
        for kind, group in both.groupby("progress_kind_corrected"):
            r = ratio.loc[group.index]
            rows.append({"rate": rate, "shape": kind, "curves": len(group),
                         "median_log_ratio": float(r.median()),
                         "iqr_log_ratio": float(r.quantile(0.75)
                                                - r.quantile(0.25)),
                         "median_pH": float(group.pH.median())})
    return pd.DataFrame(rows).set_index(["rate", "shape"])


def divergence_against_ph(block=None, reference=DIVERGENCE_REFERENCE):
    """
    Rank correlation of log(rate / reference) with pH, per candidate.

    The confound check in one number per rate. A rate that departs from the
    baseline more at high pH than at low will move every pH order by that
    much, whatever the chemistry -- so a large |rho| here is a reason to read
    that candidate's pH row with care, not a reason to reject it.
    """
    data = scope.frame(scope.archive() if block is None else block)
    live = data[data.live]
    rows = []
    for rate in RATE_CANDIDATES:
        if rate in (reference, "vmax"):
            continue
        both = live[np.isfinite(live[rate]) & (live[rate] > 0)
                    & np.isfinite(live[reference]) & (live[reference] > 0)]
        ratio = np.log(both[rate] / both[reference])
        rho = float(ratio.rank().corr(both.pH.rank()))
        rows.append({"rate": rate, "curves": len(both), "rho_pH": rho})
    return pd.DataFrame(rows).set_index("rate")


def strong_runs_by_rate(rates=tuple(RATE_CANDIDATES)):
    """
    Which two-axis runs `scope.strong_runs` would keep under each candidate.

    The strong-run set is itself defined through a rate -- see
    `scope.strong_runs` -- so replacing `vmax` would change which runs every
    "over the strong runs" result is quoted over, as well as the response.
    Returns one row per rate: the runs kept, and those gained or lost against
    `vmax`'s set.
    """
    baseline = set(scope.strong_runs(parameter="vmax"))
    rows = []
    for rate in rates:
        kept = set(scope.strong_runs(parameter=rate))
        rows.append({"rate": rate, "runs": len(kept),
                     "gained": sorted(kept - baseline),
                     "lost": sorted(baseline - kept)})
    return pd.DataFrame(rows).set_index("rate")


# THE SHAPE SPLIT. A fitted rate means a different stage of the mechanism on
# a different shape of curve, so the rates are compared WITHIN a shape, where
# each has one reading. Each class names an EARLY rate and a LATE one, and the
# per-curve ratio of the two (a `scope.frame` column), whose order is exactly
# the difference of their orders on the same curves:
#
#   burst            one phase, B < 0, tau resolved. The classical
#                    pre-steady-state burst: v(0) the fast first turnover,
#                    v_ss the slower step that then sets turnover.
#   lag              one phase, B > 0, tau resolved. v(0) the rate before the
#                    catalyst has activated, v_ss the activated steady rate.
#   lag then fall    two phases, B1 > 0 > B2, the fitted maximum INSIDE the
#                    run. v(0) before activation, the peak the activated rate
#                    before the extent-dependent decline (`product_fate`). Its
#                    v_ss lies past that decline and is not used.
#
# "mixed" (B1 < 0 < B2) is left out on purpose: its two amplitudes trade
# against each other without moving the curve, which is what put its fitted
# peak at ~3x vmax (`shape_divergence`). So are unresolved one-phase fits --
# `ProgressFit.kind` labels those from the sign of B alone, without asking
# whether tau was located, so the resolution gate is applied here.
SHAPE_CLASSES = {
    "burst": {
        "early": "v0_fit_corrected", "late": "v_ss_fit_corrected",
        "ratio": "v0_over_ss_fit_corrected",
        "reading": "fast first turnover over the slower step that sets it",
    },
    "lag": {
        "early": "v0_fit_corrected", "late": "v_ss_fit_corrected",
        "ratio": "v0_over_ss_fit_corrected",
        "reading": "rate before activation over the activated steady rate",
    },
    "lag then fall": {
        "early": "v0_fit_corrected", "late": "v_peak_corrected",
        "ratio": "v0_over_peak_corrected",
        "reading": "rate before activation over the peak before the decline",
    },
}

# The per-curve ratios: dimensionless, so never divided by [enz].
RATIO_COLUMNS = tuple(sorted({c["ratio"] for c in SHAPE_CLASSES.values()}))

# Every column a shape restriction blanks outside its class.
SHAPE_RATE_COLUMNS = tuple(RATE_CANDIDATES) + RATIO_COLUMNS


def shape_members(data, name):
    """
    Which rows of `data` belong to `SHAPE_CLASSES[name]`, as a boolean Series.

    Live curves only: a dead curve has no shape worth the name, and every
    analysis drops it anyway.
    """
    kind = data.progress_kind_corrected
    live = data.live.astype(bool)
    if name in ("burst", "lag"):
        return (live & (data.phases_corrected == 1)
                & data.v0_fit_resolved_corrected.astype(bool) & (kind == name))
    if name == "lag then fall":
        peak_at = data.v_peak_corrected_time.astype(float)
        inside = (np.isfinite(peak_at) & (peak_at > 0)
                  & (peak_at <= data.duration_s))
        return live & (data.phases_corrected == 2) & (kind == name) & inside
    raise ValueError(f"unknown shape class {name!r}")


def shape_class_counts():
    """
    How many live curves each class holds, and where.

    `two-axis runs >= 2` is the number that decides whether the two-axis
    orders are identified at all inside a class: they carry one offset per
    run, so a run contributing a single curve of the class contributes
    nothing to the slope.
    """
    data = scope.frame(scope.archive())
    places = {"archive": scope.archive(),
              "two-axis": scope.TWO_AXIS_BLOCK,
              "temperature series": scope.TEMPERATURE_SERIES}
    places.update(scope.PH_LADDERS)
    rows = []
    for name in SHAPE_CLASSES:
        members = data[shape_members(data, name)]
        row = {"class": name}
        for label, runs in places.items():
            row[label] = int(members.experiment.isin(runs).sum())
        in_block = members[members.experiment.isin(scope.TWO_AXIS_BLOCK)]
        per_run = in_block.groupby("experiment").size()
        row["two-axis runs >= 2"] = int((per_run >= 2).sum())
        row["median pH"] = float(members.pH.median())
        rows.append(row)
    return pd.DataFrame(rows).set_index("class")


def shape_selection():
    """
    Does belonging to a class depend on the composition? The bias check.

    A curve's shape is an OUTCOME of its composition, so restricting an
    order to one shape conditions on something the axes may cause. If bursts
    sit at the top of a run's peroxide ladder, a burst-only peroxide order is
    read over a truncated axis and is biased by the selection itself.

    Two numbers per class. Within the two-axis block, the correlation of
    membership with log[S] and with log[H2O2] after each is centred on its own
    run -- the same within-run contrast the orders use, so it is the
    selection those orders actually see. And across the archive, the rank
    correlation of membership with pH, which is a between-run axis. Near zero
    means the class is a slice across the design; large means read that
    class's orders as conditional on shape.
    """
    data = scope.frame(scope.archive())
    live = data[data.live]
    block = live[live.experiment.isin(scope.TWO_AXIS_BLOCK)]
    rows = []
    for name in SHAPE_CLASSES:
        member = shape_members(block, name).astype(float)
        row = {"class": name, "two-axis members": int(member.sum())}
        centred = member - member.groupby(block.experiment).transform("mean")
        for axis in ("s0", "h2o2"):
            x = np.log(block[axis].astype(float))
            x = x - x.groupby(block.experiment).transform("mean")
            ok = centred.abs().sum() > 0 and x.abs().sum() > 0
            row[f"within-run r, log {axis}"] = (
                float(np.corrcoef(centred, x)[0, 1]) if ok else np.nan)
        everywhere = shape_members(live, name).astype(float)
        row["archive rank r, pH"] = float(
            everywhere.rank().corr(live.pH.rank()))
        rows.append(row)
    return pd.DataFrame(rows).set_index("class")


def shape_split(name, reference=DIVERGENCE_REFERENCE):
    """
    Every analysis, inside one shape class, on its early rate, its late rate,
    their ratio and the reference -- ALL ON THE SAME CURVES.

    The class's members are further held to the curves where all four are
    finite and positive, so the four columns of every row share their curves.
    That is what makes the ratio's row mean something: on the same rows an
    OLS fit of log(early/late) is exactly the fit of log(early) minus the fit
    of log(late), so its estimate is the difference of the two orders and its
    error is the error OF THAT DIFFERENCE, which two separate fits sharing
    their curves cannot give. (`test_rate_choice` checks the identity.)

    Returns the long table `order_comparison` does, one row per
    (analysis, rate), with `role` saying which of the four each rate is.

    READ THE LADDER ROWS WITH THEIR chi2, AND MOSTLY DO NOT READ THEM. Inside
    one class the pH ladders keep 4 to 12 curves each, against the three
    parameters each ladder fit carries, and the pooled value can sit exactly
    at zero from three ladders that disagree at chi2 = 31 (the burst ratio
    does). The two-axis rows are the ones a class can support, and
    `shape_class_counts`' `two-axis runs >= 2` says how well.

    TWO MORE LIMITS, both measured on 2026-09-11 and both properties of the
    archive rather than of this function:

      v(0) IS NOT A POSITIVE RATE ON A DEEP LAG. It is <= 0 on 31% of the
      lag class and 34% of lag-then-fall (17 of that class's 25 two-axis
      curves), so every order of v(0) inside a lag class is read over the
      SHALLOW lags only, and the deepest -- the coldest runs, on the
      temperature series -- drop out.

      ON THE TEMPERATURE SERIES THE TWO LAG CLASSES SIT ON DISJOINT
      TEMPERATURES: one-phase lags are exactly the 15, 20 and 30 C runs and
      lag-then-fall exactly 25, 35 and 40 C. The shape a curve earns is set
      by its temperature (and run length), so each class's activation energy
      is a three-temperature fit over part of the range, and the two classes
      cannot be compared as if they spanned the same one.
    """
    spec = SHAPE_CLASSES[name]
    roles = {"early": spec["early"], "late": spec["late"],
             "ratio": spec["ratio"], "reference": reference}
    frames = _frames(common=tuple(roles.values()), shape=name)
    rows = []
    for analysis_name, extra_is, analysis in RATE_ANALYSES:
        for role, rate in roles.items():
            estimate, stderr, curves, extra = analysis(rate, frames)
            rows.append({"class": name, "analysis": analysis_name,
                         "role": role, "rate": rate,
                         "estimate": float(estimate), "stderr": float(stderr),
                         "curves": int(curves), "extra": float(extra),
                         "extra_is": extra_is})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    print("\nwhat each candidate can be read on (whole archive, live curves)")
    print(candidate_coverage().to_string())
    for title, table in (
            ("every headline result, each candidate keeping the curves it can",
             order_comparison()),
            (f"each candidate beside {DIVERGENCE_REFERENCE} on the SAME curves",
             order_comparison(common=True))):
        print(f"\n{title}   estimate +/- error (curves) [extra]")
        wide = order_comparison_wide(table)
        for analysis in wide.index:
            print(f"\n  {analysis}")
            for rate in wide.columns:
                print(f"    {rate:<20} {wide.loc[analysis, rate]}")
    print(f"\nwhere the fitted rates part from {DIVERGENCE_REFERENCE}, "
          "by shape (whole archive)")
    print(shape_divergence().round(3).to_string())
    print("\n...and whether that departure tracks pH")
    print(divergence_against_ph().round(3).to_string())
    print("\nwhich runs would be strong under each rate")
    print(strong_runs_by_rate().to_string())
    print("\n\nTHE SHAPE SPLIT: where each class lives")
    print(shape_class_counts().to_string())
    print("\n...and whether membership tracks the composition (the bias check)")
    print(shape_selection().round(3).to_string())
    for name, spec in SHAPE_CLASSES.items():
        table = shape_split(name)
        print(f"\n{name.upper()}: {spec['reading']}   "
              "estimate +/- error (curves) [extra]")
        for analysis in [a for a, _, _ in RATE_ANALYSES]:
            rows = table[table.analysis == analysis]
            if not (rows.curves > 0).any():
                continue
            print(f"\n  {analysis}")
            for r in rows.itertuples():
                text = (f"{r.estimate:+.3f} +/- {r.stderr:.3f} ({r.curves})"
                        if np.isfinite(r.estimate) else f"-- ({r.curves})")
                if np.isfinite(r.extra):
                    text += f" [{r.extra:.2f}]"
                print(f"    {r.role:<10} {r.rate:<26} {text}")
