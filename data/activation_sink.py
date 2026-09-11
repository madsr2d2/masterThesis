"""
The activation-sink form, evaluated against the form every page draws.

    python data/activation_sink.py

`summary_kinetics.fit_activation_sink` fits one rate law to every curve,

    P' = v_act + (v0 - v_act) exp(-t/tau) - k P,

the archive's own findings written as one function: the catalyst relaxing
onto its activated rate on its own clock (`induction/`), the product drained
in proportion to how much of it there is (`product_fate/`), and a rate at
mixing free to sit below zero (`early_trough/`). `scope.curve_fit` makes it
once per curve, on the rebuilt series, as `CurveFit.activation_sink`, and
`scope.frame` carries its columns. This module asks four questions of it,
without fitting anything:

  form_contest         where it holds against the drawn `fit_progress`
                       form, shape by shape -- it nests inside the two-phase
                       form (v_ss = 0), so the two-phase form's one extra
                       parameter is an F test, `asymptote_f_corrected`
  resolution_table     which of its parameters each channel resolves, and on
                       how many curves only P''(0) is determined
  sink_agreement       its k and v_act against `slowdown.sink_fit`, which
                       reads the same two numbers model-free off the rolling
                       rate's tail -- the independent check
  own_clock_plus_one   the +1 rule through the form's own (v_act, tau) pair,
                       beside `vmax_corrected` through the same clock

`rate_choice.RATE_CANDIDATES` carries `v_act_corrected` and
`v_act_where_resolved_corrected`, so every headline result on v_act is
`rate_choice.order_comparison(rates=...)`, and `main` prints it.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import functools

import induction
import rate_choice
import scope
import slowdown
from summary_kinetics import (TWO_PHASE_F, fit_activation_inhibition,
                              fit_two_step_activation)

# The form's own clock, beside the two-phase form's slow clock that the
# published +1 row reads. Same (clock, gate, windowed) shape as
# `induction.JOINT_CLOCKS`.
ACTIVATION_SINK_CLOCKS = (
    ("tau_act_corrected", "tau_act_resolved_corrected", False),
    ("tau_slow_corrected", "tau_slow_resolved_corrected", False))

# The rates `main` puts through every headline analysis: the baseline, the
# drawn form's asymptote v_act replaces, and v_act itself both ways.
ACTIVATION_SINK_RATES = ("vmax_corrected", "v_ss_fit_corrected",
                         "v_act_corrected", "v_act_where_resolved_corrected")


def _live(block):
    data = scope.frame(scope.archive() if block is None else block)
    return data[data.live]


def form_contest(block=None):
    """
    The drawn form's kind against this form's, with the misfit where it fails.

    One row per (`progress_kind_corrected`, `act_sink_kind_corrected`):
    `curves`, `free_asymptote` (curves where the two-phase form's free v_ss
    clears TWO_PHASE_F over this form -- the curve levels onto a steady rate,
    or runs through zero, where a sink can only bring it to a plateau) and
    `resid_ratio`, the median of this form's residual over the drawn form's,
    both in noise units, on those curves. On the same curves, what the drawn
    form's asymptote does instead: `through_zero` counts a negative v_ss, and
    `floor_over_peak` is the median v_ss over the fitted peak where it is
    positive -- the steady rate a sink cannot level onto.
    """
    live = _live(block)
    free = live.asymptote_f_corrected > TWO_PHASE_F
    ratio = live.act_sink_resid_corrected / live.progress_resid_corrected
    floor = live.v_ss_fit_corrected / live.v_peak_corrected
    rows = []
    for (drawn, own), group in live.groupby(
            ["progress_kind_corrected", "act_sink_kind_corrected"]):
        failing = group[free.loc[group.index]]
        positive = floor.loc[failing.index][failing.v_ss_fit_corrected > 0]
        rows.append({"drawn": drawn, "activation_sink": own,
                     "curves": len(group),
                     "free_asymptote": len(failing),
                     "resid_ratio": (float(ratio.loc[failing.index].median())
                                     if len(failing) else np.nan),
                     "through_zero": int((failing.v_ss_fit_corrected
                                          < 0).sum()),
                     "floor_over_peak": (float(positive.median())
                                         if len(positive) else np.nan)})
    return pd.DataFrame(rows).set_index(["drawn", "activation_sink"])


def resolution_table(block=None):
    """
    Which parameters resolve, by channel (substrate x catalysed).

    `v_act`, `tau`, `v0`, `accel0` are resolved counts; `accel0_only` is the
    curves where v_act is not resolved and P''(0) is -- mostly runs still
    accelerating at their last reading, which determine v_act/tau and not
    v_act -- and `neither` the rest. `sink` counts curves whose k is resolved
    (the sink earned its parameter and its interval is inside the grid), and
    `free_asymptote` as in `form_contest`.
    """
    live = _live(block).assign(catalysed=lambda d: d.e0 > 0)
    rows = []
    for (substrate, catalysed), g in live.groupby(["substrate", "catalysed"]):
        v_act = g.v_act_resolved_corrected
        accel = g.accel0_resolved_corrected
        rows.append({
            "substrate": substrate, "catalysed": bool(catalysed),
            "curves": len(g), "v_act": int(v_act.sum()),
            "tau": int(g.tau_act_resolved_corrected.sum()),
            "v0": int(g.v0_act_resolved_corrected.sum()),
            "accel0": int(accel.sum()),
            "accel0_only": int((~v_act & accel).sum()),
            "neither": int((~v_act & ~accel).sum()),
            "sink": int((g.k_sink_state_corrected == "resolved").sum()),
            "free_asymptote": int((g.asymptote_f_corrected
                                   > TWO_PHASE_F).sum())})
    return pd.DataFrame(rows).set_index(["substrate", "catalysed"])


def sink_agreement(block=None):
    """
    k and v_act against `slowdown.sink_fit`'s model-free tail regression.

    Past the rate's maximum the activation is over and the form reduces to
    P' = v_act - kP, which is exactly the line `sink_fit` regresses the rolling
    rate on the product with: its slope is -k and its intercept v_act. So on a
    curve where both exist the two should agree, and they come from different
    places -- one a fitted closed form over the whole run, the other rolling
    slopes over the tail alone.

    One row per channel over live curves where the form's sink is resolved
    and `sink_fit` returns a positive k: `curves`, the median log ratio of
    each quantity (form over tail) and the correlation of their logs.
    `sink_fit` reads the READINGS, the form the rebuilt series; on a curve
    with no detachment the two are the same.
    """
    data = scope.frame(scope.archive() if block is None else block)
    index = data.set_index(["experiment", "sample"])
    rows = []
    for curve in scope.curves(scope.archive() if block is None else block):
        row = index.loc[(curve.experiment, curve.sample)]
        if not row.live or row.k_sink_state_corrected != "resolved":
            continue
        tail = slowdown.sink_fit(curve)
        if not (np.isfinite(tail.k) and tail.k > 0 and tail.v0 > 0
                and row.v_act_corrected > 0):
            continue
        rows.append({"substrate": row.substrate, "catalysed": row.e0 > 0,
                     "k_form": np.log(row.k_sink_corrected),
                     "k_tail": np.log(tail.k),
                     "v_form": np.log(row.v_act_corrected),
                     "v_tail": np.log(tail.v0)})
    pairs = pd.DataFrame(rows)
    out = []
    for (substrate, catalysed), g in pairs.groupby(["substrate", "catalysed"]):
        out.append({
            "substrate": substrate, "catalysed": bool(catalysed),
            "curves": len(g),
            "k_log_ratio": float((g.k_form - g.k_tail).median()),
            "k_r": float(np.corrcoef(g.k_form, g.k_tail)[0, 1]),
            "v_act_log_ratio": float((g.v_form - g.v_tail).median()),
            "v_act_r": float(np.corrcoef(g.v_form, g.v_tail)[0, 1])})
    return pd.DataFrame(out).set_index(["substrate", "catalysed"])


def own_clock_plus_one(block=scope.TWO_AXIS_BLOCK):
    """
    `induction.joint_clocks` through the form's own clock, on two rates.

    The +1 of a pre-equilibrium activation is d ln v - d ln tau on the axis
    of the activating species, where v is the ACTIVATED rate and tau the
    activation's clock -- which is exactly this form's (v_act, tau) pair, off
    one fit. The published row reads `vmax_corrected` against the two-phase
    form's `tau_slow`, a rate and a clock from different statistics. Rows are
    `joint_clocks`' own, indexed (rate, clock, role).
    """
    table = scope.frame(block)
    frames = {rate: induction.joint_clocks(table, rate=rate,
                                           clocks=ACTIVATION_SINK_CLOCKS)
              for rate in ("vmax_corrected", "v_act_corrected")}
    return pd.concat(frames, names=["rate"])


def contest_by_drawn_kind(block=None):
    """`form_contest` summed over this form's kinds: one row per drawn kind."""
    live = _live(block)
    free = live[live.asymptote_f_corrected > TWO_PHASE_F]
    ratio = free.act_sink_resid_corrected / free.progress_resid_corrected
    floor = free.v_ss_fit_corrected / free.v_peak_corrected
    rows = []
    for kind, group in live.groupby("progress_kind_corrected"):
        failing = free[free.progress_kind_corrected == kind]
        positive = floor.loc[failing.index][failing.v_ss_fit_corrected > 0]
        rows.append({"drawn": kind, "curves": len(group),
                     "free_asymptote": len(failing),
                     "through_zero": int((failing.v_ss_fit_corrected
                                          < 0).sum()),
                     "floor_over_peak": (float(positive.median())
                                         if len(positive) else np.nan),
                     "resid_ratio": (float(ratio.loc[failing.index].median())
                                     if len(failing) else np.nan)})
    table = pd.DataFrame(rows).set_index("drawn")
    table.attrs["resid_ratio_quantiles"] = ratio.quantile(
        [0.25, 0.5, 0.75, 0.9]).to_dict()
    return table


# The three five-parameter forms, by the column name `alternative_forms`
# gives each one's residual sum of squares.
FIVE_PARAMETER_FORMS = ("sink", "two_step", "inhibition")


@functools.lru_cache(maxsize=4)
def _alternative_forms(block):
    rows = []
    fits = scope.fits(block)
    data = scope.frame(block)
    for row in data[data.live].itertuples():
        shape = fits[(row.experiment, row.sample)]
        times, values = shape.times, shape.corrected
        sink = shape.activation_sink
        one = shape.progress_corrected.one
        two = shape.progress_corrected.two
        step = fit_two_step_activation(times, values)
        # Started from the sink (Ki = v_act/k matches its early slope-per-
        # product) and from the one-phase form (Ki far out, which IS it).
        starts = [(one.v0, one.v_ss, one.tau, 1e3 * max(np.ptp(values), 1e-6))]
        if sink.k > 0 and sink.v_act > 0:
            starts.append((sink.v0, sink.v_act, sink.tau, sink.v_act / sink.k))
        if np.isfinite(step.tau2):
            starts.append((step.v0, step.v_act, step.tau2,
                           1e3 * max(np.ptp(values), 1e-6)))
        inhibition = fit_activation_inhibition(times, values, starts=starts)
        rows.append({
            "experiment": row.experiment, "sample": row.sample,
            "substrate": row.substrate, "catalysed": bool(row.e0 > 0),
            "drawn": row.progress_kind_corrected,
            "points": len(times), "events": len(shape.events),
            "one_phase": float(one.sse), "two_phase": float(two.sse),
            # The sink's best k > 0 cost whether or not k was earned: `sse` is
            # the SELECTED form's, which is k = 0 where the sink did not pay,
            # and a four-parameter cost read against one degree of freedom
            # overstates its failure. sink_f = (sse0 - sse5)/(sse5/(n - 5))
            # inverts exactly to sse5 = sse0 / (1 + sink_f/(n - 5)).
            "sink": float(sink.sse_no_sink
                          / (1 + sink.sink_f / max(1, len(times) - 5))
                          if np.isfinite(sink.sink_f) else sink.sse),
            "two_step": float(step.sse),
            "inhibition": float(inhibition.sse),
            "step_earned": step.step_earned,
            "tau1": step.tau1, "tau2": step.tau2, "ki": inhibition.ki,
            "inhibition_v_act": inhibition.v_act})
    return pd.DataFrame(rows)


def alternative_forms(block=None):
    """
    Every live curve's residual sum of squares under the three five-parameter
    forms -- the activation-sink, the two-step activation and the activation-
    inhibition -- beside the one- and two-phase forms, all on the rebuilt
    series. `summary_kinetics` has the three forms' chemistry.

    Adds, per form, `<form>_f`: the two-phase form's one extra parameter over
    it, ((sse_form - sse_two) / 1) / (sse_two / (n - 6)). The sink and the
    two-step forms nest inside the two-phase form, so for them this is an F
    test; the inhibition form does not, and for it the same number is a
    comparison on the same scale, not a test. `best` is the lowest-cost of the
    three and `adequate` whether it stays within TWO_PHASE_F of the two-phase
    form -- whether the curve needs the free asymptote at all once the right
    five-parameter chemistry is asked.

    MEMOISED: the inhibition form is a nonlinear fit from up to three starts,
    about 15 s over the archive.
    """
    table = _alternative_forms(scope.archive() if block is None
                               else block).copy()
    degrees = (table.points - 6).clip(lower=1)
    for form in FIVE_PARAMETER_FORMS:
        table[f"{form}_f"] = ((table[form] - table.two_phase)
                              / (table.two_phase / degrees))
    costs = table[list(FIVE_PARAMETER_FORMS)]
    table["best"] = costs.idxmin(axis=1)
    table["best_f"] = table[[f"{f}_f" for f in FIVE_PARAMETER_FORMS]].min(axis=1)
    table["adequate"] = table.best_f <= TWO_PHASE_F
    return table


def remaining_curves(block=None):
    """
    The curves the activation-sink form cannot hold, and what each one is.

    "Cannot hold" is `asymptote_f_corrected` > TWO_PHASE_F -- the two-phase
    form's free asymptote earned over the form the curve was reported in. One
    row each, with `group`:

      gas            carries a detachment, so the misfit may be the
                     reconstruction's and is not read as chemistry
      second rise    drawn "mixed" or "two lags": a rate that falls and then
                     rises again, which one relaxation cannot do
      through zero   the drawn asymptote is negative
      floor          the drawn asymptote is a positive steady rate

    and, from `alternative_forms`, which five-parameter chemistry fits it best
    and whether that one is `rescued` -- within TWO_PHASE_F of the two-phase
    form on one degree of freedom, so the free asymptote is no longer needed.
    """
    data = _live(block)
    failing = data[data.asymptote_f_corrected > TWO_PHASE_F]
    fits = scope.fits(scope.archive() if block is None else block)
    forms = alternative_forms(block).set_index(["experiment", "sample"])
    rows = []
    for row in failing.itertuples():
        key = (row.experiment, row.sample)
        shape = fits[key]
        if shape.events:
            group = "gas"
        elif row.progress_kind_corrected in ("mixed", "two lags"):
            group = "second rise"
        elif row.v_ss_fit_corrected < 0:
            group = "through zero"
        else:
            group = "floor"
        form = forms.loc[key]
        rows.append({"experiment": row.experiment, "sample": row.sample,
                     "substrate": row.substrate,
                     "catalysed": bool(row.e0 > 0),
                     "drawn": row.progress_kind_corrected, "group": group,
                     "best": form.best, "rescued": bool(form.adequate),
                     "two_step_f": form.two_step_f,
                     "inhibition_f": form.inhibition_f,
                     "sink_f": form.sink_f})
    return pd.DataFrame(rows)


def remaining_summary(block=None):
    """`remaining_curves` counted: per group, how many, and what rescues them."""
    table = remaining_curves(block)
    out = []
    for group, g in table.groupby("group"):
        out.append({"group": group, "curves": len(g),
                    "catalysed": int(g.catalysed.sum()),
                    "rescued_two_step": int(((g.best == "two_step")
                                             & g.rescued).sum()),
                    "rescued_inhibition": int(((g.best == "inhibition")
                                               & g.rescued).sum()),
                    "best_two_step": int((g.best == "two_step").sum()),
                    "best_inhibition": int((g.best == "inhibition").sum()),
                    "best_sink": int((g.best == "sink").sum())})
    return pd.DataFrame(out).set_index("group")


def second_rise_needs_the_catalyst(block=None):
    """
    How many live curves the drawn form calls "mixed" -- a rate that falls
    and then rises -- with and without the catalyst. The control
    `induction/` uses for the lag: if the pattern were mixing, thermal
    equilibration or the cell, the enzyme-free cuvettes would carry it too.
    """
    live = _live(block)
    return (live.assign(catalysed=live.e0 > 0, mixed=lambda d:
                        d.progress_kind_corrected == "mixed")
            .groupby("catalysed").agg(curves=("mixed", "size"),
                                      mixed=("mixed", "sum")))


def second_rise_clocks(block=None):
    """
    The late rise of a "mixed" curve against its own run's lags.

    If the rise is the catalyst activating after an early burst, its clock is
    the run's activation clock, which the one-phase lags beside it measure.
    One row per (mixed curve with a resolved slow clock, one-phase lag in the
    same run with a resolved clock): both clocks and their ratio. The
    reference is `lag_pairs`, the same |log ratio| between two lags of one
    run -- how far apart two readings of the SAME clock sit here.
    """
    live = _live(block)
    fits = scope.fits(scope.archive() if block is None else block)
    lags = live[(live.phases_corrected == 1)
                & (live.progress_kind_corrected == "lag")
                & live.tau_resolved_corrected]

    def clock(row):
        return fits[(row.experiment, row.sample)].progress_corrected.one.tau
    rows = []
    for row in live[live.progress_kind_corrected == "mixed"].itertuples():
        two = fits[(row.experiment, row.sample)].progress_corrected.two
        if not two.resolved:
            continue
        for lag in lags[lags.experiment == row.experiment].itertuples():
            rows.append({"experiment": row.experiment, "mixed": row.sample,
                         "lag": lag.sample, "rise_tau": float(two.tau2),
                         "lag_tau": float(clock(lag))})
    pairs = pd.DataFrame(rows)
    pairs["ratio"] = pairs.rise_tau / pairs.lag_tau
    reference = []
    for _, run in lags.groupby("experiment"):
        taus = [clock(r) for r in run.itertuples()]
        reference += [abs(np.log(a / b)) for i, a in enumerate(taus)
                      for b in taus[i + 1:]]
    pairs.attrs["lag_pairs_abs_log_ratio"] = float(np.median(reference))
    pairs.attrs["lag_pairs"] = len(reference)
    return pairs


def sink_versus_inhibition(block=None):
    """
    The sink and the inhibition form head to head, same parameter count, on
    every live curve whose sink is resolved: per channel, how many favour
    each by more than 2 in AIC (n ln(sse_sink / sse_inhibition), positive for
    inhibition). `product_fate` makes the same comparison on the tail's
    rolling rate; this is the whole curve.
    """
    forms = alternative_forms(block)
    frame = _live(block).set_index(["experiment", "sample"])
    state = frame.k_sink_state_corrected.reindex(
        pd.MultiIndex.from_frame(forms[["experiment", "sample"]])).to_numpy()
    resolved = forms[state == "resolved"].copy()
    resolved["delta_aic"] = resolved.points * np.log(resolved.sink
                                                     / resolved.inhibition)
    return resolved.groupby(["substrate", "catalysed"]).agg(
        curves=("delta_aic", "size"),
        inhibition=("delta_aic", lambda d: int((d > 2).sum())),
        sink=("delta_aic", lambda d: int((d < -2).sum())),
        median_delta_aic=("delta_aic", "median"))


def main():
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    print("\nthe drawn form's kind against the activation-sink form's "
          "(whole archive, live)")
    print(form_contest().round(2).to_string())
    by_kind = contest_by_drawn_kind()
    print(by_kind.round(2).to_string())
    print("residual ratio where the free asymptote wins, quantiles:",
          {q: round(v, 2) for q, v in by_kind.attrs["resid_ratio_quantiles"]
           .items()})
    print("\nwhat each channel resolves")
    print(resolution_table().to_string())
    print("\nk and v_act against slowdown.sink_fit's tail regression")
    print(sink_agreement().round(2).to_string())
    print("\nthe +1 rule through the form's own clock (two-axis block)")
    print(own_clock_plus_one().round(3).to_string())
    print("\nthe curves the form cannot hold, and what each one is")
    print(remaining_summary().to_string())
    print("\n...a second rise needs the catalyst")
    print(second_rise_needs_the_catalyst().to_string())
    pairs = second_rise_clocks()
    print("\n...and runs on its run's own activation clock")
    print(pairs.round(2).to_string())
    print(f"median ratio {pairs.ratio.median():.2f}; |log ratio| "
          f"{np.abs(np.log(pairs.ratio)).median():.2f} against "
          f"{pairs.attrs['lag_pairs_abs_log_ratio']:.2f} between two lags of "
          f"one run ({pairs.attrs['lag_pairs']} pairs)")
    print("\nsink against inhibition, whole curve, where the sink is resolved")
    print(sink_versus_inhibition().round(1).to_string())
    print("\nwhat each rate can be read on")
    print(rate_choice.candidate_coverage()
          .loc[list(ACTIVATION_SINK_RATES)].to_string())
    table = rate_choice.order_comparison(common=True,
                                         rates=ACTIVATION_SINK_RATES)
    wide = rate_choice.order_comparison_wide(table)
    print(f"\nevery headline result, each rate beside "
          f"{rate_choice.DIVERGENCE_REFERENCE} on the SAME curves")
    for analysis in wide.index:
        print(f"\n  {analysis}")
        for rate in wide.columns:
            print(f"    {rate:<32} {wide.loc[analysis, rate]}")


if __name__ == "__main__":
    main()
