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

import induction
import rate_choice
import scope
import slowdown
from summary_kinetics import TWO_PHASE_F

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
