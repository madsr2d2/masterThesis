"""
What pH does to the reaction: the rate and the induction clock, read the same
way across every matched pH ladder the archive holds.

`scope.PH_LADDERS` names four of them -- built 2026-09-05 for the induction
clock alone (`induction.lag_ph_ladders`, `induction.pooled_ladder`) and never
asked about the RATE. This module is that second half: the rate's order in
[HOO-], pooled the same way, with the same confound controls.

    python data/ph_role.py            the whole argument, printed
    python data/test_ph_role.py       planted recovery + the published numbers

Related: `induction.lag_ph_ladders`/`pooled_ladder` own the clock's pH order
(+0.16 to +0.33 per pH unit, chi2 0.95 on 3); `scope.ph_order`/`arm_orders`/
`hoo_consistency` own the two-axis block's own cuvette-matched reading, which
this module's cruder pooled fit is checked against; `buffer/ANALYSIS.md` S5
owns why buffer identity cannot be separated from pH anywhere in this archive.
"""
import functools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import induction
import scope

# [S] is a genuine within-run ladder on every one of the four pH ladders --
# unlike the clock, which reads one window-free number per curve, `vmax`
# varies with [S] inside every single run of every ladder. Fitting it
# alongside [HOO-] in one regression is `scope.orders`'s own case for fitting
# every axis a design moves at once: one axis at a time is a different
# regression on an L.
RATE_TERMS = ("s0", "hoo")


def _ladder_scope(experiments):
    """
    The runs a rate order should actually be read over.

    `[buf]` sits exactly fixed inside all four ladders (checked once, by
    hand, when this module was written -- 80.0, 85.0 and 75.013 mM with no
    within-ladder spread), so there is no confound there to filter for. The
    one filter that DOES apply is `scope.strong_runs()`, and only where a
    ladder is a subset of the two-axis block: `concentration_agreement`
    is built on that block's L design (both axes moving inside a run) and
    returns NaN off it, so it has nothing to say about the phosphate or
    boric ladders, whose own per-run corr([S], vmax) already runs 0.56-1.00
    with two exceptions at pH 5.6-5.9 (checked by hand: exps 9 and 10, where
    the reaction is simply slow). Nothing there rises to `strong_runs`'s
    reason for existing, which is a run whose cuvettes stop tracking their
    own composition at all.
    """
    experiments = tuple(experiments)
    if set(experiments) <= set(scope.TWO_AXIS_BLOCK):
        strong = set(scope.strong_runs())
        return tuple(e for e in experiments if e in strong)
    return experiments


def rate_ladder(experiments, parameter="vmax", terms=RATE_TERMS, frame=None):
    """
    The rate's order in [S] and [HOO-], pooled over one pH ladder.

    Every one of the four ladders holds pH fixed WITHIN a run and steps it
    only BETWEEN runs (`induction.lag_identifiability`), so there is no
    per-run offset that could carry a pH term without absorbing it --
    `scope.orders(within=False)` is the tool, not the `within=True` the
    two-axis block's own concentration orders use. That costs the day-to-day
    control a per-run offset would give the substrate axis; `schedule_
    collinearity` below is how that cost is priced instead of ignored.

    `frame` overrides the archive lookup with a caller-built table, the way
    `buffer_role.catalytic_coefficient` does, so a planted recovery test does
    not need real experiment numbers. A caller who passes it also gets no
    schedule control, since a planted frame has no `Date Collected`.
    """
    if frame is None:
        kept = _ladder_scope(experiments)
        data = scope.frame(kept)
    else:
        kept, data = None, frame
    fit = scope.orders(parameter=parameter, frame=data, terms=terms,
                       within=False)
    live = data[data.live]
    schedule_corr = np.nan
    if kept is not None:
        by_run = live.groupby("experiment").pH.first()
        dates = scope.run_dates(kept)
        common = by_run.index.intersection(dates.index[dates.order.notna()])
        if len(common) >= 3:
            schedule_corr = float(np.corrcoef(
                by_run.loc[common], dates.loc[common, "order"])[0, 1])
    return {
        "runs": int(live.experiment.nunique()) if len(live) else 0,
        "curves": int(len(live)),
        "pH_low": float(live.pH.min()) if len(live) else np.nan,
        "pH_high": float(live.pH.max()) if len(live) else np.nan,
        "hoo": {"slope": fit.get("order_hoo", np.nan),
               "stderr": fit.get("stderr_hoo", np.nan)},
        "s0": {"slope": fit.get("order_s0", np.nan),
              "stderr": fit.get("stderr_s0", np.nan)},
        "r2": fit["r2"],
        "n": fit["n"],
        "schedule_collinearity": schedule_corr,
    }


# The boric ladder does not sit flat across its own range: median vmax rises
# from 7.1e-5 at pH 8.46 to a peak of 1.3-1.4e-4 around pH 9.2-9.7 and falls
# back to 6.6e-5 by pH 10.34 (checked by hand from `scope.frame`, medians per
# experiment). A single log-log slope across a ladder that rises and then
# falls reports the two halves cancelling, which is what `rate_ladder` above
# returns for the whole nine-run set -- it is not evidence of no pH
# dependence, it is evidence of a TURNOVER `rate_ladder` cannot see because it
# only ever fits a straight line. 9.51 sits just above the ladder's own peak
# (exp 43, pH 9.50) and below the two runs where the decline is unambiguous.
BORIC_TURNOVER_SPLIT = 9.51


def boric_turnover(split=BORIC_TURNOVER_SPLIT):
    """
    The boric ladder split at its own peak: an order below it, a decline above.

    Below the split there are still only six runs and no substrate-ladder
    partner outside the two-axis block to triangulate against, so this is
    reported as what it is -- a weaker positive order than phosphate or
    pyrophosphate's, not a third confirmation of their +0.59 -- and the two
    runs above the split (pH 10.07, 10.34) are too few to fit at all, so their
    decline is reported as the medians rather than as a slope.

    `vmax_corrected` barely moves either half (checked by hand: the pH 10.07
    and 10.34 medians shift by less than 20% of themselves), which is the
    argument that the decline is not the O2 side reaction eating the reading
    -- `scope.gas_survey` puts boric's heaviest gassing in exactly this pH
    band, so the O2 explanation had to be checked and ruled out rather than
    assumed away.
    """
    data = scope.frame(scope.PH_LADDER_BORIC)
    below = tuple(sorted(int(e) for e in
                         data.loc[data.pH <= split, "experiment"].unique()))
    above = data[data.pH > split]
    by_run = data.groupby("experiment").vmax.median()
    return {
        "split": split,
        "below_runs": below,
        "below": rate_ladder(below, parameter="vmax"),
        "peak_experiment": int(by_run.idxmax()) if len(by_run) else None,
        "peak_pH": float(data.loc[data.experiment == by_run.idxmax(),
                                  "pH"].iloc[0]) if len(by_run) else None,
        "above_by_experiment": {
            int(exp): {"pH": float(g.pH.iloc[0]),
                      "median_vmax": float(g.vmax.median()),
                      "median_vmax_corrected": float(g.vmax_corrected.median())}
            for exp, g in above.groupby("experiment")},
    }


def rate_ladders(parameter="vmax", terms=RATE_TERMS, ladders=None):
    """Every pH ladder's rate order, read the same way. `scope.PH_LADDERS`."""
    ladders = scope.PH_LADDERS if ladders is None else ladders
    rows = []
    for name, experiments in ladders.items():
        result = rate_ladder(experiments, parameter=parameter, terms=terms)
        result["ladder"] = name
        rows.append(result)
    return rows


def rate_ladder_table(results=None):
    """`rate_ladders`, flattened to the table ANALYSIS.md and the figures read."""
    results = rate_ladders() if results is None else results
    rows = [{
        "ladder": row["ladder"], "runs": row["runs"], "curves": row["curves"],
        "pH_low": row["pH_low"], "pH_high": row["pH_high"],
        "order_hoo": row["hoo"]["slope"], "stderr_hoo": row["hoo"]["stderr"],
        "order_s0": row["s0"]["slope"], "stderr_s0": row["s0"]["stderr"],
        "r2": row["r2"], "n": row["n"],
        "schedule_collinearity": row["schedule_collinearity"],
    } for row in results]
    return pd.DataFrame(rows).set_index("ladder")


def pooled_rate_order(results=None, axis="hoo", drop=()):
    """
    Inverse-variance pooled rate order in [HOO-], with `chi2`/`dof` deciding
    whether pooling is honest.

    `drop` names ladders to exclude by label, so the same call answers "what
    do all four say" and "what do the ones that agree say" -- the boric
    ladder alone drives `chi2` from 0.29 on 2 (three ladders, essentially one
    number) to 174 on 3 the moment it is put back, which is the finding
    rather than a nuisance to average away. Delegates to
    `induction.pooled_ladder`, which this module does not redefine.
    """
    results = rate_ladders() if results is None else results
    kept = [row for row in results if row["ladder"] not in drop]
    return induction.pooled_ladder(kept, response=axis, controlled=False)


def clock_pooled_order(response="lag_half_s"):
    """
    The induction clock's own pooled pH order, quoted back from `induction`.

    Not recomputed: `induction.lag_ph_ladders` + `pooled_ladder` already do
    this, at a window every run in a ladder shares and with the signal
    control this module's rate side has no equivalent for (a rate does not
    have a rolling-window landmark to fail).
    """
    return induction.pooled_ladder(induction.lag_ph_ladders(),
                                   response=response, controlled=True)


def main():
    table = rate_ladder_table()
    print("the rate's order in [HOO-], per ladder")
    print(table[["runs", "pH_low", "pH_high", "order_hoo", "stderr_hoo",
                "order_s0", "stderr_s0", "schedule_collinearity"]]
         .to_string(float_format=lambda v: f"{v:.4f}"))
    pooled_all = pooled_rate_order()
    pooled_three = pooled_rate_order(drop=("boric 4OMe",))
    print(f"\npooled, all four: {pooled_all['pooled']:+.3f} +/- "
         f"{pooled_all['stderr']:.3f}, chi2 = {pooled_all['chi2']:.1f} on "
         f"{pooled_all['dof']}")
    print(f"pooled, without boric: {pooled_three['pooled']:+.3f} +/- "
         f"{pooled_three['stderr']:.3f}, chi2 = {pooled_three['chi2']:.2f} "
         f"on {pooled_three['dof']}")
    clock = clock_pooled_order()
    print(f"\nthe clock's own pooled order, for comparison: "
         f"{clock['pooled']:+.3f} +/- {clock['stderr']:.3f} per pH unit, "
         f"chi2 = {clock['chi2']:.2f} on {clock['dof']}")

    turnover = boric_turnover()
    print(f"\nthe boric ladder turns over rather than saturating: peak at "
         f"exp {turnover['peak_experiment']}, pH {turnover['peak_pH']:.2f}")
    print(f"below the peak: order_hoo = "
         f"{turnover['below']['hoo']['slope']:+.3f} +/- "
         f"{turnover['below']['hoo']['stderr']:.3f}, runs "
         f"{turnover['below_runs']}")
    for exp, row in turnover["above_by_experiment"].items():
        print(f"  exp {exp}, pH {row['pH']:.2f}: median vmax "
             f"{row['median_vmax']:.3e}, corrected "
             f"{row['median_vmax_corrected']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
