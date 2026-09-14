"""
Stage A of `PLAN_CURVES_TO_MECHANISM.md`: the curve parameters as functions of
the conditions.

Every catalysed curve in the archive carries an activation-sink fit
(`summary_kinetics.fit_activation_sink`, on the gas-corrected readings):

    dP/dt = v_act + (v0 - v_act) exp(-t/tau) - k P

Four parameters describe the curve's shape -- the rate it relaxes onto, the
rate at mixing, the relaxation clock and the product-loss constant. This module
turns those into the four element tables the plan's rate-law search is run on,
per substrate:

    v_act        log(v_act), admitted where the fit resolves it and it is
                 positive with a positive lower interval end
    k_act_lag    log(1/tau) on lag curves (v0 < v_act: the activation clock)
    k_act_burst  log(1/tau) on burst curves (v0 > v_act: the burst's decay)
    k_sink       log(k), admitted where `k_state == "resolved"`

THE ROWS. Every row of `scope.frame(scope.archive())` with `live` true and
`e0 > 0`, minus the curves `activation_sink.remaining_curves()` lists as
catalysed: on those the two-phase form's free asymptote earns over the
activation-sink form, so the form cannot hold them and no parameter of it is
read. Enzyme-free curves are not used at all -- a catalysed curve is the
catalytic INCREMENT, its reference cuvette holding the same mixture minus the
enzyme.

THE ERRORS. Each admitted row carries `se`, the log-scale standard error read
off the fit's own profile interval, (log(high) - log(low)) / (2 * 1.96), and
`weight = 1 / (se^2 + SE_FLOOR^2)`. The floor stops a few very tight curves
from deciding a table.

WHAT THIS MODULE DOES NOT DO. It does not fit anything. The rate-law models,
their health checks and the search live here too (Tasks 3-7), added in their
own commits.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import activation_sink
import scope

# The four element tables, in the order the plan lists them.
ELEMENTS = ("v_act", "k_act_lag", "k_act_burst", "k_sink")

# Why a curve is not in an element table. The first is common to every
# element: those curves were removed from the row set before any element was
# built. The other three are applied per element in this order.
REASONS = ("form cannot hold the curve", "not resolved",
           "non-positive estimate", "non-positive interval end")

# The log-scale error below which a curve's own interval stops narrowing its
# weight. A curve with se = 0.01 and one with se = 0.10 then carry the same
# weight, so the search is not decided by the few tightest interval widths.
SE_FLOOR = 0.10


def _interval_low(intervals):
    """The low end of each (low, high) profile interval, as floats."""
    return np.array([interval[0] for interval in intervals], dtype=float)


def _interval_high(intervals):
    """The high end of each (low, high) profile interval, as floats."""
    return np.array([interval[1] for interval in intervals], dtype=float)


def _log_se(interval):
    """
    The log-scale standard error of a profile interval: half its log width
    over 1.96 sigmas. On (e^-1, e^1) this is 2 / (2 * 1.96).
    """
    low, high = interval
    return (np.log(high) - np.log(low)) / (2 * 1.96)


def _curve_rows(substrate):
    """
    Every live catalysed curve of one substrate, with its sink-fit columns.

    Not yet filtered by the form's fitness -- `form_holds` marks the curves
    `activation_sink.remaining_curves()` removes, so the exclusions can be
    counted before they are applied. `family` is the sink fit's own `kind`:
    "lag" where `v0 < v_act`, "burst" where `v0 > v_act`.
    """
    frame = scope.frame(scope.archive())
    fits = scope.fits(scope.archive())
    remaining = activation_sink.remaining_curves()
    catalysed = remaining[remaining.catalysed]
    removed = set(zip(catalysed.experiment.astype(int),
                      catalysed["sample"].astype(int)))
    live = frame[frame.live & (frame.e0 > 0)
                 & (frame.substrate == substrate)]
    records = []
    for row in live.itertuples():
        fit = fits[(int(row.experiment), int(row.sample))].activation_sink
        records.append({
            "experiment": int(row.experiment),
            "sample": int(row.sample),
            "buffer": row.buffer,
            "temperature": float(row.temperature),
            "kelvin": float(row.temperature) + 273.15,
            "pH": float(row.pH),
            "s0": float(row.s0),
            "h2o2": float(row.h2o2),
            "hoo": float(row.hoo),
            "buf": float(row.buf),
            "e0": float(row.e0),
            "bubble_load": float(row.bubble_load),
            "family": "lag" if fit.kind.startswith("lag") else "burst",
            "v_act": float(fit.v_act),
            "v0": float(fit.v0),
            "tau": float(fit.tau),
            "k": float(fit.k),
            "v_act_interval": tuple(fit.v_act_interval),
            "v0_interval": tuple(fit.v0_interval),
            "tau_interval": tuple(fit.tau_interval),
            "k_interval": tuple(fit.k_interval),
            "v_act_resolved": bool(fit.v_act_resolved),
            "tau_resolved": bool(fit.tau_resolved),
            "k_state": fit.k_state,
            "sink_earned": bool(fit.sink_earned),
            "form_holds": (int(row.experiment), int(row.sample)) not in removed,
        })
    return pd.DataFrame(records)


def _element_candidates(all_rows, element):
    """
    Which rows an element could ever admit, before any resolution test.

    `v_act` and `k_sink` are properties every curve has; the two clocks are
    only defined on their own family, so a burst curve is no more a failed
    `k_act_lag` than an enzyme-free curve is a failed catalysed one.
    """
    if element in ("k_act_lag", "k_act_burst"):
        return (all_rows.family == element[len("k_act_"):]).to_numpy()
    return np.ones(len(all_rows), dtype=bool)


def _element_reasons(all_rows, element):
    """
    For every row, why it is not in `element`'s table, or None if it is.

    Returned over `all_rows`, one entry per row, so the caller can both select
    the table and count the exclusions without iterating twice.
    """
    if element == "v_act":
        resolved = all_rows.v_act_resolved.to_numpy()
        estimate = all_rows.v_act.to_numpy(dtype=float)
        low = _interval_low(all_rows.v_act_interval)
    elif element in ("k_act_lag", "k_act_burst"):
        resolved = all_rows.tau_resolved.to_numpy()
        estimate = all_rows.tau.to_numpy(dtype=float)
        low = _interval_low(all_rows.tau_interval)
    else:
        resolved = (all_rows.k_state == "resolved").to_numpy()
        estimate = all_rows.k.to_numpy(dtype=float)
        low = _interval_low(all_rows.k_interval)

    candidate = _element_candidates(all_rows, element)
    held = all_rows.form_holds.to_numpy()
    reason = np.full(len(all_rows), None, dtype=object)
    reason[candidate & ~held] = "form cannot hold the curve"
    reason[candidate & held & ~resolved] = "not resolved"
    reason[candidate & held & resolved & (estimate <= 0)] = \
        "non-positive estimate"
    reason[candidate & held & resolved & (estimate > 0)
           & ~(np.isfinite(low) & (low > 0))] = "non-positive interval end"
    return reason


def _element_column(rows, element):
    """(y, se) on the log scale for one element's admitted rows."""
    if element == "v_act":
        y = np.log(rows.v_act.to_numpy(dtype=float))
        se = np.array([_log_se(i) for i in rows.v_act_interval])
    elif element in ("k_act_lag", "k_act_burst"):
        y = -np.log(rows.tau.to_numpy(dtype=float))
        se = np.array([_log_se(i) for i in rows.tau_interval])
    else:
        y = np.log(rows.k.to_numpy(dtype=float))
        se = np.array([_log_se(i) for i in rows.k_interval])
    return y, se


def curve_parameters(substrate):
    """
    Every table Stage A fits for one substrate, with every exclusion counted.

    Returns a dict:

      rows       one row per catalysed curve the activation-sink form holds,
                 with the conditions, the family and the fit's four parameters
                 and their intervals;
      elements   {"v_act", "k_act_lag", "k_act_burst", "k_sink"} -> DataFrame,
                 each carrying `y`, `se` and `weight` on top of `rows`'
                 columns;
      excluded   (element, reason, count) over `REASONS`;
      remaining  the `activation_sink.remaining_curves()` rows for this
                 substrate that were removed.

    The expensive work is shared and memoised: `scope.frame` and `scope.fits`
    are built once for the archive, `activation_sink.remaining_curves` once.
    """
    all_rows = _curve_rows(substrate)
    rows = all_rows[all_rows.form_holds].copy()
    elements, excluded = {}, []
    for element in ELEMENTS:
        candidate = _element_candidates(all_rows, element)
        reason = _element_reasons(all_rows, element)
        admitted = candidate & (reason == None)                        # noqa: E711
        table = all_rows[admitted].copy()
        y, se = _element_column(table, element)
        table["y"] = y
        table["se"] = se
        table["weight"] = 1.0 / (se ** 2 + SE_FLOOR ** 2)
        elements[element] = table
        for name in REASONS:
            excluded.append({"element": element, "reason": name,
                             "count": int((reason == name).sum())})
    remaining = activation_sink.remaining_curves()
    remaining = remaining[remaining.catalysed
                          & (remaining.substrate == substrate)]
    return {"rows": rows, "elements": elements,
            "excluded": pd.DataFrame(excluded), "remaining": remaining}
