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

THE RATE LAWS. A model is at most one option per family:

    S   power | mm          H   power | bind | relax
    HOO power               BUF power | bind | relax | power_by_buffer
    E   power               T   arrhenius

and `None` for a family is always allowed. `law_design` (the plan's `design`,
renamed because `scope.design` holds that name) builds the weighted design at
fixed nonlinear constants, `fit_model` profiles the linear coefficients and
optimises the nonlinear constants in log10, and `rate_law_health` refuses to
let a coefficient be quoted off a fit that is at a bound, collinear, rank
deficient, carried by fewer than three runs or fitted on fewer than fifteen
curves. `term_identifiable` and `within_run_check` are the design-side guards:
what this table can separate at all, and what survives replacing the buffer
intercepts with one per run.

WHAT THIS MODULE DOES NOT DO. Nothing here searches models or fits readings;
the search (Task 4), the links (Task 5) and the global fitter (Tasks 6-7) are
added in their own commits.
"""
import itertools
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import activation_sink
import fit_dataset
import scope
import summary_kinetics

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

# The families a rate law may carry a term for, in design order.
FAMILIES = ("S", "H", "HOO", "BUF", "E", "T")

# What each element is allowed to ask of each family. "None" is always allowed
# and is not listed. The rate elements may saturate in S and H; the clocks may
# relax with them instead -- a pre-equilibrium's clock RISES with the species
# where its rate saturates.
_RATE_OPTIONS = {
    "S": ("power", "mm"),
    "H": ("power", "bind"),
    "HOO": ("power",),
    "BUF": ("power", "bind", "power_by_buffer"),
    "E": ("power",),
    "T": ("arrhenius",),
}
_CLOCK_OPTIONS = {
    "S": ("power",),
    "H": ("power", "relax"),
    "HOO": ("power",),
    "BUF": ("power", "relax", "power_by_buffer"),
    "E": ("power",),
    "T": ("arrhenius",),
}
TERM_OPTIONS = {
    "v_act": dict(_RATE_OPTIONS),
    "k_act_lag": dict(_CLOCK_OPTIONS),
    "k_act_burst": dict(_CLOCK_OPTIONS),
    "k_sink": dict(_RATE_OPTIONS),
}

# The nonlinear constant each option carries, if any. Optimised in log10.
NONLINEAR_OF = {
    ("S", "mm"): "Km",
    ("H", "bind"): "K_H",
    ("H", "relax"): "K_H",
    ("BUF", "bind"): "K_B",
    ("BUF", "relax"): "K_B",
}
NONLINEAR_BOUNDS = (-4.0, 4.0)

# The weighted design's term names, so a coefficient reads as what it is.
_T_REGRESSOR_BELOW = 1e-5         # in 1/K, for the identifiability rule
_LOG_REGRESSOR_BELOW = 0.05       # in log units, for every other family
_CORRELATION_BAR = 0.95
_GAS_CONSTANT = 8.314462618       # J/(mol K), for Ea in kJ/mol
# The floor the tie rule's log test puts under a fold score before taking its
# logarithm. Fold scores are non-negative weighted sums and a fold can score
# exactly zero, which has no logarithm.
TIE_LOG_FLOOR = 1e-12

# A numerical guard, not a constraint: the global fitter's optimiser may
# wander to coefficients whose exponent overflows before it turns back, and a
# NaN in the residual aborts the run. Clipping the LOG predictions at +/-100
# leaves every realistic rate (log AU/s from -14 to -5) untouched and only
# bounds what an overflow would otherwise do.
_LOG_RATE_CLIP = 100.0
# Where a nonlinear constant is sampled when the question is whether the
# DESIGN can identify it: if no log10 K on this grid varies the regressor, no
# K can. `_term_terms` takes log10 K, as `law_design` does.
_LOG10_K_GRID = np.linspace(-4.0, 4.0, 33)

# Where `search` saves. Never overwrites an existing save. `v2` is Amendment
# 1's rerun under the two-test tie rule, the run-clustered errors, the
# four-temperature floor and the between-run families; the 2026-09-14 saves
# sit beside it in `rate_laws/` and are not touched.
RATE_LAW_DIR = os.path.join("data", "fits", "rate_laws", "v2")


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


# --- the rate-law fitter ----------------------------------------------------


def _nonlinear_names(family, option):
    """The nonlinear constants an option carries, if any."""
    name = NONLINEAR_OF.get((family, option))
    return [name] if name else []


def _nonlinear_names_for_model(model):
    """Every nonlinear constant a model carries, sorted."""
    return sorted({name for family, option in model.items()
                   for name in _nonlinear_names(family, option)})


def _term_terms(table, family, option, nonlinear):
    """
    [(name, values, linear, family)] for one family's option.

    A LINEAR term gets a fitted coefficient. A term from a nonlinear option
    has its coefficient fixed at one -- `log([S]/(Km + [S]))` is a function of
    Km alone -- so its values become part of the offset the linear solve works
    against. Both are returned, because both are regressors a health check has
    to look at.
    """
    if family == "S":
        s0 = table.s0.to_numpy(dtype=float)
        if option == "power":
            return [("a_S", np.log(s0), True, "S")]
        if option == "mm":
            km = 10.0 ** nonlinear["Km"]
            return [("S:mm", np.log(s0 / (km + s0)), False, "S")]
    elif family == "H":
        h2o2 = table.h2o2.to_numpy(dtype=float)
        if option == "power":
            return [("a_H", np.log(h2o2), True, "H")]
        if option == "bind":
            k_h = 10.0 ** nonlinear["K_H"]
            return [("H:bind", np.log(k_h * h2o2 / (1.0 + k_h * h2o2)),
                     False, "H")]
        if option == "relax":
            k_h = 10.0 ** nonlinear["K_H"]
            return [("H:relax", np.log1p(k_h * h2o2), False, "H")]
    elif family == "HOO":
        return [("a_HOO", np.log(table.hoo.to_numpy(dtype=float)), True, "HOO")]
    elif family == "BUF":
        buf = table.buf.to_numpy(dtype=float)
        if option == "power":
            return [("a_B", np.log(buf), True, "BUF")]
        if option == "bind":
            k_b = 10.0 ** nonlinear["K_B"]
            return [("BUF:bind", np.log(k_b * buf / (1.0 + k_b * buf)),
                     False, "BUF")]
        if option == "relax":
            k_b = 10.0 ** nonlinear["K_B"]
            return [("BUF:relax", np.log1p(k_b * buf), False, "BUF")]
        if option == "power_by_buffer":
            # One slope per buffer that HAS a change of [buf]: a buffer whose
            # concentration never moves inside the table cannot carry a slope
            # distinct from its own intercept, and emitting a zero column for
            # it would make the design rank deficient.
            terms = []
            for buffer in sorted(set(table.buffer)):
                rows = table[table.buffer == buffer]
                if rows.buf.nunique() < 2:
                    continue
                mask = (table.buffer == buffer).to_numpy()
                terms.append((f"a_B[{buffer}]",
                              np.where(mask, np.log(buf), 0.0), True, "BUF"))
            return terms
    elif family == "E":
        return [("a_E", np.log(table.e0.to_numpy(dtype=float)), True, "E")]
    elif family == "T":
        inverse = 1.0 / table.kelvin.to_numpy(dtype=float)
        return [("Ea_R", -(inverse - 1.0 / 298.15), True, "T")]
    raise ValueError(f"no term for {family}:{option}")


def law_design(table, model, nonlinear):
    """
    The weighted design for one model at fixed nonlinear constants (log10).

    Named `law_design` because `scope.design` already holds `design`; the
    plan's `design(table, model, nonlinear)`. Returns:

      matrix    (n, p) linear columns: one intercept per buffer present, then
                one column per linear rate-law term, in `FAMILIES` order;
      names     the column names: `intercept[<buffer>]`, `a_S`, `a_HOO`, ...,
                `a_B[<buffer>]`, `a_E`, `Ea_R` (Ea/R, K);
      offset    the sum of the nonlinear options' terms, a quantity with no
                fitted coefficient to subtract from `y` before the solve;
      terms     (name, values, linear, family) for every non-intercept
                regressor, linear or offset, for the health checks.
    """
    buffers = sorted(set(table.buffer))
    columns = [np.where((table.buffer == buffer).to_numpy(), 1.0, 0.0)
               for buffer in buffers]
    names = [f"intercept[{buffer}]" for buffer in buffers]
    offset = np.zeros(len(table))
    terms = []
    for family in FAMILIES:
        option = model.get(family)
        if option is None:
            continue
        for name, values, linear, term_family in _term_terms(table, family,
                                                             option, nonlinear):
            terms.append((name, values, linear, term_family))
            if linear:
                columns.append(values)
                names.append(name)
            else:
                offset = offset + values
    return np.column_stack(columns), names, offset, terms


def _buffer_demeaned(values, buffers):
    """`values` with each buffer's own mean removed."""
    values = np.asarray(values, dtype=float)
    demeaned = values.copy()
    for buffer in set(buffers):
        mask = buffers == buffer
        demeaned[mask] -= values[mask].mean()
    return demeaned


def _run_demeaned(values, runs):
    """`values` with each run's own mean removed."""
    values = np.asarray(values, dtype=float)
    demeaned = values.copy()
    for run in set(runs):
        mask = runs == run
        demeaned[mask] -= values[mask].mean()
    return demeaned


def _variation_ok(values, family):
    """
    Whether a regressor varies enough to be identified: at least two distinct
    values, and an SD above the family's floor. `values` is expected already
    demeaned by whichever offsets the question carries.
    """
    values = np.asarray(values, dtype=float)
    if np.unique(values).size < 2:
        return False, "fewer than 2 distinct values"
    tolerance = _T_REGRESSOR_BELOW if family == "T" else _LOG_REGRESSOR_BELOW
    if float(np.std(values)) < tolerance:
        return False, "SD after the offsets are removed below the floor"
    return True, ""


def term_identifiable(table, family, option):
    """
    Whether this table can identify a term at all, and if not, why.

    Three rules, checked in order: the regressor must take at least two
    distinct values; its SD after removing each buffer's mean must clear 0.05
    in log units (T: 1e-5 in 1/K); and `power_by_buffer` needs at least two
    buffers that each carry two runs and a change of [buf].

    A nonlinear option is evaluated over a grid of log10 K in its own bounds,
    and passes if ANY K gives a regressor the table can identify -- the design
    is a property of the columns, and K is a free direction, so a fixed
    K = 1/mM would reject a table whose [buf] is far from it. The reason
    reported is the one from the last grid point tried.

    `T:arrhenius` carries one rule on top of the variation tests: the table
    must move at least 4 distinct temperatures. It is checked only after those
    tests pass, so a table that never moved temperature keeps its
    "fewer than 2 distinct values" reason and the temperature count is the
    reason only where the axis exists but is too short.
    """
    buffers = table.buffer.to_numpy()
    if option == "power_by_buffer":
        qualifying = 0
        for buffer in sorted(set(table.buffer)):
            rows = table[table.buffer == buffer]
            if rows.experiment.nunique() >= 2 and rows.buf.nunique() >= 2:
                qualifying += 1
        if qualifying < 2:
            return False, "fewer than 2 buffers with within-buffer variation"
    names_nonlinear = _nonlinear_names(family, option)
    if names_nonlinear:
        grids = [{name: value for name in names_nonlinear}
                 for value in _LOG10_K_GRID]
    else:
        grids = [{}]
    reason = "fewer than 2 distinct values"
    for nonlinear in grids:
        ok = True
        for _, values, _, term_family in _term_terms(table, family, option,
                                                     nonlinear):
            ok, reason = _variation_ok(_buffer_demeaned(values, buffers),
                                       term_family)
            if not ok:
                break
        if ok:
            if family == "T" and table.kelvin.nunique() < 4:
                return False, "fewer than 4 temperatures"
            return True, ""
    return False, reason


def _law_solve(table, model, nonlinear, y, weight):
    """
    The linear solve at fixed nonlinear constants: weighted lstsq of
    `y - offset` on the design, and its weighted cost.
    """
    matrix, names, offset, terms = law_design(table, model, nonlinear)
    root = np.sqrt(weight)
    target = y - offset
    beta, *_ = np.linalg.lstsq(matrix * root[:, None], target * root,
                               rcond=None)
    residual = target - matrix @ beta
    return matrix, names, offset, terms, beta, residual, float(
        np.sum(weight * residual ** 2))


def _clustered_errors(matrix, weight, error, runs):
    """
    Run-clustered standard errors of the linear coefficients.

    `B = pinv(X' diag(w) X)`, `M = sum over runs of s_g s_g'` with
    `s_g = X_g' (w_g o e_g)`, `V = G/(G - 1) B M B`, and the errors are the
    square roots of V's diagonal. The curves of one run share a day, a stock
    and a cuvette offset, so treating them as independent understates the
    error; `fit_model` keeps both and everything quoted uses this one.
    """
    bread = np.linalg.pinv(matrix.T @ (matrix * weight[:, None]))
    meat = np.zeros_like(bread)
    groups = sorted(set(runs))
    for run in groups:
        mask = runs == run
        score = matrix[mask].T @ (weight[mask] * error[mask])
        meat += np.outer(score, score)
    if len(groups) < 2:
        return np.full(matrix.shape[1], np.nan)
    covariance = len(groups) / (len(groups) - 1.0) * bread @ meat @ bread
    return np.sqrt(np.maximum(np.diag(covariance), 0.0))


def fit_model(table, model, starts=(-2.0, 0.0, 2.0)):
    """
    Weighted least squares for one model on one element table.

    The linear coefficients -- intercepts, orders, `Ea_R` -- are solved at
    fixed nonlinear constants; each nonlinear constant is optimised in log10
    by `scipy.optimize.least_squares` within (-4, 4). Every combination of
    `starts` is tried and the lowest weighted cost kept.

    Returns `model`, `coefficients` (name -> value, with nonlinear constants
    under both `log10(K)` and `K`, and `Ea` in kJ/mol beside `Ea_R` in K),
    `stderr` (the naive errors) and `stderr_clustered` (run-clustered, the one
    every comparison and quoted +/- uses; NaN on the nonlinear constants and
    on `Ea` absent), `cost` (weighted sum of squared residuals), `curves`,
    `runs`, `nonlinear` (log10 values), `nonlinear_at_bound`, `correlation` (a
    DataFrame over the linear coefficients) and `health` = {"flags", "ok"}.
    The underscore keys carry the design and the regressors the health checks
    read; they are not part of the result's contract.
    """
    y = table.y.to_numpy(dtype=float)
    weight = table.weight.to_numpy(dtype=float)
    names_nonlinear = sorted({name for family, option in model.items()
                              for name in _nonlinear_names(family, option)})
    starts = tuple(starts)
    best = None
    for start in itertools.product(starts, repeat=len(names_nonlinear)):
        if names_nonlinear:
            def objective(x):
                nonlinear = dict(zip(names_nonlinear, np.atleast_1d(x)))
                _, _, _, _, _, residual, _ = _law_solve(
                    table, model, nonlinear, y, weight)
                return np.sqrt(weight) * residual

            solution = least_squares(objective, np.asarray(start, dtype=float),
                                     bounds=NONLINEAR_BOUNDS)
            values = solution.x
            jacobian = solution.jac
        else:
            values, jacobian = np.array([]), None
        nonlinear = dict(zip(names_nonlinear, np.atleast_1d(values)))
        matrix, names, _, terms, beta, residual, cost = _law_solve(
            table, model, nonlinear, y, weight)
        if best is None or cost < best["cost"]:
            best = {"nonlinear": nonlinear, "matrix": matrix, "names": names,
                    "terms": terms, "beta": beta, "cost": cost,
                    "jacobian": jacobian}

    matrix, names, beta = best["matrix"], best["names"], best["beta"]
    cost, count = best["cost"], len(table)
    parameters = len(names) + len(names_nonlinear)
    degrees = max(1, count - parameters)
    variance = cost / degrees
    covariance = np.linalg.pinv(
        matrix.T @ (matrix * weight[:, None])) * variance
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    safe = np.where(errors > 0, errors, np.inf)
    correlation = covariance / np.outer(safe, safe)
    coefficients = {name: float(value) for name, value in zip(names, beta)}
    stderr = {name: float(value) for name, value in zip(names, errors)}
    for name, value in best["nonlinear"].items():
        coefficients[f"log10({name})"] = float(value)
        coefficients[name] = float(10.0 ** value)
        if best["jacobian"] is not None:
            nonlinear_covariance = np.linalg.pinv(
                best["jacobian"].T @ best["jacobian"]) * variance
            error = float(np.sqrt(max(np.diag(nonlinear_covariance)[
                names_nonlinear.index(name)], 0.0)))
        else:
            error = np.nan
        stderr[f"log10({name})"] = error
        stderr[name] = error * float(np.log(10.0)) * 10.0 ** value
    if "Ea_R" in coefficients:
        coefficients["Ea"] = coefficients["Ea_R"] * _GAS_CONSTANT / 1000.0
        stderr["Ea"] = stderr["Ea_R"] * _GAS_CONSTANT / 1000.0

    result = {
        "model": dict(model),
        "coefficients": coefficients,
        "stderr": stderr,
        "cost": float(cost),
        "curves": count,
        "runs": int(table.experiment.nunique()),
        "nonlinear": {name: float(value)
                      for name, value in best["nonlinear"].items()},
        "nonlinear_at_bound": [
            name for name, value in best["nonlinear"].items()
            if (abs(value - NONLINEAR_BOUNDS[0]) <= 0.05
                or abs(value - NONLINEAR_BOUNDS[1]) <= 0.05)],
        "correlation": pd.DataFrame(correlation, index=names, columns=names),
        "_matrix": matrix,
        "_names": names,
        "_terms": best["terms"],
    }
    _, prediction = _predict(table, model, result)
    clustered = _clustered_errors(matrix, weight, y - prediction,
                                  table.experiment.to_numpy(dtype=int))
    stderr_clustered = {name: float(value)
                        for name, value in zip(names, clustered)}
    if "Ea_R" in coefficients:
        stderr_clustered["Ea"] = (stderr_clustered["Ea_R"]
                                  * _GAS_CONSTANT / 1000.0)
    for name in best["nonlinear"]:
        stderr_clustered[f"log10({name})"] = np.nan
        stderr_clustered[name] = np.nan
    result["stderr_clustered"] = stderr_clustered
    flags, ok = rate_law_health(result, table)
    result["health"] = {"flags": flags, "ok": ok}
    return result


def rate_law_health(result, table):
    """
    The flags a fit must be read with, and whether it has none.

    Exactly the plan's list: "at bound: <K>" for a nonlinear constant within
    0.05 of either bound; "collinear: <a>/<b> r=+/-0.xx" for a coefficient
    pair (intercepts excluded) above |r| 0.95; "rank deficient"; "carried by
    fewer than 3 runs: <term>" where a term's regressor, after removing each
    buffer's mean, is non-zero in fewer than three runs; "too few curves" for
    fewer than 15; and "fewer than 10 runs: clustered errors unreliable" where
    the run-clustered errors have too few clusters to trust. Returns
    (flags, ok).
    """
    flags = []
    for name, value in result["nonlinear"].items():
        if (abs(value - NONLINEAR_BOUNDS[0]) <= 0.05
                or abs(value - NONLINEAR_BOUNDS[1]) <= 0.05):
            flags.append(f"at bound: {name}")
    if np.linalg.matrix_rank(result["_matrix"]) < result["_matrix"].shape[1]:
        flags.append("rank deficient")
    correlation = result["correlation"]
    columns = list(correlation.columns)
    for i, first in enumerate(columns):
        if first.startswith("intercept"):
            continue
        for second in columns[i + 1:]:
            if second.startswith("intercept"):
                continue
            r = float(correlation.loc[first, second])
            if np.isfinite(r) and abs(r) > _CORRELATION_BAR:
                flags.append(f"collinear: {first}/{second} r={r:+.2f}")
    buffers = table.buffer.to_numpy()
    experiments = table.experiment.to_numpy()
    for name, values, _, _ in result["_terms"]:
        demeaned = _buffer_demeaned(values, buffers)
        carried = {int(run) for run, value in zip(experiments, demeaned)
                   if abs(value) > 1e-9}
        if len(carried) < 3:
            flags.append(f"carried by fewer than 3 runs: {name}")
    if result["curves"] < 15:
        flags.append("too few curves")
    if result["runs"] < 10:
        flags.append("fewer than 10 runs: clustered errors unreliable")
    return flags, not flags


def _run_design(table, model, nonlinear, kept):
    """
    The design with one intercept per RUN instead of per buffer. Nonlinear
    terms among `kept` keep their offset form; every other offset family is
    dropped because its regressor is constant inside a run and the run
    intercept would absorb it.
    """
    runs = sorted(table.experiment.unique())
    columns = [np.where((table.experiment == run).to_numpy(), 1.0, 0.0)
               for run in runs]
    names = [f"run[{run}]" for run in runs]
    offset = np.zeros(len(table))
    _, _, _, terms = law_design(table, model, nonlinear)
    for name, values, linear, _ in terms:
        if name not in kept:
            continue
        if linear:
            columns.append(values)
            names.append(name)
        else:
            offset = offset + values
    return np.column_stack(columns), names, offset


def within_run_check(table, model):
    """
    Refit a model with one intercept per run, and compare the coefficients.

    A regressor constant inside every run -- pH, [enz], temperature, and every
    buffer identity as a family -- is collinear with a run intercept, so those
    terms are dropped after a test on the run-demeaned regressor. HOO, E and T
    are ALWAYS dropped, whatever their within-run variation: [HOO-] moves
    inside a run only through ionic strength, which follows [buf], so a
    within-run [HOO-] coefficient is a buffer coefficient, and [enz] and
    temperature are one value per run by design. The shared coefficients are
    then read twice: once with the buffer intercepts and once with the run
    intercepts, and both errors are run-clustered. `disagree` is
    |a - b| > 2 sqrt(se_a^2 + se_b^2), the plan's bar.

    Returns `table` (one row per shared coefficient: both estimates, both
    standard errors, `disagree`), `dropped` (each term the run design could not
    carry, with its reason, `"between-run family"` for HOO, E and T) and
    `s_buf_within_run_r`, the within-run correlation of the S and BUF
    regressors where both families are in the model.
    """
    base = fit_model(table, model)
    _, _, _, terms = law_design(table, model, base["nonlinear"])
    runs = table.experiment.to_numpy()
    kept, dropped = [], []
    for name, values, _, family in terms:
        if family in ("HOO", "E", "T"):
            dropped.append((name, "between-run family"))
            continue
        ok, reason = _variation_ok(_run_demeaned(values, runs), family)
        if ok:
            kept.append(name)
        else:
            dropped.append((name, reason))

    matrix, names, offset = _run_design(table, model, base["nonlinear"], kept)
    y = table.y.to_numpy(dtype=float)
    weight = table.weight.to_numpy(dtype=float)
    root = np.sqrt(weight)
    target = y - offset
    beta, *_ = np.linalg.lstsq(matrix * root[:, None], target * root,
                               rcond=None)
    residual = target - matrix @ beta
    errors = dict(zip(names, _clustered_errors(matrix, weight, residual, runs)))
    values = dict(zip(names, beta))

    # Nonlinear terms carry no coefficient of their own -- their constant is
    # the same in both fits -- so only the fitted coefficients are compared.
    shared = [name for name, _, linear, _ in terms
              if linear and name in kept]
    rows = []
    for name in shared:
        a, b = base["coefficients"].get(name), values.get(name)
        se_a = base["stderr_clustered"].get(name)
        se_b = errors.get(name)
        if None in (a, b) or not all(np.isfinite([a, b, se_a, se_b])):
            disagree = True
        else:
            disagree = bool(abs(a - b) > 2.0 * np.sqrt(se_a ** 2 + se_b ** 2))
        rows.append({"coefficient": name, "between_runs": a, "within_runs": b,
                     "stderr_between": se_a, "stderr_within": se_b,
                     "disagree": disagree})
    s_buf = None
    if model.get("S") and model.get("BUF"):
        first = {}
        for _, term_values, _, family in terms:
            first.setdefault(family, term_values)
        s_demeaned = _run_demeaned(first["S"], runs)
        b_demeaned = _run_demeaned(first["BUF"], runs)
        if np.std(s_demeaned) > 0 and np.std(b_demeaned) > 0:
            s_buf = float(np.corrcoef(s_demeaned, b_demeaned)[0, 1])
    return {"table": pd.DataFrame(
                rows, columns=["coefficient", "between_runs", "within_runs",
                               "stderr_between", "stderr_within", "disagree"]
            ).set_index("coefficient"),
            "dropped": tuple(dropped), "s_buf_within_run_r": s_buf}


# --- the rate-law search ----------------------------------------------------


def _model_id(model):
    """A stable string for a model: the non-None families in `FAMILIES` order."""
    return "|".join(f"{family}:{model[family]}" for family in FAMILIES
                    if model.get(family))


def enumerate_models(table, element):
    """
    Every model allowed by `TERM_OPTIONS` and `term_identifiable`.

    At most one option per family, `None` included as a choice. Returns the
    models and a `dropped` table of (family, option, reason) for every option
    this table cannot identify.
    """
    choices, dropped = {}, []
    for family in FAMILIES:
        allowed = [None]
        for option in TERM_OPTIONS[element][family]:
            ok, reason = term_identifiable(table, family, option)
            if ok:
                allowed.append(option)
            else:
                dropped.append({"family": family, "option": option,
                                "reason": reason})
        choices[family] = allowed
    models = [dict(zip(FAMILIES, combination))
              for combination in itertools.product(
                  *(choices[family] for family in FAMILIES))]
    return models, pd.DataFrame(dropped)


def _predict(table, model, fit):
    """(y, prediction) for `table` under a fit made elsewhere."""
    if len(table) == 0:
        return np.zeros(0), np.zeros(0)
    matrix, names, offset, _ = law_design(table, model, fit["nonlinear"])
    beta = np.array([fit["coefficients"][name] for name in names])
    return table.y.to_numpy(dtype=float), matrix @ beta + offset


def cross_validate(table, model):
    """
    Leave-one-run-out: fit on the other runs, predict the held-out one.

    Fold score is the held-out curves' weighted squared error
    `sum weight * (y - y_hat)^2`. A fold whose held-out run carries a buffer
    no other run in the table has is skipped and counted -- its intercept
    would be unseen. Returns `scores` (a Series over runs), `sum`, `skipped`
    and `trained` (the experiments each fold was trained on, for the test that
    no curve is ever trained on the run it is scored against).
    """
    scores, trained = {}, {}
    skipped, total = 0, 0.0
    for run in sorted(table.experiment.unique()):
        others = table[table.experiment != run]
        held = table[table.experiment == run]
        if not set(held.buffer) <= set(others.buffer):
            skipped += 1
            continue
        fit = fit_model(others, model)
        y, prediction = _predict(held, model, fit)
        score = float(np.sum(held.weight.to_numpy(dtype=float)
                             * (y - prediction) ** 2))
        scores[int(run)] = score
        trained[int(run)] = tuple(int(e) for e in others.experiment.unique())
        total += score
    return {"scores": pd.Series(scores, dtype=float).sort_index(),
            "sum": float(total), "skipped": skipped, "trained": trained}


def _cross_validate_task(payload):
    """`cross_validate` for one (table, model), for the process pool."""
    identifier, table, model = payload
    return identifier, cross_validate(table, model)


def tie_statistics(scores):
    """
    Both tie tests, per model, over the folds every model scored.

    `scores` is a DataFrame, one row per model, one column per fold. The folds
    are the columns with no NaN in any model; the best is the lowest total.
    The raw test is `d = score - score[best]`, kept if
    `mean(d) <= 2 sd(d) / sqrt(folds)`. The log test is the same inequality on
    `d = log(max(score, TIE_LOG_FLOOR)) - log(max(score[best], TIE_LOG_FLOOR))`.
    A model is tied only if BOTH keep it; the best is always tied. Returns one
    row per model in `scores.index` order, with `raw_mean`, `raw_bar`,
    `raw_kept`, `log_mean`, `log_bar`, `log_kept`, `tied` and
    `worst_fold_log_ratio` -- the largest `log(score/score[best])` over the
    folds, a reported diagnostic that changes no tie set.

    The two tests agree to within their shared bar, and both use the same
    spread as the differences they judge: a model that fails on `k` of `n`
    folds has paired `t = sqrt(k(n-1)/(n-k))` whatever the SIZE of the failure,
    so the rule cannot exclude one that fails on `k <= 4n/(n+3)` folds (3 of
    30, or 2 of 18) however badly. `worst_fold_log_ratio` is how a report
    shows where that bites.
    """
    scores = scores.dropna(axis=1, how="any")
    totals = scores.sum(axis=1)
    best = totals.idxmin()
    folds = scores.shape[1]
    raw = scores.subtract(scores.loc[best], axis=1)
    raw_mean = raw.mean(axis=1)
    raw_sd = raw.std(axis=1, ddof=1) if folds > 1 else raw_mean * 0.0
    raw_bar = 2.0 * raw_sd / np.sqrt(max(folds, 1))
    logged = np.log(scores.clip(lower=TIE_LOG_FLOOR))
    log = logged.subtract(logged.loc[best], axis=1)
    log_mean = log.mean(axis=1)
    log_sd = log.std(axis=1, ddof=1) if folds > 1 else log_mean * 0.0
    log_bar = 2.0 * log_sd / np.sqrt(max(folds, 1))
    statistics = pd.DataFrame(
        {"raw_mean": raw_mean, "raw_bar": raw_bar,
         "raw_kept": raw_mean <= raw_bar,
         "log_mean": log_mean, "log_bar": log_bar,
         "log_kept": log_mean <= log_bar,
         "worst_fold_log_ratio": log.max(axis=1)},
        index=scores.index)
    statistics["tied"] = statistics.raw_kept & statistics.log_kept
    statistics.loc[best, "tied"] = True
    return statistics


def tie_set(scores):
    """
    The models kept by BOTH tie tests, in `scores.index` order.
    """
    statistics = tie_statistics(scores)
    return list(statistics.index[statistics.tied])


def term_verdicts(tie, dropped, element):
    """
    The plan's verdict strings, per family, for a tie set of models.

    `tie` is the tied models (dicts), `dropped` the (family, option, reason)
    table from `enumerate_models`, `element` names which `TERM_OPTIONS` the
    options come from. A family whose every configured option was dropped is
    "not identifiable on this table (<reason>)".
    """
    verdicts = {}
    for family in FAMILIES:
        configured = set(TERM_OPTIONS[element][family])
        dropped_options = set(dropped[dropped.family == family].option) \
            if len(dropped) else set()
        if configured <= dropped_options:
            reasons = dropped[dropped.family == family].reason.tolist()
            verdicts[family] = (f"{family}: not identifiable on this table "
                                f"({reasons[0]})")
            continue
        options = [model.get(family) for model in tie]
        used = sorted({option for option in options if option is not None})
        if not used:
            verdicts[family] = f"{family}: absent from every tied model"
        elif len(used) == 1 and all(option == used[0] for option in options):
            verdicts[family] = f"{family}: {used[0]} in every tied model"
        elif all(option is not None for option in options):
            verdicts[family] = (f"{family}: a dependence in every tied model, "
                                f"option undecided ({', '.join(used)})")
        else:
            verdicts[family] = f"{family}: undecided"
    return verdicts


def _apply_cut(table, cut):
    """The plan's sensitivity cuts, exactly."""
    if cut in (None, "all"):
        return table
    if cut == "S-gas":
        return table[table.bubble_load <= 1.0]
    if cut == "S-pyro":
        return table[table.buffer != "Pyrophosphate"]
    if cut == "S-weak":
        weak = set(fit_dataset.TWO_AXIS_BLOCK) - set(scope.strong_runs())
        return table[~table.experiment.isin(weak)]
    raise ValueError(f"unknown cut {cut!r}")


def _search_table(table, element, workers=1):
    """
    The search's core: enumerate, cross-validate, tie, verdict, health.

    No cuts and no saving; `search` wraps this with both. Split out so the
    planted-identifiability test can run the same machinery on planted tables.
    """
    models, dropped = enumerate_models(table, element)
    identifiers = [_model_id(model) for model in models]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = dict(pool.map(
                _cross_validate_task,
                ((identifier, table, model)
                 for identifier, model in zip(identifiers, models))))
    else:
        results = {identifier: cross_validate(table, model)
                   for identifier, model in zip(identifiers, models)}
    scores = pd.DataFrame({identifier: result["scores"]
                           for identifier, result in results.items()}).T
    statistics = tie_statistics(scores)
    tied_identifiers = [identifier for identifier in identifiers
                        if bool(statistics.loc[identifier, "tied"])]
    by_identifier = dict(zip(identifiers, models))
    tied_models = [by_identifier[identifier] for identifier in tied_identifiers]
    health, disagreements = {}, {}
    for identifier in tied_identifiers:
        fit = fit_model(table, by_identifier[identifier])
        health[identifier] = fit["health"]["flags"]
        check = within_run_check(table, by_identifier[identifier])
        disagreements[identifier] = {
            name: bool(value) for name, value
            in check["table"]["disagree"].to_dict().items()}
    skipped = max((result["skipped"] for result in results.values()),
                  default=0)
    report = {
        "curves": int(len(table)),
        "runs": int(table.experiment.nunique()),
        "folds_skipped": int(skipped),
        "models": int(len(models)),
        "scores": {identifier: float(result["sum"])
                   for identifier, result in results.items()},
        "tie": list(tied_identifiers),
        "fold_scores": {
            identifier: {str(run): float(score) for run, score
                         in results[identifier]["scores"].items()}
            for identifier in identifiers},
        "tie_statistics": [
            {"model": identifier,
             "raw_mean": float(statistics.loc[identifier, "raw_mean"]),
             "raw_bar": float(statistics.loc[identifier, "raw_bar"]),
             "raw_kept": bool(statistics.loc[identifier, "raw_kept"]),
             "log_mean": float(statistics.loc[identifier, "log_mean"]),
             "log_bar": float(statistics.loc[identifier, "log_bar"]),
             "log_kept": bool(statistics.loc[identifier, "log_kept"]),
             "worst_fold_log_ratio": float(
                 statistics.loc[identifier, "worst_fold_log_ratio"]),
             "tied": bool(statistics.loc[identifier, "tied"])}
            for identifier in identifiers],
        "dropped": dropped.to_dict("records"),
        "verdicts": term_verdicts(tied_models, dropped, element),
        "health": health,
        "within_run": disagreements,
    }
    return report


def search(substrate, element, cut=None, workers=8):
    """
    Search every allowed model on one element table and save the report.

    Saves `data/fits/rate_laws/v2/<substrate>_<element>_<cut or all>.json`,
    with every model's fold scores and both tie tests; a save that already
    exists is read back and not recomputed. An element with fewer than 15
    curves after the cut is refused (the plan's floor).
    """
    os.makedirs(RATE_LAW_DIR, exist_ok=True)
    path = os.path.join(RATE_LAW_DIR,
                        f"{substrate}_{element}_{cut or 'all'}.json")
    if os.path.exists(path):
        with open(path) as handle:
            return json.load(handle)
    table = _apply_cut(curve_parameters(substrate)["elements"][element], cut)
    if len(table) < 15:
        raise ValueError(f"{substrate} {element} {cut}: too few curves "
                         f"({len(table)})")
    report = _search_table(table, element, workers=workers)
    report.update({"substrate": substrate, "element": element, "cut": cut})
    with open(path, "w") as handle:
        json.dump(report, handle, default=float, indent=1, sort_keys=True)
    return report


# --- links between the elements (Task 5) ------------------------------------


def _parse_model_id(identifier):
    """The model a saved tie identifier names."""
    model = {family: None for family in FAMILIES}
    if identifier:
        for part in identifier.split("|"):
            family, option = part.split(":", 1)
            model[family] = option
    return model


def _best_tied_model(substrate, element, family=None, option=None):
    """The lowest-CV model in a table's saved tie set, optionally with one
    family's option forced (for a link that adds the term to the law)."""
    report = search(substrate, element, workers=1)
    if not report["tie"]:
        return None
    best = min(report["tie"], key=lambda identifier: report["scores"][identifier])
    model = _parse_model_id(best)
    if family is not None:
        model[family] = option
    return model


def _other_nonlinear(table_model, name):
    """The nonlinear constants a model carries besides `name`."""
    return sorted({constant for family, option in table_model.items()
                   for constant in _nonlinear_names(family, option)
                   if constant != name})


def _binding_layout(models, name, shared):
    """[(element, local constant, global index)] for a joint binding fit."""
    layout, index = [], 0
    if shared:
        for element in range(len(models)):
            layout.append((element, name, 0))
        index = 1
    else:
        for element in range(len(models)):
            layout.append((element, name, index))
            index += 1
    for element, model in enumerate(models):
        for nuisance in _other_nonlinear(model, name):
            layout.append((element, nuisance, index))
            index += 1
    return layout


def _binding_nonlinear(layout, x, element):
    """The local nonlinear dict for one element at a global vector x."""
    return {local: float(x[global_index])
            for index, local, global_index in layout if index == element}


def _binding_residuals(x, tables, models, layout):
    """Stacked weighted residuals at a global nonlinear vector x."""
    blocks = []
    for element, (table, model) in enumerate(zip(tables, models)):
        nonlinear = _binding_nonlinear(layout, x, element)
        matrix, _, offset, _ = law_design(table, model, nonlinear)
        y = table.y.to_numpy(dtype=float)
        weight = table.weight.to_numpy(dtype=float)
        root = np.sqrt(weight)
        target = y - offset
        beta, *_ = np.linalg.lstsq(matrix * root[:, None], target * root,
                                   rcond=None)
        blocks.append((target - matrix @ beta) * root)
    return np.concatenate(blocks)


def _binding_fit(tables, models, name, shared, starts=(-2.0, 0.0, 2.0)):
    """Fit two element tables jointly with one binding constant or one each.

    Returns `fits` (one `_predict`-compatible dict per element), the global
    vector `x` (log10) and the summed weighted cost.
    """
    layout = _binding_layout(models, name, shared)
    best_x, best_cost = None, np.inf
    for start in itertools.product(starts, repeat=len(layout)):
        solution = least_squares(
            _binding_residuals, np.asarray(start, dtype=float),
            bounds=NONLINEAR_BOUNDS,
            args=(tables, models, layout))
        cost = float(2.0 * solution.cost)
        if cost < best_cost:
            best_x, best_cost = solution.x, cost
    fits = []
    for element, (table, model) in enumerate(zip(tables, models)):
        nonlinear = _binding_nonlinear(layout, best_x, element)
        matrix, names, offset, _ = law_design(table, model, nonlinear)
        y = table.y.to_numpy(dtype=float)
        weight = table.weight.to_numpy(dtype=float)
        root = np.sqrt(weight)
        target = y - offset
        beta, *_ = np.linalg.lstsq(matrix * root[:, None], target * root,
                                   rcond=None)
        fits.append({"model": model, "nonlinear": nonlinear,
                     "coefficients": {name_i: float(value)
                                      for name_i, value in zip(names, beta)}})
    return {"fits": fits, "x": np.asarray(best_x, dtype=float),
            "cost": best_cost, "layout": layout}


def _binding_cross_validate(tables, models, name, shared,
                            starts=(-2.0, 0.0, 2.0)):
    """Leave-one-run-out over the runs present in either element table."""
    runs = sorted(set().union(*(set(table.experiment) for table in tables)))
    scores, skipped = {}, 0
    for run in runs:
        train = [table[table.experiment != run] for table in tables]
        held = [table[table.experiment == run] for table in tables]
        if any(not set(h.buffer) <= set(t.buffer)
               for h, t in zip(held, train)):
            skipped += 1
            continue
        fit = _binding_fit(train, models, name, shared, starts)
        score = 0.0
        for element, (h, model) in enumerate(zip(held, models)):
            y, prediction = _predict(h, model, fit["fits"][element])
            score += float(np.sum(h.weight.to_numpy(dtype=float)
                                  * (y - prediction) ** 2))
        scores[int(run)] = score
    return scores, skipped


def _carrying_buffers(table, family):
    """
    The buffers whose own rows move the family's concentration.

    A buffer carries `BUF` or `H` if it has at least 2 curves and the log of
    the family's concentration ([buf] or [H2O2]), after removing that buffer's
    own mean, has SD 0.05 or more. A constant identified inside one buffer's
    runs cannot be separated from that buffer, so a link needs two carriers.
    """
    column = {"BUF": "buf", "H": "h2o2"}[family]
    carrying = []
    for buffer in sorted(set(table.buffer)):
        rows = table[table.buffer == buffer]
        values = np.log(rows[column].to_numpy(dtype=float))
        if len(rows) >= 2 and float(np.std(values - values.mean())) >= 0.05:
            carrying.append(buffer)
    return carrying


def link_power(substrate, family, k_rate, k_clock, seeds):
    """
    How often the link detects two different Ks at realistic noise.

    For each seed, plants `v_act` as `log(k_rate x/(1 + k_rate x))` and
    `k_act_lag` as `log(1 + k_clock x)` on the real element tables (x = [buf]
    for BUF, [H2O2] for H), with the real rows' own `sqrt(se^2 + SE_FLOOR^2)`
    as the noise SD, and runs `shared_binding_law` on the planted pair. Returns
    per seed the verdict, the shared and separate Ks, and the raw and log paired
    t of the shared-minus-separate fold scores, plus `caught`, the number of
    seeds whose verdict is `"separate K predicts better"`.

    Saves `data/fits/rate_laws/v2/link_power.json`, one entry per
    (substrate, family, k_rate, k_clock, seeds), and reads it back if it
    already holds that entry.
    """
    key = (f"{substrate}|{family}|{k_rate}|{k_clock}|"
           + ",".join(str(seed) for seed in seeds))
    path = os.path.join(RATE_LAW_DIR, "link_power.json")
    saved = {}
    if os.path.exists(path):
        with open(path) as handle:
            saved = json.load(handle)
    if key in saved:
        return saved[key]
    parameters = curve_parameters(substrate)
    rate_table = parameters["elements"]["v_act"]
    clock_table = parameters["elements"]["k_act_lag"]
    intercept = {"Boric": 0.0, "Phosphate": -0.5, "Pyrophosphate": 0.5}
    column = {"BUF": "buf", "H": "h2o2"}[family]
    detail = {}
    for seed in seeds:
        generator = np.random.default_rng(seed)
        planted = {}
        for element, table, constant in (("v_act", rate_table, k_rate),
                                         ("k_act_lag", clock_table, k_clock)):
            x = table[column].to_numpy(dtype=float)
            spread = np.sqrt(table.se.to_numpy(dtype=float) ** 2
                             + SE_FLOOR ** 2)
            trend = (np.log(constant * x / (1.0 + constant * x))
                     if element == "v_act" else np.log1p(constant * x))
            copy = table.copy()
            copy["y"] = (np.array([intercept[buffer] for buffer in copy.buffer])
                         + trend + generator.normal(0.0, spread))
            planted[element] = copy
        found = shared_binding_law(substrate, family,
                                   rate_table=planted["v_act"],
                                   clock_table=planted["k_act_lag"])
        row = {"verdict": found["verdict"], "shared_k": found["shared_k"],
               "separate_k": found["separate_k"]}
        if found["scores"] is not None:
            scores = found["scores"]
            difference = scores.loc["shared"] - scores.loc["separate"]
            logged = (np.log(scores.clip(lower=TIE_LOG_FLOOR))
                      .loc["shared"] - np.log(scores.clip(lower=TIE_LOG_FLOOR))
                      .loc["separate"])
            folds = len(difference)
            row["raw_t"] = float(difference.mean()
                                 / (difference.std(ddof=1) / np.sqrt(folds)))
            row["log_t"] = float(logged.mean()
                                 / (logged.std(ddof=1) / np.sqrt(folds)))
            row["folds"] = int(folds)
        detail[str(seed)] = row
    report = {"substrate": substrate, "family": family, "k_rate": k_rate,
              "k_clock": k_clock, "seeds": list(seeds), "detail": detail,
              "caught": int(sum(1 for row in detail.values()
                                if row["verdict"] == "separate K predicts better"))}
    saved[key] = report
    with open(path, "w") as handle:
        json.dump(saved, handle, default=float, indent=1, sort_keys=True)
    return report


def _coefficient_count(model, buffers):
    """
    The coefficients a model's law carries: every family option one,
    `power_by_buffer` one per buffer in the table, intercepts none.
    """
    count = 0
    for option in model.values():
        if option is None:
            continue
        count += len(buffers) if option == "power_by_buffer" else 1
    return count


def stage_b_candidates(substrate):
    """
    Stage B's candidate laws, from the `v2` `_all` saves, per element.

    `best` is the lowest CV total among the tied models; `simplest` the fewest
    coefficients (ties broken by the lower CV total; the empty model is
    allowed). A substrate without a fitted `k_act_burst` table gets the
    fallback `"lag law + burst_offset"` for both. `combinations` is every
    choice of one candidate per element, with `k_sink` also None.
    """
    parameters = curve_parameters(substrate)
    elements = {}
    for element in ELEMENTS:
        path = os.path.join(RATE_LAW_DIR, f"{substrate}_{element}_all.json")
        if not os.path.exists(path):
            if element != "k_act_burst":
                raise FileNotFoundError(path)
            elements[element] = {"best": "lag law + burst_offset",
                                 "simplest": "lag law + burst_offset",
                                 "candidates": ["lag law + burst_offset"]}
            continue
        report = search(substrate, element, workers=1)
        tied = report["tie"]
        buffers = sorted(set(parameters["elements"][element].buffer))
        scores = report["scores"]
        best = min(tied, key=lambda identifier: scores[identifier])
        counts = {identifier: _coefficient_count(_parse_model_id(identifier),
                                                 buffers)
                  for identifier in tied}
        simplest = min(tied, key=lambda identifier: (counts[identifier],
                                                     scores[identifier]))
        candidates = [best] if best == simplest else [best, simplest]
        elements[element] = {"best": best, "simplest": simplest,
                             "candidates": candidates}
    combinations = []
    sinks = list(elements["k_sink"]["candidates"]) + [None]
    for rate, lag, burst, sink in itertools.product(
            elements["v_act"]["candidates"], elements["k_act_lag"]["candidates"],
            elements["k_act_burst"]["candidates"], sinks):
        combinations.append({"v_act": rate, "k_act_lag": lag,
                             "k_act_burst": burst, "k_sink": sink})
    return {"elements": elements, "combinations": combinations}


def shared_binding_law(substrate, family, rate_table=None, clock_table=None):
    """
    Do `v_act`'s bind and `k_act_lag`'s relax share ONE binding constant?

    A pre-equilibrium that activates the catalyst predicts `v_act` saturating
    as K x / (1 + K x) and the clock rising as (1 + K x), with the same K in
    both. Each element takes its own best tied model, its family option forced
    to `bind` (v_act) or `relax` (k_act_lag), and the joint fit is scored by
    leave-one-run-out over the runs present in either table, one K shared
    against one K each. Three checks, in this order: if the shared K or either
    separate K sits within 0.05 of a log10 bound the verdict is
    `"not identified (K at bound)"`; if either table has fewer than two buffers
    whose own rows move the family's concentration the verdict is
    `"not identified (one buffer)"`; otherwise the verdict is the Task 4 tie
    rule on the two variants.

    Named `shared_binding_law` because `saturation.shared_binding` already
    holds the plan's `shared_binding` name. `rate_table` and `clock_table`
    override the real tables, for the planted tests.
    """
    name = {"BUF": "K_B", "H": "K_H"}[family]
    model_rate = _best_tied_model(substrate, "v_act", family, "bind")
    model_clock = _best_tied_model(substrate, "k_act_lag", family, "relax")
    if rate_table is None or clock_table is None:
        parameters = curve_parameters(substrate)
        if rate_table is None:
            rate_table = parameters["elements"]["v_act"]
        if clock_table is None:
            clock_table = parameters["elements"]["k_act_lag"]
    tables = [rate_table, clock_table]
    models = [model_rate, model_clock]
    shared_fit = _binding_fit(tables, models, name, True)
    separate_fit = _binding_fit(tables, models, name, False)
    shared_index = next(global_index for element, local, global_index
                        in shared_fit["layout"] if local == name)
    separate_indices = [global_index for element, local, global_index
                        in separate_fit["layout"] if local == name]
    health = []
    for label, fit, indices in (("shared", shared_fit, [shared_index]),
                               ("separate", separate_fit,
                                separate_indices)):
        for index in indices:
            value = float(fit["x"][index])
            if (abs(value - NONLINEAR_BOUNDS[0]) <= 0.05
                    or abs(value - NONLINEAR_BOUNDS[1]) <= 0.05):
                health.append(f"{label} at bound: {name}")
    report = {
        "shared_k": float(10.0 ** shared_fit["x"][shared_index]),
        "separate_k": [float(10.0 ** separate_fit["x"][index])
                       for index in separate_indices],
        "shared_model": model_rate,
        "clock_model": model_clock,
        "health": health,
    }
    if any(entry.endswith(f"at bound: {name}") for entry in health):
        report.update({"verdict": "not identified (K at bound)", "scores": None,
                       "shared_skipped": None, "separate_skipped": None})
        return report
    if (len(_carrying_buffers(rate_table, family)) < 2
            or len(_carrying_buffers(clock_table, family)) < 2):
        report.update({"verdict": "not identified (one buffer)", "scores": None,
                       "shared_skipped": None, "separate_skipped": None})
        return report
    shared_scores, shared_skipped = _binding_cross_validate(
        tables, models, name, True)
    separate_scores, separate_skipped = _binding_cross_validate(
        tables, models, name, False)
    scores = pd.DataFrame(
        [pd.Series(shared_scores), pd.Series(separate_scores)],
        index=["shared", "separate"]).sort_index(axis=1)
    tied = tie_set(scores)
    verdict = ("one K ties with separate K" if "shared" in tied
               else "separate K predicts better")
    report.update({"verdict": verdict, "scores": scores,
                   "shared_skipped": shared_skipped,
                   "separate_skipped": separate_skipped})
    return report


def _clock_stacked(tables, model, nonlinear, indices=None):
    """The stacked design for a shared clock law: term coefficients shared,
    one intercept per (family, buffer). `indices` names the intercept blocks
    when a single table is being predicted on its own."""
    indices = indices or list(range(len(tables)))
    sizes = [len(table) for table in tables]
    total = sum(sizes)
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    offset = np.zeros(total)
    columns, names = [], []
    for position, (index, table) in enumerate(zip(indices, tables)):
        for buffer in sorted(set(table.buffer)):
            column = np.zeros(total)
            column[starts[position]:starts[position] + len(table)] = np.where(
                (table.buffer == buffer).to_numpy(), 1.0, 0.0)
            columns.append(column)
            names.append(f"intercept[{index}][{buffer}]")
    for family in FAMILIES:
        option = model.get(family)
        if option is None:
            continue
        per_table = [_term_terms(table, family, option, nonlinear)
                     for table in tables]
        expected = [name for name, _, _, _ in per_table[0]]
        for position, name in enumerate(expected):
            if any([n for n, _, _, _ in terms][position] != name
                   for terms in per_table):
                raise ValueError(
                    "a shared clock law needs the same term names in "
                    "every table")
            stacked = np.concatenate([terms[position][1]
                                      for terms in per_table])
            if per_table[0][position][2]:
                columns.append(stacked)
                names.append(name)
            else:
                offset = offset + stacked
    return np.column_stack(columns), names, offset


def _clock_shared_fit(tables, model, starts=(-2.0, 0.0, 2.0)):
    """Fit one clock law to both clock families at once, terms shared."""
    y = np.concatenate([table.y.to_numpy(dtype=float) for table in tables])
    weight = np.concatenate([table.weight.to_numpy(dtype=float)
                             for table in tables])
    names_nonlinear = _nonlinear_names_for_model(model)

    def residuals(x):
        nonlinear = dict(zip(names_nonlinear, np.atleast_1d(x)))
        matrix, _, offset = _clock_stacked(tables, model, nonlinear)
        root = np.sqrt(weight)
        target = y - offset
        beta, *_ = np.linalg.lstsq(matrix * root[:, None], target * root,
                                   rcond=None)
        return (target - matrix @ beta) * root

    best_x, best_cost = None, np.inf
    for start in itertools.product(starts, repeat=len(names_nonlinear)):
        solution = least_squares(residuals, np.asarray(start, dtype=float),
                                 bounds=NONLINEAR_BOUNDS)
        if float(2.0 * solution.cost) < best_cost:
            best_x, best_cost = solution.x, float(2.0 * solution.cost)
    nonlinear = dict(zip(names_nonlinear, np.atleast_1d(best_x)))
    matrix, names, offset = _clock_stacked(tables, model, nonlinear)
    root = np.sqrt(weight)
    target = y - offset
    beta, *_ = np.linalg.lstsq(matrix * root[:, None], target * root,
                               rcond=None)
    return {"coefficients": {name: float(value)
                             for name, value in zip(names, beta)},
            "nonlinear": nonlinear, "cost": best_cost}


def _clock_predict(table, index, model, fit, shared=False):
    """Predict one table from a shared clock fit (or a plain fit)."""
    if len(table) == 0:
        return np.zeros(0), np.zeros(0)
    if not shared:
        return _predict(table, model, fit)
    matrix, names, offset = _clock_stacked([table], model, fit["nonlinear"],
                                           indices=[index])
    beta = np.array([fit["coefficients"][name] for name in names])
    return table.y.to_numpy(dtype=float), matrix @ beta + offset


def shared_clock_law(substrate, starts=(-2.0, 0.0, 2.0)):
    """
    Can the burst clock follow the lag clock's law, with its own intercept?

    The `k_act_lag` best tied model is fitted to both clock tables at once --
    every term coefficient shared, one intercept per (family, buffer) -- and
    scored against the two tables' own best tied models, each fitted alone.
    Leave-one-run-out over the runs present in either table; verdict by the
    Task 4 tie rule. The plan's burst-clock fallback: a substrate whose
    `k_act_burst` table is under 15 curves is reported, not fitted.
    """
    parameters = curve_parameters(substrate)
    lag_table = parameters["elements"]["k_act_lag"]
    burst_table = parameters["elements"]["k_act_burst"]
    if len(burst_table) < 15:
        return {"verdict": "too few curves"}
    lag_model = _best_tied_model(substrate, "k_act_lag")
    burst_model = _best_tied_model(substrate, "k_act_burst")
    tables = [lag_table, burst_table]
    runs = sorted(set(lag_table.experiment) | set(burst_table.experiment))
    shared_scores, separate_scores, skipped = {}, {}, 0
    for run in runs:
        train = [table[table.experiment != run] for table in tables]
        held = [table[table.experiment == run] for table in tables]
        if any(not set(h.buffer) <= set(t.buffer)
               for h, t in zip(held, train)):
            skipped += 1
            continue
        shared_fit = _clock_shared_fit(train, lag_model, starts)
        separate_fits = [fit_model(train[0], lag_model, starts=starts),
                         fit_model(train[1], burst_model, starts=starts)]
        shared_score, separate_score = 0.0, 0.0
        for index, (h, table_model) in enumerate(zip(held, [lag_model,
                                                            burst_model])):
            y, prediction = _clock_predict(h, index, lag_model, shared_fit,
                                           shared=True)
            shared_score += float(np.sum(h.weight.to_numpy(dtype=float)
                                         * (y - prediction) ** 2))
            y, prediction = _clock_predict(h, index, table_model,
                                           separate_fits[index])
            separate_score += float(np.sum(h.weight.to_numpy(dtype=float)
                                           * (y - prediction) ** 2))
        shared_scores[int(run)] = shared_score
        separate_scores[int(run)] = separate_score
    scores = pd.DataFrame([pd.Series(shared_scores), pd.Series(separate_scores)],
                          index=["shared", "separate"]).sort_index(axis=1)
    tied = tie_set(scores)
    verdict = ("one clock law ties with separate laws" if "shared" in tied
               else "separate clock laws predict better")
    return {"verdict": verdict, "scores": scores, "skipped": skipped,
            "lag_model": lag_model, "burst_model": burst_model}


def barriers(substrate="4OMe-BnOH"):
    """
    The activation energies of the elements, and whether they agree.

    Each element's best tied model containing `T:arrhenius` is fitted on its
    full table; `Ea` and its run-clustered standard error come off that fit, in
    kJ/mol. An element whose table moves fewer than 4 temperatures is left out
    and listed as `"not identified (fewer than 4 temperatures)"`; an element
    whose `Ea` is otherwise flagged (a collinear pair naming `Ea_R`, a rank
    deficiency) is left out of the comparison too. Verdict over the rest:
    "barriers equal within 2 standard errors", or
    "barriers differ by more than 2 standard errors: <pairs>", or
    "fewer than 2 elements with a barrier".
    """
    rows = {}
    tables = curve_parameters(substrate)["elements"]
    for element in ELEMENTS:
        if len(tables[element]) < 15:
            continue
        if tables[element].kelvin.nunique() < 4:
            rows[element] = {
                "model": None, "Ea": None, "stderr": None, "flags": [],
                "flagged": True,
                "status": "not identified (fewer than 4 temperatures)"}
            continue
        report = search(substrate, element, workers=1)
        candidates = [identifier for identifier in report["tie"]
                      if "T:arrhenius" in identifier.split("|")]
        if not candidates:
            continue
        best = min(candidates,
                   key=lambda identifier: report["scores"][identifier])
        fit = fit_model(tables[element], _parse_model_id(best))
        flags = fit["health"]["flags"]
        flagged = [flag for flag in flags
                   if "Ea_R" in flag or flag == "rank deficient"]
        rows[element] = {"model": best, "Ea": fit["coefficients"]["Ea"],
                         "stderr": fit["stderr_clustered"]["Ea"],
                         "flags": flags, "flagged": bool(flagged),
                         "status": None}
    usable = {element: row for element, row in rows.items()
              if not row["flagged"]}
    if len(usable) < 2:
        return {"verdict": "fewer than 2 elements with a barrier",
                "elements": rows, "pairs": []}
    pairs = []
    names = sorted(usable)
    for i, first in enumerate(names):
        for second in names[i + 1:]:
            difference = usable[first]["Ea"] - usable[second]["Ea"]
            error = np.sqrt(usable[first]["stderr"] ** 2
                            + usable[second]["stderr"] ** 2)
            if abs(difference) > 2.0 * error:
                pairs.append(f"{first}/{second}")
    verdict = ("barriers equal within 2 standard errors" if not pairs
               else "barriers differ by more than 2 standard errors: "
                    + ", ".join(pairs))
    return {"verdict": verdict, "elements": rows, "pairs": pairs}


# --- the global analytic fit (Tasks 6-7) ------------------------------------


_GLOBAL_LAWS = ("v_act", "k_act_lag", "k_act_burst", "k_sink")


def _global_curves(substrate):
    """
    The curves Stage B fits: section 3's rows with `v_act > 0`, each with its
    own times, gas-corrected readings and noise.

    Whether an element resolved on a curve does not matter here -- Stage B fits
    the readings, and a curve with an unresolved element still constrains the
    laws through its other elements and its own shape. Curves the sink form
    cannot hold are already out; the exclusion here is counted in
    `dropped_nonpositive`.
    """
    parameters = curve_parameters(substrate)
    rows = parameters["rows"]
    dropped = int((rows.v_act <= 0).sum())
    rows = rows[rows.v_act > 0].reset_index(drop=True)
    fits = scope.fits(scope.archive())
    times, values, noise = [], [], []
    for row in rows.itertuples():
        fit = fits[(int(row.experiment), int(row.sample))]
        times.append(np.asarray(fit.times, dtype=float))
        values.append(np.asarray(fit.corrected, dtype=float))
        noise.append(float(fit.noise))
    return {"substrate": substrate, "rows": rows, "times": times,
            "values": values, "noise": noise, "dropped_nonpositive": dropped}


def _global_subset(curves, keep):
    """The same curve set, restricted to a boolean mask over its rows."""
    index = np.where(keep)[0]
    return {"substrate": curves["substrate"],
            "rows": curves["rows"].iloc[index].reset_index(drop=True),
            "times": [curves["times"][i] for i in index],
            "values": [curves["values"][i] for i in index],
            "noise": [curves["noise"][i] for i in index],
            "dropped_nonpositive": 0}


def _global_layout(rows, laws):
    """
    Stage B's coefficient layout for one candidate on one curve set.

    `laws` maps each element to a model id (`""` is the intercepts-only
    model), `None` for `k_sink` only, or the string
    `"lag law + burst_offset"` where a substrate's `k_act_burst` table was too
    small. Returns the parsed models, the fallback flag, and `entries`: one
    `(key, kind, law, name)` per global coefficient -- the nonlinear constants
    first per law, the burst offset last. `key` is the coefficient's name
    everywhere downstream.
    """
    fallback = laws.get("k_act_burst") == "lag law + burst_offset"
    models, entries = {}, []
    for law in _GLOBAL_LAWS:
        if law == "k_act_burst" and fallback:
            continue
        identifier = laws.get(law)
        model = None if identifier is None else _parse_model_id(identifier)
        models[law] = model
        if model is None:
            continue
        nonlinear = _nonlinear_names_for_model(model)
        for name in nonlinear:
            entries.append((f"{law}:log10({name})", "nonlinear", law, name))
        _, names, _, _ = law_design(rows, model,
                                    {name: 0.0 for name in nonlinear})
        for name in names:
            entries.append((f"{law}:{name}", "linear", law, name))
    if fallback:
        entries.append(("burst_offset", "offset", "k_act_burst",
                        "burst_offset"))
    return {"models": models, "fallback": fallback, "entries": entries}


def _global_x(layout, values):
    """The global vector for a {key: value} dict, missing keys at 0.0."""
    return np.array([float(values.get(key, 0.0))
                     for key, _, _, _ in layout["entries"]])


def _global_terms(x, layout, rows):
    """(log v_act, log k_act, log k) per row at a global vector."""
    values = {key: float(value)
              for (key, _, _, _), value in zip(layout["entries"], x)}
    predictions = {}
    for law, model in layout["models"].items():
        if model is None:
            continue
        nonlinear = {name: values[f"{law}:log10({name})"]
                     for name in _nonlinear_names_for_model(model)}
        matrix, names, offset, _ = law_design(rows, model, nonlinear)
        beta = np.array([values[f"{law}:{name}"] for name in names])
        predictions[law] = matrix @ beta + offset
    log_v_act = predictions["v_act"]
    is_lag = (rows.family == "lag").to_numpy()
    if layout["fallback"]:
        log_k_act = predictions["k_act_lag"] + np.where(
            is_lag, 0.0, values["burst_offset"])
    else:
        log_k_act = np.where(is_lag, predictions["k_act_lag"],
                             predictions["k_act_burst"])
    if layout["models"].get("k_sink") is None:
        log_k = np.full(len(rows), -np.inf)
    else:
        log_k = np.clip(predictions["k_sink"], -_LOG_RATE_CLIP, _LOG_RATE_CLIP)
    return (np.clip(log_v_act, -_LOG_RATE_CLIP, _LOG_RATE_CLIP),
            np.clip(log_k_act, -_LOG_RATE_CLIP, _LOG_RATE_CLIP),
            log_k)


def _global_readings(curves, layout, x, c, v0):
    """The model's readings for every curve at a global vector, each curve
    with its own `c` and `v0`."""
    rows, times = curves["rows"], curves["times"]
    log_v_act, log_k_act, log_k = _global_terms(x, layout, rows)
    readings = []
    for index, curve_times in enumerate(times):
        v_act = float(np.exp(log_v_act[index]))
        tau = float(1.0 / np.exp(log_k_act[index]))
        k = float(np.exp(log_k[index])) if np.isfinite(log_k[index]) else 0.0
        h, g = summary_kinetics._activation_sink_columns(curve_times, tau, k)
        readings.append(c[index] + v0[index] * h[0] + v_act * g[0])
    return readings


def _global_stack(x, layout, curves):
    """(residuals, models, terms) at a global vector.

    Per curve the `(c, v0)` pair is solved out: a lstsq of
    `readings - v_act g` on `[1, h]`, so the stacked residual is orthogonal to
    the curve's own offset and initial rate.
    """
    rows, times, values, noise = (curves["rows"], curves["times"],
                                  curves["values"], curves["noise"])
    log_v_act, log_k_act, log_k = _global_terms(x, layout, rows)
    residuals, models = [], []
    for index, curve_times in enumerate(times):
        v_act = float(np.exp(log_v_act[index]))
        tau = float(1.0 / np.exp(log_k_act[index]))
        k = float(np.exp(log_k[index])) if np.isfinite(log_k[index]) else 0.0
        h, g = summary_kinetics._activation_sink_columns(curve_times, tau, k)
        target = values[index] - v_act * g[0]
        design = np.column_stack([np.ones_like(curve_times), h[0]])
        beta, *_ = np.linalg.lstsq(design, target, rcond=None)
        model = design @ beta + v_act * g[0]
        residuals.append((values[index] - model)
                         / (noise[index] * np.sqrt(len(curve_times))))
        models.append(model)
    return residuals, models, (log_v_act, log_k_act, log_k)


def _stage_a_starts(substrate, layout):
    """Stage A's coefficient estimates for a layout, as a {key: value} dict."""
    parameters = curve_parameters(substrate)
    starts = {}
    for law, model in layout["models"].items():
        if model is None:
            continue
        table = parameters["elements"][law]
        if len(table) < 15:
            continue
        fit = fit_model(table, model)
        for name, value in fit["coefficients"].items():
            starts[f"{law}:{name}"] = float(value)
    if layout["fallback"]:
        starts["burst_offset"] = 0.0
    return starts


def global_fit(substrate, laws, curves=None, starts=None, restarts=4):
    """
    The activation-sink form fitted to many curves at once, its parameters
    replaced by Stage A's rate laws.

    `laws` is one candidate: a model id per element, `None` for `k_sink`
    (meaning k = 0 everywhere) or the fallback string for `k_act_burst`.
    `curves` overrides the curve set (a `_global_curves` dict, readings
    replaced where planted). Each curve's `c` and `v0` are free and solved by
    lstsq inside every residual evaluation; the global coefficients -- every
    law's intercepts, orders, `Ea/R` and log10 K -- are optimised by
    `scipy.optimize.least_squares` on the stacked residuals. `starts` overrides
    Stage A's estimates; `restarts` Latin-hypercube draws within +/-1 are added
    (log10 K is bounded at (-4, 4), so a draw is clipped).

    Returns the coefficients under their law-prefixed keys, the naive and
    run-clustered standard errors, the cost, per-curve rms in units of noise,
    per-curve `net_data` and `net_model`, `converged`, `health`, and the
    underscore keys the later steps read.
    """
    data = _global_curves(substrate) if curves is None else curves
    rows = data["rows"]
    layout = _global_layout(rows, laws)
    keys = [key for key, _, _, _ in layout["entries"]]
    base = _stage_a_starts(substrate, layout) if starts is None else {
        key: float(value) for key, value in starts.items()}
    x0 = np.array([base.get(key, 0.0) for key in keys])
    kinds = [kind for _, kind, _, _ in layout["entries"]]
    lower = np.array([NONLINEAR_BOUNDS[0] if kind == "nonlinear" else -np.inf
                      for kind in kinds])
    upper = np.array([NONLINEAR_BOUNDS[1] if kind == "nonlinear" else np.inf
                      for kind in kinds])
    starts_list = [np.clip(x0, lower, upper)]
    if restarts > 0:
        sampler = qmc.LatinHypercube(d=len(keys), seed=0)
        for draw in sampler.random(restarts):
            starts_list.append(np.clip(x0 + (2.0 * draw - 1.0), lower, upper))

    def objective(x):
        residuals, _, _ = _global_stack(x, layout, data)
        return np.concatenate(residuals)

    best = None
    for start in starts_list:
        solution = least_squares(objective, start, bounds=(lower, upper),
                                 x_scale="jac")
        cost = float(2.0 * solution.cost)
        if best is None or cost < best["cost"]:
            best = {"x": solution.x, "cost": cost,
                    "success": bool(solution.success), "jacobian": solution.jac}
    x = best["x"]
    residuals, models, _ = _global_stack(x, layout, data)
    stacked = np.concatenate(residuals)
    jacobian = best["jacobian"]
    bread = np.linalg.pinv(jacobian.T @ jacobian)
    degrees = max(1, len(stacked) - len(x))
    variance = float(np.sum(stacked ** 2)) / degrees
    errors = np.sqrt(np.maximum(np.diag(bread) * variance, 0.0))
    safe = np.where(errors > 0, errors, np.inf)
    correlation = bread * variance / np.outer(safe, safe)
    spans, start = [], 0
    for index, run in enumerate(rows.experiment.to_numpy()):
        length = len(data["times"][index])
        spans.append((int(run), slice(start, start + length)))
        start += length
    groups = sorted({run for run, _ in spans})
    meat = np.zeros((len(x), len(x)))
    for run in groups:
        mask = np.zeros(len(stacked), dtype=bool)
        for curve_run, span in spans:
            if curve_run == run:
                mask[span] = True
        score = jacobian[mask].T @ stacked[mask]
        meat += np.outer(score, score)
    if len(groups) < 2:
        clustered = np.full(len(x), np.nan)
    else:
        covariance = len(groups) / (len(groups) - 1.0) * bread @ meat @ bread
        clustered = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    rms, net_data, net_model = [], [], []
    for index in range(len(rows)):
        difference = data["values"][index] - models[index]
        rms.append(float(np.sqrt(np.mean(
            (difference / data["noise"][index]) ** 2))))
        net_data.append(float(data["values"][index][-1]
                              - data["values"][index][0]))
        net_model.append(float(models[index][-1] - models[index][0]))
    coefficients = {key: float(value) for key, value in zip(keys, x)}
    result = {
        "substrate": substrate, "laws": dict(laws),
        "coefficients": coefficients,
        "stderr": {key: float(value) for key, value in zip(keys, errors)},
        "stderr_clustered": {key: float(value)
                             for key, value in zip(keys, clustered)},
        "cost": best["cost"], "converged": best["success"],
        "curves": rows, "rms": np.array(rms),
        "net_data": np.array(net_data), "net_model": np.array(net_model),
        "curves_used": int(len(rows)),
        "nonlinear_at_bound": [
            key for (key, kind, _, _), value in zip(layout["entries"], x)
            if kind == "nonlinear"
            and (abs(value - NONLINEAR_BOUNDS[0]) <= 0.05
                 or abs(value - NONLINEAR_BOUNDS[1]) <= 0.05)],
        "correlation": pd.DataFrame(correlation, index=keys, columns=keys),
        "_x": x, "_layout": layout, "_jacobian": jacobian,
        "_residuals": stacked,
    }
    result["health"] = global_fit_health(result)
    return result


def global_fit_health(result):
    """
    The flags a global fit is read with: `"not converged"`, `"at bound: <K>"`,
    and `"collinear: <a>/<b> r=+/-0.xx"` for a coefficient pair (intercepts
    excluded) above |r| 0.99. Returns `{"flags", "ok"}`.
    """
    flags = []
    if not result["converged"]:
        flags.append("not converged")
    for key in result["nonlinear_at_bound"]:
        flags.append(f"at bound: {key}")
    correlation = result["correlation"]
    columns = list(correlation.columns)
    for i, first in enumerate(columns):
        if ":intercept[" in first:
            continue
        for second in columns[i + 1:]:
            if ":intercept[" in second:
                continue
            r = float(correlation.loc[first, second])
            if np.isfinite(r) and abs(r) > 0.99:
                flags.append(f"collinear: {first}/{second} r={r:+.2f}")
    return {"flags": flags, "ok": not flags}


def global_cross_validate(substrate, laws, full=None, curves=None):
    """
    Leave-one-run-out for one Stage B candidate.

    Each fold fits on the other runs, started from the full fit's coefficients,
    predicts the held-out run from the laws, and refits each held-out curve's
    own `(c, v0)` by lstsq before scoring
    `sum ((readings - model)/(noise sqrt(n)))^2`. A fold whose held-out run
    carries a buffer no other run has is skipped and counted.
    """
    data = _global_curves(substrate) if curves is None else curves
    if full is None:
        full = global_fit(substrate, laws, curves=data)
    rows = data["rows"]
    experiments = rows.experiment.to_numpy()
    scores, skipped = {}, 0
    for run in sorted(set(experiments)):
        keep = experiments != run
        train = _global_subset(data, keep)
        held = np.where(~keep)[0]
        if not set(rows.buffer.iloc[held]) <= set(train["rows"].buffer):
            skipped += 1
            continue
        fit = global_fit(substrate, laws, curves=train,
                         starts=full["coefficients"], restarts=1)
        held_rows = rows.iloc[held].reset_index(drop=True)
        log_v_act, log_k_act, log_k = _global_terms(fit["_x"],
                                                    fit["_layout"], held_rows)
        score = 0.0
        for position, index in enumerate(held):
            curve_times = data["times"][index]
            values = data["values"][index]
            noise = data["noise"][index]
            v_act = float(np.exp(log_v_act[position]))
            tau = float(1.0 / np.exp(log_k_act[position]))
            k = (float(np.exp(log_k[position]))
                 if np.isfinite(log_k[position]) else 0.0)
            h, g = summary_kinetics._activation_sink_columns(curve_times, tau, k)
            target = values - v_act * g[0]
            design = np.column_stack([np.ones_like(curve_times), h[0]])
            beta, *_ = np.linalg.lstsq(design, target, rcond=None)
            model = design @ beta + v_act * g[0]
            score += float(np.sum(
                ((values - model) / (noise * np.sqrt(len(curve_times)))) ** 2))
        scores[int(run)] = score
    return {"scores": pd.Series(scores, dtype=float).sort_index(),
            "sum": float(sum(scores.values())), "skipped": skipped,
            "full": full}


def stage_shift(stage_a, stage_b):
    """
    Every coefficient present in both stages, shifted in units of the combined
    run-clustered errors, and the count beyond 2.

    `stage_a` is {element: `fit_model` result}; `stage_b` is a `global_fit`
    result. Keys line up as `f"{element}:{name}"` on both sides.
    """
    rows = []
    for element, fit in stage_a.items():
        for name, value in fit["coefficients"].items():
            key = f"{element}:{name}"
            if key not in stage_b["coefficients"]:
                continue
            se_a = fit["stderr_clustered"].get(name)
            se_b = stage_b["stderr_clustered"].get(key)
            if not np.isfinite([se_a, se_b]).all():
                continue
            shift = float((stage_b["coefficients"][key] - value)
                          / np.sqrt(se_a ** 2 + se_b ** 2))
            rows.append({"coefficient": key, "stage_a": float(value),
                         "stage_b": float(stage_b["coefficients"][key]),
                         "stderr_a": float(se_a), "stderr_b": float(se_b),
                         "shift": shift, "beyond_2": abs(shift) > 2.0})
    frame = pd.DataFrame(rows).set_index("coefficient") if rows else pd.DataFrame(
        columns=["stage_a", "stage_b", "stderr_a", "stderr_b", "shift",
                 "beyond_2"]).set_index(pd.Index([], name="coefficient"))
    return {"table": frame,
            "beyond": int(frame["beyond_2"].sum()) if len(frame) else 0,
            "count": int(len(frame))}


def amplitude_verdict(result, rows=None):
    """
    Whether the amplitude `net_data/net_model` tracks the conditions.

    Per curve the ratio of the readings' net rise to the model's; per run the
    median. A between-run axis is flagged at |t| > 3 with at least 4 runs, from
    an OLS of log median ratio on log[HOO-], log[buf], log[S], log[enz] and
    1/T; a within-run axis at |t| > 3 from `scope.orders` on log buf and
    log[S]. `hoo` is never put in the within-run regression. Non-positive
    ratios are dropped and counted. Verdict: `"tracks conditions: <flags>"` or
    `"no amplitude trend detected"`.
    """
    frame = (result["curves"] if rows is None else rows).copy()
    frame["ratio"] = result["net_data"] / result["net_model"]
    dropped = int(np.sum(~(frame["ratio"] > 0)))
    frame = frame[frame["ratio"] > 0].copy()
    grouped = frame.groupby("experiment")
    runs = grouped["ratio"].median()
    axes = (("hoo", "hoo", False), ("buf", "buf", False),
            ("s0", "s0", False), ("e0", "e0", False),
            ("temperature", "kelvin", True))
    flags, between = [], {}
    for name, column, inverse in axes:
        if len(runs) < 4:
            continue
        x = grouped[column].median().to_numpy(dtype=float)
        if inverse:
            x = 1.0 / x
        else:
            x = np.log(x)
        y = np.log(runs.to_numpy(dtype=float))
        design = np.column_stack([np.ones(len(runs)), x])
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        residual = y - design @ beta
        covariance = np.linalg.pinv(design.T @ design) * float(
            residual @ residual) / max(1, len(runs) - 2)
        t = float(beta[1] / np.sqrt(covariance[1, 1]))
        between[name] = t
        if abs(t) > 3.0:
            flags.append(f"{name}: t={t:+.2f}")
    within = scope.orders("ratio", frame=frame, terms=("buf", "s0"),
                          within=True, live_only=False)
    for name in ("buf", "s0"):
        order, error = within.get(f"order_{name}"), within.get(f"stderr_{name}")
        if order is not None and np.isfinite([order, error]).all() and error > 0:
            t = float(order / error)
            if abs(t) > 3.0:
                flags.append(f"{name} within runs: t={t:+.2f}")
    verdict = ("tracks conditions: " + "; ".join(flags) if flags
               else "no amplitude trend detected")
    return {"verdict": verdict, "dropped": dropped, "between": between,
            "within": within, "flags": flags}


def _laws_id(laws):
    """A file-name-safe id for one candidate's four laws."""
    return ";".join(f"{law}={laws.get(law) if laws.get(law) is not None
                             else 'None'}"
                    for law in _GLOBAL_LAWS)


def _truth_coefficients(substrate, layout):
    """The planted truth's coefficients for a layout: each element's own
    `fit_model` on its full table, and the burst offset as the median of the
    burst table's residual to the lag law."""
    parameters = curve_parameters(substrate)
    truth = {}
    for law, model in layout["models"].items():
        if model is None:
            continue
        fit = fit_model(parameters["elements"][law], model)
        for name, value in fit["coefficients"].items():
            truth[f"{law}:{name}"] = float(value)
    if layout["fallback"]:
        lag_model = layout["models"]["k_act_lag"]
        nonlinear = {name: truth[f"k_act_lag:log10({name})"]
                     for name in _nonlinear_names_for_model(lag_model)}
        table = parameters["elements"]["k_act_burst"]
        matrix, names, offset, _ = law_design(table, lag_model, nonlinear)
        beta = np.array([truth.get(f"k_act_lag:{name}", 0.0)
                         for name in names])
        residual = table.y.to_numpy(dtype=float) - (matrix @ beta + offset)
        truth["burst_offset"] = float(np.median(residual))
    return truth


def _global_save_task(payload):
    """One Stage B candidate's fit and cross-validation, for the pool."""
    substrate, laws, curves, path, cut = payload
    fit = global_fit(substrate, laws, curves=curves)
    cross = global_cross_validate(substrate, laws, full=fit, curves=curves)
    report = {
        "substrate": substrate, "identifier": _laws_id(laws),
        "laws": dict(laws), "cut": cut,
        "coefficients": fit["coefficients"], "stderr": fit["stderr"],
        "stderr_clustered": fit["stderr_clustered"],
        "cost": float(fit["cost"]), "converged": bool(fit["converged"]),
        "curves_used": int(fit["curves_used"]),
        "runs_used": int(fit["curves"].experiment.nunique()),
        "rms_median": float(np.median(fit["rms"])),
        "rms_p90": float(np.percentile(fit["rms"], 90)),
        "net_data": [float(value) for value in fit["net_data"]],
        "net_model": [float(value) for value in fit["net_model"]],
        "health": fit["health"],
        "nonlinear_at_bound": fit["nonlinear_at_bound"],
        "scores": {str(int(run)): float(value)
                   for run, value in cross["scores"].items()},
        "sum": float(cross["sum"]), "skipped": int(cross["skipped"]),
    }
    return report, path


def _global_cut(curves, cut):
    """One cut's curve set, from the rows' own cut."""
    rows = curves["rows"]
    kept = set(_apply_cut(rows, cut).index)
    return _global_subset(curves,
                          np.array([index in kept for index in rows.index]))


def global_search(substrate, cut=None, workers=8):
    """
    Every Stage B candidate for one substrate, fitted and cross-validated.

    Saves `data/fits/rate_laws/global_<substrate>_<id>[_<cut>].json`, one per
    candidate; a save that already exists is read back, the rest are fitted in
    a process pool and written as they finish.
    """
    os.makedirs(RATE_LAW_DIR, exist_ok=True)
    combinations = stage_b_candidates(substrate)["combinations"]
    data = _global_curves(substrate)
    if cut is not None:
        data = _global_cut(data, cut)
    reports, tasks = {}, []
    for laws in combinations:
        identifier = _laws_id(laws)
        suffix = f"_{cut}" if cut else ""
        path = os.path.join(
            RATE_LAW_DIR, f"global_{substrate}_{identifier}{suffix}.json")
        if os.path.exists(path):
            with open(path) as handle:
                reports[identifier] = json.load(handle)
        else:
            tasks.append((substrate, laws, data, path, cut))
    if tasks:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for report, path in pool.map(_global_save_task, tasks):
                with open(path, "w") as handle:
                    json.dump(report, handle, default=float, indent=1,
                              sort_keys=True)
                reports[report["identifier"]] = report
    return reports


def _planted_save_task(payload):
    """One candidate's planted fit and cross-validation, for the pool."""
    substrate, laws, planted, starts, path = payload
    fit = global_fit(substrate, laws, curves=planted, starts=starts)
    cross = global_cross_validate(substrate, laws, full=fit, curves=planted)
    report = {
        "identifier": _laws_id(laws),
        "coefficients": fit["coefficients"],
        "stderr_clustered": fit["stderr_clustered"],
        "health": fit["health"],
        "scores": {str(int(run)): float(value)
                   for run, value in cross["scores"].items()},
        "sum": float(cross["sum"]), "skipped": int(cross["skipped"]),
    }
    return report, path


def planted_global(substrate, seed=0, workers=8):
    """
    Stage B's realistic planted recovery, recorded and not asserted.

    The truth is each element's `best` candidate with the coefficients
    `fit_model` gives on its full table (the burst offset as the median of the
    burst table's residual to the lag law). Readings are that model at each
    curve's own times, with its own `c` and `v0`, plus Gaussian noise at that
    curve's own `noise`, drawn curve by curve in the order `global_fit`
    iterates. Every candidate from `stage_b_candidates` is fitted and
    cross-validated on the planted readings, each started from the truth + 0.5
    (never at it). Each candidate is saved as
    `data/fits/rate_laws/v2/planted_global_<substrate>_<id>.json` as it
    finishes, and the assembled record to
    `data/fits/rate_laws/v2/planted_global_<substrate>.json`.
    """
    final = os.path.join(RATE_LAW_DIR, f"planted_global_{substrate}.json")
    if os.path.exists(final):
        with open(final) as handle:
            return json.load(handle)
    chosen = stage_b_candidates(substrate)
    combinations = chosen["combinations"]
    truth_laws = {element: chosen["elements"][element]["best"]
                  for element in ELEMENTS}
    data = _global_curves(substrate)
    layout = _global_layout(data["rows"], truth_laws)
    truth = _truth_coefficients(substrate, layout)
    x_truth = _global_x(layout, truth)
    fits = scope.fits(scope.archive())
    c = np.array([fits[(int(row.experiment), int(row.sample))]
                  .activation_sink.c for row in data["rows"].itertuples()])
    v0 = np.array([fits[(int(row.experiment), int(row.sample))]
                   .activation_sink.v0 for row in data["rows"].itertuples()])
    generator = np.random.default_rng(seed)
    planted_values = [model + generator.normal(
        0.0, data["noise"][index], len(data["times"][index]))
        for index, model in enumerate(
            _global_readings(data, layout, x_truth, c, v0))]
    planted = {"substrate": substrate, "rows": data["rows"],
               "times": data["times"], "values": planted_values,
               "noise": data["noise"],
               "dropped_nonpositive": data["dropped_nonpositive"]}
    starts = {key: (truth[key] + 0.5 if key in truth else 0.0)
              for key, _, _, _ in layout["entries"]}
    reports, tasks = {}, []
    for laws in combinations:
        identifier = _laws_id(laws)
        path = os.path.join(
            RATE_LAW_DIR,
            f"planted_global_{substrate}_{identifier}.json")
        if os.path.exists(path):
            with open(path) as handle:
                reports[identifier] = json.load(handle)
        else:
            tasks.append((substrate, laws, planted, starts, path))
    if tasks:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for report, path in pool.map(_planted_save_task, tasks):
                with open(path, "w") as handle:
                    json.dump(report, handle, default=float, indent=1,
                              sort_keys=True)
                reports[report["identifier"]] = report
    scores = {identifier: report["sum"]
              for identifier, report in reports.items()}
    scores_frame = pd.DataFrame(
        {identifier: pd.Series({int(run): value for run, value
                                in report["scores"].items()})
         for identifier, report in reports.items()}).T
    statistics = tie_statistics(scores_frame)
    tied = [identifier for identifier in scores_frame.index
            if bool(statistics.loc[identifier, "tied"])]
    best = min(scores, key=scores.get)
    truth_id = _laws_id(truth_laws)
    true_fit = reports[truth_id]
    truth_record, estimate, clustered, within = {}, {}, {}, {}
    for key, kind, law, name in layout["entries"]:
        if key not in truth:
            continue
        truth_record[key] = truth[key]
        estimate[key] = true_fit["coefficients"][key]
        error = true_fit["stderr_clustered"].get(key)
        clustered[key] = error
        within[key] = bool(error is not None and np.isfinite(error)
                           and abs(estimate[key] - truth[key]) <= 2.0 * error)
    report = {
        "substrate": substrate, "seed": seed, "truth_laws": truth_laws,
        "truth": truth_record, "estimate": estimate,
        "stderr_clustered": clustered, "within_2se": within,
        "health": true_fit["health"], "tie": tied,
        "tie_statistics": [
            {"model": identifier,
             "raw_mean": float(statistics.loc[identifier, "raw_mean"]),
             "raw_bar": float(statistics.loc[identifier, "raw_bar"]),
             "raw_kept": bool(statistics.loc[identifier, "raw_kept"]),
             "log_mean": float(statistics.loc[identifier, "log_mean"]),
             "log_bar": float(statistics.loc[identifier, "log_bar"]),
             "log_kept": bool(statistics.loc[identifier, "log_kept"]),
             "worst_fold_log_ratio": float(
                 statistics.loc[identifier, "worst_fold_log_ratio"]),
             "tied": bool(statistics.loc[identifier, "tied"])}
            for identifier in scores_frame.index],
        "best": bool(best == truth_id), "tied": bool(truth_id in tied),
        "scores": {identifier: float(value)
                   for identifier, value in scores.items()},
        "skipped": {identifier: int(report["skipped"])
                    for identifier, report in reports.items()},
    }
    with open(final, "w") as handle:
        json.dump(report, handle, default=float, indent=1, sort_keys=True)
    return report
