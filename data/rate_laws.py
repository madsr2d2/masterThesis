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
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

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
# Where a nonlinear constant is sampled when the question is whether the
# DESIGN can identify it: if no log10 K on this grid varies the regressor, no
# K can. `_term_terms` takes log10 K, as `law_design` does.
_LOG10_K_GRID = np.linspace(-4.0, 4.0, 33)


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


def fit_model(table, model, starts=(-2.0, 0.0, 2.0)):
    """
    Weighted least squares for one model on one element table.

    The linear coefficients -- intercepts, orders, `Ea_R` -- are solved at
    fixed nonlinear constants; each nonlinear constant is optimised in log10
    by `scipy.optimize.least_squares` within (-4, 4). Every combination of
    `starts` is tried and the lowest weighted cost kept.

    Returns `model`, `coefficients` (name -> value, with nonlinear constants
    under both `log10(K)` and `K`, and `Ea` in kJ/mol beside `Ea_R` in K),
    `stderr`, `cost` (weighted sum of squared residuals), `curves`, `runs`,
    `nonlinear` (log10 values), `nonlinear_at_bound`, `correlation` (a
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
    fewer than 15. Returns (flags, ok).
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
    terms are dropped after a test on the run-demeaned regressor. The shared
    coefficients are then read twice: once with the buffer intercepts and once
    with the run intercepts. `disagree` is
    |a - b| > 2 sqrt(se_a^2 + se_b^2), the plan's bar.

    Returns `table` (one row per shared coefficient: both estimates, both
    standard errors, `disagree`), `dropped` (terms the run design could not
    carry) and `s_buf_within_run_r`, the within-run correlation of the S and
    BUF regressors where both families are in the model.
    """
    base = fit_model(table, model)
    _, _, _, terms = law_design(table, model, base["nonlinear"])
    runs = table.experiment.to_numpy()
    kept, dropped = [], []
    for name, values, _, family in terms:
        ok, _ = _variation_ok(_run_demeaned(values, runs), family)
        (kept if ok else dropped).append(name)

    matrix, names, offset = _run_design(table, model, base["nonlinear"], kept)
    y = table.y.to_numpy(dtype=float)
    weight = table.weight.to_numpy(dtype=float)
    root = np.sqrt(weight)
    target = y - offset
    beta, *_ = np.linalg.lstsq(matrix * root[:, None], target * root,
                               rcond=None)
    residual = target - matrix @ beta
    cost = float(np.sum(weight * residual ** 2))
    degrees = max(1, len(table) - matrix.shape[1])
    covariance = np.linalg.pinv(
        matrix.T @ (matrix * weight[:, None])) * (cost / degrees)
    errors = dict(zip(names, np.sqrt(np.maximum(np.diag(covariance), 0.0))))
    values = dict(zip(names, beta))

    # Nonlinear terms carry no coefficient of their own -- their constant is
    # the same in both fits -- so only the fitted coefficients are compared.
    shared = [name for name, _, linear, _ in terms
              if linear and name in kept]
    rows = []
    for name in shared:
        a, b = base["coefficients"].get(name), values.get(name)
        se_a, se_b = base["stderr"].get(name), errors.get(name)
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
    return {"table": pd.DataFrame(rows).set_index("coefficient"),
            "dropped": tuple(dropped), "s_buf_within_run_r": s_buf}
