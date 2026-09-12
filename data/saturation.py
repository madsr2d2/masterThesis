"""
Saturating fits, per element of the fitted curve, on both concentration axes.

    python data/saturation.py

AN ORDER IS A TANGENT WHERE THE AXIS SATURATES. `scope.orders` returns
d log(response)/d log(axis), which is a constant only if the response is a
power of the axis. This archive's are not: on the peroxide axis a fractional
order of 0.596 and a binding constant K = 0.040 /mM fit the two-axis ladder
EQUALLY well (`induction.peroxide_saturation`, SSE 10.1552 against 10.1554),
and the catalyst is then 77% bound at the archive's working 82.5 mM; on the
substrate axis the 4OMe ladders resolve a per-run Km of 2-7 mM
(`ph_role.ladder_mm_table`). An order measured across rungs that straddle K
is a range-average, and moving the rungs moves it with no chemistry changing.
So where an axis saturates, the summary is the SATURATING PARAMETER and the
order is a derived local slope.

PER ELEMENT means per fitted parameter of the curve, not just "the rate".
`summary_kinetics.fit_activation_sink` gives each curve v_act (the activated
rate, no product made), tau (the catalyst's clock), k (the product-driven
decline) and v0 (the rate at mixing), and each is a separate response with its
own concentration dependence.

THE PRE-EQUILIBRIUM PREDICTS TWO DIFFERENT FORMS, ONE CONSTANT. For
E + H2O2 <-> E.H2O2 in pre-equilibrium,

    v_act  ~  K h / (1 + K h)      the activated rate saturates
    1/tau  =  k_off (1 + K h)      the relaxation clock RISES, never saturates

and their logarithmic slopes, 1/(1 + K h) and K h/(1 + K h), sum to 1 at every
h and every K. That sum IS `induction.joint_clocks`' +1 rule, which is why the
rule survives saturation while the individual orders do not.
`shared_binding` fits both forms with ONE K against each having its own, and
the F between them is the test the +1 rule cannot make: whether the binding
that sets the activated rate is the binding that sets the clock.

The sink constant k has no such prediction here -- what destroys the product
is not identified (`product_fate`) -- so it is reported as a free power.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import induction
import scope
from induction import PROFILE_F, SATURATION_GRID, SATURATION_SPAN

# The archive's working peroxide, the concentration a bound fraction is
# quoted at so two blocks' K can be compared as one number.
WORKING_PEROXIDE = 82.5

# Each element: the column it lives in, the gate that says it is pinned, the
# form the pre-equilibrium predicts for it on the peroxide axis ("bound" is
# K h/(1 + K h), "relaxation" is (1 + K h), None asks a free power only), and
# whether it is a RATE -- which decides whether a Michaelis-Menten fit is
# even the right question on the substrate axis. v = Vmax [S]/(Km + [S])
# describes a rate; a rate CONSTANT (k_act = 1/tau, the sink's k) is not that
# shape and gets an order there instead. Fitted anyway, the clocks returned
# Km at the grid floor on 39 of 42 runs with r2 ~ 0, which is what asking the
# wrong question looks like.
SATURATION_ELEMENTS = (
    ("v_act_corrected", "v_act_resolved_corrected", "bound", True),
    ("k_act_corrected", "tau_act_resolved_corrected", "relaxation", False),
    ("k_sink_corrected", None, None, False),
    ("v0_act_corrected", "v0_act_resolved_corrected", "bound", True),
    ("vmax_corrected", None, "bound", True),
    ("v_peak_corrected", None, "bound", True),
)

# How few rungs a run may ladder and still earn a Michaelis-Menten fit. Three
# is `scope.mm_fit`'s own floor: two for (Vmax, Km) and one for a residual.
MICHAELIS_MINIMUM_RUNGS = 3


def _scheme(shape, constant, peroxide):
    """The log-offset a scheme puts on log(response), per reading."""
    if shape == "bound":
        return np.log(constant * peroxide / (1.0 + constant * peroxide))
    if shape == "relaxation":
        return np.log(1.0 + constant * peroxide)
    raise ValueError(f"unknown scheme {shape!r}")


def _levelled(y, offset, groups):
    """
    Residual sum of squares once every group is free to sit where it likes.

    The groups are one per (element, experiment): a run holds pH, buffer,
    enzyme, cell and day fixed, so its own level absorbs all of them and the
    scheme is read from the peroxide contrast inside the run. Two elements
    fitted together get separate levels for the same reason -- they are
    different quantities in different units, and only their SHAPE against
    peroxide is being asked to agree.
    """
    design = np.column_stack([(groups == g).astype(float)
                              for g in np.unique(groups)])
    beta, *_ = np.linalg.lstsq(design, y - offset, rcond=None)
    residual = y - offset - design @ beta
    return float(residual @ residual), design.shape[1]


def _element_rows(table, column, gate, block_gate=True):
    """One element's usable ladder rows: the peroxide arm, gated and positive."""
    ladder = induction.peroxide_ladder(table, column)
    if gate is not None and block_gate and len(ladder):
        ladder = ladder[ladder[gate].to_numpy(dtype=bool)]
    return ladder


def binding_by_element(block=scope.TWO_AXIS_BLOCK, elements=SATURATION_ELEMENTS,
                       gated=True, grid=SATURATION_GRID, span=SATURATION_SPAN,
                       cutoff=PROFILE_F):
    """
    Each element's own K on the peroxide axis, beside its free power.

    One row per element: the free power `order` (what `scope.orders` would
    give, and a tangent), the scheme's profiled `K` with a 95% interval, the
    bound fraction it implies at `WORKING_PEROXIDE`, and both SSEs so that
    "the scheme fits as well as a free power" can be read rather than
    asserted. `gated=True` keeps only curves where that element's own
    interval resolves it; `gated=False` is the same fit on every curve the
    element is positive on, and the gap between the two is the resolution's
    price.

    The scheme's exponent is FIXED at 1 for the reason
    `induction.peroxide_saturation` fixes it: a free exponent lets the
    saturating form buy a fit that is a curve through the points and not the
    hypothesis.
    """
    table = scope.frame(block)
    constants = np.logspace(*span, grid)
    powers = np.linspace(-1.0, 2.0, grid)
    rows = []
    for column, gate, shape, _ in elements:
        ladder = _element_rows(table, column, gate, gated)
        if len(ladder) < 10:
            rows.append({"element": column, "curves": int(len(ladder))})
            continue
        h = ladder.h2o2.to_numpy(dtype=float)
        y = np.log(ladder[column].to_numpy(dtype=float))
        groups = ladder.experiment.to_numpy()
        _, levels = _levelled(y, np.zeros(len(y)), groups)
        degrees = max(1, len(y) - levels - 1)
        power_sse = np.array([_levelled(y, a * np.log(h), groups)[0]
                              for a in powers])
        best = int(np.argmin(power_sse))
        inside = powers[power_sse
                        <= power_sse[best] * (1.0 + cutoff / degrees)]
        row = {"element": column, "scheme": shape or "-",
               "curves": int(len(ladder)),
               "runs": int(ladder.experiment.nunique()),
               "order": float(powers[best]),
               "order_low": float(inside.min()),
               "order_high": float(inside.max()),
               "power_sse": float(power_sse[best])}
        if shape is not None:
            scheme_sse = np.array([
                _levelled(y, _scheme(shape, k, h), groups)[0]
                for k in constants])
            top = int(np.argmin(scheme_sse))
            allowed = constants[scheme_sse
                                <= scheme_sse[top] * (1.0 + cutoff / degrees)]
            bound = constants[top] * WORKING_PEROXIDE
            row.update(K=float(constants[top]),
                       K_low=float(allowed.min()), K_high=float(allowed.max()),
                       scheme_sse=float(scheme_sse[top]),
                       bound_fraction=float(bound / (1 + bound)))
        rows.append(row)
    return pd.DataFrame(rows).set_index("element")


def shared_binding(block=scope.TWO_AXIS_BLOCK,
                   elements=("v_act_corrected", "k_act_corrected"),
                   gated=True, grid=SATURATION_GRID, span=SATURATION_SPAN,
                   cutoff=PROFILE_F):
    """
    Do the activated rate and the activation clock share ONE binding constant?

    The pre-equilibrium says they must: v_act ~ K h/(1 + K h) saturating and
    1/tau = k_off(1 + K h) rising, the same K in both. Each element is fitted
    alone (its own K, its own per-run levels) and then together with one K
    shared, and the F is on the single degree of freedom between them. A large
    F rejects one binding for both -- which would leave the +1 rule holding
    arithmetically while the scheme behind it does not.

    Returns the two separate K, the shared K with its interval, both SSEs and
    the F. `curves` is per element; the elements do NOT have to be measurable
    on the same curves, since each carries its own per-run levels.
    """
    table = scope.frame(block)
    shapes = {column: shape for column, _, shape, _ in SATURATION_ELEMENTS}
    gates = {column: gate for column, gate, _, _ in SATURATION_ELEMENTS}
    constants = np.logspace(*span, grid)
    pieces, separate, out = [], 0.0, {"elements": list(elements)}
    for column in elements:
        ladder = _element_rows(table, column, gates[column], gated)
        if len(ladder) < 10:
            return out | {"curves": {column: int(len(ladder))}}
        h = ladder.h2o2.to_numpy(dtype=float)
        y = np.log(ladder[column].to_numpy(dtype=float))
        groups = np.array([f"{column}/{e}" for e in ladder.experiment])
        sse = np.array([_levelled(y, _scheme(shapes[column], k, h), groups)[0]
                        for k in constants])
        best = int(np.argmin(sse))
        separate += float(sse[best])
        out[f"K_{column}"] = float(constants[best])
        out[f"curves_{column}"] = int(len(ladder))
        pieces.append((y, h, groups, shapes[column]))
    total = np.zeros(len(constants))
    for index, constant in enumerate(constants):
        total[index] = sum(_levelled(y, _scheme(shape, constant, h), groups)[0]
                           for y, h, groups, shape in pieces)
    best = int(np.argmin(total))
    points = sum(len(y) for y, _, _, _ in pieces)
    levels = sum(len(np.unique(groups)) for _, _, groups, _ in pieces)
    degrees = max(1, points - levels - len(pieces))
    allowed = constants[total <= total[best] * (1.0 + cutoff / degrees)]
    return out | {
        "K_shared": float(constants[best]),
        "K_shared_low": float(allowed.min()),
        "K_shared_high": float(allowed.max()),
        "shared_sse": float(total[best]), "separate_sse": float(separate),
        "f": (float((total[best] - separate) / (separate / degrees))
              if separate > 0 else np.nan),
        "points": int(points), "degrees": int(degrees),
        "bound_fraction": float(constants[best] * WORKING_PEROXIDE
                                / (1 + constants[best] * WORKING_PEROXIDE))}


def michaelis_by_element(block=None, elements=SATURATION_ELEMENTS,
                         gated=True, minimum=MICHAELIS_MINIMUM_RUNGS):
    """
    Per-RUN (Vmax, Km) for each element, on every run that ladders [S].

    The fit is within one experiment -- exp 48's own four cuvettes give exp
    48's own (Vmax, Km) -- and a between-run axis like pH is then read off
    those per-run pairs. Only the run's substrate arm goes in (its cuvettes
    at the run's top peroxide), for the reason `ph_role._substrate_arm`
    gives: on an L, the peroxide arm's scatter would be charged to the
    substrate curve.

    `buffer_r` is the within-run correlation of log[buf] with log[S]. In both
    4OMe ladders substrate volume displaced buffer volume, so the two move
    together at about -0.96 and a Km measured there is a Km OF THE PAIR;
    exps 135-151 hold [buf] fixed and carry no such term. Read it beside
    `km_resolved`, which says whether that run's own rungs locate a
    half-saturation point at all.
    """
    data = scope.frame(scope.archive() if block is None else block)
    live = data[data.live]
    rows = []
    for experiment, run in live.groupby("experiment"):
        arm = run[np.isclose(run.h2o2.to_numpy(dtype=float),
                             float(run.h2o2.max()))]
        if arm.s0.nunique() < minimum:
            continue
        buffer_r = np.nan
        if arm.buf.nunique() > 1 and arm.s0.nunique() > 1:
            buffer_r = float(np.corrcoef(np.log(arm.s0.to_numpy(dtype=float)),
                                         np.log(arm.buf.to_numpy(dtype=float)))[0, 1])
        for column, gate, _, is_rate in elements:
            if not is_rate:
                continue          # a Michaelis-Menten fit is for a rate
            usable = arm
            if gated and gate is not None:
                usable = arm[arm[gate].to_numpy(dtype=bool)]
            if usable.s0.nunique() < minimum:
                continue
            fit = scope.mm_fit(usable.s0.to_numpy(dtype=float),
                               usable[column].to_numpy(dtype=float))
            rows.append({
                "experiment": int(experiment), "element": column,
                "substrate": run.substrate.iloc[0],
                "buffer": run.buffer.iloc[0],
                "temperature": float(run.temperature.median()),
                "e0": float(run.e0.median()),
                "catalysed": bool(run.e0.median() > 0),
                "pH": float(run.pH.median()), "rungs": fit["n"],
                "s_low": float(usable.s0.min()), "s_high": float(usable.s0.max()),
                "vmax": fit["vmax"], "km": fit["km"],
                "km_low": fit["km_interval"][0], "km_high": fit["km_interval"][1],
                "km_resolved": bool(fit["km_resolved"]), "r2": fit["r2"],
                "buffer_r": buffer_r})
    return pd.DataFrame(rows)


def michaelis_temperature(experiments=None, elements=SATURATION_ELEMENTS,
                          resolved_only=False):
    """
    WHERE THE BARRIER SITS: in Vmax, or in Km?

    The temperature series ladders four substrate rungs at one peroxide in
    every one of its six runs, so each temperature earns its own (Vmax, Km)
    and the two can be given separate Arrhenius fits. A rate read at a fixed,
    sub-saturating [S] mixes them --

        v = Vmax [S] / (Km + [S])
        d ln v / d(1/T) = d ln Vmax / d(1/T)
                          - [Km/(Km + [S])] . d ln Km / d(1/T)

    -- so an activation energy measured on the rate alone is the Vmax barrier
    only if Km does not move with temperature. `ph/` asked the same question
    of the pH axis and found BOTH moving (`ph_role.turnover_decomposition`).

    Returns a dict per element: `vmax_kJ` (Arrhenius on Vmax/[enz], the
    turnover), `km_kJ` (the same fit applied to Km -- the enthalpy of the
    binding equilibrium if Km is one, which for Km = (k_off + k_cat)/k_on it
    need not be, so read it as apparent), and the reconstruction: for each
    substrate rung, the activation energy the fitted (Vmax, Km) pair PREDICTS
    at that rung against the one the rung's own rates give
    (`arrhenius.rung_fits`). `resolved_only` drops temperatures whose Km its
    own four cuvettes do not locate -- 35 C on this series.
    """
    import arrhenius
    experiments = arrhenius.TEMPERATURE_SERIES if experiments is None \
        else experiments
    table = michaelis_by_element(experiments, elements=elements)
    out = []
    for column, _, _, is_rate in elements:
        if not is_rate:
            continue
        rows = table[table.element == column].sort_values("temperature")
        if resolved_only:
            rows = rows[rows.km_resolved]
        if len(rows) < 3:
            continue
        kelvin = rows.temperature.to_numpy(dtype=float) + 273.15
        turnover = rows.vmax.to_numpy(dtype=float) / rows.e0.to_numpy(dtype=float)
        vmax_fit = arrhenius.arrhenius_fit(kelvin, turnover) or {}
        km_fit = arrhenius.arrhenius_fit(kelvin, rows.km.to_numpy(dtype=float)) or {}
        entry = {"element": column, "temperatures": int(len(rows)),
                 "km_resolved": int(rows.km_resolved.sum()),
                 "vmax_kJ": vmax_fit.get("activation_kJ", np.nan),
                 "vmax_stderr_kJ": vmax_fit.get("stderr_kJ", np.nan),
                 "vmax_rms": vmax_fit.get("rms", np.nan),
                 "km_kJ": km_fit.get("activation_kJ", np.nan),
                 "km_stderr_kJ": km_fit.get("stderr_kJ", np.nan),
                 "km_median": float(np.nanmedian(rows.km)),
                 "km_low": float(np.nanmin(rows.km)),
                 "km_high": float(np.nanmax(rows.km))}
        # The reconstruction, rung by rung: what the (Vmax, Km) pair predicts
        # a rate at that [S] should do with temperature, against what it does.
        measured = arrhenius.rung_fits(column, experiments=experiments)
        predicted = []
        for s0 in measured.index:
            rate = (rows.vmax.to_numpy(dtype=float) * float(s0)
                    / (rows.km.to_numpy(dtype=float) + float(s0))
                    / rows.e0.to_numpy(dtype=float))
            fit = arrhenius.arrhenius_fit(kelvin, rate) or {}
            predicted.append({
                "s0": float(s0),
                "predicted_kJ": fit.get("activation_kJ", np.nan),
                "measured_kJ": float(measured.loc[s0, "activation_kJ"]),
                "measured_stderr_kJ": float(measured.loc[s0, "stderr_kJ"])})
        entry["rungs"] = pd.DataFrame(predicted)
        out.append(entry)
    return out


def michaelis_summary(table=None):
    """`michaelis_by_element` counted per element and channel."""
    table = michaelis_by_element() if table is None else table
    return table.groupby(["substrate", "buffer", "catalysed",
                          "element"]).agg(
        runs=("km", "size"),
        resolved=("km_resolved", "sum"),
        km_median=("km", lambda k: float(np.nanmedian(k[k.notna()]))),
        r2_median=("r2", "median"),
        buffer_r=("buffer_r", "median"))


def main():
    pd.set_option("display.width", 250)
    print("\nthe peroxide axis, per element (two-axis block, gated)")
    print(binding_by_element().round(4).to_string())
    print("\n...and on every curve the element is positive on")
    print(binding_by_element(gated=False).round(4).to_string())
    print("\none binding constant for the rate and the clock?")
    shared = shared_binding()
    print({k: (round(v, 4) if isinstance(v, float) else v)
           for k, v in shared.items()})
    print("\nthe substrate axis: per-run (Vmax, Km), by element and channel")
    table = michaelis_by_element()
    print(michaelis_summary(table).round(3).to_string())
    resolved = table[table.km_resolved]
    print(f"\n{len(resolved)} of {len(table)} per-run fits resolve Km; "
          f"{resolved.experiment.nunique()} experiments")
    print("\nthe temperature series: is the barrier in Vmax or in Km?")
    for resolved_only in (False, True):
        print(f"  {'Km-resolved temperatures only' if resolved_only else 'all six temperatures'}")
        for entry in michaelis_temperature(resolved_only=resolved_only):
            print(f"    {entry['element']:18s} T {entry['temperatures']} "
                  f"(Km resolved {entry['km_resolved']})  "
                  f"Vmax {entry['vmax_kJ']:6.1f} +/- {entry['vmax_stderr_kJ']:4.1f} kJ/mol   "
                  f"Km {entry['km_kJ']:+6.1f} +/- {entry['km_stderr_kJ']:4.1f} kJ/mol   "
                  f"Km {entry['km_low']:.1f}-{entry['km_high']:.1f} mM")
    print("\n  rung by rung, what (Vmax, Km) predicts against what the rung does")
    for entry in michaelis_temperature():
        if entry["element"] != "vmax_corrected":
            continue
        print(entry["rungs"].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
