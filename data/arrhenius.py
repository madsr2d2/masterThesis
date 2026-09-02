"""
Arrhenius and Eyring fits over a temperature series, and the enzyme test.

The archive holds exactly one temperature series: exps 14-19, 4OMe-BnOH in
65 mM phosphate at pH 7.00, 82.5 mM H2O2, run at 15, 20, 25, 30, 35 and 40 C
with the SAME four-rung substrate ladder (1.850, 3.700, 5.549, 7.399 mM) in
every one. Four rungs shared across six temperatures is four independent
Arrhenius fits, not one, and that is what makes the enzyme test below possible.

WHICH RATE. `vmax`, not `v0`. Twenty-three of these 24 curves clear the 3-sigma
acceleration gate and eight exceed z = +15, so an initial rate here is the
INDUCTION rate -- `curve_metrics.peak_rate` draws exactly this distinction.
`v0_quad` is carried as a sensitivity and returns a negative rate on the 15 C
bottom rung, which drops that rung.

WHAT IT IS FOR, BEYOND THE ACTIVATION ENERGY. Exps 16 and 19 record a 11.7%
lower enzyme concentration than the other four (0.240683 against 0.272695 mM),
and exp 16 sits BETWEEN two runs at the higher value.
`verify_enzyme_stock.py` shows the two values come from two real weighings and
that exp 16 is the only experiment in the campaign whose stock interrupts a
run, so the workbooks cannot settle whether it is a restock or a copied cell.
`enzyme_hypotheses` settles it from the kinetics instead: divide the rate by
each candidate enzyme concentration and ask which makes the six points most
collinear in 1/T. See temperature_series/ANALYSIS.md.
"""
import numpy as np
import pandas as pd

import scope
from scope import TEMPERATURE_SERIES

# J/mol/K, CODATA. Only ever used to turn a fitted slope into kJ/mol.
GAS_CONSTANT = 8.314462618
# The candidate enzyme concentrations, as {experiment: mM} overrides on the
# recorded value. Named rather than inlined because the whole question is which
# of these the data prefers.
ENZYME_HYPOTHESES = {
    "recorded": {},
    "exp16 restocked": {16: 0.273},
    "exps16,19 restocked": {16: 0.273, 19: 0.273},
}


def series_frame(experiments=TEMPERATURE_SERIES):
    """The temperature series' live curves, with absolute temperature."""
    frame = scope.frame(tuple(experiments))
    frame = frame[frame.live].copy()
    frame["kelvin"] = frame.temperature + 273.15
    return frame


def _line(x, y):
    """
    Weighted-free OLS of y on x with standard errors. Returns a dict.

    Written out rather than pulled from `curve_metrics.line_fit`, which floors
    its residual variance at an ABSORBANCE quantum -- the right thing for a
    progress curve and meaningless for a fit in ln(rate) against 1/T.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        return None
    design = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    dof = max(1, len(x) - 2)
    variance = float(residual @ residual) / dof
    covariance = variance * np.linalg.pinv(design.T @ design)
    return {"intercept": float(beta[0]), "slope": float(beta[1]),
            "intercept_stderr": float(np.sqrt(covariance[0, 0])),
            "slope_stderr": float(np.sqrt(covariance[1, 1])),
            # Kept because DH and DS come from ONE line and are strongly
            # anti-correlated, so DG's error needs the covariance, not the two
            # variances. Dropping it overstates DG's error about tenfold.
            "covariance": covariance,
            "rms": float(np.sqrt(variance)), "n": int(len(x)), "dof": dof}


def arrhenius_fit(kelvin, rate):
    """
    ln(rate) against 1/T. Returns a dict, or None if under three usable points.

    Keys: activation_kJ and its stderr_kJ, intercept, rms, n. The rms is in ln
    units, so it is a relative error: 0.10 means the points scatter about 10%
    around the line.

    THE STANDARD ERROR IS THE POINT of this function existing in this form. The
    four substrate rungs returned activation energies of 94.5, 86.0, 89.0 and
    91.5 kJ/mol, and a spread of 8.5 across four fits that ought to agree looks
    like a composition dependence until the errors are computed -- see
    `rungs_agree`, which asks whether the spread exceeds what six points at this
    scatter can produce.
    """
    kelvin = np.asarray(kelvin, dtype=float)
    rate = np.asarray(rate, dtype=float)
    keep = np.isfinite(rate) & (rate > 0) & np.isfinite(kelvin)
    if keep.sum() < 3:
        return None
    fit = _line(1.0 / kelvin[keep], np.log(rate[keep]))
    if fit is None:
        return None
    fit["activation_kJ"] = -fit["slope"] * GAS_CONSTANT / 1000.0
    fit["stderr_kJ"] = fit["slope_stderr"] * GAS_CONSTANT / 1000.0
    return fit


def rung_fits(parameter="vmax", override=None, experiments=TEMPERATURE_SERIES):
    """
    One Arrhenius fit per substrate rung, on rate divided by [enzyme].

    Dividing by [enzyme] is the whole point: it is what makes runs at different
    enzyme concentrations comparable, and therefore what a wrong enzyme
    concentration corrupts. `override` is an ENZYME_HYPOTHESES entry.

    Returns a DataFrame indexed by the rung's [S], with the activation energy,
    the residual rms, and the number of temperatures that entered.
    """
    override = override or {}
    frame = series_frame(experiments)
    rows = []
    for s0 in sorted(frame.s0.unique()):
        rung = frame[np.isclose(frame.s0, s0)].sort_values("kelvin")
        enzyme = np.array([override.get(int(e), value) for e, value
                           in zip(rung.experiment, rung.e0)], dtype=float)
        fit = arrhenius_fit(rung.kelvin.to_numpy(),
                            rung[parameter].to_numpy() / enzyme)
        fit = fit or {}
        rows.append({"s0": float(s0),
                     "activation_kJ": fit.get("activation_kJ", np.nan),
                     "stderr_kJ": fit.get("stderr_kJ", np.nan),
                     "intercept": fit.get("intercept", np.nan),
                     "rms": fit.get("rms", np.nan),
                     "n": int(len(rung))})
    return pd.DataFrame(rows).set_index("s0")


def enzyme_hypotheses(parameter="vmax", hypotheses=None,
                      experiments=TEMPERATURE_SERIES):
    """
    Which enzyme concentrations make the temperature series most collinear?

    THE LOGIC. Every hypothesis divides the same rates by a different enzyme
    concentration, so they are not different data -- they are the same data
    placed differently on the vertical axis. A concentration that is wrong for
    one run displaces that run's four points together, away from the line the
    other five temperatures define. The hypothesis that leaves the smallest
    residual is the one whose concentrations are right.

    WHAT IT CANNOT DO. It cannot test a factor common to every run, which would
    move the intercept and leave the residual untouched, and it is weak against
    an error at BOTH ends of the range at once -- displacing the two endpoints
    the same way tilts nothing. So "exps16,19 restocked" is tested but is the
    hypothesis this method is least able to reject.

    Returns a DataFrame indexed by hypothesis: mean and worst rung rms, the
    mean activation energy, and how many of the rungs each hypothesis wins.
    """
    hypotheses = hypotheses or ENZYME_HYPOTHESES
    fits = {name: rung_fits(parameter, override, experiments)
            for name, override in hypotheses.items()}
    best = pd.DataFrame({name: table.rms for name, table in fits.items()})
    winner = best.idxmin(axis=1)
    rows = []
    for name, table in fits.items():
        rows.append({
            "hypothesis": name,
            "mean_rms": float(table.rms.mean()),
            "worst_rms": float(table.rms.max()),
            "activation_kJ": float(table.activation_kJ.mean()),
            "rungs_won": int((winner == name).sum()),
            "rungs": int(len(table)),
        })
    return pd.DataFrame(rows).set_index("hypothesis")


def experiment_residuals(parameter="vmax", override=None,
                         experiments=TEMPERATURE_SERIES):
    """
    Each run's mean distance from the Arrhenius line, in ln units.

    More interpretable than the pooled rms when the question is about ONE run:
    a mean of -0.06 says that run sits about 6% below the line its neighbours
    define, and the sign says which way a wrong enzyme concentration pushed it.
    """
    override = override or {}
    frame = series_frame(experiments)
    collected = {}
    for s0 in sorted(frame.s0.unique()):
        rung = frame[np.isclose(frame.s0, s0)].sort_values("kelvin")
        enzyme = np.array([override.get(int(e), value) for e, value
                           in zip(rung.experiment, rung.e0)], dtype=float)
        rate = rung[parameter].to_numpy() / enzyme
        keep = np.isfinite(rate) & (rate > 0)
        if keep.sum() < 3:
            continue
        y = np.log(rate[keep])
        design = np.column_stack([np.ones(int(keep.sum())),
                                  1.0 / rung.kelvin.to_numpy()[keep]])
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        for experiment, value in zip(rung.experiment.to_numpy()[keep],
                                     y - design @ beta):
            collected.setdefault(int(experiment), []).append(float(value))
    rows = []
    for experiment, values in sorted(collected.items()):
        run = frame[frame.experiment == experiment]
        rows.append({"experiment": experiment,
                     "temperature": float(run.temperature.iloc[0]),
                     "e0": float(override.get(experiment, run.e0.iloc[0])),
                     "mean_residual": float(np.mean(values)),
                     "rungs": len(values)})
    return pd.DataFrame(rows).set_index("experiment")


# Where the burst form is trustworthy in this block, and it is not everywhere.
# 15-30 C classify as `lag` with tau resolved on 15 of 16 curves and residuals
# 0.94-1.64x noise; 35 and 40 C flip to `burst` with tau at the top of its grid
# (resolved on 1 of 8) because a decelerating curve has no lag to measure. So
# v_ss is the asymptote at the cold end and a meaningless late rate at the hot
# end, and this is the boundary between them.
BURST_TRUSTWORTHY_BELOW_C = 32.0


def truncation_sensitivity(experiments=TEMPERATURE_SERIES):
    """
    How much does `vmax` running out of run at the cold end inflate E_a?

    THE PROBLEM the segment fit found. Every run from 15 to 30 C is a LAG
    curve, and at 15 and 20 C the steepest block sits at 70-90% of the run --
    the rate is still rising when the measurement stops, so `vmax` there is not
    a maximum, it is wherever the run ended. Under-reading the cold end tilts
    the Arrhenius line and inflates the activation energy.

    THE MEASUREMENT. On a lag curve the burst form's `v_ss` IS the rate `vmax`
    is trying to reach, and at 15 and 20 C that fit is sound -- tau resolved on
    7 of 8 curves, residuals 0.94-1.03x noise. So substituting v_ss for vmax at
    those two temperatures, and only those, bounds the bias.

    MIXING ESTIMATORS IS NORMALLY THE ERROR, not the fix, and it is done here
    ONLY to size the bias, never to produce the headline. Both numbers are
    returned so the difference is the answer rather than the substituted fit.

    Returns a dict: the E_a from vmax throughout, the E_a with v_ss at the cold
    end, their difference, and the per-temperature v_ss/vmax ratios.
    """
    frame = series_frame(experiments)
    cold = frame.temperature < 22.0
    ratios = {float(t): float((group.v_ss / group.vmax).median())
              for t, group in frame.groupby("temperature")}

    def fit(use_v_ss):
        energies = []
        for s0 in sorted(frame.s0.unique()):
            rung = frame[np.isclose(frame.s0, s0)].sort_values("kelvin")
            rate = rung.vmax.to_numpy(dtype=float)
            if use_v_ss:
                swap = (rung.temperature < 22.0).to_numpy()
                rate = np.where(swap, rung.v_ss.to_numpy(dtype=float), rate)
            fit = arrhenius_fit(rung.kelvin.to_numpy(),
                                rate / rung.e0.to_numpy(dtype=float))
            energies.append(fit["activation_kJ"] if fit else np.nan)
        return float(np.mean(energies))

    plain, corrected = fit(False), fit(True)
    return {"vmax_kJ": plain, "cold_corrected_kJ": corrected,
            "inflation_kJ": plain - corrected,
            "cold_curves": int(cold.sum()),
            "v_ss_over_vmax": ratios}


# Boltzmann over Planck, K^-1 s^-1, and its logarithm -- the Eyring intercept's
# fixed part. CODATA 2018, both constants now exact by definition.
BOLTZMANN_OVER_PLANCK = 1.380649e-23 / 6.62607015e-34
LN_BOLTZMANN_OVER_PLANCK = float(np.log(BOLTZMANN_OVER_PLANCK))
# The temperature activation parameters are quoted at, 25 C.
REFERENCE_KELVIN = 298.15


def turnover(frame, parameter="v_ss"):
    """
    A rate in AU/s turned into a pseudo-first-order constant in s^-1.

    v / (epsilon * [enz]): dividing by the sheet's own extinction coefficient
    gives mM of product per second, and dividing that by the enzyme
    concentration in mM leaves s^-1.

    WHAT THIS IS AND IS NOT. It is an apparent turnover frequency at the
    composition of that cuvette, and it is what makes an Eyring ENTROPY
    possible at all: a slope survives any constant factor but an intercept does
    not, so DH can be had from raw AU/s and DS cannot. It is NOT an elementary
    rate constant. The substrate order here is about +0.5, not 0, so the
    reaction is not saturated and this constant carries a substrate dependence;
    quoting DS from it means quoting DS for a pseudo-first-order constant at
    that [S] and [buf]. Say so wherever it is used.
    """
    return (frame[parameter].to_numpy(dtype=float)
            / frame.epsilon.to_numpy(dtype=float)
            / frame.e0.to_numpy(dtype=float))


def eyring_fit(kelvin, rate_constant):
    """
    ln(k/T) against 1/T. Returns a dict, or None if under three usable points.

    Keys: enthalpy_kJ, entropy_J (per mol per K), gibbs_kJ at REFERENCE_KELVIN,
    each with a stderr, plus rms and n.

    `rate_constant` must be an ACTUAL rate constant in s^-1 -- see `turnover`.
    Feeding it AU/s leaves the enthalpy right and the entropy meaningless,
    because only the intercept carries the units.
    """
    kelvin = np.asarray(kelvin, dtype=float)
    rate_constant = np.asarray(rate_constant, dtype=float)
    keep = (np.isfinite(rate_constant) & (rate_constant > 0)
            & np.isfinite(kelvin))
    if keep.sum() < 3:
        return None
    fit = _line(1.0 / kelvin[keep],
                np.log(rate_constant[keep] / kelvin[keep]))
    if fit is None:
        return None
    enthalpy = -fit["slope"] * GAS_CONSTANT / 1000.0
    entropy = GAS_CONSTANT * (fit["intercept"] - LN_BOLTZMANN_OVER_PLANCK)
    enthalpy_stderr = fit["slope_stderr"] * GAS_CONSTANT / 1000.0
    entropy_stderr = GAS_CONSTANT * fit["intercept_stderr"]
    fit.update({
        "enthalpy_kJ": enthalpy, "enthalpy_stderr": enthalpy_stderr,
        "entropy_J": entropy, "entropy_stderr": entropy_stderr,
        "gibbs_kJ": enthalpy - REFERENCE_KELVIN * entropy / 1000.0,
        # DH and DS are strongly anti-correlated -- both come from one line --
        # so propagating their errors independently overstates DG's by roughly
        # an order of magnitude. DG at a temperature INSIDE the measured range
        # is far better determined than either, and this is the fitted value's
        # own error at the reference temperature.
        "gibbs_stderr": _prediction_stderr(fit, 1.0 / REFERENCE_KELVIN)
                        * GAS_CONSTANT * REFERENCE_KELVIN / 1000.0,
    })
    return fit


def _prediction_stderr(fit, x):
    """
    Standard error of the fitted line's value at x, in the fit's own units.

    var(a + b x) = var(a) + x^2 var(b) + 2 x cov(a, b). The cross term is the
    whole reason this exists: at an x INSIDE the measured range the two
    coefficient errors largely cancel, which is why DG at 25 C is far better
    determined than either DH or DS alone.
    """
    covariance = fit.get("covariance")
    if covariance is None:
        return np.nan
    vector = np.array([1.0, float(x)])
    return float(np.sqrt(max(float(vector @ covariance @ vector), 0.0)))


def pooled_arrhenius(parameter="v_ss", experiments=TEMPERATURE_SERIES):
    """
    One activation energy from all four rungs at once, with a free offset each.

    WHY NOT AVERAGE THE FOUR. The rungs differ in composition -- [S] rises
    1.85 -> 7.399 mM while [buf] FALLS 80 -> 50 mM -- so they sit at different
    heights, but that is an intercept, not a slope. Fitting one slope with four
    intercepts uses all 24 points for the temperature dependence and lets each
    rung keep its own level, which is the same within-group design
    `scope.background_orders(within=True)` uses for orders.

    It is also honest about the error in a way that averaging the four separate
    fits is not: those four share the same six runs, the same six days and the
    same six cells, so their errors are correlated and averaging them would
    shrink the error by a factor that is not there.
    """
    frame = series_frame(experiments)
    rate = frame[parameter].to_numpy(dtype=float) / frame.e0.to_numpy(dtype=float)
    keep = np.isfinite(rate) & (rate > 0)
    frame, rate = frame[keep], rate[keep]
    rungs = sorted(frame.s0.unique())
    design = np.column_stack(
        [1.0 / frame.kelvin.to_numpy(dtype=float)]
        + [(np.isclose(frame.s0, s0)).astype(float) for s0 in rungs])
    y = np.log(rate)
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    dof = max(1, len(y) - design.shape[1])
    variance = float(residual @ residual) / dof
    covariance = variance * np.linalg.pinv(design.T @ design)
    return {"activation_kJ": float(-beta[0] * GAS_CONSTANT / 1000.0),
            "stderr_kJ": float(np.sqrt(covariance[0, 0])
                               * GAS_CONSTANT / 1000.0),
            "rms": float(np.sqrt(variance)), "n": int(len(y)),
            "dof": int(dof), "rungs": len(rungs)}


def rungs_agree(parameter="v_ss", experiments=TEMPERATURE_SERIES):
    """
    Is the rung-to-rung spread in activation energy real, or is it scatter?

    Four fits that should give one number gave 94.5, 86.0, 89.0 and 91.5 on
    `vmax`, and a spread of 8.5 kJ/mol reads as a composition dependence until
    the standard errors are computed -- each is about +/- 2.8. This is the
    chi-square of those four against their weighted mean.

    A caution the number cannot carry: the four rungs come from the SAME six
    runs, so their errors are correlated and the test is optimistic. It is
    sound as evidence that the spread is unremarkable, and weak as evidence
    that pooling gains precision -- which is why `pooled_arrhenius` refits
    rather than averaging.
    """
    table = rung_fits(parameter, experiments=experiments)
    values = table.activation_kJ.to_numpy(dtype=float)
    errors = table.stderr_kJ.to_numpy(dtype=float)
    keep = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    values, errors = values[keep], errors[keep]
    weights = 1.0 / errors ** 2
    mean = float((weights * values).sum() / weights.sum())
    chi2 = float((weights * (values - mean) ** 2).sum())
    dof = max(1, len(values) - 1)
    return {"weighted_mean_kJ": mean, "chi2": chi2, "dof": dof,
            "reduced_chi2": chi2 / dof, "spread_kJ": float(np.ptp(values)),
            "median_stderr_kJ": float(np.median(errors)), "n": len(values)}


# What each fitted parameter is, over which temperatures it is trustworthy, and
# why. Declared rather than chosen at the call site, because the whole lesson of
# this block is that no single estimator measures one quantity across 15-40 C.
FITTED_PARAMETERS = {
    "vmax": {
        "label": "steepest observed rate",
        "kelvin_max": None,
        "constant": "turnover",
        "note": "all six temperatures, but truncated at the cold end -- the "
                "steepest block sits at 70-90% of the run at 15 and 20 C, so "
                "this over-reads the activation energy by about 1.9 kJ/mol "
                "(`truncation_sensitivity`)",
    },
    "v_ss": {
        "label": "asymptotic rate after the induction",
        "kelvin_max": BURST_TRUSTWORTHY_BELOW_C + 273.15,
        "constant": "turnover",
        "note": "15-30 C only. Above that the burst form degenerates -- tau "
                "unresolved on 7 of 8 curves -- and v_ss becomes a late "
                "decayed rate rather than an asymptote",
    },
    "inverse_tau": {
        "label": "induction rate constant, 1/tau",
        "kelvin_max": BURST_TRUSTWORTHY_BELOW_C + 273.15,
        "constant": "raw",
        "note": "15-30 C only, same reason. This is a first-order relaxation "
                "rate in s^-1 already, so it needs neither epsilon nor [enz] "
                "-- the one parameter here whose Eyring entropy rests on no "
                "assumption about the rate law",
    },
}


def parameter_values(parameter, frame=None):
    """The parameter's value per curve, and its rate constant in s^-1."""
    frame = series_frame() if frame is None else frame
    limit = FITTED_PARAMETERS[parameter]["kelvin_max"]
    if limit is not None:
        frame = frame[frame.kelvin <= limit]
    if parameter == "inverse_tau":
        # Only where the profile actually pinned tau. An unresolved tau is the
        # optimiser's grid edge, not a measurement.
        frame = frame[frame.tau_resolved & (frame.tau > 0)]
        return frame, 1.0 / frame.tau.to_numpy(dtype=float)
    return frame, turnover(frame, parameter)


def activation_parameters(parameter="vmax", experiments=TEMPERATURE_SERIES):
    """
    Arrhenius and Eyring parameters for one fitted quantity, pooled over rungs.

    The pooled design is one slope with a free intercept per substrate rung --
    see `pooled_arrhenius` for why that rather than averaging four fits. The
    Eyring form is fitted the same way, on ln(k/T).

    Returns a dict: activation_kJ, enthalpy_kJ, entropy_J, gibbs_kJ at 25 C,
    each with a standard error, plus the rms, the curve count and the note from
    FITTED_PARAMETERS saying where the number may be trusted.
    """
    frame, constant = parameter_values(parameter, series_frame(experiments))
    keep = np.isfinite(constant) & (constant > 0)
    frame, constant = frame[keep], constant[keep]
    kelvin = frame.kelvin.to_numpy(dtype=float)
    rungs = sorted(frame.s0.unique())
    indicators = [(np.isclose(frame.s0, s0)).astype(float) for s0 in rungs]

    def pooled(y):
        design = np.column_stack([1.0 / kelvin] + indicators)
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        residual = y - design @ beta
        dof = max(1, len(y) - design.shape[1])
        variance = float(residual @ residual) / dof
        covariance = variance * np.linalg.pinv(design.T @ design)
        return beta, covariance, float(np.sqrt(variance)), dof

    beta, covariance, rms, dof = pooled(np.log(constant))
    energy = -beta[0] * GAS_CONSTANT / 1000.0
    energy_stderr = float(np.sqrt(covariance[0, 0])) * GAS_CONSTANT / 1000.0

    # Eyring, on the same design. The per-rung intercepts mean the entropy is
    # reported for the MEDIAN rung: an intercept is a level, and the four rungs
    # sit at four levels because their composition differs.
    beta_e, covariance_e, rms_e, dof_e = pooled(np.log(constant / kelvin))
    enthalpy = -beta_e[0] * GAS_CONSTANT / 1000.0
    enthalpy_stderr = float(np.sqrt(covariance_e[0, 0])) * GAS_CONSTANT / 1000.0
    # The entropy is an INTERCEPT, and the four rungs have four intercepts
    # because their composition differs, so it has to be quoted AT a
    # composition. The median rung is used and named in the returned dict --
    # never left implicit, because a DS with no composition attached is not a
    # number anyone can compare against a calculation.
    middle = 1 + len(rungs) // 2
    entropy = GAS_CONSTANT * (beta_e[middle] - LN_BOLTZMANN_OVER_PLANCK)
    entropy_stderr = GAS_CONSTANT * float(np.sqrt(covariance_e[middle, middle]))
    # DG at a temperature inside the measured range, where the anti-correlated
    # errors in DH and DS largely cancel. Propagated with the covariance rather
    # than from the two standard errors, which would overstate it several-fold.
    vector = np.zeros(len(beta_e))
    vector[0] = 1.0 / REFERENCE_KELVIN
    vector[middle] = 1.0
    predicted_stderr = float(np.sqrt(max(
        float(vector @ covariance_e @ vector), 0.0)))
    return {
        "parameter": parameter,
        "label": FITTED_PARAMETERS[parameter]["label"],
        "activation_kJ": energy, "activation_stderr": energy_stderr,
        "enthalpy_kJ": enthalpy, "enthalpy_stderr": enthalpy_stderr,
        "entropy_J": entropy, "entropy_stderr": entropy_stderr,
        "gibbs_kJ": enthalpy - REFERENCE_KELVIN * entropy / 1000.0,
        "gibbs_stderr": predicted_stderr * GAS_CONSTANT * REFERENCE_KELVIN / 1000.0,
        "at_s0": float(rungs[middle - 1]),
        "at_buf": float(frame[np.isclose(frame.s0, rungs[middle - 1])]
                        .buf.median()),
        "rms": rms, "n": int(len(constant)),
        "temperatures": int(frame.temperature.nunique()),
        "note": FITTED_PARAMETERS[parameter]["note"],
    }


def parameter_table(parameters=None, experiments=TEMPERATURE_SERIES):
    """Activation parameters for every fitted quantity, one row each."""
    parameters = parameters or list(FITTED_PARAMETERS)
    return pd.DataFrame([activation_parameters(name, experiments)
                         for name in parameters]).set_index("parameter")


# The catalysed buffer titration that lets the temperature series' substrate
# order be corrected: 4OMe-BnOH, 40 C, pH 7.00, [S] fixed at 8.251 mM, [buf]
# stepped. Exp 32 covers 50-200 mM and exp 34 covers 3.125-25 mM, and they are
# kept apart because the order is NOT constant across that range -- see
# `catalysed_buffer_order`.
BUFFER_TITRATION_HIGH = (32,)     # 50-200 mM, the range the series sits in
BUFFER_TITRATION_LOW = (34,)      # 3.125-25 mM


def catalysed_buffer_order(experiments=BUFFER_TITRATION_HIGH,
                           parameter="vmax"):
    """
    How the CATALYSED rate depends on buffer, at fixed substrate.

    The temperature series cannot measure this: its own [buf] moves only as a
    by-product of the substrate ladder. Exps 32 and 34 can, because they step
    [buf] with everything else held -- the design the BnOH set lacks entirely.

    IT IS NOT ONE NUMBER. Over 50-200 mM (exp 32) the order is about +0.40 with
    R2 = 0.99; over 3.125-25 mM (exp 34) it is about +0.80. The buffer
    dependence SATURATES, so the value to use depends on where you are, and the
    temperature series sits at 50-80 mM -- the high range.

    Use `vmax`, not `v_ss`: three of exp 34's four curves have an unresolved
    tau and return a negative v_ss, which turns the order into nonsense
    (-0.18 +/- 0.34).
    """
    return scope.background_orders(tuple(experiments), terms=("buf",),
                                   within=True, parameter=parameter)


def substrate_order(parameter="vmax", experiments=TEMPERATURE_SERIES):
    """
    The substrate order per temperature, observed and corrected for buffer.

    THE CONFOUND. The ladder is not a substrate ladder alone: [buf] falls
    80 -> 70 -> 60 -> 50 mM as [S] rises 1.850 -> 7.399 mM, because substrate
    volume displaced buffer volume in the cuvette. It is the same
    volume-displacement confound as the BnOH titrations (exps 3 and 6) and the
    4OMe background block, and it is IDENTICAL in all six runs -- which is why
    it leaves the activation energies alone (each rung is one fixed composition
    across all six temperatures) and corrupts only the substrate order.

    THE CORRECTION. If the rate goes as [S]^a [buf]^d and within the ladder
    log[buf] = g log[S] + constant, then a fit without a buffer term returns
    a' = a + d g, so a = a' - d g. Here g is about -0.339 and d comes from
    `catalysed_buffer_order` on the high range, so the correction ADDS about
    0.14 to every observed order.

    WHAT IT RESTS ON, stated because it is not testable here: that the buffer
    order measured at 40 C holds at 15 C, and that it holds at [S] = 1.85-7.4
    mM as well as at the 8.251 mM it was measured at.

    Returns a DataFrame indexed by temperature with the observed order, the
    corrected order, and the shared coupling and buffer order used.
    """
    frame = series_frame(experiments)
    # g, from the pipetted concentrations rather than from a fitted rate.
    coupling = float(np.polyfit(np.log(frame.s0.to_numpy(dtype=float)),
                                np.log(frame.buf.to_numpy(dtype=float)), 1)[0])
    buffer = catalysed_buffer_order(parameter=parameter)
    rows = []
    for temperature, group in frame.groupby("temperature"):
        group = group.sort_values("s0")
        observed = float(np.polyfit(np.log(group.s0.to_numpy(dtype=float)),
                                    np.log(group[parameter].to_numpy(dtype=float)),
                                    1)[0])
        rows.append({"temperature": float(temperature),
                     "observed": observed,
                     "corrected": observed - buffer["order_buf"] * coupling})
    table = pd.DataFrame(rows).set_index("temperature")
    table.attrs["coupling"] = coupling
    table.attrs["buffer_order"] = buffer["order_buf"]
    table.attrs["buffer_stderr"] = buffer["stderr_buf"]
    return table
