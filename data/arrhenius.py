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


def arrhenius_fit(kelvin, rate):
    """
    ln(rate) against 1/T: returns (activation energy kJ/mol, intercept, rms).

    The rms is of the residual in ln units, which is a relative error: 0.10
    means the points scatter about 10% around the line. That is the number the
    enzyme test compares, because a wrong enzyme concentration displaces one
    point vertically and nothing else.
    """
    kelvin = np.asarray(kelvin, dtype=float)
    rate = np.asarray(rate, dtype=float)
    keep = np.isfinite(rate) & (rate > 0) & np.isfinite(kelvin)
    if keep.sum() < 3:
        return np.nan, np.nan, np.nan
    y = np.log(rate[keep])
    design = np.column_stack([np.ones(int(keep.sum())), 1.0 / kelvin[keep]])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    rms = float(np.sqrt(residual @ residual / max(1, len(y) - 2)))
    return float(-beta[1] * GAS_CONSTANT / 1000.0), float(beta[0]), rms


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
        energy, intercept, rms = arrhenius_fit(
            rung.kelvin.to_numpy(), rung[parameter].to_numpy() / enzyme)
        rows.append({"s0": float(s0), "activation_kJ": energy,
                     "intercept": intercept, "rms": rms,
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
            energy, _, _ = arrhenius_fit(rung.kelvin.to_numpy(),
                                         rate / rung.e0.to_numpy(dtype=float))
            energies.append(energy)
        return float(np.mean(energies))

    plain, corrected = fit(False), fit(True)
    return {"vmax_kJ": plain, "cold_corrected_kJ": corrected,
            "inflation_kJ": plain - corrected,
            "cold_curves": int(cold.sum()),
            "v_ss_over_vmax": ratios}
