"""
Solution-chemistry quantities derived from the recorded conditions: buffer
speciation, ionic strength, and the hydroperoxide anion concentration.

These are the only columns in the dataset that are *computed* rather than
measured or read off a sheet, which makes them the easiest place for an error
to hide. Every one of the following lived in a notebook cell until now, and
every one was wrong at some point without anything noticing:

  * the buffer table was keyed 'Boric Acid' while find_buffer_type returns
    'Boric', and Carbonate was missing entirely, so 107 rows across 26
    experiments fell through to the unknown-buffer fallback and were silently
    assigned I = 0 -- precisely the high-pH runs where the correction matters
    most;
  * ionic strength counted only the buffer anions and not the sodium
    counter-ions required for electroneutrality, understating I by about a
    factor of two;
  * the Debye-Huckel term was fed I in mM while A = 0.509 is defined for
    mol/L, inflating sqrt(I) by ~31x and the pKa shift by ~9x;
  * the charge factor used z^2 = 1 rather than delta(z^2) = 2 for a neutral
    acid dissociating into two singly charged ions.

All four are fixed here and pinned by tests in test_solution_chemistry.py.
The physics is otherwise identical to the corrected notebook version, so
stamping these columns does not change any number that has been looked at.

UNITS. Concentrations throughout the dataset are in mM, and so are the inputs
and outputs here. The Debye-Huckel constants require mol/L; that conversion
happens inside effective_pka_h2o2 and nowhere else. No caller should ever pass
or receive an ionic strength in mol/L.

KNOWN LIMITATIONS -- documented rather than modelled, because the dataset does
not record what would be needed to do better:

  1. Mixed buffers. Several sheets prepare one "buffer" from two salts (exp 146
     mixes Na4P2O7*10H2O with NaH2PO4*2H2O), but the dataset carries a single
     buffer name and a single [buf]. Speciation is computed for the named
     buffer alone. For a pyrophosphate/phosphate mixture near pH 8.7 this
     misstates I by roughly 10-20%.
  2. Titrant. HCl or NaOH used to bring the buffer to its target pH adds to I
     and is not recorded anywhere.
  3. Thermodynamic pKa values. The buffer pKa values below are zero-ionic-
     strength values used at finite I, while H2O2's pKa *is* activity-
     corrected. A self-consistent treatment would correct both and iterate.
     The inconsistency is deliberate: correcting only the quantity that feeds
     the [HOO-] calculation keeps this module numerically identical to the
     results already recorded in MECHANISM.md and DATA_VERIFICATION.md.
  4. Temperature. A = 0.509 and the pKa values are 25 C values, applied across
     the dataset's 15-40 C range. The Debye-Huckel constant varies by about 5%
     over that span (~0.004 pKa units), which is negligible; the temperature
     dependence of pKa(H2O2) itself is not.
  5. Validity range. Extended Debye-Huckel is reliable to about I = 100 mM;
     70% of this dataset's live rows exceed that and the pyrophosphate runs
     reach 1069 mM, an order of magnitude beyond it. The default was therefore
     changed to Davies on 2026-08-31, which is claimed to roughly 500 mM. The
     two differ by up to 0.60 pKa units, a actor of 4 in [HOO-], at the top
     of the range.

     Davies has its own ceiling, and it is worth stating plainly rather than
     trading one silent extrapolation for another. Its -0.3I term eventually
     dominates, so the predicted shift turns around near I = 0.5 mol/L:

         I (mM)      100     400     700    1069
         pKa_eff  11.536  11.478  11.500  11.559

     Above ~500 mM the curve is no longer physical, it is merely bounded. What
     changed is that the error stops growing without limit; the highest-I rows
     are still approximate. out_of_range_fraction and the validator continue to
     report them, and a Pitzer treatment with real ion-interaction parameters
     remains the principled fix. effective_pka_h2o2(..., model="debye")
     reproduces the previous numbers exactly, so any earlier result can be
     regenerated.

Usage:
    from solution_chemistry import add_solution_columns
    df = add_solution_columns(df)          # adds 'I' (mM) and '[HOO-]' (mM)

    python data/solution_chemistry.py      # recompute and compare against the
                                           # stored columns in the dataset
"""

import numpy as np
import pandas as pd

# Dissociation constants at 25 C and zero ionic strength. Keys MUST match the
# buffer names produced by kinetics_io.find_buffer_type.
BUFFER_PKA = {
    "Phosphate": [2.15, 7.20, 12.35],
    "Pyrophosphate": [0.85, 1.96, 6.60, 9.41],
    "Boric": [9.24],
    "Carbonate": [6.35, 10.33],
}

# Charge of each species, fully protonated first. Must be one longer than the
# pKa list: n dissociation steps give n+1 species.
BUFFER_CHARGES = {
    "Phosphate": [0, -1, -2, -3],
    "Pyrophosphate": [0, -1, -2, -3, -4],
    "Boric": [0, -1],
    "Carbonate": [0, -1, -2],
}

PKA_H2O2 = 11.75  # H2O2 <-> H+ + HOO-, 25 C, zero ionic strength
DEBYE_HUCKEL_A = 0.509  # mol/L^-0.5, water at 25 C
DEBYE_HUCKEL_B = 0.328  # extended Debye-Huckel denominator constant
DELTA_Z_SQUARED = 2  # (+1)^2 + (-1)^2 - 0^2 for a neutral acid dissociating

# Ionic strength beyond which the extended Debye-Huckel equation is being used
# outside its established range. Not an error -- there is no better model here
# without activity data -- but every result that depends on it should say so.
DEBYE_HUCKEL_RELIABLE_mM = 100.0

# Davies' empirical term, the standard 0.3 in
#     log(gamma) = -A z^2 (sqrt(I)/(1 + sqrt(I)) - 0.3 I)
# which extends usable range to roughly 500 mM against Debye-Huckel's 100.
DAVIES_B = 0.3

# Which activity model effective_pka_h2o2 uses by default. Changed from
# "debye" to "davies" on 2026-08-31: 70% of the dataset's rows exceed
# DEBYE_HUCKEL_RELIABLE_mM, and the two models differ by up to a factor of 4
# in [HOO-] there. See DATA_VERIFICATION.md 2026-08-31.
ACTIVITY_MODEL = "davies"


def speciation(buffer_name, pH):
    """
    Fractional abundance of each buffer species at a given pH.

    Uses the cumulative form: the ratio of species i to the fully protonated
    form is the product of the first i dissociation constants divided by
    [H+]^i. Returned fractions sum to 1 and are ordered fully-protonated
    first, matching BUFFER_CHARGES.

    Args:
        buffer_name (str): One of the keys of BUFFER_PKA.
        pH (float): Solution pH.

    Returns:
        numpy.ndarray: Fractions, length len(BUFFER_PKA[buffer_name]) + 1.

    Raises:
        KeyError: If the buffer is not in BUFFER_PKA. This is deliberate --
            the original implementation returned a single-species fallback
            that silently produced I = 0 for every unrecognised name, which is
            how the 'Boric Acid'/'Boric' key mismatch survived unnoticed.
    """
    if buffer_name not in BUFFER_PKA:
        raise KeyError(
            f"no pKa data for buffer {buffer_name!r}; known buffers are "
            f"{sorted(BUFFER_PKA)}. Add it to BUFFER_PKA and BUFFER_CHARGES "
            f"rather than letting it fall through to a default."
        )
    proton = 10.0 ** (-float(pH))
    ratios = [1.0]
    for pKa in BUFFER_PKA[buffer_name]:
        ratios.append(ratios[-1] * (10.0 ** (-pKa) / proton))
    ratios = np.asarray(ratios, dtype=float)
    return ratios / ratios.sum()


def ionic_strength(buffer_name, pH, buffer_mM):
    """
    Ionic strength of the buffer solution, in mM.

    I = 1/2 * sum_i c_i z_i^2, summed over the buffer anions AND the sodium
    counter-ions needed for electroneutrality. The counter-ion count is taken
    as the mean anion charge at this pH, which is exact for a buffer prepared
    as a mixture of its sodium salts and titrated by adjusting that ratio --
    the preparation the sheets describe.

    Args:
        buffer_name (str): One of the keys of BUFFER_PKA.
        pH (float): Solution pH.
        buffer_mM (float): Total buffer concentration in mM.

    Returns:
        float: Ionic strength in mM.
    """
    alpha = speciation(buffer_name, pH)
    charges = np.asarray(BUFFER_CHARGES[buffer_name], dtype=float)
    if len(charges) != len(alpha):
        raise ValueError(
            f"{buffer_name}: {len(BUFFER_PKA[buffer_name])} pKa values imply "
            f"{len(alpha)} species but {len(charges)} charges are declared"
        )
    anion_term = float(np.sum(alpha * charges**2))
    counter_ion_term = float(np.sum(alpha * np.abs(charges)))
    return 0.5 * float(buffer_mM) * (anion_term + counter_ion_term)


def effective_pka_h2o2(ionic_strength_mM, model=None):
    """
    pKa of H2O2 corrected for ionic strength.

    Two models, both giving pKa_eff = pKa0 - A * delta(z^2) * f(I) with I in
    mol/L. This is the only place the mM -> mol/L conversion happens.

        davies   f = sqrt(I) / (1 + sqrt(I)) - 0.3 I      (the default)
        debye    f = sqrt(I) / (1 + B sqrt(I))

    Davies is the default because 70% of the dataset's rows sit above 100 mM,
    where extended Debye-Huckel is outside its established range, and the two
    diverge by up to 0.60 pKa units -- a factor of 4 in [HOO-] at the dataset's
    maximum of 1069 mM. Davies is itself only claimed to ~500 mM, so the
    highest-I rows remain approximate; what changes is that the approximation
    no longer degrades without bound. out_of_range_fraction still reports them.

    Args:
        ionic_strength_mM (float): Ionic strength in mM.
        model (str): "davies" or "debye"; None uses ACTIVITY_MODEL.

    Returns:
        float: Effective pKa. Equals PKA_H2O2 exactly at I = 0.

    Raises:
        ValueError: On a negative ionic strength or an unknown model.
    """
    # Resolved at call time, not bound as a default, so ACTIVITY_MODEL can be
    # switched by a caller comparing the two treatments.
    model = ACTIVITY_MODEL if model is None else model
    molar = np.asarray(ionic_strength_mM, dtype=float) / 1000.0
    if np.any(molar < 0):
        raise ValueError(f"negative ionic strength: {ionic_strength_mM} mM")
    root = np.sqrt(molar)
    if model == "davies":
        term = root / (1.0 + root) - DAVIES_B * molar
    elif model == "debye":
        term = root / (1.0 + DEBYE_HUCKEL_B * root)
    else:
        raise ValueError(f"unknown activity model: {model!r}")
    shift = DEBYE_HUCKEL_A * DELTA_Z_SQUARED * term
    result = PKA_H2O2 - shift
    return float(result) if np.ndim(ionic_strength_mM) == 0 else result


def hydroperoxide_fraction(pH, ionic_strength_mM):
    """
    Fraction of total peroxide present as HOO- (Henderson-Hasselbalch).

    Args:
        pH (float): Solution pH.
        ionic_strength_mM (float): Ionic strength in mM.

    Returns:
        float: Fraction in [0, 1].
    """
    ratio = 10.0 ** (float(pH) - effective_pka_h2o2(ionic_strength_mM))
    return ratio / (1.0 + ratio)


def hydroperoxide(pH, h2o2_mM, ionic_strength_mM):
    """
    Hydroperoxide anion concentration, in mM.

    Args:
        pH (float): Solution pH.
        h2o2_mM (float): Total peroxide concentration in mM.
        ionic_strength_mM (float): Ionic strength in mM.

    Returns:
        float: [HOO-] in mM.
    """
    return hydroperoxide_fraction(pH, ionic_strength_mM) * float(h2o2_mM)


def out_of_range_fraction(ionic_strengths_mM, limit=DEBYE_HUCKEL_RELIABLE_mM):
    """
    Fraction of the given ionic strengths that exceed the Debye-Huckel range.

    Args:
        ionic_strengths_mM (array-like): Ionic strengths in mM.
        limit (float): Threshold in mM, default DEBYE_HUCKEL_RELIABLE_mM.

    Returns:
        float: Fraction in [0, 1]. Zero for an empty input.
    """
    values = np.asarray(ionic_strengths_mM, dtype=float)
    if values.size == 0:
        return 0.0
    return float((values > limit).mean())


def add_solution_columns(
    dataframe,
    buffer_col="buffer",
    ph_col="pH",
    buffer_conc_col="[buf]",
    h2o2_col="[h2o2]",
):
    """
    Returns a copy of the dataframe with 'I' and '[HOO-]' columns added.

    Both are in mM, like every other concentration in the dataset.

    Args:
        dataframe (pd.DataFrame): Must carry the four named columns.
        buffer_col, ph_col, buffer_conc_col, h2o2_col (str): Column names.

    Returns:
        pd.DataFrame: A copy with 'I' and '[HOO-]' appended.
    """
    missing = [
        c
        for c in (buffer_col, ph_col, buffer_conc_col, h2o2_col)
        if c not in dataframe.columns
    ]
    if missing:
        raise KeyError(f"dataframe is missing required columns: {missing}")

    out = dataframe.copy()
    out["I"] = [
        ionic_strength(b, p, c)
        for b, p, c in zip(out[buffer_col], out[ph_col], out[buffer_conc_col])
    ]
    out["[HOO-]"] = [
        hydroperoxide(p, h, i) for p, h, i in zip(out[ph_col], out[h2o2_col], out["I"])
    ]
    return out


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Recompute I and [HOO-] and compare against stored columns."
    )
    parser.add_argument("--csv", default="data/experiment_data.csv")
    args = parser.parse_args()

    data = pd.read_csv(args.csv)
    recomputed = add_solution_columns(data)

    print(f"{len(recomputed)} rows, {recomputed.experiment.nunique()} experiments\n")
    summary = recomputed.groupby("buffer").agg(
        rows=("I", "size"),
        pH_min=("pH", "min"),
        pH_max=("pH", "max"),
        I_min=("I", "min"),
        I_max=("I", "max"),
        HOO_min=("[HOO-]", "min"),
        HOO_max=("[HOO-]", "max"),
    )
    print(summary.round(4).to_string())

    over = recomputed["I"] > DEBYE_HUCKEL_RELIABLE_mM
    print(
        f"\nDebye-Huckel range: {int(over.sum())} of {len(recomputed)} rows "
        f"({out_of_range_fraction(recomputed['I']):.0%}) exceed "
        f"{DEBYE_HUCKEL_RELIABLE_mM:.0f} mM, across "
        f"{recomputed.loc[over, 'experiment'].nunique()} experiments; "
        f"the maximum is {recomputed['I'].max():.0f} mM"
    )

    for column in ("I", "[HOO-]"):
        if column not in data.columns:
            print(f"\n{column}: not stored in {args.csv} (nothing to compare)")
            continue
        delta = (recomputed[column] - data[column]).abs()
        worst = delta.max()
        print(f"\n{column}: stored column present, max |difference| = {worst:.3e}")
