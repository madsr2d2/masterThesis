"""
Tests for solution_chemistry.py.

Two kinds of test live here. The first kind pins the physics against values
that can be worked out by hand, so the module has a validation gate rather
than just self-consistency. The second kind is a regression test for each of
the four bugs the notebook version carried -- each is written so that it fails
if the old behaviour ever comes back.

    python data/test_solution_chemistry.py
"""
import math
import sys

import numpy as np
import pandas as pd

from solution_chemistry import (
    BUFFER_CHARGES, BUFFER_PKA, DEBYE_HUCKEL_A, DEBYE_HUCKEL_B, PKA_H2O2,
    add_solution_columns, effective_pka_h2o2, hydroperoxide,
    hydroperoxide_fraction, ionic_strength, speciation,
)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def close(a, b, tol=1e-6):
    return abs(float(a) - float(b)) <= tol


# --- tables ---------------------------------------------------------------

def test_tables():
    print("\ntables")
    for buffer_name, pKa_values in BUFFER_PKA.items():
        check(f"{buffer_name}: charges are one longer than pKa list",
              len(BUFFER_CHARGES[buffer_name]) == len(pKa_values) + 1)
        check(f"{buffer_name}: pKa values ascend",
              list(pKa_values) == sorted(pKa_values))
        check(f"{buffer_name}: charges descend by one from zero",
              BUFFER_CHARGES[buffer_name] == list(range(0, -len(pKa_values) - 1, -1)))
    check("buffer names match the extraction's output",
          set(BUFFER_PKA) == {"Phosphate", "Pyrophosphate", "Boric", "Carbonate"})


# --- speciation -----------------------------------------------------------

def test_speciation():
    print("\nspeciation")
    for buffer_name in BUFFER_PKA:
        for pH in (5.0, 7.0, 9.0, 11.0):
            alpha = speciation(buffer_name, pH)
            check(f"{buffer_name} pH {pH}: fractions sum to 1",
                  close(alpha.sum(), 1.0, 1e-12))
            check(f"{buffer_name} pH {pH}: all fractions non-negative",
                  bool((alpha >= 0).all()))

    # A monoprotic acid at pH = pKa is exactly half dissociated.
    alpha = speciation("Boric", BUFFER_PKA["Boric"][0])
    check("Boric at pH = pKa is 50/50", close(alpha[0], 0.5, 1e-12) and close(alpha[1], 0.5, 1e-12))

    # At pH = pKa2, phosphate's two middle species are equal.
    alpha = speciation("Phosphate", BUFFER_PKA["Phosphate"][1])
    check("Phosphate at pH = pKa2 has H2PO4- == HPO4(2-)", close(alpha[1], alpha[2], 1e-12))

    # Raising the pH must shift population towards the more deprotonated end.
    low = speciation("Phosphate", 6.0)
    high = speciation("Phosphate", 8.0)
    mean_low = float(np.dot(low, np.arange(len(low))))
    mean_high = float(np.dot(high, np.arange(len(high))))
    check("higher pH means more deprotonated on average", mean_high > mean_low,
          f"{mean_low:.3f} -> {mean_high:.3f}")


# --- ionic strength -------------------------------------------------------

def test_ionic_strength():
    print("\nionic strength")
    # VALIDATION GATE, worked by hand and independent of this module:
    #   100 mM phosphate at pH 7.00, pKa2 = 7.20
    #   [HPO4(2-)]/[H2PO4-] = 10^(7.00-7.20) = 0.6310
    #   fractions 0.61313 and 0.38684 (H3PO4 and PO4(3-) negligible)
    #   Na+ needed = 0.61313 + 2*0.38684 = 1.38682 per phosphate
    #   I = 1/2 * 100 * (0.61313*1 + 0.38684*4 + 1.38682*1) = 177.37 mM
    value = ionic_strength("Phosphate", 7.00, 100.0)
    check("100 mM phosphate pH 7.00 gives I = 177.37 mM",
          close(value, 177.37, 0.02), f"got {value:.4f}")

    # Same check for the simplest case: a monoprotic buffer half dissociated.
    #   I = 1/2 * 100 * (0.5*1 + 0.5*1) = 50 mM
    value = ionic_strength("Boric", BUFFER_PKA["Boric"][0], 100.0)
    check("100 mM boric at pH = pKa gives I = 50 mM", close(value, 50.0, 1e-9),
          f"got {value:.6f}")

    check("ionic strength is linear in buffer concentration",
          close(ionic_strength("Phosphate", 7.0, 200.0),
                2 * ionic_strength("Phosphate", 7.0, 100.0), 1e-9))
    check("zero buffer gives zero ionic strength",
          close(ionic_strength("Phosphate", 7.0, 0.0), 0.0, 1e-12))

    rising = [ionic_strength("Phosphate", pH, 100.0) for pH in (5, 6, 7, 8, 9)]
    check("ionic strength rises with pH for phosphate",
          all(b > a for a, b in zip(rising, rising[1:])),
          " -> ".join(f"{v:.1f}" for v in rising))


# --- Debye-Huckel and [HOO-] ---------------------------------------------

def test_debye_huckel():
    print("\nDebye-Huckel and [HOO-]")
    check("at I = 0 the pKa is unshifted", close(effective_pka_h2o2(0.0), PKA_H2O2, 1e-12))
    check("ionic strength always lowers the pKa",
          effective_pka_h2o2(75.0) < PKA_H2O2)

    # Hand check at I = 75 mM: sqrt(0.075) = 0.273861
    #   shift = 0.509 * 2 * 0.273861 / (1 + 0.328 * 0.273861) = 0.25581
    expected = PKA_H2O2 - 0.25581
    check("pKa at I = 75 mM matches the hand calculation",
          close(effective_pka_h2o2(75.0), expected, 1e-4),
          f"got {effective_pka_h2o2(75.0):.5f}, expected {expected:.5f}")

    check("HOO- fraction is in [0, 1]",
          all(0.0 <= hydroperoxide_fraction(pH, 100.0) <= 1.0
              for pH in (4, 7, 9, 11, 12, 14)))
    check("HOO- fraction rises with pH",
          hydroperoxide_fraction(9.0, 100.0) < hydroperoxide_fraction(11.0, 100.0))
    check("at pH = pKa_eff the fraction is exactly one half",
          close(hydroperoxide_fraction(effective_pka_h2o2(75.0), 75.0), 0.5, 1e-12))
    check("[HOO-] is linear in total peroxide",
          close(hydroperoxide(9.0, 200.0, 75.0), 2 * hydroperoxide(9.0, 100.0, 75.0), 1e-9))
    check("[HOO-] never exceeds total peroxide",
          all(hydroperoxide(pH, 50.0, 75.0) <= 50.0 for pH in (5, 8, 11, 14)))


# --- regressions on the four notebook bugs -------------------------------

def test_regressions():
    print("\nregressions on the notebook bugs")

    # BUG 1: an unrecognised buffer name silently produced I = 0. The name
    # mismatch was 'Boric Acid' vs 'Boric'; the fix is to refuse, not default.
    try:
        speciation("Boric Acid", 9.0)
        check("unknown buffer name raises instead of returning I = 0", False,
              "no exception raised")
    except KeyError:
        check("unknown buffer name raises instead of returning I = 0", True)
    check("Carbonate is present in the table (it was missing entirely)",
          "Carbonate" in BUFFER_PKA and ionic_strength("Carbonate", 9.4, 100.0) > 0)

    # BUG 2: only the anions were counted, omitting the sodium counter-ions.
    alpha = speciation("Phosphate", 7.0)
    charges = np.asarray(BUFFER_CHARGES["Phosphate"], dtype=float)
    anion_only = 0.5 * 100.0 * float(np.sum(alpha * charges ** 2))
    check("counter-ions are included, not just anions",
          ionic_strength("Phosphate", 7.0, 100.0) > anion_only * 1.5,
          f"anion-only would give {anion_only:.2f} mM")

    # BUG 3: I was passed to Debye-Huckel in mM while A assumes mol/L.
    molar_mistake = (DEBYE_HUCKEL_A * 2 * math.sqrt(75.0)
                     / (1 + DEBYE_HUCKEL_B * math.sqrt(75.0)))
    actual = PKA_H2O2 - effective_pka_h2o2(75.0)
    check("ionic strength is converted mM -> mol/L before Debye-Huckel",
          actual < molar_mistake / 5,
          f"shift {actual:.4f}, the mM-fed mistake gives {molar_mistake:.4f}")

    # BUG 4: delta(z^2) is 2 for a neutral acid, not 1.
    root = math.sqrt(0.075)
    single = DEBYE_HUCKEL_A * 1 * root / (1 + DEBYE_HUCKEL_B * root)
    check("charge factor is delta(z^2) = 2, so the shift is twice the z^2 = 1 value",
          close(actual, 2 * single, 1e-9),
          f"shift {actual:.5f}, z^2=1 would give {single:.5f}")


# --- against the real dataset --------------------------------------------

def test_dataset(path="data/experiment_data.csv"):
    print("\nagainst the compiled dataset")
    try:
        data = pd.read_csv(path)
    except FileNotFoundError:
        print(f"  skip  {path} not found")
        return

    out = add_solution_columns(data)
    check("every buffer in the dataset has pKa data",
          set(data.buffer.unique()) <= set(BUFFER_PKA),
          f"unknown: {set(data.buffer.unique()) - set(BUFFER_PKA)}")
    check("I is finite and positive for every row",
          bool(np.isfinite(out["I"]).all() and (out["I"] > 0).all()))
    check("[HOO-] is finite and non-negative for every row",
          bool(np.isfinite(out["[HOO-]"]).all() and (out["[HOO-]"] >= 0).all()))
    check("[HOO-] never exceeds [h2o2]",
          bool((out["[HOO-]"] <= out["[h2o2]"] + 1e-12).all()))
    check("I is constant within each experiment where [buf] is",
          all(g["I"].nunique() == g["[buf]"].nunique()
              for _, g in out.groupby("experiment")))

    original = add_solution_columns(data)
    check("the calculation is deterministic",
          bool(np.allclose(original["I"], out["I"])
               and np.allclose(original["[HOO-]"], out["[HOO-]"])))

    lo, hi = out["I"].min(), out["I"].max()
    print(f"        I spans {lo:.1f} - {hi:.1f} mM; "
          f"[HOO-] spans {out['[HOO-]'].min():.3g} - {out['[HOO-]'].max():.3g} mM")


if __name__ == "__main__":
    test_tables()
    test_speciation()
    test_ionic_strength()
    test_debye_huckel()
    test_regressions()
    test_dataset()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
