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
    ACTIVITY_MODEL, BUFFER_CHARGES, BUFFER_PKA, DAVIES_B, DEBYE_HUCKEL_A,
    DEBYE_HUCKEL_B, PKA_H2O2,
    add_solution_columns, dominant_buffer_pair, effective_pka_h2o2,
    hydroperoxide, hydroperoxide_fraction, ionic_strength, speciation,
    OXYGEN_SOLUBILITY_mM, MOLAR_VOLUME_uL_per_umol, oxygen_budget,
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


def test_dominant_buffer_pair():
    print("\ndominant buffer pair")
    # At pH = pKa the pair is half and half, and the pair chosen is that pKa's.
    acid, base, pKa = dominant_buffer_pair("Phosphate", 7.20, 100.0)
    check("at pH = pKa2 phosphate's pair is equal", close(acid, base, 1e-12),
          f"{acid:.4f}/{base:.4f}")
    check("and the pair reported is pKa2", close(pKa, 7.20, 1e-12))
    # NOT 50.0 each. The pair is the dominant one, not the whole buffer: H3PO4
    # and PO4(3-) hold the remaining 8e-5 of it here. That is the function
    # behaving correctly -- a pair is a pair -- and a caller that needs the
    # total must use `buf`, which is why scope.frame carries both.
    check("the pair is very nearly, but not exactly, the whole buffer",
          99.9 < acid + base < 100.0, f"{acid + base:.4f} of 100")

    # The pair straddling the pH is chosen, not the first one.
    _, _, pKa_low = dominant_buffer_pair("Phosphate", 2.5, 100.0)
    check("at pH 2.5 phosphate reports pKa1", close(pKa_low, 2.15, 1e-12))

    # A pair never exceeds the total, and at a pH far from every pKa it is
    # nearly all of it -- the species that matter are the ones near the pH.
    acid, base, _ = dominant_buffer_pair("Phosphate", 8.01, 85.0)
    check("the pair does not exceed the total", acid + base <= 85.0 + 1e-9,
          f"{acid + base:.4f}")
    # 10^(8.01 - 7.20) = 6.46, which is what this must equal -- an earlier
    # draft of this test asserted > 8 and was simply wrong about the arithmetic.
    check("at pH 8.01 phosphate is mostly the basic form",
          close(base / acid, 10.0 ** (8.01 - 7.20), 1e-6),
          f"{acid:.2f} vs {base:.2f}, ratio {base / acid:.3f}")

    # Scaling: both members are linear in the total, which is the property
    # that makes the order in a species equal the order in the total at fixed
    # pH -- the degeneracy ANALYSIS.md section 6b rests on.
    one = dominant_buffer_pair("Phosphate", 6.71, 25.0)
    two = dominant_buffer_pair("Phosphate", 6.71, 85.0)
    ratio = 85.0 / 25.0
    check("both members scale with the total, so log-orders coincide",
          close(two[0] / one[0], ratio, 1e-9) and close(two[1] / one[1], ratio,
                                                        1e-9))

    # Unknown buffers must raise, for the reason speciation documents.
    try:
        dominant_buffer_pair("Tris", 8.0, 50.0)
        check("an unknown buffer raises", False)
    except KeyError:
        check("an unknown buffer raises", True)


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

    # Hand check of the DEFAULT model, Davies, at I = 75 mM:
    #   sqrt(0.075) = 0.2738613
    #   term  = 0.2738613 / 1.2738613 - 0.3 * 0.075 = 0.1924858
    #   shift = 0.509 * 2 * 0.1924858 = 0.19595
    expected = PKA_H2O2 - 0.19595
    check("Davies pKa at I = 75 mM matches the hand calculation",
          close(effective_pka_h2o2(75.0), expected, 1e-4),
          f"got {effective_pka_h2o2(75.0):.5f}, expected {expected:.5f}")

    # And of the extended Debye-Huckel model, still selectable:
    #   shift = 0.509 * 2 * 0.273861 / (1 + 0.328 * 0.273861) = 0.25581
    expected_debye = PKA_H2O2 - 0.25581
    check("Debye-Huckel pKa at I = 75 mM matches the hand calculation",
          close(effective_pka_h2o2(75.0, "debye"), expected_debye, 1e-4),
          f"got {effective_pka_h2o2(75.0, 'debye'):.5f}, expected {expected_debye:.5f}")

    check("the default model is Davies", ACTIVITY_MODEL == "davies")
    check("the two models agree at I = 0",
          close(effective_pka_h2o2(0.0, "debye"), effective_pka_h2o2(0.0, "davies"), 1e-12))
    check("Davies shifts the pKa less than Debye-Huckel wherever it matters",
          all(effective_pka_h2o2(I) > effective_pka_h2o2(I, "debye")
              for I in (50.0, 100.0, 400.0, 1069.0)))
    try:
        effective_pka_h2o2(75.0, "pitzer")
        check("an unknown activity model raises", False)
    except ValueError:
        check("an unknown activity model raises", True)

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

    # BUG 4: delta(z^2) is 2 for a neutral acid, not 1. Written against the
    # active model so switching models cannot quietly retire the check.
    root = math.sqrt(0.075)
    term = (root / (1 + root) - DAVIES_B * 0.075 if ACTIVITY_MODEL == "davies"
            else root / (1 + DEBYE_HUCKEL_B * root))
    single = DEBYE_HUCKEL_A * 1 * term
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


def test_oxygen_budget():
    """
    The budget behind two_axis/ section 5: can a side reaction make a bubble?

    Worked by hand, not against the module's own arithmetic. 2 H2O2 -> 2 H2O +
    O2 halves the peroxide; saturation is reached once twice the solubility has
    been consumed; and an ideal gas at 25 C is 24.45 L/mol, so 1 umol is 24.45
    uL.
    """
    print("\nthe oxygen budget")
    budget = oxygen_budget(73.424)
    check("the stoichiometry halves the peroxide",
          abs(budget["oxygen_mM"] - 36.712) < 1e-9,
          f"{budget['oxygen_mM']:.4f}")
    check("saturation costs twice the solubility, as a fraction of peroxide",
          abs(budget["saturation_fraction"]
              - 2 * OXYGEN_SOLUBILITY_mM / 73.424) < 1e-12)
    check("which is a few percent at the block's top peroxide",
          0.02 < budget["saturation_fraction"] < 0.05,
          f"{budget['saturation_fraction'] * 100:.1f}%")
    check("a mM of excess O2 is tens of microlitres per mL",
          abs(budget["microlitres_per_mM"] - MOLAR_VOLUME_uL_per_umol) < 1e-9)

    # More peroxide saturates the solution sooner, as a FRACTION of itself --
    # which is the direction that matters, because it says the top of the
    # block's ladder starts making gas earlier in its own turnover.
    check("more peroxide reaches saturation on a smaller fraction of itself",
          oxygen_budget(163.165)["saturation_fraction"]
          < oxygen_budget(3.671)["saturation_fraction"],
          f"{oxygen_budget(163.165)['saturation_fraction'] * 100:.2f}% against "
          f"{oxygen_budget(3.671)['saturation_fraction'] * 100:.2f}%")
    check("and the bottom of the ladder cannot saturate on a plausible "
          "conversion",
          oxygen_budget(3.671)["saturation_fraction"] > 0.5,
          f"{oxygen_budget(3.671)['saturation_fraction'] * 100:.0f}%")
    check("zero peroxide is handled rather than dividing",
          oxygen_budget(0.0)["oxygen_mM"] == 0.0)


if __name__ == "__main__":
    test_tables()
    test_speciation()
    test_dominant_buffer_pair()
    test_ionic_strength()
    test_debye_huckel()
    test_regressions()
    test_dataset()
    test_oxygen_budget()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
