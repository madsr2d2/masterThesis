"""
The saturating fits, against planted schemes.

    python data/test_saturation.py

`saturation` asks each element of the fitted curve for a binding constant
rather than an order, so the checks are recovery: a planted K comes back, a
planted SHARED K is not rejected, and two elements planted with DIFFERENT K
are. The last is the one that matters -- a test that only ever fails to reject
would call every pair of elements one binding site.

Plus the two structural rules: a Michaelis-Menten fit is for a rate and never
for a clock, and it reads one peroxide level per run.

`run_gates.py` discovers it.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import saturation
import scope

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)
    return ok


def _planted_binding(rate_k, clock_k, noise=0.02, seed=4):
    """
    A peroxide ladder carrying both schemes, one run per level of everything
    else: v_act ~ K h/(1 + K h) and 1/tau = k_off(1 + K h), each run with its
    own arbitrary level, which is what the fits' per-run offsets absorb.
    """
    generator = np.random.default_rng(seed)
    peroxide = np.array([2.5, 7.3, 24.5, 73.4])
    rows = []
    for experiment in range(1, 6):
        scale = 10 ** generator.normal(0, 0.4)
        for sample, h in enumerate(peroxide, start=1):
            wobble = lambda: float(np.exp(generator.normal(0, noise)))
            rows.append({
                "experiment": experiment, "sample": sample, "live": True,
                "s0": 5.0, "h2o2": float(h), "buf": 50.0, "e0": 0.03,
                "pH": 8.0, "substrate": "BnOH", "buffer": "Pyrophosphate",
                "v_act_corrected": scale * (rate_k * h / (1 + rate_k * h))
                * wobble(),
                "v_act_resolved_corrected": True,
                "k_act_corrected": scale * 1e-3 * (1 + clock_k * h) * wobble(),
                "tau_act_resolved_corrected": True,
                "v0_act_corrected": np.nan,
                "v0_act_resolved_corrected": False,
                "k_sink_corrected": np.nan,
                "vmax_corrected": np.nan, "v_peak_corrected": np.nan})
    return pd.DataFrame(rows)


def _with_frame(frame, function, **kwargs):
    """Run one of `saturation`'s fits against a planted frame."""
    real = scope.frame
    scope.frame = lambda block=None: frame
    try:
        return function(**kwargs)
    finally:
        scope.frame = real


def test_a_planted_binding_constant_comes_back():
    print("\na planted K, per element")
    planted = _planted_binding(0.05, 0.05)
    table = _with_frame(planted, saturation.binding_by_element)
    for element, want in (("v_act_corrected", 0.05), ("k_act_corrected", 0.05)):
        row = table.loc[element]
        check(f"{element} recovers K", abs(row.K / want - 1) < 0.35,
              f"{row.K:.4f} against {want}")
        check(f"...and its interval covers it",
              row.K_low <= want <= row.K_high,
              f"[{row.K_low:.4f}, {row.K_high:.4f}]")


def test_one_binding_constant_is_told_from_two():
    print("\nshared K against separate K")
    same = _with_frame(_planted_binding(0.05, 0.05), saturation.shared_binding)
    check("a shared K is not rejected when the two schemes share one",
          same["f"] < 4.0, f"F {same['f']:.2f}, K {same['K_shared']:.4f}")
    check("...and the shared K is the planted one",
          same["K_shared_low"] <= 0.05 <= same["K_shared_high"],
          f"[{same['K_shared_low']:.4f}, {same['K_shared_high']:.4f}]")
    apart = _with_frame(_planted_binding(0.30, 0.002), saturation.shared_binding)
    check("two different K are rejected", apart["f"] > 12.0,
          f"F {apart['f']:.2f}, separate {apart['K_v_act_corrected']:.4f} "
          f"and {apart['K_k_act_corrected']:.4f}")


def _planted_series(activation_kJ, km_enthalpy_kJ, km_25=3.0, noise=0.01,
                    seed=7):
    """
    A temperature series whose Vmax and Km each carry a KNOWN barrier.

    Six runs, four substrate rungs each, rate = Vmax [S]/(Km + [S]) with
    Vmax Arrhenius in T and Km van 't Hoff in T. `km_enthalpy_kJ` = 0 plants
    a Km that does not move, which is what the archive's own series looks
    like; a nonzero one has to be told from it, or the decomposition is
    asserting rather than measuring.
    """
    generator = np.random.default_rng(seed)
    gas = 8.314462618 / 1000.0
    reference = 298.15
    rows = []
    for experiment, celsius in zip((14, 15, 16, 17, 18, 19),
                                   (25.0, 35.0, 40.0, 30.0, 20.0, 15.0)):
        kelvin = celsius + 273.15
        step = 1 / reference - 1 / kelvin
        vmax = 1e-4 * np.exp(activation_kJ / gas * step)
        km = km_25 * np.exp(km_enthalpy_kJ / gas * step)
        for sample, s0 in enumerate((1.85, 3.7, 5.55, 7.399), start=1):
            rate = vmax * s0 / (km + s0) * float(
                np.exp(generator.normal(0, noise)))
            rows.append({
                "experiment": experiment, "sample": sample, "live": True,
                "temperature": celsius, "s0": s0, "h2o2": 82.5,
                "buf": 80.0 - 5 * sample, "e0": 0.273, "pH": 7.0,
                "substrate": "4OMe-BnOH", "buffer": "Phosphate",
                "vmax_corrected": rate, "v_peak_corrected": rate,
                "v_act_corrected": rate, "v0_act_corrected": rate,
                "v_act_resolved_corrected": True,
                "v0_act_resolved_corrected": True,
                "k_act_corrected": np.nan, "k_sink_corrected": np.nan,
                "tau_act_resolved_corrected": False})
    return pd.DataFrame(rows)


def test_the_barrier_is_split_between_vmax_and_km():
    print("\nwhere the barrier sits: Vmax or Km")
    for planted_km in (0.0, 30.0):
        frame = _planted_series(80.0, planted_km)
        entry = _with_frame(frame, saturation.michaelis_temperature,
                            experiments=(14, 15, 16, 17, 18, 19))[0]
        check(f"Vmax's barrier comes back with Km planted at "
              f"{planted_km:.0f} kJ/mol",
              abs(entry["vmax_kJ"] - 80.0) < 8.0,
              f"{entry['vmax_kJ']:.1f} +/- {entry['vmax_stderr_kJ']:.1f}")
        check(f"...and Km's own {planted_km:.0f} kJ/mol is recovered",
              abs(entry["km_kJ"] - planted_km) < 12.0,
              f"{entry['km_kJ']:+.1f} +/- {entry['km_stderr_kJ']:.1f}")
    flat = _with_frame(_planted_series(80.0, 0.0),
                       saturation.michaelis_temperature,
                       experiments=(14, 15, 16, 17, 18, 19))[0]
    moving = _with_frame(_planted_series(80.0, 30.0),
                         saturation.michaelis_temperature,
                         experiments=(14, 15, 16, 17, 18, 19))[0]
    check("a moving Km is told from a flat one",
          moving["km_kJ"] - flat["km_kJ"] > 15.0,
          f"{flat['km_kJ']:+.1f} against {moving['km_kJ']:+.1f} kJ/mol")
    # The rung reconstruction: with Km fixed every rung predicts the same
    # barrier; with Km moving they must not.
    spread = flat["rungs"].predicted_kJ.max() - flat["rungs"].predicted_kJ.min()
    moved = moving["rungs"].predicted_kJ.max() - moving["rungs"].predicted_kJ.min()
    check("a fixed Km predicts one barrier at every rung", spread < 1.0,
          f"{spread:.2f} kJ/mol across the rungs")
    check("...and a moving Km predicts a different one at each",
          moved > spread, f"{moved:.2f} kJ/mol across the rungs")


def test_michaelis_is_asked_of_rates_only():
    print("\nthe substrate axis")
    table = saturation.michaelis_by_element(scope.PH_LADDER_BORIC)
    clocks = {"k_act_corrected", "k_sink_corrected"}
    check("no clock is given a Michaelis-Menten fit",
          not (set(table.element) & clocks), f"{sorted(set(table.element))}")
    check("every run's fit reads one peroxide level",
          bool((table.rungs >= saturation.MICHAELIS_MINIMUM_RUNGS).all()),
          f"{sorted(set(table.rungs))} rungs")
    boric = table[(table.element == "v_peak_corrected") & table.km_resolved]
    check("the boric ladder resolves Km on most of its runs, as ph/ has it",
          len(boric) >= 6, f"{len(boric)} of {table.experiment.nunique()} runs")


def test_the_buffer_confound_is_reported():
    print("\nthe [S]/[buf] pair")
    table = saturation.michaelis_by_element(scope.PH_LADDER_PHOSPHATE)
    check("a 4OMe ladder's Km is flagged as a Km of the pair",
          bool((table.buffer_r < -0.9).all()),
          f"median r {table.buffer_r.median():.3f}")
    block = saturation.michaelis_by_element(scope.TWO_AXIS_BLOCK)
    check("...and the two-axis block's is not, holding [buf] fixed",
          bool(block.buffer_r.isna().all()),
          f"{block.buffer_r.notna().sum()} runs with a buffer term")


if __name__ == "__main__":
    test_a_planted_binding_constant_comes_back()
    test_one_binding_constant_is_told_from_two()
    test_the_barrier_is_split_between_vmax_and_km()
    test_michaelis_is_asked_of_rates_only()
    test_the_buffer_confound_is_reported()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
