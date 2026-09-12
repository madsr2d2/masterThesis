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
    test_michaelis_is_asked_of_rates_only()
    test_the_buffer_confound_is_reported()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
