"""
Tests for early_trough.py.

The load-bearing checks are that the three-part screen in `trough_table`
now correctly keeps exps 4.1 and 22.2 (once wrongly excluded by an overly
strict per-reading depth count -- see DATA_VERIFICATION.md), that the
archive-wide correlation the write-up quotes is reproduced exactly, that
HOO- beats total H2O2 as the archive-wide predictor, and that the rate
constant / Arrhenius / buffer-comparison figures quoted in ANALYSIS.md
match the code exactly.

    python data/test_early_trough.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import early_trough
from early_trough import (arrhenius_check, binding_rates, buffer_comparison,
                          cluster, dominance, dominance_correlation,
                          trough_table)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def test_dominance_takes_the_larger_ratio():
    """max, not min -- a curve dominated on either axis reads as dominated."""
    print("\ndominance is the larger of the two ratios")
    check("large [enz]/[HOO-], tiny [enz]/[S]",
          abs(dominance(1.0, 1000.0, 0.01) - 100.0) < 1e-9)
    check("large [enz]/[S], tiny [enz]/[HOO-]",
          abs(dominance(1.0, 0.01, 1000.0) - 100.0) < 1e-9)
    check("zero hoo does not crash, falls back to e0/s0",
          abs(dominance(1.0, 2.0, 0.0) - 0.5) < 1e-9)


def test_the_screen_now_keeps_the_two_corrected_curves():
    """
    Exps 4.1 and 22.2 clear the candidate threshold and were wrongly
    rejected by an earlier per-reading-count `sustained` test; the
    leave-one-out design (`curve_metrics.early_trough`) correctly keeps
    both, and exp 4.2 -- the strongest genuine example in the archive --
    still passes all three parts of the screen.
    """
    print("\nthe three-part screen: the two corrected curves and the "
          "strongest real example")
    table = trough_table()
    indexed = table.set_index(["experiment", "sample"])

    for exp, samp in ((4, 1), (22, 2)):
        row = indexed.loc[(exp, samp)]
        check(f"exp {exp}.{samp} is a candidate on the smoothed mean alone",
              row.candidate, f"z={row.z:.2f}")
        check(f"exp {exp}.{samp} is genuine (the corrected verdict)",
              row.genuine, f"sustained={row.sustained} genuine={row.genuine}")

    strongest = indexed.loc[(4, 2)]
    check("exp 4.2 passes all three parts of the screen", strongest.genuine,
          f"z={strongest.z:.2f} sustained={strongest.sustained} "
          f"overlap={strongest.bubble_overlap} survives={strongest.survives}")
    check("exp 4.2 is the strongest genuine curve in the archive",
          table[table.genuine].z.min() == strongest.z,
          f"{table[table.genuine].z.min():.2f} against {strongest.z:.2f}")
    check("there are currently no excluded candidates at all",
          int((table.candidate & ~table.genuine).sum()) == 0,
          f"{int((table.candidate & ~table.genuine).sum())} excluded")


def test_genuine_curves_never_overlap_a_bubble_event():
    """No curve in the genuine set owes its trough to a detected O2 event."""
    print("\nno genuine trough overlaps a detected O2 event")
    table = trough_table()
    genuine = table[table.genuine]
    check("at least one genuine curve exists", len(genuine) > 0,
          f"{len(genuine)}")
    check("none overlap a bubble event",
          not genuine.bubble_overlap.any(),
          f"{int(genuine.bubble_overlap.sum())} of {len(genuine)}")


def test_debubble_correction_does_not_explain_the_effect_away():
    """
    141.4 and 142.4 are the block's two heaviest bubblers among the genuine
    curves; if the O2 artefact were secretly responsible, correcting it
    should weaken their trough. It does the opposite.
    """
    print("\ndebubble correction does not explain the effect away")
    table = trough_table().set_index(["experiment", "sample"])
    for exp, samp in ((141, 4), (142, 4)):
        row = table.loc[(exp, samp)]
        check(f"exp {exp}.{samp}'s trough deepens after debubble correction",
              row.z_corrected < row.z,
              f"{row.z:.2f} -> {row.z_corrected:.2f}")


def test_the_dominance_correlation_by_substrate():
    """
    The write-up's headline numbers: [enz]/[HOO-] is the significant
    archive-wide predictor in BOTH substrates independently, [enz]/[S] is
    not, and -- the species test -- [enz]/[H2O2] (total, un-weighted
    peroxide) does not beat [enz]/[HOO-] either: it is not even correctly
    signed in 4OMe-BnOH.
    """
    print("\nthe archive-wide dominance correlation, by substrate")
    corr = dominance_correlation()
    check("both substrates were scanned", set(corr) == {"4OMe-BnOH", "BnOH"},
          f"{sorted(corr)}")
    for substrate in ("4OMe-BnOH", "BnOH"):
        rho_hoo, p_hoo, n_hoo = corr[substrate]["e0_hoo"]
        rho_s0, p_s0, n_s0 = corr[substrate]["e0_s0"]
        rho_h2o2, p_h2o2, n_h2o2 = corr[substrate]["e0_h2o2"]
        check(f"{substrate}: [enz]/[HOO-] correlates significantly "
              "(p < 1e-3)", p_hoo < 1e-3, f"rho={rho_hoo:+.3f} p={p_hoo:.2e}")
        check(f"{substrate}: the sign is negative (more dominance, deeper "
              "trough)", rho_hoo < 0, f"rho={rho_hoo:+.3f}")
        check(f"{substrate}: [enz]/[S] does not (p > 0.05)", p_s0 > 0.05,
              f"rho={rho_s0:+.3f} p={p_s0:.2e}")
        check(f"{substrate}: [enz]/[HOO-] beats [enz]/[H2O2]",
              abs(rho_hoo) > abs(rho_h2o2),
              f"|{rho_hoo:+.3f}| against |{rho_h2o2:+.3f}|")
    rho, p, n = corr["4OMe-BnOH"]["e0_hoo"]
    check("4OMe's own correlation is the exact figure quoted in the write-up",
          abs(rho - (-0.621)) < 0.01 and n == 147,
          f"rho={rho:+.4f} n={n}")
    rho, p, n = corr["BnOH"]["e0_hoo"]
    check("BnOH's own correlation is the exact figure quoted in the write-up",
          abs(rho - (-0.315)) < 0.01 and n == 164,
          f"rho={rho:+.4f} n={n}")
    rho_h2o2_4ome = corr["4OMe-BnOH"]["e0_h2o2"][0]
    check("[enz]/[H2O2] is not even correctly signed in 4OMe-BnOH",
          rho_h2o2_4ome > 0, f"rho={rho_h2o2_4ome:+.4f}")


def test_the_genuine_count_and_cluster_split():
    """Seventeen genuine curves, split into the two clusters by which
    ratio dominates -- this is the number every claim in ANALYSIS.md
    builds on, corrected from 15 once the sustained test's false rejections
    of exps 4.1 and 22.2 were fixed."""
    print("\nthe genuine count and the two clusters")
    table = trough_table()
    genuine = table[table.genuine]
    check("seventeen curves survive the full screen", len(genuine) == 17,
          f"{len(genuine)}")
    check("cluster is already a native column", "cluster" in table.columns)
    oxidant = int((genuine.cluster == "oxidant").sum())
    substrate = int((genuine.cluster == "substrate").sum())
    check("split fourteen oxidant-dominated, three substrate-dominated",
          oxidant == 14 and substrate == 3,
          f"oxidant={oxidant} substrate={substrate}")
    check("every substrate-cluster curve is BnOH at the block's lowest "
          "substrate rung",
          bool((genuine[genuine.cluster == "substrate"].s0 < 0.3).all()),
          f"{genuine[genuine.cluster == 'substrate'].s0.tolist()}")
    check("exps 4.1 and 22.2 both landed in the oxidant cluster (siblings "
          "of the already-confirmed 4.2 and 22.1)",
          bool((genuine.set_index(["experiment", "sample"])
               .loc[[(4, 1), (22, 2)]].cluster == "oxidant").all()))


def test_binding_rates_reproduce_the_quoted_range():
    """The k_on range and the two named extremes quoted in ANALYSIS.md."""
    print("\nbinding rates: the range quoted in the write-up")
    rates = binding_rates()
    check("seventeen rated curves", len(rates) == 17, f"{len(rates)}")
    check("k_on range is 0.51 to 32.68 M-1 s-1",
          abs(rates.k_on.min() - 0.51) < 0.02
          and abs(rates.k_on.max() - 32.68) < 0.02,
          f"{rates.k_on.min():.2f} to {rates.k_on.max():.2f}")
    check("the deepest trough (exp 4.2) is not the fastest rate",
          rates.sort_values("z").iloc[0][["experiment", "sample"]].tolist()
          != rates.sort_values("k_on", ascending=False).iloc[0][
              ["experiment", "sample"]].tolist())


def test_the_arrhenius_check_is_not_significant():
    """
    The honest result this session settled on: a proper 14-point
    regression gives an Ea whose own standard error is more than half its
    size, and the scatter among curves at one shared temperature is
    nearly as wide as the whole temperature range -- so the two-point
    estimate (85.7 kJ/mol, retracted) was never a measurement.
    """
    print("\nthe Arrhenius check: not resolved, and the code says so")
    fit = arrhenius_check()
    check("n=14 oxidant-cluster curves", fit["n"] == 14, f"{fit['n']}")
    check("Ea matches the figure quoted in the write-up",
          abs(fit["activation_kJ"] - 78.9) < 0.5, f"{fit['activation_kJ']:.1f}")
    check("its own standard error is more than half its value",
          fit["stderr_kJ"] > 0.5 * fit["activation_kJ"],
          f"{fit['activation_kJ']:.1f} +/- {fit['stderr_kJ']:.1f}")
    check("the slope is not significant at the usual 2-sigma bar",
          fit["t_statistic"] < 2.0, f"t={fit['t_statistic']:.2f}")
    check("scatter at one shared temperature is nearly as wide as the "
          "whole range (12 curves at 298.15 K)",
          fit["same_temperature_n"] == 12
          and fit["same_temperature_spread"] > 1.0,
          f"n={fit['same_temperature_n']} "
          f"spread={fit['same_temperature_spread']:.2f}")


def test_the_buffer_comparison():
    """
    Pyrophosphate's geometric-mean k_on runs well above phosphate's, even
    after normalising by [enz] identically -- the buffer-identity effect
    the write-up reads as evidence of a buffer-HOO- pre-equilibrium rather
    than a clean elementary E + HOO- step.
    """
    print("\nthe buffer comparison")
    table = buffer_comparison()
    check("both buffers present among the oxidant cluster",
          set(table.buffer) == {"Phosphate", "Pyrophosphate"},
          f"{sorted(table.buffer)}")
    phosphate = table[table.buffer == "Phosphate"].iloc[0]
    pyrophosphate = table[table.buffer == "Pyrophosphate"].iloc[0]
    check("phosphate: n=11, geometric mean matches the write-up",
          phosphate.n == 11 and abs(phosphate.geometric_mean - 2.85) < 0.05,
          f"n={phosphate.n} geo={phosphate.geometric_mean:.2f}")
    check("pyrophosphate: n=3, geometric mean matches the write-up",
          pyrophosphate.n == 3
          and abs(pyrophosphate.geometric_mean - 17.75) < 0.05,
          f"n={pyrophosphate.n} geo={pyrophosphate.geometric_mean:.2f}")
    check("pyrophosphate runs at least 5x phosphate's geometric mean",
          pyrophosphate.geometric_mean > 5 * phosphate.geometric_mean,
          f"{pyrophosphate.geometric_mean:.2f} against "
          f"{phosphate.geometric_mean:.2f}")


def test_trough_table_is_memoised_and_hands_out_a_copy():
    """
    `functools.cache` on a function returning a DataFrame is one in-place
    edit from a silent wrong answer elsewhere -- the same hazard
    `scope.frame` guards against. A caller mutating its own copy must not
    touch the next caller's.
    """
    print("\ntrough_table hands out a frame safe to mutate")
    first = trough_table()
    first.loc[first.index[0], "z"] = 999.0
    second = trough_table()
    check("mutating one caller's frame does not touch the next call's",
          second.z.iloc[0] != 999.0, f"{second.z.iloc[0]}")


if __name__ == "__main__":
    test_dominance_takes_the_larger_ratio()
    test_the_screen_now_keeps_the_two_corrected_curves()
    test_genuine_curves_never_overlap_a_bubble_event()
    test_debubble_correction_does_not_explain_the_effect_away()
    test_the_dominance_correlation_by_substrate()
    test_the_genuine_count_and_cluster_split()
    test_binding_rates_reproduce_the_quoted_range()
    test_the_arrhenius_check_is_not_significant()
    test_the_buffer_comparison()
    test_trough_table_is_memoised_and_hands_out_a_copy()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
