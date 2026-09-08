"""
Tests for early_trough.py.

The load-bearing checks are that the three-part screen in `trough_table`
actually removes the two known false positives (exps 4.1, 22.2) while
keeping the strongest real example (exp 4.2), and that the archive-wide
correlation the write-up quotes is reproduced exactly.

    python data/test_early_trough.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import early_trough
from early_trough import cluster, dominance, dominance_correlation, trough_table

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


def test_the_screen_removes_the_known_false_positives():
    """
    Exps 4.1 and 22.2 clear the candidate threshold on a smoothed mean alone
    and are the reason the `sustained` test exists; exp 4.2 is the
    strongest genuine example in the archive and must survive all three.
    """
    print("\nthe three-part screen: known false positives against the "
          "strongest real example")
    table = trough_table()
    indexed = table.set_index(["experiment", "sample"])

    for exp, samp in ((4, 1), (22, 2)):
        row = indexed.loc[(exp, samp)]
        check(f"exp {exp}.{samp} is a candidate on the smoothed mean alone",
              row.candidate, f"z={row.z:.2f}")
        check(f"exp {exp}.{samp} is rejected (not sustained)",
              not row.sustained and not row.genuine,
              f"sustained={row.sustained} genuine={row.genuine}")

    strongest = indexed.loc[(4, 2)]
    check("exp 4.2 passes all three parts of the screen", strongest.genuine,
          f"z={strongest.z:.2f} sustained={strongest.sustained} "
          f"overlap={strongest.bubble_overlap} survives={strongest.survives}")
    check("exp 4.2 is the strongest genuine curve in the archive",
          table[table.genuine].z.min() == strongest.z,
          f"{table[table.genuine].z.min():.2f} against {strongest.z:.2f}")


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
    archive-wide predictor in BOTH substrates independently, and
    [enz]/[S] is not, over every scanned live catalysed curve.
    """
    print("\nthe archive-wide dominance correlation, by substrate")
    corr = dominance_correlation()
    check("both substrates were scanned", set(corr) == {"4OMe-BnOH", "BnOH"},
          f"{sorted(corr)}")
    for substrate in ("4OMe-BnOH", "BnOH"):
        rho_hoo, p_hoo, n_hoo = corr[substrate]["e0_hoo"]
        rho_s0, p_s0, n_s0 = corr[substrate]["e0_s0"]
        check(f"{substrate}: [enz]/[HOO-] correlates significantly "
              "(p < 1e-3)", p_hoo < 1e-3, f"rho={rho_hoo:+.3f} p={p_hoo:.2e}")
        check(f"{substrate}: the sign is negative (more dominance, deeper "
              "trough)", rho_hoo < 0, f"rho={rho_hoo:+.3f}")
        check(f"{substrate}: [enz]/[S] does not (p > 0.05)", p_s0 > 0.05,
              f"rho={rho_s0:+.3f} p={p_s0:.2e}")
    rho, p, n = corr["4OMe-BnOH"]["e0_hoo"]
    check("4OMe's own correlation is the exact figure quoted in the write-up",
          abs(rho - (-0.621)) < 0.01 and n == 147,
          f"rho={rho:+.4f} n={n}")
    rho, p, n = corr["BnOH"]["e0_hoo"]
    check("BnOH's own correlation is the exact figure quoted in the write-up",
          abs(rho - (-0.315)) < 0.01 and n == 164,
          f"rho={rho:+.4f} n={n}")


def test_the_genuine_count_and_cluster_split():
    """Fifteen genuine curves, split into the two clusters by which ratio
    dominates -- this is the number every claim in ANALYSIS.md builds on."""
    print("\nthe genuine count and the two clusters")
    table = trough_table()
    genuine = table[table.genuine].copy()
    check("fifteen curves survive the full screen", len(genuine) == 15,
          f"{len(genuine)}")
    genuine["cluster"] = genuine.apply(cluster, axis=1)
    oxidant = int((genuine.cluster == "oxidant").sum())
    substrate = int((genuine.cluster == "substrate").sum())
    check("split twelve oxidant-dominated, three substrate-dominated",
          oxidant == 12 and substrate == 3,
          f"oxidant={oxidant} substrate={substrate}")
    check("every substrate-cluster curve is BnOH at the block's lowest "
          "substrate rung",
          bool((genuine[genuine.cluster == "substrate"].s0 < 0.3).all()),
          f"{genuine[genuine.cluster == 'substrate'].s0.tolist()}")


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
    test_the_screen_removes_the_known_false_positives()
    test_genuine_curves_never_overlap_a_bubble_event()
    test_debubble_correction_does_not_explain_the_effect_away()
    test_the_dominance_correlation_by_substrate()
    test_the_genuine_count_and_cluster_split()
    test_trough_table_is_memoised_and_hands_out_a_copy()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
