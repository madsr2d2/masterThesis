"""
Verifies every number quoted in ph/ANALYSIS.md against the code.

    python ph/check_numbers.py
"""
import io
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import build_manifest
import induction
import ph_role
import scope
from doc_check import Checker

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")


def main():
    doc = Checker(DOCUMENT)

    print("\nsection 1: which curves, and the hand-sorted folders")
    doc.claim("the phosphate hand-sort", "2, 4, 5, 7, 8, 9, 10, 11, 12, 20, "
              "21, 22")
    doc.check("twelve hand-sorted phosphate exps against nine in scope.py",
              scope.PH_LADDER_PHOSPHATE ==
              (11, 9, 10, 22, 14, 20, 8, 21, 12))
    doc.claim("the missing one", "Exp 14 belongs and is missing")
    doc.check("exp 14 is in the published ladder",
              14 in scope.PH_LADDER_PHOSPHATE)
    doc.check("exps 2, 4, 5, 7 are REPLICATE_RUNS, not ladder rungs",
              scope.REPLICATE_RUNS == (2, 4, 5, 7)
              and not any(e in scope.PH_LADDER_PHOSPHATE
                          for e in scope.REPLICATE_RUNS))
    doc.check("the boric ladder is 41-49, not 13 and 41-49",
              scope.PH_LADDER_BORIC == (41, 42, 46, 47, 43, 48, 45, 49, 44)
              and 13 not in scope.PH_LADDER_BORIC)
    doc.check("exp 85's ruling covers the whole hand-sorted pH-11 set",
              "86-109" in build_manifest.KNOWN_EXCLUSIONS[85])
    doc.check("PH_LADDER_BORIC_BNOH is exps 60-62",
              scope.PH_LADDER_BORIC_BNOH == (60, 61, 62))
    doc.check("exp 50 is a validated reaction-direction exclusion",
              50 in build_manifest.KNOWN_EXCLUSIONS
              and "reaction-direction" in build_manifest.KNOWN_EXCLUSIONS[50])
    doc.check("BORIC_BNOH_HIGH_ENZYME_PAIR is exps 51 and 55, not 50",
              scope.BORIC_BNOH_HIGH_ENZYME_PAIR == (51, 55))
    doc.check("both are also SUBSTRATE_PAIRS' own BnOH half",
              (42, 51) in scope.SUBSTRATE_PAIRS
              and (45, 55) in scope.SUBSTRATE_PAIRS)

    print("\nsection 2: the rate's order in [HOO-]")
    table = ph_role.rate_ladder_table()
    doc.check("four ladders", len(table) == 4, f"{len(table)}")
    for label, ladder_row in (
            ("phosphate 4OMe", "phosphate 4OMe"),
            ("two-axis low (strong)", "pyrophosphate BnOH 136-142"),
            ("two-axis high (strong)", "pyrophosphate BnOH 143-151"),
            ("boric 4OMe", "boric 4OMe")):
        row = table.loc[ladder_row]
        bold = "**" if label != "boric 4OMe" else ""
        doc.claim(f"{label}: its row",
                  f"| {row.order_s0:+.3f} ± {row.stderr_s0:.3f} | "
                  f"{bold}{row.order_hoo:+.3f} ± {row.stderr_hoo:.3f}{bold} |")
    pooled_all = ph_role.pooled_rate_order()
    pooled_three = ph_role.pooled_rate_order(drop=("boric 4OMe",))
    # NOT bolded in the document: chi2 = 174 on 3 means the +/- is not an
    # uncertainty, so the value is quoted as a contrast and not argued for.
    doc.claim("pooled over all four", f"+{pooled_all['pooled']:.3f} ± "
              f"{pooled_all['stderr']:.3f}")
    doc.claim("its chi2", f"χ² = **{pooled_all['chi2']:.0f}** on "
              f"{pooled_all['dof']}")
    doc.claim("pooled without boric", f"**+{pooled_three['pooled']:.3f} ± "
              f"{pooled_three['stderr']:.3f}**")
    doc.claim("its chi2", f"χ² = **{pooled_three['chi2']:.2f}** on "
              f"{pooled_three['dof']}")

    # ---- and which of the three can actually corroborate the block ---------
    # Two of them ARE exps 136-151. The document says so, gives the weight they
    # carry, and makes its independent comparison the phosphate ladder alone.
    independent = ph_role.independent_check()
    doc.claim("the two-axis ladders' weight in the pool",
              f"**{100 * independent['two_axis_share']:.1f}% of the weight**")
    doc.claim("and phosphate's share",
              f"the other {100 * independent['weight_share']['phosphate 4OMe']:.1f}%")
    doc.check("the two pyrophosphate ladders really are the two-axis block",
              set(scope.PH_LADDER_TWO_AXIS_LOW) <= set(scope.TWO_AXIS_BLOCK)
              and set(scope.PH_LADDER_TWO_AXIS_HIGH) <= set(scope.TWO_AXIS_BLOCK))
    doc.check("and the phosphate ladder shares no experiment with it",
              not independent["shares_experiments_with_block"])
    doc.claim("phosphate alone, the independent row",
              f"**+{independent['phosphate_order']:.3f} ± "
              f"{independent['phosphate_stderr']:.3f}**")
    doc.claim("the two-axis block's own reading",
              f"**+{independent['block_order']:.3f} ± "
              f"{independent['block_stderr']:.3f}**")
    doc.claim("and the distance between them",
              f"**{independent['sigma']:.2f}σ**")
    doc.check("that distance is under one combined standard error",
              independent["sigma"] < 1.0, f"{independent['sigma']:.3f}")

    # ---- the [S] column is an order in the [S]/[buf] pair on the 4OMe rows --
    # `_ladder_scope`'s docstring claimed [buf] was fixed inside every ladder.
    # It is fixed BETWEEN runs (which is what protects the pH axis) and not
    # within them, so this is checked both ways rather than asserted.
    for name, low, high in (("phosphate 4OMe", 50.0, 80.0),
                            ("boric 4OMe", 70.0, 85.0)):
        block = scope.frame(scope.PH_LADDERS[name])
        block = block[block.live]
        by_run = block.groupby("experiment").buf.first()
        doc.check(f"{name}: [buf] is fixed BETWEEN runs",
                  by_run.nunique() == 1, f"{sorted(by_run.unique())}")
        doc.check(f"{name}: and steps within them, {low:.0f}-{high:.0f} mM",
                  abs(block.buf.min() - low) < 0.001
                  and abs(block.buf.max() - high) < 0.001
                  and block.buf.nunique() == 4,
                  f"{sorted(block.buf.unique())}")
        pair = np.corrcoef(np.log(block.s0), np.log(block.buf))[0, 1]
        doc.claim(f"{name}: the [S]/[buf] collinearity", f"−{abs(pair):.3f}")
    for name in ("pyrophosphate BnOH 136-142", "pyrophosphate BnOH 143-151"):
        block = scope.frame(scope.PH_LADDERS[name])
        doc.check(f"{name}: [buf] is fixed outright",
                  block[block.live].buf.nunique() == 1)
    # And that saying so costs the [HOO-] column nothing.
    for name, moved in (("phosphate 4OMe", 0.002), ("boric 4OMe", 0.001)):
        block = scope.frame(ph_role._ladder_scope(scope.PH_LADDERS[name]))
        block = block[block.live]
        base = scope.orders("vmax", frame=block, terms=("s0", "hoo"),
                            within=False)
        with_buf = scope.orders("vmax", frame=block,
                                terms=("s0", "hoo", "buf"), within=False)
        doc.check(f"{name}: adding [buf] moves order_hoo by {moved:.3f}",
                  abs(base["order_hoo"] - with_buf["order_hoo"]) < moved + 5e-4,
                  f"{abs(base['order_hoo'] - with_buf['order_hoo']):.4f}")
        doc.check(f"{name}: and returns an unresolved buffer order",
                  abs(with_buf["order_buf"]) < 2 * with_buf["stderr_buf"])

    print("\nsection 3: the boric ladder turns over")
    turnover = ph_role.boric_turnover()
    doc.check("the peak is exp 43, pH 9.50",
              turnover["peak_experiment"] == 43
              and abs(turnover["peak_pH"] - 9.50) < 0.005)
    frame = scope.frame(scope.PH_LADDER_BORIC)
    medians = frame.groupby("experiment").vmax.median()
    ph_by_exp = frame.groupby("experiment").pH.first()
    doc.claim("the rise", f"7.1 × 10⁻⁵ AU/s at pH 8.46")
    doc.check("exp 41's median vmax",
              abs(medians.loc[41] - 7.1e-5) < 0.05e-5)
    doc.claim("the peak value", f"1.42 × 10⁻⁴ at pH 9.50 (exp 43)")
    doc.check("exp 43's median vmax",
              abs(medians.loc[43] - 1.42e-4) < 0.01e-4)
    doc.claim("the fall", f"8.1 × 10⁻⁵ and 6.6 × 10⁻⁵ at pH 10.07 and 10.34")
    doc.check("exps 44 and 49's medians",
              abs(medians.loc[44] - 8.1e-5) < 0.05e-5
              and abs(medians.loc[49] - 6.6e-5) < 0.05e-5)
    # NOT bolded any more: the split that produces it was chosen off the same
    # medians it is fitted to, and the sweep below is what the document argues.
    doc.claim("the refit below the peak",
              f"+{turnover['below']['hoo']['slope']:.3f} ± "
              f"{turnover['below']['hoo']['stderr']:.3f}")

    # ---- and the split it depends on, swept ---------------------------------
    sweep = ph_role.boric_split_sensitivity()
    for split, row in sweep["table"].iterrows():
        published = " *(published)*" if abs(split - ph_role.BORIC_TURNOVER_SPLIT) < 1e-9 else ""
        doc.claim(f"the split at {split:.2f}",
                  f"| {split:.2f}{published} | {int(row.runs)} | "
                  f"+{row.order_hoo:.3f} ± {row.stderr_hoo:.3f} | "
                  f"{row.sigma:.1f} |")
    doc.claim("the range the estimate spans",
              f"**+{sweep['low']:.2f} to +{sweep['high']:.2f}**")
    doc.claim("how it should be quoted",
              f"**+{sweep['published']:.2f} with a split systematic of about "
              f"+{sweep['systematic_up']:.2f}/−{sweep['systematic_down']:.2f}**")
    doc.check("the sign survives every split",
              sweep["sign_stable"])
    doc.check("and every split stays under the other ladders' shared order",
              sweep["always_weaker_than_pooled"])
    doc.check("the span really is about a factor of nine",
              8.0 < sweep["high"] / sweep["low"] < 9.5,
              f"{sweep['high'] / sweep['low']:.2f}x")
    doc.claim("vmax_corrected barely moves it",
              "8.1 → 9.4 × 10⁻⁵ and\n6.6 → 7.0 × 10⁻⁵")
    above = turnover["above_by_experiment"]
    doc.check("44's corrected median is close to 9.4e-5",
              abs(above[44]["median_vmax_corrected"] - 9.4e-5) < 0.1e-5)
    doc.check("49's corrected median is close to 7.0e-5",
              abs(above[49]["median_vmax_corrected"] - 7.0e-5) < 0.1e-5)
    doc.check("every boric run's own internal consistency is 0.91-1.00",
              all(abs(np.corrcoef(
                  frame[(frame.experiment == e) & frame.live].s0,
                  frame[(frame.experiment == e) & frame.live].vmax)[0, 1])
                  >= 0.91 for e in scope.PH_LADDER_BORIC))

    print("\nsection 3a: the turnover is not only Vmax -- and Km turns over too")
    mm_table = ph_role.ladder_mm_table(scope.PH_LADDER_BORIC,
                                       response="v_peak_corrected")
    resolved = mm_table[mm_table.km_resolved].sort_values("pH")
    doc.check("7 of 9 boric rungs resolve their own km",
              len(resolved) == 7 and set(mm_table.experiment) - set(resolved.experiment)
              == {43, 45}, f"{len(resolved)} resolved, unresolved "
              f"{sorted(set(mm_table.experiment) - set(resolved.experiment))}")
    exp43 = mm_table.set_index("experiment").loc[43]
    doc.claim("exp 43's improved fit", "**+0.65**")
    doc.check("exp 43's r2 is now positive, off the debubbled rate",
              exp43.r2 > 0.6, f"R^2 = {exp43.r2:.4f}")
    doc.check("and its km is still unresolved either way",
              not exp43.km_resolved)

    for exp, mantissa, km in (("41", "1.28", "2.96"), ("42", "3.34", "5.56"),
                             ("46", "4.45", "8.13"), ("47", "3.16", "5.18"),
                             ("48", "4.74", "10.65"), ("44", "1.90", "4.56"),
                             ("49", "2.06", "6.66")):
        row = resolved.set_index("experiment").loc[int(exp)]
        doc.claim(f"exp {exp}'s row",
                  f"| {exp} | {row.pH:.2f} | {mantissa} × 10⁻⁴ | {km} |")
        doc.check(f"exp {exp}: vmax and km against the table",
                  abs(row.vmax - float(mantissa) * 1e-4) < 0.005e-4
                  and abs(row.km - float(km)) < 0.005,
                  f"{row.vmax:.3e}, {row.km:.3f}")
    doc.check("Km and Vmax both peak at exp 48, pH 9.51",
              resolved.set_index("experiment").km.idxmax() == 48
              and resolved.set_index("experiment").vmax.idxmax() == 48)

    decomposition = ph_role.boric_vmax_km_decomposition()
    drops = decomposition["drops"]
    decomp_table = decomposition["table"].set_index("experiment")
    doc.claim("the measured statistic's own drop", f"**−47.4%**")
    doc.check("raw drop, peak to last rung, peaking at exp 46",
              abs(drops["raw"] * 100 - 47.4) < 0.1
              and decomp_table.raw.idxmax() == 46,
              f"{drops['raw'] * 100:.1f}%, peak at exp {decomp_table.raw.idxmax()}")
    doc.claim("vmax alone", f"**−56.6%**, peaking at exp 48 (pH 9.51)")
    doc.check("vmax-only drop, peaking at exp 48",
              abs(drops["vmax_only"] * 100 - 56.6) < 0.1
              and decomp_table.vmax_only.idxmax() == 48,
              f"{drops['vmax_only'] * 100:.1f}%, peak at exp {decomp_table.vmax_only.idxmax()}")
    doc.claim("km alone",
              f"**−35.3%**, peaking at exp 41, the ladder's own lowest pH")
    doc.check("km-only drop",
              abs(drops["km_only"] * 100 - 35.3) < 0.1,
              f"{drops['km_only'] * 100:.1f}%")
    doc.check("km-only's own maximum sits at exp 41, the lowest pH",
              decomp_table.km_only.idxmax() == 41, f"peak at exp {decomp_table.km_only.idxmax()}")
    doc.check("vmax alone now falls further than the measured statistic",
              drops["vmax_only"] > drops["raw"],
              f"{drops['vmax_only'] * 100:.1f}% against {drops['raw'] * 100:.1f}%")

    shared_boric = ph_role.ladder_mm_shared(scope.PH_LADDER_BORIC,
                                            response="v_peak_corrected")
    shared_diag = ph_role.km_shared_diagnostic(mm_table, shared_boric)
    doc.check("the shared km fit now resolves",
              shared_boric["km_resolved"])
    doc.claim("the shared value", "**5.87 mM (3.24–12.31 mM)**")
    doc.check("shared km and its interval",
              abs(shared_boric["km"] - 5.87) < 0.01
              and abs(shared_boric["km_interval"][0] - 3.24) < 0.01
              and abs(shared_boric["km_interval"][1] - 12.31) < 0.01,
              f"{shared_boric['km']:.3f} "
              f"({shared_boric['km_interval'][0]:.3f}-"
              f"{shared_boric['km_interval'][1]:.3f})")
    doc.claim("how many resolved rungs the shared value sits inside now",
              f"Six of the seven")
    doc.check("6 of 7 resolved rungs contain the shared km, only exp 41 outside",
              shared_diag["agree"] == 6 and shared_diag["resolved"] == 7
              and not (resolved.set_index("experiment").loc[41].km_low
                       <= shared_boric["km"]
                       <= resolved.set_index("experiment").loc[41].km_high))

    print("\nsection 3b: the boric BnOH ladder (exps 60-62) is too sparse")
    doc.check("PH_LADDER_BORIC_BNOH is exps 60-62, not in PH_LADDERS",
              scope.PH_LADDER_BORIC_BNOH == (60, 61, 62)
              and "boric BnOH" not in scope.PH_LADDERS)
    bnoh_frame = scope.frame(scope.PH_LADDER_BORIC_BNOH)
    bnoh_live = bnoh_frame[bnoh_frame.live]
    doc.claim("one fixed enzyme loading", "(0.014 mM)")
    doc.check("[enz] is fixed at 0.014 mM across the ladder",
              bnoh_live.e0.nunique() == 1
              and abs(bnoh_live.e0.iloc[0] - 0.014) < 1e-6)
    doc.claim("no detachments at all",
              "zero\n`bubble_events` on all twelve live cuvettes")
    doc.check("twelve live cuvettes, none carrying a detachment",
              len(bnoh_live) == 12 and int(bnoh_live.bubble_events.sum()) == 0)
    doc.check("v_peak_corrected equals v_peak with no gas to correct",
              bool((bnoh_live.v_peak == bnoh_live.v_peak_corrected).all()))
    bnoh_table = ph_role.ladder_mm_table(scope.PH_LADDER_BORIC_BNOH,
                                         response="v_peak_corrected")
    doc.check("only one of three runs resolves its own km",
              int(bnoh_table.km_resolved.sum()) == 1
              and bool(bnoh_table.set_index("experiment").loc[60].km_resolved))
    doc.claim("exp 60's own km", "(4.48 mM)")
    doc.check("exp 60's km against the claim",
              abs(bnoh_table.set_index("experiment").loc[60].km - 4.48) < 0.01)
    exp62 = bnoh_frame[(bnoh_frame.experiment == 62) & bnoh_frame.live]
    exp62_by_s0 = exp62.sort_values("s0")
    doc.check("exp 62's top cuvette sits an order of magnitude above its "
              "neighbour",
              exp62_by_s0.v_peak_corrected.iloc[-1]
              / exp62_by_s0.v_peak_corrected.iloc[-2] > 5,
              f"{exp62_by_s0.v_peak_corrected.iloc[-1] / exp62_by_s0.v_peak_corrected.iloc[-2]:.1f}x")

    print("\nsection 3c: the high-enzyme boric BnOH pair (exps 51, 55)")
    doc.check("the pair is not in PH_LADDERS or PH_LADDER_BORIC_BNOH",
              "boric BnOH (0.28 mM)" not in scope.PH_LADDERS
              and not set(scope.BORIC_BNOH_HIGH_ENZYME_PAIR)
              & set(scope.PH_LADDER_BORIC_BNOH))
    pair_frame = scope.frame(scope.BORIC_BNOH_HIGH_ENZYME_PAIR)
    pair_live = pair_frame[pair_frame.live]
    doc.check("both runs sit at essentially PH_LADDER_BORIC's own loading",
              bool((abs(pair_live.e0 - 0.28) < 0.01).all())
              and abs(0.28 - 0.270) < 0.02)
    pair_table = ph_role.ladder_mm_table(scope.BORIC_BNOH_HIGH_ENZYME_PAIR,
                                         response="v_peak_corrected")
    by_exp = pair_table.set_index("experiment")
    doc.claim("exp 51's own vmax", "1.24 × 10⁻⁴ AU/s")
    doc.check("exp 51's vmax and unresolved km",
              abs(by_exp.loc[51].vmax - 1.24e-4) < 0.005e-4
              and not by_exp.loc[51].km_resolved)
    doc.claim("exp 55's own vmax and km",
              "0.98 × 10⁻⁴ AU/s, with `Km` resolved at 15.8 mM\n(7.3–72.1 mM)")
    doc.check("exp 55's vmax, km and interval",
              abs(by_exp.loc[55].vmax - 0.98e-4) < 0.005e-4
              and by_exp.loc[55].km_resolved
              and abs(by_exp.loc[55].km - 15.8) < 0.1
              and abs(by_exp.loc[55].km_low - 7.3) < 0.1
              and abs(by_exp.loc[55].km_high - 72.1) < 0.1)
    exp55 = pair_live[pair_live.experiment == 55].sort_values("s0")
    doc.claim("exp 55's worst cuvette", "**25** O2\ndetachments")
    doc.check("exp 55's lowest-S cuvette carries 25 events, others at most 1",
              int(exp55.iloc[0].bubble_events) == 25
              and (exp55.iloc[1:].bubble_events <= 1).all())
    doc.check("correcting it drops that cuvette's v_peak almost sixfold",
              exp55.iloc[0].v_peak / exp55.iloc[0].v_peak_corrected > 5.9,
              f"{exp55.iloc[0].v_peak:.3e} -> {exp55.iloc[0].v_peak_corrected:.3e}")
    doc.check("the pair's own rate falls from pH 9.01 to pH 9.70",
              by_exp.loc[51].vmax > by_exp.loc[55].vmax)

    print("\nsection 4: the induction clock's own pH order")
    clock_rows = induction.lag_ph_ladders()
    by_name = {row["ladder"]: row for row in clock_rows}
    expectations = {
        "phosphate 4OMe": (-0.253, 0.159, 0.093, 0.302, -0.25, 0.87),
        "boric 4OMe": (-0.020, 0.286, 0.250, 0.316, 0.71, -0.61),
        "pyrophosphate BnOH 136-142": (0.250, 0.286, 0.427, 0.367, -0.53, 0.77),
        "pyrophosphate BnOH 143-151": (0.370, 0.144, 0.413, 0.184, -0.79, 0.85),
    }
    for name, (slope, stderr, controlled, controlled_stderr,
              schedule, signal) in expectations.items():
        row = by_name[name]
        lag = row["lag_half_s"]
        doc.check(f"{name}: clock slope {slope:+.3f} +/- {stderr:.3f}",
                  abs(lag["slope"] - slope) < 0.002
                  and abs(lag["stderr"] - stderr) < 0.002,
                  f"{lag['slope']:+.3f} +/- {lag['stderr']:.3f}")
        doc.check(f"{name}: controlled {controlled:+.3f} +/- "
                  f"{controlled_stderr:.3f}",
                  abs(lag["controlled"] - controlled) < 0.002
                  and abs(lag["controlled_stderr"] - controlled_stderr)
                  < 0.002)
        doc.check(f"{name}: schedule {schedule:+.2f}, signal {signal:+.2f}",
                  abs(row["schedule_collinearity"] - schedule) < 0.005
                  and abs(row["signal_collinearity"] - signal) < 0.005,
                  f"{row['schedule_collinearity']:+.3f}, "
                  f"{row['signal_collinearity']:+.3f}")
    clock_pooled = ph_role.clock_pooled_order()
    doc.claim("pooled clock order",
              f"**+{clock_pooled['pooled']:.3f} ± "
              f"{clock_pooled['stderr']:.3f} per pH unit**")
    doc.claim("its chi2", f"χ² = **{clock_pooled['chi2']:.2f}** on "
              f"{clock_pooled['dof']}")

    print("\nsection 5: the two confounds")
    rate_by_name = table["schedule_collinearity"]
    rate_expected = {"phosphate 4OMe": 0.32, "boric 4OMe": 0.36,
                     "pyrophosphate BnOH 136-142": 0.84,
                     "pyrophosphate BnOH 143-151": -0.74}
    for name, expected in rate_expected.items():
        doc.check(f"{name}: rate schedule collinearity {expected:+.2f}",
                  abs(rate_by_name.loc[name] - expected) < 0.01,
                  f"{rate_by_name.loc[name]:+.3f}")

    print("\nsection 6: one species, three signals")
    consistency = scope.hoo_consistency()
    doc.claim("the within/across gap",
              f"part at **3.4σ**" if False else
              f"part at **{consistency['sigma']:.1f}σ**")
    doc.claim("the two orders",
              f"+{consistency['within_order']:.3f} ± "
              f"{consistency['within_stderr']:.3f} within, "
              f"+{consistency['across_order']:.3f} ± "
              f"{consistency['across_stderr']:.3f} across")
    import early_trough
    dominance = early_trough.dominance_correlation()
    rho_4ome = dominance["4OMe-BnOH"]["e0_hoo"][0]
    p_4ome = dominance["4OMe-BnOH"]["e0_hoo"][1]
    rho_bnoh = dominance["BnOH"]["e0_hoo"][0]
    p_bnoh = dominance["BnOH"]["e0_hoo"][1]
    doc.claim("the early trough's driver",
              f"ρ = −{abs(rho_4ome):.3f} for 4OMe, −{abs(rho_bnoh):.3f} "
              f"for BnOH")
    doc.check("both p < 1e-4", p_4ome < 1e-4 and p_bnoh < 1e-4,
              f"{p_4ome:.2e}, {p_bnoh:.2e}")
    gas = scope.gas_survey()
    boric_high = gas.loc[("Boric", pd.Interval(8.5, 14.0))]
    doc.claim("the boric gas onset",
              f"24 of 64 boric curves above pH 8.5 detach gas, 1.06\n"
              f"events/hour")
    doc.check("24 of 64, 1.06/hour",
              int(boric_high.detaching) == 24 and int(boric_high.curves) == 64
              and abs(boric_high.per_hour - 1.06) < 0.01,
              f"{boric_high.detaching} of {boric_high.curves}, "
              f"{boric_high.per_hour:.3f}/hour")
    doc.check("zero detachments anywhere below pH 7.5",
              all(g.events == 0 for (buf, band), g in gas.iterrows()
                  if band.right <= 7.5))

    print("\nthe figures the document promises")
    doc.figures(os.path.join(HERE, "index.html"))

    print("\nthe curves page draws every live cuvette of all six ladders")
    page = io.open(os.path.join(HERE, "progress_curves.html"),
                   encoding="utf-8").read()
    drawn = page.count("<div class='fig panel'>")
    live = 0
    # The four `scope.PH_LADDERS` plus `PH_LADDER_BORIC_BNOH` and
    # `BORIC_BNOH_HIGH_ENZYME_PAIR`, neither in that dict on purpose (see
    # their comments in scope.py) but both still drawn in full --
    # `build_figures.CURVE_PAGE_LADDERS` is the same union, kept here as
    # scope constants rather than imported to avoid a second definition of
    # that page-only dict.
    for exps in (list(scope.PH_LADDERS.values())
                 + [scope.PH_LADDER_BORIC_BNOH,
                    scope.BORIC_BNOH_HIGH_ENZYME_PAIR]):
        block = scope.frame(exps)
        live += int(block.live.sum())
    doc.check("one panel per live cuvette across the six ladders",
              drawn == live, f"{drawn} panels, {live} live cuvettes")

    print("\nthe figures: no data point drawn outside its own frame")
    doc.unclipped(os.path.join(HERE, "index.html"),
                  os.path.join(HERE, "progress_curves.html"))
    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
