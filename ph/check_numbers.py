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
    doc.claim("pooled over all four", f"**+{pooled_all['pooled']:.3f} ± "
              f"{pooled_all['stderr']:.3f}**" if pooled_all["pooled"] > 0
              else f"**{pooled_all['pooled']:+.3f} ± {pooled_all['stderr']:.3f}**")
    doc.claim("its chi2", f"χ² = **{pooled_all['chi2']:.0f}** on "
              f"{pooled_all['dof']}")
    doc.claim("pooled without boric", f"**+{pooled_three['pooled']:.3f} ± "
              f"{pooled_three['stderr']:.3f}**")
    doc.claim("its chi2", f"χ² = **{pooled_three['chi2']:.2f}** on "
              f"{pooled_three['dof']}")
    published = scope.ph_order(parameter="vmax", scope=scope.strong_runs())
    doc.claim("the two-axis block's own reading",
              f"**+{published.loc['pooled', 'order']:.3f} ± "
              f"{published.loc['pooled', 'stderr']:.3f}**")
    doc.check("the two methods agree within one combined stderr",
              abs(pooled_three["pooled"] - published.loc["pooled", "order"])
              < (pooled_three["stderr"] + published.loc["pooled", "stderr"]))

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
    doc.claim("the refit below the peak",
              f"**+{turnover['below']['hoo']['slope']:.3f} ± "
              f"{turnover['below']['hoo']['stderr']:.3f}**")
    doc.check("six standard errors from zero",
              turnover["below"]["hoo"]["slope"]
              > 6 * turnover["below"]["hoo"]["stderr"])
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

    print("\nthe curves page draws every live cuvette of all four ladders")
    page = io.open(os.path.join(HERE, "progress_curves.html"),
                   encoding="utf-8").read()
    drawn = page.count("<div class='fig panel'>")
    live = 0
    for name, exps in scope.PH_LADDERS.items():
        block = scope.frame(exps)
        live += int(block.live.sum())
    doc.check("one panel per live cuvette across the four ladders",
              drawn == live, f"{drawn} panels, {live} live cuvettes")

    print("\nthe figures: no data point drawn outside its own frame")
    doc.unclipped(os.path.join(HERE, "index.html"),
                  os.path.join(HERE, "progress_curves.html"))
    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
