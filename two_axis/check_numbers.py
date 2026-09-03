"""
Verifies every number quoted in two_axis/ANALYSIS.md against the code.

The comparison itself is `doc_check`, shared with every sibling folder: a figure
typed into a document is a copy, and copies in this project have gone stale.

    python two_axis/check_numbers.py
"""
import io
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import induction
import scope
import slowdown
import solution_chemistry
from doc_check import Checker

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")


def main():
    doc = Checker(DOCUMENT)
    frame = scope.frame()
    summary = scope.summary()
    strong = scope.strong_runs()

    doc.section("section 1: the design")
    doc.claim("the block's size",
              f"**{summary['experiments']} runs of seven cuvettes**")
    doc.claim("its curves", f"{summary['curves']} curves, "
                            f"{summary['live_curves']} of them live")
    doc.claim("the pH span", f"{summary['pH_range'][0]:.2f} to "
                             f"{summary['pH_range'][1]:.2f}")
    doc.claim("the [HOO-] span",
              f"{summary['hoo_decades']:.1f} decades")
    doc.claim("the within-run share of [S]",
              f"**{summary['within_experiment_s0'] * 100:.1f}%**")
    doc.claim("and of [H2O2]",
              f"**{summary['within_experiment_h2o2'] * 100:.1f}%**")
    doc.claim("but of [HOO-]",
              f"**{scope.within_experiment_share('hoo') * 100:.1f}%**")
    doc.claim("the buffer is constant", f"{frame.buf.iloc[0]:g} mM")
    doc.check("and it really is, on every curve", frame.buf.nunique() == 1,
              f"{frame.buf.nunique()} values")

    design = scope.design()
    doc.claim("the run-length span",
              f"{design.duration_min.min() * 60:.0f} to "
              f"{design.duration_min.max() * 60:.0f} s")
    doc.claim("as a factor",
              f"{design.duration_min.max() / design.duration_min.min():.1f}x")
    for row in design.itertuples():
        doc.claim(f"exp {int(row.Index)}: its design row",
                  f"| {int(row.Index)} | {row.pH:.2f} | {row.hoo_mM:.3g} | "
                  f"{row.s0_ladder:.1f}x | {row.h2o2_ladder:.1f}x |")

    doc.section("section 1: the L has no interior")
    arms = induction.ladder_arms(frame)
    tops = frame.groupby("experiment").h2o2.transform("max")
    rights = frame.groupby("experiment").s0.transform("max")
    interior = frame[~np.isclose(frame.h2o2, tops)
                     & ~np.isclose(frame.s0, rights)]
    doc.check("no cuvette is off both arms", len(interior) == 0,
              f"{len(interior)} interior cuvettes")
    doc.claim("the arms' sizes",
              f"{len(arms['substrate arm'])} and "
              f"{len(arms['peroxide arm'])} live cuvettes")

    doc.section("section 2: the concentration orders")
    table = scope.order_table()
    for (parameter, fit), row in table.iterrows():
        bold = "**" if fit == "within-experiment" and parameter == "vmax" \
            else ""
        doc.claim(f"{parameter}, {fit}: the row",
                  f"| {bold}{row.order_s0:+.3f} +/- {row.stderr_s0:.3f}{bold} "
                  f"| {row.order_h2o2:+.3f} +/- {row.stderr_h2o2:.3f} | "
                  f"{int(row.n)} | {row.r2:.3f} |")
    live_curves = frame[frame.live]
    bottom = live_curves[live_curves.s0 == live_curves.s0.min()]
    doc.claim("the deepest conversion in the block",
              f"more than {live_curves.conversion.max() * 100:.1f}% of its "
              "substrate")
    doc.claim("and the bottom rung's median",
              f"the median conversion is {bottom.conversion.median() * 100:.1f}%")
    doc.check("so no curve is substrate-limited",
              live_curves.conversion.max() < 0.5,
              f"deepest {live_curves.conversion.max() * 100:.1f}%")

    doc.check("v_max has no substrate order at 2 sigma",
              abs(table.loc[("vmax", "within-experiment"), "order_s0"])
              < 2 * table.loc[("vmax", "within-experiment"), "stderr_s0"],
              f"{table.loc[('vmax', 'within-experiment'), 'order_s0']:+.3f}")
    doc.check("while v0 has one",
              abs(table.loc[("v0", "within-experiment"), "order_s0"])
              > 2 * table.loc[("v0", "within-experiment"), "stderr_s0"],
              f"{table.loc[('v0', 'within-experiment'), 'order_s0']:+.3f}")

    arm_table = scope.arm_orders()
    for (parameter, arm), row in arm_table.iterrows():
        doc.claim(f"{parameter}, {arm}: arm against joint",
                  f"| {row.arm_order:+.3f} +/- {row.arm_stderr:.3f} | "
                  f"{row.joint_order:+.3f} +/- {row.joint_stderr:.3f} | "
                  f"{row.sigma:.1f} |")
    doc.claim("the worst arm-to-joint gap",
              f"**{arm_table.sigma.max():.1f} sigma**")
    doc.check("and every arm agrees with the joint fit inside 2 sigma",
              arm_table.sigma.max() < 2.0, f"{arm_table.sigma.max():.2f}")

    induction_table = induction.induction_table(
        sorted(frame.experiment.unique()))
    saturation = induction.peroxide_saturation(induction_table)
    doc.claim("the free peroxide power",
              f"{saturation['order']:.3f} "
              f"({saturation['order_low']:.3f}-{saturation['order_high']:.3f})")
    doc.claim("first order rejected",
              f"**F = {saturation['first_order_f']:.0f}**")
    doc.claim("over what range",
              f"{saturation['peroxide_low']:g}-"
              f"{saturation['peroxide_high']:g} mM")

    doc.section("section 2: the agreement filter")
    agreement = scope.concentration_agreement()
    doc.claim("the floor", f"{scope.AGREEMENT_FLOOR:.2f}")
    doc.claim("how many runs pass", f"**{len(strong)} of "
                                    f"{summary['experiments']}**")
    weak = sorted(int(e) for e in agreement.index if e not in strong)
    doc.claim("and which fail", ", ".join(str(e) for e in weak))
    doc.check("exp 151 scores no agreement at all",
              151 not in agreement.index, f"{sorted(agreement.index)}")

    doc.section("section 3: the pH ladders")
    ladders = scope.ph_ladders(strong)
    doc.claim("how many ladders", f"**{len(ladders)} ladders**")
    for label, group in ladders.items():
        doc.claim(f"{label}: its runs",
                  f"exps {', '.join(str(int(e)) for e in sorted(group.experiment.unique()))}")
    orders = scope.ph_order("vmax", scope=strong)
    for label, row in orders.iterrows():
        bold = "**" if label == "pooled" else ""
        doc.claim(f"{label}: the pH order row",
                  f"| {int(row.runs)} | {int(row.curves)} | "
                  f"{row.pH_low:.2f}-{row.pH_high:.2f} | "
                  f"{bold}{row.order:+.3f} +/- {row.stderr:.3f}{bold} | "
                  f"{row.r2:.3f} |")
    doc.check("the pooled pH order is well below one",
              orders.loc["pooled", "order"]
              + 2 * orders.loc["pooled", "stderr"] < 1.0,
              f"{orders.loc['pooled', 'order']:+.3f}")
    doc.claim("the enzyme steps between runs",
              ", ".join(f"{value:g}" for value in sorted(frame.e0.unique())))

    unfiltered = scope.ph_order("vmax")
    doc.claim("what the weak runs do to it",
              f"{unfiltered.loc['pooled', 'order']:+.3f} +/- "
              f"{unfiltered.loc['pooled', 'stderr']:.3f}")
    doc.check("they flatten it rather than sharpen it",
              unfiltered.loc["pooled", "order"] < orders.loc["pooled", "order"],
              f"{unfiltered.loc['pooled', 'order']:+.3f} against "
              f"{orders.loc['pooled', 'order']:+.3f}")

    doc.section("section 3: the schedule control")
    schedule, verdict = scope.ph_schedule_control("vmax")
    dates = scope.run_dates()
    doc.claim("when the block was run",
              f"{dates.date.min():%-d} to "
              f"{dates.date.max():%-d September 2010}")
    doc.check("and experiment number is chronological",
              list(dates.sort_index().order) == sorted(dates.order),
              f"{list(dates.order)}")
    for label, row in schedule.iterrows():
        doc.claim(f"{label}: pH against the schedule",
                  f"{row.pH_vs_schedule:+.2f}")
    doc.check("the two ladders are opposed in the schedule",
              verdict["opposed_schedules"], f"{verdict}")
    doc.check("and agree in the sign of the pH order",
              verdict["orders_agree_in_sign"], f"{verdict}")

    doc.section("section 4: the two levers on [HOO-]")
    levers = scope.hoo_consistency("vmax")
    initial = scope.hoo_consistency("v0")
    for name, row in (("v_max", levers), ("v0", initial)):
        doc.claim(f"{name}: peroxide at fixed pH",
                  f"| {row['within_order']:+.3f} +/- "
                  f"{row['within_stderr']:.3f} | "
                  f"{row['across_order']:+.3f} +/- "
                  f"{row['across_stderr']:.3f} | "
                  f"**{row['sigma']:.1f}** |")
    doc.check("both statistics part company by more than 3 sigma",
              min(levers["sigma"], initial["sigma"]) > 3.0,
              f"{levers['sigma']:.2f}, {initial['sigma']:.2f}")
    doc.check("and peroxide is the stronger lever in both",
              levers["gap"] > 0 and initial["gap"] > 0,
              f"{levers['gap']:+.3f}, {initial['gap']:+.3f}")
    doc.claim("the levels the two contrasts sit at",
              f"{levers['within_hoo']:.4f} and "
              f"{levers['across_hoo']:.4f} mM")
    doc.check("and the higher-level contrast is the lower-order one",
              (levers["across_hoo"] > levers["within_hoo"])
              == (levers["across_order"] < levers["within_order"]),
              "which is the direction saturation predicts")

    doc.section("section 5: the gas")
    ladder = scope.bubble_ladder()
    for band, row in ladder.iterrows():
        low, high = band.left, band.right
        # The extremes are the argument -- none at the bottom, all at the top
        # -- so they carry the document's bold and the claim has to build it.
        extreme = int(row.with_drops) in (0, int(row.curves))
        mark = "**" if extreme else ""
        doc.claim(f"{low:g}-{high:g} mM: the peroxide band's row",
                  f"| {int(row.curves)} | {mark}{int(row.with_drops)}{mark} | "
                  f"{row.mean_lost:.4f}")
    doc.check("the drop count rises monotonically with peroxide",
              bool((np.diff(ladder.mean_drops.to_numpy()) > 0).all()),
              ", ".join(f"{v:.2f}" for v in ladder.mean_drops))

    runs = scope.bubble_turnover_control()
    top = runs[np.isclose(runs.top_h2o2, 73.424)]
    quiet = top[top.drops == 0]
    doc.claim("the quiet runs at the top peroxide",
              "Exps " + " and ".join(str(int(e)) for e in sorted(quiet.index)))
    doc.claim("their peroxide", f"{quiet.top_h2o2.iloc[0]:.1f} mM")
    doc.claim("their agreement",
              " and ".join(f"{v:.2f}" for v in quiet.agreement))
    doc.claim("against the noisy runs' agreement",
              f"{top[top.drops > 0].agreement.min():.2f} to "
              f"{top[top.drops > 0].agreement.max():.2f}")
    doc.check("the quiet runs really are the weakest at that peroxide",
              float(quiet.agreement.max())
              < float(top[top.drops > 0].agreement.min()))

    together = scope.bubble_synchrony()
    doc.claim("the cuvette pairs", f"{together['pairs']} cuvette pairs")
    doc.claim("coincident detachments",
              f"**{together['observed']} times against the "
              f"{together['expected']:.1f} independence predicts**")

    balance = scope.bubble_mass_balance()
    worst = balance.loc[balance.stitched.idxmax()]
    doc.claim("what stitching claims of the worst curve",
              f"**{worst.stitched:.2f}x the most absorbance its "
              f"{worst.s0:.3f} mM of substrate could ever make**")
    doc.claim("and of the next four",
              f"between {balance.nlargest(5, 'stitched').stitched.min():.2f} "
              f"and {balance.nlargest(5, 'stitched').stitched.iloc[1]:.2f}")
    doc.claim("their ramps steepen",
              f"{balance.nlargest(5, 'stitched').ramp_gain.min():.1f} to "
              f"{balance.nlargest(5, 'stitched').ramp_gain.max():.1f} times")
    doc.check("no curve as read exceeds its own substrate",
              int((balance.raw > 1.0).sum()) == 0)
    doc.check("and the subtraction leaves none exceeding it",
              int((balance.corrected > 1.0).sum()) == 0)

    asymmetry = scope.bubble_step_asymmetry()
    doc.claim("how many steps the block has",
              f"**{asymmetry['steps']}** steps")
    doc.claim("how many rise past the tail threshold",
              f"**{asymmetry['rises']} rises** beyond "
              f"{scope.BUBBLE_STEP_SIGMA:g}sigma")
    doc.claim("and how many fall", f"**{asymmetry['falls']} falls**")
    doc.claim("their ratio", f"a ratio of\n**{asymmetry['ratio']:.1f}**")
    doc.claim("the largest fall against the largest rise",
              f"({asymmetry['largest_fall']:.0f}sigma) "
              f"{abs(asymmetry['largest_fall']) / asymmetry['largest_rise']:.1f}x the\n"
              f"largest rise (+{asymmetry['largest_rise']:.0f}sigma)")
    doc.check("large steps fall far more often than they rise",
              asymmetry["falls"] > 3 * asymmetry["rises"],
              f"{asymmetry['falls']} against {asymmetry['rises']}")

    doc.claim("the solubility the budget rests on",
              f"**{solution_chemistry.OXYGEN_SOLUBILITY_mM:.2f} mM**")
    worked_h2o2 = float(frame[
        (frame.experiment == scope.BUBBLE_WORKED_EXAMPLE[0])
        & (frame["sample"] == scope.BUBBLE_WORKED_EXAMPLE[1])].h2o2.iloc[0])
    ladder_bottom = float(frame.h2o2.min())
    for label, level in (("the block's top peroxide", float(frame.h2o2.max())),
                         ("where the turnover control sits", worked_h2o2),
                         ("the bottom of the ladder", ladder_bottom)):
        budget = solution_chemistry.oxygen_budget(level)
        # The bottom of the ladder is quoted to whole percent because it is
        # past 100 and the tenth would be noise on an argument about whether a
        # bubble is possible at all.
        precision = 0 if label == "the bottom of the ladder" else 1
        doc.claim(f"the budget row at {level:.1f} mM",
                  f"| {label}, {level:.1f} mM | "
                  f"**{budget['oxygen_mM']:.1f} mM** | "
                  f"**{budget['saturation_fraction'] * 100:.{precision}f}%** |")
    doc.claim("what each further mM is",
              f"**{solution_chemistry.MOLAR_VOLUME_uL_per_umol:.2f} uL of gas "
              f"per mL of\nsolution**")
    top_budget = solution_chemistry.oxygen_budget(float(frame.h2o2.max()))
    doc.check("the top of the ladder saturates on a couple of percent",
              top_budget["saturation_fraction"] < 0.05,
              f"{top_budget['saturation_fraction'] * 100:.1f}%")
    bottom_budget = solution_chemistry.oxygen_budget(ladder_bottom)
    doc.check("and the bottom of it cannot saturate at all here",
              bottom_budget["saturation_fraction"] > 0.5,
              f"{bottom_budget['saturation_fraction'] * 100:.0f}%")
    doc.check("no live curve below 5 mM carries a detachment",
              int(frame[frame.live & (frame.h2o2 < 5.0)].bubble_drops.sum())
              == 0)

    record = scope.bubble_record(*scope.BUBBLE_WORKED_EXAMPLE)
    worked = frame[(frame.experiment == scope.BUBBLE_WORKED_EXAMPLE[0])
                   & (frame["sample"] == scope.BUBBLE_WORKED_EXAMPLE[1])].iloc[0]
    doc.claim("the worked curve's conditions",
              f"{worked.s0:g} mM substrate, {worked.h2o2:.1f} mM H2O2, "
              f"pH {worked.pH:.2f}")
    for row in record.itertuples():
        bold = "**" if row.lost == record["lost"].max() else ""
        window = "**" if row.grew_s == record.grew_s.min() else ""
        after = "**" if row.after == record.after.min() else ""
        doc.claim(f"its detachment at {row.time_s:.0f} s",
                  f"| {row.time_s:.0f} s | {bold}{row.lost:.4f}{bold} | "
                  f"{row.sigma:.1f} | {window}{row.grew_s:.0f} s{window} | "
                  f"+{row.rose:.4f} | {after}+{row.after:.4f}{after} |")
    doc.claim("what it sheds in total",
              f"**{record['lost'].sum():.4f} AU** in total against a net rise of "
              f"{worked.net:.4f}")
    doc.claim("and its load", f"load {worked.bubble_load:.2f}")
    doc.check("the largest drop follows the shortest growth window",
              record.loc[record["lost"].idxmax(), "grew_s"] == record.grew_s.min(),
              f"{record.loc[record['lost'].idxmax(), 'grew_s']:.0f} s against "
              f"{record.grew_s.min():.0f} s")
    doc.check("and the last level sits below every earlier one",
              record.after.iloc[-1] < record.after.iloc[:-1].min(),
              f"{record.after.iloc[-1]:.4f} against "
              f"{record.after.iloc[:-1].min():.4f}")
    doc.claim("the levels it climbs through",
              ", ".join(f"+{v:.4f}" for v in record.after[:-1]))
    doc.claim("and the one it falls back to",
              f"then falls to +{record.after.iloc[-1]:.4f}, below all three")

    table = scope.bubble_table()
    bands = pd.cut(table.bubble_load, [-0.001, 1e-4, 0.25, 0.5, 1.0, np.inf],
                   labels=["none", "<=0.25", "0.25-0.5", "0.5-1", ">1"])
    counts = table.groupby(bands, observed=True)
    for name, group in counts:
        # Only the band that carries no rate is argued for, so only it is bold.
        mark = "**" if str(name) == ">1" else ""
        doc.claim(f"load {name}: the row",
                  f"| {mark}{len(group)}{mark} | "
                  f"{mark}{group.h2o2.median():.1f}{mark} |")
    beyond = table[~table.repairable]
    doc.claim("how many curves carry no measurable rate", f"**{len(beyond)}**")
    doc.claim("the run that carries most of them",
              f"**all four substrate rungs of exp "
              f"{int(beyond.experiment.mode().iloc[0])}**")
    doc.claim("its peroxide", f"{beyond.h2o2.max():.0f} mM peroxide")
    doc.claim("its loads",
              f"loads {beyond[beyond.experiment == 135].bubble_load.min():.1f} "
              f"to {beyond[beyond.experiment == 135].bubble_load.max():.1f}")
    doc.claim("the other runs that carry one",
              ", ".join(str(int(e)) for e in sorted(
                  set(beyond.experiment) - {135})[:-1])
              + " and " + str(int(sorted(set(beyond.experiment) - {135})[-1])))
    doc.check("the ceiling excludes nothing -- every live curve is in the table",
              len(table) == int(frame.live.sum()),
              f"{len(table)} against {int(frame.live.sum())}")
    doc.claim("what the bound costs a typical curve",
              f"**{(1 - (table.vmax_monotone / table.vmax).median()) * 100:.1f}%**")
    doc.claim("and the worst one",
              f"**{(1 - (table.vmax_monotone / table.vmax).min()) * 100:.0f}%**")

    for emptying in (True, False):
        how = "emptying" if emptying else "partly emptying"
        recovery = scope.bubble_recovery(emptying=emptying)
        for severity, row in recovery.iterrows():
            # Bold marks what the section argues: that stitching is worse than
            # doing nothing, and that the reconstruction is close to the truth.
            doc.claim(f"{how}, the recovery row at {severity:g}x",
                      f"| {severity:g} | {row.raw:.2f} | "
                      f"**{row.stitched:.2f}** | **{row.rebuilt:.2f}** |")
        doc.check(f"{how}: stitching never beats leaving the curve alone",
                  bool((recovery.stitched >= recovery.raw - 1e-9).all()))
        doc.check(f"{how}: the reconstruction is within a tenth throughout",
                  bool(((recovery.rebuilt - 1.0).abs() < 0.10).all()),
                  f"worst {(recovery.rebuilt - 1.0).abs().max():.3f}")

    smoothness = scope.rebuild_smoothness()
    repaired = smoothness[~smoothness.clean]
    clean = smoothness[smoothness.clean]
    doc.claim("how many curves carry gas",
              f"the {len(repaired)} detaching curves")
    doc.claim("their worst fall as read",
              f"| the {len(repaired)} detaching curves, as read | "
              f"**{repaired.raw_worst.min():.1f}sigma** |")
    doc.claim("their worst fall rebuilt",
              f"| the same {len(repaired)}, rebuilt | "
              f"**{repaired.rebuilt_worst.min():.1f}sigma** |")
    doc.claim("and the clean curves' worst fall",
              f"| the {len(clean)} that never bubbled | "
              f"{clean.rebuilt_worst.min():.1f}sigma |")
    survivors = repaired[repaired.rebuilt_worst < -scope.REBUILD_STEP_SIGMA]
    doc.claim("how many repaired curves still carry one",
              f"**Exactly one repaired curve still carries a fall past "
              f"{scope.REBUILD_STEP_SIGMA:g}sigma**")
    doc.check("and it is the first-interval curve",
              [(int(r.experiment), int(r["sample"]))
               for _, r in survivors.iterrows()] == [(135, 6)],
              f"{[(int(r.experiment), int(r['sample'])) for _, r in survivors.iterrows()]}")
    for experiment, sample in ((141, 3), (144, 2)):
        row = smoothness[(smoothness.experiment == experiment)
                         & (smoothness["sample"] == sample)].iloc[0]
        doc.claim(f"exp {experiment} cuvette {sample}'s worst step",
                  f"{row.raw_worst:.1f}sigma to {row.rebuilt_worst:.1f}sigma")

    drivers = scope.gas_rate_drivers()
    doc.claim("what the gas rate does with peroxide",
              f"**first order in peroxide, {drivers['pooled_h2o2']:+.3f} +/- "
              f"{drivers['pooled_stderr_h2o2']:.3f}**")
    doc.claim("on how many curves", f"over the {drivers['n_pooled']} live curves")
    doc.claim("how far that is from zero",
              f"{drivers['pooled_h2o2'] / drivers['pooled_stderr_h2o2']:.1f}"
              "sigma from zero")
    doc.claim("and what it does with substrate",
              f"{drivers['order_s0']:+.3f} +/- {drivers['stderr_s0']:.3f}")
    doc.check("the gas rate is first order in peroxide",
              abs(drivers["pooled_h2o2"] - 1.0)
              < 2 * drivers["pooled_stderr_h2o2"])
    doc.check("and the substrate carries far less of it",
              abs(drivers["order_s0"]) < 0.5 * abs(drivers["pooled_h2o2"]))

    sensitivity = scope.bubble_sensitivity()
    published = sensitivity.loc[("vmax", "all live")]
    for (treatment, subset), row in sensitivity.iterrows():
        # The peroxide order under the reconstruction is the one number in the
        # table the section argues for, so it is the one carrying the bold.
        mark = "**" if treatment == "vmax_corrected" else ""
        doc.claim(f"{treatment} on {subset}: the order row",
                  f"| {int(row.n)} | {row.order_s0:+.3f} +/- "
                  f"{row.stderr_s0:.3f} | {mark}{row.order_h2o2:+.3f} +/- "
                  f"{row.stderr_h2o2:.3f}{mark} |")
    for (treatment, subset), row in sensitivity.iterrows():
        if (treatment, subset) == ("vmax", "all live"):
            continue
        gap = abs(row["order_s0"] - published["order_s0"])
        spread = float(np.hypot(row["stderr_s0"], published["stderr_s0"]))
        doc.check(f"{treatment}/{subset}: the substrate order does not move",
                  gap < spread, f"{gap:.3f} against {spread:.3f}")
    corrected = sensitivity.loc[("vmax_corrected", "all live")]
    peroxide_sigma = (abs(corrected["order_h2o2"] - published["order_h2o2"])
                      / float(np.hypot(corrected["stderr_h2o2"],
                                       published["stderr_h2o2"])))
    doc.check("the peroxide order comes down under the reconstruction",
              corrected["order_h2o2"] < published["order_h2o2"],
              f"{corrected['order_h2o2']:+.3f} against "
              f"{published['order_h2o2']:+.3f}")
    doc.check("and not by enough to matter", peroxide_sigma < 1.5,
              f"{peroxide_sigma:.1f} sigma")
    doc.check("the substrate order moves away from zero, not toward it",
              abs(sensitivity.loc[("vmax_monotone", "all live"), "order_s0"])
              >= abs(published["order_s0"]))
    doc.check("every live curve still carries a rate after the repair",
              int((frame[frame.live].vmax_corrected <= 0).sum()) == 0)
    strong_published = scope.orders("vmax", scope=scope.strong_runs())
    strong_rebuilt = scope.orders("vmax_corrected", scope=scope.strong_runs())
    doc.claim("the peroxide order on the strong runs",
              f"from {strong_published['order_h2o2']:+.3f} to "
              f"{strong_rebuilt['order_h2o2']:+.3f}")
    strong_sigma = (abs(strong_rebuilt["order_h2o2"]
                        - strong_published["order_h2o2"])
                    / float(np.hypot(strong_rebuilt["stderr_h2o2"],
                                     strong_published["stderr_h2o2"])))
    doc.claim("how far both shifts are",
              f"which is {peroxide_sigma:.1f}sigma and {strong_sigma:.1f}sigma")
    doc.check("the strong runs move the same way", strong_sigma < 1.5,
              f"{strong_sigma:.1f} sigma")

    doc.claim("what the worked tail curve sheds",
              f"**{smoothness[(smoothness.experiment == 149) & (smoothness['sample'] == 4)].biggest_bubble.iloc[0]:.4f} AU once**")
    ratio = repaired.gas_held / repaired.biggest_bubble
    doc.claim("the worst curve's held gas after the ceiling",
              f"no curve holds more than {ratio.max():.1f}x its own largest "
              f"bubble")
    doc.claim("and the typical one", f"the median is\n{ratio.median():.1f}x")
    negative = smoothness[smoothness.rebuilt_net < 0]
    doc.claim("the one curve that rebuilds below zero",
              f"exp {int(negative.experiment.iloc[0])} cuvette "
              f"{int(negative['sample'].iloc[0])}")
    doc.claim("what it sheds and lands at",
              f"sheds {negative.bubble_load.iloc[0]:.1f}x\n"
              f"its own net rise and lands at {negative.rebuilt_net.iloc[0]:.4f} AU")
    doc.check("only one curve rebuilds below zero", len(negative) == 1)
    doc.check("and its load already flags it",
              float(negative.bubble_load.iloc[0]) > scope.BUBBLE_LOAD_CEILING)
    doc.check("no curve holds many times its own largest bubble",
              float(ratio.max()) < 6.0, f"{ratio.max():.1f}x")

    doc.section("section 5: what the curves do")
    bands = scope.acceleration_by_ph()
    for band, row in bands.iterrows():
        doc.claim(f"{band}: the acceleration row",
                  f"| {int(row.curves)} | {int(row.accelerating)} | "
                  f"{row.share * 100:.0f}% | "
                  f"{row.median_late_over_early:+.2f} |")
    doc.check("acceleration is commoner above pH 9 than below",
              bands.loc["pH >= 9", "share"] > bands.loc["pH < 9", "share"],
              f"{bands.loc['pH >= 9', 'share']:.2f} against "
              f"{bands.loc['pH < 9', 'share']:.2f}")

    control = induction.signal_control(induction_table)
    doc.claim("the signal control",
              f"**{control['signal_slope']:+.3f} +/- "
              f"{control['signal_stderr']:.3f}**")
    doc.check("and it fails, at more than 2 sigma from zero",
              abs(control["signal_slope"]) > 2 * control["signal_stderr"],
              f"{control['signal_slope'] / control['signal_stderr']:.1f} sigma")

    live = frame[frame.live]
    signs = induction.sign_table(frame)
    signs = signs[signs.live]
    doc.claim("the phase split",
              f"{int((live.phases == 2).sum())} two-phase against "
              f"{int((live.phases == 1).sum())} one-phase")
    doc.claim("the sign split",
              f"{int(signs.lag_first.sum())} lag-first against "
              f"{int((~signs.lag_first).sum())} burst-first")

    doc.section("section 5: the burst amplitude")
    bursts = scope.burst_table()
    bounded = bursts[bursts.burst_bounded]
    doc.claim("how many curves finish their rise",
              f"**Only {len(bounded)} of the {len(bursts)} live curves finish "
              "their rise inside their run.**")
    inside = frame[frame.live].groupby("experiment").duration_s.agg(
        lambda d: d.max() - d.min())
    doc.claim("the within-run spread in length", "one 60 s reading")
    doc.check("and no run varies inside itself by more than that",
              inside.max() <= 60.0, f"{inside.max():.0f} s")
    for band, rows in bounded.groupby(np.where(bounded.pH >= 9.0,
                                               "pH >= 9", "pH < 9")):
        doc.claim(f"{band}: the turnover row",
                  f"| {band} | {len(rows)} | {rows.turnovers.median():.2f} | "
                  f"{rows.turnovers.min():.2f}-{rows.turnovers.max():.2f} |")
    doc.check("the catalyst turns over more often at high pH",
              bounded[bounded.pH >= 9].turnovers.median()
              > bounded[bounded.pH < 9].turnovers.median(),
              f"{bounded[bounded.pH >= 9].turnovers.median():.2f} against "
              f"{bounded[bounded.pH < 9].turnovers.median():.2f}")
    one_fifty = bursts[(bursts.experiment == 150) & (bursts["sample"] == 6)]
    doc.claim("exp 150 cuvette 6, the curve that prompted this",
              f"is {one_fifty.turnovers.iloc[0]:.2f}")
    doc.claim("its fitted rise",
              f"**{one_fifty.burst.iloc[0]:.4f} AU**")
    doc.claim("and where it peaks",
              f"peaks at\n{one_fifty.burst_time_s.iloc[0]:.0f} s")
    in_strong = bounded[bounded.experiment.isin(strong)]
    weakest = in_strong.loc[in_strong.turnovers.idxmin()]
    doc.claim("the weakest burst in a strong run",
              f"exp {int(weakest.experiment)} cuvette "
              f"{int(weakest['sample'])} at {weakest.turnovers:.2f}")
    doc.check("and it is the only sub-stoichiometric one there",
              int((in_strong.turnovers < 1.0).sum()) == 1,
              f"{int((in_strong.turnovers < 1.0).sum())} below one turnover")

    for label, subset in (("all live curves", scope.TWO_AXIS_BLOCK),
                          ("strong runs", strong)):
        row = scope.burst_drivers(subset)
        doc.claim(f"{label}: the burst's concentration orders",
                  f"| {label} | {row['order_s0']:+.3f} +/- "
                  f"{row['stderr_s0']:.3f} | {row['order_h2o2']:+.3f} +/- "
                  f"{row['stderr_h2o2']:.3f} | {int(row['n'])} | "
                  f"{row['r2']:.3f} |")
        doc.check(f"{label}: and the catalyst is not identified there",
                  not row["enzyme_identified"], "[enz] is constant per run")

    doc.section("section 5: the enzyme lever")
    pair, verdict = scope.enzyme_pair()
    doc.claim("which runs", f"**exps {verdict['high']} and {verdict['low']}**")
    doc.claim("the loading step",
              f"{frame[frame.experiment == verdict['high']].e0.iloc[0]:g}\n"
              f"against {frame[frame.experiment == verdict['low']].e0.iloc[0]:g} mM")
    doc.claim("the pH gap", f"{verdict['pH_gap']:.2f} pH units apart")
    doc.claim("the window", f"**{verdict['window_s']:.0f} s**")
    doc.claim("the pH correction",
              f"{verdict['pH_correction']:.3f}x against the "
              f"{verdict['expected']:.2f}x enzyme step")
    doc.claim("the result row",
              f"| **{verdict['ratio']:.2f} x/÷ {verdict['stderr_factor']:.2f}** "
              f"| {verdict['expected']:.2f} | "
              f"{verdict['sigma_first_order']:.1f} | "
              f"**{verdict['sigma_no_dependence']:.1f}** |")
    doc.check("no dependence on catalyst is excluded",
              verdict["sigma_no_dependence"] > 3.0,
              f"{verdict['sigma_no_dependence']:.2f} sigma")
    doc.check("and first order in it is not",
              verdict["sigma_first_order"] < 2.0,
              f"{verdict['sigma_first_order']:.2f} sigma")
    doc.check("all seven cuvettes pair up", verdict["cuvettes"] == 7,
              f"{verdict['cuvettes']}")

    sweep = scope.enzyme_pair_sensitivity()
    doc.claim("the window systematic",
              f"+{sweep.order.min():.2f} to +{sweep.order.max():.2f} over "
              f"{sweep.window_s.min():.0f} to {sweep.window_s.max():.0f} s")
    for _, row in sweep.iterrows():
        doc.claim(f"the sweep at {row.window_s:.0f} s",
                  f"| {row.window_s:.0f} | {row.ratio:.2f} | "
                  f"{row.order:+.2f} | {int(row.cuvettes)} |")
    doc.check("every window excludes no dependence and admits first order",
              bool((sweep.order > 0.5).all() and (sweep.order < 1.6).all()),
              f"{sweep.order.round(2).to_list()}")

    drivers = slowdown.deceleration_drivers(frame)
    doc.claim("the clock coefficient",
              f"**{drivers['span']:+.3f} +/- {drivers['span_stderr']:.3f}**")
    doc.claim("the product coefficient",
              f"{drivers['product']:+.3f} +/- "
              f"{drivers['product_stderr']:.3f}")
    doc.check("the clock carries it and the product does not",
              abs(drivers["span"]) > 2 * drivers["span_stderr"]
              and abs(drivers["product"]) < 3 * drivers["product_stderr"],
              f"{drivers['span']:+.3f}, {drivers['product']:+.3f}")

    doc.section("section 6: what the block cannot do")
    doc.check("no curve in the block is enzyme-free",
              float(frame.e0.min()) > 0, f"min [enz] {frame.e0.min():g}")
    doc.check("and the saved fits share no experiment with it",
              not (set(_fitted_experiments())
                   & set(int(e) for e in scope.TWO_AXIS_BLOCK)),
              f"{sorted(_fitted_experiments())}")

    doc.section("the figures the document promises")
    doc.figures(os.path.join(HERE, "index.html"), "ABCDEFGHIJ")
    doc.claim("the document's own count of them", "ten figures, A\nto J")

    doc.section("the curves page draws every cuvette it claims to")
    page = io.open(os.path.join(HERE, "progress_curves.html"),
                   encoding="utf-8").read()
    drawn = page.count("<div class='fig panel'>")
    doc.check("one panel per cuvette in the block",
              drawn == len(frame), f"{drawn} panels, {len(frame)} cuvettes")
    doc.check("and it draws the dead ones too",
              f"{int((~frame.live).sum())} of {len(frame)} curves" in page)
    # Every contaminated panel has to SHOW the correction, not just have one
    # computed. A page that quietly stopped drawing it would show a fit to the
    # gas again, which is exactly the state this page was in until 2026-09-03.
    chopped = frame[frame.bubble_drops > 0]
    doc.check("every contaminated panel names its detachments",
              page.count("O₂ detachment") == len(chopped),
              f"{page.count('O₂ detachment')} panels against "
              f"{len(chopped)} contaminated curves")
    doc.check("and the corrected series is drawn on each of them",
              page.count("stroke-dasharray='3 2'") == len(chopped),
              f"{page.count(chr(39).join(['stroke-dasharray=', '3 2', '']))} "
              f"against {len(chopped)}")
    beyond = frame[frame.bubble_load > scope.BUBBLE_LOAD_CEILING]
    doc.check("and the curves with no measurable rate say so on the page",
              page.count("NO MEASURABLE RATE") == len(beyond),
              f"{page.count('NO MEASURABLE RATE')} against {len(beyond)}")
    # The counts the preamble quotes are over EVERY curve the page draws, not
    # over the live ones section 5 reports -- the page draws the dead curves
    # too, and one of them carries a load past the ceiling.
    doc.claim("how many panels carry the correction",
              f"**{len(chopped)} curves carrying O₂ detachments**")
    doc.claim("and how many are flagged on their own face",
              f"**{len(beyond)} panels**")
    doc.check("but none of them is missing from the page",
              drawn == len(frame))

    doc.section("the figures: no data point drawn outside its own frame")
    doc.unclipped(os.path.join(HERE, "index.html"),
                  os.path.join(HERE, "progress_curves.html"))
    return doc.summary()


def _fitted_experiments():
    """Every experiment any saved mechanism fit was staged on."""
    import glob
    import json

    found = set()
    for path in glob.glob(os.path.join(os.path.dirname(HERE), "data", "fits",
                                       "*.json")):
        saved = json.load(io.open(path, encoding="utf-8"))
        for key, stage in saved.items():
            if isinstance(stage, dict):
                found.update(int(e) for e in stage.get("experiments", []))
    return found


if __name__ == "__main__":
    raise SystemExit(main())
