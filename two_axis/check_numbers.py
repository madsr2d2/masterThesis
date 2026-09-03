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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import induction
import scope
import slowdown
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
    doc.figures(os.path.join(HERE, "index.html"), "ABCDEF")
    doc.claim("the document's own count of them", "six figures, A to\nF")

    doc.section("the curves page draws every cuvette it claims to")
    page = io.open(os.path.join(HERE, "progress_curves.html"),
                   encoding="utf-8").read()
    drawn = page.count("<div class='fig panel'>")
    doc.check("one panel per cuvette in the block",
              drawn == len(frame), f"{drawn} panels, {len(frame)} cuvettes")
    doc.check("and it draws the dead ones too",
              f"{int((~frame.live).sum())} of {len(frame)} curves" in page)

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
