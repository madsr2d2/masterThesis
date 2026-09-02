"""
Verifies every number quoted in buffer/ANALYSIS.md against the code.

The comparison itself is `doc_check`, shared with every sibling folder: a figure
typed into a document is a copy, and copies in this project have gone stale.

    python buffer/check_numbers.py
"""
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import buffer_role
import induction
import scope
from doc_check import Checker

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")


def main():
    doc = Checker(DOCUMENT)
    frame = scope.frame(tuple(range(1, 152)))

    print("\nsection 1: what a titration at one pH can ask")
    prediction = buffer_role.species_prediction(7.00, 7.53)
    doc.claim("the pKa", f"is **{prediction['pka']:.2f}**")
    doc.claim("the two base fractions",
              f"{prediction['low_base_fraction']:.3f} at pH 7.00 and "
              f"{prediction['high_base_fraction']:.3f} at pH 7.53")
    doc.claim("what general base predicts",
              f"| **{prediction['general_base']:.2f}** |")
    doc.claim("what general acid predicts",
              f"| **{prediction['general_acid']:.2f}** |")

    print("\nsection 2: the five titrations")
    table = buffer_role.titration_table()
    doc.claim("how many curves", f"**{len(table)} curves at three pH values**"
              if False else f"**{int(table.curves.sum())} curves at three pH values**")
    for row in table.itertuples():
        doc.claim(f"exp {row.experiment}: the row",
                  f"| {row.experiment} | {row.pH:.2f} | {row.s0:.3f} | "
                  f"{row.buffer_low:g}-{row.buffer_high:g} |")
        bold = "**" if abs(row.order) > 2 * row.order_stderr else ""
        doc.claim(f"exp {row.experiment}: its order",
                  f"{bold}{row.order:+.3f} +/- {row.order_stderr:.3f}{bold} | "
                  f"{row.order_r2:.3f} |")
    doc.check("five runs, twenty curves",
              len(table) == 5 and int(table.curves.sum()) == 20,
              f"{len(table)}, {int(table.curves.sum())}")
    low = table.set_index("experiment")
    doc.claim("the coefficient at pH 7.00",
              f"is 1.86 x 10-7 at pH 7.00")
    doc.check("and it is", abs(low.loc[32].coefficient - 1.86e-7) < 5e-10,
              f"{low.loc[32].coefficient:.3g}")
    doc.claim("the coefficients at pH 7.53",
              f"3.19-3.23 x 10-7 at pH 7.53")
    doc.check("and they are",
              abs(low.loc[36].coefficient - 3.19e-7) < 5e-10
              and abs(low.loc[37].coefficient - 3.23e-7) < 5e-10,
              f"{low.loc[36].coefficient:.3g}, {low.loc[37].coefficient:.3g}")
    doc.claim("the intercepts", "from an intercept of 3.29 x 10-5 to "
                                "2.73-4.91 x 10-4")
    doc.check("and they are",
              abs(low.loc[32].intercept - 3.29e-5) < 5e-8
              and abs(low.loc[37].intercept - 2.73e-4) < 5e-7
              and abs(low.loc[36].intercept - 4.91e-4) < 5e-7,
              f"{low.loc[32].intercept:.3g}, {low.loc[37].intercept:.3g}, "
              f"{low.loc[36].intercept:.3g}")

    print("\nsection 3: which species")
    for label, drop, order in (("all four runs", (), buffer_role.SUBSTRATE_ORDER),
                               ("without exp 35", (35,),
                                buffer_role.SUBSTRATE_ORDER),
                               ("unnormalised", (), 0.0)):
        got = buffer_role.catalytic_coefficient(drop=drop,
                                                substrate_order=order)
        bold = "**" if label == "without exp 35" else ""
        doc.claim(f"{label}: the ratio",
                  f"| {bold}{got['ratio']:+.2f} +/- {got['ratio_stderr']:.2f}{bold} |")
    measured = buffer_role.catalytic_coefficient(drop=(35,))
    verdict = buffer_role.separable(measured, prediction)
    for name, label in (("general_base", "general base"),
                        ("general_acid", "general acid"),
                        ("spectator", "spectator")):
        row = verdict[name]
        doc.claim(f"{label}: its distance",
          f"| {label} | {row['predicted']:.2f} | {row['sigma']:.1f}sigma | "
                  f"survives |")
    doc.check("nothing is excluded", len(verdict["survivors"]) == 3,
              f"{verdict['survivors']}")
    doc.claim("exp 35's own R2", f"has R2 {low.loc[35].order_r2:.3f}")
    doc.check("which is why it is dropped", low.loc[35].order_r2 < 0.2)

    print("\nsection 4: the buffer as a confound")
    blocks = induction.induction_blocks(frame)
    for label, name in (("4OMe catalysed", "4OMe catalysed"),
                        ("4OMe enzyme-free", "4OMe enzyme-free"),
                        ("the temperature series", "temperature series"),
                        ("BnOH in scope", "BnOH in scope (135-151)")):
        got = induction.composition_collinearity(blocks[name])
        # The document bolds the two rows that carry the CONTRAST -- every
        # 4OMe run against the in-scope block -- and not every row whose
        # value happens to clear a threshold. Emphasis is an argument, so it
        # is named here rather than inferred from the number. This check was
        # dead until 2026-09-02: the old normaliser stripped `**` from both
        # sides before comparing, so the markers were built and discarded.
        bold = "**" if name in ("4OMe catalysed",
                                "BnOH in scope (135-151)") else ""
        # A zero collinearity is written without a sign in the document, which
        # is the right typography and needs saying here rather than there.
        shown = f"{got['median']:.2f}" if got["median"] == 0.0 \
            else f"{got['median']:+.2f}"
        doc.claim(f"{label}: runs and collinearity",
                  f"| {got['runs']} | {bold}{shown}{bold} |")
        if np.isfinite(got.get("slope", np.nan)):
            doc.claim(f"{label}: the ladder slope", f"{got['slope']:.3f} |")
    doc.check("in-scope every run holds [buf] constant",
              induction.composition_collinearity(
                  blocks["BnOH in scope (135-151)"])["constant_buffer"] == 17)

    print("\nsection 5: identity against pH")
    identity = buffer_role.identity_overlap(frame)
    for row in identity.itertuples():
        span = (f"{row.pH_low:.2f}" if row.pH_low == row.pH_high
                else f"{row.pH_low:.2f}-{row.pH_high:.2f}")
        doc.claim(f"{row.substrate} {row.channel} {row.buffer}",
                  f"| {row.curves} | {span} |")
    widest = buffer_role.overlap_width(identity)
    doc.claim("the widest overlap", f"**{widest['width']:.2f} units**")
    doc.claim("the phosphate block's peroxide",
              f"sits at exactly\n{widest['peroxide'][0][0]:g} mM")
    doc.claim("the pyrophosphate block's",
              f"at exactly {widest['peroxide'][1][0]:g} and "
              f"{widest['peroxide'][1][1]:g} mM")
    doc.check("and the two cells share no peroxide value",
              not widest["shares_peroxide"], f"{widest['peroxide']}")

    print("\nsection 6: a base, or a carrier for the peroxide")
    crossing = buffer_role.peroxide_crossing(frame=frame)
    for label, key in (("runs in the archive", "runs"),
                       ("runs that step `[buf]`", "steps_buffer"),
                       ("runs that step `[H2O2]`", "steps_peroxide")):
        doc.claim(f"the design: {label}", f"| {crossing[key]} |")
    doc.claim("and none does both", f"| **{crossing['steps_both']}** |")
    doc.check("no run crosses the two ladders", crossing["steps_both"] == 0)
    doc.claim("the titrations' one peroxide",
              f"exactly **{crossing['titration_peroxides'][0]:g} mM** peroxide")

    ladder = induction.buffer_lever(induction.induction_table(
        induction.WHOLE_ARCHIVE))
    for width in induction.BUFFER_WINDOW_SWEEP:
        rungs = induction.buffer_lever(induction.induction_table(
            induction.WHOLE_ARCHIVE), width=width)
        joint = induction.joint_buffer_order(rungs)
        control = induction.signal_control(rungs)
        passes = abs(control["signal_slope"]) < 2 * control["signal_stderr"]
        bold = "**" if passes else ""
        doc.claim(f"the joint order at {width:.0f} s",
              f"| {bold}{joint['slope']:+.3f} +/- {joint['stderr']:.3f}{bold} "
                  f"| {joint['sigma']:.1f}sigma | "
                  f"{control['signal_slope']:+.3f} +/- "
                  f"{control['signal_stderr']:.3f} "
                  f"{'passes' if passes else '**fails**'} |")
    doc.check("the two windows that pass the control straddle +1",
              all(abs(induction.joint_buffer_order(induction.buffer_lever(
                  induction.induction_table(induction.WHOLE_ARCHIVE),
                  width=width))["slope"] - 1.0) < 0.15 for width in (300.0, 450.0)))

    free = buffer_role.free_route_order(
        frame=scope.frame(buffer_role.TITRATIONS))
    doc.claim("the buffer-free level's rise",
              f"2.65 x 10-5 -> 2.09 x 10-4, **{free['level_ratio']:.2f}x**")
    doc.claim("what [HOO-] gives instead",
              f"**{free['hoo_ratio']:.2f}x** (order {free['hoo_order']:+.2f} "
              f"in [OH-]")
    doc.claim("the level's apparent order",
              f"**{free['apparent_order']:+.2f}**, not "
              f"{free['hoo_order']:+.2f}")
    excess = free["level_ratio"] / free["hoo_ratio"]
    doc.claim("the excess", f"**{excess:.2f}x more than one hydroperoxide")
    doc.check("the two runs are matched in substrate", free["matched_s0"])

    step = induction.buffer_join_step(ladder)
    doc.claim("the between-day step it has to be read against",
              f"is **{step['step']:.2f}x**")

    print("\nthe induction's buffer order, quoted back from induction/")
    order = induction.buffer_order(
        induction.buffer_lever(induction.induction_table(
            induction.WHOLE_ARCHIVE)))
    doc.claim("it", f"{order['slope']:+.3f} +/- {order['stderr']:.3f}")

    print("\nthe figures the document promises")
    doc.figures(os.path.join(HERE, "index.html"), "ABCDE")
    doc.claim("the document's own count of them", "five figures, A to\nE")

    print("\nthe curves page draws every cuvette it claims to")
    page = io.open(os.path.join(HERE, "progress_curves.html"),
                   encoding="utf-8").read()
    drawn = page.count("<div class='fig panel'>")
    live = scope.frame(buffer_role.TITRATIONS)
    doc.check("one panel per titration cuvette",
              drawn == len(live), f"{drawn} panels, {len(live)} cuvettes")
    doc.check("and it names the window it reads the landmark through",
              f"{induction.BUFFER_WINDOW:.0f} seconds" in page)

    print("\nthe figures: no data point drawn outside its own frame")
    doc.unclipped(os.path.join(HERE, "index.html"),
                  os.path.join(HERE, "progress_curves.html"))
    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
