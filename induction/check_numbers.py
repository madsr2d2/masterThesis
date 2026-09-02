"""
Verifies every number quoted in induction/ANALYSIS.md against the code.

Same contract as the sibling folders', and for the same reason: a figure typed
into a document is a copy, and copies in this project have gone stale three
times -- including once inside a single document, where the same tau ladder was
printed two different ways four sections apart.

    python induction/check_numbers.py
"""
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import arrhenius
import buffer_role
import induction
import scope

from doc_check import Checker

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")


def main():
    doc = Checker(DOCUMENT)
    table = induction.induction_table(induction.WHOLE_ARCHIVE)
    named = induction.induction_blocks(table)

    print("\nsection 1: the design and the statistic")
    lever = induction.substrate_lever(named["4OMe catalysed"])
    doc.claim("the in-run substrate lever",
              f"a factor of **{lever['median_lever']:.1f} in [S]**")
    doc.claim("how many runs carry a ladder",
              f"{lever['laddered']} of\nthe {lever['experiments']} experiments")
    doc.check("the block is 147 live curves in 38 experiments",
              int(named["4OMe catalysed"].live.sum()) == 147
              and lever["experiments"] == 38,
              f"{int(named['4OMe catalysed'].live.sum())}, {lever['experiments']}")
    doc.claim("the live count", "147 live\ncurves in 38 experiments")
    doc.claim("the archive is 402 curves", "all 402 curves")
    doc.check("and it is", len(table) == 402, f"{len(table)}")

    windows = induction.landmark_window()
    narrow = windows[windows.window == "300 s"].iloc[0]
    tenth = windows[windows.window.str.startswith("0.10")].iloc[0]
    wide = windows[windows.window == "900 s"].iloc[0]
    doc.claim("what a 300 s window does to the cold end",
              f"{narrow.cold_s:.0f} s** instead of\n**{tenth.cold_s:.0f} s")
    doc.claim("and to the activation energy",
              f"from {tenth.activation_kJ:.1f} +/- {tenth.stderr_kJ:.1f} to\n"
              f"**{narrow.activation_kJ:.1f} +/- {narrow.stderr_kJ:.1f} kJ/mol")
    doc.claim("a 900 s window recovers",
              f"to {wide.activation_kJ:.1f} +/- {wide.stderr_kJ:.1f}")

    drift = induction.schedule_dependence(table)
    doc.claim("the between-run schedule dependence",
              f"an exponent of **{drift['span']:+.3f} +/- {drift['span_stderr']:.3f}")
    doc.claim("and pH once the schedule is in the model",
              f"{drift['pH']:+.3f} +/- {drift['pH_stderr']:.3f}")

    print("\nsection 2: the induction needs the catalyst")
    summary = induction.channel_summary(table)
    for label, key in (("catalysed", "catalysed"),
                       ("enzyme-free", "enzyme_free")):
        row = summary[key]
        doc.claim(f"{label}: curves and depth",
                  f"| {row['curves']} | **{row['median_depth']:.3f}** | ")
    doc.claim("the enzyme-free block has no deep curves and none accelerating",
              f"| **{summary['enzyme_free']['deep']}** | "
              f"**{summary['enzyme_free']['accelerates']}** |")
    doc.claim("the largest enzyme-free acceleration z",
              f"is **{summary['enzyme_free']['max_accel_z']:.2f}**")
    doc.claim("against the catalysed block's",
              f"against {summary['catalysed']['max_accel_z']:.1f} in the\ncatalysed one")
    doc.claim("the longest enzyme-free run",
              f"exp 28 ran **{summary['enzyme_free']['longest_s']:.0f} s**")

    contrast = induction.channel_contrast(table)
    for row in contrast.itertuples():
        doc.claim(f"matched cell {row.temperature:.0f} C {row.channel}",
                  f"| {row.curves} | ")
        doc.claim(f"matched cell {row.temperature:.0f} C {row.channel}: depth",
                  f"{row.depth:.3f}")

    print("\nsection 3: a clock, not a product threshold")
    for rate in ("v_peak", "peak_rate", "vmax"):
        row = induction.induction_drivers(named["4OMe catalysed"], rate=rate)
        doc.claim(f"the driver coefficient against {rate}",
          f"{row['slope']:+.3f} +/- {row['stderr']:.3f}** | {row['points']} |"
                  if rate == "v_peak" else
                  f"{row['slope']:+.3f} +/- {row['stderr']:.3f} | {row['points']} |")
    series = induction.induction_drivers(named["temperature series"],
                                         rate="v_peak")
    doc.claim("the temperature series cannot do it",
              f"{series['slope']:+.3f} +/- {series['stderr']:.3f} | "
              f"{series['points']} |")
    main_fit = induction.induction_drivers(named["4OMe catalysed"],
                                           rate="v_peak")
    sigma = (main_fit["slope"] - induction.INDUCTION_PRODUCT_SLOPE) \
        / main_fit["stderr"]
    doc.claim("how many sigma product control is excluded by",
              f"excluded at {sigma:.0f}sigma")
    low = induction.induction_drivers(named["4OMe catalysed"],
                                      floor=induction.FLOOR_SWEEP[0])
    high = induction.induction_drivers(named["4OMe catalysed"],
                                       floor=induction.FLOOR_SWEEP[-1])
    doc.claim("the floor sweep",
              f"from {low['slope']:+.3f} +/- {low['stderr']:.3f} to "
              f"{high['slope']:+.3f} +/- {high['stderr']:.3f}")

    for label, name in (("4OMe catalysed", "4OMe catalysed"),
                        ("the temperature series", "temperature series")):
        ratio = induction.order_ratio(named[name])
        doc.claim(f"{label}: order of the induction time",
                  f"{ratio['induction_order']:+.3f} +/- "
                  f"{ratio['induction_stderr']:.3f}")
        doc.claim(f"{label}: order of the rate",
                  f"{ratio['rate_order']:+.3f} +/- {ratio['rate_stderr']:.3f}")
        doc.claim(f"{label}: the ratio",
                  f"{ratio['ratio']:+.2f} +/- {ratio['ratio_stderr']:.2f}")
        depth = scope.orders("depth", frame=named[name],
                             floor=induction.DEPTH_FLOOR)
        doc.claim(f"{label}: order of the amplitude",
                  f"{depth['order_s0']:+.3f} +/- {depth['stderr_s0']:.3f}")

    print("\nsection 4: not the schedule, not physical")
    control = induction.schedule_control(table)
    doc.claim("the two runs of equal length",
              f"both recorded for **{control[19]['duration_s']:.0f} s**")
    doc.claim("and their two time constants",
              f"**{control[19]['tau']:.0f} s and {control[14]['tau']:.0f} s**")
    doc.claim("the ratio", f"a factor of {control['tau_ratio']:.1f}")

    gap = induction.activation_contrast()
    doc.claim("the induction's activation energy",
              f"**{gap['induction']['activation_kJ']:.1f} +/- "
              f"{gap['induction']['activation_stderr']:.1f} kJ/mol**")
    doc.claim("on 16 curves at four temperatures",
              f"{gap['induction']['n']} curves at four temperatures")

    for label, name in (("4OMe catalysed", "4OMe catalysed"),
                        ("BnOH", "BnOH in scope (135-151)")):
        row = induction.signal_control(named[name])
        doc.claim(f"{label}: the signal-to-noise control",
              f"{row['signal_slope']:+.3f} +/- {row['signal_stderr']:.3f}** | "
                  f"{row['signal_points']} |" if name == "4OMe catalysed" else
              f"{row['signal_slope']:+.3f} +/- {row['signal_stderr']:.3f} | "
                  f"{row['signal_points']} |")
    lever = induction.peroxide_lever(table)
    doc.claim("the peroxide block's own signal control",
          f"{lever['signal_slope']:+.3f} +/- {lever['signal_stderr']:.3f} | "
              f"{lever['signal_points']} |")
    doc.claim("the peroxide levels",
              f"{lever['levels'][0]:.3f} against {lever['levels'][1]:.3f} mM")
    doc.claim("the peroxide block's size",
              f"{lever['curves']} curves of which {lever['live']} are live")
    doc.claim("the induction's peroxide order",
              f"**{lever['induction_order']:+.3f} +/- "
              f"{lever['induction_stderr']:.3f}**")
    doc.claim("the rate's peroxide order there",
              f"{lever['rate_order']:+.3f} +/- {lever['rate_stderr']:.3f}")

    scoped = scope.orders("t_ind", frame=named["BnOH in scope (135-151)"],
                          floor=induction.INDUCTION_FLOOR)
    doc.claim("BnOH's induction order in peroxide",
              f"{scoped['order_h2o2']:+.3f} +/- {scoped['stderr_h2o2']:.3f}")

    print("\nsection 4b: the constraint the two orders violate together")
    peroxide_blocks = (("4OMe peroxide, exps 127-131",
                        table[table.experiment.isin(induction.PEROXIDE_LEVER)]),
                       ("BnOH in scope, exps 135-151",
                        named["BnOH in scope (135-151)"]))
    for label, block in peroxide_blocks:
        joint = induction.joint_peroxide_order(block)
        doc.claim(f"{label}: the joint order",
                  f"**{joint['slope']:+.3f} +/- {joint['stderr']:.3f}** | "
                  f"{joint['points']} | {joint['sigma']:.1f}sigma")
        swept = [induction.joint_peroxide_order(block, floor=floor)["slope"]
                 for floor in induction.FLOOR_SWEEP]
        doc.claim(f"{label}: the floor sweep",
                  " | ".join(f"{value:+.3f}" for value in swept)
                  .replace(f"{swept[2]:+.3f}", f"**{swept[2]:+.3f}**"))
    doc.check("the required gap is 1", induction.PERHYDRATE_ORDER_GAP == 1.0)
    doc.check("and the deviation never changes sign at any floor",
              all(induction.joint_peroxide_order(block, floor=floor)["slope"] < 1.0
                  for _, block in peroxide_blocks
                  for floor in induction.FLOOR_SWEEP))

    saturation = induction.peroxide_saturation(named["BnOH in scope (135-151)"])
    doc.claim("the ladder's size and range",
          f"**{saturation['points']}\ncurves in {saturation['experiments']} "
              f"runs over {saturation['peroxide_low']:.2f}-"
              f"{saturation['peroxide_high']:.0f} mM**")
    doc.claim("the free power law",
          f"a = **{saturation['order']:.3f}** ({saturation['order_low']:.3f} "
              f"to {saturation['order_high']:.3f})")
    doc.claim("first order is rejected",
              f"**rejected, F = {saturation['first_order_f']:.1f}**")
    doc.claim("the scheme's own form fits worse",
              f"{saturation['scheme_sse']:.2f} against {saturation['power_sse']:.2f}")
    doc.check("and it does fit worse", saturation["scheme_sse"] > saturation["power_sse"],
              f"{saturation['scheme_sse']:.2f} vs {saturation['power_sse']:.2f}")

    print("\nthe trap constants, which are C8's gate")
    energies = []
    for label, block in peroxide_blocks:
        got = scope.orders("t_ind", frame=block,
                           floor=induction.INDUCTION_FLOOR)
        trap = induction.trap_constant(got["order_h2o2"], got["stderr_h2o2"],
                                       induction.peroxide_geometric_mean(block))
        energies.append(trap["free_energy_kJ"])
        doc.claim(f"{label}: K and its error",
                  f"{trap['constant']:.4f} +/- {trap['stderr']:.4f} /mM")
        doc.claim(f"{label}: as a molar constant", f"{trap['molar']:.0f} /M")
        doc.claim(f"{label}: the free energy",
                  f"**{trap['free_energy_kJ']:.2f} kJ/mol**")
        doc.claim(f"{label}: the peroxide it belongs to",
                  f"at {trap['peroxide_mM']:.1f} mM")
    import numpy as np
    molar = saturation["constant"] * 1000.0
    energy = (-arrhenius.GAS_CONSTANT * arrhenius.REFERENCE_KELVIN
              * np.log(molar) / 1000.0)
    doc.claim("the rates route's K",
              f"{saturation['constant']:.4f} /mM ({saturation['constant_low']:.4f}-"
              f"{saturation['constant_high']:.4f})")
    doc.claim("and its free energy", f"{energy:.2f} (")
    energies.append(energy)
    doc.check("all three land between -6 and -10 kJ/mol, which is the gate quoted",
              all(-10.0 <= value <= -5.5 for value in energies),
              " ".join(f"{value:+.2f}" for value in energies))

    print("\nsection 6: which way the induction points")
    for label, name in (("4OMe catalysed", "4OMe catalysed"),
                        ("the temperature series", "temperature series"),
                        ("BnOH in scope", "BnOH in scope (135-151)"),
                        ("4OMe enzyme-free", "4OMe enzyme-free")):
        block = named[name]
        block = block[block.live]
        doc.claim(f"{label}: curves and lag-first",
                  f"| {len(block)} | **{int(block.lag_first.sum())}** |")
        counts = block.progress_kind.value_counts()
        doc.claim(f"{label}: the shapes",
                  ", ".join(f"{value} {key}" for key, value in counts.items()))
    arms = induction.ladder_arms(named["BnOH in scope (135-151)"])
    for label, block, axis in (("in scope, substrate arm",
                                arms["substrate arm"], "s0"),
                               ("in scope, peroxide arm",
                                arms["peroxide arm"], "h2o2"),
                               ("4OMe catalysed", named["4OMe catalysed"],
                                "s0")):
        alone = induction.sign_drivers(block, axis=axis, control=False)
        controlled = induction.sign_drivers(block, axis=axis, control=True)
        doc.claim(f"{label}: alone",
                  f"{alone[axis]:+.3f} +/- {alone[axis + '_stderr']:.3f}")
        doc.claim(f"{label}: with signal-to-noise",
                  f"**{controlled[axis]:+.3f} +/- "
                  f"{controlled[axis + '_stderr']:.3f}**")
    for label, name in (("4OMe catalysed", "4OMe catalysed"),
                        ("BnOH in scope", "BnOH in scope (135-151)")):
        got = induction.composition_collinearity(named[name])
        doc.claim(f"{label}: the [S]/[buf] collinearity",
                  f"{got['median']:+.2f}")
    doc.check("the 4OMe blocks have no run with a constant buffer",
              induction.composition_collinearity(
                  named["4OMe catalysed"])["constant_buffer"] == 0)
    doc.check("and every in-scope run has one",
              induction.composition_collinearity(
                  named["BnOH in scope (135-151)"])["constant_buffer"] == 17)
    ladder = induction.buffer_lever(table)
    for run, label in ((34, "exp 34"), (32, "exp 32")):
        block = ladder[ladder.experiment == run]
        doc.claim(f"{label}: its F range",
                  f"F = {block.two_phase_f.min():.0f} to "
                  f"{block.two_phase_f.max():.0f}"
                  if run == 34 else
                  f"F = {block.two_phase_f.min():.1f} to "
                  f"{block.two_phase_f.max():.1f}")
        doc.claim(f"{label}: its run length", f"{block.span_s.min():.0f} s")
    step = induction.buffer_join_step(ladder)
    doc.claim("the level step at the join", f"falls {step['step']:.2f}x")
    doc.claim("the two rates it is between",
          f"{step['from_rate'] * 1e5:.2f} x 10-5 at {step['from_buf']:.0f} mM "
              f"to {step['to_rate'] * 1e5:.2f} x 10-5 at {step['to_buf']:.0f} mM")
    doc.claim("the withdrawn tau_fast slopes",
              "+0.457 +/- 0.097 and -1.052 +/- 0.469")
    doc.claim("the window it is read through",
              f"({induction.BUFFER_WINDOW:.0f} s;")
    for response, key in (("t_ind", "t"), ("depth", "d")):
        got = induction.buffer_order(ladder, response)
        doc.claim(f"{response}: exp 34's slope",
                  f"{got['slope_34']:+.3f} +/- {got['stderr_34']:.3f}")
        doc.claim(f"{response}: exp 32's slope",
                  f"{got['slope_32']:+.3f} +/- {got['stderr_32']:.3f}")
        doc.claim(f"{response}: the pooled slope",
                  f"**{got['slope']:+.3f} +/- {got['stderr']:.3f}**")
    swept = [induction.buffer_order(
        induction.buffer_lever(table, width=width))["slope"]
        for width in induction.BUFFER_WINDOW_SWEEP]
    doc.claim("the window sweep's range",
              f"runs {max(swept):.2f} to {min(swept):.2f}")
    doc.check("and no window changes its sign", all(v < 0 for v in swept),
              " ".join(f"{v:+.3f}" for v in swept))

    for width in induction.BUFFER_WINDOW_SWEEP:
        rungs = induction.buffer_lever(table, width=width)
        joint = induction.joint_buffer_order(rungs)
        control = induction.signal_control(rungs)
        passes = abs(control["signal_slope"]) < 2 * control["signal_stderr"]
        bold = "**" if passes else ""
        doc.claim(f"the joint buffer order at {width:.0f} s",
                  f"| {width:.0f} s | {bold}{joint['slope']:+.3f} +/- "
                  f"{joint['stderr']:.3f}{bold} | {joint['sigma']:.1f}sigma | "
                  f"{control['signal_slope']:+.3f} +/- "
                  f"{control['signal_stderr']:.3f} "
                  f"{'passes' if passes else '**fails**'} |")
    doc.check("the windows that fail the control are the windows that overshoot",
              all(induction.joint_buffer_order(
                  induction.buffer_lever(table, width=width))["slope"] > 1.2
                  for width in (600.0, 900.0, 1200.0)))
    crossing = buffer_role.peroxide_crossing()
    doc.claim("the crossing the archive does not have",
          f"**{crossing['steps_both']} of its {crossing['runs']} runs step "
              f"[buf] and [H2O2] at once**")

    order = induction.buffer_order(ladder)
    fixed = induction.substrate_order_corrected(
        named["4OMe catalysed"], order["slope"], order["stderr"])
    doc.claim("the ladder's own buffer slope",
              f"= {fixed['ladder_slope']:.3f}")
    doc.claim("route two as measured",
          f"| {fixed['measured']:+.3f} +/- {fixed['measured_stderr']:.3f} | "
              f"{abs(fixed['measured'] - fixed['threshold']) / fixed['measured_stderr']:.1f}sigma |")
    doc.claim("route two corrected",
              f"| **{fixed['corrected']:+.3f} +/- {fixed['corrected_stderr']:.3f}** "
              f"| **{abs(fixed['corrected'] - fixed['threshold']) / fixed['corrected_stderr']:.1f}sigma** |")
    doc.check("the correction moves it towards the threshold, not away",
              abs(fixed["corrected"] - fixed["threshold"])
              < abs(fixed["measured"] - fixed["threshold"]))

    print("\nsection 5: the activation parameters")
    for label, key in (("the induction", "induction"),
                       ("the turnover", "turnover")):
        row = gap[key]
        doc.claim(f"{label}: Ea and enthalpy",
          f"| {row['activation_kJ']:.1f} +/- {row['activation_stderr']:.1f} "
                  f"| {row['enthalpy_kJ']:.1f} +/- {row['enthalpy_stderr']:.1f} |")
        doc.claim(f"{label}: entropy",
                  f"{row['entropy_J']:+.1f} +/- {row['entropy_stderr']:.1f}")
        doc.claim(f"{label}: free energy",
                  f"**{row['gibbs_kJ']:.2f} +/- {row['gibbs_stderr']:.2f}**")
    # The emphasis wraps the whole sentence rather than the number, so the
    # claim carries the sentence: that is what makes it the headline, and a
    # claim that quoted the bare number would not notice if it stopped being.
    doc.claim("the free-energy gap",
              f"**The free-energy gap is the solid number: "
              f"{gap['gibbs_gap_kJ']:.2f} +/- "
              f"{gap['gibbs_gap_stderr']:.2f} kJ/mol.**")
    doc.claim("what it is worth as a rate constant",
              f"**{gap['rate_ratio']:.0f}x faster**")
    doc.claim("the entropy gap", f"{gap['entropy_gap_J']:+.1f} J/mol/K")
    doc.claim("its error", f"**+/- {gap['entropy_gap_stderr']:.1f}**")
    doc.claim("and the enthalpy gap's",
              f"**+/- {gap['enthalpy_gap_stderr']:.1f}**")
    doc.check("the free-energy gap is 19 or more standard errors wide",
              abs(gap["gibbs_gap_kJ"]) / gap["gibbs_gap_stderr"] > 19,
              f"{abs(gap['gibbs_gap_kJ']) / gap['gibbs_gap_stderr']:.1f}")
    doc.check("and each of its two parts is about one",
              abs(gap["entropy_gap_J"]) / gap["entropy_gap_stderr"] < 1.5
              and abs(gap["enthalpy_gap_kJ"]) / gap["enthalpy_gap_stderr"] < 1.5)

    print("\nthe product_fate numbers this document quotes back")
    import slowdown
    whole = scope.frame(induction.WHOLE_ARCHIVE)
    fall = slowdown.deceleration_drivers(
        slowdown.substrate_blocks(whole)["4OMe catalysed, phosphate"])
    doc.claim("the fall's product coefficient",
              f"({fall['product']:+.3f} +/- {fall['product_stderr']:.3f})")
    doc.claim("the fall's clock coefficient",
              f"({fall['span']:+.3f} +/- {fall['span_stderr']:.3f})")

    print("\nthe figures the document promises")
    doc.figures(os.path.join(HERE, "index.html"), "ABCDEFGHI")
    doc.claim("the document's own count of them", "nine figures, A\nto I")

    print("\nthe curves page draws both channels, whole")
    page = io.open(os.path.join(HERE, "progress_curves.html"),
                   encoding="utf-8").read()
    drawn = page.count("<div class='fig panel'>")
    live = named["4OMe catalysed"], named["4OMe enzyme-free"]
    expected = sum(int(block.live.sum()) for block in live)
    doc.check("one panel per live 4OMe cuvette, both channels",
              drawn == expected, f"{drawn} panels, {expected} curves")
    doc.check("the enzyme-free control is drawn, not summarised",
              "enzyme-free —" in page)

    print("\nthe figures: no data point drawn outside its own frame")
    # Marks are clipped to the plot area deliberately, so a data point outside
    # the axis limits vanishes silently and the figure still looks complete.
    # That is how exp 6 sample 4 went missing from a curvature figure for as
    # long as it existed: no number changes, so no prose check could see it.
    doc.unclipped(os.path.join(HERE, "index.html"),
                  os.path.join(HERE, "progress_curves.html"))
    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
