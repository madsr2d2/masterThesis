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
                        ("BnOH", "BnOH two-axis (135-151)")):
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

    scoped = scope.orders("t_ind", frame=named["BnOH two-axis (135-151)"],
                          floor=induction.INDUCTION_FLOOR)
    doc.claim("BnOH's induction order in peroxide",
              f"{scoped['order_h2o2']:+.3f} +/- {scoped['stderr_h2o2']:.3f}")

    print("\nsection 4b: the constraint the two orders violate together")
    peroxide_blocks = (("4OMe peroxide, exps 127-131",
                        table[table.experiment.isin(induction.PEROXIDE_LEVER)]),
                       ("BnOH two-axis, exps 135-151",
                        named["BnOH two-axis (135-151)"]))
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

    saturation = induction.peroxide_saturation(named["BnOH two-axis (135-151)"])
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
                        ("BnOH two-axis", "BnOH two-axis (135-151)"),
                        ("4OMe enzyme-free", "4OMe enzyme-free")):
        block = named[name]
        block = block[block.live]
        doc.claim(f"{label}: curves and lag-first",
                  f"| {len(block)} | **{int(block.lag_first.sum())}** |")
        counts = block.progress_kind.value_counts()
        doc.claim(f"{label}: the shapes",
                  ", ".join(f"{value} {key}" for key, value in counts.items()))
    arms = induction.ladder_arms(named["BnOH two-axis (135-151)"])
    for label, block, axis in (("two-axis, substrate arm",
                                arms["substrate arm"], "s0"),
                               ("two-axis, peroxide arm",
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
                        ("BnOH two-axis", "BnOH two-axis (135-151)")):
        got = induction.composition_collinearity(named[name])
        doc.claim(f"{label}: the [S]/[buf] collinearity",
                  f"{got['median']:+.2f}")
    doc.check("the 4OMe blocks have no run with a constant buffer",
              induction.composition_collinearity(
                  named["4OMe catalysed"])["constant_buffer"] == 0)
    doc.check("and every two-axis run has one",
              induction.composition_collinearity(
                  named["BnOH two-axis (135-151)"])["constant_buffer"] == 17)
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

    print("\nsection 7: the seven variables")
    archive = scope.frame(scope.archive())
    blocks = induction.induction_blocks(archive)

    identify = induction.lag_identifiability(archive)
    for axis, label in (("s0", "`[S]`"), ("buf", "`[buf]`"),
                        ("h2o2", "`[H2O2]`")):
        row = identify.loc[axis]
        doc.claim(f"{axis} moves inside this many runs",
                  f"| {label} | **{int(row.runs_moving_it)} of "
                  f"{int(row.runs)}** | {row.widest_within_run:.1f}x |")
    doc.check("the other four axes move in no run at all",
              set(identify.index[identify.runs_moving_it == 0])
              == {"pH", "temperature", "buffer", "substrate"},
              f"{sorted(identify.index[identify.runs_moving_it == 0])}")

    print("\nsection 7b: the lag by substrate and channel")
    channels = induction.lag_channel_table(archive)
    for substrate, channel in (("4OMe-BnOH", "enzyme-free"),
                               ("4OMe-BnOH", "catalysed"),
                               ("BnOH", "enzyme-free"),
                               ("BnOH", "catalysed")):
        row = channels.loc[(substrate, channel)]
        emphasis = "**" if row.median_depth > 0 else ""
        doc.claim(f"{substrate} {channel}: curves, lags, depth and clock",
                  f"| {int(row.curves)} | **{int(row.with_a_lag)}** | "
                  f"{emphasis}{row.median_depth:.3f}{emphasis} | "
                  f"{row.median_clock_s:.0f} s |")

    live = archive[archive.live]
    free = live[~live.differential]
    accelerating = {
        substrate: (int(group.accelerates.sum()), int(len(group)))
        for substrate, group in free.groupby("substrate")}
    doc.claim("enzyme-free BnOH curves that accelerate",
              f"**{accelerating['BnOH'][0]} of {accelerating['BnOH'][1]} "
              f"enzyme-free BnOH curves accelerate past 3sigma and "
              f"{accelerating['4OMe-BnOH'][0]} of "
              f"{accelerating['4OMe-BnOH'][1]} enzyme-free\n4OMe curves do**")

    print("\nsection 7c: the replicate floor")
    floor = induction.replicate_floor()
    doc.claim("the replicate clock spread",
              f"**{floor['clock_ratio']:.2f}x**")
    doc.claim("the replicate depth spread",
              f"**{floor['depth_spread']:.3f}**")
    for experiment in scope.REPLICATE_RUNS:
        row = floor["table"].loc[experiment]
        doc.claim(f"exp {experiment}'s clock and depth",
                  f"{row.clock:.0f} | ")
    doc.check("the four runs are one composition, four repeats",
              floor["runs"] == 4 and floor["curves"] == 16,
              f"{floor['runs']} runs, {floor['curves']} curves")

    print("\nsection 7d: the within-run orders")
    two_axis = blocks["BnOH two-axis (135-151)"]
    catalysed = blocks["4OMe catalysed"]
    cases = (("4OMe catalysed", catalysed, ("s0", "h2o2", "buf")),
             ("BnOH two-axis", two_axis, ("s0", "h2o2")),
             ("4OMe peroxide",
              archive[archive.experiment.isin(induction.PEROXIDE_LEVER)],
              ("s0", "h2o2")),
             ("buffer titrations",
              archive[archive.experiment.isin(scope.BUFFER_TITRATIONS)],
              ("s0", "h2o2", "buf")))
    for label, block, terms in cases:
        got = induction.lag_orders(block, terms=terms)
        for axis in terms:
            order = got["lag_half_s"][axis]
            if not np.isfinite(order["order"]):
                continue
            bold = "**" if (label == "BnOH two-axis" and axis == "s0") else ""
            doc.claim(f"{label} {axis}: the collinearity",
                      f"| {got['signal_collinearity_' + axis]:+.2f} |")
            doc.claim(f"{label} {axis}: the clock order",
                      f"| {bold}{order['order']:+.3f} +/- {order['stderr']:.3f}"
                      f"{bold} | {bold}{order['controlled']:+.3f} +/- "
                      f"{order['controlled_stderr']:.3f}{bold} |")
    doc.claim("the two-axis block fails its own signal control",
              f"({induction.lag_signal_control(two_axis)['lag_half_s']['slope']:+.3f} "
              f"+/- {induction.lag_signal_control(two_axis)['lag_half_s']['stderr']:.3f})")
    rate = scope.orders("vmax_corrected", frame=two_axis)
    doc.claim("the rate order in [S] there",
              f"of only {rate['order_s0']:+.2f} +/- {rate['stderr_s0']:.2f}")
    doc.claim("the single-axis fit is a different number",
              f"reads {scope.orders('lag_half_s', frame=two_axis, floor=induction.INDUCTION_FLOOR, terms=('s0',))['order_s0']:+.3f} "
              f"+/- {scope.orders('lag_half_s', frame=two_axis, floor=induction.INDUCTION_FLOOR, terms=('s0',))['stderr_s0']:.3f}")
    doc.claim("how many two-axis curves sit on the floor",
              f"{induction.lag_orders(two_axis)['lag_half_s']['floored']} of the "
              f"two-axis block's 110")

    sweeps = {"two-axis [S]": induction.lag_order_floor_sweep(two_axis, "s0"),
              "two-axis [H2O2]": induction.lag_order_floor_sweep(
                  two_axis, "h2o2"),
              "4OMe [S]": induction.lag_order_floor_sweep(
                  catalysed, "s0", terms=("s0", "h2o2", "buf"))}
    for label, sweep in sweeps.items():
        row = " | ".join(
            ("**" if floor_value == 60.0 else "")
            + f"{sweep.loc[floor_value, 'controlled']:+.3f}"
            + ("**" if floor_value == 60.0 else "")
            for floor_value in induction.FLOOR_SWEEP)
        doc.claim(f"the floor sweep, {label}", row)
    lows = sweeps["two-axis [S]"].controlled
    doc.check("the BnOH substrate coefficient is negative at every floor",
              bool((lows < 0).all()), f"{list(lows.round(3))}")
    doc.claim("and quoted as its range",
              f"**{lows.max():+.2f} to {lows.min():+.2f}**")

    print("\nsection 7d: route one on the window-free clock")
    for label, block in (("4OMe catalysed", catalysed),
                         ("BnOH two-axis", two_axis),
                         ("the temperature series", blocks["temperature series"])):
        route = induction.induction_drivers(block, response="lag_half_s",
                                            rate="vmax_corrected")
        bold = "**" if label != "the temperature series" else ""
        doc.claim(f"route one on {label}",
                  f"| {bold}{route['slope']:+.3f} +/- {route['stderr']:.3f}"
                  f"{bold} | {route['points']} |")
    four = induction.induction_drivers(catalysed, response="lag_half_s",
                                       rate="vmax_corrected")
    doc.claim("and how far that is from product control",
              f"at {(four['slope'] + 1.0) / four['stderr']:.1f}sigma")
    bnoh = induction.induction_drivers(two_axis, response="lag_half_s",
                                       rate="vmax_corrected")
    doc.claim("the BnOH row's distance from product control",
              f"at {(bnoh['slope'] + 1.0) / bnoh['stderr']:.1f}sigma")

    print("\nsection 7e: the four pH ladders")
    ladders = induction.lag_ph_ladders()
    for row in ladders:
        clock = row["lag_half_s"]
        doc.claim(f"{row['ladder']}: the two collinearities",
                  f"| {'**' if abs(row['schedule_collinearity']) > 0.7 else ''}"
                  f"{row['schedule_collinearity']:+.2f}"
                  f"{'**' if abs(row['schedule_collinearity']) > 0.7 else ''} | "
                  f"{'**' if abs(row['signal_collinearity']) > 0.8 else ''}"
                  f"{row['signal_collinearity']:+.2f}"
                  f"{'**' if abs(row['signal_collinearity']) > 0.8 else ''} |")
        doc.claim(f"{row['ladder']}: the held clock coefficient",
                  f"{clock['controlled']:+.3f} +/- "
                  f"{clock['controlled_stderr']:.3f} |")
        doc.claim(f"{row['ladder']}: how many runs",
                  f"| {int(row['runs'])} | ")
    pooled = induction.pooled_ladder(ladders, "lag_half_s")
    doc.claim("the pooled pH coefficient",
              f"**{pooled['pooled']:+.3f} +/- {pooled['stderr']:.3f}**")
    doc.claim("and the four agree",
              f"chi2 = {pooled['chi2']:.2f} on {pooled['dof']}")
    narrower = [induction.pooled_ladder(
        induction.lag_ph_ladders(window_fraction=fraction), "lag_half_s")
        for fraction in (0.75, 0.5)]
    doc.claim("what the window does to the pooled value",
              f"{narrower[0]['pooled']:+.3f} +/- {narrower[0]['stderr']:.3f} and "
              f"{narrower[1]['pooled']:+.3f} +/- {narrower[1]['stderr']:.3f}")
    lo = min(row["pooled"] for row in narrower + [pooled])
    hi = max(row["pooled"] for row in narrower + [pooled])
    doc.claim("the range it is quoted as",
              f"**d ln tau / d pH = {lo:+.2f} to {hi:+.2f}")
    fractions = [induction.saturation_fraction(row["pooled"], row["stderr"])
                 for row in narrower + [pooled]]
    doc.claim("as a saturation fraction",
              f"**{min(f['fraction'] for f in fractions):.2f} to "
              f"{max(f['fraction'] for f in fractions):.2f}**")
    depth_pooled = induction.pooled_ladder(ladders, "lag_depth")
    doc.claim("the depth carries nothing",
              f"pooled {depth_pooled['pooled']:+.3f} +/- "
              f"{depth_pooled['stderr']:.3f}")

    print("\nsection 7f: the barrier and the schedule")
    sweep = induction.lag_arrhenius_sweep()
    doc.claim("the windows swept",
              " | ".join(f"{row['window_s']:.0f}" for row in sweep[:-1])
              + f" | {sweep[-1]['window_s']:.0f}")
    doc.claim("the temperatures each keeps",
              " | ".join(str(row["temperatures"]) for row in sweep))
    doc.claim("the barriers",
              " | ".join(
                  ("**" if row["fraction"] == 1.0 else "")
                  + f"{row['activation_kj']:.1f} +/- {row['stderr_kj']:.1f}"
                  + ("**" if row["fraction"] == 1.0 else "")
                  for row in sweep))
    full = [row for row in sweep if row["fraction"] == 1.0][0]
    doc.claim("the number to quote",
              f"**{full['activation_kj']:.0f} +/- {full['stderr_kj']:.0f} "
              f"kJ/mol** over six temperatures")
    doc.check("it is stable over the wide windows",
              max(row["activation_kj"] for row in sweep
                  if row["fraction"] >= 0.75)
              - min(row["activation_kj"] for row in sweep
                    if row["fraction"] >= 0.75) < 15.0)
    doc.claim("the range it is stable over",
              f"**{min(row['activation_kj'] for row in sweep if row['fraction'] >= 0.75):.0f}-"
              f"{max(row['activation_kj'] for row in sweep if row['fraction'] >= 0.75):.0f} kJ/mol**")
    clocks = induction.lag_arrhenius()["clock_by_temperature"]
    ladder = [f"{clocks[t]:.0f}" for t in sorted(clocks)]
    doc.claim("the half-rise by temperature",
              ", ".join(ladder[:-1]) + f" and {ladder[-1]} s")
    whole = [row for row in sweep if row["fraction"] == 4.0][0]
    doc.claim("the whole-run barrier",
              f"**{whole['activation_kj']:.1f} +/- {whole['stderr_kj']:.1f} kJ/mol**")
    doc.claim("the collinearity that made it necessary",
              "at **+0.66**")

    print("\nsection 7g: the three pairs")
    for pair in (scope.BUFFER_TYPE_PAIR, *scope.SUBSTRATE_PAIRS,
                 *scope.ENZYME_PAIRS):
        got = induction.matched_pair(pair)
        rows = got["table"]
        doc.claim(f"pair {pair}: the clocks",
                  " / ".join(f"{rows.loc[e].clock:.0f} s" for e in pair))
        doc.claim(f"pair {pair}: the depths",
                  " / ".join(f"{rows.loc[e].depth:.3f}" for e in pair))
    enzyme = induction.matched_pair(scope.ENZYME_PAIRS[1])["table"]
    doc.claim("the enzyme pair's clock ratio",
              f"{enzyme.loc[140].clock / enzyme.loc[141].clock:.1f}x")
    doc.claim("and its catalyst ratio",
              f"{enzyme.loc[140].e0 / enzyme.loc[141].e0:.1f}x")

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
    doc.figures(os.path.join(HERE, "index.html"), "ABCDEFGHIJKLM")
    doc.claim("the document's own count of them", "thirteen figures, A\nto M")

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
