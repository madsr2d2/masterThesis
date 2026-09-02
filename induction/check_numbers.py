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
import induction
import scope

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")
FAILURES = []


def _normalise(text):
    """Fold typography onto what %-formatting emits. See the sibling modules."""
    return " ".join((text.replace("−", "-")
                         .replace("–", "-")
                         .replace("±", "+/-")
                         .replace("×", "x")
                         .replace("₀", "0").replace("₂", "2")
                         .replace("→", "->")
                         .replace("÷", "/")
                         .replace("χ²", "chi2")
                         .replace("Δ", "D").replace("‡", "")
                         .replace("τ", "tau")
                         .replace("σ", "sigma")
                         .replace("µ", "u").replace("μ", "u")
                         .replace("°", "")
                         .replace("`", "").replace("*", ""))
                    .split())


def claim(label, rendered, present=True):
    text = _normalise(io.open(DOCUMENT, encoding="utf-8").read())
    ok = (_normalise(rendered) in text) == present
    print(f"  {'pass' if ok else 'FAIL'}  {label}: {rendered!r}")
    if not ok:
        FAILURES.append(f"{label} -- {rendered!r}")


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}{': ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(f"{label} {detail}")


def main():
    table = induction.induction_table(induction.WHOLE_ARCHIVE)
    named = induction.induction_blocks(table)

    print("\nsection 1: the design and the statistic")
    lever = induction.substrate_lever(named["4OMe catalysed"])
    claim("the in-run substrate lever",
          f"a factor of {lever['median_lever']:.1f} in [S]")
    claim("how many runs carry a ladder",
          f"{lever['laddered']} of\nthe {lever['experiments']} experiments")
    check("the block is 147 live curves in 38 experiments",
          int(named["4OMe catalysed"].live.sum()) == 147
          and lever["experiments"] == 38,
          f"{int(named['4OMe catalysed'].live.sum())}, {lever['experiments']}")
    claim("the live count", "147 live\ncurves in 38 experiments")
    claim("the archive is 402 curves", "all 402 curves")
    check("and it is", len(table) == 402, f"{len(table)}")

    windows = induction.landmark_window()
    narrow = windows[windows.window == "300 s"].iloc[0]
    tenth = windows[windows.window.str.startswith("0.10")].iloc[0]
    wide = windows[windows.window == "900 s"].iloc[0]
    claim("what a 300 s window does to the cold end",
          f"{narrow.cold_s:.0f} s** instead of\n**{tenth.cold_s:.0f} s")
    claim("and to the activation energy",
          f"from {tenth.activation_kJ:.1f} +/- {tenth.stderr_kJ:.1f} to\n"
          f"**{narrow.activation_kJ:.1f} +/- {narrow.stderr_kJ:.1f} kJ/mol")
    claim("a 900 s window recovers",
          f"to {wide.activation_kJ:.1f} +/- {wide.stderr_kJ:.1f}")

    drift = induction.schedule_dependence(table)
    claim("the between-run schedule dependence",
          f"an exponent of **{drift['span']:+.3f} +/- {drift['span_stderr']:.3f}")
    claim("and pH once the schedule is in the model",
          f"{drift['pH']:+.3f} +/- {drift['pH_stderr']:.3f}")

    print("\nsection 2: the induction needs the catalyst")
    summary = induction.channel_summary(table)
    for label, key in (("catalysed", "catalysed"),
                       ("enzyme-free", "enzyme_free")):
        row = summary[key]
        claim(f"{label}: curves and depth",
              f"| {row['curves']} | **{row['median_depth']:.3f}** | ")
    claim("the enzyme-free block has no deep curves and none accelerating",
          f"| **{summary['enzyme_free']['deep']}** | "
          f"**{summary['enzyme_free']['accelerates']}** |")
    claim("the largest enzyme-free acceleration z",
          f"is **{summary['enzyme_free']['max_accel_z']:.2f}**")
    claim("against the catalysed block's",
          f"against {summary['catalysed']['max_accel_z']:.1f} in the\ncatalysed one")
    claim("the longest enzyme-free run",
          f"exp 28 ran **{summary['enzyme_free']['longest_s']:.0f} s**")

    contrast = induction.channel_contrast(table)
    for row in contrast.itertuples():
        claim(f"matched cell {row.temperature:.0f} C {row.channel}",
              f"| {row.curves} | ")
        claim(f"matched cell {row.temperature:.0f} C {row.channel}: depth",
              f"{row.depth:.3f}")

    print("\nsection 3: a clock, not a product threshold")
    for rate in ("v_peak", "peak_rate", "vmax"):
        row = induction.induction_drivers(named["4OMe catalysed"], rate=rate)
        claim(f"the driver coefficient against {rate}",
              f"{row['slope']:+.3f} +/- {row['stderr']:.3f}** | {row['points']} |"
              if rate == "v_peak" else
              f"{row['slope']:+.3f} +/- {row['stderr']:.3f} | {row['points']} |")
    series = induction.induction_drivers(named["temperature series"],
                                         rate="v_peak")
    claim("the temperature series cannot do it",
          f"{series['slope']:+.3f} +/- {series['stderr']:.3f} | "
          f"{series['points']} |")
    main_fit = induction.induction_drivers(named["4OMe catalysed"],
                                           rate="v_peak")
    sigma = (main_fit["slope"] - induction.INDUCTION_PRODUCT_SLOPE) \
        / main_fit["stderr"]
    claim("how many sigma product control is excluded by",
          f"excluded at {sigma:.0f}sigma")
    low = induction.induction_drivers(named["4OMe catalysed"],
                                      floor=induction.FLOOR_SWEEP[0])
    high = induction.induction_drivers(named["4OMe catalysed"],
                                       floor=induction.FLOOR_SWEEP[-1])
    claim("the floor sweep",
          f"from {low['slope']:+.3f} +/- {low['stderr']:.3f} to "
          f"{high['slope']:+.3f} +/- {high['stderr']:.3f}")

    for label, name in (("4OMe catalysed", "4OMe catalysed"),
                        ("the temperature series", "temperature series")):
        ratio = induction.order_ratio(named[name])
        claim(f"{label}: order of the induction time",
              f"{ratio['induction_order']:+.3f} +/- "
              f"{ratio['induction_stderr']:.3f}")
        claim(f"{label}: order of the rate",
              f"{ratio['rate_order']:+.3f} +/- {ratio['rate_stderr']:.3f}")
        claim(f"{label}: the ratio",
              f"{ratio['ratio']:+.2f} +/- {ratio['ratio_stderr']:.2f}")
        depth = scope.orders("depth", frame=named[name],
                             floor=induction.DEPTH_FLOOR)
        claim(f"{label}: order of the amplitude",
              f"{depth['order_s0']:+.3f} +/- {depth['stderr_s0']:.3f}")

    print("\nsection 4: not the schedule, not physical")
    control = induction.schedule_control(table)
    claim("the two runs of equal length",
          f"both recorded for **{control[19]['duration_s']:.0f} s**")
    claim("and their two time constants",
          f"**{control[19]['tau']:.0f} s and {control[14]['tau']:.0f} s**")
    claim("the ratio", f"a factor of {control['tau_ratio']:.1f}")

    gap = induction.activation_contrast()
    claim("the induction's activation energy",
          f"**{gap['induction']['activation_kJ']:.1f} +/- "
          f"{gap['induction']['activation_stderr']:.1f} kJ/mol**")
    claim("on 16 curves at four temperatures",
          f"{gap['induction']['n']} curves at four temperatures")

    for label, name in (("4OMe catalysed", "4OMe catalysed"),
                        ("BnOH", "BnOH in scope (135-151)")):
        row = induction.signal_control(named[name])
        claim(f"{label}: the signal-to-noise control",
              f"{row['signal_slope']:+.3f} +/- {row['signal_stderr']:.3f}** | "
              f"{row['signal_points']} |" if name == "4OMe catalysed" else
              f"{row['signal_slope']:+.3f} +/- {row['signal_stderr']:.3f} | "
              f"{row['signal_points']} |")
    lever = induction.peroxide_lever(table)
    claim("the peroxide block's own signal control",
          f"{lever['signal_slope']:+.3f} +/- {lever['signal_stderr']:.3f} | "
          f"{lever['signal_points']} |")
    claim("the peroxide levels",
          f"{lever['levels'][0]:.3f} against {lever['levels'][1]:.3f} mM")
    claim("the peroxide block's size",
          f"{lever['curves']} curves of which {lever['live']} are live")
    claim("the induction's peroxide order",
          f"**{lever['induction_order']:+.3f} +/- "
          f"{lever['induction_stderr']:.3f}**")
    claim("the rate's peroxide order there",
          f"{lever['rate_order']:+.3f} +/- {lever['rate_stderr']:.3f}")

    scoped = scope.orders("t_ind", frame=named["BnOH in scope (135-151)"],
                          floor=induction.INDUCTION_FLOOR)
    claim("BnOH's induction order in peroxide",
          f"{scoped['order_h2o2']:+.3f} +/- {scoped['stderr_h2o2']:.3f}")

    print("\nsection 5: the activation parameters")
    for label, key in (("the induction", "induction"),
                       ("the turnover", "turnover")):
        row = gap[key]
        claim(f"{label}: Ea and enthalpy",
              f"| {row['activation_kJ']:.1f} +/- {row['activation_stderr']:.1f} "
              f"| {row['enthalpy_kJ']:.1f} +/- {row['enthalpy_stderr']:.1f} |")
        claim(f"{label}: entropy",
              f"{row['entropy_J']:+.1f} +/- {row['entropy_stderr']:.1f}")
        claim(f"{label}: free energy",
              f"**{row['gibbs_kJ']:.2f} +/- {row['gibbs_stderr']:.2f}**")
    claim("the free-energy gap",
          f"**{gap['gibbs_gap_kJ']:.2f} +/- {gap['gibbs_gap_stderr']:.2f} kJ/mol")
    claim("what it is worth as a rate constant",
          f"**{gap['rate_ratio']:.0f}x faster**")
    claim("the entropy gap", f"{gap['entropy_gap_J']:+.1f} J/mol/K")
    claim("its error", f"**+/- {gap['entropy_gap_stderr']:.1f}**")
    claim("and the enthalpy gap's",
          f"**+/- {gap['enthalpy_gap_stderr']:.1f}**")
    check("the free-energy gap is 19 or more standard errors wide",
          abs(gap["gibbs_gap_kJ"]) / gap["gibbs_gap_stderr"] > 19,
          f"{abs(gap['gibbs_gap_kJ']) / gap['gibbs_gap_stderr']:.1f}")
    check("and each of its two parts is about one",
          abs(gap["entropy_gap_J"]) / gap["entropy_gap_stderr"] < 1.5
          and abs(gap["enthalpy_gap_kJ"]) / gap["enthalpy_gap_stderr"] < 1.5)

    print("\nthe product_fate numbers this document quotes back")
    import slowdown
    whole = scope.frame(induction.WHOLE_ARCHIVE)
    fall = slowdown.deceleration_drivers(
        slowdown.substrate_blocks(whole)["4OMe catalysed, phosphate"])
    claim("the fall's product coefficient",
          f"({fall['product']:+.3f} +/- {fall['product_stderr']:.3f})")
    claim("the fall's clock coefficient",
          f"({fall['span']:+.3f} +/- {fall['span_stderr']:.3f})")

    print("\nthe figures the document promises")
    letters = sorted(set(
        line.split("·")[0].strip().strip('"').split()[-1]
        for line in io.open(os.path.join(HERE, "build_figures.py"),
                            encoding="utf-8").read().splitlines()
        if "·" in line and '"' in line and line.strip().startswith('"')))
    expected = list("ABCDEFG")
    check("seven figures, A to G", letters == expected,
          f"{''.join(letters)} against {''.join(expected)}")
    claim("the document's own count of them", "seven figures, A\nto G")

    print("\nthe figures: no data point drawn outside its own frame")
    from svgplot import clipped_marks
    path = os.path.join(HERE, "index.html")
    if os.path.exists(path):
        lost = clipped_marks(io.open(path, encoding="utf-8").read())
        ok = not lost
        print(f"  {'pass' if ok else 'FAIL'}  index.html: {len(lost)} clipped")
        if not ok:
            for title, x, y, _ in lost[:4]:
                print(f"        {title!r} at ({x:.4g},{y:.4g})")
            FAILURES.append(f"index.html draws {len(lost)} point(s) outside "
                            f"the plot frame; widen the axis limits")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} MISMATCH(ES):")
        for item in FAILURES:
            print(f"  - {item}")
        return 1
    print("ANALYSIS.md agrees with the code")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
