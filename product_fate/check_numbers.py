"""
Verifies every number quoted in product_fate/ANALYSIS.md against the code.

Same contract as the sibling folders', and for the same reason: a figure typed
into a document is a copy, and copies in this project have gone stale three
times.

    python product_fate/check_numbers.py
"""
import io
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import arrhenius
import scope
import slowdown

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")
FAILURES = []


def _normalise(text):
    """Fold typography onto what %-formatting emits. See the sibling modules."""
    return " ".join((text.replace("−", "-")
                         .replace("–", "-")
                         .replace("±", "+/-")
                         .replace("×", "x")
                         .replace("₀", "0")
                         .replace("→", "->")
                         .replace("÷", "/")
                         .replace("χ²", "chi2")
                         .replace("Δ", "D").replace("‡", "")
                         .replace("τ", "tau")
                         .replace("σ", "sigma")
                         .replace("µ", "u").replace("μ", "u")
                         .replace("°", ""))
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
    print("\nsection 6: what the slowdown is")
    whole = scope.frame(tuple(range(1, 152)))
    named = slowdown.substrate_blocks(whole)
    drivers = {name: slowdown.deceleration_drivers(block)
               for name, block in named.items()}
    for label, name in (("the temperature series", "temperature series"),
                        ("4OMe catalysed, phosphate",
                         "4OMe catalysed, phosphate"),
                        ("4OMe, no enzyme at all", "4OMe enzyme-free"),
                        ("BnOH, exps 135-151", "BnOH in scope (135-151)"),
                        ("BnOH catalysed, every buffer",
                         "BnOH catalysed, all buffers")):
        row = drivers[name]
        claim(f"{label}: run length",
              f"{row['span']:+.3f} +/- {row['span_stderr']:.3f}")
        claim(f"{label}: product",
              f"{row['product']:+.3f} +/- {row['product_stderr']:.3f}")
        claim(f"{label}: curves", f"| {row['points']} |")
    controlled = slowdown.deceleration_drivers(
        named["4OMe catalysed, phosphate"], extra=("pH", "temperature"))
    claim("with pH and T in the model",
          f"{controlled['product']:+.3f} +/- "
          f"{controlled['product_stderr']:.3f}")
    claim("pH's own coefficient",
          f"{controlled['pH']:+.3f} +/- {controlled['pH_stderr']:.3f}")
    fixed = {name: slowdown.deceleration_drivers(named[name], fixed=True)
             for name in ("4OMe catalysed, phosphate",
                          "BnOH catalysed, all buffers")}
    for name in fixed:
        claim(f"one dummy per experiment: {name}",
              f"{fixed[name]['product']:+.3f} +/- "
              f"{fixed[name]['product_stderr']:.3f}")
    gap = (fixed["BnOH catalysed, all buffers"]["product"]
           - fixed["4OMe catalysed, phosphate"]["product"])
    error = np.hypot(fixed["BnOH catalysed, all buffers"]["product_stderr"],
                     fixed["4OMe catalysed, phosphate"]["product_stderr"])
    claim("the gap between the substrates", f"{gap:.2f} +/- {error:.2f}")
    claim("and how many sigma it is", f"{gap / error:.1f}sigma")
    claim("how far the two regressors correlate",
          f"{drivers['4OMe catalysed, phosphate']['collinearity']:.2f}")

    table = slowdown.sink_table(
        sorted(named["4OMe catalysed, phosphate"].experiment.unique()))
    clean = table[(table.points > 0)
                  & (table.rate_r2 > slowdown.SINK_CLEAN_R2)]
    counts = clean.prefers.value_counts()
    claim("curves deep enough to choose on", f"**{len(clean)}** catalysed")
    claim("how many favour the rate",
          f"**{int(counts.get('sink', 0))} favour the rate, "
          f"{int(counts.get('inhibition', 0))} favour the reciprocal**")
    claim("and how many tie", f"{int(counts.get('tied', 0))} tie within 1%")
    claim("the median R2 of the rate line", f"**{clean.rate_r2.median():.3f}**")
    claim("the median R2 of the reciprocal",
          f"**{clean.reciprocal_r2.median():.3f}**")

    order = slowdown.plateau_scaling(clean)
    measured = float(arrhenius.substrate_order().corrected.mean())
    claim("the rate's substrate order, from section 3", f"**{measured:+.3f}**")
    claim("the plateau's own order",
          f"is {order['order']:+.3f} +/- {order['stderr']:.3f}")
    claim("over how much substrate", f"{order['lever']:.0f}x range in [S]")
    claim("on how many curves", f"{order['points']} curves")
    picked = np.array([slowdown.selectivity(p, s, e) for p, s, e
                       in zip(clean.plateau, clean.s0, clean.epsilon)])
    picked = picked[np.isfinite(picked)]
    claim("the selectivity",
          f"median {np.median(picked):.0f}, IQR "
          f"{np.percentile(picked, 25):.0f}-{np.percentile(picked, 75):.0f}")
    kelvin = 298.15
    barrier = arrhenius.GAS_CONSTANT * kelvin * np.log(np.median(picked)) / 1000
    claim("the barrier difference it implies", f"about **{barrier:.1f} kJ/mol**")

    budget = slowdown.product_budget(whole).loc["4OMe-BnOH"]
    for experiment in scope.TEMPERATURE_SERIES:
        row = budget.loc[experiment]
        claim(f"exp {experiment}'s product budget",
              f"| {experiment} | {row.temperature:.0f} C | "
              f"{row.product_mM * 1000:.1f} uM | {row.e0 * 1000:.0f} uM | "
              f"{row.of_enzyme:.1%} | {row.late_over_early:.2f} |")
    sharpest = whole[whole.live & (whole.experiment == 21)]
    loss = 1 - float(sharpest.late_over_early.median())
    product = float((sharpest.net / sharpest.epsilon).max())
    enzyme = float(sharpest.e0.iloc[0])
    claim("the sharpest case for the arithmetic",
          f"**exp 21 loses {loss:.0%} of its rate on {product * 1000:.0f} uM "
          f"of product against {enzyme * 1000:.0f} uM of catalyst**")
    claim("what blocking that fraction would take",
          f"{loss * enzyme * 1000:.0f} uM bound")
    for label, name in (("4OMe", "4OMe catalysed, all buffers"),
                        ("BnOH", "BnOH catalysed, all buffers")):
        block = named[name]
        block = block[block.live]
        reach = float(np.max(block.net / block.epsilon))
        claim(f"how much product {label} reaches", f"{reach:.3f} mM")
    check("the deceleration is measured on scope's own statistic, not a copy",
          "late_over_early" in whole.columns)

    print("\nthe overlapping-product-range control")
    window = (0.004, 0.0772)
    for label, name in (("catalysed", "4OMe catalysed, phosphate"),
                        ("enzyme-free", "4OMe enzyme-free")):
        block = named[name]
        share = block.net / block.epsilon
        inside = block[(share >= window[0]) & (share <= window[1])]
        row = slowdown.deceleration_drivers(inside)
        claim(f"{label} inside the shared window",
              f"{row['product']:+.3f} +/- {row['product_stderr']:.3f}")
        claim(f"{label} curve count inside it", f"{row['points']} ")
    claim("the window itself",
          f"{window[0]}-{window[1]:.3f} mM")

    print("\nthe weak bound on the BnOH side")
    bnoh = slowdown.sink_table(scope.PRIMARY_SCOPE)
    good = bnoh[(bnoh.points > 0) & (bnoh.rate_r2 > slowdown.SINK_CLEAN_R2)]
    picked = np.array([slowdown.selectivity(p, s, e) for p, s, e
                       in zip(good.plateau, good.s0, good.epsilon)])
    picked = picked[np.isfinite(picked)]
    claim("the BnOH selectivity", f"**{np.median(picked):.0f}**")
    live = scope.frame(scope.PRIMARY_SCOPE)
    claim("on how few BnOH curves",
          f"{len(picked)} of {int(live.live.sum())} curves that yield a "
          "plateau at all")

    print("\nwhat the sink does to the activation parameters")
    effect = slowdown.sink_effect_on_activation()
    for label, key in (("v_peak, published", "published"),
                       ("v_prod, sink model", "corrected")):
        row = effect[key]
        claim(f"{label}: the row",
              f"| {row['activation_kJ']:.2f} +/- "
              f"{row['activation_stderr']:.2f} | {row['enthalpy_kJ']:.2f} | "
              f"{row['entropy_J']:+.1f} +/- {row['entropy_stderr']:.1f} | "
              f"{row['gibbs_kJ']:.2f} | {row['rms']:.3f} |")
    claim("how far above the published estimator the production rate sits",
          f"**{effect['lift'] - 1:.1%} above**")
    claim("and its range over the six temperatures",
          f"between {effect['lift_low'] - 1:+.1%} and "
          f"{effect['lift_high'] - 1:+.1%}")
    claim("the shift in the activation energy",
          f"**{effect['activation_shift']:+.2f} +/- "
          f"{effect['activation_shift_stderr']:.2f} kJ/mol**")
    claim("the shift in the entropy",
          f"**{effect['entropy_from_lift']:+.2f} J/mol/K**")
    claim("the sink's own activation energy",
          f"**{effect['sink_activation_kJ']:.1f} +/- "
          f"{effect['sink_stderr_kJ']:.1f} kJ/mol**")
    claim("the temperatures it was fitted to",
          ", ".join(f"{t - 273.15:.0f}" for t in
                    effect["sink_temperatures"][:2]) + " and "
          f"{effect['sink_temperatures'][-1] - 273.15:.0f} C")
    check("the corrected route is the noisier one",
          effect["corrected"]["rms"] > effect["published"]["rms"])

    print("\nthe figures: lettered once each, in order")
    # The sibling folder had J and K twice for as long as it had a section 3a,
    # because a letter is not a number and no check looked at one.
    import re
    page = io.open(os.path.join(HERE, "index.html"), encoding="utf-8").read()
    letters = re.findall(r">([A-Z]) \u00b7 ", page)
    expected = [chr(ord("A") + index) for index in range(len(letters))]
    check("every figure letter is used exactly once, A onwards",
          letters == expected,
          f"{''.join(letters)} against {''.join(expected)}")
    claim("the document's own count of them",
          "six figures, A to " + letters[-1] if len(letters) == 6
          else f"{len(letters)} figures")

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
