"""
Verifies every number quoted in buffer/ANALYSIS.md against the code.

Same contract as the sibling folders', and for the same reason: a figure typed
into a document is a copy, and copies in this project have gone stale.

    python buffer/check_numbers.py
"""
import io
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import buffer_role
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
                         .replace("Δ", "D").replace("‡", "")
                         .replace("τ", "tau")
                         .replace("σ", "sigma")
                         .replace("µ", "u").replace("μ", "u")
                         .replace("°", "")
                         .translate(str.maketrans(
                             "⁻⁰¹²³⁴"
                             "⁵⁶⁷⁸⁹",
                             "-0123456789"))
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
    frame = scope.frame(tuple(range(1, 152)))

    print("\nsection 1: what a titration at one pH can ask")
    prediction = buffer_role.species_prediction(7.00, 7.53)
    claim("the pKa", f"is **{prediction['pka']:.2f}**")
    claim("the two base fractions",
          f"{prediction['low_base_fraction']:.3f} at pH 7.00 and "
          f"{prediction['high_base_fraction']:.3f} at pH 7.53")
    claim("what general base predicts",
          f"| **{prediction['general_base']:.2f}** |")
    claim("what general acid predicts",
          f"| **{prediction['general_acid']:.2f}** |")

    print("\nsection 2: the five titrations")
    table = buffer_role.titration_table()
    claim("how many curves", f"**{len(table)} curves at three pH values**"
          if False else f"**{int(table.curves.sum())} curves at three pH values**")
    for row in table.itertuples():
        claim(f"exp {row.experiment}: the row",
              f"| {row.experiment} | {row.pH:.2f} | {row.s0:.3f} | "
              f"{row.buffer_low:g}-{row.buffer_high:g} |")
        bold = "**" if abs(row.order) > 2 * row.order_stderr else ""
        claim(f"exp {row.experiment}: its order",
              f"{bold}{row.order:+.3f} +/- {row.order_stderr:.3f}{bold} | "
              f"{row.order_r2:.3f} |")
    check("five runs, twenty curves",
          len(table) == 5 and int(table.curves.sum()) == 20,
          f"{len(table)}, {int(table.curves.sum())}")
    low = table.set_index("experiment")
    claim("the coefficient at pH 7.00",
          f"is 1.86 x 10-7 at pH 7.00")
    check("and it is", abs(low.loc[32].coefficient - 1.86e-7) < 5e-10,
          f"{low.loc[32].coefficient:.3g}")
    claim("the coefficients at pH 7.53",
          f"3.19-3.23 x 10-7 at pH 7.53")
    check("and they are",
          abs(low.loc[36].coefficient - 3.19e-7) < 5e-10
          and abs(low.loc[37].coefficient - 3.23e-7) < 5e-10,
          f"{low.loc[36].coefficient:.3g}, {low.loc[37].coefficient:.3g}")
    claim("the intercepts", "from an intercept of 3.29 x 10-5 to "
                            "2.73-4.91 x 10-4")
    check("and they are",
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
        claim(f"{label}: the ratio",
              f"| {bold}{got['ratio']:+.2f} +/- {got['ratio_stderr']:.2f}{bold} |")
    measured = buffer_role.catalytic_coefficient(drop=(35,))
    verdict = buffer_role.separable(measured, prediction)
    for name, label in (("general_base", "general base"),
                        ("general_acid", "general acid"),
                        ("spectator", "spectator")):
        row = verdict[name]
        claim(f"{label}: its distance",
              f"| {label} | {row['predicted']:.2f} | {row['sigma']:.1f}sigma | "
              f"survives |")
    check("nothing is excluded", len(verdict["survivors"]) == 3,
          f"{verdict['survivors']}")
    claim("exp 35's own R2", f"has R2 {low.loc[35].order_r2:.3f}")
    check("which is why it is dropped", low.loc[35].order_r2 < 0.2)

    print("\nsection 4: the buffer as a confound")
    blocks = induction.induction_blocks(frame)
    for label, name in (("4OMe catalysed", "4OMe catalysed"),
                        ("4OMe enzyme-free", "4OMe enzyme-free"),
                        ("the temperature series", "temperature series"),
                        ("BnOH in scope", "BnOH in scope (135-151)")):
        got = induction.composition_collinearity(blocks[name])
        bold = "**" if got["median"] == 0.0 or got["median"] < -0.95 else ""
        # A zero collinearity is written without a sign in the document, which
        # is the right typography and needs saying here rather than there.
        shown = f"{got['median']:.2f}" if got["median"] == 0.0 \
            else f"{got['median']:+.2f}"
        claim(f"{label}: runs and collinearity",
              f"| {got['runs']} | {bold}{shown}{bold} |")
        if np.isfinite(got.get("slope", np.nan)):
            claim(f"{label}: the ladder slope", f"{got['slope']:.3f} |")
    check("in-scope every run holds [buf] constant",
          induction.composition_collinearity(
              blocks["BnOH in scope (135-151)"])["constant_buffer"] == 17)

    print("\nsection 5: identity against pH")
    identity = buffer_role.identity_overlap(frame)
    for row in identity.itertuples():
        span = (f"{row.pH_low:.2f}" if row.pH_low == row.pH_high
                else f"{row.pH_low:.2f}-{row.pH_high:.2f}")
        claim(f"{row.substrate} {row.channel} {row.buffer}",
              f"| {row.curves} | {span} |")
    widest = buffer_role.overlap_width(identity)
    claim("the widest overlap", f"**{widest['width']:.2f} units**")
    claim("the phosphate block's peroxide",
          f"sits at exactly\n{widest['peroxide'][0][0]:g} mM")
    claim("the pyrophosphate block's",
          f"at exactly {widest['peroxide'][1][0]:g} and "
          f"{widest['peroxide'][1][1]:g} mM")
    check("and the two cells share no peroxide value",
          not widest["shares_peroxide"], f"{widest['peroxide']}")

    print("\nthe induction's buffer order, quoted back from induction/")
    order = induction.buffer_order(
        induction.buffer_lever(induction.induction_table(
            induction.WHOLE_ARCHIVE)))
    claim("it", f"{order['slope']:+.3f} +/- {order['stderr']:.3f}")

    print("\nthe figures the document promises")
    # Match the title STRING anywhere in the file: the titles are passed to
    # `render` and are not always at the start of a line, and a start-of-line
    # test also catches captions that happen to hold a middle dot -- figure B's
    # carries "<em>b</em>&#183;[buf]".
    letters = sorted(set(re.findall(
        r'"([A-Z]) \u00b7 ',
        io.open(os.path.join(HERE, "build_figures.py"),
                encoding="utf-8").read())))
    check("four figures, A to D", letters == list("ABCD"),
          f"{''.join(letters)}")
    claim("the document's own count of them", "four figures, A to\nD")

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
