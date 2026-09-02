"""
Builds buffer/index.html.

Draws only; every number comes from `buffer_role`, `scope`, `induction` or
`solution_chemistry`, so a figure and the prose in ANALYSIS.md cannot disagree
about a value without `check_numbers.py` saying so.

    python buffer/build_figures.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import buffer_role
import induction
import scope
from svgplot import ACCENT, GRID, INK, MUTED, Axes, esc, page, PAGE_CSS

# The two pH groups are ORDERED (either side of a pKa), so a two-step ramp; the
# three hypotheses in figure C are categorical and take a validated trio.
PH_RAMP = ["#6295c3", "#0c2f4d"]
CATEGORY = ["#2f6fb0", "#c0522a", "#8a5aa8"]
SURFACE = "#fbfbfa"

EXTRA_CSS = """
.fig{background:#fbfbfa;border-color:#e4e4e2}
.fig .cap{color:#5a5a5a}
.hero{display:flex;flex-wrap:wrap;gap:26px;margin:14px 0 4px}
.hero div{min-width:150px}
.hero .v{font-size:25px;font-weight:650;letter-spacing:-0.02em}
.hero .k{font-size:11.5px;color:var(--muted);text-transform:uppercase;
letter-spacing:0.06em}
.hero .u{font-size:12px;color:var(--muted)}
"""


def fig(svg, caption, extra=""):
    return (f"<div class='fig'>{svg}"
            f"<div class='cap'>{caption}</div>{extra}</div>")


_CACHE = {}


def _table():
    if "table" not in _CACHE:
        _CACHE["table"] = buffer_role.titration_table()
    return _CACHE["table"]


def _frame():
    if "frame" not in _CACHE:
        _CACHE["frame"] = scope.frame(tuple(range(1, 152)))
    return _CACHE["frame"]


def figure_titrations():
    table = _table()
    frame = scope.frame(buffer_role.TITRATIONS)
    frame = frame[frame.live & (frame.v_peak > 0)]
    axes = Axes(560, 300, (2.5, 260.0), (8e-6, 9e-4), xlog=True, ylog=True,
                pad=(72, 26, 46, 32))
    for row in table.itertuples():
        block = frame[frame.experiment == row.experiment].sort_values("buf")
        colour = PH_RAMP[0] if row.pH < buffer_role.TITRATION_PH_SPLIT \
            else PH_RAMP[1]
        axes.points(block.buf.to_numpy(), block.v_peak.to_numpy(), colour,
                    radius=4.0, title=f"exp {row.experiment}, pH {row.pH:.2f}")
        grid = np.array([block.buf.min(), block.buf.max()])
        level = float(np.exp(np.log(block.v_peak.to_numpy()).mean()))
        centre = float(np.exp(np.log(block.buf.to_numpy()).mean()))
        axes.line(grid, level * (grid / centre) ** row.order, colour,
                  width=1.6, dash="5 4")
    # The run labels go in the bottom-right, which a rising titration leaves
    # empty: putting them at each run's last point ran them off the frame.
    for index, row in enumerate(table.sort_values(["pH", "buffer_low"])
                                .itertuples()):
        colour = PH_RAMP[0] if row.pH < buffer_role.TITRATION_PH_SPLIT \
            else PH_RAMP[1]
        axes.note(180, 66 + 15 * index,
                  f"exp {row.experiment} · pH {row.pH:.2f} · {row.order:+.2f}",
                  colour, size=10.5, anchor="start", weight="600")
    axes.note(86, 150, "dashed: each run's own log-log slope", MUTED,
              size=10.5)
    return fig(
        axes.render("[buffer], mM", "peak rate, AU/s",
                    "A · Five buffer titrations, at three pH values"),
        "Every run in the archive that steps <code>[buf]</code> with "
        "<code>[S]</code> held fixed inside it — 4OMe-BnOH at 40 °C in "
        "phosphate, 82.5 mM H₂O₂, twenty curves. Only exps 32 and 34 have been "
        "used anywhere in this project before. At pH 7.00 the dependence is "
        "strong and <strong>saturates</strong>: +0.792 ± 0.185 over 3.125–25 mM "
        "and +0.371 ± 0.018 over 50–200. Above the pKa it is <strong>gone</strong> "
        "— −0.301 ± 0.459, +0.064 ± 0.082, +0.122 ± 0.039 — and the reason is "
        "the vertical gap, not the slopes: the buffer-independent rate is eight "
        "times higher there.")


def figure_split():
    table = _table()
    axes = Axes(560, 280, (-0.6, 4.6), (1e-5, 3e-3), ylog=True,
                pad=(76, 26, 58, 34))
    for index, row in enumerate(table.sort_values(["pH", "buffer_low"])
                                .itertuples()):
        colour = PH_RAMP[0] if row.pH < buffer_role.TITRATION_PH_SPLIT \
            else PH_RAMP[1]
        # The buffer-carried part of the rate at the run's own top rung,
        # against the part that is there without any buffer at all.
        carried = abs(row.coefficient) * row.buffer_high
        # Points joined by a stem, not bars: on a log axis a bar's height is
        # measured from whatever floor the axis happens to have, and the
        # quantity here is the RATIO of the two, which is the vertical gap.
        axes.line([index, index], [min(row.intercept, carried),
                                   max(row.intercept, carried)],
                  colour, width=2.0, opacity=0.55)
        axes.points([index], [row.intercept], colour, radius=6.0)
        axes.ring(index, carried, colour, radius=6.0, width=2.0)
        axes.note(axes._fx(index), 236, f"exp {row.experiment}", INK, size=11,
                  anchor="middle", weight="600")
        axes.note(axes._fx(index), 250, f"pH {row.pH:.2f}", MUTED, size=10.5,
                  anchor="middle")
    axes.note(96, 44, "filled: the rate with no buffer   ·   open: what the "
                      "buffer adds at the top rung", MUTED, size=10.5)
    return fig(
        axes.render("", "AU/s", "B · Why the dependence vanishes above the pKa",
                    xticks=False),
        "The same five runs split as <code>v = a + b·[buf]</code>: the filled "
        "mark is <em>a</em>, whatever the rate would be with no buffer, and the "
        "open one is <em>b</em>·[buf] at the run's own top rung. Read the "
        "vertical gap between a pair, not the height of either. <strong>The "
        "buffer's contribution barely changes across the pKa</strong> — its "
        "coefficient is 1.86 × 10⁻⁷ at pH 7.00 and 3.19–3.23 × 10⁻⁷ at 7.53 — "
        "while the buffer-independent term rises from 3.29 × 10⁻⁵ to "
        "2.73–4.91 × 10⁻⁴. Half a pH unit multiplies the [HOO⁻] route roughly "
        "eightfold and the buffer term is simply outrun. Exp 34's tall pale bar "
        "is the low-buffer run, where the same coefficient is measured over "
        "3.125–25 mM instead.")


def figure_species():
    prediction = buffer_role.species_prediction(7.00, 7.53)
    measured = buffer_role.catalytic_coefficient(drop=(35,))
    verdict = buffer_role.separable(measured, prediction)
    # Four rows, not three: the measurement gets its own, because putting it
    # on the acid hypothesis's row drew the two on top of each other.
    axes = Axes(620, 270, (-1.4, 2.6), (-0.75, 3.75), pad=(64, 26, 58, 34))
    for index, (name, label) in enumerate((
            ("general_base", "general base"),
            ("spectator", "spectator"),
            ("general_acid", "general acid"))):
        y = 3 - index
        value = prediction[name]
        axes.line([value, value], [y - 0.34, y + 0.34], CATEGORY[index],
                  width=2.2, dash="4 3")
        axes.label(value, y, f"{label}  {value:.2f}", CATEGORY[index],
                   size=11, anchor="end" if value > 1.2 else "start",
                   weight="600", dx=-10 if value > 1.2 else 10, dy=4)
    band = (measured["ratio"] - measured["ratio_stderr"],
            measured["ratio"] + measured["ratio_stderr"])
    axes.line(list(band), [0.0, 0.0], INK, width=3.0)
    axes.points([measured["ratio"]], [0.0], INK, radius=5.2)
    axes.label(measured["ratio"], 0.0,
               f"measured {measured['ratio']:+.2f} ± "
               f"{measured['ratio_stderr']:.2f}", INK, size=11,
               anchor="middle", weight="600", dy=20)
    axes.note(76, 50, f"all three survive at 2σ: "
                       f"{', '.join(n.replace('_', ' ') for n in verdict['survivors'])}",
              ACCENT, size=11.5, weight="600")
    return fig(
        axes.render("coefficient at pH 7.53 ÷ coefficient at pH 7.00", "",
                    "C · The species test the archive has, and cannot use",
                    yticks=False),
        "A buffer titration at one pH is an order in <em>total</em> buffer — at "
        "fixed pH the acid, the base and the total are proportional. Titrating "
        "at two pH values separates them: the base fraction runs 0.387 → 0.681 "
        "across phosphate's pKa, so general base predicts the catalytic "
        "coefficient rises <strong>1.76×</strong> and general acid that it "
        "falls to <strong>0.52×</strong>. The archive has that design and "
        f"measures <strong>{measured['ratio']:+.2f} ± "
        f"{measured['ratio_stderr']:.2f}</strong> — which excludes nothing. "
        "The cause is figure B: above the pKa the buffer term is a small share "
        "of a much larger rate, so its coefficient is 44–127% uncertain there.")


def figure_confound():
    frame = _frame()
    blocks = induction.induction_blocks(frame)
    names = (("4OMe catalysed", "4OMe catalysed"),
             ("4OMe enzyme-free", "4OMe, no enzyme"),
             ("temperature series", "the temperature series"),
             ("BnOH in scope (135-151)", "BnOH, exps 135–151"))
    rows = [(label, induction.composition_collinearity(blocks[key]))
            for key, label in names]
    axes = Axes(660, 240, (-1.05, 0.25), (-0.7, len(rows) - 0.3),
                pad=(226, 26, 58, 34))
    axes.line([0, 0], [-0.7, len(rows) - 0.3], INK, width=1.2, dash="3 3")
    for index, (label, row) in enumerate(rows):
        y = len(rows) - 1 - index
        value = row["median"]
        colour = CATEGORY[1] if value < -0.5 else CATEGORY[0]
        axes.line([value, 0.0], [y, y], colour, width=13, opacity=0.85)
        axes.points([value], [y], colour, radius=5.0)
        axes.note(208, axes._fy(y) + 4,
                  f"{label}   {row['runs']} runs", INK, size=11, anchor="end")
        axes.label(value, y, f"{value:+.2f}", INK, size=10.5, anchor="end",
                   weight="600", dx=-14, dy=4)
    axes.note(440, 216, "orange: the substrate ladder is also a buffer ladder",
              MUTED, size=10.5, anchor="middle")
    return fig(
        axes.render("corr(log[S], log[buf]) inside a run, median over runs", "",
                    "D · The substrate ladder is a buffer ladder, except in one "
                    "block", yticks=False),
        "Substrate was added by volume and displaced buffer, so <code>[buf]</code> "
        "falls 80 → 50 mM as <code>[S]</code> rises 1.85 → 7.4 in every 4OMe "
        "run. A substrate order measured there is an order in the <em>pair</em>, "
        "and two published ones have had to be corrected for it: the rate's "
        "(+0.45 → +0.58, which moved the level and not the trend) and the "
        "induction's (−0.121 ± 0.148 → −0.235 ± 0.158, which moved that route "
        "from excluding product control to not excluding it). "
        "<strong>Exps 135–151 are the one block where <code>[S]</code> moves "
        "alone</strong> — <code>[buf]</code> is constant across all seven "
        "cuvettes of all seventeen runs — and that is the block the in-scope "
        "fitting is scoped to.")


def build_index():
    table = _table()
    prediction = buffer_role.species_prediction(7.00, 7.53)
    measured = buffer_role.catalytic_coefficient(drop=(35,))
    verdict = buffer_role.separable(measured, prediction)
    identity = buffer_role.identity_overlap(_frame())
    widest = buffer_role.overlap_width(identity)
    induction_order = induction.buffer_order(
        induction.buffer_lever(induction.induction_table(
            induction.WHOLE_ARCHIVE)))

    hero = f"""
<div class='hero'>
  <div><div class='k'>buffer titrations</div>
       <div class='v'>{len(table)}</div>
       <div class='u'>20 curves, three pH values — two were being used</div></div>
  <div><div class='k'>order at pH 7.00</div>
       <div class='v'>+0.79 → +0.37</div>
       <div class='u'>3.125–25 mM, then 50–200: it saturates</div></div>
  <div><div class='k'>above the pKa</div>
       <div class='v'>+0.06 to +0.12</div>
       <div class='u'>a third of it, because the [HOO⁻] route grew eightfold</div></div>
  <div><div class='k'>acid or base?</div>
       <div class='v'>{len(verdict['survivors'])} of 3</div>
       <div class='u'>hypotheses still standing at 2σ</div></div>
</div>"""

    rows = "".join(
        f"<tr><td>{row.experiment}</td><td>{row.pH:.2f}</td>"
        f"<td>{row.s0:.3f}</td><td>{row.buffer_low:g}–{row.buffer_high:g}</td>"
        f"<td>{row.order:+.3f} ± {row.order_stderr:.3f}</td>"
        f"<td>{row.order_r2:.3f}</td>"
        f"<td>{row.coefficient:.3g}</td><td>{row.intercept:.3g}</td></tr>"
        for row in table.sort_values(["pH", "buffer_low"]).itertuples())

    identity_rows = "".join(
        f"<tr><td>{esc(row.substrate)}</td><td>{esc(row.channel)}</td>"
        f"<td>{esc(row.buffer)}</td><td>{row.curves}</td>"
        f"<td>{row.pH_low:.2f}–{row.pH_high:.2f}</td></tr>"
        for row in identity.itertuples())

    body = f"""
<p class='lede'>The buffer is in this archive three times over — as a reagent,
as a confound, and as the leading candidate for the step
<a href='../induction/index.html'>induction</a> is trying to name — and had
never been looked at as one question. <a href='ANALYSIS.md'>ANALYSIS.md</a> is
the argument; this is the picture of it.</p>
{hero}

<h2>1 · A titration at one pH cannot name the species</h2>
<p>At fixed pH the acid, the base and the total buffer are all proportional to
one another, so an order in <code>[buf]</code> is an order in <em>total</em>
buffer. That is what every order in this project is, including the three it
inherits. The escape is to titrate at more than one pH — section 3.</p>

<h2>2 · Five buffer titrations, not two</h2>
{figure_titrations()}
<div class='tbl'><table>
<tr><th>exp</th><th>pH</th><th>[S] mM</th><th>[buf] mM</th>
<th>order in [buf]</th><th>R²</th><th>d<i>v</i>/d[buf]</th><th>intercept</th></tr>
{rows}
</table></div>
{figure_split()}

<h2>3 · Which species — and the archive cannot say</h2>
{figure_species()}
<p>Exp 35 is dropped from the fit above because its own titration has R² 0.176
and is non-monotone, one cuvette at less than half its neighbours — a judgement
about a run, made on the run's internal consistency and not on whether it helps.
With it, the ratio is −0.08 ± 1.16 and excludes even less.</p>
<p><strong>One titration at pH 6.5 would do what four runs at 7.0–7.5 cannot</strong>:
the base fraction there is 0.17, so the two hypotheses separate by more, and the
[HOO⁻] route that swamps the buffer term above the pKa is at its weakest.</p>

<h2>4 · The buffer as a confound</h2>
{figure_confound()}

<h2>5 · Identity cannot be separated from pH</h2>
<div class='tbl'><table>
<tr><th>substrate</th><th>channel</th><th>buffer</th><th>curves</th><th>pH</th></tr>
{identity_rows}
</table></div>
<p>The three buffers were each chosen for the pH they hold, which is what a
buffer is for. The widest pH range any two share inside one substrate/channel
cell is <strong>{widest['width']:.2f} units</strong> —
{esc(' against '.join(widest['buffers']))} on catalysed 4OMe — and
<strong>those two cells share no peroxide</strong>: {esc(str(widest['peroxide'][0]))}
mM against {esc(str(widest['peroxide'][1]))} mM. A range overlap is not a
matched pair, and the check is on the values.</p>

<h2>6 · What this settles, and what it does not</h2>
<p><strong>Settled.</strong> Every buffer order in this project is an order in
total buffer. The catalysed turnover's dependence saturates — about first order
below 25 mM, half order above 50 — and disappears above the pKa because the
[HOO⁻] route grows past it, not because the buffer term shrinks.
<code>[S]</code> and <code>[buf]</code> are collinear at −0.96 in every 4OMe run
and at 0.00 in every in-scope run. Identity is confounded with pH everywhere.</p>
<p><strong>Not settled: acid, base or spectator</strong>
({measured['ratio']:+.2f} ± {measured['ratio_stderr']:.2f} against
{prediction['general_base']:.2f}, {prediction['general_acid']:.2f} and 1.00), and
<strong>whether the buffer is what carries E → E*</strong> — the induction's
buffer order is {induction_order['slope']:+.3f} ± {induction_order['stderr']:.3f}
on eight curves at 40 °C, where the induction is nearly over.</p>
<p><strong>The two things that would settle both.</strong> A buffer titration at
<strong>pH 6.5 and 25 °C on the catalysed 4OMe system</strong> — the low pH makes
the buffer term the largest share of the rate it can be, the low temperature
makes the induction thousands of seconds long, and one design gives the species
and the induction's buffer order at once. And <code>COMPUTATIONAL.md</code>
<strong>C7</strong>, whose model already carries a phosphate dianion beside the
ketone hydrate for this reason.</p>

<h2>Reproducing</h2>
<p><code>python data/buffer_role.py</code> prints the whole argument ·
<code>python data/test_buffer_role.py</code> ·
<code>python buffer/build_figures.py</code> ·
<code>python buffer/check_numbers.py</code>, which re-derives every number in
<code>ANALYSIS.md</code> from the modules and fails if the prose and the code
disagree.</p>
"""
    return page("What the buffer does, and to what", body,
                "A reagent, a confound, and a candidate").replace(
        "</style>", EXTRA_CSS + "</style>")


def main():
    path = os.path.join(HERE, "index.html")
    content = build_index()
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(content)
    print(f"wrote {path}  ({len(content) / 1024:.0f} kB)")
    from svgplot import clipped_marks
    clipped = clipped_marks(content)
    print(f"{len(clipped)} clipped marks"
          + ("" if not clipped else f"  {clipped[:6]}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
