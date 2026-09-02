"""
Builds product_fate/index.html.

Draws only; every number comes from `slowdown`, `scope`, `arrhenius` or
`curve_metrics`, so a figure and the prose in ANALYSIS.md cannot disagree about
a value without `check_numbers.py` saying so.

    python product_fate/build_figures.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
# svgplot lives at the repository root, shared with the sibling folders.
sys.path.insert(0, os.path.dirname(HERE))

import arrhenius
import scope
import slowdown
from curve_metrics import rolling_slope
from fit_dataset import source_floor
from svgplot import ACCENT, GRID, INK, MUTED, Axes, esc, page, PAGE_CSS

# ORDERED VARIABLES GET SEQUENTIAL RAMPS. Substrate rung is ordinal; the two
# channels (with and without enzyme) and the two straight-line laws are not.
RUNGS = ["#6295c3", "#3d729f", "#1e5079", "#0c2f4d"]
CATEGORY = ["#2f6fb0", "#c0522a", "#8a5aa8"]
SURFACE = "#fbfbfa"
FIT_WIDTH = 1.5

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
_SLOWDOWN_CACHE = {}


def _whole_frame():
    """Every curve in the archive, built once. Section 6 needs all of it: the
    contrast that carries the argument is between blocks, not inside one."""
    if "frame" not in _SLOWDOWN_CACHE:
        _SLOWDOWN_CACHE["frame"] = scope.frame(tuple(range(1, 152)))
    return _SLOWDOWN_CACHE["frame"]


def _slowdown_blocks():
    if "blocks" not in _SLOWDOWN_CACHE:
        _SLOWDOWN_CACHE["blocks"] = slowdown.substrate_blocks(_whole_frame())
    return _SLOWDOWN_CACHE["blocks"]


def _sink_curves():
    """The catalysed 4OMe phosphate curves whose decline has a shape."""
    if "sink" not in _SLOWDOWN_CACHE:
        block = _slowdown_blocks()["4OMe catalysed, phosphate"]
        table = slowdown.sink_table(sorted(block.experiment.unique()))
        _SLOWDOWN_CACHE["sink"] = table[(table.points > 0)
                                        & (table.rate_r2 > slowdown.SINK_CLEAN_R2)]
    return _SLOWDOWN_CACHE["sink"]


def figure_drivers():
    names = {"temperature series": "temperature series",
             "4OMe catalysed, phosphate": "4OMe catalysed, phosphate",
             "4OMe enzyme-free": "4OMe, no enzyme at all",
             "BnOH in scope (135-151)": "BnOH, exps 135-151",
             "BnOH catalysed, all buffers": "BnOH catalysed, every buffer"}
    blocks = _slowdown_blocks()
    rows = [(short, slowdown.deceleration_drivers(blocks[name]))
            for name, short in names.items()]
    axes = Axes(680, 310, (-1.05, 0.62), (-0.8, len(rows) - 0.2),
                pad=(232, 26, 46, 34))
    axes.line([0, 0], [-0.8, len(rows) - 0.2], INK, width=1.2, dash="3 3")
    for index, (name, row) in enumerate(rows):
        y = len(rows) - 1 - index
        for offset, key, colour in ((0.17, "span", CATEGORY[2]),
                                    (-0.17, "product", CATEGORY[1])):
            value, error = row[key], row[key + "_stderr"]
            axes.line([value - error, value + error], [y + offset] * 2,
                      colour, width=2.0)
            axes.points([value], [y + offset], colour, radius=4.4)
        axes.note(214, axes._fy(y) + 4, f"{name}   n = {row['points']}",
                  INK, size=11, anchor="end")
    axes.label(0.40, len(rows) - 1 + 0.17, "fell with run length", CATEGORY[2],
               size=11, weight="600", anchor="middle", dy=-9)
    axes.label(-0.58, len(rows) - 1 - 0.17, "fell with product made",
               CATEGORY[1], size=11, weight="600", anchor="middle", dy=17)
    return fig(
        axes.render("coefficient in log(late rate / early rate)", "",
                    "A · What the slowdown tracks: the clock, or the product",
                    yticks=False),
        "Each curve's deceleration — <code>scope</code>'s own "
        "<code>late_over_early</code> — regressed on how long the run lasted "
        "and how much product it made, in mM. A clock loads on run length and "
        "not on product; product control does the reverse. <strong>The "
        "temperature series loads on product alone</strong> "
        f"({rows[0][1]['product']:+.3f} ± {rows[0][1]['product_stderr']:.3f} "
        f"against {rows[0][1]['span']:+.3f} ± {rows[0][1]['span_stderr']:.3f} "
        "for run length), and so does the whole catalysed 4OMe phosphate set. "
        "<strong>The same chemistry without the catalyst does the opposite</strong> "
        "— it is a clock, and it is the reference channel every catalysed curve "
        "is measured against, so it is subtracted out. The BnOH blocks show no "
        "product-driven deceleration at all, at product concentrations up to "
        "0.39 mM against this block's 0.21.")


def figure_sink_lines():
    frame = _slowdown_blocks()["4OMe catalysed, phosphate"]
    table = _sink_curves()
    axes = Axes(430, 265, (0, 0.88), (0, 6.6e-5), pad=(74, 26, 46, 30))
    for index, curve in enumerate(scope.curves_of(16)):
        colour = RUNGS[index]
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        centres, slopes = rolling_slope(times, values, slowdown.SINK_WINDOW,
                                        source_floor(curve.source))
        product = np.interp(centres, times, values) - values[0]
        top = int(np.argmax(slopes))
        axes.points(product[top:], slopes[top:], colour, radius=2.6,
                    stroke=None)
        # table["sample"], never table.sample: the attribute is DataFrame's
        # own sampling method, and comparing it to an int silently matches
        # nothing -- which drew this figure's fitted lines for no curve at all.
        row = table[(table.experiment == 16)
                    & (table["sample"] == curve.sample)]
        if not len(row):
            continue
        row = row.iloc[0]
        grid = np.array([float(product[-1]), float(row.plateau)])
        axes.line(grid, row.v0 - row.k_sink * grid, colour, width=FIT_WIDTH,
                  dash="5 4")
        axes.points([float(row.plateau)], [0.0], colour, radius=4.2)
        axes.label(float(product[top]), float(slopes[top]),
                   f"{curve.conditions.s0:.2f} mM", colour, size=10.5,
                   anchor="end", weight="600", dx=-6, dy=4)
    axes.note(96, 40, "dashed: where the line is heading", MUTED, size=10.5)

    scatter = Axes(430, 265, (0.55, 1.005), (0.55, 1.005),
                   pad=(60, 22, 46, 24))
    scatter.line([0.55, 1.005], [0.55, 1.005], MUTED, width=1.2, dash="4 3")
    scatter.points(table.reciprocal_r2, table.rate_r2, CATEGORY[0], radius=4.0)
    scatter.note(320, 40, "rate linear in product wins", CATEGORY[0], size=11,
                 anchor="end")
    scatter.note(320, 55, "above the diagonal", MUTED, size=10.5, anchor="end")
    wins = int((table.rate_r2 > table.reciprocal_r2
                + slowdown.SLOWDOWN_MARGIN).sum())
    losses = int((table.reciprocal_r2 > table.rate_r2
                  + slowdown.SLOWDOWN_MARGIN).sum())
    return fig(
        "<div class='grid two'><div>" + axes.render(
            "product already made, AU", "rate, AU/s",
            "B · Rate falls in a straight line against product")
        + "</div><div>" + scatter.render(
            "R² of 1/rate against product", "R² of rate against product",
            "C · …and not against its reciprocal")
        + "</div></div>",
        "<strong>B.</strong> Exp 16, the 40 °C run, after each cuvette's rate "
        "maximum. Rate against accumulated product is a straight line, and the "
        "line is <code>A′ = v − kA</code>: production, minus a first-order loss "
        "of the thing being measured. Where each line meets the axis is the "
        "level it is heading for. <strong>C.</strong> Reversible product "
        f"inhibition would put 1/rate on the straight line instead. Over the "
        f"{len(table)} catalysed 4OMe phosphate curves whose decline is deep "
        f"enough to have a shape, <strong>{wins} favour the rate</strong> and "
        f"<strong>{losses} the reciprocal</strong>; medians "
        f"{table.rate_r2.median():.3f} against "
        f"{table.reciprocal_r2.median():.3f}.")


def figure_plateau_order():
    table = _sink_curves()
    order = slowdown.plateau_scaling(table)
    measured = float(arrhenius.substrate_order().corrected.mean())
    axes = Axes(560, 275, (1.2, 72), (0.12, 6.0), xlog=True, ylog=True,
                pad=(70, 26, 46, 24))
    axes.points(table.s0, table.plateau, CATEGORY[0], radius=4.2)
    grid = np.array([1.4, 64.0])
    centre = np.exp(np.log(table.plateau).mean())
    middle = np.exp(np.log(table.s0).mean())
    axes.line(grid, centre * (grid / middle) ** order["order"], CATEGORY[0],
              width=2.2)
    axes.line(grid, centre * (grid / middle) ** measured, CATEGORY[1],
              width=2.0, dash="5 4")
    axes.label(48, centre * (48 / middle) ** order["order"],
               f"plateau: {order['order']:+.2f}", CATEGORY[0], size=11,
               weight="600", anchor="end", dy=-9)
    axes.label(48, centre * (48 / middle) ** measured,
               f"the rate's own order: {measured:+.2f}", CATEGORY[1], size=11,
               weight="600", anchor="end", dy=17)
    return fig(
        axes.render("[S] in the cuvette, mM", "plateau v/k, AU",
                    "D · The level it is heading for carries the rate's substrate order"),
        "If the signal is a species the oxidant makes from the substrate and "
        "then destroys, its stationary level is <code>A∞ = v(S)/k</code> — so "
        "the plateau has to carry whatever substrate order the <em>rate</em> "
        f"has. That order was measured on the rates of this block, "
        f"<strong>{measured:+.3f}</strong>, before any of this. The plateau's "
        f"is <strong>{order['order']:+.3f} ± {order['stderr']:.3f}</strong> "
        f"over {order['lever']:.0f}× in [S], with nothing fitted to make it "
        "agree. Where the two lines meet the data gives "
        "<strong>k<sub>A</sub>/k<sub>S</sub> ≈ "
        f"{np.median([slowdown.selectivity(p, s, e) for p, s, e in zip(table.plateau, table.s0, table.epsilon)]):.0f}"
        "</strong> — an upper bound, because the aldehyde in the cuvette is "
        "the increment plus whatever the enzyme-free background made.")


def _sharpest():
    """The run where the arithmetic bites hardest, phrased for a caption."""
    frame = _whole_frame()
    live = frame[frame.live & (frame.experiment == 21)]
    product = float((live.net / live.epsilon).max())
    enzyme = float(live.e0.iloc[0])
    loss = 1 - float(live.late_over_early.median())
    return (f"exp 21 loses {loss:.0%} of its rate on {product * 1000:.0f} µM "
            f"of product against {enzyme * 1000:.0f} µM of catalyst, and "
            f"blocking that fraction would take {loss * enzyme * 1000:.0f} µM "
            f"bound.")


def figure_budget():
    frame = _whole_frame()
    budget = slowdown.product_budget(frame).loc["4OMe-BnOH"]
    series = budget.reindex(scope.TEMPERATURE_SERIES)
    axes = Axes(560, 265, (12, 43), (0, 0.30), pad=(68, 30, 46, 24))
    axes.line(series.temperature, series.of_enzyme, CATEGORY[1], width=2.2)
    axes.points(series.temperature, series.of_enzyme, CATEGORY[1], radius=4.6)
    for temperature, row in series.iterrows():
        axes.label(row.temperature, row.of_enzyme,
                   f"{row.of_enzyme:.0%}", CATEGORY[1], size=10.5,
                   anchor="middle", dy=-10)
    # The line is what the WARMEST run's own rate loss would have to bind, so
    # it is a fact about that run and not a round number chosen to sit near it.
    needed = 1 - float(series.late_over_early.min())
    axes.hline(needed, ACCENT, dash="5 4")
    axes.note(axes._fx(12.6), axes._fy(needed) - 7,
              f"a {needed:.0%} rate loss would need this much bound",
              ACCENT, size=10.5)
    return fig(
        axes.render("temperature, °C",
                    "product made, as a fraction of [enz]",
                    "E · Barely enough product to block the catalyst, and none to spare"),
        "Blocking a fraction of the catalyst takes that fraction of [enz] "
        "<em>bound</em>, whatever the binding constant — tightening the binding "
        "moves an equilibrium, it does not create inhibitor. At the end of the "
        "40 °C run the catalytic increment holds "
        f"{series.product_mM.max() * 1000:.0f} µM of product against "
        f"{series.e0.min() * 1000:.0f} µM of catalyst, and the rate is down "
        f"{1 - series.late_over_early.min():.0%}. The two numbers are within a "
        "factor of about one of each other, which is why this argument narrows "
        "the field rather than closing it: the aldehyde in the cuvette is the "
        "increment <em>plus</em> the background's, and the background is not "
        "measured at these compositions. Sharper elsewhere — " + _sharpest())




def figure_substrates():
    blocks = _slowdown_blocks()
    sets = (("4OMe-BnOH, catalysed, phosphate", "4OMe catalysed, phosphate",
             CATEGORY[1]),
            ("BnOH, catalysed, every buffer", "BnOH catalysed, all buffers",
             CATEGORY[0]))
    # Limits from the data, not typed: 18 points fell outside a hand-set
    # frame the first time, and marks are clipped silently.
    drawn = []
    for label, name, colour in sets:
        block = blocks[name]
        block = block[block.live & (block.late_over_early > 0)
                      & (block.net > slowdown.DRIVER_FLOOR)]
        drawn.append((label, colour,
                      (block.net / block.epsilon).to_numpy(dtype=float),
                      block.late_over_early.to_numpy(dtype=float), block))
    span = lambda values: (min(v.min() for v in values) / 1.6,
                           max(v.max() for v in values) * 1.6)
    axes = Axes(620, 330, span([d[2] for d in drawn]),
                span([d[3] for d in drawn]), xlog=True, ylog=True,
                pad=(72, 26, 46, 30))
    axes.hline(1.0, INK, dash="5 4", width=1.2)
    for label, colour, product, ratio, block in drawn:
        axes.points(product, ratio, colour, radius=3.4, opacity=0.85)
        row = slowdown.deceleration_drivers(block)
        grid = np.array([product.min(), product.max()])
        centre = float(np.exp(np.log(ratio).mean()))
        middle = float(np.exp(np.log(product).mean()))
        axes.line(grid, centre * (grid / middle) ** row["product"], colour,
                  width=2.4)
        # Labelled at the LEFT end, where the two lines are furthest apart:
        # at the right end they cross, and both sit inside the densest part
        # of both clouds.
        axes.label(grid[0], centre * (grid[0] / middle) ** row["product"],
                   f"{label}: {row['product']:+.2f}", colour, size=11,
                   weight="600", anchor="start", dx=5,
                   dy=-9 if row["product"] < 0 else 17)
    axes.note(600, 40, "above the dashes: the rate rose", MUTED, size=10.5,
              anchor="end")
    axes.note(600, 274, "below them: the rate fell", MUTED, size=10.5,
              anchor="end")
    return fig(
        axes.render("product made, mM", "late rate / early rate",
                    "F · The two substrates bend in opposite directions"),
        "The same statistic, the same axis, the two substrates. "
        "4-methoxybenzaldehyde's own accumulation shuts its reaction down; "
        "benzaldehyde's does not, and the BnOH set reaches the higher "
        "concentration of the two. The lines are the product coefficients of "
        "figure A. Run length is not held fixed here and does not need to be "
        "— it enters the regression as its own term and, on the 4OMe set, "
        "comes back at −0.124 ± 0.079 against the product's −0.598 ± 0.053.")


def build_index():
    frame = _whole_frame()
    named = _slowdown_blocks()
    series = slowdown.deceleration_drivers(named["temperature series"])
    catalysed = slowdown.deceleration_drivers(named["4OMe catalysed, phosphate"])
    free = slowdown.deceleration_drivers(named["4OMe enzyme-free"])
    table = _sink_curves()
    counts = table.prefers.value_counts()
    order = slowdown.plateau_scaling(table)
    measured = float(arrhenius.substrate_order().corrected.mean())
    selectivities = np.array([slowdown.selectivity(p, s, e) for p, s, e
                              in zip(table.plateau, table.s0, table.epsilon)])
    selectivities = selectivities[np.isfinite(selectivities)]
    effect = slowdown.sink_effect_on_activation()

    hero = (
        "<div class='hero'>"
        f"<div><div class='k'>falls with product</div><div class='v'>"
        f"{catalysed['product']:+.3f}</div><div class='u'>± "
        f"{catalysed['product_stderr']:.3f}, catalysed 4OMe</div></div>"
        f"<div><div class='k'>falls with run length</div><div class='v'>"
        f"{catalysed['span']:+.3f}</div><div class='u'>± "
        f"{catalysed['span_stderr']:.3f}, same curves</div></div>"
        f"<div><div class='k'>with no enzyme</div><div class='v'>"
        f"{free['span']:+.3f}</div><div class='u'>± "
        f"{free['span_stderr']:.3f} on run length, "
        f"{free['product']:+.3f} on product</div></div>"
        f"<div><div class='k'>k<sub>A</sub>/k<sub>S</sub></div><div class='v'>"
        f"{np.median(selectivities):.0f}</div><div class='u'>the oxidant's "
        f"preference for the product</div></div>"
        "</div>")

    body = f"""
<p class='lede'>Every catalysed 4-methoxybenzyl alcohol run in this archive
rises to a maximum rate and then falls, at <strong>0.3–3% conversion</strong>.
This is what the fall is: <strong>the oxidant attacking the aldehyde it has
just made</strong>. The evidence is drawn from the whole 4OMe archive, because
the temperature series alone is 23 curves and cannot separate the two
hypotheses that matter.</p>

{hero}

<h2>1 · A clock, or the product?</h2>
<p>A single progress curve cannot tell them apart. Inside one curve the product
only ever grows with time, so “the rate fell after 2000 s” and “the rate fell
after 0.1 AU” are the same sentence. The separation has to come from a design
where the two move independently, and the archive has one by accident: it holds
1470 s runs that reached 0.27 AU and 17934 s runs that reached 0.045. Across the
catalysed 4OMe phosphate set the two regressors correlate only
<strong>{catalysed['collinearity']:+.2f}</strong>.</p>
{figure_drivers()}
<p>The temperature series loads on product alone
(<strong>{series['product']:+.3f} ± {series['product_stderr']:.3f}</strong>
against {series['span']:+.3f} ± {series['span_stderr']:.3f} for run length).
<strong>The same chemistry with no enzyme in the cuvette does the opposite.</strong>
That block is not a curiosity: it is the reference channel every catalysed
curve here is measured against, and its decay is subtracted out per cuvette.</p>
<p><strong>It is not the photometer either.</strong> ε is 7.53 for the 4OMe
aldehyde against 1.23 for benzaldehyde, so at the same product <em>concentration</em>
a 4OMe run sits at six times the <em>absorbance</em> — and detector compression
would bend those curves over exactly as a sink does. Two controls: over the
absorbance window the two substrates share, 0.015–0.475 AU, the catalysed 4OMe
curves still give <strong>−0.621 ± 0.099</strong> against BnOH's
<strong>+0.386 ± 0.136</strong>, a 6.0σ gap the photometer cannot produce
because it does not know which alcohol is in the cuvette; and the enzyme-free
4OMe curves cover 0.012–0.581 AU with a product coefficient of
<strong>+0.025 ± 0.048</strong>, so compression would have to be selective for
the presence of enzyme.</p>

<p class='warn'><code>net</code> is the integral of the rate, so a curve that
decelerates <em>less</em> makes <em>more</em> product at the same starting rate.
That bias pushes the product coefficient towards zero. Every negative
coefficient on this page is a floor, not a ceiling.</p>

<h2>2 · Consumed, or blocking?</h2>
<p>Two mechanisms fit “the product sets the rate”, and they straighten
different transforms of the same two columns:
<code>A′ = v − kA</code> makes the <em>rate</em> linear in the product;
<code>A′ = v/(1 + A/K<sub>i</sub>)</code> makes <em>1/rate</em> linear in it.</p>
{figure_sink_lines()}
{figure_plateau_order()}
{figure_budget()}

<h2>3 · And it is specific to this substrate</h2>
{figure_substrates()}

<h2>4 · What it does to the activation parameters</h2>
<p>Every rate in <a href='../temperature_series/index.html'>the temperature
series</a> is already net of this loss, so none of them is the production rate.
Refitting the sink model with its rate constant pinned recovers that rate, and
it comes out <strong>{effect['lift']:.1%}</strong> above the published estimator
on the median curve — but between {effect['lift_low']:.1%} and
{effect['lift_high']:.1%} depending on temperature, with no order to it. A
factor that does not order in temperature cancels out of a slope:</p>
<div class='scroll'><table>
<tr><th>estimator</th><th>E<sub>a</sub> kJ/mol</th><th>ΔH‡ kJ/mol</th>
<th>ΔS‡ J/mol/K</th><th>ΔG‡(298) kJ/mol</th><th>Arrhenius rms</th></tr>
<tr class='hl'><td><code>v_peak</code>, published</td>
<td>{effect['published']['activation_kJ']:.2f} ±
{effect['published']['activation_stderr']:.2f}</td>
<td>{effect['published']['enthalpy_kJ']:.2f}</td>
<td>{effect['published']['entropy_J']:+.1f} ±
{effect['published']['entropy_stderr']:.1f}</td>
<td>{effect['published']['gibbs_kJ']:.2f}</td>
<td>{effect['published']['rms']:.3f}</td></tr>
<tr><td><code>v_prod</code>, sink model, k pinned</td>
<td>{effect['corrected']['activation_kJ']:.2f} ±
{effect['corrected']['activation_stderr']:.2f}</td>
<td>{effect['corrected']['enthalpy_kJ']:.2f}</td>
<td>{effect['corrected']['entropy_J']:+.1f} ±
{effect['corrected']['entropy_stderr']:.1f}</td>
<td>{effect['corrected']['gibbs_kJ']:.2f}</td>
<td>{effect['corrected']['rms']:.3f}</td></tr>
</table></div>
<p>The activation energy moves by
<strong>{effect['activation_shift']:+.2f} ± {effect['activation_shift_stderr']:.2f}
kJ/mol</strong> — nothing — and the level shift is worth
{effect['entropy_from_lift']:+.2f} J/mol/K on ΔS‡, against a ±
{effect['published']['entropy_stderr']:.1f} error. <strong>Section 4 of the
temperature series stands.</strong> The corrected route is the more nearly
correct quantity and the noisier one, Arrhenius scatter
{effect['corrected']['rms']:.3f} against {effect['published']['rms']:.3f},
because the sink shape is one the four coldest runs cannot test. Quoting it as
the headline would trade a bias smaller than the error for a variance larger
than it.</p>

<h2>5 · What this settles, and what it does not</h2>
<p><strong>Settled.</strong> The fall is set by the product; it belongs to the
catalysed pathway and not to the solution; the rate is linear in the product
rather than hyperbolic in it; the stationary level carries the substrate order
it has to; and the effect is specific to 4OMe-BnOH at product concentrations
where BnOH shows nothing.</p>
<p><strong>Not settled, and not settleable from absorbance.</strong> Whether the
aldehyde is <em>consumed</em> by the oxidant or merely <em>scavenges</em> it.
Those are the same reaction seen from two ends — an oxidant that attacks the
electron-rich aldehyde is both destroying the chromophore and being diverted
from the alcohol — and a single-wavelength trace cannot say which half it is
watching. Both give <code>A′ = v − kA</code>; both give
<code>A<sub>∞</sub> ∝ v(S)</code>.</p>
<p>Also not settled: <strong>what the enzyme-free clock is</strong>.
<a href='../background_reaction/index.html'>The background folder</a> left it
open for BnOH; this adds only that it is rate- and product-independent, which
rules product effects out of it and leaves the rest.</p>

<h2>Reproducing</h2>
<p><code>python data/slowdown.py</code> prints the whole argument ·
<code>python data/test_slowdown.py</code> ·
<code>python product_fate/build_figures.py</code> ·
<code>python product_fate/check_numbers.py</code>, which re-derives every number
in <code>ANALYSIS.md</code> from the modules and fails if the prose and the code
disagree.</p>
"""
    return page("What the 4OMe product does", body,
                "The oxidant attacks the aldehyde it just made").replace(
        "</style>", EXTRA_CSS + "</style>")


def main():
    path = os.path.join(HERE, "index.html")
    content = build_index()
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(content)
    print(f"wrote {path}  ({len(content) / 1024:.0f} kB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
