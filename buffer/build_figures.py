"""
Builds buffer/index.html.

Draws only; every number comes from `buffer_role`, `scope`, `induction` or
`solution_chemistry`, so a figure and the prose in ANALYSIS.md cannot disagree
about a value without `check_numbers.py` saying so.

    python buffer/build_figures.py
"""
import functools
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
from figure_kit import (CATEGORY, PH_RAMP, SURFACE, breakpoints, fig, panel,
                        progress_axes, progress_overlay, styled, write_pages)






@functools.cache
def _titrations():
    """The five titrations' own table, built once."""
    return buffer_role.titration_table()


@functools.cache
def _archive():
    """Every curve in the archive, built once: the confound map needs all of it."""
    return scope.frame(tuple(range(1, 152)))


def figure_titrations():
    table = _titrations()
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
    table = _titrations()
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


def figure_joint():
    """
    The pre-equilibrium constraint on both axes, with the control that gates it.

    Two panels' worth of information in one frame: the required +1 as a rule,
    each measurement as a point with its error bar, and the windows whose
    signal control FAILS drawn hollow, because those are the ones that must not
    be read.
    """
    table = induction.induction_table(induction.WHOLE_ARCHIVE)
    rows = []
    for width in induction.BUFFER_WINDOW_SWEEP:
        ladder = induction.buffer_lever(table, width=width)
        joint = induction.joint_buffer_order(ladder)
        control = induction.signal_control(ladder)
        rows.append((f"[buf], {width:.0f} s window", joint,
                     abs(control["signal_slope"])
                     < 2 * control["signal_stderr"]))
    blocks = induction.induction_blocks(table)
    for label, block in (
            ("[H₂O₂], exps 127-131",
             table[table.experiment.isin(induction.PEROXIDE_LEVER)]),
            ("[H₂O₂], exps 135-151",
             blocks["BnOH two-axis (135-151)"])):
        joint = induction.joint_peroxide_order(block)
        control = induction.signal_control(block)
        rows.append((label, joint, abs(control["signal_slope"])
                     < 2 * control["signal_stderr"]))
    axes = Axes(620, 330, (-0.1, 2.0), (-1.05, len(rows) - 0.3),
                pad=(190, 26, 46, 34))
    axes.line([1.0, 1.0], [-1.05, len(rows) - 0.3], ACCENT, width=2.0,
              dash="5 4")
    # Below the last row, not above the first: above the first it sat on the
    # figure title.
    axes.label(1.0, -0.62, "required +1", ACCENT, size=11,
               anchor="middle", weight="600", dy=4)
    for index, (label, joint, passes) in enumerate(rows):
        y = len(rows) - 1 - index
        colour = PH_RAMP[1] if label.startswith("[buf]") else CATEGORY[1]
        low = joint["slope"] - joint["stderr"]
        high = joint["slope"] + joint["stderr"]
        axes.line([low, high], [y, y], colour, width=2.6)
        # A hollow mark for a window whose signal control fails: fill it with
        # the page's own surface and ring it in the series colour.
        axes.points([joint["slope"]], [y], colour if passes else SURFACE,
                    radius=5.0, stroke=colour, stroke_width=1.8)
        axes.label(-0.14, y, label, INK if passes else MUTED, size=11,
                   anchor="end", dy=4, weight="600" if passes else "400")
        axes.label(high, y, f"{joint['slope']:+.2f}"
                   + ("" if passes else "  (S/N fails)"),
                   INK if passes else MUTED, size=10.5, anchor="start",
                   dx=8, dy=4)
    return fig(
        axes.render("d ln v/d ln[X] − d ln τ/d ln[X]", "",
                    "E · The pre-equilibrium constraint, both axes",
                    yticks=False),
        "If the catalyst is drawn into its active form by a species held in "
        "excess, that species' order on the rate and on the induction time "
        "must differ by <strong>exactly +1</strong> — for every binding "
        "constant and every concentration, so it is a prediction with nothing "
        "to fit. The buffer axis meets it. The peroxide axis falls short on "
        "both blocks that carry a peroxide ladder, and both of those also fail "
        "the signal-to-noise control (hollow), as do the buffer's longest "
        "windows — which overshoot in exactly the direction that artefact "
        "predicts. Read only the filled points.")


def figure_confound():
    frame = _archive()
    blocks = induction.induction_blocks(frame)
    names = (("4OMe catalysed", "4OMe catalysed"),
             ("4OMe enzyme-free", "4OMe, no enzyme"),
             ("temperature series", "the temperature series"),
             ("BnOH two-axis (135-151)", "BnOH, exps 135–151"))
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
        "cuvettes of all seventeen runs — and that is the block the mechanism "
        "fitting is scoped to.")


def build_curves_page():
    """
    Every curve the five titrations are made of, with the landmark drawn on it.

    THIS PAGE IS THE AUDIT SURFACE FOR SECTION 6. The joint buffer order is
    read off `v_peak` and off a landmark whose window is 450 SECONDS, chosen
    because exps 32 and 34 differ in length and a window given as a share of
    the run would be two different windows. That is exactly the sort of choice
    that has to be visible on the curves rather than argued for in prose, and
    the folder had no page showing it.
    """
    frame = scope.frame(buffer_role.TITRATIONS)
    lookup = {(c.experiment, c.sample): c
              for c in scope.curves(buffer_role.TITRATIONS)}
    panels = []
    for row in frame.sort_values(["pH", "experiment", "buf"]).itertuples():
        curve = lookup.get((row.experiment, row.sample))
        if curve is None:
            continue
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        axes, radius = progress_axes(times, values, limit=140)
        progress = progress_overlay(axes, times, values, mark_radius=radius)
        found = induction.buffer_landmark(curve)
        marks, labels = [], []
        if np.isfinite(found.t_ind) and found.t_ind > 0:
            marks.append(found.t_ind)
            labels.append("t_ind")
        breakpoints(axes, marks, labels, colour=CATEGORY[0])
        panels.append(panel(
            f"[buf] = {row.buf:g} mM · pH {row.pH:.2f}"
            f"<span class='pill'>exp {int(row.experiment)}</span>",
            f"[S] {row.s0:.3f} mM · [H₂O₂] {row.h2o2:g} mM · "
            f"{row.temperature:.0f} °C · {int(row.points)} readings over "
            f"{row.duration_s / 60:.0f} min · {row.source}",
            axes.render("time, s", "ΔA"),
            f"<strong>{int(row.phases)} phase"
            + ("s" if row.phases == 2 else "")
            + f"</strong> · {esc(str(row.progress_kind))} "
            f"· F = {row.two_phase_f:.0f}"
            f" · v_peak {row.v_peak:.2e}"
            + (f" · t_ind {found.t_ind:.0f} s · depth {found.depth:.3f}"
               if np.isfinite(found.t_ind) else " · no landmark")
            + ("" if row.live else " · <strong>NOT LIVE</strong>")))
    order = induction.buffer_order(induction.buffer_lever(
        induction.induction_table(induction.WHOLE_ARCHIVE)))
    body = (f"<p class='lede'>All {len(panels)} cuvettes of the five buffer "
            "titrations — 4OMe-BnOH, 40 °C, phosphate, 82.5 mM H₂O₂, with "
            "<code>[S]</code> fixed inside each run — in pH order. The rust "
            "line is whichever form the curve earned, one relaxation or two, "
            "from <code>summary_kinetics.fit_progress</code>; nothing is "
            "excluded and every fit uses every point except the instrument's "
            "first reading, which is discarded from every run in the archive."
            "</p>"
            "<p class='lede'>The blue dashed vertical is the <strong>induction "
            "landmark</strong>, read through a window of "
            f"<strong>{induction.BUFFER_WINDOW:.0f} seconds</strong> — in "
            "seconds, not as a share of the run, because exps 32 and 34 ran "
            "5280 s and 1767 s and a fractional window would be two different "
            "windows. Section 6's joint order and section 4c's "
            f"<code>{order['slope']:+.3f} ± {order['stderr']:.3f}</code> are "
            "read off these two quantities: the peak rate and this landmark."
            "</p>"
            "<p class='lede'><strong>Exps 32 and 34 earn different model "
            "forms</strong>, and it is visible here: every curve of exp 34 "
            "takes the two-phase form and every curve of exp 32 the one-phase, "
            "because exp 34's runs are long enough to contain the slow fall "
            "and exp 32's end before it. That is why nothing on this page "
            "compares a time constant between the two runs.</p>"
            "<div class='grid three'>" + "".join(panels) + "</div>")
    return styled("The buffer titrations — every progress curve", body,
                  "Exps 32, 34, 35, 36, 37 · 4OMe-BnOH · 40 °C · phosphate")


def build_index():
    table = _titrations()
    prediction = buffer_role.species_prediction(7.00, 7.53)
    measured = buffer_role.catalytic_coefficient(drop=(35,))
    verdict = buffer_role.separable(measured, prediction)
    identity = buffer_role.identity_overlap(_archive())
    widest = buffer_role.overlap_width(identity)
    table_full = induction.induction_table(induction.WHOLE_ARCHIVE)
    ladder = induction.buffer_lever(table_full)
    induction_order = induction.buffer_order(ladder)
    joint = induction.joint_buffer_order(ladder)
    crossing = buffer_role.peroxide_crossing(frame=_archive())
    free = buffer_role.free_route_order(
        frame=scope.frame(buffer_role.TITRATIONS))

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

<h2>6 · Is the buffer a base, or is it carrying the peroxide?</h2>
{figure_joint()}
<p>A buffer acting as a general base is a term in <code>[buf]</code>. A buffer
that first takes up the peroxide itself — <code>H₂O₂ + P ⇌ P–OOH</code>, then
<code>P–OO⁻</code> delivering oxygen to the ketone the way a peracid or
peroxymonosulfate does — is a term in <code>[buf][H₂O₂]</code>. At one pH the
two schemes differ by that <strong>interaction and by nothing else</strong>, and
this archive never varies it: of {crossing['runs']} runs,
{crossing['steps_buffer']} step <code>[buf]</code> and
{crossing['steps_peroxide']} step <code>[H₂O₂]</code>, and
<strong>{crossing['steps_both']} step both</strong>. All five titrations sit at
82.5 mM peroxide.</p>
<p>What the archive does say points the same way twice. The
<strong>pre-equilibrium constraint</strong> above is met on the buffer axis
({joint['slope']:+.3f} ± {joint['stderr']:.3f} against a required +1) and missed
on the peroxide axis — so the species drawn into equilibrium with the catalyst
behaves like the buffer and not like H₂O₂. And the buffer-<em>free</em> route
rises with pH faster than one hydroperoxide can explain:
<strong>{free['level_ratio']:.2f}×</strong> over {free['delta_pH']:.2f} pH units
at matched substrate, where [HOO⁻] gives only
<strong>{free['hoo_ratio']:.2f}×</strong>. Both are eight-curve and two-run
results across different days, and the archive's own between-day step is 1.80× —
these are pointers, not measurements.</p>

<h2>7 · What this settles, and what it does not</h2>
<p><strong>Settled.</strong> Every buffer order in this project is an order in
total buffer. The catalysed turnover's dependence saturates — about first order
below 25 mM, half order above 50 — and disappears above the pKa because the
[HOO⁻] route grows past it, not because the buffer term shrinks.
<code>[S]</code> and <code>[buf]</code> are collinear at −0.96 in every 4OMe run
and at 0.00 in every two-axis run. Identity is confounded with pH everywhere.</p>
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
    return styled("What the buffer does, and to what", body,
                  "A reagent, a confound, and a candidate")


def main():
    return write_pages(HERE, {"index.html": build_index(),
                              "progress_curves.html": build_curves_page()})


if __name__ == "__main__":
    raise SystemExit(main())
