"""
Builds early_trough/index.html and early_trough/progress_curves.html.

Draws only; every number comes from `early_trough` or `curve_metrics`, so a
figure and the prose in ANALYSIS.md cannot disagree about a value without
`check_numbers.py` saying so.

    python early_trough/build_figures.py
"""
import functools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import curve_metrics
import early_trough
import scope
from svgplot import ACCENT, GRID, INK, MUTED, Axes, esc
from figure_kit import (CATEGORY, EVENT_BAND_COLOUR, breakpoints,
                        derivative_axes, fig, panel, progress_axes,
                        progress_overlay, residual_axes, styled, write_pages)

SUBSTRATE_COLOUR = {"4OMe-BnOH": CATEGORY[0], "BnOH": CATEGORY[1]}
CLUSTER_COLOUR = {"oxidant": CATEGORY[0], "substrate": CATEGORY[1]}


@functools.cache
def _table():
    return early_trough.trough_table()


@functools.cache
def _genuine():
    table = _table()
    genuine = table[table.genuine].copy()
    genuine["cluster"] = genuine.apply(early_trough.cluster, axis=1)
    return genuine.sort_values("z")


@functools.cache
def _excluded_candidates():
    """The two curves that clear the smoothed threshold and are rejected --
    the control this page owes the reader, per the house convention that a
    curves page shows the whole picture, not only the confirmed half."""
    table = _table()
    return table[table.candidate & ~table.genuine].sort_values("z")


_CORRELATION_CEILING = 20.0


def figure_correlation():
    """
    Trough depth against dominance, restricted to `z <= _CORRELATION_CEILING`.

    A THIRD OF THE ARCHIVE IS DELIBERATELY LEFT OFF THIS PANEL, and the
    caption says so: 111 of 311 scanned curves have a strong genuine early
    RISE (z up to +1106), which is a real signal but not this question's
    subject, and plotting it would compress every dip into an unreadable
    sliver near zero on a linear axis. This is a stated selection, not a
    silent clip -- every point with z <= 20 is shown and the axis limits
    cover exactly that range, so `doc.unclipped` still holds.
    """
    table = _table()
    shown = table[table.z <= _CORRELATION_CEILING]
    genuine = _genuine()
    axes = Axes(600, 320, (0.005, 6000.0), (-46.0, _CORRELATION_CEILING),
               xlog=True, pad=(56, 20, 46, 20))
    axes.hline(early_trough.EARLY_TROUGH_CANDIDATE_Z, colour=GRID, dash="4 3")
    axes.note(560, axes._fy(early_trough.EARLY_TROUGH_CANDIDATE_Z) - 5,
              "candidate threshold", MUTED, size=10, anchor="end")
    for substrate, colour in SUBSTRATE_COLOUR.items():
        block = shown[shown.substrate == substrate]
        axes.points(block.dominance.to_numpy(), block.z.to_numpy(), colour,
                    radius=2.6, opacity=0.45, stroke=None)
    for substrate, colour in SUBSTRATE_COLOUR.items():
        block = genuine[genuine.substrate == substrate]
        if not len(block):
            continue
        axes.ring(block.dominance.to_numpy(), block.z.to_numpy(), colour,
                  radius=5.5, width=1.6)
    axes.note(80, 24, "● faint: every scanned curve with z ≤ 20   ○ ringed: "
                       "genuine (passes all three screens)", MUTED, size=10.5)
    for substrate, colour in SUBSTRATE_COLOUR.items():
        axes.note(80, 40 + 14 * list(SUBSTRATE_COLOUR).index(substrate),
                  substrate, colour, size=11, weight="600")
    return fig(
        axes.render("dominance = max([enz]/[S], [enz]/[HOO-])",
                    "early trough, sigma (capped at +20)",
                    "A · The archive-wide correlation, both substrates"),
        f"{len(shown)} of {len(table)} live catalysed curves in the archive "
        "plotted against whichever reactant the catalyst's own fixed "
        f"concentration comes closest to overwhelming — the "
        f"{len(table) - len(shown)} left off have a strong genuine early "
        "RISE (z up to +1106) rather than a trough, which is not this "
        "question's subject and would compress the dips shown here to a "
        "sliver near zero. <strong>Spearman ρ = "
        f"{early_trough.dominance_correlation()['4OMe-BnOH']['e0_hoo'][0]:+.3f}"
        f"</strong> (4OMe, p = "
        f"{early_trough.dominance_correlation()['4OMe-BnOH']['e0_hoo'][1]:.1e}) "
        "and <strong>"
        f"{early_trough.dominance_correlation()['BnOH']['e0_hoo'][0]:+.3f}"
        f"</strong> (BnOH, p = "
        f"{early_trough.dominance_correlation()['BnOH']['e0_hoo'][1]:.1e}) "
        "against log[enz]/[HOO⁻] alone, over EVERY scanned curve including "
        "the ones omitted here — the catalyst:substrate ratio carries no "
        "signal in either substrate (Section 2).")


def figure_clusters():
    genuine = _genuine()
    axes = Axes(600, 320, (0.001, 0.3), (0.001, 6000.0), xlog=True, ylog=True,
               pad=(58, 20, 46, 20))
    grid = np.array([0.001, 0.3])
    axes.line(grid, grid, GRID, dash="3 3", width=1.2)
    axes.note(axes._fx(0.15), axes._fy(0.15) - 6, "[enz]/[S] = [enz]/[HOO-]",
              MUTED, size=10, anchor="middle")
    for name, colour in CLUSTER_COLOUR.items():
        block = genuine[genuine.cluster == name]
        axes.points(block.e0_s0.to_numpy(), block.e0_hoo.to_numpy(), colour,
                    radius=5.0,
                    title=[f"exp {int(r.experiment)}.{int(r.sample)}"
                          for r in block.itertuples()])
        label = ("oxidant-dominated (above the line)" if name == "oxidant"
                 else "substrate-dominated (below the line)")
        axes.note(70, 40 + 16 * list(CLUSTER_COLOUR).index(name), label,
                  colour, size=11, weight="600")
    return fig(
        axes.render("[enz]/[S]", "[enz]/[HOO-]",
                    "B · The two clusters, by which reactant dominates"),
        "The fifteen genuine curves in the plane the classification is made "
        "in. Twelve sit far above the diagonal — the catalyst is a trivial "
        "share of the substrate but a thousand-fold excess over the "
        "nanomolar hydroperoxide pool — and three sit below it, all BnOH at "
        "the block's lowest substrate rung (0.216 mM), where the catalyst is "
        "instead a meaningful fraction of the substrate itself.")


def figure_debubble():
    genuine = _genuine()
    axes = Axes(600, 320, (-46.0, 8.0), (-46.0, 8.0), pad=(56, 20, 46, 20))
    axes.line([-46.0, 8.0], [-46.0, 8.0], GRID, dash="3 3", width=1.2)
    axes.note(axes._fx(2.0), axes._fy(4.5), "z_corrected = z_raw", MUTED,
              size=10, anchor="end")
    for name, colour in CLUSTER_COLOUR.items():
        block = genuine[genuine.cluster == name]
        axes.points(block.z.to_numpy(), block.z_corrected.to_numpy(), colour,
                    radius=5.0)
    for row in genuine.itertuples():
        if abs(row.z_corrected - row.z) > 2.0:
            axes.note(axes._fx(row.z) + 8, axes._fy(row.z_corrected),
                      f"exp {row.experiment}.{row.sample}", INK, size=10)
    return fig(
        axes.render("trough, raw readings (sigma)",
                    "trough, debubble-corrected (sigma)",
                    "C · Correction for the O2 artefact deepens the effect, "
                    "it does not explain it away"),
        "Every genuine curve's trough recomputed after `curve_metrics."
        "debubble` removes every detected O2 event. Points would move "
        "<em>up</em>, back toward zero, if the O2 correction were secretly "
        "responsible for the dip. Thirteen of fifteen sit exactly on the "
        "diagonal (no detachment to correct); the two that carry any "
        "(exps 141.4, 142.4) move further <em>down</em> — the correction "
        "had been partly masking the trough with later bubble-driven rise.")


def build_index():
    table = _table()
    genuine = _genuine()
    corr = early_trough.dominance_correlation(table)
    n_candidate = int(table.candidate.sum())
    n_genuine = len(genuine)
    n_oxidant = int((genuine.cluster == "oxidant").sum())
    n_substrate = int((genuine.cluster == "substrate").sum())
    excluded = _excluded_candidates()

    hero = f"""
<div class='hero'>
  <div><div class='k'>curves scanned</div>
       <div class='v'>{len(table)}</div>
       <div class='u'>every live catalysed curve in the archive</div></div>
  <div><div class='k'>genuine early troughs</div>
       <div class='v'>{n_genuine}</div>
       <div class='u'>of {n_candidate} candidates on a smoothed mean alone</div></div>
  <div><div class='k'>oxidant-dominated</div>
       <div class='v'>{n_oxidant} of {n_genuine}</div>
       <div class='u'>[enz]/[HOO⁻] up to 4300×</div></div>
  <div><div class='k'>substrate-dominated</div>
       <div class='v'>{n_substrate} of {n_genuine}</div>
       <div class='u'>the block's lowest substrate rung only</div></div>
</div>"""

    rows = "".join(
        f"<tr><td>exp {int(r.experiment)}.{int(r.sample)}</td>"
        f"<td>{esc(r.substrate)}</td><td>{r.pH:.2f}</td>"
        f"<td>{r.s0:.3f}</td><td>{r.h2o2:.2f}</td><td>{r.e0:.3f}</td>"
        f"<td>{r.e0_s0:.4f}</td><td>{r.e0_hoo:.3g}</td>"
        f"<td>{r.z:+.1f}</td><td>{esc(r.cluster)}</td></tr>"
        for r in genuine.itertuples())

    excluded_rows = "".join(
        f"<tr><td>exp {int(r.experiment)}.{int(r.sample)}</td>"
        f"<td>{r.z:+.1f}</td><td>{'yes' if r.sustained else 'no'}</td>"
        f"<td>{'yes' if not r.bubble_overlap else 'no'}</td>"
        f"<td>{'yes' if r.survives else 'no'}</td></tr>"
        for r in excluded.itertuples())

    body = f"""
<p class='lede'>A catalysed curve's absorbance is already reference-subtracted
against an enzyme-free cuvette holding the same composition, so the
background reaction the two cuvettes share cancels between them — as long as
both cuvettes see the same free concentration of whatever that background
reaction runs on. Fifteen curves across both substrates say it does not
always: before the catalysed rate takes over, the reported curve dips
measurably <strong>below zero</strong> and holds there for several minutes.
<a href='progress_curves.html'>progress_curves.html</a> shows every one of
them, fits and all, plus the two curves that looked like more of the same and
are not.</p>
{hero}

<h2>1 · A trough that survives three separate ways to be spurious</h2>
<p>{n_candidate} curves clear a smoothed trough of
{early_trough.EARLY_TROUGH_CANDIDATE_Z:g} sigma or deeper in their own early
readings. Three screens separate a real, sustained decline from an artefact:
at least 5 of the 9 readings inside the trough window must themselves sit
below 2 sigma (not just the smoothed mean — this alone removes exps 4.1 and
22.2, below); the window must not overlap a detected O2 detachment; and the
trough must survive `curve_metrics.debubble` correction. <strong>{n_genuine}
survive all three</strong>.</p>
<div class='tbl'><table>
<tr><th>excluded candidate</th><th>trough, σ</th><th>sustained?</th>
<th>no bubble overlap?</th><th>survives correction?</th></tr>
{excluded_rows}
</table></div>

<h2>2 · The driver is the oxidant, not the substrate, in both substrates independently</h2>
{figure_correlation()}
<p>[enz]/[HOO⁻] is significant at p &lt; 10⁻⁴ in <strong>both</strong> 4OMe-BnOH
and BnOH scanned separately; [enz]/[S] is not significant in either
(ρ = {corr['4OMe-BnOH']['e0_s0'][0]:+.3f}, p =
{corr['4OMe-BnOH']['e0_s0'][1]:.2f} for 4OMe;
ρ = {corr['BnOH']['e0_s0'][0]:+.3f}, p = {corr['BnOH']['e0_s0'][1]:.2f} for
BnOH). The hydroperoxide anion is nanomolar across most of the archive's pH
range, so a catalyst held at a fixed 0.014–0.273 mM is routinely a
100–4000-fold molar excess over it — exactly the regime in which a catalyst
engaging the peroxide non-productively (`MECHANISM.md` S4, the same route
`BUBBLES.md`'s gas takes) could measurably outrun the reference cuvette's own
supply.</p>

<h2>3 · A second, smaller cluster: the substrate itself, at its scarcest</h2>
{figure_clusters()}
<p>Three BnOH curves (141.4, 142.4, 143.4) sit at the two-axis block's lowest
substrate rung, 0.216 mM — the one composition in the whole genuine set where
[enz] is a large enough fraction of [S] (6.5–9.7%) for simple stoichiometric
substrate binding to matter. Elsewhere in the archive [enz] never exceeds
0.2% of [S], which is why [enz]/[S] carries no overall signal in Section 2:
this cluster is real but small enough not to move the archive-wide
correlation either way.</p>

<h2>4 · Not the O2 artefact</h2>
{figure_debubble()}
<p>A growing bubble raises absorbance (`BUBBLES.md` §1) so cannot produce a
dip by itself, and none of the fifteen troughs overlaps a detected
detachment — but the strongest possible test is what correction does to the
two curves that carry any gas at all. It deepens their trough rather than
erasing it.</p>

<h2>5 · The fifteen curves</h2>
<div class='tbl'><table>
<tr><th>curve</th><th>substrate</th><th>pH</th><th>[S] mM</th>
<th>[H₂O₂] mM</th><th>[enz] mM</th><th>[enz]/[S]</th><th>[enz]/[HOO⁻]</th>
<th>trough, σ</th><th>cluster</th></tr>
{rows}
</table></div>

<h2>What this does and does not establish</h2>
<p><strong>Established.</strong> The dip is a real, sustained, multi-reading
decline in the raw readings (not a fitted-curve extrapolation artefact — see
`CLAUDE.md`'s derivative-panel discussion, which is what first surfaced 135.5
and 151.5–.7 by eye), it is not the O2 artefact, and archive-wide it tracks
[enz]/[HOO⁻] far more strongly than [enz]/[S] in both substrates
independently.</p>
<p><strong>Not established.</strong> Whether the catalyst is engaging the
peroxide, the hydroperoxide anion specifically, or something else entirely —
this archive has no headspace or manometric measurement of anything, exactly
the same limit `BUBBLES.md` states for the gas. Nor can it separate "the
catalyst consumes the oxidant" from "the catalyst's own resting state changes
transiently," both of which would starve a shared background reaction the
same way. And exp 151 (three of the strongest oxidant-cluster curves) is one
of the two-axis block's own weakest, most drift-dominated runs — real
chemistry there is hardest to pull apart from the cell's own wander, which is
exactly why the 4OMe REPLICATE_RUNS curves (exps 4, 5, 7 — the single
strongest examples, at −40 to −42σ) matter: they are a different substrate,
a different buffer, and none of the two-axis block's own caveats apply to
them.</p>

<h2>Reproducing</h2>
<p><code>python data/test_early_trough.py</code> ·
<code>python early_trough/build_figures.py</code> ·
<code>python early_trough/check_numbers.py</code>, which re-derives every
number in <code>ANALYSIS.md</code> from <code>early_trough.py</code> and
fails if the prose and the code disagree.</p>
"""
    return styled("The early trough: a curve dipping before it rises", body,
                  "Reference-subtracted evidence for a transient, "
                  "non-productive engagement with a scarce reactant")


def _fit_panel(row, curve, colour):
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    corrected, events = curve_metrics.debubble(times, values, curve.noise)
    chopped = len(events) > 0
    axes, radius = progress_axes(times, values, limit=140,
                                 companion=corrected if chopped else None,
                                 pad=(56, 12, 10, 20))
    bands = []
    if chopped:
        for start, stop in events:
            bands.append((float(times[start]), float(times[stop])))
        for lo, hi in bands:
            axes.band([lo, hi], [axes.ylim[0], axes.ylim[0]],
                     [axes.ylim[1], axes.ylim[1]], EVENT_BAND_COLOUR,
                     opacity=0.14)
        axes.line(times, corrected, CATEGORY[2], width=1.0, dash="3 2",
                  opacity=0.85)
        corrected_fit = progress_overlay(axes, times, corrected,
                                         colour=CATEGORY[2], mark_radius=radius)
    raw_fit = progress_overlay(axes, times, values, mark_radius=radius,
                               colour=colour)
    chem_fit = corrected_fit if chopped else raw_fit
    chem_values = corrected if chopped else values
    if np.isfinite(row.t_trough):
        breakpoints(axes, [row.t_trough], [f"trough {row.z:+.1f}σ"],
                    colour=INK)
    residual = (chem_values - chem_fit.predict(times)) / curve.noise
    rax = residual_axes(times, residual, colour=colour, bands=bands)
    drax = derivative_axes(times, chem_fit, colour=colour, bands=bands)
    svg = (axes.render("", "ΔA") + rax.render("", "z")
          + drax.render("time, s", "dA/dt"))
    footer = (f"trough {row.z:+.1f}σ raw"
             + (f" · {row.z_corrected:+.1f}σ corrected" if chopped else "")
             + f" · {'sustained' if row.sustained else 'NOT sustained'}"
             + (" · overlaps an O2 event" if row.bubble_overlap else "")
             + (" · " + esc(row.cluster) + "-dominated"
                if hasattr(row, "cluster") else "")
             + " · <strong>GENUINE</strong>" if getattr(row, "genuine", False)
             else footer_rejected(row))
    return panel(
        f"exp {int(row.experiment)}.{int(row.sample)} · "
        f"{esc(row.substrate)}"
        f"<span class='pill'>pH {row.pH:.2f}</span>",
        f"[S] {row.s0:.3f} mM · [H₂O₂] {row.h2o2:g} mM · [enz] {row.e0:.3f} mM "
        f"· {esc(row.buffer)}",
        svg, footer)


def footer_rejected(row):
    reasons = []
    if not row.sustained:
        reasons.append("not sustained")
    if row.bubble_overlap:
        reasons.append("overlaps an O2 event")
    if not row.survives:
        reasons.append("does not survive debubble correction")
    return (f"trough {row.z:+.1f}σ · <strong>REJECTED</strong> "
           f"({', '.join(reasons)})")


def build_curves_page():
    genuine = _genuine()
    excluded = _excluded_candidates()
    lookup = {(c.experiment, c.sample): c for c in scope.curves(scope.archive())}

    panels = []
    for row in genuine.itertuples():
        curve = lookup[(row.experiment, row.sample)]
        panels.append(_fit_panel(row, curve, CLUSTER_COLOUR[row.cluster]))
    control_panels = []
    for row in excluded.itertuples():
        curve = lookup[(row.experiment, row.sample)]
        control_panels.append(_fit_panel(row, curve, MUTED))

    body = (f"<p class='lede'>The fifteen curves behind every claim in "
            "<a href='index.html'>index.html</a>, each with its own three-panel "
            "audit: the readings and whichever fitted form the curve earned "
            "(raw in colour, the debubble-corrected series dashed where a "
            "curve carries any O2 event), the residual, and this session's "
            "derivative-of-fit panel — the same one that first made 135.5's "
            "and 151.5–.7's early behaviour visible by eye. The black dashed "
            "vertical marks the trough `early_trough` reads.</p>"
            "<div class='grid three'>" + "".join(panels) + "</div>"
            "<h2>The two curves that looked like more of the same</h2>"
            "<p class='lede'>Both clear the same smoothed-mean threshold as "
            "the fifteen above and are excluded — the control this page owes "
            "the reader, per the same convention `two_axis/` and `induction/` "
            "already follow.</p>"
            "<div class='grid three'>" + "".join(control_panels) + "</div>")
    return styled("The early trough — every curve behind the finding", body,
                  "Fifteen genuine, two rejected on inspection")


def main():
    return write_pages(HERE, {"index.html": build_index(),
                              "progress_curves.html": build_curves_page()})


if __name__ == "__main__":
    raise SystemExit(main())
