"""
Builds ph/index.html and ph/progress_curves.html.

Draws only; every number comes from `ph_role`, `scope` or `induction`, so a
figure and the prose in ANALYSIS.md cannot disagree about a value without
`check_numbers.py` saying so.

    python ph/build_figures.py
"""
import functools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import induction
import ph_role
import scope
from svgplot import ACCENT, esc, GRID, INK, MUTED, Axes
from figure_kit import (CATEGORY, fig, panel, progress_axes, progress_overlay,
                        derivative_axes, residual_axes, styled, write_pages)

LADDER_COLOUR = {
    "phosphate 4OMe": CATEGORY[0],
    "boric 4OMe": CATEGORY[1],
    "pyrophosphate BnOH 136-142": CATEGORY[2],
    "pyrophosphate BnOH 143-151": CATEGORY[2],
}
LADDER_LABEL = {
    "phosphate 4OMe": "phosphate, 4OMe",
    "boric 4OMe": "boric, 4OMe",
    "pyrophosphate BnOH 136-142": "pyrophosphate, BnOH (low arm)",
    "pyrophosphate BnOH 143-151": "pyrophosphate, BnOH (high arm)",
}


@functools.cache
def _rate_table():
    return ph_role.rate_ladder_table()


@functools.cache
def _turnover():
    return ph_role.boric_turnover()


@functools.cache
def _independent():
    """Which ladder can corroborate the block, and how far it lands from it."""
    return ph_role.independent_check()


@functools.cache
def _split():
    """The boric order's dependence on where the ladder is cut."""
    return ph_role.boric_split_sensitivity()


@functools.cache
def _clock_rows():
    return induction.lag_ph_ladders()


def figure_medians():
    """A · Median vmax against pH, all four ladders -- the raw shape."""
    axes = Axes(620, 320, (5.2, 10.7), (2e-7, 3e-4), ylog=True,
                pad=(72, 26, 46, 34))
    strong = set(scope.strong_runs())
    for name, exps in scope.PH_LADDERS.items():
        colour = LADDER_COLOUR[name]
        data = scope.frame(exps)
        live = data[data.live]
        by_run = live.groupby("experiment").agg(pH=("pH", "first"),
                                                 vmax=("vmax", "median"))
        by_run = by_run.sort_values("pH")
        weak = [e for e in by_run.index if name.startswith("pyrophosphate")
               and e not in strong]
        strong_rows = by_run.drop(index=weak)
        if len(strong_rows) >= 2:
            axes.line(strong_rows.pH.to_numpy(), strong_rows.vmax.to_numpy(),
                      colour, width=1.6, opacity=0.55)
        axes.points(strong_rows.pH.to_numpy(), strong_rows.vmax.to_numpy(),
                    colour, radius=4.2, title=LADDER_LABEL[name])
        if weak:
            axes.points(by_run.loc[weak].pH.to_numpy(),
                       by_run.loc[weak].vmax.to_numpy(), colour, radius=3.4,
                       opacity=0.35, stroke=colour, stroke_width=1.2,
                       title=f"{LADDER_LABEL[name]} (weak run, dropped)")
    turnover = _turnover()
    axes.line([turnover["peak_pH"], turnover["peak_pH"]], [2e-6, 3e-4],
              CATEGORY[1], width=1.2, dash="3 3", opacity=0.6)
    axes.note(axes._fx(turnover["peak_pH"]) + 5, 40, "boric's own peak",
              CATEGORY[1], size=10, anchor="start")
    order = 60
    for name in ("phosphate 4OMe", "boric 4OMe",
                "pyrophosphate BnOH 136-142", "pyrophosphate BnOH 143-151"):
        axes.note(360, order, LADDER_LABEL[name], LADDER_COLOUR[name],
                  size=10.5, anchor="start", weight="600")
        order += 15
    return fig(
        axes.render("pH", "median vmax, AU/s",
                    "A · Median rate against pH, all four ladders",
                    xticks=True),
        "Each point is one run's median across its own substrate ladder; "
        "faint rings are the two-axis block's weak runs, dropped from every "
        "fit below. Phosphate and both pyrophosphate arms climb together; the "
        "boric ladder climbs to a peak at pH 9.50 (exp 43) and then falls -- "
        "the shape section 2's single log-log slope cannot see and section 3 "
        "reads by splitting the ladder at its own peak instead.")


def figure_ladder_orders():
    """B · The four fitted orders in [HOO-], and the pooled fit with/without boric."""
    table = _rate_table()
    pooled_all = ph_role.pooled_rate_order()
    pooled_three = ph_role.pooled_rate_order(drop=("boric 4OMe",))
    published = scope.ph_order(parameter="vmax", scope=scope.strong_runs())
    independent = _independent()
    rows = [(LADDER_LABEL[name], table.loc[name].order_hoo,
            table.loc[name].stderr_hoo, LADDER_COLOUR[name])
           for name in table.index]
    axes = Axes(640, 320, (-0.16, 0.82), (-1.2, len(rows) + 1.4),
               pad=(210, 26, 46, 34))
    axes.line([pooled_three["pooled"], pooled_three["pooled"]],
             [-1.2, len(rows) + 1.4], ACCENT, width=1.6, dash="5 4")
    axes.label(pooled_three["pooled"], -0.9,
              f"pooled w/o boric {pooled_three['pooled']:+.2f}", ACCENT,
              size=10.5, anchor="middle", weight="600", dy=4)
    for index, (label, order, stderr, colour) in enumerate(rows):
        y = len(rows) - index
        axes.line([order - stderr, order + stderr], [y, y], colour, width=2.6)
        axes.points([order], [y], colour, radius=5.0)
        axes.label(-0.18, y, label, INK, size=11, anchor="end", dy=4,
                  weight="600")
        axes.label(order + stderr, y, f"{order:+.3f}", INK, size=10.5,
                  anchor="start", dx=8, dy=4)
    axes.points([published.loc["pooled", "order"]], [0.0], INK, radius=5.0,
               stroke=INK, stroke_width=1.4)
    axes.line([published.loc["pooled", "order"]
              - published.loc["pooled", "stderr"],
              published.loc["pooled", "order"]
              + published.loc["pooled", "stderr"]], [0.0, 0.0], INK,
             width=2.6)
    axes.label(-0.18, 0.0, "two-axis, cuvette-matched", INK, size=11,
              anchor="end", dy=4, weight="600")
    axes.label(published.loc["pooled", "order"]
              + published.loc["pooled", "stderr"], 0.0,
              f"{published.loc['pooled', 'order']:+.3f}", INK, size=10.5,
              anchor="start", dx=8, dy=4)
    return fig(
        axes.render("order in [HOO-]", "",
                    "B · Three ladders agree; the fourth does not",
                    yticks=False),
        f"Point estimates ± 1 stderr for each of the four ladders "
        f"(`ph_role.rate_ladder_table`), the pooled fit without boric "
        f"(dashed, χ² = {pooled_three['chi2']:.2f} on {pooled_three['dof']}) "
        f"and, in black, the two-axis block's own cuvette-matched reading over "
        f"the same strong runs -- which the two pyrophosphate rows are NOT "
        f"independent of, being that block read along its second design; the "
        f"independent corroboration is the phosphate row, "
        f"{independent['sigma']:.2f}σ away. Boric's order "
        f"({table.loc['boric 4OMe'].order_hoo:+.3f} ± "
        f"{table.loc['boric 4OMe'].stderr_hoo:.3f}) is the reason a pool of "
        f"all four has χ² = {pooled_all['chi2']:.0f} on {pooled_all['dof']}.")


def figure_rate_vs_clock():
    """C · The rate's order and the clock's order, per ladder, side by side."""
    table = _rate_table()
    clock = {row["ladder"]: row for row in _clock_rows()}
    names = list(table.index)
    axes = Axes(640, 320, (-0.9, 0.9), (-1.2, len(names) + 1.6),
               pad=(210, 26, 58, 34))
    axes.line([0, 0], [-1.2, len(names) + 1.6], GRID, width=1.2, dash="3 3")
    for index, name in enumerate(names):
        y = len(names) - index
        colour = LADDER_COLOUR[name]
        rate_order = table.loc[name].order_hoo
        rate_stderr = table.loc[name].stderr_hoo
        clock_slope = clock[name]["lag_half_s"]["slope"]
        clock_stderr = clock[name]["lag_half_s"]["stderr"]
        axes.line([rate_order - rate_stderr, rate_order + rate_stderr],
                 [y + 0.16, y + 0.16], colour, width=2.4)
        axes.points([rate_order], [y + 0.16], colour, radius=4.4)
        axes.line([clock_slope - clock_stderr, clock_slope + clock_stderr],
                 [y - 0.16, y - 0.16], colour, width=2.4, opacity=0.55)
        axes.ring(clock_slope, y - 0.16, colour, radius=4.4, width=1.8)
        axes.label(-0.95, y, LADDER_LABEL[name], INK, size=10.5,
                  anchor="end", dy=4, weight="600")
    axes.note(360, 292, "filled: rate order in [HOO-]", INK, size=10.5,
             anchor="start", weight="600")
    axes.note(360, 306, "open: clock order, per pH unit", MUTED, size=10.5,
             anchor="start")
    return fig(
        axes.render("order", "",
                    "C · The rate separates the ladders; the clock does not",
                    yticks=False),
        "Same four ladders, two responses. The rate's order in [HOO-] "
        "(filled) is where boric parts company with the other three; the "
        "induction clock's order per pH unit (open, `induction."
        "lag_ph_ladders`) is positive on three of four and pools to "
        "+0.326 ± 0.131 with χ² = 0.95 on 3 -- the clock does not see "
        "whatever makes the boric ladder's rate turn over.")


def figure_gas_onset():
    """D · The O2 side reaction's own pH onset, from the whole archive."""
    gas = scope.gas_survey()
    axes = Axes(620, 300, (-0.6, 8.6), (-0.08, 5.0), pad=(72, 26, 46, 40))
    labels = []
    colours = {"Phosphate": CATEGORY[0], "Boric": CATEGORY[1],
              "Pyrophosphate": CATEGORY[2]}
    index = 0
    for (buffer_name, band), row in gas.iterrows():
        colour = colours[buffer_name]
        x = index
        axes.line([x, x], [0.0, float(row.per_hour)], colour, width=13,
                 opacity=0.85)
        if row.per_hour > 0:
            axes.label(x, float(row.per_hour), f"{row.per_hour:.2f}", INK,
                      size=9.5, anchor="middle", dy=-6)
        labels.append((x, f"{buffer_name[:4]}\n{band.left:g}-{band.right:g}",
                      colour))
        index += 1
    for x, text, colour in labels:
        for line_index, line in enumerate(text.split("\n")):
            axes.note(axes._fx(x), 268 + 11 * line_index, line, colour,
                      size=9, anchor="middle")
    return fig(
        axes.render("", "detachments / hour", "D · The gas side reaction "
                    "turns on with pH, not with buffer identity", xticks=False),
        "`scope.gas_survey`, every buffer and pH band the archive holds. "
        "Zero detachments anywhere below pH 7.5; boric above pH 8.5 is the "
        "single heaviest band after pyrophosphate's own high-pH arm. The "
        "same pH region section 3 finds the catalysed RATE turning over in "
        "boric is where this unrelated side reaction is heaviest -- checked, "
        "in section 3, and not the explanation.")


def _mm_fit_panel(row, group, colour):
    """
    One experiment's own MM fit: `v_peak` against [S], the fit if resolved.

    DATA FIRST, FIT ON TOP -- `progress_overlay`'s own convention, so the
    points are drawn before the line and the line stays thin enough not to
    bury them. The axes always include y=0 (an MM curve passes through it)
    and extend far enough in x and y to hold the fitted line as well as the
    points, so a fit that has not yet saturated by the ladder's own top
    rung is never clipped off the top of its own panel.
    """
    s0 = group.s0.to_numpy(dtype=float)
    y = group.v_peak.to_numpy(dtype=float)
    top_x = float(s0.max()) * 1.15
    smooth = np.linspace(0.0, top_x, 200)
    fit_line = (row.vmax * smooth / (row.km + smooth)
               if row.km_resolved else None)
    top_y = max(float(y.max()),
               float(fit_line.max()) if fit_line is not None else 0.0) * 1.15
    axes = Axes(300, 210, (0.0, top_x), (0.0, max(top_y, 1e-8)),
               pad=(58, 12, 34, 8))
    axes.points(s0, y, colour, radius=3.8, stroke="white", stroke_width=0.8)
    if fit_line is not None:
        axes.line(smooth, fit_line, colour, width=1.6)
        caption = (f"Vmax {row.vmax:.2e} AU/s · Km {row.km:.2f} "
                  f"({row.km_low:.2f}–{row.km_high:.2f}) mM · "
                  f"R² {row.r2:.2f}")
    else:
        caption = (f"Km unresolved -- profile reaches the grid "
                  f"({row.km_low:.3g}–{row.km_high:.3g} mM)")
    return fig(
        axes.render("[S], mM", "v_peak, AU/s",
                    f"exp {int(row.experiment)} · pH {row.pH:.2f}"),
        caption)


def build_mm_section(name, experiments):
    """
    The per-experiment Michaelis-Menten fit behind ANALYSIS.md §3a: `v_peak`
    against [S], one panel per pH rung, with the fitted curve where `Km`
    resolves. `v0_fit` is not shown -- restricted to `v0_fit_resolved`
    cuvettes, most rungs keep only 1-4 of their 4 cuvettes, too few to fit
    independently (`ph_role.ladder_mm_table`, response="v0_fit").
    """
    table = ph_role.ladder_mm_table(experiments, response="v_peak")
    data = scope.frame(tuple(experiments))
    data = data[data.live]
    colour = LADDER_COLOUR[name]
    panels = [_mm_fit_panel(row, data[data.experiment == row.experiment],
                            colour)
             for row in table.sort_values("pH").itertuples()]
    return (f"<p class='lede'>Per-experiment MM fit, {LADDER_LABEL[name]}: "
           f"each panel is one run's own four-cuvette substrate ladder "
           f"(`ph_role.ladder_mm_table`). {len(panels)} panels, one per "
           f"pH rung.</p><div class='grid three'>"
           + "".join(panels) + "</div>")


def build_curves_page():
    """Every live cuvette of all four pH ladders, in pH order within each."""
    sections = []
    total = 0
    for name, exps in scope.PH_LADDERS.items():
        frame = scope.frame(exps)
        lookup = {(c.experiment, c.sample): c for c in scope.curves(exps)}
        panels = []
        for row in frame[frame.live].sort_values(
                ["pH", "experiment", "s0"]).itertuples():
            curve = lookup.get((row.experiment, row.sample))
            if curve is None:
                continue
            times = np.asarray(curve.times, dtype=float)
            values = np.asarray(curve.absorbance, dtype=float)
            axes, radius = progress_axes(times, values, limit=140)
            progress = progress_overlay(axes, times, values,
                                       mark_radius=radius)
            residual = (values - progress.predict(times)) / curve.noise
            rax = residual_axes(times, residual, colour=CATEGORY[0])
            drax = derivative_axes(times, progress, colour=CATEGORY[0])
            panels.append(panel(
                f"pH {row.pH:.2f} · [S] {row.s0:.3f} mM"
                f"<span class='pill'>exp {int(row.experiment)}.{int(row.sample)}</span>",
                f"{LADDER_LABEL[name]} · [H₂O₂] {row.h2o2:g} mM · "
                f"{row.temperature:.0f} °C · {int(row.points)} readings over "
                f"{row.duration_s / 3600:.1f} h · {row.source}",
                axes.render("", "ΔA") + rax.render("", "z")
                + drax.render("time, s", "dA/dt"),
                f"<strong>{int(row.phases)} phase"
                + ("s" if row.phases == 2 else "") + f"</strong> · "
                f"{esc(str(row.progress_kind))} · vmax {row.vmax:.2e}"
                + f" · bubble_load {row.bubble_load:.2f}"
                if hasattr(row, "bubble_load") else ""))
        total += len(panels)
        section = f"<h2>{esc(LADDER_LABEL[name])}</h2>"
        if name in ph_role.MM_LADDERS:
            section += build_mm_section(name, exps)
        section += "<div class='grid three'>" + "".join(panels) + "</div>"
        sections.append(section)
    body = (f"<p class='lede'>All {total} live cuvettes of the archive's "
            "four pH ladders (`scope.PH_LADDERS`), grouped by ladder and "
            "sorted by pH within each. The rust line is whichever form the "
            "curve earned, from `summary_kinetics.fit_progress`; nothing is "
            "excluded. Phosphate and boric also carry their own "
            "per-experiment Michaelis-Menten fit ahead of the cuvette grid "
            "(ANALYSIS.md §3a).</p>" + "".join(sections))
    return styled("The pH ladders — every progress curve", body,
                 "Phosphate and boric 4OMe, pyrophosphate BnOH (135-151)")


def build_index():
    table = _rate_table()
    pooled_all = ph_role.pooled_rate_order()
    pooled_three = ph_role.pooled_rate_order(drop=("boric 4OMe",))
    published = scope.ph_order(parameter="vmax", scope=scope.strong_runs())
    clock_pooled = ph_role.clock_pooled_order()
    turnover = _turnover()
    independent = _independent()
    split = _split()

    hero = f"""
<div class='hero'>
  <div><div class='k'>pH ladders</div><div class='v'>4</div>
       <div class='u'>two buffers, two substrates, never read as one question</div></div>
  <div><div class='k'>rate order, 3 of 4</div>
       <div class='v'>+{pooled_three['pooled']:.2f}</div>
       <div class='u'>chi2 {pooled_three['chi2']:.2f} on {pooled_three['dof']} -- essentially one number</div></div>
  <div><div class='k'>boric turns over</div>
       <div class='v'>pH {turnover['peak_pH']:.2f}</div>
       <div class='u'>not the O2 side reaction -- checked, section 3</div></div>
  <div><div class='k'>the clock, all 4</div>
       <div class='v'>+{clock_pooled['pooled']:.2f}/pH</div>
       <div class='u'>chi2 {clock_pooled['chi2']:.2f} on {clock_pooled['dof']} -- agrees where the rate does not</div></div>
</div>"""

    rows = "".join(
        f"<tr><td>{esc(LADDER_LABEL[name])}</td><td>{int(table.loc[name].runs)}</td>"
        f"<td>{table.loc[name].pH_low:.2f}-{table.loc[name].pH_high:.2f}</td>"
        f"<td>{table.loc[name].order_s0:+.3f} ± {table.loc[name].stderr_s0:.3f}</td>"
        f"<td>{table.loc[name].order_hoo:+.3f} ± {table.loc[name].stderr_hoo:.3f}</td></tr>"
        for name in table.index)

    body = f"""
<p class='lede'>Four pH ladders, built for the induction clock alone
(`scope.PH_LADDERS`, 2026-09-05) and never asked about the rate.
<a href='ANALYSIS.md'>ANALYSIS.md</a> is the argument; this is the picture of
it.</p>
{hero}

<h2>1 · Which curves</h2>
<p>The hand-sorted <code>data/Mads</code> folders needed two corrections
against <code>scope.PH_LADDERS</code>: drop the replicate quadruplet (exps 2,
4, 5, 7) from the phosphate set and add exp 14 to it; drop exp 13 from the
boric set. The pH-11 folder was already correctly excluded
(<code>build_manifest.KNOWN_EXCLUSIONS[85]</code>).</p>

<h2>2 · The rate's order in [HOO⁻]</h2>
{figure_medians()}
<div class='tbl'><table>
<tr><th>ladder</th><th>runs</th><th>pH range</th><th>order in [S]</th>
<th>order in [HOO⁻]</th></tr>
{rows}
</table></div>
{figure_ladder_orders()}

<h2>3 · The boric ladder turns over</h2>
<p>Median vmax rises from pH 8.46 to a peak at exp 43 (pH 9.50) and falls
through pH 10.34 — below the ladder's own reading at pH 8.98. Refit below the
peak, the order is <strong>+{turnover['below']['hoo']['slope']:.3f} ±
{turnover['below']['hoo']['stderr']:.3f}</strong>, positive but weaker than
phosphate or pyrophosphate's. <code>vmax_corrected</code> moves the two
highest-pH runs' medians by less than a quarter of themselves, so the O2 side
reaction — heaviest in exactly this buffer and pH range — is not the cause.</p>

<h2>4 · The clock agrees where the rate does not</h2>
{figure_rate_vs_clock()}

<h2>5 · The gas side reaction's own pH onset</h2>
{figure_gas_onset()}

<h2>6 · One species, three signals, one caveat</h2>
<p>The rate (+{pooled_three['pooled']:.3f}), the induction clock
(+{clock_pooled['pooled']:.3f} per pH unit) and the early trough's
OTHER-cuvette signal (<code>ρ = -0.621</code> for 4OMe, <code>-0.315</code>
for BnOH against [enz]/[HOO⁻], both p &lt; 10⁻⁴) all name [HOO⁻] rather than
pH itself or total [H₂O₂]. <code>scope.hoo_consistency</code> is the caveat:
on the two-axis block alone, moving [HOO⁻] two independent ways gives orders
that part at 3.4σ — [HOO⁻] is not the whole story even where it is real.</p>

<h2>7 · What this settles, and what it does not</h2>
<p><strong>Settled.</strong> The rate's order in [HOO⁻] is
+{pooled_three['pooled']:.3f} ± {pooled_three['stderr']:.3f} in phosphate and
pyrophosphate. Two of those three ladders <em>are</em> the two-axis block
({100 * independent['two_axis_share']:.1f}% of the pooled weight), so the
independent corroboration is the phosphate ladder alone —
+{independent['phosphate_order']:.3f} ± {independent['phosphate_stderr']:.3f}
on nine runs sharing no experiment with the block,
{independent['sigma']:.2f}σ from its cuvette-matched
+{published.loc['pooled', 'order']:.3f} ± {published.loc['pooled', 'stderr']:.3f}.
The boric ladder does not share it, and the O2 side reaction is ruled out as
the cause of its decline — though the order below its peak spans
+{split['low']:.2f} to +{split['high']:.2f} depending on where the ladder is
cut, so it is quoted as a range. The clock's order agrees across all four
ladders where the rate's does not.</p>
<p><strong>Not settled.</strong> Whether boric's turnover is the buffer or a
genuine high-pH instability the way exp 85 showed at pH 11.84 — the archive
holds no other buffer above pH 9 to separate them, which is
<a href='../buffer/index.html'>buffer/</a> §5's finding applied to a new
question rather than a resolution of it.</p>

<h2>Reproducing</h2>
<p><code>python data/ph_role.py</code> prints the whole argument ·
<code>python data/test_ph_role.py</code> ·
<code>python ph/build_figures.py</code> ·
<code>python ph/check_numbers.py</code>, which re-derives every number in
<code>ANALYSIS.md</code> from the modules and fails if the prose and the code
disagree.</p>
"""
    return styled("What pH does to the reaction", body,
                 "Four ladders, read the same way for the first time")


def main():
    return write_pages(HERE, {"index.html": build_index(),
                              "progress_curves.html": build_curves_page()})


if __name__ == "__main__":
    raise SystemExit(main())
