"""
Builds the figures for background_reaction/.

Draws only. Every measurement comes from `scope`, `curve_metrics` or
`summary_kinetics` -- see CLAUDE.md, "Do not re-derive measurements". The only
transforms applied here are for display: centring a run's points on its own
geometric mean so that a within-experiment slope can be seen on one axis.

    python background_reaction/build_figures.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, HERE)

import scope
from curve_metrics import ACCELERATION_SIGMA, INITIAL_WINDOW
from fit_dataset import source_floor
from summary_kinetics import BURST_V0_HALFWIDTH, fit_burst_bounded
from svgplot import ACCENT, GRID, INK, MUTED, PALETTE, Axes, esc, page

# The estimator the headline numbers are quoted on, and the two reported
# beside it. v0_quad answers the "the 20% window is arbitrary" objection: it is
# the slope at t = 0 of a quadratic through EVERY point, so no window is
# chosen anywhere, and being linear in its parameters it is always identified.
HEADLINE = "v0_quad"
ESTIMATORS = ("v0_quad", "vmax", "v0", "v0_whole")
ESTIMATOR_LABEL = {
    "v0_quad": "v0 quadratic (whole curve, no window)",
    "vmax": "vmax (steepest 20% block)",
    "v0": "v0 (first 20% window)",
    "v0_whole": "slope of a straight line through the whole curve",
}
FIXED_COLOUR, TITRATION_COLOUR = PALETTE[0], ACCENT


def _geometric_centre(values):
    values = np.asarray(values, dtype=float)
    return float(np.exp(np.mean(np.log(values))))


def _within_run_points(frame, parameter):
    """
    Each cuvette as (s0 / run's geometric mean s0, rate / run's geometric mean).

    This is the WITHIN transformation the orders are fitted with, made visible:
    a per-experiment offset is exactly a per-experiment division on a log axis,
    so after it every run overlays and the common slope is the thing you see.
    Runs are not being pooled -- they are being centred.
    """
    out = []
    live = frame[frame.live & (frame[parameter] > 0)]
    for experiment, group in live.groupby("experiment"):
        cx = _geometric_centre(group.s0)
        cy = _geometric_centre(group[parameter])
        out.append((int(experiment),
                    group.s0.to_numpy(dtype=float) / cx,
                    group[parameter].to_numpy(dtype=float) / cy,
                    group.buf.to_numpy(dtype=float)))
    return out


def figure_separation(parameter=HEADLINE):
    """Panel A: the substrate order changes sign with the buffer design."""
    fixed = scope.frame(scope.BUFFER_FIXED)
    titration = scope.frame(scope.BUFFER_CONFOUNDED)
    result = scope.buffer_dependence(parameter=parameter)

    axes = Axes(560, 380, (0.14, 7.0), (0.28, 3.4), xlog=True, ylog=True)
    for series, colour in ((_within_run_points(fixed, parameter), FIXED_COLOUR),
                           (_within_run_points(titration, parameter),
                            TITRATION_COLOUR)):
        for experiment, x, y, buf in series:
            for xi, yi, bi in zip(x, y, buf):
                axes.points([xi], [yi], colour, radius=3.6,
                            title=f"exp {experiment}: [buf] {bi:.3g} mM")
    span = np.array([0.16, 6.4])
    for slope, colour, key in ((result["order_s0_fixed"], FIXED_COLOUR, "fixed"),
                               (result["order_s0_titration"], TITRATION_COLOUR,
                                "titration")):
        axes.line(span, span ** slope, colour, width=2.4)
    axes.hline(1.0, GRID, "3 3")
    axes.note(300, 44, f"[buf] held: slope {result['order_s0_fixed']:+.3f}"
                       f" ± {result['stderr_s0_fixed']:.3f}",
              FIXED_COLOUR, 11.5, weight="600")
    axes.note(300, 61, f"[buf] falling: slope {result['order_s0_titration']:+.3f}"
                       f" ± {result['stderr_s0_titration']:.3f}",
              TITRATION_COLOUR, 11.5, weight="600")
    axes.note(300, 78, "same reaction — the order changes SIGN", MUTED, 11)
    return axes.render("[BnOH] / run's geometric mean",
                       "rate / run's geometric mean",
                       "A · substrate order, within runs")


def figure_coupling():
    """Panel B: how [buf] tracks [sub] inside the titration runs."""
    titration = scope.frame(scope.BUFFER_CONFOUNDED)
    result = scope.buffer_dependence(parameter=HEADLINE)
    axes = Axes(560, 380, (0.9, 11.0), (20.0, 100.0), xlog=True, ylog=True)
    for index, (experiment, group) in enumerate(
            titration.groupby("experiment")):
        colour = PALETTE[index % len(PALETTE)]
        order = group.sort_values("s0")
        axes.line(order.s0.to_numpy(dtype=float),
                  order.buf.to_numpy(dtype=float), colour, width=1.5, dash="5 3")
        for _, row in order.iterrows():
            axes.points([row.s0], [row.buf], colour,
                        radius=4.0 if row.live else 2.6,
                        opacity=1.0 if row.live else 0.42,
                        title=f"exp {experiment} s{int(row['sample'])}"
                              f"{'' if row.live else ' (not live)'}")
        axes.note(400, 44 + 17 * index, f"exp {experiment}", colour, 11.5,
                  weight="600")
    axes.note(150, 44, f"g = dlog[buf]/dlog[sub] = {result['coupling']:+.3f}",
              INK, 11.5, weight="600")
    axes.note(150, 61, "substrate volume displaced buffer volume", MUTED, 10.5)
    return axes.render("[BnOH] / mM", "[buffer] / mM",
                       "B · the confound: buffer falls as substrate rises")


def figure_buffer_order():
    """Panel C: the recovered buffer order, every estimator, both blocks."""
    rows = []
    for estimator in ESTIMATORS:
        result = scope.buffer_dependence(parameter=estimator)
        cross = scope.background_orders(scope.FREE_4OME_40C,
                                        terms=("s0", "buf"), parameter=estimator)
        rows.append((estimator, result["order_buf"], result["stderr_buf"],
                     cross["order_buf"], cross["stderr_buf"]))

    axes = Axes(560, 380, (-0.6, len(rows) - 0.4), (-6.0, 3.0),
                pad=(62, 16, 74, 26))
    axes.hline(0.0, MUTED, "3 3", 1.4)
    axes.hline(1.0, "#3f8a5a", "6 4", 1.6)
    for index, (name, bnoh, bnoh_e, ome, ome_e) in enumerate(rows):
        axes.errorbar(index - 0.13, bnoh, bnoh - bnoh_e, bnoh + bnoh_e,
                      PALETTE[0])
        axes.errorbar(index + 0.13, ome, ome - ome_e, ome + ome_e, PALETTE[2])
        axes.note(axes._fx(index), axes.height - 52,
                  name.replace("v0_", "v0 "), MUTED, 10.5, anchor="middle")
    axes.note(axes._fx(len(rows) - 1) - 4, axes._fy(1.0) - 7,
              "first order", "#3f8a5a", 10.5, anchor="end")
    axes.note(70, 44, "BnOH 25 °C", PALETTE[0], 11.5, weight="600")
    axes.note(70, 61, "4OMe-BnOH 40 °C (independent)", PALETTE[2], 11.5,
              weight="600")
    return axes.render("", "order in [buffer]",
                       "C · buffer order, by estimator", xticks=False)


# --- the progress curves, with every fit drawn -----------------------------
WINDOW_COLOUR = PALETTE[0]     # the 20% window line
QUAD_COLOUR = "#3f8a5a"        # the whole-curve quadratic (headline)
BURST_COLOUR = "#12856a"       # the burst/lag form, as curve_dossier draws it


def curve_panel(curve, width=330, height=230):
    """One cuvette: the data, all three fits, and the numbers behind them."""
    from curve_metrics import initial_rate, quadratic_rate, whole_slope

    times = np.asarray(curve.times, dtype=float)
    times = times - times[0]
    values = np.asarray(curve.absorbance, dtype=float)
    floor = source_floor(curve.source)

    v0, v0_se, _ = initial_rate(curve.times, curve.absorbance, floor=floor)
    quad, quad_se, curvature = quadratic_rate(curve.times, curve.absorbance,
                                              floor=floor)
    whole, whole_se = whole_slope(curve.times, curve.absorbance, floor=floor)
    burst = fit_burst_bounded(curve.times, curve.absorbance)

    lo, hi = float(values.min()), float(values.max())
    margin = max((hi - lo) * 0.16, 4 * curve.noise)
    axes = Axes(width, height, (0.0, float(times[-1])), (lo - margin, hi + margin),
                pad=(58, 10, 34, 22))

    # the data
    axes.line(times, values, MUTED, width=1.0, opacity=0.55)
    step = max(1, len(times) // 90)
    axes.points(times[::step], values[::step], INK, radius=1.5, opacity=0.7)

    # the window the 20% line was fitted through, shaded
    cut = float(times[-1]) * INITIAL_WINDOW
    axes.parts.append(
        f"<rect x='{axes._fx(0):.1f}' y='{axes.top}' "
        f"width='{axes._fx(cut) - axes._fx(0):.1f}' "
        f"height='{height - axes.top - axes.bottom:.1f}' fill='{WINDOW_COLOUR}' "
        f"fill-opacity='0.09'/>")

    base = float(np.median(values[:max(1, min(5, len(values) // 10))]))
    grid = np.linspace(0, float(times[-1]), 120)
    axes.line(np.array([0, cut]), base + v0 * np.array([0, cut]),
              WINDOW_COLOUR, width=2.0)
    axes.line(grid, base + quad * grid
              + (quad_se * 0 + _quad_curvature(times, values)) * grid ** 2,
              QUAD_COLOUR, width=1.9)
    if np.isfinite(burst.tau):
        axes.line(grid, burst.predict(grid), BURST_COLOUR, width=1.6, dash="5 3")
        # the v0 profile interval, drawn as the fan of slopes it allows
        reach = float(times[-1]) * 0.22
        axes.band(np.array([0.0, reach]),
                  base + burst.v0_low * np.array([0.0, reach]),
                  base + burst.v0_high * np.array([0.0, reach]),
                  BURST_COLOUR, opacity=0.22)

    flag = "bounded" if burst.bounded else "UNBOUNDED"
    head = (f"exp {curve.experiment} s{curve.sample} · "
            f"[BnOH] {curve.conditions.s0:.3g} · [buf] {curve.buf:.4g} mM")
    axes.note(axes.left, 12, head, INK, 10.6, weight="600")
    axes.note(axes.left, height - 4,
              f"{curve.source} · noise {curve.noise:.1e} · curvature t "
              f"{curvature:+.1f}", MUTED, 9.4)
    rows = [
        (f"v0 quad {quad:.3e}", QUAD_COLOUR),
        (f"v0 win  {v0:.3e}", WINDOW_COLOUR),
        (f"burst   {burst.v0:.3e} ({flag})", BURST_COLOUR),
        (f"whole   {whole:.3e}", MUTED),
    ]
    for index, (text, colour) in enumerate(rows):
        axes.parts.append(
            f"<text x='{width - 12}' y='{26 + 12 * index}' font-size='9.3' "
            f"fill='{colour}' text-anchor='end' "
            f"font-family='ui-monospace,Menlo,monospace'>{esc(text)}</text>")
    return axes.render("time / s", "ΔA")


def _quad_curvature(times, values):
    """The quadratic's t^2 coefficient, for drawing its curve only."""
    design = np.column_stack([np.ones(len(times)), times, times ** 2])
    beta, *_ = np.linalg.lstsq(design, values, rcond=None)
    return float(beta[2])


def figure_curvature():
    """Deceleration against conversion: it is not substrate depletion."""
    frame = scope.frame(scope.FREE_BNOH_ALL)
    live = frame[frame.live]
    axes = Axes(560, 380, (0.03, 12.0), (-90.0, 15.0), xlog=True,
                pad=(62, 16, 46, 22))
    axes.hline(0.0, MUTED, "3 3", 1.3)
    axes.hline(-3.0, ACCENT, "5 3", 1.4)
    for index, (experiment, group) in enumerate(live.groupby("experiment")):
        colour = PALETTE[index % len(PALETTE)]
        for _, row in group.iterrows():
            axes.points([100 * row.conversion], [row.curvature_t], colour,
                        radius=4.0,
                        title=f"exp {experiment} s{int(row['sample'])}")
        axes.note(430, 210 + 16 * index, f"exp {experiment}", colour, 11)
    axes.note(axes._fx(0.034), axes._fy(-3.0) - 6,
              "|t| = 3: curvature is real below this line", ACCENT, 10.2)
    axes.note(90, 44, "every curve is under 8% converted,", INK, 11.5,
              weight="600")
    axes.note(90, 61, "so substrate depletion cannot bend them this far",
              INK, 11.5, weight="600")
    return axes.render("conversion at end of run / %",
                       "curvature t-statistic (negative = decelerating)",
                       "the background decelerates, and depletion does not explain it")


def figure_acceleration():
    """The in-scope curves accelerate; the background does not."""
    sets = [("enzyme-free, [buf] fixed\n(65,67,69,70)", scope.BUFFER_FIXED,
             PALETTE[0]),
            ("enzyme-free titrations\n(3,6)", scope.BUFFER_CONFOUNDED, ACCENT),
            ("in-scope catalysed\n(135-151)", scope.PRIMARY_SCOPE, PALETTE[2])]
    axes = Axes(560, 380, (-0.6, 2.6), (-12.0, 22.0), pad=(62, 16, 62, 26))
    axes.hline(ACCELERATION_SIGMA, ACCENT, "5 3", 1.6)
    axes.hline(0.0, MUTED, "3 3", 1.2)
    for index, (label, which, colour) in enumerate(sets):
        frame = scope.frame(which)
        live = frame[frame.live]
        z = np.clip(live.accel_z.to_numpy(dtype=float), -12, 22)
        jitter = np.linspace(-0.22, 0.22, len(z))
        for xi, zi in zip(jitter, z):
            axes.points([index + xi], [zi], colour, radius=3.0, opacity=0.72)
        share = int((live.accel_z > ACCELERATION_SIGMA).sum())
        axes.note(axes._fx(index), axes.height - 44,
                  label.split("\n")[0], MUTED, 10.2, anchor="middle")
        axes.note(axes._fx(index), axes.height - 32,
                  label.split("\n")[1], MUTED, 10.2, anchor="middle")
        axes.note(axes._fx(index), axes._fy(20.4),
                  f"{share}/{len(live)}", colour, 12.5, anchor="middle",
                  weight="700")
    axes.note(axes._fx(2.55), axes._fy(ACCELERATION_SIGMA) - 6,
              f"z = {ACCELERATION_SIGMA:.0f}: accelerating above this",
              ACCENT, 10.2, anchor="end")
    return axes.render("", "acceleration z (steeper later than at the start)",
                       "autocatalysis is in the catalysed increment, not the background",
                       xticks=False)


def figure_rate_law(parameter=HEADLINE):
    """The four orders, each from the design that can actually measure it."""
    fixed = scope.background_orders(scope.BUFFER_FIXED, terms=("s0",),
                                    within=True, parameter=parameter)
    pooled = scope.background_orders(terms=("s0", "h2o2", "hoo"),
                                     parameter=parameter)
    buffer = scope.buffer_dependence(parameter=parameter)
    entries = [
        ("[BnOH]", fixed["order_s0"], fixed["stderr_s0"],
         f"within runs, [buf] fixed, n={fixed['n']}"),
        ("[H2O2]", pooled["order_h2o2"], pooled["stderr_h2o2"],
         f"pooled, n={pooled['n']}"),
        ("[HOO-]", pooled["order_hoo"], pooled["stderr_hoo"],
         f"pooled, n={pooled['n']}"),
        ("[buffer]", buffer["order_buf"], buffer["stderr_buf"],
         "from the two designs' disagreement"),
    ]
    # Wide enough for the largest bar the headline estimator produces: on
    # v0_quad the [H2O2] order reaches +1.57 +/- 0.45, which a 2.0 limit clips.
    reach = max(abs(v) + e for _, v, e, _ in entries) * 1.12
    axes = Axes(560, 330, (min(-0.35, -reach * 0.2), max(2.0, reach)),
                (-0.6, len(entries) - 0.4), pad=(84, 130, 46, 26))
    axes.parts.append(
        f"<path d='M{axes._fx(0):.1f},{axes.top} L{axes._fx(0):.1f},"
        f"{axes.height - axes.bottom}' stroke='{MUTED}' stroke-width='1.2' "
        f"stroke-dasharray='3 3' fill='none'/>")
    axes.parts.append(
        f"<path d='M{axes._fx(1):.1f},{axes.top} L{axes._fx(1):.1f},"
        f"{axes.height - axes.bottom}' stroke='#3f8a5a' stroke-width='1.4' "
        f"stroke-dasharray='6 4' fill='none'/>")
    for index, (name, value, error, source) in enumerate(entries):
        y = len(entries) - 1 - index
        colour = PALETTE[index % len(PALETTE)]
        axes.parts.append(
            f"<path d='M{axes._fx(value - error):.1f},{axes._fy(y):.1f} "
            f"L{axes._fx(value + error):.1f},{axes._fy(y):.1f}' "
            f"stroke='{colour}' stroke-width='3' stroke-linecap='round'/>")
        axes.points([value], [y], colour, radius=5)
        axes.note(axes.left - 8, axes._fy(y) + 4, name, INK, 11.5,
                  anchor="end", weight="600")
        axes.note(axes.width - axes.right + 8, axes._fy(y) + 0.5,
                  f"{value:+.2f} ± {error:.2f}", colour, 11, weight="600")
        axes.note(axes.width - axes.right + 8, axes._fy(y) + 12, source,
                  MUTED, 9.4)
    axes.note(axes._fx(1.0) + 4, axes.top + 12, "first order", "#3f8a5a", 10.2)
    return axes.render("apparent reaction order", "",
                       f"the enzyme-free BnOH rate law · {ESTIMATOR_LABEL[parameter]}")


# --- pages -----------------------------------------------------------------
def _table(headers, rows, highlight=()):
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for index, row in enumerate(rows):
        cls = " class='hl'" if index in highlight else ""
        body.append(f"<tr{cls}>"
                    + "".join(f"<td>{c}</td>" for c in row) + "</tr>")
    return (f"<div class='scroll'><table><thead><tr>{head}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>")


def _fig(svg, caption):
    return f"<div class='fig'>{svg}<div class='cap'>{caption}</div></div>"


def build_index():
    buffer = scope.buffer_dependence(parameter=HEADLINE)
    rows = []
    for estimator in ESTIMATORS:
        result = scope.buffer_dependence(parameter=estimator)
        cross = scope.background_orders(scope.FREE_4OME_40C,
                                        terms=("s0", "buf"), parameter=estimator)
        rows.append([
            ESTIMATOR_LABEL[estimator],
            f"{result['order_s0_fixed']:+.3f} ± {result['stderr_s0_fixed']:.3f}",
            f"{result['order_s0_titration']:+.3f} ± {result['stderr_s0_titration']:.3f}",
            f"<b>{result['order_buf']:+.2f} ± {result['stderr_buf']:.2f}</b>",
            f"{cross['order_buf']:+.2f} ± {cross['stderr_buf']:.2f}",
        ])

    body = f"""
<p class='lede'>The uncatalysed oxidation of benzyl alcohol by H<sub>2</sub>O<sub>2</sub>,
characterised so that the catalysed progress curves of exps 135–151 — which were
recorded <em>against</em> it — can be read. The controlling question was how much of the
enzyme-free rate belongs to the buffer rather than to the substrate.</p>

<h2>The design problem</h2>
<p>No enzyme-free run in the archive varies <code>[buf]</code> at constant <code>[sub]</code>,
in any block. The five that once did — exps 32, 34–37 — were ruled <em>catalysed</em> on
2026-08-31 on their reference-channel layout. What the archive has instead is two BnOH
designs that disagree, and the disagreement is what carries the buffer order.</p>
{_fig(figure_separation(), "<b>A.</b> Every cuvette, divided by its own run's geometric "
      "mean in both axes — which is exactly what a per-experiment offset does on a log "
      "scale, so the common within-run slope is what you see. Blue: the runs that hold "
      "[buf] at 85–87.5 mM. Orange: the runs where [buf] falls as [sub] rises. The same "
      "reaction reads a positive substrate order in one design and a negative one in the "
      "other. Hover a point for its buffer concentration.")}
{_fig(figure_coupling(), "<b>B.</b> Why. Substrate was pipetted in and buffer volume shrank "
      "to compensate, so the two move together with g = dlog[buf]/dlog[sub] fitted within "
      "runs. Faded points are cuvettes with no live signal, which the fits exclude.")}
<p>If the rate goes as [S]<sup>a</sup>[buf]<sup>d</sup> and within the titrations
log[buf] = g·log[sub] + constant, then fitting the titrations without a buffer term
returns a′ = a + d·g, so <b>d = (a′ − a) / g</b>. Both a and a′ are within-run slopes,
so pH, [H<sub>2</sub>O<sub>2</sub>], cell and day are absorbed on both sides and nothing
asks one regression to separate [buf] from [sub].</p>

<h2>The answer: first order in buffer</h2>
{_fig(figure_buffer_order(), "<b>C.</b> The recovered buffer order under four different "
      "rate estimators, against an independent block. Green dashed line: first order.")}
{_table(["rate estimator", "a · [buf] fixed", "a′ · [buf] falling",
         "order in [buf] · BnOH", "cross-check · 4OMe-BnOH 40 °C"], rows,
        highlight=(0,))}
<p class='warn'>The last row is the instructive failure. A straight line through the whole
curve uses every point and chooses no window, yet it is the only estimator that gets the
answer badly wrong — because 22 of 27 curves genuinely decelerate, so it reads the average
rate, and the bias scales with run length and curvature, which differ between cuvettes.
Using the whole curve helps only if the fitted form is allowed to bend.</p>

<h2>The rate law</h2>
<p>Each order is taken from the design that can measure it: the substrate order from the
constant-buffer runs only, the peroxide and [HOO⁻] orders pooled across all six runs, the
buffer order from the two designs' disagreement. Shown under both estimators, because
neither dominates: <code>v0_quad</code> is the tighter on substrate and buffer, but the
pooled peroxide fit is better conditioned on <code>vmax</code> (R² 0.961 against 0.909),
since an extrapolated initial rate adds variance where a block slope does not.</p>
<div class='grid two'>
{_fig(figure_rate_law('v0_quad'), "Headline estimator — no window anywhere.")}
{_fig(figure_rate_law('vmax'), "The estimator the in-scope block was measured with.")}
</div>

<h2>What the background is not</h2>
{_fig(figure_curvature(), "The curves bend, and it is not substrate being consumed: every "
      "one is under 8% converted, yet most show curvature far beyond |t| = 3. Something "
      "else decays during these runs — peroxide, or the cell.")}
{_fig(figure_acceleration(), "Both sets compared here are entirely .rre, so this is not "
      "the variance-floor artefact of DATA_VERIFICATION.md 2026-09-01. The enzyme-free "
      "curves do not accelerate; several actively decelerate.")}

<h2>What this means for exps 135–151</h2>
<p><b>The confound does not reach them.</b> <code>[buf]</code> is 75.013 mM in all 119
in-scope curves across all 17 runs — zero variation — so no buffer effect can enter their
substrate order.</p>
<p><b>The background is already subtracted.</b> Every run is double-beam and the catalysed
sheets' reference omits only the enzyme, so an in-scope curve is a catalytic increment
taken against a background at identical buffer, substrate, peroxide and pH.</p>
<p><b>What is still missing is the amplitude.</b> There are 0 enzyme-free curves in the
127-curve pyrophosphate cell, so the absolute size of what was subtracted cannot be
recovered without importing it from another buffer — which MECHANISM.md's buffer section
forbids. That remains FITTING.md's F7.</p>

<p style='margin-top:34px'><a href='progress_curves.html'>→ every progress curve, with all
three fits drawn</a></p>
"""
    return page("The background reaction: BnOH + H₂O₂ without catalyst", body)


def build_curves_page():
    panels = []
    for curve in scope.curves(scope.FREE_BNOH_ALL):
        panels.append(f"<div class='fig'>{curve_panel(curve)}</div>")
    bounded = sum(fit_burst_bounded(c.times, c.absorbance).bounded
                  for c in scope.curves(scope.FREE_BNOH_ALL))
    body = f"""
<p class='lede'>All 27 enzyme-free BnOH cuvettes, each with every rate estimator drawn on
it. Nothing here is a kinetic-model fit: the mechanism has never been fitted to these
curves. These are the empirical rate measurements the analysis rests on.</p>

<h2>What is drawn</h2>
<div class='key'>
  <span><i class='sw' style='background:{QUAD_COLOUR}'></i>quadratic through every point — v0 is its slope at t = 0 <b>(headline)</b></span>
  <span><i class='sw' style='background:{WINDOW_COLOUR}'></i>least-squares line over the first {INITIAL_WINDOW:.0%} (shaded)</span>
  <span><i class='sw' style='background:{BURST_COLOUR}'></i>burst/lag form, B ≤ 0, dashed; shaded fan = v0 profile interval</span>
  <span><i class='sw' style='background:{MUTED}'></i>straight line through the whole curve</span>
</div>
<p>The <b>quadratic</b> is A = c + v<sub>0</sub>t + at². It chooses no window, uses every
point, and is linear in its parameters, so v<sub>0</sub> is always identified. Its
<code>curvature t</code> is a/se(a): |t| &gt; 3 means the curve measurably bends.</p>
<p>The <b>burst/lag</b> form is A = c + v<sub>ss</sub>t − B(1 − e<sup>−t/τ</sup>), with
v<sub>0</sub> = v<sub>ss</sub> − B/τ. It is fitted with B ≤ 0 — the lag branch is excluded
because 0 of 16 of these curves pass the acceleration test — and τ is profiled rather than
optimised. The shaded fan shows every initial slope allowed by a τ within the 95% cost
band. Where that fan is wide the curve does not determine v<sub>0</sub>, however well the
form fits: <b>{bounded} of 27</b> curves are bounded by the {BURST_V0_HALFWIDTH:.0%}
half-width criterion, and the rest are marked UNBOUNDED.</p>
<p class='warn'>Fitting well is not the same as being identified. Exp 3 sample 3 is fitted
by the unconstrained burst form to better than its own noise, yet across statistically
indistinguishable values of τ its v<sub>0</sub> ranges from −1.6e-05 to +5.8e-07 — a factor
of 28 and a change of sign. That is why v<sub>0</sub> is profiled here and why the
quadratic, not the burst form, carries the headline numbers.</p>

<h2>The curves</h2>
<div class='grid three'>{''.join(panels)}</div>
"""
    return page("Progress curves and their fits", body)


def main():
    pages = {"index.html": build_index(),
             "progress_curves.html": build_curves_page()}
    for name, content in pages.items():
        path = os.path.join(HERE, name)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(content)
        print(f"wrote {path}  ({len(content) / 1024:.0f} kB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
