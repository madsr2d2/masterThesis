"""
Builds the figures for background_reaction/.

Draws only. Every measurement comes from `scope`, `curve_metrics` or
`summary_kinetics` -- see CLAUDE.md, "Do not re-derive measurements". The only
transforms applied here are for display: centring a run's points on its own
geometric mean so that a within-experiment slope can be seen on one axis.

    python background_reaction/build_figures.py

NOTHING NUMERIC IS TYPED INTO THE PROSE OR THE CAPTIONS. The figures always
derived their numbers; the running text did not, and an audit on 2026-09-02
found four places where it had drifted -- a buffer range that still included
exp 65 after its rates were withdrawn, an experiment list in a legend that did
the same, a pooled R2 comparison whose argument had reversed sign, and a
worked identifiability example whose numbers changed when exps 3 and 6 were
recovered from the .rre. Every one of them is now an f-string over the same
call `check_numbers.py` uses.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
# svgplot lives at the repository root, not beside this file: it is shared with
# temperature_series/build_figures.py, and a second copy is exactly the
# duplication data/test_curve_metrics.py exists to forbid.
sys.path.insert(0, os.path.dirname(HERE))

import scope
from curve_metrics import (ACCELERATION_SIGMA, INITIAL_WINDOW, OUTLIER_SIGMA,
                           acceleration, isolated_outliers, local_outlier_z,
                           model_residual)
from fit_dataset import source_floor
from summary_kinetics import BURST_V0_HALFWIDTH, fit_burst_bounded
from svgplot import ACCENT, GRID, INK, MUTED, PALETTE, Axes, esc, page
from figure_kit import (decimated, fig, panel, styled,
                        write_pages)

# The estimator the headline numbers are quoted on, and the two reported
# beside it. v0_quad answers the "the 20% window is arbitrary" objection: it is
# the slope at t = 0 of a quadratic through EVERY point, so no window is
# chosen anywhere, and being linear in its parameters it is always identified.
HEADLINE = "v0_quad"
# Five now. `v0_burst` was added on 2026-09-01 when it was proposed as the
# headline: the point of this tuple is that "does the conclusion depend on
# how the rate was measured" stays a groupby rather than an argument, and
# a candidate headline has to be in it before that question can be asked.
ESTIMATORS = ("v0_quad", "v0_burst", "vmax", "v0", "v0_whole")
ESTIMATOR_LABEL = {
    "v0_quad": "v0 quadratic (whole curve, no window)",
    "v0_burst": "v0 from the burst/lag form (whole curve, one exponential)",
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

    # Limits from the data, not hardcoded: they were pinned to (-6, 3) for a
    # v0_whole value of -3.90 that turned out to be a bug in whole_slope.
    low = min(v - e for _, v, e, _, _ in rows) - 0.3
    high = max(max(v + e, o + oe) for _, v, e, o, oe in rows) + 0.3
    axes = Axes(560, 380, (-0.6, len(rows) - 0.4),
                (min(-0.4, low), max(2.2, high)), pad=(62, 16, 74, 26))
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
# Four forms are drawn on every panel, so they need four clearly separated
# hues. The first attempt used #3f8a5a for the quadratic and #12856a for the
# burst -- both green-teal, and indistinguishable on the page.
WINDOW_COLOUR = "#2f6fb0"      # the 20% window line          -- blue
QUAD_COLOUR = "#c25e00"        # the whole-curve quadratic    -- amber (headline)
BURST_COLOUR = "#7a4bb8"       # the burst/lag form           -- purple
WHOLE_COLOUR = "#3f8a5a"       # straight line, whole curve   -- green
OUTLIER_COLOUR = "#c0392b"     # ring round a suspect reading -- red


def curve_panel(curve, width=330, height=210):
    """
    One cuvette as an HTML block: the plot, then its numbers in a real table.

    The numbers were once floating text inside the SVG. They sat on top of the
    curves, were clipped by the frame, and could not be selected or searched.
    Anything textual belongs in HTML; the SVG draws only data and fits.
    """
    from curve_metrics import initial_rate, quadratic_rate, whole_slope

    times = np.asarray(curve.times, dtype=float)
    times = times - times[0]
    values = np.asarray(curve.absorbance, dtype=float)
    floor = source_floor(curve.source)

    v0, v0_se, _ = initial_rate(curve.times, curve.absorbance, floor=floor)
    quad, quad_se, curvature = quadratic_rate(curve.times, curve.absorbance,
                                              floor=floor)
    whole, whole_se = whole_slope(curve.times, curve.absorbance, floor=floor)
    # noise_floor is the SOURCE's: it reaches the acceleration z that decides
    # whether this curve's lag branch stays open. See summary_kinetics.
    burst = fit_burst_bounded(curve.times, curve.absorbance,
                              noise_floor=floor)

    lo, hi = float(values.min()), float(values.max())
    margin = max((hi - lo) * 0.16, 4 * curve.noise)
    # A little room to the left of t = 0. The first reading sits exactly there,
    # it is the one most often ringed, and t = 0 is where all four fits are
    # contested -- drawn hard against the axis its marker is half clipped.
    span = float(times[-1])
    axes = Axes(width, height, (-0.035 * span, span),
                (lo - margin, hi + margin), pad=(56, 12, 32, 12))

    cut = float(times[-1]) * INITIAL_WINDOW
    axes.parts.append(
        f"<rect x='{axes._fx(0):.1f}' y='{axes.top}' "
        f"width='{axes._fx(cut) - axes._fx(0):.1f}' "
        f"height='{height - axes.top - axes.bottom:.1f}' fill='{WINDOW_COLOUR}' "
        f"fill-opacity='0.07'/>")

    grid = np.linspace(0, float(times[-1]), 140)
    # Every line is drawn from ITS OWN fitted intercept. Displaying a fit
    # against an intercept it did not choose puts it visibly off its own data
    # and misrepresents how well it fits.
    quad_beta = _quadratic_beta(times, values)
    window = times <= cut
    window_beta = np.polyfit(times[window], values[window], 1)
    whole_beta = np.polyfit(times, values, 1)

    if np.isfinite(burst.tau):
        # The v0 profile interval, as the fan of initial slopes the 95% cost
        # band still allows. A wide fan means the curve does not determine v0,
        # however well the form fits.
        reach = float(times[-1]) * 0.25
        axes.band(np.array([0.0, reach]),
                  burst.c + burst.v0_low * np.array([0.0, reach]),
                  burst.c + burst.v0_high * np.array([0.0, reach]),
                  BURST_COLOUR, opacity=0.20)
        axes.line(grid, burst.predict(grid), BURST_COLOUR, width=1.7, dash="6 3")
    axes.line(grid, whole_beta[1] + whole_beta[0] * grid, WHOLE_COLOUR,
              width=1.5, dash="2 3")
    axes.line(np.array([0.0, cut]), window_beta[1] + v0 * np.array([0.0, cut]),
              WINDOW_COLOUR, width=2.2)
    axes.line(grid, quad_beta[0] + quad_beta[1] * grid + quad_beta[2] * grid ** 2,
              QUAD_COLOUR, width=2.2)

    # Suspect readings, ringed rather than removed. Nothing is excluded: the
    # fits above are computed on every point, and these rings say which ones a
    # reader should discount by eye. Only ISOLATED flags are ringed -- a run of
    # two or more may be real structure, and this dataset's shapes are live
    # hypotheses (curve_screen.py).
    #
    # Worked out BEFORE the readings are drawn, because the drawing is
    # decimated and the rings are not: see `shown` below.
    # Point 0 is ringed on its own z, not on isolation: a bad leading reading
    # often drags its neighbour into a run, and the pair then hides from
    # `isolated`. The case is now the leverage and the masking, NOT a raised
    # flag rate -- since the instrument's first reading is discarded, leading
    # readings are flagged on 14.7% of curves against 16.2% for last ones.
    # What the surviving rings buy is worth having: they no longer mark the
    # generic settling artefact (that is gone with the dropped reading) but
    # runs whose settling lasted LONGER than one reading, which is a small
    # nameable set -- all four of exp 65, exp 70 sample 3, exp 3 sample 6, exp
    # 6 sample 4 -- rather than a property of the archive.
    isolated, in_runs = isolated_outliers(times, values, curve.noise)
    outlier_z = local_outlier_z(times, values, curve.noise)
    ringed = sorted(set(int(i) for i in isolated) |
                    ({0} if len(outlier_z) and np.isfinite(outlier_z[0])
                     and abs(outlier_z[0]) > OUTLIER_SIGMA else set()))

    # The DATA goes on last. Drawn first it disappeared under four fit lines,
    # which inverts the point of the panel: the fits are the claim, the
    # readings are the evidence, and the evidence has to stay visible.
    # `decimated` carries the stride and the ring union.
    shown = decimated(len(times), 110, keep=ringed)
    axes.points(times[shown], values[shown], INK, radius=1.7, opacity=0.85)

    for index in ringed:
        axes.ring(times[index], values[index], OUTLIER_COLOUR,
                  title=f"suspect reading: point {index} at t={times[index]:.0f} s")

    svg = axes.render("time / s", "\u0394A")

    # The burst row names BOTH endpoints and the shape, because v0 alone is
    # ambiguous: on a burst it is the maximum rate, on a lag it is the
    # INDUCTION rate -- the reaction before it gets going. Printing one number
    # called "v0" for both puts two different quantities in one column.
    shape = {"burst": "burst, rate falls", "lag": "LAG, rate rises",
             "clamped": "lag branch shut, B→0",
             "unresolved": "unresolved"}[burst.kind]
    if not burst.shape_is_meaningful:
        # tau ran to an end of its grid, so the model has degenerated and the
        # SHAPE means nothing -- even where v0 is well determined, which is
        # most of them. Said first because it qualifies everything after it.
        shape = "τ unresolved, shape not determined"
    verdict = "bounded" if burst.bounded else f"v0 UNBOUNDED ±{burst.half_width:.0%}"
    if burst.settles_backwards:
        verdict += " · v_ss < 0"
    # HOW WELL EACH FORM FITS, in units of the curve's own noise, so the
    # reader is not left to judge it by eye. This is a different question from
    # the interval beside v0 -- that asks whether the DATA pin the parameter,
    # this asks whether the FORM describes the curve -- and the two come apart:
    # exp 65's four cuvettes report a bounded v0 on fits sitting 7-8x above
    # noise. `model_residual` defines it once for both forms so the comparison
    # is like-for-like across their different parameter counts.
    quad_fit = np.polyval(np.polyfit(times, values, 2), times)
    quad_resid = model_residual(values, quad_fit, 3, curve.noise)
    burst_resid = model_residual(
        values,
        burst.c + burst.v_ss * times - burst.B * (1 - np.exp(-times / burst.tau)),
        4, curve.noise)

    def _fit_note(residual, text):
        if not np.isfinite(residual):
            return text
        flag = " — DOES NOT FIT" if residual > 3 else ""
        return f"{text} · {residual:.2f}× noise{flag}"

    rows = [
        (QUAD_COLOUR, "v0 quadratic", f"{quad:.3e}", f"± {quad_se:.1e}",
         _fit_note(quad_resid, "whole curve, no window")),
        (WINDOW_COLOUR, "v0 window", f"{v0:.3e}", f"± {v0_se:.1e}",
         f"first {INITIAL_WINDOW:.0%}, shaded"),
        (WHOLE_COLOUR, "slope, whole", f"{whole:.3e}", f"± {whole_se:.1e}",
         "straight line, every point"),
        (BURST_COLOUR, f"burst {'v0→v_ss'}",
         f"{burst.v0:.2e} → {burst.v_ss:.2e}",
         f"[{burst.v0_low:.1e}, {burst.v0_high:.1e}]",
         _fit_note(burst_resid, f"{shape} · {verdict}")),
    ]
    body = "".join(
        f"<tr><td><i class='sw' style='background:{colour}'></i>{esc(name)}</td>"
        f"<td class='num'>{esc(value)}</td><td class='num dim'>{esc(spread)}</td>"
        f"<td class='dim'>{esc(note)}</td></tr>"
        for colour, name, value, spread, note in rows)

    # The BUFFER SALT, not just its concentration: exp 65 is boric and the rest
    # phosphate, and MECHANISM.md argues at length that the buffers are
    # chemically different reagents rather than four ways of setting pH.
    sub = (f"{curve.substrate} · pH {curve.pH:.2f} · "
           f"{curve.buffer.lower()} {curve.buf:.4g} mM · "
           f"[sub] {curve.conditions.s0:.3g} mM · "
           f"[H2O2] {curve.conditions.h2o2:.4g} mM · "
           f"{curve.temperature:.0f} °C · "
           f"{curve.source} · noise {curve.noise:.1e} AU")
    marks = []
    if ringed:
        marks.append(f"{len(ringed)} suspect reading"
                     f"{'s' if len(ringed) > 1 else ''} ringed"
                     + (" (incl. the first plotted)" if 0 in ringed else ""))
    unringed = [i for i in in_runs if i not in ringed]
    if unringed:
        marks.append(f"{len(unringed)} more in runs, not ringed")
    foot = (f"curvature t {curvature:+.1f}"
            + (f" · burst τ {burst.tau:.3g} s" if np.isfinite(burst.tau) else "")
            + (" · " + " · ".join(marks) if marks else ""))
    return panel(f"exp {curve.experiment} · sample {curve.sample}",
                 esc(sub), svg, esc(foot), body)


def _quadratic_beta(times, values):
    """
    (c, v0, a) for DRAWING the quadratic only.

    The quoted v0 and its standard error come from
    `curve_metrics.quadratic_rate`, which floors its variance by the curve's
    source. This positions the line on the page and nothing else.
    """
    design = np.column_stack([np.ones(len(times)), times, times ** 2])
    beta, *_ = np.linalg.lstsq(design, values, rcond=None)
    return [float(v) for v in beta]


def figure_curvature():
    """Deceleration against conversion: it is not substrate depletion."""
    frame = scope.frame(scope.FREE_BNOH_ALL)
    live = frame[frame.live]
    # The floor was -90 until 2026-09-02 and exp 6 sample 4 sits at -91.66, so
    # one point was drawn outside the frame and silently clipped away. Found by
    # `svgplot.clipped_marks`, which now runs in both check_numbers.
    axes = Axes(560, 380, (0.03, 12.0), (-100.0, 15.0), xlog=True,
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
    axes.note(90, 44, f"every curve is under "
                      f"{np.ceil(live.conversion.max() * 100):.0f}% converted,",
              INK, 11.5, weight="600")
    axes.note(90, 61, "so substrate depletion cannot bend them this far",
              INK, 11.5, weight="600")
    return axes.render("conversion at end of run / %",
                       "curvature t-statistic (negative = decelerating)",
                       "the background decelerates, and depletion does not explain it")


def figure_acceleration():
    """The in-scope curves accelerate; the background does not."""
    # The experiment lists come from the constants, not from the label. This
    # one read "(65,67,69,70)" for a day after exp 65's rates were withdrawn
    # from BUFFER_FIXED on 2026-09-01, which is a legend naming a curve that
    # is not on the chart.
    def _named(experiments):
        listed = sorted(int(e) for e in experiments)
        if listed == list(range(listed[0], listed[-1] + 1)) and len(listed) > 3:
            return f"({listed[0]}-{listed[-1]})"
        return "(" + ",".join(str(e) for e in listed) + ")"

    sets = [(f"enzyme-free, [buf] fixed\n{_named(scope.BUFFER_FIXED)}",
             scope.BUFFER_FIXED, PALETTE[0]),
            (f"enzyme-free titrations\n{_named(scope.BUFFER_CONFOUNDED)}",
             scope.BUFFER_CONFOUNDED, ACCENT),
            (f"in-scope catalysed\n{_named(scope.PRIMARY_SCOPE)}",
             scope.PRIMARY_SCOPE, PALETTE[2])]
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
def _html_table(headers, rows, highlight=()):
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for index, row in enumerate(rows):
        cls = " class='hl'" if index in highlight else ""
        body.append(f"<tr{cls}>"
                    + "".join(f"<td>{c}</td>" for c in row) + "</tr>")
    return (f"<div class='scroll'><table><thead><tr>{head}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>")




def build_index():
    buffer = scope.buffer_dependence(parameter=HEADLINE)
    fixed_buf = sorted(scope.frame(scope.BUFFER_FIXED).buf.unique())
    held = (f"{fixed_buf[0]:g} mM" if len(fixed_buf) == 1
            else f"{fixed_buf[0]:g}\u2013{fixed_buf[-1]:g} mM")
    pooled = {name: scope.background_orders(terms=("s0", "h2o2", "hoo"),
                                            parameter=name)
              for name in ("v0_quad", "vmax")}
    free = scope.frame(scope.FREE_BNOH_ALL)
    whole = scope.frame(tuple(range(1, 152)))
    cell = whole[(whole.buffer == "Pyrophosphate") & (whole.substrate == "BnOH")]
    in_scope = scope.frame(scope.PRIMARY_SCOPE)
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
{fig(figure_separation(), "<b>A.</b> Every cuvette, divided by its own run's geometric "
      "mean in both axes — which is exactly what a per-experiment offset does on a log "
      "scale, so the common within-run slope is what you see. Blue: the runs that hold "
      f"[buf] at {held}. Orange: the runs where [buf] falls as [sub] rises. The same "
      "reaction reads a positive substrate order in one design and a negative one in the "
      "other. Hover a point for its buffer concentration.")}
{fig(figure_coupling(), "<b>B.</b> Why. Substrate was pipetted in and buffer volume shrank "
      "to compensate, so the two move together with g = dlog[buf]/dlog[sub] fitted within "
      "runs. Faded points are cuvettes with no live signal, which the fits exclude.")}
<p>If the rate goes as [S]<sup>a</sup>[buf]<sup>d</sup> and within the titrations
log[buf] = g·log[sub] + constant, then fitting the titrations without a buffer term
returns a′ = a + d·g, so <b>d = (a′ − a) / g</b>. Both a and a′ are within-run slopes,
so pH, [H<sub>2</sub>O<sub>2</sub>], cell and day are absorbed on both sides and nothing
asks one regression to separate [buf] from [sub].</p>

<h2>The answer: first order in buffer</h2>
{fig(figure_buffer_order(), "<b>C.</b> The recovered buffer order under four different "
      "rate estimators, against an independent block. Green dashed line: first order.")}
{_html_table(["rate estimator", "a · [buf] fixed", "a′ · [buf] falling",
         "order in [buf] · BnOH", "cross-check · 4OMe-BnOH 40 °C"], rows,
        highlight=(0,))}
<p class='warn'>All four estimators agree, including the two that choose no window
anywhere. The spread across them is smaller than the distance of any one of them from
zero. An earlier version of this page reported <code>v0_whole</code> at −3.90 ± 2.25 and
argued from it that a whole-curve straight line is biased by the curves' deceleration;
that was a bug in <code>curve_metrics.whole_slope</code>, which returned the fitted
<em>intercept</em> in place of the slope. Corrected 2026-09-01.</p>

<h2>The rate law</h2>
<p>Each order is taken from the design that can measure it: the substrate order from the
constant-buffer runs only, the peroxide and [HOO⁻] orders pooled across the runs with a
usable rate, the buffer order from the two designs' disagreement. Shown under both
estimators, because neither dominates: <code>v0_quad</code> is the tighter on substrate
and buffer, and the pooled fit is a dead heat — R²
{pooled['v0_quad']['r2']:.3f} on <code>v0_quad</code> against
{pooled['vmax']['r2']:.3f} on <code>vmax</code>, on the same
{pooled['v0_quad']['n']} curves. An earlier version of this page put those at 0.961 and
0.909 and argued from the gap that an extrapolated initial rate adds variance where a
block slope does not. The gap is {abs(pooled['v0_quad']['r2'] - pooled['vmax']['r2']):.3f}
and points the other way; the argument is withdrawn.</p>
<div class='grid two'>
{fig(figure_rate_law('v0_quad'), "Headline estimator — no window anywhere.")}
{fig(figure_rate_law('vmax'), "The estimator the in-scope block was measured with.")}
</div>

<h2>What the background is not</h2>
{fig(figure_curvature(), "The curves bend, and it is not substrate being consumed: "
      f"conversion runs {free.conversion.min() * 100:.2f}\u2013"
      f"{free.conversion.max() * 100:.2f}%, yet "
      f"{int((free.curvature_t.abs() > 3).sum())} of {len(free)} show curvature beyond "
      "|t| = 3. Something else decays during these runs — peroxide, or the cell.")}
{fig(figure_acceleration(), "Both sets compared here are entirely .rre, so this is not "
      "the variance-floor artefact of DATA_VERIFICATION.md 2026-09-01. The enzyme-free "
      "curves do not accelerate; several actively decelerate.")}

<h2>What this means for exps 135–151</h2>
<p><b>The confound does not reach them.</b> <code>[buf]</code> is
{in_scope.buf.iloc[0]:.3f} mM in all {len(in_scope)} in-scope curves across all
{in_scope.experiment.nunique()} runs — zero variation — so no buffer effect can enter
their substrate order.</p>
<p><b>The background is already subtracted.</b> Every run is double-beam and the catalysed
sheets' reference omits only the enzyme, so an in-scope curve is a catalytic increment
taken against a background at identical buffer, substrate, peroxide and pH.</p>
<p><b>What is still missing is the amplitude.</b> There are
{int((~cell.differential).sum())} enzyme-free curves in the {len(cell)}-curve BnOH
pyrophosphate cell, so the absolute size of what was subtracted cannot be
recovered without importing it from another buffer — which MECHANISM.md's buffer section
forbids. That remains FITTING.md's F7.</p>

<p style='margin-top:34px'><a href='progress_curves.html'>→ every progress curve, with all
three fits drawn</a></p>
"""
    return styled("The background reaction: BnOH + H₂O₂ without catalyst",
                  body)


def build_curves_page():
    panels = [curve_panel(curve) for curve in scope.curves(scope.FREE_BNOH_ALL)]
    bounded = sum(fit_burst_bounded(c.times, c.absorbance,
                                    noise_floor=source_floor(c.source)).bounded
                  for c in scope.curves(scope.FREE_BNOH_ALL))
    # Counted here rather than written into the prose, because it is exactly
    # the kind of number that goes stale: the legend claimed "0 of 16" for
    # three commits after the gate became per-curve and after the first-reading
    # drop pushed exp 67 sample 3 across it.
    accelerating = sum(
        bool(np.isfinite(z) and z > ACCELERATION_SIGMA)
        for z, _ in (acceleration(c.times, c.absorbance,
                                  floor=source_floor(c.source))
                     for c in scope.curves(scope.FREE_BNOH_ALL)))
    curves = list(scope.curves(scope.FREE_BNOH_ALL))
    # The identifiability example, chosen by the data rather than named in the
    # prose. It used to name exp 3 sample 3 with numbers that stopped being
    # true when exps 3 and 6 were recovered from the .rre on 2026-09-01: its
    # unconstrained v0 now runs +3.2e-07 to +6.8e-07, no sign change at all.
    widest, worst = None, 0.0
    unstable = 0
    for curve in curves:
        fitted = fit_burst_bounded(curve.times, curve.absorbance,
                                   constrain=False,
                                   noise_floor=source_floor(curve.source))
        if not (np.isfinite(fitted.v0_low) and np.isfinite(fitted.v0_high)
                and fitted.v0):
            continue
        if fitted.v0_low < 0 < fitted.v0_high:
            unstable += 1
        reach = abs((fitted.v0_high - fitted.v0_low) / fitted.v0)
        if reach > worst:
            widest, worst = (curve, fitted), reach
    example, span = widest
    body = f"""
<p class='lede'>All {len(curves)} enzyme-free BnOH cuvettes, each with every rate
estimator drawn on it. Nothing here is a kinetic-model fit: the mechanism has never been fitted to these
curves. These are the empirical rate measurements the analysis rests on.</p>

<h2>What is drawn</h2>
<div class='key'>
  <span><i class='sw' style='background:{QUAD_COLOUR}'></i>quadratic through every point — v0 is its slope at t = 0 <b>(headline)</b></span>
  <span><i class='sw' style='background:{WINDOW_COLOUR}'></i>least-squares line over the first {INITIAL_WINDOW:.0%} (shaded)</span>
  <span><i class='sw' style='background:{BURST_COLOUR}'></i>burst/lag form, dashed; shaded fan = v0 profile interval</span>
  <span><i class='sw' style='background:{WHOLE_COLOUR}'></i>straight line through the whole curve</span>
  <span><i class='sw' style='background:{OUTLIER_COLOUR};height:11px;width:11px;border-radius:50%;background:none;border:2px solid {OUTLIER_COLOUR}'></i>ringed = suspect reading, still included in every fit</span>
</div>
<p><b>The instrument's first reading of every run is not plotted, fitted or scored.</b> It
is discarded in <code>fit_dataset.build_curves</code> before anything here sees a curve,
unconditionally and for every run in the archive, because it is a settling artefact rather
than a measurement: scored against the following 8 readings it carried median |z| 2.06 and
exceeded 5σ on 15.9% of curves, where the LAST reading scored the identical way carried
1.11 and 7.5%. So the leftmost point on every panel below is the instrument's <i>second</i>
reading, and where a footer says a ringed point is the first, it means the first
plotted.</p>
<p>The <b>quadratic</b> is A = c + v<sub>0</sub>t + at². It chooses no window, uses every
point, and is linear in its parameters, so v<sub>0</sub> is always identified. Its
<code>curvature t</code> is a/se(a): |t| &gt; 3 means the curve measurably bends.</p>
<p>The <b>burst/lag</b> form is A = c + v<sub>ss</sub>t − B(1 − e<sup>−t/τ</sup>), with
v<sub>0</sub> = v<sub>ss</sub> − B/τ. The lag branch (B &gt; 0) is opened or shut
<b>per curve</b>, by that curve's own <code>acceleration</code> z-score against a 3σ gate:
open on the {accelerating} of {len(curves)} that clear it, shut elsewhere. τ is profiled rather than
optimised. The shaded fan shows every initial slope allowed by a τ within the 95% cost
band. Where that fan is wide the curve does not determine v<sub>0</sub>, however well the
form fits: <b>{bounded} of {len(curves)}</b> curves are bounded by the {BURST_V0_HALFWIDTH:.0%}
half-width criterion, and the rest are marked UNBOUNDED.</p>
<p class='warn'>Fitting well is not the same as being identified. The worst case in this
set is <b>exp {example.experiment} sample {example.sample}</b>: across statistically
indistinguishable values of τ its unconstrained v<sub>0</sub> ranges from
{span.v0_low:.2e} to {span.v0_high:.2e} — a span of {worst:.0f}× its own fitted value
{'and a change of sign' if span.v0_low < 0 < span.v0_high else ''}. It is not alone:
<b>{unstable} of {len(curves)}</b> curves have an unconstrained v<sub>0</sub> interval
that straddles zero. That is why v<sub>0</sub> is profiled here and why the quadratic,
not the burst form, carries the headline numbers.</p>

<h2>The curves</h2>
<div class='grid three'>{''.join(panels)}</div>
"""
    return styled("Progress curves and their fits", body)


def main():
    return write_pages(HERE, {"index.html": build_index(),
                              "progress_curves.html": build_curves_page()})


if __name__ == "__main__":
    raise SystemExit(main())
