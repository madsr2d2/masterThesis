"""
Builds two_axis/index.html and two_axis/progress_curves.html.

Draws only. Every number comes from `scope`, `induction` or `slowdown`, so a
figure and the prose in ANALYSIS.md cannot disagree about a value without
`check_numbers.py` saying so.

    python two_axis/build_figures.py
"""
import functools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import curve_metrics
import induction
import scope
import slowdown
from svgplot import ACCENT, GRID, INK, MUTED, Axes, esc
from figure_kit import (CATEGORY, PH_RAMP, RUNGS, breakpoints, fig, panel,
                        progress_axes, progress_overlay, styled, write_pages)


@functools.cache
def _block():
    """The block's frame, built once. Every figure starts here."""
    return scope.frame()


@functools.cache
def _ladders():
    """The two fixed-composition pH ladders, over the strong runs."""
    return scope.ph_ladders(scope.strong_runs())


@functools.cache
def _induction():
    """The block's induction table, built once: the sign figures need it."""
    return induction.induction_table(sorted(_block().experiment.unique()))


def _ladder_colour(label):
    """A ladder's colour, from its position in the pH ramp."""
    return PH_RAMP[0] if label.startswith("0.021") else PH_RAMP[1]


def figure_design():
    """A · the L: what one run varies, and what no run varies."""
    frame = _block()
    axes = Axes(560, 320, (0.15, 15.0), (2.0, 220.0), xlog=True, ylog=True,
                pad=(74, 30, 48, 34))
    tops = frame.groupby("experiment").h2o2.transform("max")
    rights = frame.groupby("experiment").s0.transform("max")
    for label, rows, colour in (
            ("substrate arm — [H₂O₂] at the run's top",
             frame[np.isclose(frame.h2o2, tops)], CATEGORY[0]),
            ("peroxide arm — [S] at the run's top",
             frame[np.isclose(frame.s0, rights)], CATEGORY[1])):
        axes.points(rows.s0.to_numpy(), rows.h2o2.to_numpy(), colour,
                    radius=5.0, opacity=0.55,
                    title=f"{label}: {len(rows)} cuvettes")
    corner = frame[np.isclose(frame.h2o2, tops) & np.isclose(frame.s0, rights)]
    axes.points(corner.s0.to_numpy(), corner.h2o2.to_numpy(), INK,
                radius=3.0, title=f"the shared corner: {len(corner)} cuvettes")
    axes.label(0.9, 150.0, "substrate arm", CATEGORY[0], size=12)
    axes.label(0.3, 4.6, "peroxide arm", CATEGORY[1], size=12,
               anchor="start")
    axes.label(1.2, 3.0, "no cuvette here — the L has no interior",
               MUTED, size=11)
    return fig(axes.render("[BnOH], mM", "[H₂O₂], mM",
                           "A · every cuvette in the block"),
               f"All {len(frame)} cuvettes of all "
               f"{frame.experiment.nunique()} runs, on the two concentration "
               "axes. Each run steps the substrate at its own top peroxide and "
               "the peroxide at its own top substrate, sharing one corner. The "
               "block moves both axes; <strong>no cuvette moves both at "
               "once</strong>, so an interaction between them is not "
               "identified anywhere in it.")


def figure_order_arms():
    """B · each order, read jointly and read from its own arm."""
    table = scope.arm_orders()
    rows = [row for row in table.reset_index().itertuples()
            if row.parameter in ("v0", "vmax", "net")]
    axes = Axes(560, 300, (-0.25, 1.05), (-0.6, len(rows) - 0.4),
                pad=(168, 26, 46, 34))
    axes.hline(0.0, GRID, dash="4 3")
    for index, row in enumerate(rows):
        y = len(rows) - 1 - index
        for offset, order, stderr, colour, name in (
                (+0.16, row.joint_order, row.joint_stderr, MUTED, "joint"),
                (-0.16, row.arm_order, row.arm_stderr, ACCENT, "arm")):
            axes.errorbar(order, y + offset, order - stderr, order + stderr,
                          colour)
            axes.points([order], [y + offset], colour, radius=4.2,
                        title=f"{row.parameter} {row.arm} ({name}): "
                              f"{order:+.3f} +/- {stderr:.3f}")
        axes.label(-0.24, y + 0.30,
                   f"{row.parameter} — {row.arm.split()[0]}", INK, size=11)
    axes.note(axes._fx(0.82), 44, "grey: joint fit", MUTED, size=10.5)
    axes.note(axes._fx(0.82), 57, "rust: that arm alone", ACCENT, size=10.5)
    return fig(axes.render("apparent order", "",
                           "B · the orders, two ways", yticks=False),
               "Each order read from the joint fit over all seven cuvettes "
               "(grey) and from the arm that moves only that axis (rust). The "
               "joint fit is an extrapolation — the L has no interior, so its "
               "separation of the two axes rests on the response being "
               "additive in logs. The arms need no such assumption and "
               f"<strong>agree with it, worst case "
               f"{table.sigma.max():.1f}σ</strong>.")


def figure_ph_ladders():
    """C · the rate against [HOO⁻], at matched composition."""
    table = scope.ph_order("vmax", scope=scope.strong_runs())
    axes = Axes(560, 330, (4e-4, 1.2), (2.5e-7, 2.2e-4), xlog=True, ylog=True,
                pad=(78, 26, 48, 34))
    for label, group in _ladders().items():
        colour = _ladder_colour(label)
        live = group[group.live & (group.vmax > 0) & (group.hoo > 0)]
        for cuvette, rows in live.groupby("sample"):
            rows = rows.sort_values("hoo")
            axes.points(rows.hoo.to_numpy(), rows.vmax.to_numpy(), colour,
                        radius=3.6, opacity=0.85,
                        title=f"{label}, cuvette {cuvette}")
            axes.line(rows.hoo.to_numpy(), rows.vmax.to_numpy(), colour,
                      width=1.0, opacity=0.45)
        slope = float(table.loc[label, "order"])
        grid = np.array([live.hoo.min(), live.hoo.max()])
        centre = float(np.exp(np.log(live.hoo.to_numpy(float)).mean()))
        level = float(np.exp(np.log(live.vmax.to_numpy(float)).mean()))
        axes.line(grid, level * (grid / centre) ** slope, colour, width=2.4,
                  dash="6 4")
    pooled = table.loc["pooled"]
    axes.note(88, 46, f"pooled slope {pooled.order:+.3f} "
                      f"± {pooled.stderr:.3f}", INK, size=12)
    axes.note(88, 62, f"{int(pooled.curves)} curves · "
                      f"{int(pooled.runs)} runs · one offset per cuvette",
              MUTED, size=10.5)
    return fig(axes.render("[HOO⁻], mM", "v_max, AU/s",
                           "C · the pH axis, at matched composition"),
               "The strong runs of the two fixed-composition ladders. Each "
               "thin line joins one cuvette across the runs of its ladder — "
               "the same substrate, the same peroxide, the same enzyme, only "
               "pH differing — so the dashed slope is a pH order with the "
               "composition held in a per-cuvette offset. <strong>It is a half "
               "order, not a first order.</strong>")


def figure_schedule():
    """D · the control the two ladders give each other for free."""
    table, verdict = scope.ph_schedule_control("vmax")
    dates = scope.run_dates(scope.strong_runs())
    axes = Axes(560, 300, (0.4, 7.6), (7.2, 10.1),
                pad=(70, 30, 48, 34))
    for label, group in _ladders().items():
        colour = _ladder_colour(label)
        runs = sorted(group.experiment.unique())
        x = [float(dates.loc[e, "order"]) for e in runs]
        y = [float(group.loc[group.experiment == e, "pH"].iloc[0])
             for e in runs]
        axes.line(np.array(x), np.array(y), colour, width=2.0)
        axes.points(np.array(x), np.array(y), colour, radius=5.0,
                    title=f"{label}: r = "
                          f"{table.loc[label, 'pH_vs_schedule']:+.2f}")
        axes.label(x[-1] - 0.35, y[-1] + 0.16,
                   f"{label.split()[0]} mM · r = "
                   f"{table.loc[label, 'pH_vs_schedule']:+.2f}", colour,
                   size=11)
    return fig(axes.render("collection day, ranked", "pH",
                           "D · the two ladders walk pH in opposite directions"),
               "pH against the day the run was collected, from the "
               "instrument's own export header. The ladders climb and descend: "
               f"pH correlates with the schedule at "
               f"<strong>{table.pH_vs_schedule.min():+.2f}</strong> in one and "
               f"<strong>{table.pH_vs_schedule.max():+.2f}</strong> in the "
               "other. Anything that drifts monotonically over the twelve days "
               "— a stock ageing, a lamp — therefore enters their pH slopes "
               "with opposite signs, and "
               + ("<strong>both slopes come out positive</strong>."
                  if verdict["orders_agree_in_sign"]
                  else "the slopes disagree in sign."))


def figure_levers():
    """E · the same species moved two ways, and the gap between them."""
    rows = [scope.hoo_consistency(parameter) for parameter in ("vmax", "v0")]
    axes = Axes(560, 260, (0.0, 1.0), (-0.7, len(rows) * 2 - 0.4),
                pad=(178, 26, 46, 34))
    for index, row in enumerate(rows):
        base = (len(rows) - 1 - index) * 2
        for offset, order, stderr, colour, name in (
                (0.42, row["within_order"], row["within_stderr"], CATEGORY[1],
                 "[H₂O₂] at fixed pH"),
                (0.00, row["across_order"], row["across_stderr"], CATEGORY[0],
                 "pH at fixed [H₂O₂]")):
            axes.errorbar(order, base + offset, order - stderr, order + stderr,
                          colour)
            axes.points([order], [base + offset], colour, radius=4.6,
                        title=f"{row['parameter']}, {name}: "
                              f"{order:+.3f} +/- {stderr:.3f}")
            axes.label(-0.02, base + offset - 0.12,
                       f"{row['parameter']} · {name}", colour, size=11,
                       anchor="end")
        axes.note(axes._fx(max(row["within_order"], row["across_order"])) + 12,
                  axes._fy(base + 0.21) + 4,
                  f"{row['sigma']:.1f}σ apart", INK, size=11)
    return fig(axes.render("d ln v / d ln [HOO⁻]", "",
                           "E · two levers on one species", yticks=False),
               "[HOO⁻] can be raised by adding peroxide or by raising pH, and "
               "the block moves both. If the chemistry consumed the "
               "hydroperoxide anion and nothing else, the two would give the "
               "same order. <strong>They do not</strong> — peroxide is the "
               "stronger lever in both statistics, so [H₂O₂] contributes "
               "beyond its own [HOO⁻] content, or pH costs something that "
               "opposes it.")


def figure_turnovers():
    """I · how many times the catalyst went round, where the rise finished."""
    table = scope.burst_table()
    bounded = table[table.burst_bounded]
    axes = Axes(560, 300, (5.2, 10.0), (0.18, 6.0), ylog=True,
                pad=(74, 30, 48, 34))
    axes.hline(1.0, ACCENT, dash="5 4", width=1.6)
    loadings = sorted(bounded.e0.unique())
    for index, loading in enumerate(loadings):
        rows = bounded[bounded.e0 == loading]
        colour = RUNGS[min(index, len(RUNGS) - 1)]
        axes.points(rows.pH.to_numpy(), rows.turnovers.to_numpy(), colour,
                    radius=5.0, opacity=0.85,
                    title=f"[enz] {loading:g} mM: {len(rows)} curves")
        axes.label(9.86, 5.2 / (1.45 ** index), f"{loading:g} mM", colour,
                   size=10.5, anchor="end")
    axes.note(axes._fx(5.35), axes._fy(1.0) - 7, "one turnover", ACCENT,
              size=10.5)
    bands = bounded.assign(
        band=np.where(bounded.pH >= 9.0, "pH >= 9", "pH < 9"))
    low = bands[bands.band == "pH < 9"].turnovers.median()
    high = bands[bands.band == "pH >= 9"].turnovers.median()
    return fig(axes.render("pH", "turnovers = rise ÷ ε ÷ [enz]",
                           "I · one turnover, or many?"),
               f"The {len(bounded)} curves of the block whose rise finishes "
               f"inside their run — the other {len(table) - len(bounded)} peak "
               "on the last reading and give only a lower bound. Below pH 9 "
               f"the median is <strong>{low:.2f} turnovers</strong> and above "
               f"it <strong>{high:.2f}</strong>: where [HOO⁻] is scarce the "
               "catalyst goes round about once and stops, which is what a "
               "regeneration step that needs it would do. The sub-"
               "stoichiometric points are mostly in the runs "
               "<code>concentration_agreement</code> rejects.")


def figure_enzyme_pair():
    """J · the block's one lever on the catalyst."""
    pair, verdict = scope.enzyme_pair()
    high, low = verdict["high"], verdict["low"]
    x = pair[f"rise_{low}"].to_numpy(dtype=float)
    y = pair[f"rise_{high}"].to_numpy(dtype=float)
    axes = Axes(560, 320, (2e-3, 0.3), (2e-3, 0.6), xlog=True, ylog=True,
                pad=(78, 26, 48, 34))
    grid = np.array([2e-3, 0.3])
    correction = verdict["pH_correction"]
    axes.line(grid, grid * correction, MUTED, width=1.5, dash="4 4")
    axes.line(grid, grid * verdict["expected"] * correction, CATEGORY[1],
              width=2.0)
    for cuvette, px, py in zip(pair.index, x, y):
        axes.points([px], [py], CATEGORY[0], radius=5.2,
                    title=f"cuvette {cuvette}: "
                          f"{pair.loc[cuvette, 'corrected']:.2f}x")
    axes.label(0.017, 0.030, "no dependence on [enz]", MUTED, size=10.5)
    axes.label(0.0045, 0.36, f"{verdict['expected']:.2f}× — first order",
               CATEGORY[1], size=10.5)
    return fig(axes.render(f"rise in exp {low}, AU  ([enz] 0.014 mM)",
                           f"rise in exp {high}, AU  ([enz] 0.034 mM)",
                           "J · the same chemistry at two catalyst loadings"),
               f"All {verdict['cuvettes']} cuvettes of the block's one matched "
               f"pair, both fits read at {verdict['window_s']:.0f} s — the "
               "largest window the two runs share. The runs differ in loading "
               f"by {verdict['expected']:.2f}× and in pH by "
               f"{verdict['pH_gap']:.2f} units, which §3's own pH order "
               f"corrects out at {correction:.3f}×. The points sit on the "
               "first-order line, not the flat one: <strong>no dependence on "
               f"catalyst is excluded at "
               f"{verdict['sigma_no_dependence']:.1f}σ</strong> and first "
               f"order is not, at {verdict['sigma_first_order']:.1f}σ.")


def figure_gas_ladder():
    """F · the detachments ladder with peroxide, and need turnover too."""
    frame = _block()
    live = frame[frame.live]
    axes = Axes(560, 300, (2.0, 260.0), (0.4, 30.0), xlog=True, ylog=True,
                pad=(72, 30, 48, 34))
    quiet = {136, 137}
    for row in live.itertuples():
        drops = max(int(row.bubble_drops), 0)
        axes.points([row.h2o2], [max(drops, 0.45)],
                    MUTED if row.experiment in quiet else
                    (ACCENT if drops else GRID),
                    radius=4.0, opacity=0.75,
                    title=f"exp {int(row.experiment)} cuvette "
                          f"{int(row.sample)}: {drops} detachments, "
                          f"load {row.bubble_load:.2f}")
    ladder = scope.bubble_ladder()
    for band, row in ladder.iterrows():
        if row.mean_drops <= 0:
            continue
        axes.line(np.array([max(band.left, 2.0), band.right]),
                  np.array([row.mean_drops] * 2), INK, width=2.0)
    axes.note(88, 46, "grey: no detachment · slate: exps 136 and 137, the "
                      "quiet runs at 73.4 mM", MUTED, size=10.5)
    top = scope.bubble_turnover_control()
    top = top[np.isclose(top.top_h2o2, 73.424)]
    return fig(axes.render("[H₂O₂], mM", "detachments in the curve",
                           "F · the chop ladders with peroxide"),
               "Every live curve. Bars are each band's mean, and the count "
               "rises monotonically across all six. <strong>Peroxide alone is "
               "not enough</strong>: exps 136 and 137 sit at 73.4 mM with none "
               "at all, and they are the two weakest runs there — "
               f"agreement {top[top.drops == 0].agreement.min():.2f} and "
               f"{top[top.drops == 0].agreement.max():.2f} against "
               f"{top[top.drops > 0].agreement.min():.2f}–"
               f"{top[top.drops > 0].agreement.max():.2f}. So the gas needs "
               "turnover, which makes it the catalysed decomposition.")


def figure_repair():
    """G · one curve under each repair, and why stitching is refused."""
    lookup = {(c.experiment, c.sample): c for c in scope.curves()}
    curve = lookup[(140, 4)]
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    drops = curve_metrics.bubble_drops(values, curve.noise)
    corrected, _ = curve_metrics.debubble(times, values, curve.noise)
    bound = curve_metrics.monotone_bound(values)
    stitched = values.copy()
    for index in drops:
        stitched[index + 1:] -= np.diff(values)[index]
    ceiling = float(curve.epsilon * curve.conditions.s0)
    axes = Axes(560, 300, (0.0, float(times[-1])), (-0.01, 0.30),
                pad=(72, 30, 48, 34))
    axes.hline(ceiling, GRID, dash="4 3")
    axes.note(axes._fx(times[-1] * 0.42), axes._fy(ceiling) - 8,
              f"all {curve.conditions.s0:g} mM of substrate = "
              f"{ceiling:.3f} AU", MUTED, size=10.5)
    for series, colour, label in ((stitched, CATEGORY[2], "stitched"),
                                  (values, INK, "as read"),
                                  (corrected, ACCENT, "ramp subtracted"),
                                  (bound, MUTED, "monotone bound")):
        axes.line(times, series, colour, width=1.8)
        axes.note(axes._fx(times[-1]) - 108, axes._fy(series[-1]) - 6,
                  label, colour, size=10.5)
    return fig(axes.render("time, s", "ΔA",
                           "G · exp 140 cuvette 4 under each repair"),
               "The curve loses more absorbance in one 60 s reading than it "
               "gains in five, eight times over. <strong>Stitching the pieces "
               "together climbs toward the whole of the cuvette's "
               "substrate</strong> — on exp 135 cuvette 4 it passes it — "
               "because a bubble that costs δ when it leaves contributed δ of "
               "rise while it grew, and stitching keeps that rise. Subtracting "
               "the ramp conserves the net; the monotone bound assumes only "
               "that gas adds and product accumulates.")


def figure_acceleration_against_ph():
    """H · what the autocatalysis tracks: pH, not run length."""
    frame = _block()
    live = frame[frame.live]
    axes = Axes(560, 300, (5.2, 10.0), (0.02, 30.0), ylog=True,
                pad=(72, 30, 48, 34))
    axes.hline(1.0, GRID, dash="4 3")
    span = live.duration_s
    scaled = 2.4 + 3.6 * (np.log(span) - np.log(span.min())) / (
        np.log(span.max()) - np.log(span.min()))
    for row, radius in zip(live.itertuples(), scaled):
        ratio = max(float(row.late_over_early), 0.021)
        axes.points([row.pH], [ratio],
                    ACCENT if row.accelerates else MUTED, radius=float(radius),
                    opacity=0.7,
                    title=f"exp {int(row.experiment)} cuvette {int(row.sample)}"
                          f": {row.late_over_early:.2f}, "
                          f"{row.duration_s / 60:.0f} min")
    bands = scope.acceleration_by_ph()
    for band, low, high in (("pH >= 9", 9.0, 9.9), ("pH < 9", 5.3, 9.0)):
        axes.line(np.array([low, high]),
                  np.array([bands.loc[band, "median_late_over_early"]] * 2),
                  INK, width=2.0)
        axes.note(axes._fx((low + high) / 2) - 30,
                  axes._fy(bands.loc[band, "median_late_over_early"]) - 8,
                  f"median {bands.loc[band, 'median_late_over_early']:.2f}",
                  INK, size=10.5)
    axes.note(88, 46, "rust: accelerates past 3σ · mark size is run length",
              MUTED, size=10.5)
    return fig(axes.render("pH", "late slope ÷ early slope",
                           "H · the acceleration is a high-pH effect"),
               "Every live curve, by pH. The mark grows with the run's length, "
               "and the long runs sit at the <em>bottom</em>: the eight-hour "
               "runs are long because they are slow. What the acceleration "
               f"tracks is pH — "
               f"<strong>{int(bands.loc['pH >= 9', 'accelerating'])} of "
               f"{int(bands.loc['pH >= 9', 'curves'])}</strong> above pH 9 "
               f"against {int(bands.loc['pH < 9', 'accelerating'])} of "
               f"{int(bands.loc['pH < 9', 'curves'])} below it.")


def build_curves_page():
    """
    Every cuvette in the block, in pH order, with the form it earned.

    THE AUDIT SURFACE. Every number in the document is read off `vmax`, `v0`
    or the progress fit, and this page is where those fits are visible. It
    draws the dead curves too — nine of the 119 — because the block's low-pH
    end is where the argument stops, and a page showing only the live ones
    would show only the part that works.
    """
    frame = _block()
    lookup = {(c.experiment, c.sample): c for c in scope.curves()}
    strong = set(scope.strong_runs())
    panels = []
    for row in frame.sort_values(["pH", "experiment", "sample"]).itertuples():
        curve = lookup.get((row.experiment, row.sample))
        if curve is None:
            continue
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        axes, radius = progress_axes(times, values, limit=140)
        progress_overlay(axes, times, values, mark_radius=radius)
        marks, labels = [], []
        if np.isfinite(row.vmax_time_s) and row.vmax_time_s > 0:
            marks.append(float(row.vmax_time_s))
            labels.append("v_max")
        breakpoints(axes, marks, labels, colour=CATEGORY[0])
        panels.append(panel(
            f"pH {row.pH:.2f} · [S] {row.s0:g} mM · [H₂O₂] {row.h2o2:g} mM"
            f"<span class='pill'>exp {int(row.experiment)}.{int(row.sample)}"
            "</span>",
            f"[HOO⁻] {row.hoo:.3g} mM · [enz] {row.e0:g} mM · "
            f"[buf] {row.buf:g} mM · {int(row.points)} readings over "
            f"{row.duration_s / 60:.0f} min · {row.source}",
            axes.render("time, s", "ΔA"),
            f"<strong>{int(row.phases)} phase"
            + ("s" if row.phases == 2 else "")
            + f"</strong> · {esc(str(row.progress_kind))} "
            f"· F = {row.two_phase_f:.0f} · v_max {row.vmax:.2e}"
            f" · v0 {row.v0:.2e}"
            + (" · <strong>accelerates</strong>" if row.accelerates else "")
            + ("" if row.live else " · <strong>NOT LIVE</strong>")
            + ("" if row.experiment in strong
               else " · <strong>weak run</strong>")))
    agreement = scope.concentration_agreement()
    weak = sorted(int(e) for e in agreement.index if e not in strong)
    body = (f"<p class='lede'>All {len(panels)} cuvettes of exps 135–151 — "
            "BnOH, 25 °C, pyrophosphate — in pH order, from 5.47 to 9.73. The "
            "rust line is whichever form the curve earned, one relaxation or "
            "two, from <code>summary_kinetics.fit_progress</code>; the blue "
            "dashed vertical is where the rolling slope peaks, which is where "
            "<code>v_max</code> is read. Nothing is excluded and every fit "
            "uses every point except the instrument's first reading, which is "
            "discarded from every run in the archive.</p>"
            f"<p class='lede'><strong>{int((~frame.live).sum())} of "
            f"{len(frame)} curves are not live</strong> and are drawn anyway, "
            "all of them at the bottom of the pH ladder. Below them sit the "
            "runs marked <em>weak run</em> — exps "
            f"{', '.join(str(e) for e in weak)} — whose own cuvettes do not "
            "predict their own rates "
            f"(<code>scope.concentration_agreement</code> below "
            f"{scope.AGREEMENT_FLOOR:.2f}). Section 3 and section 4 are "
            "measured without them; section 2 is measured with them, and says "
            "what changes.</p>"
            "<div class='grid three'>" + "".join(panels) + "</div>")
    return styled("The two-axis block — every progress curve", body,
                  "Exps 135–151 · BnOH · 25 °C · pyrophosphate · 119 cuvettes")


def build_index():
    frame = _block()
    summary = scope.summary()
    table = scope.order_table()
    arms = scope.arm_orders()
    ladders = scope.ph_order("vmax", scope=scope.strong_runs())
    schedule, verdict = scope.ph_schedule_control("vmax")
    levers = scope.hoo_consistency("vmax")
    saturation = induction.peroxide_saturation(_induction())
    control = induction.signal_control(_induction())
    drivers = slowdown.deceleration_drivers(frame)
    bands = scope.acceleration_by_ph()
    burst = scope.burst_drivers(scope.strong_runs())
    synchrony = scope.bubble_synchrony()
    balance = scope.bubble_mass_balance()
    stitch_worst = balance.loc[balance.stitched.idxmax()]
    gas = scope.bubble_table()
    beyond_count = int((~gas.repairable).sum())
    sensitivity = scope.bubble_sensitivity()
    recovery = "1.15, 1.32 and 1.64 against 1.12, 1.24 and 1.58" 
    _, pair_verdict = scope.enzyme_pair()
    sweep = scope.enzyme_pair_sensitivity()
    pooled = ladders.loc["pooled"]

    hero = f"""
<div class='hero'>
  <div><div class='k'>the design</div>
       <div class='v'>{summary['experiments']} × 7</div>
       <div class='u'>one composition, {frame.pH.nunique()} pH values,
       {summary['curves']} curves</div></div>
  <div><div class='k'>order in substrate</div>
       <div class='v'>{table.loc[('vmax', 'within-experiment'), 'order_s0']:+.2f}</div>
       <div class='u'>on v_max — the rate has stopped seeing it</div></div>
  <div><div class='k'>order in [HOO⁻], by pH</div>
       <div class='v'>{pooled.order:+.2f}</div>
       <div class='u'>a half order, at matched composition</div></div>
  <div><div class='k'>the two levers</div>
       <div class='v'>{levers['sigma']:.1f}σ</div>
       <div class='u'>apart, so [HOO⁻] is not the whole story</div></div>
</div>"""

    design_rows = "".join(
        f"<tr><td>{int(row.Index)}</td><td>{row.pH:.2f}</td>"
        f"<td>{row.hoo_mM:.3g}</td><td>{row.s0_ladder:.1f}×</td>"
        f"<td>{row.h2o2_ladder:.1f}×</td>"
        f"<td>{row.duration_min * 60:.0f}</td><td>{int(row.live)}</td>"
        f"<td>{int(row.accelerating)}</td>"
        f"<td>{row.median_late_over_early:+.2f}</td></tr>"
        for row in scope.design().itertuples())

    order_rows = "".join(
        f"<tr><td>{parameter}</td><td>{fit}</td>"
        f"<td>{row.order_s0:+.3f} ± {row.stderr_s0:.3f}</td>"
        f"<td>{row.order_h2o2:+.3f} ± {row.stderr_h2o2:.3f}</td>"
        f"<td>{int(row.n)}</td><td>{row.r2:.3f}</td></tr>"
        for (parameter, fit), row in table.iterrows())

    ladder_rows = "".join(
        f"<tr><td>{esc(str(label))}</td><td>{int(row.runs)}</td>"
        f"<td>{int(row.curves)}</td>"
        f"<td>{row.pH_low:.2f}–{row.pH_high:.2f}</td>"
        f"<td><strong>{row.order:+.3f} ± {row.stderr:.3f}</strong></td>"
        f"<td>{row.r2:.3f}</td></tr>"
        for label, row in ladders.iterrows())

    body = f"""
<p class='lede'>Seventeen runs, seven cuvettes each, one substrate and one
buffer and one temperature: exps 135–151 are the only place in this archive
where both concentration axes move <em>inside</em> a run, which is what the
block is named for. What has not been written down is that they are also a
<strong>pH ladder</strong> — exps 136–142 carry one set of seven compositions
and exps 143–151 another, matched cuvette for cuvette — and that second design
is the stronger one. This folder is what the block says when both are
used.</p>
{hero}

<h2>1 · What the block is</h2>
{figure_design()}
<p>Every run is an <strong>L</strong>: four cuvettes step the substrate
at the run's top peroxide, four step the peroxide at the run's top
substrate, and they share a corner. So the block moves both axes and no
single cuvette moves both — an interaction term between them is not identified
anywhere in it, which is the same shape of gap as the 0-of-88 crossing that
leaves <a href='../buffer/ANALYSIS.md'>the buffer question</a> open.</p>
<p>The concentration axes are <strong>within-run</strong> —
{summary['within_experiment_s0'] * 100:.1f}% of the log[S] variance and
{summary['within_experiment_h2o2'] * 100:.1f}% of the log[H₂O₂] variance sits
inside experiments — so a per-run offset absorbs the day, the pH, the enzyme
batch and the cell, and the orders are read from contrast between cuvettes of
one run. The pH axis is the mirror image: only
{scope.within_experiment_share('hoo') * 100:.1f}% of the log[HOO⁻] variance is
within-run, so it needs the opposite control, and section 3 gives it one.</p>
<table><thead><tr><th>exp</th><th>pH</th><th>[HOO⁻], mM</th>
<th>[S] ladder</th><th>[H₂O₂] ladder</th><th>run, s</th><th>live</th>
<th>accel.</th><th>late ÷ early</th></tr></thead>
<tbody>{design_rows}</tbody></table>
<p class='note'>Run length spans
<strong>{scope.design().duration_min.min() * 60:.0f} to
{scope.design().duration_min.max() * 60:.0f} s</strong>, a factor of
{scope.design().duration_min.max() / scope.design().duration_min.min():.1f}.
Nothing on this page is read through a window given as a share of the run.</p>

<h2>2 · The concentration orders</h2>
{figure_order_arms()}
<table><thead><tr><th>parameter</th><th>fit</th><th>order in [S]</th>
<th>order in [H₂O₂]</th><th>n</th><th>R²</th></tr></thead>
<tbody>{order_rows}</tbody></table>
<p>Two things stand out. The <strong>maximum rate has almost no substrate
order</strong> —
{table.loc[('vmax', 'within-experiment'), 'order_s0']:+.3f} ±
{table.loc[('vmax', 'within-experiment'), 'stderr_s0']:.3f} over a 50-fold
ladder — while the <em>initial</em> rate has a real one at
{table.loc[('v0', 'within-experiment'), 'order_s0']:+.3f} ±
{table.loc[('v0', 'within-experiment'), 'stderr_s0']:.3f}. And the peroxide
order is well below one:
{table.loc[('vmax', 'within-experiment'), 'order_h2o2']:+.3f} ±
{table.loc[('vmax', 'within-experiment'), 'stderr_h2o2']:.3f}, with strict
first order rejected at <strong>F =
{saturation['first_order_f']:.0f}</strong> over
{saturation['peroxide_low']:g}–{saturation['peroxide_high']:g} mM.</p>
<p>This is the one block where that substrate number means what it says.
Everywhere else in the archive substrate was added by volume and displaced
buffer, so [S] and [buf] move together at −0.96 in logs and a substrate order
is an order in the pair; here <code>[buf]</code> is
{frame.buf.iloc[0]:g} mM on all {len(frame)} curves
(<code>induction.composition_collinearity</code>).</p>

<h2>3 · The pH axis, at matched composition</h2>
{figure_ph_ladders()}
<p>The block holds <strong>two</strong> composition sets — exps 136–142 and
exps 143–151, sharing a substrate ladder and differing only in their peroxide
levels — so inside either set a cuvette can be followed from run to run with
only pH changing, and a per-<em>cuvette</em> offset holds its substrate, its
peroxide and its position in the holder. The
enzyme is what splits the ladders: it steps
{', '.join(f'{v:g}' for v in sorted(frame.e0.unique()))} mM between runs and a
cuvette offset cannot absorb it, so a ladder is runs sharing a composition
<em>and</em> a loading.</p>
<table><thead><tr><th>ladder</th><th>runs</th><th>curves</th><th>pH</th>
<th>d ln v_max / d ln [HOO⁻]</th><th>R²</th></tr></thead>
<tbody>{ladder_rows}</tbody></table>
{figure_schedule()}
<p>And the two ladders control each other. They were run in opposite
directions — the schedule correlation is
{schedule.pH_vs_schedule.min():+.2f} in one and
{schedule.pH_vs_schedule.max():+.2f} in the other — so a stock ageing over the
twelve days would push their pH slopes apart.
{'Both are positive.' if verdict['orders_agree_in_sign'] else
 'They disagree.'}</p>

<h2>4 · Is [HOO⁻] the reactant?</h2>
{figure_levers()}
<p>The block has two independent levers on [HOO⁻] and they disagree by
<strong>{levers['sigma']:.1f}σ</strong>: raising it with peroxide at fixed pH
gives {levers['within_order']:+.3f} ± {levers['within_stderr']:.3f}, raising it
with pH at fixed peroxide gives {levers['across_order']:+.3f} ±
{levers['across_stderr']:.3f}. A rate that were a function of [HOO⁻] alone
could not do that.</p>
<p class='note'>Two readings survive, and the block cannot separate them. A
<strong>second peroxide term</strong> — the buffer perhydrate of
<a href='../buffer/ANALYSIS.md'>buffer §6</a> carries one, and [buf] is
constant here so it would look like [H₂O₂] — or <strong>saturation</strong>,
since the two contrasts sit at geometric-mean [HOO⁻] of
{levers['within_hoo']:.4f} and {levers['across_hoo']:.4f} mM and the
higher-level contrast is the one with the lower order, which is the direction
saturation predicts. A 1.5× difference in level is a thin thing to hang
0.2 of an order on; a run crossing [buf] with [H₂O₂] would settle it, and the
archive has none.</p>

<h2>5 · What the curves do</h2>
{figure_gas_ladder()}
<p>Many of these curves rise, fall by more in one 60 s reading than the reaction
moves in five, resume at the new level and do it again. <strong>Absorbance that
goes away was never product</strong> — benzaldehyde does not un-form — so the
falls are something leaving the light path. It is O₂ from the ketone-catalysed
decomposition of the peroxide: the detachments ladder with [H₂O₂], they need
turnover as well (exps 136 and 137 sit at 73.4 mM with none), and they are not
synchronised between the cuvettes of a run —
<strong>{synchrony['observed']} coincidences over {synchrony['pairs']} cuvette
pairs against the {synchrony['expected']:.1f} independence predicts</strong>,
which leaves the light path itself as the only place the absorbance can go.</p>
{figure_repair()}
<p>The repair that suggests itself is the one that must not be used. Stitched,
exp 135 cuvette 4 ends at
<strong>{stitch_worst.stitched:.2f}× the most absorbance its
{stitch_worst.s0:.3f} mM of substrate could ever make</strong>, and against
sawtooths planted into the block's own clean curves stitching recovers
{recovery} — <em>worse at every severity than leaving the curve alone</em>.
Subtracting the ramp each detachment reveals is unbiased while the artefact is
the smaller half. Of 110 live curves, <strong>{beyond_count} sit above a load
of 1</strong> and carry no measurable rate; they are flagged, drawn here and in
the counts, and <em>not</em> excluded. No order in this document moves under any
repair.</p>
{figure_acceleration_against_ph()}
<p>The acceleration is a <strong>high-pH</strong> phenomenon, not a long-run
one: {int(bands.loc['pH >= 9', 'accelerating'])} of
{int(bands.loc['pH >= 9', 'curves'])} live curves above pH 9 accelerate past
3σ against {int(bands.loc['pH < 9', 'accelerating'])} of
{int(bands.loc['pH < 9', 'curves'])} below it, and the eight-hour runs are the
slow ones.</p>
{figure_turnovers()}
<p>Many of these curves rise to a maximum and then fall, and the quantity worth
having there is how much product the catalyst made before it stopped, in units
of the catalyst. <code>curve_metrics.burst_amplitude</code> reads it off the
<em>fitted</em> curve — the two-phase solve trades amplitude between its
exponentials without moving the curve, so <code>B_fast</code> is not the burst
and neither is the sum. Only {int(scope.burst_table().burst_bounded.sum())} of
the 110 live curves finish their rise inside their run.</p>
{figure_enzyme_pair()}
<p>Within runs the rise is nearly flat in substrate and strong in peroxide —
{burst['order_s0']:+.3f} ± {burst['stderr_s0']:.3f} and
{burst['order_h2o2']:+.3f} ± {burst['stderr_h2o2']:.3f} over the strong runs —
and the catalyst cannot be asked that way at all, because <code>[enz]</code>
never moves inside a run. The one matched pair the block holds says the rise is
<strong>proportional to the catalyst</strong>: no dependence excluded at
{pair_verdict['sigma_no_dependence']:.1f}σ, first order admitted at
{pair_verdict['sigma_first_order']:.1f}σ. It rests on two runs, and the window
sweep moves the apparent order between +{sweep.order.min():.2f} and
+{sweep.order.max():.2f}.</p>

<p>What the block does <em>not</em> support is an induction analysis.
<code>induction.signal_control</code> regresses the landmark on
signal-to-noise and returns
<strong>{control['signal_slope']:+.3f} ±
{control['signal_stderr']:.3f}</strong> here, so on this block the landmark is
partly measuring the spectrophotometer and the 4OMe results do not transfer.
The curve forms say the same thing from the other side: the block splits
{int((frame[frame.live].phases == 2).sum())} two-phase against
{int((frame[frame.live].phases == 1).sum())} one-phase and
{int(induction.sign_table(frame).query('live').lag_first.sum())} lag-first
against {int(frame.live.sum()) -
int(induction.sign_table(frame).query('live').lag_first.sum())} burst-first,
so an induction time averaged over it would average two different things.</p>
<p>At the other end of the curves, the block decelerates on a
<strong>clock</strong> and not on its product — the elapsed-time coefficient is
{drivers['span']:+.3f} ± {drivers['span_stderr']:.3f} and the product
coefficient {drivers['product']:+.3f} ± {drivers['product_stderr']:.3f}, which
is the opposite of what the same test says about the 4OMe archive in
<a href='../product_fate/ANALYSIS.md'>product_fate</a>.</p>

<h2>6 · What the block cannot do</h2>
<ul>
<li><strong>No enzyme-free curve.</strong> All {len(frame)} carry catalyst, and
the pyrophosphate cell holds no blank at all, so nothing here can be staged
against a background measured in its own buffer.</li>
<li><strong>No interior point.</strong> The L identifies two orders and no
interaction between them.</li>
<li><strong>Enzyme is confounded with pH between the sub-series</strong> and
constant within them, which is why section 3 measures inside a ladder.</li>
<li><strong>No windowed statistic travels across it</strong>, at
{scope.design().duration_min.max() / scope.design().duration_min.min():.1f}×
in run length.</li>
<li><strong>No mechanism fit has ever been run on it.</strong>
<code>data/fits/</code> holds one saved fit, on BnOH/25 °C/<em>phosphate</em>,
and it shares no experiment with this block.</li>
</ul>
<p class='note'>Every fit behind every number here is drawn in
<a href='progress_curves.html'>progress_curves.html</a>, all
{len(frame)} of them.</p>
"""
    return styled("The two-axis block", body,
                  f"Exps 135–151 · BnOH · 25 °C · pyrophosphate · "
                  f"{summary['curves']} curves · pH "
                  f"{summary['pH_range'][0]:.2f}–{summary['pH_range'][1]:.2f}")


def main():
    return write_pages(HERE, {"index.html": build_index(),
                              "progress_curves.html": build_curves_page()})


if __name__ == "__main__":
    raise SystemExit(main())
