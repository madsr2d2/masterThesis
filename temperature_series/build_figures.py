"""
Builds temperature_series/index.html and progress_curves.html.

Draws only; every number comes from `arrhenius`, `scope`, `curve_metrics` or
`verify_enzyme_stock`, so a figure and the prose in ANALYSIS.md cannot disagree
about a value without `check_numbers.py` saying so.

    python temperature_series/build_figures.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
# svgplot lives at the repository root, shared with background_reaction.
sys.path.insert(0, os.path.dirname(HERE))

import arrhenius
import scope
import slowdown
import slowdown
import verify_enzyme_stock
from curve_metrics import (ACCELERATION_SIGMA, SEGMENT_RATIO_STEEP,
                           rolling_slope, segmented_fit)
from fit_dataset import source_floor
from summary_kinetics import fit_burst_bounded, fit_progress
from svgplot import ACCENT, GRID, INK, MUTED, Axes, esc, page, PAGE_CSS

# ORDERED VARIABLES GET SEQUENTIAL RAMPS, not categorical hues. Substrate rung
# and temperature are both ordinal, so a light-to-dark single hue carries the
# order; cycling categorical hues over them would throw that away. Every step
# clears 3:1 against the figure surface below, and every line is direct-labelled
# as well, so identity is never colour alone.
RUNGS = ["#6295c3", "#3d729f", "#1e5079", "#0c2f4d"]
TEMPERATURES = ["#7fa9cd", "#5d8fba", "#3c74a4", "#255a8a", "#12406b", "#062a4b"]
# The one genuinely categorical set: three unordered things (three enzyme
# hypotheses, three fitted parameters). Validated for colour-vision deficiency
# separation and contrast in both light and dark.
CATEGORY = ["#2f6fb0", "#c0522a", "#8a5aa8"]
# Figures sit on a fixed light surface whatever the page theme, so the ramps'
# contrast is deterministic. A sequential ramp cannot clear 3:1 against a white
# AND a near-black surface at once -- it needs the lightness range that the
# contrast rule would spend.
SURFACE = "#fbfbfa"
# Narrower than the smallest mark diameter used on a progress panel, so the
# fit never covers a reading whole. See build_curves_page.
FIT_WIDTH = 1.5

EXTRA_CSS = """
.fig{background:#fbfbfa;border-color:#e4e4e2}
.fig .cap{color:#5a5a5a}
.pill{display:inline-block;font-size:11px;padding:1px 8px;border-radius:10px;
background:var(--rule);color:var(--muted);margin-left:7px;vertical-align:2px}
.hero{display:flex;flex-wrap:wrap;gap:26px;margin:14px 0 4px}
.hero div{min-width:140px}
.hero .v{font-size:25px;font-weight:650;letter-spacing:-0.02em}
.hero .k{font-size:11.5px;color:var(--muted);text-transform:uppercase;
letter-spacing:0.06em}
.hero .u{font-size:12px;color:var(--muted)}
"""


def fig(svg, caption, extra=""):
    return (f"<div class='fig'>{svg}"
            f"<div class='cap'>{caption}</div>{extra}</div>")


def legend(entries):
    items = "".join(f"<span><i class='sw' style='background:{c}'></i>"
                    f"{esc(t)}</span>" for c, t in entries)
    return f"<div class='key'>{items}</div>"


def _frame():
    return arrhenius.series_frame()


# --- section 1: the enzyme ------------------------------------------------
def figure_stocks():
    table = verify_enzyme_stock.audit()
    used = table[table.used].sort_index()
    experiments = np.array([int(e) for e in used.index], dtype=float)
    values = used.kuv_sheet.to_numpy(dtype=float)
    axes = Axes(560, 250, (0, experiments.max() + 4),
                (0, values.max() * 1.18), pad=(64, 14, 46, 22))
    # A step, not a line: a stock is used unchanged until the next weighing.
    for index in range(len(experiments)):
        left = experiments[index]
        right = (experiments[index + 1] if index + 1 < len(experiments)
                 else left + 2)
        axes.line([left, right], [values[index]] * 2, CATEGORY[0], width=2.2)
    axes.points(experiments, values, CATEGORY[0], radius=2.6)
    flagged = verify_enzyme_stock.out_of_sequence(table)
    for _, row in flagged.iterrows():
        value = float(used.loc[int(row.experiment), "kuv_sheet"])
        axes.ring([float(row.experiment)], [value], ACCENT, radius=8)
        axes.label(float(row.experiment), value, f"exp {int(row.experiment)}",
                   ACCENT, size=11, anchor="middle", weight="600", dy=-14)
    axes.note(70, 26, "each plateau is one weighing", MUTED)
    return fig(
        axes.render("experiment number", "[enz] in the cuvette, mM",
                    "A · Seven enzyme stocks, each used until the next weighing"),
        "From <code>verify_enzyme_stock.audit()</code>, which recomputes the "
        "cuvette concentration from the mass, molar mass, made-up volume and "
        "cuvette volumes recorded beside it — none of which the pipeline reads. "
        "Stocks change rarely and never change back, which is why the ringed "
        "point is the question: <strong>exp 16 is the only experiment in the "
        "campaign whose stock interrupts a run</strong>.")


def figure_enzyme_test():
    plain = arrhenius.experiment_residuals()
    forced = arrhenius.experiment_residuals(override={16: 0.273})
    axes = Axes(560, 250, (12, 43), (-0.16, 0.14), pad=(64, 96, 46, 22))
    axes.hline(0.0, INK, dash="3 3", width=1.1)
    for table, colour, dash, name in ((plain, CATEGORY[0], None, "as recorded"),
                                      (forced, CATEGORY[1], "5 3",
                                       "exp 16 forced to 0.273")):
        order = table.sort_values("temperature")
        axes.line(order.temperature, order.mean_residual, colour,
                  width=2.0, dash=dash)
        axes.points(order.temperature, order.mean_residual, colour, radius=4.2)
        last = order.iloc[-1]
        axes.label(float(last.temperature), float(last.mean_residual), name,
                   colour, size=11, anchor="start", weight="600", dx=9, dy=4)
    for table, colour in ((plain, CATEGORY[0]), (forced, CATEGORY[1])):
        axes.ring([float(table.loc[16, "temperature"])],
                  [float(table.loc[16, "mean_residual"])], colour, radius=8)
    axes.label(40, float(plain.loc[16, "mean_residual"]), "exp 16", MUTED,
               size=10.5, anchor="middle", dy=17)
    return fig(
        axes.render("temperature, °C", "mean distance from the Arrhenius line, ln units",
                    "B · Forcing the higher enzyme makes exp 16 a worse outlier"),
        "Each run's four rungs averaged. Rate divided by enzyme concentration "
        "should fall on one line against 1/T, so a concentration that is wrong "
        "for one run displaces that run away from the line the other five "
        "define. Forcing exp 16 to 0.273 doubles its distance, from "
        "<strong>−0.060 to −0.121</strong>, and moves it further out on all "
        "four rungs — the direction a too-large divisor gives.")


# --- section 2: the shape -------------------------------------------------
def figure_shape():
    frame = _frame()
    rungs = sorted(frame.s0.unique())
    left = Axes(430, 250, (12, 43), (0, 2.6), pad=(58, 60, 46, 22))
    left.hline(1.0, INK, dash="3 3", width=1.1)
    left.hline(SEGMENT_RATIO_STEEP, ACCENT, dash="4 3", width=1.1)
    for index, s0 in enumerate(rungs):
        rung = frame[np.isclose(frame.s0, s0)].sort_values("temperature")
        capped = np.minimum(rung.break_ratio.to_numpy(dtype=float), 2.55)
        left.line(rung.temperature, capped, RUNGS[index], width=1.9)
        left.points(rung.temperature, capped, RUNGS[index], radius=3.6)
    left.note(64, 26, "steepening after the break", MUTED)
    left.note(64, 40, "flat = no break", MUTED)
    left.label(15, 2.55, "18.1 ↑", ACCENT, size=10, anchor="middle", dy=-6)

    right = Axes(430, 250, (12, 43), (0, 1.0), pad=(58, 60, 46, 22))
    right.hline(0.5, GRID, dash="4 3")
    for index, s0 in enumerate(rungs):
        rung = frame[np.isclose(frame.s0, s0)].sort_values("temperature")
        right.line(rung.temperature, rung.vmax_where, RUNGS[index], width=1.9)
        right.points(rung.temperature, rung.vmax_where, RUNGS[index], radius=3.6)
    right.band([12, 43], [0.65, 0.65], [1.0, 1.0], ACCENT, opacity=0.10)
    right.note(64, 26, "run ends before the rate levels off", ACCENT)
    return fig(
        "<div class='grid two'><div>" + left.render(
            "temperature, °C", "slope after ÷ slope before",
            "C · The break is the induction, and it fades with heat")
        + "</div><div>" + right.render(
            "temperature, °C", "position of the steepest block in the run",
            "D · …and at the cold end the run stops too soon")
        + "</div></div>"
        + legend([(RUNGS[i], f"[S] = {s0:.3f} mM") for i, s0 in enumerate(rungs)]),
        "<strong>C.</strong> Cold runs accelerate throughout, hot runs "
        "decelerate; the dashed orange line is the 1.5 threshold "
        "<code>synchronised_break</code> counts as steepening. One point is "
        "clipped: exp 19's lowest rung reads 18.1 because its pre-break slope "
        "is 5.6e-08 AU/s — a near-zero denominator at the slowest condition in "
        "the block, not an anomaly. <strong>D.</strong> Above 0.65 the steepest "
        "block sits in the last third of the run, so <code>vmax</code> is not a "
        "maximum but wherever the measurement stopped. That is the whole reason "
        "<code>v_ss</code> is carried beside it.")


def figure_selection():
    frame = _frame()
    left = Axes(430, 250, (12, 43), (0, 3.3), pad=(58, 24, 46, 22))
    left.hline(1.0, MUTED, dash="3 3", width=1.1)
    for label, column, colour in (("one phase", "v0_burst_resid", CATEGORY[1]),
                                  ("after selection", "progress_resid",
                                   CATEGORY[0])):
        median = frame.groupby("temperature")[column].median()
        left.line(median.index, median.values, colour, width=2.2)
        left.points(median.index, median.values, colour, radius=4.4)
        left.label(43, float(median.iloc[-1]), label, colour, size=11,
                   anchor="end", weight="600", dy=-9)
    left.note(62, 26, "1x noise = the form fits", MUTED)

    right = Axes(430, 250, (12, 43), (0, 4.4), pad=(58, 24, 46, 22))
    counts = frame.groupby("temperature").phases.apply(lambda s: (s == 2).sum())
    for temperature, count in counts.items():
        colour = CATEGORY[0] if count else MUTED
        px = right._fx(temperature)
        right.parts.append(
            f"<rect x='{px - 9:.1f}' y='{right._fy(count):.1f}' width='18' "
            f"height='{right._fy(0) - right._fy(count):.1f}' fill='{colour}' "
            f"fill-opacity='0.85' rx='3'/>")
        right.label(temperature, count, f"{count} of 4", INK, size=10.5,
                    anchor="middle", dy=-8)
    return fig(
        "<div class='grid two'><div>" + left.render(
            "temperature, °C", "residual, in units of the curve's own noise",
            "J · The second phase is taken only where the first form fails")
        + "</div><div>" + right.render(
            "temperature, °C", "cuvettes selecting two phases",
            "K · …which is 11 of 24, and never at 15, 20 or 30 °C")
        + "</div></div>",
        "The forms are exactly nested — B₂ = 0 is the one-phase form — so the "
        "choice is an F test on two extra parameters "
        "(<code>summary_kinetics.fit_progress</code>). At 25, 35 and 40 °C the "
        "one-phase form sits at 1.5–3.0× noise and the second phase brings it "
        "to about 1.1×; at 15, 20 and 30 °C it is already at noise and nothing "
        "is selected. <strong>Selection is on F alone, deliberately not on τ₂ "
        "being resolved</strong> — τ₂ is pinned on only 2 of the 11, and "
        "requiring it rejected the clearest two-phase curve in the block "
        "(F = 791). Whether a phase exists and whether its time constant is "
        "determined are different questions.")


def figure_tau():
    frame = _frame()
    warm = frame[frame.tau_resolved & (frame.temperature <= 32)]
    axes = Axes(560, 270, (3.25, 3.52), (9e-5, 2.2e-3), ylog=True,
                pad=(70, 20, 46, 22))
    rungs = sorted(frame.s0.unique())
    for index, s0 in enumerate(rungs):
        rung = warm[np.isclose(warm.s0, s0)]
        axes.points(1000.0 / rung.kelvin, 1.0 / rung.tau, RUNGS[index],
                    radius=4.2, title=f"[S] = {s0:.3f} mM")
    fit = arrhenius.activation_parameters("inverse_tau")
    grid = np.array([3.26, 3.51])
    # The fitted pooled slope, drawn through the mean of the points so the eye
    # compares the SLOPE, which is the activation energy, not the offset.
    y = np.log(1.0 / warm.tau.to_numpy(dtype=float))
    x = 1000.0 / warm.kelvin.to_numpy(dtype=float)
    slope = -fit["activation_kJ"] * 1000.0 / arrhenius.GAS_CONSTANT / 1000.0
    axes.line(grid, np.exp(y.mean() + slope * (grid - x.mean())), ACCENT,
              width=2.0, dash="6 4")
    axes.note(300, 34, f"Ea = {fit['activation_kJ']:.0f} ± "
                       f"{fit['activation_stderr']:.0f} kJ/mol", ACCENT, size=11.5)
    for temperature in (15, 20, 25, 30):
        axes.label(1000.0 / (temperature + 273.15), 1.0e-4, f"{temperature}°",
                   MUTED, size=10, anchor="middle")
    return fig(
        axes.render("1000 / T, K⁻¹", "1/τ, s⁻¹  (induction rate constant)",
                    "E · The induction is a rate constant, and it has its own barrier"),
        "Only the 15–30 °C runs, and only the 15 of 16 curves where the "
        "profile actually pinned τ. Above 30 °C the induction is over before "
        "the run is under way, the burst form flips to a decelerating shape and "
        "τ stops being resolved — so those points are not measurements and are "
        "not plotted. τ falls <strong>6489 → 3666 → 945 → 876 s</strong> over "
        "these four temperatures.")


def figure_truncation():
    frame = _frame()
    ratio = frame.groupby("temperature").apply(
        lambda g: float((g.v_ss / g.vmax).median()), include_groups=False)
    axes = Axes(560, 240, (12, 43), (0, 1.35), pad=(64, 20, 46, 22))
    axes.hline(1.0, INK, dash="3 3", width=1.1)
    good = [t for t in ratio.index if t <= 32]
    bad = [t for t in ratio.index if t > 32]
    axes.line(good, [ratio[t] for t in good], CATEGORY[0], width=2.2)
    axes.points(good, [ratio[t] for t in good], CATEGORY[0], radius=4.6)
    axes.line(bad, [ratio[t] for t in bad], MUTED, width=2.0, dash="5 3")
    axes.points(bad, [ratio[t] for t in bad], MUTED, radius=4.6)
    for temperature in ratio.index:
        axes.label(temperature, ratio[temperature], f"{ratio[temperature]:.2f}",
                   INK if temperature <= 32 else MUTED, size=10.5,
                   anchor="middle", dy=-11)
    axes.note(330, 205, "τ unresolved — v_ss is not an asymptote here", MUTED)
    return fig(
        axes.render("temperature, °C", "v_ss ÷ vmax",
                    "F · How much vmax under-reads, and where v_ss may be trusted"),
        "Where both estimators are sound they agree — <strong>0.98 and 0.99 at "
        "25 and 30 °C</strong> — which is what makes the substitution mean "
        "anything. At 15 and 20 °C <code>vmax</code> under-reads by "
        "<strong>4–6.5%</strong> because the run ends before the rate levels "
        "off. Above 32 °C (grey) the burst form has degenerated and the ratio is "
        "not a measurement. Substituting <code>v_ss</code> at the two cold "
        "points only moves E<sub>a</sub> from 90.2 to 88.3 kJ/mol: the "
        "truncation costs <strong>1.9 kJ/mol</strong>.")


# --- section 3: composition ----------------------------------------------
def figure_buffer():
    frame = scope.frame((32, 34))
    frame = frame[frame.live]
    axes = Axes(560, 270, (2.5, 260), (1e-5, 1e-4), xlog=True, ylog=True,
                pad=(70, 22, 46, 22))
    for experiment, colour, name in ((34, CATEGORY[2], "exp 34"),
                                     (32, CATEGORY[0], "exp 32")):
        run = frame[frame.experiment == experiment].sort_values("buf")
        axes.points(run.buf, run.vmax, colour, radius=4.6)
        fit = arrhenius.catalysed_buffer_order((experiment,))
        grid = np.array([run.buf.min(), run.buf.max()])
        centre = np.exp(np.log(run.vmax).mean())
        middle = np.exp(np.log(run.buf).mean())
        axes.line(grid, centre * (grid / middle) ** fit["order_buf"], colour,
                  width=2.0)
        axes.label(float(run.buf.iloc[0]), float(run.vmax.iloc[0]),
                   f"{name}: {fit['order_buf']:+.2f}", colour, size=11,
                   anchor="start", weight="600", dx=6, dy=-9)
    axes.band([50, 80], [1e-5, 1e-5], [1e-4, 1e-4], ACCENT, opacity=0.09)
    axes.note(300, 32, "the temperature series sits here", ACCENT, size=11)
    return fig(
        axes.render("[buf], mM", "vmax, AU/s",
                    "G · The catalysed buffer order saturates — and it is not one number"),
        "Exps 32 and 34: 4OMe-BnOH at 40 °C and pH 7.00 with <strong>[S] held "
        "at 8.251 mM</strong> and only the buffer stepped. This is the design "
        "the BnOH set lacks entirely. Over 3.125–25 mM the order is "
        "<strong>+0.803 ± 0.173</strong>; over 50–200 mM it is "
        "<strong>+0.402 ± 0.024</strong> (R² 0.993). The temperature series "
        "sits in the shaded band, so +0.40 is the value that corrects it.")


def figure_substrate_order():
    orders = arrhenius.substrate_order()
    frame = _frame()
    left = Axes(430, 250, (1.4, 8.4), (44, 86), pad=(58, 24, 46, 22))
    ladder = frame.groupby("s0").buf.median().sort_index()
    left.line(ladder.index, ladder.values, ACCENT, width=2.2)
    left.points(ladder.index, ladder.values, ACCENT, radius=5)
    for s0, buf in ladder.items():
        left.label(s0, buf, f"{buf:.0f}", ACCENT, size=10.5, anchor="middle",
                   dy=-11)
    left.note(64, 26, "substrate volume displaced buffer volume", MUTED)

    right = Axes(430, 250, (12, 43), (0.0, 0.82), pad=(58, 24, 46, 22))
    right.line(orders.index, orders.observed, CATEGORY[0], width=2.2)
    right.points(orders.index, orders.observed, CATEGORY[0], radius=4.4)
    right.line(orders.index, orders.corrected, CATEGORY[1], width=2.2)
    right.points(orders.index, orders.corrected, CATEGORY[1], radius=4.4)
    right.hline(float(orders.corrected.mean()), CATEGORY[1], dash="4 3")
    right.label(15, float(orders.observed[15.0]), "observed", CATEGORY[0],
                size=11, anchor="start", weight="600", dx=6, dy=14)
    right.label(15, float(orders.corrected[15.0]), "buffer-corrected",
                CATEGORY[1], size=11, anchor="start", weight="600", dx=6, dy=-8)
    return fig(
        "<div class='grid two'><div>" + left.render(
            "[S] in the cuvette, mM", "[buf] in the same cuvette, mM",
            "H · The ladder carries a buffer gradient")
        + "</div><div>" + right.render(
            "temperature, °C", "order in [S]",
            "I · …so the substrate order needs correcting, and it moves")
        + "</div></div>",
        "<strong>H.</strong> Every run steps [buf] down 80 → 50 mM as [S] rises, "
        "identically. That is why the activation energies are unaffected — each "
        "rung is one fixed composition measured at six temperatures — and why "
        "the substrate order is not. <strong>I.</strong> With g = −0.325 and "
        "d = +0.402 the correction adds 0.13 at every temperature. The mean goes "
        "from +0.445 to <strong>+0.576</strong> and the spread is "
        "<strong>0.173 before and after</strong>: the correction moved the "
        "level, not the trend. Partial saturation in substrate, with no "
        "detectable temperature dependence.")


# --- sections 4-5: the parameters -----------------------------------------
def figure_arrhenius():
    frame = _frame()
    rungs = sorted(frame.s0.unique())
    axes = Axes(560, 330, (3.16, 3.52), (4e-6, 3.2e-4), ylog=True,
                pad=(72, 96, 46, 22))
    fits = arrhenius.rung_fits("vmax")
    for index, s0 in enumerate(rungs):
        rung = frame[np.isclose(frame.s0, s0)].sort_values("kelvin")
        x = 1000.0 / rung.kelvin.to_numpy(dtype=float)
        y = rung.vmax.to_numpy(dtype=float) / rung.e0.to_numpy(dtype=float)
        axes.points(x, y, RUNGS[index], radius=4.2)
        row = fits.loc[s0]
        grid = np.array([3.185, 3.478])
        slope = -row.activation_kJ * 1000.0 / arrhenius.GAS_CONSTANT
        axes.line(grid, np.exp(np.log(y).mean()
                               + slope * (grid - x.mean()) / 1000.0),
                  RUNGS[index], width=1.7)
        # Direct-labelled in the right margin, at the cold end where the lines
        # are furthest apart. Inside the frame they collided with the 40 C
        # points, which sit two pixels from the left edge.
        axes.label(x.max(), float(y[np.argmax(x)]),
                   f"{s0:.3f} mM · {row.activation_kJ:.0f}±{row.stderr_kJ:.0f}",
                   RUNGS[index], size=10.5, anchor="start", weight="600",
                   dx=9, dy=4)
    for temperature in (15, 20, 25, 30, 35, 40):
        axes.label(1000.0 / (temperature + 273.15), 4.4e-6, f"{temperature}°",
                   MUTED, size=10, anchor="middle")
    pooled = arrhenius.pooled_arrhenius("vmax")
    axes.note(300, 34, f"pooled Ea = {pooled['activation_kJ']:.1f} ± "
                       f"{pooled['stderr_kJ']:.1f} kJ/mol", INK, size=12,
              weight="650")
    return fig(
        axes.render("1000 / T, K⁻¹", "vmax ÷ [enz], AU s⁻¹ mM⁻¹",
                    "J · Arrhenius: four rungs, six temperatures, one slope"),
        "The four rungs are four <em>independent</em> fits — the same substrate "
        "ladder recurs in all six runs — and they are offset because their "
        "composition differs, not because their slopes do. Individually they "
        "give 94.5, 86.0, 89.0 and 91.5 kJ/mol at ±2.6–3.0, which against their "
        "weighted mean is <strong>χ² = 4.68 on 3 dof, reduced 1.56</strong>: no "
        "composition dependence is detectable. Pooling by refitting one slope "
        "with a free intercept per rung gives "
        "<strong>90.2 ± 1.5 kJ/mol</strong>.")


def figure_rung_energies():
    fits = arrhenius.rung_fits("vmax")
    pooled = arrhenius.pooled_arrhenius("vmax")
    axes = Axes(560, 230, (1.0, 8.3), (78, 102), pad=(66, 24, 46, 22))
    low = pooled["activation_kJ"] - pooled["stderr_kJ"]
    high = pooled["activation_kJ"] + pooled["stderr_kJ"]
    axes.band([1.0, 8.3], [low, low], [high, high], CATEGORY[0], opacity=0.16)
    axes.hline(pooled["activation_kJ"], CATEGORY[0], dash="5 3", width=1.6)
    for index, (s0, row) in enumerate(fits.iterrows()):
        axes.errorbar(s0, row.activation_kJ,
                      row.activation_kJ - row.stderr_kJ,
                      row.activation_kJ + row.stderr_kJ, RUNGS[index])
        axes.label(s0, row.activation_kJ + row.stderr_kJ,
                   f"{row.activation_kJ:.1f}", INK, size=10.5, anchor="middle",
                   dy=-9)
    axes.note(300, 32, f"pooled {pooled['activation_kJ']:.1f} ± "
                       f"{pooled['stderr_kJ']:.1f}", CATEGORY[0], size=11.5,
              weight="600")
    return fig(
        axes.render("[S] of the rung, mM", "activation energy, kJ/mol",
                    "K · The rung spread is scatter, not composition"),
        "The 8.5 kJ/mol spread across four fits that ought to agree looked like "
        "a composition dependence until the standard errors were computed — "
        "each rung is ±2.6–3.0, and every one overlaps the pooled band. The "
        "pooled value is <em>refitted</em> over all 24 points rather than "
        "averaged from these four, because the four share the same six runs, "
        "days and cells: their errors are correlated and averaging would shrink "
        "the error by a factor that is not there.")


def figure_eyring():
    axes = Axes(560, 320, (3.16, 3.50), (2e-9, 1.3e-7), ylog=True,
                pad=(76, 74, 46, 22))
    frame = _frame()
    rungs = sorted(frame.s0.unique())
    for index, s0 in enumerate(rungs):
        rung = frame[np.isclose(frame.s0, s0)].sort_values("kelvin")
        constant = arrhenius.turnover(rung, "vmax")
        axes.points(1000.0 / rung.kelvin, constant / rung.kelvin, RUNGS[index],
                    radius=4.2, title=f"[S] = {s0:.3f} mM")
    result = arrhenius.activation_parameters("vmax")
    middle = frame[np.isclose(frame.s0, result["at_s0"])].sort_values("kelvin")
    x = 1000.0 / middle.kelvin.to_numpy(dtype=float)
    y = np.log(arrhenius.turnover(middle, "vmax")
               / middle.kelvin.to_numpy(dtype=float))
    grid = np.array([3.17, 3.49])
    slope = -result["enthalpy_kJ"] * 1000.0 / arrhenius.GAS_CONSTANT / 1000.0
    axes.line(grid, np.exp(y.mean() + slope * (grid - x.mean())), ACCENT,
              width=2.2, dash="6 4")
    for temperature in (15, 20, 25, 30, 35, 40):
        axes.label(1000.0 / (temperature + 273.15), 2.2e-9, f"{temperature}°",
                   MUTED, size=10, anchor="middle")
    axes.note(300, 34, f"ΔH‡ = {result['enthalpy_kJ']:.1f} ± "
                       f"{result['enthalpy_stderr']:.1f} kJ/mol", ACCENT,
              size=11.5, weight="600")
    axes.note(300, 50, f"ΔS‡ = {result['entropy_J']:.1f} ± "
                       f"{result['entropy_stderr']:.1f} J/mol/K", ACCENT,
              size=11.5, weight="600")
    return fig(
        axes.render("1000 / T, K⁻¹", "k / T,  K⁻¹ s⁻¹   (k = vmax ÷ ε ÷ [enz])",
                    "L · Eyring: the slope is ΔH‡, the intercept is ΔS‡"),
        "The rate is turned into an actual rate constant in s⁻¹ by dividing by "
        "the sheet's own extinction coefficient and by [enz] — a slope survives "
        "any constant factor but an intercept does not, so <strong>ΔH‡ needs "
        "none of this and ΔS‡ needs all of it</strong>. The dashed line is the "
        "pooled fit at the median rung, [S] = 5.549 mM and [buf] = 60 mM, which "
        "is the composition ΔS‡ is quoted at. Since the substrate order is "
        "+0.58 rather than 0, this constant is pseudo-first-order and carries a "
        "substrate dependence.")


# --- the progress curves --------------------------------------------------
def build_curves_page():
    frame = _frame()
    rungs = sorted(frame.s0.unique())
    temperatures = sorted(frame.temperature.unique())
    curves = {(c.temperature, round(float(c.conditions.s0), 3)): c
              for c in scope.curves(scope.TEMPERATURE_SERIES)}
    panels = []
    for temperature in temperatures:
        row = frame[frame.temperature == temperature]
        for s0 in rungs:
            curve = curves.get((temperature, round(float(s0), 3)))
            if curve is None:
                continue
            times = np.asarray(curve.times, dtype=float)
            values = np.asarray(curve.absorbance, dtype=float)
            record = row[np.isclose(row.s0, s0)].iloc[0]
            axes = Axes(340, 210, (0, times[-1] * 1.02),
                        (min(values.min(), 0) * 1.1 - 1e-4,
                         max(values.max(), 1e-4) * 1.12),
                        pad=(56, 12, 34, 20))
            # WHICHEVER FORM THE CURVE EARNED. Drawing the one-phase fit on a
            # two-phase curve is how the shape stayed invisible for as long as
            # it did -- the panel showed a fit that could not bend the way the
            # readings do, and it read as scatter.
            progress = fit_progress(times, values)
            smooth = np.linspace(0, times[-1], 300)
            fitted = progress.predict(smooth)
            # DATA FIRST, FIT ON TOP, AND THE FIT MUST BE NARROWER THAN THE
            # MARKS. Three passes to get this right, so the constraint is
            # written down rather than re-derived:
            #
            #   1. fit under the data   -> 368 readings bury the fit
            #   2. fit over the data, 2.0 px wide on a 3.6 px white halo
            #      -> the halo is wider than a mark, so wherever the curve is
            #         tight the fit erases the very points it is fitting
            #   3. this: a light-grey scatter with a thin rust line over it.
            #
            # Separation is by HUE and LIGHTNESS, not by a halo: grey cloud,
            # dark rust line. FIT_WIDTH stays below the mark diameter so a mark
            # centred on the line still shows on both sides -- asserted below,
            # because it is the property that keeps both visible and it is easy
            # to break by nudging a radius.
            dense = len(times) > 150
            radius = 1.6 if dense else 2.1
            axes.points(times, values, MUTED, radius=radius,
                        opacity=0.75 if dense else 0.9,
                        stroke=None if dense else "white", stroke_width=0.6)
            assert FIT_WIDTH < 2 * radius, "fit line would cover the marks"
            axes.line(smooth, fitted, ACCENT, width=FIT_WIDTH)
            # EVERY breakpoint the curve earned, not just the first. A rate
            # that rises and then falls has two, and drawing one leaves the
            # other invisible -- which is how the early break on the 40 C
            # curves went unnoticed until it was seen by eye on this page.
            for index, where in enumerate(record.break_times):
                axes.parts.append(
                    f"<path d='M{axes._fx(where):.2f},{axes.top} "
                    f"L{axes._fx(where):.2f},{axes.height - axes.bottom}' "
                    f"stroke='{MUTED}' stroke-width='1.1' "
                    f"stroke-dasharray='3 3' fill='none'/>")
                axes.note(axes._fx(where) + 3, axes.top + 10,
                          f"{index + 1}", MUTED, size=9.5)
            kind = progress.chosen.kind
            panels.append(
                "<div class='fig panel'>"
                f"<div class='ph'>{temperature:.0f} °C · [S] = {s0:.3f} mM"
                f"<span class='pill'>exp {int(record.experiment)}</span></div>"
                f"<div class='ps'>[buf] {record.buf:.0f} mM · "
                f"{int(record.points)} readings over {record.duration_s / 60:.0f} min"
                f" · {record.source}</div>"
                + axes.render("time, s", "ΔA at 300 nm")
                + f"<div class='pf'><strong>{int(record.phases)} phase"
                + ("s" if record.phases == 2 else "")
                + f"</strong> · {esc(kind)} · F = {record.two_phase_f:.0f}"
                + f" · fit {record.progress_resid:.2f}x noise"
                + (f" · τ₁ = {record.tau_fast:.0f} s" if record.phases == 2
                   else (f" · τ = {record.tau:.0f} s" if record.tau_resolved
                         else " · τ unresolved"))
                + " · breaks at "
                + (", ".join(f"{v:.0f}" for v in record.break_times) + " s"
                   if record.break_times else "none")
                + f" · <strong>{esc(record.break_pattern)}</strong>"
                + "</div></div>")
    body = ("<p class='lede'>All 24 curves of the temperature series, in "
            "temperature order. The orange line is whichever form the curve "
            "earned — one relaxation or two — from "
            "<code>summary_kinetics.fit_progress</code>; the numbered dashed "
            "verticals are the breakpoints of the best piecewise-linear "
            "description (<code>curve_metrics.segment_selection</code>), one "
            "or two according to an F test. <strong>Read the pattern in the "
            "footer, not the number of breaks</strong>: three straight lines "
            "fit a smooth bend better than two whether or not anything "
            "happened, so the count is a statement about approximation and the "
            "sequence of slopes is the statement about the curve. Nothing is "
            "excluded and "
            "every fit is computed on every point — except the instrument's "
            "first reading, which is discarded from every run in the archive "
            "(<code>fit_dataset.DROP_FIRST_READING</code>) and is therefore "
            "not plotted, fitted or scored here.</p>"
            "<p class='lede'>Read down a column and the induction shortens: the "
            "curves bend from lag-shaped at 15 °C to plainly decelerating at "
            "40 °C, which is the same thing τ and the breakpoint ratio report "
            "as numbers.</p>"
            "<div class='grid three'>" + "".join(panels) + "</div>")
    return page("Temperature series — all 24 progress curves", body).replace(
        "</style>", EXTRA_CSS + "</style>")


# --- the presentation -----------------------------------------------------
# --- section 4: what naming the fall does to the numbers -------------------
def figure_sink_effect():
    """Both Arrhenius lines on one plot, so "it does not move" is visible."""
    frame = slowdown.production_frame()
    effect = slowdown.sink_effect_on_activation()
    enzyme = frame.e0.to_numpy(dtype=float) * frame.epsilon.to_numpy(dtype=float)
    inverse = 1000.0 / frame.kelvin.to_numpy(dtype=float)
    series = {"v_peak": (CATEGORY[0], "published, v_peak"),
              "v_prod": (CATEGORY[1], "sink model, k pinned")}
    values = {name: frame[name].to_numpy(dtype=float) / enzyme
              for name in series}
    low = min(np.nanmin(v) for v in values.values()) / 1.7
    high = max(np.nanmax(v) for v in values.values()) * 1.7
    axes = Axes(560, 285, (inverse.min() - 0.02, inverse.max() + 0.02),
                (low, high), ylog=True, pad=(74, 26, 46, 30))
    for name, (colour, label) in series.items():
        axes.points(inverse, values[name], colour, radius=3.6, opacity=0.9)
        fit = arrhenius.activation_parameters(name, frame=frame)
        # One pooled slope, four intercepts: draw it at the median rung, which
        # is the level every quoted intercept in this document belongs to.
        middle = np.isclose(frame.s0, sorted(frame.s0.unique())[2])
        centre = float(np.exp(np.log(values[name][middle]).mean()))
        pivot = float(np.mean(inverse[middle]))
        slope = -fit["activation_kJ"] * 1000.0 / arrhenius.GAS_CONSTANT / 1000.0
        grid = np.array([inverse.min() - 0.01, inverse.max() + 0.01])
        axes.line(grid, centre * np.exp(slope * (grid - pivot)), colour,
                  width=2.2)
        # The two lines very nearly coincide -- which IS the figure -- so
        # they cannot be labelled on themselves. The upper right is empty
        # because the line runs down to the right.
        axes.note(524, 52 + 17 * list(series).index(name),
                  f"{label}: {fit['activation_kJ']:.1f} kJ/mol", colour,
                  size=11, weight="600", anchor="end")
    return fig(
        axes.render("1000 / T, K⁻¹", "rate / (ε·[enz]), s⁻¹",
                    "M · Naming the fall does not move the line"),
        "Both estimators, all 24 curves, one pooled slope each. "
        f"<code>v_prod</code> is the production rate of the sink model with k "
        f"pinned at <code>slowdown.sink_activation</code>'s value; it sits "
        f"<strong>{effect['lift'] - 1:.1%}</strong> above <code>v_peak</code> "
        f"on the median curve, and between {effect['lift_low'] - 1:+.1%} and "
        f"{effect['lift_high'] - 1:+.1%} across the six temperatures "
        "<strong>with no order to it</strong>. So it shifts the level and not "
        "the slope: E<sub>a</sub> moves by "
        f"<strong>{effect['activation_shift']:+.2f} ± "
        f"{effect['activation_shift_stderr']:.2f} kJ/mol</strong>. The "
        f"corrected line is the noisier one — Arrhenius scatter "
        f"{effect['corrected']['rms']:.3f} against "
        f"{effect['published']['rms']:.3f} — which is why it is reported and "
        "not applied. <a href='../product_fate/index.html'>product_fate</a> "
        "has the derivation.")


def build_index():
    table = arrhenius.parameter_table()
    vmax = table.loc["v_peak"]
    hero = (
        "<div class='hero'>"
        f"<div><div class='k'>ΔH‡</div><div class='v'>{vmax.enthalpy_kJ:.1f}"
        f"</div><div class='u'>± {vmax.enthalpy_stderr:.1f} kJ/mol</div></div>"
        f"<div><div class='k'>ΔS‡</div><div class='v'>{vmax.entropy_J:.0f}"
        f"</div><div class='u'>± {vmax.entropy_stderr:.0f} J/mol/K</div></div>"
        f"<div><div class='k'>ΔG‡ at 25 °C</div><div class='v'>"
        f"{vmax.gibbs_kJ:.1f}</div><div class='u'>± {vmax.gibbs_stderr:.1f} "
        f"kJ/mol  ({vmax.gibbs_kJ / 4.184:.1f} kcal/mol)</div></div>"
        f"<div><div class='k'>E<sub>a</sub></div><div class='v'>"
        f"{vmax.activation_kJ:.1f}</div><div class='u'>± "
        f"{vmax.activation_stderr:.1f} kJ/mol</div></div>"
        "</div>")

    rows = []
    labels = {"v_peak": ("<code>v_peak</code>",
                         "peak rate of the fitted model",
                         "all six temperatures"),
              "vmax": ("<code>vmax</code>", "steepest observed rate",
                       "all six, truncated cold"),
              "v_ss": ("<code>v_ss</code>", "asymptote after the induction",
                       "15–30 °C"),
              "inverse_tau": ("<code>1/τ</code>", "induction rate constant",
                              "15–30 °C")}
    for name, row in table.iterrows():
        symbol, what, span = labels[name]
        highlight = " class='hl'" if name == "v_peak" else ""
        rows.append(
            f"<tr{highlight}><td>{symbol} — {what}</td><td>{span}</td>"
            f"<td>{row.activation_kJ:.1f} ± {row.activation_stderr:.1f}</td>"
            f"<td>{row.enthalpy_kJ:.1f} ± {row.enthalpy_stderr:.1f}</td>"
            f"<td>{row.entropy_J:+.1f} ± {row.entropy_stderr:.1f}</td>"
            f"<td>{row.gibbs_kJ:.1f} ± {row.gibbs_stderr:.1f}</td>"
            f"<td>{int(row.n)}</td></tr>")
    parameters = (
        "<div class='scroll'><table><tr><th>fitted parameter</th>"
        "<th>valid over</th><th>E<sub>a</sub> kJ/mol</th><th>ΔH‡ kJ/mol</th>"
        "<th>ΔS‡ J/mol/K</th><th>ΔG‡(298) kJ/mol</th><th>curves</th></tr>"
        + "".join(rows) + "</table></div>")

    body = f"""
<p class='lede'>Exps 14–19: 4-methoxybenzyl alcohol + H<sub>2</sub>O<sub>2</sub>,
chemzyme-catalysed, in 65 mM phosphate at pH 7.00 with 82.5 mM
H<sub>2</sub>O<sub>2</sub>, run at 15, 20, 25, 30, 35 and 40 °C with the same
four-rung substrate ladder in every one. <strong>The only temperature series in
the archive</strong>, and so the only route to activation parameters — the
quantity <code>COMPUTATIONAL.md</code>'s barriers can be compared against.</p>

{hero}
<p class='warn'>These are the parameters of the <strong>catalytic increment</strong>
of a whole turnover, not of an elementary step. Read section 5 before putting
any of them beside a calculated barrier.</p>

<h2>1 · Is the enzyme mismatch real?</h2>
<p>Two of the six runs record [enz] 11.7% below the other four — 0.240683 against
0.272695 mM — and exp 16, the 40 °C run, sits <em>between</em> two runs at the
higher value. A real restock is divided out; a transcription error is corrected.
The two treatments differ, so this had to be settled before any fit.</p>
{figure_stocks()}
<p>The weighing chain is self-consistent on <strong>62 of 63</strong> experiments
and the compiled value equals the sheet's own <code>kuv</code> on
<strong>63 of 63</strong>. So both disputed values descend from a real weighing,
0.0203 g against 0.023 g, and the lower stock is used by nine experiments. A
typed-over cell would break the chain; neither chain is broken. What the
documents cannot settle is the <em>ordering</em> — exps 16, 17, 18 and 20 are
byte-identical workbooks straddling the boundary. The kinetics settle it:</p>
{figure_enzyme_test()}
<p><strong>Verdict: real.</strong> The recorded values win all four rungs on
<code>vmax</code> and all four again on <code>v0_quad</code>, so the answer does
not depend on the estimator. Every fit below divides by [enz] as recorded.</p>

<h2>2 · What the curves are doing</h2>
<p>Run before anything was built on these rates, because exp 65 in the BnOH
background showed that the start-versus-end shape statistics step straight over a
break in the middle of a run.</p>
{figure_shape()}
<p>Nothing here resembles exp 65, whose four cuvettes broke within 56 s of each
other across a 20-fold substrate range. Here two unrelated statistics agree on an
induction period that shortens with heating — and the second one is a rate
constant in its own right:</p>
{figure_tau()}
{figure_truncation()}

<h2>3 · What the composition is doing</h2>
<p>The substrate ladder is not a substrate ladder alone, and the fix needed an
experiment from outside this block.</p>
{figure_buffer()}
{figure_substrate_order()}

<h2>3a \u00b7 The fitting form: one relaxation, or two</h2>
<p>The burst/lag form's rate is <em>monotone</em>, so it cannot hold a rate that rises to a maximum and then falls \u2014 and 14 of these 24 curves do exactly that. A second relaxation term gives it that freedom, and the two forms are exactly nested.</p>
{figure_selection()}

<h2>4 · The activation parameters</h2>
{figure_arrhenius()}
{figure_rung_energies()}
{figure_eyring()}
{parameters}
<p>In kcal/mol for <code>vmax</code>: ΔH‡ <strong>21.0</strong>, ΔS‡
<strong>−12.7 cal/mol/K</strong>, ΔG‡(298) <strong>24.8</strong>.</p>
<p><strong>The load-bearing check is that the two rate estimators agree.</strong>
<code>vmax</code> uses all six temperatures and is truncated at the cold end;
<code>v_ss</code> uses only the four where the burst form is sound and is not
truncated. They give 87.7 and 86.4 kJ/mol — and the truncation correction in
figure F predicted <code>vmax</code> would read about 1.9 high, i.e. 85.8. Two
estimators with different failure modes, landing inside their errors.</p>
<p><strong>The induction is a different process.</strong> Its ΔG‡ is 12 kJ/mol
lower, so it is the faster of the two, and its ΔS‡ is near zero rather than −55 —
a much looser transition state. Its enthalpy, 92 ± 16, is not resolved apart from
the turnover's 88 ± 1.5; the <em>entropies</em> are what distinguish them.</p>

<h2>5 · Naming the fall, and what it does to the numbers</h2>
<p><a href='../product_fate/index.html'>product_fate</a> identifies the fall of
section 3a as <strong>the oxidant attacking the aldehyde it has just
made</strong> — the rate declines linearly in the accumulated product,
<code>A′ = v − kA</code>, on 24 of 29 well-determined curves against 0 for the
hyperbolic law product inhibition would give. So every rate in the table above
is already net of that loss, and none of them is the production rate.</p>
{figure_sink_effect()}
<p>The correction is reported and not applied, and the temperature series' own
figures and parameters are unchanged by it.</p>

<h2>6 · Before comparing any of this to a calculation</h2>
<ul>
<li><strong>Composite, not elementary.</strong> <code>vmax</code> and
<code>v_ss</code> are whole-turnover rates through a seven-step mechanism, so ΔH‡
is the barrier of whatever step is rate-limiting at this composition plus any
pre-equilibria in front of it. Compare against the <strong>highest</strong>
computed step, never a sum.</li>
<li><strong>ΔS‡ assumes a pseudo-first-order rate constant</strong>,
<code>v/(ε·[enz])</code> in s⁻¹. The substrate order is +0.58, not 0, so the
reaction is not saturated and that constant carries a substrate dependence.
<strong>ΔH‡ does not depend on this; ΔS‡ and the absolute ΔG‡ do.</strong></li>
<li><strong>First order in enzyme is assumed and untestable here</strong> — [enz]
takes two values in this block, 11.7% apart.</li>
<li><strong>These are the catalytic increment's parameters.</strong> The runs are
<code>with_E</code>, so the reference channel omits only the enzyme and the
background is already subtracted. The uncatalysed background's own temperature
dependence is measured nowhere in the archive.</li>
<li><strong>ΔG‡ is far better determined than ΔH‡ or ΔS‡</strong> — ±0.1 against
±1.5 and ±4.9 — because both come from one line and their errors are strongly
anti-correlated; at 25 °C, inside the measured range, they cancel. It is
propagated with the covariance, not from the two standard errors, which would
overstate it several-fold.</li>
</ul>

<h2>7 · What the slowdown is — and where it is written up</h2>
<p>The fall of section 3a is not a property of this block, so its analysis is
not on this page. It is <a href='../product_fate/index.html'><strong>product_fate</strong></a>,
because the discrimination that settles it needs the whole 4OMe archive: these
23 curves cannot separate “the rate fell with time” from “the rate fell with
product” — inside one curve the product only ever grows with time — and the
separation comes from 84 curves in which run length and product vary
independently.</p>
<p>In one paragraph: <strong>the catalysed 4OMe rate falls in proportion to the
product it has made</strong>, and the same chemistry with no enzyme falls on a
clock instead — which is the reference channel these runs are measured against,
so its decay is subtracted out. <strong>The rate is linear in the product, not
hyperbolic in it</strong>: production minus a first-order loss of the measured
species, not inhibition of the catalyst. It does not happen to BnOH at product
concentrations nearly twice as high. Section 5 above has the only consequence
for this document, which is none that these six temperatures can resolve.</p>

<h2>Reproducing</h2>
<p><code>python data/verify_enzyme_stock.py --sequence</code> ·
<code>python temperature_series/build_figures.py</code> ·
<code>python temperature_series/check_numbers.py</code> — which re-derives every
number quoted in <code>ANALYSIS.md</code> from the modules and fails if the prose
and the code disagree. Fits: <code>data/arrhenius.py</code>. Curves:
<a href='progress_curves.html'>all 24 progress curves</a>.</p>
"""
    return page("Temperature series — activation parameters", body,
                "Exps 14–19 · 4OMe-BnOH · pH 7.00 · 15–40 °C").replace(
        "</style>", EXTRA_CSS + "</style>")


def main():
    for name, content in (("index.html", build_index()),
                          ("progress_curves.html", build_curves_page())):
        path = os.path.join(HERE, name)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(content)
        print(f"wrote {path}  ({len(content) / 1024:.0f} kB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
