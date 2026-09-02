"""
Builds induction/index.html.

Draws only; every number comes from `induction`, `scope`, `arrhenius` or
`curve_metrics`, so a figure and the prose in ANALYSIS.md cannot disagree about
a value without `check_numbers.py` saying so.

    python induction/build_figures.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
# svgplot lives at the repository root, shared with the sibling folders.
sys.path.insert(0, os.path.dirname(HERE))

import arrhenius
import induction
import scope
from curve_metrics import LAG_WINDOW, rolling_slope
from fit_dataset import source_floor
from svgplot import ACCENT, GRID, INK, MUTED, Axes, esc, page, PAGE_CSS

# ORDERED VARIABLES GET SEQUENTIAL RAMPS. Temperature and substrate rung are
# ordinal; the two channels and the two competing laws are not.
RUNGS = ["#6295c3", "#3d729f", "#1e5079", "#0c2f4d"]
TEMPERATURES = ["#9fc0dd", "#6295c3", "#3d729f", "#1e5079", "#123c60", "#0c2f4d"]
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


_CACHE = {}


def _table():
    """Every curve in the archive with its landmark, built once."""
    if "table" not in _CACHE:
        _CACHE["table"] = induction.induction_table(induction.WHOLE_ARCHIVE)
    return _CACHE["table"]


def _blocks():
    if "blocks" not in _CACHE:
        _CACHE["blocks"] = induction.induction_blocks(_table())
    return _CACHE["blocks"]


def _rate_track(curve, window=LAG_WINDOW):
    """The rolling rate this module's landmark is read off, for drawing."""
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    centres, slopes = rolling_slope(times, values, window,
                                    source_floor(curve.source))
    return centres - centres[0], slopes


# The matched pair drawn in figure A: exp 14 is the temperature series' 25 C
# run and exp 40 the enzyme-free run at the same substrate, buffer, peroxide
# and pH. `induction.channel_contrast` is what says they are matched.
SHOWN_CATALYSED = 14
SHOWN_FREE = 40


def figure_two_channels():
    # ylim to 1.35 so the legend band above 1.05 is one no curve can enter:
    # every track is scaled to its own maximum and so cannot exceed 1.
    axes = Axes(560, 290, (0, 1.0), (0, 1.35), pad=(66, 24, 46, 32))
    shown = []
    for experiment, colour, label in ((SHOWN_CATALYSED, CATEGORY[0],
                                       "with the chemzyme"),
                                      (SHOWN_FREE, CATEGORY[1],
                                       "no chemzyme")):
        for index, curve in enumerate(scope.curves_of(experiment)):
            centres, slopes = _rate_track(curve)
            top = float(np.max(slopes))
            if not top > 0:
                continue
            span = float(centres[-1]) or 1.0
            axes.line(centres / span, slopes / top, colour, width=1.4,
                      opacity=0.85 if index else 1.0)
        shown.append((colour, label))
    axes.hline(1.0, GRID, dash="3 3")
    axes.note(80, 44, f"catalysed: exp {SHOWN_CATALYSED} · enzyme-free: "
                      f"exp {SHOWN_FREE} · both 25 °C, pH ≈ 7, 82.5 mM H₂O₂",
              MUTED, size=10.5)
    for offset, (colour, label) in enumerate(shown):
        axes.note(80 + 190 * offset, 61, label, colour, size=11.5,
                  weight="600")
    return fig(
        axes.render("time, as a fraction of the run", "rate / its own maximum",
                    "A · The induction is there with the catalyst and absent "
                    "without it"),
        "The catalysed cuvettes climb for the first third of the run before they reach "
        "their own maximum rate; the enzyme-free cuvettes are at their maximum "
        "in the first window and never come back to it. Scaled per curve, so "
        "the comparison is of <em>shape</em> and not of size — the catalysed "
        "runs are four times faster in absolute terms. Across the whole "
        "archive this is 91 accelerating curves of 151 with the catalyst and "
        "<strong>0 of 49 without it</strong>, the largest acceleration "
        "z-score in the enzyme-free block being 2.25. Both channels' tracks are "
        "drawn from the rolling slope this folder's landmark is read off, and "
        "the dashed line is each curve's own maximum.")


def figure_channel_depth():
    table = _table()
    four = table[table.substrate == "4OMe-BnOH"]
    axes = Axes(560, 300, (12, 43), (-0.04, 1.32), pad=(66, 26, 46, 32))
    for offset, (differential, colour, label) in enumerate(
            ((True, CATEGORY[0], "with the chemzyme"),
             (False, CATEGORY[1], "no chemzyme"))):
        block = four[four.differential == differential]
        # Jitter WITHIN each temperature column: spreading one ramp over the
        # whole block gives the 15 C column four adjacent offsets and the 25 C
        # column the whole width, which reads as a trend that is not there.
        for temperature, column in block.groupby("temperature"):
            spread = np.linspace(-0.85, 0.85, len(column))
            axes.points(temperature + spread,
                        np.clip(column.depth.to_numpy(), 0.0, 1.0), colour,
                        radius=3.0, opacity=0.7, stroke=None)
        medians = block.groupby("temperature").depth.median()
        axes.line(medians.index.to_numpy(), np.clip(medians.to_numpy(), 0, 1),
                  colour, width=2.2)
        axes.points(medians.index.to_numpy(),
                    np.clip(medians.to_numpy(), 0, 1), colour, radius=5.0)
        axes.note(300 + 130 * offset, 46, label, colour, size=11.5,
                  weight="600")
    axes.note(80, 46, "one point per curve, medians joined", MUTED, size=10.5)
    return fig(
        axes.render("temperature, °C", "induction depth  (1 − start / peak)",
                    "B · How deep the induction is, and when it disappears"),
        "<code>depth</code> is the fraction of the peak rate missing from the "
        "first window. With the catalyst it falls smoothly from "
        "<strong>0.79 at 15 °C to 0.06 at 40 °C</strong> — the induction is "
        "still there when it is warm, it is just over before the run is "
        "properly under way. Without the catalyst it is <strong>zero at every "
        "temperature</strong>, including 37 curves at 40 °C and a five-hour "
        "run (exp 28). Depth is read through a window a tenth of the run wide, "
        "so it understates the true amplitude, and understates it most where "
        "the induction is fastest.")


def figure_clock_or_product():
    blocks = _blocks()
    rows = []
    for name, short in (("4OMe catalysed", "4OMe catalysed, 38 runs"),
                        ("4OMe catalysed, 25 C", "…the 25 °C runs alone"),
                        ("temperature series", "the temperature series")):
        rows.append((short, induction.induction_drivers(blocks[name],
                                                        rate="v_peak")))
    axes = Axes(680, 250, (-1.35, 0.75), (-0.7, len(rows) - 0.3),
                pad=(268, 26, 46, 34))
    for value, colour, label in ((induction.INDUCTION_PRODUCT_SLOPE,
                                  CATEGORY[1], "product threshold"),
                                 (induction.INDUCTION_CLOCK_SLOPE,
                                  CATEGORY[2], "a clock")):
        axes.line([value, value], [-0.7, len(rows) - 0.3], colour, width=1.6,
                  dash="4 3")
        axes.label(value, len(rows) - 1, label, colour, size=11,
                   weight="600", anchor="middle", dy=-14)
    for index, (name, row) in enumerate(rows):
        y = len(rows) - 1 - index
        axes.line([row["slope"] - row["stderr"], row["slope"] + row["stderr"]],
                  [y, y], CATEGORY[0], width=2.4)
        axes.points([row["slope"]], [y], CATEGORY[0], radius=4.6)
        axes.note(250, axes._fy(y) + 4, f"{name}   n = {row['points']}", INK,
                  size=11, anchor="end")
    return fig(
        axes.render("d log(induction time) / d log(rate)", "",
                    "C · The induction runs on a clock, not on product",
                    yticks=False),
        "Each curve's induction time regressed on its own peak rate, with one "
        "offset per experiment — so temperature, pH, buffer, catalyst, cell, "
        "day and <em>the run length</em> are absorbed, and only the "
        "four-fold substrate ladder inside each run carries the fit. A clock "
        "predicts 0; if the induction ended when enough product had "
        "accumulated, a curve twice as fast would get there in half the time "
        "and the coefficient would be −1. On the whole catalysed 4OMe archive "
        f"it is <strong>{rows[0][1]['slope']:+.3f} ± "
        f"{rows[0][1]['stderr']:.3f}</strong>, which excludes product control "
        f"by {(rows[0][1]['slope'] + 1) / rows[0][1]['stderr']:.0f} standard "
        "errors. The temperature series on its own cannot do it: 24 curves and "
        "a two-fold lever give ±0.44.")


def figure_orders():
    blocks = _blocks()
    rows = []
    for name, short in (("4OMe catalysed", "4OMe catalysed"),
                        ("temperature series", "the temperature series")):
        block = blocks[name]
        rows.append((short,
                     scope.orders("t_ind", frame=block,
                                  floor=induction.INDUCTION_FLOOR),
                     scope.orders("depth", frame=block,
                                  floor=induction.DEPTH_FLOOR),
                     scope.orders("v_peak", frame=block)))
    axes = Axes(680, 272, (-0.75, 0.75), (-0.7, len(rows) - 0.3),
                pad=(214, 26, 60, 40))
    axes.line([0, 0], [-0.7, len(rows) - 0.3], INK, width=1.2, dash="3 3")
    legend = (("the rate", CATEGORY[0], 0.26), ("the induction time",
                                                CATEGORY[1], 0.0),
              ("the induction depth", CATEGORY[2], -0.26))
    for index, (name, timing, depth, rate) in enumerate(rows):
        y = len(rows) - 1 - index
        for (label, colour, offset), row in zip(legend, (rate, timing, depth)):
            value, error = row["order_s0"], row["stderr_s0"]
            axes.line([value - error, value + error], [y + offset] * 2, colour,
                      width=2.2)
            axes.points([value], [y + offset], colour, radius=4.2)
        axes.note(196, axes._fy(y) + 4, name, INK, size=11, anchor="end")
    for position, (label, colour, offset) in enumerate(legend):
        axes.note(230 + 150 * position, 248, label, colour, size=11,
                  weight="600")
    return fig(
        axes.render("order in [S], within experiments", "",
                    "D · The substrate moves the rate and nothing else",
                    yticks=False),
        "The same log-log fit <code>scope.orders</code> uses everywhere, on "
        "three quantities. <strong>The rate carries a substrate order of "
        f"{rows[0][3]['order_s0']:+.3f} ± {rows[0][3]['stderr_s0']:.3f}</strong>"
        " — a real lever, half an order over a four-fold ladder — while the "
        f"induction time carries {rows[0][1]['order_s0']:+.3f} ± "
        f"{rows[0][1]['stderr_s0']:.3f} and its depth "
        f"{rows[0][2]['order_s0']:+.3f} ± {rows[0][2]['stderr_s0']:.3f}. "
        "This is figure C's answer again with the regressor replaced by the "
        "composition the operator set, which removes the one bias figure C "
        "cannot: noise in a measured rate attenuates its coefficient towards "
        "zero, and zero is the answer being defended.")


def figure_arrhenius():
    frame = arrhenius.series_frame()
    table = _table()
    resolved = frame[frame.tau_resolved & (frame.tau > 0)
                     & (frame.kelvin <= arrhenius.BURST_TRUSTWORTHY_BELOW_C
                        + 273.15)]
    fitted = arrhenius.activation_parameters("inverse_tau")
    x = 1000.0 / resolved.kelvin.to_numpy(dtype=float)
    y = 1.0 / resolved.tau.to_numpy(dtype=float)
    axes = Axes(560, 285, (float(x.min()) - 0.02, float(x.max()) + 0.02),
                (0.6 * float(y.min()), 1.7 * float(y.max())), ylog=True,
                pad=(72, 26, 46, 32))
    rungs = sorted(resolved.s0.unique())
    for index, rung in enumerate(rungs):
        keep = np.isclose(resolved.s0.to_numpy(), rung)
        axes.points(x[keep], y[keep], RUNGS[index], radius=4.2,
                    title=f"[S] = {rung:.3f} mM")
    grid = np.array([float(x.min()) - 0.01, float(x.max()) + 0.01])
    centre = float(np.exp(np.mean(np.log(y))))
    slope = -fitted["activation_kJ"] * 1000.0 / arrhenius.GAS_CONSTANT / 1000.0
    axes.line(grid, centre * np.exp(slope * (grid - float(np.mean(x)))),
              ACCENT, width=2.0, dash="6 4")
    axes.note(100, 186, f"Ea = {fitted['activation_kJ']:.0f} ± "
                        f"{fitted['activation_stderr']:.0f} kJ/mol", ACCENT,
              size=11.5, weight="600")
    axes.note(100, 202, "dissolution and diffusion run at 15–25 kJ/mol",
              MUTED, size=10.5)
    matched = resolved[resolved.experiment.isin(induction.SCHEDULE_PAIR)]
    for row in matched.itertuples():
        axes.ring(1000.0 / row.kelvin, 1.0 / row.tau, ACCENT, radius=8.0,
                  width=1.5)
    axes.note(100, 218, "ringed: every curve of exps 19 and 14, both "
                        "recorded for 17934 s", ACCENT, size=10.5)
    return fig(
        axes.render("1000 / T, K⁻¹", "1 / τ, s⁻¹",
                    "E · The induction has a barrier, and it is a chemical one"),
        "The induction's own rate constant, from the burst form's τ, on the "
        "four temperatures where the profile pins it. <strong>95 ± 16 "
        "kJ/mol</strong> is four times too large for dissolution, diffusion or "
        "a cuvette reaching temperature: whatever the catalyst is doing, it is "
        "making or breaking a bond. The two ringed temperatures are the "
        "schedule control — exps 19 and 14 were recorded for the same 17934 s "
        "and their τ differ <strong>6.9-fold</strong>, so the trend is not the "
        "operator's stopwatch.")


def figure_activation_gap():
    gap = induction.activation_contrast()
    rows = (("the induction, 1/τ", gap["induction"], CATEGORY[2]),
            ("turnover, v_peak", gap["turnover"], CATEGORY[0]))
    axes = Axes(560, 250, (-0.7, 1.7), (-70, 135), pad=(70, 24, 48, 34))
    axes.hline(0.0, GRID, dash=None, width=1.0)
    width = 0.20
    for index, (name, row, colour) in enumerate(rows):
        pieces = ((row["enthalpy_kJ"], "ΔH‡", 0.75),
                  (-arrhenius.REFERENCE_KELVIN * row["entropy_J"] / 1000.0,
                   "−TΔS‡", 0.42),
                  (row["gibbs_kJ"], "ΔG‡", 1.0))
        for position, (value, label, opacity) in enumerate(pieces):
            x = index + (position - 1) * 0.30
            axes.line([x, x], [0.0, value], colour, width=15, opacity=opacity)
            axes.label(x, value, f"{value:+.0f}", INK, size=10.5,
                       anchor="middle", dy=-6 if value > 0 else 14)
            if not index:
                axes.label(x, -58, label, MUTED, size=11, anchor="middle")
        # In pixels, below the frame: at these limits a data-coordinate label
        # for the column name lands on the axis line itself.
        axes.note(axes._fx(index), 226, name, colour, size=11.5,
                  anchor="middle", weight="600")
    axes.note(300, 182, f"ΔG‡ gap {gap['gibbs_gap_kJ']:+.2f} ± "
                        f"{gap['gibbs_gap_stderr']:.2f} kJ/mol  →  "
                        f"{gap['rate_ratio']:.0f}× faster at 298 K", ACCENT,
              size=11.5, weight="600", anchor="middle")
    return fig(
        axes.render("", "kJ/mol at 298 K",
                    "F · The induction is the faster step, as it has to be",
                    xticks=False),
        "Both columns are <code>arrhenius.activation_parameters</code>. "
        "<strong>The free-energy gap is the solid number</strong>: "
        f"{gap['gibbs_gap_kJ']:+.2f} ± {gap['gibbs_gap_stderr']:.2f} kJ/mol, "
        f"which makes the induction step {gap['rate_ratio']:.0f}× faster as a "
        "rate constant — required of anything that has to finish before "
        "turnover can begin, and predicted before it was looked at. "
        "<strong>Its decomposition is not solid.</strong> The point estimates "
        "put the whole gap in the entropy, but that gap carries "
        f"±{gap['entropy_gap_stderr']:.0f} J/mol/K and the enthalpy gap "
        f"±{gap['enthalpy_gap_stderr']:.0f} kJ/mol, so each is about one "
        "standard error from zero. The claim that the induction is "
        "unimolecular rests on figure D, not on this.")


def figure_signal_control():
    blocks = _blocks()
    table = _table()
    rows = [("4OMe catalysed", induction.signal_control(blocks["4OMe catalysed"])),
            ("BnOH, exps 135–151",
             induction.signal_control(blocks["BnOH in scope (135-151)"])),
            ("4OMe peroxide, exps 127–131",
             induction.signal_control(
                 table[table.experiment.isin(induction.PEROXIDE_LEVER)]))]
    axes = Axes(680, 235, (-0.45, 1.15), (-0.7, len(rows) - 0.3),
                pad=(272, 26, 46, 34))
    axes.line([0, 0], [-0.7, len(rows) - 0.3], INK, width=1.2, dash="3 3")
    for index, (name, row) in enumerate(rows):
        y = len(rows) - 1 - index
        clean = abs(row["signal_slope"]) < 2 * row["signal_stderr"]
        colour = CATEGORY[0] if clean else CATEGORY[1]
        axes.line([row["signal_slope"] - row["signal_stderr"],
                   row["signal_slope"] + row["signal_stderr"]], [y, y],
                  colour, width=2.4)
        axes.points([row["signal_slope"]], [y], colour, radius=4.6)
        axes.note(254, axes._fy(y) + 4,
                  f"{name}   n = {row['signal_points']}", INK, size=11,
                  anchor="end")
    axes.note(465, 186, "blue: passes · orange: the landmark is reading the "
                        "spectrophotometer", MUTED, size=10.5, anchor="middle")
    return fig(
        axes.render("d log(induction time) / d log(signal-to-noise)", "",
                    "G · The control that decides which blocks may be used",
                    yticks=False),
        "The landmark is the first crossing of half the <em>largest</em> "
        "rolling slope, so on a curve with no signal the largest rolling slope "
        "is a noise excursion and the landmark measures the instrument. "
        "Regressed on each curve's own <code>net/noise</code>, same offsets. "
        "<strong>The block this analysis rests on passes</strong> at "
        f"{rows[0][1]['signal_slope']:+.3f} ± "
        f"{rows[0][1]['signal_stderr']:.3f}. The two blocks that carry a "
        "peroxide ladder — the only design that could test whether the "
        "induction is the catalyst binding H₂O₂ — <strong>fail it</strong>, "
        "because [H₂O₂] is also what sets the signal. That question stays open "
        "and goes to <code>COMPUTATIONAL.md</code>.")


def figure_peroxide():
    table = _table()
    ladder = induction.peroxide_ladder(_blocks()["BnOH in scope (135-151)"])
    fitted = induction.peroxide_saturation(_blocks()["BnOH in scope (135-151)"])
    # Each run sits at its own level, which is what the per-experiment offsets
    # in the fit absorb; dividing every run by its own geometric mean puts them
    # on one panel without changing a single slope.
    scaled = ladder.copy()
    centre = scaled.groupby("experiment").vmax.transform(
        lambda column: float(np.exp(np.log(column).mean())))
    scaled["relative"] = scaled.vmax / centre
    axes = Axes(560, 300, (2.0, 200.0), (0.06, 9.0), xlog=True, ylog=True,
                pad=(66, 26, 46, 32))
    axes.points(scaled.h2o2.to_numpy(), scaled.relative.to_numpy(),
                CATEGORY[0], radius=3.4, opacity=0.8)
    grid = np.array([2.2, 190.0])
    middle = float(np.exp(np.log(scaled.h2o2.to_numpy()).mean()))
    # The two laws cross in the middle of the cloud, so their labels go in the
    # bottom-right corner, which is the one region a rising ladder leaves empty.
    for offset, (exponent, colour, label, dash) in enumerate((
            (fitted["order"], CATEGORY[0],
             f"a = {fitted['order']:.3f}, what the ladder says", "6 4"),
            (1.0, CATEGORY[1], "a = 1, a single unsaturated adduct", "2 4"))):
        axes.line(grid, (grid / middle) ** exponent, colour, width=2.0,
                  dash=dash)
        axes.note(520, 206 + 18 * offset, label, colour, size=11,
                  anchor="end", weight="600")
    axes.note(84, 44, f"{fitted['points']} cuvettes in "
                      f"{fitted['experiments']} runs, each divided by its own "
                      f"level", MUTED, size=10.5)
    joint = {name: induction.joint_peroxide_order(block) for name, block in
             (("4OMe", table[table.experiment.isin(induction.PEROXIDE_LEVER)]),
              ("BnOH", _blocks()["BnOH in scope (135-151)"]))}
    return fig(
        axes.render("[H₂O₂], mM", "rate / the run's own level",
                    "H · The rate is not first order in peroxide"),
        "The peroxide arm of each in-scope run, over a 67-fold range. "
        f"<strong>Strict first order is rejected at F = "
        f"{fitted['first_order_f']:.0f}</strong>; the free power law is "
        f"a = {fitted['order']:.3f}. It matters because first order is the "
        "<em>unsaturated</em> limit of <code>K + H₂O₂ ⇌ KP</code>, not a "
        "general consequence of it — and the order stays fractional all the "
        "way across rather than falling from 1 towards 0 the way one binding "
        "equilibrium saturating would make it (forcing that form on these "
        f"points fits worse, {fitted['scheme_sse']:.2f} against "
        f"{fitted['power_sse']:.2f}). The same scheme fixes the difference of "
        "the two orders at exactly 1 whatever K is; measured as one regression "
        f"on <code>log(v/t_ind)</code> it is "
        f"{joint['4OMe']['slope']:+.3f} ± {joint['4OMe']['stderr']:.3f} on the "
        f"4OMe block ({joint['4OMe']['sigma']:.1f}σ short) and "
        f"{joint['BnOH']['slope']:+.3f} ± {joint['BnOH']['stderr']:.3f} here "
        f"({joint['BnOH']['sigma']:.1f}σ short). Neither is clean — figure G — "
        "but both fall the same way.")


def build_index():
    table = _table()
    summary = induction.channel_summary(table)
    drivers = induction.induction_drivers(_blocks()["4OMe catalysed"],
                                          rate="v_peak")
    gap = induction.activation_contrast()
    contrast = induction.channel_contrast(table)

    hero = f"""
<div class='hero'>
  <div><div class='k'>with the catalyst</div>
       <div class='v'>{summary['catalysed']['accelerates']} / {summary['catalysed']['curves']}</div>
       <div class='u'>4OMe curves that accelerate</div></div>
  <div><div class='k'>without it</div>
       <div class='v'>{summary['enzyme_free']['accelerates']} / {summary['enzyme_free']['curves']}</div>
       <div class='u'>largest z-score 2.25</div></div>
  <div><div class='k'>clock or product</div>
       <div class='v'>{drivers['slope']:+.3f}</div>
       <div class='u'>± {drivers['stderr']:.3f}, where product control needs −1</div></div>
  <div><div class='k'>faster than turnover</div>
       <div class='v'>{gap['rate_ratio']:.0f}×</div>
       <div class='u'>ΔG‡ gap {gap['gibbs_gap_kJ']:+.2f} ± {gap['gibbs_gap_stderr']:.2f} kJ/mol</div></div>
</div>"""

    matched = "".join(
        f"<tr><td>{row.temperature:.0f} °C</td><td>{row.channel}</td>"
        f"<td>{row.curves}</td><td>{row.depth:.3f}</td>"
        f"<td>{row.accelerates}</td><td>{row.longest_s:.0f} s</td></tr>"
        for row in contrast.itertuples())

    windows = induction.landmark_window()
    window_rows = "".join(
        f"<tr><td>{esc(row.window)}</td><td>{row.cold_s:.0f} s</td>"
        f"<td>{row.warm_s:.0f} s</td>"
        f"<td>{row.activation_kJ:.1f} ± {row.stderr_kJ:.1f}</td></tr>"
        for row in windows.itertuples())

    body = f"""
<p class='lede'>Every catalysed 4-methoxybenzyl alcohol run from 15 to 30 °C
begins slowly and takes thousands of seconds to reach the rate the temperature
series puts on an Arrhenius plot. <a href='ANALYSIS.md'>ANALYSIS.md</a> is the
argument; this is the picture of it.</p>
{hero}

<h2>1 · Four candidates, and what each one waits for</h2>
<p>Seeding of the catalyst-free Cannizzaro loop and a scavenger in the reagents
both end when enough <strong>product</strong> or <strong>turnover</strong> has
happened, and both run in a cuvette with no catalyst in it. Catalyst activation
ends after a fixed <strong>time</strong> and cannot happen without the catalyst.
One curve cannot separate them — inside a curve the product only grows with
time — so the separation comes from the substrate ladder, which moves the rate
four-fold while the schedule, the peroxide, the catalyst and the temperature
stay where they are.</p>
{figure_two_channels()}
{figure_channel_depth()}

<h2>2 · The induction needs the catalyst</h2>
<div class='tbl'><table>
<tr><th>T</th><th>channel</th><th>curves</th><th>median depth</th>
<th>accelerating</th><th>longest run</th></tr>
{matched}
</table></div>
<p>Matched on substrate, buffer, 82.5 mM H₂O₂, pH and temperature.
<strong>Not one of the 49 enzyme-free 4OMe curves has an induction</strong>, and
the objection that they were too short does not hold: exp 28 ran 17934 s at
40 °C, the longest run in the archive, and is flat from its first window.</p>

<h2>3 · It runs on a clock, not on product</h2>
{figure_clock_or_product()}
{figure_orders()}
<p><strong>This is the exact mirror of the fall.</strong> In the same block
<a href='../product_fate/index.html'>product_fate</a> finds the deceleration
carried by the product (−0.598 ± 0.053) and not by the clock (−0.124 ± 0.079).
The rise and the fall of one curve are driven by different things, and the
intuition that would assign them the other way round — a rise that waits for
product, a fall that runs down a clock — is wrong at both ends.</p>

<h2>4 · Not the schedule, not the instrument, not physical</h2>
{figure_arrhenius()}
<div class='tbl'><table>
<tr><th>window</th><th>induction at 15 °C</th><th>at 40 °C</th>
<th>E<sub>a</sub> of 1/t<sub>ind</sub>, kJ/mol</th></tr>
{window_rows}
</table></div>
<p>The landmark's window is a fraction of the run, and it has to be: a 15 °C run
is 0.04 AU over five hours, so a 300 s window reads its slope through noise, the
maximum lands on an early excursion and the induction collapses from 4289 s to
529 s — taking the activation energy with it. A window that is a fraction of the
run is safe <em>within</em> a run, where every cuvette shares the schedule, and
unsafe between runs, where it is not. <strong>Every concentration order on this
page is measured within experiments.</strong></p>
{figure_signal_control()}

<h2>4a · What an adduct with H₂O₂ would require, and does not get</h2>
<p>With <code>h = [H₂O₂]</code> in 100–6000× excess, <code>K + H₂O₂ ⇌ KP</code>
gives <code>1/τ = k_f·h + k_r</code> and <code>[KP]/E₀ = Kh/(1+Kh)</code>.
<strong><code>1/τ</code> increases with h whatever the constants are</strong> —
no K and no concentration make the approach slower — so a <em>positive</em>
order on the induction time is not a statement about the regime, it falsifies
the scheme. The only 4OMe block that moves [H₂O₂] gives
<strong>+0.302 ± 0.092</strong>.</p>
<p>The two orders are also locked: their difference is
<strong>1 identically</strong>, for every K and every h, so the constraint can
be tested without knowing where on the saturation curve the design sits.</p>
{figure_peroxide()}
<p>If the sign survives, the perhydrate is a <strong>trap</strong> rather than
the activation — <code>1/τ = k_act/(1+Kh)</code>, which is positive, bounded by
1, and keeps everything section 3 establishes. Inverting the orders gives K =
11 /M and 29 /M on the two blocks, and the saturating fit to the rates gives
29 /M (13–54): <strong>ΔG° between −6 and −10 kJ/mol</strong>. Three routes
within a factor of three — which is <em>not</em> evidence, since two of them
share a confound, but is a target. A computed ΔG° in that range would
corroborate the trap from a direction with no spectrophotometer in it.</p>
<h2>5 · What the activation parameters say</h2>
{figure_activation_gap()}

<h2>6 · What this settles, and what it does not</h2>
<p><strong>Settled.</strong> The induction is a property of the catalysed
reaction; it ends on a clock and not at a product threshold; its amplitude is a
fraction of the eventual rate and does not scale with the substrate; its barrier
is chemical rather than physical; and it is faster than turnover by the factor
it must be.</p>
<p><strong>Not settled: which step it is.</strong> "The catalyst becomes active"
is a shape, not a mechanism. Everything that survives section 3 is unimolecular
in what the cuvette holds — the ketone's gem-diol hydrate dehydrating, the
perhydrate collapsing to the dioxirane, or a conformational change of the
cyclodextrin — and absorbance at one wavelength cannot choose between them.</p>
<p><strong>Not settled: whether the peroxide is involved at all.</strong>
Section 4a: three things point the same way and none is clean — the induction's
peroxide order has the wrong sign for an adduct, the joint constraint the scheme
puts on both orders at once falls short by 2.6σ and 3.7σ, and the rate is not
first order in peroxide either. Both induction orders come from blocks that fail
figure G, so their agreement is also what one shared artefact looks like.</p>
<p><strong>What would finish it.</strong> <code>COMPUTATIONAL.md</code>
<strong>C7</strong> — the hydration equilibrium and dehydration barrier of the
chemzyme's ketone against the barrier for adding H₂O₂ to it — and
<strong>C8</strong>, the profile from the perhydrate to the dioxirane. And the
experiment that will not be run: pre-incubate the catalyst with H₂O₂, then add
substrate. If the induction is catalyst activation it disappears. One
cuvette.</p>

<h2>Reproducing</h2>
<p><code>python data/induction.py</code> prints the whole argument ·
<code>python data/test_induction.py</code> ·
<code>python induction/build_figures.py</code> ·
<code>python induction/check_numbers.py</code>, which re-derives every number in
<code>ANALYSIS.md</code> from the modules and fails if the prose and the code
disagree.</p>
"""
    return page("What happens before the 4OMe curves start", body,
                "The catalyst wakes up on its own clock").replace(
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
