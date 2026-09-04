"""
One self-contained HTML page per experiment, for judging curves by eye.

`build_dossier.py` answers "is this experiment's METADATA right?" and puts all
hundred experiments on one 6 MB page. This answers a different question -- "is
this CUVETTE's curve any good?" -- and puts each experiment in its own small
file, because that is the unit a curation judgement is actually made in.

The screens in `summary_kinetics.py` can flag a curve but cannot tell you why
it is odd. A curve reading 25x below its neighbours might be a failed cuvette
or the slow bottom rung of a titration, and the only reliable way to tell them
apart is to look. So every page shows the automatic verdict AND the evidence
behind it, and never hides a curve because a rule disliked it.

Deliberately shows MORE than the fitting dataset does: curves excluded by
`build_manifest` appear too, marked with the rule that removed them, so the
exclusions themselves can be reviewed rather than taken on trust.

Plots are hand-generated inline SVG -- no matplotlib, no JavaScript, no
external files. A page opens instantly, zooms without pixelating, and survives
being emailed or committed.

    python data/curve_dossier.py                    # every experiment
    python data/curve_dossier.py --experiments 25,34
    python data/curve_dossier.py --out dossier
"""
import argparse
import html
import os

import numpy as np
import pandas as pd

from build_manifest import (EXCLUDED_BUFFERS, KNOWN_EXCLUSIONS,
                            KNOWN_SAMPLE_EXCLUSIONS)
from curve_metrics import (ACCELERATION_SIGMA, QUANTISATION_SIGMA,
                           acceleration, peak_rate)
from fit_dataset import (BASELINE_POINTS, curve_noise, read_all_curves,
                         source_floor)
from solution_chemistry import add_solution_columns
from curve_screen import MINIMUM_WINDOW_QUANTA, eligibility
from summary_kinetics import (BURST_TAU_CAP, DEAD_CURVE_SNR, INITIAL_WINDOW,
                              OUTLIER_FACTOR, fit_burst, initial_rate, line_fit,
                              slope_ratio, window_quanta, window_size)

DATASET_PATH = "data/experiment_data.csv"
CURVE_DIRECTORY = "data/data"
OUTPUT_DIRECTORY = "dossier"

# NOT svgplot.PALETTE, which is EIGHT colours and diverges from this one from
# its sixth on (#4a9ab0/#a03a5a/#6a7a3a against #4a8a9a/#a83f5a/#6a6a6a). Both
# were called PALETTE, in three modules, until 2026-09-04. This is the nine
# the two dossier pages cycle over their cuvettes; `build_dossier` imports it
# from here rather than keeping the third copy.
DOSSIER_PALETTE = ["#2f6fb0", "#c0522a", "#3f8a5a", "#8a5aa8", "#b08a2f",
                   "#4a8a9a", "#a83f5a", "#6a6a6a", "#5a7a3f"]

LINE_COLOUR = "#d0342c"    # the straight line v0 comes from
# NOT figure_kit.BURST_COLOUR, which is the purple the analysis folders
# draw this form in. This page is a different set of marks on a different
# ground and picked teal; the two shared the bare name until 2026-09-04.
DOSSIER_BURST_COLOUR = "#12856a"  # the burst/lag form, never used for a rate

# A downward excursion larger than this many noise sigma is a curve going
# backwards, not a curve wobbling. Seven curves in the dataset move 0.03-0.10 AU
# yet start with a negative slope; an initial-rate window measures the transient
# on those, not the reaction, so they must be visible at a glance.
BACKTRACK_SIGMA = 5.0


# --- tiny SVG plotting ----------------------------------------------------

def _nice_ticks(low, high, count=5):
    """Round tick values spanning [low, high], at a 1/2/5 x 10^n step."""
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        return [low] if np.isfinite(low) else [0.0]
    raw = (high - low) / max(count, 1)
    magnitude = 10 ** np.floor(np.log10(raw))
    step = min((s for s in (1, 2, 5, 10) if s * magnitude >= raw),
               default=10) * magnitude
    first = np.ceil(low / step) * step
    return [first + i * step for i in range(int((high - first) / step) + 1)]


def _format_tick(value):
    if value == 0:
        return "0"
    if 1e-3 <= abs(value) < 1e5:
        text = f"{value:,.6g}"
    else:
        text = f"{value:.0e}".replace("e-0", "e-").replace("e+0", "e")
    return text


def svg_chart(series, width=780, height=330, xlabel="", ylabel="",
              logx=False, logy=False, legend=True, margin=None):
    """
    A line/scatter chart as standalone inline SVG.

    `series` is a list of dicts with keys x, y, label, colour, and optional
    dash (bool) and marker (bool). Written by hand rather than through
    matplotlib so a page stays a few tens of kB and needs nothing to render.
    """
    margin = margin or dict(left=68, right=14, top=12, bottom=44)
    if legend:
        margin = dict(margin, right=margin["right"] + 96)
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]

    def usable(entry):
        x = np.asarray(entry["x"], dtype=float)
        y = np.asarray(entry["y"], dtype=float)
        keep = np.isfinite(x) & np.isfinite(y)
        if logx:
            keep &= x > 0
        if logy:
            keep &= y > 0
        return x[keep], y[keep]

    points = [usable(entry) for entry in series]
    xs = np.concatenate([p[0] for p in points if len(p[0])]) if any(
        len(p[0]) for p in points) else np.array([0.0, 1.0])
    ys = np.concatenate([p[1] for p in points if len(p[1])]) if any(
        len(p[1]) for p in points) else np.array([0.0, 1.0])

    def limits(values, log):
        low, high = float(np.min(values)), float(np.max(values))
        if log:
            low, high = np.log10(low), np.log10(high)
        if high == low:
            high, low = high + 0.5, low - 0.5
        pad = (high - low) * 0.06
        return low - pad, high + pad

    x_low, x_high = limits(xs, logx)
    y_low, y_high = limits(ys, logy)

    def to_px(value, low, high, size, log, flip=False):
        value = np.log10(value) if log else value
        fraction = (value - low) / (high - low)
        return (1 - fraction) * size if flip else fraction * size

    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" '
             f'preserveAspectRatio="xMidYMid meet" class="chart">',
             f'<g transform="translate({margin["left"]},{margin["top"]})">']

    # grid and ticks
    for value in _nice_ticks(y_low, y_high):
        py = to_px(10 ** value if logy else value, y_low, y_high, plot_h, logy, True)
        if not 0 <= py <= plot_h:
            continue
        parts.append(f'<line class="grid" x1="0" y1="{py:.1f}" x2="{plot_w}" y2="{py:.1f}"/>')
        label = _format_tick(10 ** value if logy else value)
        parts.append(f'<text class="tick" x="-8" y="{py + 4:.1f}" text-anchor="end">{label}</text>')
    for value in _nice_ticks(x_low, x_high):
        px = to_px(10 ** value if logx else value, x_low, x_high, plot_w, logx)
        if not 0 <= px <= plot_w:
            continue
        parts.append(f'<line class="grid" x1="{px:.1f}" y1="0" x2="{px:.1f}" y2="{plot_h}"/>')
        label = _format_tick(10 ** value if logx else value)
        parts.append(f'<text class="tick" x="{px:.1f}" y="{plot_h + 18}" '
                     f'text-anchor="middle">{label}</text>')
    parts.append(f'<rect class="frame" x="0" y="0" width="{plot_w}" height="{plot_h}"/>')

    for entry, (x, y) in zip(series, points):
        if not len(x):
            continue
        colour = entry.get("colour", "#333")
        coords = " ".join(
            f"{to_px(a, x_low, x_high, plot_w, logx):.1f},"
            f"{to_px(b, y_low, y_high, plot_h, logy, True):.1f}"
            for a, b in zip(x, y))
        dash = ' stroke-dasharray="5 3"' if entry.get("dash") else ""
        if entry.get("marker"):
            for a, b in zip(x, y):
                parts.append(
                    f'<circle cx="{to_px(a, x_low, x_high, plot_w, logx):.1f}" '
                    f'cy="{to_px(b, y_low, y_high, plot_h, logy, True):.1f}" '
                    f'r="3.2" fill="{colour}"/>')
        if len(x) > 1 or not entry.get("marker"):
            parts.append(f'<polyline class="trace" points="{coords}" '
                         f'stroke="{colour}"{dash}/>')

    parts.append(f'<text class="axis" x="{plot_w / 2:.0f}" y="{plot_h + 38}" '
                 f'text-anchor="middle">{html.escape(xlabel)}</text>')
    parts.append(f'<text class="axis" transform="translate(-52,{plot_h / 2:.0f}) rotate(-90)" '
                 f'text-anchor="middle">{html.escape(ylabel)}</text>')

    if legend:
        for index, entry in enumerate(series):
            if not entry.get("label"):
                continue
            y = 6 + index * 17
            parts.append(f'<line x1="{plot_w + 10}" y1="{y}" x2="{plot_w + 30}" y2="{y}" '
                         f'stroke="{entry.get("colour", "#333")}" stroke-width="2.4"/>')
            parts.append(f'<text class="legend" x="{plot_w + 35}" y="{y + 4}">'
                         f'{html.escape(entry["label"])}</text>')
    parts.append("</g></svg>")
    return "".join(parts)


# --- per-cuvette diagnostics ---------------------------------------------

def backtrack(values, noise):
    """
    Largest downward excursion from the running maximum, in noise sigma.

    A progress curve should not go back down. When it does by more than a few
    sigma the signal is not a single rising reaction, and an initial-rate
    window fitted to it is measuring something else.
    """
    values = np.asarray(values, dtype=float)
    if len(values) < 3 or noise <= 0:
        return 0.0
    drop = np.maximum.accumulate(values) - values
    return float(drop.max() / noise)


def describe(times, values, noise, epsilon, s0, floor=QUANTISATION_SIGMA):
    """Everything the page reports about one cuvette."""
    count = window_size(len(times), INITIAL_WINDOW)
    # line_fit, not initial_rate: the page draws the fitted line, so it needs
    # the intercept the fit actually chose. Reconstructing it from the baseline
    # put the drawn line up to 124 sigma away from the real one on some curves.
    intercept, v0, stderr, rms = line_fit(times[:count], values[:count], floor)
    amplitude = float(values.max() - values.min())
    # v0 is the rate before any catalyst has built up. On an accelerating
    # curve that is the induction period, not the reaction, so the page
    # also carries the steepest block's rate and says which is which.
    accel_z, accel_where = acceleration(times, values, INITIAL_WINDOW, floor)
    vmax, vmax_stderr, vmax_where = peak_rate(times, values, INITIAL_WINDOW,
                                              floor)
    return dict(
        accel_z=accel_z, accel_where=accel_where,
        vmax=vmax, vmax_stderr=vmax_stderr, vmax_where=vmax_where,
        accelerates=bool(np.isfinite(accel_z) and accel_z > ACCELERATION_SIGMA),
        intercept=intercept,
        burst=fit_burst(times, values),
        window_quanta=window_quanta(times, values, INITIAL_WINDOW),
        v0=v0, v0_stderr=stderr, window_rms=rms,
        window_points=count,
        amplitude=amplitude,
        snr=amplitude / noise if noise > 0 else np.inf,
        noise=noise,
        ratio=slope_ratio(times, values, INITIAL_WINDOW),
        backtrack=backtrack(values, noise),
        conversion=amplitude / (epsilon * s0) if epsilon > 0 and s0 > 0 else np.nan,
        points=len(times), duration=float(times[-1]),
    )


def verdicts(row, block):
    """
    Short flags for one cuvette, each with the evidence that raised it.

    Every flag is a reason to LOOK, never a decision. The wording says what was
    measured, so a flag can be dismissed from the page without rerunning
    anything.
    """
    notes = []
    if row["experiment"] in KNOWN_EXCLUSIONS:
        notes.append(("excluded", f"experiment excluded: {KNOWN_EXCLUSIONS[row['experiment']]}"))
    pair = (row["experiment"], row["sample"])
    if pair in KNOWN_SAMPLE_EXCLUSIONS:
        notes.append(("excluded", f"cuvette excluded: {KNOWN_SAMPLE_EXCLUSIONS[pair]}"))
    if row["buffer"] in EXCLUDED_BUFFERS:
        notes.append(("excluded", f"{row['buffer']} buffer excluded wholesale"))

    if not np.isfinite(row["v0"]) or row["v0"] <= 0:
        notes.append(("bad", f"initial slope is not positive ({row['v0']:+.2e} AU/s)"))
    if row["backtrack"] > BACKTRACK_SIGMA:
        notes.append(("bad", f"goes backwards by {row['backtrack']:.0f} sigma "
                             f"({row['backtrack'] * row['noise']:.3f} AU)"))
    if row["snr"] < DEAD_CURVE_SNR:
        notes.append(("weak", f"whole curve rose only {row['snr']:.0f}x its own noise"))

    rates = np.array([r["v0"] for r in block if np.isfinite(r["v0"]) and r["v0"] > 0])
    if len(rates) >= 3 and np.isfinite(row["v0"]) and row["v0"] > 0:
        centre = float(np.median(rates))
        if row["v0"] * OUTLIER_FACTOR < centre:
            notes.append(("weak", f"{centre / row['v0']:.0f}x below this experiment's "
                                  f"median rate ({centre:.2e} AU/s)"))
    if row["window_rms"] > 2.5 * row["noise"]:
        notes.append(("shape", f"the initial-rate line misfits its own window by "
                               f"{row['window_rms'] / row['noise']:.1f}x noise"))
    # The verdict is the sigma statistic, not the late/early ratio. A ratio
    # compares two windows without asking whether the difference clears
    # their own noise, and it reads a finished sigmoid -- lag, burst, then
    # plateau -- as a DEceleration, which is how exp 135 s3 and s4 came to
    # look flat while accelerating at 11 and 28 sigma.
    if row["accelerates"]:
        notes.append(("shape",
                      f"accelerates: steepest at {row['accel_where']:.0%} of the "
                      f"run, {row['vmax'] / row['v0']:.1f}x the initial slope "
                      f"({row['accel_z']:+.0f} sigma)"
                      if np.isfinite(row["v0"]) and row["v0"] > 0 else
                      f"accelerates: steepest at {row['accel_where']:.0%} of the "
                      f"run ({row['accel_z']:+.0f} sigma)"))
    elif np.isfinite(row["accel_z"]) and row["accel_z"] < -ACCELERATION_SIGMA:
        notes.append(("shape", f"decelerates: steepest block is the first "
                               f"({row['accel_z']:+.0f} sigma)"))
    return notes


# --- page ------------------------------------------------------------------

# NOT build_dossier.STYLE, which is a different sheet for a different page.
CURVE_DOSSIER_STYLE = """
:root{--ink:#1b1e23;--dim:#5d646e;--faint:#8a919b;--rule:#d8dce2;--bg:#ffffff;
      --panel:#f6f7f9;--bad:#b3261e;--weak:#8a5a00;--shape:#2f5fa0;--excl:#6a3d9a;
      --power:#3f7a6a;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
     font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;}
.wrap{max-width:1080px;margin:0 auto;padding:28px 22px 60px}
h1{font-size:22px;margin:0 0 2px;letter-spacing:-.01em}
h2{font-size:13px;text-transform:uppercase;letter-spacing:.09em;color:var(--dim);
   margin:34px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--rule)}
.sub{color:var(--dim);font-size:14px;margin:0 0 18px}
.sub b{color:var(--ink);font-weight:600}
.chart{display:block;background:var(--bg)}
.grid{stroke:var(--rule);stroke-width:1;shape-rendering:crispEdges}
.frame{fill:none;stroke:var(--faint);stroke-width:1;shape-rendering:crispEdges}
.trace{fill:none;stroke-width:1.8;stroke-linejoin:round;stroke-linecap:round}
.tick{font-size:11px;fill:var(--dim)}
.axis{font-size:12px;fill:var(--dim)}
.legend{font-size:11.5px;fill:var(--ink)}
table{border-collapse:collapse;width:100%;font-size:13px;
      font-variant-numeric:tabular-nums}
th{text-align:right;font-weight:600;color:var(--dim);font-size:11px;
   text-transform:uppercase;letter-spacing:.05em;padding:6px 9px;
   border-bottom:1px solid var(--rule);white-space:nowrap}
th:first-child,td:first-child{text-align:left}
td{text-align:right;padding:6px 9px;border-bottom:1px solid var(--rule)}
tr.flagged td{background:#fdf6f4}
.swatch{display:inline-block;width:11px;height:11px;border-radius:2px;
        margin-right:7px;vertical-align:-1px}
.scroll{overflow-x:auto}
.small{display:grid;grid-template-columns:repeat(auto-fill,minmax(238px,1fr));gap:14px}
.card{border:1px solid var(--rule);border-radius:6px;padding:10px 10px 4px;
      background:var(--panel)}
.card h3{margin:0 0 2px;font-size:13px;font-weight:600}
.card p{margin:0 0 6px;font-size:11.5px;color:var(--dim);
        font-variant-numeric:tabular-nums}
p.accel{margin:2px 0 6px;font-size:11px;color:var(--muted)}
.notes{list-style:none;padding:0;margin:0}
.notes li{padding:7px 11px;border-left:3px solid var(--rule);margin-bottom:6px;
          background:var(--panel);font-size:13px;border-radius:0 4px 4px 0}
.notes li b{font-weight:600}
.bad{border-left-color:var(--bad)} .bad b{color:var(--bad)}
.weak{border-left-color:var(--weak)} .weak b{color:var(--weak)}
.shape{border-left-color:var(--shape)} .shape b{color:var(--shape)}
.excluded{border-left-color:var(--excl)} .excluded b{color:var(--excl)}
.power{border-left-color:var(--power)} .power b{color:var(--power)}
.none{color:var(--dim);font-style:italic}
.foot{margin-top:44px;padding-top:14px;border-top:1px solid var(--rule);
      color:var(--faint);font-size:11.5px}
.foot a{color:var(--faint)}
@media (prefers-color-scheme:dark){
 :root{--ink:#e6e8ec;--dim:#9aa2ad;--faint:#6e7681;--rule:#2f343c;--bg:#14171b;
       --panel:#1b1f25;--bad:#ff6b5e;--weak:#e0a341;--shape:#79a9ee;--excl:#c08ae8;
       --power:#5fc0a5;}
 tr.flagged td{background:#241b1a}
}
"""


def _cell(value, spec=".3g", dash="--"):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return dash
    return format(value, spec)


def render(number, rows, curves_by_sample):
    """The HTML page for one experiment."""
    first = rows[0]
    title = (f"Experiment {number} &middot; {html.escape(str(first['substrate']))} "
             f"&middot; {html.escape(str(first['buffer']))} "
             f"&middot; {first['T']:.0f} &deg;C")
    catalysed = any(r["[enz]"] > 0 for r in rows)
    total_points = sum(r["points"] for r in rows)
    longest = max(r["duration"] for r in rows)

    def colour(index):
        return DOSSIER_PALETTE[index % len(DOSSIER_PALETTE)]

    # --- overlay of every cuvette
    overlay = []
    for index, row in enumerate(rows):
        times, values = curves_by_sample[row["sample"]]
        overlay.append(dict(x=times, y=values, colour=colour(index),
                            label=f"{row['sample']}: [S]={row['[sub]']:.3g} mM"))
    overlay_svg = svg_chart(overlay, xlabel="time (s)",
                            ylabel="absorbance (baseline-subtracted)")

    # --- one small panel per cuvette, with the window the rate came from
    cards = []
    for index, row in enumerate(rows):
        times, values = curves_by_sample[row["sample"]]
        count = row["window_points"]
        fit_x = times[:count]
        fit_y = row["intercept"] + row["v0"] * (fit_x - fit_x[0])
        layers = [dict(x=times, y=values, colour=colour(index), label="")]
        burst = row["burst"]
        if np.isfinite(burst.tau):
            layers.append(dict(x=times, y=burst.predict(times),
                               colour=DOSSIER_BURST_COLOUR, label=""))
        layers.append(dict(x=fit_x, y=fit_y, colour=LINE_COLOUR, dash=True, label=""))
        panel = svg_chart(
            layers, width=250, height=150, legend=False,
            margin=dict(left=44, right=8, top=8, bottom=28),
            xlabel="t (s)", ylabel="AU")
        tone = ("bad" if any(k == "bad" for k, _ in row["notes"])
                else "weak" if any(k == "weak" for k, _ in row["notes"]) else "")
        cards.append(
            f'<div class="card"><h3><span class="swatch" style="background:{colour(index)}">'
            f'</span>sample {row["sample"]}{"  &#9888;" if tone else ""}</h3>'
            f'<p>v<sub>0</sub> = {_cell(row["v0"], ".3e")} AU/s &middot; '
            f'v<sub>max</sub> = {_cell(row["vmax"], ".3e")} &middot; '
            f'SNR {_cell(row["snr"], ".0f")}</p>'
            f'<p class="accel">{"accelerates " + format(row["accel_z"], "+.0f") + " sigma, peak at " + format(row["accel_where"], ".0%") if row["accelerates"] else "no acceleration (" + _cell(row["accel_z"], "+.0f") + " sigma)"}</p>'
            f'{panel}</div>')

    # --- the ladder: rate against whatever this experiment varied
    varied = [(name, key) for name, key in
              (("[S] (mM)", "[sub]"), ("[H2O2] (mM)", "[h2o2]"),
               ("[HOO-] (mM)", "[HOO-]"), ("[buffer] (mM)", "[buf]"))
              if len({round(r[key], 12) for r in rows}) > 1]
    ladder_svg = ""
    if varied and len(rows) >= 3:
        name, key = varied[0]
        good = [r for r in rows if np.isfinite(r["v0"]) and r["v0"] > 0 and r[key] > 0]
        flagged = [r for r in good if r["defects"]]
        clean = [r for r in good if not r["defects"]]
        ladder = []
        if clean:
            ladder.append(dict(x=[r[key] for r in clean], y=[r["v0"] for r in clean],
                               colour="#2f6fb0", marker=True, label="unflagged"))
        if flagged:
            ladder.append(dict(x=[r[key] for r in flagged], y=[r["v0"] for r in flagged],
                               colour="#c0522a", marker=True, label="flagged"))
        ladder_svg = (f'<h2>Rate against {html.escape(name)}</h2>'
                      f'<p class="sub">Log-log. A broken rung shows up here as a point off '
                      f'the ladder, which is the quickest way to separate a failed cuvette '
                      f'from a genuinely slow one.</p>'
                      + svg_chart(ladder, width=560, height=300, logx=True, logy=True,
                                  xlabel=name, ylabel="initial rate (AU/s)"))

    # --- burst / lag
    burst_columns = [("sample", None, None), ("shape", None, None),
                     ("v0 = v_ss - B/tau", "v0", ".3e"), ("v_ss", "v_ss", ".3e"),
                     ("B (AU)", "B", "+.4f"), ("tau (s)", "tau", ".0f"),
                     ("tau 95%", None, None), ("lag = B/v_ss (s)", "lag_time", ".0f"),
                     ("rms/noise", None, None)]
    head = "".join(f"<th>{html.escape(label)}</th>" for label, _, _ in burst_columns)
    body, any_resolved = "", False
    for index, row in enumerate(rows):
        fit = row["burst"]
        any_resolved |= fit.resolved
        low, high = fit.tau_interval
        interval = ("--" if not np.isfinite(low)
                    else f"{low:.0f} &ndash; {high:.0f}")
        badge = ({"lag": "lag", "burst": "burst", "linear": "linear"}
                 .get(fit.kind, "<b>unresolved</b>"))
        body += (
            f'<tr class="{"" if fit.resolved else "flagged"}">'
            f'<td><span class="swatch" style="background:{colour(index)}"></span>'
            f'{row["sample"]}</td>'
            f'<td>{badge}</td>'
            f'<td>{_cell(fit.v0, ".3e")}</td><td>{_cell(fit.v_ss, ".3e")}</td>'
            f'<td>{_cell(fit.B, "+.4f")}</td><td>{_cell(fit.tau, ".0f")}</td>'
            f'<td>{interval}</td><td>{_cell(fit.lag_time, ".0f")}</td>'
            f'<td>{_cell(fit.rms / row["noise"] if row["noise"] > 0 else np.nan, ".2f")}</td>'
            f"</tr>")
    caveat = ("" if any_resolved else
              " <b>On this experiment &tau; is unresolved on every cuvette</b>, so none of "
              "v<sub>0</sub>, B or &tau; above is a measurement &mdash; they are one point "
              "on a flat valley.")
    burst_section = f"""<h2>Burst / lag fit</h2>
<p class="sub">A = c + v<sub>ss</sub>t &minus; B(1 &minus; e<sup>&minus;t/&tau;</sup>),
&tau; profiled on a log grid from {1 / 300:.4f} to {BURST_TAU_CAP:.0f} times the run
length. B &gt; 0 is a lag (rate rising), B &lt; 0 a burst (rate falling).
A row is highlighted when the 95% profile interval on &tau; reaches an end of that
grid, which means the data do not locate the transient: as &tau; grows past the run
length the term B(1 &minus; e<sup>&minus;t/&tau;</sup>) becomes
&minus;(B/&tau;)t, exactly collinear with v<sub>ss</sub>t, and the two slide along a
valley. As &tau; goes to zero the term becomes a step, B is absorbed into c, and
v<sub>0</sub> = v<sub>ss</sub> &minus; B/&tau; diverges.{{caveat}}</p>
<div class="scroll"><table><tr>{{head}}</tr>{{body}}</table></div>"""
    burst_section = burst_section.format(head=head, body=body, caveat=caveat)

    # --- tables
    condition_columns = [("sample", "sample", "d"), ("[S] mM", "[sub]", ".4g"),
                         ("[H2O2] mM", "[h2o2]", ".4g"), ("[buf] mM", "[buf]", ".4g"),
                         ("pH", "pH", ".2f"), ("[HOO-] mM", "[HOO-]", ".3e"),
                         ("[enz] mM", "[enz]", ".4g"), ("eps", "e", ".4g")]
    head = "".join(f"<th>{html.escape(label)}</th>" for label, _, _ in condition_columns)
    body = ""
    for index, row in enumerate(rows):
        cells = ""
        for label, key, spec in condition_columns:
            value = row[key]
            if key == "sample":
                cells += (f'<td><span class="swatch" style="background:{colour(index)}">'
                          f'</span>{value}</td>')
            else:
                cells += f"<td>{_cell(value, spec)}</td>"
        body += f'<tr class="{"flagged" if row["defects"] else ""}">{cells}</tr>'
    conditions_table = f'<div class="scroll"><table><tr>{head}</tr>{body}</table></div>'

    measure_columns = [("sample", "sample", "d"), ("v0 AU/s", "v0", ".3e"),
                       ("+/-", "v0_stderr", ".1e"),
                       ("vmax AU/s", "vmax", ".3e"), ("+/-", "vmax_stderr", ".1e"),
                       ("vmax at", "vmax_where", ".0%"),
                       ("accel sigma", "accel_z", "+.1f"),
                       ("window rise (quanta)",
                        "window_quanta", ".1f"),
                       ("amplitude AU", "amplitude", ".4f"),
                       ("noise AU", "noise", ".4f"), ("SNR", "snr", ".0f"),
                       ("late/early", "ratio", ".2f"), ("backtrack sigma", "backtrack", ".0f"),
                       ("conversion", "conversion", ".2%"), ("points", "points", "d"),
                       ("duration s", "duration", ".0f")]
    head = ("".join(f"<th>{html.escape(label)}</th>" for label, _, _ in measure_columns)
            + "<th>usable for</th>")
    body = ""
    for row in rows:
        cells = "".join(f"<td>{_cell(row[key], spec)}</td>"
                        for _, key, spec in measure_columns)
        use = row["eligibility"].summary
        tone = "" if use == "rate+shape" else ' style="font-weight:600"'
        cells += f"<td{tone}>{html.escape(use)}</td>"
        body += f'<tr class="{"flagged" if row["defects"] else ""}">{cells}</tr>'
    measures_table = f'<div class="scroll"><table><tr>{head}</tr>{body}</table></div>'

    # --- the verdict list
    items = []
    for row in rows:
        for kind, text in row["notes"]:
            items.append(f'<li class="{kind}"><b>sample {row["sample"]}</b> &mdash; '
                         f'{html.escape(text)}</li>')
        # Eligibility is not a defect and never appears as one. It says what a
        # curve can be USED for -- a slow rung loses its rate and keeps its
        # shape, and deleting it would take the K_M information with it.
        for use, text in row["eligibility"].reasons:
            items.append(f'<li class="power"><b>sample {row["sample"]}</b> &mdash; '
                         f'not usable for <b>{use}</b>: {html.escape(text)}</li>')
    notes_html = (f'<ul class="notes">{"".join(items)}</ul>' if items else
                  '<p class="none">Nothing flagged. Every cuvette rose monotonically, '
                  'cleared its noise, and sits on its experiment\'s ladder.</p>')

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Experiment {number} curves</title>
<style>{CURVE_DOSSIER_STYLE}</style></head><body><div class="wrap">
<h1>{title}</h1>
<p class="sub"><b>{len(rows)}</b> cuvettes &middot; <b>{total_points:,}</b> points &middot;
longest run <b>{longest / 60:.0f}</b> min &middot;
{"catalysed" if catalysed else "enzyme-free"} &middot;
<b>{sum(1 for r in rows if r["defects"])}</b> with defects &middot;
<b>{sum(1 for r in rows if any(k == "shape" for k, _ in r["notes"]))}</b> with shape notes</p>

<h2>All cuvettes</h2>
{overlay_svg}

<h2>Each cuvette, with both fits</h2>
<p class="sub"><b style="color:{LINE_COLOUR}">Dashed red</b> is the straight line
fitted to the first {INITIAL_WINDOW:.0%} of the run &mdash; the one that produces
v<sub>0</sub>. If it does not lie along the start of the curve, the rate for that
cuvette is not trustworthy.
<b style="color:{DOSSIER_BURST_COLOUR}">Solid green</b> is the burst/lag form
A = c + v<sub>ss</sub>t &minus; B(1 &minus; e<sup>&minus;t/&tau;</sup>) fitted to the
<em>whole</em> curve. It is shown as a description of shape only &mdash; no rate on
this page comes from it, for the reason set out under the table below.</p>
<div class="small">{"".join(cards)}</div>

{ladder_svg}

{burst_section}

<h2>What was in each cuvette</h2>
{conditions_table}

<h2>What was measured</h2>
{measures_table}

<h2>Flags</h2>
<p class="sub">Reasons to look, not decisions. Nothing here has been excluded
from any fit unless a purple note says so.
<b style="color:var(--bad)">Red</b> means the curve is not a clean rising
reaction. <b style="color:var(--weak)">Amber</b> means it is faint or far below
its neighbours &mdash; which is what a failed cuvette AND the slow bottom rung
of a titration both look like, so this one always needs the eye.
<b style="color:var(--shape)">Blue</b> is a description of the curve's shape,
not a fault. <b style="color:var(--power)">Green</b> says what a curve cannot be
USED for &mdash; a window rising fewer than {MINIMUM_WINDOW_QUANTA:.0f} of the
instrument's 0.001&nbsp;AU steps cannot carry a rate, whether the cuvette worked
or not. That is a statement about resolution, never about the cuvette, and the
same curve usually remains perfectly good for reading shape.</p>
{notes_html}

<p class="foot">Generated by <code>data/curve_dossier.py</code>.
Absorbance is baseline-subtracted (median of the first {BASELINE_POINTS} points).
Noise is the median absolute second difference, floored at the 0.001 AU
quantisation. Flag thresholds: SNR &lt; {DEAD_CURVE_SNR:.0f}, backtrack &gt;
{BACKTRACK_SIGMA:.0f} sigma, rate &gt; {OUTLIER_FACTOR:.0f}x below the experiment median.<br>
v<sub>0</sub> is the slope of the first {INITIAL_WINDOW:.0%} of the run and
v<sub>max</sub> the slope of the steepest such block; a curve is called
accelerating when the two differ by more than {ACCELERATION_SIGMA:.0f} of their
combined standard errors. That test replaces the old late/early ratio, which
read a finished sigmoid as a deceleration.</p>
</div></body></html>"""


def build(experiments=None, out_directory=OUTPUT_DIRECTORY,
          dataset_path=DATASET_PATH, curve_directory=CURVE_DIRECTORY):
    """Write one page per experiment. Returns the paths written."""
    data = add_solution_columns(pd.read_csv(dataset_path))
    exports = read_all_curves(curve_directory)
    os.makedirs(out_directory, exist_ok=True)

    written = []
    for number in sorted(data.experiment.unique()):
        number = int(number)
        if experiments and number not in experiments:
            continue
        if number not in exports:
            continue
        block = data[data.experiment == number].to_dict("records")
        samples = exports[number]

        rows, curves_by_sample = [], {}
        for row in block:
            sample = int(row["sample"])
            if not 1 <= sample <= len(samples):
                continue
            times, values, source = samples[sample - 1]
            if len(times) < 3:
                continue
            baseline = float(np.median(values[:max(1, min(BASELINE_POINTS, len(values) // 10))]))
            times, values = times - times[0], values - baseline
            # The noise floor is the source's, not a constant: a .rre curve
            # floored at the .txt export's quantisation reports 2.4x its own
            # noise, and every flag on this page is scaled by that number.
            # The same floor goes on to the rates, whose standard errors
            # line_fit floors by it -- see fit_dataset.source_floor.
            floor = source_floor(source)
            noise = curve_noise(values, floor)
            entry = dict(row)
            entry["sample"] = sample
            entry["experiment"] = number
            entry.update(describe(times, values, noise,
                                  float(row["e"] or 0.0), float(row["[sub]"]),
                                  floor=floor))
            rows.append(entry)
            curves_by_sample[sample] = (times, values)
        if not rows:
            continue
        for row in rows:
            row["eligibility"] = eligibility(row)
            row["notes"] = verdicts(row, rows)
            # A shape note describes the curve; it is not a reason to distrust
            # it. Keeping the two apart matters: experiment 26 is the only true
            # replicate set in the enzyme-free data and every one of its
            # cuvettes is sound, but two of them accelerate mildly. Counting
            # that as a fault would train the reader to ignore the flags.
            row["defects"] = [n for n in row["notes"] if n[0] != "shape"]

        path = os.path.join(out_directory, f"exp_{number:03d}.html")
        with open(path, "w") as handle:
            handle.write(render(number, rows, curves_by_sample))
        written.append(path)
    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--experiments", default=None,
                        help="comma-separated experiment numbers; default is all")
    parser.add_argument("--out", default=OUTPUT_DIRECTORY)
    arguments = parser.parse_args()
    wanted = ({int(t) for t in arguments.experiments.split(",")}
              if arguments.experiments else None)
    written = build(wanted, arguments.out)
    total = sum(os.path.getsize(p) for p in written)
    print(f"{len(written)} page(s) -> {arguments.out}/  ({total / 1024:.0f} kB total)")
    for path in written[:10]:
        print(f"  {path}  ({os.path.getsize(path) / 1024:.0f} kB)")
    if len(written) > 10:
        print(f"  ... and {len(written) - 10} more")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
