"""
The style every analysis folder's `build_figures.py` draws in.

`svgplot` is the drawing primitives -- axes, marks, lines, the page shell. This
is the layer above: the palettes, the figure wrapper and the page writer that
were declared separately in each of the five folders. As with `doc_check`, the
copies had already drifted where it mattered least visibly: `TEMPERATURES` was
a different six-step ramp in `induction/` than in `temperature_series/`, so the
same ordinal variable was drawn two ways in two folders of one thesis, and only
one of the two carried the contrast claim.

    from figure_kit import CATEGORY, EXTRA_CSS, RUNGS, SURFACE, fig, write_pages

Nothing here computes anything. A number reaches a figure from `scope`,
`curve_metrics` or an analysis module, never from here.
"""
import os

import numpy as np

from svgplot import ACCENT, MUTED


# ORDERED VARIABLES GET SEQUENTIAL RAMPS, not categorical hues. Substrate rung
# and temperature are both ordinal, so a light-to-dark single hue carries the
# order; cycling categorical hues over them would throw that away. Every step
# clears 3:1 against the figure surface below, and every line is direct-labelled
# as well, so identity is never colour alone.
RUNGS = ["#6295c3", "#3d729f", "#1e5079", "#0c2f4d"]
TEMPERATURES = ["#7fa9cd", "#5d8fba", "#3c74a4", "#255a8a", "#12406b", "#062a4b"]

# Two steps, for a variable split either side of a threshold -- the buffer
# titrations above and below phosphate's pKa. The same ramp as RUNGS at its
# ends, so a reader moving between folders reads light-to-dark the same way.
PH_RAMP = ["#6295c3", "#0c2f4d"]

# The one genuinely categorical set: three unordered things (three enzyme
# hypotheses, three fitted parameters, three species tests). Validated for
# colour-vision deficiency separation and contrast in both light and dark.
CATEGORY = ["#2f6fb0", "#c0522a", "#8a5aa8"]

# Figures sit on a fixed light surface whatever the page theme, so the ramps'
# contrast is deterministic. A sequential ramp cannot clear 3:1 against a white
# AND a near-black surface at once -- it needs the lightness range the contrast
# rule would spend.
SURFACE = "#fbfbfa"

# Narrower than the smallest mark diameter used on a progress panel, so the fit
# never covers a reading whole. `temperature_series/check_numbers.py` asserts
# the margin on every panel it draws, because it is easy to break by nudging a
# radius and invisible once broken.
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
    """One figure: the drawing, its caption, and anything that follows it."""
    return (f"<div class='fig'>{svg}"
            f"<div class='cap'>{caption}</div>{extra}</div>")


def decimated(count, limit, keep=()):
    """
    The indices to draw for a long run, with `keep` forced in whatever the stride.

    A 400-reading curve drawn whole is a solid bar, so long runs are thinned to
    about `limit` marks. RINGS ARE NOT THINNED, so anything ringed has to be
    unioned back in or a ring lands where no mark was drawn and reads as a flag
    on nothing -- which it did: one curve ringed reading 85 and another reading
    225, both odd, both skipped by a stride of 2.
    """
    step = max(1, count // limit)
    return sorted(set(range(0, count, step)) | set(int(i) for i in keep))


def progress_overlay(axes, times, values, colour=ACCENT, width=FIT_WIDTH,
                     mark_radius=None, samples=300):
    """
    Draw whichever form the curve earned, and return the fit.

    WHICHEVER FORM IT EARNED. Drawing the one-phase fit on a two-phase curve is
    how the shape stayed invisible for as long as it did -- the panel showed a
    fit that could not bend the way the readings do, and it read as scatter.
    `summary_kinetics.fit_progress` fits both and returns the one an F test
    chooses.

    DATA FIRST, FIT ON TOP, AND THE FIT NARROWER THAN THE MARKS. Three passes
    to get this right, so the constraint is written down rather than re-derived:

      1. fit under the data   -> 368 readings bury the fit
      2. fit over the data, 2.0 px wide on a 3.6 px white halo -> the halo is
         wider than a mark, so wherever the curve is tight the fit erases the
         very points it is fitting
      3. this: a light scatter with a thin rust line over it.

    Separation is by HUE and LIGHTNESS, not by a halo. Pass `mark_radius` and
    the constraint is asserted rather than trusted; it is easy to break by
    nudging a radius and invisible once broken.
    """
    from summary_kinetics import fit_progress          # local: heavy import
    progress = fit_progress(times, values)
    smooth = np.linspace(0, float(times[-1]), samples)
    if mark_radius is not None:
        assert width < 2 * mark_radius, "fit line would cover the marks"
    axes.line(smooth, progress.predict(smooth), colour, width=width)
    return progress


def breakpoints(axes, where, labels=None, colour=MUTED):
    """
    EVERY landmark the curve earned, labelled, not just the first.

    A rate that rises and then falls has two breakpoints, and drawing one
    leaves the other invisible -- which is how the early break on the 40 C
    curves went unnoticed until it was seen by eye on a page like this. The
    same drawing serves an induction landmark or a window edge; `labels`
    defaults to the ordinal.
    """
    for index, cut in enumerate(where):
        x = axes._fx(cut)
        axes.parts.append(
            f"<path d='M{x:.2f},{axes.top} "
            f"L{x:.2f},{axes.height - axes.bottom}' "
            f"stroke='{colour}' stroke-width='1.1' "
            f"stroke-dasharray='3 3' fill='none'/>")
        text = f"{index + 1}" if labels is None else labels[index]
        axes.note(x + 3, axes.top + 10, text, colour, size=9.5)


def progress_axes(times, values, width=340, height=210, limit=None,
                  colour=MUTED, pad=(56, 12, 34, 20), companion=None):
    """
    Axes over one progress curve with its readings already drawn on.

    The limits are the convention the two existing curve pages settled on: from
    t = 0, and a y range that always includes zero so a curve that starts
    negative is not cropped into looking like it starts at its own minimum.

    `limit` thins a long run to about that many marks (see `decimated`); left
    None every reading is drawn, which is right for a 24-panel page and wrong
    for a 200-panel one. Returns the axes and the mark radius, so the caller
    can pass the radius to `progress_overlay` and have the width asserted.

    `companion` is a SECOND series the caller will draw itself, given here only
    so the limits hold it. A panel that draws a transformed copy of the curve
    beside the readings -- `two_axis/` draws the O2-corrected series beside the
    raw one -- must not let that copy fall off the frame, and a mark drawn
    outside the limits vanishes silently while the figure still looks finished.
    Nothing is drawn from it; only `values` gets marks.
    """
    from svgplot import Axes
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    extent = (values if companion is None
              else np.concatenate([values, np.asarray(companion, dtype=float)]))
    axes = Axes(width, height, (0, float(times[-1]) * 1.02),
                (min(float(extent.min()), 0.0) * 1.1 - 1e-4,
                 max(float(extent.max()), 1e-4) * 1.12), pad=pad)
    dense = len(times) > 150
    radius = 1.6 if dense else 2.1
    shown = (slice(None) if limit is None
             else decimated(len(times), limit))
    axes.points(times[shown], values[shown], colour, radius=radius,
                opacity=0.75 if dense else 0.9,
                stroke=None if dense else "white", stroke_width=0.6)
    return axes, radius


def panel(header, subhead, svg, footer="", table=""):
    """
    One cuvette as an HTML block: heading, conditions, the plot, its numbers.

    The numbers were once floating text inside the SVG. They sat on top of the
    curves, were clipped by the frame, and could not be selected or searched.
    Anything textual belongs in HTML; the SVG draws only data and fits.
    """
    return (f"<div class='fig panel'>"
            f"<div class='ph'>{header}</div>"
            f"<div class='ps'>{subhead}</div>{svg}"
            + (f"<table class='nums'>{table}</table>" if table else "")
            + (f"<div class='pf'>{footer}</div>" if footer else "")
            + "</div>")


def styled(title, body, subtitle=""):
    """`svgplot.page` with the shared figure CSS appended."""
    from svgplot import page
    return page(title, body, subtitle).replace("</style>",
                                               EXTRA_CSS + "</style>")


def write_pages(directory, pages):
    """
    Write each page and report what it drew outside its own frame.

    THE CLIP REPORT IS THE POINT. Marks are clipped to the plot area
    deliberately, so a data point outside the axis limits vanishes silently and
    the figure still looks complete -- which is how a reading went missing from
    a curvature figure for as long as that figure existed. Two of the five
    folders reported it at build time and three did not, so in three folders the
    only warning came from `check_numbers.py` minutes later, if at all.

    `pages` maps a file name to its content. Returns an exit code, non-zero if
    anything was clipped, so a builder can be run as a gate.
    """
    from svgplot import clipped_marks
    clipped_total = 0
    for name, content in pages.items():
        path = os.path.join(directory, name)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(content)
        clipped = clipped_marks(content)
        clipped_total += len(clipped)
        print(f"wrote {path}  ({len(content) / 1024:.0f} kB)  "
              f"{len(clipped)} clipped marks"
              + ("" if not clipped else f"  {clipped[:4]}"))
    return 1 if clipped_total else 0
