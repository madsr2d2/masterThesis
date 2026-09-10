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

from svgplot import ACCENT, GRID, MUTED


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

# The forms a progress curve can earn, drawn together on one panel and so
# needing four clearly separated hues rather than four steps of a ramp. Not
# CATEGORY, which is a validated trio. The first attempt used #3f8a5a for the
# quadratic and #12856a for the burst -- both green-teal, and indistinguishable
# on the page. Declared here and not in the folder that draws them, because
# NEVER DECLARE A COLOUR IN A FOLDER: `background_reaction/build_figures.py`
# held all five until 2026-09-04, one of them under a name `curve_dossier` also
# used, for a different colour.
WINDOW_COLOUR = "#2f6fb0"      # the 20% window line          -- blue
QUAD_COLOUR = "#c25e00"        # the whole-curve quadratic    -- amber (headline)
BURST_COLOUR = "#7a4bb8"       # the burst/lag form           -- purple
WHOLE_COLOUR = "#3f8a5a"       # straight line, whole curve   -- green
OUTLIER_COLOUR = "#c0392b"     # ring round a suspect reading -- red

# The window a correction was read through, shaded behind the marks it
# covers -- the two-axis curves page's gas growth-windows, and anything else
# that needs to show WHERE a statistic was measured rather than just quoting
# it in prose. Amber, so it reads as "this span mattered" without competing
# with WINDOW_COLOUR/QUAD_COLOUR/BURST_COLOUR, which are all fit lines.
EVENT_BAND_COLOUR = "#e0a530"

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
.fig svg{display:block}
.fig svg + svg{margin-top:2px}
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


def residual_axes(times, residual, width=340, height=72, pad=(56, 12, 30, 8),
                  colour=ACCENT, bands=(), band_colour=EVENT_BAND_COLOUR):
    """
    A thin strip of (data - fit) / noise against time, to pair beneath a
    `progress_axes` panel.

    Shares `progress_axes`' default left/right padding, so the two plot AREAS
    line up when the two SVGs are stacked in one panel div even though they
    are independent drawings -- a reader should be able to look straight down
    from a feature in the residual to the reading that caused it.

    `bands` is a sequence of `(start, stop)` time pairs shaded the full height
    of the strip -- the window a correction was read through, so a residual
    spike that lines up with one is explained rather than mysterious. Pass the
    same pairs to shade the panel above for the same reason.

    Y-limits are `max(|residual|) * 1.15` each side, symmetric about zero
    because a residual has no natural floor the way a reading does.
    """
    from svgplot import Axes
    times = np.asarray(times, dtype=float)
    residual = np.asarray(residual, dtype=float)
    finite = residual[np.isfinite(residual)]
    span = max(float(np.abs(finite).max()), 1.0) if len(finite) else 1.0
    axes = Axes(width, height, (0, float(times[-1]) * 1.02),
               (-span * 1.15, span * 1.15), pad=pad)
    for start, stop in bands:
        axes.band([start, stop], [-span * 1.15, -span * 1.15],
                 [span * 1.15, span * 1.15], band_colour, opacity=0.16)
    axes.hline(0.0, colour=GRID, dash="2 2", width=1.0)
    dense = len(times) > 150
    axes.line(times, residual, colour, width=0.9, opacity=0.5)
    axes.points(times, residual, colour, radius=1.4 if dense else 1.8,
               opacity=0.8, stroke=None)
    return axes


def derivative_axes(times, progress, width=340, height=72, pad=(56, 12, 30, 8),
                    colour=ACCENT, bands=(), band_colour=EVENT_BAND_COLOUR,
                    samples=300):
    """
    A third strip beneath `residual_axes`: dA/dt of the FITTED curve.

    The rate, not the readings differenced. `progress.rate(times)` is
    `summary_kinetics.ProgressFit`'s own analytic derivative of whichever form
    the curve earned -- a steady rate minus decaying exponentials, so the
    slope is exact and needs no finite differencing of noisy points. That
    matters here specifically: differencing 402 curves' worth of readings
    would draw the noise floor, not the reaction, and would need its own
    smoothing choice on every panel. The fit already made that choice; this
    panel reads it off rather than re-deriving it.

    Evaluated on a SMOOTH GRID (`samples` points from 0 to the run's own
    length), the same way `progress_overlay` draws the curve itself -- not at
    the reading times, which would draw a jagged line wherever readings are
    sparse and hide that the underlying rate is smooth.

    Y-limits include zero, the same convention `progress_axes` uses for the
    readings above it, because a rate that goes negative (a two-phase curve's
    decaying second term can do this briefly) should visibly cross the axis
    rather than be cropped into looking like it never does.

    `bands` and `pad`'s left/right match `residual_axes` exactly, so all
    three panels' plot AREAS line up when stacked in one panel div -- a
    reader should be able to look straight down from a reading, through its
    residual, to the rate the fit assigns it at that instant.
    """
    from svgplot import Axes
    times = np.asarray(times, dtype=float)
    grid = np.linspace(0.0, float(times[-1]), samples)
    rate = np.asarray(progress.rate(grid), dtype=float)
    finite = rate[np.isfinite(rate)]
    lo = min(float(finite.min()), 0.0) if len(finite) else 0.0
    hi = max(float(finite.max()), 1e-9) if len(finite) else 1e-9
    margin = max((hi - lo) * 0.12, 1e-12)
    axes = Axes(width, height, (0, float(times[-1]) * 1.02),
               (lo - margin, hi + margin), pad=pad)
    for start, stop in bands:
        axes.band([start, stop], [lo - margin, lo - margin],
                 [hi + margin, hi + margin], band_colour, opacity=0.16)
    axes.hline(0.0, colour=GRID, dash="2 2", width=1.0)
    axes.line(grid, rate, colour, width=1.3)
    return axes


def breakpoints(axes, where, labels=None, colour=MUTED, row=0):
    """
    EVERY landmark the curve earned, labelled, not just the first.

    A rate that rises and then falls has two breakpoints, and drawing one
    leaves the other invisible -- which is how the early break on the 40 C
    curves went unnoticed until it was seen by eye on a page like this. The
    same drawing serves an induction landmark or a window edge; `labels`
    defaults to the ordinal.

    `row` drops the LABELS one line, leaving the rules where they are. A panel
    that marks several different parameters calls this once per parameter, in
    each one's own colour, and two landmarks a minute apart would otherwise
    write their names on top of each other -- which is the failure mode of
    drawing more than one thing and the reason a panel used to draw only
    `v_max`.

    LABELS THAT CARRY A VALUE CAN OVERFLOW THE RIGHT EDGE, where a bare
    ordinal or "v_max" never did -- "v_max* 8.67e-05" is four times the
    width. Past 60% of the plot a label is right-anchored and grows back
    towards the line instead of off the panel; `note`'s `anchor` was already
    there for this, just never driven by the mark's own position.
    """
    x0, x1 = axes.left, axes.width - axes.right
    for index, cut in enumerate(where):
        x = axes._fx(cut)
        axes.parts.append(
            f"<path d='M{x:.2f},{axes.top} "
            f"L{x:.2f},{axes.height - axes.bottom}' "
            f"stroke='{colour}' stroke-width='1.1' "
            f"stroke-dasharray='3 3' fill='none'/>")
        text = f"{index + 1}" if labels is None else labels[index]
        if not text:
            continue
        late = (x - x0) > 0.6 * (x1 - x0)
        axes.note(x - 3 if late else x + 3, axes.top + 10 + 11 * row, text,
                  colour, size=9.5, anchor="end" if late else "start")


def stacked_landmarks(axes_list, where, labels=None, colour=MUTED, row=0):
    """
    The same vertical landmark(s), drawn on every axes in a stack.

    A panel showing readings, residual and derivative is three independent
    SVGs, each with its own `breakpoints` call -- so a dashed rule marking
    `v_max` or `tau` used to stop at the bottom of the top panel, where the
    reader's eye keeps going. `progress_axes`/`residual_axes`/
    `derivative_axes` all promise the same left/right padding so their plot
    AREAS line up when stacked; drawing the same rule at the same `where` in
    each one lines the columns up into what reads as one continuous line
    running through the whole panel, even though `breakpoints` runs once per
    axes.

    Only the FIRST axes in `axes_list` gets the label -- `breakpoints`
    already skips its own label when passed an empty string, so the rest
    draw the bare rule only.
    """
    blank = [""] * len(where)
    for index, axes in enumerate(axes_list):
        breakpoints(axes, where, labels if index == 0 else blank,
                   colour=colour, row=row)


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
    from svgplot import clipped_marks, colliding_clips
    faults = 0
    for name, content in pages.items():
        path = os.path.join(directory, name)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(content)
        clipped = clipped_marks(content)
        # The same failure one level up: a clip id bound to two rectangles
        # sends a figure's marks through another figure's frame, and the marks
        # it drops are invisible in exactly the way a clipped mark is. See
        # svgplot.colliding_clips -- this was live on five of the six index
        # pages, and `clipped_marks` is structurally unable to see it.
        colliding = colliding_clips(content)
        faults += len(clipped) + len(colliding)
        print(f"wrote {path}  ({len(content) / 1024:.0f} kB)  "
              f"{len(clipped)} clipped marks"
              + ("" if not clipped else f"  {clipped[:4]}")
              + ("" if not colliding
                 else f"  {len(colliding)} COLLIDING CLIP IDS {colliding[:2]}"))
    return 1 if faults else 0
