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
import dataclasses
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

# The same device for the OTHER direction the gas moves. A detachment span and
# an arrival span are both "this stretch is the artefact, not the chemistry",
# so an arrival gets a band too rather than a bare rule -- but the two have to
# be told apart at a glance, and at band opacity a hue difference reads far
# better than a lightness one. Cool against the amber's warm.
#
# BLUE AND NOT VIOLET, SINCE 2026-09-12. It was #8a5aa8, which is CATEGORY[2]
# to the byte -- the colour the reconstruction's own line and fit are drawn
# in on every curves page. Nothing caught it because no folder drew an arrival
# at all until `progress_panel`; the constant sat in this file, used only by
# `scratch/build_detection_review.py`, for as long as it was wrong. A band and
# a line of one hue on one panel is exactly the collision the four form
# colours above were separated to avoid.
ARRIVAL_BAND_COLOUR = "#2f6fb0"

# The two series a progress panel draws, and the two fits over them. Named
# here rather than reached for as CATEGORY[n] in a folder, because these are
# not three unordered categories -- they are "the readings" and "the readings
# with the gas taken out", and a reader has to carry the pairing between
# panels and between folders.
FIT_COLOUR = ACCENT                # the fit to the readings as recorded
CHEMISTRY_COLOUR = "#8a5aa8"       # the rebuilt series, and the fit to it

# A quiet tail -- the stretch after the last detachment, where `debubble`
# subtracts nothing because nothing was watched to leave. Past about one
# shedding interval the run may simply have ended mid-bubble, and there the
# reconstruction keeps gas it should have removed. Grey, and the lightest
# wash on the panel: it marks where the correction STOPS, which is a weaker
# statement than either band.
QUIET_TAIL_COLOUR = MUTED
QUIET_TAIL_INTERVALS = 1.0

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
.ph .detail{font-weight:400}
/* The panel sits on a FIXED light background whatever the page theme (the
   same reason the SVGs draw in fixed hex, not theme tokens) -- so its own
   text needs fixed colours too. `.ph`/`.ps`/`.pf` used to inherit
   var(--ink)/var(--muted), which flip light in dark mode and then sit on
   this always-light panel: pale grey on off-white, unreadable. `.cap` was
   already fixed for the same reason; these three were the gap. */
.fig .ph{color:#1a1a1a}
.fig .ps,.fig .pf{color:#5a5a5a}
/* Panel text wraps to the SAME width the SVG renders at (340px, the
   progress-panel default), not the wider grid column it sits in -- a
   `.fig.panel` used to size to its column and let the header wrap wider
   than the 340px chart beneath it. */
.fig.panel{max-width:340px}
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


def draw_progress(axes, times, progress, colour=FIT_COLOUR, width=FIT_WIDTH,
                  mark_radius=None, samples=300):
    """
    Draw a fit that has ALREADY BEEN MADE. This function does not fit.

    IT TAKES THE FIT RATHER THAN MAKING ONE, since 2026-09-12, and that is the
    whole point of it. It was `progress_overlay(axes, times, values)` and it
    called `summary_kinetics.fit_progress` itself -- so every panel in the
    repository refitted a curve `scope._frame` had already fitted, and nothing
    could assert that the line a panel DREW came from the numbers its footer
    QUOTED. The gap was live: panels marked a vertical labelled "tau" off
    `fit_burst_bounded`, the ONE-phase form, on 36 of the two-axis block's 110
    live curves whose drawn line was the TWO-phase form. Pass
    `scope.curve_fit(curve).progress` (or `.progress_corrected`) and the two
    cannot come apart.

    WHICHEVER FORM THE CURVE EARNED. Drawing the one-phase fit on a two-phase
    curve is how the shape stayed invisible for as long as it did -- the panel
    showed a fit that could not bend the way the readings do, and it read as
    scatter. `fit_progress` fits both and returns the one an F test chooses.

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


def breakpoints(axes, where, labels=None, colour=MUTED, row=0, dash="3 3"):
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

    `dash` is the rule's own stroke pattern, `None` for a solid line. It is
    how two KINDS of landmark share one colour: the gas marks draw arrivals
    solid and detachments dashed, so a reader tells "the beam took gas on"
    from "the beam shed it" by line style and does not have to hold a second
    hue in mind to do it.

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
            + (f"stroke-dasharray='{dash}' " if dash else "")
            + f"fill='none'/>")
        text = f"{index + 1}" if labels is None else labels[index]
        if not text:
            continue
        late = (x - x0) > 0.6 * (x1 - x0)
        axes.note(x - 3 if late else x + 3, axes.top + 10 + 11 * row, text,
                  colour, size=9.5, anchor="end" if late else "start")


def stacked_landmarks(axes_list, where, labels=None, colour=MUTED, row=0,
                      dash="3 3"):
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
                   colour=colour, row=row, dash=dash)


# The two forms `summary_kinetics.fit_progress` chooses between, as a reader
# sees them in a panel's own header -- named once so a panel never states
# "the fitted curve" without saying which model that is. Keys are
# `row.phases`.
FITTED_FUNCTION = {
    1: "A(t) = c + v_ss·t − B(1−e^(−t/τ))",
    2: "A(t) = c + v_ss·t − B₁(1−e^(−t/τ₁)) − B₂(1−e^(−t/τ₂))",
}


def panel_header(experiment, sample, pH, conditions, phases):
    """
    One line, one weight, only the experiment label bold: "Exp. X.Y : pH =
    ..., cond, cond, ..., <fitted function>".

    `conditions` is a sequence of already-formatted "[name] = value unit"
    strings, folder-specific (a pH ladder's substrate/peroxide/buffer, the
    two-axis block's own). `.ph .detail` (figure_kit's own CSS) carries the
    normal weight; the surrounding `.ph` class is bold by default, which is
    what leaves only "Exp. X.Y" bold.
    """
    detail = ", ".join([f"pH = {pH:.2f}", *conditions,
                        FITTED_FUNCTION[int(phases)]])
    return (f"<strong>Exp. {int(experiment)}.{int(sample)}</strong> "
           f"<span class='detail'>: {detail}</span>")


def bubble_labels(count, held):
    """
    One "bubble" label per O2 detachment, for `stacked_landmarks`.

    `held` (`row.terminal_gas > 0`) is a DIFFERENT fact from a detachment --
    the run ended before the beam's own last bubble released, so `debubble`'s
    correction past that point is a lower bound -- and rides the LAST
    detachment's line rather than replacing its label, since it is not
    itself a detachment.
    """
    labels = ["bubble"] * count
    if held and count:
        labels[-1] = "bubble · gas held"
    return labels


def progress_axes(times, values, width=340, height=210, limit=None,
                  colour=MUTED, pad=(56, 12, 34, 20), companion=None,
                  keep=()):
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
             else decimated(len(times), limit, keep=keep))
    axes.points(times[shown], values[shown], colour, radius=radius,
                opacity=0.75 if dense else 0.9,
                stroke=None if dense else "white", stroke_width=0.6)
    return axes, radius


@dataclasses.dataclass(frozen=True)
class Mark:
    """
    A landmark that is NOT a parameter of the fit the panel draws.

    `v_max` is the steepest 20% block slope of the raw readings
    (`curve_metrics.peak_rate`); `t_ind` is the first crossing of half the
    largest rolling slope through a window a tenth of the run wide; a
    breakpoint is a piecewise-linear segmentation; `tail_start` is a rolling
    rate's maximum. Not one of them appears in
    `A(t) = c + v_ss.t - sum B_i(1 - e^(-t/tau_i))`, which is the curve the
    panel actually draws.

    THEY ARE STILL WORTH DRAWING -- each is a quantity some folder's document
    is read off, and a curves page is the audit surface for exactly that. What
    they may not do is arrive UNLABELLED beside the fit's own clocks and read
    as if the fit produced them. So every one carries `source`, the estimator
    it came from, and `label` renders as "v_max (20% block)". The fitted
    parameters need no such tag: they are named in the header, in the fitted
    function the header prints.
    """
    time: float
    label: str
    source: str
    colour: str = None

    @property
    def text(self):
        return f"{self.label} ({self.source})"


@dataclasses.dataclass(frozen=True)
class Estimator:
    """
    ANOTHER estimator of the same curve, drawn over it for comparison.

    `Mark` is a landmark in TIME and this is a whole drawn series -- the
    quadratic through every point, the straight line through every point, the
    least-squares line over the first 20%, the burst form. `background_reaction/`
    draws four of them, because that folder's question is precisely "does the
    conclusion depend on how the rate was measured", and the answer has to be
    visible on the curve and not only in a table.

    IT CARRIES `source` FOR THE SAME REASON `Mark` DOES. None of these is the
    relaxation form the header prints, and four unlabelled lines on one panel
    is how a reader comes to believe the fit is whichever line is nearest.

    `low`/`high` draw a band instead of nothing where an estimator has an
    interval -- the burst form's v0 profile fan, which is the whole of that
    folder's identifiability argument: where the fan is wide the curve does
    not determine v0, however well the form fits.
    """
    times: np.ndarray
    values: np.ndarray
    colour: str
    label: str
    source: str
    dash: str = None
    width: float = 1.6
    low: np.ndarray = None
    high: np.ndarray = None
    band_times: np.ndarray = None

    @property
    def text(self):
        return f"{self.label} ({self.source})"


@dataclasses.dataclass(frozen=True)
class Region:
    """
    A span of the run a folder's statistic was read THROUGH, shaded.

    `slowdown.sink_fit` reads its decline over a tail starting at the rolling
    rate's maximum; `induction.buffer_landmark` reads through a 450 s window
    chosen because two runs differ in length. CLAUDE.md's rule is that the
    window a statistic was read through gets drawn and gets taken FROM THE
    FIT -- `sink_fit.tail_start` exists precisely so a page cannot guess a
    tail and then draw a different window from the one that was fitted.

    Same device as the gas bands, deliberately: both say "this stretch is
    where something other than the plain reaction is being read".
    """
    start: float
    stop: float
    colour: str
    opacity: float = 0.09


def _asymptote(fit, times):
    """
    The straight line the fitted curve settles onto: `c - sum B + v_ss.t`.

    THREE OF THE FOUR FITTED PARAMETERS, IN ONE LINE. Its slope is `v_ss`, and
    its vertical distance from the fit at t = 0 is exactly `sum B_i` -- so a
    reader gets the steady rate and the amplitude off the panel by looking,
    without either needing a mark of its own. A vertical rule was never the
    right device for `v_ss`, which is a slope and not a time; that is half of
    why the panels ended up marking `v_max` instead, which IS a time and is
    also not in the model.

    Returns (y0, y1) at the ends of `times`, or None where the form did not
    resolve.
    """
    chosen = fit.chosen
    amplitudes = [b for b in fit.amplitudes if np.isfinite(b)]
    if not np.isfinite(chosen.c) or not np.isfinite(chosen.v_ss):
        return None
    if not amplitudes:
        return None
    intercept = float(chosen.c) - float(sum(amplitudes))
    v_ss = float(chosen.v_ss)
    return (intercept, intercept + v_ss * float(times[-1]))


def _clocks(fit):
    """
    The drawn fit's own time constants, as (time, label) pairs.

    ONE PHASE GIVES ONE CLOCK AND TWO GIVE TWO, and the labels are the symbols
    the header's fitted function uses -- so tau_1 on the panel is tau_1 in
    `A(t) = c + v_ss.t - B_1(1 - e^(-t/tau_1)) - B_2(1 - e^(-t/tau_2))` above
    it, and a reader can pair them without being told.

    tau_2 is drawn only where `TwoPhaseFit.resolved` says its 95% profile
    interval stays inside the grid. That flag says whether tau_2 is QUOTABLE,
    not whether the second phase is real -- the two are different questions
    and confusing them is what `model_residual` exists to prevent -- but a
    rule drawn at an unlocated time constant asserts a location, so
    resolution is the right gate for DRAWING it.
    """
    if fit.phases == 2:
        two = fit.two
        found = []
        if np.isfinite(two.tau1):
            found.append((float(two.tau1), "τ₁"))
        if two.resolved and np.isfinite(two.tau2):
            found.append((float(two.tau2), "τ₂"))
        return found
    one = fit.one
    if np.isfinite(one.tau) and one.tau > 0:
        return [(float(one.tau), "τ")]
    return []


def progress_panel(fit, row, conditions, subhead, footer="", marks=(),
                   regions=(), estimators=(), rings=(), table="",
                   limit=140, width=340):
    """
    ONE progress curve, drawn the ONE way. Every folder's curves page uses it.

    There were eight of these, ~890 lines, and they marked eight different
    vocabularies on the same object: `two_axis` and `ph` drew v_max/v_max*/tau,
    `buffer` and `induction` drew t_ind, `product_fate` drew a rolling rate's
    maximum, `early_trough` a trough, `temperature_series` piecewise-linear
    breakpoints, and `background_reaction` four estimator lines on axes it
    built itself. `two_axis` and `ph` were ~70 lines of byte-comparable code
    differing in line wrapping. None of them marked a parameter of the function
    they drew.

    WHAT IS DRAWN, IN ORDER:

      the readings              grey marks, thinned to `limit`
      the gas                   detachment spans amber, arrival spans blue,
                                a quiet tail greyed -- on all three strips
      the rebuilt series        dashed, where there is gas to take out
      both fits                 the readings' in rust, the rebuilt one's in
                                purple; the rebuilt fit is the CHEMISTRY and
                                the residual strip is read against it
      the fit's own parameters  tau (and tau_2 where resolved) as rules; v_ss
                                and sum B as the asymptote they define
      v_peak                    on the DERIVATIVE strip, at its maximum, where
                                it is exact rather than guessed at
      `estimators`              other whole series drawn for comparison, UNDER
                                the fit, each naming itself
      `regions`                 a window a folder's statistic was read through
      `rings`                   suspect readings, ringed and never removed
      `marks`                   anything else, each naming its own estimator

    `fit` is a `scope.CurveFit` and `row` its frame row. They are the same
    computation -- see `scope.curve_fit` -- so the rule this draws and the
    number the footer prints cannot disagree.
    """
    times, values = fit.times, fit.values
    chem_values, chem_fit = fit.chemistry
    chopped = bool(fit.events)
    chem_colour = CHEMISTRY_COLOUR if chopped else FIT_COLOUR

    # RINGS ARE NOT THINNED, so anything ringed is forced back into the drawn
    # stride -- a ring landing where no mark was drawn reads as a flag on
    # nothing, which is what `decimated`'s `keep` exists for.
    rings = sorted(int(i) for i in rings)
    axes, radius = progress_axes(times, values, limit=limit, width=width,
                                 companion=fit.corrected if chopped else None,
                                 pad=(56, 12, 10, 20), keep=rings)
    # THE GAS FIRST, so every later mark sits on top of it. Bands are spans and
    # not rules because these are events with a WIDTH: one bubble can take
    # several 60 s readings to leave, and a rule at one reading says it was
    # instantaneous. Detachments and arrivals both, since 2026-09-12 -- the
    # arrival half of the detection layer was drawn on no published page at
    # all, though it is the correction that moved `joint_clocks`' headline
    # row to +0.757 +/- 0.289.
    detach_bands = [(float(times[a]), float(times[b])) for a, b in fit.events]
    arrive_bands = [(float(times[a]), float(times[b])) for a, b in fit.arrivals]
    quiet_band = []
    if fit.events and fit.quiet_intervals > QUIET_TAIL_INTERVALS:
        quiet_band = [(float(times[fit.events[-1][1]]), float(times[-1]))]
    washes = [(lo, hi, QUIET_TAIL_COLOUR, 0.07) for lo, hi in quiet_band]
    # A folder's own window goes under the gas, which is the artefact and has
    # to stay the most visible thing on the panel.
    washes += [(r.start, r.stop, r.colour, r.opacity) for r in regions
               if np.isfinite(r.start) and np.isfinite(r.stop)]
    washes += [(lo, hi, EVENT_BAND_COLOUR, 0.14) for lo, hi in detach_bands]
    washes += [(lo, hi, ARRIVAL_BAND_COLOUR, 0.12) for lo, hi in arrive_bands]
    for lo, hi, colour, opacity in washes:
        axes.band([lo, hi], [axes.ylim[0], axes.ylim[0]],
                  [axes.ylim[1], axes.ylim[1]], colour, opacity=opacity)

    # OTHER ESTIMATORS UNDER THE FIT, never over it. The panel's own claim is
    # the relaxation form its header prints; a comparison line drawn on top of
    # it would make whichever estimator was passed last look like the fit.
    for other in estimators:
        if other.low is not None and other.high is not None:
            axes.band(other.band_times if other.band_times is not None
                      else other.times,
                      other.low, other.high, other.colour, opacity=0.20)
        axes.line(other.times, other.values, other.colour, width=other.width,
                  dash=other.dash)

    if chopped:
        axes.line(times, fit.corrected, CHEMISTRY_COLOUR, width=1.0,
                  dash="3 2", opacity=0.85)
        draw_progress(axes, times, fit.progress_corrected,
                      colour=CHEMISTRY_COLOUR, mark_radius=radius)
    draw_progress(axes, times, fit.progress, colour=FIT_COLOUR,
                  mark_radius=radius)

    # SUSPECT READINGS ARE RINGED, NEVER REMOVED. Every fit above is computed
    # on every point; a ring says which ones a reader should discount by eye.
    # `curve_metrics.isolated_outliers` nominates them, and it is also what
    # nominates the excursions the detection layer REJECTED as gas -- those
    # stay in the curve on purpose (BUBBLES.md), so a page that draws the gas
    # and not these is showing only the half that was acted on.
    for index in rings:
        axes.ring(times[index], values[index], OUTLIER_COLOUR,
                  title=f"suspect reading: point {index} at "
                        f"t={times[index]:.0f} s")

    residual = (chem_values - chem_fit.predict(times)) / fit.noise
    bands = detach_bands + arrive_bands
    rax = residual_axes(times, residual, colour=chem_colour, bands=bands)
    drax = derivative_axes(times, chem_fit, colour=chem_colour, bands=bands)
    stack = [axes, rax, drax]

    # THE ASYMPTOTE: v_ss and sum B, read off one line. Drawn under the rules
    # so a clock's label is never crossed by it, and clipped by the plot frame
    # the way any extrapolation leaving the axes is -- on a deep lag it starts
    # well below the readings, which is the shape being visible rather than a
    # fault.
    ends = _asymptote(chem_fit, times)
    if ends is not None:
        axes.line([0.0, float(times[-1])], list(ends), chem_colour,
                  width=1.0, dash="1 3", opacity=0.7)

    # THE FIT'S OWN CLOCKS, in the colour of the fit they belong to. Nothing
    # here is labelled with its estimator, because the estimator is the fitted
    # function printed in the header.
    clocks = _clocks(chem_fit)
    for index, (when, name) in enumerate(clocks):
        if 0 < when < float(times[-1]):
            stacked_landmarks(stack, [when], [f"{name} {when:.0f} s"],
                              colour=chem_colour, row=index)
    used_rows = len(clocks)

    # v_peak ON THE DERIVATIVE STRIP, which is the only panel it is exact on:
    # it IS the maximum of that curve. Marked there rather than as a vertical
    # through the readings, where a rate has no position to point at.
    peak, peak_at = chem_fit.peak_rate
    if np.isfinite(peak) and np.isfinite(peak_at) and 0 < peak_at < times[-1]:
        drax.points([peak_at], [peak], chem_colour, radius=3.0, stroke="white",
                    stroke_width=0.8)
        drax.note(drax._fx(peak_at) + 4, drax._fy(peak) - 4,
                  f"v_peak {peak:.2e}", chem_colour, size=9.0)

    # THE GAS, LABELLED. One label per detachment and one per arrival, so the
    # two directions read apart at a glance: a detachment is dashed, an
    # arrival solid, which is the convention `breakpoints`' own `dash`
    # argument exists for.
    if fit.events:
        stacked_landmarks(
            stack, [float(times[a]) for a, _ in fit.events],
            bubble_labels(len(fit.events), fit.terminal_held > 0),
            colour=GRID, row=used_rows)
        used_rows += 1
    if fit.arrivals:
        stacked_landmarks(
            stack, [float(times[a]) for a, _ in fit.arrivals],
            ["gas in"] * len(fit.arrivals), colour=ARRIVAL_BAND_COLOUR,
            row=used_rows, dash=None)
        used_rows += 1

    # AND ANYTHING ELSE THE FOLDER READS ITS DOCUMENT OFF, each naming the
    # estimator it came from. See `Mark`.
    for offset, mark in enumerate(marks):
        if not np.isfinite(mark.time) or not 0 < mark.time < times[-1]:
            continue
        stacked_landmarks(stack, [float(mark.time)], [mark.text],
                          colour=mark.colour or MUTED, row=used_rows + offset)

    return panel(
        panel_header(fit.experiment, fit.sample, row.pH, conditions,
                     chem_fit.phases),
        subhead,
        axes.render("", "ΔA", xticks=False)
        + rax.render("", "z", xticks=False)
        + drax.render("time, s", "dA/dt"),
        footer, table)


def curves_page_title(block):
    """
    Every curves page's title, one shape: "<block> — every progress curve".

    Four of the eight already read that way and four did not -- "Progress
    curves and their fits", "Temperature series — all 24 progress curves",
    "The 4OMe curves and their induction landmarks", "The early trough —
    every curve behind the finding". A reader moving between folders should
    not have to work out whether a page is the audit surface or an argument;
    `index.html` presents the argument and this page shows the fits it is read
    off, in every folder, and the title is where that promise is made. The
    count belongs in the subtitle, where it can change without the title
    changing.
    """
    return f"{block} — every progress curve"


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
