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
