"""
A small inline-SVG plotting helper, in the style data/curve_dossier.py uses.

No matplotlib, no JavaScript, no external files: a page opens instantly, zooms
without pixelating, and survives being emailed or committed. This module draws;
it never computes a measurement. Every number that reaches a figure comes from
`scope`, `curve_metrics` or `summary_kinetics`.
"""
import hashlib
import html

import numpy as np

PALETTE = ["#2f6fb0", "#c0522a", "#3f8a5a", "#8a5aa8", "#b08a2f",
           "#4a9ab0", "#a03a5a", "#6a7a3a"]
INK = "#1a1a1a"
MUTED = "#6b6b6b"
GRID = "#d8d8d8"
ACCENT = "#c0522a"


def esc(text):
    return html.escape(str(text))


class Axes:
    """A linear-or-log axis pair mapping data coordinates to SVG pixels."""

    def __init__(self, width, height, xlim, ylim, xlog=False, ylog=False,
                 pad=(62, 16, 46, 14)):
        self.width, self.height = width, height
        self.left, self.right, self.bottom, self.top = pad
        self.xlog, self.ylog = xlog, ylog
        self.xlim = tuple(np.log10(v) if xlog else v for v in xlim)
        self.ylim = tuple(np.log10(v) if ylog else v for v in ylim)
        self.parts = []      # marks: clipped to the plot area
        self.overlay = []    # text: never clipped, always drawn on top

    def _fx(self, value):
        value = np.log10(value) if self.xlog else np.asarray(value, dtype=float)
        span = self.xlim[1] - self.xlim[0] or 1.0
        inner = self.width - self.left - self.right
        return self.left + (value - self.xlim[0]) / span * inner

    def _fy(self, value):
        value = np.log10(value) if self.ylog else np.asarray(value, dtype=float)
        span = self.ylim[1] - self.ylim[0] or 1.0
        inner = self.height - self.top - self.bottom
        return self.height - self.bottom - (value - self.ylim[0]) / span * inner

    # --- marks -----------------------------------------------------------
    def points(self, x, y, colour, radius=3.2, opacity=1.0, title=None,
               stroke="white", stroke_width=0.8):
        """
        `stroke=None` for a DENSE series. The white ring separates overlapping
        marks and is right for a dozen points; on a progress curve of 368 it
        inverts -- each mark's ring covers its neighbour's fill, so a run of
        closely spaced readings turns into a band of white with dots in it.
        Exps 14 and 19 sample 368 readings into 270 px, about 0.7 px apart.
        """
        ring = (f"stroke='{stroke}' stroke-width='{stroke_width}'"
                if stroke else "stroke='none'")
        for xi, yi in zip(np.atleast_1d(x), np.atleast_1d(y)):
            if not (np.isfinite(xi) and np.isfinite(yi)):
                continue
            tip = f"<title>{esc(title)}</title>" if title else ""
            self.parts.append(
                f"<circle cx='{self._fx(xi):.2f}' cy='{self._fy(yi):.2f}' "
                f"r='{radius}' fill='{colour}' fill-opacity='{opacity}' "
                f"{ring}>{tip}</circle>")
        return self

    def ring(self, x, y, colour, radius=6.0, width=1.7, title=None):
        """An open circle, for marking a point without hiding its value."""
        for xi, yi in zip(np.atleast_1d(x), np.atleast_1d(y)):
            if not (np.isfinite(xi) and np.isfinite(yi)):
                continue
            tip = f"<title>{esc(title)}</title>" if title else ""
            self.parts.append(
                f"<circle cx='{self._fx(xi):.2f}' cy='{self._fy(yi):.2f}' "
                f"r='{radius}' fill='none' stroke='{colour}' "
                f"stroke-width='{width}'>{tip}</circle>")
        return self

    def line(self, x, y, colour, width=1.8, dash=None, opacity=1.0):
        pairs = [(self._fx(xi), self._fy(yi)) for xi, yi in zip(x, y)
                 if np.isfinite(xi) and np.isfinite(yi)]
        if len(pairs) < 2:
            return self
        d = " ".join(("M" if i == 0 else "L") + f"{px:.2f},{py:.2f}"
                     for i, (px, py) in enumerate(pairs))
        style = f" stroke-dasharray='{dash}'" if dash else ""
        self.parts.append(
            f"<path d='{d}' fill='none' stroke='{colour}' "
            f"stroke-width='{width}' stroke-opacity='{opacity}'{style}/>")
        return self

    def band(self, x, low, high, colour, opacity=0.16):
        top = [(self._fx(xi), self._fy(yi)) for xi, yi in zip(x, high)
               if np.isfinite(xi) and np.isfinite(yi)]
        bot = [(self._fx(xi), self._fy(yi)) for xi, yi in zip(x, low)
               if np.isfinite(xi) and np.isfinite(yi)]
        if len(top) < 2 or len(bot) < 2:
            return self
        d = ("M" + " L".join(f"{px:.2f},{py:.2f}" for px, py in top)
             + " L" + " L".join(f"{px:.2f},{py:.2f}" for px, py in reversed(bot))
             + " Z")
        self.parts.append(f"<path d='{d}' fill='{colour}' "
                          f"fill-opacity='{opacity}' stroke='none'/>")
        return self

    def errorbar(self, x, y, low, high, colour, cap=4):
        px, ylo, yhi = self._fx(x), self._fy(low), self._fy(high)
        self.parts.append(
            f"<path d='M{px:.2f},{ylo:.2f} L{px:.2f},{yhi:.2f} "
            f"M{px-cap:.2f},{ylo:.2f} L{px+cap:.2f},{ylo:.2f} "
            f"M{px-cap:.2f},{yhi:.2f} L{px+cap:.2f},{yhi:.2f}' "
            f"fill='none' stroke='{colour}' stroke-width='1.6'/>")
        self.points([x], [y], colour, radius=4)
        return self

    def hline(self, y, colour=GRID, dash="4 3", width=1.2):
        py = self._fy(y)
        self.parts.append(
            f"<path d='M{self.left},{py:.2f} L{self.width - self.right},{py:.2f}' "
            f"stroke='{colour}' stroke-width='{width}' stroke-dasharray='{dash}' "
            f"fill='none'/>")
        return self

    def label(self, x, y, text, colour=INK, size=11, anchor="start",
              weight="normal", dx=0, dy=0):
        self.overlay.append(
            f"<text x='{self._fx(x) + dx:.2f}' y='{self._fy(y) + dy:.2f}' "
            f"font-size='{size}' fill='{colour}' text-anchor='{anchor}' "
            f"font-weight='{weight}'>{esc(text)}</text>")
        return self

    def note(self, px, py, text, colour=MUTED, size=10.5, anchor="start",
             weight="normal"):
        """Place text in PIXEL coordinates, for annotations outside the data."""
        self.overlay.append(
            f"<text x='{px:.2f}' y='{py:.2f}' font-size='{size}' "
            f"fill='{colour}' text-anchor='{anchor}' "
            f"font-weight='{weight}'>{esc(text)}</text>")
        return self

    # --- frame -----------------------------------------------------------
    def _ticks(self, lo, hi, log):
        if log:
            ticks = []
            for power in range(int(np.floor(lo)), int(np.ceil(hi)) + 1):
                if lo - 1e-9 <= power <= hi + 1e-9:
                    ticks.append((power, _power_label(power)))
            if len(ticks) < 2:
                ticks = [(lo, _power_label(lo)), (hi, _power_label(hi))]
            return ticks
        span = hi - lo
        if span <= 0:
            return [(lo, f"{lo:g}")]
        step = 10 ** np.floor(np.log10(span / 4.0))
        for multiple in (1, 2, 2.5, 5, 10):
            if span / (step * multiple) <= 6:
                step *= multiple
                break
        start = np.ceil(lo / step) * step
        values = np.arange(start, hi + step * 0.5, step)
        # arange's endpoint is computed in floating point and overshoots, which
        # puts a gridline and a tick label outside the frame.
        values = values[(values <= hi + 1e-9) & (values >= lo - 1e-9)]
        # `f"{-1e-16:g}"` renders as "-0", which reads as a distinct value.
        return [(v, f"{0.0 if abs(v) < step * 1e-6 else v:g}") for v in values]

    def render(self, xlabel="", ylabel="", title="", xticks=True,
               yticks=True):
        """`xticks=False` for a categorical x-axis, whose labels the caller
        draws itself -- numeric ticks under named categories are noise. Same
        for `yticks` on a horizontal categorical chart, where the y coordinate
        is a row index and printing it invites the reader to read meaning into
        it."""
        out = [f"<svg viewBox='0 0 {self.width} {self.height}' "
               f"width='100%' style='max-width:{self.width}px' "
               f"font-family=\"system-ui,-apple-system,'Segoe UI',sans-serif\">"]
        x0, x1 = self.left, self.width - self.right
        y0, y1 = self.height - self.bottom, self.top
        for value, text in (self._ticks(*self.ylim, self.ylog)
                            if yticks else []):
            py = self._fy(10 ** value if self.ylog else value)
            out.append(f"<path d='M{x0},{py:.2f} L{x1},{py:.2f}' stroke='{GRID}' "
                       f"stroke-width='1' fill='none'/>")
            out.append(f"<text x='{x0 - 7}' y='{py + 3.5:.2f}' font-size='10.5' "
                       f"fill='{MUTED}' text-anchor='end'>{esc(text)}</text>")
        for value, text in (self._ticks(*self.xlim, self.xlog)
                            if xticks else []):
            px = self._fx(10 ** value if self.xlog else value)
            out.append(f"<path d='M{px:.2f},{y0} L{px:.2f},{y1}' stroke='{GRID}' "
                       f"stroke-width='1' fill='none' stroke-opacity='0.55'/>")
            out.append(f"<text x='{px:.2f}' y='{y0 + 15}' font-size='10.5' "
                       f"fill='{MUTED}' text-anchor='middle'>{esc(text)}</text>")
        out.append(f"<path d='M{x0},{y1} L{x0},{y0} L{x1},{y0}' fill='none' "
                   f"stroke='{INK}' stroke-width='1.3'/>")
        # Marks are clipped to the plot area. A fit whose extrapolation leaves
        # the frame -- the burst fan on the curves where v0 is unbounded runs
        # to 28x the line rate -- should visibly run off the top, not be drawn
        # across the axis labels.
        # THE ID NAMES THE REGION, NOT THE OBJECT. It was `id(self)` until
        # 2026-09-04 -- a memory address, which CPython reuses the moment a
        # panel's Axes is freed, so a page of many figures drew far fewer
        # distinct ids than it had clips. `two_axis/progress_curves.html` had
        # 119 clipPaths sharing 5 ids. Ids are DOCUMENT-scoped in HTML, so
        # `url(#c)` resolves to the FIRST clipPath of that name on the page,
        # not to the one beside the marks -- and on five of the six index
        # pages the collisions bound one id to conflicting rectangles, which
        # clips a figure's marks to a DIFFERENT figure's plot area and deletes
        # whatever falls outside it. Silently: the mark is in the file.
        # Hashing the geometry fixes both halves. Two clips of the same region
        # may still share an id, which is what we want -- the referent is
        # identical, so resolving to the first is correct -- and two clips of
        # different regions can no longer collide. It is also deterministic,
        # so a rebuild that changes nothing now produces a byte-identical
        # page; with an address in the markup no such check could ever pass.
        rect = (f"<rect x='{x0}' y='{y1}' "
                f"width='{x1 - x0}' height='{y0 - y1}'/>")
        clip = "clip" + hashlib.blake2s(rect.encode(), digest_size=6).hexdigest()
        out.append(f"<defs><clipPath id='{clip}'>{rect}</clipPath></defs>")
        out.append(f"<g clip-path='url(#{clip})'>")
        out.extend(self.parts)
        out.append("</g>")
        # Annotations sit OUTSIDE the clip. Axis category labels live below the
        # frame and legends above it, and clipping them deleted both when this
        # clip was introduced -- figure C lost its estimator names entirely.
        out.extend(self.overlay)
        if xlabel:
            out.append(f"<text x='{(x0 + x1) / 2:.1f}' y='{self.height - 3}' "
                       f"font-size='11.5' fill='{INK}' text-anchor='middle'>"
                       f"{esc(xlabel)}</text>")
        if ylabel:
            out.append(f"<text transform='translate(13,{(y0 + y1) / 2:.1f}) "
                       f"rotate(-90)' font-size='11.5' fill='{INK}' "
                       f"text-anchor='middle'>{esc(ylabel)}</text>")
        if title:
            out.append(f"<text x='{x0}' y='{y1 - 3}' font-size='12' "
                       f"fill='{INK}' font-weight='600'>{esc(title)}</text>")
        out.append("</svg>")
        return "".join(out)


def _power_label(power):
    power = int(round(power))
    return f"10<tspan dy='-4' font-size='8'>{power}</tspan>" if False else f"1e{power}"


PAGE_CSS = """
:root{--ink:#1a1a1a;--muted:#6b6b6b;--rule:#e2e2e2;--bg:#ffffff;--panel:#fafafa;
--accent:#c0522a;--good:#3f8a5a}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
--ink:#e8e8e8;--muted:#a0a0a0;--rule:#333;--bg:#141414;--panel:#1c1c1c}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.6 system-ui,-apple-system,'Segoe UI',sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:32px 22px 72px}
h1{font-size:23px;margin:0 0 6px;letter-spacing:-0.01em}
h2{font-size:16px;margin:38px 0 10px;padding-top:16px;border-top:1px solid var(--rule)}
h3{font-size:13.5px;margin:22px 0 6px;color:var(--muted);font-weight:600;
text-transform:uppercase;letter-spacing:0.06em}
p{margin:9px 0;max-width:78ch}
.lede{color:var(--muted);font-size:14.5px;max-width:82ch}
.fig{background:var(--panel);border:1px solid var(--rule);border-radius:7px;
padding:15px 15px 9px;margin:16px 0}
.cap{font-size:12.2px;color:var(--muted);margin:8px 2px 2px;max-width:88ch}
.grid{display:grid;gap:14px}
@media(min-width:900px){.grid.two{grid-template-columns:1fr 1fr}
.grid.three{grid-template-columns:repeat(3,1fr)}}
table{border-collapse:collapse;font-size:12.7px;margin:12px 0;width:100%}
th,td{padding:5px 10px;border-bottom:1px solid var(--rule);text-align:right}
th:first-child,td:first-child{text-align:left}
th{font-weight:600;color:var(--muted);font-size:11.5px;text-transform:uppercase;
letter-spacing:0.05em}
tr.hl td{background:rgba(192,82,42,0.09);font-weight:600}
code{font:12.3px ui-monospace,SFMono-Regular,Menlo,monospace;
background:var(--rule);padding:1px 5px;border-radius:3px}
.key{display:flex;flex-wrap:wrap;gap:7px 20px;font-size:12.2px;
color:var(--muted);margin:6px 2px 10px}
.key span{display:flex;align-items:center;gap:6px}
.sw{width:15px;height:3px;border-radius:2px;display:inline-block}
.scroll{overflow-x:auto}
.warn{border-left:3px solid var(--accent);padding:2px 0 2px 13px;
color:var(--muted);font-size:12.8px;margin:12px 0;max-width:82ch}
.panel{padding:12px 12px 8px}
.ph{font-size:12.5px;font-weight:650;letter-spacing:-0.005em}
.ps{font-size:10.9px;color:var(--muted);margin:1px 0 6px}
.pf{font-size:10.6px;color:var(--muted);margin-top:5px;
border-top:1px solid var(--rule);padding-top:5px}
table.nums{margin:7px 0 0;font-size:11.2px;width:100%}
table.nums td{padding:2.5px 5px;border-bottom:none;text-align:left;
white-space:nowrap}
table.nums td.num{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
text-align:right;font-size:11px}
table.nums td.dim{color:var(--muted);font-size:10.4px}
table.nums tr:nth-child(odd){background:rgba(128,128,128,0.055)}
.sw{width:13px;height:3px;border-radius:2px;display:inline-block;
margin-right:6px;vertical-align:middle}
"""


def page(title, body, subtitle=""):
    return (f"<!doctype html><html><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{esc(title)}</title><style>{PAGE_CSS}</style></head><body>"
            f"<div class='wrap'><h1>{esc(title)}</h1>"
            + (f"<p class='lede'>{subtitle}</p>" if subtitle else "")
            + body + "</div></body></html>")


def colliding_clips(html_text):
    """
    Every clip id on the page that names more than one rectangle.

    SVG ids are DOCUMENT-scoped, so `clip-path='url(#c)'` resolves to the FIRST
    clipPath called `c` anywhere on the page -- not to the one sitting three
    lines above the marks that reference it. Bind one id to two different
    rectangles and a figure's marks are clipped to a DIFFERENT figure's plot
    area, and whatever falls outside that rectangle is gone from the rendered
    page while remaining present in the file.

    That happened, on five of the six index pages, for as long as the pages
    existed. The id was `id(self)` -- a CPython memory address, reused the
    moment a panel's Axes is freed -- so `two_axis/progress_curves.html` drew
    119 clipPaths under 5 distinct ids. `clipped_marks` could not see it: it
    reads each figure's OWN first clipPath, which is the one the browser
    ignores.

    Two clips of the same region sharing an id is fine and is not reported --
    the referent is identical, so resolving to the first is correct.

    Returns a list of (clip id, [rects]) for the conflicts, empty when clean.
    Pass the rendered HTML of a whole page.
    """
    import collections
    import re
    regions = collections.defaultdict(list)
    for cid, rect in re.findall(r"<clipPath id='([^']+)'><rect ([^/]*)/>",
                                html_text):
        rect = rect.strip()
        if rect not in regions[cid]:
            regions[cid].append(rect)
    return [(cid, rects) for cid, rects in regions.items() if len(rects) > 1]


def clipped_marks(html_text):
    """
    Every drawn point that falls outside its own figure's plot frame.

    Marks are clipped to the plot area on purpose -- a fit whose extrapolation
    leaves the frame should visibly run off rather than be drawn over the axis
    labels. The cost is that a DATA point outside the axis limits disappears in
    exactly the same way, silently, and the figure looks complete.

    That happened: the curvature figure's y-axis floor was -90 while exp 6
    sample 4 sits at -91.66, so a live curve was missing from a published page
    for as long as the page existed. Nothing caught it because every check in
    this project read the prose, and a missing point changes no number.

    Returns a list of (figure title, x, y, frame) tuples, empty when clean.
    Pass the rendered HTML of a whole page.
    """
    import re
    found = []
    for svg in re.findall(r"<svg .*?</svg>", html_text, re.S):
        frame = re.search(r"<clipPath id='[^']*'><rect x='([\d.]+)' "
                          r"y='([\d.]+)' width='([\d.]+)' height='([\d.]+)'",
                          svg)
        if not frame or "clip-path=" not in svg:
            continue
        x0, y0, width, height = (float(v) for v in frame.groups())
        title = re.search(r"font-weight='600'>([^<]*)</text>", svg)
        body = svg.split("clip-path=", 1)[1]
        for cx, cy in re.findall(r"<circle cx='(-?[\d.]+)' cy='(-?[\d.]+)'",
                                 body):
            cx, cy = float(cx), float(cy)
            # Half a pixel of slack: a point exactly on the axis is drawn on
            # the frame line and is not lost.
            if not (x0 - 0.5 <= cx <= x0 + width + 0.5
                    and y0 - 0.5 <= cy <= y0 + height + 0.5):
                found.append((title.group(1) if title else "(untitled)",
                              cx, cy, (x0, y0, width, height)))
    return found
