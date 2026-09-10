"""
Eyeball page for the 2026-09-11 segmentation rewrite: every curve it moved.

    python scratch/build_detection_review.py            # against HEAD~1
    python scratch/build_detection_review.py <ref>      # against any ref

WHAT THIS IS FOR. `data/bubble_cases.py` says whether the layer agrees with
the twenty-five curves somebody argued about, and `run_gates.py` says whether
anything published moved. Neither shows what the detector now DRAWS, and the
three faults the old layer had were all found by eye rather than by a test.
This is the page for looking.

It renders BOTH versions of `curve_metrics` on the same axes -- the shipped
one from `data/`, and the pre-rewrite one loaded straight out of a git
worktree -- so every panel is a before/after on the same readings rather than
an assertion that something improved. Nothing here is a gate and nothing
imports it; it is a scratch page, and the numbers in it come from the modules
like everywhere else.

Amber bands are the beam SHEDDING gas (a detachment's own span), violet bands
are the beam TAKING IT ON (an arrival's). Both are `figure_kit`'s, the same
device `two_axis/progress_curves.html` uses, so a band means the same thing
here as it does there.
"""
import importlib.util
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPOSITORY = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPOSITORY, "data"))
sys.path.insert(0, REPOSITORY)

import curve_metrics
import scope
from figure_kit import (ARRIVAL_BAND_COLOUR, EVENT_BAND_COLOUR, MUTED,
                        panel, progress_axes, styled, write_pages)

DOCUMENT = "detection_review.html"


def _previous(ref):
    """
    The pre-rewrite `curve_metrics`, loaded from a throwaway worktree.

    A COPY OF THE OLD CODE WOULD ROT. The whole point of the comparison is
    that the left-hand column is what the package actually did, so it is read
    out of git rather than pasted here and left to drift out of agreement with
    the history it claims to show. The module imports nothing but numpy, so it
    loads standalone.
    """
    directory = tempfile.mkdtemp(prefix="detection-review-")
    subprocess.run(["git", "worktree", "add", "--detach", directory, ref],
                   cwd=REPOSITORY, check=True, capture_output=True)
    path = os.path.join(directory, "data", "curve_metrics.py")
    spec = importlib.util.spec_from_file_location("curve_metrics_before", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module._worktree = directory
    return module


def _discard(module):
    subprocess.run(["git", "worktree", "remove", "--force", module._worktree],
                   cwd=REPOSITORY, check=False, capture_output=True)


def _read(module, curve):
    """One version's verdict on one curve: events and arrivals."""
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    events = module.detachments(values, curve.noise)
    arrivals = module.bubble_arrivals(times, values, curve.noise, events)
    return events, arrivals


def _draw(curve, events, arrivals, window=None):
    """
    One version's reading of one curve, banded.

    `window` is a (low, high) span in READINGS. A 480-reading run drawn whole
    puts a two-reading detachment inside one pixel, which is exactly the
    detail this page exists to show, so a panel about a specific stretch is
    drawn over that stretch and says so in its own caption.
    """
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    if window is not None:
        low, high = window
        low, high = max(0, low), min(len(times), high)
        times, values = times[low:high], values[low:high]
        shift = low
    else:
        shift = 0
    axes, _ = progress_axes(times, values, limit=140, width=330, height=170,
                            pad=(56, 12, 30, 20))

    def band(lo_index, hi_index, colour):
        # OVERLAP FIRST, THEN CLAMP, and in that order for a reason: clamping
        # an event that lies wholly outside the window lands both its edges on
        # the same boundary reading, which draws a zero-width band at the frame
        # rather than nothing at all. On a curve with seventeen events and a
        # twenty-reading window that is fifteen smears down the edges, and it
        # reads as detections the panel is not making.
        if hi_index < shift or lo_index > shift + len(times) - 1:
            return
        lo = float(times[max(0, lo_index - shift)])
        hi = float(times[min(len(times) - 1, hi_index - shift)])
        axes.band([lo, max(hi, lo + 1e-9)],
                  [axes.ylim[0], axes.ylim[0]], [axes.ylim[1], axes.ylim[1]],
                  colour, opacity=0.16)

    for start, stop in events:
        band(start, stop, EVENT_BAND_COLOUR)
    for index, _gain in arrivals:
        band(index - 1, index, ARRIVAL_BAND_COLOUR)
    axes.line(times, values, MUTED, width=0.8, opacity=0.55)
    return axes.render("time, s", "absorbance")


def _pair(curve, before, after, window, heading, note):
    """A before/after row for one stretch of one curve."""
    old_events, old_arrivals = _read(before, curve)
    new_events, new_arrivals = _read(after, curve)
    inside = lambda pairs: [p for p in pairs
                            if window[0] <= p[0] <= window[1]]
    left = panel(
        f"before &mdash; {len(old_events)} detachments, "
        f"{len(old_arrivals)} arrivals",
        f"{len(inside(old_events))} and "
        f"{len([a for a in old_arrivals if window[0] <= a[0] <= window[1]])} "
        f"in this window",
        _draw(curve, old_events, old_arrivals, window))
    right = panel(
        f"after &mdash; {len(new_events)} detachments, "
        f"{len(new_arrivals)} arrivals",
        f"{len(inside(new_events))} and "
        f"{len([a for a in new_arrivals if window[0] <= a[0] <= window[1]])} "
        f"in this window",
        _draw(curve, new_events, new_arrivals, window))
    return (f"<h3>{heading}</h3>"
            f"<div class='grid two'>{left}{right}</div>"
            f"<div class='cap'>{note}</div>")


def _curves():
    return {(c.experiment, c.sample): c for c in scope.curves(scope.archive())}


def _classify(before, after, lookup):
    """
    Every curve the rewrite moved, sorted into what moved on it.

    Returns four lists of `(key, detail)`. A curve can appear in more than
    one: exp 135.1 both merges releases and gains a full-length arrival.
    """
    merged, credited, dropped, unchanged = [], [], [], []
    for key, curve in sorted(lookup.items()):
        old_events, old_arrivals = _read(before, curve)
        new_events, new_arrivals = _read(after, curve)
        if not (old_events or old_arrivals or new_events or new_arrivals):
            continue
        old_falls = {i for a, b in old_events for i in range(a, b)}
        new_falls = {i for a, b in new_events for i in range(a, b)}
        old_landings = {i for i, _ in old_arrivals}
        new_landings = {i for i, _ in new_arrivals}
        if (len(new_events) < len(old_events)) and new_falls >= old_falls:
            merged.append((key, (len(old_events), len(new_events))))
        gained = sorted(new_landings - old_landings)
        if gained:
            credited.append((key, gained))
        lost = sorted(old_landings - new_landings)
        if lost and not gained:
            dropped.append((key, lost))
        if (old_events == new_events and old_landings == new_landings):
            unchanged.append((key, len(new_events)))
    return merged, credited, dropped, unchanged


def main(ref="HEAD~1"):
    before = _previous(ref)
    try:
        after = curve_metrics
        lookup = _curves()
        merged, credited, dropped, unchanged = _classify(before, after, lookup)
        body = []

        body.append(
            "<p class='lede'>Every curve the 2026-09-11 segmentation rewrite "
            "moved, drawn both ways on the same readings. The left panel is "
            f"<code>{ref}</code>, the right is the working tree. Amber is the "
            "beam shedding gas, violet is the beam taking it on.</p>")
        body.append(
            "<div class='key'>"
            f"<span><i class='sw' style='background:{EVENT_BAND_COLOUR}'></i>"
            "detachment span</span>"
            f"<span><i class='sw' style='background:{ARRIVAL_BAND_COLOUR}'></i>"
            "arrival span</span></div>")

        body.append("<h2>1. Stuttering releases, now one event</h2>")
        body.append(
            "<p>A bubble sliding out of the beam can pause or tick up for a "
            "single reading. Strict-consecutive grouping cut those releases "
            "into fragments and read the tick between two of them as gas "
            "<em>arriving</em>. Same falls either way &mdash; fewer, longer "
            "events, and the false arrival gone.</p>")
        for key, window, heading, note in MERGE_CASES:
            body.append(_pair(lookup[key], before, after, window, heading, note))

        body.append("<h2>2. Multi-reading arrivals, credited in full</h2>")
        body.append(
            "<p>One bubble growing over several readings. Scoring every step "
            "alone credited its largest single step and left the rest for the "
            "production rate to invent.</p>")
        for key, window, heading, note in CREDIT_CASES:
            body.append(_pair(lookup[key], before, after, window, heading, note))

        body.append("<h2>3. Rises that were never anomalous</h2>")
        body.append(
            "<p>Forty-four of the old eighty-nine arrivals. Scored on absolute "
            "&sigma; they ran 6.5 to 160&sigma; and looked convincing; scored "
            "against their own neighbourhoods they run 1.00 to 1.91, meaning "
            "the curve was already climbing that fast. These are the panels "
            "to disagree with if any of it is wrong.</p>")
        for key, window, heading, note in ORDINARY_CASES:
            body.append(_pair(lookup[key], before, after, window, heading, note))

        body.append("<h2>4. The pinned cases, unmoved</h2>")
        body.append(
            "<p>The curves the old thresholds were pinned to, which a rewrite "
            "is most likely to break. <code>data/bubble_cases.py</code> "
            "asserts all of these; here they are drawn.</p>")
        for key, window, heading, note in PINNED_CASES:
            body.append(_pair(lookup[key], before, after, window, heading, note))

        body.append("<h2>5. What moved, archive-wide</h2>")
        body.append(_summary(before, after, lookup, merged, credited,
                             dropped, unchanged))

        page = styled("Detection Layer Review", "".join(body),
                      "segmentation against the pointwise tests it replaced")
        return write_pages(HERE, {DOCUMENT: page})
    finally:
        _discard(before)


def _summary(before, after, lookup, merged, credited, dropped, unchanged):
    """The archive-wide tally, recomputed rather than quoted."""
    old_e = old_a = new_e = new_a = 0
    old_falls = new_falls = new_span = 0
    for curve in lookup.values():
        oe, oa = _read(before, curve)
        ne, na = _read(after, curve)
        old_e += len(oe); old_a += len(oa)
        new_e += len(ne); new_a += len(na)
        old_falls += sum(b - a for a, b in oe)
        steps = np.diff(np.asarray(curve.absorbance, dtype=float))
        new_span += sum(b - a for a, b in ne)
        new_falls += sum(1 for a, b in ne for i in range(a, b) if steps[i] < 0)
    rows = [
        ("detachment events", old_e, new_e),
        ("falling readings inside them", old_falls, new_falls),
        ("readings those events span", old_falls, new_span),
        ("arrivals", old_a, new_a),
    ]
    table = "".join(f"<tr><td>{name}</td><td>{was}</td>"
                    f"<td><strong>{now}</strong></td></tr>"
                    for name, was, now in rows)
    counts = (f"<p>{len(merged)} curves merge a stuttering release, "
              f"{len(credited)} gain a fully-credited arrival, "
              f"{len(dropped)} lose an arrival that was never anomalous, and "
              f"{len(unchanged)} are untouched.</p>")
    return (counts + "<div class='scroll'><table>"
            "<tr><th>quantity</th><th>before</th><th>after</th></tr>"
            + table + "</table></div>")


# The stretches worth looking at, one per finding. Windows are in READINGS and
# are deliberately tight: the point of each panel is a two-to-four reading
# feature that a whole-run view puts inside one pixel.
MERGE_CASES = (
    ((131, 1), (205, 230),
     "Exp 131 cuvette 1, readings 215&ndash;218",
     "Two events become one over exactly the same falls, and reading 216 "
     "&mdash; the stutter itself &mdash; is gained rather than lost. This "
     "curve's count goes 18 &rarr; 17; its <code>bubble_cases</code> row "
     "moved with it, and the row's claim (a heavy bubbler the SNR floor must "
     "leave alone) did not."),
    ((131, 2), (58, 80),
     "Exp 131 cuvette 2, readings 66&ndash;69",
     "The same shape again, twice on this curve (readings 66&ndash;69 and "
     "143&ndash;146), 19 &rarr; 17 events."),
    ((43, 1), (48, 66),
     "Exp 43 cuvette 1, readings 54&ndash;57",
     "The fragmentation artefact at its clearest: readings 55&ndash;58 are "
     "one stuttering release, and the +9.6&sigma; uptick between two "
     "fragments of it scored as an arrival of 0.0033 AU. Now swallowed by "
     "the release it sits inside. The curve's total stays at seven events "
     "because these two merging (&minus;1) is offset by a real fall at "
     "reading 114 that the old excursion test had rejected (+1)."),
    ((135, 2), (150, 168),
     "Exp 135 cuvette 2, readings 153&ndash;158",
     "9180&ndash;9540 s is one messy release cut into three events, with the "
     "uptick between two of them credited as gas arriving. Three amber bands "
     "become one and the violet band inside them goes. The curve's own "
     "arrival count is unchanged at five only because a different, genuine "
     "one at reading 147 is now credited over three readings instead of one "
     "&mdash; that is section 2's business, and the window here starts after "
     "it."),
)

CREDIT_CASES = (
    ((44, 1), (80, 100),
     "Exp 44 cuvette 1, readings 87&ndash;90",
     "One bubble growing over three readings &mdash; +69, +96 and "
     "+75&sigma;, totalling 0.0667 AU against the 0.0663 the next detachment "
     "sheds. Credited 0.0231 AU at reading 89 before; 0.0569 at reading 90 "
     "now. The near-perfect cancellation with the fall that follows is the "
     "bubble lifecycle, which is why <code>_excursions</code> may not reject "
     "an <code>accumulate &rarr; release</code> pair on size alone."),
    ((138, 4), (385, 405),
     "Exp 138 cuvette 4, readings 392&ndash;398",
     "The pair the old code got wrong in both halves: it rejected the real "
     "2-reading detachment because a rise follows, and would then have "
     "refused the real 4-reading arrival for following a rejected fall. The "
     "rise does not <em>return</em> &mdash; it overshoots the pre-fall level "
     "by 0.0176 AU and holds, so the pair does not cancel and both are real."),
    ((135, 1), (265, 285),
     "Exp 135 cuvette 1, readings 273&ndash;277",
     "One arrival that a run rule tolerating no interruption split in two: "
     "+0.0031, +0.0192, a &minus;0.0029 wobble, +0.0117, then flat. The "
     "wobble is 6.6&sigma; and was scored as a detachment of its own."),
    ((49, 1), (30, 50),
     "Exp 49 cuvette 1, readings 38&ndash;39",
     "Credited 0.0210 AU before and 0.0714 now &mdash; the largest single "
     "step against the whole of what the beam took on."),
)

ORDINARY_CASES = (
    ((13, 4), (55, 80),
     "Exp 13 cuvette 4, reading 64 &mdash; 43.8&sigma;, and ordinary",
     "The clearest of the forty-four. The step is 18.4 thousandths; the six "
     "around it are 11.8, 11.8, 10.3, 10.9, 14.6 and 12.7. It scores 1.57 "
     "against its own neighbourhood and had been credited 0.0067 AU of gas "
     "for rising 1.6&times; an ordinary step."),
    ((45, 3), (40, 62),
     "Exp 45 cuvette 3, reading 50 &mdash; 61.9&sigma;, scoring 1.91",
     "The closest of the forty-four to the bar, and still under it. The curve "
     "is climbing at 10&ndash;15 thousandths a reading throughout."),
    ((144, 2), (25, 48),
     "Exp 144 cuvette 2, readings 29&ndash;42 &mdash; the case that proves it",
     "Fourteen consecutive readings of real, smooth 7&ndash;12&sigma; "
     "acceleration, every one scoring about 1.0 locally. The old rule kept "
     "this out only by refusing to merge <em>any</em> rise ever, which "
     "protected this curve and left the same mistake standing on 27 others. "
     "Scored locally it needs no special rule, and the real detachment at "
     "readings 43&ndash;44 is untouched."),
    ((142, 4), (50, 72),
     "Exp 142 cuvette 4, readings 58 and 67",
     "Two arrivals dropped on a heavy bubbler whose detachments are all "
     "unchanged &mdash; the rises scored 1.72 and 1.41 against a curve "
     "already stepping 4&ndash;5 thousandths."),
)

PINNED_CASES = (
    ((149, 5), (0, 20),
     "Exp 149 cuvette 5, readings 8&ndash;9 &mdash; still rejected",
     "The curve that forced the recovery clause. It falls 0.00206 and the "
     "next reading climbs 0.00222 straight back: a release phase and the "
     "accumulate phase beside it, cancelling. All four of its falls are "
     "rejected and the curve is returned untouched, as before."),
    ((135, 1), (215, 235),
     "Exp 135 cuvette 1, readings 222&ndash;223 &mdash; still kept",
     "41.3&sigma;, sitting right after four readings of fast, real "
     "acceleration. Extending a recovery test backwards flagged this one; "
     "segmentation reaches it the other way, because an ordinary reading "
     "separates the rise from the fall so the two never pair."),
    ((135, 1), (262, 280),
     "Exp 135 cuvette 1, readings 272&ndash;273 &mdash; kept by margin now",
     "6.2&sigma;, and the curve climbs 0.0312 AU over the four readings after "
     "a fall that cost 0.0027. The old test kept it on a 2% margin "
     "(3.155e&minus;3 against a bar of 3.210e&minus;3); a two-sided "
     "cancellation test keeps it by a factor of ten, because a counter-move "
     "that overshoots does not cancel."),
    ((130, 2), (46, 68),
     "Exp 130 cuvette 2, readings 54&ndash;55 and 60&ndash;61 &mdash; kept",
     "12.4&sigma; and 16.6&sigma; on a curve climbing 2.5 thousandths a "
     "reading. Their neighbours score 1.0&times; and 1.4&times; locally, so "
     "they are not phases at all and there is nothing for either fall to "
     "pair off against."),
)


if __name__ == "__main__":
    raise SystemExit(main(*sys.argv[1:]))
