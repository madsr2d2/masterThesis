"""
The bubble machinery's accumulated calibration, as one table.

WHY THIS EXISTS. `detachments`, `bubble_arrivals` and the constants they read
are not derived from anything -- they are the residue of about a dozen
specific curves that were looked at by eye, argued over, and then written into
a threshold. That knowledge is the real asset in this corner of the package,
and until this file existed it was scattered across four docstrings, three
constant comments, several test bodies and a run of DATA_VERIFICATION
entries. Nothing could be restructured with confidence because nobody could
say, in one place, what the current behaviour is *for*.

AND THEN IT WAS. On 2026-09-11 the detection layer was rewritten from
pointwise tests to segmentation (`bubble_segments`, `_excursions`), and this
table is what made that judgeable instead of arguable: all seven OPEN rows
flipped to agreeing, every PINNED row held but two, and the two that moved
moved for a reason the table itself could state. Four `source` fields below
name functions and constants that no longer exist; they are left pointing at
the DATA_VERIFICATION entry that carries the argument instead.

Each row names a curve, an event, the verdict the evidence supports, and
where that evidence is written down. Two statuses:

  PINNED  the code agrees, and a change that moves it is a regression --
          these are the cases the thresholds were pinned to.
  OPEN    the code DISAGREES and is believed wrong. Found by eye, argued in
          DATA_VERIFICATION.md, not yet fixed. `test_bubble_cases` asserts
          these still fail, so whoever fixes one is told to promote it.

THIS IS A FIXTURE, NOT A MEASUREMENT. Nothing here is computed; every field
is a transcription of something already written down, and the `source` field
says where. If a row and its source disagree, the source wins.
"""

# kind        what the row is about
#   "fall"    a candidate detachment at `where` = (start, stop)
#   "arrival" a confirmed rise landing at `where` = index
#   "curve"   a property of the whole curve; `where` is None
#
# verdict     what the evidence supports
#   fall      "kept" | "rejected"
#   arrival   "released" | "unreleased" | "none"
#   curve     "no events" | "untouched" | ("kept", n)
#
# `gain` is optional and only meaningful on an arrival row. Where it is
# present the CREDITED AMOUNT is the point, not merely that something was
# detected there: a multi-reading arrival credited with its largest single
# step is present at roughly the right reading and still wrong.
BUBBLE_CASES = (
    # ---- exp 149 cuvette 5: the curve that forced the recovery clause -----
    # An unfiltered repair set a production rate off these four and then
    # removed 0.0097 AU from a curve that rose 0.0262, flattening a real
    # early rise into a straight line.
    dict(experiment=149, sample=5, kind="fall", where=(8, 9),
         verdict="rejected", status="pinned",
         why="falls 0.00206 and the next reading climbs 0.00222 straight "
             "back -- 7.2x the local baseline step. A release phase and the "
             "accumulate phase beside it, cancelling",
         source="curve_metrics._excursions"),
    dict(experiment=149, sample=5, kind="fall", where=(56, 57),
         verdict="rejected", status="pinned",
         why="departs from a reading that is itself an isolated spike, 2.1x "
             "the local baseline -- one interval up and straight back down, "
             "which is the only shape accumulate -> release may be rejected on",
         source="curve_metrics._excursions"),
    dict(experiment=149, sample=5, kind="fall", where=(36, 37),
         verdict="rejected", status="pinned",
         why="a candidate only since the 2026-09-07 threshold change; its "
             "recovery lands one or two readings late",
         source="DATA_VERIFICATION.md 2026-09-07"),
    dict(experiment=149, sample=5, kind="fall", where=(148, 149),
         verdict="rejected", status="pinned",
         why="the same delayed-recovery shape as (36, 37)",
         source="DATA_VERIFICATION.md 2026-09-07"),

    # ---- exp 135 cuvette 1: two falls that must NOT be read as spikes ----
    dict(experiment=135, sample=1, kind="fall", where=(272, 273),
         verdict="kept", status="pinned",
         why="6.2 sigma, and the curve accelerates hard in the two readings "
             "after it -- 0.0312 AU over four, against the 0.0027 the fall "
             "cost. A recovery test read that acceleration as the fall's own "
             "reversal and survived only on a 2% margin; a two-sided "
             "cancellation test keeps it by a factor of ten, because a "
             "counter-move that overshoots does not cancel",
         source="curve_metrics._excursions"),
    dict(experiment=135, sample=1, kind="fall", where=(222, 223),
         verdict="kept", status="pinned",
         why="41.3 sigma, sitting right after four readings of fast, real "
             "acceleration. Extending the recovery reach to the reading "
             "BEFORE a fall as well flagged this one. Segmentation reaches "
             "it the other way: the rise is a phase of its own, separated "
             "from the fall by an ordinary reading, so the two never pair",
         source="curve_metrics._excursions"),

    # ---- exp 130 cuvette 2: real gas on a fast-rising curve --------------
    dict(experiment=130, sample=2, kind="fall", where=(54, 55),
         verdict="kept", status="pinned",
         why="12.4 sigma. Its neighbouring steps are unremarkable against "
             "the ~0.0025 AU/reading the curve climbs at, but exceeded half "
             "of this fall's own modest size -- the recovery test alone "
             "rejected it until the local step scale was added. Those "
             "neighbours score 1.0x and 1.4x, under ANOMALY_BAR, so they are "
             "not phases at all and there is nothing for the fall to pair off "
             "against",
         source="curve_metrics.local_step_scale"),
    dict(experiment=130, sample=2, kind="fall", where=(60, 61),
         verdict="kept", status="pinned",
         why="16.6 sigma, rejected the same way and for the same reason",
         source="curve_metrics.local_step_scale"),

    # ---- the curve-level SNR gate ----------------------------------------
    dict(experiment=150, sample=1, kind="curve", where=None,
         verdict="no events", status="pinned",
         why="net/noise 20.7, under DETACHMENT_SNR_FLOOR. Local noise in the "
             "stretches its candidates sit in runs up to 5x the curve's "
             "global estimate, and no per-event test could separate its "
             "falls from real gas elsewhere in the block",
         source="curve_metrics.DETACHMENT_SNR_FLOOR"),
    dict(experiment=151, sample=6, kind="curve", where=None,
         verdict="no events", status="pinned",
         why="net/noise 13.9. Reached by the floor now; before that it was "
             "the 2026-09-07 depth extension that took both its candidates",
         source="curve_metrics.DETACHMENT_SNR_FLOOR"),
    dict(experiment=66, sample=3, kind="curve", where=None,
         verdict="no events", status="pinned",
         why="a dead curve whose one candidate is a mixing transient at the "
             "very start of the run",
         source="curve_metrics.DETACHMENT_SNR_FLOOR"),
    dict(experiment=131, sample=1, kind="curve", where=None,
         verdict=("kept", 17), status="pinned",
         why="net/noise 44.6 -- a genuine heavy bubbler the floor must leave "
             "alone. Its bubble_load, 6.5-8.3, is as high as exp 150.1's, so "
             "load alone cannot tell the two apart. WAS 18 UNTIL 2026-09-11, "
             "on the same falls: readings 215-218 are one stuttering release "
             "that strict-consecutive grouping cut in two. No falling reading "
             "was lost -- 216 was gained",
         source="curve_metrics.DETACHMENT_SNR_FLOOR"),
    dict(experiment=131, sample=2, kind="curve", where=None,
         verdict=("kept", 17), status="pinned",
         why="net/noise 36.8, the other side of the gap the floor sits in. "
             "WAS 19 until 2026-09-11 and for the same reason as its sibling: "
             "readings 66-69 and 143-146 are each one stuttering release, and "
             "readings 67 and 144 were gained rather than lost",
         source="curve_metrics.DETACHMENT_SNR_FLOOR"),

    # ---- the one curve no rate can pay for -------------------------------
    dict(experiment=135, sample=6, kind="curve", where=None,
         verdict="untouched", status="pinned",
         why="its fall is in the FIRST reading interval, so no rate explains "
             "a bubble grown before the run began. bubble_rate returns inf "
             "and debubble hands the readings back unchanged",
         source="curve_metrics.bubble_rate"),

    # ---- arrivals: the cases the rise tests were built against ------------
    dict(experiment=144, sample=2, kind="no arrival in", where=(29, 43),
         verdict="none", status="pinned",
         why="readings 29-42 climb 20-30 sigma a step for fourteen "
             "consecutive readings, real and smooth. Merged the way a fall's "
             "candidates are, this scored as one ~0.034 AU jump -- larger "
             "than any real gain in the block. Not one step here may be an "
             "arrival",
         source="curve_metrics.bubble_arrivals"),
    dict(experiment=135, sample=5, kind="arrival", where=163,
         verdict="released", status="pinned",
         why="the jump at 9780 s the kink test was built for: the reading "
             "before scores -8.4 and the one it lands on +10.0",
         source="curve_metrics.bubble_arrivals"),
    dict(experiment=146, sample=4, kind="arrival", where=449,
         verdict="unreleased", status="pinned",
         why="a bubble that arrived near the end and never left -- the only "
             "shape for which apply_gains' permanent shift is the right "
             "operator. The curve carries no detachment at all",
         source="curve_metrics.split_arrivals"),
    dict(experiment=139, sample=2, kind="arrival", where=71,
         verdict="released", status="pinned",
         why="a 17.1 sigma step two readings before a confirmed detachment. "
             "Rejected until 2026-09-10 because the detachment that released "
             "it read, to the recovery test, as a spike reverting",
         source="DATA_VERIFICATION.md 2026-09-10"),

    # ---- PROMOTED 2026-09-11: open until segmentation fixed all seven -----
    # These seven were found by eye and argued out in DATA_VERIFICATION.md
    # while the code disagreed with every one of them. They are what the
    # rewrite was FOR, and they are pinned now because it delivered them --
    # not one at a time, but as a consequence of representing the phase
    # structure that all seven were symptoms of the absence of.
    dict(experiment=44, sample=1, kind="arrival", where=90, gain=0.0569,
         verdict="released", status="pinned",
         why="ONE bubble growing over three readings -- +69, +96 and +75 "
             "sigma, totalling 0.0667 AU, against the 0.0663 the next "
             "detachment sheds. The rise tests never merge, so only the "
             "largest single step is credited: 0.0231, about a third of it",
         source="DATA_VERIFICATION.md 2026-09-10 (second entry)"),
    dict(experiment=43, sample=1, kind="arrival", where=57,
         verdict="none", status="pinned",
         why="not an arrival at all. Readings 55-58 are one stuttering "
             "release; detachments groups only strictly consecutive falls, "
             "so the uptick between two fragments of it scores as gas "
             "arriving",
         source="DATA_VERIFICATION.md 2026-09-10 (second entry)"),
    dict(experiment=135, sample=2, kind="arrival", where=158,
         verdict="none", status="pinned",
         why="the same fragmentation: 9180-9540 s is one messy release cut "
             "into three events, and this uptick sits between two of them",
         source="DATA_VERIFICATION.md 2026-09-10 (second entry)"),
    dict(experiment=135, sample=1, kind="arrival", where=277, gain=0.0257,
         verdict="released", status="pinned",
         why="readings 274-277 are ONE arrival -- +0.0031, +0.0192, a "
             "-0.0029 wobble, +0.0117, then flat. A run rule that tolerates "
             "no interruption splits it in two, and the code credits 0.0103 "
             "here: present at the right reading, and still the wrong event",
         source="DATA_VERIFICATION.md 2026-09-10 (second entry)"),
    dict(experiment=135, sample=1, kind="arrival", where=275,
         verdict="none", status="pinned",
         why="the first half of that same split. Merged into the run ending "
             "at 277 it is not an event of its own",
         source="DATA_VERIFICATION.md 2026-09-10 (second entry)"),
    dict(experiment=138, sample=4, kind="fall", where=(392, 394),
         verdict="kept", status="pinned",
         why="a real 2-reading detachment rejected as an excursion because a "
             "rise follows it. That rise does not RETURN -- it overshoots "
             "the pre-fall level by 0.0176 AU and then goes flat, which is a "
             "new bubble growing, not a spike reverting",
         source="DATA_VERIFICATION.md 2026-09-10 (second entry)"),
    dict(experiment=138, sample=4, kind="arrival", where=398, gain=0.0296,
         verdict="released", status="pinned",
         why="the 4-reading arrival that grows after that rejected "
             "detachment, 0.0296 AU over readings 395-398",
         source="DATA_VERIFICATION.md 2026-09-10 (second entry)"),
)


def pinned():
    """The rows the current code is expected to satisfy."""
    return tuple(case for case in BUBBLE_CASES if case["status"] == "pinned")


def open_cases():
    """The rows the current code is expected to FAIL."""
    return tuple(case for case in BUBBLE_CASES if case["status"] == "open")
