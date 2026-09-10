"""
The contract every folder's progress_curves.html is held to.

    python test_progress_panels.py

THERE WERE EIGHT WAYS TO DRAW A PROGRESS CURVE UNTIL 2026-09-12, about 890
lines of them, and they marked eight different vocabularies on the same
object. `two_axis` and `ph` drew v_max/v_max*/tau and were ~70 lines of
byte-comparable code differing in line wrapping; `buffer` and `induction` drew
an induction landmark through two different windows; `product_fate` drew a
rolling rate's maximum; `early_trough` a windowed trough; `temperature_series`
piecewise-linear breakpoints on axes it built itself; `background_reaction`
four estimator lines on axes it built itself and no progress fit at all.

NOT ONE OF THEM MARKED A PARAMETER OF THE FUNCTION IT DREW, and two of them
marked something worse. `scope.frame`'s `tau` is `fit_burst_bounded`'s clock --
the ONE-phase form -- while the line those panels drew is `fit_progress`'s,
and on 36 of the two-axis block's 110 live curves that line was the TWO-phase
form, whose own clocks are `tau_fast` and `tau_slow`. Exp 137.5 carried a rule
labelled "tau 179 s" over a drawn fit whose clocks are 431 s and 6562 s. Exp
138.3's read 9067 s against 3298 s and 3779 s. The panels were not wrong about
`tau`; they were drawing a quantity from a model that curve had rejected, with
nothing on the page to say so, because the fit was made TWICE -- once in
`scope._frame` for the numbers and again in the builder for the line.

The gates below are the two halves of not going back:

    the STRUCTURAL half   no folder may draw a progress panel itself
    the CONTENT half      what is on the page is what the fit says

`run_gates.py` discovers this file by its name, like any other gate.
"""
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "data"))
sys.path.insert(0, HERE)

import scope
from doc_check import Checker
from figure_kit import FITTED_FUNCTION

# Every folder that publishes a curves page. Discovered, not listed, for the
# reason `run_gates.gate_paths` globs: a ninth folder written and never
# checked is the failure this file exists to prevent one level up.
PAGES = sorted(glob.glob(os.path.join(HERE, "*", "progress_curves.html")))

# The primitives that DRAW a progress panel. `figure_kit.progress_panel` is
# the one caller; a folder reaching past it is a ninth vocabulary starting.
# `scratch/` is exempt -- an eyeball page is allowed its own axes, and
# `build_detection_review.py` says in its own docstring why it needs them.
PANEL_PRIMITIVES = ("progress_axes", "residual_axes", "derivative_axes",
                    "draw_progress", "stacked_landmarks", "breakpoints")

# The clock symbols each fitted form is allowed to draw, keyed by phase count.
# A one-phase fit has one relaxation and a two-phase fit has two, and the
# symbols are the ones `FITTED_FUNCTION` prints in the header directly above
# them -- so a reader pairs the rule with the parameter without being told.
CLOCK_SYMBOLS = {1: {"τ"}, 2: {"τ₁", "τ₂"}}

FAILURES = []


def _panels(page_text):
    """Each panel of a page as (header, body)."""
    found = []
    for chunk in page_text.split("<div class='fig panel'>")[1:]:
        header = chunk.split("</div>")[0]
        found.append((header, chunk))
    return found


def test_no_folder_draws_its_own_panel(doc):
    """
    THE STRUCTURAL GATE. Only `figure_kit` may call the panel primitives.

    This is `test_curve_metrics.test_no_duplicate_definitions` applied to
    drawing rather than to naming, and it catches what that guard cannot: the
    eight copies were never eight FUNCTIONS with one name, they were eight
    inlined bodies inside eight `build_curves_page`s, where nothing could see
    them. A duplicate guard that reads top-level definitions is structurally
    blind to a copy that lives inside a function, which is exactly where this
    one lived for as long as it did.
    """
    offenders = []
    for path in sorted(glob.glob(os.path.join(HERE, "*", "build_figures.py"))):
        source = open(path, encoding="utf-8").read()
        # Comments and docstrings may NAME a primitive -- several explain why
        # the folder no longer calls one -- so only a call counts.
        for name in PANEL_PRIMITIVES:
            for hit in re.finditer(rf"\b{name}\s*\(", source):
                line = source[:hit.start()].count("\n") + 1
                offenders.append(f"{os.path.basename(os.path.dirname(path))}"
                                 f"/build_figures.py:{line} {name}(")
    doc.check("no folder calls a progress-panel primitive itself",
              not offenders,
              "clean" if not offenders else
              f"{len(offenders)} calls: {offenders[:4]}")


def test_every_page_is_titled_the_same_way(doc):
    """
    THE TITLE IS A PROMISE and every page makes the same one.

    `index.html` presents the argument and `progress_curves.html` shows the
    fits the argument is read off -- in every folder. Four of the eight titles
    said so and four said something else ("Progress curves and their fits",
    "The 4OMe curves and their induction landmarks"), so a reader had to open
    a page to learn which kind it was. The count belongs in the SUBTITLE,
    where it can change without the title changing; "all 24 progress curves"
    was in a title.
    """
    wrong = []
    for path in PAGES:
        text = open(path, encoding="utf-8").read()
        title = re.search(r"<title>(.*?)</title>", text, re.S)
        name = os.path.basename(os.path.dirname(path))
        if not title or not title.group(1).strip().endswith(
                "— every progress curve"):
            wrong.append(f"{name}: {title.group(1) if title else '(none)'}")
    doc.check("every curves page is '<block> — every progress curve'",
              not wrong, "all eight" if not wrong else str(wrong))


def test_every_panel_prints_its_fitted_function(doc):
    """
    A PANEL SAYS WHICH MODEL IT DREW, in the header, as the function itself.

    Not "the fitted curve" and not a phase count: the algebra, so the clocks
    marked below it have somewhere to be read from.
    """
    forms = set(FITTED_FUNCTION.values())
    missing = []
    total = 0
    for path in PAGES:
        name = os.path.basename(os.path.dirname(path))
        for header, _ in _panels(open(path, encoding="utf-8").read()):
            total += 1
            if not any(form in header for form in forms):
                missing.append(f"{name}: {header[:60]}")
    doc.check("every panel's header prints the fitted function it drew",
              not missing and total > 0,
              f"{total} panels" if not missing else
              f"{len(missing)} without: {missing[:3]}")


def test_every_clock_belongs_to_the_form_above_it(doc):
    """
    THE GATE THE tau DEFECT WOULD HAVE FAILED.

    A panel whose header prints the two-phase function may draw tau_1 and
    tau_2 and may NOT draw a bare tau, because there is no bare tau in that
    function -- and the panels drew one, off `fit_burst_bounded`, on 36 of
    110 curves in `two_axis/` alone. The converse too: a one-phase header may
    not carry tau_1.
    """
    wrong = []
    total = 0
    for path in PAGES:
        name = os.path.basename(os.path.dirname(path))
        for header, body in _panels(open(path, encoding="utf-8").read()):
            phases = next((n for n, form in FITTED_FUNCTION.items()
                           if form in header), None)
            if phases is None:
                continue
            drawn = set(re.findall(r">(τ[₁₂]?) [\d.]+ s<", body))
            total += len(drawn)
            stray = drawn - CLOCK_SYMBOLS[phases]
            if stray:
                wrong.append(f"{name} {header[:24]}: {sorted(stray)} "
                             f"on a {phases}-phase fit")
    doc.check("every clock drawn is a clock of the form its header prints",
              not wrong and total > 0,
              f"{total} clocks drawn" if not wrong else
              f"{len(wrong)}: {wrong[:3]}")


def test_every_clock_is_the_fits_own_value(doc):
    """
    AND IT CARRIES THE FIT'S OWN NUMBER, not a number from somewhere else.

    Checked against `scope.curve_fit` directly rather than against the frame,
    because the frame is what the footer prints and this is about the
    DRAWING. Run over `two_axis/`, which is the block the defect was found on
    and the only page whose every cuvette is one scope.
    """
    path = os.path.join(HERE, "two_axis", "progress_curves.html")
    shapes = scope.fits()
    checked, wrong = 0, []
    for header, body in _panels(open(path, encoding="utf-8").read()):
        label = re.search(r"Exp\. (\d+)\.(\d+)", header)
        if not label:
            continue
        key = (int(label.group(1)), int(label.group(2)))
        fit = shapes.get(key)
        if fit is None:
            continue
        # The chemistry fit is the one drawn and marked: the rebuilt series'
        # where there is gas to take out, the readings' otherwise.
        _, chosen = fit.chemistry
        allowed = {"τ": None, "τ₁": None, "τ₂": None}
        if chosen.phases == 2:
            allowed["τ₁"] = float(chosen.two.tau1)
            allowed["τ₂"] = (float(chosen.two.tau2) if chosen.two.resolved
                             else None)
        else:
            allowed["τ"] = float(chosen.one.tau)
        for symbol, value in re.findall(r">(τ[₁₂]?) ([\d.]+) s<", body):
            checked += 1
            want = allowed.get(symbol)
            if want is None or abs(float(value) - want) > 1.0:
                wrong.append(f"{key} {symbol} drawn {value} s, fit says "
                             f"{want}")
    doc.check("every drawn clock is the drawn fit's own, to a second",
              not wrong and checked > 0,
              f"{checked} clocks on two_axis/" if not wrong else
              f"{len(wrong)}: {wrong[:3]}")


def test_a_non_fitted_landmark_names_its_estimator(doc):
    """
    ANYTHING THAT IS NOT A FITTED PARAMETER SAYS WHERE IT CAME FROM.

    `v_max` is `curve_metrics.peak_rate`, the steepest 20% block slope of the
    readings; `t_ind` is a rolling-window crossing; a breakpoint is a
    piecewise-linear segmentation. None appears in
    `A(t) = c + v_ss.t - sum B_i(1 - e^(-t/tau_i))`, and drawn unlabelled
    beside tau they read as if the fit had produced them. `figure_kit.Mark`
    renders them as "v_max 8.7e-05 (20% block)"; this is the check that none
    arrives bare.

    The vocabulary is CLOSED. A landmark whose source is not in this set is a
    ninth vocabulary starting, so adding one here is a deliberate act.
    """
    known = {"20% block", "20% block, rebuilt", "rolling rate", "windowed z",
             "run/10 window", "piecewise fit",
             f"{450:.0f} s window"}
    bare, unknown = [], set()
    total = 0
    for path in PAGES:
        name = os.path.basename(os.path.dirname(path))
        for _, body in _panels(open(path, encoding="utf-8").read()):
            # THE DRAWING ONLY. A footer is prose and may say "t_ind 402 s"
            # in a sentence; what this is about is the LABEL on the rule,
            # where a bare name reads as the fit's own.
            drawn = "".join(re.findall(r"<svg .*?</svg>", body, re.S))
            for text in re.findall(r"<text [^>]*>([^<>]*(?:v_max|t_ind|"
                                   r"trough|break |rate max)[^<>]*)</text>",
                                   drawn):
                if "(" not in text:
                    bare.append(f"{name}: {text[:40]}")
                    continue
                total += 1
                source = text.rsplit("(", 1)[1].rstrip(")")
                if source not in known:
                    unknown.add(f"{name}: {source}")
    doc.check("every non-fitted landmark names its estimator in brackets",
              not bare, f"{total} labelled" if not bare else
              f"{len(bare)} bare: {bare[:3]}")
    doc.check("...and its estimator is one of the known few",
              not unknown, f"{len(known)} known sources" if not unknown
              else str(sorted(unknown)))


def test_the_gas_is_drawn_in_both_directions(doc):
    """
    ARRIVALS TOO, WHICH NO PUBLISHED PAGE DREW AT ALL UNTIL 2026-09-12.

    `bubble_arrivals` and `split_arrivals` are the correction that moved
    `joint_clocks`' `tau_slow` row to +0.757 +/- 0.289 (CLAUDE.md,
    DATA_VERIFICATION.md 2026-09-10), and the only place the archive's 80
    arrivals were ever drawn was `scratch/build_detection_review.py`, an
    eyeball page. The audit surfaces showed the falls and not the rises --
    half of a detection layer, on the pages that exist to audit it.
    `figure_kit.ARRIVAL_BAND_COLOUR` was meanwhile #8a5aa8, byte-identical to
    the reconstruction's own line colour, and nothing caught that either
    because nothing used it.
    """
    from figure_kit import ARRIVAL_BAND_COLOUR, EVENT_BAND_COLOUR, \
        CHEMISTRY_COLOUR
    doc.check("the arrival band is not the reconstruction's own colour",
              ARRIVAL_BAND_COLOUR != CHEMISTRY_COLOUR,
              f"{ARRIVAL_BAND_COLOUR} vs {CHEMISTRY_COLOUR}")
    path = os.path.join(HERE, "two_axis", "progress_curves.html")
    shapes = scope.fits()
    missing_falls, missing_rises = [], []
    falls = rises = 0
    for header, body in _panels(open(path, encoding="utf-8").read()):
        label = re.search(r"Exp\. (\d+)\.(\d+)", header)
        if not label:
            continue
        fit = shapes.get((int(label.group(1)), int(label.group(2))))
        if fit is None:
            continue
        if fit.events:
            falls += 1
            if EVENT_BAND_COLOUR not in body:
                missing_falls.append(label.group(0))
        if fit.arrivals:
            rises += 1
            if ARRIVAL_BAND_COLOUR not in body:
                missing_rises.append(label.group(0))
    doc.check("every curve with a detachment shades it",
              not missing_falls and falls > 0,
              f"{falls} curves" if not missing_falls else str(missing_falls))
    doc.check("every curve with an arrival shades it too",
              not missing_rises and rises > 0,
              f"{rises} curves" if not missing_rises else str(missing_rises))


def test_the_frame_and_the_panel_share_one_fit(doc):
    """
    ONE COMPUTATION, WHICH IS WHAT MAKES THE CHECKS ABOVE POSSIBLE.

    `scope.curve_fit` is memoised and `scope._frame` reads it, so the object a
    panel draws IS the object the frame's columns were flattened from. Asserted
    by identity, not by value: two equal fits made twice is the arrangement
    that let the clocks drift apart in the first place.
    """
    curve = scope.curves()[0]
    first, second = scope.curve_fit(curve), scope.curve_fit(curve)
    doc.check("curve_fit is memoised, not recomputed", first is second,
              "same object")
    doc.check("...and hands out read-only arrays",
              not first.times.flags.writeable
              and not first.values.flags.writeable
              and not first.corrected.flags.writeable,
              "times, values and corrected are frozen")
    block = scope.fits()
    table = scope.frame()
    doc.check("fits() and frame() cover the same curves",
              set(block) == set(zip(table["experiment"], table["sample"])),
              f"{len(block)} curves")


def main():
    doc = Checker(os.path.join(HERE, "CLAUDE.md"),
                  label="the curves pages",
                  document_label="the shared progress panel")
    doc.section("The shared progress panel")
    # CALLED BY NAME, one line each, and not through a tuple of function
    # objects: `data/test_curve_metrics.test_every_test_is_actually_run` reads
    # this file's AST for a Call whose func is a Name, so a runner that
    # iterates a tuple defines eight tests and calls none of them as far as
    # that guard can see -- which is exactly the "defined but never called"
    # failure it exists to catch, arriving through the back door.
    test_no_folder_draws_its_own_panel(doc)
    test_every_page_is_titled_the_same_way(doc)
    test_every_panel_prints_its_fitted_function(doc)
    test_every_clock_belongs_to_the_form_above_it(doc)
    test_every_clock_is_the_fits_own_value(doc)
    test_a_non_fitted_landmark_names_its_estimator(doc)
    test_the_gas_is_drawn_in_both_directions(doc)
    test_the_frame_and_the_panel_share_one_fit(doc)
    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
