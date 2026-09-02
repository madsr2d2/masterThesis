"""
Tests for doc_check.py, the contract every folder's check_numbers.py runs on.

This one is worth having because the thing it guards is invisible when it
breaks. A checker that folds too much passes on a document that disagrees with
the code, and prints "agrees with the code" while doing it -- which is exactly
what happened: three folders stripped `**` before comparing while their claims
built it, so the emphasis was asserted and then discarded, and the runs stayed
green for as long as that was true.

So the tests below are in two halves. Half one: every fold that SHOULD happen
does. Half two -- the half that matters -- the folds that should NOT happen do
not, planted as a document that is wrong in exactly one way.

    python test_doc_check.py
"""
import contextlib
import io
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from doc_check import Checker, normalise

FAILURES = []


def check(name, condition, detail=""):
    print(f"  {'pass' if condition else 'FAIL'}  {name}"
          f"{': ' + detail if detail else ''}")
    if not condition:
        FAILURES.append(f"{name} {detail}")


@contextlib.contextmanager
def _quiet():
    """
    Swallow a planted Checker's own output.

    Its lines say FAIL by design -- that is what is being tested -- and a log
    full of expected failures is a log nobody reads carefully.
    """
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def _document(text):
    """A throwaway document, so a Checker can be pointed at real content."""
    handle = tempfile.NamedTemporaryFile("w", suffix=".md", delete=False,
                                         encoding="utf-8")
    handle.write(text)
    handle.close()
    return handle.name


def test_typography_folds_onto_what_a_format_string_emits():
    """Every difference between how a number is typed and how it is printed."""
    print("\nwhat folds")
    for label, typed, emitted in (
            ("the minus sign", "−0.433", "-0.433"),
            ("the en dash in a range", "50–200 mM", "50-200 mM"),
            ("plus-minus", "±0.201", "+/-0.201"),
            ("the multiplication sign", "1.80×", "1.80x"),
            ("a subscript zero", "v₀", "v0"),
            ("a subscript two", "H₂O₂", "H2O2"),
            ("subscript ss", "vₛₛ", "vss"),
            ("a superscript exponent", "3.29 × 10⁻⁵", "3.29 x 10-5"),
            ("an arrow", "6489 → 3190", "6489 -> 3190"),
            ("tau", "τ_fast", "tau_fast"),
            ("sigma", "3.7σ", "3.7sigma"),
            ("chi-squared", "χ² = 4", "chi2 = 4"),
            ("a degree sign", "40 °C", "40 C"),
            ("micro against mu", "µM", "μM"),
            ("delta and the double dagger", "ΔG‡", "DG"),
            ("a rewrapped line", "five figures, A to\nE",
             "five figures, A to E"),
            ("code formatting", "`[buf]`", "[buf]")):
        check(label, normalise(typed) == normalise(emitted),
              f"{normalise(typed)!r} vs {normalise(emitted)!r}")


def test_emphasis_survives_because_it_is_an_argument():
    """
    The one thing that must NOT fold, and the reason this module exists.

    In these documents bold marks the number a section argues for. Folding it
    away turns `bold = "**" if ... else ""` into an assertion that cannot fail.
    """
    print("\nwhat does not fold")
    check("bold is not stripped",
          normalise("**+1.094**") != normalise("+1.094"),
          f"{normalise('**+1.094**')!r}")
    document = _document("The pooled slope is **−0.433 ± 0.201** on eight "
                         "curves, and +0.371 ± 0.018 in exp 32.\n")
    with _quiet():
        doc = Checker(document, label="planted")
        doc.claim("bolded, claimed with its emphasis", "**-0.433 +/- 0.201**")
        doc.claim("unbolded, claimed without", "+0.371 +/- 0.018")
    check("a claim that matches the document's emphasis passes",
          not doc.failures, f"{doc.failures}")

    with _quiet():
        loud = Checker(document, label="planted")
        loud.claim("unbolded in the document, claimed as bold",
                   "**+0.371 +/- 0.018**")
    check("a claim that asserts emphasis the document does not carry FAILS",
          len(loud.failures) == 1, f"{loud.failures}")

    # THE CONTRACT IS ASYMMETRIC, and it is worth writing down rather than
    # discovering later. `**` sits OUTSIDE the number, so a claim that omits it
    # is still a substring of a bolded document and still matches. What that
    # buys is the direction that matters: the CODE declares which number is the
    # headline, by building the markers, and the DOCUMENT has to carry them.
    # The reverse -- a document that bolds something the code did not call a
    # headline -- is not caught, and no folder needs it to be.
    with _quiet():
        lenient = Checker(document, label="planted")
        lenient.claim("bolded, claimed without its emphasis",
                      "-0.433 +/- 0.201")
    check("and a claim that merely omits emphasis still matches",
          not lenient.failures,
          "the markers are outside the number, so this is a substring")
    os.unlink(document)


def test_a_wrong_number_is_caught_and_says_where():
    """The base case, and that the failure names the document it missed."""
    print("\na planted disagreement")
    document = _document("The gap is **−11.99 ± 0.62 kJ/mol**.\n")
    with _quiet():
        doc = Checker(document, label="planted")
        doc.claim("as published", "**-11.99 +/- 0.62 kJ/mol**")
    check("the right number passes", not doc.failures)
    with _quiet():
        doc.claim("drifted by one digit", "**-11.98 +/- 0.62 kJ/mol**")
        code = doc.summary()
    check("the wrong number fails", len(doc.failures) == 1)
    check("and the failure names the file it was not found in",
          os.path.basename(document) in doc.failures[0], doc.failures[0])
    check("the summary returns a non-zero exit code", code == 1)
    os.unlink(document)


def test_absence_can_be_asserted_too():
    """`present=False` is how a withdrawn claim is kept withdrawn."""
    print("\nasserting absence")
    document = _document("The argument is withdrawn.\n")
    with _quiet():
        doc = Checker(document, label="planted")
        doc.claim("the withdrawn number is gone", "0.961", present=False)
    check("absence passes when the text is absent", not doc.failures)
    with _quiet():
        doc.claim("a phrase that IS present, asserted absent",
                  "withdrawn", present=False)
    check("and fails when it is not", len(doc.failures) == 1,
          f"{doc.failures}")
    check("the message says it is still present",
          "still present in" in doc.failures[0], doc.failures[0])
    os.unlink(document)


def test_the_figure_letters_catch_a_repeat_and_an_order():
    """
    Read from the rendered page, which is the only place a repeat is visible.

    Two folders read the BUILDER source and one of them sorted the result, so
    neither a duplicate letter nor a letter out of order could be seen. The
    induction page drew A B C D E G H F I when this was written.
    """
    print("\nfigure letters")
    good = _document("<h3>A · one</h3><h3>B · two</h3><h3>C · three</h3>")
    with _quiet():
        doc = Checker(good, label="planted")
        doc.figures(good, "ABC")
    check("A B C in order passes", not doc.failures, f"{doc.failures}")

    swapped = _document("<h3>A · one</h3><h3>C · three</h3><h3>B · two</h3>")
    with _quiet():
        out_of_order = Checker(swapped, label="planted")
        out_of_order.figures(swapped)
    check("A C B is caught even with no expected set given",
          len(out_of_order.failures) == 1, f"{out_of_order.failures}")

    repeated = _document("<h3>A · one</h3><h3>B · two</h3><h3>B · again</h3>")
    with _quiet():
        twice = Checker(repeated, label="planted")
        twice.figures(repeated)
    check("a letter used twice is caught", len(twice.failures) == 1,
          f"{twice.failures}")
    for path in (good, swapped, repeated):
        os.unlink(path)


def test_every_folder_runs_this_contract_and_no_other():
    """No folder may keep a private copy of the comparison."""
    print("\nno folder has its own copy")
    here = os.path.dirname(os.path.abspath(__file__))
    folders = sorted(name for name in os.listdir(here)
                     if os.path.exists(os.path.join(here, name,
                                                    "check_numbers.py")))
    check("all five analysis folders are present", len(folders) == 5,
          f"{folders}")
    for folder in folders:
        source = io.open(os.path.join(here, folder, "check_numbers.py"),
                         encoding="utf-8").read()
        check(f"{folder} imports the shared contract",
              "from doc_check import" in source)
        check(f"{folder} defines no comparison of its own",
              "def _normalise" not in source and "def claim(" not in source
              and "\nFAILURES" not in source)


if __name__ == "__main__":
    test_typography_folds_onto_what_a_format_string_emits()
    test_emphasis_survives_because_it_is_an_argument()
    test_a_wrong_number_is_caught_and_says_where()
    test_absence_can_be_asserted_too()
    test_the_figure_letters_catch_a_repeat_and_an_order()
    test_every_folder_runs_this_contract_and_no_other()
    print(f"\n{len(FAILURES)} failures")
    sys.exit(1 if FAILURES else 0)
