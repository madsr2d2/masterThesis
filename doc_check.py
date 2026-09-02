"""
The shared contract behind every folder's `check_numbers.py`.

Each analysis folder verifies that the numbers in its `ANALYSIS.md` are the ones
its modules actually produce. That check is only as good as the comparison
underneath it, and until 2026-09-02 there were FIVE copies of that comparison
and no two were the same:

    background_reaction   7 substitutions, folded neither tau nor sigma
    temperature_series   14, and did not strip markup
    product_fate         16
    buffer               17, stripped `**`
    induction            17, stripped `**`, folded every superscript digit

So five documents were being held to five different standards, and the folder
with the most numbers in it ran the weakest one. Two consequences were live:

  * `background_reaction` could not compare a Greek letter in prose against the
    `tau`/`sigma` a format string emits, so any claim carrying one had to be
    written around the checker.
  * `buffer`, `induction` and `product_fate` STRIPPED `**` while their claims
    built it -- `f"| {bold}{value:+.2f}{bold} |"` with `bold = "**" if ...`.
    The markers were constructed and then removed before the comparison, so the
    emphasis was never checked. Dead assertions that read like live ones.

This module is the union of the five, with one deliberate asymmetry.

WHAT IS FOLDED AWAY, AND WHY. Everything here is a difference between how a
number is TYPED and how `%`-formatting EMITS it: U+2212 against a hyphen, `±`
against `+/-`, `τ` against `tau`, superscripts against `10-5`, `H₂O₂` against
`H2O2`. Whitespace runs collapse too, so a claim does not fail because the prose
rewrapped. None of that is a number, and a checker that failed on it would fail
on typography while a genuinely wrong value typed with a hyphen went through.

WHAT IS NOT FOLDED AWAY. **Emphasis.** In these documents bold is not
decoration -- it marks the number a section is arguing for, and the convention
is followed closely enough to be worth enforcing. So `**` survives
normalisation and a claim has to carry the document's own emphasis. Backticks do
not survive: code formatting is noise, `[buf]` and `` `[buf]` `` are the same
claim.

    from doc_check import Checker
    doc = Checker(os.path.join(HERE, "ANALYSIS.md"))
    doc.claim("the pooled slope", f"**{fit['slope']:+.3f} +/- {fit['stderr']:.3f}**")
    doc.check("and it is negative", fit["slope"] < 0)
    raise SystemExit(doc.summary())
"""
import io
import os
import re

# Single characters that fold one-for-one. Applied before the multi-character
# rules below, none of which they can interfere with.
CHARACTERS = {
    "−": "-",      # minus sign, which reads better in a table than a hyphen
    "–": "-",      # en dash
    "±": "+/-",    # plus-minus
    "×": "x",      # multiplication sign
    "₀": "0",      # subscript zero, as in v0
    "₂": "2",      # subscript two, as in H2O2
    "→": "->",     # arrow, as in tau falling
    "÷": "/",
    "Δ": "D",      # capital delta, as in dG
    "‡": "",       # double dagger, as in dG-doubledagger
    "τ": "tau",
    "σ": "sigma",
    "µ": "u",      # micro sign
    "μ": "u",      # greek mu, which is a different codepoint
    "°": "",       # degree, so "40 C" matches "40 °C"
    "`": "",            # code formatting is noise; emphasis is not -- see above
}

# Longer than one character, so they run first: folding the superscript two of
# "chi-squared" to a digit would leave "chi2" unreachable.
SEQUENCES = (
    ("χ²", "chi2"),       # chi-square, as typed
    ("ₛₛ", "ss"),         # subscript ss, as in v_ss
)

# Superscript digits and minus, so "10⁻⁵" folds onto the "10-5" a
# format string emits.
SUPERSCRIPTS = str.maketrans(
    "⁻⁰¹²³⁴⁵⁶⁷⁸⁹",
    "-0123456789")


def normalise(text):
    """Fold typography onto what %-formatting emits. The contract is the docstring."""
    for source, target in SEQUENCES:
        text = text.replace(source, target)
    for source, target in CHARACTERS.items():
        text = text.replace(source, target)
    return " ".join(text.translate(SUPERSCRIPTS).split())


class Checker:
    """
    One folder's document check: assertions, their failures, and an exit code.

    A class rather than five module-level functions with a shared `FAILURES`
    list, because `background_reaction` checks two documents and the list has to
    know which one a claim missed. `default` is the document a bare `claim` goes
    to; pass `document=` for any other.
    """

    def __init__(self, default, label=None):
        self.default = default
        self.label = label or os.path.basename(os.path.dirname(
            os.path.abspath(default)))
        self.failures = []
        self.claims = 0
        self._cache = {}

    def text(self, document=None):
        """The normalised document, read once per path."""
        path = document or self.default
        if path not in self._cache:
            self._cache[path] = normalise(
                io.open(path, encoding="utf-8").read())
        return self._cache[path]

    def claim(self, label, rendered, present=True, document=None):
        """Require `rendered` to appear in the document, as the code renders it."""
        path = document or self.default
        ok = (normalise(rendered) in self.text(path)) == present
        self.claims += 1
        print(f"  {'pass' if ok else 'FAIL'}  {label}: {rendered!r}")
        if not ok:
            self.failures.append(
                f"{label} -- {rendered!r} "
                f"{'not found in' if present else 'still present in'} "
                f"{os.path.basename(path)}")
        return ok

    def check(self, label, ok, detail=""):
        """Require a condition that is not a quoted string."""
        ok = bool(ok)
        print(f"  {'pass' if ok else 'FAIL'}  {label}"
              f"{': ' + detail if detail else ''}")
        if not ok:
            self.failures.append(f"{label} {detail}")
        return ok

    def fail(self, message):
        """
        Record a failure whose line was printed by the caller.

        For a guard that wants a different sentence in the summary than the one
        it printed -- "section 6b's degeneracy argument needs rewriting" rather
        than the value that triggered it. `check` covers everything else.
        """
        self.failures.append(message)

    def section(self, title):
        print(f"\n{title}")

    def figures(self, page, expected=None):
        """
        The figure letters a page actually draws: once each, A onwards, in order.

        READ FROM THE RENDERED PAGE, not from the builder's source. Three of the
        five folders did this and two read the source instead, which cannot see
        a letter used twice -- and one folder carried J and K twice for as long
        as it had a section 3a, because a letter is not a number and no check
        looked at one.

        `expected` names the letters the document promises, when the document
        says how many there are. Without it this still catches a gap, a repeat
        and a letter out of order.
        """
        letters = re.findall(r">([A-Z]) \u00b7 ",
                             io.open(page, encoding="utf-8").read())
        ordered = [chr(ord("A") + index) for index in range(len(letters))]
        self.check("every figure letter is used exactly once, A onwards",
                   letters == ordered,
                   f"{''.join(letters)} against {''.join(ordered)}")
        if expected is not None:
            wanted = list(expected)
            self.check(f"{len(wanted)} figures, {wanted[0]} to {wanted[-1]}",
                       letters == wanted, f"{''.join(letters)}")
        return letters

    def unclipped(self, *pages):
        """No data point drawn outside its own frame, on each page that exists."""
        from svgplot import clipped_marks
        for path in pages:
            if not os.path.exists(path):
                continue
            lost = clipped_marks(io.open(path, encoding="utf-8").read())
            name = os.path.basename(path)
            ok = not lost
            print(f"  {'pass' if ok else 'FAIL'}  {name}: {len(lost)} clipped")
            if not ok:
                for title, x, y, *_ in lost[:4]:
                    print(f"        {title!r} at ({x:.4g},{y:.4g})")
                self.failures.append(
                    f"{name} draws {len(lost)} point(s) outside the plot "
                    f"frame; widen the axis limits")

    def summary(self):
        """Print the verdict and return the exit code."""
        print()
        if self.failures:
            print(f"{len(self.failures)} MISMATCH(ES):")
            for item in self.failures:
                print(f"  - {item}")
            return 1
        print(f"{self.label}/ANALYSIS.md agrees with the code "
              f"({self.claims} claims)")
        return 0
