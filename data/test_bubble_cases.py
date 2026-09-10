"""
Run `bubble_cases.BUBBLE_CASES` against the code.

Two assertions, not one. The PINNED rows must pass -- they are what the
thresholds were pinned to, and a change that moves one is a regression. The
OPEN rows must still FAIL: they record behaviour believed wrong, so if one
starts passing, somebody has fixed it and the row needs promoting rather than
quietly agreeing. A fixture that only checked the happy half would let the
known-wrong list rot.

`python data/test_bubble_cases.py`, and `run_gates.py` discovers it.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import bubble_cases
import curve_metrics
import scope

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)
    return ok


def verdict_of(case):
    """What the code says about this row, in the row's own vocabulary."""
    curve = next((c for c in scope.curves_of(case["experiment"])
                  if c.sample == case["sample"]), None)
    if curve is None:
        return "curve missing"
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    events = curve_metrics.detachments(values, curve.noise)
    kind = case["kind"]
    if kind == "fall":
        start, stop = case["where"]
        return "kept" if (start, stop) in events else "rejected"
    if kind == "arrival":
        arrivals = curve_metrics.bubble_arrivals(times, values, curve.noise,
                                                 events)
        released, unreleased = curve_metrics.split_arrivals(arrivals, events)
        for name, group in (("released", released), ("unreleased", unreleased)):
            for index, gain in group:
                if index != case["where"]:
                    continue
                # Where the row states a gain, the AMOUNT is the claim: a
                # multi-reading arrival credited with one step of itself sits
                # at the right reading and is still the wrong event.
                if "gain" in case and abs(gain - case["gain"]) > 0.002:
                    return f"{name} but {gain:.4f} AU"
                return name
        return "none"
    if kind == "no arrival in":
        low, high = case["where"]
        arrivals = curve_metrics.bubble_arrivals(times, values, curve.noise,
                                                 events)
        return ("none" if not any(low <= index <= high for index, _ in arrivals)
                else "some")
    if kind == "curve":
        rate = curve_metrics.bubble_rate(times, values, events)
        if not np.isfinite(rate):
            return "untouched"
        if not events:
            return "no events"
        return ("kept", len(events))
    raise ValueError(f"unknown kind {kind!r}")


def _describe(case):
    where = case["where"]
    place = "" if where is None else f" at {where}"
    return (f"exp {case['experiment']}.{case['sample']} "
            f"{case['kind']}{place}")


def test_the_pinned_cases_still_hold():
    """
    Every curve a threshold in this module was pinned to, checked at once.

    These are not arbitrary regression values: each row is a curve somebody
    looked at, and the constant that decides it was chosen so that this row
    comes out the way it does. Moving any of them means the calibration
    argument has changed and DATA_VERIFICATION.md needs an entry saying so.
    """
    print("\nthe pinned bubble cases")
    for case in bubble_cases.pinned():
        got = verdict_of(case)
        check(_describe(case) + f" -> {case['verdict']}",
              got == case["verdict"], f"got {got}")


def test_the_open_cases_still_fail():
    """
    The known-wrong list, asserted to still be wrong.

    Backwards on purpose. Each of these was found by eye and argued out in
    DATA_VERIFICATION.md, and the code has not been changed to match. If one
    starts agreeing, that is progress and the row's `status` should become
    "pinned" -- but it has to be noticed, and a fixture nobody re-reads is
    exactly how a known-wrong list turns into a forgotten one.
    """
    print("\nthe open bubble cases, still open")
    for case in bubble_cases.open_cases():
        got = verdict_of(case)
        check(_describe(case) + f" still disagrees (wants {case['verdict']})",
              got != case["verdict"],
              f"got {got} -- if this is right now, promote the row to pinned")


def test_every_case_names_its_evidence():
    """
    A row with no source is a number somebody remembered, which is the thing
    this file exists to stop.
    """
    print("\nevery case cites where its evidence is written down")
    for case in bubble_cases.BUBBLE_CASES:
        check(_describe(case) + " cites a source",
              bool(case.get("source")) and bool(case.get("why")))
    kinds = {case["kind"] for case in bubble_cases.BUBBLE_CASES}
    check("no unknown kinds crept in",
          kinds <= {"fall", "arrival", "curve", "no arrival in"}, f"{kinds}")


if __name__ == "__main__":
    test_the_pinned_cases_still_hold()
    test_the_open_cases_still_fail()
    test_every_case_names_its_evidence()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
