"""
The activation-sink form, as `scope` carries it, on real curves.

    python data/test_activation_sink.py

`test_summary_kinetics.test_activation_sink_form` holds the fit to planted
curves: it solves its rate law, recovers what was planted, and resolves only
what a run determines. This holds the WIRING on the two-axis block's real
curves, where the checks a planted curve cannot reach live:

  the frame is the object     every column is the `CurveFit.activation_sink`
                              the curves page would draw, not a second fit
  the form is the earned one  an unearned sink is reported at k = 0, with the
                              one-phase column's cost; an earned one never
                              costs more than k = 0
  the convention holds        k <= 1/tau on every curve, the faster
                              relaxation being the catalyst's
  the blank is the flag       `v_act_where_resolved_corrected` is blank
                              exactly where `v_act_resolved_corrected` is False
  the contest adds up         `form_contest` covers every live curve once

`run_gates.py` discovers it.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import activation_sink
import scope

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)
    return ok


def test_the_frame_is_the_fit():
    print("\nthe frame's activation-sink columns are the CurveFit's")
    data = scope.frame()
    fits = scope.fits()
    worst = 0.0
    for row in data.itertuples():
        fit = fits[(row.experiment, row.sample)].activation_sink
        for column, value in (("v_act_corrected", fit.v_act),
                              ("k_sink_corrected", fit.k),
                              ("tau_act_corrected", fit.tau),
                              ("v0_act_corrected", fit.v0)):
            got = getattr(row, column)
            if np.isfinite(value) or np.isfinite(got):
                worst = max(worst, abs(got - value))
    check("every column equals the object the page would draw", worst == 0.0,
          f"worst difference {worst:.1e}")


def test_the_form_is_the_earned_one():
    print("\nthe sink is reported only where it is earned")
    fits = [f.activation_sink for f in scope.fits().values()]
    unearned = [f for f in fits if np.isfinite(f.v_act) and not f.sink_earned]
    earned = [f for f in fits if f.sink_earned]
    check("an unearned sink is reported at k = 0",
          all(f.k == 0 for f in unearned), f"{len(unearned)} curves")
    check("with the one-phase column's own cost",
          all(f.sse == f.sse_no_sink for f in unearned))
    check("an earned sink never costs more than k = 0",
          all(f.sse <= f.sse_no_sink for f in earned), f"{len(earned)} curves")
    check("the faster relaxation is the catalyst's on every curve",
          all(f.k <= (1 / f.tau) * (1 + 1e-9) for f in fits
              if np.isfinite(f.tau)))


def test_the_blank_is_the_flag():
    print("\nv_act_where_resolved_corrected")
    data = scope.frame()
    blank = ~np.isfinite(data.v_act_where_resolved_corrected)
    check("is blank exactly where v_act is unresolved",
          (blank == ~data.v_act_resolved_corrected).all())


def test_the_contest_adds_up():
    print("\nform_contest")
    data = scope.frame()
    table = activation_sink.form_contest(scope.TWO_AXIS_BLOCK)
    check("covers every live curve once",
          int(table.curves.sum()) == int(data.live.sum()),
          f"{int(table.curves.sum())} of {int(data.live.sum())}")
    check("never counts more failures than curves",
          (table.free_asymptote <= table.curves).all())


if __name__ == "__main__":
    test_the_frame_is_the_fit()
    test_the_form_is_the_earned_one()
    test_the_blank_is_the_flag()
    test_the_contest_adds_up()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
