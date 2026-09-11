"""
The rate comparison is only evidence if its baseline IS the published analysis.

    python data/test_rate_choice.py

`rate_choice` puts every headline result on six candidate rates side by side,
so that replacing `vmax` can be decided on evidence. That is worth nothing if
the `vmax` column is a re-implementation that happens to agree today: the
comparison would then be between a copy of the analysis and five other rates,
and the copy would drift. So the gate here is not "the numbers are these" --
the folders' own `check_numbers.py` already hold the published values -- but
"the comparison's baseline cell IS what the published function returns when
called the published way".

Plus the three things the module's plumbing added, each checked against the
code it must agree with:

  the fitted rates    are the drawn fit's own, not a second fit
  the pairwise mask   is a no-op where a candidate is positive everywhere
  peroxide_ladder     now filters on the rate it is asked about, and gives
                      the same ladder as before for `vmax`

`python data/test_rate_choice.py`, and `run_gates.py` discovers it.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import arrhenius
import induction
import ph_role
import rate_choice
import scope

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)
    return ok


def _comparison_cell(table, analysis, rate):
    row = table[(table.analysis == analysis) & (table.rate == rate)]
    return row.iloc[0]


def test_every_candidate_is_a_frame_column():
    """A candidate the frame does not carry would be a column computed here."""
    print("\nevery candidate rate is a column scope.frame builds")
    columns = set(scope.frame(scope.archive()).columns)
    for rate in rate_choice.RATE_CANDIDATES:
        check(f"{rate} is in the frame", rate in columns)
    check("the reference is one of the candidates",
          rate_choice.DIVERGENCE_REFERENCE in rate_choice.RATE_CANDIDATES)


def test_the_baseline_is_the_published_machinery():
    """
    THE vmax COLUMN IS THE PUBLISHED ANALYSIS, CALLED THE PUBLISHED WAY.

    Equality to a direct call, not to a literal: the literals are the folders'
    business and are gated there. What this proves is that the comparison
    routes through the same function with the same arguments, so any other
    column is that analysis on a different rate and nothing else.
    """
    print("\nthe baseline column is the published analysis itself")
    table = rate_choice.order_comparison(rates=("vmax", "vmax_corrected"))
    block = scope.frame()

    fit = scope.orders("vmax_corrected", within=True)
    cell = _comparison_cell(table, "two-axis [H2O2] order, all live", "vmax_corrected")
    check("the two-axis peroxide order is scope.orders' own",
          np.isclose(cell.estimate, fit["order_h2o2"])
          and np.isclose(cell.stderr, fit["stderr_h2o2"]),
          f"{cell.estimate:+.3f} against {fit['order_h2o2']:+.3f}")

    strong = scope.orders("vmax_corrected", scope=scope.strong_runs(),
                          within=True)
    cell = _comparison_cell(table, "two-axis [H2O2] order, strong runs", "vmax_corrected")
    check("...and over the strong runs",
          np.isclose(cell.estimate, strong["order_h2o2"]),
          f"{cell.estimate:+.3f} against {strong['order_h2o2']:+.3f}")

    ph = scope.ph_order("vmax", scope=scope.strong_runs()).loc["pooled"]
    cell = _comparison_cell(table, "two-axis pH order, strong runs", "vmax")
    check("the pH order is scope.ph_order's pooled row",
          np.isclose(cell.estimate, ph.order)
          and np.isclose(cell.stderr, ph.stderr),
          f"{cell.estimate:+.3f} against {ph.order:+.3f}")

    pooled = ph_role.pooled_rate_order(ph_role.rate_ladders("vmax"),
                                       drop=("boric 4OMe",))
    cell = _comparison_cell(table, "[HOO-] order, three ladders pooled", "vmax")
    check("the three-ladder order is ph_role's",
          np.isclose(cell.estimate, pooled["pooled"])
          and np.isclose(cell.extra, pooled["chi2"]),
          f"{cell.estimate:+.3f} chi2 {cell.extra:.2f} against "
          f"{pooled['pooled']:+.3f} chi2 {pooled['chi2']:.2f}")

    energy = arrhenius.pooled_arrhenius("vmax")
    cell = _comparison_cell(table, "activation energy, kJ/mol", "vmax")
    check("the activation energy is arrhenius.pooled_arrhenius'",
          np.isclose(cell.estimate, energy["activation_kJ"]),
          f"{cell.estimate:.3f} against {energy['activation_kJ']:.3f}")

    clocks = induction.joint_clocks(block).loc[("tau_slow_corrected", "axis")]
    cell = _comparison_cell(table, "+1 rule through tau_slow", "vmax_corrected")
    check("the +1 row is joint_clocks' default, which reads vmax_corrected",
          np.isclose(cell.estimate, clocks.order)
          and np.isclose(cell.extra, clocks.sigma),
          f"{cell.estimate:+.3f} ({cell.extra:.2f} sigma) against "
          f"{clocks.order:+.3f} ({clocks.sigma:.2f} sigma)")

    saturation = induction.peroxide_saturation(block)
    cell = _comparison_cell(table, "peroxide saturation order (+/- half profile)", "vmax")
    check("the saturation F is peroxide_saturation's default",
          np.isclose(cell.extra, saturation["first_order_f"]),
          f"F {cell.extra:.2f} against {saturation['first_order_f']:.2f}")


def test_an_injected_ladder_frame_changes_nothing_for_vmax():
    """
    `rate_ladders(frame=...)` is new plumbing; handed the archive's own frame
    it must return what the default path does, ladder by ladder.
    """
    print("\nrate_ladders with an injected frame agrees with its default")
    runs = tuple(sorted(set().union(*scope.PH_LADDERS.values())))
    default = {r["ladder"]: r for r in ph_role.rate_ladders("vmax")}
    injected = {r["ladder"]: r for r in
                ph_role.rate_ladders("vmax", frame=scope.frame(runs))}
    for name, row in default.items():
        other = injected[name]
        check(f"{name}: same [HOO-] order and same curves",
              np.isclose(row["hoo"]["slope"], other["hoo"]["slope"])
              and row["n"] == other["n"],
              f"{row['hoo']['slope']:+.4f} ({row['n']}) against "
              f"{other['hoo']['slope']:+.4f} ({other['n']})")


def test_the_fitted_rates_are_the_drawn_fit():
    """
    The three fitted rates come off `scope.CurveFit.progress_corrected` -- the
    fit the curves pages draw -- and not off a second fit.

    And one identity that ties two of them together: on a one-phase LAG the
    fitted rate rises monotonically to v_ss, so its maximum IS v_ss, reached
    at t = infinity (`ProgressFit.peak_rate`).
    """
    print("\nthe fitted candidates are the drawn fit's own rates")
    table = scope.frame()
    fits = scope.fits()
    worst_ss = worst_v0 = 0.0
    for row in table.itertuples():
        fit = fits[(row.experiment, row.sample)].progress_corrected
        if np.isfinite(fit.v_ss):
            worst_ss = max(worst_ss, abs(row.v_ss_fit_corrected - fit.v_ss))
        rate = float(fit.rate(np.array([0.0]))[0])
        if np.isfinite(rate):
            worst_v0 = max(worst_v0, abs(row.v0_fit_corrected - rate))
    check("v_ss_fit_corrected is progress_corrected.v_ss", worst_ss == 0.0,
          f"worst gap {worst_ss:.1e}")
    check("v0_fit_corrected is progress_corrected.rate(0)", worst_v0 == 0.0,
          f"worst gap {worst_v0:.1e}")
    lags = table[np.isinf(table.v_peak_corrected_time.astype(float))]
    check("on every one-phase lag the fitted peak IS v_ss",
          len(lags) > 0 and np.allclose(lags.v_peak_corrected,
                                        lags.v_ss_fit_corrected),
          f"{len(lags)} curves")


def test_the_pairwise_mask_is_a_pair():
    """
    `order_comparison(common=True)` holds a candidate and the reference to the
    candidate's own usable curves -- and on a candidate positive on every live
    curve that is a no-op, so its row must equal the standard one exactly.
    """
    print("\nthe pairwise mask")
    rates = ("v_peak_corrected", "v0_fit_corrected")
    standard = rate_choice.order_comparison(rates=rates)
    paired = rate_choice.order_comparison(common=True, rates=rates)
    everywhere = paired[paired.rate == "v_peak_corrected"]
    same = standard[standard.rate == "v_peak_corrected"]
    check("a candidate positive everywhere is unmoved by the mask",
          np.allclose(everywhere.estimate.to_numpy(), same.estimate.to_numpy(),
                      equal_nan=True))
    v0 = paired[paired.rate == "v0_fit_corrected"]
    check("the reference is reported beside every paired row",
          v0.reference_estimate.notna().all(), f"{len(v0)} rows")
    joint = rate_choice._common(scope.frame(), ("v0_fit_corrected",
                                                "vmax_corrected"))
    check("and the mask blanks both columns on exactly the same curves",
          (joint.v0_fit_corrected.isna() == joint.vmax_corrected.isna()).all())


def test_peroxide_ladder_honours_its_parameter():
    """
    `peroxide_ladder` filtered on `vmax` whatever rate it was asked about,
    which crashed the first caller to vary it. For `vmax` the ladder must be
    exactly the one it always was.
    """
    print("\nperoxide_ladder filters on the rate it is asked about")
    block = scope.frame()
    before = induction.peroxide_ladder(block)
    # THE LADDER AS IT WAS BUILT BEFORE 2026-09-11, restated here as the pin:
    # the filter was `live & vmax > 0 & h2o2 > 0`, then the top-substrate arm,
    # then runs stepping at least two peroxides.
    old = block[block.live & (block.vmax > 0) & (block.h2o2 > 0)]
    old = old[np.isclose(old.s0, old.groupby("experiment").s0.transform("max"))]
    old = old[old.groupby("experiment").h2o2.transform("nunique") >= 2]
    check("for vmax the ladder is exactly the one the old filter built",
          before.index.equals(old.index), f"{len(before)} curves")
    ladder = induction.peroxide_ladder(block, "v0_fit_corrected")
    check("for a rate that goes negative, only its positive curves are kept",
          len(ladder) > 0 and (ladder.v0_fit_corrected > 0).all(),
          f"{len(ladder)} of {len(before)}")


def test_strong_runs_default_is_unchanged():
    """`strong_runs` grew a `parameter`; its default must be what it was."""
    print("\nstrong_runs' default is still vmax")
    check("strong_runs() is strong_runs(parameter='vmax')",
          scope.strong_runs() == scope.strong_runs(parameter="vmax"),
          f"{scope.strong_runs()}")


if __name__ == "__main__":
    test_every_candidate_is_a_frame_column()
    test_the_baseline_is_the_published_machinery()
    test_an_injected_ladder_frame_changes_nothing_for_vmax()
    test_the_fitted_rates_are_the_drawn_fit()
    test_the_pairwise_mask_is_a_pair()
    test_peroxide_ladder_honours_its_parameter()
    test_strong_runs_default_is_unchanged()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
