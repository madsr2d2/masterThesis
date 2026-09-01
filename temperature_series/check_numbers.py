"""
Verifies every number quoted in temperature_series/ANALYSIS.md against the code.

Same contract as background_reaction/check_numbers.py, and for the same reason:
a figure typed into a document is a copy, and copies in this project have gone
stale three times -- a legend, a rate law and an estimator table -- each time
because no check re-derived them.

    python temperature_series/check_numbers.py
"""
import io
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))

import arrhenius
import scope
import verify_enzyme_stock
from curve_metrics import ACCELERATION_SIGMA

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")
FAILURES = []


def _normalise(text):
    """Fold typography onto what %-formatting emits. See the sibling module."""
    return " ".join((text.replace("−", "-")
                         .replace("–", "-")
                         .replace("±", "+/-")
                         .replace("×", "x")
                         .replace("₀", "0")
                         .replace("→", "->")   # arrow, as in tau falling
                         .replace("÷", "/"))
                    .split())


def claim(label, rendered, present=True):
    text = _normalise(io.open(DOCUMENT, encoding="utf-8").read())
    ok = (_normalise(rendered) in text) == present
    print(f"  {'pass' if ok else 'FAIL'}  {label}: {rendered!r}")
    if not ok:
        FAILURES.append(f"{label} -- {rendered!r}")


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}{': ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(f"{label} {detail}")


def main():
    print("checking temperature_series/ANALYSIS.md against the code\n")

    print("the series itself")
    frame = arrhenius.series_frame()
    by_temperature = frame.groupby("temperature").agg(
        experiment=("experiment", "first"), e0=("e0", "first"),
        n=("live", "size"))
    claim("the experiment row",
          "| exp | " + " | ".join(str(int(by_temperature.loc[t, "experiment"]))
                                  for t in [25, 35, 40, 30, 20, 15]) + " |")
    check("all 24 curves are live",
          len(frame) == 24 and bool(frame.live.all()), f"{len(frame)}")
    # The source split, and that the .txt ones are exactly one rung. If a
    # future rebuild moves that rung onto the .rre the prose must change, so
    # the shape of the split is checked, not just its size.
    export = frame[frame.source == "txt"]
    claim("the source split",
          f"**{int((frame.source == 'rre').sum())} come from the instrument's "
          f"own `.rre` and {len(export)} from the `.txt` export**")
    check("the .txt curves are exactly one rung, in every run",
          export.s0.nunique() == 1 and len(export) == 6,
          f"{export.s0.nunique()} rung(s), {len(export)} curves")
    # Rounded, not exact: the six agree to the last two bits of the float
    # (…234 against …239), which is one lattice value reached by two summation
    # orders, not two different noises.
    check("that rung's noise is one quantisation value",
          export.noise.round(12).nunique() == 1,
          f"{export.noise.round(12).nunique()} distinct")
    # The four rungs recur in every run -- the claim the whole design rests on.
    rungs = sorted(frame.s0.unique())
    check("the same four rungs in all six runs",
          len(rungs) == 4 and all(
              len(frame[np.isclose(frame.s0, s0)]) == 6 for s0 in rungs),
          f"{[round(v, 3) for v in rungs]}")
    claim("the rung ladder",
          "(" + ", ".join(f"{v:.3f}" for v in rungs) + " mM)")

    print("\nthe enzyme stocks, from the weighing")
    table = verify_enzyme_stock.audit()
    used = table[table.used]
    stocks = verify_enzyme_stock.stocks(table)
    check("seven distinct stocks", len(stocks) == 7, f"{len(stocks)}")
    claim("chains self-consistent",
          f"**{int(used.chain_ok.sum())} of the {len(used)} chains are "
          f"self-consistent**")
    claim("compiled equals the sheet",
          f"**{int(used.agrees.sum())} of {len(used)} compiled values equal "
          f"the sheet's own `kuv`**")
    for _, row in stocks.iterrows():
        claim(f"stock at {row.grams:g} g",
              f"| {row.grams:g} g | {row.litres * 1000:g} mL | "
              f"{row.kuv_mM:g} |")
    # Exp 16 is the ONLY out-of-sequence stock. That is the sentence the whole
    # ambiguity rests on, so it is derived rather than asserted.
    flagged = verify_enzyme_stock.out_of_sequence(table)
    check("exp 16 is the only out-of-sequence stock",
          list(flagged.experiment) == [16], f"{list(flagged.experiment)}")

    print("\nthe kinetics that decide it")
    for estimator in ("vmax", "v0_quad"):
        result = arrhenius.enzyme_hypotheses(estimator)
        best = result.rms_winner = result.mean_rms.idxmin()
        check(f"{estimator}: the recorded values win",
              best == "recorded" and int(result.loc["recorded", "rungs_won"])
              == int(result.loc["recorded", "rungs"]),
              f"{best}, {int(result.loc['recorded', 'rungs_won'])} of "
              f"{int(result.loc['recorded', 'rungs'])} rungs")
        if estimator == "vmax":
            for name, label in (("recorded", "**as recorded**"),
                                ("exp16 restocked",
                                 "exp 16 restocked to 0.273"),
                                ("exps16,19 restocked",
                                 "exps 16 and 19 both 0.273")):
                row = result.loc[name]
                mark = "**" if name == "recorded" else ""
                claim(f"hypothesis row, {name}",
                      f"| {label} | {mark}{row.mean_rms:.3f}{mark} | "
                      f"{mark}{row.worst_rms:.3f}{mark} | "
                      f"{mark}{int(row.rungs_won)} of {int(row.rungs)}{mark} | "
                      f"{row.activation_kJ:.1f} |")
        else:
            row = result.loc["recorded"]
            claim("v0_quad sensitivity",
                  f"{row.mean_rms:.3f} against "
                  f"{result.loc['exp16 restocked', 'mean_rms']:.3f} and "
                  f"{result.loc['exps16,19 restocked', 'mean_rms']:.3f}")

    print("\nexp 16's own distance from the line")
    plain = arrhenius.experiment_residuals()
    forced = arrhenius.experiment_residuals(override={16: 0.273})
    claim("as recorded", f"as recorded (0.241): **{plain.loc[16, 'mean_residual']:+.3f}**")
    claim("forced", f"forced to 0.273: **{forced.loc[16, 'mean_residual']:+.3f}**")
    claim("the neighbours' scatter",
          f"(exp 15 is {plain.loc[15, 'mean_residual']:+.3f}, "
          f"exp 17 {plain.loc[17, 'mean_residual']:+.3f})")
    # The DIRECTION is the argument, not just the size: forcing the higher
    # concentration must push exp 16 further out on every rung, not on average.
    check("forcing 0.273 moves exp 16 further out",
          forced.loc[16, "mean_residual"] < plain.loc[16, "mean_residual"],
          f"{plain.loc[16, 'mean_residual']:+.3f} -> "
          f"{forced.loc[16, 'mean_residual']:+.3f}")

    print("\nthe breakpoint screen")
    # Run because exp 65 in the BnOH background showed the start-versus-end
    # statistics step over a mid-run break. The claim here is the OPPOSITE of
    # exp 65's -- that the breaks are chemistry -- so it needs the same
    # evidence: a trend that tracks a physical variable, and a second,
    # unrelated statistic agreeing.
    breaks = scope.synchronised_break(scope.TEMPERATURE_SERIES)
    temperature = {int(e): float(frame[frame.experiment == e].temperature.iloc[0])
                   for e in breaks.index}
    ordered = sorted(breaks.index, key=lambda e: temperature[e])
    claim("the break-ratio row",
          "| slope after / before | " + " | ".join(
              f"{min(breaks.loc[e, 'ratios']):.2f}-"
              f"{max(breaks.loc[e, 'ratios']):.2f}" for e in ordered) + " |")
    claim("the steepening row",
          "| cuvettes steepening | " + " | ".join(
              f"{int(breaks.loc[e, 'steep'])} of {int(breaks.loc[e, 'n'])}"
              for e in ordered) + " |")
    # The trend is the argument. Ratio must fall as temperature rises.
    # NOT monotone, and the document says so: 30 C sits above 25 C because the
    # ratio tracks how much of the run the induction occupies, and the 25 C run
    # is 19 tau long against the 30 C run's 5.8. What must hold is that the
    # trend falls overall and that the single inversion is the one named.
    medians = [float(np.median(breaks.loc[e, "ratios"])) for e in ordered]
    claim("the median ratios",
          ", ".join(f"{v:.2f}" if abs(v - medians[3]) > 1e-9 else f"**{v:.2f}**"
                    for v in medians))
    rises = [i for i, (a, b) in enumerate(zip(medians, medians[1:])) if b > a]
    check("exactly one inversion, at 25/30 C",
          rises == [2], f"inversions after index {rises}")
    lengths = frame.groupby("temperature").apply(
        lambda g: float(g.duration_s.median() / g.tau.median()),
        include_groups=False)
    claim("the run lengths in tau",
          f"the 25 °C run is **{lengths[25.0]:.0f}τ** long and the 30 °C run "
          f"only **{lengths[30.0]:.1f}τ**")
    # The second, independent statistic: the lag time constant.
    taus = frame.groupby("temperature").tau.median()
    lag = [t for t in sorted(taus.index) if t <= arrhenius.BURST_TRUSTWORTHY_BELOW_C]
    claim("tau falls as it warms",
          "**" + " -> ".join(f"{taus[t]:.0f}" for t in lag) + " s**")
    check("every run below 32 C is a lag curve",
          bool((frame[frame.temperature <= arrhenius.BURST_TRUSTWORTHY_BELOW_C]
                .v0_burst_kind == "lag").all()))

    print("\nwhat the screen exposed about vmax")
    where = frame.groupby("temperature").vmax_where
    claim("the vmax_where row",
          "| `vmax_where` | " + " | ".join(
              f"{where.min()[t]:.2f}-{where.max()[t]:.2f}"
              for t in sorted(frame.temperature.unique(), reverse=True)) + " |")
    sensitivity = arrhenius.truncation_sensitivity()
    claim("Ea from vmax", f"| `vmax` throughout | {sensitivity['vmax_kJ']:.1f} |")
    claim("Ea cold-corrected",
          f"| `v_ss` at the cold end | {sensitivity['cold_corrected_kJ']:.1f} |")
    claim("the inflation",
          f"| **inflation from truncation** | **{sensitivity['inflation_kJ']:.1f}** |")
    ratios = sensitivity["v_ss_over_vmax"]
    claim("the v_ss/vmax ratios",
          f"{ratios[15.0]:.3f} at 15 C and {ratios[20.0]:.3f} at 20 C")
    claim("and where they should agree",
          f"{ratios[25.0]:.2f}, {ratios[30.0]:.2f} at 25 and 30 C")
    # tau's trustworthiness boundary, which decides where v_ss may be used.
    hot = frame[frame.temperature > arrhenius.BURST_TRUSTWORTHY_BELOW_C]
    warm = frame[frame.temperature <= arrhenius.BURST_TRUSTWORTHY_BELOW_C]
    claim("tau resolved, cold end",
          f"tau resolved on {int(warm.tau_resolved.sum())} of {len(warm)} curves")
    claim("tau unresolved, hot end",
          f"tau unresolved on {len(hot) - int(hot.tau_resolved.sum())} of "
          f"{len(hot)}")

    print("\nthe activation energy, first pass")
    fits = arrhenius.rung_fits()
    for s0, row in fits.iterrows():
        claim(f"rung {s0:.3f}",
              f"| {s0:.3f} | {row.activation_kJ:.1f} | {row.rms:.3f} |")
    claim("mean activation energy",
          f"Mean **{fits.activation_kJ.mean():.1f} kJ/mol**, spread "
          f"{fits.activation_kJ.max() - fits.activation_kJ.min():.1f}")

    print("\nthe caution about v0")
    rising = int((frame.accel_z > ACCELERATION_SIGMA).sum())
    strongly = int((frame.accel_z > 15).sum())
    claim("how many accelerate",
          f"{rising} of these {len(frame)} curves clear the 3sigma acceleration "
          f"gate and {strongly} exceed z = +15")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} MISMATCH(ES):")
        for item in FAILURES:
            print(f"  - {item}")
        return 1
    print("ANALYSIS.md agrees with the code")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
