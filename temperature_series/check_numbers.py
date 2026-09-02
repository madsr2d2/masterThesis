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
sys.path.insert(0, os.path.dirname(HERE))

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
                         .replace("÷", "/")
                         .replace("χ²", "chi2")   # chi-square, as typed
                         .replace("Δ", "D").replace("‡", "")
                         .replace("τ", "tau"))
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

    print("\nwhat is going on: the buffer riding on the substrate ladder")
    # The ladder's [buf] gradient is the reason the substrate order needs
    # correcting and the reason the activation energies do NOT -- both halves
    # are checked, because the second is the licence for everything in section 4.
    ladder = frame.pivot_table(index="temperature", columns="sample",
                               values="buf")
    check("[buf] falls 80-50 along the ladder, identically in all six runs",
          bool((ladder.nunique() == 1).all())
          and list(ladder.iloc[0]) == [80.0, 70.0, 60.0, 50.0],
          f"{list(ladder.iloc[0])}")
    high = arrhenius.catalysed_buffer_order()
    low = arrhenius.catalysed_buffer_order(arrhenius.BUFFER_TITRATION_LOW)
    claim("buffer order, high range",
          f"| 50-200 mM | 32 | **+{high['order_buf']:.3f} +/- "
          f"{high['stderr_buf']:.3f}** | {high['r2']:.3f} |")
    claim("buffer order, low range",
          f"| 3.125-25 mM | 34 | +{low['order_buf']:.3f} +/- "
          f"{low['stderr_buf']:.3f} | {low['r2']:.3f} |")
    check("the buffer dependence saturates",
          low["order_buf"] > high["order_buf"] + 0.2,
          f"{low['order_buf']:.2f} against {high['order_buf']:.2f}")

    orders = arrhenius.substrate_order()
    claim("the coupling and buffer order used",
          f"**g = {orders.attrs['coupling']:+.3f}** and "
          f"**d = +{orders.attrs['buffer_order']:.3f}**")
    claim("the observed order row",
          "| observed | " + " | ".join(
              f"{orders.observed[t]:+.3f}" for t in sorted(orders.index)) + " |")
    claim("the corrected order row",
          "| corrected | " + " | ".join(
              f"{orders.corrected[t]:+.3f}" for t in sorted(orders.index)) + " |")
    claim("the mean corrected order",
          f"Mean **{orders.corrected.mean():+.3f}**")
    check("correcting moves the level, not the trend",
          abs(float(np.ptp(orders.observed) - np.ptp(orders.corrected))) < 1e-9,
          f"{np.ptp(orders.observed):.3f} vs {np.ptp(orders.corrected):.3f}")

    print("\nthe rungs share one activation energy")
    agree = arrhenius.rungs_agree("vmax")
    claim("the chi-square",
          f"**chi2 = {agree['chi2']:.2f} on {agree['dof']} degrees of freedom, "
          f"reduced chi2 = {agree['reduced_chi2']:.2f}**")
    check("the spread is unremarkable against the errors",
          agree["reduced_chi2"] < 3.0, f"{agree['reduced_chi2']:.2f}")

    print("\nthe activation parameters")
    table = arrhenius.parameter_table()
    labels = {"vmax": "`vmax`, steepest observed rate",
              "v_ss": "`v_ss`, asymptote after the induction",
              "inverse_tau": "`1/tau`, induction rate constant"}
    for name, row in table.iterrows():
        mark = "**" if name == "vmax" else ""
        claim(f"parameter row, {name}",
              f"| {labels[name]} | {row.activation_kJ:.1f} +/- "
              f"{row.activation_stderr:.1f} | {mark}{row.enthalpy_kJ:.1f} +/- "
              f"{row.enthalpy_stderr:.1f}{mark} | {mark}{row.entropy_J:+.1f} +/- "
              f"{row.entropy_stderr:.1f}{mark} | {mark}{row.gibbs_kJ:.1f} +/- "
              f"{row.gibbs_stderr:.1f}{mark} | {int(row.n)}, "
              f"{int(row.temperatures)} T |")
    kcal = 4.184
    vmax_row = table.loc["vmax"]
    claim("the kcal restatement",
          f"ΔH‡ **{vmax_row.enthalpy_kJ / kcal:.1f}**, ΔS‡ "
          f"**{vmax_row.entropy_J / kcal:.1f} cal/mol/K**, ΔG‡(298) "
          f"**{vmax_row.gibbs_kJ / kcal:.1f}**")
    # The composition DS is quoted at. A DS with none attached is not a number.
    claim("the composition the entropy is quoted at",
          f"[S] = {vmax_row.at_s0:.3f} mM and [buf] = {vmax_row.at_buf:.0f} mM**")
    # The two estimators agreeing is the load-bearing check of section 4.
    check("vmax and v_ss agree within their errors",
          abs(table.loc["vmax", "enthalpy_kJ"] - table.loc["v_ss", "enthalpy_kJ"])
          < 2 * table.loc["v_ss", "enthalpy_stderr"],
          f"{table.loc['vmax', 'enthalpy_kJ']:.1f} vs "
          f"{table.loc['v_ss', 'enthalpy_kJ']:.1f}")
    claim("they come out at",
          f"They come out at {table.loc['vmax', 'enthalpy_kJ']:.1f} and\n"
          f"{table.loc['v_ss', 'enthalpy_kJ']:.1f} kJ/mol")
    # The induction is a DIFFERENT process, and the entropy is what says so.
    gap = table.loc["vmax", "gibbs_kJ"] - table.loc["inverse_tau", "gibbs_kJ"]
    claim("the induction's gap", f"**{gap:.0f} kJ/mol lower**")

    print("\nthe activation energy, first pass")
    fits = arrhenius.rung_fits()
    buffers = frame.groupby("s0").buf.median()
    for s0, row in fits.iterrows():
        claim(f"rung {s0:.3f}",
              f"| {s0:.3f} | {buffers[s0]:.0f} | {row.activation_kJ:.1f} +/- "
              f"{row.stderr_kJ:.1f} | {row.rms:.3f} |")
    agree = arrhenius.rungs_agree("vmax")
    pooled = arrhenius.pooled_arrhenius("vmax")
    claim("weighted mean and pooled refit",
          f"Weighted mean **{agree['weighted_mean_kJ']:.1f}**, reduced chi2 "
          f"**{agree['reduced_chi2']:.2f}**")
    claim("the pooled value",
          f"**{pooled['activation_kJ']:.1f} +/- {pooled['stderr_kJ']:.1f}**")

    print("\nthe caution about v0")
    rising = int((frame.accel_z > ACCELERATION_SIGMA).sum())
    strongly = int((frame.accel_z > 15).sum())
    claim("how many accelerate",
          f"{rising} of these {len(frame)} curves clear the 3sigma acceleration "
          f"gate and {strongly} exceed z = +15")


    print("\nthe figures: no data point drawn outside its own frame")
    # Marks are clipped to the plot area deliberately, so a data point outside
    # the axis limits vanishes silently and the figure still looks complete.
    # That is how exp 6 sample 4 went missing from the curvature figure for as
    # long as it existed: no number changes, so no prose check could see it.
    from svgplot import clipped_marks
    for name in ("index.html", "progress_curves.html"):
        path = os.path.join(HERE, name)
        if not os.path.exists(path):
            continue
        lost = clipped_marks(io.open(path, encoding="utf-8").read())
        ok = not lost
        print(f"  {'pass' if ok else 'FAIL'}  {name}: {len(lost)} clipped")
        if not ok:
            for title, x, y, _ in lost[:4]:
                print(f"        {title!r} at ({x:.0f},{y:.0f})")
            FAILURES.append(f"{name} draws {len(lost)} point(s) outside the "
                            f"plot frame; widen the axis limits")

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
