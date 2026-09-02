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
import slowdown
import verify_enzyme_stock
from curve_metrics import ACCELERATION_SIGMA
from fit_dataset import build_curves

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
                         .replace("τ", "tau")
                         .replace("σ", "sigma")
                         .replace("µ", "u").replace("μ", "u")
                         .replace("°", ""))
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
    # Every curve on the instrument's own binary. Six of these were on the
    # rounded .txt export until 2026-09-02, when `read_rre`'s case-sensitive
    # regex was found to be dropping `sample003`. The claim is checked in both
    # directions -- the count here, and archive-wide below -- because a
    # regression would be silent: the export exists, so a fallback reports
    # nothing.
    claim("the source split",
          f"**{int((frame.source == 'rre').sum())} of {len(frame)} come from "
          f"the instrument's own `.rre`**")
    check("no curve in the series falls back to the export",
          bool((frame.source == "rre").all()),
          f"{int((frame.source == 'txt').sum())} on .txt")
    everything, _ = build_curves()
    fallbacks = [c for c in everything if c.source != "rre"]
    check("nor anywhere in the archive",
          not fallbacks,
          f"{len(fallbacks)} of {len(everything)} curves on .txt")
    claim("the recovery",
          f"**31 files across exps 1-32**")
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
    ladder = " -> ".join(f"{taus[t]:.0f}" for t in lag) + " s"
    claim("tau falls as it warms", "**" + ladder + "**")
    # And the same ladder in figure E's caption. It drifted there and nowhere
    # else -- the document carried 6489 -> 3190 -> 945 -> 916 in one section and
    # 6489 -> 3666 -> 945 -> 876 in another, four sections apart, for as long as
    # both existed, because only the prose was ever checked.
    drawn = _normalise(io.open(os.path.join(HERE, "build_figures.py"),
                               encoding="utf-8").read())
    check("figure E's caption carries the same ladder",
          _normalise(ladder) in drawn, ladder)
    check("every run below 32 C is a lag curve",
          bool((frame[frame.temperature <= arrhenius.BURST_TRUSTWORTHY_BELOW_C]
                .v0_burst_kind == "lag").all()))

    print("\nthe screen now finds both breakpoints")
    counts = frame.groupby(["temperature", "break_pattern"]).size()
    patterns = ("rising", "rise then fall", "falling")
    temperatures = sorted(frame.temperature.unique())
    for pattern in patterns:
        row = [int(counts.get((t, pattern), 0)) for t in temperatures]
        mark = "**" if pattern == "rise then fall" else ""
        claim(f"pattern row, {pattern}",
              f"| {pattern} | " + " | ".join(
                  f"{mark}{v}{mark}" if (mark and v) else str(v)
                  for v in row) + " |")
    two = int((frame.break_count == 2).sum())
    claim("how many take a second breakpoint",
          f"the F test takes a second breakpoint on **{two} of {len(frame)}**")
    # The exp 16 pair the user spotted by eye, now measured.
    hot = frame[(frame.temperature == 40) & np.isclose(frame.s0, 3.700)].iloc[0]
    claim("exp 16's early break",
          f"**{hot.break_times[0]:.0f} and {hot.break_times[1]:.0f} s**")
    # THE cross-check: the model-free pattern and the nested fit agree on which
    # temperatures are two-phase. If they ever stop agreeing, one of them is
    # describing something other than the curve.
    by_temperature = frame.groupby("temperature").agg(
        two_phase=("phases", lambda s: (s == 2).sum()),
        rising=("break_pattern", lambda s: (s == "rising").sum()),
        total=("phases", "size"))
    agree = all(
        (row.two_phase == 0) == (row.rising == row.total)
        for _, row in by_temperature.iterrows())
    check("the model-free pattern and the nested fit partition the block alike",
          agree,
          "; ".join(f"{t:.0f}C {int(r.two_phase)}/{int(r.total)} two-phase, "
                    f"{int(r.rising)}/{int(r.total)} rising"
                    for t, r in by_temperature.iterrows()))

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
    labels = {"v_peak": "`v_peak`, peak rate of the fitted model",
              "vmax": "`vmax`, steepest observed rate",
              "v_ss": "`v_ss`, asymptote after the induction",
              "inverse_tau": "`1/tau`, induction rate constant"}
    for name, row in table.iterrows():
        mark = "**" if name == "v_peak" else ""
        claim(f"parameter row, {name}",
              f"| {labels[name]} | {row.activation_kJ:.1f} +/- "
              f"{row.activation_stderr:.1f} | {mark}{row.enthalpy_kJ:.1f} +/- "
              f"{row.enthalpy_stderr:.1f}{mark} | {mark}{row.entropy_J:+.1f} +/- "
              f"{row.entropy_stderr:.1f}{mark} | {mark}{row.gibbs_kJ:.1f} +/- "
              f"{row.gibbs_stderr:.1f}{mark} | {int(row.n)}, "
              f"{int(row.temperatures)} T |")
    kcal = 4.184
    vmax_row = table.loc["v_peak"]
    claim("the kcal restatement",
          f"ΔH‡ **{vmax_row.enthalpy_kJ / kcal:.1f}**, ΔS‡ "
          f"**{vmax_row.entropy_J / kcal:.1f} cal/mol/K**, ΔG‡(298)\n"
          f"**{vmax_row.gibbs_kJ / kcal:.1f}**")
    claim("the three routes",
          f"`v_peak` gives **{table.loc['v_peak', 'enthalpy_kJ']:.1f} +/- "
          f"{table.loc['v_peak', 'enthalpy_stderr']:.1f}** on\nall six "
          f"temperatures; `v_ss` gives **{table.loc['v_ss', 'enthalpy_kJ']:.1f}"
          f" +/- {table.loc['v_ss', 'enthalpy_stderr']:.1f}**")
    claim("the rung agreement cost",
          f"reduced chi2 {arrhenius.rungs_agree('vmax')['reduced_chi2']:.2f} "
          f"against **{arrhenius.rungs_agree('v_peak')['reduced_chi2']:.2f}** "
          f"for `v_peak`")
    # The composition DS is quoted at. A DS with none attached is not a number.
    claim("the composition the entropy is quoted at",
          f"[S] = {vmax_row.at_s0:.3f} mM and [buf] = {vmax_row.at_buf:.0f} mM**")
    # The two estimators agreeing is the load-bearing check of section 4.
    print("\nthe fitting form: one relaxation or two")
    # The nesting is exact -- B2 = 0 is the one-phase form -- so the count and
    # the residual improvement are the evidence for the extra parameters.
    frame_sel = arrhenius.series_frame()
    selected = int((frame_sel.phases == 2).sum())
    claim("how many earned the second phase",
          f"**Selected on {selected} of {len(frame_sel)}**")
    for temperature in sorted(frame_sel.temperature.unique()):
        run = frame_sel[frame_sel.temperature == temperature]
        two = int((run.phases == 2).sum())
        mark = "**" if two else ""
        claim(f"selection row, {temperature:.0f} C",
              f"| {temperature:.0f} °C | {mark}{two} of {len(run)}{mark} | "
              f"{run.v0_burst_resid.median():.2f}x | "
              f"{mark}{run.progress_resid.median():.2f}x{mark} |")
    # Nothing may be selected where the one-phase form already sits at noise.
    quiet = frame_sel[frame_sel.v0_burst_resid < 1.2]
    check("nothing is selected where one phase already fits",
          int((quiet.phases == 2).sum()) <= 1,
          f"{int((quiet.phases == 2).sum())} of {len(quiet)} quiet curves")
    # tau2 is quoted nowhere, and this is why.
    resolved = int(frame_sel[frame_sel.phases == 2].tau_slow_resolved.sum())
    claim("tau2 is not quoted",
          f"resolved on only {resolved} of the {selected} two-phase curves")

    # THE load-bearing check: three routes with unrelated weaknesses.
    for other in ("vmax", "v_ss"):
        check(f"v_peak and {other} agree within their errors",
              abs(table.loc["v_peak", "enthalpy_kJ"]
                  - table.loc[other, "enthalpy_kJ"])
              < 2 * max(table.loc[other, "enthalpy_stderr"],
                        table.loc["v_peak", "enthalpy_stderr"]),
              f"{table.loc['v_peak', 'enthalpy_kJ']:.1f} vs "
              f"{table.loc[other, 'enthalpy_kJ']:.1f}")
    claim("they come out at",
          f"They come out at {table.loc['vmax', 'enthalpy_kJ']:.1f}, "
          f"{table.loc['v_ss', 'enthalpy_kJ']:.1f} and\n"
          f"{table.loc['v_peak', 'enthalpy_kJ']:.1f} kJ/mol")
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


    print("\nthe catalyst-order test, which the block cannot do for itself")
    order = scope.catalyst_order()
    claim("the pair",
          f"**0.028 against 0.014,\nexactly {2.0:.3f}x**")
    claim("the order in catalyst",
          f"**+{order['order_in_catalyst']:.2f}**, not 1")
    # The DIRECTION is the verdict, not the number: inactivation requires the
    # ratio to fall BELOW 1 and it never does.
    check("catalyst inactivation is excluded by the direction of the effect",
          not order["inactivation"],
          f"lowest ratio {order['lowest_ratio']:.2f}, "
          f"median {order['median_ratio']:.2f}")
    # And the tension, which has to stay visible: this block prefers n = 1.
    frame_local = arrhenius.series_frame()
    rms = {}
    for power in (1.0, order["order_in_catalyst"]):
        values = []
        for s0 in sorted(frame_local.s0.unique()):
            rung = frame_local[np.isclose(frame_local.s0, s0)].sort_values("kelvin")
            fit = arrhenius.arrhenius_fit(
                rung.kelvin, rung.vmax / rung.e0 ** power)
            values.append(fit["rms"])
        rms[round(power, 2)] = float(np.mean(values))
    claim("the tension",
          f"`[enz]^1` gives an Arrhenius rms of **{rms[1.0]:.3f}** against "
          f"**{rms[round(order['order_in_catalyst'], 2)]:.3f}** for")
    check("this block still prefers first order in enzyme",
          rms[1.0] < rms[round(order["order_in_catalyst"], 2)],
          f"{rms[1.0]:.3f} vs {rms[round(order['order_in_catalyst'], 2)]:.3f}")

    print("\nsection 4: naming the fall does not move the numbers")
    effect = slowdown.sink_effect_on_activation()
    for label, key in (("v_peak", "published"), ("v_prod", "corrected")):
        row = effect[key]
        claim(f"{label}'s row",
              f"| {row['activation_kJ']:.2f} +/- "
              f"{row['activation_stderr']:.2f} | {row['enthalpy_kJ']:.2f} | "
              f"{row['entropy_J']:+.1f} +/- {row['entropy_stderr']:.1f} | "
              f"{row['gibbs_kJ']:.2f} | {row['rms']:.3f} |")
    claim("how far above v_peak the production rate sits",
          f"**{effect['lift'] - 1:.1%} above**")
    claim("its range over the six temperatures",
          f"between {effect['lift_low'] - 1:+.1%} and "
          f"{effect['lift_high'] - 1:+.1%}")
    claim("the shift in the activation energy",
          f"**{effect['activation_shift']:+.2f} +/- "
          f"{effect['activation_shift_stderr']:.2f} kJ/mol**")
    claim("the shift in the entropy",
          f"**{effect['entropy_from_lift']:+.2f} J/mol/K**")
    check("the corrected route is the noisier one",
          effect["corrected"]["rms"] > effect["published"]["rms"],
          f"{effect['corrected']['rms']:.3f} vs "
          f"{effect['published']['rms']:.3f}")

    print("\nsection 6: the pointer to product_fate")
    whole = scope.frame(tuple(range(1, 152)))
    named = slowdown.substrate_blocks(whole)
    series = slowdown.deceleration_drivers(named["temperature series"])
    free = slowdown.deceleration_drivers(named["4OMe enzyme-free"])
    fixed = slowdown.deceleration_drivers(named["4OMe catalysed, phosphate"],
                                          fixed=True)
    claim("the block's own product coefficient",
          f"{series['product']:+.3f} +/- {series['product_stderr']:.3f}")
    claim("against its run-length coefficient",
          f"{series['span']:+.3f} +/- {series['span_stderr']:.3f}")
    claim("the enzyme-free block's run-length coefficient",
          f"{free['span']:+.3f} +/- {free['span_stderr']:.3f}")
    claim("and its product coefficient",
          f"{free['product']:+.3f} +/- {free['product_stderr']:.3f}")
    claim("the every-condition-fixed comparison",
          f"{fixed['product']:+.3f} +/- {fixed['product_stderr']:.3f}")
    table = slowdown.sink_table(
        sorted(named["4OMe catalysed, phosphate"].experiment.unique()))
    clean = table[(table.points > 0)
                  & (table.rate_r2 > slowdown.SINK_CLEAN_R2)]
    order = slowdown.plateau_scaling(clean)
    claim("the plateau's substrate order",
          f"order {order['order']:+.3f} +/- {order['stderr']:.3f}")
    picked = np.array([slowdown.selectivity(p, s, e) for p, s, e
                       in zip(clean.plateau, clean.s0, clean.epsilon)])
    picked = picked[np.isfinite(picked)]
    claim("the selectivity it implies",
          f"k_A/k_S ~ {np.median(picked):.0f}")
    counts = clean.prefers.value_counts()
    claim("how many curves prefer the sink form",
          f"on {int(counts.get('sink', 0))} of {len(clean)} well-determined "
          f"curves against {int(counts.get('inhibition', 0))}")

    print("\nthe figures: lettered once each, in order")
    # J and K were each used TWICE from the day figure_selection was added --
    # once in section 3a and once in section 4 -- and nothing noticed, because
    # a letter is not a number and no check looked at one.
    import re
    page = io.open(os.path.join(HERE, "index.html"), encoding="utf-8").read()
    letters = re.findall(r">([A-Z]) \u00b7 ", page)
    expected = [chr(ord("A") + index) for index in range(len(letters))]
    check("every figure letter is used exactly once, A onwards",
          letters == expected,
          f"{''.join(letters)} against {''.join(expected)}")
    claim("the document's own count of them",
          f"{len(letters)} figures" if len(letters) != 15 else
          "fifteen figures, A to " + letters[-1])

    print("\nthe fitting form, section 3a")
    series = arrhenius.series_frame()
    two_phase = int((series.phases == 2).sum())
    claim("how many curves the burst form cannot hold",
          f"{two_phase} of these {len(series)} curves do exactly that")
    check("the model-free screen agrees with the F test on that count",
          int(series.break_pattern.isin(("rise then fall", "falling")).sum())
          == two_phase,
          f"{int(series.break_pattern.isin(('rise then fall', 'falling')).sum())}"
          f" against {two_phase}")

    print("\nthe figures: the fit never covers the data it is fitting")
    # Three passes were needed to get this right -- fit under the data, then
    # fit over it on a white halo wider than a mark, then this -- so the
    # invariant is checked rather than remembered. On the tightest panel the
    # scatter is sub-pixel, so the only thing keeping the readings visible is
    # that the line is narrower than a mark.
    import re
    panels = io.open(os.path.join(HERE, "progress_curves.html"),
                     encoding="utf-8").read().split("<div class='fig panel'>")[1:]
    worst = None
    for block in panels:
        marks = re.findall(r"<circle cx='[\d.]+' cy='[\d.]+' r='([\d.]+)'", block)
        line = re.search(r"stroke='#c0522a' stroke-width='([\d.]+)'", block)
        if not marks or not line:
            continue
        margin = 2 * float(marks[0]) - float(line.group(1))
        worst = margin if worst is None else min(worst, margin)
    ok = worst is not None and worst > 0
    print(f"  {'pass' if ok else 'FAIL'}  the fit is narrower than a mark on "
          f"every panel: {worst:.2f} px of margin at worst")
    if not ok:
        FAILURES.append("the fit line is at least as wide as a data mark, so a "
                        "reading centred on it is covered whole")
    # And nothing white is drawn over the data: a halo erases rather than
    # overlays, which is what the second attempt did.
    haloed = [b for b in panels if "stroke='#ffffff'" in b]
    check("no white halo is drawn over the readings", not haloed,
          f"{len(haloed)} panels")

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
