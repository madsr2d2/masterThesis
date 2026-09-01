"""
Verifies every number quoted in ANALYSIS.md against the code that produces it.

The point is drift. A figure typed into a document is a copy, and copies in this
project have gone stale before -- the whole reason `curve_metrics` exists. So
each claim below is re-derived from the modules AND its rendered form is
required to appear verbatim in ANALYSIS.md. Change the code and this fails;
edit the prose inconsistently and this fails.

    python background_reaction/check_numbers.py
"""
import io
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))

import scope
from curve_metrics import (ACCELERATION_SIGMA, OUTLIER_SIGMA, acceleration,
                           isolated_outliers, local_outlier_z)
from fit_dataset import build_curves, source_floor
from summary_kinetics import fit_burst, fit_burst_bounded

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")
# The RENDERED page, checked as well as the prose. The legend claimed the
# burst fit was blanket "B <= 0 ... 0 of 16" for three commits after the gate
# became per-curve, and nothing noticed, because every check here read
# ANALYSIS.md and the stale sentence lived in build_figures.py.
CURVES_PAGE = os.path.join(HERE, "progress_curves.html")
# MECHANISM.md carries the same boric-probe numbers, in ASCII. A figure quoted
# in two documents is two copies, and this project has had a number go stale in
# one file while the other stayed right.
MECHANISM_DOC = os.path.join(os.path.dirname(HERE), "MECHANISM.md")
ESTIMATORS_ALL = ("v0_quad", "v0_burst", "vmax", "v0", "v0_whole")
FAILURES = []


def _normalise(text):
    """
    Fold the typographic characters prose uses onto the ones %-formatting emits.

    U+2212 MINUS SIGN reads better than a hyphen in a document and is what the
    tables are written with, but `f"{value:+.3f}"` produces U+002D. Comparing
    glyphs would make this checker fail on typography while a genuinely wrong
    number that happened to be typed with a hyphen went through.

    Runs of whitespace are collapsed for the same reason: a claim must not fail
    because the prose rewrapped or because a list item is indented. That is
    typography, not a number.
    """
    return " ".join((text.replace("\u2212", "-")     # minus sign
                         .replace("\u2013", "-")     # en dash
                         .replace("\u00b1", "+/-")   # plus-minus
                         .replace("\u00d7", "x")     # multiplication sign
                         .replace("\u2080", "0")     # subscript zero, as in v0
                         .replace("\u207b", "-")     # superscript minus, [HOO-]
                         .replace("\u209b\u209b", "ss"))  # v_ss
                    .split())


def claim(label, rendered, present=True, document=None):
    """Require `rendered` to appear in the document."""
    path = document or DOCUMENT
    text = _normalise(io.open(path, encoding="utf-8").read())
    ok = (_normalise(rendered) in text) == present
    print(f"  {'pass' if ok else 'FAIL'}  {label}: {rendered!r}")
    if not ok:
        FAILURES.append(f"{label} -- {rendered!r} not found in "
                        f"{os.path.basename(path)}")


def main():
    print("checking ANALYSIS.md against the code\n")

    print("the enzyme-free census")
    import pandas as pd
    data = pd.read_csv(os.path.join(os.path.dirname(HERE), "data",
                                    "experiment_data.csv"))
    free = data[data["[enz]"] == 0]
    claim("enzyme-free curves", f"{len(free)} enzyme-free curves over "
                                f"{free.experiment.nunique()} experiments")
    fixed = scope.frame(scope.BUFFER_FIXED)
    titration = scope.frame(scope.BUFFER_CONFOUNDED)
    claim("BUFFER_FIXED size", f"| 65, 67, 69, 70 | {len(fixed)} |")
    claim("BUFFER_CONFOUNDED size",
          f"| 3, 6 | {len(titration)} ({int(titration.live.sum())} live) |")

    print("\nthe separation")
    quad = scope.buffer_dependence(parameter="v0_quad")
    vmax = scope.buffer_dependence(parameter="vmax")
    claim("a, v0_quad", f"**{quad['order_s0_fixed']:+.3f} ± "
                        f"{quad['stderr_s0_fixed']:.3f}** (n = {quad['n_fixed']})")
    claim("a', v0_quad", f"**{quad['order_s0_titration']:+.3f} ± "
                         f"{quad['stderr_s0_titration']:.3f}** "
                         f"(n = {quad['n_titration']})")
    claim("a, vmax", f"{vmax['order_s0_fixed']:+.3f} ± "
                     f"{vmax['stderr_s0_fixed']:.3f} (n = {vmax['n_fixed']})")
    claim("coupling g", f"g = {quad['coupling']:.3f}")
    claim("d, v0_quad", f"**d = {quad['order_buf']:+.2f} ± "
                        f"{quad['stderr_buf']:.2f}** on `v0_quad`")

    print("\nthe result table")
    for estimator, marker in (("v0_quad", "**"), ("vmax", ""),
                              ("v0", ""), ("v0_whole", "")):
        result = scope.buffer_dependence(parameter=estimator)
        cross = scope.background_orders(scope.FREE_4OME_40C,
                                        terms=("s0", "buf"), parameter=estimator)
        claim(f"{estimator} row",
              f"| {marker}{result['order_buf']:+.2f} ± "
              f"{result['stderr_buf']:.2f}{marker} | "
              f"{marker}{cross['order_buf']:+.2f} ± "
              f"{cross['stderr_buf']:.2f}{marker} |")

    print("\nthe joint fit that cannot carry it")
    joint = scope.background_orders(scope.BUFFER_CONFOUNDED,
                                    terms=("s0", "buf"), within=True,
                                    parameter="v0_quad")
    # The live count is derived, not spelled: it went 8 -> 10 on 2026-09-01
    # when exps 3 and 6 moved onto the instrument's own readings and two more
    # of exp 3's cuvettes cleared the live-signal threshold.
    spelled = {8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}
    claim("joint VIF", f"VIF {joint['vif_s0']:.1f} and {joint['vif_buf']:.1f} "
                       f"on {spelled.get(joint['n'], joint['n'])} live curves; "
                       f"returns {joint['order_buf']:+.2f} ± "
                       f"{joint['stderr_buf']:.2f}")

    print("\nthe burst diagnostics")
    curves = scope.curves(scope.FREE_BNOH_ALL)
    plain = [fit_burst(c.times, c.absorbance) for c in curves]
    bounded_free = [fit_burst_bounded(c.times, c.absorbance, constrain=False,
                                      noise_floor=source_floor(c.source))
                    for c in curves]
    bounded = [fit_burst_bounded(c.times, c.absorbance,
                                 noise_floor=source_floor(c.source))
               for c in curves]
    negative = sum(1 for f in plain if np.isfinite(f.v0) and f.v0 < 0)
    claim("negative v0 count",
          f"**{negative} of {len(curves)} curves return a negative v₀**")
    always = [fit_burst_bounded(c.times, c.absorbance, constrain=True,
                                noise_floor=source_floor(c.source))
              for c in curves]
    claim("bounded count, lag branch shut everywhere",
          f"from **{sum(f.bounded for f in bounded_free)} to "
          f"{sum(f.bounded for f in always)} of {len(curves)}**")
    # `bounded` is the DEFAULT, constrain="auto": the branch is gated on each
    # curve's own acceleration rather than shut for the whole block.
    claim("bounded count, lag branch gated per curve",
          f"That bounds **{sum(f.bounded for f in bounded)} of {len(curves)}**")

    print("\nshape and scope")
    live_fixed = fixed[fixed.live]
    inscope = scope.frame()
    live_inscope = inscope[inscope.live]
    claim("background acceleration",
          f"| **{int((live_fixed.accel_z > ACCELERATION_SIGMA).sum())} of "
          f"{len(live_fixed)}** |")
    claim("in-scope acceleration",
          f"| **{int((live_inscope.accel_z > ACCELERATION_SIGMA).sum())} of "
          f"{len(live_inscope)} live** |")
    buffers = sorted(inscope.buf.round(3).unique())
    claim("in-scope buffer constant",
          f"`[buf]` = {buffers[0]:.3f} mM in all {len(inscope)} in-scope curves")
    everything = scope.frame(scope.FREE_BNOH_ALL)
    curved = int((everything.curvature_t.abs() > 3).sum())
    claim("curvature count", f"{curved}\nof {len(everything)} curves show "
                             f"curvature at |t| > 3")

    print("\nthe estimator description table")
    # This table carries no fitted numbers, so nothing else here would notice
    # if it were clobbered -- and a bulk regex did clobber it in c41f459..
    # 4b344ac, replacing the descriptions with buffer orders for three
    # commits. Pinned by content now.
    for phrase in ("| `vmax` | steepest of four 20% blocks | no | straight line |",
                   "| `v0` | first 20% | no | straight line |",
                   "| `v0_whole` | none | yes | straight line |"):
        claim("description row", phrase)

    print("\nthe accelerating-curve sensitivity")
    for estimator in ESTIMATORS_ALL:
        dropped = scope.buffer_dependence(parameter=estimator,
                                          drop_accelerating=True)
        claim(f"{estimator} dropped-accelerating",
              f"| {dropped['order_buf']:+.2f} ± {dropped['stderr_buf']:.2f} |")

    print("\nthe in-scope substrate order")
    orders = scope.orders("vmax", scope.PRIMARY_SCOPE, within=True)
    strong = scope.orders("vmax", scope.strong_runs(), within=True)
    claim("in-scope order", f"{orders['order_s0']:+.3f} ± "
                            f"{orders['stderr_s0']:.3f} over all "
                            f"{orders['n']} live in-scope curves")
    claim("strong-run order", f"({strong['order_s0']:+.3f} ± "
                              f"{strong['stderr_s0']:.3f} on the "
                              f"{len(scope.strong_runs())} strong runs)")

    print("\nthe rate law, both scopes")
    # Pinned at last. The line in section 6 was the PRE-first-reading-drop
    # value from 2026-09-01 to b144b9e -- [H2O2]^+1.78 where the code said
    # +1.87, [HOO-]^+0.68 where it said +0.63 -- because no check re-derived
    # it. Same failure as the legend and the estimator table before it.
    for estimator in ("v0_quad", "vmax"):
        for label, block, anchor in (
                ("all", scope.FREE_BNOH_ALL, scope.BUFFER_FIXED),
                ("phosphate", scope.FREE_BNOH_PHOSPHATE, (67, 69, 70))):
            law = scope.background_orders(block, terms=("s0", "h2o2", "hoo"),
                                          parameter=estimator)
            buf = scope.buffer_dependence(anchor=anchor, parameter=estimator)
            claim(f"rate law, {estimator}, {label}",
                  f"[S]^{buf['order_s0_fixed']:+.2f}  "
                  f"[H2O2]^{law['order_h2o2']:+.2f}  "
                  f"[HOO-]^{law['order_hoo']:+.2f}  "
                  f"[buf]^{buf['order_buf']:+.2f}     ({estimator})")
            claim(f"pooled R2, {estimator}, {label}",
                  f"{law['r2']:.3f}")

    print("\nthe species/total degeneracy behind the buffer order")
    # Section 6b claims the design cannot separate a buffer SPECIES from the
    # total. That claim is the reason the peroxophosphate hypothesis is left
    # open rather than tested, so it is checked rather than asserted.
    live = scope.frame(scope.FREE_BNOH_PHOSPHATE)
    live = live[live.live]
    titration = live[live.experiment.isin(scope.BUFFER_CONFOUNDED)]
    logs = {name: np.log(titration[name].to_numpy())
            for name in ("buf", "buf_acid", "buf_base")}
    worst = min(float(np.corrcoef(logs["buf"], logs[name])[0, 1])
                for name in ("buf_acid", "buf_base"))
    claim("species and total are one variable in the titration",
          f"correlation **{worst:.6f}**")
    inflation_total = scope.variance_inflation(
        live, ("s0", "h2o2", "hoo", "buf"))["buf"]
    inflation_base = scope.variance_inflation(
        live, ("s0", "h2o2", "hoo", "buf_base"))["buf_base"]
    claim("substituting the basic form inflates the variance",
          f"factor to **{inflation_base:.1f}**, against **{inflation_total:.1f}** for the")
    # The two pH levels, and that there are only two of them.
    levels = sorted(live.pH.unique())
    ok = len(levels) == 2
    print(f"  {'pass' if ok else 'FAIL'}  the phosphate set has two pH levels: "
          f"{levels}")
    if not ok:
        FAILURES.append("the phosphate set no longer has exactly two pH levels; "
                        "section 6b's degeneracy argument needs rewriting")
    span = live.groupby("pH").agg(buf=("buf", "median"),
                                  acid=("buf_acid", "median"),
                                  base=("buf_base", "median"),
                                  hoo=("hoo", "median"))
    claim("what moves between the two pH levels",
          f"`[buf]` {span.buf.iloc[0]:.0f} -> {span.buf.iloc[1]:.0f} mM, "
          f"`[H2PO4-]` {span.acid.iloc[0]:.1f} -> {span.acid.iloc[1]:.1f}, "
          f"`[HPO4^2-]` {span.base.iloc[0]:.1f} -> {span.base.iloc[1]:.1f}, "
          f"`[HOO-]` {span.hoo.iloc[0]:.4f} -> {span.hoo.iloc[1]:.3f}")

    print("\nthe boric probe: does a peroxo-forming buffer run fast?")
    # Section 6b's answer to "why first order in buffer" now rests on this
    # table, in ANALYSIS.md AND in MECHANISM.md. Both are checked.
    peroxo = scope.peroxo_buffer_test()
    for estimator, row in peroxo.iterrows():
        note = " (2 of 4 cuvettes)" if int(row["n"]) < 4 else ""
        claim(f"boric probe, {estimator}",
              f"| `{estimator}` | {row['predicted']:.2f}x | "
              f"{row['observed']:.2f}x | **{row['excess']:.2f}x**{note} |")
        claim(f"boric probe in MECHANISM.md, {estimator}",
              f"{row['excess']:.2f}x ({estimator})", document=MECHANISM_DOC)
    # The claim that three of four fall BELOW the law, which is the sentence
    # the conclusion is written on.
    below = int((peroxo.excess < 1).sum())
    ok = below == 3
    print(f"  {'pass' if ok else 'FAIL'}  three of four estimators put boric "
          f"below the phosphate law: {below} of {len(peroxo)}")
    if not ok:
        FAILURES.append("the boric probe no longer puts 3 of 4 estimators below "
                        "the phosphate law; section 6b's conclusion needs rewriting")
    # And that the pair really is matched, which is what lets the substrate and
    # peroxide orders cancel out of the prediction.
    pair = scope.frame(tuple(sorted(scope.PEROXO_PAIR)))
    pair = pair[pair.live]
    matched = (pair.groupby("experiment").h2o2.nunique().eq(1).all()
               and pair.h2o2.nunique() == 1
               and len(set(pair[pair.experiment == scope.PEROXO_PAIR[0]].s0)
                       ^ set(pair[pair.experiment == scope.PEROXO_PAIR[1]].s0)) == 0)
    print(f"  {'pass' if matched else 'FAIL'}  exps {scope.PEROXO_PAIR} are "
          f"matched in [S] and [H2O2]")
    if not matched:
        FAILURES.append("exps 65 and 67 no longer match in [S] and [H2O2]; the "
                        "boric probe's prediction needs the other two orders")

    print("\nthe mid-run break that weakens the boric probe")
    # Found from the plots, not by any statistic in the package -- every other
    # shape column compares a curve's start to its end and steps over a break
    # in the middle. Section 6b's downgrade of the boric probe rests on this
    # table, so it is re-derived rather than typed.
    block = tuple(sorted(set(scope.FREE_BNOH_ALL) | {59, 60, 61, 62}))
    breaks = scope.synchronised_break(block)
    exp65 = breaks.loc[scope.PEROXO_PAIR[0]]
    times = [str(v) for v in sorted(exp65["breaks"])]
    claim("exp 65's break times",
          f"{', '.join(times[:-1])} and {times[-1]} s")
    ratios = [f"{v:.2f}" for v in sorted(exp65["ratios"])]
    claim("exp 65's steepening",
          f"**{', '.join(ratios[:-1])} and {ratios[-1]}x**")
    claim("the break span", f"a span of **{exp65['span']:.0f} s**")
    # The claim the argument turns on: exp 65 is the ONLY run that steepens.
    others = breaks.drop(index=scope.PEROXO_PAIR[0])
    # "All of them" is the claim, not "any of them": exp 3 has one steep
    # cuvette of six, which is scatter, not a run-wide event.
    ok = (int(exp65["steep"]) == int(exp65["n"])
          and not (others.steep == others.n).any())
    print(f"  {'pass' if ok else 'FAIL'}  exp 65 is the only run whose cuvettes "
          f"all steepen: {int(exp65['steep'])} of {int(exp65['n'])}, "
          f"next best {int(others.steep.max())} of {int(others.n[others.steep.idxmax()])}")
    if not ok:
        FAILURES.append("exp 65 is no longer the only run that steepens after "
                        "its break; section 6b's downgrade needs rewriting")
    # And that no OTHER boric run does it, which is what rules the shape out as
    # borate chemistry.
    boric_others = breaks.loc[[59, 60, 61, 62]]
    ok = int(boric_others.steep.sum()) == 0
    print(f"  {'pass' if ok else 'FAIL'}  no other boric run steepens "
          f"(exps 59-62, {int(boric_others.n.sum())} curves): "
          f"{int(boric_others.steep.sum())}")
    if not ok:
        FAILURES.append("a boric run other than 65 now steepens; the shape may "
                        "be borate chemistry after all")

    print("\nthe second probe, which does not use exp 65")
    catalysed = scope.peroxo_buffer_test(pair=scope.CATALYSED_PEROXO_PAIR,
                                         orders_scope=None)
    for estimator, row in catalysed.iterrows():
        claim(f"catalysed probe, {estimator}",
              f"| `{estimator}` | {row['observed']:.2f}x |")
    # Quoted in MECHANISM.md as one line; checked as one string so it cannot
    # collide with the enzyme-free probe's numbers, which share a value.
    claim("catalysed probe in MECHANISM.md",
          "boric/phosphate = " + ", ".join(
              f"{row['observed']:.2f}x ({estimator})"
              for estimator, row in catalysed.iterrows()),
          document=MECHANISM_DOC)
    pair = scope.frame(scope.CATALYSED_PEROXO_PAIR)
    pair = pair[pair.live]
    hoo = float(pair[pair.experiment == scope.CATALYSED_PEROXO_PAIR[0]].hoo.median()
                / pair[pair.experiment == scope.CATALYSED_PEROXO_PAIR[1]].hoo.median())
    claim("the hydroperoxide it carries", f"**{hoo:.2f}x the [HOO-]**")

    print("\nthe buffer order is a phosphate number")
    # buffer_dependence's default anchor includes the boric run. Section 6b
    # says the order survives dropping it; that is checked, not asserted.
    for estimator, key in (("vmax", "vmax"), ("v0_quad", "v0_quad")):
        without = scope.buffer_dependence(
            anchor=tuple(e for e in scope.BUFFER_FIXED
                         if e not in scope.BORIC_BUFFER),
            parameter=estimator)
        claim(f"buffer order without boric, {estimator}",
              f"**+{without['order_buf']:.2f} +/- {without['stderr_buf']:.2f}**")

    print("\nwhat the boric run carries")
    spread = scope.boric_spread()
    for term, name in (("s0", "`[S]`"), ("h2o2", "`[H2O2]`"),
                       ("hoo", "`[HOO-]`"), ("buf", "`[buf]`")):
        row = spread.loc[term]
        highlight = "**" if row["phosphate"] < row["all"] else ""
        claim(f"estimator spread in {term}",
              f"| {name} | {row['all']:.3f} | "
              f"{highlight}{row['phosphate']:.3f}{highlight} |")
    # The chemistry reason has to stay attached to the number.
    claim("the borate pH window", "maximum at **pH 8.4-9**")
    claim("boric is a scope, not a deletion",
          "0.089 mM against 0.041 and 0.0012")

    print("\nthe burst form as a candidate headline")
    # Every number in the "why not the burst/lag v0" section, re-derived. The
    # case for keeping v0_quad rests on these, so drift here would leave a
    # decision standing on figures that no longer hold.
    frame = scope.frame(scope.FREE_BNOH_ALL)
    difference = frame.v0_burst_resid - frame.v0_quad_resid
    claim("burst clearly better",
          f"| clearly better (> 0.1x noise) | {int((difference > 0.1).sum())} "
          f"curves | {int((difference < -0.1).sum())} curves |")
    claim("indistinguishable",
          f"| indistinguishable | {int((difference.abs() <= 0.1).sum())} curves | |")
    claim("median residuals",
          f"| median residual | {frame.v0_quad_resid.median():.2f}x noise | "
          f"{frame.v0_burst_resid.median():.2f}x noise |")
    claim("bounded burst v0",
          f"| v0 defined on | {len(frame)} of {len(frame)} | "
          f"{int(frame.v0_burst_bounded.sum())} of {len(frame)} bounded |")
    exp65 = frame[frame.experiment == 65].sort_values("sample")
    claim("exp 65 burst residuals",
          ", ".join(f"{v:.1f}" for v in exp65.v0_burst_resid[:3])
          + f" and {exp65.v0_burst_resid.iloc[3]:.1f}x noise")
    claim("exp 65 quadratic residuals",
          "against the quadratic's "
          + ", ".join(f"{v:.1f}" for v in exp65.v0_quad_resid[:3])
          + f" and {exp65.v0_quad_resid.iloc[3]:.1f}x")
    # The degeneracy that makes `bounded` misleading on those four.
    degenerate = exp65[(exp65.v0_burst_kind == "clamped")
                       & (~exp65.tau_resolved) & (exp65.v0_burst_bounded)]
    ok = len(degenerate) == 4
    print(f"  {'pass' if ok else 'FAIL'}  exp 65 burst fits are degenerate yet "
          f"'bounded': {len(degenerate)} of 4")
    if not ok:
        FAILURES.append("exp 65's burst fits are no longer clamped-and-bounded; "
                        "the argument for keeping v0_quad rests on that")
    dropped_boric = scope.buffer_dependence(anchor=(67, 69, 70),
                                            parameter="v0_burst")
    claim("dropping the boric block",
          f"moves it to +{dropped_boric['order_buf']:.2f} ± "
          f"{dropped_boric['stderr_buf']:.2f}")

    print("\nthe outlier rings, and what the first-reading drop left them doing")
    curves = list(scope.curves(scope.FREE_BNOH_ALL))
    rings = zeros = 0
    for curve in curves:
        z = local_outlier_z(curve.times, curve.absorbance, curve.noise)
        isolated, _ = isolated_outliers(curve.times, curve.absorbance,
                                        curve.noise)
        ringed = set(int(i) for i in isolated)
        if len(z) and np.isfinite(z[0]) and abs(z[0]) > OUTLIER_SIGMA:
            ringed.add(0)
        rings += len(ringed)
        zeros += 0 in ringed
    claim("rings over the 27 panels", f"from 16 to {rings}")
    claim("point-0 rings", f"point-0 rings from 14 to {zeros}")

    # The archive-wide split, which moved when the first reading went.
    all_curves, _ = build_curves()
    split = [0, 0]
    for curve in all_curves:
        isolated, in_runs = isolated_outliers(curve.times, curve.absorbance,
                                              curve.noise)
        split[0] += len(isolated)
        split[1] += len(in_runs)
    claim("isolated vs runs", f"**{split[0]}\nisolated against {split[1]} in runs**")

    # The leading reading is no longer the outlier class it was, and the
    # document is required to say so with the current numbers.
    flag = [0, 0]
    for curve in all_curves:
        z = local_outlier_z(curve.times, curve.absorbance, curve.noise)
        if not len(z):
            continue
        flag[0] += bool(np.isfinite(z[0]) and abs(z[0]) > OUTLIER_SIGMA)
        flag[1] += bool(np.isfinite(z[-1]) and abs(z[-1]) > OUTLIER_SIGMA)
    total = len(all_curves)
    claim("leading vs last flag rate",
          f"**{100 * flag[0] / total:.1f}% of leading\nreadings are flagged "
          f"against {100 * flag[1] / total:.1f}% of last ones**")

    print("\nevery ring encloses a drawn reading")
    # Rings are not decimated and the readings are, so a ring could land where
    # no mark was drawn -- it did, on exp 3 samples 1 and 6, and it reads as a
    # flag on nothing. Checked against the rendered SVG rather than the code
    # that emits it, because the two sets are built separately.
    page_text = io.open(CURVES_PAGE, encoding="utf-8").read()
    orphans = []
    for panel in re.split(r"<div class='ph'>", page_text)[1:]:
        name = panel.split("<")[0]
        marks = set(re.findall(
            r"<circle cx='([-\d.]+)' cy='([-\d.]+)' r='1.7' fill='#", panel))
        for spot in re.findall(
                r"<circle cx='([-\d.]+)' cy='([-\d.]+)' r='6.0' fill='none'",
                panel):
            if spot not in marks:
                orphans.append(f"{name} at {spot}")
    ok = not orphans
    print(f"  {'pass' if ok else 'FAIL'}  no ring without a reading inside it"
          + ("" if ok else f": {orphans}"))
    if not ok:
        FAILURES.append(f"rings drawn over nothing: {orphans}")

    print("\nthe rendered page")
    accelerating = sum(
        bool(np.isfinite(z) and z > ACCELERATION_SIGMA)
        for z, _ in (acceleration(c.times, c.absorbance,
                                  floor=source_floor(c.source))
                     for c in curves))
    claim("legend: per-curve gate", f"open on the {accelerating} of 27 that clear it",
          document=CURVES_PAGE)
    claim("legend: the discard is disclosed",
          "first reading of every run is not plotted, fitted or scored",
          document=CURVES_PAGE)
    claim("the constant-buffer accelerating count",
          f"| enzyme-free BnOH, `[buf]` fixed | **{sum(
              bool(np.isfinite(z) and z > ACCELERATION_SIGMA)
              for z, _ in (acceleration(c.times, c.absorbance,
                                        floor=source_floor(c.source))
                           for c in scope.curves(scope.BUFFER_FIXED)))} of 16** |")

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
