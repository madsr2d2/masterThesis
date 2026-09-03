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
sys.path.insert(0, os.path.dirname(HERE))

import scope
from curve_metrics import (ACCELERATION_SIGMA, OUTLIER_SIGMA, acceleration,
                           isolated_outliers, local_outlier_z)
from fit_dataset import build_curves, source_floor
from summary_kinetics import fit_burst, fit_burst_bounded

from doc_check import Checker

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")
# The RENDERED page, checked as well as the prose. The legend claimed the
# burst fit was blanket "B <= 0 ... 0 of 16" for three commits after the gate
# became per-curve, and nothing noticed, because every check here read
# ANALYSIS.md and the stale sentence lived in build_figures.py.
CURVES_PAGE = os.path.join(HERE, "progress_curves.html")
# And the index, for the same reason. Its prose was typed rather than derived
# until 2026-09-02 and four of its numbers had drifted -- see build_figures.py.
INDEX_PAGE = os.path.join(HERE, "index.html")
# MECHANISM.md carries the same boric-probe numbers, in ASCII. A figure quoted
# in two documents is two copies, and this project has had a number go stale in
# one file while the other stayed right.
MECHANISM_DOC = os.path.join(os.path.dirname(HERE), "MECHANISM.md")
ESTIMATORS_ALL = ("v0_quad", "v0_burst", "vmax", "v0", "v0_whole")




def main():
    doc = Checker(DOCUMENT)
    print("checking ANALYSIS.md against the code\n")

    print("the enzyme-free census")
    import pandas as pd
    data = pd.read_csv(os.path.join(os.path.dirname(HERE), "data",
                                    "experiment_data.csv"))
    free = data[data["[enz]"] == 0]
    doc.claim("enzyme-free curves", f"{len(free)} enzyme-free curves over "
                                    f"{free.experiment.nunique()} experiments")
    fixed = scope.frame(scope.BUFFER_FIXED)
    titration = scope.frame(scope.BUFFER_CONFOUNDED)
    doc.claim("BUFFER_FIXED size",
              f"| {', '.join(str(e) for e in scope.BUFFER_FIXED)} | {len(fixed)} |")
    doc.claim("BUFFER_CONFOUNDED size",
              f"| 3, 6 | {len(titration)} ({int(titration.live.sum())} live) |")

    print("\nthe separation")
    quad = scope.buffer_dependence(parameter="v0_quad")
    vmax = scope.buffer_dependence(parameter="vmax")
    doc.claim("a, v0_quad", f"**{quad['order_s0_fixed']:+.3f} ± "
                            f"{quad['stderr_s0_fixed']:.3f}** (n = {quad['n_fixed']})")
    doc.claim("a', v0_quad", f"**{quad['order_s0_titration']:+.3f} ± "
                             f"{quad['stderr_s0_titration']:.3f}** "
                             f"(n = {quad['n_titration']})")
    doc.claim("a, vmax", f"{vmax['order_s0_fixed']:+.3f} ± "
                         f"{vmax['stderr_s0_fixed']:.3f} (n = {vmax['n_fixed']})")
    doc.claim("coupling g", f"g = {quad['coupling']:.3f}")
    doc.claim("d, v0_quad", f"**d = {quad['order_buf']:+.2f} ± "
                            f"{quad['stderr_buf']:.2f}** on `v0_quad`")

    print("\nthe result table")
    for estimator, marker in (("v0_quad", "**"), ("vmax", ""),
                              ("v0", ""), ("v0_whole", "")):
        result = scope.buffer_dependence(parameter=estimator)
        cross = scope.background_orders(scope.FREE_4OME_40C,
                                        terms=("s0", "buf"), parameter=estimator)
        doc.claim(f"{estimator} row",
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
    doc.claim("joint VIF", f"VIF {joint['vif_s0']:.1f} and {joint['vif_buf']:.1f} "
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
    doc.claim("negative v0 count",
              f"**{negative} of {len(curves)} curves return a negative v₀**")
    always = [fit_burst_bounded(c.times, c.absorbance, constrain=True,
                                noise_floor=source_floor(c.source))
              for c in curves]
    doc.claim("bounded count, lag branch shut everywhere",
              f"from **{sum(f.bounded for f in bounded_free)} to "
              f"{sum(f.bounded for f in always)} of {len(curves)}**")
    # `bounded` is the DEFAULT, constrain="auto": the branch is gated on each
    # curve's own acceleration rather than shut for the whole block.
    doc.claim("bounded count, lag branch gated per curve",
              f"That bounds **{sum(f.bounded for f in bounded)} of {len(curves)}**")

    print("\nshape and scope")
    live_fixed = fixed[fixed.live]
    inscope = scope.frame()
    live_inscope = inscope[inscope.live]
    doc.claim("background acceleration",
              f"| **{int((live_fixed.accel_z > ACCELERATION_SIGMA).sum())} of "
              f"{len(live_fixed)}** |")
    doc.claim("two-axis acceleration",
          f"| **{int((live_inscope.accel_z > ACCELERATION_SIGMA).sum())} of "
              f"{len(live_inscope)} live** |")
    buffers = sorted(inscope.buf.round(3).unique())
    doc.claim("two-axis buffer constant",
              f"`[buf]` = {buffers[0]:.3f} mM in all {len(inscope)} two-axis curves")
    everything = scope.frame(scope.FREE_BNOH_ALL)
    curved = int((everything.curvature_t.abs() > 3).sum())
    doc.claim("curvature count", f"{curved}\nof {len(everything)} curves show "
                                 f"curvature at |t| > 3")

    print("\nthe estimator description table")
    # This table carries no fitted numbers, so nothing else here would notice
    # if it were clobbered -- and a bulk regex did clobber it in c41f459..
    # 4b344ac, replacing the descriptions with buffer orders for three
    # commits. Pinned by content now.
    for phrase in ("| `vmax` | steepest of four 20% blocks | no | straight line |",
                   "| `v0` | first 20% | no | straight line |",
                   "| `v0_whole` | none | yes | straight line |"):
        doc.claim("description row", phrase)

    print("\nthe accelerating-curve sensitivity")
    for estimator in ESTIMATORS_ALL:
        dropped = scope.buffer_dependence(parameter=estimator,
                                          drop_accelerating=True)
        doc.claim(f"{estimator} dropped-accelerating",
                  f"| {dropped['order_buf']:+.2f} ± {dropped['stderr_buf']:.2f} |")

    print("\nthe two-axis substrate order")
    orders = scope.orders("vmax", scope.TWO_AXIS_BLOCK, within=True)
    strong = scope.orders("vmax", scope.strong_runs(), within=True)
    doc.claim("two-axis order", f"{orders['order_s0']:+.3f} ± "
                                f"{orders['stderr_s0']:.3f} over all "
                                f"{orders['n']} live two-axis curves")
    doc.claim("strong-run order", f"({strong['order_s0']:+.3f} ± "
                                  f"{strong['stderr_s0']:.3f} on the "
                                  f"{len(scope.strong_runs())} strong runs)")

    print("\nthe rate law")
    # Pinned at last. The line in section 6 was the PRE-first-reading-drop
    # value from 2026-09-01 to b144b9e -- [H2O2]^+1.78 where the code said
    # +1.87, [HOO-]^+0.68 where it said +0.63 -- because no check re-derived
    # it. Same failure as the legend and the estimator table before it.
    #
    # ONE SCOPE, not two. The "all six runs" law was withdrawn on 2026-09-01
    # when exp 65's rates were ruled unusable, and this loop checked it until
    # then. A withdrawn number must stop being checked, or the check keeps it
    # alive in the document.
    for estimator in ("v0_quad", "vmax"):
        law = scope.background_orders(scope.FREE_BNOH_PHOSPHATE,
                                      terms=("s0", "h2o2", "hoo"),
                                      parameter=estimator)
        buf = scope.buffer_dependence(parameter=estimator)
        doc.claim(f"rate law, {estimator}",
                  f"[S]^{buf['order_s0_fixed']:+.2f}  "
                  f"[H2O2]^{law['order_h2o2']:+.2f}  "
                  f"[HOO-]^{law['order_hoo']:+.2f}  "
                  f"[buf]^{buf['order_buf']:+.2f}     ({estimator})")
        doc.claim(f"pooled R2, {estimator}", f"{law['r2']:.3f}")
    # And that the withdrawn version is really gone from the document.
    doc.claim("no all-six rate law survives", "All six runs, boric included",
              present=False)

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
    doc.claim("species and total are one variable in the titration",
              f"correlation **{worst:.6f}**")
    inflation_total = scope.variance_inflation(
        live, ("s0", "h2o2", "hoo", "buf"))["buf"]
    inflation_base = scope.variance_inflation(
        live, ("s0", "h2o2", "hoo", "buf_base"))["buf_base"]
    doc.claim("substituting the basic form inflates the variance",
              f"factor to **{inflation_base:.1f}**, against **{inflation_total:.1f}** for the")
    # The two pH levels, and that there are only two of them.
    levels = sorted(live.pH.unique())
    ok = len(levels) == 2
    print(f"  {'pass' if ok else 'FAIL'}  the phosphate set has two pH levels: "
          f"{levels}")
    if not ok:
        doc.fail("the phosphate set no longer has exactly two pH levels; "
                            "section 6b's degeneracy argument needs rewriting")
    span = live.groupby("pH").agg(buf=("buf", "median"),
                                  acid=("buf_acid", "median"),
                                  base=("buf_base", "median"),
                                  hoo=("hoo", "median"))
    doc.claim("what moves between the two pH levels",
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
        doc.claim(f"boric probe, {estimator}",
                  f"| `{estimator}` | {row['predicted']:.2f}x | "
                  f"{row['observed']:.2f}x | **{row['excess']:.2f}x**{note} |")
        doc.claim(f"boric probe in MECHANISM.md, {estimator}",
                  f"{row['excess']:.2f}x ({estimator})", document=MECHANISM_DOC)
    # The claim that three of four fall BELOW the law, which is the sentence
    # the conclusion is written on.
    below = int((peroxo.excess < 1).sum())
    ok = below == 3
    print(f"  {'pass' if ok else 'FAIL'}  three of four estimators put boric "
          f"below the phosphate law: {below} of {len(peroxo)}")
    if not ok:
        doc.fail("the boric probe no longer puts 3 of 4 estimators below "
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
        doc.fail("exps 65 and 67 no longer match in [S] and [H2O2]; the "
                            "boric probe's prediction needs the other two orders")

    print("\nthe mid-run break that weakens the boric probe")
    # Found from the plots, not by any statistic in the package -- every other
    # shape column compares a curve's start to its end and steps over a break
    # in the middle. Section 6b's downgrade of the boric probe rests on this,
    # so it is re-derived rather than typed.
    #
    # THE POPULATION IS THE ENZYME-FREE RUNS AND ONLY THOSE. The first version
    # of this check passed exps 59-62 as controls; they are catalysed, their
    # reference channel already subtracted the background, and they cannot
    # show a background shape. `synchronised_break` now raises on the mix.
    census = pd.read_csv(os.path.join(os.path.dirname(HERE), "data",
                                      "experiment_data.csv"))
    # EVERY row must be enzyme-free, not merely some: exp 128 is a catalysed
    # run whose reference row 5 carries [enz] = 0, and filtering on rows put it
    # in the background population until `synchronised_break` refused it.
    highest = census.groupby("experiment")["[enz]"].max()
    background = tuple(sorted(highest[highest == 0].index))
    breaks = scope.synchronised_break(background)
    exp65 = breaks.loc[scope.PEROXO_PAIR[0]]
    times = [str(v) for v in sorted(exp65["breaks"])]
    doc.claim("exp 65's break times",
              f"{', '.join(times[:-1])} and {times[-1]} s")
    ratios = [f"{v:.2f}" for v in sorted(exp65["ratios"])]
    doc.claim("exp 65's steepening",
              f"**{', '.join(ratios[:-1])} and {ratios[-1]}x**")
    doc.claim("the break span", f"a span of **{exp65['span']:.0f} s**")

    # The claim the whole section now turns on: among the runs where a
    # background shape is VISIBLE, exp 65 is the only boric one and the only
    # one that breaks. Both halves are checked, because either failing would
    # change the conclusion in a different direction.
    salt = census.groupby("experiment").buffer.first()
    boric = [e for e in breaks.index if salt.get(e) == "Boric"]
    ok = boric == [scope.PEROXO_PAIR[0]]
    print(f"  {'pass' if ok else 'FAIL'}  exactly one boric background run has "
          f"live curves: {boric}")
    if not ok:
        doc.fail("the boric background population is no longer exp 65 "
                            "alone; section 6b's 1-of-1 count needs rewriting")
    whole = breaks[breaks.steep == breaks.n]
    ok = list(whole.index) == [scope.PEROXO_PAIR[0]]
    print(f"  {'pass' if ok else 'FAIL'}  and it is the only background run "
          f"whose cuvettes all steepen: {list(whole.index)} of "
          f"{len(breaks)} runs")
    if not ok:
        doc.fail("exp 65 is no longer the only background run that "
                            "steepens throughout; section 6b needs rewriting")
    doc.claim("the count both ways",
              f"**1 boric background run of 1 showing the break, against\n"
              f"{len(breaks) - 1} of {len(breaks) - 1} not**")
    # Exp 64, the only other boric background run, stops BEFORE exp 65 breaks,
    # which is why it cannot serve as the control it was first used as.
    from read_rre import read_rre
    last = max(float(np.asarray(t)[-1]) for _, t, _ in
               read_rre(os.path.join(os.path.dirname(HERE), "data", "Mads",
                                     "rate064.rre")))
    ok = last < float(exp65["breaks"][0])
    print(f"  {'pass' if ok else 'FAIL'}  exp 64 ends at {last:.0f} s, before "
          f"exp 65's first break at {exp65['breaks'][0]} s")
    if not ok:
        doc.fail("exp 64 now runs past exp 65's break, so it could test "
                            "it after all -- section 6b says it cannot")

    print("\nthe withdrawn second probe, kept as a catalysed comparison")
    catalysed = scope.peroxo_buffer_test(pair=scope.CATALYSED_PEROXO_PAIR,
                                         orders_scope=None)
    for estimator, row in catalysed.iterrows():
        doc.claim(f"catalysed probe, {estimator}",
                  f"| `{estimator}` | {row['observed']:.2f}x |")
    # Quoted in MECHANISM.md as one line; checked as one string so it cannot
    # collide with the enzyme-free probe's numbers, which share a value.
    doc.claim("catalysed probe in MECHANISM.md",
              "boric/phosphate = " + ", ".join(
                  f"{row['observed']:.2f}x ({estimator})"
                  for estimator, row in catalysed.iterrows()),
              document=MECHANISM_DOC)
    pair = scope.frame(scope.CATALYSED_PEROXO_PAIR)
    pair = pair[pair.live]
    hoo = float(pair[pair.experiment == scope.CATALYSED_PEROXO_PAIR[0]].hoo.median()
                / pair[pair.experiment == scope.CATALYSED_PEROXO_PAIR[1]].hoo.median())
    doc.claim("the hydroperoxide it carries", f"**{hoo:.2f}x the [HOO-]**")

    print("\nthe buffer order is a phosphate number")
    # buffer_dependence's default anchor includes the boric run. Section 6b
    # says the order survives dropping it; that is checked, not asserted.
    for estimator, key in (("vmax", "vmax"), ("v0_quad", "v0_quad")):
        without = scope.buffer_dependence(
            anchor=tuple(e for e in scope.BUFFER_FIXED
                         if e not in scope.BORIC_BUFFER),
            parameter=estimator)
        doc.claim(f"buffer order without boric, {estimator}",
                  f"**+{without['order_buf']:.2f} +/- {without['stderr_buf']:.2f}**")

    print("\nwhat the boric run carries")
    spread = scope.boric_spread()
    for term, name in (("s0", "`[S]`"), ("h2o2", "`[H2O2]`"),
                       ("hoo", "`[HOO-]`"), ("buf", "`[buf]`")):
        row = spread.loc[term]
        highlight = "**" if row["phosphate"] < row["all"] else ""
        doc.claim(f"estimator spread in {term}",
                  f"| {name} | {row['all']:.3f} | "
                  f"{highlight}{row['phosphate']:.3f}{highlight} |")
    # The chemistry reason has to stay attached to the number.
    doc.claim("the borate pH window", "maximum at **pH 8.4-9**")
    doc.claim("boric is a scope, not a deletion",
              "0.089 mM against 0.041 and 0.0012")

    print("\nthe burst form as a candidate headline")
    # Every number in the "why not the burst/lag v0" section, re-derived. The
    # case for keeping v0_quad rests on these, so drift here would leave a
    # decision standing on figures that no longer hold.
    frame = scope.frame(scope.FREE_BNOH_ALL)
    difference = frame.v0_burst_resid - frame.v0_quad_resid
    doc.claim("burst clearly better",
          f"| clearly better (> 0.1x noise) | {int((difference > 0.1).sum())} "
              f"curves | {int((difference < -0.1).sum())} curves |")
    doc.claim("indistinguishable",
              f"| indistinguishable | {int((difference.abs() <= 0.1).sum())} curves | |")
    doc.claim("median residuals",
          f"| median residual | {frame.v0_quad_resid.median():.2f}x noise | "
              f"{frame.v0_burst_resid.median():.2f}x noise |")
    doc.claim("bounded burst v0",
              f"| v0 defined on | {len(frame)} of {len(frame)} | "
              f"{int(frame.v0_burst_bounded.sum())} of {len(frame)} bounded |")
    exp65 = frame[frame.experiment == 65].sort_values("sample")
    doc.claim("exp 65 burst residuals",
              ", ".join(f"{v:.1f}" for v in exp65.v0_burst_resid[:3])
              + f" and {exp65.v0_burst_resid.iloc[3]:.1f}x noise")
    doc.claim("exp 65 quadratic residuals",
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
        doc.fail("exp 65's burst fits are no longer clamped-and-bounded; "
                            "the argument for keeping v0_quad rests on that")
    dropped_boric = scope.buffer_dependence(anchor=(67, 69, 70),
                                            parameter="v0_burst")
    doc.claim("dropping the boric block",
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
    doc.claim("rings over the 27 panels", f"from 16 to {rings}")
    doc.claim("point-0 rings", f"point-0 rings from 14 to {zeros}")

    # The archive-wide split, which moved when the first reading went.
    all_curves, _ = build_curves()
    split = [0, 0]
    for curve in all_curves:
        isolated, in_runs = isolated_outliers(curve.times, curve.absorbance,
                                              curve.noise)
        split[0] += len(isolated)
        split[1] += len(in_runs)
    doc.claim("isolated vs runs",
              f"**{split[0]}\nisolated against {split[1]} in runs**")

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
    doc.claim("leading vs last flag rate",
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
        doc.fail(f"rings drawn over nothing: {orphans}")

    print("\nthe rendered page")
    accelerating = sum(
        bool(np.isfinite(z) and z > ACCELERATION_SIGMA)
        for z, _ in (acceleration(c.times, c.absorbance,
                                  floor=source_floor(c.source))
                     for c in curves))
    doc.claim("legend: per-curve gate", f"open on the {accelerating} of 27 that clear it",
              document=CURVES_PAGE)
    doc.claim("legend: the discard is disclosed",
              "first reading of every run is not plotted, fitted or scored",
              document=CURVES_PAGE)
    # The DENOMINATOR is derived too. It was hardcoded at 16 and survived
    # BUFFER_FIXED dropping exp 65 on 2026-09-01, so the check went on
    # comparing a live numerator against a dead total.
    anchor = list(scope.curves(scope.BUFFER_FIXED))
    rising = sum(bool(np.isfinite(z) and z > ACCELERATION_SIGMA)
                 for z, _ in (acceleration(c.times, c.absorbance,
                                           floor=source_floor(c.source))
                              for c in anchor))
    doc.claim("the constant-buffer accelerating count",
              f"| enzyme-free BnOH, `[buf]` fixed | **{rising} of {len(anchor)}** |")


    print("\nsection 7: what the background is not")
    # Audited 2026-09-02. The conversion range read 0.04-7.9% and the exp 69
    # claim "roughly halves"; the range was stale and the verb understated a
    # curve that keeps an eighth of its rate. Neither had ever been checked.
    free = scope.frame(scope.FREE_BNOH_ALL)
    doc.claim("the conversion range",
              f"{free.conversion.min() * 100:.2f}-{free.conversion.max() * 100:.2f}%")
    doc.claim("how many curves bend",
          f"{int((free.curvature_t.abs() > 3).sum())} of {len(free)} curves "
              "show curvature at |t| > 3")
    worst = free[(free.experiment == 69) & (free["sample"] == 3)].iloc[0]
    doc.claim("what exp 69 sample 3 does",
              f"ends at **{worst.late_over_early:.0%} of its early rate** at "
              f"{worst.conversion * 100:.2f}% conversion")

    print("\nwhat product_fate narrowed the open question to")
    import slowdown
    whole_frame = scope.frame(tuple(range(1, 152)))
    free4 = slowdown.substrate_blocks(whole_frame)["4OMe enzyme-free"]
    row = slowdown.deceleration_drivers(free4)
    doc.claim("the clock loads on run length",
              f"{row['span']:+.3f} +/- {row['span_stderr']:.3f} on run length")
    doc.claim("and not on product",
              f"{row['product']:+.3f} +/- {row['product_stderr']:.3f} on\nproduct")
    live4 = free4[free4.live]
    doc.claim("over what absorbances",
              f"over {live4.net.min():.3f}-{live4.net.max():.3f} AU")

    print("\nthe pyrophosphate cell the amplitude is missing from")
    whole = scope.frame(tuple(range(1, 152)))
    cell = whole[(whole.buffer == "Pyrophosphate") & (whole.substrate == "BnOH")]
    doc.claim("its size, and that it holds no enzyme-free curve",
              f"{int((~cell.differential).sum())} enzyme-free curves in the "
              f"{len(cell)}-curve BnOH pyrophosphate cell")

    print("\nthe page: every number in its prose is derived, not typed")
    # The index and the curves page quote figures the document also quotes.
    # They used to be typed and four of them had drifted; they are f-strings
    # now, and these checks read them back OUT of the rendered HTML so that a
    # future edit cannot quietly retype one.
    page = io.open(INDEX_PAGE, encoding="utf-8").read()
    pooled = {name: scope.background_orders(terms=("s0", "h2o2", "hoo"),
                                            parameter=name)
              for name in ("v0_quad", "vmax")}
    held = sorted(scope.frame(scope.BUFFER_FIXED).buf.unique())
    in_scope = scope.frame(scope.TWO_AXIS_BLOCK)
    doc.claim("the buffer the fixed design holds", f"[buf] at {held[0]:g} mM",
              document=INDEX_PAGE)
    doc.claim("the pooled R2 on v0_quad",
              f"R\u00b2 {pooled['v0_quad']['r2']:.3f} on <code>v0_quad</code>",
              document=INDEX_PAGE)
    doc.claim("the pooled R2 on vmax",
              f"{pooled['vmax']['r2']:.3f} on <code>vmax</code>",
              document=INDEX_PAGE)
    doc.claim("the two-axis buffer",
              f"{in_scope.buf.iloc[0]:.3f} mM in all {len(in_scope)}",
              document=INDEX_PAGE)
    doc.claim("the pyrophosphate cell on the page", f"{len(cell)}-curve BnOH",
              document=INDEX_PAGE)
    doc.check("the withdrawn 0.961/0.909 comparison is named as withdrawn",
              "argument is withdrawn" in page or "0.961" not in page)

    print("\nthe figures: lettered once each, in order")
    doc.figures(os.path.join(HERE, "index.html"))

    print("\nthe figures: no data point drawn outside its own frame")
    # Marks are clipped to the plot area deliberately, so a data point outside
    # the axis limits vanishes silently and the figure still looks complete.
    # That is how exp 6 sample 4 went missing from the curvature figure for as
    # long as it existed: no number changes, so no prose check could see it.
    doc.unclipped(os.path.join(HERE, "index.html"),
                  os.path.join(HERE, "progress_curves.html"))

    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
