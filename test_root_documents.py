"""
The document contract, applied to the documents that had no gate at all.

    python test_root_documents.py

Every folder's `check_numbers.py` re-derives its own `ANALYSIS.md` from the
modules -- 883 claims across eight documents, and that machinery has worked
well. It stopped at the folder boundary. Until 2026-09-08 the ROOT documents
were checked barely at all: `background_reaction` and `two_axis` each verify a
handful of strings in `MECHANISM.md` and `FITTING.md`, and CLAUDE.md,
BUBBLES.md and the `analyse-kinetics` skill had NOTHING pointed at them --
between them about 480 numeric tokens, in the three documents that are read
FIRST and that steer every future analysis.

That gap is where the drift went, exactly as the duplicate guard's own history
predicted. An end-to-end review on 2026-09-08 found, in ungated prose only:

  * the skill file quoting a lag fraction of 37.6% (151/402) when the code
    says 160/402 and `test_lag_statistic`'s docstring says never to quote it
    to three digits at all -- two revisions stale;
  * MECHANISM.md quoting the same superseded 151/402, in bold;
  * the skill's `bubble_synchrony` at 17-against-16.0 where it is 23 over 357
    cuvette pairs against 21.3, its gas-rate order at +1.417 where it is
    +1.477, its `bubble_load` count at thirteen where it is fourteen, its
    `tau_slow` resolution at 33 where it is 34, and its joint-clock sigma at
    1.4 where the gains-corrected value is 0.3;
  * BUBBLES.md calling 27-of-243 an ARCHIVE-WIDE excursion count when it is
    the two-axis block's -- archive-wide the same pipeline sees 40 of 414.

None of those was wrong when written. Each was correct, was superseded, and
had no check to notice. This module is that check. It is deliberately NOT a
full re-derivation of every number in the root documents -- it covers the
load-bearing ones and the ones that have actually drifted, which is what a
gate is for.

`doc_check.Checker` is the same comparison every folder runs on, so the
typography rules are identical: hyphens fold against U+2212, backticks and
rewrapping are noise, and `**` survives because bold marks the number a
passage argues for.
"""
import os
import sys

import numpy as np

REPOSITORY = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(REPOSITORY, "data"))
sys.path.insert(0, REPOSITORY)

import curve_metrics
import fit_dataset
import induction
import ph_role
import scope
from doc_check import Checker

CLAUDE = os.path.join(REPOSITORY, "CLAUDE.md")
BUBBLES = os.path.join(REPOSITORY, "BUBBLES.md")
MECHANISM = os.path.join(REPOSITORY, "MECHANISM.md")
FITTING = os.path.join(REPOSITORY, "FITTING.md")
SKILL = os.path.join(REPOSITORY, ".claude", "skills", "analyse-kinetics",
                     "SKILL.md")

# The documents that carry numbers and now have a gate. `README.md` and
# `COMPUTATIONAL.md` are deliberately absent: the first quotes almost nothing,
# and the second is a task register whose numbers are TARGETS for calculations
# not yet run, so there is nothing in the modules to check them against.
GUARDED_DOCUMENTS = (CLAUDE, BUBBLES, MECHANISM, FITTING, SKILL)


def _excursion_tally(curves):
    """
    Candidate falls, rejections and emptied curves, over a set of curves.

    `detachments` returns only the events it KEEPS, so the rejection count --
    which BUBBLES.md and CLAUDE.md both quote -- cannot be read back off it.
    This runs the same three stages the function does (SNR gate, `bubble_drops`,
    group consecutive, drop excursions) and counts what falls out between them.
    """
    candidates = rejected = 0
    emptied = []
    for curve in sorted(curves, key=lambda c: (c.experiment, c.sample)):
        values = np.asarray(curve.absorbance, dtype=float)
        noise = curve.noise
        if len(values) < 2 or not np.isfinite(noise) or noise <= 0:
            continue
        if float(values[-1] - values[0]) / noise < curve_metrics.DETACHMENT_SNR_FLOOR:
            continue
        drops = curve_metrics.bubble_drops(values, noise)
        if not len(drops):
            continue
        events, start = [], int(drops[0])
        previous = start
        for index in drops[1:]:
            index = int(index)
            if index == previous + 1:
                previous = index
            else:
                events.append((start, previous + 1))
                start = previous = index
        events.append((start, previous + 1))
        kept = [e for e in events
                if not curve_metrics._is_excursion(values, e)]
        candidates += len(events)
        rejected += len(events) - len(kept)
        if events and not kept:
            emptied.append((int(curve.experiment), int(curve.sample)))
    return {"candidates": candidates, "rejected": rejected,
            "emptied": emptied}


def main():
    doc = Checker(CLAUDE, label="the repository root",
                  document_label="CLAUDE.md + BUBBLES.md + MECHANISM.md + "
                                 "FITTING.md + SKILL.md")

    doc.section("the lag statistic, in every document that quotes it")
    # The number that drifted furthest, and the one the code explicitly refuses
    # to have quoted precisely. It has moved four times without the statistic
    # changing once, so every document has to say "about 40%".
    curves, _ = fit_dataset.build_curves()
    positions = np.array([curve_metrics.peak_position(c.absorbance, c.times)
                          for c in curves], dtype=float)
    lagging = int(np.nansum(positions > curve_metrics.LAG_THRESHOLD))
    doc.check("402 fittable curves", len(curves) == 402, f"{len(curves)}")
    doc.check("160 of them lag", lagging == 160, f"{lagging}")
    for document in (FITTING, MECHANISM, SKILL):
        doc.claim(f"{os.path.basename(document)}: the count",
                  f"{lagging}/402", document=document)
        doc.claim(f"{os.path.basename(document)}: quoted as about 40%",
                  "about 40%", document=document)
    # And that the two superseded values are no longer asserted as current.
    # 151/402 survives in FITTING.md's history table and in MECHANISM.md's
    # and the skill's account of WHY it moved, so this checks the phrasing
    # that made it a live claim rather than the digits.
    doc.claim("MECHANISM.md no longer argues the superseded 37.6%",
              "**37.6%** — 151/402", present=False, document=MECHANISM)
    doc.claim("...and neither does the skill",
              "now **37.6%** (151/402)", present=False, document=SKILL)

    doc.section("the O2 artefact: the numbers BUBBLES.md and CLAUDE.md share")
    block_curves = [c for c in curves
                    if int(c.experiment) in scope.TWO_AXIS_BLOCK]
    block = _excursion_tally(block_curves)
    archive = _excursion_tally(curves)
    doc.check("the block's candidate falls", block["candidates"] == 243,
              f"{block['candidates']}")
    doc.check("the block's rejections", block["rejected"] == 27,
              f"{block['rejected']}")
    doc.check("two block curves lose all of theirs",
              len(block["emptied"]) == 2, f"{block['emptied']}")
    # The distinction BUBBLES.md got wrong: this is a BLOCK count, and the
    # archive-wide one is half again as large.
    doc.check("archive-wide it is a different number",
              archive["candidates"] == 414 and archive["rejected"] == 40
              and len(archive["emptied"]) == 5,
              f"{archive['candidates']} candidates, {archive['rejected']} "
              f"rejected, {len(archive['emptied'])} emptied")
    for document in (CLAUDE, BUBBLES, SKILL):
        doc.claim(f"{os.path.basename(document)}: the block's tally",
                  f"{block['rejected']} of {block['candidates']}",
                  document=document)
    doc.claim("BUBBLES.md no longer calls the block count archive-wide",
              "of 243 archive-wide", present=False, document=BUBBLES)
    doc.claim("BUBBLES.md gives the archive-wide pair instead",
              f"{archive['rejected']} of {archive['candidates']}",
              document=BUBBLES)
    doc.claim("and the skill does too", f"{archive['rejected']} of "
              f"{archive['candidates']}", document=SKILL)

    synchrony = scope.bubble_synchrony()
    doc.claim("CLAUDE.md: the synchrony control",
              f"{synchrony['observed']} coincidences over "
              f"{synchrony['pairs']} cuvette pairs against "
              f"{synchrony['expected']:.1f} expected")
    doc.claim("the skill: the same control, same numbers",
              f"{synchrony['observed']} coincidences over "
              f"{synchrony['pairs']} cuvette", document=SKILL)

    gas = scope.gas_rate_drivers()
    for document in (CLAUDE, BUBBLES, SKILL):
        doc.claim(f"{os.path.basename(document)}: the fitted gas rate's orders",
                  f"+{gas['pooled_h2o2']:.3f} ± "
                  f"{gas['pooled_stderr_h2o2']:.3f} in peroxide",
                  document=document)
    doc.claim("CLAUDE.md: and its substrate order",
              f"{gas['order_s0']:+.3f} ± {gas['stderr_s0']:.3f} in substrate")

    live = scope.frame()
    live = live[live.live]
    loaded = int((scope.bubble_table().bubble_load > 1).sum())
    doc.check("14 live curves carry bubble_load > 1", loaded == 14,
              f"{loaded}")
    doc.claim("CLAUDE.md: the flagged count", f"{loaded} of {len(live)} live")
    doc.claim("the skill: the same count, spelled out",
              "Fourteen of 110", document=SKILL)
    # The list of runs they sit in -- the skill omitted exp 149 for a while.
    flagged = sorted({int(e) for e in
                      scope.bubble_table().query("bubble_load > 1").experiment})
    doc.check("they sit in exps 135, 138, 140, 141, 142, 149, 150",
              flagged == [135, 138, 140, 141, 142, 149, 150], f"{flagged}")
    for document in (CLAUDE, SKILL):
        doc.claim(f"{os.path.basename(document)}: and it names all of them",
                  "138, 140, 141, 142, 149 and 150", document=document)

    doc.section("the clocks, corrected and not")
    resolved = {name: int(live[name].sum()) for name in
                ("tau_resolved", "tau_resolved_corrected",
                 "tau_slow_resolved", "tau_slow_resolved_corrected")}
    doc.check("the correction buys resolution: 62 to 69, 25 to 34",
              resolved == {"tau_resolved": 62, "tau_resolved_corrected": 69,
                           "tau_slow_resolved": 25,
                           "tau_slow_resolved_corrected": 34}, f"{resolved}")
    doc.claim("CLAUDE.md: the resolution it buys",
              f"({resolved['tau_resolved']} to "
              f"{resolved['tau_resolved_corrected']} and "
              f"{resolved['tau_slow_resolved']} to "
              f"{resolved['tau_slow_resolved_corrected']} curves)")
    doc.claim("the skill: the resolved count it must quote",
              f"resolved on {resolved['tau_slow_resolved_corrected']} of "
              f"{len(live)} curves", document=SKILL)

    clocks = induction.joint_clocks(scope.frame())
    slow = clocks.loc[("tau_slow_corrected", "axis")]
    fast = clocks.loc[("tau_corrected", "axis")]
    doc.check("the fully-corrected tau_slow row sits 0.8 sigma from +1",
              abs(slow.sigma - 0.84) < 0.05, f"{slow.sigma:.3f}")
    doc.check("and the tau row 2.3", abs(fast.sigma - 2.30) < 0.05,
              f"{fast.sigma:.3f}")
    for document in (CLAUDE, SKILL):
        doc.claim(f"{os.path.basename(document)}: both clock sigmas",
                  f"{fast.sigma:.1f}σ and {slow.sigma:.1f}σ"
                  if document is CLAUDE else
                  f"{fast.sigma:.1f} and {slow.sigma:.1f}", document=document)
    # How many curves the correction actually moves -- the count that says the
    # repair is doing something even where it does not change a conclusion.
    moved = {
        "tau": int((~np.isclose(live.tau.fillna(-1),
                                live.tau_corrected.fillna(-1))).sum()),
        "tau_slow": int((~np.isclose(live.tau_slow.fillna(-1),
                                     live.tau_slow_corrected.fillna(-1))).sum()),
    }
    doc.check("it moves 38 curves' tau and 34 curves' tau_slow",
              moved == {"tau": 38, "tau_slow": 34}, f"{moved}")
    for document in (CLAUDE, SKILL):
        doc.claim(f"{os.path.basename(document)}: the curves it moves",
                  f"{moved['tau_slow']} of {len(live)} live curves",
                  document=document)

    doc.section("the two-axis block's own headline orders")
    # CLAUDE.md is written in ASCII throughout -- "+/-" and "sigma", where the
    # folder documents use "±" and "σ". `normalise` folds chi2 and superscripts
    # but not those, so the claims are rendered the way each document writes.
    published = scope.ph_order(parameter="vmax", scope=scope.strong_runs())
    pooled = published.loc["pooled"]
    doc.claim("CLAUDE.md: the pH order",
              f"+{pooled.order:.3f} +/- {pooled.stderr:.3f}")
    doc.claim("...and that it is called a half order",
              "a HALF order in [HOO-]")
    consistency = scope.hoo_consistency()
    doc.claim("CLAUDE.md: the [HOO-] within/across gap",
              f"{consistency['sigma']:.1f} sigma")
    doc.check("that gap really is 3.4 sigma",
              abs(consistency["sigma"] - 3.4) < 0.05,
              f"{consistency['sigma']:.3f}")

    doc.section("the pH ladders, and what may be quoted off them")
    three = ph_role.pooled_rate_order(drop=("boric 4OMe",))
    doc.claim("CLAUDE.md: three ladders agree",
              f"+{three['pooled']:.3f} +/- {three['stderr']:.3f}")
    doc.claim("CLAUDE.md: and the chi2 that says so",
              f"chi2 = {three['chi2']:.2f} on {three['dof']}")
    doc.claim("the skill: the same pooled order",
              f"+{three['pooled']:.3f} +/- {three['stderr']:.3f}",
              document=SKILL)
    independent = ph_role.independent_check()
    # The overlap that ph/ANALYSIS.md overstated. CLAUDE.md and the skill both
    # have to say which ladder actually corroborates the block.
    doc.claim("the skill: phosphate is the independent one",
              f"+{independent['phosphate_order']:.3f} +/- "
              f"{independent['phosphate_stderr']:.3f}", document=SKILL)
    doc.claim("the skill: and how far it lands from the block",
              f"{independent['sigma']:.2f} sigma", document=SKILL)
    doc.check("the two pyrophosphate ladders ARE the block",
              set(scope.PH_LADDER_TWO_AXIS_LOW) <= set(scope.TWO_AXIS_BLOCK)
              and set(scope.PH_LADDER_TWO_AXIS_HIGH)
              <= set(scope.TWO_AXIS_BLOCK))
    doc.claim("the skill says so rather than calling them independent",
              "ARE the two-axis block", document=SKILL)
    sweep = ph_role.boric_split_sensitivity()
    doc.check("the boric split systematic is larger than its error",
              sweep["span"] > 5 * sweep["published_stderr"],
              f"span {sweep['span']:.3f} against "
              f"{sweep['published_stderr']:.3f}")

    doc.section("the concentration-agreement screen, as the skill states it")
    agreement = scope.concentration_agreement()
    strong = scope.strong_runs()
    doc.check("exp 151 gets no row -- fewer than four live cuvettes",
              151 not in agreement.index)
    doc.check("and is therefore not a strong run", 151 not in strong)
    doc.check("the weak runs go negative, to -0.612 on exp 150",
              abs(agreement.agreement.min() + 0.612) < 0.005,
              f"{agreement.agreement.min():.3f}")
    doc.claim("the skill quotes the negative floor, not 0.005",
              "-0.612 on exp 150", document=SKILL)
    doc.claim("...and no longer the superseded 0.005 floor",
              "as low as 0.005", present=False, document=SKILL)
    doc.check("exps 135/138/139/140/142 run 0.92 to 0.97",
              abs(agreement.loc[[135, 138, 139, 140, 142]].agreement.min()
                  - 0.922) < 0.005
              and abs(agreement.loc[[135, 138, 139, 140, 142]].agreement.max()
                      - 0.974) < 0.005,
              f"{agreement.loc[[135, 138, 139, 140, 142]].agreement.min():.3f}"
              f" to "
              f"{agreement.loc[[135, 138, 139, 140, 142]].agreement.max():.3f}")
    doc.claim("the skill: that range", "0.92 to\n  0.97", document=SKILL)

    doc.section("every guarded document is real and was actually read")
    # The guard is worth nothing if a path stopped resolving: a missing file
    # would make every `present=False` claim pass by finding nothing.
    for document in GUARDED_DOCUMENTS:
        doc.check(f"{os.path.relpath(document, REPOSITORY)} exists and is read",
                  os.path.exists(document) and len(doc.text(document)) > 2000,
                  f"{len(doc.text(document)) if os.path.exists(document) else 0}"
                  " chars")
    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
