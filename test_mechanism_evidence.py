"""
The gate for `MECHANISM_EVIDENCE.md`: every number re-derived from the modules.

    python test_mechanism_evidence.py

The document is the evidence register for the mechanism question -- what the
kinetics archive decides, what it cannot, and what would decide the rest. It
quotes results from three places (`rate_laws`, `mechanism_discrimination` and
`buffer_role`), and a register whose numbers drift is worse than no register,
so every one of them is re-derived here through the function that owns it
rather than restated as a literal.

`doc_check.Checker` is the same comparison every folder's `check_numbers.py`
runs on, so the typography rules are identical and `**` survives, because bold
marks the number a passage argues for.
"""
import glob
import json
import os
import sys

REPOSITORY = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(REPOSITORY, "data"))
sys.path.insert(0, REPOSITORY)

import activation_sink
import buffer_role
import mechanism_discrimination
import rate_laws
from doc_check import Checker

DOCUMENT = os.path.join(REPOSITORY, "MECHANISM_EVIDENCE.md")

_SUBSTRATES = ("4OMe-BnOH", "BnOH")
_FIVE = ("C0", "C1", "C2", "C3", "C4")


def _stage_b_best(substrate):
    """The Stage B report the document quotes: the lowest uncut CV total."""
    reports = []
    for path in glob.glob(os.path.join(
            rate_laws.RATE_LAW_DIR, f"global_{substrate}_*.json")):
        with open(path) as handle:
            report = json.load(handle)
        if not report.get("cut"):
            reports.append(report)
    return min(reports, key=lambda report: report["sum"])


def _element_range(column):
    """(min, max) of a `law_scatter` column over the four element rows."""
    values = []
    for substrate in _SUBSTRATES:
        table = rate_laws.law_scatter(substrate)
        for element in rate_laws.ELEMENTS:
            value = table.loc[element, column]
            if value == value:                       # not NaN
                values.append(float(value))
    return min(values), max(values)


def main():
    doc = Checker(DOCUMENT, label="the mechanism evidence register",
                  document_label="MECHANISM_EVIDENCE.md")

    doc.section("the material")
    rows = {s: mechanism_discrimination.summary_table(s)["rows"]
            for s in _SUBSTRATES}
    parameters = {s: rate_laws.curve_parameters(s) for s in _SUBSTRATES}
    doc.claim("curves per substrate", "| live catalysed curves | "
              + " | ".join(str(len(rows[s])) for s in _SUBSTRATES) + " |")
    doc.claim("runs per substrate", "| runs | "
              + " | ".join(str(int(rows[s].experiment.nunique()))
                           for s in _SUBSTRATES) + " |")
    doc.claim("curves the form cannot hold",
              "| curves the activation-sink form cannot hold | "
              + " | ".join(str(len(parameters[s]["remaining"]))
                           for s in _SUBSTRATES) + " |")
    live = sum(len(rows[s]) for s in _SUBSTRATES)
    remaining = activation_sink.remaining_curves()
    missed = int(remaining.catalysed.sum())
    doc.claim("curves with a fit",
              f"**{live - missed} of the {live} live catalysed curves carry "
              "an activation-sink fit**")
    doc.claim("curves without one", f"and {missed} do not")

    doc.section("the reproducibility ceiling")
    replicate = rate_laws.replicate_parameter_scatter()
    doc.claim("v_act replicate SD",
              f"**{replicate.loc['v_act', 'y_sd']:.3f}**")
    doc.claim("k_act_lag replicate SD",
              f"**{replicate.loc['k_act_lag', 'y_sd']:.3f}**")
    doc.claim("the lag clock's own replicate SD",
              f"{replicate.loc['v_act', 'lag_half_s_sd']:.3f} / "
              f"{replicate.loc['k_act_lag', 'lag_half_s_sd']:.3f}")
    doc.claim("vmax's replicate SD",
              f"{replicate.loc['v_act', 'vmax_corrected_sd']:.3f}")
    scatter = mechanism_discrimination.summary_scatter("4OMe-BnOH")
    doc.claim("the summaries' replicate SDs",
              f"{scatter.loc['L', 'replicate_sd']:.3f} (level), "
              f"{scatter.loc['E', 'replicate_sd']:.3f} (early shape) and "
              f"{scatter.loc['D', 'replicate_sd']:.3f} (late shape)")
    low, high = _element_range("total_sd")
    doc.claim("the law residuals", f"residuals of {low:.3f} to {high:.3f}")
    low, high = _element_range("median_se")
    doc.claim("the parameters' own errors",
              f"own errors are {low:.3f} to {high:.3f}")
    low, high = _element_range("within_run_sd")
    doc.claim("the within-run share", f"({low:.3f} to {high:.3f})")

    doc.section("what rate laws bought")
    baselines = {s: rate_laws.law_free_baselines(s) for s in _SUBSTRATES}
    costs = {s: baselines[s]["rows"].set_index("name")["cost"]
             for s in _SUBSTRATES}
    for label, name in (
            ("the per-curve fit", "own activation-sink fit"),
            ("the per-curve quadratic", "own quadratic"),
            ("the laws, cross-validated", "Stage B best, cross-validation"),
            ("the law-free baseline", "own line and one global sink"),
            ("Stage A unchanged", "Stage A laws unchanged")):
        doc.claim(label, "| " + " | ".join(
            f"{costs[s][name]:,.0f}" for s in _SUBSTRATES) + " |")
    for substrate in _SUBSTRATES:
        doc.claim(f"{substrate}'s baseline verdict",
                  f'"{baselines[substrate]["verdict"]}"')
        best = _stage_b_best(substrate)
        verdict = rate_laws.stage_b_degeneracy(
            substrate, best["laws"], best["coefficients"])["verdict"]
        doc.claim(f"{substrate}'s Stage B degeneracy",
                  f'"{verdict.split(";")[0]}')

    doc.section("what mechanism discrimination decided")
    result = mechanism_discrimination.candidate_verdicts("4OMe-BnOH")
    tie = result["tie"]
    doc.check("no fold was skipped", result["skipped"] == 0,
              f"{result['skipped']}")
    for candidate in tie.index:
        doc.claim(f"{candidate}'s fold total and status",
                  f"| {candidate} | {tie.loc[candidate, 'total']:.1f} |")
        doc.check(f"{candidate}'s status is in the document",
                  tie.loc[candidate, "status"] in ("best", "tied", "excluded"),
                  tie.loc[candidate, "status"])
    excluded = result["verdicts"]["C3"].split(" (")[0]
    doc.claim("C3's verdict", f'"{excluded}"')
    refused = mechanism_discrimination.candidate_verdicts(
        "BnOH", candidates=_FIVE)["verdicts"].iloc[0]
    doc.claim("BnOH's refusal", f'"{refused}"')
    design = mechanism_discrimination.discrimination_table("4OMe-BnOH")
    pairs = design["pairs"]
    for pair in ("C2 vs C3", "C2 vs C4", "C2 vs C5"):
        doc.check(f"{pair} is not distinguishable",
                  pairs[pair] == "not distinguishable", pairs[pair])
    doc.check("C1 vs C4 is distinguishable",
              pairs["C1 vs C4"] == "distinguishable", pairs["C1 vs C4"])
    doc.claim("the pair that moved with the candidate set",
              f'"{pairs["C0 vs C1"]}"')

    doc.section("the live lead and the experiments")
    trend = mechanism_discrimination.residual_trends(
        "4OMe-BnOH", "C1")["verdict"]
    doc.claim("C1's residual trend", f'"{trend}"')
    crossing = buffer_role.peroxide_crossing()
    doc.claim("no run crosses the two ladders",
              f"Of {crossing['runs']} runs, {crossing['steps_buffer']} step")
    doc.claim("the peroxide arm",
              f"{crossing['steps_peroxide']} step")
    doc.claim("and none crosses",
              f"**{crossing['steps_both']} step both**")

    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
