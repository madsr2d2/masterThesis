"""
Checks the compiled substrate concentrations against the instrument's own
export header.

Every other concentration check in this repository reads the spreadsheet. The
volume route (recompute_concentrations.py), the dilution chain
(verify_dilutions.py) and the buffer recovery (verify_buffer.py) all re-derive
a cuvette's concentration from cells in the same workbook the value was
extracted from. They catch arithmetic mistakes; they cannot catch a workbook
that was copied forward from the previous run and only partly updated, because
a copied workbook is perfectly self-consistent.

The .txt exports carry a source that is not the workbook. The Evolution 600
writes a "Substrate Conc." line per sample, typed into the instrument's method
at the bench, and it is exported alongside the readings:

    Substrate Conc.  7.3100 mmol/L   Substrate Conc.  3.6550 mmol/L   ...

75 of the 100 exports carry it, and in 60 it agrees with the compiled [sub] to
within 2%. That is 60 experiments whose substrate concentrations now rest on
two independent records instead of one.

WHAT THE HEADER IS NOT. It is not an authority. It is an operator-typed method
field carried between runs, and it goes stale exactly the way a filename does
-- which is why the precedence rule is the same one CLAUDE.md already states
for filenames: THE SHEET WINS. Exp 72 is the proof. Its header repeats exp
71's first three values and omits the fourth, while its own workbook and exp
71's workbook compile different ladders; a header that can be a truncated copy
of the previous run's is not evidence about this run.

So a disagreement here is a tripwire, not a verdict: it says one of the two
records is stale and the experiment needs a look. That is worth having --
running it for the first time on 2026-08-31 found exps 69 and 70, whose
headers carry exp 68's ladder and which nothing in the manifest had recorded.

Checks performed:

  instconc    the compiled [sub] disagrees with the exported header and
              nothing accounts for it -- graded as an error
  instopen    it disagrees for a reason recorded but not yet explained --
              graded as a warning (exp 36 is the only one)
  instnote    it disagrees for a reason already adjudicated -- graded as a
              note, so the evidence stays visible

Reports only; changes nothing.

Usage:
    python data/verify_instrument.py
    python data/verify_instrument.py --verbose
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

EXPORT_DIR = "data/data"
DATASET_PATH = "data/experiment_data.csv"
MANIFEST_PATH = "data/manifest.csv"

# A cuvette agrees if it is within this fraction of the exported value. The
# header is written to four decimal places and the compiled value to three, so
# the floor is rounding, not measurement.
CONCENTRATION_TOLERANCE = 0.02

# The unit is written both "mmol/L" and "mmol/l" across the archive, so the
# match is case-insensitive. It was not at first, which silently dropped 40 of
# the 80 headers and made this check look half as useful as it is.
HEADER_PATTERN = re.compile(r"([-\d.]+)\s*mmol/l", re.IGNORECASE)

# The disagreements adjudicated on 2026-08-31, with the evidence. All are ruled
# to the SHEET; none changes a compiled number. Recorded here rather than in
# build_manifest.RULINGS because no manifest field is in dispute -- the
# manifest and the dataset already agree, and it is the instrument's header
# that is wrong. DATA_VERIFICATION.md carries the same ruling with the working.
ADJUDICATED = {
    69: "header is exp 68's ladder verbatim; exps 69-71 are one low-ladder set "
        "run 6/8/2010 and exp 71's header was updated while 69 and 70's were "
        "not. Ruled to the sheet.",
    70: "header is exp 68's ladder verbatim; see exp 69. Ruled to the sheet.",
    72: "header is exp 71's first three values with the fourth missing -- a "
        "truncated copy of the previous run. Ruled to the sheet.",
    83: "header is exp 62's ladder (0.60/2.98/5.37/7.75), a template block that "
        "also appears stale inside exp 66's workbook. Ruled to the sheet.",
    57: "samples 2-4 agree; sample 1's header value 1.5286 belongs to no "
        "dilution in the sheet's own series. Ruled to the sheet.",
    58: "as exp 57, and the header is exp 57's verbatim. Ruled to the sheet.",
    82: "sample 1 agrees; samples 2-4 disagree by 4-12%. Header and sheet "
        "hold different dilution series off the same starting value -- "
        "header 1, 4/5, 2/3, 4/9 against the sheet's 1, 5/6, 5/7, 1/2. "
        "Ruled to the sheet, whose volumes produce its series.",
    84: "single cuvette; the header disagrees by 29% and the substrate itself "
        "was ruled 4OMe-BnOH on 2026-08-30 against a filename saying BnOH. "
        "Ruled to the sheet.",
    127: "header holds 102.093 mmol/L in two cuvettes and 1.276 in two, against "
         "a constant compiled 9.47 -- the field was reused for something that "
         "is not the substrate. Not a concentration conflict.",
    129: "as exp 127, at 195.882 / 3.879 mmol/L.",
    130: "as exp 127; header is exp 129's verbatim.",
    131: "as exp 127; header is exp 129's verbatim.",
    30: "sample 3's header reads 2856.0000 mmol/l where the sheet has 0.2856 "
        "-- a misplaced decimal typed at the instrument. The other three "
        "cuvettes agree to four figures. Ruled to the sheet.",
    40: "header is exp 39's ladder verbatim and exp 40 compiles exactly half "
        "of it at every rung, which is the dilution exp 40 adds. Ruled to the "
        "sheet.",
}

# Disagreements that are real and NOT yet explained. These grade as warnings,
# not errors: the sheet still wins by default, but nothing here accounts for
# the gap and it is recorded rather than waved through. Moving an experiment
# out of here means either finding the cause (move it to ADJUDICATED with the
# evidence) or changing a number -- never deleting the line.
OPEN_QUESTIONS = {
    36: "header reads 59.8997 mmol/l in all four cuvettes against a compiled "
        "57.90, a uniform 3.5% gap. It is not stale -- exps 34, 35 and 41 all "
        "declare something else -- so neither record is obviously the copy. "
        "Raised 2026-08-31; see DATA_VERIFICATION.md.",
}


def export_path(experiment, directory=EXPORT_DIR):
    """The .txt export for an experiment, or None when it has none."""
    path = os.path.join(directory, f"data{experiment}.txt")
    return path if os.path.exists(path) else None


def declared_concentrations(path):
    """
    The substrate concentrations the instrument's export header declares.

    Returns them in sample order, or None when the export carries no such
    header -- 20 of the 100 do not.
    """
    with open(path, encoding="latin-1", errors="replace") as handle:
        for line in handle:
            if "Substrate Conc" in line:
                return [float(v) for v in HEADER_PATTERN.findall(line)]
    return None


def compare_header(experiment, compiled, directory=EXPORT_DIR):
    """
    Compares one experiment's compiled [sub] against its export header.

    Returns (status, message) where status is "absent", "agrees" or "differs".
    """
    path = export_path(experiment, directory)
    if path is None:
        return "absent", "no .txt export"
    declared = declared_concentrations(path)
    if not declared:
        return "absent", "export carries no Substrate Conc. header"

    compiled = np.asarray(compiled, dtype=float)
    declared = np.asarray(declared, dtype=float)

    # A header may cover fewer cuvettes than were run -- the operator stopped
    # typing. That is not a disagreement about the cuvettes it does cover, so
    # compare the common prefix and mention the shortfall only if the values
    # agree; if they disagree the values are the story, not the length.
    shared = min(len(declared), len(compiled))
    truncated = len(declared) < len(compiled)
    declared, compiled_shared = declared[:shared], compiled[:shared]
    relative = np.abs(declared - compiled_shared) / np.maximum(np.abs(compiled_shared), 1e-12)
    off = np.nonzero(relative > CONCENTRATION_TOLERANCE)[0]
    if not len(off):
        note = (f" (header covers {shared} of {len(compiled)} cuvettes)"
                if truncated else "")
        return "agrees", (f"{shared} cuvettes agree within "
                          f"{CONCENTRATION_TOLERANCE:.0%}{note}")
    detail = ", ".join(
        f"sample {i + 1}: header {declared[i]:.4g} vs dataset "
        f"{compiled_shared[i]:.4g}" for i in off)
    return "differs", f"{len(off)} of {shared} cuvettes disagree -- {detail}"


def analyse(dataset_path=DATASET_PATH, manifest_path=MANIFEST_PATH,
            directory=EXPORT_DIR):
    """
    Compares every compiled experiment against its export header.

    Returns (findings, summary). Findings are (experiment, check, message);
    the check is "instnote" where ADJUDICATED explains the disagreement and
    "instconc" where nothing does.
    """
    data = pd.read_csv(dataset_path)
    findings = []
    summary = {"checked": 0, "agreeing": 0, "differing": 0, "absent": 0}
    for experiment in sorted(data.experiment.unique()):
        rows = data[data.experiment == experiment].sort_values("sample")
        status, message = compare_header(int(experiment),
                                         rows["[sub]"].to_numpy(dtype=float),
                                         directory)
        if status == "absent":
            summary["absent"] += 1
            continue
        summary["checked"] += 1
        if status == "agrees":
            summary["agreeing"] += 1
            continue
        summary["differing"] += 1
        number = int(experiment)
        if number in ADJUDICATED:
            findings.append((number, "instnote",
                             f"{message}  ({ADJUDICATED[number]})"))
        elif number in OPEN_QUESTIONS:
            findings.append((number, "instopen",
                             f"{message}  ({OPEN_QUESTIONS[number]})"))
        else:
            findings.append((number, "instconc", message))
    return findings, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--csv", default=DATASET_PATH)
    parser.add_argument("--manifest", default=MANIFEST_PATH)
    parser.add_argument("--verbose", action="store_true",
                        help="also list the experiments that agree")
    arguments = parser.parse_args()

    findings, summary = analyse(arguments.csv, arguments.manifest)
    print(f"{summary['checked']} experiments carry a Substrate Conc. header "
          f"({summary['absent']} do not)")
    print(f"{summary['agreeing']} agree with the compiled [sub] within "
          f"{CONCENTRATION_TOLERANCE:.0%}, {summary['differing']} do not\n")
    marks = {"instnote": "NOTE ", "instopen": "WARN ", "instconc": "OPEN "}
    unadjudicated = 0
    for experiment, check, message in sorted(findings):
        unadjudicated += check == "instconc"
        print(f"  {marks[check]} exp {experiment:>3}  {message}")
    if unadjudicated:
        print(f"\n{unadjudicated} disagreement(s) not yet adjudicated. The "
              f"sheet wins by default -- see this module's docstring -- but "
              f"each needs recording.")
    return 1 if unadjudicated else 0


if __name__ == "__main__":
    raise SystemExit(main())
