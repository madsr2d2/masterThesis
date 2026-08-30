"""
Validates the compiled dataset against data/manifest.csv, the declared ground
truth for this dataset.

The extraction in kinetics_io.py works by scanning heterogeneous spreadsheets
for anything that looks like a pH, a temperature, a buffer name or a
concentration table. That approach cannot fail loudly: when it guesses wrong it
returns a plausible number, and historically every such error was found by
accident during unrelated analysis rather than by the pipeline itself. This
module exists to remove that failure mode. The manifest declares what each
experiment IS; this compares the extraction against it and reports every
disagreement.

Checks performed:
  coverage    every compiled experiment has a manifest row, and vice versa
  metadata    substrate / buffer / pH / T / enzyme-presence agree
  structure   sample counts agree; excluded experiments are flagged
  optics      the monitoring wavelength and extinction coefficient follow the
              substrate (BnOH 285 nm / 4OMe-BnOH 300 nm) and are uniform across
              the dataset; sheets that declare something else are reported as
              notes, since their cells are working notes, not measurements
  invariants  extinction coefficient matches the substrate, [enz] is zero
              exactly when the manifest says no enzyme, concentrations are
              finite and non-negative, required fields are non-null

With --deep, additionally re-derives every concentration from the sheets'
volume tables (recompute_concentrations.py), traces the serially diluted ones
back through the recorded dilution chain (verify_dilutions.py), and traces
[enz] back to the weighed mass of catalyst (verify_enzyme.py), and cross-checks
[enz] against the experimenter's own Rate(pH).xls analysis
(verify_rate_workbook.py). This is slower
-- it reopens all 98 spreadsheets several times -- so it is opt-in.

Findings are graded:
  ERROR    a disagreement nothing has accounted for -- the pipeline has drifted
  WARN     a disagreement already recorded in the manifest's open_questions
           column, i.e. known and awaiting adjudication
  NOTE     a deviation that has been adjudicated and recorded as a ruling;
           shown so the evidence stays visible, but not a defect

Exit status is non-zero if any ERROR is raised, so this can gate a rebuild.

Usage:
    python data/validate_dataset.py
    python data/validate_dataset.py --csv data/experiment_data.csv
    python data/validate_dataset.py --quiet     # findings only
"""
import argparse
import sys

import numpy as np
import pandas as pd

from kinetics_io import SUBSTRATE_PROPERTIES

MANIFEST_PATH = "data/manifest.csv"
DATASET_PATH = "data/experiment_data.csv"

REQUIRED_COLUMNS = ["experiment", "sample", "substrate", "abs", "e", "buffer",
                    "pH", "T", "[enz]", "[buf]", "[h2o2]", "[sub]"]
CONCENTRATION_COLUMNS = ["[enz]", "[buf]", "[h2o2]", "[sub]"]


class Findings:
    """Collects graded findings and renders them as a report."""

    def __init__(self):
        self.items = []

    def add(self, level, check, experiment, message):
        self.items.append({"level": level, "check": check,
                           "experiment": experiment, "message": message})

    def error(self, check, experiment, message):
        self.add("ERROR", check, experiment, message)

    def warn(self, check, experiment, message):
        self.add("WARN", check, experiment, message)

    def note(self, check, experiment, message):
        """Records something worth seeing that is not a defect -- an adjudicated
        deviation, say. Notes never gate a rebuild."""
        self.add("NOTE", check, experiment, message)

    @property
    def errors(self):
        return [i for i in self.items if i["level"] == "ERROR"]

    @property
    def warnings(self):
        return [i for i in self.items if i["level"] == "WARN"]

    @property
    def notes(self):
        return [i for i in self.items if i["level"] == "NOTE"]

    def report(self, quiet=False):
        """Prints the findings and returns the number of errors."""
        if not self.items and not quiet:
            print("no findings -- dataset agrees with the manifest")
        for level in ("ERROR", "WARN", "NOTE"):
            group = [i for i in self.items if i["level"] == level]
            if not group:
                continue
            print(f"\n{level}  ({len(group)})")
            for item in sorted(group, key=lambda i: (i["check"], i["experiment"] or 0)):
                where = f"exp {item['experiment']}" if item["experiment"] else "dataset"
                print(f"  [{item['check']:10s}] {where:9s} {item['message']}")
        return len(self.errors)


def _open_fields(open_questions):
    """Returns the set of field names the manifest records as unadjudicated."""
    if not isinstance(open_questions, str) or not open_questions.strip():
        return set()
    return {part.split(":", 1)[0].strip()
            for part in open_questions.split("|") if ":" in part}


def _agrees(a, b, tolerance=0.02):
    """Compares two metadata values, tolerating float noise and nulls."""
    a_null = a is None or (isinstance(a, float) and np.isnan(a))
    b_null = b is None or (isinstance(b, float) and np.isnan(b))
    if a_null or b_null:
        return a_null and b_null
    if isinstance(a, (bool, np.bool_)) or isinstance(b, (bool, np.bool_)):
        return bool(a) == bool(b)
    if isinstance(a, (int, float, np.number)) and isinstance(b, (int, float, np.number)):
        return abs(float(a) - float(b)) < tolerance
    return str(a) == str(b)


def validate(dataset_path=DATASET_PATH, manifest_path=MANIFEST_PATH):
    """
    Compares a compiled dataset against the manifest.

    Args:
        dataset_path (str): Path to the compiled per-sample CSV.
        manifest_path (str): Path to the declared per-experiment manifest.

    Returns:
        Findings: The graded findings.
    """
    findings = Findings()
    data = pd.read_csv(dataset_path)
    manifest = pd.read_csv(manifest_path).set_index("experiment")

    missing_columns = [c for c in REQUIRED_COLUMNS if c not in data.columns]
    if missing_columns:
        findings.error("schema", None, f"dataset is missing columns: {missing_columns}")
        return findings

    # --- coverage -----------------------------------------------------------
    in_data = set(data.experiment.unique())
    in_manifest = set(manifest.index)
    for number in sorted(in_data - in_manifest):
        findings.error("coverage", number, "compiled but absent from the manifest")
    for number in sorted(in_manifest - in_data):
        findings.error("coverage", number, "in the manifest but not compiled")

    # --- dataset-level optics invariant ------------------------------------
    # e and the wavelength are conventions applied per substrate. Nothing
    # asserted that they are actually uniform, which is the property the whole
    # dataset's comparability rests on: within a substrate a wrong e is one
    # global scale factor absorbed into the fitted constants, but a per-
    # experiment e would silently rescale individual runs against each other.
    for substrate, group in data.groupby("substrate"):
        for column, label in (("abs", "monitoring wavelength"),
                              ("e", "extinction coefficient")):
            values = sorted(group[column].dropna().unique())
            if len(values) > 1:
                findings.error("optics", None,
                               f"{substrate}: {label} is not uniform across the "
                               f"dataset -- {values} appear in experiments "
                               f"{sorted(group.loc[group[column] == values[-1], 'experiment'].unique())[:6]} "
                               f"and others")

    # --- per-experiment agreement ------------------------------------------
    for number in sorted(in_data & in_manifest):
        rows = data[data.experiment == number]
        declared = manifest.loc[number]
        disputed = _open_fields(declared.get("open_questions"))

        observed = {
            "substrate": rows.substrate.iloc[0],
            "buffer": rows.buffer.iloc[0],
            "pH": rows.pH.iloc[0],
            "T": rows["T"].iloc[0],
            "has_enzyme": bool((rows["[enz]"].fillna(0) > 0).any()),
        }
        for field, value in observed.items():
            expected = declared[field]
            if field == "has_enzyme":
                expected = bool(expected)
            if _agrees(expected, value):
                continue
            message = f"{field}: manifest says {expected!s}, dataset has {value!s}"
            if field in disputed:
                findings.warn("metadata", number, message + "  (open question)")
            else:
                findings.error("metadata", number, message)

        if not _agrees(declared["n_samples"], len(rows), tolerance=0.5):
            findings.error("structure", number,
                           f"manifest declares {declared['n_samples']} samples, "
                           f"dataset has {len(rows)}")

        if declared["status"] == "exclude":
            findings.warn("status", number,
                          f"manifest marks this excluded ({declared['exclude_reason']}); "
                          f"the raw compile retains it by design, so this is expected "
                          f"unless the consumer forgets to drop it")

        # Surface anything the manifest itself flags as unadjudicated. These are
        # disagreements between the manifest and an independent source (usually
        # the filename), not between the manifest and the dataset, so no other
        # check can see them -- but they are exactly what needs a human.
        questions = declared.get("open_questions")
        if isinstance(questions, str) and questions.strip():
            for question in questions.split("|"):
                findings.warn("unresolved", number, question.strip())

        # --- optics ---------------------------------------------------------
        # The monitoring wavelength is a property of the SUBSTRATE, not of the
        # individual run: BnOH is read at 285 nm and 4OMe-BnOH at 300 nm
        # throughout. And e is a conversion factor applied uniformly at analysis
        # time, not a per-experiment measurement -- so the authority for both is
        # SUBSTRATE_PROPERTIES, and a sheet's own cells are notes rather than
        # ground truth. Nine of the 98 sheets disagree, every one of them a 4OMe
        # workbook carrying the BnOH template's 285 nm; those are recorded as
        # rulings in the manifest and reported below as notes, not defects.
        # See DATA_VERIFICATION.md 2026-08-30.
        expected = SUBSTRATE_PROPERTIES.get(observed["substrate"], {})
        actual_e = rows.e.iloc[0]
        actual_nm = rows["abs"].iloc[0]
        if expected.get("e") is not None and not _agrees(expected["e"], actual_e, 1e-6):
            findings.error("optics", number,
                           f"extinction coefficient {actual_e} does not match "
                           f"{observed['substrate']} (expected {expected['e']})")
        if expected.get("abs") is not None and not _agrees(expected["abs"], actual_nm, 0.5):
            findings.error("optics", number,
                           f"monitoring wavelength {actual_nm:g} nm does not match "
                           f"{observed['substrate']} (expected {expected['abs']:g} nm)")

        # The sheet's own cells, kept visible so the deviations are not lost.
        # These come from the *_sheet columns, which a ruling never overwrites.
        declared_nm = declared.get("abs_nm_sheet")
        declared_e = declared.get("e_sheet")
        deviations = []
        if declared_nm is not None and not pd.isna(declared_nm) \
                and not _agrees(declared_nm, actual_nm, 0.5):
            deviations.append(f"{declared_nm:g} nm vs {actual_nm:g} nm")
        if declared_e is not None and not pd.isna(declared_e) \
                and not _agrees(declared_e, actual_e, 1e-6):
            deviations.append(f"e = {declared_e} vs {actual_e}")
        if deviations:
            ruled = "abs_nm=ruling" in str(declared.get("provenance") or "")
            message = ("sheet declares " + " and ".join(deviations)
                       + ("; adjudicated, see the manifest's notes" if ruled
                          else "; NOT adjudicated"))
            (findings.note if ruled else findings.warn)("optics", number, message)

        enzyme_present = bool((rows["[enz]"].fillna(0) > 0).any())
        if bool(declared["has_enzyme"]) != enzyme_present and "has_enzyme" not in disputed:
            findings.error("invariant", number,
                           "[enz] is "
                           f"{'non-zero' if enzyme_present else 'zero'} but the manifest "
                           f"says enzyme was {'present' if declared['has_enzyme'] else 'absent'}")

        for column in CONCENTRATION_COLUMNS:
            values = rows[column]
            if values.isna().any():
                findings.error("invariant", number, f"{column} contains nulls")
            elif not np.isfinite(values.astype(float)).all():
                findings.error("invariant", number, f"{column} contains non-finite values")
            elif (values < 0).any():
                findings.error("invariant", number, f"{column} contains negative values")

        for column in ("pH", "T", "substrate", "buffer"):
            if rows[column].nunique(dropna=False) > 1:
                findings.error("structure", number,
                               f"{column} varies between samples within one experiment")

    return findings


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", default=DATASET_PATH)
    parser.add_argument("--manifest", default=MANIFEST_PATH)
    parser.add_argument("--quiet", action="store_true", help="Print findings only.")
    parser.add_argument("--deep", action="store_true",
                        help="Also re-derive concentrations from the volume tables.")
    args = parser.parse_args()

    findings = validate(args.csv, args.manifest)

    if args.deep:
        from recompute_concentrations import analyse
        deep_findings, _, verified, unverifiable = analyse(args.csv, args.manifest)
        manifest = pd.read_csv(args.manifest).set_index("experiment")
        for number, check, message in deep_findings:
            accepted = ""
            if number in manifest.index:
                accepted = str(manifest.loc[number].get("accepted_deviations") or "")
            label = f"deep:{check}"[:10]
            if check in [a.strip() for a in accepted.split(";") if a.strip()]:
                findings.warn(label, number, message + "  (accepted deviation)")
            else:
                findings.error(label, number, message)
        from verify_enzyme import analyse as analyse_enzyme
        enzyme_findings, enzyme_summary = analyse_enzyme(args.csv, args.manifest)
        for number, check, message in enzyme_findings:
            accepted = ""
            if number in manifest.index:
                accepted = str(manifest.loc[number].get("accepted_deviations") or "")
            label = f"deep:{check}"[:10]
            if check in [a.strip() for a in accepted.split(";") if a.strip()]:
                findings.warn(label, number, message + "  (accepted deviation)")
            else:
                findings.error(label, number, message)

        from verify_rate_workbook import analyse as analyse_workbook
        workbook_findings, workbook_summary = analyse_workbook(args.csv, args.manifest)
        for number, check, message in workbook_findings:
            # ratemap is informational: it records which pH the experimenter's
            # own analysis used, not a disagreement about a value.
            if check == "ratemap":
                findings.note("deep:rate", number, message)
                continue
            accepted = ""
            if number in manifest.index:
                accepted = str(manifest.loc[number].get("accepted_deviations") or "")
            if check in [a.strip() for a in accepted.split(";") if a.strip()]:
                findings.warn("deep:rate", number, message + "  (accepted deviation)")
            else:
                findings.error("deep:rate", number, message)

        from verify_dilutions import analyse as analyse_dilutions
        dilution_findings, dilution_summary = analyse_dilutions(args.csv, args.manifest)
        for number, check, message in dilution_findings:
            accepted = ""
            if number in manifest.index:
                accepted = str(manifest.loc[number].get("accepted_deviations") or "")
            label = f"deep:{check}"[:10]
            if check in [a.strip() for a in accepted.split(";") if a.strip()]:
                findings.warn(label, number, message + "  (accepted deviation)")
            else:
                findings.error(label, number, message)

        if not args.quiet:
            with_blocks = dilution_summary[dilution_summary.blocks > 0]
            traced = int(with_blocks[["sub", "h2o2"]].fillna(0).to_numpy().sum())
            traced_enzyme = int(enzyme_summary.from_mass.notna().sum())
            print(f"deep check: {traced_enzyme} experiment(s) traced [enz] back to a "
                  f"weighed mass of catalyst")
            confirmed = int(workbook_summary.matched.apply(bool).sum()) if len(workbook_summary) else 0
            print(f"deep check: {confirmed} of {len(workbook_summary)} pH points in "
                  f"Rate(pH).xls confirm the dataset's [enz] independently")
            print(f"deep check: {len(verified)} concentration column(s) confirmed from the "
                  f"volume tables; the remaining {len(unverifiable)} were serially diluted "
                  f"and {traced} cuvette values were traced back to a recorded dilution")
    if not args.quiet:
        data = pd.read_csv(args.csv)
        print(f"validating {args.csv} ({len(data)} rows, "
              f"{data.experiment.nunique()} experiments) against {args.manifest}")
    error_count = findings.report(quiet=args.quiet)
    print(f"\n{len(findings.errors)} error(s), {len(findings.warnings)} warning(s), "
          f"{len(findings.notes)} note(s)")
    sys.exit(1 if error_count else 0)
