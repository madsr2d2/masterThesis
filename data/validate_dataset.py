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
  optics      the monitoring wavelength and extinction coefficient the sheet
              declares match the ones the dataset uses
  invariants  extinction coefficient matches the substrate, [enz] is zero
              exactly when the manifest says no enzyme, concentrations are
              finite and non-negative, required fields are non-null

With --deep, additionally re-derives every concentration from the sheets'
volume tables (see recompute_concentrations.py) and checks that the compiled
values match the cuvettes that were actually measured. This is slower -- it
reopens all 98 spreadsheets -- so it is opt-in rather than the default.

Findings are graded:
  ERROR    a disagreement nothing has accounted for -- the pipeline has drifted
  WARN     a disagreement already recorded in the manifest's open_questions
           column, i.e. known and awaiting adjudication

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

    @property
    def errors(self):
        return [i for i in self.items if i["level"] == "ERROR"]

    @property
    def warnings(self):
        return [i for i in self.items if i["level"] == "WARN"]

    def report(self, quiet=False):
        """Prints the findings and returns the number of errors."""
        if not self.items and not quiet:
            print("no findings -- dataset agrees with the manifest")
        for level in ("ERROR", "WARN"):
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

        # --- invariants ----------------------------------------------------
        expected_e = SUBSTRATE_PROPERTIES.get(observed["substrate"], {}).get("e")
        actual_e = rows.e.iloc[0]
        if expected_e is not None and not _agrees(expected_e, actual_e, tolerance=1e-6):
            findings.error("invariant", number,
                           f"extinction coefficient {actual_e} does not match "
                           f"{observed['substrate']} (expected {expected_e})")

        # kinetics_io hardcodes e per SUBSTRATE in SUBSTRATE_PROPERTIES, but each
        # sheet declares its own alongside the monitoring wavelength. Where they
        # diverge the compiled concentrations are scaled wrong by exactly that
        # ratio, so this is checked against the sheet, not against the hardcode.
        # The monitoring wavelength is declared per experiment but the dataset
        # takes it, like e, from the SUBSTRATE_PROPERTIES hardcode. Until this
        # check existed the manifest's abs_nm was read only to interpolate into
        # an error message, so a new drift would have passed silently.
        declared_nm = declared.get("abs_nm")
        actual_nm = rows["abs"].iloc[0]
        if declared_nm is not None and not pd.isna(declared_nm):
            if not _agrees(declared_nm, actual_nm, tolerance=0.5):
                message = (f"manifest declares {declared_nm:g} nm but the dataset "
                           f"uses {actual_nm:g} nm")
                if "abs_nm" in disputed:
                    findings.warn("optics", number, message + "  (open question)")
                else:
                    findings.error("optics", number, message)

        declared_e = declared.get("e_declared")
        if declared_e is not None and not pd.isna(declared_e):
            if not _agrees(declared_e, actual_e, tolerance=1e-6):
                message = (f"sheet declares e = {declared_e} at "
                           f"{declared.get('abs_nm')} nm but the dataset uses "
                           f"{actual_e} (all concentrations scaled by "
                           f"{actual_e / declared_e:.3g}x)")
                if "e_declared" in disputed:
                    findings.warn("optics", number, message + "  (open question)")
                else:
                    findings.error("optics", number, message)

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
            print(f"deep check: {len(verified)} concentration column(s) confirmed from the "
                  f"volume tables; the remaining {len(unverifiable)} were serially diluted "
                  f"and {traced} cuvette values were traced back to a recorded dilution")
    if not args.quiet:
        data = pd.read_csv(args.csv)
        print(f"validating {args.csv} ({len(data)} rows, "
              f"{data.experiment.nunique()} experiments) against {args.manifest}")
    error_count = findings.report(quiet=args.quiet)
    print(f"\n{len(findings.errors)} error(s), {len(findings.warnings)} warning(s)")
    sys.exit(1 if error_count else 0)
