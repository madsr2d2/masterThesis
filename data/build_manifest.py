"""
Bootstraps data/manifest.csv -- the declared ground truth for the dataset --
from the independent sources that exist outside the .xls sheets:

  1. Filenames.       They encode pH, temperature, substrate, buffer and
                      whether enzyme was present. Confirmed as ground truth by
                      the experimentalist; recovers 75-92% of fields.
  2. data/Mads/.      Hand-sorted folders: "No enzyme" states enzyme status
                      outright, "bad data"/"bad data pH ca. 11" record which
                      runs were judged unusable.
  3. The extraction.  data/kinetics_io.py, used only to fill fields the two
                      sources above are silent on, and to supply sample counts.

Where a filename and the extraction disagree, the filename wins and the
disagreement is reported for adjudication -- those rows are exactly the ones
worth a human's attention.

This is a BOOTSTRAP, not part of the normal pipeline. Run it once to propose a
manifest, hand-check the reported conflicts, then edit data/manifest.csv
directly from then on. Routine checking is data/validate_dataset.py's job.

Usage:
    python data/build_manifest.py                  # report only, writes nothing
    python data/build_manifest.py --write          # write data/manifest.csv
    python data/build_manifest.py --write -o path  # write elsewhere
"""
import argparse
import os
import re

import pandas as pd

from kinetics_io import load_experiment

MADS_DIR = "data/Mads"

# Experiments removed by clean_experiment_dataframe, with the reasons recorded
# in DATA_VERIFICATION.md. Seeded here so status/exclusion becomes declared data
# rather than a literal buried in a notebook cell.
KNOWN_EXCLUSIONS = {
    58: "reaction-direction: whole experiment runs backwards",
    72: "flat progress curves in every sample; same-day sibling 71 normal",
    77: "reaction-direction: whole experiment runs backwards",
    78: "reaction-direction: whole experiment runs backwards",
    79: "reaction-direction: whole experiment runs backwards",
    82: "flat progress curves in every sample; same-day sibling 83 normal",
    84: "hand-sorted into data/Mads/'bad data'",
}

# Substrate keys are checked in order: '4OMe-BnOH' spellings must be tested
# before the bare 'bnoh' that they all contain.
SUBSTRATE_PATTERNS = [
    ("4ome-bnoh", "4OMe-BnOH"),
    ("pmeobnoh", "4OMe-BnOH"),
    ("4-meoh-bnoh", "4OMe-BnOH"),
    ("meobnoh", "4OMe-BnOH"),
    ("4-brom-bnoh", "4Br-BnOH"),
    ("bnoh", "BnOH"),
]

# Likewise 'pyrophosphat' must precede 'phosphat', which it contains.
BUFFER_PATTERNS = [
    ("pyrophosphat", "Pyrophosphate"),
    ("carbonate", "Carbonate"),
    ("co3", "Carbonate"),
    ("boric", "Boric"),
    ("b(oh)3", "Boric"),
    ("na2hpo4", "Phosphate"),
    ("phosphat", "Phosphate"),
]


def parse_filename(filename):
    """
    Extracts whatever metadata the filename declares.

    Args:
        filename (str): The .xls filename for one experiment.

    Returns:
        dict: Any of "pH", "T", "substrate", "buffer", "has_enzyme" that the
            filename states. Keys are simply absent when it is silent.
    """
    parsed = {}
    lowered = filename.lower()

    match = re.search(r"ph[=_]?(\d+[.,]\d+|\d+)", lowered)
    if match:
        parsed["pH"] = float(match.group(1).replace(",", "."))

    match = re.search(r"_t=(\d+)", lowered)
    if match:
        parsed["T"] = float(match.group(1))

    for pattern, value in SUBSTRATE_PATTERNS:
        if pattern in lowered:
            parsed["substrate"] = value
            break

    for pattern, value in BUFFER_PATTERNS:
        if pattern in lowered:
            parsed["buffer"] = value
            break

    # "with_NO_E" lowercases to "with_no_e", which does not contain "with_e",
    # so the no-enzyme test is safe to run first.
    if "no_e" in lowered:
        parsed["has_enzyme"] = False
    elif "with_e" in lowered:
        parsed["has_enzyme"] = True

    return parsed


def scan_mads_folders(mads_dir=MADS_DIR):
    """
    Reads the hand-sorted folders under data/Mads for curation ground truth.

    Args:
        mads_dir (str): Path to the data/Mads directory.

    Returns:
        dict: {experiment_number: {"has_enzyme": bool, "folder_flag": str}}
            carrying only what the folder placement actually asserts.
    """
    claims = {}
    if not os.path.isdir(mads_dir):
        return claims

    folder_meaning = {
        "No enzyme": {"has_enzyme": False},
        "bad data": {"folder_flag": "bad data"},
        "bad data pH ca. 11": {"folder_flag": "bad data pH ca. 11"},
        "good data BnOH": {"folder_flag": "good data BnOH"},
    }

    for folder, meaning in folder_meaning.items():
        path = os.path.join(mads_dir, folder)
        if not os.path.isdir(path):
            continue
        for entry in os.listdir(path):
            match = (re.search(r"mads_t0*(\d+)[_.]", entry)
                     or re.match(r"data(\d+)\.txt$", entry)
                     or re.match(r"rate0*(\d+)\.rre$", entry)
                     or re.search(r"#(\d+)\.", entry))
            if match:
                claims.setdefault(int(match.group(1)), {}).update(meaning)
    return claims


def build(directory="data/data"):
    """
    Proposes a manifest row per experiment and collects every disagreement
    between the declared sources and the extraction.

    Args:
        directory (str): Directory holding the .txt/.xls experiment files.

    Returns:
        tuple: (pd.DataFrame manifest, list of conflict dicts)
    """
    folder_claims = scan_mads_folders()
    rows, conflicts = [], []

    txt_numbers = []
    for entry in sorted(os.listdir(directory)):
        match = re.match(r"data(\d+)\.txt$", entry)
        if match:
            txt_numbers.append(int(match.group(1)))

    for number in sorted(txt_numbers):
        experiment = load_experiment(number, directory=directory)
        if experiment is None:
            continue

        from_name = parse_filename(str(experiment["xls_file"]))
        from_folder = folder_claims.get(number, {})
        extracted = {
            "pH": experiment["pH"],
            "T": experiment["T"],
            "substrate": experiment["substrate"],
            "buffer": experiment["buffer"],
            "has_enzyme": any(s["[enz]"] for s in experiment["samples"]),
        }

        resolved, provenance, open_questions = {}, {}, []
        for field in ("pH", "T", "substrate", "buffer", "has_enzyme"):
            declared = from_folder.get(field, from_name.get(field))
            declared_by = ("folder" if field in from_folder
                           else "filename" if field in from_name else None)
            if declared is None:
                resolved[field] = extracted[field]
                provenance[field] = "extracted"
                continue

            if _agrees(declared, extracted[field]):
                resolved[field] = declared
                provenance[field] = declared_by
                continue

            # Declared and extracted disagree. Keep the CURRENT (extracted)
            # value so that adopting the manifest changes no data, and record
            # the disagreement as an open question for a human to rule on.
            # Adjudicating means editing the value here and clearing the note.
            resolved[field] = extracted[field]
            provenance[field] = "extracted (disputed)"
            open_questions.append(
                f"{field}: {declared_by} says {declared!s}, extraction says {extracted[field]!s}")
            conflicts.append({
                "experiment": number,
                "field": field,
                "declared": declared,
                "declared_by": declared_by,
                "extracted": extracted[field],
                "file": str(experiment["xls_file"]),
            })

        flag = from_folder.get("folder_flag")
        excluded = number in KNOWN_EXCLUSIONS or flag in ("bad data", "bad data pH ca. 11")
        rows.append({
            "experiment": number,
            "substrate": resolved["substrate"],
            "buffer": resolved["buffer"],
            "pH": resolved["pH"],
            "T": resolved["T"],
            "has_enzyme": bool(resolved["has_enzyme"]),
            "n_samples": len(experiment["samples"]),
            "status": "exclude" if excluded else "use",
            "exclude_reason": KNOWN_EXCLUSIONS.get(number, flag if excluded else ""),
            "provenance": ";".join(f"{k}={v}" for k, v in provenance.items()),
            "open_questions": " | ".join(open_questions),
            "xls_file": experiment["xls_file"],
            "notes": "",
        })

    return pd.DataFrame(rows), conflicts


def _agrees(a, b):
    """Compares two field values, tolerating float noise and None."""
    if a is None or b is None:
        return a is b
    if isinstance(a, bool) or isinstance(b, bool):
        return bool(a) == bool(b)
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(float(a) - float(b)) < 0.02
    return str(a) == str(b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--directory", default="data/data")
    parser.add_argument("--write", action="store_true",
                        help="Write the manifest (otherwise report only).")
    parser.add_argument("-o", "--output", default="data/manifest.csv")
    args = parser.parse_args()

    manifest, conflicts = build(args.directory)

    print(f"experiments: {len(manifest)}   "
          f"use: {(manifest.status == 'use').sum()}   "
          f"exclude: {(manifest.status == 'exclude').sum()}")
    print(f"\nconflicts between declared sources and the extraction: {len(conflicts)}")
    if conflicts:
        print(f"\n{'exp':>5} {'field':10s} {'declared':>14} {'by':9s} {'extracted':>14}  file")
        for c in conflicts:
            print(f"{c['experiment']:>5} {c['field']:10s} {str(c['declared']):>14} "
                  f"{c['declared_by']:9s} {str(c['extracted']):>14}  {c['file'][:46]}")

    if args.write:
        manifest.to_csv(args.output, index=False)
        print(f"\nwrote {args.output}")
    else:
        print("\n(report only -- pass --write to create the manifest)")
