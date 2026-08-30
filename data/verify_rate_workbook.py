"""
Checks [enz] against Rate(pH).xls, a third and fully independent record.

verify_enzyme.py traces [enz] through the sheet that ran the experiment: its
weighed catalyst, its stock, its cuvette volumes. That is one document. If the
enzyme block of a sheet were itself wrong -- copied from another run, say, which
is exactly what happened to exps 57/58's substrate block -- the whole chain
would agree with itself and still be wrong.

Rate(pH).xls sits in data/Mads and is the experimenter's own pH-dependence
analysis of the 4OMe-BnOH series. It tabulates, per pH, the enzyme concentration
that was used:

    pH     [E] [mmol/l]
    5.64   0.17533
    5.87   0.272695
    6.50   0.240683
    ...

Those numbers were written down for a different purpose, at a different time,
from the analyst's side rather than the bench's. Agreement between them and the
compiled dataset is therefore worth more than either document alone.

The workbook labels its rows by pH, not by experiment number, so a row is
matched to every 4OMe with-enzyme experiment at that pH; the row passes if any
of them carries the tabulated [E]. Where the matched experiment's pH differs
from the label the difference is reported, since it says which pH the
experimenter's own analysis used.

Usage:
    python data/verify_rate_workbook.py
"""
import argparse
import os

import numpy as np
import pandas as pd

MANIFEST_PATH = "data/manifest.csv"
DATASET_PATH = "data/experiment_data.csv"
WORKBOOK_PATH = "data/Mads/Rate(pH).xls"

SUBSTRATE = "4OMe-BnOH"
PH_WINDOW = 0.05        # a row is matched to experiments within this of its label
TOLERANCE = 0.02        # relative, on [E]


def read_workbook(path=WORKBOOK_PATH):
    """
    Reads the (pH, [E]) table from Rate(pH).xls.

    Located by its header labels rather than by fixed coordinates, so a shifted
    sheet fails loudly instead of reading the wrong columns.

    Args:
        path (str): Path to the workbook.

    Returns:
        list: (pH, enzyme_mM) pairs, or [] if the table is not found.
    """
    if not os.path.exists(path):
        return []
    sheet = pd.read_excel(path, sheet_name="Sheet1", header=None)
    text = sheet.map(lambda v: "" if pd.isna(v) else str(v).strip().lower())

    header_row = ph_column = enzyme_column = None
    for row in range(min(30, sheet.shape[0])):
        columns = {text.iat[row, c]: c for c in range(sheet.shape[1])}
        if "ph" in columns and any(k.startswith("[e]") for k in columns):
            header_row = row
            ph_column = columns["ph"]
            enzyme_column = next(c for k, c in columns.items() if k.startswith("[e]"))
            break
    if header_row is None:
        return []

    pairs = []
    for row in range(header_row + 1, min(header_row + 20, sheet.shape[0])):
        pH = pd.to_numeric(pd.Series([sheet.iat[row, ph_column]]), errors="coerce").iloc[0]
        enzyme = pd.to_numeric(pd.Series([sheet.iat[row, enzyme_column]]),
                               errors="coerce").iloc[0]
        if pd.isna(pH):
            break
        if pd.notna(enzyme):
            pairs.append((float(pH), float(enzyme)))
    return pairs


def analyse(dataset_path=DATASET_PATH, manifest_path=MANIFEST_PATH,
            workbook_path=WORKBOOK_PATH):
    """
    Compares the workbook's tabulated [E] against the compiled dataset.

    Returns:
        tuple: (findings, summary). findings is a list of
            (experiment or None, check, message); summary is a DataFrame.
    """
    pairs = read_workbook(workbook_path)
    findings, rows = [], []
    if not pairs:
        findings.append((None, "ratebook",
                         f"could not read the (pH, [E]) table from {workbook_path}"))
        return findings, pd.DataFrame(rows)

    data = pd.read_csv(dataset_path)
    manifest = pd.read_csv(manifest_path).set_index("experiment")

    for pH, enzyme in pairs:
        candidates = []
        for number in manifest.index:
            if manifest.loc[number, "substrate"] != SUBSTRATE:
                continue
            if not bool(manifest.loc[number, "has_enzyme"]):
                continue
            if abs(float(manifest.loc[number, "pH"]) - pH) > PH_WINDOW:
                continue
            values = sorted(set(np.round(data[data.experiment == number]["[enz]"], 5)))
            candidates.append((number, values, float(manifest.loc[number, "pH"])))

        matched = [c for c in candidates
                   if any(abs(v - enzyme) <= max(0.0006, TOLERANCE * enzyme)
                          for v in c[1])]
        rows.append({"pH": pH, "workbook_E": enzyme,
                     "candidates": [c[0] for c in candidates],
                     "matched": [c[0] for c in matched]})

        if not candidates:
            findings.append((None, "ratebook",
                             f"Rate(pH).xls tabulates pH {pH} with [E] = {enzyme:g} mM, "
                             f"but no {SUBSTRATE} enzyme experiment sits within "
                             f"{PH_WINDOW} of that pH"))
        elif not matched:
            detail = "; ".join(f"exp {n} has {v}" for n, v, _ in candidates)
            findings.append((candidates[0][0], "ratee",
                             f"Rate(pH).xls tabulates [E] = {enzyme:g} mM at pH {pH}, "
                             f"but no experiment at that pH carries it ({detail})"))
        else:
            for number, _, experiment_pH in matched:
                if abs(experiment_pH - pH) > 1e-9:
                    findings.append((number, "ratemap",
                                     f"matched on [E] = {enzyme:g} mM, but Rate(pH).xls "
                                     f"labels the row pH {pH} where the dataset has "
                                     f"{experiment_pH} -- the experimenter's own analysis "
                                     f"used {pH}"))
    return findings, pd.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", default=DATASET_PATH)
    parser.add_argument("--manifest", default=MANIFEST_PATH)
    parser.add_argument("--workbook", default=WORKBOOK_PATH)
    args = parser.parse_args()

    findings, summary = analyse(args.csv, args.manifest, args.workbook)
    print(f"{len(summary)} tabulated pH points; "
          f"{int(summary.matched.apply(bool).sum())} confirmed against the dataset\n")
    print(summary.to_string(index=False))
    print(f"\n{len(findings)} finding(s)")
    for number, check, message in findings:
        where = f"exp {number}" if number else "workbook"
        print(f"  [{check:8s}] {where:<9s} {message}")
