"""
Traces [enz] back to the weighed catalyst, the one concentration with no second
source.

recompute_concentrations checks [buf], [h2o2] and [sub] against the cuvette
volume tables, and verify_dilutions traces the serially diluted ones back
through the recorded dilution chain to weighed grams. [enz] has had neither: it
is read from a single cell of the cuvette table, and when that cell is wrong
nothing notices. Exps 79 and 80 are the demonstration -- their "[Enz] mmol/l"
column holds 0.000001, roughly 14,000x too low, and the error survived until it
was spotted by eye.

The sheets do carry a second route. Every sheet that uses enzyme prepares it in
a header block that records the molar mass, the mass weighed, the volume, and
the resulting stock concentration, then a "kuv" row giving the concentration in
the cuvette:

    g/mol   1054.29        <- the catalyst
    g       0.0236
    l       0.004
    mol     0.0000224
    mol/l   0.005596
    mmol/l  5.596183       <- stock
    kuv     0.139905  mmol/l   <- in the cuvette

That is a complete chain from a weighed mass to the number the dataset uses, and
this module checks every link of it:

  stock    the block's own arithmetic: g / (g/mol) / l == the declared mol/l
  kuv      the cuvette value: stock * V_enz / V_total == the declared kuv
  compiled the dataset agrees with kuv

Usage:
    python data/verify_enzyme.py
"""
import argparse
import contextlib
import io
import os

import numpy as np
import pandas as pd

from build_manifest import read_sheet_optics  # noqa: F401  (keeps import parity)
from kinetics_io import find_header_row

MANIFEST_PATH = "data/manifest.csv"
DATASET_PATH = "data/experiment_data.csv"
SHEET_DIR = "data/data"

# Molar mass of the catalyst, as every sheet that weighs it records.
CATALYST_MW = 1054.29
TOLERANCE = 0.02        # relative, on concentrations
LABELS = ("g/mol", "g", "l", "ml", "mol", "mol/l", "mmol/l", "kuv")


def read_enzyme_block(sheet):
    """
    Reads the enzyme preparation block a sheet declares in its header.

    Anchored on the catalyst's molar mass, 1054.29 g/mol, which appears in no
    other block: the substrate blocks carry 108.14, 138.17 or 187.03 and the
    buffer blocks their own. Anchoring on the "kuv" label instead missed the
    17 sheets of the later series (exps 135-151), which lay the block out
    differently -- volume in ml rather than l, and no kuv row at all, the
    cuvette concentration living in the table instead.

    Args:
        sheet (pd.DataFrame): A raw, header-less sheet.

    Returns:
        dict or None: the labelled values found, with "l" always in litres.
    """
    text = sheet.map(lambda v: "" if pd.isna(v) else str(v).strip().lower())
    for row in range(min(60, sheet.shape[0])):
        for column in range(sheet.shape[1] - 1):
            if text.iat[row, column] != "g/mol":
                continue
            molar_mass = pd.to_numeric(pd.Series([sheet.iat[row, column + 1]]),
                                       errors="coerce").iloc[0]
            if pd.isna(molar_mass) or abs(float(molar_mass) - CATALYST_MW) > 0.5:
                continue
            block = {}
            for scan in range(row, min(row + 10, sheet.shape[0])):
                label = text.iat[scan, column]
                if label not in LABELS or label in block:
                    continue
                found = pd.to_numeric(pd.Series([sheet.iat[scan, column + 1]]),
                                      errors="coerce").iloc[0]
                if pd.notna(found):
                    block[label] = float(found)
            if "ml" in block and "l" not in block:
                block["l"] = block.pop("ml") / 1000.0     # the later series
            return block or None
    return None


def read_enzyme_volumes(sheet):
    """
    Returns (V_enz, V_total) for the measured cuvettes, or (None, None).

    Reference rows are skipped: they are the paired no-enzyme channels, so
    their enzyme volume is zero by design.
    """
    with contextlib.redirect_stdout(io.StringIO()):
        header_row = find_header_row(sheet)
    if header_row is None:
        return None, None

    labels = []
    for column in range(sheet.shape[1]):
        parts = [str(sheet.iat[row, column])
                 for row in range(header_row, min(header_row + 2, len(sheet)))
                 if str(sheet.iat[row, column]) != "nan"]
        labels.append(" ".join(parts).strip().lower())

    enzyme_column = next((i for i, l in enumerate(labels)
                          if l.startswith("enz") and "ml" in l), None)
    total_column = next((i for i, l in enumerate(labels) if l.startswith("vol")), None)
    first_column = next((i for i, l in enumerate(labels) if l), None)
    if enzyme_column is None or total_column is None:
        return None, None

    volumes, totals = [], []
    for row in range(header_row + 1, min(header_row + 26, len(sheet))):
        label = str(sheet.iat[row, first_column]).strip().lower()
        if label.startswith(("ref", "sum")) or label == "nan":
            continue
        enzyme = pd.to_numeric(pd.Series([sheet.iat[row, enzyme_column]]),
                               errors="coerce").iloc[0]
        total = pd.to_numeric(pd.Series([sheet.iat[row, total_column]]),
                              errors="coerce").iloc[0]
        if pd.notna(enzyme) and pd.notna(total) and total > 0:
            volumes.append(float(enzyme))
            totals.append(float(total))
    if not volumes:
        return None, None
    return volumes, totals


def _off(a, b):
    """Relative difference, guarding against a zero reference."""
    if a is None or b is None:
        return None
    if abs(b) < 1e-12:
        return None if abs(a) < 1e-12 else float("inf")
    return abs(a - b) / abs(b)


def analyse(dataset_path=DATASET_PATH, manifest_path=MANIFEST_PATH,
            directory=SHEET_DIR):
    """
    Checks every link of the enzyme chain for every experiment.

    Args:
        dataset_path, manifest_path (str): Inputs.
        directory (str): Where the .xls sheets live.

    Returns:
        tuple: (findings, summary). findings is a list of
            (experiment, check, message); summary is a DataFrame.
    """
    data = pd.read_csv(dataset_path)
    manifest = pd.read_csv(manifest_path).set_index("experiment")
    findings, rows = [], []

    for number in sorted(manifest.index):
        sheet = pd.read_excel(os.path.join(directory, manifest.loc[number, "xls_file"]),
                              sheet_name="Sheet1", header=None)
        block = read_enzyme_block(sheet)
        compiled = sorted(set(np.round(data[data.experiment == number]["[enz]"], 5)))

        row = {"experiment": number, "block": bool(block), "compiled": compiled[0]
               if len(compiled) == 1 else None, "kuv": None, "stock": None,
               "from_mass": None, "from_volumes": None}
        if block is None:
            rows.append(row)
            continue

        row["kuv"] = block.get("kuv")
        row["stock"] = block.get("mmol/l")

        # link 1 -- the block's own arithmetic, from the weighed mass
        mass, volume, molar_mass = block.get("g"), block.get("l"), block.get("g/mol")
        if mass and volume and molar_mass:
            row["from_mass"] = mass / molar_mass / volume * 1000.0
            if abs(molar_mass - CATALYST_MW) > 0.5:
                findings.append((number, "enzmw",
                                 f"enzyme block declares M = {molar_mass} g/mol, "
                                 f"not the catalyst's {CATALYST_MW}"))
            drift = _off(row["from_mass"], row["stock"])
            if drift is not None and drift > TOLERANCE:
                findings.append((number, "enzstock",
                                 f"stock: {mass} g / {molar_mass} g/mol / {volume} l "
                                 f"= {row['from_mass']:.4f} mM, but the sheet declares "
                                 f"{row['stock']:.4f} mM"))

        # link 2 -- stock to cuvette, through the recorded volumes
        volumes, totals = read_enzyme_volumes(sheet)
        implied = None
        if volumes and row["stock"]:
            candidates = sorted({round(row["stock"] * v / t, 6)
                                 for v, t in zip(volumes, totals)})
            if len(candidates) == 1:
                implied = candidates[0]
                row["from_volumes"] = implied
                if row["kuv"] is not None:
                    drift = _off(implied, row["kuv"])
                    if drift is not None and drift > TOLERANCE:
                        findings.append((number, "enzkuv",
                                         f"cuvette: stock {row['stock']:.4f} mM x "
                                         f"{volumes[0]}/{totals[0]} ml = {implied:.5f} mM, "
                                         f"but the sheet's kuv says {row['kuv']:.5f} mM"))

        # the sheets of the later series carry no kuv row, so the volume route
        # is the only statement of the cuvette concentration there
        expected = row["kuv"] if row["kuv"] is not None else implied

        # link 3 -- cuvette to the compiled dataset
        if row["compiled"] is not None and expected is not None:
            drift = _off(row["compiled"], expected)
            if drift is not None and drift > max(TOLERANCE, 0.0005 / max(expected, 1e-9)):
                findings.append((number, "enzuse",
                                 f"dataset uses [enz] = {row['compiled']} mM but the "
                                 f"sheet's enzyme block gives {expected:.5f} mM"))
        rows.append(row)

    return findings, pd.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", default=DATASET_PATH)
    parser.add_argument("--manifest", default=MANIFEST_PATH)
    args = parser.parse_args()

    findings, summary = analyse(args.csv, args.manifest)
    with_block = summary[summary.block]
    print(f"{len(summary)} experiments; {len(with_block)} declare an enzyme block")
    traced = with_block[with_block.from_mass.notna()]
    print(f"{len(traced)} trace [enz] back to a weighed mass of catalyst")
    print(f"\n{len(findings)} finding(s)")
    for number, check, message in findings:
        print(f"  [{check:9s}] exp {number:<5d} {message}")
