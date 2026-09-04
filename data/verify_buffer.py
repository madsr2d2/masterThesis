"""
Checks [buf] against the buffer stock each experiment actually used.

[buf] is the least-verified concentration in the dataset. [sub], [h2o2] and
[enz] are all traceable to a weighed mass through a recorded dilution chain;
[buf] is not, because the buffer was made up once and used across many runs.
The extraction takes it from one of two places:

  a declared column     the sheet carries a '[buffer] mmol/l' column, and the
                        value is read straight off it
  the volume fallback   no such column, so [buf] is computed as
                        (V_buf / 2 ml) * 100 mM -- which silently hardcodes
                        both a 2 ml cuvette and a 0.1 M stock

Thirty-three experiments rest on that fallback: the phosphate campaign of
2010-04-08 to 2010-05-06, which predates the '0.1M' filename convention that
starts at t041. I and [HOO-] scale linearly with any error in it.

This module verifies the fallback rather than trusting it, in three ways.

**The stock is recoverable wherever [buf] is declared.** Dividing the declared
concentration by the volume ratio gives the stock, and it must come out the same
for every cuvette in an experiment -- the buffer came from one bottle. A
disagreement means the declared column is not what it claims to be.

**Where the filename states a molarity, the recovered stock must match it.**
Eighteen of the twenty-one experiments whose filename says '0.1M' recover
exactly 100 mM. The three that do not recover exactly 33 mM, which is reported
rather than assumed away.

**The fallback's 0.1 M is evidenced, not assumed, for one experiment.** Exp 14's
workbook carries the experimenter's own arithmetic in a scratch column beside
the cuvette table -- in the .ods copy under data/Mads/Variable Temperature it is
still a live formula:

    P27  =(C20/H20)*0.1   = 0.08     C = Buf [ml], H = Vol [ml]
    P28  =(C21/H21)*0.1   = 0.07
    P29  =(C22/H22)*0.1   = 0.06
    P30  =(C23/H23)*0.1   = 0.05

That is (V_buf / V_tot) x 0.1 mol/L, which is the fallback exactly. Exp 14 sits
in the middle of the campaign (2010-04-18) and twenty-one of the thirty-three
share its cuvette recipe to the millilitre, so the check pins the one piece of
direct evidence the campaign has. See DATA_VERIFICATION.md 2026-08-30.

Usage:
    python data/verify_buffer.py
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

MANIFEST_PATH = "data/manifest.csv"
DATASET_PATH = "data/experiment_data.csv"
SHEET_DIR = "data/data"

TOLERANCE = 0.02              # relative, on a recovered stock
FALLBACK_STOCK_mM = 100.0     # what the extraction assumes where nothing is declared
CUVETTE_SPREAD = 1e-6         # a stock must be identical across cuvettes

# The experimenter's own buffer arithmetic, and where to find it. Keyed by
# experiment -> (sheet, first row (1-based), column (0-based), stock in mol/L).
SHEET_ARITHMETIC = {14: ("Sheet1", 27, 15, 0.1)}

# NOT build_manifest.FILENAME_MOLARITY, a looser pattern for the same job --
# this one anchors on the separators either side. They shared the bare name
# until 2026-09-04 and were never the same regex.
BUFFER_FILENAME_MOLARITY = re.compile(r"[_ ](\d+(?:[.,]\d+)?)\s*M[_.]", re.IGNORECASE)

# How far down the sheet the header's own buffer block can sit.
HEADER_DEPTH = 14


def header_stock(sheet):
    """
    Returns (buffer name, stock in mol/L) from the sheet's own header block, or
    (None, None).

    Twenty-one sheets state the buffer outright -- a name beside a 'Buffer'
    label and a molarity under it:

        N3 Buffer   O3 hexametaphosphate (pyrophosphate)
        N4 mol/l    O4 0.033

    This is the experimenter recording what went into the cuvette, and the
    sheet's [buffer] column is computed from it. It outranks the filename,
    which is copied forward between runs and only partly updated -- exps 75, 76
    and 78 all say '0.1M' in the filename over a sheet that says 0.033.
    """
    for row in range(min(HEADER_DEPTH, len(sheet))):
        for column in range(sheet.shape[1]):
            if str(sheet.iat[row, column]).strip().lower() != "buffer":
                continue
            name = None
            for right in range(column + 1, min(column + 3, sheet.shape[1])):
                value = sheet.iat[row, right]
                if str(value) != "nan" and not isinstance(value, (int, float)):
                    name = str(value).strip()
                    break
            stock = None
            for below in range(row + 1, min(row + 4, len(sheet))):
                if str(sheet.iat[below, column]).strip().lower() != "mol/l":
                    continue
                for right in range(column + 1, min(column + 3, sheet.shape[1])):
                    try:
                        candidate = float(sheet.iat[below, right])
                    except (TypeError, ValueError):
                        continue
                    if np.isfinite(candidate) and 0 < candidate < 2:
                        stock = candidate
                        break
            if name or stock:
                return name, stock
    return None, None


def filename_stock_mM(filename):
    """Returns the stock a filename declares, in mM, or None."""
    match = BUFFER_FILENAME_MOLARITY.search(str(filename))
    if match is None:
        return None
    return float(match.group(1).replace(",", ".")) * 1000.0


def recipe(number, manifest):
    """
    Returns (buffer volumes, total volumes, declared [buf]) for one experiment's
    measured cuvettes, read from its sheet. Any element may be empty.
    """
    from build_dossier import (column_labels, header_index, measured_rows,
                               read_sheet, recipe_rows)

    filename = str(manifest.loc[number, "xls_file"])
    if not os.path.exists(os.path.join(SHEET_DIR, filename)):
        return [], [], []
    sheet = read_sheet(filename)
    header = header_index(sheet)
    if header is None:
        return [], [], []
    labels = column_labels(sheet, header)

    buffer_key = total_key = declared_key = None
    for label in labels:
        text = str(label).lower().replace(" ", "")
        if "buf" in text and "ml" in text and buffer_key is None:
            buffer_key = label
        elif "vol" in text and "ml" in text and total_key is None:
            total_key = label
        elif ("buffer" in text or "[buf" in text) and "ml" not in text:
            declared_key = label
    if buffer_key is None or total_key is None:
        return [], [], []

    volumes, totals, declared = [], [], []
    for row in measured_rows(recipe_rows(sheet, header, labels)):
        try:
            buffer_volume = float(row.get(buffer_key))
            total = float(row.get(total_key))
        except (TypeError, ValueError):
            continue
        # A trailing 'Sum:' row can carry a column total in both cells and so
        # survive recipe_rows' blank-total guard. Exp 128 sums to 6.88 ml of
        # buffer in a 0.32 ml cuvette, which no real cuvette can be.
        if not (np.isfinite(buffer_volume) and np.isfinite(total)
                and buffer_volume > 0 and total > 0
                and buffer_volume <= total):
            continue
        volumes.append(buffer_volume)
        totals.append(total)
        value = np.nan
        if declared_key is not None:
            try:
                value = float(row.get(declared_key))
            except (TypeError, ValueError):
                pass
        declared.append(value)
    return volumes, totals, declared


def read_arithmetic(number, manifest):
    """
    Returns the stock in mol/L implied by the experimenter's own scratch
    calculation, or None if the expected cells do not hold it.
    """
    if number not in SHEET_ARITHMETIC:
        return None
    sheet_name, first_row, column, _ = SHEET_ARITHMETIC[number]
    filename = str(manifest.loc[number, "xls_file"])
    path = os.path.join(SHEET_DIR, filename)
    if not os.path.exists(path):
        return None
    grid = pd.read_excel(path, sheet_name=sheet_name, header=None)
    volumes, totals, _ = recipe(number, manifest)
    if not volumes or column >= grid.shape[1]:
        return None

    ratios = np.array(volumes) / np.array(totals)
    values = pd.to_numeric(
        grid.iloc[first_row - 1:first_row - 1 + len(ratios), column],
        errors="coerce").values
    usable = np.isfinite(values)
    if usable.sum() < 3:
        return None
    stocks = values[usable] / ratios[:len(values)][usable]
    if stocks.std() > CUVETTE_SPREAD * max(abs(stocks.mean()), 1.0):
        return None
    return float(stocks.mean())


def analyse(dataset_path=DATASET_PATH, manifest_path=MANIFEST_PATH):
    """
    Returns:
        tuple: (findings, summary). findings is a list of
            (experiment or None, check, message); summary is a DataFrame.
    """
    from build_dossier import read_sheet

    data = pd.read_csv(dataset_path)
    manifest = pd.read_csv(manifest_path).set_index("experiment")
    findings, rows = [], []

    declared_stock = {}
    for number in manifest.index:
        filename = str(manifest.loc[number, "xls_file"])
        if os.path.exists(os.path.join(SHEET_DIR, filename)):
            declared_stock[number] = header_stock(read_sheet(filename))

    for number in manifest.index:
        volumes, totals, declared = recipe(number, manifest)
        if not volumes:
            continue
        ratios = np.array(volumes) / np.array(totals)
        declared = np.array(declared, dtype=float)
        row = {"experiment": number, "cuvettes": len(volumes),
               "provenance": manifest.loc[number, "buf_provenance"],
               "stock_mM": np.nan, "source": "", "buffer_name": None}

        # The volume fallback divides by a hardcoded 2 ml cuvette, so a sheet
        # whose cuvettes are not 2 ml silently gets the wrong [buf].
        if manifest.loc[number, "buf_provenance"] in ("assumed", "filename"):
            odd = sorted({t for t in totals if abs(t - 2.0) > 1e-6})
            if odd:
                findings.append((number, "bufvolume",
                                 f"[buf] here comes from the volume fallback, which "
                                 f"assumes a 2 ml cuvette, but the sheet's cuvettes are "
                                 f"{', '.join(f'{t:g}' for t in odd)} ml -- so every "
                                 f"[buf] in this experiment is off by that ratio"))

        # 1. the stock recovered from a declared column must be one number
        usable = np.isfinite(declared) & (declared > 0)
        if usable.any():
            stocks = declared[usable] / ratios[usable]
            spread = float(stocks.max() - stocks.min())
            row["stock_mM"] = float(np.median(stocks))
            row["source"] = "declared"
            if spread > CUVETTE_SPREAD * max(abs(stocks.mean()), 1.0):
                findings.append((number, "bufstock",
                                 f"the declared [buf] column implies a buffer stock that "
                                 f"changes between cuvettes, {stocks.min():.4g} to "
                                 f"{stocks.max():.4g} mM. One experiment draws from one "
                                 f"bottle, so the column is not [buf] = stock * "
                                 f"V_buf/V_total"))

            # 2. where the sheet states the stock outright, the [buf] column
            # must agree with it -- two places in one document, so a
            # disagreement means one of them was left behind by an edit.
            recovered = float(np.median(stocks))
            name, stated_stock = declared_stock.get(number, (None, None))
            if stated_stock is not None:
                row["buffer_name"] = name
                if abs(recovered - stated_stock * 1000.0) > TOLERANCE * recovered:
                    findings.append((number, "bufdeclared",
                                     f"the sheet's header states a {stated_stock:g} M "
                                     f"buffer but its own [buf] column implies "
                                     f"{recovered:.4g} mM"))

            # 3. and the filename, which ranks below both
            stated = filename_stock_mM(manifest.loc[number, "xls_file"])
            if stated is not None and abs(recovered - stated) > TOLERANCE * stated:
                if stated_stock is not None:
                    findings.append((number, "bufstale",
                                     f"the filename says {stated / 1000:g} M, but the "
                                     f"sheet states {stated_stock:g} M in its header and "
                                     f"computes [buf] from it. The filename is stale; the "
                                     f"sheet is the measurement"
                                     + (f" (buffer: {name})" if name else "")))
                else:
                    findings.append((number, "buflabel",
                                     f"the filename declares a {stated / 1000:g} M buffer "
                                     f"but the sheet's own [buf] column implies "
                                     f"{recovered:.4g} mM, and the sheet states no stock "
                                     f"in its header to adjudicate between them"))

        # 3. the fallback's 0.1 M must still be evidenced where it can be
        if number in SHEET_ARITHMETIC:
            expected = SHEET_ARITHMETIC[number][3]
            recovered = read_arithmetic(number, manifest)
            if recovered is None:
                findings.append((number, "bufarith",
                                 f"the experimenter's own buffer arithmetic is no longer "
                                 f"readable at {SHEET_ARITHMETIC[number][0]} "
                                 f"R{SHEET_ARITHMETIC[number][1]}"
                                 f"C{SHEET_ARITHMETIC[number][2]}. It is the only direct "
                                 f"evidence for the 0.1 M stock behind 33 experiments"))
            elif abs(recovered - expected) > TOLERANCE * expected:
                findings.append((number, "bufarith",
                                 f"the experimenter's own buffer arithmetic implies a "
                                 f"{recovered:g} mol/L stock, not the {expected:g} mol/L "
                                 f"the extraction assumes for this campaign"))
            else:
                row["source"] = "sheet-arithmetic"
                row["stock_mM"] = recovered * 1000.0

        # Where no stock is recoverable, fall back to the campaign's 0.1 M so
        # the compiled column is still checked against the sheet's volumes.
        # For those experiments this is a consistency check on the CSV, not
        # independent evidence: the extraction assumed the same number.
        if not np.isfinite(row["stock_mM"]):
            row["stock_mM"] = FALLBACK_STOCK_mM
            row["source"] = "fallback"

        # 4. every compiled [buf] must be one of the values that stock gives
        # when diluted by a cuvette's own volumes. Comparing against the SET
        # rather than pairwise, because a sheet lays out twice as many rows as
        # the run has samples -- the second half is the reference channel -- so
        # the rows do not line up with the measured cuvettes.
        compiled = data[data.experiment == number].sort_values("sample")["[buf]"].values
        expected = np.unique(np.round(row["stock_mM"] * ratios, 6))
        # sheet-text experiments are buffer titrations that vary the STOCK at a
        # fixed volume -- exps 32 and 34-37 step 0.1/0.2/0.3/0.4 M through
        # identical cuvettes. One stock per experiment is the wrong model for
        # them, so the single-stock comparison does not apply.
        if manifest.loc[number, "buf_provenance"] == "sheet-text":
            expected = np.array([])
        if len(compiled) and len(expected):
            bad = [i for i, value in enumerate(compiled)
                   if not np.any(np.abs(expected - value)
                                 <= np.maximum(0.01, TOLERANCE * np.abs(expected)))]
            if bad:
                where = ", ".join(f"cuvette {i + 1}: {compiled[i]:.4g}" for i in bad[:4])
                findings.append((number, "bufcompiled",
                                 f"the compiled [buf] is not what a "
                                 f"{row['stock_mM']:.4g} mM stock gives when diluted by "
                                 f"the sheet's own cuvette volumes "
                                 f"({{{', '.join(f'{v:.4g}' for v in expected)}}}): "
                                 f"{where}"))
        rows.append(row)

    return findings, pd.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", default=DATASET_PATH)
    parser.add_argument("--manifest", default=MANIFEST_PATH)
    args = parser.parse_args()

    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    problems, summary = analyse(args.csv, args.manifest)
    for number, check, message in problems:
        label = f"exp {number}" if number is not None else "dataset"
        print(f"  [{check:12}] {label:9} {message}")

    recovered = summary[summary["stock_mM"].notna()]
    print(f"\n{len(recovered)} experiment(s) with a recoverable buffer stock:")
    for stock, group in recovered.groupby(recovered["stock_mM"].round(1)):
        print(f"   {stock:8.1f} mM   {len(group):2} experiment(s)   "
              f"[{', '.join(sorted(set(group['source'])))}]")
    print(f"\n{len(problems)} problem(s)")
    raise SystemExit(1 if problems else 0)
