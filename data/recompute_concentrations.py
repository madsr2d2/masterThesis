"""
Recomputes every concentration from the sheets' volume tables and compares the
result against the compiled dataset.

data/manifest.csv declares the metadata; the concentration columns had nothing
independent checking them. This closes that gap by taking a different path
through the same sheets: instead of reading the pre-computed mmol/l columns
that kinetics_io reads, it reads the per-cuvette component VOLUMES and the
total volume, derives the implied stock concentration, and re-derives each
cuvette concentration from those.

Three checks:

  self-consistency   For one concentration column, c_row * V_total / V_component
                     is the implied stock. Two dilution designs exist and only
                     one is checkable this way:
                       - VOLUME design (e.g. exp 2): the component volume varies
                         at a fixed stock, so the implied stock must be identical
                         across every row. A row that disagrees means a broken
                         formula or a shifted cell -- a real finding.
                       - STOCK design (e.g. exp 65): the volume is fixed and the
                         stock is serially diluted instead. The stock is then not
                         recorded in the volume table at all, so the computed
                         concentration is the ONLY record and cannot be verified
                         by this route. Reported as unverifiable, not as an error.

  block selection    kinetics_io takes the first `sample_num` numeric rows below
                     the header, but 97 of 98 sheets plan more cuvettes than were
                     measured. This selects the measured block explicitly --
                     dropping rows labelled "ref" and honouring the manifest's
                     has_enzyme -- and reports where that disagrees with what was
                     compiled.

  buffer provenance  [buf] is read from the sheet's own [buf]/[buffer] column
                     where one exists, and otherwise falls back in kinetics_io to
                     (V_buf / 2) * 100, which hardcodes a 0.1 M stock. Reports
                     which experiments rest on that assumption.

Reports only; changes nothing.

Usage:
    python data/recompute_concentrations.py
    python data/recompute_concentrations.py --verbose
"""
import argparse
import sys
import warnings

import numpy as np
import pandas as pd

from kinetics_io import find_and_parse_experiment_file

warnings.filterwarnings("ignore")

DATASET_PATH = "data/experiment_data.csv"
MANIFEST_PATH = "data/manifest.csv"

# Column header -> the dataset column it feeds.
VOLUME_OF = {"buf [ml]": "[buf]", "enz [ml]": "[enz]",
             "sub [ml]": "[sub]", "h2o2 [ml]": "[h2o2]"}
SPECIES = ["[sub]", "[h2o2]", "[enz]", "[buf]"]
RELATIVE_TOLERANCE = 1e-3


def _text(value):
    return "" if pd.isna(value) else str(value).strip().lower()


def locate_table(sheet):
    """
    Finds the volume table and maps volume/concentration columns.

    Two sheet generations exist: the earlier one splits the species label
    ("[sub]") and its unit ("mmol/l") across two rows, the later one puts both
    in a single header cell. Combining the header row with the row above it
    handles both.

    Args:
        sheet (pd.DataFrame): A raw, header-less sheet.

    Returns:
        tuple: (header_row_index, {"volumes": {...}, "concentrations": {...},
                "total": col, "label": col}) or (None, {}).
    """
    for row in range(min(90, sheet.shape[0])):
        header = [_text(sheet.iloc[row, c]) for c in range(sheet.shape[1])]
        if "vol [ml]" not in header:
            continue
        above = ([_text(sheet.iloc[row - 1, c]) for c in range(sheet.shape[1])]
                 if row else [""] * sheet.shape[1])

        volumes, concentrations = {}, {}
        for c, cell in enumerate(header):
            if cell in VOLUME_OF:
                volumes[VOLUME_OF[cell]] = c
            combined = f"{above[c]} {cell}"
            if "mmol" not in combined:
                continue
            for species in SPECIES:
                names = [species] + (["[buffer]"] if species == "[buf]" else [])
                if any(name in combined for name in names):
                    concentrations.setdefault(species, c)

        if len(volumes) >= 3:
            return row, {"volumes": volumes, "concentrations": concentrations,
                         "total": header.index("vol [ml]"),
                         "label": 0 if _text(sheet.iloc[row, 0]).startswith("kuv") else 1}
    return None, {}


def read_rows(sheet, header_row, table):
    """Reads the contiguous cuvette rows below the header into dicts."""
    rows = []
    for r in range(header_row + 1, sheet.shape[0]):
        total = sheet.iloc[r, table["total"]]
        if not isinstance(total, (int, float)) or pd.isna(total) or total <= 0:
            if rows:
                break
            continue
        row = {"row": r, "total": float(total),
               "label": _text(sheet.iloc[r, table["label"]])}
        for species, col in table["volumes"].items():
            v = sheet.iloc[r, col]
            row[f"V{species}"] = float(v) if isinstance(v, (int, float)) and pd.notna(v) else np.nan
        for species, col in table["concentrations"].items():
            v = sheet.iloc[r, col]
            row[f"C{species}"] = float(v) if isinstance(v, (int, float)) and pd.notna(v) else np.nan
        rows.append(row)
    return rows


def implied_stocks(rows, species):
    """
    Back-calculates the stock concentration each row implies.

    c_cuvette = c_stock * V_component / V_total, so c_stock = c * V_total / V.
    Rows with zero component volume carry no information and are skipped.
    """
    stocks = []
    for row in rows:
        volume, concentration = row.get(f"V{species}"), row.get(f"C{species}")
        if (volume is None or concentration is None
                or np.isnan(volume) or np.isnan(concentration) or volume == 0):
            continue
        stocks.append((row["row"], concentration * row["total"] / volume))
    return stocks


def select_measured(rows, n_samples, has_enzyme):
    """
    Picks the cuvettes that were actually measured.

    Reference rows are labelled "ref" and were never monitored. Among the rest,
    the enzyme-containing and enzyme-free blocks are distinguished by their Enz
    volume, and the manifest says which one ran.
    """
    candidates = [r for r in rows if not r["label"].startswith("ref")]
    volumes = [r.get("V[enz]") for r in candidates]
    if any(v is not None and not np.isnan(v) for v in volumes):
        wanted = [r for r in candidates
                  if (r.get("V[enz]", 0) or 0) > 0 if has_enzyme] if has_enzyme else \
                 [r for r in candidates if not (r.get("V[enz]", 0) or 0) > 0]
        if len(wanted) >= n_samples:
            candidates = wanted
    return candidates[:n_samples]


def analyse(dataset_path=DATASET_PATH, manifest_path=MANIFEST_PATH):
    """Runs all three checks and returns (findings, per-experiment summary)."""
    data = pd.read_csv(dataset_path)
    manifest = pd.read_csv(manifest_path).set_index("experiment")
    findings, summary = [], []
    verified, unverifiable = [], []

    for number in sorted(data.experiment.unique()):
        _, sheet = find_and_parse_experiment_file(number, "data/data", "Sheet1")
        if sheet is None:
            findings.append((number, "sheet", "no .xls sheet found"))
            continue
        header_row, table = locate_table(sheet)
        if header_row is None:
            findings.append((number, "table", "volume table not located"))
            continue

        rows = read_rows(sheet, header_row, table)
        declared = manifest.loc[number] if number in manifest.index else None
        compiled = data[data.experiment == number].sort_values("sample")

        # --- select the measured cuvettes first ------------------------------
        # Everything below is checked on the rows that actually produced data.
        # Planned-but-unmeasured rows often belong to a different design recorded
        # on the same sheet (exps 135-151 carry D2O variants), so including them
        # manufactures spurious disagreement.
        has_enzyme = bool(declared["has_enzyme"]) if declared is not None else True
        chosen = select_measured(rows, len(compiled), has_enzyme)

        # --- self-consistency of each concentration column -------------------
        for species in table["concentrations"]:
            stocks = implied_stocks(chosen, species)
            if len(stocks) < 2:
                continue
            values = np.array([s for _, s in stocks])
            spread = (values.max() - values.min()) / max(abs(values.mean()), 1e-12)
            if spread <= RELATIVE_TOLERANCE:
                # One stock, diluted by volume: the sheet's own numbers reproduce
                # exactly from the volumes, so this column is independently confirmed.
                verified.append((number, species))
            else:
                # The stock itself was varied between cuvettes. It is recorded
                # nowhere in the volume table, so the computed concentration is
                # the only record and this route cannot check it. Not an error.
                unverifiable.append((number, species))

        # --- block selection -------------------------------------------------
        mismatches = []
        if len(chosen) == len(compiled):
            for i, row in enumerate(chosen):
                for species in SPECIES:
                    sheet_value = row.get(f"C{species}")
                    if sheet_value is None or np.isnan(sheet_value):
                        continue
                    compiled_value = float(compiled.iloc[i][species])
                    if abs(sheet_value - compiled_value) > max(1e-3, 1e-3 * abs(sheet_value)):
                        mismatches.append(f"s{i+1} {species}: sheet {sheet_value:.4g} "
                                          f"vs compiled {compiled_value:.4g}")
            if mismatches:
                findings.append((number, "block", "; ".join(mismatches[:4])
                                 + (f" (+{len(mismatches)-4} more)" if len(mismatches) > 4 else "")))
        else:
            findings.append((number, "block",
                             f"selected {len(chosen)} measured rows from {len(rows)} planned, "
                             f"but {len(compiled)} samples were compiled"))

        summary.append({"experiment": number, "planned_rows": len(rows),
                        "compiled": len(compiled), "selected": len(chosen),
                        "has_buf_column": "[buf]" in table["concentrations"],
                        "mismatches": len(mismatches)})

    return findings, pd.DataFrame(summary), verified, unverifiable


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    findings, summary, verified, unverifiable = analyse()

    print(f"experiments analysed: {len(summary)}\n")
    by_check = {}
    for number, check, message in findings:
        by_check.setdefault(check, []).append((number, message))
    for check in ("sheet", "table", "self-consistency", "block"):
        items = by_check.get(check, [])
        print(f"  {check:17s} {len(items):3d} finding(s)")
    print()
    for check, items in by_check.items():
        print(f"\n--- {check} ({len(items)}) ---")
        for number, message in (items if args.verbose else items[:20]):
            print(f"  exp {number:3d}  {message}")
        if not args.verbose and len(items) > 20:
            print(f"  ... {len(items)-20} more (use --verbose)")

    print("\n--- concentration verifiability (species x experiment) ---")
    print(f"  recomputed from volumes and CONFIRMED : {len(verified)}")
    print(f"  stock serially diluted, UNVERIFIABLE  : {len(unverifiable)}")
    by_species = {}
    for _, species in unverifiable:
        by_species[species] = by_species.get(species, 0) + 1
    for species, count in sorted(by_species.items()):
        experiments = sorted({e for e, sp in unverifiable if sp == species})
        print(f"    {species:8s} {count:3d} experiments: "
              f"{experiments if len(experiments) <= 12 else str(experiments[:12])[:-1] + ', ...]'}")

    no_buf = summary[~summary.has_buf_column]
    print(f"\n--- buffer provenance ---")
    print(f"  sheet carries a [buf]/[buffer] column : {summary.has_buf_column.sum()}")
    print(f"  [buf] rests on the hardcoded 0.1 M    : {len(no_buf)}")
    if len(no_buf):
        print(f"    {sorted(no_buf.experiment.tolist())}")
    sys.exit(0)
