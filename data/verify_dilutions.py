"""
Verifies the serially-diluted concentrations that recompute_concentrations.py
cannot reach.

Two dilution designs exist in this dataset. The volume design (exps 2-62) steps
a component's volume at a fixed stock, so its concentrations re-derive from the
volume table alone -- recompute_concentrations.py handles those. The stock
design (exps 63 onward) holds the volume fixed and serially dilutes the stock
instead, which the volume table does not record.

But the sheets do record it. Every such sheet carries a dilution-series table
("Fortyndingsraekke" / "oplosning") listing, per dilution, the stock volume
taken, the final volume made up to, and the resulting concentration -- and above
it the master stock traced back to a weighed mass and a molar mass. So the whole
chain is checkable:

    mass / molar mass / volume            -> master stock
    master stock * V_stock / V_final      -> each serial dilution
    dilution * V_component / V_cuvette    -> the compiled concentration

Checks performed:

  internal    every row of a dilution table implies the same master stock
  chain       a measured cuvette's concentration equals one of the recorded
              dilutions scaled by its own volume ratio

The chain check runs ONLY on the (experiment, species) pairs that
recompute_concentrations.py could not confirm. Anything the volume route
already verified needs no second opinion, and checking it here produces noise:
in several experiments H2O2 was taken straight from the 30% stock rather than
serially diluted, so it appears in no dilution table while being perfectly
verified by the volume route.

Reports only; changes nothing.

Usage:
    python data/verify_dilutions.py
    python data/verify_dilutions.py --verbose
"""
import argparse
import warnings

import numpy as np
import pandas as pd

from kinetics_io import _text, find_and_parse_experiment_file
from recompute_concentrations import (analyse as analyse_volumes, locate_table,
                                      read_rows, select_measured)

warnings.filterwarnings("ignore")

DATASET_PATH = "data/experiment_data.csv"
MANIFEST_PATH = "data/manifest.csv"
RELATIVE_TOLERANCE = 1e-4


def _number(value):
    return float(value) if isinstance(value, (int, float)) and pd.notna(value) else None


def find_dilution_blocks(sheet):
    """
    Locates every dilution-series table in a sheet.

    A dilution header is a row carrying a stock-volume column, a final-volume
    column and a concentration unit. Both sheet generations use this shape; only
    the unit differs (mol/L in the earlier one, mmol/L in the later).

    Args:
        sheet (pd.DataFrame): A raw, header-less sheet.

    Returns:
        list: One dict per block, each {"row", "entries": [(label, mM)], "unit"}.
    """
    blocks = []
    for row in range(min(120, sheet.shape[0])):
        cells = [_text(sheet.iloc[row, c]) for c in range(min(20, sheet.shape[1]))]
        stock_col = next((c for c, v in enumerate(cells) if "stam" in v), None)
        final_col = next((c for c, v in enumerate(cells) if "final vol" in v), None)
        conc_col = next((c for c, v in enumerate(cells) if v in ("mol/l", "mmol/l")), None)
        if stock_col is None or final_col is None or conc_col is None:
            continue

        scale = 1000.0 if cells[conc_col] == "mol/l" else 1.0
        entries = []
        for r in range(row + 1, min(row + 12, sheet.shape[0])):
            concentration = _number(sheet.iloc[r, conc_col])
            if concentration is None:
                break
            entries.append({
                "label": _text(sheet.iloc[r, 0]),
                "v_stock": _number(sheet.iloc[r, stock_col]),
                "v_final": _number(sheet.iloc[r, final_col]),
                "mM": concentration * scale,
            })
        if len(entries) >= 2:
            blocks.append({"row": row, "entries": entries, "unit": cells[conc_col]})
    return blocks


def find_declared_stocks(sheet):
    """
    Collects stock concentrations declared outside a dilution table.

    Not every stock is made by serial dilution. Several sheets state one
    directly -- a cell reading "[H2O2]" with a "mmol/l" unit and its value in
    the neighbouring rows (exps 127-131 record both peroxide stocks this way).
    Restricted to species-labelled declarations so the candidate pool stays
    small enough for a numeric match to mean something.

    Args:
        sheet (pd.DataFrame): A raw, header-less sheet.

    Returns:
        list: Candidate concentrations in mM.
    """
    declared = []
    rows, cols = min(120, sheet.shape[0]), min(20, sheet.shape[1])
    for r in range(rows):
        for c in range(cols):
            label = _text(sheet.iloc[r, c])
            if label not in ("[h2o2]", "[sub]", "[h202]"):
                continue
            # The unit and value sit within a couple of cells of the label.
            for dr in range(0, 4):
                for dc in range(-1, 3):
                    ur, uc = r + dr, c + dc
                    if not (0 <= ur < rows and 0 <= uc < cols):
                        continue
                    unit = _text(sheet.iloc[ur, uc])
                    if unit not in ("mol/l", "mmol/l"):
                        continue
                    value = _number(sheet.iloc[ur, uc + 1]) if uc + 1 < cols else None
                    if value:
                        declared.append(value * (1000.0 if unit == "mol/l" else 1.0))
    return declared


def check_internal(block):
    """
    Confirms every row of a dilution table implies the same master stock.

    Returns:
        tuple: (ok, implied_master_mM, relative_spread) -- ok is None when the
            table lacks the volumes needed to check.
    """
    implied = []
    for entry in block["entries"]:
        if not entry["v_stock"] or not entry["v_final"]:
            continue
        implied.append(entry["mM"] * entry["v_final"] / entry["v_stock"])
    if len(implied) < 2:
        return None, None, None
    implied = np.array(implied)
    spread = (implied.max() - implied.min()) / max(abs(implied.mean()), 1e-12)
    return spread <= RELATIVE_TOLERANCE, float(implied.mean()), float(spread)


def analyse(dataset_path=DATASET_PATH, manifest_path=MANIFEST_PATH):
    """Runs both checks across the dataset. Returns (findings, summary rows)."""
    data = pd.read_csv(dataset_path)
    manifest = pd.read_csv(manifest_path).set_index("experiment")
    findings, summary = [], []

    # The gap this module exists to fill: what the volume route could not confirm.
    _, _, _, unverifiable = analyse_volumes(dataset_path, manifest_path)
    gap = set(unverifiable)

    for number in sorted(data.experiment.unique()):
        _, sheet = find_and_parse_experiment_file(number, "data/data", "Sheet1")
        if sheet is None:
            continue
        blocks = find_dilution_blocks(sheet)
        if not blocks:
            summary.append({"experiment": number, "blocks": 0,
                            "internal_ok": None, "sub": None, "h2o2": None})
            continue

        internal_ok = True
        for block in blocks:
            ok, master, spread = check_internal(block)
            if ok is False:
                internal_ok = False
                findings.append((number, "internal",
                                 f"dilution table at row {block['row']} implies "
                                 f"inconsistent master stock ({100*spread:.2f}% spread)"))

        # Every concentration the sheet records: serial dilutions plus any
        # stock declared directly.
        available = ([e["mM"] for b in blocks for e in b["entries"]]
                     + find_declared_stocks(sheet))

        header_row, table = locate_table(sheet)
        if header_row is None:
            continue
        rows = read_rows(sheet, header_row, table)
        compiled = data[data.experiment == number].sort_values("sample")
        has_enzyme = (bool(manifest.loc[number]["has_enzyme"])
                      if number in manifest.index else True)
        chosen = select_measured(rows, len(compiled), has_enzyme)

        matched = {}
        for species in ("[sub]", "[h2o2]"):
            if len(chosen) != len(compiled) or (number, species) not in gap:
                continue
            hits = 0
            for i, row in enumerate(chosen):
                volume, total = row.get(f"V{species}"), row["total"]
                if volume is None or np.isnan(volume) or volume == 0:
                    continue
                target = float(compiled.iloc[i][species])
                if any(abs(candidate * volume / total - target)
                       <= max(1e-3, 1e-3 * abs(target)) for candidate in available):
                    hits += 1
                else:
                    findings.append((number, "chain",
                                     f"s{i+1} {species}={target:.4g} matches no recorded "
                                     f"dilution scaled by {volume}/{total}"))
            matched[species] = hits

        summary.append({"experiment": number, "blocks": len(blocks),
                        "internal_ok": internal_ok,
                        "sub": matched.get("[sub]"), "h2o2": matched.get("[h2o2]")})

    return findings, pd.DataFrame(summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    findings, summary = analyse()
    with_blocks = summary[summary.blocks > 0]

    print(f"experiments with a dilution series : {len(with_blocks)} / {len(summary)}")
    print(f"dilution tables internally consistent : "
          f"{(with_blocks.internal_ok == True).sum()} / {len(with_blocks)}")
    print(f"\ncuvette concentrations traced back to a recorded dilution")
    print(f"(restricted to what the volume route could NOT confirm):")
    for species in ("sub", "h2o2"):
        traced = with_blocks[species].fillna(0).sum()
        print(f"  {species:5s} {int(traced):4d} cuvettes across "
              f"{int((with_blocks[species].fillna(0) > 0).sum())} experiments")

    by_check = {}
    for number, check, message in findings:
        by_check.setdefault(check, []).append((number, message))
    print(f"\nfindings: {len(findings)}")
    for check, items in by_check.items():
        print(f"\n--- {check} ({len(items)}) ---")
        for number, message in (items if args.verbose else items[:15]):
            print(f"  exp {number:3d}  {message}")
        if not args.verbose and len(items) > 15:
            print(f"  ... {len(items)-15} more (use --verbose)")
