"""
Verifies every experiment's [enz] against the WEIGHING it was made from.

WHY THIS EXISTS. The 4OMe-BnOH temperature series (exps 14-19) carries two
different enzyme concentrations -- 0.272695 mM on exps 14, 15, 17, 18 and
0.240683 mM on exps 16 and 19 -- and exp 16 sits between two runs that use the
higher value. An 11.7% step at one interior point of a six-point Arrhenius
series is either a real change of stock or a transcription error, and the two
call for opposite treatments: a real change is divided out, an error is
corrected. Nothing in the pipeline could tell them apart, because `[enz]`
reaches `experiment_data.csv` from a single cell.

WHAT MAKES THIS AN INDEPENDENT CHECK. The compiled `[enz]` is the sheet's `kuv`
cell. Beside it every sheet also records the preparation that produced it:

    mass / molar mass / volume          -> stock concentration, mmol/l
    stock * Enz[ml] / Vol[ml]           -> the cuvette concentration, `kuv`

Four numbers written at four different moments, none of them the one the
pipeline reads. Recomputing `kuv` from the mass and the cuvette volumes
therefore confirms or refutes it without consulting it -- and a transcription
error into `kuv` would leave the chain inconsistent, while a real restock
changes the mass and keeps the chain intact.

WHAT IT FOUND (2026-09-01). Four stock preparations across the campaign, each
arithmetically self-consistent to six figures and each used by a contiguous run
of experiments:

    0.0122 g / 3.3 ml  ->  3.5066  mM  ->  kuv 0.17533   exps 2-9
    0.023  g / 4.0 ml  ->  5.45391 mM  ->  kuv 0.272695  exps 10-15, 17, 18
    0.0203 g / 4.0 ml  ->  4.81367 mM  ->  kuv 0.240683  exps 16, 19-22, 32-36
    0.0228 g / 4.0 ml  ->  5.40648 mM  ->  kuv 0.270324  exps 37, 41-49

Exp 16 is the only experiment whose stock is out of sequence with its
neighbours, which is what `--sequence` reports. It is NOT evidence of an error:
see temperature_series/ANALYSIS.md, where the kinetics decide it independently.

Reports only; changes nothing.

Usage:
    python data/verify_enzyme_stock.py
    python data/verify_enzyme_stock.py --sequence
"""
import argparse
import glob
import os
import re
import warnings

import numpy as np
import pandas as pd

SHEET_DIR = "data/data"
COMPILED = "data/experiment_data.csv"
# The label that marks the enzyme's own preparation block. The sheets name the
# compound rather than "enzyme", and the name is stable across the campaign.
ENZYME_LABEL = "a-diesterketon"
# The cuvette table's column headers, used to find the volume ratio.
VOLUME_HEADERS = ("Enz [ml]", "Vol [ml]")
# A recomputed kuv must match the sheet's own to this relative tolerance. The
# sheet stores a spreadsheet formula's output, so the two agree to machine
# precision when the chain is intact; 1e-4 is loose enough for a rounded cell
# and far tighter than the 11.7% step this module was written to adjudicate.
CHAIN_TOLERANCE = 1e-4
# `experiment_data.csv` stores [enz] to three decimals, so 0.17533 is compiled
# as 0.175. Comparing at full precision reports 35 of 83 experiments as
# disagreeing when every one of them is the same number rounded. The
# extraction is checked at the precision it actually stores.
COMPILED_DECIMALS = 3


def sheet_paths(directory=SHEET_DIR):
    """Every compiled experiment's workbook, keyed by experiment number."""
    found = {}
    for path in sorted(glob.glob(os.path.join(directory, "mads_t*.xls"))):
        match = re.search(r"mads_t(\d+)", os.path.basename(path))
        if match:
            found.setdefault(int(match.group(1)), path)
    return found


def read_preparation(path):
    """
    The enzyme preparation recorded in one workbook.

    Returns a dict with molar_mass, grams, litres, stock_mM, kuv_sheet and the
    cuvette volumes, or None where the sheet has no enzyme block -- which is
    what a genuine background run looks like and is not a fault.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            sheet = pd.read_excel(path, header=None, sheet_name="Sheet1")
        except (ValueError, KeyError):
            return None

    column = None
    for index in range(sheet.shape[1]):
        if (sheet[index].astype(str).str.strip() == ENZYME_LABEL).any():
            column = index
    if column is None or column < 1:
        return None
    labels = sheet[column - 1].astype(str).str.strip()

    def value(label):
        rows = sheet.index[labels == label]
        if not len(rows):
            return np.nan
        try:
            return float(sheet.iloc[rows[0], column])
        except (TypeError, ValueError):
            return np.nan

    # The cuvette volumes come from the measured rows of the volume table, not
    # from the preparation block: the ratio that turns a stock into a cuvette
    # concentration is a property of how the cuvette was filled.
    header = None
    for index in range(sheet.shape[0]):
        row = sheet.iloc[index].astype(str).str.strip()
        if all((row == name).any() for name in VOLUME_HEADERS):
            header = index
            break
    enzyme_ml = total_ml = np.nan
    if header is not None:
        names = sheet.iloc[header].astype(str).str.strip()
        enzyme_column = int(np.argmax((names == "Enz [ml]").to_numpy()))
        total_column = int(np.argmax((names == "Vol [ml]").to_numpy()))
        body = sheet.iloc[header + 1:header + 5]
        enzyme_ml = pd.to_numeric(body[enzyme_column], errors="coerce").max()
        total_ml = pd.to_numeric(body[total_column], errors="coerce").median()

    return {
        "molar_mass": value("g/mol"), "grams": value("g"),
        "litres": value("l"), "stock_mM": value("mmol/l"),
        "kuv_sheet": value("kuv"),
        "enzyme_ml": float(enzyme_ml) if np.isfinite(enzyme_ml) else np.nan,
        "total_ml": float(total_ml) if np.isfinite(total_ml) else np.nan,
    }


def audit(directory=SHEET_DIR, compiled=COMPILED):
    """
    Recompute every experiment's [enz] from its weighing and compare.

    Returns a DataFrame indexed by experiment with the recorded value, the
    value the mass implies, their ratio, and whether the sheet's own chain is
    internally consistent.
    """
    census = pd.read_csv(compiled)
    recorded = census.groupby("experiment")["[enz]"].max()
    rows = []
    for experiment, path in sorted(sheet_paths(directory).items()):
        if experiment not in recorded.index:
            continue
        preparation = read_preparation(path)
        if preparation is None:
            continue
        molar, grams = preparation["molar_mass"], preparation["grams"]
        litres = preparation["litres"]
        stock = (grams / molar / litres * 1000.0
                 if molar and litres and np.isfinite(grams) else np.nan)
        ratio = (preparation["enzyme_ml"] / preparation["total_ml"]
                 if preparation["total_ml"] else np.nan)
        implied = stock * ratio
        sheet_kuv = preparation["kuv_sheet"]
        # Two separate questions, kept apart on purpose. `chain_ok` asks
        # whether the SHEET is self-consistent -- mass, volumes and kuv telling
        # one story. `agrees` asks whether the pipeline read that kuv
        # correctly. A transcription error into kuv breaks the first; a broken
        # extraction breaks only the second.
        chain_ok = bool(np.isfinite(implied) and np.isfinite(sheet_kuv)
                        and abs(implied - sheet_kuv)
                        <= CHAIN_TOLERANCE * max(abs(sheet_kuv), 1e-12))
        rows.append({
            "experiment": experiment,
            "grams": grams, "litres": litres,
            "stock_mM": stock, "enzyme_ml": preparation["enzyme_ml"],
            "total_ml": preparation["total_ml"],
            "kuv_sheet": sheet_kuv, "implied_mM": implied,
            "recorded_mM": float(recorded[experiment]),
            "chain_ok": chain_ok,
        })
    table = pd.DataFrame(rows).set_index("experiment")
    table["agrees"] = (table.recorded_mM.round(COMPILED_DECIMALS)
                       == table.kuv_sheet.round(COMPILED_DECIMALS))
    # An enzyme-free run's sheet still carries a preparation block -- the
    # workbook was copied from a catalysed run -- but its cuvettes took 0 ml of
    # it. That block describes a stock the run did not use, so it must not
    # enter the stock census or the sequence check, and it is not a fault.
    table["used"] = table.recorded_mM > 0
    return table


def stocks(table=None):
    """
    The distinct preparations, and which experiments used each.

    A stock is identified by its weighed mass and made-up volume -- the two
    numbers a restock changes -- rather than by the resulting concentration, so
    two preparations that happen to land on the same concentration stay
    distinct.
    """
    if table is None:
        table = audit()
    live = table[(table.grams > 0) & table.used]
    grouped = live.groupby(["grams", "litres"])
    rows = []
    for (grams, litres), group in grouped:
        experiments = sorted(int(e) for e in group.index)
        rows.append({
            "grams": grams, "litres": litres,
            "stock_mM": float(group.stock_mM.iloc[0]),
            "kuv_mM": float(group.kuv_sheet.median()),
            "n": len(experiments), "first": experiments[0],
            "last": experiments[-1], "experiments": experiments,
        })
    return pd.DataFrame(rows).sort_values("first").reset_index(drop=True)


def out_of_sequence(table=None):
    """
    Experiments whose stock is not the one their neighbours used.

    Stocks are prepared once and used until exhausted, so in experiment order
    the mass should change rarely and never change back. An experiment whose
    mass differs from BOTH the experiment before it and the one after it
    interrupts a run, and is the shape a copied workbook would leave.

    IT IS A FLAG, NOT A FAULT. Being out of sequence says only that the
    documentary evidence is ambiguous; deciding it needs a source outside the
    workbooks. For exp 16 that source is the Arrhenius fit -- see
    temperature_series/ANALYSIS.md.
    """
    if table is None:
        table = audit()
    live = table[(table.grams > 0) & table.used].sort_index()
    masses = live.grams.to_numpy()
    experiments = [int(e) for e in live.index]
    flagged = []
    for index in range(1, len(masses) - 1):
        if masses[index] != masses[index - 1] and masses[index] != masses[index + 1]:
            flagged.append({
                "experiment": experiments[index],
                "grams": masses[index],
                "before": f"exp {experiments[index - 1]} at {masses[index - 1]:g} g",
                "after": f"exp {experiments[index + 1]} at {masses[index + 1]:g} g",
            })
    return pd.DataFrame(flagged)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--sequence", action="store_true",
                        help="also list experiments whose stock interrupts a run")
    arguments = parser.parse_args()

    table = audit()
    # Only runs that actually used the enzyme: an enzyme-free run's inherited
    # block is not a chain anyone relied on.
    table = table[table.used] if len(table) else table
    broken = table[~table.chain_ok]
    disagree = table[~table.agrees]
    print(f"{len(table)} experiments carry an enzyme preparation\n")
    print("the distinct stocks, by weighed mass:")
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(stocks(table).to_string(index=False,
                                      float_format=lambda v: f"{v:g}"))

    print(f"\n  sheet chain self-consistent (mass -> stock -> kuv): "
          f"{int(table.chain_ok.sum())} of {len(table)}")
    if len(broken):
        print("  BROKEN CHAINS -- a kuv its own preparation does not produce:")
        print(broken[["grams", "litres", "implied_mM", "kuv_sheet"]]
              .to_string(float_format=lambda v: f"{v:.6g}"))
    print(f"  compiled [enz] equals the sheet's kuv: "
          f"{int(table.agrees.sum())} of {len(table)}")
    if len(disagree):
        print("  EXTRACTION DISAGREES WITH THE SHEET:")
        print(disagree[["kuv_sheet", "recorded_mM"]]
              .to_string(float_format=lambda v: f"{v:.6g}"))
    # A broken chain on an already-excluded experiment is recorded, not failed
    # on: exp 58 runs backwards and no result depends on its [enz]. A broken
    # chain anywhere else is a real fault and should stop a gate.
    from build_manifest import KNOWN_EXCLUSIONS
    live_faults = ([e for e in broken.index if e not in KNOWN_EXCLUSIONS]
                   + [e for e in disagree.index if e not in KNOWN_EXCLUSIONS])
    if len(broken) or len(disagree):
        excluded = sorted(set(list(broken.index) + list(disagree.index))
                          & set(KNOWN_EXCLUSIONS))
        if excluded:
            print(f"  ({', '.join(f'exp {e}' for e in excluded)}: already "
                  f"excluded, so not treated as a failure)")

    if arguments.sequence:
        flagged = out_of_sequence(table)
        print(f"\nstock out of sequence with its neighbours: {len(flagged)}")
        if len(flagged):
            print(flagged.to_string(index=False))
        print("  A flag, not a fault: it says the workbooks alone cannot decide,\n"
              "  not that the value is wrong. See temperature_series/ANALYSIS.md.")

    return 1 if live_faults else 0


if __name__ == "__main__":
    raise SystemExit(main())
