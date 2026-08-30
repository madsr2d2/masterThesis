"""
Bootstraps data/manifest.csv -- the declared ground truth for the dataset --
from the independent sources that exist outside the .xls sheets:

  0. The sheet header. Each sheet declares the monitoring wavelength and the
                      extinction coefficient (`abs [nm]`, `e [U/mM]`).
                      kinetics_io ignores both, hardcoding e per substrate in
                      SUBSTRATE_PROPERTIES, so recording them here makes any
                      divergence visible.
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
import contextlib
import io
import os
import re

import pandas as pd

from kinetics_io import find_header_row, find_and_parse_experiment_file, load_experiment

MADS_DIR = "data/Mads"

# Experiments removed by clean_experiment_dataframe, with the reasons recorded
# in DATA_VERIFICATION.md. Seeded here so status/exclusion becomes declared data
# rather than a literal buried in a notebook cell.
KNOWN_EXCLUSIONS = {
    64: ("aborted run: 7 minutes at dt = 28 s, and three of the four curves are "
         "flat or backwards (net +0.000, -0.001, -0.007 against sample 1's "
         "+0.006). Same criterion as exps 72 and 82. Its session was troubled "
         "throughout -- the sibling sheet mads_t063_..._no_E_94 is named "
         "NO_DATA_FILE and t063 has no instrument export at all. Compiled "
         "anyway, on 2026-08-30, so the archive-to-dataset mapping is complete "
         "and the exclusion is recorded rather than silent"),
    57: ("[sub] is not recoverable: the workbook was copied from t056, a "
         "4-brom-BnOH run, and the whole substrate stock block came with it. "
         "The compound was 4OMe-BnOH (ruled 2026-08-30), so a fresh stock must "
         "have been weighed, and the sheet records t056's instead -- 0.3511 g "
         "in 0.1 L, matching t056 to four digits, which no independent weighing "
         "would. No earlier 4OMe preparation matches those numbers either. "
         "Rescaling by 187.03/138.17 assumes only the molar mass was stale, "
         "which is the one assumption the evidence argues against, so the "
         "concentration is unknown rather than merely wrong"),
    50: ("reaction-direction: all four curves descend 0.075-0.105 AU, and the "
         "descent has no ordering by substrate at all -- rho(|net| vs [sub]) = "
         "0.00 against +1.00 for every healthy run. Nothing that depends on "
         "substrate produced that signal, and negating the curves does not "
         "help, since a sign flip preserves the ordering. The repeat under "
         "identical conditions (exp 55: same buffer, pH, ladder and [enz]) "
         "rises +0.145 to +0.476 monotonically, and the same-day sibling exp "
         "51 is clean -- the criterion that excluded exps 72 and 82. Ruled "
         "2026-08-30. It survived earlier passes only because it is "
         "hand-sorted into data/Mads/'good data BnOH'"),
    58: "reaction-direction: whole experiment runs backwards",
    72: "flat progress curves in every sample; same-day sibling 71 normal",
    77: "reaction-direction: whole experiment runs backwards",
    78: "reaction-direction: whole experiment runs backwards",
    79: "reaction-direction: whole experiment runs backwards",
    82: "flat progress curves in every sample; same-day sibling 83 normal",
    85: ("reaction-direction: all seven curves crash 0.30-0.63 AU, most of it "
         "in the first 15 of 60 minutes, which is a decay rather than a "
         "reaction. The drop is ANTI-correlated with substrate -- 14.30 mM "
         "falls -0.298 where 0.51 mM falls -0.604 -- so the observable is not "
         "product formation. Substrate inhibition was considered and does not "
         "explain it: inhibition drives the rate toward zero, never below, and "
         "it predicts the low-[sub] cuvette to be the largest POSITIVE. The "
         "sheet gives pH 11.84 in Na2CO3, 1.5 units above carbonate's pKa2, so "
         "the run is barely buffered; exps 86-109 at pH 11 were hand-sorted "
         "into 'bad data pH ca. 11' and exp 85 escaped only because it sits in "
         "'carbonate buffer'. Ruled 2026-08-30. It carries the widest substrate "
         "range in the archive, 0.51-14.30 mM, and may be worth citing as "
         "evidence of catalyst instability at pH 11.8 -- but not as kinetics"),
    84: "hand-sorted into data/Mads/'bad data'",
}

# Deviations that are understood and accepted, so the deep checks warn rather
# than error on them. Keyed by experiment -> (check names, reason).
_PLANNED_ENZYME_ROWS = (
    "enzuse",
    "the sheet's enzyme block describes the four PLANNED with-enzyme cuvettes, "
    "which were never run. These are enzyme-free buffer titrations -- the "
    "filenames say with_NO_E and all five sit in data/Mads/'No enzyme'/ -- and "
    "only cuvettes 5-8, the no-enzyme half of the eight-row plan, were measured. "
    "[enz] = 0 is correct and the block's 0.24-0.27 mM is the plan. See "
    "DATA_VERIFICATION.md 2026-08-30."
)

_BROKEN_ENZ_COLUMN = (
    "block",
    "[enz] is deliberately taken from the sheet's header block rather than from "
    "its cuvette table. The table's '[Enz] mmol/l' column holds 0.000001 in every "
    "measured row -- a broken formula roughly 14,000x too low -- while the header "
    "block's 'kuv' row gives 0.01399 (exp 79) / 0.027981 (exp 80), which agrees "
    "with the recorded volumes and the enzyme stock. The deep check reads the "
    "table, so it necessarily disagrees here. Ruled 2026-08-30, see "
    "DATA_VERIFICATION.md."
)

KNOWN_ACCEPTED_DEVIATIONS = {
    32: _PLANNED_ENZYME_ROWS,
    34: _PLANNED_ENZYME_ROWS,
    35: _PLANNED_ENZYME_ROWS,
    36: _PLANNED_ENZYME_ROWS,
    37: _PLANNED_ENZYME_ROWS,
    58: ("bufvolume;bufcompiled",
         "the cuvettes here are 2.10-2.11 ml, not the 2 ml the volume fallback "
         "assumes, so every [buf] is about 5% high (90.0 where the sheet's own "
         "volumes give 85.71). Found 2026-08-30 by verify_buffer.py. The "
         "experiment is excluded for running backwards, so no result depends on "
         "it; the deviation is recorded rather than corrected because correcting "
         "an excluded run would change the compiled CSV for no benefit. It is "
         "the only experiment in the archive where the hardcoded 2 ml is wrong"),
    79: _BROKEN_ENZ_COLUMN,
    80: _BROKEN_ENZ_COLUMN,
    128: ("block;chain",
          'sample 5 is the reference row "ref 5" (a matched no-enzyme blank) with its own '
          "cuvette volume, not a 5th titration condition, so both its concentrations and "
          "its volume ratio differ from the measured cuvettes. Documented in "
          "DATA_VERIFICATION.md Round 4; the sample is already handled downstream."),
}

# Questions raised outside the automatic filename/extraction comparison, seeded
# so they survive a rebuild. Listing a field here also tells validate_dataset.py
# to warn rather than error on it, since it is known and awaiting a ruling.
KNOWN_OPEN_QUESTIONS = {}
# Adjudicated questions: a value the sheet states that we have decided against,
# with the evidence. A ruling overrides the extracted value, marks the field's
# provenance as "ruling" so it is never mistaken for something read off a file,
# and records what the sheet said in notes so the evidence is not erased.
# Rulings are settled; open questions are not. Both survive a rebuild.
_ABS_285_RULING = (
    300.0,
    "sheet declares 285 nm; ruled 300 nm on 2026-08-30. The early series "
    "alternates BnOH and 4OMe runs (t001 BnOH 285/1.23, t002 4OMe 285/7.53, "
    "t003 BnOH 285/1.23, ...), so the 4OMe workbooks were copied from the BnOH "
    "one with the substrate and e changed and the abs [nm] cell left at 285. "
    "The cell is corrected to 300 at t011 and stays 300 for every later 4OMe "
    "run, while e stays 7.53 across that boundary -- had the instrument really "
    "been retuned, e would have had to change with it, which is exactly what "
    "exps 57/58 do (285 nm paired with e = 1.59). No concentration is affected "
    "either way here, since sheet and dataset agree on e = 7.53. Not directly "
    "confirmable: the .txt instrument export records batch, date, data mode and "
    "smoothing but not the wavelength, and conversion stays under 1%, so no "
    "internal consistency test has any power."
)
RULINGS = {number: {"abs_nm": _ABS_285_RULING} for number in (2, 4, 5, 7, 8, 9, 10)}

# Two experiments disagree with their filename on pH, both by hundredths, and
# both were ruled to the sheet on 2026-08-30: the sheet carries the reading taken
# on the day and the filename is typed from it. Neither changes a number -- the
# dataset already held the sheet value in both cases -- but the provenance is now
# stated rather than left as an unadjudicated conflict.
_PH_FROM_SHEET = (
    "filename says {filename}; ruled to the sheet's {sheet} on 2026-08-30. The "
    "sheet carries the reading taken on the day and the filename is typed from "
    "it. The dataset already held the sheet value, so no number changes."
)
for _exp, _sheet_pH, _filename_pH in ((9, 5.67, "5.64"), (38, 7.00, "6.97")):
    RULINGS.setdefault(_exp, {})["pH"] = (
        _sheet_pH, _PH_FROM_SHEET.format(filename=_filename_pH, sheet=_sheet_pH))

# Exps 79 and 80 are enzyme runs. Their "[Enz] mmol/l" column holds 0.000001 in
# every measured row -- a broken formula about 14,000x too low -- which the
# extraction read faithfully and rounding turned into zero, so the manifest
# inherited has_enzyme = False against filenames that say "with_E". The right
# value is in the header block's "kuv" row and agrees with the volumes (stock
# 0.559618 mM; 0.05/2 ml -> 0.01399 for exp 79, 0.1/2 ml -> 0.027981 for exp 80).
# Corrected in kinetics_io by EXPERIMENT_CORRECTIONS. Ruled 2026-08-30.
_ENZYME_RUN = (
    True,
    "extraction said no enzyme; ruled an enzyme run on 2026-08-30. The cuvette "
    "table carries only 'Enz [ml]' and no '[Enz] mmol/l' column, so the "
    "extraction fell through to its zero default. The sheet's header block "
    "declares the cuvette concentration on the 'kuv' row and it agrees with the "
    "volumes. Corrected in kinetics_io.EXPERIMENT_CORRECTIONS."
)
for _exp in (79, 80):
    RULINGS.setdefault(_exp, {})["has_enzyme"] = _ENZYME_RUN

# Exps 84 and 85 are 4OMe-BnOH runs whose FILENAMES say BnOH. Ruled 2026-08-30
# on four independent grounds inside the sheets, against the filename alone:
#   the stock-solution block is labelled "Stamopløsning 4 / 4-MeO-BnOH [g]";
#   the molar mass it computes from is 138.17 (4-methoxybenzyl alcohol;
#     benzyl alcohol is 108.14), so every concentration already assumes 4OMe;
#   the method string reads "1h_kcat(4OMeBnOH)_7cuv";
#   the optics are 300 nm and e = 7.53, the 4OMe convention.
# Both sheets also carry "83" as their experiment number, i.e. the workbook was
# copied from t083 (a genuine BnOH run) and the header never updated -- which is
# almost certainly where the filename's "BnOH" came from too. The dataset
# already carries 4OMe-BnOH, so this ruling changes no number.
_SUBSTRATE_84_85 = (
    "4OMe-BnOH",
    "filename says BnOH; ruled 4OMe-BnOH on 2026-08-30. The sheet's stock block "
    "is labelled '4-MeO-BnOH [g]' and computes from M = 138.17 g/mol "
    "(4-methoxybenzyl alcohol; benzyl alcohol is 108.14), so every concentration "
    "already assumes 4OMe; the method string reads '1h_kcat(4OMeBnOH)_7cuv'; and "
    "the optics are the 4OMe convention. The sheet's own experiment number is "
    "83, so the workbook was copied from the BnOH run t083 and the header never "
    "updated -- the likely source of the filename error. No number changes."
)
for _exp in (84, 85):
    RULINGS.setdefault(_exp, {})["substrate"] = _SUBSTRATE_84_85

# Exps 57 and 58 are 4-methoxybenzyl alcohol runs (ruled 2026-08-30) whose
# workbook was copied wholesale from t056, a 4-bromobenzyl alcohol run: the
# substrate label in cell C4 and the filename were updated, the entire stock
# block and the optics were not. So the sheet's 285 nm, its e = 1.59 and its
# M = 187.03 are all t056's, and the dataset's 300 nm / e = 7.53 are correct by
# the substrate convention. The molar mass is the consequential one -- the stock
# of 0.3511 g in 0.1 L was divided by 187.03 instead of 138.17, so every [sub]
# in these two experiments was low by that ratio. Corrected in kinetics_io by
# CONCENTRATION_RESCALINGS. See DATA_VERIFICATION.md 2026-08-30.
_COPIED_FROM_T056 = (
    "sheet declares {declared}; ruled to the 4OMe-BnOH convention on 2026-08-30. "
    "The workbook was copied from t056, a 4-brom-BnOH run, with the substrate "
    "label and the filename changed and the stock block and optics left behind: "
    "the stock (0.3511 g / 0.1 L) and the molar mass (187.03) are byte-identical "
    "to t056's, and e = 1.59 is the 4-bromo convention. The stale M = 187.03 "
    "made every [sub] low by 187.03/138.17 = 1.3536x, corrected in kinetics_io "
    "by CONCENTRATION_RESCALINGS."
)
for _exp in (57, 58):
    RULINGS.setdefault(_exp, {}).update({
        "abs_nm": (300.0, _COPIED_FROM_T056.format(
            declared="285 nm, paired with e = 1.59")),
        "e_declared": (7.53, _COPIED_FROM_T056.format(declared="e = 1.59")),
    })

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


def read_sheet_claims(sheet):
    """
    Returns the fields the SHEET itself supplies, independently of any filename.

    kinetics_io's find_* helpers are already sheet-first, taking the filename
    only as a fallback, so calling them with filename=None isolates what the
    sheet actually declares. That distinction is what the provenance column
    needs: stamping "filename" whenever a filename agreed made the sheet's
    corroboration invisible, and understated the evidence for 83 of the 100 pH
    values, which the sheets carry beside a 'pH' label.

    Args:
        sheet (pd.DataFrame): A raw, header-less sheet.

    Returns:
        dict: Any of "pH", "T", "buffer", "substrate" the sheet declares.
    """
    from kinetics_io import (find_buffer_type, find_pH_value_in_range,
                             find_substrate_type, find_temperature_value_in_range)

    readers = {"pH": find_pH_value_in_range, "T": find_temperature_value_in_range,
               "buffer": find_buffer_type, "substrate": find_substrate_type}
    claims = {}
    with contextlib.redirect_stdout(io.StringIO()):
        for field, reader in readers.items():
            try:
                value = reader(sheet, (0, 100), (0, 100), filename=None)
            except Exception:
                value = None
            if value is not None and not (isinstance(value, float) and pd.isna(value)):
                claims[field] = value
    return claims


def read_sheet_optics(sheet):
    """
    Reads the monitoring wavelength and extinction coefficient a sheet declares.

    Args:
        sheet (pd.DataFrame): A raw, header-less sheet.

    Returns:
        tuple: (wavelength_nm, e_per_mM), either possibly None.
    """
    wavelength = extinction = None
    for row in range(min(40, sheet.shape[0])):
        for col in range(min(16, sheet.shape[1])):
            cell = ("" if pd.isna(sheet.iloc[row, col])
                    else str(sheet.iloc[row, col]).strip().lower())
            if cell == "nm" or cell.startswith("abs ["):
                for offset in (-1, 1):
                    if 0 <= col + offset < sheet.shape[1]:
                        value = sheet.iloc[row, col + offset]
                        if isinstance(value, (int, float)) and pd.notna(value) and 200 < value < 800:
                            wavelength = float(value)
            if cell == "e" or cell.startswith("e ["):
                for offset in (1, 2):
                    if col + offset < sheet.shape[1]:
                        value = sheet.iloc[row, col + offset]
                        if isinstance(value, (int, float)) and pd.notna(value) and 0 < value < 1000:
                            extinction = float(value)
    return wavelength, extinction


# A buffer-preparation sentence written into a sheet as free text. Exactly one
# exists in the whole archive -- exp 13's -- and it is the only direct evidence
# for the 0.1 M buffer stock that 32 further experiments rely on by convention:
#
#   "The Buffer was prepeared by mixing 0.0699g NaBH4 with 0.1964 g B(OH)3 in
#    0.05 l water and the pH was adjusted with a NaOH solution."
#
#   B(OH)3  0.1964 g / 61.83 g/mol = 3.176 mmol
#   NaBH4   0.0699 g / 37.83 g/mol = 1.848 mmol
#                            total = 5.024 mmol boron in 0.05 l = 100.5 mM
#
# i.e. 0.1 M to within half a percent. See DATA_VERIFICATION.md 2026-08-30 for
# the caveats: NaBH4 is almost certainly a mis-transcription (it is a reducing
# agent, and would not survive contact with H2O2), though every plausible
# substitute lands at 64-85 mM and only the literal reading gives 0.1 M.
RECIPE_PATTERN = re.compile(r"buffer\s+was\s+prep", re.I)

# A stock molarity written as a text label beside a cuvette, e.g. "1 (0.1M)".
LABEL_MOLARITY = re.compile(r"\(\s*(\d+[.,]?\d*)\s*M\s*\)", re.I)
# A stock molarity in a filename, e.g. "phosphate_0.1M".
FILENAME_MOLARITY = re.compile(r"(\d+[.,]?\d*)\s*M(?![a-z])")


def find_recipe(sheet):
    """
    Returns the buffer-preparation sentence a sheet carries, or None.

    Args:
        sheet (pd.DataFrame): A raw, header-less sheet.

    Returns:
        str or None: The sentence, whitespace-collapsed.
    """
    for value in sheet.values.ravel():
        if isinstance(value, str) and RECIPE_PATTERN.search(value):
            return " ".join(value.split())
    return None


def classify_buffer(sheet, filename):
    """
    Classifies how each experiment's [buf] is known, and how it was titrated.

    Neither is stated as such in any file, and both are needed before fitting:
    the design says whether [buf] is expected to vary across the cuvettes at
    all, and the provenance says how far the recorded value can be trusted.

    design
        'stock'   fixed volumes throughout; concentration varied by using
                  different stock solutions, so [buf] is constant
        'volume'  buffer volume traded against substrate volume, so [buf]
                  genuinely varies across the cuvettes
    buf_provenance, in descending order of authority
        'declared'     the sheet states [buf] per cuvette
        'sheet-recipe' the sheet carries a weighed buffer recipe in free text
        'filename'     the stock molarity is in the filename
        'sheet-text'   it appears only as a label beside a cuvette
        'assumed'      stated nowhere; (V_buf/2)*100 assumes a 0.1 M stock

    Reference channels are excluded from the design test: several sheets run
    the paired no-enzyme cuvettes at different volumes than the samples (exp
    146 uses 0.5 ml buffer in 1 ml for samples and 1.0 ml in 2 ml for
    references), which makes a fixed-volume experiment look like a titration.

    Args:
        sheet (pd.DataFrame): A raw, header-less sheet.
        filename (str): The experiment's .xls filename.

    Returns:
        tuple: (design, buf_provenance, recipe or None)
    """
    header_row = None
    with contextlib.redirect_stdout(io.StringIO()):
        header_row = find_header_row(sheet)

    labels, measured = [], []
    if header_row is not None:
        for column in range(sheet.shape[1]):
            parts = [str(sheet.iat[row, column])
                     for row in range(header_row, min(header_row + 2, len(sheet)))
                     if str(sheet.iat[row, column]) != "nan"]
            labels.append(" ".join(parts).strip())
        volume_column = next((i for i, l in enumerate(labels)
                              if l.lower().startswith("buf") and "ml" in l.lower()), None)
        total_column = next((i for i, l in enumerate(labels)
                             if l.lower().startswith("vol")), None)
        first_column = next((i for i, l in enumerate(labels) if l), None)
        if volume_column is not None:
            for row in range(header_row + 1, min(header_row + 26, len(sheet))):
                label = str(sheet.iat[row, first_column]).strip().lower()
                if label.startswith(("ref", "sum")) or label == "nan":
                    continue
                if total_column is not None and pd.isna(sheet.iat[row, total_column]):
                    continue
                value = pd.to_numeric(pd.Series([sheet.iat[row, volume_column]]),
                                      errors="coerce").iloc[0]
                if pd.notna(value):
                    measured.append(float(value))

    declared = any(l.lower().startswith("[buf") for l in labels)
    design = ("stock" if declared and not measured
              else "volume" if len(set(measured)) > 1
              else "stock" if measured else "unknown")

    recipe = find_recipe(sheet)
    if declared:
        provenance = "declared"
    elif recipe:
        provenance = "sheet-recipe"
    elif FILENAME_MOLARITY.search(filename):
        provenance = "filename"
    elif any(isinstance(v, str) and LABEL_MOLARITY.search(v)
             for v in sheet.values.ravel()):
        provenance = "sheet-text"
    else:
        provenance = "assumed"
    return design, provenance, recipe


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

        _, sheet = find_and_parse_experiment_file(number, directory, "Sheet1")
        wavelength, extinction = read_sheet_optics(sheet) if sheet is not None else (None, None)

        from_name = parse_filename(str(experiment["xls_file"]))
        from_folder = folder_claims.get(number, {})
        extracted = {
            "pH": experiment["pH"],
            "T": experiment["T"],
            "substrate": experiment["substrate"],
            "buffer": experiment["buffer"],
            "has_enzyme": any(s["[enz]"] for s in experiment["samples"]),
        }

        rulings = RULINGS.get(number, {})
        from_sheet = read_sheet_claims(sheet) if sheet is not None else {}
        resolved, provenance, open_questions = {}, {}, []
        for field in ("pH", "T", "substrate", "buffer", "has_enzyme"):
            if field in rulings:
                resolved[field] = rulings[field][0]
                provenance[field] = "ruling"
                continue
            declared = from_folder.get(field, from_name.get(field))
            declared_by = ("folder" if field in from_folder
                           else "filename" if field in from_name else None)
            if declared is None:
                resolved[field] = extracted[field]
                provenance[field] = "sheet" if field in from_sheet else "extracted"
                continue

            if _agrees(declared, extracted[field]):
                resolved[field] = declared
                # Two independent sources agreeing is the strongest state the
                # metadata reaches, and it deserves its own label: stamping
                # "filename" here hid that the sheet says the same thing.
                provenance[field] = (f"sheet+{declared_by}"
                                     if field in from_sheet and
                                     _agrees(from_sheet[field], declared)
                                     else declared_by)
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

        # What the sheet itself says, kept unmodified so a ruling never erases
        # the evidence it was made against. abs_nm/e_declared below carry the
        # adjudicated values; these two carry the observation.
        sheet_wavelength, sheet_extinction = wavelength, extinction
        design, buf_provenance, recipe = (
            classify_buffer(sheet, str(experiment["xls_file"]))
            if sheet is not None else ("unknown", "assumed", None))

        notes = []
        for field, (value, reason) in rulings.items():
            if field == "abs_nm":
                wavelength = value
            elif field == "e_declared":
                extinction = value
            provenance[field] = "ruling"
            notes.append(f"{field}: {reason}")

        rows.append({
            "experiment": number,
            "substrate": resolved["substrate"],
            "buffer": resolved["buffer"],
            "pH": resolved["pH"],
            "T": resolved["T"],
            "has_enzyme": bool(resolved["has_enzyme"]),
            "n_samples": len(experiment["samples"]),
            "abs_nm": wavelength,
            "e_declared": extinction,
            "abs_nm_sheet": sheet_wavelength,
            "e_sheet": sheet_extinction,
            "design": design,
            "buf_provenance": buf_provenance,
            "status": "exclude" if excluded else "use",
            "exclude_reason": KNOWN_EXCLUSIONS.get(number, flag if excluded else ""),
            "accepted_deviations": KNOWN_ACCEPTED_DEVIATIONS.get(number, ("", ""))[0],
            "accepted_deviation_reason": KNOWN_ACCEPTED_DEVIATIONS.get(number, ("", ""))[1],
            "provenance": ";".join(f"{k}={v}" for k, v in provenance.items()),
            "open_questions": " | ".join(open_questions + KNOWN_OPEN_QUESTIONS.get(number, [])),
            "xls_file": experiment["xls_file"],
            "notes": " | ".join(notes + ([f"buffer recipe: {recipe}"] if recipe else [])),
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
