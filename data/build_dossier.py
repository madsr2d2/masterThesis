"""
Builds a review dossier: one HTML page presenting every source of truth for
each experiment side by side, so the per-experiment metadata can be curated by
eye instead of by heuristic.

The compiled dataset's per-cuvette concentrations are the best-verified part of
this project -- 248 columns checked against volume proportionality and 341
cuvette values traced through the recorded dilution chains back to weighed
grams. The per-experiment metadata is the opposite: guessed by a scanner that
cannot fail loudly. This page is aimed at that second layer.

For each experiment it shows, in one place:

  sources     what the filename says, what the sheet declares, and what the
              compiled dataset ended up with -- disagreements highlighted
  recipe      the sheet's own cuvette table, verbatim
  compiled    the rows that reached experiment_data.csv
  curves      every progress curve, with dead and backwards ones marked
  inferred    design (volume vs stock) and buffer-stock provenance, proposed
              by rule for confirmation rather than asserted
  ruling      the fields to fill in on the manifest

Nothing here writes to the manifest. The dossier is for reading; the rulings go
back by hand.

Usage:
    python data/build_dossier.py                        # all experiments
    python data/build_dossier.py --experiments 2,26,146  # a sample
    python data/build_dossier.py --out /tmp/dossier.html
"""
import argparse
import base64
import contextlib
import html
import io
import os
import re
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from build_manifest import read_sheet_optics
from kinetics_io import (EXPERIMENT_CORRECTIONS, find_header_row,
                         parse_experiment_data)
from solution_chemistry import add_solution_columns

MANIFEST_PATH = "data/manifest.csv"
DATASET_PATH = "data/experiment_data.csv"
SHEET_DIR = "data/data"

PALETTE = ["#2f6fb0", "#c0522a", "#3f8a5a", "#8a5aa8", "#b08a2f",
           "#4a8a9a", "#a83f5a", "#6a6a6a", "#5a7a3f"]


# --- reading the sheet ----------------------------------------------------

def read_sheet(xls_file):
    """Returns Sheet1 of one experiment's workbook, headerless."""
    return pd.read_excel(os.path.join(SHEET_DIR, xls_file),
                         sheet_name="Sheet1", header=None)


def header_index(sheet):
    """Returns the cuvette-table header row index, or None."""
    with contextlib.redirect_stdout(io.StringIO()):
        return find_header_row(sheet)


def column_labels(sheet, header_row):
    """
    Joins the header row with the one below it, since several sheets split a
    column label across two rows ('[sub]' over 'mmol/l').
    """
    labels = []
    for column in range(sheet.shape[1]):
        parts = [str(sheet.iat[row, column])
                 for row in range(header_row, min(header_row + 2, len(sheet)))
                 if str(sheet.iat[row, column]) != "nan"]
        labels.append(" ".join(parts).strip())
    return labels


def recipe_rows(sheet, header_row, labels, limit=24):
    """
    Returns the cuvette table as a list of dicts.

    Stops at the first fully blank row, and skips the trailing 'Sum:' rows
    several sheets carry. Those sum rows are why the design classification
    used to read exp 26 -- four identical cuvettes -- as a volume design: the
    sum row contributes a different buffer volume and makes the column look
    like it varies. A row only counts as a cuvette if its total-volume cell is
    filled in, which the sum rows leave blank.
    """
    keep = [i for i, label in enumerate(labels) if label]
    volume_column = next((i for i, label in enumerate(labels)
                          if label.lower().startswith("vol")), None)
    rows = []
    start = header_row + (2 if any(
        str(sheet.iat[header_row + 1, i]) != "nan" and
        not _is_number(sheet.iat[header_row + 1, i]) for i in keep) else 1)
    for row in range(start, min(start + limit, len(sheet))):
        values = [sheet.iat[row, i] for i in keep]
        if all(pd.isna(v) for v in values):
            break
        first = str(sheet.iat[row, keep[0]]).strip().lower()
        if first.startswith("sum"):
            continue
        if volume_column is not None and pd.isna(sheet.iat[row, volume_column]):
            continue
        entry = {labels[i]: sheet.iat[row, i] for i in keep}
        entry["_reference"] = first.startswith("ref")
        rows.append(entry)
    return rows


def measured_rows(recipe):
    """
    The cuvettes that were actually read, excluding the paired no-enzyme
    reference channels.

    The distinction matters for classification: exp 146 runs every sample
    cuvette at 0.5 ml buffer in 1 ml total and every reference at 1.0 ml in
    2 ml, so counting the references makes a fixed-volume experiment look like
    a volume titration.
    """
    return [row for row in recipe if not row.get("_reference")]


def _is_number(value):
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


SHEET_NOTE_PATTERN = re.compile(
    r"^(ph\b.*|t|temp.*|dt \[s\]|method:?|d2o ?%|substrate|exp\. ?#|buf ?[12]?\.?)$")


def sheet_notes(sheet, limit_rows=45):
    """
    Label/value pairs the sheet states in its header block.

    Surfaces things no column check can: several sheets record the pH measured
    before AND after the run, which is the only evidence in the files that the
    single pH in the dataset is an approximation.
    """
    notes = []
    text = sheet.map(lambda v: "" if pd.isna(v) else str(v).strip())
    for row in range(min(limit_rows, sheet.shape[0])):
        for column in range(sheet.shape[1] - 1):
            label = text.iat[row, column]
            if not label or not SHEET_NOTE_PATTERN.match(label.lower()):
                continue
            value = text.iat[row, column + 1]
            if not value and row + 1 < sheet.shape[0]:
                value = text.iat[row + 1, column]   # several sheets stack the
                                                    # value under the label,
                                                    # e.g. 'pH before run.'
            if not value or value == label:
                continue
            pair = (label.rstrip(":"), value)
            if pair not in notes:
                notes.append(pair)
    return notes


def sheet_optics(sheet):
    """Returns (wavelength_nm, e) as the sheet's header declares them."""
    found = {}
    text = sheet.map(lambda v: str(v).strip().lower())
    for row in range(sheet.shape[0]):
        for column in range(sheet.shape[1] - 1):
            label = text.iat[row, column]
            value = pd.to_numeric(pd.Series([sheet.iat[row, column + 1]]),
                                  errors="coerce").iloc[0]
            if pd.isna(value):
                continue
            if label.startswith("abs") and "nm" in label:
                found.setdefault("abs_nm", float(value))
            elif label.startswith("e ") or label.startswith("e["):
                found.setdefault("e", float(value))
    return found.get("abs_nm"), found.get("e")


# --- classification -------------------------------------------------------

def classify(sheet, header_row, labels, xls_file, recipe):
    """
    Proposes a design and a buffer-stock provenance for confirmation.

    design      'stock'  fixed volumes, concentration varied by using
                         different stock solutions -- [buf] is constant
                'volume' buffer volume traded against substrate volume --
                         [buf] genuinely varies across the cuvettes

    buf_provenance
                'declared'   the sheet states [buf] per cuvette
                'filename'   the stock molarity is in the filename
                'sheet-text' it appears only as free text in the sheet
                'assumed'    stated nowhere; (V_buf/2)*100 assumes 0.1 M
    """
    declared = any(label.lower().startswith("[buf") for label in labels)
    volume_label = next((label for label in labels
                         if label.lower().startswith("buf")
                         and "ml" in label.lower()), None)

    design = "unknown"
    measured = measured_rows(recipe)
    if volume_label and measured:
        values = pd.to_numeric(pd.Series([row.get(volume_label) for row in measured]),
                               errors="coerce").dropna()
        design = "volume" if values.nunique() > 1 else "stock"
    elif declared:
        design = "stock"

    if declared:
        provenance = "declared"
    elif re.search(r"(\d+[.,]?\d*)\s*M(?![a-z])", xls_file):
        provenance = "filename"
    elif any(isinstance(v, str) and re.search(r"\(\s*\d+[.,]?\d*\s*M\s*\)", v)
             for v in sheet.values.ravel()):
        provenance = "sheet-text"
    else:
        provenance = "assumed"
    return design, provenance


def curve_diagnostics(experiment_number):
    """Returns a DataFrame of per-curve diagnostics, or None if no txt found."""
    for name in sorted(os.listdir(SHEET_DIR)):
        if not name.endswith(".txt"):
            continue
        parsed = parse_experiment_data(os.path.join(SHEET_DIR, name))
        if parsed is None or parsed.get("num") != experiment_number:
            continue
        rows = []
        for index, (sample_name, sample) in enumerate(parsed["samples"].items(), start=1):
            time = np.asarray(sample["time"], dtype=float)
            values = np.asarray(sample["values"], dtype=float)
            rows.append({
                "sample": index, "name": sample_name, "points": len(time),
                "dt": float(np.median(np.diff(time))) if len(time) > 1 else np.nan,
                "duration_min": float(time[-1] - time[0]) / 60.0 if len(time) else np.nan,
                "start": float(values[0]), "end": float(values[-1]),
                "max": float(values.max()), "net": float(values[-1] - values[0]),
                "time": time, "values": values,
            })
        return pd.DataFrame(rows), parsed.get("date")
    return None, None


def curve_flags(row):
    """Returns the list of problems this curve has, if any."""
    flags = []
    if row["net"] < -1e-9:
        flags.append(f"ends {row['net']:+.3f} below start")
    if row["max"] < 0:
        flags.append("never rises above baseline")
    if row["points"] < 20:
        flags.append(f"only {row['points']} points")
    return flags


# --- rendering ------------------------------------------------------------

def plot_curves(diagnostics):
    """Renders the progress curves to a base64 PNG."""
    figure, axis = plt.subplots(figsize=(7.2, 3.2), dpi=110)
    for index, row in diagnostics.iterrows():
        flagged = bool(curve_flags(row))
        axis.plot(np.asarray(row["time"]) / 60.0, row["values"],
                  color=PALETTE[index % len(PALETTE)],
                  linestyle="--" if flagged else "-",
                  linewidth=1.1 if flagged else 1.5,
                  alpha=0.75 if flagged else 1.0,
                  label=f"{row['sample']}" + (" ⚠" if flagged else ""))
    axis.axhline(0, color="#999", linewidth=0.6, zorder=0)
    axis.set_xlabel("time (min)")
    axis.set_ylabel("absorbance")
    axis.legend(fontsize=7, ncol=min(len(diagnostics), 8), frameon=False,
                loc="upper left")
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png")
    plt.close(figure)
    return base64.b64encode(buffer.getvalue()).decode("ascii")


def table(rows, columns=None, classes=""):
    """Renders a list of dicts as an HTML table."""
    if not rows:
        return "<p class='none'>nothing to show</p>"
    columns = columns or list(rows[0].keys())
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in columns)
    body = []
    for row in rows:
        cells = []
        for column in columns:
            value = row.get(column, "")
            if isinstance(value, float):
                value = "" if pd.isna(value) else f"{value:g}"
            elif pd.isna(value) if not isinstance(value, (list, str)) else False:
                value = ""
            cells.append(f"<td>{html.escape(str(value))}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return (f"<table class='{classes}'><thead><tr>{head}</tr></thead>"
            f"<tbody>{''.join(body)}</tbody></table>")


def parse_filename(name):
    """Pulls what the filename asserts. Filenames are ground truth here."""
    facts = {}
    match = re.search(r"pH[=_ ]?(\d+[.,]\d+|\d+)", name, re.I)
    if match:
        facts["pH"] = float(match.group(1).replace(",", "."))
    match = re.search(r"t=(\d+)", name, re.I)
    if match:
        facts["T"] = float(match.group(1))
    if re.search(r"4OMe|pMeOBnOH|4MeOBnOH", name, re.I):
        facts["substrate"] = "4OMe-BnOH"
    elif re.search(r"4-brom", name, re.I):
        facts["substrate"] = "4Br-BnOH"
    elif re.search(r"BnOH|BnO", name, re.I):
        facts["substrate"] = "BnOH"
    if re.search(r"no[_ ]E|NO[_ ]E", name):
        facts["has_enzyme"] = False
    elif re.search(r"with[_ ]E", name, re.I):
        facts["has_enzyme"] = True
    match = re.search(r"(\d+[.,]?\d*)\s*M(?![a-z])", name)
    if match:
        facts["buffer_stock_M"] = match.group(1).replace(",", ".")
    return facts


def render_experiment(number, manifest, dataset):
    """Renders one experiment's section."""
    declared = manifest.loc[number]
    xls_file = declared["xls_file"]
    rows = dataset[dataset.experiment == number]

    sheet = read_sheet(xls_file)
    header_row = header_index(sheet)
    labels = column_labels(sheet, header_row) if header_row is not None else []
    recipe = recipe_rows(sheet, header_row, labels) if header_row is not None else []
    design, provenance = classify(sheet, header_row, labels, xls_file, recipe)
    sheet_nm, sheet_e = read_sheet_optics(sheet)
    notes = sheet_notes(sheet)
    from_filename = parse_filename(xls_file)

    # --- sources side by side
    source_rows = []
    for field, filename_value, sheet_value, dataset_value in [
        ("substrate", from_filename.get("substrate"), None, rows.substrate.iloc[0]),
        ("pH", from_filename.get("pH"), None, rows.pH.iloc[0]),
        ("T", from_filename.get("T"), None, rows["T"].iloc[0]),
        ("has_enzyme", from_filename.get("has_enzyme"), None,
         bool((rows["[enz]"] > 0).any())),
        ("buffer", None, None, rows.buffer.iloc[0]),
        ("wavelength nm", None, sheet_nm, rows["abs"].iloc[0]),
        ("e mM-1cm-1", None, sheet_e, rows.e.iloc[0]),
        ("buffer stock M", from_filename.get("buffer_stock_M"), None,
         "0.1 (assumed)" if provenance == "assumed" else None),
    ]:
        known = [v for v in (filename_value, sheet_value, dataset_value) if v is not None]
        conflict = len({str(v) for v in known}) > 1 and len(known) > 1
        source_rows.append({
            "field": field,
            "filename says": "" if filename_value is None else filename_value,
            "sheet says": "" if sheet_value is None else sheet_value,
            "dataset has": "" if dataset_value is None else dataset_value,
            "_conflict": conflict,
        })

    source_html = ["<table class='sources'><thead><tr>"
                   "<th>field</th><th>filename says</th><th>sheet says</th>"
                   "<th>dataset has</th></tr></thead><tbody>"]
    for row in source_rows:
        klass = " class='conflict'" if row["_conflict"] else ""
        source_html.append(
            f"<tr{klass}><td>{html.escape(row['field'])}</td>"
            f"<td>{html.escape(str(row['filename says']))}</td>"
            f"<td>{html.escape(str(row['sheet says']))}</td>"
            f"<td>{html.escape(str(row['dataset has']))}</td></tr>")
    source_html.append("</tbody></table>")

    # --- recipe as the sheet states it
    display = []
    for row in recipe:
        entry = {k: v for k, v in row.items() if k != "_reference"}
        entry["role"] = "reference" if row.get("_reference") else "measured"
        display.append(entry)
    recipe_html = table(display, classes="recipe") if display else \
        "<p class='none'>no cuvette table found on Sheet1</p>"

    # --- compiled rows
    compiled = add_solution_columns(rows)
    compiled_columns = ["sample", "[enz]", "[buf]", "[h2o2]", "[sub]", "I", "[HOO-]"]
    compiled_html = table(compiled[compiled_columns].round(4).to_dict("records"),
                          classes="compiled")

    # --- curves
    diagnostics, collected = curve_diagnostics(number)
    if diagnostics is not None and len(diagnostics):
        image = plot_curves(diagnostics)
        curve_table = []
        for _, row in diagnostics.iterrows():
            flags = curve_flags(row)
            curve_table.append({
                "sample": row["sample"], "points": row["points"],
                "dt s": row["dt"], "duration min": round(row["duration_min"], 1),
                "start": round(row["start"], 4), "end": round(row["end"], 4),
                "net": round(row["net"], 4),
                "flags": "; ".join(flags),
            })
        curves_html = (f"<img alt='progress curves for experiment {number}' "
                       f"src='data:image/png;base64,{image}'>"
                       + table(curve_table, classes="curves"))
        bad = sum(1 for _, row in diagnostics.iterrows() if curve_flags(row))
    else:
        curves_html = "<p class='none'>no time series found</p>"
        bad = 0

    # --- flags
    flags = []
    if declared["status"] == "exclude":
        flags.append(("excluded", declared["exclude_reason"]))
    if isinstance(declared.get("open_questions"), str) and declared["open_questions"].strip():
        for question in declared["open_questions"].split("|"):
            flags.append(("open question", question.strip()))
    if number in EXPERIMENT_CORRECTIONS:
        flags.append(("corrected", "values overridden by kinetics_io.EXPERIMENT_CORRECTIONS: "
                      + ", ".join(EXPERIMENT_CORRECTIONS[number])))
    if provenance == "assumed":
        flags.append(("assumption", "[buf] computed as (V_buf/2)*100; the 0.1 M stock "
                                    "is stated nowhere in the filename or the sheet"))
    if bad:
        flags.append(("curves", f"{bad} of {len(diagnostics)} curve(s) flagged below"))
    over = compiled["I"] > 100
    if over.any():
        flags.append(("ionic strength", f"I = {compiled['I'].max():.0f} mM exceeds the "
                                        f"Debye-Huckel range, so [HOO-] carries a "
                                        f"systematic error here"))

    flag_html = "".join(
        f"<div class='flag'><span class='tag'>{html.escape(kind)}</span>"
        f"{html.escape(str(text))}</div>" for kind, text in flags
    ) or "<div class='flag ok'><span class='tag ok'>clean</span>nothing flagged</div>"

    notes_html = ("<table class='notes'><tbody>" + "".join(
        f"<tr><td>{html.escape(label)}</td><td>{html.escape(value)}</td></tr>"
        for label, value in notes) + "</tbody></table>") if notes else \
        "<p class='none'>no labelled header values found</p>"

    if provenance == "declared":
        stock_note = "n/a -- [buf] read per cuvette from the sheet"
    elif provenance == "filename":
        stock_note = str(from_filename.get("buffer_stock_M", "?"))
    elif provenance == "sheet-text":
        found = sorted({m.group(1) for v in sheet.values.ravel()
                        if isinstance(v, str)
                        for m in [re.search(r"\(\s*(\d+[.,]?\d*)\s*M\s*\)", v)] if m})
        stock_note = ", ".join(found) + "  (per cuvette, from sheet text)"
    else:
        stock_note = "0.1 assumed -- stated nowhere"

    status_class = "exclude" if declared["status"] == "exclude" else "use"
    return f"""
<section id="exp{number}">
  <div class="gutter">{number}<small>exp</small></div>
  <div class="record">
  <h2>{html.escape(str(declared['substrate']))} in {html.escape(str(declared['buffer']))},
      pH {html.escape(str(declared['pH']))}, {html.escape(str(declared['T']))}&deg;C</h2>
  <p class="file">{html.escape(str(xls_file))}
     {'&middot; collected ' + html.escape(str(collected)) if collected else ''}</p>
  <div class="badges">
      <span class="badge {status_class}">{html.escape(str(declared['status']))}</span>
      <span class="badge design">{html.escape(design)} design</span>
      <span class="badge prov-{provenance}">[buf] {html.escape(provenance)}</span>
      <span class="badge">{len(rows)} cuvettes</span>
  </div>

  {flag_html}

  <div class="grid">
    <div>
      <h3>Sources</h3>
      {''.join(source_html)}
      <h3>Compiled rows</h3>
      {compiled_html}
      <h3>What else the sheet states</h3>
      {notes_html}
    </div>
    <div>
      <h3>Curves</h3>
      {curves_html}
    </div>
  </div>

  <h3>The sheet's own cuvette table</h3>
  <div class="scroll">{recipe_html}</div>

  <h3>Ruling</h3>
  <div class="ruling"><table><tbody>
    <tr><td>design</td><td class="proposed">{html.escape(design)}</td><td class="blank"></td></tr>
    <tr><td>buf_provenance</td><td class="proposed">{html.escape(provenance)}</td><td class="blank"></td></tr>
    <tr><td>buffer stock (M)</td><td class="proposed">{html.escape(stock_note)}</td><td class="blank"></td></tr>
    <tr><td>status</td><td class="proposed">{html.escape(str(declared['status']))}</td><td class="blank"></td></tr>
    <tr><td>open questions resolved?</td><td class="proposed">{
        len([f for f in flags if f[0] == 'open question'])} outstanding</td><td class="blank"></td></tr>
  </tbody></table></div>
  </div>
</section>
"""


STYLE = """
/* Palette grounded in the subject: a UV spectrophotometer readout on lab
   paper. Neutrals carry a cool bias toward the 285 nm accent rather than
   sitting at pure grey; the warn hue is the burnt orange of a flagged trace. */
:root {
  --ground:#f6f6f3; --panel:#edeee9; --sunk:#e4e5df;
  --ink:#16181c; --body:#33363d; --dim:#71757e;
  --rule:#d6d7d0; --hair:#e6e7e1;
  --accent:#33459b; --accent-soft:#e7e9f4;
  --warn:#9d4419; --warn-soft:#f6ebe3;
  --ok:#3a6a48; --ok-soft:#e8f0e9;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --ground:#111317; --panel:#191c21; --sunk:#1f232a;
    --ink:#eceef1; --body:#c3c7cf; --dim:#858b96;
    --rule:#2b3038; --hair:#22262d;
    --accent:#8f9fe4; --accent-soft:#1c2140;
    --warn:#e08a52; --warn-soft:#2c1f16;
    --ok:#7cb28c; --ok-soft:#16241b;
  }
}
:root[data-theme="dark"] {
  --ground:#111317; --panel:#191c21; --sunk:#1f232a;
  --ink:#eceef1; --body:#c3c7cf; --dim:#858b96;
  --rule:#2b3038; --hair:#22262d;
  --accent:#8f9fe4; --accent-soft:#1c2140;
  --warn:#e08a52; --warn-soft:#2c1f16;
  --ok:#7cb28c; --ok-soft:#16241b;
}

body { background:var(--ground); color:var(--body);
       font:400 15px/1.6 "IBM Plex Serif", Georgia, serif;
       margin:0; padding:0 28px 96px; }
.wrap { max-width:1220px; margin:0 auto; }

h1 { font:600 30px/1.15 "IBM Plex Sans Condensed", "IBM Plex Sans",
     system-ui, sans-serif; color:var(--ink); letter-spacing:-.01em;
     margin:40px 0 10px; text-wrap:balance; }
p.lede { max-width:66ch; color:var(--body); margin:0 0 22px; }
p.lede em { color:var(--ink); font-style:italic; }

/* Section = one specimen sheet: a number in the gutter, the record beside it */
section { display:grid; grid-template-columns:88px minmax(0,1fr); gap:0 22px;
          border-top:2px solid var(--ink); padding-top:18px; margin:0 0 52px; }
@media (max-width:820px){ section { grid-template-columns:minmax(0,1fr); } }
.gutter { font:600 40px/1 "IBM Plex Sans Condensed", system-ui, sans-serif;
          color:var(--ink); font-variant-numeric:tabular-nums;
          letter-spacing:-.02em; position:sticky; top:12px; align-self:start; }
.gutter small { display:block; font:600 10px/1.4 "IBM Plex Mono", monospace;
                color:var(--dim); text-transform:uppercase; letter-spacing:.12em;
                margin-top:4px; }
.record { min-width:0; }

h2 { font:600 19px/1.25 "IBM Plex Sans Condensed", system-ui, sans-serif;
     color:var(--ink); margin:0 0 3px; }
h3 { font:600 10.5px/1 "IBM Plex Mono", monospace; text-transform:uppercase;
     letter-spacing:.13em; color:var(--dim); margin:22px 0 7px; }
p.file { font:400 11.5px/1.5 "IBM Plex Mono", monospace; color:var(--dim);
         margin:0 0 14px; word-break:break-all; }

/* Badges encode the three classifications, not decoration */
.badges { display:flex; flex-wrap:wrap; gap:6px; margin:8px 0 12px; }
.badge { font:600 10px/1.7 "IBM Plex Mono", monospace; text-transform:uppercase;
         letter-spacing:.09em; padding:1px 9px; border-radius:2px;
         background:var(--sunk); color:var(--dim); }
.badge.use { background:var(--accent-soft); color:var(--accent); }
.badge.exclude, .badge.prov-assumed { background:var(--warn-soft); color:var(--warn); }

.flag { display:grid; grid-template-columns:112px minmax(0,1fr); gap:10px;
        background:var(--panel); border-left:3px solid var(--warn);
        padding:6px 12px; margin:0 0 3px; font-size:13.5px; }
.flag.ok { border-left-color:var(--ok); background:var(--ok-soft); }
.tag { font:600 10px/1.9 "IBM Plex Mono", monospace; text-transform:uppercase;
       letter-spacing:.08em; color:var(--warn); }
.flag.ok .tag { color:var(--ok); }

.grid { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1fr);
        gap:26px; align-items:start; margin-top:4px; }
@media (max-width:900px){ .grid { grid-template-columns:minmax(0,1fr); } }

table { border-collapse:collapse; width:100%;
        font:400 12.5px/1.55 "IBM Plex Mono", ui-monospace, monospace;
        font-variant-numeric:tabular-nums; }
th, td { padding:3px 12px 3px 0; text-align:left; white-space:nowrap;
         border-bottom:1px solid var(--hair); }
th { font-weight:600; font-size:10px; text-transform:uppercase;
     letter-spacing:.09em; color:var(--dim); border-bottom:1px solid var(--rule); }
td { color:var(--body); }
tbody tr:last-child td { border-bottom:1px solid var(--rule); }

/* A conflict is marked, not washed: the stripe carries it */
tr.conflict td { color:var(--warn); font-weight:600; background:var(--warn-soft); }
tr.conflict td:first-child { box-shadow:inset 3px 0 0 var(--warn); padding-left:9px; }

.notes td:first-child, .sources td:first-child { color:var(--dim); }
.scroll { overflow-x:auto; padding-bottom:2px; }

img { display:block; width:100%; height:auto; margin:0 0 10px;
      border:1px solid var(--rule); background:#fff; }

/* The ruling block is a form: proposed value, then a rule to write on */
.ruling { background:var(--panel); padding:12px 14px; margin-top:6px; }
.ruling table { font-size:12.5px; }
.ruling td { border-bottom:none; padding:6px 12px 6px 0; white-space:normal; }
.ruling td:first-child { color:var(--dim); width:170px; }
.proposed { color:var(--ink); font-weight:600; width:46%; }
.blank { border-bottom:1px solid var(--rule); min-width:150px; }
.none { color:var(--dim); font-style:italic; font-size:13px; }

nav { position:sticky; top:0; z-index:5; background:var(--ground);
      border-bottom:1px solid var(--rule); padding:9px 0 8px; margin:0 0 30px;
      font:400 12px/1.7 "IBM Plex Mono", monospace; }
nav strong { color:var(--dim); font-weight:600; text-transform:uppercase;
             letter-spacing:.09em; font-size:10px; margin-right:10px; }
nav a { color:var(--accent); text-decoration:none; padding:1px 7px;
        border:1px solid var(--rule); border-radius:2px; margin-right:4px;
        display:inline-block; }
nav a:hover, nav a:focus-visible { background:var(--accent-soft);
        border-color:var(--accent); outline:none; }
"""


def build(experiments=None, out_path="dossier.html",
          dataset_path=DATASET_PATH, manifest_path=MANIFEST_PATH):
    """
    Writes the dossier and returns its path.

    Args:
        experiments (list[int] or None): Which experiments; None means all.
        out_path (str): Where to write the HTML.
        dataset_path, manifest_path (str): Inputs.

    Returns:
        str: out_path
    """
    manifest = pd.read_csv(manifest_path).set_index("experiment")
    dataset = pd.read_csv(dataset_path)
    numbers = sorted(experiments) if experiments else sorted(manifest.index)

    sections = []
    for number in numbers:
        print(f"  experiment {number}", flush=True)
        sections.append(render_experiment(number, manifest, dataset))

    links = "".join(f'<a href="#exp{n}">{n}</a>' for n in numbers)
    page = f"""<title>Kinetics Experiment Dossier</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans+Condensed:wght@600&family=IBM+Plex+Serif:ital,wght@0,400;1,400&display=swap">
<style>{STYLE}</style>
<div class="wrap">
<h1>Experiment review dossier</h1>
<p class="lede">Every source of truth for each experiment, side by side. The
per-cuvette concentrations below are already verified against the volume tables
and the recorded dilution chains; what needs a human is the per-experiment
metadata &mdash; the design, the buffer stock, and whether the run should be
used at all. Conflicts between sources are highlighted. Values under
<em>Ruling</em> are proposed by rule, not asserted.</p>
<nav><strong>{len(numbers)} experiments</strong>{links}</nav>
{''.join(sections)}
</div>
"""
    with open(out_path, "w") as handle:
        handle.write(page)
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--experiments", default=None,
                        help="Comma-separated experiment numbers; default all.")
    parser.add_argument("--out", default="dossier.html")
    parser.add_argument("--csv", default=DATASET_PATH)
    parser.add_argument("--manifest", default=MANIFEST_PATH)
    args = parser.parse_args()

    chosen = ([int(x) for x in args.experiments.split(",")]
              if args.experiments else None)
    path = build(chosen, args.out, args.csv, args.manifest)
    print(f"\nwrote {path} ({os.path.getsize(path) / 1e6:.2f} MB)")
