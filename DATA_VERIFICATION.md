# Data Verification Log

Record of checks performed on the raw kinetics dataset (`data/data/*.txt` + `*.xls`)
and its derived artifacts, and any corrections made as a result. New entries go at
the top.

---

## 2026-08-29 — Raw [P] time-series audit

Audited the raw progress-curve data (`data/data/*.txt`, parsed by
`parse_experiment_data`) that feeds the `[P]` column, independently of the
metadata audit above.

### 1. Structural integrity — clean

Reparsed all 98 `.txt` files (443 sample-series) standalone. No parse failures,
no empty series, no backwards or duplicate timestamps, uniform sampling interval
within every series. Point counts range 10–481, all plausible.

### 2. Reaction-direction check — confirms existing fixes, no new cases

Flagged every sample series with a net-negative absorbance trend (signal
decreasing over time, which is wrong for product formation). This reproduced
*exactly* the set already special-cased in `clean_experiment_dataframe`:
experiments 50 and 85 (already sign-flipped), 58/77/78/79/84 (already removed),
experiment 128 samples 2 and 3 (already removed), and experiment 80 (already
dropped by the "remove all Carbonate buffer" filter). No new sign-flip
candidates found.

### 3. Physical-plausibility check — clean

Converted each curve's amplitude to Δ[P] via its extinction coefficient (`e`)
and confirmed it never exceeds the sample's initial `[sub]` (product can't
exceed starting substrate) across all 443 samples. Zero violations.

### 4. Experiments 32 and 82 are dead runs — removed

Experiments 72 and 82 show essentially zero signal (amplitude 0.001–0.008
absorbance units, vs. a dataset median of 0.044) in *every* sample, despite
normal, nonzero `[enz]`/`[h2o2]`/`[sub]` per the (trusted) xls values. Compared
against same-day, same-instrument sibling runs:
- Exp 72 (6/8/2010, 7:13pm) vs exp 71 (same day, 4:51pm): exp 71 shows real
  signal (amplitude 0.017–0.031) at *lower* substrate concentrations than exp 72
  used. Exp 72 is flat despite higher substrate.
- Exp 82 (6/14/2010, 1:55pm) vs exp 83 (same day, 5:03pm): exp 83 shows strong,
  clean signal (amplitude 0.031–0.043) at comparable/lower substrate. Exp 82 is
  flat despite comparable/higher substrate.

Same-day neighbors with comparable reagents produced normal curves, so this
isn't a bad-reagent-stock day or a data-entry problem — the working hypothesis
is a failed individual run (e.g. an inactive enzyme aliquot), the same failure
mode presumed for the already-excluded 58/77/78/79/84.

**Fix applied:** added 72 and 82 to `experiments_to_remove` in
`clean_experiment_dataframe`, alongside 58/77/78/79/84. See
`masterThesis.ipynb`, cell 16.

### 5. Open item, not chased: `.txt` "Substrate Conc." line disagrees with xls `[sub]`

Each `.txt` file has its own `Substrate Conc.` annotation, independent of the
`.xls` parsing and **not read by the pipeline at all**. Cross-checked it against
the xls-derived `[sub]` anyway as a bonus sanity check: 51 of 313 comparable
rows disagree, with no single clean pattern (includes a clear raw transcription
typo — exp 30 sample 3 literally reads `2856.0000 mmol/l` against neighbors
`0.0953, 0.1906, _, 0.3812`, almost certainly meant `0.2856`).

**Disposition: the `.xls` files are the source of truth for concentrations —
this `.txt` field was not reliably filled in and should not be trusted.** Not
worth further investigation since the pipeline never reads it; logged here only
in case it resurfaces.

---

## 2026-08-29 — Metadata sweep + buffer-titration fix

### 1. `data/experiment_data.csv` verified against the raw `.xls` files

Re-implemented the notebook's own extraction pipeline (`parse_experiment_data`,
`find_and_parse_experiment_file`, `find_header_row`, `find_numeric_values_below_header`,
`find_pH_value_in_range`, `find_temperature_value_in_range`, `find_buffer_type`,
`find_substrate_type` — notebook cells 3, 5–11) as a standalone script and ran it
fresh against all 98 experiment `.txt`/`.xls` pairs in `data/data/`.

**Result: all 443 rows match `experiment_data.csv` exactly, field-by-field
(`substrate`, `abs`, `e`, `buffer`, `pH`, `T`, `[enz]`, `[buf]`, `[h2o2]`, `[sub]`).
Zero mismatches.** Spot-checked the raw cell layout of experiment 10 by hand to
confirm the extraction logic itself (not just internal consistency) is reading the
right cells.

No warnings or extraction failures were emitted across any of the 98 files.

### 2. `data/data.toml` is stale and unreliable — do not use

- Does not parse as TOML: 8 experiments have `pH = [6.71, 6.71, 6.71. 6.71]` (a
  stray period instead of a comma), and `experiment_022.references` has
  `pH = [.50, ...]` (missing leading digit).
- Not referenced anywhere in `masterThesis.ipynb` — the notebook builds its
  dataframe directly from `data/data/`, never from `data.toml` or
  `experiment_data.csv`.
- Where it disagrees with the raw `.xls`/filename ground truth, it's the one
  that's wrong:
  - exp 2: toml says pH 7.0; raw filename and csv both say 6.71 (correct).
  - exp 38: toml lists 4 samples; the raw `.txt` file only has 3
    (`Sample001–003`) — toml has a phantom 4th sample.
- Decision: **not fixing this file for now.** It's slated for a full rewrite, so
  effort goes into the source data instead. Do not treat it as authoritative in
  the meantime.

### 3. Buffer-titration experiments were silently flattened — fixed

Experiments 32, 34, 35, 36, 37 are buffer-concentration titrations: their `.xls`
files label each sample `1 (0.1M)`, `2 (0.2M)`, `3 (0.3M)`, `4 (0.4M)` (confirmed
with the user: these are the actual buffer **stock** concentrations used). But
that molarity only exists as text in the label column — there's no numeric `[buf]`
column like there is for `[sub]`/`[h2o2]`. The generic `[buf]` extraction fallback
(`find_numeric_values_below_header`, notebook cell 7) instead derives `[buf]` from
sample volume assuming a fixed 100 mM buffer:
```python
adjusted_value = round((value / 2) * 100, 3)  # assumes the buffer is 100 mM
```
Since the pipetted volume is constant across all 4 samples in these experiments
(only the stock concentration changes), this produces an identical, wrong
`[buf] = 50.0` for every sample, regardless of which stock was actually used.

`clean_experiment_dataframe` (notebook cell 16) already had a hardcoded override
for this exact bug on experiments 34 and 35/36 — but 32 and 37 were missed and
still carried the flat wrong value.

**Confirmed these were meant to vary, not a mislabeled flat-concentration run.**
This experiment type is a *buffer-concentration titration*, distinct from the
bulk of the dataset (which titrate substrate at a nominal fixed buffer level).
Evidence:
- Filenames encode a list of four distinct concentrations, not a single value —
  `Phosphat_0.1_0.2_0.3_0.4` (32, 35, 36, 37) and
  `Phosphat_0.05_0.025_0.0125_0.00625` (34, a 2-fold serial dilution) — unlike
  ordinary experiment filenames, which name just one buffer.
- The raw sheet labels each sample individually: `1 (0.1M)`, `2 (0.2M)`,
  `3 (0.3M)`, `4 (0.4M)`.
- `Sub [ml]` and `Buf [ml]` volumes are both held constant across all 4 samples
  (e.g. `Sub [ml] = 0.8`, `Buf [ml] = 1` for every sample) — the opposite of a
  normal substrate titration (e.g. exp 10), where `Sub [ml]` climbs
  (`0.2, 0.4, 0.6, 0.8`) and `Buf [ml]` shrinks to compensate
  (`1.6, 1.4, 1.2, 1.0`). Here substrate is pinned and only the buffer stock
  pipetted in differs sample-to-sample — the signature of buffer concentration
  being the deliberately titrated variable.
- User confirmed 0.1M/0.2M/0.3M/0.4M as the actual buffer stock concentrations
  used.

**Fix applied:** extended the existing override in `clean_experiment_dataframe`
to cover experiments 32 and 37 as well, using the same `[100, 200, 300, 400]` mM
convention already established for 35/36. See `masterThesis.ipynb`, cell 16.

### 4. Full sweep for the same class of bug across all 98 experiments

Checked every experiment for any concentration field (`[enz]`, `[buf]`, `[h2o2]`,
`[sub]`) that's flat across all samples, then checked whether the filename/raw
label implies it should vary (the tell-tale sign of this bug class).

- **Experiment 26**: all four fields flat. Checked the raw sheet directly —
  genuinely a 4-way replicate design at fixed conditions (no enzyme; comparing
  "with H2O2" vs "without H2O2" blanks), not a titration. Not a bug.
- **Experiments 32, 34, 35, 36, 37**: buffer titrations — fixed (see §3).
- **Experiments 127–131** ("sub H2O2" filenames): `[sub]` flat, `[h2o2]` varies —
  correct by design, these substitute an H2O2 concentration range for the usual
  substrate range.
- No other experiment showed an unexplained flat field. No further metadata bugs
  found in this pass.

### 5. Experiment 32 has a stray duplicate `.xls` file

`data/data/` contains two files for experiment 32:
`mads_t032_..._Phosphat_0.1_0.2_0.xls` (looks like a truncated/interrupted save)
and `mads_t032_..._0.1_0.2_0.3_0.4_..._with_NO_E.xls` (the complete one). Checked
both — they contain byte-identical concentration tables, so whichever one the
pipeline's glob picks doesn't affect results. Not fixed/removed, just noted here
in case the stray file causes confusion later.

---

## Open items

- `data.toml` still doesn't parse and isn't reconciled with the verified csv
  (deprioritized — planned rewrite of the toml-generating script).
- Metadata (pH, T, buffer, substrate, `[enz]`/`[buf]`/`[h2o2]`/`[sub]`) and the
  raw `[P]` time-series data have both now been audited (see entries above).
  Remaining known-unreliable field: the `.txt` files' own `Substrate Conc.`
  line — not used by the pipeline, don't trust it if it resurfaces.
