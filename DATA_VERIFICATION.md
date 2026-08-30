# Data Verification Log

Record of checks performed on the raw kinetics dataset (`data/data/*.txt` + `*.xls`)
and its derived artifacts, and any corrections made as a result. New entries go at
the top.

---

## 2026-08-30 — Mechanism research: consequences for how this data can be fitted

Three literature research passes were run to pin down the reaction mechanism
before fitting (full write-up, reasoning and 51 references in `MECHANISM.md`).
Several findings bear directly on the data itself and on what can and cannot be
concluded from it, so they are recorded here too.

### 1. The four buffers are chemically different reagents, not just pH setters

This is the most important finding for data interpretation. It also gives a
mechanistic explanation for the buffer-dependent behaviour noticed by eye in
the dataset.

- **Boric buffer points should be treated as suspect.** Borate does three
  separate things to this chemistry: it forms **peroxoborate** with H2O2
  (significant above pH ~7.7), whose anions are much faster oxidants than H2O2
  and "deliver the hydroperoxide anion at a lower pH than when H2O2 is used";
  it generates **dioxaborirane**, a highly reactive cyclic peroxide that is a
  competing oxidant unrelated to the catalyst; and — most damaging — **boric
  acid catalyses peroxyacid hydrolysis ~12-fold, with a maximum at pH 8.4–9**.
  The dataset's Boric experiments span pH 8.46–10.34, i.e. straight through
  that maximum. If the proposed mechanism's peracid intermediate is real,
  borate buffer is actively destroying it.
- **Carbonate buffer points should also be treated as suspect.** Bicarbonate +
  H2O2 forms **peroxymonocarbonate (HCO4-)**, a two-electron oxidant ~300x
  faster than H2O2 for sulfide oxidation, formed within minutes near neutral
  pH. In the dataset's 7 Carbonate experiments (pH 9.40–11.84) the effective
  oxidant is partly HCO4-, not H2O2.
- **Phosphate catalyses the first step of the mechanism directly.** H2O2
  addition to a carbonyl is subject to **both general acid and general base
  catalysis** (Sander & Jencks 1968), so buffer concentration is a genuine
  kinetic variable, not a nuisance parameter.
- **Pyrophosphate is probably chelating trace metals** (reasoning, not
  sourced) — trace Fe/Cu catalyse H2O2 decomposition, which is why dioxirane
  papers routinely add EDTA. Pyrophosphate chelates; phosphate and carbonate
  do not, to the same degree.

**Consequence:** rate constants from different buffer systems at the same
nominal pH are **not directly comparable**. Any pH-rate profile built by
pooling across buffers is confounded. This should be stated explicitly
alongside any such plot.

### 2. `[buf]` and `[sub]` are collinear within every titration experiment — a buffer-concentration effect cannot be isolated from the existing data

Checked directly: in every titration experiment, `[buf]` decreases in lockstep
as `[sub]` increases across samples 1→4 (substrate stock dilutes the buffer).
Example, experiment 2: `[buf]` 80→70→60→50 mM while `[sub]` 1.52→3.04→4.56→6.08
mM. 50 of 98 experiments vary `[buf]` within the experiment, and **not one of
them holds `[sub]` fixed while doing so**.

So any within-experiment "buffer effect" is perfectly confounded with the
substrate effect. Separating them requires either a multi-variable regression
of v0 against `[sub]`, `[h2o2]`, `[enz]`, `[buf]` jointly across experiments,
or (cleanly, but this needs new bench work) a buffer-dilution series at fixed
pH, fixed `[sub]` and fixed ionic strength.

Buffer *type* is the cleaner comparison available now: the pH ranges overlap
across buffers (Phosphate 5.64–8.95, Pyrophosphate 5.47–9.73, Boric 8.46–10.34,
Carbonate 9.40–11.84), so different buffer species can be compared at matched
pH — the classic diagnostic for general acid/base catalysis. Bear §1 in mind
when interpreting the result, though: borate and carbonate bring their own
chemistry, so a buffer-type difference is not automatically general acid/base
catalysis.

### 3. The pH range may not be able to discriminate the mechanism's key branch

The mechanism predicts that dioxirane formation (its central catalytic step) is
pH-controlled, because it needs an anionic peroxide species. But the same pH
dependence arises trivially from HOO- (pKa 11.6) being the nucleophile in the
*addition* step. Both predict rate rising across pH 5.5→11.8. The discriminator
is **where the inflection sits**: near 11.6 points to peroxyanion
nucleophilicity; well below 10 would point to ionization of the tetrahedral
adduct. Worth fitting the pH profile carefully enough to locate the inflection
— but note that no pKa has ever been measured for the relevant adduct, so this
is an open question, not a calibrated test.

### 4. Reminder: the differential-measurement design still governs everything

Nothing in the mechanism work changes the round-2 §4 finding — with-enzyme
`[P]` is already background-subtracted at the source. But it is now clearer
*why* the background is substantial and accelerating: the proposed
catalyst-free loop (aldehyde + HOO- → peracid → oxidises more alcohol) is
autocatalytic on its own, with no catalyst required. Whether that specific
ionic mechanism is right is unresolved (it has no literature precedent — see
`MECHANISM.md`), and a radical/O2-chain alternative is equally consistent with
the data and cannot be excluded without dark/anaerobic controls.

## 2026-08-29 — Round 3: row-level block-structure audit (all 98 experiments)

Round 2 confirmed `.txt`↔`.xls` sample alignment indirectly (structurally on
two examples, statistically via rate-vs-concentration correlation on 12
experiments). This round checks it **directly, on every experiment**, now
that the paired-reference-cuvette design (round 2 §4) is understood: for
each of the 98 experiments, located the `[Enz]` concentration column in the
raw `.xls`, read its full numeric run (not just the truncated
`sample_num`-length slice the pipeline keeps), and split it into the
"leading" block (what the pipeline actually extracts) vs. whatever follows
(the reference block, if recorded).

Checked for two failure modes that would be serious silent labeling errors:
the leading block containing a mix of zero/nonzero `[Enz]` (would mean a
reference row leaked into the real samples), and an all-zero leading block
sitting in front of a nonzero trailing block (would mean the block order is
reversed and the pipeline is reading the reference cuvette instead of the
real one).

**Result: 97 of 98 experiments clean.** No reversed blocks, no reference
row ever extracted in place of a real sample. Two apparent flags turned out
to be benign and are noted for completeness:

- **Exp 6, exp 65** ("no_E" designs): these have no separate `[Enz]`
  concentration column at all (enzyme volume is 0 throughout, so no
  concentration is computed), or only a single shared reference row instead
  of one per sample. Both fall back correctly to `[enz] = 0` — already
  understood from round 2, not a new issue.

- **Exp 128 — real finding, not benign.** Its leading (extracted) block is
  `[0.032, 0.032, 0.032, 0.032, 0.0]` — 4 real enzyme-containing values
  followed by a **0** still inside the block the pipeline keeps, because
  this `.txt` file declares **5** samples while its sibling experiments in
  the same series (127, 129, 130, 131 — same `sub_H2O2` design, same day)
  all declare 4. Read the raw sheet directly (`mads_t128_..._sub H2O2_
  pyrophosphate_02.xls`, rows 32–39): it has the standard 4 "kuv" (real)
  rows followed by 4 "ref" (matched, no-enzyme) rows — but unlike other
  experiments, where the reference block is written to the sheet purely as
  a record and never separately monitored over time, **this run's 5th
  measured channel (`Sample005`) is the first reference row itself
  (`ref 5`)**, not a 5th titration condition. Its `[enz]=0` and its
  `[buf]`/`[h2o2]`/`[sub]` (193.6 / 195.9 / 9.47) exactly match `ref 5`,
  which is the matched blank for `kuv 1`/`kuv 2` (H2O2 = 195.9 mM).

  Current `data/experiment_data.csv` (rows 155–159) already reflects this:
  sample 5 has `[enz] = 0.000` while samples 1–4 have `[enz] = 0.032`.
  Samples 2 and 3 of this experiment were already excluded by the existing
  `clean_experiment_dataframe` (pre-dating this verification effort, for an
  unrelated reason — a backwards/negative time-series trend, see the
  reaction-direction check above). **Resolved in the follow-up investigation
  below** — see "Round 4" for the final disposition of all five samples.

---

## 2026-08-29 — Round 4: plotted experiment 128 directly, resolved samples 2/3/4/5

Built two small reusable modules in `data/` for this and future
investigations, rather than another one-off scratchpad script:

- `data/kinetics_io.py` — the notebook's extraction functions (cells 3,
  5–12 of `masterThesis.ipynb`), copied verbatim so there is one canonical
  implementation, plus a new `load_experiment(exp_num)` convenience wrapper
  that returns one experiment's metadata *and* raw time series together
  (the whole-dataset `populate_experimental_data_from_directory` returns
  metadata only).
- `data/plot_kinetics.py` — `plot_experiment(exp_num, mark_samples={...})`,
  raw absorbance and Beer-Lambert-converted `[P]` side by side, one fixed
  color per sample index, optional dashed styling + legend annotation for
  samples under investigation. Runnable as `python data/plot_kinetics.py
  <exp_num>` from the repo root, or imported directly.

Plotting experiment 128 (`python data/plot_kinetics.py 128`) showed all 5
samples together for the first time and changed the picture from round 3:

| sample | condition | shape | current status |
|---|---|---|---|
| 1 | [enz]=0.032, H2O2=196mM | clean rise, 0.030→0.047 | kept |
| 2 | [enz]=0.032, H2O2=196mM (same as 1) | flat/noisy, net trend −0.008 | **excluded** |
| 3 | [enz]=0.032, H2O2=3.88mM | flat/noisy, net trend −0.004 | was excluded, **reinstated** |
| 4 | [enz]=0.032, H2O2=3.88mM (same as 3) | flat/noisy, net trend +0.003 | kept |
| 5 | [enz]=0, H2O2=196mM (ref channel, see round 3 §4) | flat, net trend +0.002 | **excluded** |

The key observation: samples 3 and 4 are two *independent replicates* of
the same low-H2O2 condition, and both go flat — reproducible, not a
one-off. Meanwhile at the high-H2O2 condition, one replicate (1) works
cleanly and the other (2) doesn't, on the same enzyme aliquot — that's the
signature of a single failed measurement, not a real effect (the aliquot
demonstrably works, since sample 1 proves it).

**Decision:** treat 3.88 mM H2O2 as a genuine, reproducible near-zero-rate
condition rather than two failed measurements. Sample 3's original
exclusion (under the general "any net-negative time-series trend is
unphysical" rule from the earlier reaction-direction check) was conflating
a small, noise-scale negative trend with the large, clearly-wrong-sign
trends that rule was designed to catch (e.g. exp 50/85 pre-correction, or
58/77/78/79/84). Reinstating sample 3 alongside 4 avoids cherry-picking
whichever replicate's noise happened to land on the "nice" side of zero.
Sample 2 stays excluded — unlike 3/4, it disagrees with its own matched
replicate (1) at a condition proven capable of strong signal, which is a
single-measurement failure, not a low-rate result.

**Fix applied** to `clean_experiment_dataframe` (`masterThesis.ipynb`,
cell 16): `exp_num_sample_num` changed from `[[128, 2], [128, 3]]` to
`[[128, 2], [128, 5]]` — drops the sample-3 exclusion, adds the sample-5
(reference channel) exclusion from round 3. Net effect on experiment 128:
samples 1, 3, 4 feed the fit as real data points; 2 and 5 are excluded.

**Open question flagged, not yet chased:** the same "small negative trend
excluded under the blanket backwards-trend rule" pattern that applied to
128,3 may exist elsewhere in the dataset — worth a follow-up pass
distinguishing large-magnitude (real sign-flip) exclusions from
small/noise-scale ones before finalizing the excluded-sample list.

---

## 2026-08-29 — Round 5: broad pass on all "backwards trend" exclusions — closed, one mischaracterization fixed

Followed up on round 4's open question: checked every sample in
experiments 58, 77, 78, 79, 84 (the ones grouped together as pre-existing
"backwards trend" removals) plus 50/85 (sign-flip corrected, as a sanity
check) with a proper significance test instead of eyeballing the sign of
the net trend. For each sample, ran linear regression (`scipy.stats.
linregress`) of raw absorbance vs. time and compared the fitted slope to
its standard error (t-test for slope ≠ 0), plus amplitude vs. the residual
noise level. The question: is each "backwards" trend a real, high-SNR
inverted signal (genuine error, correctly excluded) or a small one
indistinguishable from noise (candidate to reinstate, per the 128,3
reasoning)?

**Result: no new reinstatement candidates.** 128,3 was the exception, not
part of a wider pattern:

- **77, 78:** every sample, large and highly significant negative slopes
  (t down to −95, amplitude 0.02–0.06 vs. noise ~0.001–0.003). Real,
  high-confidence inverted signal — correctly excluded. (Also Carbonate
  buffer, so excluded via that filter regardless of this list.)
- **58:** mixed — 3 samples strongly negative, 1 strongly positive, all
  large real signal (not noise), internally inconsistent between samples
  (possibly a cuvette-position mixup). Correctly excluded; also Carbonate.
- **79:** all 4 samples tiny amplitude (0.001–0.003, at the noise floor) —
  but unlike 128, there's no sibling sample in the group with clean signal
  to contrast against, so this reads as a group dead-run (same flavor as
  72/82) rather than "noise obscuring a real, localized effect." Left
  excluded; also Carbonate.
- **50, 85** (sign-flip corrected, not excluded): every sample in both,
  large amplitude (0.09–0.63), highly significant, consistently inverted.
  Confirms the existing sign-flip fix is well-founded, not a noise
  artifact.

**One mischaracterization found and fixed here:** round 1's original audit
grouped **84** in with the "backwards trend" removals, and round 4 repeated
that framing. It's wrong — checked directly, 84's one recorded sample has
a small but *correctly-signed, positive* trend (t = +22), not a backwards
one. Plotting it (`python data/plot_kinetics.py 84`) shows why it's still
correctly excluded, for an entirely different reason: its `.xls` shows a
planned 6-sample substrate titration (14.3→1.0 mM `[sub]`), but the `.txt`
only ever recorded **one** of those six samples, and that lone curve is a
coarse, quantized 5-step staircase over just 385 s (56 points, dt=7s —
every other experiment in this set uses dt=28–31s). This is an incomplete/
truncated run, not usable kinetic data, regardless of its trend direction.

No code changes from this round — every disposition above already matches
what `clean_experiment_dataframe` currently does. This closes the open
question from round 4.

---

## 2026-08-29 — Round 2: concentration/curve cross-checks

A second verification pass, this time checking things that need *both* data
sources together (round 1 checked the csv against the xls, and the raw
curves against themselves, but not whether a given `.txt` sample and its
`.xls` row are actually the same physical sample).

### 1. `.txt` sample ↔ `.xls` row alignment — verified, no mismatches found

The extraction pipeline assumes the *n*-th `Sample00n` in a `.txt` file is
the *n*-th data row under the concentration header in the matching `.xls`
sheet. This was never directly checked before. Two lines of evidence:

- **Structural**: many `.xls` sheets (e.g. exp 5, exp 135–151) list each
  condition twice — a block of "Kuv." (cuvette, the real monitored sample)
  rows immediately followed by an equal-sized block of "Ref." (no-enzyme
  reference/blank) rows in the *same* concentration columns. Hand-checked
  exp 5 and exp 135 directly against the raw sheet: the `.txt` file only
  declares as many samples as there are Kuv rows, and the pipeline's
  truncate-to-sample-count logic grabs exactly that leading Kuv block, in
  order — it does not spill into the Ref block or reorder anything.
- **Statistical**: for every experiment with exactly one concentration field
  varying across its samples (a clean titration), computed each sample's
  initial reaction rate (linear fit, first 30% of the curve) and checked
  its Spearman correlation with the varied concentration. If `.txt` order
  and `.xls` row order were ever out of sync for a given experiment, this
  would show up as a scrambled or inverted (rate falls as substrate rises)
  correlation. Result, 12 substrate-titration experiments (65, 66, 68, 69,
  70, 71, 73, 74, 75, 76, 83 — excluding the already-known no-enzyme
  controls) all show positive correlation (rho 0.6–1.0), consistent with
  expected Michaelis-Menten behavior. No inverted or scrambled cases.

  A separate "numeric run length" check was tried first (count how many
  numeric cells actually sit under each concentration header, independent
  of sample count, to catch a mismatch directly) but it produces mostly
  false positives, because of the legitimate Kuv+Ref double-block layout
  above — abandoned in favor of the two checks described here.

### 2. Full range/outlier sweep — clean

Checked all 405 non-removed sample rows: pH in [2, 12], T in [0, 80] °C, and
no negative `[enz]`/`[buf]`/`[h2o2]`/`[sub]` values. Zero violations.

The only rows flagged were 72 samples across experiments 6, 23–31, 38–40,
52, 65, 67, 69, 70, 128 (`[enz] = 0` for every sample in the experiment).
Spot-checked several filenames directly — e.g.
`mads_t065_..._BnOH_no_E_H2O2_122.xls` — confirming these are deliberate
no-enzyme blank/background-control runs, not an extraction failure.

### 3. Baseline check (`[P]` near zero at t=0) — clean

Converted each curve to `[P]` and checked the first point sits near the
curve's minimum. 3 borderline flags out of 405, all in samples with a tiny
absolute amplitude (near the noise floor at the lowest substrate
concentration in a titration) — not real calibration offsets, just a small
denominator inflating the fractional metric.

### 4. Explained: `[P]` is a *differential* (reference-subtracted) signal, not an absolute one — this is why no-enzyme "blank" runs can rate as high as or higher than their catalyzed pair

Initially flagged as an anomaly (no-enzyme blanks 65/67/69/70 showing rates
at or above their enzyme-containing pairs 66/68/71) and logged as an open
item. Root cause identified by reading the raw cuvette layout directly, not
just the extracted concentration columns — **every kinetics run in this
dataset is a two-channel differential measurement, not a single absolute
one.** Each `.xls` sheet lists twice as many cuvette rows as the matching
`.txt` file has samples: a leading block of "real" sample rows (labelled
`kuv`/`Kuv.`/`prøve`), immediately followed by an equal-sized block of
reference rows (labelled `ref`/`Ref.`). The instrument (or the operator,
manually) reads each sample cuvette *against* its paired reference cuvette,
and it is that difference — not the sample's raw absorbance — that ends up
in the `.txt` progress curve. The extraction pipeline only ever reads the
leading `sample_num` rows (the real samples); the trailing reference block
is present in the spreadsheet purely as a record of what the reference
cuvette contained, and is correctly never pulled into `[enz]`/`[buf]`/
`[h2o2]`/`[sub]`.

**What differs between experiment types is what the reference cuvette omits:**

- **"with_E" experiments (the bulk of the dataset, e.g. exp 10, exp 66):**
  the reference row for each sample has **identical `Sub [ml]` and
  `H2O2 [ml]` to that sample**, and only `Enz [ml]` set to 0. Example, exp 10
  (`mads_t010_..._with_E.xls`), rows 19–26:

  | kuv | Buf[ml] | H2O[ml] | Enz[ml] | Sub[ml] | H2O2[ml] | role |
  |---|---|---|---|---|---|---|
  | 1–4 | 1.6→1.0 | 0 | 0.1 | 0.2→0.8 | 0.1 | real sample (enzyme present) |
  | 5–8 | 1.6→1.0 | 0.1 | **0** | 0.2→0.8 | 0.1 | matched reference (no enzyme, same Sub/H2O2) |

  So the reported `[P]` for a with_E sample is **already net of the
  non-enzymatic background** — the enzyme-free reaction is running in the
  paired reference cuvette under identical Sub/H2O2 conditions and is
  subtracted out by construction. This is the design in exp 5, exp 10, exp
  66, exp 135–151, and appears to be the general pattern for the with-enzyme
  half of the dataset.

- **"no_E" experiments (65, 67, 69, 70 — standalone background-characterization
  runs, not part of a with/without-enzyme pair in the same sheet):** the
  reference cuvette is missing a *different* reagent, not the enzyme (there
  is no enzyme in either channel):
  - Exp 65 (Boric, pH 8.51): reference has `Sub [ml] = 0` (water makes up the
    volume), but **H2O2 is still present** (`H2O2 = 122.4 mM`, same as the
    sample). So the reference isolates "does H2O2 alone drift?" — it does
    **not** subtract the substrate+H2O2 background reaction, since the
    reference has no substrate for that reaction to run on.
  - Exp 67/69/70 (Phosphate, pH 8.01): reference has `H2O2 = 0`, but
    **substrate is still present** at the same concentration as the sample.
    So here the reference isolates "does substrate alone drift?" — it does
    **not** subtract the substrate+H2O2 background reaction either, since
    the reference has no H2O2 to drive it.

  Either way, the `.txt` curve for a no_E experiment is the **raw,
  undiminished non-enzymatic (substrate + H2O2) background reaction** —
  nothing chemically equivalent to it is subtracted. This background is
  expected to be substantial in its own right and likely autocatalytic
  (accelerating over the course of the reaction rather than a small linear
  drift), which further explains why its measured rate can look as large
  as, or larger than, the *net* (background-subtracted) rate reported for
  a with_E experiment.

**Practical implications for rate-expression fitting:**

- `[P]` for with_E samples throughout the main dataset should be treated as
  **already background-corrected** — do **not** additionally subtract a
  no_E curve from a with_E curve; that would double-subtract, since the
  with_E curve's own paired reference already performed that subtraction
  at the instrument level.
- The four no_E experiments (65, 67, 69, 70) are **not on the same
  reference basis** as the rest of the dataset (no matched enzyme-free
  reference of their own) and are not directly comparable, sample-for-sample,
  to any with_E experiment's raw rate. They exist to characterize the size
  of the background reaction on its own terms, not to be arithmetically
  combined with catalyzed runs.
- This was previously logged as an unresolved "open item" — it no longer
  is. No further action needed on this point.

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
- Metadata (pH, T, buffer, substrate, `[enz]`/`[buf]`/`[h2o2]`/`[sub]`), the
  raw `[P]` time-series data, and the cross-linkage between the two (sample
  alignment, range sanity, baseline, rate-vs-concentration monotonicity)
  have now all been audited (see entries above). Remaining known-unreliable
  field: the `.txt` files' own `Substrate Conc.` line — not used by the
  pipeline, don't trust it if it resurfaces.
- **Read this before fitting rate expressions:** `[P]` for with-enzyme
  samples throughout the dataset is a *differential* signal — each sample
  cuvette is measured against a paired, matched, no-enzyme reference
  cuvette (same `[sub]`/`[h2o2]`), so the non-enzymatic background is
  already subtracted at the source. Do not subtract a no-enzyme experiment
  (65, 67, 69, 70) from a with-enzyme one — that double-subtracts. See
  round 2 §4 for the full explanation and raw-cell evidence.
- Row-level block structure (which rows in each `.xls` the pipeline
  actually extracts vs. the paired reference rows it correctly leaves
  behind) has now been checked directly for all 98 experiments — see
  round 3. 97/98 clean; experiment 128 was the one exception, resolved in
  round 4 (samples 2 and 5 excluded, 3 and 4 kept as a reproducible
  near-zero-rate replicate pair).
- ~~Round 4 open question about other "backwards trend" exclusions~~ —
  **closed in round 5.** Checked every sample in 58/77/78/79/84 (plus 50/85
  as a sanity check) with a proper significance test; 128,3 was the only
  noise-scale case. One mischaracterization fixed: 84 was never actually a
  backwards-trend case — it's an incomplete run (1 of a planned 6-sample
  titration was ever recorded, and that one is a coarse quantized-staircase
  artifact) — but it's still correctly excluded either way. No code changes.
- `data/kinetics_io.py` and `data/plot_kinetics.py` now exist as reusable,
  single-source-of-truth tools for loading/plotting one experiment at a
  time (`python data/plot_kinetics.py <exp_num>`) — use these instead of
  writing new one-off extraction scripts for future investigations.
- **Buffer systems are not interchangeable** — borate, carbonate, phosphate
  and pyrophosphate each bring their own chemistry to an H2O2 reaction (see
  the 2026-08-30 mechanism entry, §1). Boric and Carbonate points in
  particular should be flagged as suspect in any pooled analysis. Rate
  constants from different buffers at the same pH are not directly
  comparable.
- **A buffer-concentration effect cannot be isolated from the current data**
  — `[buf]` and `[sub]` are perfectly collinear within every titration
  experiment (50 of 98 experiments vary `[buf]`; none hold `[sub]` fixed
  while doing so). Buffer *type* at matched pH is the comparison that is
  available; buffer *concentration* needs either a joint multi-variable
  regression or new bench work.
- `MECHANISM.md` now holds the proposed 7-step mechanism, its per-step
  literature support (51 references), the competing side reactions that
  belong in the ODE model as sink terms, and the open questions. Read it
  before building the kinetic model.
