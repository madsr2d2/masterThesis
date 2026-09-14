# Step 3, revised — the rate law first, then a module search

Handover plan. This version was written 2026-09-14 after reviewing R0.0 (commit
2a2dab4). It replaces every earlier version of this file (see git history) and
stages 3.2-3.3 of `PLAN_MECHANISM_NEXT_STEPS.md` section 4.

## Why this plan is written the way it is

Earlier rounds went wrong where the plan left a judgement to the executor: a
reading table turned into a mechanistic verdict, a test was made to pass by
starting at the truth, an F value was read without checking whether the two
fits were comparable, a constant's collapse was read without asking whether
another constant had absorbed it. So in this plan:

- **Verdicts are computed by functions, with fixed thresholds, and pinned by
  tests.** Reports quote the verdict strings those functions return.
- **Every task has the same parts:** Goal, Changes, Tests, Run, Report, Gate.
  Do them in order; do not skip a Report.
- **Anchors** are numbers you must reproduce before trusting a new function.
- **Assumptions** (section 4) each have a trigger. When a trigger fires, stop
  and report; do not adapt the plan.

---

## Contents

0. Rules for the executor
1. Where Step 3 stands
2. Glossary — what each parameter is
3. The curves, and what each block can identify
4. Assumptions and their stop triggers
5. Task list and stops
6. T1 — the reference-design parser, and exps 3 and 6
7. T2 — fit health, constant shifts and effective seed rates
8. T3 — the amplitude diagnosis on the R0.0 fits — STOP 1
9. T4 — model terms and fitter machinery
10. T5 — the search machinery, and a planted search — STOP 2
11. T6 — the background search
12. T7 — the catalysed rate law — STOP 3
13. T8 — the catalysed mechanism search
14. T9 — robustness, T10 — joint fit, T11 — acceptance — STOP 4
15. Report templates

---

## 0. Rules for the executor

**Repository rules** (`CLAUDE.md`, `PLAN_MECHANISM_NEXT_STEPS.md` section 0):
- `.venv/bin/python` from the repository root, never `python`.
- The duplicate guard covers every top-level name, private names and constants
  included. Check each new name with `grep -rn` before adding it.
- `.venv/bin/python run_gates.py` must print `0 failed` before every commit.
  `run_gates.py --all` (adds the ~9-minute optimiser suite) once per task that
  touches `fit_kinetics.py`, `kinetic_model.py` or `test_fit_kinetics.py`.
- Commit straight to master. Subject `area: what changed`, a paragraph, then
  exactly the two attribution lines from `PLAN_MECHANISM_NEXT_STEPS.md`
  section 5.
- A new `DATA_VERIFICATION.md` entry at the top for every task that produces a
  result, using the templates in section 15. End each with "Nothing is
  adopted."
- Never edit `FITTING.md`, `MECHANISM.md`, `CLAUDE.md`, `BUBBLES.md` or the
  skill. Propose text to the user in the stop message.
- Long fits run in the background, save as they go, and are resumable.
  Parallel workers set `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` and
  `MKL_NUM_THREADS=1` BEFORE importing numpy.
- Never overwrite an existing file in `data/fits/`.

**Rules that exist because they were broken before:**
1. **Tests.** Never loosen an assertion, threshold or tolerance in this plan to
   make a test pass. Never start a recovery fit at, or within 0.7 decades of,
   the planted truth. Never slice a planted design list; build the exact
   design the plan gives. If a test fails, stop and report the failing values.
2. **Comparisons.** Never compare two fits' costs, or read a constant's change
   between two fits, without `fit_health` on both and `constant_shift_report`
   on the pair (T2).
3. **Constants.** Never quote a constant as a value if `fit_health` flags it.
   Quote it as "not identified (reason)".
4. **Words.** Use the glossary names (section 2). Forbidden unless a function
   in this plan returns it as a verdict string:
   - that a term or module "earns", "is established", "is required", "is
     rejected" or "is the mechanism";
   - that a correction is "small", "negligible" or "does not matter";
   - that a result "rescues" or "confirms" a finding;
   - that the data "prefer" a species (H2O2 or HOO-, a buffer form).
5. **Numbers in documents** come from a named function in the repository,
   never from a throwaway script.
6. **Stops.** At a STOP, commit, post the stop message (section 15), and do
   nothing further until the user replies.

## 1. Where Step 3 stands

### 1.1 Done

- **8954989:** `kinetic_model.py`'s extension terms (`k_sink`, `K4`, `km_s`,
  `k_act_r`, `K_act`, `km_s_background`), each exactly OFF at its default.
- **da21765:** `fit_kinetics.ladder_f_test` and `ladder_checks`, and
  `data/test_fit_ladder.py`. The M0-M4b ladder was re-read; nothing in it is
  established.
- **2a2dab4 (R0.0):** `kinetic_model.increment`; `fit_dataset.Curve.reference_omits`;
  `fit_kinetics.observed_signal`, which fits a catalysed curve as the INCREMENT
  (full signal minus the same mixture at `e0 = 0`) where the reference omitted
  the enzyme; `_assert_observation_design`; `increment_comparison`; the saves
  `data/fits/<block>_R_M0.json`.

### 1.2 What the review of R0.0 found (the facts this plan builds on)

**F-A. The observable is now right.** Stage 1 of `_R_M0` matches the old M0 to
the last bit, and 31 gates pass.

**F-B. M0 cannot show whether the increment correction mattered.**
- Each block has ONE catalyst loading: 0.028 mM on BnOH; 0.241 mM on 4OMe
  (0.270 on exp 37). So M0's background seed `k0·[H2O2][S]` and catalysed seed
  `k5·E0·[H2O2][S]` have the same form under the catalysed curves.
- Removing the background from the observable therefore moves it into `k5`.
  BnOH `k5` went 2.24e-7 → 3.92e-7, against old `k5 + k0/E0` = 3.55e-7. On
  4OMe `k5` went 1.66e-8 → 4.24e-7, and `k6` (the catalysed peracid loop)
  7.68 → 4.5e-6.
- So the small cost change (−2.6% BnOH, −6.0% 4OMe) is uninformative. It must
  not be read as "the correction is small".
- Consequence for every later fit: **nothing here can test the order in E0.**
  `k5`, `k6` and `k_sink_e` are each identified only as that constant times
  the block's one E0.

**F-C. The remaining misfit is an amplitude error that tracks the conditions.**
Per run, the median over its curves of `net_data / net_model`, from the
`per_curve` table of `<block>_R_M0.json`. These were measured in the review
with a scratch script; T3's function must reproduce them (anchors):

| save, stage | run: median `net_data/net_model` |
|---|---|
| 4OMe `_R_M0`, stage 2 | 16: 0.658 · 32: 0.352 · 34: 0.402 · 35: 1.561 · 36: 0.933 · 37: 1.985 |
| BnOH `_R_M0`, stage 1 | 3: 0.567 · 6: 1.304 · 67: 14.808 · 69: 47.625 · 70: 64.174 |
| BnOH `_R_M0`, stage 2 | 68: 4.944 · 71: 8.398 · 73: 2.278 · 74: 7.275 · 83: 0.861 |

- **pH:** the pH 7.5 runs sit above the pH 7.0 runs on 4OMe, and the pH 8.0
  runs above the pH 6.7-7.5 runs on BnOH.
- **Buffer:** inside 4OMe catalysed runs 16, 32, 34 and 37, the ratio rises
  with [buf] (per-run correlation of log[buf] with log ratio 0.92-0.99).
- **Substrate:** the high-[S] runs (4OMe 36, BnOH 83) sit below their
  pH-mates.
- The reduced model's two seed steps have no [HOO-] and no [buf] dependence.
  The descriptive work measured both: a rate order of about +0.55 in [HOO-]
  (`ph/`), buffer orders (`buffer/`, Step 2) and substrate saturation
  (`saturation`). **Mechanism modules must not be searched until the rate law
  gets the amplitude right**; otherwise they absorb the amplitude error, as
  `km_s_background` did in 8954989.

**F-D. The two-axis block's "other" design is a parser bug.** On exps 135-151,
`verify_enzyme.reference_design` reads the `Sum:` and `Sum*9:` rows beneath the
cuvette table as cuvettes, because they carry a total volume. That shifts the
halves. The `Ref.` rows carry `Enz` 0.000, so the reference omits only the
enzyme. T1 fixes this. The block is NOT fitted in this plan.

## 2. Glossary — what each parameter is

Use these names and meanings. Do not rename. "OFF" is the value that removes
the term.

| parameter | meaning | fitted on | OFF |
|---|---|---|---|
| `k0` | uncatalysed seed: alcohol + peroxide → aldehyde | background | 0 |
| `k_can`, `k3` | steps 1-3: aldehyde → peracid (needs HOO-), peracid oxidises alcohol | background | 0 |
| `r` | ε(benzoate)/ε(aldehyde) at the assay wavelength; spectroscopy says 0.08-0.33 | background | — |
| `km_s_background` | the uncatalysed seed saturates in [S] | background | ∞ |
| `n_buf_background` | order of the uncatalysed seed in [buf] | background | 0 |
| `n_hoo_background` | order of the uncatalysed seed in [HOO-] | background | 0 |
| `k_sink` | first-order aldehyde loss in BOTH cuvettes (background sink) | control only | 0 |
| `k5` | catalysed seed (identified only as `k5·E0`, F-B) | catalysed | 0 |
| `n_hoo_catalysed` | order of the catalysed seed in [HOO-] | catalysed | 0 |
| `n_buf_catalysed` | order of the catalysed seed in [buf] | catalysed | 0 |
| `km_s` | the catalysed seed saturates in [S] | catalysed | ∞ |
| `k6` | steps 6-7: the catalysed PERACID LOOP. Not a sink. | catalysed | 0 |
| `k_sink_e` | catalyst-mediated aldehyde loss, `k_sink_e·E0·[A]` (the catalysed product sink) | catalysed | 0 |
| `k_act_r` | activation relaxation rate: catalyst activates as `1 − exp(−t·(k_act_r + k_act_f·[buf]))` | catalysed | ∞ |
| `k_act_f` | the activation rate's increase per mM buffer | catalysed | 0 |
| `K4`, `K_act` | held at 0 in this plan (not identified; see 3) | — | 0 |

Sign conventions: a clock τ = 1/k. An order of τ is MINUS the order of k.
State which one a number is.

## 3. The curves, and what each block can identify

Two separate global fits, 104 curves: every curve `build_curves` assigns to
the two groups. Reference designs are from `verify_enzyme`.

### BnOH, 25 °C, phosphate — 43 curves

| exp | kind | cuv | pH | [HOO-], mM | [S], mM | [H2O2], mM | [buf], mM | [enz], mM | run, s |
|---|---|---|---|---|---|---|---|---|---|
| 3 | background (design unruled) | 7 | 6.71 | — | 1.28–8.98 | 82.5 | 25–85 | 0 | 11025 |
| 6 | background (design unruled) | 4 | 6.71 | — | 3.67–9.16 | 165 | 40–70 | 0 | 17934 |
| 67 | background (omits H2O2) | 4 | 8.01 | — | 0.37–7.31 | 122.4 | 85 | 0 | 1178 |
| 69, 70 | background (omits H2O2) | 4 each | 8.01 | — | 0.21–2.11 | 122.4 | 85 | 0 | 1116–1767 |
| 68 | increment (omits enzyme) | 4 | 8.01 | 0.0405 | 0.36–7.31 | 122.4 | 85 | 0.028 | 1829 |
| 71 | increment | 4 | 8.01 | 0.0405 | 0.21–2.11 | 122.4 | 85 | 0.028 | 1225 |
| 73 | increment | 4 | 7.50 | 0.0122 | 0.32–6.33 | 122.4 | 75 | 0.028 | 2666 |
| 74 | increment | 4 | 8.00 | 0.0391 | 0.32–6.33 | 122.4 | 75 | 0.028 | 1984 |
| 83 | increment | 4 | 7.50 | 0.0122 | 6.36–14.31 | 122.4 | 75 | 0.028 | 3565 |

### 4-methoxybenzyl alcohol, 40 °C, phosphate — 61 curves

| exp | kind | cuv | pH | [HOO-], mM | [S], mM | [H2O2], mM | [buf], mM | [enz], mM | run, s |
|---|---|---|---|---|---|---|---|---|---|
| 23, 24, 28 | background (omits H2O2) | 4 each | 7.00 | — | 2.06–8.25 | 82.5 | 50–80 | 0 | 1764–17934 |
| 25 | background | 2 | 7.00 | — | 2.06, 6.19 | 82.5 | 60, 80 | 0 | 1078 |
| 26 | background | 4 identical | 7.00 | — | 4.13 | 82.5 | 75 | 0 | 392 |
| 27 | background | 4 | 7.00 | — | 1.03–4.13 | 82.5 | 75–90 | 0 | 1764 |
| 29 | background | 4 | 7.00 | — | 0.38–1.53 | 82.5 | 75–90 | 0 | 1008 |
| 30, 39 | background | 4 each | 7.00, 6.97 | — | 0.10–0.38 | 82.5 | 75–90 | 0 | 729–853 |
| 38 | background | 3 | 7.00 | — | 0.19–0.57 | 82.5 | 65–85 | 0 | 7161 |
| 16 | increment (omits enzyme) | 4 | 7.00 | 0.0025 | 1.85–7.40 | 82.5 | 50–80 | 0.241 | 6542 |
| 32 | increment | 4 | 7.00 | 0.0024 | 8.25 | 82.5 | 50–200 | 0.241 | 1767 |
| 34 | increment | 4 | 7.00 | 0.0021 | 8.25 | 82.5 | 3.1–25 | 0.241 | 5280 |
| 35 | increment | 4 | 7.50 | 0.0078 | 8.25 | 82.5 | 50–200 | 0.241 | 7920 |
| 36 | increment | 4 | 7.53 | 0.0083 | 57.9 | 82.5 | 50–200 | 0.241 | 2100 |
| 37 | increment | 4 | 7.53 | 0.0083 | 12.23 | 82.5 | 50–200 | 0.270 | 3210 |

The filenames of exps 32 and 34-37 say `with_NO_E`. They are catalysed
(`build_manifest.RULINGS`). Trust `[enz]`, never a filename.

### Which parameters each block can identify, and from what contrast

The search library (T5) enforces this table in code. Do not add a module to a
block this table excludes.

| parameter | BnOH | 4OMe | contrast / reason |
|---|---|---|---|
| `km_s_background` | yes | yes | [S] within background runs |
| `n_buf_background` | weak | confounded with [S] | BnOH: within exps 3 and 6 only. 4OMe: [buf] moves against [S] in all 9 ladders (median r = −0.974); `km_s_background` and `n_buf_background` may not separate — the search reports whether they do |
| `n_hoo_background` | confounded | **excluded** | BnOH: pH 6.71 vs 8.01 runs also differ in [H2O2], buffer range and run length (11-18 ks vs 1.1-1.8 ks). 4OMe: one pH |
| `n_hoo_catalysed` | yes, between runs | yes, between runs | BnOH 0.0122 vs 0.039-0.041 mM; 4OMe 0.0021-0.0025 vs 0.0078-0.0083 mM |
| `n_buf_catalysed` | **excluded** | yes, within runs | BnOH catalysed [buf] is 75 or 85 mM, one value per run |
| `km_s` | yes, within runs | yes, between runs | 4OMe: exps 36 vs 37 share pH, [buf] and [HOO-] at [S] 57.9 vs 12.2 (and exp 16's ladder) |
| `k6` (loop) | yes | yes | needs the background's `k_can`, `k3` on — no peracid otherwise |
| `k_sink_e` | yes | yes | late-run shape; needs long runs |
| `k_act_r` | yes | yes | early-run shape |
| `k_act_f` | **excluded** | yes, within runs | BnOH has no within-run [buf] |
| `K4` | **excluded** | **excluded** | one [H2O2] per catalysed block |

**Within-run [HOO-] is not an axis.** [HOO-] shifts slightly with [buf] inside
a run (ionic strength), so a within-run regression on [HOO-] would read the
buffer. [HOO-] is used only between runs.

## 4. Assumptions and their stop triggers

| id | assumption | trigger → STOP and report |
|---|---|---|
| A1 | every catalysed curve's reference omits the enzyme; every background curve's omits the H2O2 (exps 3, 6 pending T1) | T1 classifies any curve in either block otherwise |
| A2 | T1's parser fix changes classifications ONLY for exps 135-151 | any other experiment's design changes |
| A3 | the F-C anchors reproduce | T3's per-run medians differ from the anchors by more than 1e-3 relative |
| A4 | the amplitude error tracks conditions (F-C) | T3's `amplitude_verdict` on 4OMe stage 2 does not flag within-run `buf`, or on BnOH stage 1 does not flag between-run `hoo` |
| A5 | a correct rate law removes most of the amplitude error | T7's best rate-law model still gets any axis flagged by `amplitude_verdict` |
| A6 | the fitter can recover what the search reports | T4's recovery test or T5's planted search fails |
| A7 | a larger model never fits worse than a model it contains | `nesting_violation` still true after the one allowed refit |
| A8 | the search is affordable | T5's pilot projects more than 48 h wall time, even for the greedy fallback |

## 5. Task list and stops

| task | content | ends with |
|---|---|---|
| T1 | parser fix; exps 3 and 6 ruling | commit |
| T2 | `fit_health`, `constant_shift_report`, `effective_seed_rates` | commit |
| T3 | `amplitude_by_condition`, `amplitude_verdict`, `residual_thirds` on the `_R_M0` saves | **STOP 1** |
| T4 | new model terms; warm starts; frozen background; bound overrides; recovery test | commit |
| T5 | `data/mechanism_search.py`; planted search; pilot | **STOP 2** |
| T6 | background search | commit |
| T7 | catalysed rate-law search; amplitude re-check; matched compositions; background leverage | **STOP 3** |
| T8 | catalysed mechanism search on the rate-law base | commit |
| T9 | robustness cuts | commit |
| T10 | joint fit | commit |
| T11 | `acceptance_report` | **STOP 4** |

---

## 6. T1 — the reference-design parser, and exps 3 and 6

**Goal.** Make `verify_enzyme.reference_design` read the two-axis block's
sheets correctly (F-D), and rule exps 3 and 6's designs.

**Changes.**
1. In `reference_design`, stop collecting cuvette rows at the first row whose
   first-column label, lower-cased and stripped, starts with `sum`. Do not
   change anything else in the function.
2. Exps 3 and 6: open each sheet (`manifest` → `xls_file`, `Sheet1`) and write
   down every cuvette row's `Enz`, `H2O2`, `Sub`, `Buf` and total volume. Say
   which rows are measured and which are references, and what the references
   omit.
   - If the layout is a readable variant, extend `reference_design` to read it.
   - If it is not, add a `RULINGS` entry in `data/build_manifest.py` stating
     the design and the evidence, and make
     `verify_enzyme.reference_omits_by_experiment` use it.
   - Remove the ruled numbers from `fit_dataset.REFERENCE_OMITS_UNRULED`.

**Tests** (`data/test_fit_ladder.py`):
- `test_the_two_axis_sheets_read_as_enzyme_references`: exps 135, 140 and 151
  classify as "enzyme".
- `test_the_sum_rows_do_not_change_other_designs`: over every experiment in
  the manifest, the classification before and after the change differs only
  for experiments in `fit_dataset.TWO_AXIS_BLOCK`. Compute "before" by calling
  the function with the new stop disabled (add a keyword
  `stop_at_sum=True`; the test passes `False` for "before").
- `test_every_fitted_curve_has_the_design_its_observable_assumes`: remove exps
  3 and 6 from its exception list once ruled.

**Run.** `.venv/bin/python data/verify_enzyme.py`; record the design counts.

**Report** (template T1). The design counts before and after; the exps 3 and 6
tables; the ruling.

**Gate.** A1 or A2 fires → STOP. If exp 3 or 6 is not a raw background, drop
it from stage 1 by a `RULINGS` entry, and say so. Otherwise commit and go on.

## 7. T2 — fit health, constant shifts and effective seed rates

**Goal.** Put the checks that were skipped before into code, so that no later
report can omit them.

**Changes** (`data/fit_kinetics.py`):

1. `fit_health(save_stage, sub_stages=(), recovered=None) -> dict` on one
   saved stage dict. Returns the list `flags` of strings, and `ok`
   (`flags == []`). Flags, exactly:
   - `"not converged"` if `converged` is false.
   - `"no correlation matrix"` if `correlation` is None.
   - `"correlated: a/b r=±0.xxx"` for every free pair with |r| > 0.99.
   - `"at bound: name side"` for every `at_bound` entry whose side is not the
     parameter's OFF end (OFF ends from the glossary: lower for parameters OFF
     at 0 in log10; upper for parameters OFF at ∞; for linear orders OFF is
     0, which is not a bound, so any bound hit is flagged).
   - `"nested worse than sub-model M"` for any stage in `sub_stages` (same
     observation, same frozen background) whose cost is below this one's by
     more than `NESTING_TOLERANCE`.
   - `"not recovered: name"` for every free parameter absent from `recovered`
     when `recovered` (a set of names) is given.
   And `unidentified`: the set of parameter names named in a `correlated`,
   `at bound` or `not recovered` flag.
2. `effective_seed_rates(save, curves) -> DataFrame`. One row per catalysed run
   and per background run: the run's conditions, and at the save's constants
   the initial seed rates the model uses, in mM/s at t = 0:
   `background_seed = k0 · (all background factors) · [H2O2] · [S]` and
   `catalysed_seed = k5 · E0 · (all catalysed factors) · [H2O2] · [S]`, with
   the median over the run's curves. This is what makes F-B visible: two fits
   whose constants differ can have the same effective rates.
3. `constant_shift_report(save_a, save_b, curves) -> dict` with:
   - `constants`: per parameter, a, b and the shift in log10 (or linear for
     `r` and orders), with each fit's `fit_health` flags;
   - `effective`: `effective_seed_rates` for both, merged, with the log10
     ratio `catalysed_seed_b / catalysed_seed_a` per run;
   - `verdict`, exactly one of:
     - `"constants moved, effective rates did not"`: some free constant moved
       by more than 0.3 decades AND every run's |log10 effective ratio| < 0.1;
     - `"effective rates moved"`: any run's |log10 effective ratio| ≥ 0.1;
     - `"no material change"`: neither.

**Tests** (`data/test_fit_ladder.py`, fast, on hand-built stage dicts and a
3-curve planted group):
- `test_fit_health_flags`: one planted dict per flag; each produces exactly
  that flag.
- `test_a_constant_swap_is_named`: two saves with (`k0` = 1e-8, `k5` = 1e-7)
  and (`k0` = 0, `k5` = 1e-7 + 1e-8/E0), same everything else → verdict
  "constants moved, effective rates did not".
- `test_effective_seed_rates_anchor`: on the real `BnOH_25C_Phosphate_M0.json`
  and `_R_M0.json`, BnOH run 68's catalysed seed at old constants plus its
  background seed equals, within 15%, its catalysed seed at the new constants.
  (That is F-B in numbers. If it fails, stop: F-B is wrong.)

**Run.** `constant_shift_report` on `<block>_M0.json` against
`<block>_R_M0.json`, both blocks.

**Report** (template T2). Both verdicts, both effective-rate tables.

**Gate.** Commit and go on.

## 8. T3 — the amplitude diagnosis on the R0.0 fits — STOP 1

**Goal.** Measure, by function, whether the misfit is an amplitude error, and
whether it tracks pH, buffer or substrate.

**Changes** (`data/fit_kinetics.py`):

1. `amplitude_by_condition(block, save, stage, curves=None) -> dict`:
   - `per_curve`: the save's `per_curve` rows merged with each curve's
     `pH`, `hoo`, `buf`, `s0`, and `ratio = net_data / net_model`. Rows with
     `net_model <= 0` or `net_data <= 0` are dropped and COUNTED
     (`dropped`).
   - `per_run`: per experiment, the median ratio, its min and max, and the
     run's share of the stage's summed `rms²`.
   - `between`: for each axis in (`hoo`, `buf`, `s0`), an OLS of the per-run
     log median ratio on the per-run log median of the axis: slope, stderr,
     t = slope/stderr, number of runs. Skip an axis that is constant across
     runs. Also the pairwise correlation of the per-run log axes.
   - `within`: `scope.orders("ratio", frame=per_curve, terms=("buf", "s0"),
     within=True, live_only=False)` — orders, stderrs and t. NEVER include
     `hoo` here (section 3). Also `induction.composition_collinearity` on the
     same rows (it filters on a `live` column: add `live = True` to
     `per_curve` first).
2. `amplitude_verdict(result) -> dict` with `flags` (a list) and `verdict` (a
   string). Thresholds, fixed:
   - Between-run axis flagged if |t| > 3 with at least 4 runs.
   - Within-run axis flagged if |t| > 3.
   - If two flagged between-run axes correlate at |r| > 0.8 across runs, they
     are reported as ONE flag `"between: a+b (confounded)"`.
   - If `buf` and `s0` are both flagged within runs and
     `composition_collinearity`'s |median| > 0.9, they are reported as ONE flag
     `"within: buf+s0 (confounded)"`.
   - BnOH stage 1 `hoo` is always suffixed `"(confounded with [H2O2], buffer
     range, run length)"`.
   - `verdict`:
     - `"tracks conditions: <flags>"` if any flag;
     - else `"uniform amplitude offset"` if the median of per-run medians is
       outside [0.8, 1.25] and the SD of per-run log medians is below 0.2;
     - else `"concentrated in runs: <runs>"` if the two largest `rms²`
       shares sum above 0.6;
     - else `"no amplitude error detected"`.
3. `residual_thirds(block, save, stage, curves=None) -> DataFrame`: per curve,
   the mean residual (data − model, the save's observation) in the first,
   middle and last third of the run, in AU and in units of the curve's noise.

**Tests** (`data/test_fit_ladder.py`, fast):
- `test_amplitude_anchors`: the per-run medians of section 1.2 F-C, each to
  1e-3 relative, on the real `_R_M0` saves.
- `test_amplitude_verdict_rules`: planted per-curve tables:
  - ratio ∝ [buf]^0.5 within runs → flags within `buf`;
  - ratio ∝ [HOO-]^1 between 5 runs → flags between `hoo`;
  - ratio constant 3.0 → "uniform amplitude offset";
  - one run at ratio 10, others at 1 → "concentrated in runs: <that run>";
  - [buf] and [S] collinear at r = −0.97, both driving → the single
    confounded flag.
- `test_residual_thirds_are_zero_on_a_planted_fit`.

**Run.** Both functions and `residual_thirds` on both `_R_M0` saves, both
stages.

**Report** (template T3). Four `amplitude_verdict` results; four per-run tables;
the `between` and `within` tables; the ten curves with the largest rms and
their thirds.

**Gate.** A3 or A4 fires → STOP and report (do not change thresholds). Else
**STOP 1** (template STOP 1).

## 9. T4 — model terms and fitter machinery

**Goal.** Everything the searches need, each piece tested. No real-data fits in
this task except the recovery test.

**Changes.**

1. **`kinetic_model.RateConstants` gains**, each OFF at its default:
   `k_sink_e = 0.0`, `k_act_f = 0.0`, `n_buf_background = 0.0`,
   `n_hoo_background = 0.0`, `n_hoo_catalysed = 0.0`, `n_buf_catalysed = 0.0`.
   Module constants `BUF_REFERENCE = 75.0` (mM) and `HOO_REFERENCE`: the
   geometric mean of [HOO-] over the catalysed curves of both blocks, computed
   once and written as a literal with a comment saying how.
2. **`rates`** becomes, with every new factor exactly 1 at its default:
   ```
   bg_factor  = sat_background · (buf/BUF_REFERENCE)^n_buf_background · (hoo/HOO_REFERENCE)^n_hoo_background
   cat_factor = phi · free · sat_s · (buf/BUF_REFERENCE)^n_buf_catalysed · (hoo/HOO_REFERENCE)^n_hoo_catalysed
   v_seed     = (k0 · bg_factor + k5 · E0 · cat_factor) · [H2O2] · [S]
   v6         = k6 · E0 · phi · free · [PBA]
   v_sink     = (k_sink + k_sink_e · E0) · [A]
   phi        = 1 − exp(−t · (k_act_r · (1 + K_act·[buf]) + k_act_f·[buf]))   (1 when k_act_r = ∞)
   ```
   A power factor is 1 when its order is 0 or its concentration is NaN. The
   catalysed orders multiply the SEED only, not `k6`: the descriptive orders
   were measured on the initial rate, which the seed sets.
3. `PARAMETER_NAMES`, `BOUNDS`, `EXTENDED_INITIAL`, `UNITS` gain all six. The
   four orders are LINEAR (not in `LOG_PARAMETERS`), bounds (−2, 2), start 0.
   `k_sink_e` log10 bounds (−10, 2), start −3. `k_act_f` log10 bounds
   (−10, 0), start −5.
4. **Warm starts.** `fit_group(..., extra_starts=(), bounds=None)`:
   `extra_starts` are full `RateConstants` whose free parameters are always
   run as starts; `bounds` overrides `BOUNDS` for this fit only.
   `nested_start(sub_constants, free_names)`: the sub-model's constants with
   every parameter it does not free set just inside its OFF end, such that the
   curve differs from the sub-model's by less than 1e-6 AU (verify by
   simulation; step further in if not).
   `nesting_violation(cost, sub_costs)`: cost exceeds min(sub_costs) by more
   than `NESTING_TOLERANCE` (relative).
5. **Frozen background.** `sequential_fit(..., background=None)`: a
   `RateConstants`; when given, stage 1 is skipped and those constants are
   used. Saves record `"background"` (an id string) and `"observation"`.

**Tests.**

`data/test_kinetic_model.py`:
- the existing `FROZEN` arrays still reproduce to 1e-10 at defaults;
- `test_each_new_order_is_its_log_derivative`: for each of the four orders set
  to 0.7 alone, the numerical d ln(initial seed)/d ln(concentration) equals 0.7
  to 1%;
- `test_the_catalysed_sink_leaves_the_background_alone`: with `k_sink_e` on,
  `observable` at `e0 = 0` is unchanged to 1e-12, and an aldehyde pulse at
  `e0 = 0.1` with every seed off decays as `exp(−k_sink_e·0.1·t)` to 1e-6
  relative;
- `test_forward_activation_is_the_same_relaxation`: `(k_act_r=k, K_act=K)` and
  `(k_act_r=k, K_act=0, k_act_f=k·K)` give the same observable to 1e-10.

`data/test_fit_ladder.py` (fast):
- `test_a_nested_start_reproduces_the_sub_model` (1e-6 AU);
- `test_a_warm_started_larger_model_never_fits_worse` (3 planted curves, 20
  points, `restarts=0`);
- `test_bound_overrides_do_not_leak` (`BOUNDS` unchanged after the fit).

`data/test_fit_kinetics.py` — **replace** `test_extended_parameter_recovery`
(keep the name). Build it exactly like this:
- **Curves:** every curve of the 4OMe/40 °C/phosphate group from
  `build_curves()`, with its own `times`, `conditions`, `noise` and
  `reference_omits`. Do not subset.
- **Truth:** `k_can = 5.0`, `k3 = 2e-2`, `k0 = 2e-8`, `r = 0.2`,
  `n_buf_background = 0.8`; `k5 = 4e-7`, `n_hoo_catalysed = 0.6`,
  `n_buf_catalysed = 0.4`, `km_s = 5.0`, `k_sink_e = 5e-3`, `k_act_r = 1e-3`,
  `k_act_f = 1e-5`; `k6 = 0`. Assert, inside the test, that at least 5 of the
  24 catalysed curves have `peak_model > LAG_PEAK_FRACTION`. If not, stop:
  do not change the truth to make it pass.
- **Data:** `baseline_like_data(epsilon · observed_signal)` plus Gaussian noise
  at each curve's `noise`, seeds 0, 1, 2; plus one noiseless copy.
- **Starts:** a start table built in the test, each log10 parameter at truth
  +0.8 decades for even positions in `PARAMETER_NAMES` and −0.8 for odd ones,
  orders at 0, `r` at 0.5. Assert every start is ≥ 0.7 decades from truth.
  `restarts=4`. Stage 1 is fitted, then stage 2 on the fitted stage 1.
- **Assertions:** noiseless, every free parameter within 0.1 decade (orders
  and `r` within 0.05). Noisy, a parameter is "recovered" if all three seeds put
  it within 0.3 decade (orders and `r` within 0.15). Write the recovered set to
  `data/fits/search/recovered.json`. The test passes if the noiseless check
  passes, and fails listing the unrecovered parameters otherwise; the noisy
  result is recorded, not asserted, and T5 uses it.
- If it runs over 20 minutes, report the time; do not cut the design.

**Run.** `run_gates.py --all`.

**Report** (template T4). The recovered set; the recovery test's time.

**Gate.** A6 fires (noiseless recovery fails) → STOP. Otherwise commit. The
parameters NOT in `recovered.json` are reported as "not identified" wherever
they appear later.

## 10. T5 — the search machinery, and a planted search — STOP 2

**Goal.** `data/mechanism_search.py` and its gate `data/test_mechanism_search.py`.
Check that the module name and every top-level name in it are unused.

### T5.1 The library, in code

Each module is a record: `name`, `phase`, `frees` (parameter names), `fixes`
(dict), `bounds` (dict), `requires` (module names), `blocks` (which blocks).
A model is a frozenset of module names. `enumerate_models(block, phase,
base=frozenset(), candidates=None)` returns every model that respects
`requires` and `blocks`, and always contains `base`.

**Background phase** (base frees `k0`, `r`; fixes `k_can = k3 = 0`):

| module | frees / sets | blocks |
|---|---|---|
| `bg_can` | frees `k_can`, `k3` | both |
| `bg_sat` | frees `km_s_background` | both |
| `bg_buf` | frees `n_buf_background` | both |
| `bg_hoo` | frees `n_hoo_background` | BnOH only |
| `bg_rspec` | bounds `r` to (0.08, 0.33) | both |

**Catalysed rate-law tier** (base frees `k5`; fixes `k6 = 0`, `k_act_r = ∞`,
`k_sink_e = 0`):

| module | frees | blocks |
|---|---|---|
| `rl_hoo` | `n_hoo_catalysed` | both |
| `rl_buf` | `n_buf_catalysed` | 4OMe only |
| `rl_sat` | `km_s` | both |

**Catalysed mechanism tier** (base = a rate-law model from T7):

| module | frees | requires | blocks |
|---|---|---|---|
| `mech_loop` | `k6` | the background model contains `bg_can` | both |
| `mech_sink` | `k_sink_e` | — | both |
| `mech_act` | `k_act_r` | — | both |
| `mech_actbuf` | `k_act_f` | `mech_act` | 4OMe only |

Counts the test must find: background 32 (BnOH) and 16 (4OMe); rate-law 4
(BnOH) and 8 (4OMe); mechanism, on a background with `bg_can`, 16 (BnOH) and
12 (4OMe).

### T5.2 Fitting and leave-one-run-out

- `fit_model(block, phase, model, curves, background=None, holdout=None,
  sub_saves=())`: frees the base's and every module's parameters; applies
  `fixes` and `bounds`; `extra_starts` = `nested_start` from each sub-model
  save in `sub_saves`; `observation="design"`; returns the save dict, which
  also holds `modules`, `phase`, `background`, `holdout` and `fit_health`.
- Full-data fits go in order of model size. After each, `nesting_violation`
  against every fitted sub-model; if true, refit once with `restarts` doubled;
  if still true, keep the flag (A7 → STOP after the phase finishes).
- `cross_validate(...)`: for each run of the phase (background runs in the
  background phase; catalysed runs otherwise), refit on the other runs
  warm-started from the full-data fit (`restarts=1`), predict the held-out run
  with its own conditions and observation, and score it: the sum over its
  curves of the fitter's weighted cost (`weighting="curve"`). CV score = sum
  over folds. A failed integration scores `inf`; the model is listed, not
  ranked.
- In catalysed phases the background stays frozen at the background model's
  full-data fit.
- Saves: `data/fits/search/<block>/<model_id>.json` and
  `<model_id>__without_<exp>.json`. `model_id` is `phase:` plus the sorted
  module names joined by `+`, plus `@` and the background id for catalysed
  phases. Skip any fit whose save exists.
- `search(block, phase, ..., workers=8)`: runs everything with
  `ProcessPoolExecutor`.

### T5.3 Ranking and verdicts, in code

- `tie_set(table)`: model `m` is tied with the best if its per-fold
  differences `d` satisfy `mean(d) <= 2 · sd(d) / sqrt(folds)`.
- `module_verdicts(table)`: for each candidate module, among the tie set:
  `"in every tied model"`, `"in no tied model"`, or `"in some tied models"`.
  These three strings are the only module verdicts reports may use.
- `promote(verdicts)`: modules "in every tied model" become part of the next
  tier's base; "in some tied models" stay candidates; "in no tied model" are
  dropped.
- A tie-set model with `fit_health` flags stays in the tie set; its flagged
  constants are reported as "not identified".

### T5.4 Tests

`data/test_mechanism_search.py` (fast):
- `test_the_library_counts_and_constraints` (the counts in T5.1; no
  `mech_loop` without `bg_can`; no `mech_actbuf` without `mech_act`; no BnOH
  `bg_hoo`-excluded modules on 4OMe and vice versa);
- `test_model_ids_round_trip`;
- `test_every_run_is_held_out_once_and_never_trained_on`;
- `test_the_tie_rule` (hand-built fold scores, one tied, one not);
- `test_module_verdicts_and_promotion` (hand-built tie set);
- `test_a_finished_fit_is_not_refitted` (monkeypatched fitter never called).

`data/test_fit_kinetics.py` (slow) — `test_the_search_finds_a_planted_mechanism`,
exactly:
- Curves: the BnOH/25 °C/phosphate group from `build_curves()`, all of them.
- Truth: background M0 constants from `BnOH_25C_Phosphate_R_M0.json` stage 1
  (`bg_can` on, nothing else); catalysed `k5` = 4e-7, `n_hoo_catalysed` = 0.6,
  `k_sink_e` = 2e-2, `k_act_r` = 2e-3; `k6 = 0`, `km_s = ∞`.
- Data: the observed signal plus noise at each curve's `noise`, seed 0.
- Run the rate-law tier on the true background, promote, then the mechanism
  tier.
- Assert: `rl_hoo` "in every tied model"; `rl_sat` not "in every tied model";
  `mech_sink` and `mech_act` "in every tied model"; `mech_loop` not "in every
  tied model".

### T5.5 Pilot

Time one full-data fit and one fold refit of the LARGEST model of each phase on
each block. Project the wall time: sum over phases of
`models × (folds + 1) × fit time / workers`. Use these model counts: background
32 + 16; rate-law 4 × 3 + 8 × 3 (up to 3 backgrounds each); mechanism 16 × 3 +
12 × 3 (up to 3 rate-law bases each).

If a phase projects more than 48 h, it uses the greedy search: from the base,
add the candidate module that lowers the CV score most, while the result is not
tied with the current best; then remove modules while the result stays tied.
Run a second greedy pass from the full model downwards. Report both end points.

**Report** (template T5). The planted search's module verdicts; the pilot
table; full or greedy per phase.

**Gate.** A6 (planted search) or A8 fires → STOP. Else **STOP 2**.

## 11. T6 — the background search

**Run.** `search` on the background phase for both blocks. Also, for BnOH,
repeat the whole phase without exps 3 and 6 if T1 kept them.

**Report** (template T6). The CV ranking with per-fold scores; the tie set;
`module_verdicts`; `fit_health` flags of tied models; `amplitude_verdict` on the
best tied model's stage 1; for 4OMe, whether `bg_sat` and `bg_buf` are each "in
every", "in no" or "in some" tied models.

**Backgrounds carried forward.** The tie set, capped at the 3 best. On 4OMe,
if the best model containing `bg_sat` or the best containing `bg_buf` is
outside the cap, add it, whatever its rank. The catalysed runs reach [buf]
3.1-200 mM, and a background measured at 50-90 mM is extrapolated there.

**Gate.** A7 → STOP. Else commit and go on.

## 12. T7 — the catalysed rate law — STOP 3

**Run.** The rate-law tier on each carried background, both blocks.

**Then, on the best tied (background, rate-law) model of each block:**
1. `amplitude_verdict` on stage 2.
2. `constant_shift_report` against `<block>_R_M0.json`.
3. `matched_compositions(block, save)` (BnOH only), a new function:
   - For each matched cuvette pair — exp 67 with exp 68; exps 69 and 70 with
     exp 71 — interpolate both measured curves onto the shorter run's times.
     Form `absolute_data = background + increment`.
   - Set it beside `observable(E0)` and `increment(E0)` at the save's
     constants. Report net rises at the common end time and the rms of each
     difference.
   - Also report the exp 69 against exp 70 difference, the between-run
     spread at one composition.
   - Test it on a planted pair made from one set of constants
     (`absolute_data == absolute_model` to 1e-10).
4. `background_leverage(block, save)`, a new function:
   - Per catalysed curve, the net increment at the save's constants, and the
     same with `k0 = k_can = k3 = 0`.
   - The share of the increment that difference makes up; flag curves above
     0.10.
   - Report the block median share.
   - Test on a planted curve with `k0 = k_can = k3 = 0` (share 0 to 1e-12).

**Report** (template T7).

**Gate.** A5 (any axis still flagged) → STOP and report which axis. Otherwise
**STOP 3**, and do not start T8 until the user replies.

## 13. T8 — the catalysed mechanism search

**Run.** For each carried background, `promote` the rate-law verdicts into the
base. Then run the mechanism tier on up to 3 (background, rate-law base)
combinations, chosen by CV score. Candidates are the mechanism modules plus
any rate-law modules "in some tied models".

**Report** (template T8). Ranking with per-fold scores; the tie set;
`module_verdicts`; `fit_health` flags; `ladder_checks` (substrate order, lag
counts) and `amplitude_verdict` on every tied model.

**Gate.** A7 → STOP. Else commit and go on.

## 14. T9 — robustness, T10 — joint fit, T11 — acceptance — STOP 4

### T9 — robustness

On every model in each block's final tie set (capped at 3), rerun the
catalysed phase (or, for S4, the background phase):

| cut | drops | reason |
|---|---|---|
| S1 | 4OMe exp 36; BnOH exp 83 | [S] beyond the background's measured range |
| S2 | 4OMe exps 35, 36, 37; BnOH exps 73, 83 | pH outside the background's measured pH |
| S3 | 4OMe exps 32, 34 | [buf] outside 50-90 mM |
| S4 | background without BnOH exps 3 and 6 (if T1 kept them) | only pH 6.71 background |
| B-sink | each tied background plus free `k_sink` | a background sink that lowers the CV score contradicts `product_fate`'s enzyme-free clock |

Report every constant's shift in standard errors, and every module verdict that
changes. S3 removes most of 4OMe's buffer axis, so `rl_buf` and `mech_actbuf`
moving to "in some tied models" there is expected.

### T10 — joint fit

For each final tied model (capped at 3), fit all curves of the block at once,
all parameters free, each curve with its design's observation, with
`extra_starts` from the sequential fit. Report rms per kind of curve,
`fit_health`, and `constant_shift_report` against the sequential fit.

### T11 — acceptance

`acceptance_report(block, tie_saves) -> DataFrame` with columns `item`,
`value`, `threshold`, `status` ("met" / "not met" / "not identified"). One row
per item, per tied model:

1. every catalysed curve fitted as an increment;
2. the mechanism tie set holds fewer than half the models searched;
3. `amplitude_verdict` is not "tracks conditions";
4. substrate order within 2 sigma of the data's (`ladder_checks`), per kind of
   curve;
5. lag counts within 3 curves of the data's, per kind of curve;
6. activation clock: over catalysed curves whose `act_sink_kind_corrected` in
   `scope.frame(tuple(experiments))` is "lag" or "lag then sink", the median of
   `(1/(k_act_r + k_act_f·[buf])) / tau_act_corrected` lies within [0.5, 2];
   runs whose curves mix kinds are excluded and counted;
7. `r` within 0.08-0.33 at T10;
8. if `mech_sink` is in the model, `k_sink_e·E0` against
   `slowdown.sink_constants` for 4OMe at 40 °C, within its error;
9. matched compositions (BnOH): the model's absolute signal within twice the
   exp 69 vs exp 70 spread;
10. no reported constant shifts beyond 2 standard errors under S1-S4.

**Report** (template T11). **STOP 4.**

## 15. Report templates

Fill the brackets. Tables come from the named functions. Add no verdict
sentences beyond the ones shown.

### DATA_VERIFICATION.md entry (every task that produces a result)

```
## <date> — Step 3 <task>: <the verdict string(s), verbatim>

What was asked: <one sentence from this plan's Goal>.
What was built: <function names>.
Bugs found: <before/after numbers, or "none">.

<the task's tables>

Verdicts (from <function names>): <verdict strings, verbatim>.
Health: <every fit_health flag of every fit quoted, or "no flags">.
Not identified: <parameter names and reasons, or "none">.
Could overturn this: <the task's list below>.

Nothing is adopted.
```

"Could overturn this" lists:
- T1: a sheet layout the parser still misreads.
- T2: none.
- T3: the per-curve ratio is a net-rise proxy, not an initial rate; the
  between-run tests have 5-6 runs.
- T4: the recovery design is one block.
- T5: the planted search is one block and one truth.
- T6: 5 (BnOH) and 10 (4OMe) folds; the [S]/[buf] confound on 4OMe.
- T7: E0 is one value per block (F-B); the matched compositions are
  different days.
- T8-T11: all of the above; the tie rule is generous at 5-6 folds.

### Stop messages (post in chat, nothing else)

```
STOP <n> — <task>

Done:
- <at most five bullets, each naming a function or file>

<the task's tables>

Verdicts: <verbatim strings>
Health flags: <list or "none">
Could overturn this: <the task's list>

Question for you: <the stop's question below>
```

Questions:
- **STOP 1:** "The amplitude verdicts are above. Proceed to T4-T5 (new terms,
  fitter machinery, search machinery and a pilot)?"
- **STOP 2:** "The planted search and pilot are above; projected wall time
  <h>. Run the background search and the catalysed rate law (T6-T7)?"
- **STOP 3:** "The catalysed rate law is above. Run the mechanism search
  (T8-T11)? Separately: FITTING.md F5/F6 were measured on the absolute
  observable. Proposed note: <text>. Add it?"
- **STOP 4:** "Acceptance is above. Proposed FITTING.md / MECHANISM.md text:
  <text>. Adopt any of it?"
