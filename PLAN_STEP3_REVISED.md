# Step 3, revised — a ladder that can answer the question

Handover plan, written 2026-09-14 after the review in `DATA_VERIFICATION.md`
2026-09-14 (second entry), and extended the same day once the curves in each
fit were listed. It REPLACES stages 3.2-3.3 and the acceptance list of
`PLAN_MECHANISM_NEXT_STEPS.md` section 4. Stages 3.0-3.1 of that plan are done
(commit 8954989) and are not repeated.

**Read section 1.2 before anything else.** The fitter has compared the
catalysed curves with the wrong quantity since it was written, and every
stage-2 result — including FITTING.md F5 and F6 — was measured that way.

## Contents

0. Read first
1. Where Step 3 stands
2. The curves in each fit, and what the design limits
3. What Steps 1 and 2 can and cannot feed in
4. Stage R0 — fix the machinery before fitting anything
5. Stage R1 — look before adding terms
6. Stage R2 — the background, with its confound
7. Stage R3 — the catalysed ladder on one background
8. Stage R4 — the joint fit
9. Acceptance
10. Stop conditions
11. Deliverables, reporting, committing

---

## 0. Read first

- **Every rule in `PLAN_MECHANISM_NEXT_STEPS.md` section 0 still applies**:
  `.venv/bin/python` from the repository root, never `python`; the duplicate
  guard covers private names and constants; `run_gates.py` green before every
  commit; commit straight to master with the two attribution lines; new
  `DATA_VERIFICATION.md` entries at the top; ASK the user before editing
  `FITTING.md` or `MECHANISM.md`.
- Read `DATA_VERIFICATION.md` 2026-09-14 (second entry) IN FULL, and the
  2026-08-29 entry "every kinetics run in this dataset is a two-channel
  differential measurement".
- **Fits take minutes to tens of minutes.** Run them in the background and
  save every one with `--save`. Never overwrite `data/fits/BnOH_25C_Phosphate.json`
  or the `_M*.json` saves; this round's saves are named `_R<rung>.json`.
- `data/test_fit_kinetics.py` is the slow optimiser suite (~9 min). Run it
  when you change the fitter and once at the end (`run_gates.py --all`), not
  after every edit.

## 1. Where Step 3 stands

### 1.1 What stands and what does not

**Stands, do not redo:**
- `kinetic_model.py`'s six extension terms, each exactly OFF at its default,
  and the six tests in `data/test_kinetic_model.py`.
- `fit_kinetics.MODEL_STAGES`, `--model`, `EXTENDED_INITIAL`.
- `fit_kinetics.ladder_f_test` (reads each rung against `MODEL_PARENTS`, with
  verdicts "not nested", "not nested: stage 1 differs", "nothing added",
  "optimiser failure", "not converged", "earns", "does not earn") and
  `fit_kinetics.ladder_checks` (convergence, correlations, data-vs-model
  substrate order, lag counts, [S]/[buf] collinearity). **Use both on every
  ladder you fit.**

**Does not stand:** the 8954989 entry's verdict, the recovery test as a
recovery test, every F from a pair where the larger model fits worse — and,
per 1.2, every stage-2 fit in `data/fits/`.

### 1.2 The catalysed curves are increments, and the fitter models them as absolutes

Every run is double-beam, and what the reference cuvette omits is what the
recorded curve is net of (`DATA_VERIFICATION.md` 2026-08-29 and 2026-08-31;
`verify_enzyme.reference_design`; `build_manifest.RULINGS` for exps 32, 34-37):

| the reference omits | the recorded curve is | runs in these fits |
|---|---|---|
| the **enzyme** (same S, H2O2, buffer, pH) | the catalytic INCREMENT: sample with catalyst minus the same mixture without | all 11 catalysed runs |
| the **H2O2** | the raw background reaction | 15 of the 17 enzyme-free runs (exps 3 and 6 unclassified, R0.6) |

`fit_kinetics.residuals` compares every curve with
`observable(constants, conditions)` — the FULL signal, background chemistry
included (`k0`, `k_can`, `k3`). That is right for the enzyme-free curves and
wrong for the catalysed ones, whose data have the background removed. The
catalysed model should be

    increment(t) = observable(constants, conditions) - observable(constants, conditions with e0 = 0)

It is not small. In the M0 save's per-curve table, matched compositions on
BnOH: exp 67.1's background rises +0.0332 AU and exp 68.1's increment +0.0329;
exp 69.1's background +0.0450 against exp 71.1's increment +0.0307. The
background is as large as the increment the model is being asked to fit, so
stage 2 has been fitting "catalyst plus a background that was already
subtracted". F5's catalysed misfit, F6's inert loop and every stage-2
constant in `data/fits/` were measured on that observable. They may survive
the correction; nobody knows yet. R0.0 is therefore the first thing done, and
it ends with a stop.

The increment is not just the catalysed chemistry on its own. The catalyst
also acts on aldehyde and peracid the background makes in the SAMPLE cuvette
(steps 6-7 consume background-made peracid), so the difference of two
simulations is the correct observable and "k0 = k_can = k3 = 0" is not.

## 2. The curves in each fit, and what the design limits

Two separate global fits, 104 curves, every curve `build_curves` assigns to the
two blocks (none dropped). Taken from the M0 saves and `build_curves`.

### BnOH, 25 °C, phosphate — 43 curves

| exp | kind (reference omits) | cuvettes | pH | [S], mM | [H2O2], mM | [buf], mM | [enz], mM | run, s |
|---|---|---|---|---|---|---|---|---|
| 3 | background (unclassified, R0.6) | 7 | 6.71 | 1.28–8.98 | 82.5 | 25–85 | 0 | 11025 |
| 6 | background (unclassified, R0.6) | 4 | 6.71 | 3.67–9.16 | 165 | 40–70 | 0 | 17934 |
| 67 | background (H2O2) | 4 | 8.01 | 0.37–7.31 | 122.4 | 85 | 0 | 1178 |
| 69 | background (H2O2) | 4 | 8.01 | 0.21–2.11 | 122.4 | 85 | 0 | 1767 |
| 70 | background (H2O2) | 4 | 8.01 | 0.21–2.11 | 122.4 | 85 | 0 | 1116 |
| 68 | increment (enzyme) | 4 | 8.01 | 0.37–7.31 | 122.4 | 85 | 0.028 | 1829 |
| 71 | increment (enzyme) | 4 | 8.01 | 0.21–2.11 | 122.4 | 85 | 0.028 | 1225 |
| 73 | increment (enzyme) | 4 | 7.50 | 0.32–6.33 | 122.4 | 75 | 0.028 | 2666 |
| 74 | increment (enzyme) | 4 | 8.00 | 0.32–6.33 | 122.4 | 75 | 0.028 | 1984 |
| 83 | increment (enzyme) | 4 | 7.50 | 6.36–14.31 | 122.4 | 75 | 0.028 | 3565 |

### 4-methoxybenzyl alcohol, 40 °C, phosphate — 61 curves

| exp | kind (reference omits) | cuvettes | pH | [S], mM | [H2O2], mM | [buf], mM | [enz], mM | run, s |
|---|---|---|---|---|---|---|---|---|
| 23, 24, 28 | background (H2O2) | 4 each | 7.00 | 2.06–8.25 | 82.5 | 50–80 | 0 | 1764–17934 |
| 25 | background (H2O2) | 2 (cuv 1, 3) | 7.00 | 2.06, 6.19 | 82.5 | 60, 80 | 0 | 1078 |
| 26 | background (H2O2) | 4 identical | 7.00 | 4.13 | 82.5 | 75 | 0 | 392 |
| 27 | background (H2O2) | 4 | 7.00 | 1.03–4.13 | 82.5 | 75–90 | 0 | 1764 |
| 29 | background (H2O2) | 4 | 7.00 | 0.38–1.53 | 82.5 | 75–90 | 0 | 1008 |
| 30, 39 | background (H2O2) | 4 each | 7.00, 6.97 | 0.10–0.38 | 82.5 | 75–90 | 0 | 729–853 |
| 38 | background (H2O2) | 3 | 7.00 | 0.19–0.57 | 82.5 | 65–85 | 0 | 7161 |
| 16 | increment (enzyme) | 4 | 7.00 | 1.85–7.40 | 82.5 | 50–80 | 0.241 | 6542 |
| 32 | increment (enzyme) | 4 | 7.00 | 8.25 | 82.5 | 50–200 | 0.241 | 1767 |
| 34 | increment (enzyme) | 4 | 7.00 | 8.25 | 82.5 | 3.1–25 | 0.241 | 5280 |
| 35 | increment (enzyme) | 4 | 7.50 | 8.25 | 82.5 | 50–200 | 0.241 | 7920 |
| 36 | increment (enzyme) | 4 | 7.53 | 57.9 | 82.5 | 50–200 | 0.241 | 2100 |
| 37 | increment (enzyme) | 4 | 7.53 | 12.23 | 82.5 | 50–200 | 0.270 | 3210 |

**The filenames of exps 32 and 34-37 say `with_NO_E`. They are catalysed**
(ruled to the sheet, `build_manifest.RULINGS`); trust `[enz]`, not the name.

### What that design limits

1. **Peroxide binding cannot be tested.** Every catalysed curve in each block
   sits at ONE [H2O2] (122.4 mM on BnOH, 82.5 mM on 4OMe), so `K4` cannot be
   separated from `k5`. Drop it this round (section 3).
2. **The background is carried to conditions it was never measured at.**
   - pH: the BnOH background is measured at pH 6.71 and 8.01 only, and the
     catalysed exps 73 and 83 sit at 7.50; the 4OMe background is pH 7.00
     only, and the catalysed exps 35-37 sit at 7.50-7.53. The model's
     uncatalysed seed `k0` has NO pH dependence.
   - [S]: exp 36 runs at 57.9 mM, seven times the largest enzyme-free [S]
     (8.25 mM); BnOH exp 83 reaches 14.3 mM against 9.16.
   - [buf]: exp 34 runs at 3.1-25 mM against the enzyme-free 50-90 mM, and
     exps 32 and 35-37 reach 200 mM.
   Once R0.0 makes the catalysed observable an increment, the background
   enters those curves only through the catalyst acting on what the background
   makes (1.2), so how much the extrapolation matters is a number to measure
   (R1.3), not a given.
3. **Only two catalysed compositions have a background run at the same
   composition:** BnOH 67 ↔ 68 and 69/70 ↔ 71. They are the only place the
   model's ABSOLUTE catalysed signal can be checked (R1.2).
4. **On 4OMe the catalysed substrate order rests on one run**, exp 16 (4
   curves); the other five catalysed runs step [buf] at fixed [S]. `km_s` is
   not identified on 4OMe.
5. **Run length spans 46x**, 392 s (exp 26) to 17934 s (exps 6, 28). Only the
   long runs can show a late fall; the short ones barely show a lag. The fit
   weights each curve equally (`weighting="curve"`), not each reading.
6. **[buf] moves with [S] in the background runs** — rising with it in BnOH
   exps 3 and 6, falling against it in all nine 4OMe background ladders
   (median r = -0.974, `ladder_checks`). That is R2's confound.

## 3. What Steps 1 and 2 can and cannot feed in

- **Peroxide binding (`K4`): drop it from this round** (section 2, limit 1).
  Leave `species` as total [H2O2] and `K4` at 0.
- **The activation (`k_act_r`, `K_act`): no prior and no fixed-constant run.**
  The review of Step 2 found the buffer design cannot tell an activating buffer
  from an inhibiting one at the real clock scatter, and that runs 35 and 37 mix
  fit forms across their buffer rungs. The old plan's "fit with Step 2's
  constants FIXED" is dropped for that reason. The activation is tested here
  on its own evidence, and cross-checked against the progress fits' own
  `tau_act_corrected` (Acceptance 6).
- **The product sink:** `slowdown.sink_constants` is the comparison value for
  a CATALYSED 4OMe sink (Acceptance 8).

## 4. Stage R0 — fix the machinery before fitting anything

All in `data/kinetic_model.py`, `data/fit_kinetics.py` and `data/fit_dataset.py`,
each with tests. R0.0 is the one sanctioned change to existing behaviour; every
other R0 change must leave the default path computing what it computes today.

### R0.0 The catalysed observable is an increment — FIRST, and then STOP

1. **`kinetic_model.increment(constants, conditions, times, **kwargs)`**:
   `observable(constants, conditions, times)` minus
   `observable(constants, dataclasses.replace(conditions, e0=0.0), times)`;
   None if either integration fails. Check the name is unique first.
2. **`fit_dataset.Curve` gains `reference_omits`**: "enzyme", "h2o2",
   "substrate", "other" or None. Populate it in `build_curves` from
   `verify_enzyme.analyse()`'s `design` column (read once per build, keyed by
   experiment — do not re-read a sheet per curve). Exps 3 and 6 stay None
   until R0.6 rules them.
3. **`fit_kinetics.residuals` and `_per_curve`** use `increment` when
   `curve.reference_omits == "enzyme"`, `observable` otherwise. Add
   `observation="design"` to `fit_group` / `sequential_fit` with a second
   value `"absolute"` that reproduces today's behaviour exactly, so the old
   fits can be reproduced for comparison. `plot_fit.py` must use the same
   choice (it calls `baseline_like_data`; find every caller of `observable` in
   `data/` and `*/` and route each one).
4. **Guard:** `sequential_fit` raises if a curve with `e0 > 0` has
   `reference_omits` other than "enzyme", or a curve with `e0 == 0` has
   "enzyme". A catalysed curve whose design is unknown is not fitted as
   either.
5. **Tests** — `data/test_kinetic_model.py`:
   - `test_the_increment_vanishes_without_catalyst`: at `e0 = 0` it is zero
     to 1e-12.
   - `test_the_increment_is_the_whole_signal_without_a_background`: with
     `k0 = k_can = k3 = 0` it equals `observable` to 1e-10.
   - `test_the_increment_keeps_the_coupling`: with every constant on, it
     differs from the `k0 = k_can = k3 = 0` observable (the catalyst acting
     on background-made peracid) — assert the difference exceeds 1e-6 mM on a
     BnOH-like planting.
   `data/test_fit_ladder.py` (fast):
   - `test_every_fitted_curve_has_the_design_its_observable_assumes`: over
     both blocks' `build_curves` curves, `e0 > 0` ⇔ `reference_omits ==
     "enzyme"`, and every `e0 == 0` curve is "h2o2" except exps 3 and 6
     (listed by number until R0.6).
   `data/test_fit_kinetics.py`: every synthetic catalysed planting must now
   plant the INCREMENT (update `_synthetic`, `_planted_extended` and the stage-2
   recovery tests), and add `test_absolute_observation_reproduces_the_old_fit`
   on one small planted group.
6. **Refit M0 on both blocks with the increment** (`--model M0`, saved as
   `<block>_R_M0.json`) and report beside the existing M0 save: per stage,
   rms in AU, `k5`, `k6`, their correlation, peak [PBA] on three catalysed
   curves (F6), and `ladder_checks`' substrate order and lag counts. Stage 1
   must be byte-identical to the old M0 stage 1 — the change touches only
   catalysed curves; if it is not, that is a bug.
7. **Commit, write the `DATA_VERIFICATION.md` entry, tell the user, and WAIT.**
   This changes FITTING.md F5 and F6's evidence base; the user decides whether
   FITTING.md is revised before anything else is fitted.

### R0.1 Warm start: a child can never fit worse than its parent

- Add `extra_starts=()` to `fit_group`: a sequence of full `RateConstants`
  whose free parameters are appended to `_guaranteed_points`' list, so they
  always run.
- Add `parent=None` to `sequential_fit`: a loaded save (the dict `to_dict`
  writes, both stages). For each stage, build a start from the parent's
  fitted constants with the child's NEW parameters at a value just inside
  their OFF bound (the value that reproduces the parent's curve to better than
  1e-6 AU), and pass it in `extra_starts`.
- After each stage, if `child.cost > parent.cost * (1 + NESTING_TOLERANCE)`,
  raise `RuntimeError` naming the model and stage. That is a bug, not a
  result.
- Tests, in `data/test_fit_ladder.py` (fast — use a tiny planted problem, 3
  curves, 20 points, `restarts=0`):
  - `test_a_warm_started_child_never_fits_worse`.
  - `test_a_child_above_its_parent_raises` (monkeypatch the cost).

### R0.2 One background for every stage-2 rung

- Add `--background PATH` to `main`: load stage 1's constants from a save and
  skip stage 1 (`sequential_fit(..., background=RateConstants)`), so every
  stage-2 rung sits on the identical frozen base.
- Save the background's path inside each stage-2 save (`"background": PATH`),
  and extend `ladder_f_test`: two stage-2 rows are comparable only if their
  `background` fields match (in addition to the constants check it makes).
- Save the observation (`"observation": "design"` or `"absolute"`) too, and
  make `ladder_f_test` refuse to compare rows whose observations differ.
- Test: `test_stage_two_rungs_share_the_named_background` on planted saves.

### R0.3 A `parents` argument to the ladder readers

`ladder_f_test` and `ladder_checks` take `models=`; add `parents=MODEL_PARENTS`
so the R-ladder can pass its own map. Test on a planted R-ladder.

### R0.4 The terms this round needs, each OFF at its default

In `kinetic_model.py`, extending `RateConstants` exactly as 8954989 did:

| field | default (off) | meaning |
|---|---|---|
| `k_act_f` | 0.0 | forward activation rate on [buf], 1/(mM s) |
| `n_buf_background` | 0.0 | order of the uncatalysed seed in [buf]: `k0 [H2O2][S] ([buf]/BUF_REFERENCE)^n` |
| `n_hoo_background` | 0.0 | order of the uncatalysed seed in [HOO-]: `... ([HOO-]/HOO_REFERENCE)^m` |

- The activation exponent becomes
  `t * (k_act_r * (1 + K_act*[buf]) + k_act_f*[buf])`, and `phi` stays exactly
  1 while `k_act_r` is infinite. **Why:** `k_act_r` and `K_act` correlate at
  1.000 on 4OMe because only their product is identified when the fit wants
  `k_r -> 0`. Fitting `(k_act_r, k_act_f)` with `K_act` held at 0 is the same
  relaxation, `1/tau = k_r + k_f[buf]`, without the degenerate pair.
- `BUF_REFERENCE` (e.g. 75.0 mM) and `HOO_REFERENCE` (the geometric mean of
  the two blocks' [HOO-], computed once and written as a literal with a
  comment saying so): module constants, so `k0` keeps its units and scale.
  Check both names are unique.
- Both background factors apply to the UNCATALYSED seed only, and are 1 when
  their order is 0 or the concentration is NaN.
- Add all three to `PARAMETER_NAMES`, `BOUNDS`, `EXTENDED_INITIAL`, `UNITS`.
  The two orders are fitted LINEARLY (not in `LOG_PARAMETERS`), bounds
  (-2, 2), like `r`. `k_act_f` in log10, bounds (-10, 0).
- Tests in `data/test_kinetic_model.py`:
  - defaults still reproduce `FROZEN` to 1e-10;
  - `test_forward_activation_is_the_same_relaxation`: `(k_act_r=k, K_act=K)`
    and `(k_act_r=k, K_act=0, k_act_f=k*K)` give the same observable to 1e-10;
  - `test_the_background_buffer_order` and `test_the_background_hoo_order`:
    the numerical log-derivative of the initial seed equals the order to 1%.

### R0.5 Rewrite the recovery test so it tests recovery

Replace `test_extended_parameter_recovery` in `data/test_fit_kinetics.py`
(keep its name so the gate list does not change):

- **Real conditions**: the 4OMe/40 C/phosphate curves from `build_curves()`,
  their own `times`, `conditions` and `reference_omits`; replace `absorbance`
  with the planted OBSERVATION (increment for catalysed curves) plus Gaussian
  noise at each curve's own `noise`.
- **Starts away from the truth**: `EXTENDED_INITIAL` at least 0.7 decades
  from every planted log10 value; assert that in the test. `restarts=4`.
  Stage 1 is FITTED, not given.
- **Truth in the regime the ladder will report**: the model lags on some
  catalysed curves (check `peak_model`).
- Noiseless first (every planted parameter to 0.1 decade), then **three** noise
  seeds; "recovered" means all three within 0.3 decade. Assert `k_act_r` and
  `k_act_f` separately. Anything not recovered is later reported as "not
  identified", never as a value.
- Runs under `--all` only. If it exceeds 15 minutes, cut to one run per
  experiment and say so in the docstring.

### R0.6 Rule exps 3 and 6's reference design

`verify_enzyme.reference_design` cannot classify exp 3 (the cuvette table does
not read as two halves) or exp 6 ("other"), though the 2026-08-31 entry lists
both as omitting the H2O2.
- Open both sheets (`manifest` → `xls_file`, `Sheet1`) and read the cuvette
  table by eye: which rows are the measured cuvettes, which the references,
  and what the references omit.
- Record the ruling twice, as the repository requires: a `RULINGS` entry in
  `data/build_manifest.py` with the reasoning, and a dated
  `DATA_VERIFICATION.md` entry with the evidence. Extend `reference_design`
  if the layout is a readable variant; do not hand-patch the result.
- If either run's reference omits anything other than the H2O2, it is not a
  raw background curve: drop it from stage 1 and say so. Fit stage 1 with and
  without exps 3 and 6 either way (R2), because they are also the only pH
  6.71, only 11-18 ks, only moving-buffer BnOH background runs.

## 5. Stage R1 — look before adding terms

After R0 is committed and the user has seen R0.0's M0 refit.

### R1.1 The residual shape

- Add `fit_kinetics.residual_shape(block, model, stage=2, directory="data/fits")`.
  For each curve of that save re-simulate with the saved constants and the
  save's observation, and report the mean residual (data − model) in the
  first, middle and last third of the run, in AU and in units of the curve's
  noise, plus `net_data/net_model`. One row per curve, sorted by rms.
- Print it for the R0.0 M0 refit of both blocks, both stages. For the top ten
  curves answer, in the entry:
  1. AMPLITUDE (`net_data/net_model` far from 1, one sign in every third)?
     Then it is an extinction / `r` / concentration question, not a term.
  2. EARLY (first third): a lag the model lacks or overshoots?
  3. LATE (last third): a fall the model lacks (the product sink)?
  4. Concentrated in a few experiments? Name them and check `KNOWN_EXCLUSIONS`
     / `RULINGS` in `data/build_manifest.py`.
- **STOP and report** if the answer is (1) or (4) for most of the cost.
- Test: `test_residual_shape_is_zero_on_a_planted_fit`.

### R1.2 The matched compositions: the only absolute check

BnOH exp 68 has the composition of background exp 67, and exp 71 that of 69
and 70 — same pH, [S] per cuvette, [H2O2] and [buf].
- Add `fit_kinetics.matched_compositions(block, save)`. For each matched
  cuvette pair: interpolate both measured curves onto the shorter run's times
  (after each is baseline-subtracted as recorded) and form
  `absolute_data = background + increment`; set it beside
  `absolute_model = observable(E0)` and `increment_model = increment(E0)` at
  the save's constants. Report net rises at the common end time and the rms of
  each difference.
- What it can show: whether the model gets the ABSOLUTE catalysed signal right
  as well as the increment. What it cannot: the two runs are different days
  and cuvettes, and the background run's reference (no H2O2) is not the
  increment run's (no enzyme); static absorbance cancels in the baseline
  subtraction, drifts do not. Say so. Treat a mismatch beyond twice the
  exp 69-against-70 difference (two background runs at one composition) as
  real.
- Test: on a planted pair made from one set of constants,
  `absolute_data == absolute_model` to 1e-10.

### R1.3 How much the background still decides on each catalysed curve

- Add `fit_kinetics.background_leverage(block, save)`. For each catalysed
  curve: the net increment at the save's constants, and the same with
  `k0 = k_can = k3 = 0` (r kept). Their difference is what the frozen
  background contributes to the increment through the catalyst.
- Report the share per curve and flag every curve above 10%. Those are the
  curves where limit 2's extrapolations (pH, [S], [buf]) reach the catalysed
  constants; list them with their experiment's pH, [S] and [buf] against the
  background's measured range.

## 6. Stage R2 — the background, with its confound

Stage 1 only, both blocks, each rung warm-started from its parent. Stage 1 is
untouched by R0.0, but refit it here with R0's warm start.

| rung | frees, added to STAGE_ONE | question |
|---|---|---|
| B0 | — | baseline |
| B1 | `km_s_background` | a background saturating in [S] |
| B2 | `n_buf_background` | a background with an order in [buf] |
| B3 | `km_s_background`, `n_buf_background` | both |
| B5 | winner + `n_hoo_background` | **BnOH only**: does `k0` depend on pH? |

- B1 and B2 are NOT nested in each other; same parameter count, so compare
  their costs directly. B3 against each is an F.
- Read with `ladder_f_test(parents={"B1": "B0", "B2": "B0", "B3": "B1"})` (and
  again with B3 against B2) and `ladder_checks`.
- **B5 is confounded and must say so**: the BnOH pH 6.71 runs (exps 3, 6)
  also differ from the pH 8.01 runs (67, 69, 70) in [H2O2] (82.5/165 against
  122.4), in buffer range and in run length (11-18 ks against 1.1-1.8 ks). A
  pH order there is an order in all of that. 4OMe's background is one pH, so
  B5 is not identified there; carry `k0` unchanged to the pH 7.5 catalysed
  runs and let R3's sensitivity cut (below) price it.
- On BnOH, where 3 of 5 background ladders hold [buf] fixed, the fixed-buffer
  runs are the discriminating ones: also report `ladder_checks`' substrate
  order on exps 67, 69, 70 alone, data against B1 and B2.
- Fit every BnOH rung with and without exps 3 and 6 (R0.6).

| outcome | reading |
|---|---|
| B1 and B2 within a factor 1.05 in cost, B3 does not earn | **not decided** — carry B1 AND B2 into R3 and report both |
| B2 clearly lower, and its substrate order matches the data within 2 sigma | the background has a buffer order, not a saturation |
| B1 clearly lower, and its substrate order matches the data within 2 sigma | the background saturates in [S] |
| the winner's substrate order misses the data's by more than 2 sigma | it wins by bending something else; not accepted |

- Add `k_sink` to stage 1 as rung **B4** (parent = the R2 winner) only as a
  CONTROL: a background sink that earns would contradict `product_fate`'s
  enzyme-free clock, and must be reported, not carried.
- **On 4OMe, run R3 on BOTH B1 and B2 regardless of R2's outcome**: exps 32
  and 34-37 take [buf] from 3.1 to 200 mM against a background measured at
  50-90, so the two backgrounds diverge exactly there.

## 7. Stage R3 — the catalysed ladder on one background

Stage 2 only, observation "design" (increments), `--background` = the R2
choice, each rung warm-started from its parent.

| rung | parent | frees, added to STAGE_TWO | question |
|---|---|---|---|
| C0 | — | — | catalysed baseline on this background |
| C1 | C0 | `km_s` | catalysed substrate saturation (**BnOH only**; 4OMe: "not identified") |
| C2 | best of C0/C1 | `k_sink` | the catalysed product sink |
| C3 | C2 | `k_act_r` | activation on one clock (every block) |
| C4 | C3 | `k_act_f` | the clock moving with [buf] (**4OMe only**: no BnOH catalysed run moves [buf]) |
| C5 | best so far | `r` | does `r` move once a lag is supplied? |

- C5 frees `r` in stage 2. `r` also shapes the background, so report the
  enzyme-free rms at C5's `r` (re-simulate stage 1 with it). If the background
  worsens by more than 10% in rms, `r` is being bought by the catalysed curves
  at the enzyme-free curves' expense — say so and defer to R4.
- "Best of" means the lower cost among COMPARABLE rows whose verdict is
  "earns"; a term that does not earn is not carried.
- Read every rung with `ladder_f_test(parents=...)`, `ladder_checks` and
  `background_leverage`.

**Sensitivity cuts on the final R3 model** (refit stage 2 only, same
background, warm-started from the full fit):

| cut | drops | why |
|---|---|---|
| S1 | 4OMe exp 36; BnOH exp 83 | [S] beyond the background's measured range |
| S2 | 4OMe exps 35, 36, 37; BnOH exps 73, 83 | pH outside the background's measured pH |
| S3 | 4OMe exps 32, 34 | [buf] outside 50-90 mM (4OMe only) |

Report every constant's shift in units of its standard error. A shift beyond
2 is a constant decided by extrapolated background; say which. S3 removes
most of 4OMe's buffer axis, so `k_act_f` is expected to go unidentified there —
that is the point of the cut, not a failure.

## 8. Stage R4 — the joint fit

Take the final R3 model and fit ALL curves of the block at once, both stages'
parameters free, each curve with its design's observation, starting from the
sequential optimum (`extra_starts`). This is the only place `r` is fitted to
both kinds of curve together. Report rms per kind of curve, and whether any
constant moved by more than its standard error from the sequential value.

## 9. Acceptance

State each as **met** or **not met**, with its number, in the entry. Nothing
counts from a row whose verdict is not "earns", from an unconverged fit, or
from a fit with observation "absolute".

1. **The observable matches the design**: every fitted catalysed curve used
   the increment (the R0.0 guard and test).
2. **No larger rung fits worse than its parent**, and every comparable row
   converged.
3. **rms in AU** against C0 on the same background, per kind of curve.
4. **Substrate order**: the final model's (`ladder_checks`) within 2 sigma of
   the data's, per kind of curve. On 4OMe enzyme-free, state R2's reading.
5. **Lag counts**: model within 3 curves of the data per kind of curve.
6. **Activation clock**: for catalysed curves whose `act_sink_kind_corrected`
   is a lag kind ("lag" or "lag then sink") in
   `scope.frame(tuple(save["stage_2"]["experiments"]))`, the median ratio of the
   model's `1/(k_act_r + k_act_f*[buf])` to that curve's `tau_act_corrected`
   within a factor of 2. Exclude, and count, runs whose rungs mix fit forms.
7. **r**: at C5 or R4, inside 0.08-0.33, or state where it went.
8. **Sink**: C2's `k_sink` against `slowdown.sink_constants` for 4OMe at 40 C,
   within its error.
9. **Machinery (F6)**: peak [PBA] at the final constants on three
   representative catalysed curves.
10. **Matched compositions** (R1.2): the final model's absolute catalysed
    signal on BnOH 68 and 71 within the exp 69-against-70 spread.
11. **Sensitivity**: no reported constant shifts by more than 2 standard
    errors under S1-S3, or each one that does is named as extrapolation-bound.

## 10. Stop conditions

Stop and report to the user, without continuing, if:
- R0.0 is done (always a stop: the user decides about FITTING.md);
- R0.0's stage 1 is not byte-identical to the old M0 stage 1;
- any existing test fails after an R0 change other than the synthetic
  plantings R0.0 says to update;
- R0.5's recovery fails for a parameter a later stage would report (then do
  not report it; if it is `k_act_r` or `k_act_f`, stop);
- R0.6 finds exp 3 or 6 is not a raw background curve (rule it first);
- R1.1 finds the misfit is an amplitude error or a handful of runs;
- R0.1's guard raises on a real fit;
- a fitted constant sits at a bound that is not its OFF end and you cannot
  explain why.

## 11. Deliverables, reporting, committing

- [ ] **R0.0**: `increment`, `Curve.reference_omits`, the observation switch
      and guard, updated plantings, the tests, the M0 refit of both blocks,
      a `DATA_VERIFICATION.md` entry. **One commit. Tell the user and WAIT.**
- [ ] **R0.1-R0.6**: warm start, `--background`, `parents=`, the three terms,
      the rewritten recovery test, the exps 3/6 ruling, and their tests.
      `run_gates.py` green; `run_gates.py --all` once. **One commit.**
- [ ] **R1**: `residual_shape`, `matched_compositions`,
      `background_leverage`, their tests, and an entry answering R1.1-R1.3.
      **One commit. Tell the user and wait.**
- [ ] **R2-R4**: the saves `data/fits/<block>_R<rung>.json`, and one entry
      with, per block: the R2 table and reading, the R3 `ladder_f_test`,
      `ladder_checks` and `background_leverage` tables, the S1-S3 shifts, the
      R4 comparison, and the acceptance list. **One commit.**
- Entries follow `PLAN_MECHANISM_NEXT_STEPS.md` section 5 and end "Nothing is
  adopted".
- Propose any `FITTING.md` / `MECHANISM.md` text to the user; do not commit it
  unasked.
- After each commit, tell the user in plain words: the verdict, the one or two
  numbers that decide it, and what could still overturn it.
