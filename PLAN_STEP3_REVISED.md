# Step 3, revised — a search that can answer the question

Handover plan, written 2026-09-14 after the review in `DATA_VERIFICATION.md`
2026-09-14 (second entry), and extended the same day twice: once the curves in
each fit were listed (section 2, R0.0), and again to replace the step-by-step
ladder with a bounded search over mechanism modules, ranked by leave-one-run-out
cross-validation (R2). It REPLACES stages 3.2-3.3 and the acceptance list of
`PLAN_MECHANISM_NEXT_STEPS.md` section 4. Stages 3.0-3.1 of that plan are done
(commit 8954989) and are not repeated.

**Read section 1.2 before anything else.** The fitter has compared the
catalysed curves with the wrong quantity since it was written, and every
stage-2 result — including FITTING.md F5 and F6 — was measured that way.

## Contents

0. Read first
1. Where Step 3 stands, and why a search
2. The curves in each fit, and what the design limits
3. What Steps 1 and 2 can and cannot feed in
4. Stage R0 — fix the machinery before fitting anything
5. Stage R1 — look before searching
6. Stage R2 — the module search
7. Stage R3 — robustness of the tie set
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
- **Fits take minutes to tens of minutes, and R2 runs hundreds of them.** Run
  them in the background, save every one, and make every long job RESUMABLE
  (skip a fit whose save exists). Never overwrite `data/fits/BnOH_25C_Phosphate.json`
  or the `_M*.json` saves.
- **Parallel fits must pin BLAS to one thread** (`OMP_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` set before numpy is imported
  in each worker). Several processes each spawning a full BLAS thread pool
  oversubscribes the machine and has doubled gate time in this repository
  before.
- `data/test_fit_kinetics.py` is the slow optimiser suite (~9 min). Run it
  when you change the fitter and once at the end (`run_gates.py --all`), not
  after every edit.

## 1. Where Step 3 stands, and why a search

### 1.1 What stands and what does not

**Stands, do not redo:**
- `kinetic_model.py`'s six extension terms, each exactly OFF at its default,
  and the six tests in `data/test_kinetic_model.py`.
- `fit_kinetics.MODEL_STAGES`, `--model`, `EXTENDED_INITIAL`.
- `fit_kinetics.ladder_checks` (convergence, correlations, data-vs-model
  substrate order, lag counts, [S]/[buf] collinearity) — R2 runs it on every
  model in the final tie set. `fit_kinetics.ladder_f_test` stays for reading
  the old `_M*` saves; R2 does not use F tests.

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

### 1.3 Why a search, and why a bounded one

The first attempt added terms one at a time and let an F test decide. Three
things make that the wrong instrument here:
- **The answer depends on the order terms are added in**, and two of the
  questions are not nested at all (a background saturating in [S] against one
  with an order in [buf]).
- **F > 12 cannot discriminate at this misfit.** The fits sit 23-650x the
  noise over thousands of serially correlated readings; one term cleared the
  bar at F = 259,201 by bending the substrate order to zero.
- **A term that absorbs an artefact in the fitted runs looks like chemistry.**
  The only test that catches it is predicting runs the fit did not see.

A fully agnostic search over arbitrary reactions is not the fix: the
observable is one absorbance channel, the peracid and every catalyst state are
unseen, and the design lacks the contrasts to separate most networks, so an
open search returns a crowd of equivalent networks. R2 instead enumerates
COMBINATIONS OF A SMALL LIBRARY OF CHEMICAL MODULES, every one already a
switchable term in `kinetic_model.py`, scores each combination by how well it
predicts held-out runs, and reports the SET of combinations the data cannot
tell apart. That set, and which modules are in all of it, in none of it, or
undecided, is the mechanistic result.

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
   separated from `k5`. It is not in R2's library.
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
   curves); the other five catalysed runs step [buf] at fixed [S]. Catalysed
   substrate saturation is not in 4OMe's library.
5. **Run length spans 46x**, 392 s (exp 26) to 17934 s (exps 6, 28). Only the
   long runs can show a late fall; the short ones barely show a lag. The fit
   weights each curve equally (`weighting="curve"`), not each reading.
6. **[buf] moves with [S] in the background runs** — rising with it in BnOH
   exps 3 and 6, falling against it in all nine 4OMe background ladders
   (median r = -0.974, `ladder_checks`). R2's background library carries both
   readings so the search can say whether they separate.
7. **Leave-one-run-out has few folds**: 5 and 10 background runs, 5 and 6
   catalysed runs. Its tie rule (R2.3) is therefore generous, and some held-out
   runs are extrapolations by construction (exp 34's buffer, exp 36's [S]).
   That is a feature — a model that extrapolates badly should lose — but the
   per-fold table must be reported so a single run's weight is visible.

## 3. What Steps 1 and 2 can and cannot feed in

- **Peroxide binding (`K4`): out of this round** (section 2, limit 1). Leave
  `species` as total [H2O2] and `K4` at 0.
- **The activation (`k_act_r`, `K_act`): no prior and no fixed-constant run.**
  The review of Step 2 found the buffer design cannot tell an activating buffer
  from an inhibiting one at the real clock scatter, and that runs 35 and 37 mix
  fit forms across their buffer rungs. The activation is tested here on its own
  evidence, and cross-checked against the progress fits' own
  `tau_act_corrected` (Acceptance 7).
- **The product sink:** `slowdown.sink_constants` is the comparison value for
  a CATALYSED 4OMe sink (Acceptance 9).

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
   choice (find every caller of `observable` in `data/` and `*/` and route
   each one).
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
     differs from the `k0 = k_can = k3 = 0` observable by more than 1e-6 mM on
     a BnOH-like planting.
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

### R0.1 Warm starts, and a guard that a larger model never fits worse

- Add `extra_starts=()` to `fit_group`: a sequence of full `RateConstants`
  whose free parameters are appended to `_guaranteed_points`' list, so they
  always run.
- Add `nested_start(fitted, free_names)` in `fit_kinetics`: from a fitted
  smaller model's constants, the start for a larger one with every added
  parameter just inside its OFF bound (the value that reproduces the smaller
  model's curve to better than 1e-6 AU). R2 uses it to warm-start every model
  from each of its already-fitted sub-models.
- Add `nesting_violation(cost, sub_costs)`: True if `cost` exceeds the
  smallest fitted sub-model cost by more than `NESTING_TOLERANCE`. R2 refits
  such a model once with `restarts` doubled, and flags it if it still
  violates.
- Tests, in `data/test_fit_ladder.py` (fast — a tiny planted problem, 3
  curves, 20 points, `restarts=0`):
  - `test_a_warm_started_larger_model_never_fits_worse`.
  - `test_a_nesting_violation_is_detected` (a hand-built cost pair).

### R0.2 A frozen background for the catalysed search

- `sequential_fit(..., background=RateConstants)` skips stage 1 and uses those
  constants, so every catalysed model sits on the identical frozen base; add
  `--background PATH` to `main` for the same from a save.
- Every save records `"background"` (a path or a background model id) and
  `"observation"` (`"design"` or `"absolute"`).
- Test: `test_a_frozen_background_is_not_refitted` on a planted group.

### R0.3 Per-fit bound overrides

`fit_group(..., bounds=None)`: a dict overriding `BOUNDS` for named parameters
in this fit only. R2 needs it to hold `r` inside the spectroscopic range as a
module. Test that an override is respected and `BOUNDS` itself is unchanged.

### R0.4 The terms the library needs, each OFF at its default

In `kinetic_model.py`, extending `RateConstants` exactly as 8954989 did:

| field | default (off) | meaning |
|---|---|---|
| `k_sink_e` | 0.0 | catalyst-mediated product sink: `v_sink = (k_sink + k_sink_e·E0)·[A]`, 1/(mM s) |
| `k_act_f` | 0.0 | forward activation rate on [buf], 1/(mM s) |
| `n_buf_background` | 0.0 | order of the uncatalysed seed in [buf]: `k0 [H2O2][S] ([buf]/BUF_REFERENCE)^n` |
| `n_hoo_background` | 0.0 | order of the uncatalysed seed in [HOO-]: `... ([HOO-]/HOO_REFERENCE)^m` |

- **Why `k_sink_e` and not `k_sink`:** `k_sink` is one constant shared by
  both cuvettes, so in the increment observable it drains the reference too,
  and a background fitted without it becomes inconsistent. `product_fate`'s
  sink is a CATALYSED finding — the enzyme-free curves decline on a clock
  instead — so the catalysed module is a sink proportional to the catalyst.
  `k_sink` stays, as the background sink control (R3).
- The activation exponent becomes
  `t * (k_act_r * (1 + K_act*[buf]) + k_act_f*[buf])`, and `phi` stays exactly
  1 while `k_act_r` is infinite. `k_act_r` and `K_act` correlate at 1.000 on
  4OMe because only their product is identified when the fit wants
  `k_r -> 0`; `(k_act_r, k_act_f)` with `K_act` held at 0 is the same
  relaxation, `1/tau = k_r + k_f[buf]`, without the degenerate pair.
- `BUF_REFERENCE` (e.g. 75.0 mM) and `HOO_REFERENCE` (the geometric mean of
  the two blocks' [HOO-], computed once and written as a literal with a
  comment saying so): module constants. Check every new name is unique.
- Both background factors apply to the UNCATALYSED seed only, and are 1 when
  their order is 0 or the concentration is NaN.
- Add all four to `PARAMETER_NAMES`, `BOUNDS`, `EXTENDED_INITIAL`, `UNITS`.
  The two orders are fitted LINEARLY (not in `LOG_PARAMETERS`), bounds
  (-2, 2), like `r`. `k_sink_e` in log10, bounds (-10, 2); `k_act_f` in
  log10, bounds (-10, 0).
- Tests in `data/test_kinetic_model.py`:
  - defaults still reproduce `FROZEN` to 1e-10;
  - `test_the_catalysed_sink_leaves_the_background_alone`: with only
    `k_sink_e` on, `observable` at `e0 = 0` is unchanged to 1e-12, and a pulse
    of aldehyde at `e0 > 0` decays as `exp(-k_sink_e*e0*t)`;
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
- **Truth in the regime the search will report**: the planting switches on
  the catalysed sink and the activation, and the model lags on some catalysed
  curves (check `peak_model`).
- Noiseless first (every planted parameter to 0.1 decade), then **three** noise
  seeds; "recovered" means all three within 0.3 decade. Assert `k_act_r`,
  `k_act_f` and `k_sink_e` separately. Anything not recovered is later
  reported as "not identified", never as a value.
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
  raw background curve: drop it from stage 1 and say so. Either way R3 refits
  the tie set without exps 3 and 6, because they are also the only pH 6.71,
  only 11-18 ks, only moving-buffer BnOH background runs.

## 5. Stage R1 — look before searching

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
     Then it is an extinction / `r` / concentration question, not a module.
  2. EARLY (first third): a lag the model lacks or overshoots?
  3. LATE (last third): a fall the model lacks (the product sink)?
  4. Concentrated in a few experiments? Name them and check `KNOWN_EXCLUSIONS`
     / `RULINGS` in `data/build_manifest.py`.
- **STOP and report** if the answer is (1) or (4) for most of the cost. A
  search over reaction modules cannot fix an amplitude or a bad run; it will
  pick whichever module best disguises it.
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
- **Also report the block-wide median share.** If it is below 10% on both
  blocks, say in the entry that the two-axis pyrophosphate block — excluded
  so far because it has no background runs (FITTING.md F7) — may become
  fittable as increments. Do not fit it; that is the user's decision after
  this plan.

## 6. Stage R2 — the module search

New module `data/mechanism_search.py` (check the name and every top-level
name in it are unique), gate `data/test_mechanism_search.py`. It runs in two
phases on each block: the BACKGROUND search on the enzyme-free curves, then the
CATALYSED search on the increments, on each background the first phase could
not rule out.

### R2.1 The library

Each module is a set of parameters it frees, parameters it fixes, bound
overrides, which modules it requires, and which blocks it applies to. A model
is a set of modules; its free parameters are the base's plus every module's.

**Background phase** — base frees `k0`, `r`, with `k_can = k3 = 0` fixed:

| module | frees / sets | chemistry | blocks |
|---|---|---|---|
| `can` | frees `k_can`, `k3` | steps 1-3, the uncatalysed peracid loop | both |
| `sat` | frees `km_s_background` | the seed saturates in [S] | both |
| `buf` | frees `n_buf_background` | the seed has an order in [buf] | both |
| `hoo` | frees `n_hoo_background` | the seed has an order in [HOO-] | **BnOH only** (4OMe's background is one pH) |
| `rspec` | bounds `r` to (0.08, 0.33) | `r` held to the spectroscopic range | both |

BnOH: 2^5 = **32** background models. 4OMe: 2^4 = **16**.
`hoo` is confounded and the entry must say so: the BnOH pH 6.71 runs (exps 3,
6) also differ from the pH 8.01 runs (67, 69, 70) in [H2O2] (82.5/165 against
122.4), in buffer range and in run length (11-18 ks against 1.1-1.8 ks).

**Catalysed phase** — base frees `k5`, with `k6 = 0` fixed, on a frozen
background:

| module | frees | chemistry | requires | blocks |
|---|---|---|---|---|
| `loop` | `k6` | steps 6-7, the catalysed peracid loop | the background has `can` | both |
| `sat` | `km_s` | catalysed substrate saturation | — | **BnOH only** (limit 4) |
| `sink` | `k_sink_e` | the catalysed product sink | — | both |
| `act` | `k_act_r` | activation on one clock | — | both |
| `actbuf` | `k_act_f` | the clock moving with [buf] | `act` | **4OMe only** (no BnOH catalysed run moves [buf]) |

`loop` on a background without `can` makes no peracid and does nothing, so
`enumerate_models` skips it there (and the test checks that). BnOH: up to
2^4 = **16** catalysed models per background; 4OMe: 2 × 2 × 3 = **12**.

### R2.2 Fitting and cross-validation

- `enumerate_models(block, phase, background_modules=None)` → the models,
  respecting `requires` and `blocks`. `model_id(phase, modules)` → a stable
  string such as `bg:buf+can` or `cat:act+sink@bg:buf+can`.
- **Full-data fit** of every model, in order of increasing size, each
  warm-started (`nested_start`) from every already-fitted model it contains.
  Background phase: stage 1 on all enzyme-free curves. Catalysed phase: stage
  2 on all catalysed curves, background frozen (R0.2) at that background
  model's full-data constants. After fitting, apply `nesting_violation` (R0.1).
- **Leave-one-run-out**: for each run of the phase (background: 5 BnOH, 10
  4OMe; catalysed: 5 BnOH, 6 4OMe), refit the model on the other runs
  warm-started from its full-data fit (`restarts=1`), then predict the
  held-out run's curves with that run's own conditions and the design's
  observation. In the catalysed phase the background stays frozen at its
  full-data fit — the folds hold out CATALYSED runs, which the background
  never saw.
- **Score** = the held-out run's cost in the fitter's own weighting
  (`weighting="curve"`, residual / (noise·sqrt(n))), summed over its curves;
  a model's CV score is the sum over folds. Save the per-fold scores and the
  held-out rms in AU.
- A fit whose integration fails (`FAILURE_RESIDUAL`) scores infinity for that
  fold; a model with any failed fold is listed but not ranked.
- **Saves**: `data/fits/search/<block>/<model_id>.json` and
  `.../<model_id>__without_<exp>.json`, `to_dict` plus `modules`, `phase`,
  `background`, `observation`, `held_out`, `cv_score`. `search(...)` skips any
  fit whose save exists, so the job is resumable.
- **Parallelism**: `search(..., workers=N)` with `concurrent.futures.
  ProcessPoolExecutor`, BLAS pinned to one thread per worker (section 0). Fold
  refits of one model are independent and parallelise; a model waits for its
  sub-models' full-data fits.

### R2.3 Ranking, the tie set, and the module verdicts

- `rank(table)`: models by CV score.
- `tie_set(table)`: every model `m` whose per-fold differences against the
  best, `d_f = score(m, f) − score(best, f)`, have
  `mean(d) <= 2 * sd(d) / sqrt(folds)`. With 5-10 folds this is generous, and
  it is meant to be: the entry reports the tie set, not a winner.
- `module_verdicts(table)`: for each module, the share of tie-set models that
  contain it — **required** (in all), **rejected** (in none), **undecided**
  (in some). That table is the mechanistic result.
- **Identifiability flags** on every tie-set model's full-data fit: any
  |correlation| > 0.99 pair, any constant at a bound that is not its OFF end,
  any parameter R0.5 did not recover. A flagged model stays in the tie set but
  its constants are reported as "not identified".
- Run `ladder_checks` (substrate order, lag counts) and `background_leverage`
  on every tie-set model.
- **Which backgrounds go on to the catalysed phase**: the background tie set,
  capped at the 3 best. On 4OMe, always include the best background containing
  `sat` and the best containing `buf` if either is outside the cap: exps 32
  and 34-37 take [buf] from 3.1 to 200 mM against a background measured at
  50-90, so those two readings diverge exactly under the catalysed curves.

### R2.4 Budget: pilot first, and a greedy fallback

- **R2.0 pilot** before any search: time one full-data fit and one fold refit
  of the LARGEST model in each phase on each block. Project the total,
  `sum over phases of models × (folds + 1) × fit time / workers`, and write it
  in the entry. The full enumeration is about 900 fits.
- If the projection exceeds **48 hours of wall time**, switch that phase to a
  **greedy search**: start from the base; repeatedly add the module (respecting
  `requires`) that lowers the CV score most, stopping when no addition leaves
  the tie set of the current best; then repeatedly remove modules while the
  result stays in the tie set. Run it a second time starting from the FULL
  model and removing. Report both end points; if they differ, the difference
  is itself a result.

### R2.5 Tests — `data/test_mechanism_search.py` (fast)

- `test_the_library_respects_its_constraints`: model counts 32/16 background
  and, on a background with `can`, 16/12 catalysed; no `loop` without `can`;
  no `actbuf` without `act`; no `hoo` or `actbuf` on the wrong block.
- `test_model_ids_round_trip` and every model's free parameters ⊆ `BOUNDS`.
- `test_every_run_is_held_out_once_and_never_trained_on`.
- `test_the_tie_rule`: on hand-built fold scores, a model within the rule is
  tied and one beyond it is not.
- `test_module_verdicts`: hand-built tie set → required / rejected /
  undecided as constructed.
- `test_a_finished_fit_is_not_refitted`: with a save present, `search` makes
  no fit call (monkeypatch the fitter).

And in `data/test_fit_kinetics.py` (slow): `test_the_search_finds_a_planted_mechanism`
— plant a small catalysed design with `sink` and `act` on and `loop` off,
run the catalysed phase on it, and assert `sink` and `act` are **required**
and `loop` is not.

## 7. Stage R3 — robustness of the tie set

On every model in each block's final (catalysed) tie set, capped at the 3 best:

| check | refit | why |
|---|---|---|
| S1 | catalysed phase without 4OMe exp 36 / BnOH exp 83 | [S] beyond the background's measured range |
| S2 | without 4OMe exps 35, 36, 37 / BnOH exps 73, 83 | pH outside the background's measured pH |
| S3 | without 4OMe exps 32, 34 | [buf] outside 50-90 mM (4OMe only) |
| S4 | background phase without BnOH exps 3 and 6 | R0.6, and they are the only pH 6.71 background |
| B-sink | each tie-set background plus `k_sink` | a background sink that improves the CV score contradicts `product_fate`'s enzyme-free clock; report, never carry |

Report every constant's shift in units of its standard error, and whether each
module verdict changes. A shift beyond 2, or a verdict that flips, is decided
by extrapolated background or by those runs; say which. S3 removes most of
4OMe's buffer axis, so `actbuf` is expected to become undecided there — that is
the point of the cut, not a failure.

## 8. Stage R4 — the joint fit

For each model in the final tie set (capped at 3), fit ALL curves of the
block at once, both phases' parameters free, each curve with its design's
observation, starting from the sequential optimum (`extra_starts`). This is
the only place `r` is fitted to both kinds of curve together. Report rms per
kind of curve, and whether any constant moved by more than its standard error
from the sequential value.

## 9. Acceptance

State each as **met** or **not met**, with its number, in the entry, for every
model in each block's final tie set. Nothing counts from an unconverged fit, a
nesting violation that survived its refit, or a fit with observation
"absolute".

1. **The observable matches the design**: every fitted catalysed curve used
   the increment (the R0.0 guard and test).
2. **The search discriminates**: the catalysed tie set holds fewer than half
   of the models searched. If it holds more, the data cannot rank these
   modules at all — report that as the result.
3. **The module verdict table** is reported per block, with S1-S4's changes.
4. **Held-out rms in AU**: the best model against the base model of each phase.
5. **Substrate order**: the model's (`ladder_checks`) within 2 sigma of the
   data's, per kind of curve. On 4OMe enzyme-free, state whether `sat` and
   `buf` separated.
6. **Lag counts**: model within 3 curves of the data per kind of curve.
7. **Activation clock**: for catalysed curves whose `act_sink_kind_corrected`
   is a lag kind ("lag" or "lag then sink") in
   `scope.frame(tuple(save["stage_2"]["experiments"]))`, the median ratio of the
   model's `1/(k_act_r + k_act_f*[buf])` to that curve's `tau_act_corrected`
   within a factor of 2. Exclude, and count, runs whose rungs mix fit forms.
8. **r**: whether `rspec` is required, rejected or undecided, and `r` at R4.
9. **Sink**: if `sink` is required, `k_sink_e·[enz]` against
   `slowdown.sink_constants` for 4OMe at 40 C, within its error.
10. **Machinery (F6)**: whether `loop` is required; peak [PBA] at the final
    constants on three representative catalysed curves.
11. **Matched compositions** (R1.2): the model's absolute catalysed signal on
    BnOH 68 and 71 within the exp 69-against-70 spread.
12. **Sensitivity**: no reported constant shifts by more than 2 standard
    errors under S1-S4, or each one that does is named as extrapolation-bound.

## 10. Stop conditions

Stop and report to the user, without continuing, if:
- R0.0 is done (always a stop: the user decides about FITTING.md);
- R0.0's stage 1 is not byte-identical to the old M0 stage 1;
- any existing test fails after an R0 change other than the synthetic
  plantings R0.0 says to update;
- R0.5's recovery fails for a parameter the search would report (then that
  parameter is "not identified"; if it is `k_act_r`, `k_act_f` or `k_sink_e`,
  stop);
- R0.6 finds exp 3 or 6 is not a raw background curve (rule it first);
- R1.1 finds the misfit is an amplitude error or a handful of runs;
- R2.0's pilot projects more than 48 hours even for the greedy search;
- `test_the_search_finds_a_planted_mechanism` fails — then the search cannot
  find what is there, and nothing it reports on real data means anything;
- a fitted constant sits at a bound that is not its OFF end and you cannot
  explain why.

## 11. Deliverables, reporting, committing

- [ ] **R0.0**: `increment`, `Curve.reference_omits`, the observation switch
      and guard, updated plantings, the tests, the M0 refit of both blocks,
      a `DATA_VERIFICATION.md` entry. **One commit. Tell the user and WAIT.**
- [ ] **R0.1-R0.6**: warm starts and the nesting guard, the frozen
      background, bound overrides, the four terms, the rewritten recovery
      test, the exps 3/6 ruling, and their tests. `run_gates.py` green;
      `run_gates.py --all` once. **One commit.**
- [ ] **R1**: `residual_shape`, `matched_compositions`,
      `background_leverage`, their tests, and an entry answering R1.1-R1.3.
      **One commit. Tell the user and wait.**
- [ ] **R2.0-R2.5**: `mechanism_search.py`, its fast gate and the slow planted
      search, the pilot's projection in an entry. **One commit. Tell the user
      the projected wall time and wait for the go-ahead to run it.**
- [ ] **R2-R4 results**: the saves under `data/fits/search/`, and one entry
      with, per block: the background and catalysed rankings with per-fold
      scores, the tie sets, the module verdict table, `ladder_checks` and
      `background_leverage` on the tie set, R3's shifts and verdict changes,
      the R4 comparison, and the acceptance list. **One commit.**
- Entries follow `PLAN_MECHANISM_NEXT_STEPS.md` section 5 and end "Nothing is
  adopted".
- Propose any `FITTING.md` / `MECHANISM.md` text to the user; do not commit it
  unasked.
- After each commit, tell the user in plain words: the verdict, the one or two
  numbers that decide it, and what could still overturn it.
