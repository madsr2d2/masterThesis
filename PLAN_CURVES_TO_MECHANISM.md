# From curve parameters to a mechanism

Handover plan, written 2026-09-14. It replaces `PLAN_STEP3_REVISED.md` (removed;
see git history) and everything after stage 3.1 of `PLAN_MECHANISM_NEXT_STEPS.md`.

> **AMENDMENT 2 (2026-09-15). Task 5a is done (commit 5b3d0ac). Resume at
> section 11b (Task 5b), then continue to Task 6.** Amendment 2 adds link
> power, a one-buffer rule for the links, an exact definition of Stage B's
> candidates and a realistic planted recovery for Stage B. It overrides the
> items of Tasks 5-7 it names.
>
> **AMENDMENT 1 (2026-09-15), done.** Section 11a changed the tie rule, the
> standard errors, the temperature floor and `within_run_check`, reran Tasks 4
> and 5 into `data/fits/rate_laws/v2`, and overrides the items of Tasks 3-6 it
> names. Those overrides still apply.

## The goal, in one paragraph

Find the set of equations that best PREDICTS the progress curves, where "best"
means best at predicting runs the fit did not see. Get there in three stages,
cheap to expensive:

- **Stage A — rate laws for the curve parameters.** Every catalysed curve
  already carries a fitted activation-sink form with four parameters. Write
  each parameter as a function of the conditions ([S], [H2O2], [HOO-], [buf],
  [enz], buffer identity, temperature), and find which functions predict held-
  out runs. Archive-wide, no ODEs, fast.
- **Stage B — one global analytic fit.** Fit the activation-sink form to the
  readings of all those curves at once, with its parameters replaced by
  Stage A's rate laws. Still no ODEs.
- **Stage C — the chemical mechanism.** Only for what the analytic form cannot
  represent. **Not in this plan.** It is written after STOP 2, from Stage B's
  result.

## Why this plan is written the way it is

Earlier rounds went wrong where the plan left a judgement to the executor:
- a reading table became a mechanistic verdict;
- a test passed because it started at the truth;
- an F value was read without checking the two fits were comparable;
- a constant's change was read without checking whether another constant had
  absorbed it.

So in this plan:
- **Verdicts are strings returned by functions**, with fixed thresholds, pinned
  by tests. Reports quote them and add no verdict of their own.
- **Every task has the same parts:** Goal, Changes, Tests, Run, Report, Gate.
- **Anchors** are numbers you must reproduce before trusting new code.
- **Assumptions** (section 5) each have a trigger. When a trigger fires, stop
  and report. Do not adapt the plan.

---

## Contents

0. Rules for the executor
1. What is already true
2. Glossary
3. The data and the exclusions
4. What the data can and cannot identify
5. Assumptions and their stop triggers
6. Tasks and stops
7. Task 1 — the reference-design parser
8. Task 2 — the curve-parameter tables
9. Task 3 — the rate-law fitter, its health checks and planted tests
10. Task 4 — the rate-law search
11. Task 5 — links between the parameters (halted; see 11a)
11a. Amendment 1 — Task 5a: the tie rule, clustered errors, the rerun — STOP 1 (done)
11b. Amendment 2 — Task 5b: link power, the one-buffer rule, Stage B's candidates and its realistic planting
12. Task 6 — the global analytic fitter and its pilot
13. Task 7 — the global fits — STOP 2
14. Report templates

---

## 0. Rules for the executor

**Repository rules** (from `CLAUDE.md`):
- `.venv/bin/python`, from the repository root. Never `python`, never from
  inside a folder.
- Before adding any top-level name (function, class, constant, private or
  not), `grep -rn "NAME" --include=*.py .` must find nothing. The duplicate
  guard in `data/test_curve_metrics.py` fails otherwise.
- `.venv/bin/python run_gates.py` must print `0 failed` before every commit.
- Commit straight to master. Subject `area: what changed`, one paragraph, then
  exactly these two lines:
  ```
  Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01Vht919cHDPUdJWghusZFdS
  ```
- Every task that produces a result gets a `DATA_VERIFICATION.md` entry at the
  top of the file, using the template in section 14.
- Never edit `FITTING.md`, `MECHANISM.md`, `CLAUDE.md`, `BUBBLES.md` or
  `.claude/skills/`. Propose text in a stop message instead.
- Never overwrite a file in `data/fits/`. This plan saves under
  `data/fits/rate_laws/`.
- Parallel worker processes set `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`
  and `MKL_NUM_THREADS=1` in the environment BEFORE numpy is imported.
- Long jobs run in the background, save each fit as it finishes, and skip any
  fit whose save already exists.
- A number that goes into a document comes from a function in the repository,
  never from a throwaway script.

**Rules that exist because they were broken before:**
1. **Tests.** Never loosen an assertion, threshold or tolerance written in this
   plan. Never start a recovery fit at the planted truth. Build planted designs
   exactly as written; never slice or subsample them. If a test fails, stop and
   report the failing values.
2. **Comparisons.** Never compare two fits, or read a coefficient's change
   between two fits, without the health check of section 9 on both.
3. **Coefficients.** Never quote a coefficient that its health check flags.
   Write "not identified (<flag>)".
4. **Words.** Use the glossary. These are forbidden unless a function in this
   plan returns them as a verdict string:
   - that a term "is required", "is established", "is the mechanism", "earns";
   - that a result "confirms", "rescues" or "rules out" anything;
   - that a correction or difference is "small", "negligible" or "does not
     matter";
   - that the data "prefer" a species (H2O2 or HOO-, a buffer form).
5. **Stops.** At a STOP: commit, post the stop message (section 14), and do
   nothing more until the user replies.

## 1. What is already true

**1.1 The activation-sink form is a small mechanism.** Every live curve in
the archive has a fit (`summary_kinetics.fit_activation_sink`) of

    dP/dt = v_act + (v0 − v_act)·e^(−t/τ) − k·P

on its gas-corrected readings. That is exactly: product formed at a rate that
relaxes from v0 to v_act with time constant τ, and lost by a first-order sink
with rate constant k. Its solution is analytic:
`P(t) = c + v0·h(t) + v_act·g(t)`, where `(h, g) =
summary_kinetics._activation_sink_columns(times, τ, k)`.

**1.2 Where the fits live.** `scope.fits(scope.archive())` returns
`{(experiment, sample): CurveFit}`. Each `CurveFit` has `times`, `corrected`
(the gas-corrected readings), `noise`, `events` and `activation_sink`, an
`ActivationSinkFit` with:
- estimates `c`, `v0`, `v_act`, `tau`, `k`, `acceleration`;
- profile intervals `v0_interval`, `v_act_interval`, `tau_interval`,
  `k_interval` (each `(low, high)`);
- flags `v_act_resolved`, `v0_resolved`, `tau_resolved`, `k_state`
  ("none", "resolved", "unbounded", "unresolved"), `sink_earned`, `kind`.

`scope.frame(scope.archive())` is the same fits flattened, one row per curve
(`v_act_corrected`, `tau_act_corrected`, `k_sink_corrected`, their resolved
flags, `act_sink_kind_corrected`, and the conditions `substrate`, `buffer`,
`temperature`, `pH`, `s0`, `h2o2`, `hoo`, `buf`, `e0`, `live`, `bubble_load`).
It takes about 22 s to build.

**1.3 A catalysed curve is a catalytic INCREMENT.** Its reference cuvette holds
the same mixture minus the enzyme (`DATA_VERIFICATION.md` 2026-08-29). So every
parameter in this plan describes the catalyst's contribution, with the
uncatalysed background already removed. Enzyme-free curves are not used in
Stages A and B.

**1.4 The two meanings of τ.** `kind` is "lag" when v0 < v_act and "burst"
when v0 > v_act. On a lag curve τ is how long the catalyst takes to activate.
On a burst curve it is how long an early fast rate takes to decay. These are
different processes, so τ is fitted separately for the two families, and
Task 5 tests whether they share one rate law. Lag and burst curves occur in
the same run (see anchors), so the family is a per-curve label.

**1.5 What the empirical analyses already found** (the source of the candidate
terms; none is assumed true here):
- the rate saturates in [S] (`saturation.michaelis_by_element`);
- the rate saturates in peroxide, about 77% bound at 82.5 mM
  (`saturation.binding_by_element`);
- the rate has a pH order of about +0.55 in [HOO-] (`ph/`, `scope.ph_order`);
- the rate has a positive buffer order (`buffer/`, `saturation.activation_by_buffer`);
- the activation clock does not respond to peroxide
  (`saturation.binding_by_element`, `k_act_corrected`);
- the barrier sits in the saturated rate, not in Km
  (`saturation.michaelis_temperature`);
- 4OMe's catalysed rate falls linearly in the product made (`product_fate/`).

**1.6 Open bug, fixed in Task 1.** `verify_enzyme.reference_design` reads the
`Sum:` and `Sum*9:` rows beneath the two-axis block's cuvette tables as
cuvettes, and classifies exps 135-151 as "other". Their `Ref.` rows have `Enz`
0.000: the reference omits only the enzyme.

## 2. Glossary

### Curve parameters

| name | meaning | units |
|---|---|---|
| `v_act` | the production rate the curve relaxes onto (lag: the activated rate; burst: the late rate) | AU/s |
| `v0` | the production rate at t = 0 | AU/s |
| `k_act` | 1/τ. Lag family: the activation rate constant. Burst family: the burst's decay rate constant | 1/s |
| `k_sink` | first-order product loss rate constant `k` | 1/s |
| family | "lag" or "burst", from `kind` | — |

Rates are in AU/s. The rate laws are fitted per substrate, so a substrate's
extinction coefficient is absorbed into the intercepts.

### Rate-law terms

Each rate law is `log(parameter) = intercept[buffer] + Σ terms`. The terms:

| term | contribution to the log | coefficients |
|---|---|---|
| `S:power` | `a_S · log[S]` | order `a_S` |
| `S:mm` | `log([S]/(Km + [S]))` | `Km` (mM) |
| `H:power` | `a_H · log[H2O2]` | order |
| `H:bind` | `log(K_H·[H2O2]/(1 + K_H·[H2O2]))` | `K_H` (1/mM) |
| `H:relax` | `log(1 + K_H·[H2O2])` (a pre-equilibrium clock) | `K_H` |
| `HOO:power` | `a_HOO · log[HOO-]` | order |
| `BUF:power` | `a_B · log[buf]` | order |
| `BUF:bind` | `log(K_B·[buf]/(1 + K_B·[buf]))` | `K_B` (1/mM) |
| `BUF:relax` | `log(1 + K_B·[buf])` | `K_B` |
| `BUF:power_by_buffer` | a separate `a_B` per buffer identity | orders |
| `E:power` | `a_E · log[enz]` | order |
| `T:arrhenius` | `−(Ea/R)·(1/T − 1/298.15 K)` | `Ea` (kJ/mol) |

A "family" is S, H, HOO, BUF, E or T. A model picks at most one option per
family, or none.

### Verdict strings

Only these may describe a term, and only as a function returns them:
- `"<family>: <option> in every tied model"`
- `"<family>: a dependence in every tied model, option undecided (<options>)"`
- `"<family>: absent from every tied model"`
- `"<family>: undecided"`
- `"<family>: not identifiable on this table (<reason>)"`

## 3. The data and the exclusions

**The table rows**: every row of `scope.frame(scope.archive())` with `live`
true and `e0 > 0`, minus the curves `activation_sink.remaining_curves()` lists
(the form cannot hold them; their `catalysed` column is true).

**Per element**, a curve is in the element's table if:
- `v_act`: `v_act_resolved` and `v_act > 0` and `v_act_interval[0] > 0`;
- `k_act`: `tau_resolved` and `tau_interval[0] > 0`;
- `k_sink`: `k_state == "resolved"` and `k_interval[0] > 0`.

Every excluded curve is COUNTED by reason in the report.

**Sensitivity cuts** (Tasks 4, 7), defined exactly:
- `S-gas`: drop curves with `bubble_load > 1`.
- `S-weak`: drop the two-axis block's weak runs, the experiments in
  `fit_dataset.TWO_AXIS_BLOCK` not in `scope.strong_runs()`.
- `S-pyro`: drop every curve with `buffer == "Pyrophosphate"`. The pyrophosphate
  two-axis block is most of BnOH, so this checks it does not decide BnOH alone.

## 4. What the data can and cannot identify

Checked in code (Task 3, `term_identifiable`). Know these before reading any
result:

- **pH is one value per run everywhere.** `HOO:power` is identified only
  between runs, where day, enzyme batch and run length also differ. It is
  never a within-run axis: [HOO-] shifts slightly with [buf] inside a run.
- **Temperature moves only in 4OMe phosphate** (15-40 °C, one run per
  temperature). `T:arrhenius` is dropped from BnOH tables.
- **[enz] moves only between runs**: 4OMe 0.032-0.273 mM, BnOH 0.014-0.28 mM.
- **Buffer identity and pH are confounded.** Every boric run in these tables
  sits at pH 8.46 or above (to 10.34). The
  buffer intercepts absorb part of any pH effect, and `HOO:power` absorbs
  part of any buffer-identity effect.
- **[S] and [buf] move together inside many 4OMe runs** (substrate volume
  displaced buffer volume). Where both families are in a model, the code
  reports their within-run correlation.
- **BnOH is mostly the pyrophosphate two-axis block**, the only runs moving
  [S] and [H2O2] within a run. `S-pyro` checks it.
- **Run-level leave-out folds:** 34 (4OMe) and 31 (BnOH) before element
  exclusions. A run whose buffer has no other run in the table cannot be held
  out (its intercept would be unseen); the code skips that fold and counts it.

## 5. Assumptions and their stop triggers

| id | assumption | trigger → STOP and report |
|---|---|---|
| A1 | the anchors of Task 2 reproduce | any Task 2 anchor differs |
| A2 | Task 1 changes reference designs only inside the two-axis block | any other experiment's design changes |
| A3 | the rate-law fitter recovers a strong planted model exactly | Task 3's strong planted test fails |
| A4 | the global fitter recovers a strong planted model exactly | Task 6's planted test fails |
| A5 | each stage is affordable | Task 4's or Task 6's pilot projects more than 12 hours of wall time |
| A6 | Stage B's best model fits | its health check still flags "not converged" after one refit with doubled starts |
| A7 | Stage A's coefficients are not an artefact of keeping only resolved curves | more than half of the coefficients shared by Stage A and Stage B shift by more than 2 standard errors |
| A8 | Amendment 1's rerun reproduces the re-scored tie sets | any tie size or verdict change in 11a's A8 table differs |
| A9 | Amendment 1's clustered barrier errors reproduce | any value in 11a's A9 list differs by more than 0.2, or a listed verdict differs |
| A10 | Amendment 2's link verdicts, link power and Stage B candidates reproduce | any item in 11b's A10 list differs |

An element with fewer than 15 curves after exclusions is not fitted in Tasks
4-5. Report it as "too few curves"; that is not a stop. Today this applies to
exactly one table, 4OMe-BnOH `k_act_burst` (11 curves).

**The burst-clock fallback.** When a substrate's `k_act_burst` table is too
small, Task 6 gives its burst curves the `k_act_lag` law plus one free extra
intercept for the burst family (`burst_offset`, fitted in Stage B). Task 5's
`shared_clock_law` is then not run for that substrate; report it as "too few
curves".

## 6. Tasks and stops

| task | content | ends with |
|---|---|---|
| 1 | reference-design parser fix | commit |
| 2 | `data/rate_laws.py`: the curve-parameter tables, anchors | commit |
| 3 | rate-law fitter, health, identifiability, planted tests | commit |
| 4 | rate-law search per element and substrate, cuts | commit |
| 5 | links: shared binding constant, lag/burst clocks, barriers | halted on its own test; continued by 5a |
| 5a | Amendment 1: two-test tie rule, run-clustered errors, temperature floor, within-run families; rerun Tasks 4-5 into `v2` | **STOP 1** (done, 5b3d0ac) |
| 5b | Amendment 2: one-buffer rule and link power; Stage B candidates defined; realistic planted recovery for Stage B | commit, then straight on to Task 6 |
| 6 | global analytic fitter, planted test, pilot | commit (STOP if A4 or A5) |
| 7 | global fits, comparison with Stage A, amplitude check | **STOP 2** |

---

## 7. Task 1 — the reference-design parser

**Goal.** `verify_enzyme.reference_design` reads the two-axis block's sheets as
"enzyme" (section 1.6).

**Changes.** In `reference_design`, add a keyword `stop_at_sum=True`. When
true, stop collecting cuvette rows at the first row whose first-column value,
as a lower-cased stripped string, starts with `sum`. Change nothing else.

**Tests** (`data/test_fit_ladder.py`):
- `test_the_two_axis_sheets_read_as_enzyme_references`: exps 135, 140 and 151
  classify as "enzyme".
- `test_the_sum_rows_change_only_the_two_axis_designs`: for every experiment
  in the manifest, `reference_design(sheet, stop_at_sum=False)` equals
  `reference_design(sheet)` unless the experiment is in
  `fit_dataset.TWO_AXIS_BLOCK`.

**Run.** `.venv/bin/python data/verify_enzyme.py`. Record the design counts
before (`git stash` the change, or call with `stop_at_sum=False`) and after.

**Report.** Template E. Table: design → count before → count after.

**Gate.** A2 → STOP. Otherwise commit.

## 8. Task 2 — the curve-parameter tables

**Goal.** One function that builds every table Stage A fits, with every
exclusion counted.

**Changes.** New module `data/rate_laws.py`, with its gate
`data/test_rate_laws.py`. Check the module name and every top-level name are
unused.

`curve_parameters(substrate) -> dict` returns:
- `rows`: one row per catalysed curve of that substrate (section 3's table rows),
  with `experiment`, `sample`, `buffer`, `temperature` (°C), `kelvin`, `pH`,
  `s0`, `h2o2`, `hoo`, `buf`, `e0`, `bubble_load`, `family` ("lag"/"burst"
  from `kind`), and, from the `ActivationSinkFit`, `v_act`, `v0`, `tau`, `k`,
  their intervals, `v_act_resolved`, `tau_resolved`, `k_state`, `sink_earned`.
- `elements`: a dict of DataFrames `"v_act"`, `"k_act_lag"`, `"k_act_burst"`,
  `"k_sink"`. Each holds the rows admitted by section 3, plus
  - `y`: `log(v_act)`, `log(1/tau)` or `log(k)`;
  - `se`: the log-scale standard error from the interval,
    `(log(high) − log(low)) / (2 × 1.96)`. For `k_act`, use `tau_interval`
    (the log-width is the same for 1/τ);
  - `weight`: `1 / (se² + SE_FLOOR²)`, with module constant
    `SE_FLOOR = 0.10`. The floor stops a few very tight curves from deciding
    a table.
- `excluded`: a DataFrame of (element, reason, count). Reasons: "form cannot
  hold the curve", "not resolved", "non-positive estimate", "non-positive
  interval end".
- `remaining`: the `remaining_curves()` rows for this substrate that were
  removed.

**Tests** (`data/test_rate_laws.py`) — the anchors, exact:
- `test_the_archive_anchors`, over both substrates together: 311 live
  catalysed curves in the frame; 94 of them removed as "form cannot hold the
  curve"; 217 remain.
- `test_the_substrate_anchors`:
  - 4OMe-BnOH: 105 curves over 34 runs; `v_act` resolved on 73, `tau` on 81,
    `k_state == "resolved"` on 44; `bubble_load > 1` on 1.
  - BnOH: 112 curves over 31 runs; `v_act` resolved on 65, `tau` on 74,
    `k_state == "resolved"` on 27; `bubble_load > 1` on 8.
  These are counts BEFORE the positivity exclusions.
- `test_the_element_table_anchors`, AFTER every exclusion of section 3
  (curves / runs):

  | element | 4OMe-BnOH | BnOH |
  |---|---|---|
  | `v_act` | 73 / 32 | 65 / 29 |
  | `k_act_lag` | 70 / 29 | 28 / 18 |
  | `k_act_burst` | 11 / 5 | 46 / 24 |
  | `k_sink` | 44 / 21 | 27 / 18 |

  and the curves whose own `v_act <= 0`: 4 (4OMe-BnOH), 10 (BnOH).
  4OMe-BnOH's `k_act_burst` table is below 15 curves, so Tasks 4-5 do not fit
  it (section 5), and Task 6 uses the fallback written there.
- `test_lag_and_burst_share_runs`: the 4OMe runs containing both families
  include 16, 21, 42, 46, 47, 48 and 130; the BnOH runs include 51, 66, 71,
  75, 76 and 83.
- `test_se_is_the_interval_width`: on a hand-built row with interval
  (e^−1, e^1), `se` = 2/(2 × 1.96).

**Report.** Template E. Per substrate and element: curves, runs, buffers and
their curve counts, and the `excluded` table.

**Gate.** A1 → STOP. Otherwise commit.

## 9. Task 3 — the rate-law fitter, its health checks and planted tests

**Goal.** A fitter for one model on one element table, with the checks that
were skipped in earlier rounds built in.

**Changes** (`data/rate_laws.py`). **Amendment 1 (section 11a) overrides
items 3, 5 and 6 below**: `fit_model` also returns run-clustered errors,
`T:arrhenius` needs 4 temperatures, and `within_run_check` always drops HOO, E
and T.

1. `TERM_OPTIONS`: a dict per element of the options each family may take:
   - `v_act` and `k_sink`: S {power, mm}; H {power, bind}; HOO {power};
     BUF {power, bind, power_by_buffer}; E {power}; T {arrhenius}.
   - `k_act_lag` and `k_act_burst`: S {power}; H {power, relax};
     HOO {power}; BUF {power, relax, power_by_buffer}; E {power};
     T {arrhenius}.
   "None" is always allowed and is not listed.

2. `design(table, model, nonlinear)`: the design matrix for a model (a dict
   family → option) at given nonlinear constants (`Km`, `K_H`, `K_B`, in
   log10). One intercept column per buffer present. Returns the matrix and
   the column names.

3. `fit_model(table, model, starts=(-2.0, 0.0, 2.0))`. Weighted least squares
   by variable projection:
   - the linear coefficients (intercepts, orders, `Ea/R`) come from a weighted
     `lstsq` at given nonlinear constants;
   - each nonlinear constant is optimised in log10 by
     `scipy.optimize.least_squares` within bounds (−4, 4). Try every
     combination of `starts` and keep the lowest cost.

   Returns a dict: `model`, `coefficients` (name → value), `stderr` (from the
   weighted normal equations at the optimum, scaled by the residual variance),
   `cost` (weighted sum of squared residuals), `curves`, `runs`,
   `nonlinear_at_bound` (names), `correlation` (matrix of the coefficients),
   `health`.

4. `rate_law_health(result, table)` → `flags` (list of strings) and `ok`.
   Flags, exactly:
   - `"at bound: <K>"` for a nonlinear constant within 0.05 of either bound;
   - `"collinear: <a>/<b> r=±0.xx"` for any coefficient pair with
     |correlation| > 0.95, intercepts excluded;
   - `"rank deficient"` if the design's rank is below its column count;
   - `"carried by fewer than 3 runs: <term>"` if a term's regressor, after
     removing each buffer's mean, is non-zero in fewer than 3 runs;
   - `"too few curves"` if `curves < 15`.

5. `term_identifiable(table, family, option) -> (bool, reason)`. False if:
   - the regressor takes fewer than 2 distinct values;
   - the SD of the regressor after removing each buffer's mean is below
     0.05 in log units (T: 1e-5 in 1/K);
   - the option is `power_by_buffer` and fewer than 2 buffers each have
     2 runs with within-buffer variation.
   The reason string names which.

6. `within_run_check(table, model)`: refit the model with one intercept per
   RUN instead of per buffer, dropping every family `term_identifiable` flags
   on the demeaned-by-run regressor (always HOO, E and T). Returns, per shared
   coefficient, both estimates and
   `disagree = |a − b| > 2·sqrt(se_a² + se_b²)`. Also returns, for S and BUF
   when both are in the model, the within-run correlation of their regressors
   (each demeaned by run).

**Tests** (`data/test_rate_laws.py`):
- `test_a_strong_planted_model_is_recovered` — **the machinery test (A3)**.
  - Rows: the real 4OMe-BnOH `v_act` element table, all of it, with its real
    conditions and `weight`.
  - Replace `y` with `intercept[buffer] + 0.8·log[S] + 0.5·log[HOO-] +
    log(0.02·[buf]/(1 + 0.02·[buf]))`, intercepts 0, −0.5 and 0.5 for Boric,
    Phosphate and Pyrophosphate, plus Gaussian noise of SD 0.01, seed 0.
  - Fit the TRUE model, with `starts=(-2.0, 0.0, 2.0)`.
  - Assert: `a_S` within 0.03 of 0.8; `a_HOO` within 0.03 of 0.5; log10 `K_B`
    within 0.1 of log10 0.02; health flags contain no "at bound" and no
    "rank deficient".
- `test_health_flags`: one hand-built case per flag, each producing exactly
  that flag.
- `test_identifiability_rules`: a table with one buffer at constant [S] makes
  `S:power` unidentifiable, with the reason "fewer than 2 distinct values".
- `test_within_run_check_catches_a_between_run_confound`: planted rows where
  a per-run offset correlates with [S] across runs, but [S] has no effect
  within runs. The between-run `a_S` is non-zero, the within-run one near
  zero, and `disagree` is true.

**Run.** Tests only.

**Report.** Template E. The planted recovery values.

**Gate.** A3 → STOP. Otherwise commit.

## 10. Task 4 — the rate-law search

**Goal.** For each substrate and element, every model allowed by
`TERM_OPTIONS` and `term_identifiable`, scored by leave-one-run-out.

**Changes** (`data/rate_laws.py`). **Amendment 1 (section 11a) overrides
items 3 and 5 below**: `tie_set` needs both the raw and the log test, and
`search` saves every model's fold scores into `data/fits/rate_laws/v2`.

1. `enumerate_models(table, element)`: every combination of at most one
   option per family, including none, keeping only options
   `term_identifiable` allows. Returns the models and a `dropped` table
   (family, option, reason).
2. `cross_validate(table, model)`: for each run, fit on the other runs and
   predict the held-out run's curves. Fold score is
   `Σ weight·(y − ŷ)²` over the held-out curves. Skip a fold whose held-out
   run has a buffer no other run in the table has, and count it. Returns
   per-fold scores and their sum.
3. `tie_set(scores)`: models whose per-fold differences `d` against the best
   satisfy `mean(d) <= 2·sd(d)/sqrt(folds)`, over the folds every model
   scored.
4. `term_verdicts(tie, dropped)`: the verdict strings of section 2, per
   family.
   - "an option in every tied model": every tied model has that option.
   - "a dependence in every tied model, option undecided": every tied model
     has an option, but not the same one.
   - "absent from every tied model": no tied model has any option.
   - "undecided": otherwise.
   - "not identifiable on this table": the family had no option left after
     `term_identifiable`.
5. `search(substrate, element, cut=None, workers=8)`: runs 1-4, applies
   `rate_law_health` and `within_run_check` to every tied model, and saves
   `data/fits/rate_laws/<substrate>_<element>_<cut or all>.json`. Skip if the
   save exists.

**Tests** (`data/test_rate_laws.py`):
- `test_every_run_is_held_out_once_and_never_trained_on`.
- `test_the_tie_rule` (hand-built fold scores: one tied, one not).
- `test_term_verdicts` (a hand-built tie set for each of the five verdicts).
- `test_the_search_finds_a_realistic_planted_model`, recorded and not asserted:
  - Rows: the real 4OMe-BnOH `v_act` table.
  - Truth: `S:mm` with Km = 2 mM, `HOO:power` 0.5, `BUF:bind` with
    K_B = 0.02 /mM, the Task 3 intercepts.
  - Noise: Gaussian with each row's own `sqrt(se² + SE_FLOOR²)`, seeds 0, 1
    and 2.
  - Run the search each time and write the three verdict tables to
    `data/fits/rate_laws/planted_identifiability.json`.
  - This is how well the real design identifies each family. Task 4's report
    quotes it beside every real verdict.

**Pilot.** Before the real search, time `cross_validate` on the largest model
of the largest table. Project the total:
`Σ over tables of models × folds × fold time / workers`. A5 → STOP.

**Run.** `search` for every substrate × element, cut `None`. Then, for every
table, the cuts `S-gas`, `S-weak` (BnOH only; the weak runs are all BnOH) and
`S-pyro`.

**Report.** Template E, one entry for the whole task. Per substrate × element:
- curves, runs, folds skipped;
- the ten best models with CV scores, and the tie-set size;
- `term_verdicts`, beside the planted identifiability for the same family;
- health flags of every tied model;
- `within_run_check` disagreements;
- for each cut, only the verdicts that change.

**Gate.** Commit.

## 11. Task 5 — links between the parameters — STOP 1

**Goal.** Three tests that couple elements. Mechanisms predict couplings that
single-element rate laws cannot see.

**Amendment 1 (section 11a) overrides items 1 and 3 below, and this task's
Report and STOP**: `shared_binding_law` returns "not identified (K at bound)"
when a K sits at a bound, `barriers` uses clustered errors and needs 4
temperatures, and STOP 1 happens at the end of Task 5a, not here.

**Changes** (`data/rate_laws.py`). Each test fits two elements jointly: the
summed weighted cost, one coefficient set per element, with the stated
coefficients shared. Each is scored by leave-one-run-out over the runs present
in either table.

1. `shared_binding(substrate, family)`, for family in BUF and H:
   - `v_act` takes the family's `bind` option, `k_act_lag` its `relax` option,
     and every other term is its own best tied model's.
   - Two variants: one shared K, and a separate K per element.
   - Verdict, exactly one of `"one K ties with separate K"` and `"separate K
     predicts better"` (the tie rule of Task 4).
   - Why: a pre-equilibrium in which the family's species activates the
     catalyst predicts one K for both.
2. `shared_clock_law(substrate)`:
   - `k_act_lag` and `k_act_burst` share every term coefficient, with a
     separate intercept per family; against fully separate laws.
   - Verdict, exactly one of `"one clock law ties with separate laws"` and
     `"separate clock laws predict better"`.
3. `barriers(substrate="4OMe-BnOH")`:
   - `Ea` with its standard error, from each element's best tied model
     containing `T:arrhenius`.
   - Verdict, exactly one of `"barriers equal within 2 standard errors"` and
     `"barriers differ by more than 2 standard errors: <pairs>"`, over the
     elements whose `Ea` is not flagged.

**Tests** (`data/test_rate_laws.py`):
- `test_shared_binding_on_a_planted_pre_equilibrium`:
  - Rows: the real 4OMe-BnOH `v_act` and `k_act_lag` tables.
  - Plant `v_act ∝ K·buf/(1 + K·buf)` and `k_act ∝ (1 + K·buf)` with
    K = 0.03 /mM in both, noise SD 0.01, seed 0.
  - Assert the verdict "one K ties with separate K" and the shared K within
    0.1 decade of 0.03.
- `test_separate_binding_is_detected`: the same with K = 0.003 and 0.3 /mM.
  Assert "separate K predicts better".

**Run.** All three per substrate (`barriers` on 4OMe only).

**Report.** Template E, then **STOP 1** (template S).

**Question at STOP 1:** "Stage A's rate laws and links are above. Proceed to
Stage B (Tasks 6-7), using as candidates each element's best tied model and its
simplest tied model?"

## 11a. Amendment 1 — Task 5a: the tie rule, clustered errors, the rerun — STOP 1

Added 2026-09-15. **This section overrides** Task 3 items 3, 5 and 6, Task 4
items 3 and 5, Task 5 items 1 and 3, and Task 6's standard errors. Where their
text differs from this section, this section wins.

### Why

- **The tie rule failed the plan's own test.** Task 5 halted on
  `test_separate_binding_is_detected`. The planted constants, K = 0.003 and
  0.3, were recovered separately as 0.0031 and 0.296. Yet the tie rule called
  one shared K tied: mean difference 31.20, sd 110.95, 33 folds, bar 38.63.
  The rule compares RAW fold scores, and a few runs at extreme conditions
  dominate their spread, so a model that is worse on nearly every run can
  still tie.
- **The standard errors treat the curves of one run as independent.** They
  share a day, a stock and a cuvette offset. Clustered by run, the 4OMe-BnOH
  `k_sink` activation energy is 54.7 ± 20.7 kJ/mol, not ± 11.8.
- **A temperature term was fitted on 3 temperatures.** That is the 4OMe-BnOH
  `k_sink` table.
- **`within_run_check` keeps [HOO-].** Task 3 item 6 described dropping it
  instead of requiring it. [HOO-] moves inside a run only through ionic
  strength, which follows [buf], so a within-run [HOO-] coefficient is a
  buffer coefficient.

**Keep Task 5's uncommitted code in the working tree and adapt it. Do not
revert it.** Do not delete or edit anything in `data/fits/rate_laws/`.

### Goal

Fix the four flaws, rerun Tasks 4 and 5 under the fixes into a new directory,
and stop at STOP 1 again.

### Changes (`data/rate_laws.py`; check every new name is unused first)

1. **`tie_set(scores)` keeps a model only if BOTH tests keep it.**
   - Folds: the columns with no NaN in any model, as now. Best: the lowest
     total, as now.
   - Raw test, as now: `d = score − score[best]`; kept if
     `mean(d) <= 2·sd(d)/sqrt(folds)`.
   - Log test: `d = log(max(score, TIE_LOG_FLOOR)) −
     log(max(score[best], TIE_LOG_FLOOR))`, with module constant
     `TIE_LOG_FLOOR = 1e-12`. Kept by the same inequality.
   - Tied = kept by both. The best is always tied. Return the models in
     `scores.index` order.
   - Add `tie_statistics(scores)`: one row per model, with `raw_mean`,
     `raw_bar`, `raw_kept`, `log_mean`, `log_bar`, `log_kept`, `tied`, and
     `worst_fold_log_ratio` = the largest `log(score/score[best])` over the
     folds (added 2026-09-15). `worst_fold_log_ratio` is a DIAGNOSTIC: it
     changes no tie set, verdict or anchor. Report it for every tied model.
     It exists because of the limit stated under Tests.
2. **`fit_model` adds run-clustered standard errors,**
   `result["stderr_clustered"]`:
   - `e = y − prediction` at the optimum (from `_predict`). `X` is the
     result's `_matrix`; `w` is the table's weights; `G` is the number of runs.
   - `B = pinv(Xᵀ·diag(w)·X)`. `M = Σ over runs g of s_g·s_gᵀ`, with
     `s_g = X_gᵀ·(w_g ∘ e_g)`. `V = G/(G − 1)·B·M·B`.
   - `stderr_clustered[name] = sqrt(V[i, i])` for every linear coefficient.
     `stderr_clustered["Ea"] = stderr_clustered["Ea_R"] × 8.314462618/1000`.
     Nonlinear constants get NaN.
   - `stderr` (the naive errors) stays exactly as it is.
   - **From now on, every verdict, comparison and quoted ± in this plan uses
     `stderr_clustered`.**
   - `rate_law_health` adds the flag
     `"fewer than 10 runs: clustered errors unreliable"` when `G < 10`.
3. **`term_identifiable`:** `T:arrhenius` also needs at least 4 distinct
   temperatures in the table. Reason string: `"fewer than 4 temperatures"`.
4. **`within_run_check`:**
   - Families HOO, E and T are ALWAYS dropped from the run design, whatever
     their within-run variation. List each in `dropped` with the reason
     `"between-run family"`.
   - `disagree` uses `stderr_clustered` for the between-run estimate. For the
     within-run estimate, use the same clustered formula applied to the run
     design.
5. **`shared_binding_law`:** if the shared K or either separate K is within
   0.05 of a log10 bound, the verdict is exactly
   `"not identified (K at bound)"`. Check this BEFORE the tie rule.
6. **`barriers`:**
   - Use `stderr_clustered`.
   - An element whose table has fewer than 4 distinct temperatures is left out
     and listed as `"not identified (fewer than 4 temperatures)"`.
   - The verdict compares only the remaining elements. If fewer than 2 remain,
     the verdict is `"fewer than 2 elements with a barrier"`.
7. **Saves:**
   - `RATE_LAW_DIR` becomes `data/fits/rate_laws/v2`.
   - `search` saves `fold_scores` for EVERY model, not only tied ones, and
     `tie_statistics`.
   - `_best_tied_model` reads `v2`.
8. **For Task 6, when you get there:** `global_fit` also reports
   `stderr_clustered`, computed the same way on the stacked problem:
   - `J` is `least_squares`'s Jacobian of the weighted residuals with respect
     to the global coefficients, and `r` the weighted residuals, both at the
     optimum.
   - `B = pinv(JᵀJ)`, `M = Σ over runs g of (J_gᵀ·r_g)(J_gᵀ·r_g)ᵀ`,
     `V = G/(G − 1)·B·M·B`.
   - Treat each curve's `c` and `v0` as fixed.
   - `stage_shift` uses `stderr_clustered` on both sides.

### Tests (`data/test_rate_laws.py`)

Every existing test stays unchanged. These three must pass under the new
rule: `test_the_tie_rule`,
`test_shared_binding_on_a_planted_pre_equilibrium` and
`test_separate_binding_is_detected`. Add:

- **`test_the_tie_rule_needs_both_tests`.** Hand-built, 30 folds, base fold
  scores `b_f = 10**(f/5)` for f = 0..29, so the folds span six decades.
  Models:
  - `"best"`: `b_f`.
  - `"consistent"`: `1.10·b_f`. Assert NOT tied (the log test catches it).
  - `"concentrated"`: `b_f`, plus `40·max(b)` added to folds 0 through 5
    (6 of 30). Assert NOT tied, and assert `raw_kept` is False for it (the raw
    test catches it: raw t = 2.69, log t = 2.66). **Corrected 2026-09-15 at the
    executor's stop:** the first version added the outlier to folds 0 and 1
    only, which no paired test can catch (see the next test).
  - `"equal"`: `b_f·(1 + 0.02·(−1)**f)`. Assert tied; assert `"best"` is also
    tied.
- **`test_the_tie_rule_cannot_see_three_fold_failures`** (added 2026-09-15).
  The same `b_f` over 30 folds, with a model `"hidden"` = `b_f` plus
  `40·max(b)` on folds 0, 1 and 2 only. Assert `"hidden"` IS tied, and assert
  its `worst_fold_log_ratio` is above 10.
  - **Why this test exists.** A model that differs from the best on k of n
    folds has paired t = sqrt(k·(n − 1)/(n − k)), whatever the size of the
    difference, because the difference inflates the spread by the same factor.
  - So the rule cannot exclude a model that fails on k <= 4n/(n + 3) folds,
    however badly: 3 of 30, or 2 of 18.
  - The test pins that limit, and `worst_fold_log_ratio` is how a report shows
    where it bites.
- **`test_run_clustered_errors_cover_a_between_run_null`.** 200 planted tables,
  seeds 0..199. Each has 12 runs × 4 curves and one buffer. Per run, draw
  `x_run ~ N(0, 1)`; set `hoo = exp(x_run)` on every curve of that run; set
  every other condition column to 1.
  - `y = run_offset + noise`, with `run_offset ~ N(0, 0.3)` per run and
    `noise ~ N(0, 0.02)` per curve. `weight = 1`. The true HOO order is 0.
  - Fit the model `{"HOO": "power"}`.
  - Assert: the share of seeds with `|a_HOO / stderr| > 2` is above 0.22, and
    the share with `|a_HOO / stderr_clustered| > 2` is below 0.17.
- **`test_temperature_needs_four_temperatures`.** A table whose temperatures
  are 25, 30 and 35 °C makes `T:arrhenius` unidentifiable, with the reason
  `"fewer than 4 temperatures"`.
- **`test_within_run_check_never_keeps_between_run_families`.** A planted table
  where `hoo` varies by 20% inside every run. The within-run design has no HOO
  coefficient, and `dropped` lists it with `"between-run family"`.
- **`test_shared_binding_at_bound_is_not_identified`.**
  - Rows: the real 4OMe-BnOH `v_act` and `k_act_lag` tables.
  - Plant `v_act` with no buffer dependence (its intercepts plus noise SD
    0.01, seed 0), and `k_act` proportional to [buf] (noise SD 0.01).
  - Assert the verdict `"not identified (K at bound)"`.

### Run, in order

1. `.venv/bin/python data/test_rate_laws.py`. Every test passes.
2. `search` for every substrate × element × cut that Task 4 ran (the same 25),
   into `v2`. Then the realistic planted identifiability search, into `v2`.
3. **Anchors A8** — the number of tied models in the `v2` saves must be
   exactly:

   | save | tied |
   |---|---|
   | 4OMe-BnOH `v_act`: all / S-gas / S-pyro | 32 / 32 / 13 |
   | 4OMe-BnOH `k_act_lag`: all / S-gas / S-pyro | 24 / 24 / 10 |
   | BnOH `v_act`: all / S-gas / S-weak / S-pyro | 33 / 27 / 78 / 89 |
   | BnOH `k_act_lag`: all / S-gas / S-weak | 4 / 10 / 24 |
   | BnOH `k_act_burst`: all / S-gas / S-weak / S-pyro | 70 / 72 / 68 / 72 |
   | BnOH `k_sink`: all / S-gas / S-weak | 104 / 104 / 100 |

   and these verdicts must change from the old saves exactly as follows, with
   no other change in these tables:
   - 4OMe-BnOH `k_act_lag` (all, S-gas, S-pyro): T becomes
     `"T: arrhenius in every tied model"`.
   - 4OMe-BnOH `v_act` (all, S-gas, S-pyro): S becomes
     `"S: a dependence in every tied model, option undecided (mm, power)"`.
   - BnOH `k_act_lag` (all): H becomes
     `"H: a dependence in every tied model, option undecided (power, relax)"`.

   4OMe-BnOH `k_sink` is not anchored: change 3 removes its temperature term.
4. Task 5's three links, under `v2`.
5. **Anchors A9** — `barriers("4OMe-BnOH")` must give:
   - `v_act`: Ea 90.4, clustered se 7.9;
   - `k_act_lag`: Ea 84.9, clustered se 9.8 (each value ±0.2);
   - `k_sink`: `"not identified (fewer than 4 temperatures)"`;
   - verdict: `"barriers equal within 2 standard errors"`.

### Report

Two `DATA_VERIFICATION.md` entries, template E:

1. **"Amendment 1: the tie rule, clustered errors, the temperature floor and
   the within-run families."** Contents:
   - the Why above, with its numbers;
   - the changes;
   - one table per save: tied models in the old save → in `v2`, and every
     verdict that changed;
   - the A8 and A9 anchors, reproduced;
   - a sentence stating that the tie sizes and verdicts of the 2026-09-14
     "rate-law search" entry are superseded by this one.
2. **"Task 5: links between the parameters"**, with the `v2` results and
   clustered errors.

### Gate

A8 or A9 → STOP and report. Otherwise:
1. `run_gates.py` prints `0 failed`.
2. Commit: Task 5's code, the amendment's code, the tests, the `v2` saves and
   both entries.
3. **STOP 1**, template S.

**Question at STOP 1** (replaces Task 5's question): "Stage A under Amendment 1
is above. Proceed to Stage B (Tasks 6-7), with each element's best tied model
and simplest tied model from the `v2` saves as candidates?"

## 11b. Amendment 2 — Task 5b: link power, the one-buffer rule, Stage B's candidates and its realistic planting

Added 2026-09-15, after STOP 1. The user approved Stage B with these additions.

**This section overrides** the parts of the tasks named below, wherever their
text differs:
- Task 5 item 1: the order of the verdicts.
- Task 6: the candidate set, the tests and the pilot.
- Task 7: the Run and the Report.

### Why

- **Stage A's link tests were planted at noise 0.01.** Real rows carry
  sqrt(se² + 0.10²). A review re-ran `shared_binding_law` on 4OMe-BnOH BUF at
  realistic noise:
  - K = 0.003 (v_act) and 0.3 (k_act_lag): "separate K predicts better" in 3
    of 3 seeds. The raw t was about 1.6; the log t was 3.1 to 4.7.
  - K = 0.03 in both: "one K ties with separate K" in 3 of 3 seeds.
  - K = 0.0999 and 0.00439, the real separate estimates: caught in 4 of 5
    seeds; seed 2 tied.

  So the real "one K ties" on 4OMe-BnOH BUF is informative, with about a 1 in
  5 chance of a tie even if the two Ks really were that far apart. That result
  exists only in a scratch script; it must become a recorded function.
- **The 4OMe-BnOH H link rests on one buffer.** Only 4 pyrophosphate runs (exps
  127, 129, 130 and 131) have a [H2O2] other than 82.5 mM. A constant
  identified inside one buffer's runs cannot be separated from that buffer.
- **Stage B has only a strong planted test** (noise 0.1 × each curve's noise).
  Stage A's realistic planting showed what a strong planting hides. At real
  noise it recovered:
  - the HOO order in 3 of 3 seeds;
  - an S dependence in 2 of 3;
  - a buffer dependence in 1 of 3;
  - the option (mm or power, bind or power) in none.
- **"Simplest tied model" was never defined precisely.**

### Goal

Record the link power, add the one-buffer rule, define Stage B's candidates
exactly, and add a realistic planted recovery to Stage B. Then continue to
Task 6 without stopping.

### Changes (`data/rate_laws.py`; check every new name is unused first)

1. **`shared_binding_law` checks three things, in this order:**
   1. K at a bound → `"not identified (K at bound)"` (as Amendment 1).
   2. **One buffer** → `"not identified (one buffer)"`. In each of the two
      tables, a buffer "carries" the family if its own rows have at least 2
      curves and the log of the family's concentration ([buf] for BUF,
      [H2O2] for H) has SD ≥ 0.05 after removing that buffer's mean. If either
      table has fewer than 2 carrying buffers, return this verdict.
   3. The tie rule.
2. **`link_power(substrate, family, k_rate, k_clock, seeds)`:**
   - **Rows:** the substrate's real `v_act` and `k_act_lag` element tables.
   - **For each seed:**
     - Make one generator, `numpy.random.default_rng(seed)`.
     - Build the `v_act` table first, then the `k_act_lag` table, in that
       order. For each: `y = intercept[buffer] + term +
       generator.normal(0.0, spread)`, where `spread` is the vector
       `sqrt(se² + SE_FLOOR²)` over the table's rows, drawn in one call.
     - Intercepts: Boric 0.0, Phosphate −0.5, Pyrophosphate 0.5.
     - Terms: `log(k_rate·x/(1 + k_rate·x))` for `v_act` and
       `log(1 + k_clock·x)` for `k_act_lag`, with x = [buf] for BUF and
       [H2O2] for H.
     - Call `shared_binding_law(substrate, family, rate_table=...,
       clock_table=...)`.
   - **Returns**, per seed: the verdict, the shared K, the separate Ks, and the
     raw t and log t of the shared-minus-separate fold scores, each
     `t = mean/(sd/sqrt(folds))`. Also `caught`, the number of seeds whose
     verdict is "separate K predicts better".
   - **Saves** `data/fits/rate_laws/v2/link_power.json`. Reads that file back if
     it exists.
3. **`stage_b_candidates(substrate)`**, from the `v2` `_all` saves, per element:
   - `best`: the lowest CV total among the tied models.
   - `simplest`: the fewest coefficients among the tied models.
     - Each family option counts 1; `power_by_buffer` counts the number of
       buffers in that element's table; intercepts do not count.
     - Ties are broken by the lower CV total.
     - The empty model (intercepts only) is allowed.
   - If `best` and `simplest` are the same model, the element has one
     candidate.
   - For a substrate without a fitted `k_act_burst` table, the burst candidate
     is `"lag law + burst_offset"`.
   - The combinations are every choice of one candidate per element, each also
     with `k_sink` None. These are exactly Task 6's candidates.
4. **`planted_global(substrate, seed=0)`**, recorded and not asserted:
   - **Truth:** each element's `best` candidate, with the coefficients
     `fit_model` gives on that element's full table.
     - Where the burst fallback applies, the truth `burst_offset` is the median
       over that substrate's `k_act_burst` element rows of
       `y − (lag-law prediction at their conditions)`.
   - **Curves:** exactly the curves `global_fit` uses for that substrate.
     - Each curve keeps its own family label.
     - Each curve's `c` and `v0` are its own `ActivationSinkFit.c` and `.v0`.
   - **Readings:** the model at each curve's times, plus `normal(0,
     curve.noise)` from `default_rng(seed)`, drawn curve by curve in the order
     `global_fit` iterates.
   - **Starts:** never the truth. Use the truth + 0.5 on every linear
     coefficient and +0.5 decade on every log10 K, plus Task 6's restarts.
   - **Run:** `global_fit` and `global_cross_validate` for EVERY combination
     from `stage_b_candidates`, the same set Task 7 fits.
   - **Record**, in `data/fits/rate_laws/v2/planted_global_<substrate>.json`:
     - for the true combination: per coefficient, the truth, the estimate,
       `stderr_clustered`, and `within_2se`;
     - its health flags;
     - the candidates' two-test tie set, with `tie_statistics` including
       `worst_fold_log_ratio`;
     - whether the true combination is the best, and whether it is tied.
5. **Task 7's Report.**
   - Beside every Stage B coefficient, give its planted twin's `within_2se`.
   - A Stage B coefficient whose planted twin has `within_2se` False is quoted
     as `"not identified (planted recovery failed)"`.
   - Beside the Stage B tie set, give whether the planted truth was tied.
6. **Task 6's pilot** projects the real fits AND `planted_global` together, at
   about twice the real fits. A5's 12-hour bar applies to that total.

### Tests (`data/test_rate_laws.py`)

Every existing test stays unchanged. Add:
- **`test_one_buffer_links_are_not_identified`:** `shared_binding_law(
  "4OMe-BnOH", "H")` on the real tables returns
  `"not identified (one buffer)"`.
- **`test_link_power_anchors`:** reads `link_power.json` (produced by Run
  step 3) and asserts the A10 link-power anchors below.
- **`test_stage_b_candidates_anchors`:** asserts the A10 candidate table below.
- **`test_planted_global_is_recorded`:** after Task 7, both
  `planted_global_<substrate>.json` files exist and contain every key listed in
  change 4. No value is asserted.

### Run, in order

1. `.venv/bin/python data/test_rate_laws.py`. Every test except
   `test_planted_global_is_recorded` passes; that one passes after Task 7.
2. Rerun `shared_binding_law` for both substrates, for BUF and for H.
3. `link_power("4OMe-BnOH", "BUF", ...)` three times:
   - `(0.003, 0.3, seeds=(0, 1, 2))`;
   - `(0.03, 0.03, seeds=(0, 1, 2))`;
   - `(0.0999, 0.00439, seeds=(0, 1, 2, 3, 4))`.
4. Continue to Task 6 with change 6, then Task 7 with changes 4 and 5. Run
   `planted_global` before the real global fits.

### Anchors A10

- **`shared_binding_law` verdicts:**

  | substrate, family | verdict |
  |---|---|
  | 4OMe-BnOH, BUF | "one K ties with separate K" |
  | 4OMe-BnOH, H | "not identified (one buffer)" |
  | BnOH, BUF | "not identified (K at bound)" |
  | BnOH, H | "one K ties with separate K" |

- **`link_power`:**
  - (0.003, 0.3): caught 3 of 3.
  - (0.03, 0.03): caught 0 of 3.
  - (0.0999, 0.00439): caught 4 of 5, the one tie at seed 2.
- **`stage_b_candidates`** (model ids as in the saves; `''` is the empty model):

  | substrate | element | best | simplest |
  |---|---|---|---|
  | 4OMe-BnOH | v_act | `S:power\|HOO:power\|BUF:bind\|T:arrhenius` | `S:power\|HOO:power\|T:arrhenius` |
  | 4OMe-BnOH | k_act_lag | `S:power\|E:power\|T:arrhenius` | `E:power\|T:arrhenius` |
  | 4OMe-BnOH | k_act_burst | `lag law + burst_offset` | same |
  | 4OMe-BnOH | k_sink | `S:power` | `''` |
  | BnOH | v_act | `S:mm\|H:power\|HOO:power\|BUF:power` | `HOO:power` |
  | BnOH | k_act_lag | `H:relax\|E:power` | `H:power` |
  | BnOH | k_act_burst | `S:power\|H:power\|HOO:power` | `''` |
  | BnOH | k_sink | `BUF:power` | `''` |

### Report

One `DATA_VERIFICATION.md` entry, template E, for Task 5b:
- the Why, with its numbers;
- the changes;
- the A10 anchors, reproduced;
- a link table giving each verdict beside its power, e.g. 4OMe-BnOH BUF: "one
  K ties with separate K" with "caught 4 of 5 at the observed separation";
- a sentence stating that the 2026-09-15 Task 5 entry's 4OMe-BnOH H verdict is
  superseded.

### Gate

A10 → STOP and report. Otherwise:
1. `run_gates.py` prints `0 failed`. `test_planted_global_is_recorded` is
   skipped until Task 7: make it print `"skipped until Task 7"` and return
   when the files are absent.
2. Commit.
3. Continue to Task 6. **There is no stop here.**

## 12. Task 6 — the global analytic fitter and its pilot

**Amendment 2 (section 11b) overrides this task's candidate set, tests and
pilot**: the candidates are exactly `stage_b_candidates(substrate)`, and the
pilot projects the real fits and `planted_global` together. Amendment 1
change 8 (`stderr_clustered` in `global_fit`, used by `stage_shift`) also
applies.

**Goal.** Fit the activation-sink form to the readings of many curves at once,
with v_act, τ and k given by rate laws, and with c and v0 free per curve.

**Changes** (`data/rate_laws.py`):

1. `global_fit(substrate, laws, curves=None, starts=None, restarts=4)`:
   - `laws`: one model per element (`"v_act"`, `"k_act_lag"`,
     `"k_act_burst"`, `"k_sink"`), or `None` for `"k_sink"`, meaning k = 0
     on every curve. For a substrate whose `k_act_burst` table was too small
     (section 5), `"k_act_burst"` is the string `"lag law + burst_offset"`:
     burst curves take the `k_act_lag` law plus one free global intercept.
     That substrate's candidate count is then at most 2 × 2 × 2 = 8.
   - Curves: section 3's table rows for the substrate, whether or not any
     element is resolved, minus curves whose own `v_act <= 0` (counted).
     Readings are `CurveFit.corrected` and `CurveFit.times`; each curve's
     residuals are divided by `noise·sqrt(n)`.
   - For each curve: `v_act`, `k_act` (by the curve's family) and `k` come
     from the laws at its conditions; τ = 1/k_act;
     `(h, g) = _activation_sink_columns(times, τ, k)`. Then `(c, v0)` is a
     weighted `lstsq` of `readings − v_act·g` on the columns `[1, h[0]]`.
   - The global coefficients (every law's intercepts, orders, `Ea/R` and
     log10 K) are optimised by `scipy.optimize.least_squares` on the stacked
     residuals.
   - Starts: Stage A's estimates, plus `restarts` Latin-hypercube draws within
     ±1 of each (±1 decade for log10 K). Keep the best.
   - Returns coefficients, standard errors (from the Jacobian), cost,
     per-curve rms in units of noise, per-curve `net_data` and `net_model`,
     `converged`, and `health` (below).
2. `global_fit_health(result)`: flags `"not converged"`, `"at bound: <K>"`,
   `"collinear: <a>/<b> r=±0.xx"` (|r| > 0.99).
3. `global_cross_validate(substrate, laws)`: for each run, fit on the other
   runs (starts = the full fit's coefficients, `restarts=1`), then predict the
   held-out run's curves. For each held-out curve, refit only its own `(c, v0)`
   by lstsq, then score `Σ ((readings − model)/(noise·sqrt(n)))²`.
   Skip-and-count as in Task 4.
4. `stage_shift(stage_a, stage_b)`: for every coefficient present in both, the
   shift in units of `sqrt(se_a² + se_b²)`, and the count beyond 2.
5. `amplitude_verdict(result, rows)`:
   - Per curve, `ratio = net_data/net_model` (non-positive ratios dropped and
     counted); per run, the median ratio.
   - Between runs: OLS of log median ratio on each of log[HOO-], log[buf],
     log[S], log[enz], 1/T, over runs, with t = slope/se.
   - Within runs: `scope.orders("ratio", frame=..., terms=("buf", "s0"),
     within=True, live_only=False)`. Never put `hoo` here.
   - Flag a between-run axis if |t| > 3 with at least 4 runs, and a
     within-run axis if |t| > 3.
   - Verdict: `"tracks conditions: <flags>"` if any flag, else `"no amplitude
     trend detected"`.

**Tests** (`data/test_rate_laws.py`):
- `test_a_strong_planted_global_model_is_recovered` — **A4**:
  - Curves: every 4OMe-BnOH table curve, with its real `times` and `noise`.
  - Truth: `v_act` law `S:power` 0.8, `BUF:power` 0.5, intercepts −13,
    −12.5, −12 (Boric, Phosphate, Pyrophosphate). Over the table's [S]
    (1.5-58 mM) and [buf] (3.1-200 mM) that puts v_act between about 1e-5
    and 1e-3 AU/s, the real range. `k_act_lag` law intercepts −6.5,
    `BUF:power` 0.3; `k_act_burst` law intercepts −7.5, no terms;
    `k_sink` None.
  - Each curve's family is its real one; `c = 0`, `v0 = 0.3·v_act` for lag,
    `3·v_act` for burst.
  - Readings = model + Gaussian noise at 0.1 × each curve's noise, seed 0.
  - Starts: the truth +0.5 in every coefficient.
  - Assert every coefficient within 0.05 of truth, and health has no flag.
- `test_c_and_v0_are_per_curve`: two planted curves with different c and v0
  under one law are both fitted to 1e-8.

**Pilot.** Time `global_fit` and one `global_cross_validate` fold with the
largest candidate laws, per substrate. The candidate set per substrate is,
for each of the three rate elements, the best tied model and the simplest tied
model (fewest coefficients), in all combinations; plus each combination with
`"k_sink"` None. That is at most 2 × 2 × 2 × 2 = 16 per substrate. Project
`candidates × (folds + 1) × time / workers`.

**Gate.** A4 or A5 → STOP. Otherwise commit and go straight on to Task 7.

## 13. Task 7 — the global fits — STOP 2

**Amendment 2 (section 11b) overrides this task's Run and Report where they
differ**:
- run `planted_global(substrate)` before the real fits, over the same
  candidates;
- beside every Stage B coefficient, report its planted twin's `within_2se`;
- quote a coefficient whose planted twin was not recovered as
  "not identified (planted recovery failed)".

**Run.** For each substrate, `global_fit` and `global_cross_validate` on every
candidate, saving `data/fits/rate_laws/global_<substrate>_<id>.json`. Then:
1. `tie_set` over the candidates;
2. `global_fit_health` on every tied candidate (A6);
3. `stage_shift` of the best tied candidate against Stage A (A7);
4. `amplitude_verdict` on the best tied candidate;
5. the same fits under `S-gas` and `S-pyro`, reporting only verdicts that
   change.

**Report.** Template E. Per substrate:
- curves used, and how many were not in any Stage A table;
- the candidates' CV scores and the tie set;
- the best tied candidate's coefficients with standard errors, unflagged ones
  only;
- the per-curve rms distribution in units of noise (median, 90th percentile);
- the `stage_shift` table;
- `amplitude_verdict`;
- the cut changes.

Then **STOP 2** (template S).

**Question at STOP 2:** "Stage B is above. Stage C — the ODE mechanism for what
the analytic form cannot represent (the 94 curves it misses, substrate
depletion, the peracid loop, the background) — needs its own plan written from
this result. Write it next?"

## 14. Report templates

Fill the brackets. Tables come from named functions. Add no sentences beyond
the template.

### Template E — the DATA_VERIFICATION.md entry

```
## <date> — <task name>: <verdict strings, verbatim, or "machinery only">

What was asked: <the task's Goal, one sentence>.
What was built: <function names>.
Bugs found: <before/after numbers, or "none">.

<the task's tables>

Verdicts (from <function names>): <verbatim>.
Health flags: <every flag, or "none">.
Not identified: <coefficients and reasons, or "none">.
Could overturn this: <the task's list below>.

Nothing is adopted.
```

"Could overturn this", per task:
- Task 1: a sheet layout the parser still misreads.
- Task 2: `remaining_curves` and the resolution flags decide which curves
  enter; a different form would admit different curves.
- Task 3: none.
- Task 4: pH, [enz] and temperature are between-run axes, confounded with day
  and batch; buffer identity is confounded with pH; the planted identifiability
  shows what this design can separate.
- Task 5: the links use only curves resolved in both elements.
- Task 5a: the two-test tie rule excludes a truly equal model somewhat more
  often than one test would; clustered errors are unreliable below 10 runs
  (flagged); the A8 anchors were re-scored from the old saves, whose fold
  scores covered only tied models; the tie rule cannot exclude a model that
  fails on 3 or fewer of 30 held-out runs, however badly —
  `worst_fold_log_ratio` shows where that happens.
- Task 5b: link power was measured with only buffer terms planted, while the
  real link fits carry each element's other terms, so the real power may be
  lower; a "ties" verdict at 4 of 5 power still has about a 1 in 5 chance of
  hiding a real difference.
- Tasks 6-7: the planted recovery is one seed per substrate, with the best
  candidate as its truth; a coefficient recovered there may not be recovered
  under a different truth.
- Tasks 6-7: the 94 curves the form cannot hold are not in the fit; Stage A's
  candidate set limits Stage B's.

### Template S — the stop message (in chat; nothing else)

```
STOP <n> — <task>

Done:
- <at most five bullets, each naming a function or a file>

<the task's tables>

Verdicts: <verbatim>
Health flags: <list, or "none">
Could overturn this: <the task's list>

Question: <the stop's question>
```
