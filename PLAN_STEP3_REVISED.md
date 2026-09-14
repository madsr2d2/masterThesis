# Step 3, revised — a ladder that can answer the question

Handover plan, written 2026-09-14 after the review in `DATA_VERIFICATION.md`
2026-09-14 (second entry). It REPLACES stages 3.2-3.3 and the acceptance list
of `PLAN_MECHANISM_NEXT_STEPS.md` section 4. Stages 3.0-3.1 of that plan are
done (commit 8954989) and are not repeated.

## Contents

0. Read first
1. Where Step 3 stands
2. What Steps 1 and 2 can and cannot feed in
3. Stage R0 — fix the machinery before fitting anything
4. Stage R1 — look at the 4OMe catalysed residuals first
5. Stage R2 — the background, with its confound
6. Stage R3 — the catalysed ladder on one background
7. Stage R4 — the joint fit
8. Acceptance
9. Stop conditions
10. Deliverables, reporting, committing

---

## 0. Read first

- **Every rule in `PLAN_MECHANISM_NEXT_STEPS.md` section 0 still applies**:
  `.venv/bin/python` from the repository root, never `python`; the duplicate
  guard covers private names and constants; `run_gates.py` green before every
  commit; commit straight to master with the two attribution lines; new
  `DATA_VERIFICATION.md` entries at the top; ASK the user before editing
  `FITTING.md` or `MECHANISM.md`.
- Read `DATA_VERIFICATION.md` 2026-09-14 (second entry) IN FULL. It lists what
  went wrong last time; most of this plan exists to stop it recurring.
- **Fits take minutes to tens of minutes.** Run them in the background and
  save every one with `--save`. Never overwrite `data/fits/BnOH_25C_Phosphate.json`
  or the `_M*.json` saves; this round's saves are named `_R<rung>.json`.
- `data/test_fit_kinetics.py` is the slow optimiser suite (~9 min). Run it
  when you change the fitter and once at the end (`run_gates.py --all`), not
  after every edit.

## 1. Where Step 3 stands

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

**Does not stand:** the entry's verdict, the recovery test as a recovery test,
and every F from a pair where the larger model fits worse.

**The four facts the new ladder has to respect:**

| fact | number | consequence |
|---|---|---|
| a child fitted from scratch can land above its parent | 4OMe M1, M2 stage 2 | warm-start every child from its parent (R0.1) |
| stage-2 costs on different backgrounds are not comparable | M1b, M3 | one shared background for every stage-2 rung (R0.2) |
| 4OMe enzyme-free ladders move [buf] against [S] | median r -0.974, 9 of 9 runs | fit a buffer-order background beside the saturating one (R2) |
| 4OMe catalysed misfit is 378-650x noise under every model | `ladder_f_test` rms | inspect residuals before adding terms (R1) |

## 2. What Steps 1 and 2 can and cannot feed in

- **Peroxide binding (`K4`): drop it from this round.** Both phosphate blocks
  hold one [H2O2] per run, so `K4` is not identified whatever it binds (BnOH
  K4/k5 correlation 1.000). Step 1's exponent measures how the effective
  binding constant moves with pH on the two-axis block; it does not name a
  species these blocks could use. Leave `species` as total [H2O2] and `K4` at 0.
- **The activation (`k_act_r`, `K_act`): no prior and no fixed-constant run.**
  The review of Step 2 found the buffer design cannot tell an activating buffer
  from an inhibiting one at the real clock scatter, and that runs 35 and 37 mix
  fit forms across their buffer rungs. So Step 2's K is not a constraint.
  The old plan's "fit with Step 2's constants FIXED" is dropped for that
  reason. The activation is tested here on its own evidence, and cross-checked
  against the progress fits' own `tau_act_corrected` (Acceptance 5).
- **The product sink:** `slowdown.sink_constants` is the comparison value for
  a CATALYSED 4OMe sink (Acceptance 7).

## 3. Stage R0 — fix the machinery before fitting anything

All in `data/fit_kinetics.py` (and `data/kinetic_model.py` for R0.4), each
with tests. **Do not change what any existing default path computes**: the
existing tests must pass unchanged.

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
  - `test_a_warm_started_child_never_fits_worse`: plant M0-shaped data, fit
    M0, then fit M1 with `parent=` the M0 result; assert
    `M1.cost <= M0.cost * (1 + NESTING_TOLERANCE)`.
  - `test_a_child_above_its_parent_raises`: monkeypatch the fitted cost above
    the parent's and assert the `RuntimeError`.

### R0.2 One background for every stage-2 rung

- Add `--background PATH` to `main`: load stage 1's constants from a save and
  skip stage 1 (`sequential_fit(..., background=RateConstants)`), so every
  stage-2 rung sits on the identical frozen base.
- Save the background's path inside each stage-2 save (`"background": PATH`),
  and extend `ladder_f_test` to check it: two stage-2 rows are comparable only
  if their `background` fields match (in addition to the constants check it
  already makes).
- Test: `test_stage_two_rungs_share_the_named_background` on planted saves.

### R0.3 A `parents` argument to the ladder readers

`ladder_f_test` and `ladder_checks` take `models=`; add `parents=MODEL_PARENTS`
so the R-ladder below can pass its own map. Test on a planted R-ladder.

### R0.4 The terms this round needs, each OFF at its default

In `kinetic_model.py`, extending `RateConstants` exactly as 8954989 did:

| field | default (off) | meaning |
|---|---|---|
| `k_act_f` | 0.0 | forward activation rate on [buf], 1/(mM s) |
| `n_buf_background` | 0.0 | order of the uncatalysed seed in [buf]: `k0 [H2O2][S] ([buf]/BUF_REFERENCE)^n` |

- The activation exponent becomes
  `t * (k_act_r * (1 + K_act*[buf]) + k_act_f*[buf])`, and `phi` stays exactly
  1 while `k_act_r` is infinite. **Why:** `k_act_r` and `K_act` correlate at
  1.000 on 4OMe because only their product is identified when the fit wants
  `k_r -> 0`. Fitting `(k_act_r, k_act_f)` with `K_act` held at 0 is the same
  relaxation, `1/tau = k_r + k_f[buf]`, without the degenerate pair.
- `BUF_REFERENCE`: one module constant, e.g. 75.0 mM, the middle of the
  archive's buffer range, so `k0` keeps its units and scale. Check the name is
  unique.
- The buffer factor applies to the UNCATALYSED seed only, and is 1 when
  `n_buf_background` is 0 or `[buf]` is NaN.
- Add both to `PARAMETER_NAMES`, `BOUNDS`, `EXTENDED_INITIAL`, `UNITS`.
  `n_buf_background` is fitted LINEARLY (not in `LOG_PARAMETERS`), bounds
  (-2, 2), like `r`. `k_act_f` in log10, bounds (-10, 0).
- Tests in `data/test_kinetic_model.py`:
  - extend `FROZEN`'s check: defaults still reproduce the old arrays to 1e-10;
  - `test_forward_activation_is_the_same_relaxation`: `(k_act_r=k, K_act=K)`
    and `(k_act_r=k, K_act=0, k_act_f=k*K)` give the same observable to 1e-10;
  - `test_the_background_buffer_order`: the numerical d ln(initial seed)/
    d ln[buf] with only `n_buf_background` on equals it to 1%.

### R0.5 Rewrite the recovery test so it tests recovery

Replace `test_extended_parameter_recovery` in `data/test_fit_kinetics.py`
(keep its name so the gate list does not change):

- **Real conditions**: take the 4OMe/40 C/phosphate curves from
  `build_curves()` (`c.group == ("4OMe-BnOH", 40.0, "Phosphate")`), keep their
  own `times` and `conditions`, and replace `absorbance` with the planted
  observable plus Gaussian noise at each curve's own `noise`.
- **Starts away from the truth**: `EXTENDED_INITIAL` must be at least 0.7
  decades from every planted log10 value; assert that in the test itself.
  `restarts=4`. Stage 1 is FITTED, not given.
- **Truth in the regime the ladder will report**: pick the planted values so
  the model lags on some catalysed curves (check `peak_model`).
- Noiseless first (every planted parameter to 0.1 decade), then **three** noise
  seeds. A parameter is "recovered" only if all three seeds put it within 0.3
  decade. Assert `k_act_r` and `k_act_f` separately. Anything not recovered is
  never reported as a value later — say "not identified" instead.
- This runs under `--all` only. Time it; if it exceeds 15 minutes, cut the
  curve count to one run per experiment and say so in the docstring.

## 4. Stage R1 — look at the 4OMe catalysed residuals first

The 4OMe catalysed block sits at 378-650x the noise under every model. A term
cannot be judged against a misfit that large, so find out what it is.

- Add `fit_kinetics.residual_shape(block, model, stage=2, directory="data/fits")`.
  For each curve of that save: re-simulate with the saved constants
  (`observable`, the curve's own times, `baseline_like_data` as `_per_curve`
  does), and report the mean residual (data - model) in the first, middle and
  last third of the run, in AU and in units of the curve's noise, plus
  `net_data/net_model`. Return one row per curve, sorted by rms.
- Print it for `4OMe-BnOH_40C_Phosphate_M4b` and `_M1b`. Read the top ten
  curves and answer, in the entry:
  1. Is the misfit an AMPLITUDE error (`net_data/net_model` far from 1, the
     same sign in every third)? Then it is an extinction/`r`/epsilon or
     concentration question, not a missing term.
  2. Is it EARLY (first third) — a lag the model lacks or overshoots?
  3. Is it LATE (last third) — a fall the model lacks (the product sink)?
  4. Is it concentrated in a few experiments? Name them and check
     `KNOWN_EXCLUSIONS` / `RULINGS` in `data/build_manifest.py` for them.
- **STOP after R1 and report to the user** if the answer is (1) or (4) for
  most of the cost. Do not start R2 on a block whose misfit is an amplitude or
  a handful of runs.
- Test: `test_residual_shape_is_zero_on_a_planted_fit` (planted curves fitted
  exactly give thirds within 1e-8 AU).

## 5. Stage R2 — the background, with its confound

Stage 1 only, both blocks, each rung warm-started from its parent.

| rung | frees, added to STAGE_ONE | question |
|---|---|---|
| B0 | — | baseline (reuse M0's stage 1 as `parent` only; refit) |
| B1 | `km_s_background` | a background saturating in [S] |
| B2 | `n_buf_background` | a background with an order in [buf] |
| B3 | `km_s_background`, `n_buf_background` | both |

- B1 and B2 are NOT nested in each other; they have the same parameter count,
  so compare their costs directly. B3 against each is an F.
- Read each with `ladder_f_test(parents={"B1": "B0", "B2": "B0", "B3": "B1"})`
  (and again with B3 against B2) and `ladder_checks`.
- **On 4OMe the reading is in the table below.** On BnOH, where 3 of 5
  enzyme-free ladders hold [buf] fixed, the fixed-buffer runs are the
  discriminating ones: also report `ladder_checks`' substrate order on those
  three runs alone (exps 67, 69, 70), data against B1 and B2.

| outcome | reading |
|---|---|
| B1 and B2 within a factor 1.05 in cost, B3 does not earn | **not decided** — [S] saturation and a buffer order cannot be told apart on this block; carry B2 AND B1 into R3 and report both |
| B2 clearly lower, and its substrate order matches the data within 2 sigma | the background has a buffer order, not a saturation |
| B1 clearly lower, and its substrate order matches the data within 2 sigma | the background saturates in [S] |
| the winner's substrate order misses the data's by more than 2 sigma | it wins by bending something else; not accepted |

- Add `k_sink` to stage 1 as rung **B4** (`B4` parent = the R2 winner) only as
  a CONTROL: a background sink that earns would contradict `product_fate`'s
  enzyme-free clock, and must be reported, not carried.
- The chosen background's save is the `--background` for all of R3.

## 6. Stage R3 — the catalysed ladder on one background

Stage 2 only, `--background` = the R2 choice (run it twice if R2 is not
decided), each rung warm-started from its parent.

| rung | parent | frees, added to STAGE_TWO | question |
|---|---|---|---|
| C0 | — | — | catalysed baseline on this background |
| C1 | C0 | `km_s` | catalysed substrate saturation (**BnOH only**: on 4OMe only exp 16 steps [S]; say "not identified" and skip) |
| C2 | best of C0/C1 | `k_sink` | the catalysed product sink |
| C3 | C2 | `k_act_r` | activation on one clock (every block) |
| C4 | C3 | `k_act_f` | the clock moving with [buf] (**4OMe only**: no BnOH catalysed run moves [buf]) |
| C5 | best so far | `r` | does `r` move once a lag is supplied? |

- C5 frees `r` in stage 2. Because `r` also shapes the background, report the
  enzyme-free rms at C5's `r` as well (re-simulate stage 1 with it). If the
  background worsens by more than 10% in rms, `r` is being bought by the
  catalysed curves at the enzyme-free curves' expense — say so, and defer the
  question to R4.
- "Best of" means the lower cost among COMPARABLE rows whose verdict is
  "earns"; a term that does not earn is not carried.
- Read every rung with `ladder_f_test(parents=...)` and `ladder_checks`.

## 7. Stage R4 — the joint fit

Take the final R3 model and fit ALL curves of the block at once, both stages'
parameters free, starting from the sequential optimum (`extra_starts`). This
is the only place `r` is fitted to both kinds of curve together. Report rms
per kind of curve, and whether any constant moved by more than its standard
error from the sequential value. A large move means the sequential split was
hiding a trade-off; report it.

## 8. Acceptance

State each as **met** or **not met**, with its number, in the entry. Nothing
counts from a row whose verdict is not "earns" or from an unconverged fit.

1. **No larger rung fits worse than its parent**, and every comparable row
   converged (R0.1 guarantees the first; if it ever raises, that is a stop).
2. **rms in AU** against C0 on the same background, per kind of curve.
3. **Substrate order**: the final model's order (`ladder_checks`) within 2
   sigma of the data's, per kind of curve. On 4OMe enzyme-free, state which
   background reading R2 reached.
4. **Lag counts**: model within 3 curves of the data per kind of curve.
5. **Activation clock**: for catalysed curves whose `act_sink_kind_corrected`
   is a lag kind ("lag" or "lag then sink") in
   `scope.frame(tuple(save["stage_2"]["experiments"]))`, the median ratio of the
   model's `1/(k_act_r + k_act_f*[buf])` to that curve's `tau_act_corrected`
   is within a factor of 2. Exclude, and count, runs whose rungs mix fit forms.
6. **r**: at C5 or R4, inside 0.08-0.33, or state where it went.
7. **Sink**: C2's `k_sink` against `slowdown.sink_constants` for 4OMe at 40 C,
   within its error.
8. **Machinery (F6)**: report the peak [PBA] at the final constants on three
   representative catalysed curves.

## 9. Stop conditions

Stop and report to the user, without continuing, if:
- any existing test fails after an R0 change;
- R0.5's recovery fails for a parameter a later stage would report (then do
  not report it; if it is `k_act_r` or `k_act_f`, stop);
- R1 finds the 4OMe misfit is an amplitude error or a handful of runs;
- R0.1's guard raises on a real fit;
- a fitted constant sits at a bound that is not its OFF end and you cannot
  explain why.

## 10. Deliverables, reporting, committing

- [ ] R0: `extra_starts`, `parent=`, the nesting guard, `--background`,
      `parents=`, `k_act_f`, `n_buf_background`, the rewritten recovery test,
      and their tests. `run_gates.py` green; `run_gates.py --all` once.
      **One commit.**
- [ ] R1: `residual_shape`, its test, and a `DATA_VERIFICATION.md` entry with
      the top-ten table and the answer. **One commit. Then tell the user and
      wait.**
- [ ] R2-R4: the saves `data/fits/<block>_R<rung>.json`, and one
      `DATA_VERIFICATION.md` entry with, per block: the R2 table and reading,
      the R3 `ladder_f_test` and `ladder_checks` tables, the R4 comparison,
      and the acceptance list. **One commit.**
- Entries follow `PLAN_MECHANISM_NEXT_STEPS.md` section 5 and end "Nothing is
  adopted".
- Propose any `FITTING.md` / `MECHANISM.md` text to the user; do not commit it
  unasked.
- After each commit, tell the user in plain words: the verdict, the one or two
  numbers that decide it, and what could still overturn it.
