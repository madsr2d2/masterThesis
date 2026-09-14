# Plan: three next steps towards the mechanism

Written 2026-09-13 for an agent picking this work up cold. Read the whole file
before touching anything. Do the steps **in order** — each one feeds the next.

This is a working plan, not a findings document. It is **not** read by any
number gate, so the numbers below are **baselines to reproduce**, not facts to
copy into other documents. Every number you put in any other document must come
from a function in a module that you can re-run (see Rule 6).

---

## Contents

0. [Rules of this repository — read first](#0-rules-of-this-repository--read-first)
1. [What is already established — do not re-derive](#1-what-is-already-established--do-not-re-derive)
2. [Step 1 — which peroxide species saturates the catalyst: H2O2 or HOO-?](#2-step-1--which-peroxide-species-saturates-the-catalyst)
3. [Step 2 — what activates the catalyst: the buffer, through the new clock](#3-step-2--what-activates-the-catalyst)
4. [Step 3 — extend the mechanism model with what the archive established, then fit](#4-step-3--extend-the-mechanism-model-then-fit)
5. [Reporting and committing, for every step](#5-reporting-and-committing-for-every-step)

---

## 0. Rules of this repository — read first

These are not style preferences. Every one was learned from a real mistake.

1. **Use `.venv/bin/python`, never `python`.** A bare `python3` has no pandas.
2. **Always run from the repository root** (`/home/madsr2d2/masterThesis`).
   Never `cd data`. The scripts open `data/experiment_data.csv` by a
   root-relative path and fail with `FileNotFoundError` from anywhere else. Run
   modules as `.venv/bin/python data/saturation.py`.
3. **Read `CLAUDE.md` in full before starting.** It is the rulebook. The
   `analyse-kinetics` skill (`.claude/skills/analyse-kinetics/SKILL.md`) is the
   short version.
4. **Gates must be green before every commit:**
   `.venv/bin/python run_gates.py` (about 130 s, 30 gates). The slow optimiser
   suite (`data/test_fit_kinetics.py`, about 9 minutes) runs only with
   `.venv/bin/python run_gates.py --all` — **run it only in Step 3, once, at the
   end**, not routinely.
5. **The duplicate guard** (`data/test_curve_metrics.py`) fails if any name —
   function, class, **constant, or private `_helper`** — is defined at top level
   in two modules. Before adding any name, run
   `grep -rn "def <name>\|^<NAME> =" --include=*.py .` and pick a unique name
   if it exists. This session hit it with `_planted` (renamed `_planted_binding`).
6. **Never compute a reported number inline in a throwaway script.** Put the
   computation in a function in a module (`data/saturation.py` for Steps 1–2,
   `data/kinetic_model.py` / `data/fit_kinetics.py` for Step 3) and print it from
   that module's `main()`. Throwaway numbers end up in documents and go stale.
7. **Never fit a progress curve yourself.** Fits are made once, in
   `scope.curve_fit`; `scope.frame(block)` is that object flattened into one row
   per curve. Read columns off the frame. Do not call `fit_progress` or
   `fit_activation_sink` on a curve in an analysis.
8. **Use the gas-corrected columns** (`*_corrected`). The O2 bubbles in the
   sample beam distort readings, worst at high pH and high peroxide.
9. **Test on planted data before real data.** Every new fit gets a test that
   plants a known answer on the REAL experimental design and recovers it — and a
   second planting with a DIFFERENT answer that the test must tell apart. A test
   that can only pass is not a test.
10. **In tight numerical loops use `np.einsum`, not `@`**, for matrix-vector
    products over many nodes. `@` goes to multithreaded BLAS, and with
    `run_gates.py` running 8 processes at once that doubled the suite's runtime.
11. **Report negative and inconclusive results as such.** "Not decided" is a
    valid result. Do not claim what an interval does not exclude.
12. **Do not edit numbers in `CLAUDE.md`, `BUBBLES.md`, `MECHANISM.md`,
    `FITTING.md` or the skill.** Those five are guarded by
    `test_root_documents.py`; a number added there needs a claim added to that
    gate. Record results in `DATA_VERIFICATION.md` instead (Section 5) and ask
    the user before editing a guarded document.
13. **Commit straight to `master`**, no branches. Commit message format is in
    Section 5.

---

## 1. What is already established — do not re-derive

Pointers, so you know where things live. Re-run the named function if you need
the number.

**The curve model.** Every curve is fitted with the *activation-sink form*
(`summary_kinetics.fit_activation_sink`), the product P = A − c obeying

    P' = v_act + (v0 − v_act)·e^(−t/τ) − k·P

`scope.frame` columns off it (all on the gas-corrected series):

| column | meaning | resolution flag |
|---|---|---|
| `v_act_corrected` | activated rate, no product made yet | `v_act_resolved_corrected` |
| `tau_act_corrected` | the catalyst's activation clock, s | `tau_act_resolved_corrected` |
| `k_act_corrected` | 1/τ, the activation rate constant, 1/s | use `tau_act_resolved_corrected` |
| `k_sink_corrected` | product-driven decline constant, 1/s | `k_sink_state_corrected == "resolved"` |
| `v0_act_corrected` | rate at mixing | `v0_act_resolved_corrected` |

Plus the model-free rates: `vmax_corrected` (steepest 20% block, a linear fit
over 5 consecutive blocks, max taken — robust, slightly upward-biased by
selection) and `v_peak_corrected` (analytic maximum of the fitted rate).

**Saturation machinery** — `data/saturation.py` (built this session):
- `binding_by_element(block, elements, gated)` — each element's binding constant
  K on the peroxide axis, per-run offsets, on the peroxide arm of the L
  (`induction.peroxide_ladder`). Forms: `"bound"` is K·h/(1+K·h), `"relaxation"`
  is (1+K·h).
- `shared_binding(block, elements, gated)` — one K shared by two elements vs a
  K each; F on 1 degree of freedom.
- `michaelis_by_element(block)` — per-RUN (Vmax, Km) for rate-like elements, on
  each run's substrate arm; carries `buffer_r`, the [S]/[buf] collinearity.
- `michaelis_temperature()` — Vmax vs Km barrier on the temperature series.
- `_scheme(shape, constant, x)` and `_levelled(y, offset, groups)` are the two
  helpers every new fit in Steps 1–2 should reuse.

**Chemistry findings you are building on** (details in `CLAUDE.md`,
`MECHANISM.md`, `DATA_VERIFICATION.md` 2026-09-11/12 entries):
- The induction is the catalyst activating on its own clock (`induction/`).
- The decline is a first-order sink on the product (`product_fate/`).
- The rate saturates in peroxide: K ≈ 0.04 /mM, catalyst ≈ 77% bound at 82.5 mM
  (`induction.peroxide_saturation`).
- The rate saturates in substrate: per-run Km resolved on 115 of 248 fits,
  median 2.3 mM on catalysed 4OMe phosphate.
- The temperature series' barrier is in Vmax, not Km.
- **The activation clock has no peroxide dependence** (k_act order +0.012,
  −0.24 to +0.26) while the activated rate carries +0.548; one shared binding
  constant for both is rejected at F = 4.87. This is the observation Steps 1–2
  follow up.
- The **buffer** axis meets the +1 pre-equilibrium rule through the landmark:
  `induction.joint_buffer_order` = +1.094 ± 0.150 — but on **only 2 runs,
  8 curves** (`induction.BUFFER_LEVER = (34, 32)`).

**The two schemes the +1 rule and the bounds come from** (read
`induction.joint_buffer_order`'s docstring and `MECHANISM.md` step 0):
- A species X in pre-equilibrium that draws the catalyst active:
  `1/τ = k_f[X] + k_r` (the "relaxation" form, `1/τ = k_r(1 + K[X])`), and the
  active fraction `K[X]/(1 + K[X])`. Log-slopes of rate and clock sum to +1.
- Bounds, rate-constant free: an **activating** species gives
  `d ln τ/d ln[X] ∈ (−1, 0)`; a species holding the catalyst **off** the path
  gives `(0, +1)`. A coefficient outside its bound falsifies the scheme.

---

## 2. Step 1 — which peroxide species saturates the catalyst

### 2.1 The question, in one paragraph

The rate saturates in peroxide. Is the species doing the binding **H2O2** or its
anion **HOO-**? Three existing results point at HOO- (the half order in [HOO-]
across three pH ladders; `early_trough`'s [enz]/[HOO-] driver; the gas's pH
onset) but none of them can name the *bound* species. This test can.

### 2.2 Why the design can answer it

The two-axis block's peroxide arm (`induction.peroxide_ladder(scope.frame(),
"vmax_corrected")`) holds **17 runs** spanning **pH 5.47 to 9.73**. Inside each
run pH is fixed and [H2O2] moves over 4 rungs (2–3 on exps 136, 150, 151). The
HOO- fraction `hoo/h2o2` is constant inside a run and varies by a factor of
**~17,600** between runs. (Reproduce these first — Checkpoint 1.)

With one free level per run:
- **H2O2 binds** → saturation sets in at the **same [H2O2]** in every run.
- **HOO- binds** → saturation sets in at a [H2O2] that shifts by
  `10^(pKa − pH)` between runs: high-pH runs are already saturated at their
  lowest rung, low-pH runs are straight lines.

The per-run levels absorb each run's height, but **not** where its curvature
sits, so the two hypotheses make different predictions the offsets cannot hide.

### 2.3 The statistic: a species exponent α, not a two-way contest

Do not just compare two SSEs. Define, per curve,

    x = h2o2 · (hoo / h2o2)^α

- α = 0 → x = [H2O2] (H2O2 binds)
- α = 1 → x = [HOO-] (HOO- binds)

For each α on a grid, profile K on the existing log grid and record the minimum
SSE. That gives a **profile over α**: an estimate α̂ with a 95% interval, using
the same convention the module already uses (`sse ≤ best · (1 + PROFILE_F /
degrees)`). This turns "which one wins" into a measurement with an error bar,
and it shows honestly when the data cannot decide.

### 2.4 Implementation — exact changes

All in `data/saturation.py`. Do not change the behaviour of any existing
function; `run_gates.py` must stay green with no edits to existing tests.

1. **Add a module constant** (check the name is unique first — Rule 5):

   ```python
   # The species exponent's grid: x = h2o2 * (hoo/h2o2)**alpha, alpha = 0 is
   # H2O2 binding and alpha = 1 is HOO- binding. The range reaches past both
   # so an estimate outside [0, 1] can be seen rather than clipped.
   SPECIES_ALPHAS = np.linspace(-0.5, 1.5, 81)
   ```

2. **Add `binding_species(block=scope.TWO_AXIS_BLOCK, element="vmax_corrected",
   gated=True, alphas=SPECIES_ALPHAS, grid=SATURATION_GRID,
   span=SATURATION_SPAN, cutoff=PROFILE_F, runs=None, exclude_bubbles=False)`.**
   Steps inside:
   - `table = scope.frame(block)`; look up `gate` and `shape` for `element` from
     `SATURATION_ELEMENTS`; `ladder = _element_rows(table, element, gate, gated)`.
   - If `runs` is given, keep `ladder[ladder.experiment.isin(runs)]`.
   - If `exclude_bubbles`, keep `ladder[ladder.bubble_load <= 1]`.
   - Drop rows with non-finite or non-positive `hoo`. If fewer than 10 rows
     remain, return `{"curves": n}` and stop.
   - `h = ladder.h2o2`, `fraction = ladder.hoo / ladder.h2o2`,
     `y = np.log(ladder[element])`, `groups = ladder.experiment.to_numpy()`.
   - `degrees = len(y) − (number of runs) − 2` (per-run levels, K, α). Guard
     `max(1, ...)`, **and** if it is below 1 return NaN — an exactly determined
     fit is not a measurement (DATA_VERIFICATION 2026-09-11, third entry).
   - For each α: `x = h * fraction**alpha`; for each K in
     `np.logspace(*span, grid)`: `sse = _levelled(y, _scheme(shape, K, x),
     groups)[0]`. Keep `profile[α] = min over K` and the K that achieved it.
   - `best = argmin(profile)`; interval = α values with
     `profile ≤ profile[best] · (1 + cutoff/degrees)`.
   - Also record SSE and best K at the grid nodes nearest α = 0 and α = 1, and
     `delta_aic = n · ln(sse_alpha0 / sse_alpha1)` (positive favours HOO-).
   - Return a dict: `element, curves, runs, alpha, alpha_low, alpha_high,
     alpha_at_edge` (True if the interval touches either end of the grid),
     `K_best, sse_best, sse_alpha0, K_alpha0, sse_alpha1, K_alpha1, delta_aic,
     degrees`.
   - **Cost:** 81 α × 600 K × one small `lstsq` each is about 49,000 solves per
     call. If that is too slow (> 60 s), reduce the K grid to 200 for the α
     profile only, then refine K at α̂ on the full grid. Measure before
     optimising.

3. **Add `species_table(block=scope.TWO_AXIS_BLOCK)`** that runs
   `binding_species` over the control matrix in 2.6 and returns one DataFrame,
   one row per (element, cut). Print it from `main()` under a heading
   `"which species saturates the catalyst"`.

### 2.5 Tests first — `data/test_saturation.py`

Write these **before** running on real data. Reuse `_with_frame` (already in
that file) to feed a planted frame.

1. **`_planted_species(alpha_true, K_true, noise, seed)`**: build a frame from
   the REAL peroxide arm — take
   `induction.peroxide_ladder(scope.frame(), "vmax_corrected")` and keep its
   `experiment, sample, h2o2, hoo, s0, pH, live` columns and its rows. Then
   overwrite `vmax_corrected` with
   `level_run · x/(1/K_true + x) · exp(normal(0, noise))`, where
   `x = h2o2 · (hoo/h2o2)**alpha_true` and `level_run` is a random factor per
   run (`10**normal(0, 0.4)`). Add any other columns the function reads, e.g.
   `bubble_load = 0.0`.
   - `K_true` must be chosen **per α** so the curvature sits inside the design,
     or neither planting is identifiable: for α = 0 use K ≈ 0.04 /mM (the real
     fitted value); for α = 1 choose K so `K · median(hoo) ≈ 1` —
     compute it from the frame, do not guess.
   - `noise`: use the real residual scatter. Run the fit once on real data,
     take `sqrt(sse_best / degrees)`, and plant at that value. Also run at half
     of it.
2. **`test_a_planted_species_comes_back`**: for α_true in (0, 1), the returned
   interval contains α_true.
3. **`test_the_two_species_are_told_apart`**: the α = 0 planting's interval
   excludes 1, and the α = 1 planting's interval excludes 0. **If this fails at
   the real noise level, the design cannot answer the question — stop, report
   that, and do not run the real-data interpretation.**
4. **`test_linear_data_leaves_alpha_unidentified`**: plant with K so small that
   every run is linear (K · max(x) < 0.05). The interval must span most of the
   grid (`alpha_at_edge` True). This protects against reading an unidentified α
   as a result.

Run `.venv/bin/python data/test_saturation.py` until all pass.

### 2.6 Run on real data — the control matrix

Checkpoint 1 before anything else:
`binding_by_element()` must still reproduce, for `vmax_corrected` gated:
order ≈ +0.598, K ≈ 0.0397, bound fraction ≈ 0.766. If not, stop — something
upstream changed.

Then `binding_species` for every row of this matrix:

| # | element | cut | why |
|---|---|---|---|
| A | `vmax_corrected` | all 17 runs | the headline |
| B | `vmax` (raw readings) | all runs | does the gas repair drive it? |
| C | `v_peak_corrected` | all runs | model-based rate |
| D | `v_act_corrected` (gated) | all runs | the activated rate |
| E | `vmax_corrected` | `runs=scope.strong_runs()` | weak runs flatten ladders |
| F | `vmax_corrected` | `exclude_bubbles=True` | the gas is heaviest at high pH |
| G | `vmax_corrected` | leave-one-run-out, 17 fits | is one run deciding it? report min/max α̂ |
| H | `vmax_corrected` | `runs=scope.PH_LADDER_TWO_AXIS_LOW` | exps 136–142, one composition set |
| I | `vmax_corrected` | `runs=scope.PH_LADDER_TWO_AXIS_HIGH` | exps 143–151, the other set |

Note: `SATURATION_ELEMENTS` does not list raw `vmax`. For row B add an optional
`shape` argument to `binding_species` (default: look it up) and pass
`shape="bound", gate=None` for `vmax`.

**The gas warning, specifically.** Gas is made from peroxide and rises steeply
with pH. Uncorrected, it adds signal at high [H2O2] in high-pH runs, which
*straightens* exactly the runs HOO- binding says should be most saturated. That
biases α̂ **towards 0**. So if rows A and B disagree, or A and F disagree, the
gas repair matters and the result must be reported with that caveat.

### 2.7 How to read the result

| α̂ interval (row A), consistent across E, F, G, H, I | verdict |
|---|---|
| contains 1, excludes 0 | **HOO- is the saturating species** |
| contains 0, excludes 1 | **H2O2 is the saturating species** — then the pH dependence of the rate must come from another step (e.g. general base on activation) |
| contains both | **not decided** — say so; carry both into Step 3 |
| entirely outside [0, 1] | **neither species alone** — likely a pH effect on K itself (e.g. the catalyst's own protonation). Report, do not force |
| `alpha_at_edge` True | **not identified** — the runs are too linear or too saturated |

If rows H and I give non-overlapping intervals, the two composition sets
disagree — report that as the main result, not the pooled one.

### 2.8 Deliverables for Step 1

- [ ] `binding_species` and `species_table` in `data/saturation.py`, printed by
      `main()`.
- [ ] The four tests in `data/test_saturation.py`, passing.
- [ ] `.venv/bin/python run_gates.py` green.
- [ ] A `DATA_VERIFICATION.md` entry (Section 5) with the full control matrix
      table and the verdict from 2.7.
- [ ] One commit.

---

## 3. Step 2 — what activates the catalyst

### 3.1 The question

The activation clock does not respond to peroxide. The buffer meets the +1 rule
through the landmark, on 2 runs. **Does the new, window-free clock
(`k_act_corrected = 1/τ`) confirm the buffer as the activating species, on all
five buffer titrations — and with what binding constant?**

### 3.2 The block

`scope.BUFFER_TITRATIONS = (32, 34, 35, 36, 37)`. Verified this session:

| exp | pH | [buf] rungs, mM | [S], mM | [enz], mM | run length, s |
|---|---|---|---|---|---|
| 32 | 7.00 | 50–200 (4) | 8.251 | 0.241 | 1767 |
| 34 | 7.00 | 3.125–25 (4) | 8.251 | 0.241 | 5280 |
| 35 | 7.50 | 50–200 (4) | 8.251 | 0.241 | 7920 |
| 36 | 7.53 | 50–200 (4) | 57.9 | 0.241 | 2100 |
| 37 | 7.53 | 50–200 (4) | 12.228 | 0.270 | 3210 |

All 4OMe-BnOH, phosphate, 40 °C, catalysed. 20 live curves; `tau_act` resolved
on **20 of 20**, `v_act` on 19 of 20. Inside every run [S] and [H2O2] are fixed
and only [buf] moves — so, unlike the 4OMe substrate ladders, **there is no
[S]/[buf] collinearity here**. pH, [S] and [enz] differ between runs and are
absorbed by the per-run levels. Reproduce this table first (Checkpoint 2).

### 3.3 What to measure, and the predictions

For the responses `k_act_corrected` (gate `tau_act_resolved_corrected`),
`v_act_corrected` (gate `v_act_resolved_corrected`), `vmax_corrected`,
`v_peak_corrected`, all with one free level per run:

1. **Free power** — the order in [buf], with interval.
2. **The bound check on the clock** (MECHANISM.md step 0): the order of **τ** in
   [buf] is **minus** the order of `k_act`. An activating species requires it
   in **(−1, 0)**; an inhibiting one in **(0, +1)**. Report which, and whether
   the interval sits inside the bound.
3. **Two forms for the clock**, compared by SSE (same parameter count):
   - relaxation: `1/τ = k_r · (1 + K[buf])` — `_scheme("relaxation", K, buf)`
   - saturating activation: `1/τ = k_max · K[buf]/(1 + K[buf])` —
     `_scheme("bound", K, buf)`
4. **The rates** in the bound form `K[buf]/(1 + K[buf])`.
5. **One shared K** for `v_act` (bound) and `k_act` (relaxation) against a K
   each — the saturating form of the +1 rule, as `shared_binding` does for
   peroxide.
6. **The +1 rule itself, through the new clock**, with
   `induction.joint_order(table, axis="buf", rate="v_act_corrected",
   timescale="tau_act_corrected", gate="tau_act_resolved_corrected",
   minimum=4, live_only=True)` — and again with `rate="vmax_corrected"`.
   Read `joint_order`'s docstring first. Pass **`gate=`, never `floor=`**, for a
   fitted clock (CLAUDE.md explains why).

### 3.4 Implementation — exact changes

In `data/saturation.py`:

1. **Generalise the axis.** `binding_by_element` and `shared_binding` hard-code
   `ladder.h2o2`. Add an `axis="h2o2"` argument to both, used as
   `x = ladder[axis]`, and a `rows=None` argument: a function
   `(table, column, gate, gated) -> DataFrame` defaulting to the existing
   `_element_rows`. With defaults, behaviour must be **byte-identical** —
   re-run Checkpoint 1 to prove it.
2. **Add `_buffer_rows(table, column, gate, gated=True)`**: live rows with
   finite, positive `column` and positive `buf`, gated if `gated`, restricted to
   runs where `buf` takes at least 3 distinct values (a buffer ladder). Check
   the name is unique.
3. **Add `activation_by_buffer(block=scope.BUFFER_TITRATIONS, gated=True,
   runs=None)`** returning one DataFrame with, per response: `curves, runs,
   order, order_low, order_high`; for `k_act_corrected` also
   `tau_order = −order` and `bound` ∈ {"activating (−1,0)", "inhibiting (0,+1)",
   "outside both", "straddles"} judged on the interval; `K_relaxation,
   sse_relaxation, K_bound, sse_bound`; for the rates `K, K_low, K_high`.
   Build it by calling the generalised `binding_by_element(block, elements,
   gated, axis="buf", rows=_buffer_rows)`, plus a second pass for
   `k_act_corrected` with `shape="bound"` to get the saturating-activation SSE.
4. **Add `activation_plus_one(block=scope.BUFFER_TITRATIONS)`** returning the
   `joint_order` results from 3.3 item 6 for both rates, on all five runs and on
   `induction.BUFFER_LEVER` alone.
5. Print both from `main()` under `"what activates the catalyst"`.

### 3.5 Tests first — `data/test_saturation.py`

1. **`_planted_activation(K_true, scheme, seed)`** on the REAL titration design
   (take `scope.frame(scope.BUFFER_TITRATIONS)` rows and overwrite the response
   columns): `k_act_corrected = level_run · k_r · (1 + K_true·buf)` for
   `scheme="relaxation"`, or `level_run · K·buf/(1 + K·buf)` for
   `scheme="bound"`; `v_act_corrected = level_run · K·buf/(1 + K·buf)`; set the
   gate columns True. Noise ~2% lognormal.
2. **`test_a_planted_buffer_K_comes_back`** — relaxation planting: K recovered
   within its interval, relaxation SSE < bound SSE.
3. **`test_an_activator_is_told_from_an_inhibitor`** — plant an **inhibitor**,
   `k_act = level · k/(1 + K·buf)`: the τ-order interval must land in
   (0, +1) and `bound` must say "inhibiting". Plant an activator: it must say
   "activating".
4. **`test_the_peroxide_axis_is_unchanged`** — `binding_by_element()` with all
   defaults returns exactly what it returned before your edit (freeze the
   current output as expected values in the test, computed once from the
   *unchanged* code before you edit anything).

### 3.6 Run on real data — controls

Checkpoint 2: reproduce the 3.2 table from `scope.frame(scope.BUFFER_TITRATIONS)`.

| # | cut | why |
|---|---|---|
| A | all five runs | the headline |
| B | `runs=induction.BUFFER_LEVER` (34, 32) | the same two runs the published +1 result used — compare like with like first |
| C | leave-one-run-out (5 fits) | one run should not decide it |
| D | without exp 34 | exp 34 is the ONLY run below 50 mM buffer; without it the curvature (K) may be unresolved — report whether K survives |
| E | `gated=False` | the price of resolution |

**Specific cautions for this block:**
- **Run length differs 4.5×** (1767 to 7920 s). Within-run contrast is safe,
  but the fitted clock's grid is capped at 2× the run length. Count curves
  with `tau_act_corrected > 1.8 · duration_s` and report them; if the
  shortest run (32) has capped clocks, repeat A without it.
- **pH 7.50–7.53 is near the gas onset** (none below 7.5 archive-wide). Check
  `bubble_load` on exps 35–37 and report the count above 0.
- **Only total buffer is measurable.** pH is one value per run, so acid and
  base forms are proportional inside a run. Do not claim which buffer species
  (`solution_chemistry.dominant_buffer_pair` docstring explains).
- **General base vs buffer perhydrate (MECHANISM.md step 6b) cannot be
  separated here:** they differ by a [buf]·[H2O2] term and no run moves both.
  Say so in the entry.

### 3.7 How to read the result

| finding | reading |
|---|---|
| τ-order in (−1, 0), relaxation fits ≤ bound form, +1 rule met through `tau_act`, shared K not rejected | **the buffer is the activating species, in pre-equilibrium**, with constant K — the strongest possible outcome |
| τ-order in (−1, 0) but +1 not met, or shared K rejected | buffer activates, but not by the simple pre-equilibrium scheme |
| τ-order in (0, +1) | buffer holds the catalyst **off** the path — contradicts the landmark result; check row B against the published +1.094 before believing either |
| τ-order interval straddles 0 | **not decided** on this block |
| row B disagrees with the published +1.094 ± 0.150 beyond errors | the landmark and the fitted clock measure different things on these runs — that is itself the result to report |

### 3.8 Deliverables for Step 2

- [ ] Generalised `binding_by_element` / `shared_binding` (defaults unchanged),
      `_buffer_rows`, `activation_by_buffer`, `activation_plus_one` in
      `data/saturation.py`.
- [ ] The four tests, passing; Checkpoint 1 still reproduces.
- [ ] Gates green.
- [ ] `DATA_VERIFICATION.md` entry with the control table and the 3.7 verdict.
- [ ] One commit.

---

## 4. Step 3 — extend the mechanism model, then fit

> **Superseded 2026-09-14.** Stages 3.0-3.1 below were done (commit 8954989).
> Stages 3.2-3.3 and the acceptance list are replaced by
> `PLAN_CURVES_TO_MECHANISM.md` (which replaced `PLAN_STEP3_REVISED.md`),
> after the review in `DATA_VERIFICATION.md`
> 2026-09-14 (second entry). Its verdict was built on F values read against
> the wrong models and on fits that did not reach their optimum.

### 4.1 The question, and why the existing fit cannot answer it

The archive has one mechanism fit (`data/fits/BnOH_25C_Phosphate.json`), made by
`data/fit_kinetics.py` on the reduced model in `data/kinetic_model.py`. Read
**`FITTING.md` sections F1–F7, "What to do next" and "What this does not license
saying"** in full before starting. The short version:

- **F1**: the reduced model is exactly first order in [S]; the data is not
  (+0.30 on BnOH).
- **F3**: without a lag mechanism, the model can only produce a lag if the
  extinction ratio r > 1, against the spectroscopic 0.08–0.33.
- **F5**: it misfits by **24.0×** (enzyme-free) and **20.6×** (catalysed) the
  noise on BnOH/25 °C/phosphate.
- **F6**: at its best fit the autocatalytic machinery is switched off.

And the reduction in `MECHANISM.md` ("Four-stage reduction") made assumptions the
archive has since shown to be wrong:

| reduction assumption | what the archive now says |
|---|---|
| stage 4: `K4[H2O2] ≪ 1`, catalyst fully free | rate ≈ 77% saturated in peroxide at 82.5 mM |
| no saturable [S] term | per-run Km resolved, 2–7 mM on 4OMe |
| catalyst fully active at t = 0 | activation on its own clock (step 0); `MECHANISM.md` itself says "the reduction below has not been redone under it" |
| aldehyde consumed only by the second-order Cannizzaro steps | a **first-order** sink on the product (`product_fate`) |
| no buffer route | step 6b proposed; buffer activates (Step 2) |

So Step 3 is **not** "run the fitter". It is: add the established terms to the
model one at a time, keeping the old model as the exact limit, prove each on
planted data, then fit and let each term earn its place.

### 4.2 Stop conditions — read before starting

Stop and report to the user, without continuing, if any of these happens:
- Steps 1 and 2 are not both committed with verdicts.
- The extended model at its default parameters does **not** reproduce the old
  model to < 1e-10 (4.4, test 1).
- Planted-parameter recovery fails for any parameter you intend to report.
- A fit hits a parameter bound you cannot justify.
- Any change makes an existing test in `data/test_kinetic_model.py` or
  `data/test_fit_kinetics.py` fail.

### 4.3 Stage 3.0 — reproduce the baseline (no code changes)

From the repository root:

```bash
.venv/bin/python data/fit_kinetics.py --list --scope all
.venv/bin/python data/fit_kinetics.py --substrate BnOH --temperature 25 --buffer Phosphate --scope all --save /tmp/v1_bnoh.json
.venv/bin/python data/fit_kinetics.py --substrate 4OMe-BnOH --temperature 40 --buffer Phosphate --scope all --save /tmp/v1_4ome40.json
```

Expected from `--list` (verified 2026-09-13): **BnOH 25 Phosphate 23 enzyme-free
/ 20 catalysed; 4OMe-BnOH 40 Phosphate 37 / 24**. These are the only two blocks
the sequential fit supports (F7). Expect the BnOH fit to reproduce F5's rms of
about 24× and 20.6× noise. Record both runs' `rms_sigma`, fitted constants and
correlations. **These are the numbers every later stage is compared against.**
The fits take minutes; use `run_in_background` if your tool supports it.

Do **not** use the two-axis block (BnOH 25 Pyrophosphate, 127 catalysed, 0
enzyme-free): it has no background data, and importing background constants
from another buffer is exactly what F7 forbids.

### 4.4 Stage 3.1 — extend `data/kinetic_model.py`, old model as the exact limit

**Do not create a second model module.** Every name in `kinetic_model.py`
(`RateConstants`, `Conditions`, `rates`, `rhs`, `simulate`, `observable`, …)
would be a duplicate and fail the guard. Extend in place, with every new
parameter defaulting to the value that switches it **off**.

1. **`Conditions`** — add `buf: float = float("nan")` and
   `species: float = float("nan")`, the concentration the catalyst binds
   (set from Step 1: `hoo` if HOO- binds, `h2o2` if H2O2). In
   `data/fit_dataset.py` (around the `Conditions(...)` call, ~line 352) pass
   `buf=float(row["[buf]"])` — `Curve` already carries `buf` — and `species`
   per the Step 1 verdict. Keep `hoo` exactly as it is.
2. **`RateConstants`** — add, each with an OFF default:

   | field | default (off) | meaning | source of prior |
   |---|---|---|---|
   | `k_sink` | 0.0 | first-order aldehyde loss, 1/s | `slowdown.sink_constants` |
   | `K4` | 0.0 | catalyst binding of the peroxide species, 1/mM | Step 1 K |
   | `km_s` | `inf` | substrate half-saturation, mM | `saturation.michaelis_by_element` |
   | `k_act_r` | `inf` | activation relaxation rate k_r, 1/s (inf = active at t = 0) | Step 2 |
   | `K_act` | 0.0 | activation constant on [buf], 1/mM | Step 2 K |

3. **`rates`** — change only the catalysed terms, and only by factors that
   equal 1 at the defaults:
   - active fraction `phi(t) = 1 − exp(−t · k_act_r · (1 + K_act·buf))` —
     exactly 1 when `k_act_r = inf`. `rates` has no time argument today; add
     one (`rhs` already receives `_time`), and update every caller.
     **If Step 2 selected the saturating-activation form instead of
     relaxation, use `1/τ = k_act_r · K_act·buf/(1 + K_act·buf)` and choose
     defaults that make phi ≡ 1.**
   - free catalyst fraction `free = 1 / (1 + K4 · species)` — 1 at `K4 = 0`.
   - substrate factor `sat_s = 1 / (1 + S/km_s)` — 1 at `km_s = inf`; use it
     so that the catalysed seed is `k5·E0·phi·free·[H2O2]·S·sat_s`.
   - `v6 = k6 · E0 · phi · free · [PBA]`.
   - `v_sink = k_sink · A`, and `dA/dt` gains `− v_sink` (the product goes to
     benzoate: add `+ v_sink` to BA through the aryl conservation — check
     `aryl_residual` still returns ~0).
   - **Do not add a term for the gas (S4).** [H2O2] consumption is < 1.5% on
     these blocks (`solution_chemistry.oxygen_budget`), so `[H2O2] ≈ const`
     still holds; the gas's competition for catalyst is absorbed into the
     effective constants. Say so in a comment.
   - Whether `sat_s` also applies to the **uncatalysed** seed `k0` is a model
     choice FITTING.md F1 raises (enzyme-free orders +0.28–0.49). Implement it
     as a separate switch `km_s_background` (default `inf`) and fit both
     variants in Stage 3.3.
4. **`PARAMETER_NAMES` / `LOG_PARAMETERS`** — append the new names.
   `K4`, `K_act` and `k_sink` can be exactly 0 (off), so they cannot be fitted
   in log10 as they are. Either fit them in log10 with a lower bound well below
   anything physical (e.g. `−8`) and treat "at the lower bound" as off, or add a
   separate linear-parameter list like `r`. Pick one, write down why in a
   comment.

### 4.5 Stage 3.1 tests — `data/test_kinetic_model.py`

Mirror the existing tests (`test_conservation`, `test_linear_in_e0`,
`test_enzyme_free_without_seed_is_frozen`, `test_acceleration_requires_r_above_one`, …).
**All existing tests must pass unchanged.**

1. **`test_defaults_are_the_old_model`**: for several `Conditions` and
   `RateConstants` (take the fitted constants from
   `data/fits/BnOH_25C_Phosphate.json`), the new `observable` at default new
   parameters equals a frozen copy of the old `observable` output to < 1e-10.
   Freeze those expected arrays by running the **unmodified** code once before
   editing, and save them in the test.
2. **`test_conservation_survives_the_extension`**: with every new term on,
   `aryl_residual` < 1e-8.
3. **`test_activation_makes_a_lag_with_r_below_one`**: with `r = 0.2` and a
   finite `k_act_r`, the observable accelerates (use
   `curve_metrics.acceleration`, pass `floor=` explicitly). This is the
   counterpart of `test_acceleration_requires_r_above_one`, which must still
   pass at the defaults.
4. **`test_substrate_saturation_order`**: the numerical d ln(initial rate)/d ln S
   with only `km_s` finite equals `km_s/(km_s + S)` to 1%.
5. **`test_peroxide_binding_saturates`**: the initial catalysed rate against
   `species` follows `x/(1 + K4·x)` to 1%.
6. **`test_sink_is_first_order`**: with only `k_sink` on and the seed off after
   a pulse of aldehyde, `A` decays as `exp(−k_sink·t)`.

### 4.6 Stage 3.2 — make the fitter carry the new parameters

In `data/fit_kinetics.py`:

1. Add each new parameter to `BOUNDS` and `INITIAL` with a comment giving the
   prior it came from (Steps 1–2, `slowdown.sink_constants`,
   `saturation.michaelis_by_element`).
2. Stage split: **stage 1 (E0 = 0)** may free `k_sink` and `km_s_background`;
   **stage 2 (E0 > 0)** frees `K4`, `km_s`, `k_act_r`, `K_act`. Leave
   `STAGE_ONE` / `STAGE_TWO` as they are and add **new** tuples (e.g.
   `STAGE_ONE_EXTENDED`), selected by a new `--extended` flag, so the old
   path is untouched.
3. **Synthetic recovery first** — add to `data/test_fit_kinetics.py`,
   following `test_parameter_recovery` and `_synthetic`: plant the extended
   model at the REAL conditions of the 4OMe/40 °C/phosphate block with the
   block's real noise; the fitter must recover every parameter you intend to
   report. Noiseless first (must be exact), then 3 noise seeds. Any parameter
   with |correlation| > 0.95 to another must be **profiled**, not reported as a
   value (FITTING.md F4 is the precedent: `k3` and `k5'` are lower bounds only).
4. The fitter reads raw readings (`curve.absorbance`). On these two blocks
   (pH 6.71–8.01 BnOH; pH 7.00–7.53 4OMe) the gas is minor, so keep raw
   readings for now and **state it**. Do not attempt the pyrophosphate block.

### 4.7 Stage 3.3 — fit, adding one term at a time

For **each** of the two blocks, fit this ladder of nested models, recording for
each: `rms_sigma` per stage, cost, every constant with its standard error and
`at_bound`, the correlation matrix, and the per-curve table (`_per_curve`:
`peak_data` vs `peak_model`, `net_data` vs `net_model`).

| model | adds | the question it answers |
|---|---|---|
| M0 | nothing (Stage 3.0 baseline) | — |
| M1 | `km_s` | F1: does substrate saturation fix the order? |
| M1b | `km_s` + `km_s_background` | does the background saturate too? |
| M2 | M1 + `K4` | does peroxide binding matter within one block? (weak lever: one [H2O2] per run on these blocks — expect K4 poorly determined; say so) |
| M3 | M2 + `k_sink` | does the product sink fix the late curve? |
| M4 | M3 + `k_act_r`, `K_act` | does activation supply the lag — and does `r` then drop below 1? |

**Each term must earn its parameters.** Use the weighted cost: compute
`F = ((cost_small − cost_big)/Δp) / (cost_big/(N − p_big))` and require F > 12
(the package's `TWO_PHASE_F` bar, which allows for serially correlated
residuals). Report terms that do not earn as "not supported", and do not carry
them into the next model.

**Also run each model with Step 2's constants FIXED** (K_act from Step 2,
`k_act_r` from Step 2's k_r) and compare with the model where they are free.
Agreement is a strong cross-check; disagreement is a result.

### 4.8 Acceptance — what "the extended model works" means

Compare against M0. State each explicitly as met or not met:

1. `rms_sigma` falls substantially below M0's (~24× / 20.6× on BnOH).
2. The **fitted r** moves into or towards the spectroscopic 0.08–0.33 once M4
   supplies a lag (FITTING.md F3).
3. The model's **lag fraction** (curves with `peak_model > 0.15`) matches the
   data's (`peak_data > 0.15`) within a few curves.
4. The model's **substrate order**, computed from simulated initial rates
   across each run's rungs, reproduces the data's: +0.30 on BnOH/25 °C
   (FITTING.md F1), +0.60 ± 0.25 on 4OMe/40 °C.
5. Fitted `k_sink` is consistent with `slowdown.sink_constants` for that
   temperature and substrate.
6. Fitted `K_act` is consistent with Step 2's buffer K (4OMe/40 °C block
   contains the titrations, so this is the same data — consistency is required).
7. The autocatalytic machinery is **not** switched off (F6): report `[PBA]` at
   the fitted constants.

### 4.9 Deliverables for Step 3

- [ ] Extended `kinetic_model.py` with OFF defaults; the 6 new tests passing and
      every existing test unchanged.
- [ ] `--extended` path in `fit_kinetics.py`; synthetic recovery test.
- [ ] Fits M0–M4 on both blocks, saved with `--save` into `data/fits/` with new
      filenames (e.g. `4OMe-BnOH_40C_Phosphate_M4.json`); do not overwrite
      `BnOH_25C_Phosphate.json`.
- [ ] `.venv/bin/python run_gates.py --all` green (the one time the slow suite
      runs).
- [ ] A `DATA_VERIFICATION.md` entry with the M0–M4 table per block and the
      acceptance checklist.
- [ ] **Ask the user** before editing `FITTING.md` or `MECHANISM.md` (guarded
      documents, Rule 12). Propose the text; do not commit it unasked.

---

## 5. Reporting and committing, for every step

**`DATA_VERIFICATION.md` entry** — new entries go at the **top**, under the `---`
after the header. Heading format:

    ## 2026-MM-DD — <one line saying what was found>

(Add "(second entry)" etc. if one already exists for that date.) Contents, in
this order:
1. What was asked and why (2–3 sentences).
2. What was built (function names).
3. **Any bug found on the way**, with the before/after numbers.
4. The results table(s), each number from a function you name.
5. The verdict, using the reading tables in this plan — including
   "not decided" if that is what it is.
6. Cautions (gas, weak runs, collinearity, anything that could have produced
   the result).
7. The tests added.
8. The words **"Nothing is adopted"** unless the user has approved a change to a
   published result.

**Before committing:**

```bash
.venv/bin/python run_gates.py        # must say "0 failed"
git status --short                   # only your intended files
```

**Commit** straight to master. Message: a short `area: what changed` subject
line, a paragraph on what was found, then these two lines exactly:

    Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
    Claude-Session: https://claude.ai/code/session_01Vht919cHDPUdJWghusZFdS

**Tell the user** after each step, in plain words: the verdict, the one or two
numbers that decide it, and what could still overturn it. Then wait for their
go-ahead before starting the next step.
