# Which mechanisms survive the progress curves

Handover plan, written 2026-09-15. It continues `PLAN_CURVES_TO_MECHANISM.md`,
whose Stages A and B closed at STOP 3 (commit 031bccf). This is its Stage C.
Read this file from the top. You do not need the earlier plan, except where a
section names a function from it.

> **AMENDMENT 2 (2026-09-16). Stage C is done through STOP 2 (commits 0577792,
> 5fa9641). Read section 13 and run it; it is the last task in this plan.**
> Amendment 2 adds ONE candidate, C5, because `residual_trends` found a
> within-run peroxide dependence in the late shape that every candidate misses.
> The write-up is NOT yours: do not start one.

> **AMENDMENT 1 (2026-09-16). Tasks 1 and 2 are done (commits 8211a58,
> 6f87773). Task 3's fits are all on disk and its table is fixed; read section
> 12 BEFORE finishing Task 3.** Amendment 1 records a bug a review fixed in
> your working tree, adds the discrimination anchors, and stops a substrate
> whose planted truth is not recovered from being read in Task 4.

## The goal, in one paragraph

Take five hand-written candidate mechanisms for the catalysed curves, score
each on how well it PREDICTS runs it did not see, and report which candidates
survive and which are excluded. Before looking at the real result, measure on
planted data which candidates this archive can tell apart at all. **The
deliverable is a verdict table, not a best-fit set of equations.** No
candidate is "the mechanism".

## Why this plan is written the way it is

The previous plan failed in ways this one is built to avoid.

- **Stage B let every curve keep its own free starting rate**, in the fit and
  in the cross-validation. The laws were never tested on how big a curve is,
  and the optimiser pushed the clocks outside the run window. So here:
  - no curve and no run gets a free amplitude;
  - run-to-run variation is a random effect whose size is estimated;
  - a held-out run is predicted, never refitted.
- **The activation-sink parameters scatter twice as much between repeat runs
  as window-free summaries of the same curves** (ln SD 0.64-0.67 against
  0.25-0.37). So candidates are scored on three summaries read identically
  off the readings and off the candidate's curve (section 3), never on fitted
  curve parameters.
- **A planted test passed while the real fit was degenerate**, because the
  planted truth sat in a well-behaved region. So:
  - the planted truths here are each candidate's own fit to the REAL data;
  - the planted noise is that fit's own estimated scatter;
  - every fit, planted or real, gets the same degeneracy flags.
- **Earlier rounds went wrong where a plan left a judgement to the executor.**
  So:
  - every verdict is a string returned by a function, with fixed thresholds;
  - anchors must reproduce before new code is trusted;
  - assumptions have stop triggers;
  - the words for results are fixed (section 2).

---

## Contents

0. Rules for the executor
1. What is already true
2. Glossary and verdict strings
3. The data and the three summaries
4. The five candidates
5. Assumptions and their stop triggers
6. Tasks and stops
7. Task 1 — the summaries and their scatter
8. Task 2 — candidates, likelihood, fitter, cross-validation, pilot
9. Task 3 — real full fits and the planted discrimination — STOP 1
10. Task 4 — real cross-validation and the verdicts — STOP 2
11. Report templates
12. Amendment 1 — the fixed table, its anchors, and what BnOH may not be read for
13. Amendment 2 — C5, the oxidant-depletion candidate — STOP 3

---

## 0. Rules for the executor

**Repository rules** (from `CLAUDE.md`):
- **Interpreter.** Use `.venv/bin/python`, from the repository root. Never
  `python`, never from inside a folder.
- **New names.** Before adding any top-level name (function, class or
  constant, private or not), `grep -rn "NAME" --include=*.py .` must find
  nothing. The duplicate guard in `data/test_curve_metrics.py` fails
  otherwise. Every name this plan introduces was checked unused on
  2026-09-15; check again anyway.
- **Gates.** `.venv/bin/python run_gates.py` must print `0 failed` before
  every commit.
- **Commits.** Commit straight to master. Subject `area: what changed`, one
  paragraph, then exactly these two lines:
  ```
  Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01Vht919cHDPUdJWghusZFdS
  ```
- **Log.** Every task that produces a result gets a `DATA_VERIFICATION.md`
  entry at the top of the file, using template E (section 11). Never edit an
  earlier entry. Correct one by writing a new entry.
- **Protected documents.** Never edit `FITTING.md`, `MECHANISM.md`,
  `CLAUDE.md`, `BUBBLES.md`, `.claude/skills/` or `PLAN_CURVES_TO_MECHANISM.md`.
  Propose text in a stop message instead.
- **Saves.** Never overwrite a file in `data/fits/`. This plan saves under
  `data/fits/mechanism_discrimination/` only.
- **Worker processes.** Parallel workers set `OMP_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1` in the environment BEFORE
  numpy is imported.
- **Long jobs** run in the background, save each fit as it finishes, and skip
  any fit whose save already exists.
- **Numbers.** A number that goes into a document comes from a function in
  the repository, never from a throwaway script.

**Rules that exist because they were broken before:**
1. **Tests.** Never loosen an assertion, threshold or tolerance written in
   this plan. Never start a fit at a planted truth. Build planted data exactly
   as written. If a test fails, stop and report the failing values.
2. **Fits.** Never read a parameter from a fit that `candidate_degeneracy`
   flags. Write "not read (<flag>)".
3. **Words.** Use section 2. Forbidden unless a function in this plan returns
   them as a verdict string:
   - that a candidate "is the mechanism", "is established", "is required",
     "is correct" or "is wrong";
   - that a result "confirms", "proves", "rescues" or "rules out" anything;
   - that a difference is "small", "negligible" or "does not matter";
   - that the data "prefer" a species or a candidate.
4. **Stops.** At a STOP: commit, post the stop message (template S), and do
   nothing more until the user replies.
5. **No improvising.** Do not add a candidate, a summary, a parameter, a bound
   or a start. Do not change a window or a threshold. If the plan seems wrong,
   stop and say why.

## 1. What is already true

**1.1 The readings.** `scope.fits(scope.archive())` returns
`{(experiment, sample): CurveFit}`. A `CurveFit` has:
- `times`: seconds, and `times[0]` is 0.0 on every curve used here;
- `corrected`: the gas-corrected readings, in AU;
- `noise`: the per-reading noise, in AU.

`scope.frame(scope.archive())` has one row per curve, with `experiment`,
`sample`, `substrate`, `buffer`, `temperature` (°C), `pH`, `s0`, `h2o2`,
`hoo`, `buf`, `e0` (all mM), and `live`. The frame takes minutes to build the
first time in a process, and it is memoised after that.

**1.2 A catalysed curve is a catalytic INCREMENT.** Its reference cuvette holds
the same mixture minus the enzyme. The candidates below therefore describe
only the catalyst's contribution. The uncatalysed background, including BnOH's
catalyst-free autocatalytic loop (steps 1-3 of `MECHANISM.md`), is assumed to
cancel in the subtraction. That is an approximation, and it is on the
"Could overturn this" list.

**1.3 Design facts** (checked 2026-09-15):
- **Temperature** moves only in 4OMe-BnOH (15-40 °C). Every BnOH curve is at
  25 °C, so BnOH candidates carry no barrier parameters.
- **pH** is one value per run everywhere. **[enz]** moves only between runs.
- **Buffer identity is confounded with pH.**
- **No run moves [buf] and [H2O2] together.** Candidates C1 and C2 may
  therefore be impossible to separate. Task 3 measures that.
- **[S] and [buf] move together inside many 4OMe-BnOH runs.**

**1.4 What earlier work found** (the source of the candidates; none is assumed
true here):
- **Activation is on a clock.** The catalysed curves begin with the catalyst
  activating on its own clock, which needs the catalyst.
- **The clock slows with pH on four ladders**, the sign a base-held trap
  requires.
- **The buffer axis meets the pre-equilibrium "+1" rule**, and the peroxide
  axis does not.
- **The rate rises with [HOO-]** and saturates in peroxide and in [S].
- **4OMe's catalysed rate falls with the product made.**

## 2. Glossary and verdict strings

| word | meaning here |
|---|---|
| summary | one of L, E, D (section 3), read off a curve |
| candidate | one of C0-C4 (section 4) |
| full fit | a candidate fitted to every curve of a substrate |
| fold | one run held out; the candidate is fitted to the other runs and scored on that run |
| fold score | the held-out run's negative log-likelihood under the fold's fit (lower is better) |
| best | the candidate with the lowest total fold score |
| tied | kept by `candidate_tie`: its fold scores are not worse than the best's beyond 2 standard errors |
| excluded | not tied |
| distinguishable | on planted data, each of the two candidates excludes the other when it is the truth |

**Verdict strings.** Only these may describe a result, and only as a function
returns them:
- `candidate_tie`: `"best"`, `"tied"`, `"excluded"`;
- `discrimination_table`, per pair: `"distinguishable"`,
  `"one-way (<A> as truth excludes <B>)"`, `"not distinguishable"`; per
  truth: `"truth recovered"`, `"truth not recovered"`;
- `candidate_verdicts`, per candidate (section 10);
- `candidate_degeneracy`: `"not degenerate"` or `"degenerate: <flags>"`;
- `residual_trends`: `"no residual trend detected"` or
  `"unmodelled dependence: <flags>"`.

## 3. The data and the three summaries

**The curves.** For each substrate ("4OMe-BnOH", "BnOH"), take every row of
`scope.frame(scope.archive())` with `live` true, `e0 > 0` and that substrate,
in frame order, with the index reset. No other exclusion. Curves the
activation-sink form cannot hold ARE included.

**The windows.** For a curve with times `t` and duration `T = t[-1] - t[0]`,
a window `(a, b)` holds the readings with
`t[0] + a*T <= t <= t[0] + b*T` (both ends inclusive). The four windows are:

| name | (a, b) |
|---|---|
| `Q3` | (0.50, 0.75) |
| `H1` | (0.00, 0.50) |
| `H2` | (0.50, 1.00) |
| `Q4` | (0.75, 1.00) |

**A window's slope** is the ordinary least-squares slope of `values` on `t`
inside the window: with `tc = t - mean(t)`, `slope = sum(tc*(y - mean(y))) /
sum(tc*tc)`. Its error is `noise / sqrt(sum(tc*tc))`.

**The summaries**, from the observed `corrected` readings and `noise`:

| summary | value | admitted when | se |
|---|---|---|---|
| `L` (level) | `ln(slope_Q3)` | `slope_Q3 > 2*se_Q3` | `se_Q3/|slope_Q3|` |
| `E` (early shape) | `ln(slope_H1/slope_H2)` | `slope_H1 > 2*se_H1` and `slope_H2 > 2*se_H2` | `hypot(se_H1/|slope_H1|, se_H2/|slope_H2|)` |
| `D` (late shape) | `ln(slope_Q4/slope_Q3)` | `slope_Q4 > 2*se_Q4` and `slope_Q3 > 2*se_Q3` | `hypot(se_Q4/|slope_Q4|, se_Q3/|slope_Q3|)` |

- A summary that is not admitted is NaN and is never scored.
- `se` is computed for every curve, admitted or not.
- E is negative on a lag-first curve and positive on a burst-first curve. D is
  negative where the rate declines late in the run.

**A candidate's predicted summaries** use the same windows and the same slope
formula on the candidate's curve `A(t)` at the curve's own `times` (`noise`
irrelevant). If a slope that a summary takes the log of is ≤ 0, the
prediction is `-50.0`.

## 4. The five candidates

**4.1 The shared structure.** Each candidate is a small mechanism for the
catalyst's active fraction θ(t) and the product A(t) it makes. For one curve:

    dθ/dt = q·(1 − θ) − k_r·θ,          θ(0) = θ0
    C0-C3:  dA/dt = V·θ − k·A,           A(0) = 0
    C4:     dA/dt = V·θ·exp(−k·t),       A(0) = 0

- θ0 is the active fraction the catalyst stock brings in, shared by every
  curve of a substrate.
- q and k_r are the forward and back rate constants of activation.
- V is the turnover rate at full activation.
- k is the product loss constant (C0-C3), or the catalyst loss constant (C4).

**4.2 The rate expressions**, for a curve with conditions `s0, h2o2, hoo,
buf, e0, pH`, temperature `T` (°C) and buffer `b`:

    invT  = 1/(T + 273.15) − 1/298.15            R = 8.314462618e-3 kJ/(mol K)
    a     = 1/(1 + 10^(pH − pKa))                 the base-held trap
    q     = 10^lk_f · X · a · exp(−Ea_act·invT/R)
    k_r   = 10^lk_r · exp(−Ea_act·invT/R)
    τ     = 1/(q + k_r)          θss = q/(q + k_r)
    Y     = z/(10^lK_O + z)
    V     = 10^lk_cat[b] · e0 · s0/(10^lK_S + s0) · Y · exp(−Ea_cat·invT/R)
    k     = 10^lk_s · exp(−Ea_s·invT/R)          (C4: 10^lk_d · exp(−Ea_d·invT/R))

On BnOH every `Ea_*` is absent and every `exp(...)` factor is 1.

**4.3 What distinguishes the candidates:**

| id | X (what draws the catalyst active) | z (what the rate saturates in) | loss | reading |
|---|---|---|---|---|
| C0 | 1 | `hoo` | product, `lk_s` | activation unimolecular, held back by base |
| C1 | `buf` | `hoo` | product, `lk_s` | the buffer draws the catalyst active (E + buf ⇌ E*) |
| C2 | `buf` | `buf*hoo` | product, `lk_s` | as C1, and the oxidant is delivered through a buffer perhydrate (step 6b) |
| C3 | `hoo` | `hoo` | product, `lk_s` | HOO- draws the catalyst active (step 4 as the entry) |
| C4 | `buf` | `hoo` | catalyst, `lk_d` | as C1, but the catalyst is lost on a clock instead of the product |

Common to all five:
- one turnover constant per buffer;
- saturation in [S];
- a rate proportional to [enz];
- an activation clock independent of [enz];
- one θ0 per substrate.

These are shared assumptions, not results. `residual_trends` (Task 4) reports
any dependence they miss.

**4.4 The curve, in closed form** (verified against direct ODE integration to
5e-12 on 2026-09-15). For one curve, with `times` as `t`:
- **C0-C3:** `(h, g) = summary_kinetics._activation_sink_columns(t, τ, k)`;
  `A = V*θ0*h[0] + V*θss*g[0]`.
- **C4:** `r = k + 1/τ`; `A = V*(θss*(-expm1(-k*t))/k + (θ0 - θss)*(-expm1(-r*t))/r)`.

**4.5 The parameter vector**, in exactly this order:

1. `lk_cat[<buffer>]`, one per buffer present in the substrate's curves,
   buffers sorted alphabetically;
2. `lK_S`, `lK_O`, `lk_f`, `lk_r`, `pKa`, `theta0`;
3. `lk_s` (C0-C3) or `lk_d` (C4);
4. 4OMe-BnOH only: `Ea_cat`, `Ea_act`, then `Ea_s` (C0-C3) or `Ea_d` (C4);
5. `ls_w_L`, `ls_w_E`, `ls_w_D`, `ls_b_L`, `ls_b_E`, `ls_b_D`.

That is 19 parameters on 4OMe-BnOH and 16 on BnOH.

**Bounds** (L-BFGS-B box bounds):

| parameter | lower | upper |
|---|---|---|
| `lk_cat[*]` | −8 | 2 |
| `lK_S` | −2 | 4 |
| `lK_O` | −8 | 2 (C2: 4) |
| `lk_f` | −12 | 6 |
| `lk_r` | −10 | −1 |
| `pKa` | 4 | 14 |
| `theta0` | 0 | 1 |
| `lk_s`, `lk_d` | −9 | −2 |
| `Ea_cat`, `Ea_act`, `Ea_s`, `Ea_d` (kJ/mol) | 0 | 200 |
| `ls_w_*`, `ls_b_*` | −3 | 1 |

`sigma_w_q = 10^ls_w_q` is summary q's scatter between curves within a run.
`sigma_b_q = 10^ls_b_q` is its scatter between runs.

**4.6 The likelihood.** For summary q, let `e = observed − predicted` over the
admitted curves being scored. For each run g among them:

    d_i   = se_i² + sigma_w²
    A_g   = Σ 1/d_i        B_g = Σ e_i/d_i        C_g = Σ e_i²/d_i
    NLL_g = ½ [ C_g − sigma_b²·B_g²/(1 + sigma_b²·A_g)
                + Σ ln d_i + ln(1 + sigma_b²·A_g) + n_g·ln(2π) ]

This is the exact negative log-density of a normal with a shared random
offset per run (verified against `scipy.stats.multivariate_normal` on
2026-09-15). The total is the sum over q in (L, E, D) and over runs. If the
total is not finite, use `1e12`.

## 5. Assumptions and their stop triggers

| id | assumption | trigger → STOP and report |
|---|---|---|
| A1 | Task 1's anchors reproduce | any Task 1 anchor differs |
| A2 | the candidate code is the plan's | a Task 2 test fails, or any midpoint anchor differs by more than 1e-6 |
| A3 | Stage C is affordable | Task 2's pilot projects more than 12 hours of wall time |
| A4 | the fitter reaches the known optimum on real data | the real 4OMe-BnOH full fit of C1 or of C0 has NLL above 125.0 |
| A5 | real full fits converge | a real full fit has `success` False after its one refit |

## 6. Tasks and stops

| task | content | ends with |
|---|---|---|
| 1 | `data/mechanism_discrimination.py`: summaries, scatter, anchors | commit |
| 2 | candidates, likelihood, fitter, cross-validation, tie rule, degeneracy, tests, pilot | commit (STOP if A2 or A3) |
| 3 | real full fits; planted data from each; planted fits and cross-validation; the discrimination table | **STOP 1** |
| 4 | real cross-validation; tie sets; verdicts; residual trends | **STOP 2** |

---

## 7. Task 1 — the summaries and their scatter

**Goal.** Build the three summaries for every live catalysed curve and
measure their scatter, before any candidate exists.

**Changes** (new module `data/mechanism_discrimination.py`; its docstring
summarises sections 3 and 4 of this plan):
1. **`SUMMARY_WINDOWS`**: the constant `(("Q3", 0.50, 0.75), ("H1", 0.00,
   0.50), ("H2", 0.50, 1.00), ("Q4", 0.75, 1.00))`.
2. **`window_slopes(times, values, noise)`**: `{name: (slope, se)}` for the
   four windows, as section 3.
3. **`curve_summaries(times, values, noise)`**: `{"L", "E", "D", "L_se",
   "E_se", "D_se"}` as section 3, with NaN for a summary not admitted.
4. **`summary_table(substrate)`**: returns
   `{"rows": DataFrame, "times": [array per row]}`.
   - The rows are section 3's curves, in frame order, with columns
     `experiment, sample, buffer, temperature, pH, s0, h2o2, hoo, buf, e0`
     plus the six summary columns.
   - The summaries are read from `CurveFit.corrected` and `CurveFit.noise`.
5. **`summary_scatter(substrate)`**: a DataFrame indexed by `L`, `E`, `D`.
   - **Columns:** `admitted`, `curves`, `q10`, `q50`, `q90` (`numpy.nanquantile`
     over admitted values), `median_se` (over admitted curves), `within_run_sd`,
     `within_run_df`, `replicate_sd`, `replicate_df`, `replicate_n`.
   - **`within_run_sd`:** pooled over runs with ≥ 2 admitted curves,
     `sqrt(Σ_runs Σ (value − run mean)² / Σ_runs (n_run − 1))`, with the
     denominator as `within_run_df`. It still contains every design effect,
     so it is a RAW spread, not a noise estimate.
   - **`replicate_*`:** the same pooling over the rows whose experiment is in
     `scope.REPLICATE_RUNS`, grouped by `sample` instead of by run. NaN when
     there are no such rows (BnOH).

**Tests** (new file `data/test_mechanism_discrimination.py`). Copy the
`FAILURES`/`check` harness and the `__main__` list pattern of
`data/test_rate_laws.py`.
- **`test_window_slopes_on_a_line`:**
  - `t = numpy.arange(0, 3600.0 + 1, 60.0)`, `y = 2e-5*t + 0.01`,
    `noise = 1e-3`;
  - every slope equals 2e-5 to a relative 1e-10;
  - `se_Q3` equals `1e-3/sqrt(sum(tc*tc))` over the 16 readings with
    `1800 <= t <= 2700`, to a relative 1e-12.
- **`test_summaries_are_nan_when_not_admitted`:** the same `t` with `y = 0.01`
  (flat): L, E and D are all NaN.
- **`test_summary_anchors`:** A1 below.

**Anchors A1:**

| substrate | curves | runs | Boric (curves/runs) | Phosphate | Pyrophosphate |
|---|---|---|---|---|---|
| 4OMe-BnOH | 147 | 38 | 40/10 | 92/23 | 15/5 |
| BnOH | 164 | 31 | 27/7 | 19/5 | 118/19 |

`summary_scatter` (q10/q50/q90 to 0.001, SDs to 0.001, counts exact):

| substrate | summary | admitted | q10 / q50 / q90 | median_se | within_run_sd (df) |
|---|---|---|---|---|---|
| 4OMe-BnOH | L | 147 | −13.227 / −10.576 / −8.727 | 0.008 | 0.739 (109) |
| 4OMe-BnOH | E | 145 | −0.623 / −0.025 / 0.539 | 0.004 | 0.279 (107) |
| 4OMe-BnOH | D | 146 | −0.308 / −0.024 / 0.204 | 0.011 | 0.129 (108) |
| BnOH | L | 158 | −14.627 / −11.978 / −10.212 | 0.024 | 0.919 (127) |
| BnOH | E | 158 | −0.867 / 0.224 / 0.781 | 0.012 | 0.810 (127) |
| BnOH | D | 156 | −0.416 / −0.095 / 0.435 | 0.034 | 0.447 (125) |

Replicate (4OMe-BnOH), `replicate_sd (df, n)`: L 0.236 (12, 16); E 0.705 (12,
16); D 0.132 (12, 16).

(`median_se` is anchored to 0.001; the table shows the three-decimal value.)

**Report.** Template E. The tables above as the functions return them.
Verdict: "machinery only".

**Gate.** A1 → STOP. Otherwise `run_gates.py` prints `0 failed`, then commit,
then go straight on to Task 2.

---

## 8. Task 2 — candidates, likelihood, fitter, cross-validation, pilot

**Goal.** Implement section 4 exactly, the fitter, the cross-validation, the
tie rule and the degeneracy flags. Prove them on fixed anchors, and measure
the cost.

**Changes** (`data/mechanism_discrimination.py`):
1. **`CANDIDATES`:** the constant `("C0", "C1", "C2", "C3", "C4")`.
2. **`candidate_names(candidate, substrate, buffers)`:** the section 4.5 list.
3. **`candidate_bounds(candidate, substrate, buffers)`:** `(lower, upper)`
   arrays in that order, from section 4.5.
4. **`candidate_curves(candidate, parameters, rows, times)`:**
   - `parameters` is a `{name: value}` dict; `rows` and `times` come from a
     `summary_table`;
   - returns `{"A": [array per curve], "tau", "theta_ss", "V", "k"}`, each an
     array over curves, from section 4.2 and 4.4.
5. **`candidate_predict(candidate, parameters, rows, times)`:** the predicted
   `{"L", "E", "D"}` arrays (section 3, `-50.0` rule).
6. **`run_likelihood(errors, se, sigma_w, sigma_b, runs)`:** a
   `pandas.Series` of `NLL_g` indexed by run, from section 4.6. It takes
   arrays over the curves being scored.
7. **`candidate_nll(candidate, x, table, keep=None)`:**
   - `table` is a `summary_table` dict and `keep` is a boolean mask over its
     rows (None = all);
   - for each q in (L, E, D), score the rows that are admitted AND kept, with
     `sigma_w = 10**ls_w_q` and `sigma_b = 10**ls_b_q`;
   - returns the total (1e12 if not finite);
   - `buffers` are always the full table's sorted buffers, so the vector
     layout never changes between folds.
8. **`fit_candidate(candidate, table, substrate, keep=None, start=None,
   starts=16, seed=0, maxiter=1000)`:**
   - **One start or many.** If `start` is given, run one
     `scipy.optimize.minimize(..., method="L-BFGS-B", bounds=...,
     options={"maxiter": maxiter})` from it. Otherwise run one from each row
     of `lower + (upper - lower) * qmc.LatinHypercube(d=len(names),
     seed=seed).random(starts)`.
   - **Returns** `{"names", "x", "nll", "success", "nit", "starts": [(nll,
     success) per start], "flags"}` for the lowest NLL, with
     `flags = candidate_degeneracy(...)`.
   - **Refit.** If the lowest start has `success` False, refit ONCE from its
     `x` with `maxiter=3000` and keep the refit, recording `"refit": True`.
9. **`candidate_degeneracy(candidate, fit, table, substrate)`:** the flags,
   in this order:
   - `"not converged"` if `success` is False;
   - `"at bound: <name>"` for each parameter within `0.01*(upper − lower)` of
     either bound;
   - `"clocks outside the run window on <n> of <N> curves"` if more than half
     the curves have `τ > 10*T` or `τ < T/300` (T is the curve's duration).

   Returns `"not degenerate"` or `"degenerate: " + "; ".join(flags)`, plus the
   flag list.
10. **`cross_validate_candidate(candidate, table, substrate, full)`:** for
    each run in sorted experiment order:
    - **skip** the fold and count it if the held-out run has a buffer no
      other run has;
    - otherwise **fit** with `keep = (experiment != run)`,
      `start = full["x"]`, `maxiter=1000`;
    - **score** the held-out run: the `candidate_nll` terms of that run's rows
      only, at the fold's `x`.

    Returns `{"folds": {run: score}, "skipped": count, "fold_success":
    {run: bool}}`.
11. **`candidate_tie(scores)`:**
    - `scores` is a DataFrame, candidates × folds; drop any fold with a NaN;
    - `best` is the lowest total;
    - for each candidate, `d = scores − scores.loc[best]`,
      `mean = d.mean()`, `bar = 2*d.std(ddof=1)/sqrt(folds)`;
    - `status` is `"best"` for the best, `"tied"` if `mean <= bar`, else
      `"excluded"`;
    - also report `worst_fold_difference = d.max()`;
    - returns one row per candidate.

    Like the earlier tie rule, it cannot exclude a candidate that fails on
    only a few folds, however badly. `worst_fold_difference` shows where that
    happens.

**Tests** (`data/test_mechanism_discrimination.py`):
- **`test_candidate_curves_match_the_ode`:** integrate section 4.1 with
  `scipy.integrate.solve_ivp` (`rtol=1e-11`, `atol=1e-16`) on
  `t = numpy.linspace(0, 6000, 101)`, and compare A(t) with section 4.4 at
  t > 0, to a relative 1e-8, for:
  - sink form: V=2e-5, τ=1500, θ0=0.1, θss=0.6, k=2e-4;
  - sink form: V=2e-5, τ=900, θ0=0.9, θss=0.2, k=1e-9;
  - decay form: V=2e-5, τ=1500, θ0=0.8, θss=0.3, k=2e-4.

  In the ODE, write `dθ/dt = (θss − θ)/τ`, which is section 4.1 rearranged.
- **`test_run_likelihood_matches_the_dense_normal`:**
  - `rng = numpy.random.default_rng(0)`; `e = rng.normal(size=5)`, then
    `se = rng.uniform(0.1, 0.3, 5)`, in that order;
  - `sigma_w = 0.4`, `sigma_b = 0.7`, `runs = [1, 1, 1, 2, 2]`;
  - the Series' sum equals the negative `logpdf` of a zero-mean
    `scipy.stats.multivariate_normal` whose covariance is
    `diag(se² + sigma_w²) + sigma_b²·[runs_i == runs_j]`, to an absolute
    1e-10;
  - both equal 3.9404970364222 to 1e-10.
- **`test_midpoint_likelihood_anchors`:** `candidate_nll` at the midpoint of
  every bound (`(lower + upper)/2`), all rows kept, equals A2 to an absolute
  1e-6.
- **`test_candidate_tie_rule`:**
  - folds 0-9; `A = [1..10]`; `B = A + 1.0`; `C = A + [0.6, −0.5]` repeated
    five times;
  - A is `"best"`, B is `"excluded"`, C is `"tied"`.
- **`test_degeneracy_flags`:** on 4OMe-BnOH C1, a fake fit with `x` at the
  midpoint except `Ea_act = 199.0`, and `success` True:
  - the verdict contains `"at bound: Ea_act"`;
  - with `lk_f = −12.0` and `lk_r = −10.0` as well, it contains
    `"clocks outside the run window on 147 of 147 curves"`.

**Anchors A2** (midpoint `candidate_nll`):

| candidate | 4OMe-BnOH | BnOH |
|---|---|---|
| C0 | 4873.435372 | 15622.845649 |
| C1 | 4730.711800 | 15005.891228 |
| C2 | 6297.876950 | 15975.468150 |
| C3 | 3082.772937 | 14998.498241 |
| C4 | 4730.544204 | 15005.932644 |

**Pilot** (not a gate test):
1. Time `fit_candidate("C1", table, "4OMe-BnOH")` with `starts=16`.
2. Time one `cross_validate_candidate` fold for it, the first run not
   skipped.
3. Do the same on BnOH.
4. Project the wall time at 16 worker processes:
   `(10 real full fits + 10 real cross-validations + 50 planted full fits +
   50 planted cross-validations) × measured time / 16`, using each
   substrate's own times.
5. Record the projection.

A review ran a four-start version of this pilot on 2026-09-15: about 60-90 s
per start on 4OMe-BnOH, and best NLLs of 117.9 (C1) and 115.0 (C0).

**Report.** Template E: the changes, the A2 table reproduced, the pilot
timings and the projection. Verdict: "machinery only".

**Gate.** A2 or A3 → STOP. Otherwise `run_gates.py` prints `0 failed`, then
commit, then go straight on to Task 3.

---

## 9. Task 3 — real full fits and the planted discrimination — STOP 1

**Goal.** Fit every candidate to the real summaries. Then, with each fit as a
planted truth, find which candidates the scoring can tell apart on this
archive's own design and scatter.

**Changes** (`data/mechanism_discrimination.py`):
1. **`DISCRIMINATION_DIR`:** the constant
   `os.path.join("data", "fits", "mechanism_discrimination")`.
2. **`real_full_fits(substrate, workers=16)`:**
   - `fit_candidate` for every candidate on `summary_table(substrate)`;
   - saves `DISCRIMINATION_DIR/real_<substrate>_<candidate>_full.json`
     (`names`, `x`, `nll`, `success`, `nit`, `starts`, `refit`, `flags`),
     each as it finishes;
   - reads an existing save back instead of refitting.
3. **`planted_table(substrate, truth, seed=0)`:** a copy of the real
   `summary_table` whose L, E, D are replaced by planted values.
   - `generator = numpy.random.default_rng(seed)`; the truth parameters come
     from the truth's real full-fit save; `pred` is `candidate_predict` at
     them.
   - For q in (L, E, D), in that order:
     - `runs` = the sorted unique experiments among rows where the REAL q is
       admitted;
     - `u = generator.normal(0.0, 10**ls_b_q, len(runs))`, one per run in that
       order;
     - over the admitted rows in table order, `w = generator.normal(0.0, 1.0,
       count) * sqrt(10**(2*ls_w_q) + se_q**2)`;
     - planted q = `pred_q + u[run of the row] + w` on admitted rows, NaN
       elsewhere.
   - `se` columns are unchanged.
4. **`planted_fits(substrate, truth, workers=16)`:** for every candidate on
   `planted_table(substrate, truth)`:
   - `fit_candidate` (16 Latin-hypercube starts; NEVER the truth's `x`), then
     `cross_validate_candidate`;
   - saves `DISCRIMINATION_DIR/planted_<substrate>_<truth>_<candidate>.json`
     (the full fit's fields, plus `folds`, `skipped`, `fold_success`), each as
     it finishes;
   - reads existing saves back.

   Workers receive the table in their payload and never call `scope`.
5. **`discrimination_table(substrate)`:** from the planted saves, run
   `candidate_tie` per truth. Returns:
   - `status`: truth × candidate → `"best"` / `"tied"` / `"excluded"`;
   - `recovered`: per truth, `"truth recovered"` if the truth's own status is
     best or tied, else `"truth not recovered"`;
   - `pairs`, per unordered pair (A, B):
     - `"distinguishable"` if B is excluded under truth A AND A is excluded
       under truth B;
     - `"one-way (<A> as truth excludes <B>)"` if only one of the two holds;
     - `"not distinguishable"` if neither holds;
   - `flags`: truth × candidate → the planted full fit's degeneracy verdict.

**Run, in order:**
1. `.venv/bin/python data/test_mechanism_discrimination.py`: every test
   passes.
2. `real_full_fits` for both substrates. Check A4 and A5.
3. `planted_fits` for every (substrate, truth), background, resumable.
4. `discrimination_table` for both substrates.

**Report.** Template E. Per substrate:
- the real full fits: NLL, `success`, `refit`, degeneracy verdict. **No
  parameter values**; they are read in Task 4 only.
- the status table (truth rows × candidate columns);
- `recovered` per truth;
- the pair table;
- the planted fits' degeneracy verdicts where not `"not degenerate"`.

**Gate.** A4 or A5 → STOP and report. Otherwise `run_gates.py` prints
`0 failed`, then commit, then **STOP 1** (template S).

**Question at STOP 1:** "The discrimination table is above. Run Task 4, the
real cross-validation and verdicts?"

---

## 10. Task 4 — real cross-validation and the verdicts — STOP 2

**Goal.** Score every candidate on held-out real runs. Read each result
against what Task 3 says the design can separate. Report which dependences
the best candidate still misses.

**Changes** (`data/mechanism_discrimination.py`):
1. **`real_cross_validation(substrate, workers=16)`:**
   - `cross_validate_candidate` for every candidate from its real full-fit
     save;
   - saves `DISCRIMINATION_DIR/real_<substrate>_<candidate>_folds.json`, and
     never rewrites the `_full` save;
   - reads existing saves back.
2. **`candidate_verdicts(substrate)`:**
   - `candidate_tie` on the real folds; `A` = the real best; the Task 3 tables
     from `discrimination_table(substrate)`.
   - For each candidate C:
     - C is A: `"best"`;
     - C tied, and under truth A the planted status of C is `"excluded"`:
       `"tied with <A>; the design can separate them"`;
     - C tied otherwise: `"tied with <A>; the design cannot separate them"`;
     - C excluded, and under truth A the planted status of C is `"excluded"`:
       `"excluded against <A>; the planted design reproduces this separation"`;
     - C excluded otherwise: `"excluded against <A>; the planted design does
       not reproduce this separation"`.
   - If C's real full fit is degenerate, append `" (<its degeneracy
     verdict>)"`.
   - Returns the verdict per candidate, plus the tie table with
     `worst_fold_difference`.
3. **`residual_trends(substrate, candidate)`:** on that candidate's real full
   fit, `e = observed − predicted` over admitted rows, per summary q.
   - **Within runs:**
     - regress e on `log s0`, `log buf`, `log h2o2`, after subtracting each
       run's mean from e and from each regressor;
     - drop a regressor whose demeaned SD is below 0.05;
     - OLS with `sigma² = RSS/(n − runs − regressors)` and
       `se = sqrt(sigma²·diag((XᵀX)⁻¹))`;
     - flag `"<q> within runs: <regressor> t=<t:+.2f>"` if `|t| > 3`.
   - **Between runs:**
     - take each run's mean e, and each run's median `log hoo`, `log e0`,
       `pH`, and (4OMe-BnOH only) `invT`;
     - a separate OLS of run-mean e on each axis, with an intercept;
     - flag `"<q> between runs: <axis> t=<t:+.2f>"` if `|t| > 3` and there
       are at least 4 runs.
   - The verdict is `"unmodelled dependence: " + "; ".join(flags)` or
     `"no residual trend detected"`.

**Run, in order:**
1. `real_cross_validation` for both substrates (background, resumable).
2. `candidate_verdicts` for both substrates.
3. `residual_trends` on each substrate's real best, and on every tied
   candidate.

**Report.** Template E. Per substrate:
- the tie table (total fold score, mean difference, bar, status,
  `worst_fold_difference`, skipped folds);
- the verdict per candidate, verbatim;
- the best candidate's parameters, only if its degeneracy verdict is
  `"not degenerate"`, one line per parameter, with no uncertainty. This plan
  computes no standard errors, so quote none;
- the `residual_trends` verdicts.

**Gate.** `run_gates.py` prints `0 failed`, then commit, then **STOP 2**
(template S).

**Question at STOP 2:** "Stage C's verdicts are above. Next: uncertainty for
the surviving candidates' parameters, a revised candidate set built from the
residual trends, or a list of experiments that would separate the survivors?"

---

## 11. Report templates

Fill the brackets. Tables come from named functions. Add no sentences beyond
the template.

### Template E — the DATA_VERIFICATION.md entry

```
## <date> — Stage C Task <n>: <verdict strings, verbatim, or "machinery only">

`PLAN_MECHANISM_DISCRIMINATION.md` Task <n>.

What was asked: <the task's Goal, one sentence>.
What was built: <function names>.
Bugs found: <before/after numbers, or "none">.

<the task's tables>

Verdicts (from <function names>): <verbatim>.
Degeneracy: <every verdict that is not "not degenerate", or "none">.
Could overturn this: <the task's list below>.

Nothing is adopted.
```

"Could overturn this", per task:
- Task 1: the summaries read windows that are fixed fractions of each run, so
  a run's length changes what its windows see; the replicate spread rests on
  one composition, one buffer and one pH, with 12 degrees of freedom.
- Task 2: none.
- Task 3:
  - one planted seed per truth;
  - each truth is a candidate's own fit, so a mechanism outside the five is
    not represented;
  - the tie rule cannot exclude a candidate that fails on a few folds only.
- Task 4:
  - every candidate shares the assumptions listed in section 4.3, and the
    background is assumed to cancel in the reference subtraction (section
    1.2);
  - buffer identity is confounded with pH, and [enz] and temperature move
    only between runs;
  - no run moves [buf] and [H2O2] together;
  - no parameter uncertainty is computed.

### Template S — the stop message (in chat; nothing else)

```
STOP <n> — <task>

Done:
- <at most five bullets, each naming a function or a file>

<the task's tables>

Verdicts: <verbatim>
Degeneracy: <list, or "none">
Could overturn this: <the task's list>

Question: <the stop's question>
```

---

## 12. Amendment 1 — the fixed table, its anchors, and what BnOH may not be read for

Added 2026-09-16, after a review of Task 3's saves. **This section overrides
Task 3's Report and Task 4 wherever they differ.**

### What happened

Task 3's background run finished every fit -- the 10 real full fits and all 50
planted fits are saved -- and then crashed in `discrimination_table` with
`KeyError: 'C0'`. `pd.DataFrame(columns)` builds folds x candidates, while
`candidate_tie` takes candidates x folds, so the tie rule chose a run number
as the best candidate.

**The review fixed it in your working tree** (`data/mechanism_discrimination.py`,
in `discrimination_table`):

    scores = pd.DataFrame(columns).T.sort_index(axis=1)

Do not revert it and do not refit anything. Record it in Task 3's entry under
"Bugs found", with the wording: "`discrimination_table` passed
`candidate_tie` the transpose of its scores table; every fit was already
saved, so only the table was rebuilt."

### Changes

1. **`candidate_verdicts` refuses a substrate whose planted truth was not
   recovered.** Before anything else, read `discrimination_table(substrate)`.
   If any truth's `recovered` is `"truth not recovered"`, return for EVERY
   candidate of that substrate the single string
   `"not readable (planted truth <X> not recovered)"`, naming the first such
   truth, and return no tie table. A scoring that rejects a mechanism on that
   mechanism's own data cannot be read on real data.
2. **Task 4 runs only on a substrate that passes that check.** Skip
   `real_cross_validation` for a substrate that fails it; its folds are not
   worth the compute. On today's saves that means **Task 4 runs on 4OMe-BnOH
   only**.
3. **No parameters are reported.** All ten real full fits carry a degeneracy
   flag, so section 10's "the best candidate's parameters" line is struck.
   Report `"not read (<its degeneracy verdict>)"` in its place.
4. **Task 4's report adds one line per candidate**: its verdict beside the
   pair verdict from `discrimination_table` for that candidate against the
   best.
5. **This plan ends at STOP 2.** The write-up is a new plan, written after
   it. Do not begin one.

### Anchors A6 — the discrimination tables

Rows are the planted truth, columns the fitted candidate. These were
reproduced from the saves with the repository's own `candidate_tie` on
2026-09-16. Any difference → STOP.

**4OMe-BnOH** (every truth recovered):

| truth | C0 | C1 | C2 | C3 | C4 |
|---|---|---|---|---|---|
| C0 | best | tied | excluded | excluded | excluded |
| C1 | tied | best | excluded | excluded | excluded |
| C2 | tied | tied | best | tied | tied |
| C3 | tied | tied | tied | best | tied |
| C4 | tied | excluded | tied | excluded | best |

Pairs: C1 vs C4 `"distinguishable"`; C0 vs C2, C0 vs C3, C0 vs C4, C1 vs C2,
C1 vs C3 one-way with the first as truth; C3 vs C4 one-way with C4 as truth;
C0 vs C1, C2 vs C3, C2 vs C4 `"not distinguishable"`.

**BnOH** (C0 is `"truth not recovered"`; C1-C4 recovered):

| truth | C0 | C1 | C2 | C3 | C4 |
|---|---|---|---|---|---|
| C0 | excluded | tied | tied | tied | best |
| C1 | tied | tied | best | tied | tied |
| C2 | tied | tied | best | tied | tied |
| C3 | tied | tied | best | tied | tied |
| C4 | tied | tied | best | tied | tied |

All ten BnOH pairs are `"not distinguishable"`.

### Run, in order

1. Rerun `discrimination_table` for both substrates and check A6.
2. Write Task 3's entry (template E), including the bug above and the
   `recovered` line for BnOH.
3. `run_gates.py` prints `0 failed`; commit; post **STOP 1**.
4. On the user's go-ahead, Task 4 with changes 1-4.

### Could overturn this

- The tie rule cannot exclude a candidate that fails on a few folds only, so
  `"not distinguishable"` is a statement about this rule at this fold count.
- One planted seed per truth, and each truth is a candidate's own degenerate
  fit, so a different seed or a non-degenerate truth could separate more.

---

## 13. Amendment 2 — C5, the oxidant-depletion candidate — STOP 3

Added 2026-09-16, after STOP 2. **This section overrides sections 4, 9 and 10
where they differ.** It is the last task in this plan.

### Why

- `residual_trends` reports the SAME missing dependence under every candidate,
  the best included: the late shape `D` tracks `log h2o2` WITHIN runs at
  t = +3.2. Within a run is not day-to-day scatter; it is structure. Its sign
  says more peroxide, less late decline.
- Every candidate so far explains the late decline by product loss (C0-C3) or
  catalyst loss (C4). Neither depends on peroxide.
- **The design can separate late-decline mechanisms:** C1 against C4 came back
  `"distinguishable"` in Task 3. So this test has power, which is why it is
  worth one more round.
- **C5 is the minimum test.** It is C1 with the product sink replaced by the
  oxidant draining, so C1 against C5 isolates the late decline and nothing
  else. The catalyst decomposing peroxide is `BUBBLES.md`'s gas reaction.

### The candidate

C5 keeps C1's activation exactly (`X = buf`, the base trap, the per-buffer
turnover constant, the [S] saturation, one `theta0`). It has NO product sink.
Instead the oxidant falls from the start of the run:

    z(t)  = hoo * exp(-k_ox * t)
    Y(t)  = z(t)/(10^lK_O + z(t))
    k_ox  = 10^lk_ox * exp(-Ea_ox*invT/R)
    V_base = 10^lk_cat[b] * e0 * s0/(10^lK_S + s0) * exp(-Ea_cat*invT/R)
    theta(t) = theta_ss + (theta0 - theta_ss)*exp(-t/tau)
    dA/dt = V_base * Y(t) * theta(t),      A(0) = 0

`tau` and `theta_ss` are C1's, unchanged. Note `V_base` EXCLUDES `Y`, which is
now inside the integral.

**The curve** is the running integral on the curve's own times, exactly:

    rate = V_base * Y(t) * theta(t)
    A = numpy.concatenate([[0.0],
        numpy.cumsum(numpy.diff(t) * (rate[1:] + rate[:-1]) / 2.0)])

### Changes (`data/mechanism_discrimination.py`)

1. **`CANDIDATES`** becomes `("C0", "C1", "C2", "C3", "C4", "C5")`.
2. **`candidate_names`** for C5: the C1 list with `lk_s` replaced by `lk_ox`,
   and on 4OMe-BnOH `Ea_s` replaced by `Ea_ox`. **`candidate_bounds`:**
   `lk_ox` is (-9, -2) and `Ea_ox` is (0, 200), as the constants they replace.
3. **`candidate_curves`** gains the C5 branch above. Its returned `k` is
   `k_ox`, and `tau`, `theta_ss` and `V` are as C1 (`V` being `V_base`).
4. **`discrimination_table(substrate, candidates=None)`** and
   **`candidate_verdicts(substrate, candidates=None)`** take the candidate
   tuple to use, defaulting to `CANDIDATES`.
5. **`_real_scores` takes the same candidate tuple**, so a call restricted to
   the five never looks for a BnOH C5 folds save.
6. **C5 is fitted on 4OMe-BnOH ONLY.** Every BnOH call passes
   `("C0", "C1", "C2", "C3", "C4")`. BnOH's refusal from Amendment 1 stands
   unchanged.
7. Nothing else changes: the same tie rule, the same degeneracy flags, the
   same refusal to report parameters from a degenerate fit.

### Tests (`data/test_mechanism_discrimination.py`)

- **`test_c5_matches_the_ode`:** `t = numpy.linspace(0, 6000, 601)`,
  `V_base = 2e-5`, `tau = 1500.0`, `theta0 = 0.1`, `theta_ss = 0.6`,
  `k_ox = 2e-4`, `lK_O = -3.0`, `hoo = 1e-2`. Integrate `dA/dt` above with
  `scipy.integrate.solve_ivp` (`rtol=1e-11`, `atol=1e-16`) and compare with
  the running integral at t > 0, to a relative **1e-4**. (A review measured
  1.9e-5 on this grid on 2026-09-16; the trapezoid is the model, so do not
  tighten this.)
- **`test_midpoint_likelihood_anchors`** gains C5 (A7 below).

### Anchors A7

`candidate_nll` at the midpoint of every bound, all rows kept, to 1e-6:

| candidate | 4OMe-BnOH | BnOH |
|---|---|---|
| C5 | 4829.060949 | 15031.666318 |

**A6 still applies**, checked with the five original candidates:
`discrimination_table(substrate, candidates=("C0", "C1", "C2", "C3", "C4"))`
must reproduce section 12's tables on BOTH substrates. The six-candidate table
is new output and has no anchor.

### Run, in order

1. `.venv/bin/python data/test_mechanism_discrimination.py`: every test passes,
   A7 and A6 included.
2. `real_full_fits("4OMe-BnOH")`. The five existing saves are read back; only
   C5 is fitted. If C5's fit has `success` False after its one refit → STOP
   (A5).
3. `planted_fits("4OMe-BnOH", truth)` for truths `C0`-`C4`: each adds the
   single C5 fit and reads the rest back. Then `planted_fits("4OMe-BnOH",
   "C5")`, which needs step 2's save as its truth and fits all six. Eleven new
   fits in all. Background, resumable. Project the wall time first; A3's
   12-hour bar applies.
4. `discrimination_table("4OMe-BnOH")` with all six, and the A6 re-check with
   the five.
5. `real_cross_validation("4OMe-BnOH")`: only C5's folds are new.
6. `candidate_verdicts("4OMe-BnOH")` with all six, then `residual_trends` on
   the new best and on every tied candidate.

### Report

Template E, plus these, which the entry must state explicitly:
- the six-candidate tie table and the verdict per candidate;
- C5's own verdict, verbatim;
- the pair verdict for **C1 against C5**, verbatim, since that pair is the
  question this amendment asks;
- every new pair verdict involving C5;
- `residual_trends` for the new best: whether `D within runs: log h2o2`
  survives;
- C5's degeneracy verdict.

### Gate

A3, A5, A6 or A7 → STOP and report. Otherwise `run_gates.py` prints
`0 failed`, commit, then **STOP 3** (template S).

**Question at STOP 3:** "C5's result is above. This plan is finished."

**After STOP 3 the plan is closed.** Do not write an analysis folder, an
`ANALYSIS.md`, a figure or any summary document: the write-up is handled
elsewhere. Stop and wait.

### Could overturn this

- C5 keeps every shared assumption of section 4.3, and adds one: the oxidant
  falls as a first-order decay from t = 0, at a rate that does not depend on
  the catalyst loading or on the substrate.
- The archive never measured the peroxide concentration during a run, and
  never measured the gas, so the decay is inferred from the curve shape alone.
- One planted seed per truth, as before.
