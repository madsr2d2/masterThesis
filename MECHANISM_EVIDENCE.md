# What the progress curves decide about the mechanism

The evidence register for the mechanism question. `MECHANISM.md` holds the
proposed chemistry and `FITTING.md` the record of fitting it; this document
says which parts of it the **kinetics archive can and cannot decide**, and what
would decide the rest. Written 2026-09-16, after the three-stage effort of
`PLAN_CURVES_TO_MECHANISM.md` and `PLAN_MECHANISM_DISCRIMINATION.md`.

`test_mechanism_evidence.py` re-derives every number below from the modules.

## The four statements

1. **One mechanistic candidate is excluded by the curves**: peroxide drawing
   the catalyst into its active form. It agrees with the pH evidence in
   `induction/`, reached independently.
2. **Nothing else is separated.** Five other candidates tie, on a design that
   demonstrably cannot tell most of them apart.
3. **The ceiling is reproducibility, not statistics.** Repeat runs of one
   composition differ by more than the differences between mechanisms.
4. **The strongest mechanistic statements in this project need no fit at
   all** — they are the sign and bound tests in `induction/`, `ph/` and
   `buffer/`.

## 1. The material

| | 4OMe-BnOH | BnOH |
|---|---|---|
| live catalysed curves | 147 | 164 |
| runs | 38 | 31 |
| curves the activation-sink form cannot hold | 42 | 52 |

**217 of the 311 live catalysed curves carry an activation-sink fit** and 94
do not (`activation_sink.remaining_curves`). A catalysed curve is the
catalytic INCREMENT: its reference cuvette holds the same mixture minus the
enzyme, so everything here describes the catalyst's own contribution.

## 2. The reproducibility ceiling

`scope.REPLICATE_RUNS` (exps 2, 4, 5, 7) is the archive's only four-fold
repeat of one composition. Pooled by cuvette, in natural logs
(`rate_laws.replicate_parameter_scatter`):

| quantity | SD |
|---|---|
| the fitted form's `v_act` | **0.642** |
| the fitted form's `k_act_lag` | **0.665** |
| `lag_half_s`, the same curves | 0.250 / 0.374 |
| `vmax_corrected`, the same curves | 0.279 |

Two consequences, and they govern everything below.

- **A quantity read through the fitted form scatters about twice as much as
  the same curves' window-free summaries.** Part of what looks like chemistry
  in a fitted parameter is the form trading its parameters against each other.
  This rests on 6-7 degrees of freedom at one composition, one buffer and one
  pH.
- **Nothing can predict a curve better than the curve repeats.** On the
  summaries Stage C scores (`mechanism_discrimination.summary_scatter`), the
  replicate SDs are 0.236 (level), 0.705 (early shape) and 0.132 (late shape).
  The early shape — which is where an activation clock shows itself — is the
  least reproducible of the three.

Against that, the rate laws fitted to the curve parameters leave residuals of
0.596 to 1.200 in log units (`rate_laws.law_scatter`), where those parameters'
own errors are 0.014 to 0.212, and much of the residual is WITHIN runs (0.371
to 0.978). The conditions do not explain the curve-to-curve variation.

## 3. What rate laws bought, at the readings

Stage B fitted the activation-sink form to the readings of every curve at
once, with its parameters replaced by rate laws. Costs are summed squared
residuals in units of the noise, about 1 per curve at the noise floor
(`rate_laws.law_free_baselines`):

| model | 4OMe-BnOH | BnOH |
|---|---|---|
| each curve's own activation-sink fit | 3,709 | 763 |
| each curve's own quadratic, no laws | 13,117 | 1,387 |
| the rate laws, cross-validated | 56,225 | 6,080 |
| each curve's own line and ONE global sink constant, no conditions | 58,959 | 6,195 |
| Stage A's laws plugged in unchanged | 2,499,344 | 48,071 |

The verdicts are **"laws add less than 10% over the law-free baseline
(4.6%)"** and **"laws add less than 10% over the law-free baseline (1.9%)"**.

And the fits that produced them are not readable as rate laws:
`rate_laws.stage_b_degeneracy` returns **"degenerate: clocks outside the run
window on 86 of 101 curves"** on 4OMe-BnOH and **"degenerate: clocks outside
the run window on 66 of 102 curves"** on BnOH, each with flat coefficients
listed. The optimiser switched activation off rather than place a clock per
curve.

**Why this is a result and not a failure.** It measures how much of a progress
curve the conditions actually determine, once the curve's own size is not
given away for free.

## 4. What mechanism discrimination decided

Six candidate mechanisms, each a small ODE for the catalyst's active fraction
and the product, scored on three window-free summaries with a per-run random
offset, and cross-validated by leaving out whole runs
(`PLAN_MECHANISM_DISCRIMINATION.md` section 4):

| id | what it says |
|---|---|
| C0 | activation is unimolecular, held back by base |
| C1 | the buffer draws the catalyst into its active form |
| C2 | as C1, with the oxidant delivered through a buffer perhydrate |
| C3 | HOO⁻ draws the catalyst into its active form |
| C4 | as C1, with the catalyst lost on a clock instead of the product |
| C5 | as C1, with the oxidant draining during the run |

On 4OMe-BnOH (`mechanism_discrimination.candidate_verdicts`, 38 folds, none
skipped):

| candidate | fold total | status |
|---|---|---|
| C1 | 208.9 | best |
| C0 | 211.5 | tied |
| C4 | 237.3 | tied |
| C5 | 237.8 | tied |
| C3 | 308.7 | **excluded** |
| C2 | 340.5 | tied |

**C3's verdict is "excluded against C1; the planted design reproduces this
separation".** That is the one mechanistic exclusion the curves earned, and
`induction/` reaches the same conclusion from the clock's pH sign without
fitting anything.

**BnOH says nothing.** Its verdict is **"not readable (planted truth C0 not
recovered)"**: on data generated by C0 itself, C0 was excluded. A scoring that
rejects a mechanism on that mechanism's own data cannot be read on real data.

**Most pairs cannot be separated even in principle here.** Planting each
candidate in turn and scoring all six gives, among others, **"not
distinguishable"** for C2 against C3, C2 against C4 and C2 against C5, and
only C1 against C4 comes back **"distinguishable"**. Two cautions:

- the planted truth is only ever *tied* on its own data, never best, and a
  wrong candidate often takes "best" — these mechanisms mimic one another
  through the summaries;
- a pair's verdict depends on which other candidates are in the comparison.
  C0 against C1 read "not distinguishable" with five candidates and "one-way
  (C1 as truth excludes C0)" with six.

**No parameter values are quoted from any of this.** Every real fit carries a
degeneracy flag, most of them clocks outside the run window.

## 5. The live lead

Under every candidate, including the best, the late shape's residual tracks
peroxide WITHIN runs: `mechanism_discrimination.residual_trends` returns
**"unmodelled dependence: D within runs: log h2o2 t=+3.22"** for C1. Within a
run is not day-to-day scatter.

Neither product loss (C0-C3), nor catalyst decay (C4), nor first-order oxidant
depletion (C5) accounts for it — C5 was built for exactly this and ties rather
than wins, and the trend survives in its own residuals. The sign says more
peroxide, less late decline. `BUBBLES.md`'s gas reaction and `MECHANISM.md`'s
S4 are where to look next.

## 6. What would settle the open questions

Nothing in the archive separates C0 from C1, or C1 from C2. Four experiments
and two calculations would:

1. **A run crossing [buf] with [H2O2] at one pH.** Of 88 runs, 53 step
   `[buf]`, 20 step `[H2O2]` and **0 step both**
   (`induction.peroxide_crossing`), which is why a general base cannot be told
   from a buffer perhydrate.
2. **An enzyme-free control in pyrophosphate.** The archive's largest
   catalysed block has none, so its background constants cannot be fitted
   (`FITTING.md` F7).
3. **A catalyst ladder inside one run.** `[enz]` moves only between runs, so
   the prediction that the activation clock does not depend on it is untested.
4. **Replicates at more than one composition.** The whole ceiling in section 2
   rests on four runs at one composition.
5. **`COMPUTATIONAL.md` C7** — the hydrate pKa and dehydration barrier — is
   the non-experimental route to C0 against C1, and it has a measured barrier
   to reproduce.
6. **`COMPUTATIONAL.md` C9** — phosphate against pyrophosphate perhydrate —
   is the same for C1 against C2.

## What this does not license saying

- **Not** that C1 is the mechanism, or that the buffer activates the catalyst.
  C1 predicts held-out runs best; C0, C2, C4 and C5 tie with it.
- **Not** that a tie supports a candidate. Most ties here are the design
  failing to separate, which the planted tables measure directly.
- **Not** any rate constant, barrier or binding constant from Stage B or
  Stage C. Every fit behind them is flagged.
- **Not** anything about BnOH from Stage C.
- **Not** that the activation clock is absent. It is not identifiable from
  these curves, which is a different statement.
