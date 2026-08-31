# Kinetic model fitting

What has been fitted, what the fits established, and what they cannot decide.
Companion to `MECHANISM.md` (the chemistry and its reduction),
`DATA_VERIFICATION.md` (the data and every ruling made on it) and
`COMPUTATIONAL.md` (calculations that would settle what the data cannot).

`MECHANISM.md` reduces the 7-step mechanism to 3 ODEs and 4 rate constants and
derives the observation equation. This document covers the implementation of
that reduction, the first fits, and the fact that **the reduced model does not
fit** — together with a reasonably precise account of which layer is at fault.

## Conventions

- **Confidence**: `CERTAIN` (follows from the rate laws or from data alone, no
  fit involved) · `STRONG` (measured, one block) · `PROVISIONAL` (one fit, not
  yet reproduced elsewhere).
- A finding earns a place here only if it **changes what to do next**. Numbers
  that merely describe a fit belong in the saved results file.
- Rate constants are quoted only where they are identifiable. Two of the six are
  not (see [F4](#f4--two-of-six-constants-are-lower-bounds-not-values)), and
  quoting those as values would be a mistake, not a rounding.
- Units throughout: concentrations mM, time s. Hence `k_can` mM⁻² s⁻¹, `k3`,
  `k0` and `k6` mM⁻¹ s⁻¹, `k5'` mM⁻² s⁻¹, `r` dimensionless.

## The code

| module | what it is |
|---|---|
| `data/kinetic_model.py` | the reduced system as ODEs. No I/O, no fitting, no data |
| `data/fit_dataset.py` | which curves are fittable, and what each one is |
| `data/fit_kinetics.py` | the sequential fit, its diagnostics, and the `r` profile |
| `data/plot_fit.py` | a saved fit drawn against the curves it was fitted to |
| `data/test_kinetic_model.py` | 24 checks on the chemistry and the solver |
| `data/test_fit_kinetics.py` | selection counts, and parameter recovery |
| `data/fits/*.json` | saved results — a fit costs ~30 min, the plots seconds |

```bash
python data/fit_kinetics.py --list                          # blocks with both stages
python data/fit_kinetics.py --substrate BnOH --save data/fits/BnOH_25C_Phosphate.json
python data/plot_fit.py data/fits/BnOH_25C_Phosphate.json   # -> figures/
```

## Findings register

| id | finding | confidence | what it changes |
|---|---|---|---|
| [F1](#f1--the-model-is-first-order-in-substrate-and-the-data-is-not) | the model is first order in `[S]`, the data is ~half order | **CERTAIN** | the reduction needs a saturable substrate term. Highest leverage |
| [F2](#f2--the-enzyme-free-limit-is-a-fixed-point) | the `E0 = 0` limit needs a third constant, not two | **CERTAIN** | `k0` added; stage 1 fits four parameters |
| [F3](#f3--the-observation-equation-only-rescues-the-lag-if-r--1) | `signal = A + r·BA` can only lag if `r > 1` | **CERTAIN** | the observation equation inherits the falsification |
| [F4](#f4--two-of-six-constants-are-lower-bounds-not-values) | `k3` and `k5'` saturate out of the observable | **CERTAIN** | never quote them as values |
| [F5](#f5--the-reduced-model-misfits-by-20-24x-the-noise) | 20–24× the curves' own noise on the fitted block | **STRONG** | the reduction, not necessarily the chemistry |
| [F6](#f6--at-its-own-best-fit-the-mechanisms-machinery-is-inert) | `[PBA] ≈ 10⁻⁹ mM` at the best fit | **PROVISIONAL** | the fitted model is just a linear seed |
| [F7](#f7--only-two-of-eleven-blocks-support-the-sequential-fit) | 2 of 11 blocks can be fitted sequentially | **CERTAIN** | names the missing experiment |

---

## F1 — the model is first order in substrate, and the data is not

**This is the finding with the most leverage, and the only one that needs no fit
at all.**

At `t = 0` with `A = PBA = 0`, every term in the reduced system vanishes except
the seed:

```
d(signal)/dt |_{t=0}  =  v_seed(0)  =  (k0 + k5'·E0)·[H2O2]·[S]
```

which is *exactly* proportional to `[S]`. No choice of the six rate constants
changes that — `v_can` carries no `[S]`, `v3` and `v_seed` carry it linearly, and
nothing divides by it. Measured across the fitted block the model's effective
order is **+1.01**, which is just the algebra confirming itself.

The data does not do this. Measured as the slope of log rate against log `[S]`
*within* each experiment, so pH, buffer, temperature and `[H2O2]` are all held
fixed:

| | experiments | data | model |
|---|---|---|---|
| `[buf]` held constant | 8 | **+0.30** (range 0.08–0.49) | **+1.01** |
| `[buf]` varied along the ladder | 2 (exps 3, 6) | −0.23 | +1.03 |

**The second row must not be read as a substrate order.** Exps 3 and 6 are
buffer titrations in which `[buf]` falls 85 → 25 mM as `[sub]` rises 1.28 → 8.98
mM — the `[buf]`/`[sub]` collinearity recorded in `DATA_VERIFICATION.md` since
2026-08-29. Exp 3's rate *falling* monotonically with substrate is a buffer
effect wearing substrate's clothes. The two series are reported separately and
never pooled; `plot_fit.py` hatches them.

On the eight clean experiments a 20× substrate range produces roughly a 3×
signal change — exp 67 runs 0.0115 → 0.0345 AU net for `[S]` 0.36 → 7.31 mM,
monotonic and clean, an order of 0.37.

Three things worth stating about this:

- **It is not substrate depletion.** Conversion at the lowest `[S]` in exp 67 is
  2.6%; across the whole archive the median is 0.9%.
- **It holds for the enzyme-free runs too** (exps 67, 69, 70: +0.49, +0.37,
  +0.28), so it is not only the catalyst that saturates. Whatever the background
  reaction is, it is not first order in alcohol either.
- **Saturation is expected here on independent grounds.** The Bols group reports
  `Km = 1.25 mM` for benzyl alcohol (`MECHANISM.md` item 4) and these ladders
  span 0.21–14.31 mM, straddling it. An order near +0.3 is what
  Michaelis–Menten gives across such a range.

The reduction has no saturable `[S]` term anywhere: stage 3 makes the catalyst
states algebraic and stage 4 drops the remaining denominator. **This is a
structural deficiency of the reduced model, not a bad fit.** It is the first
thing to change.

## F2 — the enzyme-free limit is a fixed point

`MECHANISM.md` writes the `E0 = 0` limit as collapsing to "2 ODEs and 2
parameters (`k_can`, `k3`)". That is true of the algebra and false of the
trajectory: with `A = PBA = 0` at `t = 0`, `v_can` goes as `[A]²` and `v3` as
`[PBA]`, so both are zero and the system never starts. Integrated to `t = 10⁵` s
it returns `A = PBA = 0` exactly.

Seeding with a trace of aldehyde does not rescue it. Steps 1–2 destroy two
aldehydes per peracid and step 3 returns one, so the catalyst-free loop is a net
aldehyde **sink**: `A(0) = 10⁻⁴` mM decays, reaching `[PBA] = 1.2 × 10⁻¹²` and
`[BA] = 1.0 × 10⁻⁹` mM after 10⁵ s. This is `MECHANISM.md`'s own "step 5 is the
only net source of aldehyde" followed to its conclusion.

Since the enzyme-free controls demonstrably react, an **E0-independent source**
is required: the uncatalysed `S + H2O2 → A`, i.e. the `E0 → 0` limit of step 5.
The seed is written `(k0 + k5'·E0)·[H2O2]·[S]`, which keeps the system exactly
linear in `E0` as the reduction promises, with an intercept rather than through
the origin. **Stage 1 fits three rate constants plus `r`, not two.**

## F3 — the observation equation only rescues the lag if r > 1

`MECHANISM.md` derives `dA/dt ≤ v5(0)` — the pure-aldehyde reading cannot show a
lag — and proposes `signal = [A] + r·[BA]` as the fix. The same bound survives
the fix for every `r ≤ 1`:

```
d(signal)/dt = v_seed + (1 + r)(v3 + v6) − 2·v_can
dPBA/dt      = v_can − (v3 + v6) ≥ 0        while peracid accumulates
⟹  d(signal)/dt ≤ v_seed + (r − 1)·v_can ≤ v_seed(0)
```

and `v_seed` only ever decreases, since `S` and `H2O2` only deplete. Checked
numerically as well: a random search over **227 parameter sets** spanning eight
orders of magnitude in every rate constant, at both `E0 = 0` and `E0 > 0`,
produced **not one** accelerating curve at `r ≤ 1`. At `r > 1` acceleration
appears at once, up to 4.7× the initial slope.

So reproducing any induction period requires `eps(benzoate) > eps(benzaldehyde)`
at 285/300 nm, against `MECHANISM.md`'s band-shape bracket of `r ~ 0.08–0.33`.

**And the fitted model errs the other way.** At its best-fit `r = 1.52` it lags
in 19 of 23 enzyme-free curves and 19 of 20 catalysed ones, where the block's
data lags in 7 of 43 (16%). Both directions are closed off:

| `r` | model lags in | data |
|---|---|---|
| `≤ 1` | no curve, at any rate constants | — |
| `1.52` (best fit) | 19/23 and 19/20 | 7/43 |

Whether some intermediate `r` reproduces 16% is what `--profile-r` exists to
answer and has not been run.

*Re-measured statistic:* `MECHANISM.md` reports 52% of curves reaching peak slope
past 15% into the run (n = 326). Over the 404 curves the fitting code selects, by
the same smoothed method, it is **34%** (136/404). The selections differ — the
n = 326 predates the carbonate rule and the exclusions of exps 50, 64 and 85.

## F4 — two of six constants are lower bounds, not values

Both `k3` and `k5'` stop affecting the observable above a threshold, so the
curves bound them from below and say nothing above.

**`k3` (step 3, uncatalysed peracid oxidation).** Step 3 consumes the peracid
steps 1–2 make. Once `k3` is large enough that `[PBA]` reaches quasi-steady
state, `v3 → v_can`; benzoate is produced at `v3`, so the observable is set by
`v_can` alone and `k3` has dropped out of it. Raising `k3` further only lowers
standing `[PBA]` in exact inverse proportion, which nothing measures.

| change to `k3` | effect on the signal | `[PBA]` |
|---|---|---|
| ×1000 | **0.39%** | ÷1000 (to within 2%) |
| ÷100 | **25%** | ×~65 |

**`k5'` (step 5, the catalysed seed).** The same phenomenon one step round the
loop: once the seed is fast enough that starting the loop is not what limits the
observable, making it faster changes almost nothing. Cost against the true value:

| `k5'` | ×0.1 | ×1 | ×100 |
|---|---|---|---|
| cost | 1.3 × 10³ | 1.5 × 10⁻⁵ | 7.6 |

Steep below, flat above. This is consistent with `MECHANISM.md`'s instinct that
`k5` be "treated as free/unconstrained", and extends the same verdict to `k3`. It
also explains the fitted block's `corr(k_can, k3) = −0.979`: near the
quasi-steady-state threshold the two trade off.

**So only `k_can`, `k0`, `r` and `k6` are determined.**

## F5 — the reduced model misfits by 20-24x the noise

BnOH / 25 °C / phosphate, the one block where the sequential strategy is fully
supported. Both stages converged with no parameter at a bound.

| stage | curves | fitted | rms residual |
|---|---|---|---|
| 1, enzyme-free | 23 / 5 exps | `k_can` 10.86, `k3` 23.8, `k0` 2.14 × 10⁻⁹, `r` 1.523 | 0.0147 AU = **24.0×** noise |
| 2, catalysed | 20 / 5 exps | `k5'` 1.4 × 10⁻⁷, `k6` 4.4 × 10⁻³ | 0.0140 AU = **20.6×** noise |

Per-curve noise here is 0.0003–0.0006 AU, so this is a decisive failure rather
than a rough fit. On the worst curves the model produces a fraction of the
observed signal: exp 69 sample 1 rises +0.047 AU where the model gives +0.001.

The plots make the pattern legible. On exp 6 — 300 min at pH 6.71 — the model
tracks the data well. On exps 67/69/70 at pH 8.01 over 20–30 min it is
essentially flat while the data climbs to 0.02–0.047 AU: the model's only route
to signal is the pH-independent seed `k0`, so it has no fast pH-dependent path
and cannot produce anything in twenty minutes.

Stage-1 correlations, which is where the remaining degeneracy sits:

| | k_can | k3 | k0 | r |
|---|---|---|---|---|
| **k_can** | 1.000 | **−0.979** | −0.341 | 0.048 |
| **k3** | −0.979 | 1.000 | 0.376 | −0.115 |
| **k0** | −0.341 | 0.376 | 1.000 | **−0.902** |
| **r** | 0.048 | −0.115 | −0.902 | 1.000 |

Condition number 523. Note `k_can` and `r` are *not* degenerate here. On
synthetic curves at a **single** pH they come out −0.999 correlated at condition
number 6 × 10⁷, because they enter the observable almost purely as the product
`k_can·r`; pH enters only through `[HOO⁻]`, which multiplies `k_can` and not `r`,
so varying pH separates them. The real block spans pH 6.71–8.01 and does.

## F6 — at its own best fit, the mechanism's machinery is inert

At the fitted constants `[PBA] ≈ 10⁻⁹ mM` — nine orders below the aldehyde pool.
The Cannizzaro/peracid apparatus contributes essentially nothing, and the model
reduces to linear aldehyde accumulation from the uncatalysed seed `k0`, with the
autocatalytic loop switched off. Given F2 this is what it must do: the loop is a
net aldehyde sink, so the optimiser turns it down as far as it can.

Stage 2 tells the same story from the other side — it drives `k5' → 1.4 × 10⁻⁷`
and carries the catalysed route entirely on `k6`. Given F4 the honest reading is
that these runs give **no evidence for a seed step at all**, not that the seed
constant is small.

## F7 — only two of eleven blocks support the sequential fit

Rate constants may be pooled only within one (substrate, temperature, buffer)
cell: temperature moves every constant through Arrhenius, the two substrates are
different molecules, and `MECHANISM.md`'s buffer section argues at length that
the four buffers are chemically different reagents rather than four ways of
setting pH.

| cell | `E0 = 0` | `E0 > 0` | sequential fit |
|---|---|---|---|
| BnOH, 25 °C, phosphate | 23 | 20 | **yes** — the block fitted here |
| 4OMe-BnOH, 40 °C, phosphate | 59 | 4 | yes, but almost no catalysed data |
| 4OMe-BnOH, 25 °C, phosphate | 8 | 52 | weak background |
| BnOH, 25 °C, boric | 4 | 28 | weak background |
| BnOH, 25 °C, pyrophosphate | 0 | **127** | none at all |
| the other six | 0 | 4–40 | none at all |

The archive's largest catalysed block — 127 BnOH pyrophosphate curves — has no
enzyme-free counterpart, so its background constants could only be imported from
a different buffer, which is exactly what the buffer section forbids.
**The missing experiment is an enzyme-free control in pyrophosphate**, and it is
cheap.

## Why the fitter is believable

A fitter that has never recovered a known answer is evidence about the
optimiser, not about the mechanism. Synthetic curves are generated from known
constants at the real experimental conditions, with the block's real noise
(0.0006 AU against curves of 0.004–0.12 AU), and the fitter has to find them
back:

| | `k_can` | `k3` | `k0` | `r` |
|---|---|---|---|---|
| noiseless | 0.00 dec | 0.00 dec | 0.00 dec | 0.00 |
| at instrument noise, seeds 7/11/23 | −0.02/−0.10/+0.16 | +0.10/+1.06/−0.97 | 0.00/+0.01/+0.01 | +0.03/+0.01/−0.13 |

`k3` swings by two decades, which is F4 showing up rather than a fitter fault.

**Four bugs were found this way, all of which would have biased every fit**, and
none of which would have been visible from real data alone — on real curves each
just looks like model error:

| bug | how it showed | caught by |
|---|---|---|
| model not baseline-subtracted like the data | cost **2.03** at the *true* parameters, not 0 | noiseless recovery |
| optimiser slid along the `k_can`/`r` valley | `k_can` three decades out, `r` pinned at its bound | synthetic recovery |
| blind multi-start reported a bad local minimum as **converged** | cost 4.07 × 10⁷ where 7.67 × 10³ was available — 5300× | hand-checking the cost at `k5, k6 → 0` |
| screening then discarded the good starts | noiseless recovery went from exact to +4.47 decades | the same recovery test |

The third is the most dangerous: "converged, no parameter at a bound" is exactly
what a reader would trust. The fix was to screen starting points — a 64-point
Latin hypercube over the bounded box, ranked by one cheap residual evaluation
each — and the fourth bug was that fix discarding the nominal start and the `r`
ladder whenever random points happened to score better. **A low cost *at* a point
is a poor predictor of where an optimiser starting there converges.** Both
mechanisms now compose: guaranteed starts always run, screening only adds.

Suites: `test_kinetic_model` 24/24, `test_fit_kinetics` 39/39.

## What to do next

1. **Add a saturable substrate term and refit.** F1 is the highest-leverage
   finding, is independent of any fit, and is expected on independent grounds
   from Bols's `Km = 1.25 mM`. A small change to `kinetic_model.rhs`.
2. **Fit the 4OMe-BnOH / 40 °C block.** 59 enzyme-free curves, more background
   data than the block fitted here. If the substrate-order gap reproduces there,
   F1 stops being a claim about five experiments.
3. **Run `--profile-r`.** Cheap, and it settles whether any `r` reproduces the
   observed lag fraction or whether F3 closes the door completely.
4. **Measure `eps(benzoate)` at 285 nm** — `COMPUTATIONAL.md` task C1, or one UV
   spectrum. `r` is fitted at 1.52 and the spectroscopy says 0.08–0.33; that
   contradiction is currently absorbed into a fitted parameter.
5. **Consider whether the observable includes `[PBA]`.** Perbenzoic acid is the
   one species in the mechanism whose concentration genuinely accelerates, and
   no extinction coefficient for it has been sought. One line in
   `kinetic_model.observable`.

## What this does not license saying

The evidence bears on the **reduction** and the **observation equation** far more
than on the chemistry:

| layer | verdict |
|---|---|
| the reduction (stages 1–4 of `MECHANISM.md`) | **clearly deficient** — no saturable `[S]` term, and the enzyme-free limit is a fixed point as written |
| the observation equation | **suspect** — `r ≤ 1` cannot lag, fitted `r > 1` over-lags |
| the 7-step chemistry | **not yet convicted** |

Stage 4 of the reduction explicitly discarded the saturation denominators, so F1
is at least as likely to be an artefact of the reduction as a fault in the
mechanism. One block has been fitted. The 4OMe block, the `r` profile, and a
saturable variant should all be done before the chemistry is blamed.

## Log

### 2026-08-31 — the fitter, and the first fits

Implemented the reduction as `kinetic_model.py` and fitted BnOH/25 °C/phosphate
sequentially. Findings F1–F7 above; full working, including the four bugs and
the two corrections to earlier claims, in `DATA_VERIFICATION.md` under the same
date. Figures rebuild from `data/fits/BnOH_25C_Phosphate.json` in seconds.

Two claims made and withdrawn during the work, recorded so they are not
resurrected: that `k_can` and `r` are non-identifiable (true only at a single
pH — F5), and that the observed lag fraction is 52% (34% on this selection —
F3).
