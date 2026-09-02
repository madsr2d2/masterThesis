# The temperature series: exps 14-19

4OMe-BnOH + H<sub>2</sub>O<sub>2</sub>, chemzyme-catalysed, in 65 mM phosphate
at pH 7.00 with 82.5 mM H<sub>2</sub>O<sub>2</sub>, run at **15, 20, 25, 30, 35
and 40 °C**. The same four-rung substrate ladder (1.850, 3.700, 5.549,
7.399 mM) in every one, so the block is four independent Arrhenius fits over
six temperatures, not one fit over 24 points.

**It is the only temperature series in the archive.** Nothing in the BnOH set
varies temperature at all, so this is the sole route to activation parameters
and the only kinetic quantity that can be compared against the ORCA barriers in
`COMPUTATIONAL.md`.

    scope.TEMPERATURE_SERIES        the six experiments, in temperature order
    data/arrhenius.py               the fits
    python data/verify_enzyme_stock.py --sequence
    python temperature_series/build_figures.py
    python temperature_series/check_numbers.py

**Figures**: `index.html` is the presentation — twelve figures, A to L, one per
claim in this document. `progress_curves.html` carries all 24 curves with the
burst/lag fit and the breakpoint drawn on each. Both are rebuilt by
`build_figures.py`, which computes nothing: every number in a figure comes from
`arrhenius`, `scope`, `curve_metrics` or `verify_enzyme_stock`, so a figure and
this document cannot disagree without `check_numbers.py` saying so.

| exp | 14 | 15 | 16 | 17 | 18 | 19 |
|---|---|---|---|---|---|---|
| temperature | 25 °C | 35 °C | **40 °C** | 30 °C | 20 °C | 15 °C |
| `[enz]` mM | 0.273 | 0.273 | **0.241** | 0.273 | 0.273 | **0.241** |
| curves, live | 4 | 4 | 4 | 4 | 4 | 4 |

All 24 curves are live and **24 of 24 come from the instrument's own `.rre`**.

That was not true until 2026-09-02. Six of them — the 5.549 mM rung in every one
of the six runs — were being read from the `.txt` export at 1096× the noise
floor, and the cause was a case-sensitive regex: the instrument wrote sample 3's
label as `sample003`, lowercase, in **31 files across exps 1–32**, and
`read_rre` matched only `Sample00\d`. Nothing reported a missing curve because
the export exists and was silently used instead. It was found from a plot — one
panel per run looked coarser than its neighbours. All 28 recovered cuvettes
agree with their own export within its rounding step; none failed. See
`DATA_VERIFICATION.md` 2026-09-02.

## 1. The enzyme mismatch is real, not a typing mistake

Two of the six runs record an enzyme concentration **11.7% lower** than the
other four, and one of them — exp 16, the 40 °C run — sits *between* two runs at
the higher value. That had to be settled before any fit, because the two
readings call for opposite treatments: a real restock is divided out, a
transcription error is corrected.

**Three independent lines agree that the values are real.**

### The weighing chain, which the pipeline never reads

`[enz]` reaches `experiment_data.csv` from the sheet's `kuv` cell. Beside it
each sheet also records the preparation that produced that cell:

    mass / molar mass / volume        ->  stock concentration
    stock x Enz[ml] / Vol[ml]         ->  kuv

Four numbers written at four different moments, none of them the one the
pipeline reads. `data/verify_enzyme_stock.py` recomputes `kuv` from the mass and
the cuvette volumes for every experiment in the archive:

| stock | weighed | made up to | cuvette mM | experiments |
|---|---|---|---|---|
| 1 | 0.0122 g | 3.3 mL | 0.17533 | 2-9 |
| 2 | 0.023 g | 4 mL | 0.272695 | 10-15, **17, 18** |
| 3 | 0.0203 g | 4 mL | 0.240683 | **16**, 19-22, 32, 34-36 |
| 4 | 0.0228 g | 4 mL | 0.270324 | 37, 41-49 |
| 5 | 0.0236 g | 4 mL | 0.167885 | 50-59, 85 |
| 6 | 0.0236 g | 40 mL | 0.0279809 | 60-84 |
| 7 | 0.0027 g | 4 mL | 0.0320121 | 127-131 |

**62 of the 63 chains are self-consistent** and **63 of 63 compiled values equal
the sheet's own `kuv`**. The single break is exp 58, which is already excluded
for running backwards, and its `kuv` misses its own preparation by 0.5%.

So both values in question descend from a real weighing — 0.0203 g against
0.023 g — each consistent to six figures, and stock 3 is used by **nine**
experiments. A transcription error into one cell would break the chain. Neither
chain is broken.

### But the ordering is genuinely odd, and the workbooks cannot settle it

Stocks are made once and used until exhausted, so in experiment order the mass
should change rarely and never change back. **Exp 16 is the only experiment in
the entire campaign whose stock interrupts a run** (`--sequence`): stock 3
between exp 15 and exp 17, both stock 2.

That is exactly the shape a copied workbook leaves, and the file evidence is
ambiguous. Exps **16, 17, 18 and 20 are byte-identical in size (92672) and
structure** — one lineage — and it straddles the boundary: 16 and 20 carry
stock 3, 17 and 18 carry stock 2. Whichever direction the copying went, someone
edited an enzyme block; the documents cannot say who was edited into what.

The instrument files add nothing here. The `.rre` binaries store an internal
path, and it does catch two copies elsewhere in the archive — exp 003's file is
internally named `mads_t006.rre` and exp 029's is `mads_t023.rre` — but exp
016's is named `mads_t016.rre`. The `161008` field is a firmware build code,
identical in all 143 files, not a date.

### The kinetics decide it

Rate divided by enzyme concentration should fall on a straight line against
1/T. A concentration that is wrong for one run displaces that run's four points
together, off the line the other five temperatures define. So: divide by each
candidate and ask which is most collinear (`arrhenius.enzyme_hypotheses`).

On `vmax`, residual rms in ln units — a relative scatter:

| hypothesis | mean rms | worst rung | rungs won | E<sub>a</sub> kJ/mol |
|---|---|---|---|---|
| **as recorded** | **0.078** | **0.085** | **4 of 4** | 90.1 |
| exp 16 restocked to 0.273 | 0.108 | 0.115 | 0 of 4 | 87.5 |
| exps 16 and 19 both 0.273 | 0.130 | 0.143 | 0 of 4 | 90.2 |

**The recorded values win every rung.** On `v0_quad` as a sensitivity they win
every rung again — 0.236 against 0.262 and 0.275 — so the verdict does not
depend on the estimator, and the four rungs are independent fits rather than
one result counted four times.

Per run it is clearer still. Exp 16's mean distance from the line:

- as recorded (0.241): **−0.060** — about 6% low, inside the others' scatter
  (exp 15 is +0.048, exp 17 +0.029)
- forced to 0.273: **−0.121** — twice as far out, and further out on all four
  rungs

Forcing the higher concentration makes exp 16 a **worse** outlier, in the
direction a too-large divisor produces.

### Verdict

**Real. Use the recorded values and divide by them.** Exp 16 and exp 19 were run
with a genuinely different enzyme stock, and the fits below normalise by
`[enz]` per run.

**What this does not establish**: *why* exp 16 sits out of sequence. Most likely
it was re-run later with the 19-22 session and kept its design number, but
nothing in the archive records that, and it does not change the treatment. Note
also that the method is weakest against an error common to *both* stock-3 runs —
displacing 15 °C and 40 °C the same way tilts the line hardly at all — so
"exps 16 and 19 both wrong" is rejected less firmly than "exp 16 alone", though
it loses every rung too.

## 2. The breakpoint screen: the breaks are the induction, not an artefact

Run first, because exp 65 in the BnOH background showed that the start-versus-end
shape statistics step straight over a mid-run break, and every conclusion here
rests on a single rate per curve. `scope.synchronised_break`:

| T | 15 °C | 20 °C | 25 °C | 30 °C | 35 °C | 40 °C |
|---|---|---|---|---|---|---|
| slope after / before | 1.50-18.08 | 1.33-2.27 | 0.99-1.81 | 1.27-1.36 | 0.89-1.00 | 0.70-0.86 |
| cuvettes steepening | 3 of 4 | 2 of 4 | 1 of 4 | 0 of 4 | 0 of 4 | 0 of 4 |

**Nothing here looks like exp 65.** That run's four cuvettes broke within 56 s of
each other regardless of a 20-fold substrate range, which is the signature of
something that is not the reaction. Here the ratio falls with temperature — cold
runs accelerate throughout, hot runs decelerate — and a second, unrelated
statistic says the same thing: every run from 15 to 30 °C classifies as a **lag**
curve, with τ falling **6489 → 3190 → 945 → 916 s** as it warms, and 35 and 40 °C
flip to `burst`. Two independent statistics agreeing on an induction period that
shortens with heating is chemistry, not an artefact.

**The fall is not quite monotone, and the exception is instructive.** Median
ratios run 1.74, 1.56, 1.19, **1.33**, 0.91, 0.82 — 30 °C sits above 25 °C. The
break ratio does not track temperature directly, it tracks **how much of the run
the induction occupies**, and the runs are not the same length: measured in units
of its own τ, the 25 °C run is **19τ** long and the 30 °C run only **5.5τ**. The
quantity that *is* monotone in temperature is τ itself.

**A limitation of the statistic, seen here for the first time.** Exp 19's
lowest-substrate cuvette reports a ratio of **18.08**, far outside everything
else. It is a near-zero denominator: that curve is flat at 5.6e-08 AU/s for the
first 6076 s — the coldest temperature at the lowest substrate, the slowest
condition in the block — and then rises at 1.0e-06. The reading is "the
pre-break slope is indistinguishable from flat", not "the post-break slope is
remarkable". `break_ratio` is unstable wherever `slope_before` approaches zero,
and a lag curve at the slowest condition is exactly where that happens.

### The screen itself looked for one break, and there are two

Added 2026-09-02, after the same observation that produced §3a: a rate that
rises and then falls has **two** breakpoints, and `segmented_fit` searches for
one. It lands on whichever change is stronger and never reports the other —
which is exactly how the early break on the 40 °C curves went unseen while the
late one was written up.

`curve_metrics.segment_breaks` now searches one **or** two, exhaustively, on
straight-line errors precomputed from prefix sums so a two-break search costs
O(1) per candidate pair rather than refitting lines 68 000 times.
`segment_selection` chooses between them on an F test — three extra parameters,
one breakpoint and one line — with the same nested discipline as the fitting
forms.

Exp 16's rise-then-fall rungs now report the pair: **490 and 3626 s** (slopes
2.39 → 3.31 → 2.81 ×10⁻⁵) and **539 and 3871 s** (2.93 → 4.19 → 3.53).

**Read the pattern, not the count.** A progress curve is smooth, and three
straight lines fit a smooth bend better than two whether or not anything
happened, so the F test takes a second breakpoint on **20 of 24** — including
the 15 °C curves whose slopes simply rise 0.08 → 0.15 → 0.17. The *count* is a
statement about piecewise-linear approximation; the *sequence of slopes* is the
statement about the curve, and only "rise then fall" needs a maximum in the
middle.

| T | 15 °C | 20 °C | 25 °C | 30 °C | 35 °C | 40 °C |
|---|---|---|---|---|---|---|
| rising | 4 | 4 | 1 | 4 | 0 | 0 |
| rise then fall | 0 | 0 | **3** | 0 | **2** | **2** |
| falling | 0 | 0 | 0 | 0 | 2 | 2 |

**This is an independent check on §3a.** The model-free description says
"rising" at exactly the temperatures where the nested fit selects one phase
(15, 20, 30 °C) and "rise then fall" or "falling" at exactly those where it
selects two (25, 35, 40 °C). Two unrelated procedures — piecewise lines and a
nested exponential fit — partition the block the same way.

### But it exposed a real problem with `vmax`

`vmax_where` — where the steepest block sits, as a fraction of the run:

| T | 40 °C | 35 °C | 30 °C | 25 °C | 20 °C | 15 °C |
|---|---|---|---|---|---|---|
| `vmax_where` | 0.10-0.30 | 0.29-0.49 | 0.48-0.87 | 0.30-0.50 | 0.70-0.90 | 0.70-0.90 |

At 15 and 20 °C the rate is **still rising when the run ends**, so `vmax` there
is not a maximum — it is wherever the measurement stopped. Under-reading the
cold end tilts the Arrhenius line and inflates E<sub>a</sub>.

**The two estimators fail at opposite ends**, which is the awkward part. On a lag
curve the burst form's `v_ss` is the rate `vmax` is trying to reach, and at
15-30 °C that fit is sound (tau resolved on 16 of 16 curves, residuals 0.94-1.64x
noise). At 35–40 °C it degenerates — τ unresolved on 7 of 8, the form flipping
to `burst`, tau unresolved on 7 of 8, because a decelerating curve has no lag to
measure — so `v_ss` there
is a meaningless late rate, and `vmax` is the good one.

**How big is the bias?** `arrhenius.truncation_sensitivity` substitutes `v_ss`
for `vmax` at 15 and 20 °C only, purely to size it:

| | E<sub>a</sub> kJ/mol |
|---|---|
| `vmax` throughout | 90.1 |
| `v_ss` at the cold end | 88.4 |
| **inflation from truncation** | **1.7** |

`v_ss`/`vmax` is 1.042 at 15 C and 1.035 at 20 C — a 3.5–4% under-read — and
0.98, 0.99 at 25 and 30 C where the two should agree, which is the check that
the substitution means anything. (It is 0.67 and 0.29 at 35 and 40 °C, the
degenerate end, and is not used there.)

**So truncation costs about 1.9 kJ/mol, against a spread of 8.5 across the four
substrate rungs that should agree.** It is real, it is in a known direction, and
it is *not* the limiting problem. Mixing estimators is normally the error rather
than the fix, and it is used here only to size the bias, never for a headline.

## 3. What is actually going on in these runs

Three things, and they have to be separated before any activation energy means
anything.

### The reaction has an induction period, and it is what the breaks are

Every run from 15 to 30 °C is a **lag** curve and τ falls monotonically as it
warms: **6489 → 3666 → 945 → 876 s**. By 35 °C the induction is over before the
measurement is properly under way and the curves decelerate instead, which is
why the burst form flips to `burst` there and τ stops being resolved. So the
block is not six measurements of one shape; it is a shape *changing with
temperature*, and which end of a curve is informative changes with it.

### The substrate ladder is not a substrate ladder

`[buf]` falls **80 → 70 → 60 → 50 mM** as `[S]` rises **1.850 → 7.399 mM**,
because substrate volume displaced buffer volume in the cuvette. It is the same
volume-displacement confound as the BnOH titrations (exps 3 and 6) and the 4OMe
background block.

**It leaves the activation energies alone.** The gradient is identical in all
six runs, so each rung is one *fixed composition* measured at six temperatures —
which is exactly what an Arrhenius fit needs. What it corrupts is the substrate
order.

### The catalysed buffer order, and it is not one number

The temperature series cannot measure it — its own `[buf]` moves only as a
by-product. Exps **32 and 34** can: 4OMe-BnOH at 40 °C and pH 7.00 with `[S]`
held at 8.251 mM and `[buf]` stepped. This is the design the BnOH set lacks
entirely.

| range | experiment | order in `[buf]` | R² |
|---|---|---|---|
| 50–200 mM | 32 | **+0.400 ± 0.028** | 0.990 |
| 3.125–25 mM | 34 | +0.803 ± 0.173 | 0.915 |

**The buffer dependence saturates** — about first order at low buffer, about
half order above 50 mM. The temperature series sits at 50–80 mM, so +0.40 is
the value that applies to it. (`v_ss` cannot be used here: three of exp 34's
four curves have an unresolved τ and return a negative `v_ss`, giving a
meaningless −0.18 ± 0.34.)

### So the substrate order is about +0.58, not +0.45

With `[S]^a [buf]^d` and `log[buf] = g·log[S] + c`, a fit without a buffer term
returns `a' = a + d·g`. Here **g = −0.325** and **d = +0.400**, so the
correction adds 0.13:

| T | 15 °C | 20 °C | 25 °C | 30 °C | 35 °C | 40 °C |
|---|---|---|---|---|---|---|
| observed | +0.460 | +0.521 | +0.490 | +0.349 | +0.438 | +0.420 |
| corrected | +0.590 | +0.651 | +0.620 | +0.479 | +0.569 | +0.550 |

Mean **+0.577**, and it does not drift with temperature — the spread, 0.173, is
the same before and after correcting, so the correction moved the level and not
the trend. **Genuine partial saturation in substrate, with no detectable
temperature dependence**, which is the same as saying K<sub>M</sub> is roughly
constant over 15–40 °C.

*What this rests on*: that the buffer order measured at 40 °C holds at 15 °C,
and at 1.85–7.4 mM substrate as well as at the 8.251 mM it was measured at.
Neither is testable here.

### And the four rungs share one activation energy

The rung-to-rung spread of 8.5 kJ/mol looked like a composition dependence until
the standard errors were computed — each rung is **± 2.6–3.0**. Against their
weighted mean, **χ² = 4.91 on 3 degrees of freedom, reduced χ² = 1.64**. No
composition dependence is detectable, so the rungs can be pooled.

They are pooled by refitting — one slope with a free intercept per rung, over
all 24 points — rather than by averaging the four fits, because those four share
the same six runs, the same six days and the same six cells, so their errors are
correlated and averaging would shrink the error by a factor that is not there.

## 3a. The fitting form: one relaxation, or two

The burst/lag form is `A = c + v_ss·t − B(1 − e^(−t/τ))`, so its rate is
`v_ss − (B/τ)e^(−t/τ)` — **monotone**. It approaches `v_ss` from below (lag,
B > 0) or from above (burst, B < 0) and never turns. **It cannot represent a
rate that rises to a maximum and then falls**, and 14 of these 24 curves do
exactly that: exp 16's upper rungs peak at 1666–2303 s and end about 30% below
their peak.

`summary_kinetics.fit_two_phase` adds a second relaxation:

    A = c + v_ss·t − B₁(1 − e^(−t/τ₁)) − B₂(1 − e^(−t/τ₂))

**B₂ = 0 recovers the one-phase form exactly**, so the two are properly nested
and the choice is a plain F test on two extra parameters
(`summary_kinetics.fit_progress`). Given both time constants the model is linear
in the other four, so the pair is profiled on a log grid rather than handed to
an optimiser — the same reason `fit_burst` profiles τ, and it removes local
minima from the search.

**Selection is on F alone, and deliberately not on τ₂ being resolved.** The
first version required both and rejected the clearest two-phase curves in the
block — 40 °C at 5.549 mM, F = 791, residual falling from 3.81× noise to
1.05× — because τ₂ is 7562 s in a 6517 s run and anything above 2545 s fits as
well. Whether a phase **exists** and whether its time constant is **pinned** are
different questions; answering the first with the second is the error exp 65
taught. `TwoPhaseFit.resolved` still reports the second, and τ₂ is resolved on
only 2 of the 11 two-phase curves, so **τ₂ is not quoted anywhere**.

**Selected on 11 of 24**, and on exactly the curves the one-phase form was
failing:

Residuals below are medians over the run's four cuvettes, in units of each
curve's own noise.

| T | two-phase | one-phase | after selection |
|---|---|---|---|
| 15 °C | 0 of 4 | 1.03× | 1.03× |
| 20 °C | 0 of 4 | 0.98× | 0.98× |
| 25 °C | **4 of 4** | 1.99× | **1.15×** |
| 30 °C | 0 of 4 | 1.09× | 1.09× |
| 35 °C | **4 of 4** | 1.49× | **1.13×** |
| 40 °C | **3 of 4** | 2.97× | **1.04×** |

Nothing is selected where the one-phase form already sits at noise. Every
selected fit is `lag then fall` bar one.

**The second phase is descriptive and must stay that way**, even now that §6
has named the process behind it. τ₂ is not resolved on 9 of the 11 curves, and
`fit_two_phase`'s second exponential is not the sink's rate constant in any
case: the sink is a loss term on the accumulated signal, `A′ = v − kA`, and its
k is read off the rate-against-product line in §6, not off τ₂. An activation
energy taken from τ₂ would be an activation energy of a fitting parameter.

### And it changes which rate to put on an Arrhenius plot

Not `v_ss`. In the two-phase form that is the t → ∞ asymptote, extrapolated far
outside the data and unconstrained once a decay is in the fit: it comes out
**negative** on two of the 35 °C curves and gives an Arrhenius scatter of
0.19–1.18 in ln units against 0.08 for the alternative.

**`v_peak`** — the maximum of the *fitted* rate — is the observable. It is
defined the same way on both forms: `v_ss` on a one-phase lag, where the rate
rises monotonically to it; the interior maximum on a two-phase curve. It is what
`vmax` measures off the raw readings **with the truncation removed**, because
the fit knows the shape and does not have to stop where the readings do.

## 4. Activation parameters

`arrhenius.parameter_table()`. Three fitted quantities, each pooled over the
four rungs with a free intercept each, Eyring fitted on the same design.

| parameter | E<sub>a</sub> kJ/mol | ΔH‡ kJ/mol | ΔS‡ J/mol/K | ΔG‡(298) kJ/mol | curves |
|---|---|---|---|---|---|
| `v_peak`, peak rate of the fitted model | 88.8 ± 1.8 | **86.3 ± 1.8** | **−57.7 ± 5.9** | **103.5 ± 0.1** | 24, 6 T |
| `vmax`, steepest observed rate | 90.1 ± 1.5 | 87.6 ± 1.5 | −53.5 ± 5.0 | 103.6 ± 0.1 | 24, 6 T |
| `v_ss`, asymptote after the induction | 89.0 ± 2.6 | 86.6 ± 2.6 | −56.5 ± 9.0 | 103.4 ± 0.1 | 16, 4 T |
| `1/τ`, induction rate constant | 95.0 ± 15.7 | 92.6 ± 15.7 | +3.7 ± 53.2 | 91.5 ± 0.6 | 16, 4 T |

In kcal/mol for `v_peak`: ΔH‡ **20.6**, ΔS‡ **−13.8 cal/mol/K**, ΔG‡(298)
**24.7**.

**Three routes now agree on ΔH‡ ≈ 86 kJ/mol.** `v_peak` gives **86.3 ± 1.8** on
all six temperatures; `v_ss` gives **86.6 ± 2.6** on the four where the
one-phase form is sound; and `vmax` gives 87.6 ± 1.5, which §2 measured as
running about 1.7 high from truncation — i.e. **85.9** corrected. The two-phase
form reaches the same answer *without* the estimator substitution, which is what
it was built for.

`vmax` keeps the tighter error and the better-agreeing rungs (reduced χ² 1.64
against **3.99** for `v_peak`), because a model-based peak inherits the fit's
scatter and the τ grid is coarse. So: **quote `v_peak`, and note that `vmax`
agrees once its known bias is removed.**

**The load-bearing check is that estimators with different failure modes
agree.** `vmax` uses all six temperatures and is truncated at the cold end;
`v_ss` uses only the four where the one-phase form is sound and is not
truncated; `v_peak` uses all six and is not truncated because it is read off a
fitted shape rather than off the readings. They come out at 87.6, 86.6 and
86.3 kJ/mol, and the truncation correction predicted `vmax` would read 1.7 high.
Three routes with unrelated weaknesses landing within their errors is the best
internal evidence this block offers.

**ΔS‡ ≈ −55 J/mol/K is strongly negative** — a well-ordered transition state, as
an associative step should be. Note what it is quoted *at*: an intercept is a
level, the four rungs sit at four levels, and this is the **median rung,
[S] = 5.549 mM and [buf] = 60 mM**. A ΔS‡ with no composition attached is not
comparable to anything.

**The induction is a different process.** Its ΔG‡ is **12 kJ/mol lower** than the
turnover's, so it is the faster of the two, and its ΔS‡ is near zero rather than
−55 — a much looser transition state. Its enthalpy, 92 ± 16, is indistinguishable
from the turnover's 88 ± 1.5, but with an error that wide that is not evidence
of anything. Treat the induction's numbers as an order-of-magnitude statement:
the *entropies* differ, the enthalpies are not resolved apart.

### The fall has a name now, and it does not move these numbers

`product_fate/` identifies the fall of §3a as **the oxidant attacking the
aldehyde it has just made**: the rate declines linearly in the accumulated
product, `A′ = v − kA`, on 24 of 29 well-determined curves against 0 for the
hyperbolic law product inhibition would give. That has an immediate consequence
here — **every rate in the table above is already net of that loss**, so none of
them is the production rate.

Recovering it is not as simple as adding `k·A` back at the peak. On the four
coldest runs the fitted rate peaks at the *end* of the run rather than inside
it, so that recipe quietly swaps `v_peak` — an asymptote, and truncation-free by
construction — for a value read at the last reading, which is exactly the
cold-end truncation `v_peak` exists to avoid. It moves E<sub>a</sub> by
**+2.4 kJ/mol of pure bookkeeping**. What is done instead is to fit the sink
model itself with its rate constant **pinned** at the value
`slowdown.sink_activation` gives for the temperature, so the production rate is
a parameter rather than a reconstruction (`slowdown.production_frame`).

| estimator | E<sub>a</sub> kJ/mol | ΔH‡ kJ/mol | ΔS‡ J/mol/K | ΔG‡(298) kJ/mol | rms |
|---|---|---|---|---|---|
| `v_peak`, the table above | 88.77 ± 1.77 | 86.28 | −57.7 ± 5.9 | 103.47 | 0.099 |
| `v_prod`, sink model, k pinned | 85.82 ± 3.40 | 83.32 | −67.2 ± 11.3 | 103.34 | 0.189 |

The production rate is **6.7% above** `v_peak` on the median curve and between
−6.6% and +12.7% across the six temperatures, **with no order to it** — and a
factor that does not order in temperature cancels out of a slope. E<sub>a</sub>
moves by **−2.96 ± 3.83 kJ/mol**, which is nothing, and the level shift is worth
**+0.54 J/mol/K** on ΔS‡ against its ±5.9 error.

**So the table stands**, and now with a bound on this rather than a hope. It is
reported and not applied: `v_prod` is the more nearly correct quantity and the
noisier one — Arrhenius scatter 0.189 against 0.099, because the sink shape is
one the four coldest runs cannot test — so quoting it would trade a bias
smaller than the error for a variance larger than it. `product_fate/` §5 has
the derivation and what it rests on.

### Before comparing any of this to a calculation

- **These are composite, not elementary.** `vmax` and `v_ss` are whole-turnover
  rates through a seven-step mechanism; ΔH‡ is the barrier of whatever step is
  rate-limiting at that composition, plus any pre-equilibria in front of it.
  `COMPUTATIONAL.md`'s barriers are for single steps. The comparison to make is
  against the **highest** step, not against a sum.
- **ΔS‡ assumes a pseudo-first-order rate constant.** The rate constant used is
  `v / (ε·[enz])`, in s⁻¹, which is what makes an entropy possible at all — a
  slope survives any constant factor, an intercept does not. But the substrate
  order is +0.58, not 0, so the reaction is not saturated and this constant
  carries a substrate dependence. ΔH‡ does **not** depend on this; ΔS‡ and the
  absolute ΔG‡ do.
- **It also assumes first order in enzyme, and the archive's one test of that
  disagrees.** This block cannot test its own assumption — `[enz]` takes two
  values here, 11.7% apart. Exps 59 and 60 can: the only pair in the archive
  differing in `[enz]` and nothing else that matters, **0.028 against 0.014,
  exactly 2.000×**. A Selwyn test on them
  (`scope.selwyn_test`, `scope.catalyst_order`) gives an order in catalyst of
  **+0.34**, not 1.

  If that transferred here, **ΔH‡ would not move at all** — it is a slope, and
  no constant factor touches it — but ΔS‡ would go **−53.5 → −46.4 J/mol/K** and
  ΔG‡(298) **103.6 → 101.4 kJ/mol**. It may well not transfer: exps 59/60 are
  BnOH in boric at pH 8.51, and this block's own normalisation prefers the
  opposite — `[enz]^1` gives an Arrhenius rms of **0.078** against **0.109** for
  `[enz]^0.34`, a 40% worse fit. The lever is only 11.7%, so that is weak, but
  it points the other way. **ΔH‡ is the number to put beside a calculation;
  ΔS‡ and the absolute ΔG‡ carry an unresolved ±7 J/mol/K and ∓2 kJ/mol on top
  of their fit errors.**
- **The curves are background-subtracted** (`with_E`, reference omits the
  enzyme), so these are the activation parameters of the *catalytic increment*,
  not of the overall reaction. The uncatalysed background has its own
  temperature dependence and it is not measured anywhere in the archive.
- **ΔG‡ is much better determined than either ΔH‡ or ΔS‡** — ± 0.1 against
  ± 1.5 and ± 4.9 — because both come from one line and their errors are
  strongly anti-correlated. That is real, not a mistake: at 25 °C, inside the
  measured range, the errors cancel. Quote ΔG‡ with confidence and the split
  with care.

## 5. Activation energy per rung

From `arrhenius.rung_fits()` on `vmax`, kept because the agreement between them
is the evidence that pooling is allowed:

| `[S]` mM | `[buf]` mM | E<sub>a</sub> kJ/mol | rms |
|---|---|---|---|
| 1.850 | 80 | 94.5 ± 2.8 | 0.078 |
| 3.700 | 70 | 86.0 ± 3.0 | 0.085 |
| 5.549 | 60 | 88.4 ± 2.7 | 0.074 |
| 7.399 | 50 | 91.5 ± 2.7 | 0.076 |

Weighted mean **90.2**, reduced χ² **1.64** — see §3. The pooled refit over all
24 points gives **90.1 ± 1.5**.

**A caution for whatever comes next.** 20 of these 24 curves clear the 3sigma
acceleration gate and 8 exceed z = +15, so `v0` is an *induction* rate
here and must not be used — `vmax` is the estimator throughout. And none of
these curves has been through `curve_metrics.segmented_fit`; exp 65 in the BnOH
background taught us that the start-versus-end statistics step straight over a
mid-run break.

## 6. What the slowdown is — and where it is written up

The fall of §3a is not a property of this block, so its analysis is not in this
document. It is **[`product_fate/`](../product_fate/ANALYSIS.md)**, because the
discrimination that settles it needs the whole 4OMe archive: these 23 curves
cannot separate "the rate fell with time" from "the rate fell with product" —
inside one curve the product only ever grows with time — and the separation
comes from 84 curves in which run length and product vary independently.

What it concludes, in one paragraph. **The catalysed 4OMe rate falls in
proportion to the product it has made** (−0.525 ± 0.071 here on product against
+0.118 ± 0.122 on run length; −0.919 ± 0.161 on the comparison that holds every
condition fixed). **The same chemistry with no enzyme falls on a clock instead**
(−0.361 ± 0.047 on run length, +0.025 ± 0.048 on product) — and that is the
reference channel these runs are measured against, so its decay is subtracted
out. **The rate is linear in the product, not hyperbolic in it**, which is
`A′ = v − kA`: production minus a first-order loss of the measured species, not
inhibition of the catalyst. The stationary level it implies carries substrate
order +0.610 ± 0.067 against the +0.577 §3 measured on the rates, and the
implied selectivity is `k_A/k_S ~ 54`.

And it does not happen to BnOH, at product concentrations nearly twice as high —
`MECHANISM.md`'s S2, switched on by the methoxy group. The effect on the
activation parameters is in §4 above: none that this block can resolve.
