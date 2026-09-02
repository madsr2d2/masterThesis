# The background reaction: BnOH + H<sub>2</sub>O<sub>2</sub> without catalyst

What the uncatalysed oxidation of benzyl alcohol by hydrogen peroxide does, how
much of its rate belongs to the **buffer** rather than to the substrate, and what
that means for reading the catalysed progress curves of exps 135–151.

Figures: [`index.html`](index.html) (the argument),
[`progress_curves.html`](progress_curves.html) (all 27 cuvettes and every fit).
Rebuild both with `python background_reaction/build_figures.py`.

Every number below is re-derived and checked by
`python background_reaction/check_numbers.py`, which fails if this document and
the code disagree.

---

## 1. Why this was needed

The in-scope curves were recorded *against* the background, so the background has
to be understood before they can be read. The first obstacle was that many
enzyme-free runs varied buffer concentration at the same time as substrate, which
confounds the two.

## 2. What the archive actually holds

83 enzyme-free curves over 21 experiments. On BnOH there are six usable runs, in
two designs that do different things:

| set | experiments | curves | `[buf]` | `[sub]` |
|---|---|---|---|---|
| `scope.BUFFER_FIXED` | 67, 69, 70 | 12 | held at 85 mM | ladder |
| `scope.BUFFER_CONFOUNDED` | 3, 6 | 11 (10 live) | falls 85 → 25 mM | ladder, rising |

All but one of these 27 cuvettes now read from the instrument's own `.rre`. Until
2026-09-01 exps 3 and 6 were entirely on the `.txt` export, whose 0.001 AU
rounding is 1096× coarser — see `DATA_VERIFICATION.md`. The buffer order moved
by 0.03 when they were corrected, which is the best evidence available that it
does not rest on the readings' resolution.

Exp 64 is excluded by the manifest (aborted run, 7 minutes at dt = 28 s).

**No enzyme-free run in the archive varies `[buf]` at constant `[sub]`, in any
block.** The five that once did — exps 32, 34–37, a 3.125–200 mM titration at
fixed substrate — were ruled *catalysed* on 2026-08-31 from their
reference-channel layout, so their curves are catalytic increments, not
backgrounds. `DATA_VERIFICATION.md` records that withdrawal.

Two facts about the instrument matter throughout. Every run is double-beam, and
**what the reference channel omits is what the reported curve is net of**
(`kinetics_io.py`, `verify_enzyme.py`):

- enzyme-free sheets — reference omits the **H<sub>2</sub>O<sub>2</sub>** → the curve is the raw background;
- catalysed sheets — reference omits the **enzyme** → the curve is a catalytic *increment*, already net of the background at matched buffer, substrate, peroxide and pH.

## 3. Choosing a rate estimator

This was settled empirically, not by preference, because the answer had to be
shown not to depend on it.

| estimator | window? | whole curve? | form |
|---|---|---|---|
| `vmax` | steepest of four 20% blocks | no | straight line |
| `v0` | first 20% | no | straight line |
| **`v0_quad`** | **none** | **yes** | **A = c + v₀t + at²** |
| `v0_whole` | none | yes | straight line |
| `v0_burst` | none | yes | A = c + v_ss t − B(1 − e^(−t/τ)) |

`v0_quad` is the headline. It answers the objection that `INITIAL_WINDOW = 0.20`
is arbitrary — it chooses no window, uses every point, and being linear in its
parameters is always identified.

### Why not the burst/lag v₀, which is the better *shape*

Considered on 2026-09-01 and declined, on measurement rather than taste. The
burst form is the only one here that is a shape rather than a polynomial, and
where the deceleration really is exponential it is decisively better — exp 67
sample 2 goes from 3.06× noise to 0.94×, exp 69 sample 2 from 3.47× to 1.48×.
But it is not better everywhere, and it fails worst exactly where the quadratic
was objected to:

| | quadratic | burst/lag |
|---|---|---|
| clearly better (> 0.1× noise) | 5 curves | 6 curves |
| indistinguishable | 16 curves | |
| median residual | 1.18× noise | 1.10× noise |
| v₀ defined on | 27 of 27 | 24 of 27 bounded |

**Exp 65 is the case that decides it.** All four of its cuvettes are fitted
badly by both forms, and *worse* by the burst: 6.9, 8.3, 8.0 and 6.9× noise
against the quadratic's 5.2, 5.7, 5.2 and 4.8×. The reason is worth stating,
because the panel does not show it unless you read the numbers: on those four
the burst fit has **collapsed to a straight line** — B = 0 exactly, τ pinned at
the floor of its grid (≈3 s), `tau_resolved` False — while still reporting
`bounded = True`. Fitting with `constrain="none"` returns the identical fit, so
the 3σ lag gate is not the cause. A four-parameter form that has degenerated
into a two-parameter one, and reports its v₀ confidently, is the worst thing to
put in a headline column. On exp 65 samples 3 and 4 it returns +2.6e−5 and
+1.6e−5 — the whole-curve average slope — where the quadratic returns −2.0e−6
and −5.3e−6, which is the curve's genuine flat start being reported as flat.

**And the result does not depend on the choice**, which is the point of the
table in §5: `v0_burst` gives **+1.19 ± 0.25** against `v0_quad`'s +1.30 ± 0.25,
inside one standard error, and the two agree exactly at +1.37 once the
accelerating curves are dropped. Dropping the whole boric block (exp 65) moves
it to +1.20 ± 0.28. Nothing about first order in buffer turns on it.

So `v0_burst` is carried as a fifth estimator in the robustness tables, and
every panel now prints both forms' residuals in units of the curve's own noise
with `— DOES NOT FIT` above 3×. The real finding behind the objection is not
that the wrong estimator is headlined; it is that **exp 65 is described by
neither form**, which is now stated on its panels rather than left to the eye.

**The burst/lag form was tested and rejected as a rate source.** Fitting
A = c + v<sub>ss</sub>t − B(1 − e<sup>−t/τ</sup>) and taking v₀ = v<sub>ss</sub> − B/τ
is the natural whole-curve alternative, and it is what `summary_kinetics.fit_burst`
implements. On these curves it is not identified:

- Unconstrained, **4 of 27 curves return a negative v₀** where the line fit is firmly positive. Exp 67 sample 3 gives −2.07e-04 against a line's +3.35e-06, with a profile interval, [−3.72e-04, −1.34e-04], that never reaches zero.
- The existing `resolved` flag does not catch it: exp 67 sample 3 is `resolved=True`. `resolved` asks whether **τ** is located, which is a different question from whether **v₀** is — and it errs both ways. Exp 6's four cuvettes are `resolved=False` with v₀ pinned to a 0.00 half-width, because as τ → ∞ the curve is a straight line, which kills τ and B but leaves v₀ → v_ss exact.
- **A close fit does not license the extrapolation.** Exp 3 sample 7 is fitted to within its own noise (rms/noise 1.00), yet across values of τ that are statistically indistinguishable — the rms does not move at two decimal places — its v₀ ranges from **−1.06e-05 to +5.20e-07**, a change of sign, against a line fit of +4.01e-07.
- Shutting the lag branch (**B ≤ 0**) removes every negative v₀ and raises the bounded count from **16 to 25 of 27**. It does not solve the problem: it trades the τ → ∞ degeneracy for τ → 0, and exps 65 samples 1 and 2 come back at **10× and 28×** their line rate.
- **A blanket B ≤ 0 is not justified, and was corrected on 2026-09-01.** It rested on "0 of 16 curves accelerate", which was measured on the constant-buffer runs and generalised to all 27 without checking (that count is now 1 of 16, exp 67 sample 3 having crossed the gate when the first reading was dropped). Exps 3 and 6 hold **four curves that do accelerate**, two at z = +8.4 and +11.8, and the blanket rule bound on two of them — forcing a decelerating shape onto curves whose own z-score says they rise. `fit_burst_bounded(constrain="auto")` now asks each curve: the branch stays open where `acceleration` clears 3σ and is shut elsewhere. That bounds **24 of 27**, admits no negative v₀, and imposes no shape a curve's own statistic contradicts.

`summary_kinetics.fit_burst_bounded` implements the fit and profiles **both** v₀
and τ. The two are different questions and the second is the harsher: across
these 27 curves **v₀ is bounded on 19, τ is resolved on only 8, and both hold on 5**.
Where τ runs to an end of its grid the model has degenerated — to a step at the
floor, to a straight line at the cap — so v₀ → v_ss becomes exactly determined
while the *burst* means nothing. That is why exps 69 s3 and 70 s4 can report a
bounded v₀ alongside a negative `v_ss`. The panels say `τ unresolved, shape not
determined` whenever that holds. It is drawn on every curve panel as
a diagnostic — the shaded fan is the range of initial slopes the data still
allows — but it carries no reported number.

## 4. Separating substrate from buffer

Because no run moves `[buf]` at fixed `[sub]`, the buffer order is recovered from
the **disagreement between the two designs**. Read within runs, so that pH,
`[H2O2]`, cell and day are absorbed by a per-experiment offset on both sides:

| where `[buf]` is | order in `[sub]`, `v0_quad` | order in `[sub]`, `vmax` |
|---|---|---|
| held constant (67, 69, 70) | **+0.321 ± 0.056** (n = 12) | +0.343 ± 0.071 (n = 12) |
| falling as `[sub]` rises (3, 6) | **−0.306 ± 0.111** (n = 10) | −0.210 ± 0.092 (n = 8) |

**The apparent substrate order changes sign with the buffer design.** If the rate
goes as [S]<sup>a</sup>[buf]<sup>d</sup> and within the titrations
log[buf] = g·log[sub] + constant, then fitting the titrations without a buffer
term returns a′ = a + d·g, so

    d = (a' - a) / g,     g = -0.487

giving **d = +1.29 ± 0.25** on `v0_quad`. Both inputs are within-run contrast, so
nothing here asks one regression to separate two collinear predictors, and
nothing lets `[buf]` stand in for pH.

### Two routes that look better and are not

- **Fitting `[sub]` and `[buf]` jointly on exps 3 and 6.** VIF 9.7 and 8.3 on ten live curves; returns −0.09 ± 0.65 The design cannot carry it.
- **Fitting `[buf]` as a fourth pooled term across all six runs.** Returns a tight +0.95 ± 0.44, but `[buf]` is 85+ mM in every pH 8.0–8.5 run and sweeps only in the pH 6.71 ones, so it is partly a label for pH. The tell is that the [HOO⁻] order falls +0.836 → +0.739 when the buffer term is added — the buffer term stealing the pH effect.

## 5. Result: first order in buffer

| rate estimator | order in `[buf]`, BnOH 25 °C | cross-check, 4OMe-BnOH 40 °C |
|---|---|---|
| **`v0_quad`** | **+1.29 ± 0.25** | **+0.91 ± 0.38** |
| `v0_burst` | +1.20 ± 0.28 | +0.69 ± 0.24 |
| `vmax` | +1.31 ± 0.23 | +0.84 ± 0.29 |
| `v0` | +1.33 ± 0.32 | +0.86 ± 0.28 |
| `v0_whole` | +1.60 ± 0.23 | +1.38 ± 0.60 |

The anchor is exps 67, 69 and 70 — **phosphate only**. Exp 65 was in it until
2026-09-01 and is not any more, because its curves have no usable rate at all
(§6a). Every number in this table moved by at most 0.12 when it left.

The cross-check is independent in substrate, temperature and design geometry: its
buffer contrast lies *between* experiments at fixed pH (6.97–7.00) and fixed
`[H2O2]` (82.5 mM), so it takes no offsets and has no pH for `[buf]` to proxy.

**First order in buffer survives every reasonable estimator.** That is what
general acid/base catalysis of the H<sub>2</sub>O<sub>2</sub>/carbonyl addition
predicts (Sander & Jencks; `MECHANISM.md`).

**All five estimators agree**, including the two that choose no window at all.
They span **+1.20 to +1.60**, a spread narrower than the smallest of them, and
the weakest sits 4.3σ from zero.

> **Corrected 2026-09-01.** An earlier version of this document reported
> `v0_whole` at −3.90 ± 2.25 and built an argument on it: that a straight line
> through the whole curve is biased by the curves' deceleration and so is the
> one estimator that fails. **That was a bug, not a finding.**
> `curve_metrics.whole_slope` returned `line_fit(...)[:2]`, and `line_fit`
> returns `(intercept, slope, stderr, rms)` — so it handed back the *intercept*,
> an absorbance, as if it were a rate. Fixed; `v0_whole` agrees with the rest.
> The deceleration in section 7 is measured by `curvature_t` and is unaffected.

### A caveat on "initial rate", and what it costs

Four of the 27 curves are still **accelerating** where the initial rate is
read, and all four are in the titrations (exp 3 samples 2 and 3; exp 6 samples
1 and 2). On such a curve *every* initial-rate estimator — the quadratic, the
20% window, the burst form's v₀ — measures the **induction** rate rather than
the reaction at the stated concentrations. It is the same distinction
`curve_metrics.peak_rate` draws between `v0` and `vmax`, and it is why the
curve panels now name the fitted shape and both endpoints (`v₀ → v_ss`) instead
of printing one number called v₀.

They are 4 of the 10 live curves the titration arm rests on, so this is worth
measuring rather than asserting. `buffer_dependence(drop_accelerating=True)`:

| rate estimator | all live curves | accelerating dropped (n = 10 → 6) |
|---|---|---|
| `v0_quad` | +1.29 ± 0.25 | +1.36 ± 0.33 |
| `v0_burst` | +1.20 ± 0.28 | +1.38 ± 0.31 |
| `vmax` | +1.31 ± 0.23 | +1.34 ± 0.25 |
| `v0` | +1.33 ± 0.32 | +1.41 ± 0.35 |
| `v0_whole` | +1.60 ± 0.23 | +1.64 ± 0.27 |

Every estimator moves **up** slightly and every shift is inside its own
standard error, so the first-order reading does not depend on them. The
headline keeps all ten, because dropping a curve for the shape of its transient
is a stronger claim than the evidence needs.

## 6. The rate law

**Phosphate buffer only** — exps 3, 6, 67, 69, 70 (`scope.FREE_BNOH_PHOSPHATE`,
now the default everywhere in `scope` that fits a rate):

    v_background  ~  [S]^+0.32  [H2O2]^+1.49  [HOO-]^+0.82  [buf]^+1.29     (v0_quad)
    v_background  ~  [S]^+0.34  [H2O2]^+1.14  [HOO-]^+0.84  [buf]^+1.31     (vmax)

Sub-linear in substrate, first order in buffer, and about 0.8 order in [HOO⁻] —
so the background climbs nearly tenfold per pH unit, and a background measured
at pH 8 says little about one at pH 7.

**There is no "all six runs" version any more.** This document reported one
until 2026-09-01, labelled as the wider-in-pH alternative. It is withdrawn, not
demoted: exp 65 has no usable rate (§6a), so a law fitted with it in is not a
second opinion, it is the same law plus four numbers that do not mean anything.
The two estimators differ by 0.02 in [S], 0.35 in [H₂O₂] and 0.02 in [HOO⁻];
`v0_quad`'s pooled R² is 0.968 against `vmax`'s 0.966.

**On "first order in peroxide".** That names the `[H2O2]` coefficient with
[HOO⁻] carried separately. Since [HOO⁻] = f(pH)·[H₂O₂], the dependence on TOTAL
peroxide at fixed pH is the sum of the two, ≈ **+2.3** — second order, not
first. Whether that is two peroxide-derived species in one rate-determining
step or an artefact of splitting one variable in two is not settled here; the
two terms are separable in this design (VIF 1.3 and 3.1) but the split is a
modelling choice, not a measurement.

### 6b. Why first order in buffer? Probably catalysis; the boric evidence is real but confounded

The obvious reading of a first-order buffer dependence is that the buffer anion
is *making an oxidant* — phosphate + H₂O₂ → **peroxomonophosphate**, the way
borate gives peroxoborate and carbonate gives peroxymonocarbonate. If that were
happening, the buffer would not be catalysing the reaction; it would be a
reagent in it, and "background" would mean something different.

**Two mechanisms predict exactly this dependence:**

1. **General acid/base catalysis.** Sander & Jencks (MECHANISM.md item 45)
   showed H₂O₂ addition to carbonyls is subject to *both* general acid and
   general base catalysis. Buffer species accelerate adduct formation
   independently of pH, first order in the catalysing species. Already the
   repo's standing explanation.
2. **A buffer-derived peroxo oxidant.** A pre-equilibrium
   HPO₄²⁻ + H₂O₂ ⇌ HPO₅²⁻ + H₂O gives a rate first order in phosphate and
   first order in H₂O₂ while the equilibrium is far to the left.

**On the chemistry, (2) is much weaker than for borate or carbonate.**
Peroxomonophosphoric acid is a real and competent oxidant, but it is not made
by mixing inorganic phosphate with H₂O₂ in water — it is prepared from P₄O₁₀ or
from concentrated H₂O₂ with strong acid, and in dilute near-neutral solution it
*hydrolyses* to phosphate + H₂O₂. That is the thermodynamically downhill
direction. Borate is different because B(OH)₃ is a Lewis acid that adds HOO⁻
directly; carbonate is different because its carbon is electrophilic
(K ≈ 0.3 M⁻¹, formed in minutes). Phosphate at pH 7–8 is H₂PO₄⁻/HPO₄²⁻ — an
anion, so attack by HOO⁻ at a tetrahedral phosphorus is both electrostatically
disfavoured and slow. Expect a very small equilibrium constant, not a
kinetically useful one.

That is an argument, not a measurement, and the possibility is not zero: a tiny
equilibrium concentration of a *much* faster oxidant can still carry a rate.

**The phosphate runs cannot separate them.** Both mechanisms are first order
in a buffer *species*, and the design cannot resolve a species from the total.
Within the titration runs (exps 3 and 6, the only sweep of `[buf]` at fixed pH)
the pH is constant, so log[buf], log[H₂PO₄⁻] and log[HPO₄²⁻] are the same
variable — correlation **1.000000**, identical ranges. The measured +1.29 is
simultaneously an order in the total, in the acid form and in the base form.

Across pH there are only two phosphate levels (6.71 and 8.01), and everything
moves at once — medians over the live curves: `[buf]` 58 -> 85 mM, `[H2PO4-]`
43.4 -> 11.4, `[HPO4^2-]` 14.1 -> 73.6, `[HOO-]` 0.0013 -> 0.041 mM. (The
pH 6.71 runs are the titration, so their `[buf]` sweeps 85 → 25 mM; 58 is its
median, not a setpoint.) Substituting `buf_base` for `buf` in the rate
law drives the variance inflation factor to **30.2**, against **2.8** for the
total — above the package's own threshold of 10, where "the coefficient is
arithmetic, not evidence".

**The experiments that would settle it**, cheapest first:

- **³¹P NMR of the buffer under run conditions.** Peroxomonophosphate is a
  distinct resonance. Present or absent, one spectrum, no kinetics required.
  This is the decisive test and it should be done before the interpretation is
  written up either way.
- **Saturation in [H₂O₂].** A pre-equilibrium peroxo adduct saturates when
  buffer is limiting; general catalysis does not. Our peroxide order is ≥ 1
  everywhere measured, which is *not* what saturation looks like — weak
  evidence against (2), and confounded.
- **A Brønsted plot**: several buffers at one pH. General catalysis predicts
  the catalytic constant tracks the buffer's pKa; a peroxo route does not.
- **`[buf]` swept at a third pH**, which would break the species/total
  degeneracy that stops this dataset from answering at all.

`scope.frame` now carries `buf_acid`, `buf_base` and `buf_pka`
(`solution_chemistry.dominant_buffer_pair`) so all of this is a query rather
than an argument.

**But the boric run can, and it answers (2) no.** Added 2026-09-01,
`scope.peroxo_buffer_test`, `python data/scope.py --buffer`.

Borate is the buffer where the peroxo route is not a hypothesis: MECHANISM.md
item 39 has B(OH)₃ + H₂O₂ → peroxoborate with K = 2.0e-8, significant above
pH ≈ 7.7, and the anionic peroxoborates are *much faster oxidants than H₂O₂
itself*. Exp 65 is boric buffer at **pH 8.51** — 0.8 units above that
threshold, at 122 mM H₂O₂. If a buffer-derived peroxo oxidant is what carries a
first-order buffer term, exp 65 has to run far above a rate law fitted without
one.

The archive happens to hold the controlled comparison. **Exps 65 and 67 are
matched cuvette for cuvette**: the same substrate ladder (7.310 / 3.655 /
1.827 / 0.365 mM), the same [H₂O₂] (122.426 mM), the same temperature, the same
instrument, the same `.rre` source, 87.5 against 85.0 mM buffer. Only the salt
and the pH differ. Because [S] and [H₂O₂] match exactly, the predicted ratio
depends only on the `[buf]` and `[HOO⁻]` orders — the substrate and peroxide
orders, the two worst determined here, drop out — and the law is fitted on
`FREE_BNOH_PHOSPHATE`, so exp 65 is out of sample.

| estimator | predicted | observed | excess |
|---|---|---|---|
| `v0` | 1.89× | 1.25× | **0.66×** |
| `vmax` | 1.85× | 2.50× | **1.35×** |
| `v0_whole` | 1.73× | 1.44× | **0.83×** |
| `v0_quad` | 1.83× | 0.36× | **0.20×** (2 of 4 cuvettes) |

**No excess.** Three of four estimators put boric *below* the phosphate law and
one puts it 1.35× above. And the first-order buffer term is itself a phosphate
result — dropping exp 65 from `buffer_dependence`'s anchor as well leaves
**+1.31 ± 0.23** (`vmax`), **+1.29 ± 0.25** (`v0_quad`).

**But this probe is weaker than it first looked, and the reason is exp 65's
curve shape.** Raised from the plots on 2026-09-01: all four of its cuvettes
run bumpy to ~450 s, break sharply upward, and plateau. Measured
(`scope.synchronised_break`, `curve_metrics.segmented_fit`), the four break at
504, 532, 560 and 560 s — a span of **56 s**, two sampling intervals — and
across the break each one *steepens*, by **1.82, 2.04, 5.59 and 15.94×**.

| run | buffer | pH | enzyme | break span | slope after ÷ before | steepening |
|---|---|---|---|---|---|---|
| **65** | **boric** | **8.51** | **none** | **56 s** | **1.82–15.94** | **4 of 4** |
| 59–62 | boric | 8.51–9.00 | yes | 112–448 s | 0.65–0.95 | 0 of 16 |
| 67, 69, 70 | phosphate | 8.01 | none | 310–434 s | 0.22–1.13 | 0 of 12 |

Every other run decelerates after its break, as a reaction followed to a few
percent conversion should. Exp 65 is the only one that accelerates, and its
cuvette 4 ranks 6th of 402 archive-wide. **A single rate number therefore does
not describe exp 65**: `vmax` reads the post-break stretch and `v0` the
pre-break one. The excesses in the table above are a comparison between two
different processes.

They still say borate is **nowhere fast**, which is the direction that matters.
They do not say borate matches the law.

**The control set first used here was wrong, and the correction matters.** The
first version of this section argued the break was *not* borate chemistry
because exps 59–62 — boric, same and higher pH — are smooth. That does not
follow. **Every run is double-beam, and what the reference channel omits
decides what the curve means:** an enzyme run's reference omits the *enzyme*,
so the background sits in both beams and **cancels**, and the curve is a
catalytic increment. Exps 59–62 and 66 are all `with_E`. They cannot show a
background shape at all, so their smoothness is not evidence about one.

Restricted to the runs where a background feature is even visible — the 20
background experiments, 18 of which have live curves (`frame`'s
`differential == False`):

| buffer | experiments | any run whose cuvettes all steepen |
|---|---|---|
| phosphate | 17 | **0** (ratios 0.22–1.23; one isolated cuvette of six at 2.44) |
| **boric** | **1 (exp 65)** | **1 — all four cuvettes** |

So the honest count is **1 boric background run of 1 showing the break, against
17 of 17 not**. That is consistent with borate chemistry and it rests on a
single run. The other boric background run, **exp 64, was aborted at 448 s —
before exp 65's break at 504–560 s** — and is dead besides, so it cannot test
it either.

**What survives from the earlier version:** not the substrate (20-fold range,
break time unmoved, steepening *largest at the lowest* substrate); not
conversion (0.3–3.0% at the plateau). And the first bump *is* identified as
instrument settling — exp 64, flat and dead in the same session, shows the
identical 0–112 s transient with no chemistry at all. The break and the
plateau remain unexplained, and borate chemistry is now among the live
explanations rather than excluded.

### The second probe was withdrawn: it is blind, not merely confounded

Exps **66 and 68** were added here as a probe that avoids exp 65 — the same
design with enzyme, boric vs phosphate, both 0.028 mM enzyme, 85.0 mM buffer,
122.426 mM H₂O₂, 2.741 mM BnOH, 25 °C, both running smooth. **It cannot do the
job.** Both runs are catalysed, so both references omit the enzyme, so a
buffer-made oxidant — which acts on the *background* — is present in both beams
of each run and cancels. The comparison is structurally incapable of detecting
the effect it was built for.

It is kept as what it actually is, a catalysed boric-vs-phosphate comparison.
Uncorrected, cuvette for cuvette (`scope.CATALYSED_PEROXO_PAIR`):

| estimator | boric ÷ phosphate |
|---|---|
| `v0` | 0.66× |
| `vmax` | 0.65× |
| `v0_whole` | 0.54× |
| `v0_quad` | 1.11× |

Boric carries **2.19× the [HOO⁻]** and lands at or below phosphate — the
*catalysed* reaction is 2–4× slower in borate than pH alone would give it,
which is what §6a's peroxyacid-hydrolysis argument predicts (~12-fold, maximum
at pH 8.4–9, exactly where exp 66 sits). It says nothing about the background.

**Net position, and it is weaker than this section first claimed.** The peroxo
question is **open**, and the one piece of direct evidence — exp 65 — now
points mildly *toward* a buffer-made oxidant rather than against it, because
the only boric background run in the archive is the only background run in the
archive that breaks. Against that: the break's cause is genuinely unknown, n is
1, and exp 65 carries six other complaints.

The buffer term is still most likely **catalysis** on the strength of the
chemistry (§6b above: phosphate is an anion, HOO⁻ attack at tetrahedral
phosphorus is disfavoured, and peroxomonophosphate hydrolyses in the downhill
direction). But that is an argument, not a measurement, and the kinetics have
not closed it.

**Two experiments would.** ³¹P NMR of the buffer under run conditions settles
peroxomonophosphate outright, one spectrum, no kinetics. And **a repeat of exp
65 — enzyme-free, boric, run long — settles whether the break is real**; exp 64
was that experiment and it died at 448 s.

**What neither probe settles.** pH is matched in neither pair, so both lean on
a `[HOO⁻]` correction across 8.01 → 8.51. Each is one run of four cuvettes
against one run of four, on days and cells that are not controlled. And
strictly they show the peroxo mechanism failing to appear *where it certainly
operates*, which is evidence against it for phosphate but not a measurement on
phosphate.

Exp 65 now carries six independent complaints — the break, neither rate form
fitting it, four degenerate-yet-`bounded` burst fits, a negative `v0_quad` on
two cuvettes, the settling that outlasts the discarded reading, and a session
in which exp 63 has no data file and exp 64 is dead. It remains a named scope
rather than an exclusion; whether that should change is a ruling for the
author, recorded as open in `DATA_VERIFICATION.md`.

**So what IS the buffer catalysing? Not settled, and the standing citation has
a gap.** Sander & Jencks is the right *physics* — a buffer species in the
transition state of the rate-determining step, which is what a first-order
buffer term means and all it means — but their catalysis is of H₂O₂ addition to
a **carbonyl**, and the enzyme-free background has no carbonyl at t = 0. Two
ways out:

- **Benzaldehyde supplies it as it accumulates.** That route is autocatalytic,
  and this background does not accelerate — 1 of 12 at the fixed-buffer anchor
  against 50 of 110 catalysed increments, with exps 67, 69 and 70 actively
  decelerating (§7). Disfavoured by the data.
- **The buffer catalyses the oxidant's formation, before the substrate is
  involved.** The rate law has that shape: total peroxide ≈ **+2.3** at fixed
  pH against a substrate order of **+0.32**. Two peroxide-derived species in
  the rate-determining step, the substrate barely in it — making the oxidant
  is slow and oxidising BnOH is fast. A general base assisting the proton
  transfer that makes it is first order in buffer, which is what is measured.

The second is a hypothesis consistent with the orders, not a result, and it is
recorded as one in `MECHANISM.md`. What **is** established: the order is
**+1.19 ± 0.21** (`vmax`) and **+1.30 ± 0.25** (`v0_quad`) — 5σ from zero,
under 1.2σ from exactly one — independently **+0.83 ± 0.29** on the
4OMe-BnOH / 40 °C block; it survives dropping boric from both arms; and it is
catalysis rather than oxidant-making.

**A consequence worth stating plainly.** An order of 1 with no sign of
flattening means the buffer-free term is small over 25–85 mM: the background
is not a property of BnOH and H₂O₂, it is roughly *proportional to how much
buffer was in the cuvette*. Two "uncatalysed rates" measured in different
buffer concentrations are not comparable — and adding a buffer term to
`scope.literature_comparison` does move the median excess over the literature's
uncatalysed rate, from **34× to 25×** (§9).

### 6a. Exp 65 has no usable rate — ruled 2026-09-01

Exp 65 is the only boric-buffer run in the set. `MECHANISM.md` already says to
treat boric points as suspect and gives three reasons, every one of which bites
hardest exactly where exp 65 sits at pH 8.51: borate forms **peroxoborate**
with H₂O₂ (significant above pH ≈ 7.7), a much faster oxidant than H₂O₂ itself;
it generates **dioxaborirane**, a competing electrophilic oxidant; and boric
acid **catalyses peroxyacid hydrolysis** ~12-fold with a maximum at **pH 8.4–9**
— destroying the intermediate the mechanism runs through, right where this run
was measured.

The data said the same thing before that was consulted. Exp 65 is the one run
that neither rate form fits (§3a), and its samples 3 and 4 return a negative
`v0_quad`, so they silently leave any log-log fit — `v0_quad`'s "all runs"
numbers above rest on two of its four cuvettes.

**The ruling, and why it is not a sensitivity.** Until 2026-09-01 this section
reported the fit both ways and called excluding exp 65 a defensible choice.
That was too weak. The run's curves cannot be reduced to a rate at all:

- All four cuvettes **break upward mid-run**, at 504–560 s, steepening by
  1.82–15.94× across the break. So `v0` measures the pre-break stretch and
  `vmax` the post-break one; they are not two estimates of one quantity that
  disagree, they are **two different quantities**.
- `v0_quad` returns a **negative** rate on two of the four, which silently drops
  them from any log-log fit — so the withdrawn "all six runs" law rested on two
  of exp 65's four cuvettes without saying so.
- **Neither rate form fits it**: 6.9–8.3× noise for the burst, 4.8–5.7× for the
  quadratic, against a median 1.08–1.18× over the rest of the block.
- All four burst fits collapse to B = 0 with τ at its grid floor while still
  reporting `bounded = True`.

There is no defensible choice of estimator, so there is no defensible number,
and reporting the fit both ways implies one of them is right. Exp 65 is
therefore excluded from **every** default in `scope` that fits a rate
(`BORIC_RATE_UNUSABLE`): `background_orders`, `buffer_dependence`'s anchor,
`literature_comparison`, `background_model`, `predicted_enhancement`.

**It is deliberately not a `KNOWN_EXCLUSIONS` entry.** That would drop exp 65
from the dataset, and its *shape* is the most informative thing in the boric
block — it is the only background run in the archive that breaks, which is the
only direct evidence bearing on §6b's question. **The shape is evidence; the
rate is not.** `FREE_BNOH_ALL` still exists and still contains exp 65, for shape
work and for the record below.

**What the exclusion bought**, kept as the record of the decision rather than a
live alternative (`scope.FREE_BNOH_PHOSPHATE`, or `--scope free-bnoh-phosphate`):

    v_background  ~  [S]^+0.32  [H2O2]^+1.49  [HOO-]^+0.82  [buf]^+1.29     (v0_quad)
    v_background  ~  [S]^+0.34  [H2O2]^+1.14  [HOO-]^+0.84  [buf]^+1.31     (vmax)

**The estimators stop disagreeing.** The [HOO⁻] order was the worst of it,
0.63 against 0.84; without boric it is 0.82 against 0.84. `v0_quad`'s pooled R²
goes 0.903 → 0.968, overtaking `vmax`'s 0.966, so the sentence above about an
extrapolated rate being worse conditioned turns out to have been a statement
about exp 65 rather than about extrapolation.

Across all five estimators (`scope.boric_spread()`), the max-minus-min spread:

| order | with boric | phosphate only |
|---|---|---|
| `[S]` | 0.241 | **0.060** |
| `[H2O2]` | 0.719 | **0.349** |
| `[HOO-]` | 0.210 | **0.148** |
| `[buf]` | 0.329 | 0.393 |

Three of the four tighten, the substrate order fourfold. The buffer order is the
exception and does not need the help: it barely moves for any estimator
(+1.30 → +1.29 on the headline, largest move +0.12 on `vmax`, every one inside
its own standard error). **The headline result is indifferent to boric**, which
is why §5 is unchanged.

**What excluding it costs, and it is real.** Exp 65 is the only run at pH 8.51
and carried the top of the [HOO⁻] range by itself — 0.089 mM against 0.041 and
0.0012. Without it the pH axis has **two levels**, so the [HOO⁻] order is a
two-point contrast that cannot be checked for curvature, and at those two levels
pH, [buf] (58 vs 85 mM) and run design all move together.

That cost is accepted rather than traded off. A rate law fitted with exp 65 in
it does not *have* the pH coverage it appears to have: it has two good levels
plus a third whose rate has no defined value. **Missing coverage is preferable
to coverage that is wrong**, and the honest statement of the range is now in the
caveats. The experiment that restores it is a repeat of exp 65 — enzyme-free,
boric, run long — which §6b already wants for a different reason.

## 7. What the background is *not*

**It does not accelerate.** Source-matched — both sets entirely `.rre`, so this is
not the variance-floor artefact of `DATA_VERIFICATION.md` 2026-09-01:

| set | accelerating, > 3σ |
|---|---|
| enzyme-free BnOH, `[buf]` fixed | **1 of 12** |
| in-scope catalysed increments | **50 of 110 live** |

Exps 67, 69 and 70 actively *decelerate* (median `accel_z` −1.67, −5.15, −2.26).
**The autocatalytic signature in the in-scope block is not inherited from the
background.**

**Its deceleration is not substrate depletion.** Conversion runs 0.04–7.9%, yet 21
of 27 curves show curvature at |t| > 3, and exp 69 sample 3 roughly halves its
rate at 3.9% conversion. Something else decays during these runs — peroxide, or
the cell. This is an open question, not a settled finding.

## 8. What it means for exps 135–151

1. **The confound does not reach them.** `[buf]` = 75.013 mM in all 119 in-scope curves across all 17 runs — zero variation — so no buffer effect can enter their substrate order. The correction derived here is needed for exps 3 and 6, not for the scope.
2. **The background is already subtracted**, per cuvette, at matched conditions, because the catalysed reference omits only the enzyme.
3. **The amplitude is still missing.** There are 0 enzyme-free curves in the 127-curve pyrophosphate cell, so the absolute size of what was subtracted cannot be recovered without importing it from another buffer, which `MECHANISM.md` forbids. This remains `FITTING.md` F7: *the missing experiment is an enzyme-free control in pyrophosphate.*
4. **For scale**, on `vmax` within runs, the recorded increment is *flatter* in substrate than the background it was measured against: +0.091 ± 0.052 over all 110 live in−scope curves (+0.004 ± 0.046 on the 11 strong runs) against +0.343 ± 0.071 for the background. Relevant to `FITTING.md` F1, which is not restated here.

## 8a. Suspect readings are flagged, never removed

Each reading is scored by its **leave-one-out residual against a local
quadratic through its 8 neighbours**, in units of the curve's own noise
(`curve_metrics.local_outlier_z`). Readings beyond 5σ are ringed in red on
`progress_curves.html`. **Nothing is excluded** — every fit on every panel is
computed on every point, and the rings say which readings to discount by eye.

Two rules make that safe:

- **The local fit is quadratic, not linear.** A line is wrong wherever the
  curve bends, so it reads real kinetics as outliers: on exp 65 sample 3 a
  local line flags the genuine flat-to-rise transition at −6.5σ, where the
  quadratic scores it −2.7σ and leaves it alone.
- **Only *isolated* flags are ringed.** At 30–60 s sampling nothing chemical
  moves in one interval and reverts in the next, so a single reading out of
  line with both neighbours is an artefact; two or more consecutive ones are
  not separated from chemistry at all. Across the archive as it is now read
  that splits **445 isolated against 1429 in runs**, the longest run being 16
  readings. Runs are counted in the panel footer and never ringed —
  `curve_screen.py` is explicit that curve shape is never a defect.

The leading reading is the exception and gets its own flag. Note what that flag
now means. It was added when the instrument's first reading was the worst-
behaved point in the archive — **15.9% beyond 5σ against 7.5% of last
readings** on the identical extrapolation test — but that reading has since
been discarded from every run (`fit_dataset.DROP_FIRST_READING`), and on the
curves the flag now sees the excess is largely gone: **14.7% of leading
readings are flagged against 16.2% of last ones**, so point 0 is no longer a
distinct outlier class. The flag is kept for two reasons that do survive: v₀ is
an extrapolation *to* t = 0, so that point carries leverage no interior reading
has; and a bad leading point drags its neighbour past the threshold, so the
pair reads as a run and hides from `isolated` — on 21 of the 86 flagged curves
before the drop and 8 of the 56 after it.

**What the surviving rings are for.** Across the 27 panels the drop took the
rings from 16 to 11 and the point-0 rings from 14 to 7, and the ones it cleared
are exactly the ones it should have: exps 67 (samples 1–3), 69 (1, 4), 70
(1, 2) and 3 (sample 2) were ringed on their first reading alone and are now
clean. What is left is not the generic settling artefact — that went with the
discarded reading — but the runs whose settling lasted *longer than one
reading*: **all four cuvettes of exp 65** (leading z = +31.0, +19.6, +17.9,
+6.2, with a third reading now flagged on samples 1 and 2), **exp 70 sample 3**
(+6.5, which its own bad first reading had been masking), **exp 3 sample 6**
(+5.3) and **exp 6 sample 4** (−8.5). That is a short, named list of runs to
distrust at the start, which is a far more useful thing than a flag that fired
on one curve in five.

This is what exp 70 sample 2 turned out to be. Its first reading sits 11.5σ
below the trend while every later step is smooth; removing it takes the burst
v₀ profile from a half-width of 0.327 (unbounded) to 0.068 (bounded) and the
fit from 1.51× its own noise to 0.93×. It remains in every fit, and is ringed.

Known limitations, none fatal for an advisory flag: adjacent spikes mask each
other; a bad first reading can drag its neighbour so the pair reads as a run;
and a genuinely instantaneous kink would be flagged. All three are pinned by
`test_curve_metrics.test_outlier_flagging`.

In the enzyme-free BnOH set this ringed 10 isolated readings across 27
curves. In-scope it is 222 over 119 curves, 23 of them a first reading.

## 9. Caveats

- The quadratic does not fit within noise where deceleration is strongest (rms/noise 2–7 on exps 65, 67, 69), so its v₀ is an extrapolation from an approximate form. This is why three estimators are reported rather than one.
- The spread in the buffer order across anchors (+0.65 to +1.40 on `vmax`) is driven by which runs supply the substrate order; boric (exp 65) reads a much flatter substrate order than phosphate. One more reason to treat boric as suspect.
- Recovering the buffer order assumes the substrate order is the same at the titrations' pH 6.71 as at the anchor's pH 8.01–8.51. The archive has no run that could test this.
- Including a buffer term moves the literature comparison's median excess from **34× to 25×**. The default in `scope.literature_comparison` is deliberately left un-buffered so a figure quoted in four files does not move silently; the alternative is `orders_terms=("s0","h2o2","hoo","buf")`.

## 10. Reproducing

```bash
python data/scope.py --buffer                  # the separation result
python data/scope.py --scope free-bnoh-all     # the enzyme-free BnOH block
python background_reaction/build_figures.py    # both HTML pages
python background_reaction/check_numbers.py    # verify this document
```

Code added for this analysis: `curve_metrics.quadratic_rate` and
`curve_metrics.whole_slope`; `summary_kinetics.fit_burst_bounded`;
`scope.buffer_dependence`, `scope.buffer_cross_check`, `scope.blocks`,
`scope.variance_inflation`, and `terms`/`within`/`parameter` arguments on
`scope.background_orders`. `fit_dataset.Curve` gained `buf`.
