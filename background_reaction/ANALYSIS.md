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
| `scope.BUFFER_FIXED` | 65, 67, 69, 70 | 16 | held at 85–87.5 mM | ladder |
| `scope.BUFFER_CONFOUNDED` | 3, 6 | 11 (8 live) | falls 85 → 25 mM | ladder, rising |

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

`v0_quad` is the headline. It answers the objection that `INITIAL_WINDOW = 0.20`
is arbitrary — it chooses no window, uses every point, and being linear in its
parameters is always identified.

**The burst/lag form was tested and rejected as a rate source.** Fitting
A = c + v<sub>ss</sub>t − B(1 − e<sup>−t/τ</sup>) and taking v₀ = v<sub>ss</sub> − B/τ
is the natural whole-curve alternative, and it is what `summary_kinetics.fit_burst`
implements. On these curves it is not identified:

- Unconstrained, **5 of 27 curves return a negative v₀** where the line fit is firmly positive. Exp 67 sample 3 gives −2.07e-04 against a line's +3.35e-06.
- The existing `resolved` flag does not catch them: exp 3 sample 3 and exp 67 sample 3 are both `resolved=True`. `resolved` asks whether **τ** is located, which is a different question from whether **v₀** is.
- **A close fit does not license the extrapolation.** Exp 3 sample 3 is fitted to *better than its own noise* (rms/noise 0.73), yet across values of τ that are statistically indistinguishable (rms/noise 0.73 → 0.74) its v₀ ranges from **−1.62e-05 to +5.79e-07** — a factor of 28 and a change of sign.
- Constraining **B ≤ 0** — justified because 0 of 16 of these curves pass the acceleration test — removes every negative v₀ and raises the bounded count from **13 to 21 of 27**. It does not solve it: it trades the τ → ∞ degeneracy for τ → 0, and exps 65 samples 1 and 2 come back at **10× and 28×** their line rate.

`summary_kinetics.fit_burst_bounded` implements the constrained fit and profiles
v₀ directly (`v0_low`, `v0_high`, `bounded`). It is drawn on every curve panel as
a diagnostic — the shaded fan is the range of initial slopes the data still
allows — but it carries no reported number.

## 4. Separating substrate from buffer

Because no run moves `[buf]` at fixed `[sub]`, the buffer order is recovered from
the **disagreement between the two designs**. Read within runs, so that pH,
`[H2O2]`, cell and day are absorbed by a per-experiment offset on both sides:

| where `[buf]` is | order in `[sub]`, `v0_quad` | order in `[sub]`, `vmax` |
|---|---|---|
| held constant (65, 67, 69, 70) | **+0.330 ± 0.066** (n = 14) | +0.297 ± 0.080 (n = 16) |
| falling as `[sub]` rises (3, 6) | **−0.244 ± 0.084** (n = 8) | −0.210 ± 0.092 (n = 8) |

**The apparent substrate order changes sign with the buffer design.** If the rate
goes as [S]<sup>a</sup>[buf]<sup>d</sup> and within the titrations
log[buf] = g·log[sub] + constant, then fitting the titrations without a buffer
term returns a′ = a + d·g, so

    d = (a' - a) / g,     g = -0.432

giving **d = +1.33 ± 0.25** on `v0_quad`. Both inputs are within-run contrast, so
nothing here asks one regression to separate two collinear predictors, and
nothing lets `[buf]` stand in for pH.

### Two routes that look better and are not

- **Fitting `[sub]` and `[buf]` jointly on exps 3 and 6.** VIF 14.2 and 11.4 on eight live curves; returns +0.18 ± 0.66. The design cannot carry it.
- **Fitting `[buf]` as a fourth pooled term across all six runs.** Returns a tight +0.95 ± 0.44, but `[buf]` is 85+ mM in every pH 8.0–8.5 run and sweeps only in the pH 6.71 ones, so it is partly a label for pH. The tell is that the [HOO⁻] order falls +0.836 → +0.739 when the buffer term is added — the buffer term stealing the pH effect.

## 5. Result: first order in buffer

| rate estimator | order in `[buf]`, BnOH 25 °C | cross-check, 4OMe-BnOH 40 °C |
|---|---|---|
| **`v0_quad`** | **+1.33 ± 0.25** | **+0.87 ± 0.38** |
| `vmax` | +1.17 ± 0.28 | +0.83 ± 0.27 |
| `v0` | +1.67 ± 0.52 | +0.77 ± 0.26 |
| `v0_whole` | +1.62 ± 0.28 | +1.38 ± 0.60 |

The cross-check is independent in substrate, temperature and design geometry: its
buffer contrast lies *between* experiments at fixed pH (6.97–7.00) and fixed
`[H2O2]` (82.5 mM), so it takes no offsets and has no pH for `[buf]` to proxy.

**First order in buffer survives every reasonable estimator.** That is what
general acid/base catalysis of the H<sub>2</sub>O<sub>2</sub>/carbonyl addition
predicts (Sander & Jencks; `MECHANISM.md`).

**All four estimators agree**, including the two that choose no window at all.
The spread across them, +0.77 to +1.67, is smaller than the distance of any of
them from zero.

> **Corrected 2026-09-01.** An earlier version of this document reported
> `v0_whole` at −3.90 ± 2.25 and built an argument on it: that a straight line
> through the whole curve is biased by the curves' deceleration and so is the
> one estimator that fails. **That was a bug, not a finding.**
> `curve_metrics.whole_slope` returned `line_fit(...)[:2]`, and `line_fit`
> returns `(intercept, slope, stderr, rms)` — so it handed back the *intercept*,
> an absorbance, as if it were a rate. Fixed; `v0_whole` agrees with the rest.
> The deceleration in section 7 is measured by `curvature_t` and is unaffected.

## 6. The rate law

    v_background  ~  [S]^+0.33  [H2O2]^+1.57  [HOO-]^+0.67  [buf]^+1.33     (v0_quad)
    v_background  ~  [S]^+0.30  [H2O2]^+0.96  [HOO-]^+0.84  [buf]^+1.17     (vmax)

Roughly first order in peroxide and in buffer, sub-linear in substrate, and about
0.7–0.8 order in [HOO⁻] — so the background climbs nearly tenfold per pH unit, and
a background measured at pH 8 says little about one at pH 7. The peroxide and pH
orders are better conditioned on `vmax` (pooled R² 0.961 against 0.909): an
extrapolated initial rate adds variance where a block slope does not.

## 7. What the background is *not*

**It does not accelerate.** Source-matched — both sets entirely `.rre`, so this is
not the variance-floor artefact of `DATA_VERIFICATION.md` 2026-09-01:

| set | accelerating, > 3σ |
|---|---|
| enzyme-free BnOH, `[buf]` fixed | **0 of 16** |
| in-scope catalysed increments | **51 of 110 live** |

Exps 67, 69 and 70 actively *decelerate* (median `accel_z` −1.67, −5.15, −2.26).
**The autocatalytic signature in the in-scope block is not inherited from the
background.**

**Its deceleration is not substrate depletion.** Conversion runs 0.05–8%, yet 22
of 27 curves show curvature at |t| > 3, and exp 69 sample 3 roughly halves its
rate at 3.9% conversion. Something else decays during these runs — peroxide, or
the cell. This is an open question, not a settled finding.

## 8. What it means for exps 135–151

1. **The confound does not reach them.** `[buf]` = 75.013 mM in all 119 in-scope curves across all 17 runs — zero variation — so no buffer effect can enter their substrate order. The correction derived here is needed for exps 3 and 6, not for the scope.
2. **The background is already subtracted**, per cuvette, at matched conditions, because the catalysed reference omits only the enzyme.
3. **The amplitude is still missing.** There are 0 enzyme-free curves in the 127-curve pyrophosphate cell, so the absolute size of what was subtracted cannot be recovered without importing it from another buffer, which `MECHANISM.md` forbids. This remains `FITTING.md` F7: *the missing experiment is an enzyme-free control in pyrophosphate.*
4. **For scale**, on `vmax` within runs, the recorded increment is *flatter* in substrate than the background it was measured against: +0.097 ± 0.052 over all 110 live in-scope curves (+0.012 ± 0.044 on the 11 strong runs) against +0.297 ± 0.080 for the background. Relevant to `FITTING.md` F1, which is not restated here.

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
