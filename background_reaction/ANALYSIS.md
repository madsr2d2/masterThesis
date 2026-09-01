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
| median residual | 1.18× noise | 1.08× noise |
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
| held constant (65, 67, 69, 70) | **+0.328 ± 0.054** (n = 14) | +0.282 ± 0.057 (n = 16) |
| falling as `[sub]` rises (3, 6) | **−0.306 ± 0.110** (n = 10) | −0.210 ± 0.092 (n = 8) |

**The apparent substrate order changes sign with the buffer design.** If the rate
goes as [S]<sup>a</sup>[buf]<sup>d</sup> and within the titrations
log[buf] = g·log[sub] + constant, then fitting the titrations without a buffer
term returns a′ = a + d·g, so

    d = (a' - a) / g,     g = -0.487

giving **d = +1.30 ± 0.25** on `v0_quad`. Both inputs are within-run contrast, so
nothing here asks one regression to separate two collinear predictors, and
nothing lets `[buf]` stand in for pH.

### Two routes that look better and are not

- **Fitting `[sub]` and `[buf]` jointly on exps 3 and 6.** VIF 9.7 and 8.3 on ten live curves; returns −0.09 ± 0.65 The design cannot carry it.
- **Fitting `[buf]` as a fourth pooled term across all six runs.** Returns a tight +0.95 ± 0.44, but `[buf]` is 85+ mM in every pH 8.0–8.5 run and sweeps only in the pH 6.71 ones, so it is partly a label for pH. The tell is that the [HOO⁻] order falls +0.836 → +0.739 when the buffer term is added — the buffer term stealing the pH effect.

## 5. Result: first order in buffer

| rate estimator | order in `[buf]`, BnOH 25 °C | cross-check, 4OMe-BnOH 40 °C |
|---|---|---|
| **`v0_quad`** | **+1.30 ± 0.25** | **+0.91 ± 0.38** |
| `v0_burst` | +1.19 ± 0.25 | +0.69 ± 0.24 |
| `vmax` | +1.19 ± 0.21 | +0.83 ± 0.29 |
| `v0` | +1.30 ± 0.30 | +0.85 ± 0.28 |
| `v0_whole` | +1.52 ± 0.22 | +1.38 ± 0.60 |

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
| `v0_quad` | +1.30 ± 0.25 | +1.37 ± 0.33 |
| `v0_burst` | +1.19 ± 0.25 | +1.37 ± 0.28 |
| `vmax` | +1.19 ± 0.21 | +1.22 ± 0.23 |
| `v0` | +1.30 ± 0.30 | +1.38 ± 0.34 |
| `v0_whole` | +1.52 ± 0.22 | +1.56 ± 0.27 |

Every estimator moves **up** slightly and every shift is inside its own
standard error, so the first-order reading does not depend on them. The
headline keeps all ten, because dropping a curve for the shape of its transient
is a stronger claim than the evidence needs.

## 6. The rate law

    v_background  ~  [S]^+0.33  [H2O2]^+1.87  [HOO-]^+0.63  [buf]^+1.30     (v0_quad)
    v_background  ~  [S]^+0.28  [H2O2]^+1.15  [HOO-]^+0.84  [buf]^+1.19     (vmax)

Roughly first order in peroxide and in buffer, sub-linear in substrate, and about
0.6–0.8 order in [HOO⁻] — so the background climbs nearly tenfold per pH unit, and
a background measured at pH 8 says little about one at pH 7.

The two estimators disagree on the peroxide and pH orders, and `vmax` is the
better conditioned of them (pooled R² 0.974 against 0.903). **That disagreement
is one experiment**, and §6a is about which.

### 6a. What the boric run carries

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

Excluding it (`scope.FREE_BNOH_PHOSPHATE`, or `--scope free-bnoh-phosphate`):

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
| `[H2O2]` | 0.718 | **0.348** |
| `[HOO-]` | 0.210 | **0.148** |
| `[buf]` | 0.329 | 0.392 |

Three of the four tighten, the substrate order fourfold. The buffer order is the
exception and does not need the help: it barely moves for any estimator
(+1.30 → +1.29 on the headline, largest move +0.12 on `vmax`, every one inside
its own standard error). **The headline result is indifferent to boric**, which
is why §5 is unchanged.

**What excluding it costs, and why this is a scope and not a deletion.** Exp 65
is the only run at pH 8.51 and carries the top of the [HOO⁻] range by itself —
0.089 mM against 0.041 and 0.0012. Without it the pH axis has **two levels**, so
the [HOO⁻] order becomes a two-point contrast that cannot be checked for
curvature, and at those two levels pH, [buf] (55 vs 85 mM) and run design all
move together. A tighter number is not automatically a better-determined one.
Both versions are reported for that reason.

## 7. What the background is *not*

**It does not accelerate.** Source-matched — both sets entirely `.rre`, so this is
not the variance-floor artefact of `DATA_VERIFICATION.md` 2026-09-01:

| set | accelerating, > 3σ |
|---|---|
| enzyme-free BnOH, `[buf]` fixed | **1 of 16** |
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
4. **For scale**, on `vmax` within runs, the recorded increment is *flatter* in substrate than the background it was measured against: +0.091 ± 0.052 over all 110 live in−scope curves (+0.004 ± 0.046 on the 11 strong runs) against +0.297 ± 0.080 for the background. Relevant to `FITTING.md` F1, which is not restated here.

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
  that splits **442 isolated against 1429 in runs**, the longest run being 16
  readings. Runs are counted in the panel footer and never ringed —
  `curve_screen.py` is explicit that curve shape is never a defect.

The leading reading is the exception and gets its own flag. Note what that flag
now means. It was added when the instrument's first reading was the worst-
behaved point in the archive — **15.9% beyond 5σ against 7.5% of last
readings** on the identical extrapolation test — but that reading has since
been discarded from every run (`fit_dataset.DROP_FIRST_READING`), and on the
curves the flag now sees the excess is largely gone: **13.9% of leading
readings are flagged against 15.2% of last ones**, so point 0 is no longer a
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
