# What pH does to the reaction

The archive holds four pH ladders and this project had used exactly one of
them for exactly one question. `PH_LADDER_PHOSPHATE`, `PH_LADDER_BORIC` and
the two-axis block's own `PH_LADDER_TWO_AXIS_LOW`/`_HIGH` were named in
`scope.py` on 2026-09-05 for the induction clock alone
(`induction.lag_ph_ladders`) and never asked about the **rate**. This folder
is that second half, read the same way, plus what the two questions say when
put next to each other.

    data/ph_role.py               the rate's order in [HOO-], per ladder and pooled
    python data/ph_role.py        the whole argument, printed
    python data/test_ph_role.py   planted recovery + the published numbers
    python ph/build_figures.py
    python ph/check_numbers.py

**Figures**: [`index.html`](index.html) is the presentation.
[`progress_curves.html`](progress_curves.html) carries every live cuvette of
all four ladders, which is the audit surface for every order below.

Related: [`../induction/`](../induction/ANALYSIS.md) owns the clock's own pH
order and the machinery this folder reuses rather than redefining;
[`../two_axis/`](../two_axis/ANALYSIS.md) owns the two-axis block's own
cuvette-matched pH reading, which section 2 checks this folder's cruder pooled
fit against; [`../buffer/`](../buffer/ANALYSIS.md) §5 owns why buffer identity
cannot be separated from pH anywhere in this archive, which is this folder's
main caveat rather than a settled question it can redo;
[`../early_trough/`](../early_trough/ANALYSIS.md) owns [HOO⁻] as the driver of
a different signal entirely, and section 6 is where the two meet.

## 1. Which curves, and a correction to how they had been collected by hand

Three folders under `data/Mads/` were hand-sorted for exactly this question
before this analysis existed: `VARIABLE pH (first run) phosphate`,
`done/Variable pH borich acid (4MeOBnOH)` and `bad data pH ca. 11`. Checked
against `scope.py`'s own selections:

- **The phosphate folder over-includes and under-includes.** It holds exps 2,
  4, 5, 7, 8, 9, 10, 11, 12, 20, 21, 22 — twelve runs.
  `scope.PH_LADDER_PHOSPHATE` holds nine: 8, 9, 10, 11, 12, **14**, 20, 21, 22.
  Exps 2, 4, 5 and 7 are `scope.REPLICATE_RUNS`, four repeats of ONE
  composition fixed at pH 6.71 — not four ladder rungs, and pooling them in
  would quadruple-weight that one pH value against every other point. Exp 14
  belongs and is missing: same phosphate/pH 7.00 design as the rest of the
  ladder, filed instead under `Variable Temperature/` because it also serves
  as the temperature series' 25 °C rung.
- **The boric folder over-includes by one.** It holds exp 13 alongside 41–49;
  `scope.PH_LADDER_BORIC` is 41–49 only. Exp 13 shares the cell and buffer
  type but carries a different substrate ladder and 65 mM buffer against the
  other nine's 70–85 mM — not matched to the rest.
- **The pH-11 folder is correctly excluded**, and independently so:
  `build_manifest.KNOWN_EXCLUSIONS[85]` rules the whole hand-sorted set
  unusable for kinetics on its own evidence — pH 11.84 in carbonate is 1.5
  units past carbonate's pKa2, and what the curves show runs backwards
  (anti-correlated with substrate), consistent with catalyst decay rather than
  turnover.
- **Neither folder touches the two-axis block's own pH ladders**
  (`PH_LADDER_TWO_AXIS_LOW`/`_HIGH`, exps 136–151), never hand-sorted because
  that design was found by structure rather than by folder. They are the
  best-powered pH result in the archive and this folder's own section 2 lives
  or dies by whether the other two ladders agree with them.

**This folder uses `scope.PH_LADDERS`, all four**, not any hand-sorted folder.

| ladder | substrate | buffer | runs | pH range | design |
|---|---|---|---|---|---|
| phosphate 4OMe | 4OMe-BnOH | phosphate | 9 | 5.64–8.95 | one substrate ladder per run |
| boric 4OMe | 4OMe-BnOH | boric | 9 | 8.46–10.34 | one substrate ladder per run |
| two-axis low | BnOH | pyrophosphate | 7 (5 strong) | 6.95–9.43 | seven cuvettes matched across runs |
| two-axis high | BnOH | pyrophosphate | 9 (5 strong) | 5.47–9.73 | seven cuvettes matched across runs |

## 2. The rate's order in [HOO⁻]: three ladders agree, one does not

`[buf]` sits exactly fixed inside every one of the four ladders (80.0, 85.0
and 75.013 mM, checked by hand — no within-ladder spread), and pH is one value
per run everywhere in this archive, so there is no per-run offset that can
carry a pH term without absorbing it. `ph_role.rate_ladder` fits
`log(vmax) = a·log[S] + b·log[HOO⁻] + intercept` pooled, with no offset —
every axis the design moves, in one regression, the way `scope.orders`
already argues for on the two-axis block's own L. On the two-axis ladders the
scope is filtered to `scope.strong_runs()` first, for the reason
`MECHANISM.md` gives: the weak runs sit at the cell's own drift and flatten
whatever they are pooled into.

| ladder | runs | pH range | order in [S] | order in [HOO⁻] |
|---|---|---|---|---|
| phosphate 4OMe | 9 | 5.64–8.95 | +0.369 ± 0.169 | **+0.596 ± 0.032** |
| two-axis low (strong) | 5 | 8.16–9.43 | −0.220 ± 0.072 | **+0.616 ± 0.066** |
| two-axis high (strong) | 5 | 7.60–9.73 | +0.042 ± 0.069 | **+0.570 ± 0.057** |
| boric 4OMe | 9 | 8.46–10.34 | +0.611 ± 0.097 | −0.034 ± 0.040 |

Pooled inverse-variance over all four: **+0.411 ± 0.022**, χ² = **174** on 3
degrees of freedom — nowhere close to consistent. Drop the boric row and the
other three collapse onto one number: **+0.594 ± 0.026**, χ² = **0.29** on 2 —
three ladders in two buffers, on two substrates, effectively reading the same
order.

**This checks against the two-axis block's own, better-powered reading.**
`scope.ph_order`, which matches cuvettes one-for-one across the seven-
composition ladder instead of pooling without an offset, already puts the
pooled two-axis order at **+0.554 ± 0.040** over the same strong runs
(`../two_axis/ANALYSIS.md`, `MECHANISM.md`). This folder's cruder fit —
+0.594 ± 0.026 pooled with phosphate — lands within one combined standard
error of it. The two independent methods agree; the outlier is the fourth
ladder, not the method.

**A half order in [HOO⁻] is what the archive's other reading of this axis
already found**, and this folder's phosphate and pyrophosphate numbers are a
second, independent confirmation of it from ladders `../two_axis/` does not
touch at all.

## 3. The boric ladder turns over, and it is not the O2 side reaction

A single log-log slope over a ladder that rises and then falls reports the two
halves cancelling, not the absence of a pH dependence — and that is what
section 2's boric row is. Read by experiment, median `vmax` **rises** from
7.1 × 10⁻⁵ AU/s at pH 8.46 to a peak of 1.42 × 10⁻⁴ at pH 9.50 (exp 43) and
**falls** to 8.1 × 10⁻⁵ and 6.6 × 10⁻⁵ at pH 10.07 and 10.34 — below every run
from pH 8.98 on.

Refitting the six runs at or below the peak (`ph_role.boric_turnover`) gives
**+0.222 ± 0.037** in [HOO⁻] — positive, six standard errors from zero, and
still markedly weaker than phosphate or pyrophosphate's +0.57 to +0.62. The
two runs above the peak are too few to fit a slope at all, so their decline is
reported as what it is rather than forced into one.

**This is not the O2 side reaction reading as a rate decline.** Boric buffer
is exactly where the archive's gas artefact is heaviest at high pH
(`scope.gas_survey`: 24 of 64 boric curves above pH 8.5 detach gas, 1.06
events/hour, against a hard floor of zero detachments anywhere in the archive
below pH 7.5) — so the possibility that
`debubble`'s correction would restore the missing rate had to be checked
rather than assumed away. It does not: `vmax_corrected` moves the pH 10.07 and
10.34 medians by less than a quarter of themselves (8.1 → 9.4 × 10⁻⁵ and
6.6 → 7.0 × 10⁻⁵), in the direction that makes the decline slightly smaller,
not the direction that would erase it.

**What is left is a genuine turnover in the catalysed rate above about
pH 9.5–10, in boric buffer.** The archive has one other observation in this
direction: exp 85, ruled unusable for kinetics at pH 11.84
(`build_manifest.KNOWN_EXCLUSIONS`), was flagged there as possible evidence of
"catalyst instability at pH 11.8." This ladder's own decline starts a full 1.8
pH units below that, in a different buffer, substrate and enzyme loading, and
its curves are otherwise perfectly well-behaved — every one of the nine boric
runs has corr([S], vmax) between 0.91 and 1.00 within itself
(checked by hand from `scope.frame`), so this is not a data-quality problem
the way exp 85's was. **Whether the same cause is responsible is not
established**; see section 7.

## 4. The induction clock's own pH order, quoted back

`induction.lag_ph_ladders` already reads all four ladders' clocks at a window
every run in the ladder shares, and `induction.pooled_ladder` combines them —
built 2026-09-05, before this folder existed, and not recomputed here.

| ladder | slope, per pH unit | schedule collinearity | signal collinearity |
|---|---|---|---|
| phosphate 4OMe | −0.253 ± 0.159 (controlled +0.093 ± 0.302) | −0.25 | +0.87 |
| boric 4OMe | −0.020 ± 0.286 (controlled +0.250 ± 0.316) | +0.71 | −0.61 |
| two-axis low | +0.250 ± 0.286 (controlled +0.427 ± 0.367) | −0.53 | +0.77 |
| two-axis high | +0.370 ± 0.144 (controlled +0.413 ± 0.184) | −0.79 | +0.85 |

Pooled with the signal held (`clock_pooled_order`, `controlled=True`):
**+0.326 ± 0.131 per pH unit**, χ² = **0.95** on 3 — the four ladders agree,
unlike the rate's. **More alkaline, longer induction**, on all four, in three
buffers and on both substrates.

## 5. The two confounds run in different directions on the two questions

Every pH ladder in this archive is also, unavoidably, a schedule and a
signal-to-noise ladder (`induction.lag_ph_ladders`'s own `schedule_
collinearity`/`signal_collinearity`, section 4's table; `ph_role.rate_ladder`
computes the schedule figure separately for the rate, below). Neither question
can be answered by one ladder alone for exactly this reason.

| ladder | pH vs. schedule (rate reading) | pH vs. schedule (clock reading) |
|---|---|---|
| phosphate 4OMe | +0.32 | −0.25 |
| boric 4OMe | +0.36 | +0.71 |
| two-axis low | +0.84 | −0.53 |
| two-axis high | −0.74 | −0.79 |

**On the rate**, phosphate and boric share a weak positive schedule confound
while the two two-axis ladders are strongly opposed to each other (+0.84,
−0.74) — so a coefficient positive on both of the opposed pair cannot be a
monotone drift, and it is (+0.616, +0.570). **On the clock**, three of four
run negative and boric alone is positive; the pooled clock order is still
positive, in the direction none of the three negatively-confounded ladders'
own drift could produce by itself. The two questions are protected by
different ladders playing the control role, which is the reason to read all
four rather than the two `../two_axis/` already covers.

## 6. One species, three independent signals — with a caveat the archive cannot close

[HOO⁻] rather than pH itself, and rather than total [H2O2], is the quantity
that keeps coming up as the operative one:

- **The rate.** Section 2: +0.594 ± 0.026 pooled over phosphate and both
  two-axis ladders, checked against `scope.ph_order`'s own +0.554 ± 0.040.
- **The induction clock.** Section 4: +0.326 ± 0.131 per pH unit, pooled and
  signal-controlled, χ² = 0.95 on 3.
- **The early trough.** `../early_trough/ANALYSIS.md`'s driver of the
  reference-cuvette dip is [enz]/[HOO⁻], not [enz]/[S] and not total [H2O2],
  at p < 10⁻⁴ in both substrates independently (ρ = −0.621 for 4OMe, −0.315
  for BnOH) — a signal that has nothing to do with the catalysed rate or the
  induction clock, measured on the OTHER cuvette of the pair, and it names the
  same species.
- **The O2 side reaction.** `scope.gas_survey`: zero detachments below pH 7.5
  in 270 hours of catalysed running across every buffer, rising steeply above
  it, in three buffers and on both substrates. A catalyst engaging peroxide
  non-productively should have the same pH signature as one engaging it
  productively, and it does.

**The caveat.** `scope.hoo_consistency`, on the two-axis block alone, moves
[HOO⁻] two independent ways — within a run (via pH, holding [H2O2] fixed) and
across runs (via [H2O2], holding pH close to fixed) — and the two orders part
at **3.4σ** (+0.762 ± 0.046 within, +0.554 ± 0.040 across). [HOO⁻] is not the
whole story even where the rate's dependence on it is real and reproduced
three times over; either `[H2O2]` contributes something beyond its
hydroperoxide content, or the rate saturates in a way [HOO⁻] alone does not
capture. Section 3's boric turnover is a second, independent way the simple
picture breaks — a NEW way, since it happens at fixed [H2O2] and moves with pH
alone.

## 7. What this settles, and what it does not

**Settled.**

- The archive's hand-sorted pH folders needed two corrections (§1): drop the
  replicate quadruplet from the phosphate set, add exp 14 to it, drop exp 13
  from the boric set. The pH-11 folder was already correctly excluded.
- The catalysed rate's order in [HOO⁻] is **+0.594 ± 0.026** in phosphate and
  pyrophosphate (χ² = 0.29 on 2, three ladders, two buffers, two substrates),
  checked against and consistent with the two-axis block's own independent,
  cuvette-matched +0.554 ± 0.040.
- The boric ladder does **not** share that order: it turns over near pH 9.5,
  and `vmax_corrected` rules out the O2 side reaction as the cause of the
  decline above it.
- The induction clock's order in pH, +0.326 ± 0.131 per pH unit, agrees across
  all four ladders (χ² = 0.95 on 3) where the rate's does not.
- Three independent measurements — the rate, the clock, and the early
  trough's OTHER-cuvette signal — name [HOO⁻] as the species that matters, and
  a fourth (the O2 side reaction's pH onset) is consistent with the same
  reading.

**Not settled.**

- **Why boric turns over.** A buffer-specific effect on the catalyst (borate
  is a known Lewis-acidic complexing agent) and a genuine high-pH instability
  of the kind exp 85 already showed at pH 11.84 are both consistent with one
  ladder alone; the archive holds no other buffer above pH 9 to separate them,
  which is exactly `../buffer/ANALYSIS.md` §5's finding that buffer identity
  and pH cannot be told apart anywhere in this archive.
- **Whether [HOO⁻] is the whole story even where its order is real.**
  `hoo_consistency`'s 3.4σ gap on the two-axis block says no, on the block
  where it is best measured.
- **What in the catalyst engages HOO⁻.** `MECHANISM.md`'s own caveat on the
  gas applies here without modification: nothing in this archive has directly
  measured the species, and a cyclodextrin or a trace transition metal would
  give the same peroxide-anion order the ketone would.

**What would settle the open one.** A fifth pH ladder in a buffer that holds
pH 9–10.5 without boric's own chemistry — a carbonate or a phosphate/borate
mixed system buffered past phosphate's own working range — run at the same
substrate ladder and enzyme loading as this folder's other three, would tell
the turnover apart from the buffer. Short of a new experiment, the closest
existing evidence is `background_reaction/ANALYSIS.md` §6b's finding that
boric is separately unusable for the ENZYME-FREE rate (`scope.
BORIC_RATE_UNUSABLE`) for reasons that are about that buffer specifically, not
about pH — worth reading before assuming this ladder's turnover is chemistry
rather than the buffer.
