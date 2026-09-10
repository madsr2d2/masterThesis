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

**Two more BnOH-in-boric sources exist, both too small for `PH_LADDERS` but
not too small to show.** `data/Mads/done/Boric acid buffer BnOH/` holds exps
60–62 at 0.014 mM chemzyme (`scope.PH_LADDER_BORIC_BNOH`, §3b) — a genuine
BnOH pH ladder in boric buffer that an earlier version of a comment in
`scope.py` wrongly asserted could not exist. `data/Mads/good data BnOH/`
holds exp 50 alongside the two-axis block's own exps 135–151: it is BnOH in
boric buffer too, at 0.28 mM chemzyme — the SAME loading `PH_LADDER_BORIC`
itself uses — but `build_manifest.KNOWN_EXCLUSIONS[50]` rules it a
reaction-direction failure (all four curves descend with no ordering by
substrate at all), and its same-composition, same-pH repeat, exp 55, is
clean. Exp 55 and its own composition partner exp 51 (pH 9.01) are
`scope.BORIC_BNOH_HIGH_ENZYME_PAIR` (§3c) — two pH points, not three, and
already named elsewhere in `scope.py` as half of `SUBSTRATE_PAIRS`, so
neither curve is new to the project, only to this folder's own pH-vs-rate
reading of them. Both sets are on the curves page in full.

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

**The `[S]` column is an order in the [S]/[buf] pair on the two 4OMe rows.**
`[buf]` is fixed between the runs of every ladder — 80.0, 85.0 and 75.013 mM —
which is what leaves the pH axis clean, since pH is a between-run axis here.
But it steps *within* the runs of the phosphate and boric ladders (50–80 and
70–85 mM, four values each), collinear with `log[S]` at −0.960 and −0.974,
because substrate volume displaced buffer volume. That is the archive-wide
pairing `induction.composition_collinearity` measures, and it is why this
folder concludes from the [HOO⁻] column and not the `[S]` one. Adding `buf` as
a third term moves `order_hoo` by 0.002 and 0.001 and returns an unresolved
buffer order either way.

Pooled inverse-variance over all four: +0.411 ± 0.022, χ² = **174** on 3
degrees of freedom — nowhere close to consistent, so that ± is not a real
uncertainty and the value is quoted only as the contrast it makes. Drop the
boric row and the other three collapse onto one number: **+0.594 ± 0.026**,
χ² = **0.29** on 2 — three ladders in two buffers, on two substrates,
effectively reading the same order.

**Two of those three ladders are not independent of `../two_axis/`, and the
one that is carries the corroboration.** `PH_LADDER_TWO_AXIS_LOW` and `_HIGH`
are exps 136–142 and 143–151 — the two-axis block itself, read along its
second design. They contribute **35.9% of the weight** in the pooled
+0.594 ± 0.026 (phosphate carries the other 64.1%), so that pooled number and
the block's own reading share curves and cannot check each other.

The comparison that *is* independent is the phosphate ladder alone:

| | order in [HOO⁻] | shares experiments with the block? |
|---|---|---|
| phosphate 4OMe, this folder | **+0.596 ± 0.032** | no — exps 8–22, 4OMe, phosphate |
| `scope.ph_order`, strong runs | **+0.554 ± 0.040** | it *is* the block |

The two differ by 0.042 against a combined standard error of 0.051 — **0.83σ**.
Nine phosphate runs on the other substrate, in another buffer, sharing not one
experiment with exps 135–151, land on the order the block measures by matching
cuvettes one-for-one. That is the second reading of this axis, and the two
pyrophosphate rows are better read as a consistency check that the pooled
method reproduces the cuvette-matched one on the same curves (it does) than as
new evidence.

**A half order in [HOO⁻] is therefore the archive's answer on two independent
substrate/buffer systems**, not one — with the outlier being the fourth ladder,
not the method.

## 3. The boric ladder turns over, and it is not the O2 side reaction

A single log-log slope over a ladder that rises and then falls reports the two
halves cancelling, not the absence of a pH dependence — and that is what
section 2's boric row is. Read by experiment, median `vmax` **rises** from
7.1 × 10⁻⁵ AU/s at pH 8.46 to a peak of 1.42 × 10⁻⁴ at pH 9.50 (exp 43) and
**falls** to 8.1 × 10⁻⁵ and 6.6 × 10⁻⁵ at pH 10.07 and 10.34 — below every run
from pH 8.98 on.

Refitting the six runs at or below the peak (`ph_role.boric_turnover`) gives
+0.222 ± 0.037 in [HOO⁻] — but **that ± is not the uncertainty that matters,
and the split it comes from was chosen by looking at the medians above.**
Where the ladder is cut is an analyst's choice, and the order below the cut
depends on it far more than on its own standard error
(`ph_role.boric_split_sensitivity`):

| split | runs below | order in [HOO⁻] | σ from zero |
|---|---|---|---|
| 9.01 | 2 | +0.481 ± 0.023 | 21.2 |
| 9.24 | 3 | +0.353 ± 0.044 | 8.1 |
| 9.41 | 4 | +0.259 ± 0.044 | 5.8 |
| 9.51 *(published)* | 6 | +0.222 ± 0.037 | 6.0 |
| 9.71 | 7 | +0.186 ± 0.035 | 5.3 |
| 10.08 | 8 | +0.056 ± 0.043 | 1.3 |

The estimate spans **+0.06 to +0.48** — a factor of nine — while its quoted
error stays near ±0.04. Quote it as **+0.22 with a split systematic of about
+0.26/−0.17**, never as +0.222 ± 0.037; this is the same discipline
`slowdown.sink_window_sensitivity` imposes on the sink's 72 kJ/mol and
`induction.lag_window_sweep` on the clock's pooled order.

**What survives the choice is the sign and the comparison**, and both are what
this section actually argues: the order below the peak is positive at every
split, and weaker than the +0.594 the other three ladders share at every split
— including the two splits chosen without reference to the peak at all
(9.24 is boric's own pKa). The two runs above the peak are too few to fit a
slope, so their decline is reported as medians rather than forced into one.

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

### 3a. The turnover is not only `Vmax` — and `Km` turns over too, not just up

Section 3's "median vmax" is read at whatever `[S]` each cuvette happened to
carry, not at saturation, so a decline in it is not automatically a decline in
`Vmax`. Fitting each boric run's own four-cuvette substrate ladder to
Michaelis-Menten (`ph_role.ladder_mm_table`, profiling `Km` the way
`summary_kinetics.profile_km` profiles the same parameter in its own log-log
form) resolves a per-experiment `Km` on **7 of the 9 runs** — exps 43
(pH 9.50) and 45 (pH 9.70) do not (`km_resolved` false, the profile reaching
both edges of its own grid) and are excluded from what follows.

**Read off `v_peak_corrected`, not `v_peak`.** The fitted rate and the fit
itself are read off the debubbled curve (`scope.py`'s `v_peak_corrected`
column, off the exact fit `tau_corrected` already uses), not the raw
readings: exps 44 and 49, the ladder's two highest-pH runs, are also two of
the archive's heaviest boric gassers (`scope.gas_survey`), and leaving their
own O2 in moves both `Vmax` and `Km` at exactly those rungs. Exp 43's own fit
improves once its own worst-gassing cuvette is corrected too — its R² turns
positive, **+0.65**, against a negative value uncorrected — though its `Km`
stays unresolved either way.

| experiment | pH | `Vmax` (AU/s) | `Km` (mM) |
|---|---|---|---|
| 41 | 8.46 | 1.28 × 10⁻⁴ | 2.96 |
| 42 | 8.98 | 3.34 × 10⁻⁴ | 5.56 |
| 46 | 9.23 | 4.45 × 10⁻⁴ | 8.13 |
| 47 | 9.40 | 3.16 × 10⁻⁴ | 5.18 |
| 48 | 9.51 | 4.74 × 10⁻⁴ | 10.65 |
| 44 | 10.07 | 2.28 × 10⁻⁴ | 5.98 |
| 49 | 10.34 | 2.06 × 10⁻⁴ | 6.43 |

**`Km` does not simply climb with pH — it turns over, at the same rung
`Vmax` does.** Both peak together at exp 48 (pH 9.51: `Km` 10.65 mM, `Vmax`
4.74 × 10⁻⁴ AU/s) and both fall over the ladder's top two rungs. This revises
the uncorrected reading: there, exps 44 and 49's own uncorrected gas inflated
their fitted `Km` past exp 48's (11.9 and 14.2 mM against 10.5), which read
as a `Km` still climbing to the ladder's own top rung. Debubbled, it does
not — the two rungs' own gas was the reason it looked like it did.

`ph_role.boric_vmax_km_decomposition` prices the two effects apart the same
way as before: holding ONE parameter at exp 41's own value (the lowest-pH
resolved rung) and letting the OTHER move as fitted, evaluated at each run's
own cuvette concentrations — every boric run shares the same substrate
ladder, so this is not comparing different `[S]` sets.

| series | what moves | peak → last rung |
|---|---|---|
| debubbled `v_peak`, as measured | both | **−47.4%**, peaking at exp 46 (pH 9.23) |
| `Vmax` alone | `Km` held at exp 41's value | **−56.6%**, peaking at exp 48 (pH 9.51) |
| `Km` alone | `Vmax` held at exp 41's value | **−35.3%**, peaking at exp 41, the ladder's own lowest pH |

**`Vmax` alone now falls FURTHER than the measured statistic does — the
opposite of the uncorrected reading.** `Km`'s own fall over the top two
rungs is why: once `Km` is past its own peak, a FALLING `Km` at fixed `Vmax`
*raises* the rate read at a fixed, sub-saturating `[S]`, so `Km`'s decline
there partly OFFSETS `Vmax`'s rather than compounding it — which is why the
measured decline (−47.4%) is milder than `Vmax` acting alone would produce
(−56.6%). `Km` alone's own counterfactual maximum still sits at exp 41, the
ladder's lowest pH, exactly as it did uncorrected: that counterfactual holds
`Vmax` fixed and lets `Km` climb, and a climbing `Km` can only ever shrink
the rate at fixed sub-saturating `[S]`, so its own peak has to sit where
`Km` is smallest regardless of what the real ladder does above the split.

**The ladder now supports one `Km` held fixed across it, where it did not
before.** `ph_role.ladder_mm_shared` fits one `Km` shared across all nine
runs (`Vmax` free per run) and now resolves it — **5.87 mM (3.24–12.31 mM)**
— against the grid-floor, unresolved fit the uncorrected readings gave. Six
of the seven resolved rungs' own 95% intervals contain it; only exp 41, the
ladder's lowest-pH rung and its most tightly resolved `Km` (2.66–3.36 mM),
sits outside it. That single disagreement — not a ladder-wide climb — is the
one `Km`-vs-pH signal this ladder still supports once the gas is removed.

### 3b. The other substrate's own boric ladder (exps 60-62) is too sparse to extend this

BnOH's own pH ladder in boric buffer — `scope.PH_LADDER_BORIC_BNOH`, exps 60
(pH 8.51), 61 (pH 8.75) and 62 (pH 9.00), one fixed `[enz]` (0.014 mM) and the
same four-point substrate ladder across all three runs — sat outside this
folder until now: an earlier version of a comment in `scope.py` wrongly
claimed the two-axis block's own pyrophosphate ladders were "the only pH
ladders in the archive on the other substrate." They were not, and this is
the ladder that comment missed. It carries no O2 detachments at all (zero
`bubble_events` on all twelve live cuvettes), consistent with
`scope.gas_survey` needing pH above about 8.5 for the side reaction and this
ladder barely reaching it — `v_peak_corrected` is identical to `v_peak` here.

Only exp 60 (pH 8.51) resolves its own `Km` (4.48 mM); exp 61's profile
reaches the grid's own ceiling, and so does exp 62's, whose top cuvette sits
an order of magnitude above its own neighbours' trend for a reason these
readings alone cannot diagnose. Three runs with one resolved `Km` is not
enough to ask whether `Vmax` or `Km` moves with pH here — a shared-`Km` fit
or a decomposition both need more than one resolved rung to mean anything, so
neither is attempted, and it is why this ladder is kept out of
`scope.PH_LADDERS`: three points cannot be pooled into the four-ladder order
without a size the fit cannot support. Its own progress curves and
per-experiment MM panels are on the curves page in full.

### 3c. BnOH in boric buffer at `PH_LADDER_BORIC`'s own enzyme loading (exps 51, 55) — two points, one of them very gassy

The other BnOH-in-boric enzyme tier, `scope.BORIC_BNOH_HIGH_ENZYME_PAIR`,
sits at ~0.28 mM chemzyme — essentially the SAME loading as
`PH_LADDER_BORIC` (0.270 mM), not a third one — at two pH points: exp 51
(pH 9.01) and exp 55 (pH 9.70, the clean repeat of the excluded exp 50; §1).
Each resolves its own per-experiment fit: exp 51's `Vmax` is
1.24 × 10⁻⁴ AU/s (`Km` unresolved, the profile reaching the grid's own
ceiling); exp 55's is 0.98 × 10⁻⁴ AU/s, with `Km` resolved at 15.8 mM
(7.3–72.1 mM). Exp 55's own lowest-`[S]` cuvette carries **25** O2
detachments — the heaviest single cuvette this whole folder's ladders carry
— and correcting it drops that cuvette's own `v_peak` almost sixfold
(7.5 × 10⁻⁵ to 1.2 × 10⁻⁵ AU/s); its other three cuvettes carry at most one
event each.

Two points cannot support an order or a shared-`Km` fit any more than
exps 60–62's three could, so none is attempted. What the pair DOES show,
for what a two-point comparison is worth: the debubbled rate falls from
pH 9.01 to pH 9.70, the same direction `PH_LADDER_BORIC` itself falls in
above its own pH 9.5 peak (§3) — a pointer in the same direction as the
4OMe finding, at the SAME enzyme loading, on the OTHER substrate, and
nothing stronger than that. `early_trough.arrhenius_check`'s own retraction
of a two-point activation energy is the standing reason not to read more
into two points than they can carry.

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
**+0.318 ± 0.131 per pH unit**, χ² = **1.06** on 3 — the four ladders agree,
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
- **The induction clock.** Section 4: +0.318 ± 0.131 per pH unit, pooled and
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
  pyrophosphate (χ² = 0.29 on 2, three ladders, two buffers, two substrates).
  Two of those three ladders **are** the two-axis block, so the independent
  corroboration is the phosphate ladder alone — +0.596 ± 0.032 on nine runs
  sharing no experiment with the block, **0.83σ** from its cuvette-matched
  +0.554 ± 0.040 (§2).
- The boric ladder does **not** share that order: it turns over near pH 9.5,
  and `vmax_corrected` rules out the O2 side reaction as the cause of the
  decline above it. The order below the peak is positive and weaker than the
  other three ladders' at every split tried, but its magnitude is set by where
  the ladder is cut (+0.06 to +0.48), so it is quoted as a range (§3).
- **The turnover is a mixture, and debubbling the fit flips which half
  dominates.** Per-experiment Michaelis-Menten fits, read off the debubbled
  rate (`v_peak_corrected`), resolve `Km` on 7 of the 9 boric runs. `Km` and
  `Vmax` turn over TOGETHER, peaking at the same rung (exp 48, pH 9.51),
  rather than `Km` climbing to the ladder's own top rung the way the
  uncorrected fit suggested — the two rungs the uncorrected fit was leaning
  on (exps 44, 49) are the ladder's two heaviest boric gassers. Priced apart, `Vmax` alone
  now falls FURTHER than the measured statistic (−56.6% against −47.4%, peak
  to last rung), because `Km`'s own fall past its peak partly offsets
  `Vmax`'s decline rather than compounding it. The ladder now supports one
  shared `Km` (5.87 mM, resolved) where the uncorrected fit could not — six
  of the seven resolved rungs' own intervals contain it, and only the
  lowest-pH rung (exp 41) disagrees (§3a).
- The induction clock's order in pH, +0.318 ± 0.131 per pH unit, agrees across
  all four ladders (χ² = 0.95 on 3) where the rate's does not.
- Three independent measurements — the rate, the clock, and the early
  trough's OTHER-cuvette signal — name [HOO⁻] as the species that matters, and
  a fourth (the O2 side reaction's pH onset) is consistent with the same
  reading.

**Not settled.**

- **Why boric turns over** is now two questions, not one, and §3a only
  separates them — it does not answer either. A buffer-specific effect on the
  catalyst (borate is a known Lewis-acidic complexing agent) and a genuine
  high-pH instability of the kind exp 85 already showed at pH 11.84 are both
  consistent with one ladder alone; the archive holds no other buffer above
  pH 9 to separate them, which is exactly `../buffer/ANALYSIS.md` §5's finding
  that buffer identity and pH cannot be told apart anywhere in this archive.
- **What `Km`'s own turnover means** is open in a different direction than it
  was. A borate–benzyl-alcohol ester lowering the *free* substrate
  concentration at high pH was floated as a candidate for a `Km` that kept
  climbing; debubbled, `Km` peaks and falls together with `Vmax` instead,
  which a simple free-substrate-sequestration story does not predict on its
  own. Nothing in this archive measures free versus bound alcohol, so this is
  a candidate revised, not a finding either way.
- **The other substrate's own boric ladder (§3b, §3c) cannot yet be compared
  to this.** Exps 60-62 show no O2 artefact to correct for and share the
  design this section reads, but two of its three runs do not resolve their
  own `Km` at all; the second enzyme tier (exps 51, 55) is two points, one
  of them the folder's own heaviest single-cuvette gasser. Whether BnOH's
  `Vmax`/`Km` turn over the same way 4OMe's do in the same buffer is open,
  not merely unmeasured to a lower precision — though the pair's own rate
  falling from pH 9.01 to 9.70 points the same direction 4OMe's turnover
  does past its peak.
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
