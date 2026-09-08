# The bubble artefact, and the machinery that removes it

Every catalysed curve in the two-axis block (`two_axis/`, exps 135-151) rides
on a second process the reaction itself does not produce: a gas that grows in
the light path and periodically detaches, chopping the trace with falls that
have no chemical explanation and, occasionally, jumps that don't either. This
document is the physics of why that happens, the evidence for what the gas
is, and the machinery in `data/curve_metrics.py` and `data/scope.py` that
separates it from the chemistry underneath. Everything quoted here is
reproduced by `two_axis/check_numbers.py` from the live code, or is in
`DATA_VERIFICATION.md`; nothing here is a source of truth in its own right.

## 1. What the instrument actually measures

A UV-Vis spectrophotometer does not measure absorption. It measures one
thing: how much light reaches the detector, expressed as a ratio against a
reference,

```
A = -log10(I / I0)
```

and reports that ratio as "absorbance," on the assumption written into the
Beer-Lambert law that the only thing standing between `I0` and `I` is
molecules absorbing photons via an electronic transition. The instrument has
no way to check that assumption. Anything else that keeps a photon from
reaching the detector — reflection, refraction, obstruction — is invisible to
it and gets folded into the same number, indistinguishable from real
absorption.

A bubble is exactly such a thing. It is a gas-liquid interface with a large
refractive-index mismatch (n ≈ 1.0 for the gas against ≈ 1.33 for the
aqueous solution), and at the size these bubbles reach it behaves as a
crude, uncontrolled lens rather than a small-particle (Rayleigh/Mie)
scatterer: geometric refraction and partial reflection at a curved surface,
not diffraction. The detector only accepts light travelling within a narrow
cone along the optical axis — in a clear, bubble-free cuvette essentially
every photon that enters on-axis exits on-axis and lands on the detector. A
bubble sitting in that path bends a fraction of those rays off-axis. Some of
the deflected light misses the detector's acceptance cone entirely: the
photons are not destroyed, they simply land somewhere other than the sensor.
Fewer photons arrive, `I` falls, and `A = -log10(I/I0)` rises — with no
chromophore having done anything. This is the general **turbidity/scattering
artefact** every cuvette-based photometric method is vulnerable to; it is not
specific to this instrument or this reaction.

The shape of the artefact follows directly:

- **A growing bubble gives a gradual climb.** As it grows it occupies a
  growing fraction of the beam's cross-section, intercepting more rays as it
  goes, so the apparent absorbance ramps up roughly in step with its size.
- **A detachment gives a sudden fall.** Once the bubble breaks free of its
  nucleation site and leaves the beam (rising out of it, or drifting off the
  optical axis), the full ray bundle is restored inside one 60 s reading
  interval, so `I` — and the readings — jump back up in one step.
- **An arrival gives a sudden rise** (§6): a bubble that is *already formed*
  entering the beam abruptly, rather than nucleating and growing inside it,
  occludes the beam in one reading instead of ramping in.

This is a dual-beam instrument (`CLAUDE.md`'s "sample beam"/"reference beam"
language), and that geometry is itself part of the evidence in §2: the
reference cuvette sees the same peroxide but not the same enzyme, so if the
gas came from peroxide alone, both beams would bubble, and they do not.

## 2. The evidence, and its limits

Nothing published claims to have *measured* a gas — there is no headspace
analysis, no manometry, no electrode anywhere in this archive. Everything
below is inferred from the shape and statistics of the optical artefact
itself, and the inference has a stated ceiling (end of this section).

**It scales with peroxide.** `scope.bubble_ladder` bands the block by
`[H2O2]`: 0 of 7 curves below 5 mM detach, 5 of 5 above 80 mM do, monotone
across six bands in between.

**It needs turnover, not just peroxide standing in solution.** Exps 136 and
137 sit at 73.4 mM — near the top of the ladder — and carry no detachment at
all (`bubble_turnover_control`). This reads as a control for peroxide alone
not being enough, but it is confounded with pH: those two runs also sit at
pH 6.95 and 7.53, and no run anywhere in the archive detaches below pH 7.5
regardless of turnover (§ below). `turnover_control_confound` checks whether
pH alone already predicts their silence — scaling the matched-peroxide run
exp 138 (73.4 mM, pH 8.16, 0.50 events/h) by first order in [HOO⁻] (the
shallowest pH dependence the survey supports) predicts 0.68 and 1.48 events
over exps 136 and 137's own run lengths. Zero is consistent with pH alone, so
this control does not independently establish that turnover is required —
only that its own silence is unsurprising either way.

**It is synchronised across neither instrument nor lamp.** `bubble_synchrony`
finds 23 coincident detachment times over 357 cuvette pairs against 21.3
expected by chance — indistinguishable from chance, so the events are not an
instrument or lamp artefact shared across cuvettes on a run.

**It needs enough peroxide to actually saturate the solution.**
`solution_chemistry.oxygen_budget` answers the "too small to see" objection
directly: at the top of the ladder the solution saturates on only 1.5% of
the peroxide present, and at the bottom of the ladder it cannot saturate at
all — so a real side reaction consuming a small fraction of the peroxide is
enough to produce a visible bubble at high [H2O2] and none at low [H2O2],
exactly the pattern `bubble_ladder` shows.

**It is in the enzyme-containing cuvette, and every run is its own
±enzyme control.** `bubble_step_asymmetry` reads the *sign* of the large
steps: a bubble growing in the sample beam climbs slowly and falls suddenly
(large falls should dominate); a bubble in the reference beam would do the
same to the reference, inverting the *sign* of what the ratio reports (large
rises would dominate). At 20σ the block shows 122 falls against 23 rises —
5.3:1, decisively the sample-beam signature. The reference cuvette omits only
the enzyme (it holds the same peroxide), so if the gas came from peroxide
alone in solution it would bubble too, and by this measurement it does not:
`gas_enzyme_control` checks the archive's directly enzyme-free runs at
matched peroxide and pH (exps 65, 67, 69, 70 — 16 curves at 122.4 mM, pH ≥
8) and finds none detaching, against roughly one event expected from the
same buffer's catalysed rate (p = 0.4 for observing zero) — consistent with,
though not by itself proof of, the enzyme's necessity.

**It is not confined to the two-axis block.** `scope.gas_curves` runs the
identical detachment test over all 402 curves of 88 experiments, keeping the
non-bubbling ones as the control population. It appears with both substrates
— 27 of 58 4OMe curves against 23 of 68 BnOH curves, catalysed, above 40 mM
and pH 8, in three different buffers (`gas_substrate_control`) — which is
exactly what a catalyst decomposing peroxide predicts, since that reaction
involves no alcohol at all. And it is strongly pH-gated: the archive's median
[H2O2] is 82.5 mM (278 of 402 curves sit above 80 mM), of which only 36
chop, and inside every buffer the detachment rate climbs with pH from a hard
floor of **zero** detachments in 270 hours of run time below pH 7.5, over 23
experiments (`gas_survey`). **pH is the trigger, not peroxide concentration
alone** — a curve at high [H2O2] and low pH is quiet, and a curve at
moderate [H2O2] and high pH is not.

**CO2 is close to excluded.** Carbonic acid's first pKa is 6.35, so above pH
8 essentially all dissolved carbonate sits as bicarbonate, which does not
leave solution as gas — and the detachment rate gets *more* common as pH
rises, which is backwards for a carbonate-release mechanism.

**What this does and does not establish.** The beam asymmetry says the gas
forms in the cuvette that holds the enzyme; it does not say *what* in that
cuvette produces it — a cyclodextrin host effect or a trace transition metal
contaminant would give the same peroxide order and the same pH dependence as
a ketone-catalysed disproportionation. **The ketone is not established as the
catalyst, and the gas is not established as O2.** The correct statement is
"a non-condensable gas, made in the enzyme-containing cuvette, growing with
peroxide turnover and rising steeply with pH" — O2 (from 2 H2O2 → 2 H2O +
O2, the catalysed disproportionation of the oxidant) is the *reading*, the
chemically obvious candidate given everything above, not a finding in its
own right. In `MECHANISM.md` this is **S4**, a fourth competing sink
alongside the three the seven-step scheme already carries: an unproductive,
catalyst-dependent drain on the oxidant that is fastest exactly where the
productive chemistry is strongest, and accounted for by no rate constant in
the mechanism as fitted. Where it sits in the catalytic cycle — the
peroxo-complex shedding O2 directly, a second peroxide equivalent attacking
that complex, or the reduced catalyst being reoxidised before it can turn
over the substrate — is open; no run in the archive moves the productive and
unproductive routes against each other, so the archive cannot separate them
(`COMPUTATIONAL.md` C10).

## 3. Why it cannot be stitched back in

The repair that suggests itself — add each detachment's size back onto every
reading after it — is exactly wrong, and it is worth being precise about why.
A bubble that costs `Δ` in absorbance when it leaves *contributed* `Δ` of
apparent rise while it was growing: that upward slope was never chemistry
either. Adding `Δ` back at the detachment restores the artefact's whole
upward half and removes none of it. On exp 135 cuvette 4, stitching ends the
curve at 126% of the most absorbance its own 0.219 mM of substrate could
ever produce (`bubble_mass_balance` — the mass-balance test this repair
fails outright, since a fraction above 1.0 of the diagnostic substrate
ceiling is not chemically possible). Tested against planted sawtooth
artefacts of increasing severity, stitching never beats *doing nothing at
all*: at severities 0.25×/0.5×/1×/2× of the curve's own rise, raw readings
recover a planted `vmax` at 1.14/1.26/1.62/2.27× the truth and stitching is
*worse*, at 1.16/1.33/1.66/2.34× (`scope.bubble_recovery`, bubbles emptying
completely and production still running at the last reading — the other
three combinations of emptying/ends-holding tell the same story). Stitching
is, in fact, a special case of
the model in §4 below: it is that model with its one rate parameter fixed at
zero, i.e. every bubble springing into being full-sized at the instant it
leaves rather than growing over time.

## 4. The falls model: `curve_metrics.debubble`

The model splits each reading into

```
A_obs(t) = f(t) + b(t)
```

a non-decreasing chemistry `f` and a non-negative gas `b` held in the beam,
made at a steady rate and released in whole detachments. `f = A_obs - b` is
what `debubble` returns.

### 4.1 Finding candidate falls: `bubble_drops`

`bubble_drops(values, noise, sigma=BUBBLE_DROP_SIGMA)` (6.0σ) flags every
reading-to-reading step that falls by more than `sigma` of the curve's own
noise. This needs only an amplitude test, and that is a genuine asymmetry
with the rise side (§6): the reaction is monotonic by construction
(benzaldehyde does not un-form), so *any* fall this large is already outside
what chemistry can produce, full stop.

**The noise floor matters more than it looks.** Every function here takes
`noise` as an argument rather than recomputing it, because the floor differs
by 1096× between the instrument's own `.rre` readings and the `.txt` export's
0.001 AU quantisation (`fit_dataset.source_floor`), and every one of these
functions divides by it somewhere. Passing the export's floor on `.rre` data
— which is all 402 curves in this archive since 2026-08-31 — silently
under-counts detachments and suppresses every downstream z-score.

### 4.2 Rejecting spikes: `_is_excursion`

Not every candidate fall is gas. A fall that a single adjacent reading undoes
by a comparable amount is a spike — noise, or a mixing transient — because
gas that leaves the beam does not come back, and a bubble cannot grow half
its own size again in one 60 s interval. `_is_excursion` tests this with
three clauses, in increasing order of how much history they were added to
handle:

1. **The original test**: the reading immediately before or after the fall
   recovers more than `BUBBLE_RECOVERY_FRACTION` (0.5) of the drop's own
   size, *and* that recovery is itself anomalous against
   `_local_step_scale` — the median local step size in an 8-reading window
   either side of the event, excluding the event itself
   (`EXCURSION_LOCAL_WINDOW`), scaled by `EXCURSION_LOCAL_SIGMA` (2.0). The
   second half of that clause exists because on a curve rising fast enough
   that ordinary steps are themselves a sizeable fraction of a modest,
   genuine detachment, the recovery-fraction test alone fires on real gas —
   it rejected two genuine 12.4σ and 16.6σ falls on exp 130 cuvette 2 before
   this half was added.
2. **The depth extension** (added 2026-09-07): on the block's two weakest,
   most drift-dominated curves (exps 150.1 and 151.6), a spike's recovery
   landed one or two readings later rather than at the very next one, and
   the single-reading test could not see it. `_is_excursion` now looks up to
   `EXCURSION_RECOVERY_DEPTH` (3) readings past a fall, crediting recovery
   only up to `EXCURSION_RECOVERY_CEILING` (1.0) × the drop's own size — capped
   so a genuine acceleration right after a real detachment (exp 135 cuvette
   1's 6.2σ fall at readings 272/273) is not read as that fall reversing
   itself. The reach is deliberately one-directional: extending it to the
   reading *before* a fall as well reads real pre-fall acceleration as a
   spike, which is what wrongly rejected exp 135 cuvette 1's genuine 41.3σ
   detachment during development.
3. **`local_outlier_z` cannot do this job.** Its fitting window spans the
   fall itself, so a genuine step change flags itself as an outlier — exp
   135 cuvette 2's real 0.1196 AU detachment scores +130σ under it.
   `_is_excursion` only ever looks at the readings on either side of an
   event, never across it.

### 4.3 The curve-level gate: `DETACHMENT_SNR_FLOOR`

Even with all three clauses, some curves are too weak for any per-event test
to resolve. Exp 150.1 (net/noise 20.7, barely over the `live` threshold of
20) has four candidate falls that survived every per-event refinement tried
against them — a second-difference local-noise estimate excluding every
other candidate on the curve, gated on having enough clean points to trust —
because on a curve this weak, real detachments elsewhere in the block (exp
140.4's fall at readings 55/56, exp 142.4's densely-packed run) score *below*
several confirmed excursions on the identical metric: no per-event threshold
separates them here, because the curve's own noise characterisation cannot
resolve "anomalous" at this signal level.

`DETACHMENT_SNR_FLOOR = 30.0` excludes a curve's detachments entirely below
that net/noise ratio, sitting in a genuine gap in the archive: nothing with a
candidate fall sits between exp 150.1's 20.7 and exp 131's two cuvettes at
36.8-44.6, which are genuine heavy bubblers (18 and 19 real detachments each)
the floor leaves untouched — and whose own `bubble_load` (6.5-8.3, §5) is as
high as exp 150.1's (5.4), so load alone cannot substitute for the SNR gate.
The floor also incidentally catches one dead curve's single mixing-transient
candidate (exp 66.3).

### 4.4 Grouping into events: `detachments`

`detachments(values, noise, sigma, recovery, floor)` runs the full pipeline:
gate the curve by `floor`, find raw candidates with `bubble_drops`, group
*consecutive* candidate indices into one event (a bubble sliding out of the
beam is not obliged to finish inside one 60 s reading — treating a
two-reading fall as two separate events dated one of them from a zero-second
growth window and left the second one's mass entirely uncorrected, in an
earlier version of this model), then drop any event `_is_excursion` rejects.
27 of 243 archive-wide candidate falls are rejected as excursions; two
curves lose every one of theirs and are returned untouched.

**Grouping consecutive candidates is safe for falls specifically because
chemistry cannot produce a multi-reading run of large falls at all** — any
such run is unambiguously gas, regardless of length. This assumption does
*not* transfer to rises (§6).

### 4.5 The mass balance: `unreleased_gas`, `bubble_profile`, `bubble_rate`

Given a set of falls events and one rate parameter, `bubble_profile` builds
`b(t)`, the gas held in the beam at every reading, under three clauses:

1. **It may not outrun the curve it rides on.** `b` grows by at most what
   the reading itself gained (`room = max(diff(values), 0)`), so
   `f = A_obs - b` can never fall across an *ordinary* step — this is the
   entirety of `f' ≥ 0` outside a detachment.
2. **It may not go negative.** A detachment takes `b` down by exactly the
   size of the fall and no further.
3. **It may not exceed `unreleased_gas`**: the total of the detachments
   still to come at that reading. A reading after the last fall therefore
   carries no gas at all in this model — the beam may hold at most what a
   *later, watched* release proves it was holding, never what an
   extrapolated rate merely predicts it holds. This is the clause that
   prevents a rate fixed by early drops from being carried unchecked across
   hours in which nothing detached (§5 has the cost of what this clause
   deliberately leaves uncorrected).

The gas does not reset to zero at each detachment — it is *reduced* by the
detachment's own size and carries the remainder forward, because one
detachment can empty one bubble of several growing at once:
`bubble_record(141, 3)` sheds its *largest* drop after its *shortest* growth
window (r = -0.91 between window length and drop size across its own
detachments), which a single bubble emptying each time cannot produce.

`bubble_rate` fixes the one free parameter by bisection (`BISECTION_ROUNDS`
= 80 halvings) to the *least* rate that pays for every detachment in full
(`bubble_shortfall(...) ≤ 0`) — an upper bound on what continuous production
this model can justify, never an inflated one. This fit is entirely
per-curve and entirely blind to composition: it sees only that one curve's
own detachment timing and sizes, never its `[S]`, `[H2O2]` or anything else
about what was in the cuvette.

**That blindness is what makes `gas_rate_drivers` an independent check, and
it works by pooling the fitted rate across many runs, not by fitting it
across them.** Each of the 46 detaching, rate-bearing curves in the block
already carries its own `bubble_rate` — one number, fit in isolation from
that curve's readings alone. `gas_rate_drivers` takes that whole set of 46
independently-fit numbers and regresses *them*, as data, against each
curve's own `[S]` and `[H2O2]` — the composition step happens entirely
*after* the per-curve fitting, on quantities the fitting itself never saw.
If the block's artefact really is the catalysed decomposition of the
peroxide, a rate measured this way — one curve's optical evidence at a time,
concentration-blind — should still come out first order in peroxide and
flat or negative in substrate once regressed across the runs that carry it.
It does: **+1.477 ± 0.258 in peroxide, -0.344 ± 0.093 in substrate** —
first-order-ish in the oxidant, mildly *negative* in substrate, exactly the
signature of a side reaction competing for the same catalyst rather than one
that consumes the alcohol. This is the strongest independent support the gas
argument has, because nothing about how each curve's own rate was fit could
have produced this pattern by construction — the concentrations were never
in the room.

### 4.6 `debubble` itself, and what it guarantees

`debubble(times, values, noise, sigma)` returns `(A_obs - b, events)` where
`b` is built from the least-paying rate. Its two hard guarantees, checked by
`scope.rebuild_smoothness`:

- **`worst_at_event` is zero or above on every one of 110 live curves except
  one.** Every one of the block's 216 detachments is corrected in full — the
  reconstruction never falls at a point the model itself calls gas. The one
  exception is exp 135 cuvette 6, whose fall sits in the *first* interval: no
  rate can be fit to a bubble that finished growing before the run started
  (there is no preceding rise to date it from), so `bubble_rate` returns
  `inf` and `debubble` returns that one curve untouched, worst step exactly
  -9.6σ.
- **`rebuilt_worst` is *not* the guarantee**, and reads -61.1σ somewhere in
  the block — that number belongs to the 27 rejected excursions, which are
  deliberately left in the curve, because laundering an instrument spike
  through a gas model removes real chemistry (§4.2's exp 149 cuvette 5 case
  lost a third of a genuine early rise this way, in an earlier version).

**On a curve with no detachment, `debubble` returns the readings
unchanged** — checked exactly, not merely closely, over the curves that
carry no fall and no gain (§6) at all.

**Against a planted truth** — real chemistry and real noise borrowed from
the block's own clean curves (curves carrying neither a fall nor a gain,
§6.3), a synthetic sawtooth artefact added on top — `scope.bubble_recovery`
sweeps all four combinations of how the planted bubbles behave. Where
production has stopped by the last reading (`ends_holding=False`, nothing
left un-shed) the repair recovers a planted `vmax` at 0.99/0.99/0.98/0.97×
the truth as the artefact grows from 0.25× to 2× the underlying chemistry
when each bubble empties completely, and 1.00/1.00/1.02/1.04× when bubbles
empty only partially (40-100%, matching `bubble_record(141,3)`'s own
behaviour) — both close to exact. Where production is still running at the
last reading (`ends_holding=True`, §5 has the reason), the repair instead
reads 1.00/1.01/1.15/1.67× and 1.02/1.02/1.16/1.80× respectively — the gap
between the two `ends_holding` values at each severity *is* the systematic
priced in §5. Stitching, over the same severity range, gives 1.16/1.33/
1.66/2.34×.

## 5. What the falls model cannot do: bracketed, not cured

**The bubble that never left is priced, not hidden.** `unreleased_gas`'s
third clause means a run still actively producing gas at its final reading
keeps the whole of that un-shed bubble in the reconstruction — the model
subtracts only gas it was *watched* to release. `curve_metrics.terminal_gas`
bounds what the beam may still hold at the last reading (the fitted rate over
the quiet tail, capped by how much the tail itself rose), and
`scope.terminal_bubbles` is the resulting table: 41 of 110 live curves carry
some terminal gas, 12 of them above a fifth of their own net rise.

That bound cannot by itself distinguish a run that ended *mid-bubble* from
one that had genuinely stopped making gas — it only asks the fitted rate, not
the readings. `curve_metrics.tail_excess` (the tail's own slope minus the
body's) can: exp 140.4 runs +7.4×10⁻⁵ AU/s faster in its tail than its body
(1.24× its own fitted gas rate — still making gas), while exp 149.4 runs
slower, and 8 of the 9 curves that ran past two full shedding intervals
without a detachment do too (consistent with production having stopped). The
comparison is one-sided by construction: an accelerating curve with no gas
in it ends steeper for reasons that have nothing to do with a bubble, so a
*positive* excess is evidence of ongoing production and a *negative* one is
not evidence against it. `vmax_terminal` is the far end of this bracket —
every bit of terminal gas charged to the artefact rather than the chemistry
— and nothing published lives inside the resulting range: across the whole
bracket the two concentration orders move by only +0.026 in substrate and
+0.002 in peroxide, against their own standard errors of 0.047 and 0.071.

**Some curves carry more artefact than chemistry, and are flagged rather
than excluded.** `bubble_load` divides the absorbance lost to detachments by
the curve's own net rise; above `BUBBLE_LOAD_CEILING` (1.0) more absorbance
left the beam than the reaction ever produced net, and no repair here
recovers a trustworthy rate from that curve. 14 of 110 live curves sit above
this ceiling — all four substrate rungs of exp 135, plus inner rungs of
exps 138, 140, 141, 142, 149 and 150 — and they stay in the frame, the live
counts, and the curves page; `bubble_sensitivity` confirms nothing published
moves under any repair variant because of them.

## 6. The rise model: `curve_metrics.bubble_gains` (added 2026-09-08)

A bubble does not only leave the beam abruptly — occasionally one *arrives*
abruptly: already formed, entering the optical path in one reading rather
than nucleating and growing inside it. Exp 135 cuvette 5's jump at 9780 s is
the case that surfaced this: readings 162 and 163 climb 0.00529 AU (21.2σ)
in a single interval, immediately after a reading that itself sits 8.4σ
*below* the local trend, and the curve resumes its exact pre-jump slope
right after — the mirror image of a detachment's shape, sign flipped.

**A rise cannot be read the way a fall is**, and this is the whole reason
`bubble_gains` is not `bubble_drops` with a sign flip. §4.1 rests on
chemistry being unable to produce a large *fall* at all; there is no
equivalent restriction on rises — real accelerating kinetics produces large
rises constantly. At `BUBBLE_DROP_SIGMA`'s own 6σ threshold, the two-axis
block shows **809 rises against 303 falls** — the *opposite* of the
122-against-23 asymmetry `bubble_step_asymmetry` reports at its own much
stricter 20σ threshold (§2). Most large single-reading rises in this archive
are the reaction, not gas, so a rise has to clear two tests a fall does not
before it counts as an arrival.

**Test one: recovery, reused rather than reinvented.** Negating the curve
turns a rise into a fall, so `_is_excursion(-values, event)` asks exactly
the question a rise needs: does the very next reading undo a comparable
amount? If so it is a spike — real fast chemistry, or noise — not gas that
arrived and stayed. This matters because the leave-one-out fit in test two
cannot tell the two apart on its own: a one-reading spike that reverts
immediately pulls its neighbours' local polynomial fits exactly the same way
a persistent jump does.

**Test two: the kink, past recovery.** `local_outlier_z` — the same
leave-one-out local quadratic fit `isolated_outliers` uses, explicitly
excluding the point being scored — flags the reading just before a level
jump as anomalously *low* and the reading the jump lands on as anomalously
*high*, because the fit is pulled between the two levels it straddles. A
genuine acceleration builds curvature over several readings and does not do
this *at any interior step*, even though the region as a whole would score
as a kink if it were tested as a single span (which is exactly why events
are never grouped — see below). `_is_excursion`'s own docstring calls this
same property the reason `local_outlier_z` "cannot be used" for a fall —
there, a step's anomalousness was never in question, only whether it
reverses; for a rise, past the recovery test, anomalousness is the only
question left, and it is exactly what this test answers.

**Only the excess over an ordinary step counts as gas.**
`_local_step_scale` — the same median-local-step estimator `_is_excursion`
uses for its own recovery test — gives what a step there looks like with
nothing arriving, so `gain = max(step - _local_step_scale(...), 0)` removes
only what is anomalous beyond the curve's own ordinary movement.

### 6.1 Why merging consecutive candidates — safe for falls — is wrong here

An early implementation copied `detachments`'s grouping strategy exactly:
merge consecutive rise-candidates into one span, test the span's own
boundary. On exp 144 cuvette 2, readings 29-42 climb 20-30σ a step for
*fourteen consecutive readings* — real, smooth, fast kinetics. Merged into
one span, its two endpoints still score as a level jump, because reading 42
is genuinely anomalously high relative to a local fit still anchored near
reading 29's level — the same signature a true one-reading jump gives. That
merged "event" would have counted as a single ~0.034 AU gain, larger than
any genuine gain found anywhere else in the block, straight out of a real
acceleration. `bubble_gains` scores every candidate step *individually,
never merged*: no gain now falls inside that curve's real 14-reading climb,
and every other confirmed gain elsewhere in the archive is unchanged.
`data/test_curve_metrics.py::test_bubble_gains` checks both the true
positive (exp 135.5, exp 146.4) and this specific false-positive case by
name.

### 6.2 `apply_gains`, and how it folds into `debubble`

A confirmed gain's size is read directly off the jump — net of the curve's
own local step scale — so unlike a fall it needs no rate fit and no
bisection: `apply_gains` is a permanent level shift, lowering every reading
from the gain's own index onward by exactly that amount. It is computed
independently of the falls model (different candidates, different tests, no
shared state), and `debubble` applies both:

```
reconstructed = apply_gains(A_obs - b, gains)
```

**Why the shift is permanent, with no release mechanism.** A fall has a
natural place to stop being owed — the detachment that pays for it. A gain
has no equivalent: there is no future event to date a release from, so
unlike the falls model's `unreleased_gas` cap, a confirmed gain simply stays
subtracted for the rest of the curve. This is the direct analogue of a
curve that "ends holding" an un-shed bubble in §5 — the mirror image, priced
the same way.

**This is provably additive, not just observed to be.** Because gains are
computed from the raw readings independently of the falls model's own rate
and events, `gas_at_end` (what a reconstruction still holds at its last
reading) equals a curve's own `gain_total` *exactly*, checked over all 110
live curves with zero mismatches — proving the falls component's own
end-of-run promise (`unreleased_gas`'s guarantee that nothing is held past
the last watched release) is completely untouched by adding gains, and the
gains component is purely additive on top of it. `worst_at_event` and
`rebuilt_worst` (§4.6) are both unmoved by this addition once the merge bug
(§6.1) was fixed.

### 6.3 Archive-wide scope

Restricted to the two-axis block: **15 curves, 26 confirmed gains** — exps
135.1 (2), 135.2 (4), 135.3 (3), 135.4 (1), 135.5 (1), 138.2 (1), 138.4 (2),
140.4 (1), 141.3 (2), 141.4 (2), 142.4 (3), 143.2 (1), 144.2 (1), 146.2 (1),
146.4 (1). All but one already carry confirmed detachments — the heavy
bubblers, unsurprisingly, since gas arriving and gas leaving are the same
underlying process seen from opposite ends. The exception is **exp 146
cuvette 4**, which carries no detachment at all: a bubble that arrived near
the end of the run and never left before the recording stopped — the "ends
holding" case, but for an arrival rather than a departure. `bubble_gains` is
gated by the identical `DETACHMENT_SNR_FLOOR` as `detachments`: exp 150.1
carries no gain, for the same reason it carries no fall.

### 6.4 One known limitation, bounded rather than hidden

A synthetic 120σ fall placed in the very first reading interval (an existing
test for the falls model's own "no affordable rate" edge case) distorts
`local_outlier_z`'s local fit for a few readings after it — the same
"masking" limitation `isolated_outliers` already documents for two adjacent
real spikes — producing a spurious 0.0006 AU gain there. That is four orders
of magnitude under the fall that caused it, and under anything found
anywhere on a real curve; the archive-wide sweep in `test_bubble_gains`
confirms no real curve triggers it. The affected test assertion was widened
from exact equality to a documented, explicit bound rather than silently
weakened.

## 7. Downstream consequences

`vmax_corrected`, `tau_corrected` and `tau_slow_corrected` are the rate and
the progress fit's two time constants, read off the fully reconstructed
curve rather than the raw readings — computed in `scope.frame` from
`debubble`'s output, so §6 reaches all three automatically.

The most notable single consequence is on the induction "+1" test
(`induction.joint_clocks` — see `CLAUDE.md`'s "What the catalyst does
first"): a pre-equilibrium activation step predicts
`d ln v / d ln h - d ln τ / d ln h = 1` for the activating species, and the
two-axis block's peroxide axis, asked through `tau_slow_corrected`, had sat
1.4σ from that prediction under the falls-only correction. Correcting gains
as well moves it to **0.3σ from +1**, matching what the *raw, uncorrected*
readings already said (0.3σ) almost exactly. This is not the two corrections
cancelling by coincidence: the falls-only correction itself had been part of
what pushed this row away from +1, because gains change *which* curves'
`tau_slow` resolves (exp 135.4 loses its resolved value once its own gain is
also removed; exps 138.2, 141.4 and 146.4 gain one — the resolved count is
unchanged at 34 of 110, but its composition is not), and the fitted
regression order depends on which curves are in that set, not only how many.
The correction still tightens the fast clock's fitted error (0.196 → 0.146)
but does *not* tighten the slow clock's (0.261 → 0.366) —
`data/test_scope.py::test_the_clocks_are_corrected_like_the_rate` asserts
this directly, as the true finding, rather than a threshold ("the correction
moves the axis by more than a rounding") that stopped being true once gains
were added.

Every other published number touching these 15 curves moves by a
correspondingly small amount and none reverses direction: `vmax_corrected`'s
own peroxide order (+0.794 → +0.696, was → +0.706), the saturation test's F
statistic (46, was 44), the enzyme-pair clock ratio, the substrate-pair
clocks, several pH-ladder and floor-sweep numbers in `induction/`. The full,
itemised list — with before/after values for every one of them — is in
`DATA_VERIFICATION.md`'s 2026-09-08 entry; it is not repeated here because
none of it changes the argument this document makes, only its third decimal
place.

## 8. Where this lives

| Function | Module | What it does |
|---|---|---|
| `bubble_drops` | `curve_metrics` | Raw candidate falls, amplitude test only |
| `_local_step_scale` | `curve_metrics` | Median local step size, excluding one event |
| `_is_excursion` | `curve_metrics` | Recovery test: is a fall (or, negated, a rise) a spike? |
| `detachments` | `curve_metrics` | Grouped, excursion-filtered, SNR-gated falls |
| `unreleased_gas` | `curve_metrics` | The cap: gas still to be watched leaving |
| `bubble_profile` | `curve_metrics` | Builds `b(t)` from events and a rate |
| `bubble_rate` | `curve_metrics` | The least rate that pays for every fall |
| `debubble` | `curve_metrics` | `(A_obs - b) with gains applied, events)` |
| `bubble_gains` | `curve_metrics` | Confirmed rises: recovery test + kink test, unmerged |
| `apply_gains` | `curve_metrics` | The permanent level shift a gain applies |
| `terminal_gas`, `tail_excess` | `curve_metrics` | The bracket on gas never watched to leave |
| `bubble_load` | `curve_metrics` | Fraction of net rise lost to detachments |
| `bubble_ladder`, `bubble_turnover_control`, `bubble_synchrony`, `bubble_step_asymmetry` | `scope` | The evidence in §2 |
| `oxygen_budget` | `solution_chemistry` | Is a small side reaction enough to see? |
| `rebuild_smoothness` | `scope` | The test a repair has to pass |
| `bubble_recovery` | `scope` | Recovery against a planted, known truth |
| `terminal_bubbles` | `scope` | The bracket, per curve |
| `gas_rate_drivers` | `scope` | What the fitted rate depends on |

## 9. What remains open

- **The gas has never been directly measured.** No headspace analysis, no
  manometry, no electrode. O2 is the chemically motivated reading of a
  non-condensable, pH-and-peroxide-dependent, catalyst-cuvette-localised gas
  — not a directly confirmed identity.
- **The ketone is not established as the catalytic species behind it.** A
  cyclodextrin host effect or a trace transition-metal contaminant would
  produce the same peroxide order and the same pH dependence.
- **Where S4 sits in the seven-step catalytic cycle is unresolved.** No run
  in the archive moves the productive route (turning over substrate) and the
  unproductive route (S4) against each other, so the archive cannot
  distinguish the peroxo-complex shedding O2 directly from a second
  peroxide equivalent attacking it, or from the reduced catalyst being
  reoxidised before it can turn over.
- **A gain event's arrival mechanism is not itself explained** — whether a
  bubble nucleated elsewhere in the cuvette and drifted into the beam, or
  detached from one site and immediately re-settled in the optical path, is
  not distinguishable from the absorbance trace alone.
