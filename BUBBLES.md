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
what chemistry can produce, full stop. It is the *only* place absolute σ is
the right currency, and §4.2 is where everything else asks the question
instead.

**The noise floor matters more than it looks.** Every function here takes
`noise` as an argument rather than recomputing it, because the floor differs
by 1096× between the instrument's own `.rre` readings and the `.txt` export's
0.001 AU quantisation (`fit_dataset.source_floor`), and every one of these
functions divides by it somewhere. Passing the export's floor on `.rre` data
— which is all 402 curves in this archive since 2026-08-31 — silently
under-counts detachments and suppresses every downstream z-score.

### 4.2 One score, and phases: `step_anomaly` and `bubble_segments`

**Rewritten 2026-09-11.** Until then this layer detected events *pointwise* —
a fall here, a rise there, each judged by its own local veto — and then tried
to recover the physical picture by stitching those points together
afterwards. The physics is **phases**: the beam alternates between
accumulating gas and releasing it, over runs of readings, and nothing in the
code represented one. Every defect found on 2026-09-10 was a case where the
missing phase structure had to be patched around, and the patches had reached
four rules and three constants, one of them already withdrawn.

**`step_anomaly(values, noise)`** is the single place the question "is this
step unusual" is now decided: one signed score per interval, the step in units
of `local_step_scale` — the median |step| in an 8-reading window either side
(`ANOMALY_WINDOW`), which is what a step there looks like with nothing
happening. Before this there were two currencies for one question,
`BUBBLE_DROP_SIGMA` in units of the curve's global noise and
`EXCURSION_LOCAL_SIGMA` in units of the local step, and nothing could compare
them. `noise` floors the scale, so a stretch the instrument reports as
identical readings cannot make its neighbours infinitely anomalous; that floor
binds on about a third of all intervals and changes no result in the archive.

**`bubble_segments(values, noise)`** returns the phases: maximal runs of
consistently-signed anomalous intervals, each an `accumulate` or a `release`.
The two kinds are nominated differently, and that asymmetry is the physics of
§4.1 rather than a fitted choice — a fall counts on amplitude alone
(`sigma` of the curve's noise, because chemistry cannot fall at all), a rise
only where it is anomalous against its own neighbourhood (`ANOMALY_BAR`, 2.0,
because chemistry rises constantly).

**A phase survives one reading of interruption** (`BUBBLE_SEGMENT_GAP`),
provided that reading gives back less than `BUBBLE_RECOVERY_FRACTION` of what
the phase has moved so far. A bubble sliding out of the beam can stutter, and
a rule that tolerated nothing cut those releases into fragments and read the
tick between two fragments as gas *arriving* — exps 43.1 and 135.2 are that
false arrival, and exp 135.1's readings 274–277 are the mirror, one arrival
split in two by a 6.6σ wobble in the middle of it. A bridged reading is
consumed by the phase that bridged it and cannot start one of its own.

The gap width is a **continuum, not a break**, and is stated as one: across
gaps of 0/1/2/3 readings the archive's detachment events run 381/364/352/343
and the readings those events span run 431/446/476/507. Nothing picks 1 out of
that curve. What picks it is the timescale argument the rest of this module
already rests on — one reading is the shortest interruption there is — plus the
cost of going wider being one-directional: past one reading an event swallows
ordinary readings faster than it gains falling ones, drawing a longer box
around the same gas.

`ANOMALY_BAR` is better placed. Among the 89 arrivals the archive carried
under the old absolute-σ rule, the 44 this bar throws out score **1.00 to
1.91** and the 42 it keeps score **2.06** upward, with nothing in between —
though that is a gap among rises that had already passed the kink test, not
among all rises, and swept over the whole archive the arrival count still
runs 75/58/44/43/33/32/24 at bars of 1.25 through 3.0. Both readings are
true; the second is the caveat. The value itself is unmoved from
`EXCURSION_LOCAL_SIGMA`, pinned between exp 149.5's two genuine excursions
(counter-moves at 7.2× and 2.1× the local step) and exp 130.2's two real
detachments (flanked by ordinary steps at 1.0× and 1.4×).

**`local_outlier_z` still cannot do a fall's job.** Its fitting window spans
the fall itself, so a genuine step change flags itself as an outlier — exp
135.2's real 0.1196 AU detachment scores +130σ under it. Nothing here ever
looks across an event; the rise side uses it only at a phase's own ends,
where a level jump is exactly what it is good at (§6).

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

### 4.4 Events, and rejecting spikes: `detachments` and `_excursions`

`detachments(values, noise, sigma, recovery, floor)` runs the pipeline: gate
the curve by `floor`, build the phases with `bubble_segments`, drop the ones
`_excursions` pairs off, and keep the `release` phases big enough to be an
event. **22 of 233 candidate falls in the two-axis block** are rejected as
excursions and two curves (144.7, 149.5) lose every one of theirs and are
returned untouched. Archive-wide it is **33 of 397** and five curves, adding
exps 3.1, 45.1 and 65.1.

**A spike is not one event that comes back — it is two adjacent phases that
cancel**, and saying it that way is what collapsed three rules into one. The
pair cancels when its smaller half is at least `BUBBLE_RECOVERY_FRACTION` of
its larger: the trace went somewhere and came back to about where it was.

**The order of the pair carries the physics, and it is the only asymmetry.**

| pair | reading | verdict |
|---|---|---|
| `release` → `accumulate` | gas cannot leave before it arrived, so a fall that is given back is a spike however long the give-back runs | both artefact |
| `accumulate` → `release`, rise grew over ≥ 2 readings | the ordinary bubble lifecycle — gas arrives and then all of it leaves, which cancels *exactly* | both real |
| `accumulate` → `release`, rise is one interval | a bubble that reaches full size inside one 60 s reading and detaches inside the next is not a bubble | both artefact |
| either order, does not cancel | two real events | both real |

Exp 44.1 is why the second row cannot be decided on size: it grows one bubble
over three readings, 0.0667 AU, and the next detachment sheds 0.0663 of it —
a near-perfect cancellation that is the lifecycle working. Exp 149.5's
readings 55–57 are the third row: 0.00096 AU up in one reading and 0.00183
back down in the next.

**A pair that does *not* cancel is two real events, and that is the case the
old code had no way to express.** Exp 138.4's real 2-reading detachment is
followed immediately by a real 4-reading arrival that overshoots the pre-fall
level by 0.0176 AU and holds. `_is_excursion` rejected the fall for being
followed by a rise, and would then have refused the rise for following a
rejected fall. One rule, correctly scoped, decides both.

**Two constants went with `_is_excursion`.** `EXCURSION_RECOVERY_DEPTH` (3)
reached past the adjacent reading for recoveries that landed one or two
readings late on exps 150.1 and 151.6 — both of which `DETACHMENT_SNR_FLOOR`
(§4.3) now excludes whole, so it had no live case left.
`EXCURSION_RECOVERY_CEILING` (1.0) capped a recovery at the drop's own size so
that a genuine acceleration just after a real fall was not read as that fall
reversing; a two-sided cancellation test needs no cap, because a counter-move
that *overshoots* does not cancel. Exp 135.1's 6.2σ fall at readings 272/273
— the curve accelerates 0.0312 AU over the four readings after a fall that
cost 0.0027 — survived the old test on a 2% margin and survives this one by a
factor of ten.

**Grouping consecutive falls was always safe** because chemistry cannot
produce a multi-reading run of large falls at all, so any such run is
unambiguously gas regardless of length. Since 2026-09-11 the grouping also
bridges a single reading of interruption (§4.2). Archive-wide that merges 14
pairs of adjacent events — exp 131 cuvettes 1 and 2, the heaviest bubblers
here, go from 18 and 19 events to **17 each over exactly the same falls** —
and costs exactly one falling reading anywhere, on exp 44.1 at reading 10,
which is not the grouping at all but a 16.7σ rise and a 14.8σ fall in
consecutive readings on a curve stepping 2.4σ. The old test kept that one on
a 3% margin. **This assumption does not transfer to rises** (§6).

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
It does: **+1.356 ± 0.256 in peroxide, -0.312 ± 0.089 in substrate** —
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

## 6. The rise model: `curve_metrics.bubble_arrivals` (added 2026-09-08)

A bubble does not only leave the beam abruptly — occasionally one *arrives*
abruptly: already formed, entering the optical path in one reading rather
than nucleating and growing inside it, or nucleating and growing fast over
two or three. Exp 135 cuvette 5's jump at 9780 s is the case that surfaced
this: readings 162 and 163 climb 0.00529 AU (21.2σ) in a single interval,
immediately after a reading that itself sits 8.4σ *below* the local trend,
and the curve resumes its exact pre-jump slope right after — the mirror image
of a detachment's shape, sign flipped.

**A rise cannot be read the way a fall is**, and this is the whole reason
`bubble_arrivals` is not `bubble_drops` with a sign flip. §4.1 rests on
chemistry being unable to produce a large *fall* at all; there is no
equivalent restriction on rises — real accelerating kinetics produces large
rises constantly. At `BUBBLE_DROP_SIGMA`'s own 6σ threshold, the two-axis
block shows **809 rises against 303 falls** — the *opposite* of the
122-against-23 asymmetry `bubble_step_asymmetry` reports at its own much
stricter 20σ threshold (§2). Most large single-reading rises in this archive
are the reaction, not gas.

**So an arrival is an `accumulate` phase (§4.2) that clears three bars.** It
must survive `_excursions`; it must move at least `sigma` of the curve's own
noise, so it is an event at all; and it must pass the kink test below.

**The kink test.** `local_outlier_z` — the same leave-one-out local quadratic
fit `isolated_outliers` uses, explicitly excluding the point being scored —
flags the reading a level jump departs from as anomalously *low* and the one
it lands on as anomalously *high*, because the fit is pulled between the two
levels it straddles. A genuine acceleration builds curvature over several
readings and does not do this. It is read at the phase's own ends, and only
across **the part of the phase that moved** (`_moved_span`): a phase can pick
up an ordinary reading or two at its edges through the gap rule, and reading
the kink there scores the wrong reading. Exp 146.4 is the case — its 35.7σ
jump is bridged to a 3.5σ step and a 1.2σ one before it, whose leading
reading scores only −1.5 against a bar of −5.

**Only the excess over an ordinary step counts as gas.** `local_step_scale`
gives what a step there looks like with nothing arriving, so
`gain = rise − n × local_step_scale(...)` removes only what is anomalous
beyond the curve's own ordinary movement — over the whole phase, not merely
its largest step.

### 6.1 Scoring rises on absolute σ was wrong, and it was wrong widely

Until 2026-09-11 candidate rises were nominated at 6σ of the curve's *global*
noise and scored one step at a time, never merged. Both halves of that were
wrong, and the archive went from **89 arrivals to 43**. The four ways an old
one could go are worth separating rather than summing:

| | what happened |
|---|---|
| **44** | **not anomalous at all** — no phase there |
| **31** | kept, unchanged |
| **11** | **merged** into a longer arrival, replaced by the 12 landings those arrivals actually end on |
| **3** | **inside a release phase** — not arrivals but upticks between two fragments of one stuttering detachment (exps 43.1, 135.2) |

**The 44 are the important number.** Scored on absolute σ they ran 6.5 to
160σ and looked convincing; scored against their own neighbourhoods they run
**1.00 to 1.91** — the curve was already climbing that fast. Exp 13 cuvette
4's jump at reading 64 is 43.8σ, and the six steps around it are 11.8, 11.8,
10.3, 10.9, 14.6 and 12.7 thousandths against its 18.4. It was credited
0.0067 AU of gas for rising 1.6× an ordinary step.

This is **exp 144 cuvette 2's documented failure mode turning out to be
general**. That curve's readings 29–42 climb 20–30σ a step for fourteen
consecutive readings, real and smooth; merged into one span its endpoints
still score as a level jump, and the "event" would have counted as a single
~0.034 AU gain, larger than any genuine gain in the block. The old rule kept
it out by refusing to merge *any* rise, ever — which protected that one curve
and left the same mistake standing on 27 others, one step at a time. Scored
locally, not one of its fourteen steps is anomalous (every one is about 1.0),
so it needs no special rule and merging becomes safe.

**The 11 are the other half of the same repair.** Exp 44 cuvette 1 grows one
bubble over three readings — +69, +96 and +75σ, totalling 0.0667 AU against
the 0.0663 the next detachment sheds — and the unmerged rule credited only
its largest single step, 0.0231, at reading 89. It is now 0.0569 at reading
90: present at the right reading, and finally the right event.

### 6.2 `split_arrivals`, `apply_gains`, and how they fold into `debubble`

A confirmed arrival's size is read directly off the jump — net of the curve's
own local step scale — so unlike a fall it needs no rate fit and no
bisection. **Which half it falls in decides which operator is correct**, and
getting that wrong was the whole of the 2026-09-10 repair:

- **released** — a detachment departs at or after the landing. The gas was
  held for a while and then went, so it belongs in `bubble_profile`'s `b(t)`,
  where the release takes it back out, and `bubble_rate` must see it or the
  rate is bid up to pay for a fall the arrival already covers.
- **unreleased** — no detachment follows. The gas was never watched to leave,
  so it is still in the beam at the last reading and a permanent shift from
  its own index on (`apply_gains`) is right. This is the direct analogue of a
  curve that "ends holding" an un-shed bubble in §5.

Until 2026-09-10 every arrival took the permanent shift, on 69 of the 80 the
archive then carried, so those curves were lowered for the whole of their
remaining length on account of gas they had demonstrably shed — on top of the
rate model already paying for the same fall. Exp 49 cuvette 1 carried four
such shifts. `debubble` is what routes each half to the operator that fits it.

### 6.3 Archive-wide scope

**43 arrivals over the whole archive**, on 24 curves. Restricted to the
two-axis block it is 15 curves. Almost all already carry confirmed
detachments — the heavy bubblers, unsurprisingly, since gas arriving and gas
leaving are the same underlying process seen from opposite ends. The
exception is **exp 146 cuvette 4**, which carries no detachment at all: a
bubble that arrived near the end of the run and never left before the
recording stopped — the "ends holding" case, but for an arrival rather than a
departure. `bubble_arrivals` is gated by the identical `DETACHMENT_SNR_FLOOR`
as `detachments`: exp 150.1 carries no arrival, for the same reason it
carries no fall.

### 6.4 One known limitation, bounded rather than hidden

A synthetic 120σ fall placed in the very first reading interval (an existing
test for the falls model's own "no affordable rate" edge case) distorts
`local_outlier_z`'s local fit for a few readings after it — the same
"masking" limitation `isolated_outliers` already documents for two adjacent
real spikes — producing a spurious 0.0006 AU gain there. That is four orders
of magnitude under the fall that caused it, and under anything found
anywhere on a real curve; the archive-wide sweep in `test_bubble_arrivals`
confirms no real curve triggers it. The affected test assertion was widened
from exact equality to a documented, explicit bound rather than silently
weakened.

### 6.5 The arrival set has a soft edge, and the bar is inherited

Almost every candidate rise is decided by **one number**: `z_before`, the
kink test's score on the reading the jump departs from. `z_after` clears the
bar comfortably on all of them, so the verdict turns on that single score
and the threshold it is compared against.

**That threshold does not sit in a gap.** `scope.arrival_margins` is the
table. Sorted, the deciding scores run smoothly through the bar: the closest
admitted sits at **−5.25σ** and the closest rejected at **−4.83σ**, a gap of
**0.41σ**, with candidates densely either side. Compare
`DETACHMENT_SNR_FLOOR`, which is defensible exactly because nothing in the
archive falls between 20.7 and 36.8 — there is no such break here.

**And the bar is inherited rather than calibrated for this question.**
`bubble_arrivals` takes its `kink_sigma` from `OUTLIER_SIGMA`, pinned for
`isolated_outliers`, which asks whether a single reading is suspect — not
whether a level stepped. Nothing has ever been calibrated against the
arrival question itself, because there is no break in the distribution to
calibrate against.

So a handful of admitted arrivals sit just past an arbitrary line, and a
handful of rejected ones just short of it. Nothing published rests on where
it falls — `bubble_sensitivity` moves no order under any repair — but **an
individual jump argued about one curve at a time has to be read against
this**, not against the bar alone. Exp 139 cuvette 2's rise at 3720 s is the
worked example: a 13.8σ step that misses admission at −4.63σ, argued over at
length, and unresolvable on that curve's own evidence in either direction.

### 6.6 The pair test measures cancellation in absolute absorbance

`_excursions` asks whether two adjacent phases cancel, comparing raw
absorbance while the chemistry climbs underneath. That comparison is biased,
and **in opposite directions for the two orders of pair**: for
`release → accumulate` the reaction *adds* to the apparent give-back, so a
real detachment on a fast-rising curve looks more like a spike; for
`accumulate → release` it adds to the rise and subtracts from the fall.

**Most of that bias is now absorbed upstream rather than compensated.** A
counter-move only becomes a phase at all if it is anomalous against its own
neighbourhood (`ANOMALY_BAR`, §4.2), and an ordinary trend-driven step is not
— which is exactly what keeps exp 130.2's two real 12.4σ and 16.6σ
detachments, whose neighbours score 1.0× and 1.4×. Before the rewrite that
job was done after the fact, by `_local_step_scale` inside the recovery test,
and only in one of the two directions.

**What is left was measured on 2026-09-11 rather than left as a worry, and
the measurement is again why nothing was changed.** Detrending the pair test
— comparing each phase's *gain* over an ordinary step there, instead of its
raw total — gains 4 detachments net on 5 curves and one arrival. It moves
nothing published. But it re-admits **exp 149 cuvette 5**, whose falls are the
documented instrument excursions that forced the recovery clause in the first
place (§4.4) — the curve where an unfiltered repair removed 0.0097 AU from a
trace that rose 0.0262, flattening a real early rise into a straight line.

The trade is the same one the old code faced and much smaller: 15 falls
against that 1 before the rewrite, 4 against it now. The bias is real,
measured, bounded, and deliberately left in place, because it buys nothing
published and cannot be separated from a known false positive.
DATA_VERIFICATION.md 2026-09-11 has the sweep.

### 6.7 One case the rewrite could have decided by accident, and did not

Exp 139 cuvette 2's rise at 3720 s (§6.5) sits two readings before a −4.1σ
fall. Under segmentation a fall that size could in principle have paired with
it and changed its verdict without anyone deciding to. It does not: 4.1σ is
under `BUBBLE_DROP_SIGMA`, so no release phase forms there, no pair exists,
and the jump is still rejected on the identical kink score it always had,
**−4.628σ** against the −5.0 bar. The case remains exactly as unresolvable as
§6.5 says it is, on the same evidence.

## 7. Downstream consequences

`vmax_corrected`, `tau_corrected` and `tau_slow_corrected` are the rate and
the progress fit's two time constants, read off the fully reconstructed
curve rather than the raw readings — computed in `scope.frame` from
`debubble`'s output, so §6 reaches all three automatically.

The most notable single consequence is on the induction "+1" test
(`induction.joint_clocks` — see `CLAUDE.md`'s "What the catalyst does
first"): a pre-equilibrium activation step predicts
`d ln v / d ln h - d ln τ / d ln h = 1` for the activating species, and the
two-axis block's peroxide axis is asked through `tau_slow_corrected`.

**This row has been misread twice, and the history is the point.** Asked of
the raw readings it sat 0.3σ from +1. The falls-only correction moved it to
1.4σ. Adding arrivals appeared to bring it back to 0.3σ — and that reading
stood for two days on a bug: an arrival a later detachment had already shed
was still being applied as a *permanent* level shift, depressing those
curves' tails and lengthening the clocks read off them (§6.2). Routed through
`split_arrivals` the row sits at **+0.774 ± 0.291, 0.8σ below +1** — between
the two earlier readings, nowhere near rejecting +1, and now moving in the
direction the artefact argument requires, since taking peroxide-made gas out
has to move `d ln v − d ln τ` *away* from the +1 it flattered.

The correction still moves individual curves — `tau_slow` differs from
`tau_slow_corrected` on 33 of 110 live curves and `tau` from `tau_corrected`
on 38 — and it buys resolution rather than costing it, 62 → 68 curves for
`tau` and 25 → 33 for `tau_slow`, because the artefact was what those fits
could not pin. It tightens the fast clock's fitted error (0.196 → 0.144) but
not the slow clock's (0.261 → 0.291), because arrivals change *which* curves
resolve, not only how many.
`data/test_scope.py::test_the_clocks_are_corrected_like_the_rate` asserts
that directly, as the true finding, rather than a threshold ("the correction
moves the axis by more than a rounding") that stopped being true once
arrivals were added.

Every other published number touching these curves moves by a
correspondingly small amount and none reverses direction. The 2026-09-11
segmentation rewrite (§4.2, §6.1) moved them again and by less:
`vmax_corrected`'s own peroxide order +0.704 → **+0.703**, the fitted gas
rate +1.343 → **+1.356 ± 0.256**, the `tau_slow` row +0.757 → **+0.774**, all
far inside their own errors. The itemised before/after lists are in
`DATA_VERIFICATION.md`'s 2026-09-08 and 2026-09-11 entries; they are not
repeated here because none of it changes the argument this document makes,
only its third decimal place.

**One number in this document was wrong for longer than that.** §7 quoted the
peroxide-saturation test's F statistic as "46, was 44"; it is **32.0**
(`induction.peroxide_saturation`, `first_order_f`), and was 32.0 before the
rewrite too. It still rejects `a = 1` decisively — that is what the test
says — but the figure had drifted in prose with nothing watching, in exactly
the way §8's own history describes. `test_root_documents.py` now carries a
claim on it.

## 8. Where this lives

| Function | Module | What it does |
|---|---|---|
| `bubble_drops` | `curve_metrics` | Raw candidate falls, amplitude test only |
| `local_step_scale` | `curve_metrics` | Median local step size, excluding one event |
| `step_anomaly` | `curve_metrics` | Every step in units of what a step looks like there |
| `bubble_segments` | `curve_metrics` | The phases: runs of consistently-signed anomalous steps |
| `_excursions` | `curve_metrics` | Which adjacent phases pair off as an instrument spike |
| `detachments` | `curve_metrics` | Grouped, excursion-filtered, SNR-gated falls |
| `unreleased_gas` | `curve_metrics` | The cap: gas still to be watched leaving |
| `bubble_profile` | `curve_metrics` | Builds `b(t)` from events and a rate |
| `bubble_rate` | `curve_metrics` | The least rate that pays for every fall |
| `debubble` | `curve_metrics` | `(A_obs - b) with gains applied, events)` |
| `bubble_arrivals` | `curve_metrics` | Confirmed rises: surviving accumulate phases, sized, kinked |
| `split_arrivals` | `curve_metrics` | Released against unreleased: which operator each needs |
| `apply_gains` | `curve_metrics` | The permanent level shift an UNRELEASED arrival applies |
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
