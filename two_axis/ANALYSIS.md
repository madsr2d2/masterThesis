# The two-axis block

Seventeen runs, seven cuvettes each, one substrate and one buffer and one
temperature. Exps 135–151 are the only place in this archive where **both**
concentration axes move *inside* a run, which is what the block is named for
and why the mechanism fitting was scoped to it.

What has not been written down is that they are also a **pH ladder**. Exps
136–142 carry one set of seven compositions at seven pH values and exps 143–151
another set at nine, so inside each group a cuvette can be followed from run to
run with only pH changing. That second design is the stronger one, and using
both is what this folder is.

    data/scope.py                     the block, its orders, its pH ladders
    python data/scope.py              the block in one paragraph
    python data/scope.py --design     the per-experiment design table
    python two_axis/build_figures.py
    python two_axis/check_numbers.py

**Figures**: [`index.html`](index.html) is the presentation — ten figures, A
to J, one per claim below. [`progress_curves.html`](progress_curves.html) carries
all 119 cuvettes with the form each earned, which is the audit surface for
every number here. The **46 curves carrying O₂ detachments** are drawn twice
there — the raw readings and the reconstruction, with a fit to each and a rule
at every detachment — so a fit to the gas cannot pass for a fit to the reaction.
A curve whose only falls were instrument excursions is drawn once, because
there is nothing there the gas model will touch.
The **14 panels** whose load clears 1 say so on their own face.

Related: `FITTING.md` owns the mechanism fitting and why the block was chosen
for it; [`../buffer/`](../buffer/ANALYSIS.md) §6 owns the perhydrate question
§4 runs into; [`../product_fate/`](../product_fate/ANALYSIS.md) owns the 4OMe
deceleration §5 contrasts with; [`../induction/`](../induction/ANALYSIS.md)
owns the induction §5 declines to import.

## 1. What the block is

**17 runs of seven cuvettes**: 119 curves, 110 of them live, BnOH at 25 °C in
pyrophosphate, pH 5.47 to 9.73 — 5.1 decades of [HOO⁻].

| exp | pH | [HOO⁻], mM | [S] ladder | [H₂O₂] ladder | run, s | live | accel. | late ÷ early |
|---|---|---|---|---|---|---|---|---|
| 135 | 7.70 | 0.0272 | 40.0x | 66.7x | 18780 | 7 | 4 | +0.09 |
| 136 | 6.95 | 0.00218 | 50.1x | 20.0x | 11280 | 5 | 3 | +1.16 |
| 137 | 7.53 | 0.00827 | 50.1x | 20.0x | 6540 | 7 | 3 | +1.03 |
| 138 | 8.16 | 0.0352 | 50.1x | 20.0x | 28740 | 7 | 2 | +0.50 |
| 139 | 8.50 | 0.0768 | 50.1x | 20.0x | 4560 | 7 | 4 | +0.64 |
| 140 | 9.22 | 0.396 | 50.1x | 20.0x | 3780 | 7 | 5 | +0.30 |
| 141 | 9.15 | 0.338 | 50.1x | 20.0x | 4680 | 7 | 4 | +0.53 |
| 142 | 9.43 | 0.635 | 50.1x | 20.0x | 5220 | 7 | 5 | +1.14 |
| 143 | 9.73 | 0.596 | 50.1x | 6.9x | 3000 | 7 | 5 | +2.32 |
| 144 | 9.32 | 0.238 | 50.1x | 6.9x | 4380 | 7 | 3 | +0.56 |
| 145 | 9.04 | 0.126 | 50.1x | 6.9x | 4380 | 7 | 5 | +1.46 |
| 146 | 8.66 | 0.0532 | 50.1x | 6.9x | 28740 | 7 | 1 | +0.39 |
| 147 | 8.02 | 0.0122 | 50.1x | 6.9x | 8280 | 7 | 2 | +0.49 |
| 148 | 7.60 | 0.00466 | 50.1x | 6.9x | 28740 | 7 | 2 | +0.59 |
| 149 | 7.10 | 0.00148 | 50.1x | 6.9x | 28740 | 7 | 0 | +0.08 |
| 150 | 6.26 | 0.000212 | 50.1x | 6.9x | 28740 | 4 | 0 | -0.01 |
| 151 | 5.47 | 3.38e-05 | 50.1x | 6.9x | 28740 | 3 | 2 | +1.37 |

### The run is an L, not a grid

Four cuvettes step the substrate at the run's top peroxide; four step the
peroxide at the run's top substrate; they share a corner — 63 and 63 live
cuvettes, and no cuvette off both arms. So the block moves both axes and **no
single cuvette moves both at once**. An interaction between them is not
identified anywhere in it, which is the same shape of gap that leaves the
perhydrate question open in `buffer/` §6: there, 0 of 88 runs cross `[buf]`
with `[H₂O₂]`; here, 0 of 119 cuvettes cross `[S]` with `[H₂O₂]`.

### The axes point in opposite directions

The concentration axes are **within-run**: **100.0%** of the log[S] variance
and **94.1%** of the log[H₂O₂] variance sits inside experiments, so a per-run
offset absorbs the day, the pH, the enzyme batch and the cell, and the orders
of §2 are read from contrast between cuvettes of one run.

The pH axis is the mirror image. Only **11.5%** of the log[HOO⁻] variance is
within-run, so it needs the opposite control — a per-*cuvette* offset — and §3
gives it one. Nothing else in the block varies: `[buf]` is 75.013 mM on every
one of the 119 curves.

Run length spans 3000 to 28740 s, a factor of 9.6x. **Nothing on this page is
read through a window given as a share of the run.**

## 2. The concentration orders

| parameter | fit | order in [S] | order in [H₂O₂] | n | R² |
|---|---|---|---|---|---|
| v0 | within-experiment | +0.374 ± 0.066 | +0.663 ± 0.093 | 106 | 0.677 |
| v0 | pooled | +0.384 ± 0.091 | +0.680 ± 0.124 | 106 | 0.246 |
| v_max | within-experiment | **+0.091 ± 0.052** | +0.794 ± 0.077 | 110 | 0.811 |
| v_max | pooled | +0.079 ± 0.093 | +0.768 ± 0.133 | 110 | 0.264 |
| net | within-experiment | +0.199 ± 0.041 | +0.722 ± 0.061 | 110 | 0.784 |
| net | pooled | +0.179 ± 0.066 | +0.683 ± 0.094 | 110 | 0.335 |
| gain | within-experiment | -0.228 ± 0.052 | +0.135 ± 0.072 | 106 | 0.498 |
| gain | pooled | -0.207 ± 0.059 | +0.100 ± 0.080 | 106 | 0.185 |

Two things stand out.

**The maximum rate has no substrate order.** Over a fifty-fold ladder it is
+0.091 ± 0.052 — inside two standard errors of zero — while the *initial* rate
over the same cuvettes has a real one at +0.374 ± 0.066. The rate the run
climbs to has stopped seeing the substrate; the rate it starts at has not.

The obvious objection is depletion — a cuvette that runs out of substrate
cannot show an order in it. It does not apply here. **No live curve in the
block converts more than 37.6% of its substrate**, and at the bottom rung, the
one that would empty first, the median conversion is 6.7%. The rate is flat in
substrate while the substrate is still there.

This is the one block where that number means what it says. Everywhere else in
the archive substrate was added by volume and displaced buffer, so `[S]` and
`[buf]` move together at −0.96 in logs and a substrate order measured there is
an order in the pair (`induction.composition_collinearity`). Here `[buf]` is
fixed and `[S]` moves alone.

**And the peroxide order is well below one**: +0.794 ± 0.077 on `v_max`. Fitted
as a free power over the peroxide ladder it is 0.654 (0.534-0.775), and strict
first order is rejected at **F = 32** over 2.447-163.165 mM. "First order in
H₂O₂" is not what this block says.

### Read from the arm that moves only its own axis

The joint fit above reads both coefficients from all seven cuvettes at once.
That is correct and it is also an extrapolation: the L has no interior, so the
fit's separation of the two axes rests on the response being additive in logs.
Each arm holds the other axis fixed within a run and needs no such assumption.

| parameter | arm | that arm alone | joint fit | σ apart |
|---|---|---|---|---|
| v0 | substrate arm | +0.463 ± 0.078 | +0.374 ± 0.066 | 0.9 |
| v0 | peroxide arm | +0.526 ± 0.059 | +0.663 ± 0.093 | 1.3 |
| v_max | substrate arm | +0.154 ± 0.054 | +0.091 ± 0.052 | 0.8 |
| v_max | peroxide arm | +0.655 ± 0.061 | +0.794 ± 0.077 | 1.4 |
| net | substrate arm | +0.241 ± 0.042 | +0.199 ± 0.041 | 0.7 |
| net | peroxide arm | +0.624 ± 0.050 | +0.722 ± 0.061 | 1.2 |
| gain | substrate arm | -0.264 ± 0.052 | -0.228 ± 0.052 | 0.5 |
| gain | peroxide arm | +0.130 ± 0.053 | +0.135 ± 0.072 | 0.1 |

The worst disagreement is **1.4 sigma**. The additivity the joint fit assumes
is not costing anything measurable, and the sub-first-order peroxide result
survives on the arm that never touches the substrate.

### Which runs may carry an argument

`scope.concentration_agreement` correlates each run's observed rates against
the rates its own cuvette concentrations imply. Above the floor of 0.70 a run
predicts itself; below it, its cuvettes differ by something other than what was
put in them. **11 of 17** pass. The five that fail are 136, 137, 147, 149, 150
— and exp 151 does not appear at all, its seven cuvettes scattering with two
negative rates. §2's orders are measured over all of them, which is the
conservative direction; §3 and §4 are not, and §3 says what that changes.

## 3. The pH axis, at matched composition

The block holds **two** composition sets, not one: exps 136–142 and exps
143–151. They share the substrate ladder exactly — 0.216, 0.865, 3.028 and
10.816 mM — and differ only in their four peroxide levels, 3.671 to 73.424 mM
against 5.140 to 35.244. Inside either set a cuvette can be followed across the
runs of a series with only pH changing, and a per-cuvette offset holds its
substrate, its peroxide and its position in the holder.

The **enzyme** is what splits the series. It is not constant over the block —
0.014, 0.021, 0.034, 0.069 mM — and it steps between runs, exactly where pH
does, so a per-cuvette offset cannot absorb it. Exps 141 and 142 sit at
0.014 mM against exps 136–140's 0.034 on the same composition and at the top of
that group's pH range, so pooling them would suppress a rate that scales with
enzyme precisely at high pH and bias the slope down with it. A ladder is
therefore runs sharing a composition **and** a loading, which leaves
**2 ladders** over the strong runs: exps 143, 144, 145, 146, 148 and exps 138,
139, 140.

| ladder | runs | curves | pH | d ln v_max / d ln [HOO⁻] | R² |
|---|---|---|---|---|---|
| 0.021 mM chemzyme | 5 | 35 | 7.60-9.73 | +0.522 ± 0.047 | 0.877 |
| 0.034 mM chemzyme | 3 | 21 | 8.16-9.22 | +0.697 ± 0.068 | 0.958 |
| pooled | 8 | 56 | 7.60-9.73 | **+0.554 ± 0.040** | 0.911 |

**A half order in [HOO⁻], not a first order.** Doubling the hydroperoxide anion
by moving pH buys about 1.5× the rate, not 2×.

Over all seventeen runs rather than the strong eleven the pooled slope is
+0.400 ± 0.034, and the two ladders disagree by 5σ instead of 2σ. The weak runs
**flatten** the ladder rather than sharpening it, and they do it at the end
where they live: exps 149, 150 and 151 sit at the bottom of the pH ladder, and
what they measure there is the cell's own wander rather than a rate, so the
ladder stops falling when the chemistry does.

### The control the two ladders give each other

The block was collected 3 to 14 September 2010, and experiment number is
chronological. The two ladders were **not run in the same direction**: exps
138–140 climb in pH as the days pass, exps 143–148 descend. pH correlates with
the schedule at +0.75 in one and -0.74 in the other.

So anything that drifts monotonically over those twelve days — an enzyme stock
losing activity, a lamp, a peroxide stock decomposing — enters the two ladders'
pH slopes with **opposite signs**. The control needs no assumption about what
drifts or how fast; it only needs both slopes to have the same sign.

Both are positive. No monotone schedule effect produced the pH order.

## 4. Is [HOO⁻] the reactant?

[HOO⁻] = [H₂O₂]·Kₐ/([H⁺] + Kₐ), so there are two levers on it and the block
moves both — peroxide within a run at fixed pH, pH between runs at fixed
peroxide. Nothing is shared between the two measurements: different contrast,
different offsets, different runs carrying the signal. **If the chemistry
consumed the hydroperoxide anion and nothing else, they would give the same
number.**

| statistic | [H₂O₂] at fixed pH | pH at fixed [H₂O₂] | σ apart |
|---|---|---|---|
| v_max | +0.762 ± 0.046 | +0.554 ± 0.040 | **3.4** |
| v0 | +0.569 ± 0.045 | +0.291 ± 0.055 | **3.9** |

They do not. Peroxide is the stronger lever in both statistics, by 3.4σ and
3.9σ. Adding H₂O₂ buys more rate than raising pH does at the same [HOO⁻], so
either **[H₂O₂] contributes beyond its own hydroperoxide content**, or **pH
costs something that opposes it**.

Two readings survive and the block cannot separate them.

- **A second peroxide term.** The buffer perhydrate of `buffer/` §6 is one:
  its rate carries `[buf][H₂O₂]`, and `[buf]` is constant here, so it would
  appear on this block as an extra dependence on [H₂O₂] alone with a pH
  dependence set by the *perhydrate's* pKₐ rather than by H₂O₂'s. That is
  exactly the shape of the gap.
- **Saturation.** The two contrasts sit at geometric-mean [HOO⁻] of 0.0329 and
  0.0499 mM, and the higher-level contrast is the one with the lower order —
  the direction saturation predicts. A 1.5× difference in level is a thin thing
  to hang 0.2 of an order on, but it is not nothing.

What would settle it is a run that crosses `[buf]` with `[H₂O₂]`, and the
archive has none — the same missing experiment `buffer/` §6 ends on.

## 5. What the curves do

### The chop is gas, and it is measured rather than excluded

Many of the block's curves rise, fall by more in one 60 s reading than the
reaction moves in five, resume at the new level, and do it again. Exps 143.3,
140.4 and 144.2 are ordinary examples and exp 135 does it twenty times a
cuvette. **Absorbance that goes away was never product** — benzaldehyde does not
un-form — so the falls are something leaving the light path, and the rises
before them are the same thing arriving.

It is O₂, from the ketone-catalysed decomposition of the peroxide
(`MECHANISM.md` refs 34–35), growing on the window and detaching. Three things
say so and no two of them could be arranged by the same accident.

**It ladders with peroxide.** `scope.bubble_ladder`:

| [H₂O₂] mM | live curves | carrying a detachment | mean absorbance lost |
|---|---|---|---|
| 0–5 | 7 | **0** | 0.0000 AU |
| 5–10 | 14 | 2 | 0.0005 |
| 10–25 | 18 | 3 | 0.0007 |
| 25–40 | 38 | 21 | 0.0071 |
| 40–80 | 28 | 19 | 0.0469 |
| 80–200 | 5 | **5** | 0.2916 AU |

**But peroxide alone will not do it.** Exps 136 and 137 sit at the same 73.4 mM
as exps 138–142 and carry **no detachment at all**, and they are the two weakest
runs at that peroxide — `concentration_agreement` 0.25 and 0.21 against 0.83 to
0.97. `scope.bubble_turnover_control`.

**That control is confounded with pH, and this document read it too strongly
until 2026-09-04.** Those two runs also sit at **pH 6.95 and 7.53**, third and
fifth lowest in the block, and across the whole archive *nothing* detaches below
pH 7.5 — 270 hours of catalysed, high-peroxide running over 23
experiments (`scope.gas_survey`). Priced at first order in [HOO⁻] from exp 138,
which carries the same 73.4 mM at pH 8.16, their silence is worth **0.68 and
1.48 expected events** against zero seen (`scope.turnover_control_confound`).
pH alone predicts the same zero, so this pair does **not** establish that the
gas needs turnover. What does is the beam asymmetry below.

**And it is not the instrument.** A lamp, a shutter or the carousel would take
all seven cuvettes of a run at the same reading. Over the block's 357 cuvette
pairs, detachments coincide **17 times against the 16.0 independence predicts**
(`scope.bubble_synchrony`), so each detachment belongs to one window.

### The curves may not be stitched back together

The repair that suggests itself — add each step back to everything after it —
is the one repair guaranteed to make things worse, and the reason is
conservation. **A bubble that costs δ when it leaves must have contributed δ of
rise while it grew.** Stitching removes the artefact's downward half and keeps
all of its upward half, so it inflates the curve by exactly the sum of the
drops.

`scope.bubble_mass_balance` is that stated against the substrate. Stitched, exp
135 cuvette 4 ends at **1.26× the most absorbance its 0.219 mM of substrate
could ever make**, and four more land between 0.49 and 0.86 — while their ramps
*steepen*, at 1.1 to 23.5 times their starting slope. A reaction that has spent
half its substrate cannot be running twenty times faster than it started.

**Take the gas out instead.** Write the readings as a sum of two things and
give each one the only property it cannot lack: `A_obs = f + b`, where `f` is
the chemistry and never decreases, and `b` is the gas and is never negative.
Gas is made steadily — the peroxide is in enough excess that its decomposition
does not slow over a run — and leaves in whole detachments, and two clauses
close the model:

- **the gas may not outrun the curve it rides on.** Over an ordinary step `b`
  grows by at most what the reading itself gained, which is exactly `f' ≥ 0`.
- **a bubble cannot shed gas that was never made.** The production rate is set
  to the least that pays for every detachment, so `b ≥ 0` throughout.
- **the gas may not exceed what leaves.** The beam holds at most the total of
  the detachments still to come, so a reading after the last fall carries no
  gas at all — `curve_metrics.unreleased_gas` and the section below.

That leaves one free parameter, the production rate, and it is *pinned* rather
than fitted — the least rate consistent with what was seen to leave, which
makes the result an upper bound on the chemistry, the same direction as the
monotone bound. `curve_metrics.bubble_profile` and `curve_metrics.bubble_rate`.

**Stitching is this model at a rate of zero** — every bubble springing into
being full-sized at the instant it leaves. That is not a caricature of it, it
is the same arithmetic with the one parameter set to the one value the
mass balance excludes.

Against sawtooths planted into the block's own clean curves, where the rate
before the planting is the truth. Planted twice, because the difference between
the plantings is the difference between the two repairs: with each bubble
emptying at its detachment, and with each emptying only partly and carrying the
rest over, which is what one curve below actually does. The last column is the
same planting stopped at its final release, and the gap between the last two
columns is the third clause's whole cost — see below.

| artefact ÷ chemistry | left alone | **stitched** | rebuilt, ends holding | **rebuilt, shed out** |
|---|---|---|---|---|
| 0.25 | 1.12 | **1.15** | 1.00 | **0.99** |
| 0.5 | 1.24 | **1.32** | 1.00 | **0.98** |
| 1 | 1.58 | **1.64** | 1.15 | **1.00** |
| 2 | 2.26 | **2.34** | 1.65 | **0.99** |

*(each bubble empties; recovered `vmax` ÷ true `vmax`, 1.00 exact)*

| artefact ÷ chemistry | left alone | **stitched** | rebuilt, ends holding | **rebuilt, shed out** |
|---|---|---|---|---|
| 0.25 | 1.12 | **1.12** | 1.03 | **1.01** |
| 0.5 | 1.26 | **1.31** | 1.04 | **1.02** |
| 1 | 1.51 | **1.64** | 1.17 | **1.02** |
| 2 | 2.22 | **2.32** | 1.79 | **1.07** |

*(each bubble empties only partly)*

Stitching never beats doing nothing, and from 0.5x on it is strictly worse
under both plantings; at the smallest artefact with the bubbles emptying
only partly the two are equal to four decimals, because barely a drop clears
the detector and there is nothing for stitching to add back.

**On gas that was seen to leave the reconstruction is exact** — within
**0.07 of the truth at every severity under both plantings** — where the
segment ramp it replaced held to 1.01 while the artefact was small and reached
1.60 by 2×. What separates the last two columns is entirely the bubble that
never detached, and that is the third clause stating its price rather than
assuming it away. `scope.bubble_recovery` is all four columns and
`data/test_scope.py::test_the_correction_recovers_a_planted_rate` asserts them.

**On a curve with no detachment it returns the readings unchanged** — not
close, identical — so nothing here can move a clean curve, and 60 of the
block's 110 live curves are clean.

### Every detachment is corrected, and nothing else is touched

A repair is only a repair if the falls it claims to explain are gone from what
comes out. That is the test the segment ramp failed, and it failed it two ways.

It dated each bubble from the previous detachment and made it reach the full
size of the drop, so a large drop after a short window implied a bubble growing
faster than the curve: across exp 141 cuvette 3's last six readings it
subtracted **+0.0061 per reading from a trace rising +0.0018**, and the
"corrected" curve fell at every one of them. And it treated a fall spanning two
readings as two bubbles, gave the second a growth window of **zero seconds**,
skipped it, and left the whole of it in place — −0.0165 AU at 60σ in exp 144
cuvette 2.

Both are gone. Over the block's **180 detachments**, every one is corrected in
full: `scope.rebuild_smoothness`'s `worst_at_event` is zero or above on all 44
repairable curves. The single exception is exp 135 cuvette 6, whose fall is in
the *first interval* — a bubble grown before the run began leaves no rise to
date it from, no rate explains it, and `debubble` returns that curve untouched
rather than guessing.

**`rebuilt_worst` is not the guarantee, and it should not be.** The
reconstructions still fall by up to −61.1σ somewhere, and the curves doing that
are behaving correctly: those are the **34 falls rejected as instrument
excursions**, which the gas model deliberately leaves where they are. The next
section is why.

### A fall that comes straight back is not gas

Exps 149.1, 149.2 and 149.5 came out pulled far too low, and the cause was not
the repair but the *detector*: it was correcting falls that were never gas.

Exp 149 cuvette 5 has two, at 9.3σ and 8.2σ. The first falls 0.00206 AU and the
**next reading climbs 0.00222 straight back** — readings 8 and 9 are a spike up
followed by a spike down, and reading 10 resumes the trend exactly. The second
falls off a reading that is an isolated spike. Between them they set a
production rate the curve has no business carrying, and the repair then removed
**0.0021 AU
from a curve that rose 0.0262** — flattening a real early rise,
while staying perfectly monotone and passing every test there was. When the
fault was found on 2026-09-03 it was worse still, 0.0097 AU or a third of the
curve, because the tail rule then in force let that rate run on to the end of
the run; the clause below has since taken most of it back. **Neither number is
gas**, and that is the point: a repair that launders an instrument spike
through a gas model is wrong at any size.

Absorbance that goes away because gas left the beam does not come back, and at
60 s sampling a bubble cannot grow half its size in one interval. So a fall
that a single adjacent reading undoes by more than half is an excursion, not a
detachment — the same timescale argument `curve_metrics.isolated_outliers`
rests on, applied to the other artefact. The block separates cleanly on it: its
ten largest detachments have the readings either side moving −0.001 to 0.110 of
the fall, while the excursions sit at 0.8 to 1.66.

**`local_outlier_z` cannot do this job**, though it is the obvious tool. Its
window spans the fall, so a genuine step change flags itself — exp 135 cuvette
2's 0.1196 AU detachment scores **+130**, and using it would have rejected 200
of 221 falls including every large one. The excursion test looks only at the
two readings immediately either side, which no step change can make anomalous.

Of 214 candidate falls, **34 are rejected** and 180 kept. **5 curves** lose all
of theirs and are returned exactly as they were read, exps 149.1 and 149.5
among them. Nothing is deleted: an excursion stays in the readings, visible as
the instrument problem it is, and `isolated_outliers` is what nominates those.

### Only gas that was watched to leave is subtracted

Monotone is necessary and not sufficient, and the fault this cures was
invisible to the test above: pulling a curve down by a smooth ramp leaves it
smooth. The production rate is fixed by the detachments, and where those are
early and the run is long it was extrapolated across hours in which nothing
detached.

Exp 149 cuvette 4 sheds **0.0031 AU once**, 1920 s into an 8 h run, and nothing
else ever leaves. An unbounded profile keeps making gas at the rate that one
detachment implies, so by the end it has **0.0273 AU** sitting in the beam —
8.8× the largest bubble that curve ever shed — and the reconstruction reads
**−0.0209 against a raw rise of +0.0064**. Twelve of the 49 detaching curves
held more than twice their own largest bubble, one at 26.6×, and three finished
below zero. All twelve passed the smoothness test.

Capping the tail at the most the beam had carried *bounded* that without curing
it, and the residue was visible by eye on this page: exp 149 cuvette 3 sheds
0.0022 AU once at 5100 s, runs on for another **23 580 s with nothing leaving**,
and the reconstruction sat that flat 0.0022 AU below its own readings for 82% of
the run — under a last arm that shows no sign of a bubble at all. Exp 150
cuvette 1 sat below its readings by **99% of everything it rose**.

**The readings themselves refute the gas that was being subtracted.** If exp 149
cuvette 4 had gone on making gas at its fitted rate, the tail would have carried
0.0736 AU of it; that tail rises 0.0041 in total. On **18 of the 44** repairable
curves the rate would have made more gas over the quiet tail than the trace rose
altogether, and a bubble that keeps growing detaches — **11 of 44** ran more
than a full shedding interval past their last detachment without one
(`curve_metrics.quiet_tail`; exps 149.4 at 14.4, 150.1 at 7.7, 146.1 at 6.9,
149.2 at 5.7). Those runs stopped making gas, and the last thing they did was
shed, so their beams are empty.

So the third clause is evidence and not extrapolation: **the beam holds at most
the total of the detachments still to come**, and after the last one, nothing.
`curve_metrics.unreleased_gas`. The cap never blocks a detachment — at the
reading before the *j*th fall it is *d<sub>j</sub>* plus everything later, which
already exceeds what *d<sub>j</sub>* has to pay — so no rate changed and no
detachment went uncorrected. Every reconstruction now lands back **on** the
readings at the last reading, which is checkable by eye on the page opposite.

The price is stated rather than assumed. A run still making gas when the
recording stopped keeps the bubble it never shed, and that is the whole of the
gap between the last two columns of the recovery tables above: exact on gas that
went, blind to gas that did not. `quiet_tail` says which curves could be paying
it, and on this block the tails are long.
`data/test_scope.py::test_the_gas_may_not_outlast_the_evidence` names the curves
so it cannot come back quietly.

### The bubble that never left, and what keeping it costs

The clause above has a price, and it is visible on the page opposite rather
than only in a recovery table. A run that was still growing a bubble when the
recording stopped is handed its tail back untouched, so the reconstruction can
climb at one rate for hours and then, from the last detachment onward, at
another. **Exp 140 cuvette 4 is the one the eye catches**: corrected down for
3180 s, then lying on its own readings for the last 600, over which it rises
faster than at any earlier point in the run.

`scope.terminal_bubbles` measures what is left in. The bound is what the fitted
rate would have made since the last detachment, capped by what the tail
actually rose — gas may not be taken off faster than the trace climbed
(`curve_metrics.terminal_gas`). **Thirty-eight of the 110 live curves** carry
one, **13** of them more than a fifth of everything they rose.

| curve | load | bound, AU | of its rise | tail excess, AU/s | ÷ its own gas rate |
|---|---|---|---|---|---|
| 149.4 | 0.49 | 0.0041 | **64%** | −4.9e-07 | −0.18 |
| 135.4 | 8.67 | 0.0179 | 51% | +1.4e-05 | 0.27 |
| 140.4 | 1.24 | 0.0322 | 32% | +4.4e-05 | **0.83** |
| 142.4 | 2.96 | 0.0128 | 32% | +3.4e-05 | 0.75 |
| 142.2 | 1.27 | 0.0118 | 16% | +3.0e-05 | **1.08** |

**The bound alone cannot tell the two endings apart**, because it asks the
fitted rate and not the readings — it charges a run that stopped making gas
exactly as it charges one that ended mid-bubble, which is planted in
`test_curve_metrics`. The tail's own slope can. A run still growing a bubble
hands back a stretch **steeper** than the body it was corrected against, and a
run that stopped hands back one that is not: exp 140 cuvette 4's tail runs
4.4e-05 AU/s faster, which is 83% of that curve's own fitted gas rate, while
exp 149 cuvette 4's runs **slower**. Of the 24 curves with a positive excess
the median is **1.15** of their own rate — the statistic recovers the gas rate
it never saw. And the long silences line up with the quiet ones: of the **8**
curves that ran more than two shedding intervals without a detachment, **7**
come back negative, which is `unreleased_gas` being simply right on them.

`curve_metrics.tail_excess`, and it is **one-sided**. A curve's own curvature
is in it: 51 of these 110 accelerate past 3σ and an accelerating curve ends
steeper than it began with no gas in it at all, while a decelerating one hides
a bubble it is still growing — both are planted. So a positive excess is
evidence that gas is still being made; a negative one is not evidence that it
is not.

**So this is a bracket and not a further repair.** `vmax_corrected` calls the
whole tail chemistry and `vmax_terminal` calls the whole of it gas, and the
answer to "are these curves too pathological to correct" is that the question
is already settled elsewhere: the two orders move by **+0.016** in substrate
and **−0.005** in peroxide across the whole bracket, against their own standard
errors of 0.047 and 0.071, so **nothing this document reports lives inside it**.
Exps 140.4 and 142.4 — the two the bracket bites hardest on — sit at loads of
1.24 and 2.96 and already carry no rate under any repair. **No curve is
excluded for ending mid-bubble**, and none needs to be.

One thing the pair does say about the older statistic: **`bubble_load` is not a
complete severity axis.** It divides the absorbance that *left* by the curve's
rise, so it is blind by construction to gas that never left. **4 curves the
load calls unbiased** carry a bracket wider than 5% of their own rate with a
positive tail excess to support it, the two widest being 145.4 and 140.5 at
loads of 0.34 and 0.28, whose `v_max` falls to 0.73× and 0.78×. They are not
excluded either, and nothing rests on them individually; but a curve read one
at a time should be read with both numbers.

### The rate the model fits is the peroxide decomposing

The production rate is read off the timing and size of the detachments alone —
**the fit never sees a concentration** — so what it correlates with is a
prediction the gas argument makes rather than a parameter it was given.

It is **first order in peroxide, +1.389 ± 0.251** over the 44 live curves that
carry one, 5.5σ from zero and within 2σ of exactly first order. That is the
catalysed decomposition of H₂O₂, measured a second and independent way: the
ladder above shows detachments get *more common* with peroxide, and this shows
the gas is made *faster*, in proportion.

The substrate carries −0.307 ± 0.103 with one offset per experiment — a weak
negative, 3.0σ, of the sign a substrate competing for the same catalyst would
give. What matters for the diagnosis is that it is nothing like the substrate
order of a rate: **the alcohol is not what is turning into gas.**
`scope.gas_rate_drivers`.

**What the model still cannot do.** The last bubble on a curve never detaches,
so nothing reveals how much gas was sitting in the beam when the run ended, and
the reconstruction deliberately leaves it there; the rate is a lower bound
whenever a bubble leaves only partly, which is why the recovery sits at
1.01–1.07 rather than 1.00 under that planting; and the detector cannot see a
bubble that grows and leaves inside one reading. The load
below is a ceiling on use rather than a threshold to tune.

### Which beam the bubble is in, and why it raises the absorbance

The trace says something accumulates and then releases. Two questions decide
what: **which cuvette** the gas is in, and **which of two opposite optical
effects** a bubble has. Every run here is double-beam, and the reference omits
only the *enzyme* — it holds the same peroxide, substrate and buffer — so
`ΔA = A_sample − A_reference` and a bubble's sign flips with the beam.

A bubble in the beam does two opposite things. It **displaces absorbing
solution**, which lowers that beam's absorbance, and it **scatters and refracts
light out of the collection aperture**, which raises it. Crossing that against
the two cuvettes gives four possibilities, and only two of them produce the
slow-rise/sudden-fall shape the curves actually have:

| bubble grows in | if displacement wins | if scattering wins |
|---|---|---|
| sample | ΔA falls slowly, jumps up on release | **ΔA rises slowly, drops on release** |
| reference | ΔA rises slowly, drops on release | ΔA falls slowly, jumps up on release |

The two survivors are *scattering in the sample* and *displacement in the
reference*. Three things pick the first.

**Displacement cannot make a step that large.** It is capped by the absorbance
the displaced solution actually carried, and the largest step in exp 141
cuvette 3 is **0.0303 AU**. The exports record no wavelength, but ε = 1.23
mM⁻¹cm⁻¹ places these runs in the benzaldehyde n→π\* region near 285–300 nm,
where 73 mM H₂O₂ contributes only a few hundredths of an AU across the entire
1 cm path. A bubble would have to swallow a large fraction of the cuvette —
millimetres — and a bubble that size does not sit still for ten minutes.
Scattering has no such ceiling.

**The gas is made where the catalyst is, and this is the ± enzyme control the
block otherwise lacks.** Both cuvettes hold the same peroxide, substrate,
buffer and pH and differ only in the enzyme, so peroxide standing in solution
would bubble both beams alike and the shape would be symmetric. It is not — see
the asymmetry below. Every run is therefore its own matched control, which is
worth more than the archive's four enzyme-free runs: those are short, and
against their own buffers' catalysed rates the whole enzyme-free set is worth
about **one expected event**, so their zero decides nothing
(`scope.gas_enzyme_control`). This argument needs no wavelength and no optics.

**And the large steps are lopsided the right way.** `scope.bubble_step_asymmetry`
over the block's **28827** steps finds **23 rises** beyond 20σ against
**122 falls**, a ratio of **5.3**, with the largest fall (−260σ) 4.6× the
largest rise (+57σ). Sample-beam events dominate, and the minority of sudden
rises is what a stray reference-beam bubble — or one sliding *into* the sample
beam — would contribute.

**What that does and does not establish.** It puts the gas in the cuvette that
holds the enzyme, which is the catalyst-dependence S4 needs. It does **not**
identify the catalyst as the chemzyme's *ketone*: anything arriving with the
enzyme stock would do the same job — the cyclodextrin scaffold, or a trace
transition metal, which is catalase-like and would show exactly this pH
profile. The archive holds no run with the cyclodextrin alone, the ketone
alone, or a chelator added, so the attribution to the ketone is a hypothesis
consistent with the data and not a result. `MECHANISM.md` S4,
`COMPUTATIONAL.md` C10.

### The gas budget makes this unavoidable, not lucky

The objection to reading a few percent of a side reaction as visible bubbles is
that it sounds too small. It is not, by a wide margin
(`solution_chemistry.oxygen_budget`). 2 H₂O₂ → 2 H₂O + O₂, so the ceiling is
half the peroxide, and water at 25 °C holds only **1.25 mM** O₂:

| | O₂ it can make | decomposed before the solution saturates |
|---|---|---|
| the block's top peroxide, 163.2 mM | **81.6 mM** | **1.5%** |
| where the turnover control sits, 73.4 mM | **36.7 mM** | **3.4%** |
| the bottom of the ladder, 2.4 mM | **1.2 mM** | **102%** |

Past saturation every further mM of O₂ is **24.45 µL of gas per mL of
solution**. One percent of net decomposition beyond that point is millimetres
of bubble in a 1 cm cuvette. Given a ketone known to decompose peroxide to O₂
(`MECHANISM.md` refs 34–35), a curve at 73 mM and high pH that did *not* chop
would be the thing needing explanation.

The bottom row is the other half of the ladder in §5's first table, and it is
emphatic: at 2.4 mM, decomposing **all** of the peroxide would still leave the
solution able to dissolve every molecule of O₂ it made. A bubble there is not
merely unlikely, it is unavailable — which is why **no curve below 5 mM carries
a detachment at all**.

### The gas is not a property of this block — the archive says so

Every claim above was measured inside exps 135–151, and the prediction S4 makes
reaches further: a catalyst that decomposes the peroxide involves **no alcohol
at all**, so the gas should appear wherever the peroxide is high and the pH
sufficient, whatever the substrate. `scope.gas_curves` runs the same detachment
test over all **402 curves of 88 experiments**, keeping the ones that did not
bubble, because those are the control.

**It appears with both substrates** (`scope.gas_substrate_control`, catalysed,
≥ 40 mM, pH ≥ 8):

| substrate | curves | detaching | experiments | buffers | events/h |
|---|---|---|---|---|---|
| 4OMe-BnOH | 58 | **20** | 16 | 3 | 0.72 |
| BnOH | 68 | **22** | 17 | 3 | 1.05 |

So the chop is not benzyl alcohol's business and not this block's.

**And peroxide is not the trigger this block made it look like.** The archive's
median [H₂O₂] is **82.5 mM** and **278 of 402** curves sit above 80 — if
peroxide alone did it, most of the archive would chop, and only **28 of those
278** do. The discriminator is pH, and it holds *inside each buffer*
(`scope.gas_survey`, catalysed, ≥ 40 mM, detachment events per hour of run):

| buffer | pH ≤ 7.5 | 7.5–8.5 | > 8.5 |
|---|---|---|---|
| phosphate | **0 in 250.8 h** (80 curves, 20 exps) | 0 in 33.1 h | 0.251 |
| pyrophosphate | 0 in 18.9 h | 1.357 | 3.548 |
| boric | — | 0 in 2.9 h | 0.551 |

**Zero detachments in 270 hours of catalysed, high-peroxide running below
pH 7.5**, over 23 experiments, is the strongest single statement this archive
makes about the gas,
and it is what a reaction consuming HOO⁻ predicts. It is also what leaves the
turnover control above with nothing to add.

### Is it O₂? Nothing in this archive measures the gas

There is **no headspace analysis, no manometry and no oxygen electrode**
anywhere in this project. The identification is an inference, and it is worth
stating what carries it.

*For O₂:* the budget above; the fitted production rate is **first order in
peroxide, +1.389 ± 0.251**, from a fit that never saw a concentration
(`scope.gas_rate_drivers`); and the same fit puts it **negative in substrate**,
−0.307 ± 0.103.

*Against CO₂ — and this is the stronger half:* CO₂ generated in solution at
these pH values **is not volatile**. H₂CO₃/HCO₃⁻ has pKa₁ 6.35, so above pH 8
essentially all dissolved inorganic carbon is bicarbonate and stays in
solution. And the gas becomes *more* common as pH rises, which is precisely
backwards for carbonate, which is driven out of solution by acid. That also
disposes of the one CO₂ route the composition cannot exclude — oxidative
degradation of the cyclodextrin itself, which would be enzyme-, peroxide- and
pH-dependent like the real thing, but would still be captured as bicarbonate at
pH 9.

So the defensible statement is: **a non-condensable gas, made in the
enzyme-containing cuvette, at a rate first order in peroxide and rising steeply
with pH.** O₂ is by a distance the best candidate and CO₂ is close to excluded;
neither has been measured.

### One curve, read detachment by detachment

`scope.bubble_record(141, 3)` — 0.865 mM substrate, 73.4 mM H₂O₂, pH 9.15:

| detached at | drop | σ | grew for | rose over that window | level after |
|---|---|---|---|---|---|
| 2460 s | 0.0023 | −8.4 | 2460 s | +0.0490 | +0.0472 |
| 3120 s | 0.0214 | −77.7 | 600 s | +0.0300 | +0.0558 |
| 3780 s | 0.0153 | −55.6 | 600 s | +0.0242 | +0.0647 |
| 4140 s | **0.0303** | −109.9 | **300 s** | +0.0086 | **+0.0430** |

It sheds **0.0694 AU** in total against a net rise of 0.0551 — load 1.26, past
the ceiling.

**And it is not one bubble.** A single site growing steadily would shed in
proportion to the time it had to grow; here the *largest* drop follows the
*shortest* window. The level after each release says the same from the other
end: it climbs +0.0472, +0.0558, +0.0647 — consistent with product
accumulating underneath — and then falls to +0.0430, below all three. A
monotone chemistry beneath a single bubble cannot do that. Something detached
at 4140 s that had been growing since before the release at 3780 s, so **at
least two bubbles were on different sites at once** and the trace is their sum.

That is the honest limit of the repair. `debubble` charges each drop to growth
since the previous drop, which is right on average — hence the unbiased
recovery above — and wrong curve by curve. It is why the correction is quoted
with the monotone bound beside it rather than on its own.

### What each curve is entitled to

`bubble_load` is the absorbance lost to detachments over the curve's net rise.
`scope.bubble_table`:

| load | live curves | median [H₂O₂], mM | what the curve carries |
|---|---|---|---|
| none | 60 | 25.0 | nothing to correct |
| ≤ 0.25 | 19 | 35.2 | the correction is unbiased |
| 0.25–0.5 | 12 | 36.0 | the correction is unbiased |
| 0.5–1 | 6 | 73.4 | a rate, with about a tenth of systematic |
| **> 1** | **13** | **73.4** | **no rate this package can measure** |

The thirteen are not scattered: **all four substrate rungs of exp 135** (163 mM
peroxide, the archive's highest, loads 4.0 to 8.7), and inner rungs of 138, 140,
141, 142 and 150. They are **flagged and not excluded** — every one stays in the
frame, in the live counts and on `progress_curves.html`. Curve shape is never a
defect in this project, and these are not excluded for their shape; what they
lack is a measurable rate.

Beside the correction sits the assumption-free bracket. Gas only adds absorbance
and product only accumulates, so the greatest non-decreasing function under the
readings is an upper bound on the chemistry with no bubble model in it at all
(`curve_metrics.monotone_bound`). It costs a median **0.6%** of `vmax` across the
block and **83%** on the worst curve — that spread *is* the contamination.
Quote the corrected rate with the gap to the bound as its systematic, the way
`product_fate/` quotes the sink's activation energy.

### The substrate order does not rest on the gas, and the peroxide order moves

This is the check that had to pass before §2 could stand, because a
substrate-blind, peroxide-driven additive artefact is precisely what would
manufacture a flat substrate order. `scope.bubble_sensitivity`:

| `vmax` from | n | order in [S] | order in [H₂O₂] |
|---|---|---|---|
| the readings | 110 | +0.091 ± 0.052 | +0.794 ± 0.077 |
| the reconstruction | 110 | +0.093 ± 0.047 | **+0.700 ± 0.070** |
| the monotone bound | 110 | +0.141 ± 0.049 | +0.770 ± 0.072 |
| the reconstruction, terminal bubble charged to gas | 110 | +0.110 ± 0.051 | +0.696 ± 0.077 |
| readings, load ≤ 1 only | 97 | +0.139 ± 0.059 | +0.767 ± 0.086 |

**The substrate order does not move** under any repair — by less than the two
estimates' errors combined every time — and under the monotone bound it moves
*away* from zero rather than toward it, which is the direction that matters
since the artefact could only have flattened it. §2 stands.

**The peroxide order does move.** The reconstruction takes it from +0.794 to
+0.700, and on the strong runs alone from +0.871 to +0.765 — about 0.10 either
way, which is 0.9σ and 1.2σ of the two estimates' errors combined and the
largest shift any repair here produces. That is not the repair failing, it is
the repair working: the gas is *made from peroxide*, so an uncorrected artefact
has to inflate the apparent peroxide order, and taking the gas out has to bring
it down. Neither reading changes any conclusion in this document — but it is
the one number here that the choice of repair visibly touches, and the earlier
claim that *no* order moved was made with a repair too weak to move it.

Every live curve still carries a rate after the repair.

### The acceleration is a high-pH effect, not a long-run one

| band | live curves | accelerating | share | median late ÷ early |
|---|---|---|---|---|
| pH >= 9 | 42 | 27 | 64% | +1.17 |
| pH < 9 | 68 | 23 | 34% | +0.49 |

The eight-hour runs are the *slow* ones — exps 138, 146 and 148–151 are long
because they are slow, and they decelerate — so "the acceleration builds up
over a long run" reads the design backwards. What it tracks is pH.

### The block will not carry a *landmark* induction analysis

`induction.signal_control` regresses the landmark on each curve's own
signal-to-noise. On the catalysed 4OMe block it passes; here it returns
**+0.619 ± 0.228**, two and a half standard errors from zero. On this block the landmark is partly measuring the
spectrophotometer, so the induction results of `induction/` do not transfer to
it and nothing here is read off one.

The curve forms say the same from the other side. The block splits 56 two-phase
against 54 one-phase, so `tau_fast` is τ₁ of one fit on half the curves and τ
of another on the rest; and it splits 46 lag-first against 64 burst-first, so
an induction time averaged over it would average two different things.

### But the +1 rule does not need a landmark

The constraint above is about a *statistic*, not about the block. A catalyst
drawn into its active form by a species X held in excess gives

    E + X ⇌ E*     1/τ = k_f[X] + k_r,     [E*]/E₀ = K[X]/(1 + K[X])

so `d ln v/d ln[X] = 1/(1 + K[X])` and `d ln τ/d ln[X] = −K[X]/(1 + K[X])`, and
their **difference is 1 for every K and every [X]** — a parameter-free
prediction, and an *exact* one rather than a local slope. It says nothing about
which clock τ is. Every use of it in this project has nonetheless gone through
`t_ind`, and `t_ind` is the rolling window this block cannot use.

`tau` and `tau_slow` come from the **progress fit**. They carry no window, so
neither the 9.6× run-length span nor the failed signal control touches them, and
the identity applies to them unchanged. `induction.joint_clocks` asks it through
each clock in turn, as one regression rather than two differenced by hand — the
two orders are fitted to the same curves on the same design, so their errors are
correlated and a hand-difference would quote an error that is not theirs.

**Both sides of the ratio are read off the rebuilt curves** (§5), and on this
axis that is not housekeeping. The O₂ is *made from peroxide*, so leaving it in
inflates the rate's peroxide order and shortens the apparent clock, and both
push `d ln v − d ln τ` towards the +1 being tested. Asked of the readings the
`tau_slow` row sits 0.3σ from +1; asked of the rebuilt curves, 1.4σ. Taking the
gas out costs no resolution either — it buys some, because the artefact was what
some of those fits could not pin: 62 to 65 curves for `tau` and 25 to 32 for
`tau_slow`, with smaller errors throughout.

| clock | window | curves | order in [H₂O₂] | from +1 | control, [S] |
|---|---|---|---|---|---|
| `t_ind` | a tenth of the run | 110 | +0.304 ± 0.188 | 3.7σ | 6.4σ |
| `tau` | none, from the fit | 65 | +0.707 ± 0.158 | 1.9σ | 8.5σ |
| `tau_slow` | none, from the fit | 32 | **+0.669 ± 0.242** | **1.4σ** | 6.0σ |

*(all 110 live curves, rebuilt; the fitted clocks are gated on being resolved.
The `t_ind` row is **not** gas-corrected and cannot be — it is the comparator
this block rejects on other grounds, not a like-for-like line.)*

**The two routes disagree, and that is the result here.** Through the landmark
the peroxide axis falls 3.7σ short of the +1; through the clocks that carry no
window it falls 1.9σ and 1.4σ short. Over the strong runs alone `tau` gives
+0.768 ± 0.180 on 40 curves and `tau_slow` +0.850 ± 0.283 on 20.

**The substrate axis is the control, and it is the reason any of this counts.**
The +1 belongs to the species that draws the catalyst into its active form; the
alcohol does not, and the clock carries no substrate order, so that axis must
*miss*. It does, in every cut — by 4.2σ to 8.5σ. A reading of this table where
both axes met +1 would be a regression that had stopped discriminating.

**No conclusion is drawn from it, and the reason is the count.** `tau_slow` is
resolved on 32 of 110 live curves and 20 of the 77 strong ones, and across
those cuts the peroxide-axis estimate moves from **+0.67 to +0.85** — short of
+1 everywhere, but nowhere by enough to reject it. What has
changed is not the answer but the question's availability: the landmark's
failure closed a statistic, not the block, and the route that stays open is the
one the block's design actually supports. An enzyme-free arm, or enough
resolved slow constants to halve that error, would settle it.

### How big is the early rise, and is it one turnover or many?

Many of the block's curves rise to a maximum and then fall — exp 150 cuvette 6
is the plainest: its fit climbs **0.0129 AU** above its start, peaks at
6900 s, and then declines steadily for six hours. The chemically interesting
quantity there is **how much product the catalyst made before it stopped**, in
units of the catalyst itself.

`curve_metrics.burst_amplitude` measures it off the **fitted** curve rather
than the readings, and the reason is not noise. When the two exponentials are
nearly degenerate the linear solve trades enormous opposite amplitudes between
them — exp 135 cuvette 3 returns `B_fast` = -241 against `B_slow` = +303 on a
curve that moves 0.06 AU — so neither `B_fast` nor their sum is the burst. The
trade leaves the fitted *curve* alone and moves only the split, so the value of
the curve at a time inside the window is what survives it.

**Only 17 of the 110 live curves finish their rise inside their run.** The rest
peak on the last reading, and their amplitude is a lower bound; `burst_bounded`
carries the distinction, and it matters *between* runs and not within them —
inside a run every cuvette shares the length to within one 60 s reading, while
between runs it spans 9.6x.

| band | bounded curves | median turnovers | range |
|---|---|---|---|
| pH < 9 | 11 | 1.41 | 0.28-3.62 |
| pH >= 9 | 6 | 3.31 | 1.98-4.22 |

So at the bottom of the pH ladder the catalyst gets **about one turnover and
stops** — exp 150 cuvette 6 is 0.50 — and at the top it gets three or four.
That is what a catalyst whose regeneration needs HOO⁻ should do. It is also
where the caveat bites hardest: the sub-stoichiometric curves are in exps 149,
150 and 151, which are the runs `concentration_agreement` rejects, and the one
sub-stoichiometric burst in a strong run is exp 135 cuvette 4 at 0.58.

### And the rise scales with the catalyst, not with the substrate

Within runs, where both concentration axes move, the rise is nearly flat in
substrate and strong in peroxide:

| set | order in [S] | order in [H₂O₂] | n | R² |
|---|---|---|---|---|
| all live curves | +0.201 ± 0.041 | +0.757 ± 0.062 | 110 | 0.795 |
| strong runs | +0.136 ± 0.045 | +0.793 ± 0.063 | 77 | 0.805 |

The catalyst cannot be asked that way: `[enz]` is constant across every cuvette
of a run and moves only between runs, so the per-experiment offsets absorb it
and `scope.orders` refuses it. The block offers exactly one comparison instead
— **exps 140 and 141**, which share a composition, step the loading 0.034
against 0.014 mM, and sit 0.07 pH units apart, the smallest gap in the block.
Both fits are read at **3780 s**, the largest window the two runs share, because
their own maxima fall at different times; and the residual pH gap is corrected
with §3's own pH order, worth 1.065x against the 2.43x enzyme step.

| corrected rise ratio | expected if first order in [enz] | σ from first order | σ from no dependence |
|---|---|---|---|
| **1.99 ×/÷ 1.21** | 2.43 | 1.0 | **3.6** |

**The early rise is proportional to the catalyst.** No dependence on it is
excluded at 3.6σ; first order in it is not distinguishable at 1.0σ. Sweeping
the window gives an apparent order of +0.77 to +1.35 over 945 to 3780 s — a
systematic about the size of the statistical error, and every window in the
sweep says the same two things.

| window, s | ratio | apparent order in [enz] | cuvettes |
|---|---|---|---|
| 945 | 2.58 | +1.07 | 5 |
| 1890 | 3.30 | +1.35 | 7 |
| 2835 | 2.47 | +1.02 | 7 |
| 3780 | 1.99 | +0.77 | 7 |

It rests on **one pair of runs**, which is all the block has, and it should be
read as such. What it does say is that the thing accumulating early is counted
by the catalyst rather than by the substrate — the substrate order of the same
quantity is +0.14 ± 0.05 over a fifty-fold ladder.

### It decelerates on a clock, not on its product

| coefficient | value | what it means |
|---|---|---|
| elapsed time | **-0.697 ± 0.136** | the rate falls with the clock |
| product made | +0.283 ± 0.129 | and not with what it has made |

This is the *opposite* of what the same test says about the 4OMe archive, where
the rate declines linearly in the product and the clock coefficient is the one
that vanishes (`product_fate/`). One progress curve cannot tell time from
product — inside a curve the product only grows with time — and the separation
exists only across curves, which is why it takes a block to ask.

## 6. What the block cannot do

- **No enzyme-free curve.** All 119 carry catalyst, and the pyrophosphate cell
  holds no blank at all, so nothing here can be staged against a background
  measured in its own buffer. `scope.PAIRED_CONTROLS` is the nearest thing and
  it is boric and phosphate.
- **No interior point.** The L identifies two orders and no interaction between
  them, which is what §4 runs into.
- **Enzyme is confounded with pH** between the sub-series and constant within
  them, which is why §3 measures inside a ladder rather than across the block.
- **No windowed statistic travels across it**, at 9.6× in run length — which
  closes the induction landmark and, with it, `induction/`'s results. It does
  not close the +1 rule, which §5 asks through the progress fit instead.
- **Thirteen curves carry no measurable rate**, their O₂ detachments having
  moved more absorbance than the reaction did (§5). They are flagged, drawn
  and counted; they are not excluded, and nothing here is read off one.
- **No mechanism fit has ever been run on it.** `data/fits/` holds one saved
  fit, on BnOH / 25 °C / *phosphate*, and it shares no experiment with this
  block: stage 1 on exps 3, 6, 67, 69, 70 and stage 2 on exps 68, 71, 73, 74,
  83. The block with the strongest design is the one the fitter has never
  seen — and §2's substrate order and §4's gap are both things a fit would have
  to reproduce.
