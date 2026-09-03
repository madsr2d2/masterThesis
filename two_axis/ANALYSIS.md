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

**Figures**: [`index.html`](index.html) is the presentation — six figures, A to
F, one per claim below. [`progress_curves.html`](progress_curves.html) carries
all 119 cuvettes with the form each earned, which is the audit surface for
every number here.

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

### The acceleration is a high-pH effect, not a long-run one

| band | live curves | accelerating | share | median late ÷ early |
|---|---|---|---|---|
| pH >= 9 | 42 | 27 | 64% | +1.17 |
| pH < 9 | 68 | 23 | 34% | +0.49 |

The eight-hour runs are the *slow* ones — exps 138, 146 and 148–151 are long
because they are slow, and they decelerate — so "the acceleration builds up
over a long run" reads the design backwards. What it tracks is pH.

### The block will not carry an induction analysis

`induction.signal_control` regresses the landmark on each curve's own
signal-to-noise. On the catalysed 4OMe block it passes; here it returns
**+0.619 ± 0.228**, two and a half standard errors from zero. On this block the landmark is partly measuring the
spectrophotometer, so the induction results of `induction/` do not transfer to
it and nothing here is read off one.

The curve forms say the same from the other side. The block splits 56 two-phase
against 54 one-phase, so `tau_fast` is τ₁ of one fit on half the curves and τ
of another on the rest; and it splits 46 lag-first against 64 burst-first, so
an induction time averaged over it would average two different things.

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
- **No windowed statistic travels across it**, at 9.6x in run length.
- **No mechanism fit has ever been run on it.** `data/fits/` holds one saved
  fit, on BnOH / 25 °C / *phosphate*, and it shares no experiment with this
  block: stage 1 on exps 3, 6, 67, 69, 70 and stage 2 on exps 68, 71, 73, 74,
  83. The block with the strongest design is the one the fitter has never
  seen — and §2's substrate order and §4's gap are both things a fit would have
  to reproduce.
