# What happens before the 4OMe curves start

Every catalysed 4-methoxybenzyl alcohol run from 15 to 30 °C begins slowly and
takes **thousands of seconds** to reach the rate `../temperature_series/` puts
on an Arrhenius plot. That document names the induction and times it. This one
asks what it *is*.

The answer the archive supports: **the catalyst is not in its active form when
the run starts, and it converts into it on its own clock.** The induction needs
the catalyst, it does not wait for product, its size is a fraction rather than a
concentration, and its barrier is 95 kJ/mol — a covalent step, not a physical
one. What the archive cannot say is *which* covalent step, and in particular it
cannot test the obvious candidate, the ketone binding the peroxide: the one
experiment that moves `[H2O2]` on this substrate is confounded with
signal-to-noise, and §4a shows the confounding rather than hiding it.

The evidence is drawn from the **whole catalysed 4OMe archive** — 147 live
curves in 38 experiments — not from the temperature series alone. The
temperature series is 24 curves and its substrate ladder moves the rate by only
a factor of two; the discrimination in §3 needs that lever repeated across
enough runs to have an error bar worth quoting.

    data/induction.py               the landmark, the regressions, the controls
    python data/induction.py        the whole argument, printed
    python data/test_induction.py   the statistic, its floor, the numbers below
    python induction/build_figures.py
    python induction/check_numbers.py

**Figures**: [`index.html`](index.html) is the presentation — seven figures, A
to G, one per claim below. It is rebuilt by `build_figures.py`, which computes
nothing.

Related: [`../temperature_series/`](../temperature_series/ANALYSIS.md) for the
activation parameters, [`../product_fate/`](../product_fate/ANALYSIS.md) for the
other end of the same curves, and `COMPUTATIONAL.md` C7 and C8 for the two
calculations that would finish this.

## 1. Four candidates, and what each one waits for

| | ends when | needs the catalyst |
|---|---|---|
| **seeding** — the catalyst-free Cannizzaro loop (`MECHANISM.md` 1–3) needs product before it can run | enough **A** has been made | no |
| **scavenger** — something in the reagents consumes the oxidant until it is used up | enough **turnover** has happened | no |
| **activation** — the catalyst converts into its active form | enough **time** has passed | yes |
| **schedule** — the induction is an artefact of when recording started and stopped | — | no |

The first two end on the **product**, the third on a **clock**. That is the same
discrimination `product_fate` makes for the fall at the other end of these
curves, and it has to be made the same way: within one curve the product only
grows with time, so "the rate rose after 2000 s" and "the rate rose after
0.02 AU" are the same sentence. Only a design that moves the rate while holding
the schedule can separate them — and the substrate ladder inside every 4OMe run
does exactly that, a factor of **4.0 in `[S]`** inside the median run — 28 of
the 38 experiments carry a ladder — at fixed peroxide, catalyst, buffer, pH,
temperature and run length.

### The statistic, and why a landmark is allowed here

`slowdown` withdrew a landmark — when does the rate fall to three quarters of
its peak — because a curve whose rate never falls that far has none, and
dropping those curves biases the answer. **The rise has no such problem, for a
structural reason rather than a lucky one**: the landmark is a fraction of the
curve's *own* maximum, and every curve reaches its own maximum.
`induction.induction_point` is therefore defined for all 402 curves in the
archive, and a curve that is fastest in its first window returns **zero**, which
is a measurement and not a gap. Nothing here is censored.

It is `curve_metrics.lag_time` without that function's acceleration gate, and
the difference matters: the gate is passed more often by fast curves, so gating
would put exactly the correlation this section is testing for into the data by
hand.

What it costs instead is a **window**. The rolling slope is read through a tenth
of the run, because a slow curve needs a wide window to have a slope at all — at
15 °C the whole run is 0.04 AU — and a fixed window in seconds turns the cold
curves into noise: a 300 s window puts the 15 °C induction at **529 s** instead of
**4289 s** and drags the block's activation energy from 86.5 ± 10.7 to
**24.0 ± 19.1 kJ/mol**. Widen it to 900 s and it recovers, to 80.3 ± 11.9 —
what fails is a window too narrow for the coldest curves, not the idea of an
absolute window. A
window that is a fraction of the run is safe **within** a run, where every
cuvette shares the schedule, and unsafe **between** runs, where it is not:
across the 4OMe archive at 25 °C the induction time regresses on run length with
an exponent of **+0.437 ± 0.181**, and on pH — once run length is in the model —
with **−0.004 ± 0.123**, so the schedule carries what pH does not.

**So every concentration order below is measured within experiments.** The one
between-run quantity used anywhere here is the temperature dependence, and that
is taken from `arrhenius`'s fitted `inverse_tau`, which is not windowed.

## 2. The induction needs the catalyst

`differential` is the structural classification — what the reference channel
omitted — not the filename. A catalysed run's reference omits the *enzyme*, so
the curve is the catalytic increment; a background run's omits the *H₂O₂*, so
the curve is the raw reaction.

| | curves | median depth | depth > 0.25 | accelerating | longest run |
|---|---|---|---|---|---|
| **4OMe catalysed** | 151 | **0.289** | 83 | 91 | 17940 s |
| **4OMe, no enzyme** | 49 | **0.000** | **0** | **0** | 17934 s |

`depth` is the fraction of the peak rate that is missing at the start. **Not one
of the 49 enzyme-free curves has an induction** — their largest rolling slope is
the first one — and none of them accelerates past 3σ; the largest acceleration
z-score in the whole enzyme-free block is **2.25**, against 49.6 in the
catalysed one.

The obvious objection is that the enzyme-free runs were too short to show
anything. They were not: exp 28 ran **17934 s** at 40 °C, the longest run in the
archive, and is flat from its first window.

The second objection is that the two channels are not comparable because they
sit at different compositions. Matched, they are:

| T | channel | curves | experiments | median depth | accelerating | longest |
|---|---|---|---|---|---|---|
| 25 °C | enzyme-free | 8 | 40, 52 | 0.039 | 0 | 3565 s |
| 25 °C | **catalysed** | 8 | 14, 20 | **0.291** | **8** | 17934 s |
| 40 °C | enzyme-free | 37 | 23–30, 38, 39 | 0.000 | 0 | 17934 s |
| 40 °C | **catalysed** | 12 | 16, 32, 34 | **0.231** | **9** | 6566 s |

Same substrate, same buffer, same 82.5 mM H₂O₂, same pH to within 0.03–0.5, same
temperature. **Seeding and a reagent-borne scavenger both run in a cuvette with
no catalyst in it**, and neither shows.

## 3. It runs on a clock, not on product

Two routes to the same coefficient, because each has a weakness the other does
not.

**Route one — regress the induction time on the rate**, with one offset per
experiment so that temperature, pH, buffer, enzyme, cell, day *and the run
length* are all absorbed and only the in-run substrate ladder carries the fit.

> A clock predicts **0**: the induction takes the time it takes.
> Product control predicts **−1**: a curve twice as fast reaches the threshold
> in half the time.

| block | rate used | slope | curves |
|---|---|---|---|
| **4OMe catalysed** | `v_peak` | **−0.025 ± 0.109** | 147 |
| 4OMe catalysed | `peak_rate` | +0.138 ± 0.122 | 147 |
| 4OMe catalysed | `vmax` | +0.148 ± 0.123 | 147 |
| the temperature series alone | `v_peak` | −0.132 ± 0.442 | 24 |

**Product control is excluded at 9σ** on the whole catalysed block. The
temperature series on its own cannot do it — 24 curves and a two-fold lever give
±0.44 — which is why this section is not written from that block.

The answer does not come from the floor that the zeros are placed at. Sweeping
it from 1 s to 300 s moves the coefficient from +0.181 ± 0.242 to +0.094 ± 0.083,
never near −1.

**Route two — divide the two substrate orders.** Route one puts a rate measured
off the same curve on the right-hand side, and noise in a regressor attenuates
its coefficient *towards zero*, which is the clock's answer. That bias cannot be
argued away, so the same question is asked again with the regressor replaced by
the composition the operator set:

| block | order of `t_ind` in [S] | order of `v_peak` in [S] | ratio |
|---|---|---|---|
| **4OMe catalysed** | −0.121 ± 0.148 | +0.471 ± 0.084 | **−0.26 ± 0.32** |
| the temperature series | −0.008 ± 0.220 | +0.468 ± 0.039 | −0.02 ± 0.47 |

Less sharp, and it has no bias towards the answer. It agrees.

**And the amplitude is a fraction, not a concentration.** `depth`'s substrate
order is **−0.114 ± 0.169** over the catalysed block and −0.211 ± 0.387 over the
temperature series: how far below its eventual speed the run starts does not
depend on how fast that eventual speed is. A catalyst that begins partly in an
inactive form gives a fraction; an intermediate accumulating to a
substrate-dependent steady state does not.

**This is the exact mirror of the fall.** In the same block `product_fate` finds
the deceleration carried by the product (−0.598 ± 0.053) and not by the clock
(−0.124 ± 0.079). The rise and the fall of one curve are driven by different
things, and the intuition that would assign them the other way round — a rise
that waits for product, a fall that runs down a clock — is wrong at both ends.

## 4. It is not the schedule, and it is not physical

**Not the schedule.** Exps 19 and 14 were both recorded for **17934 s**, and
their induction time constants are **6489 s and 945 s** — a factor of 6.9 at the
same run length, 15 °C against 25 °C. Whatever the schedule does to τ it does
equally to both.

**Not the instrument warming up.** τ at 15 °C is a hundred minutes. A 1 cm
cuvette in a thermostatted holder equilibrates in minutes, and the effect would
not scale with the catalyst.

**Not dissolution or aggregation of the cyclodextrin.** The induction's
activation energy is **95.0 ± 15.7 kJ/mol** (`arrhenius`, `inverse_tau`, 16
curves at four temperatures). Diffusion-limited and dissolution processes in
water run at 15–25 kJ/mol. A barrier of 95 kJ/mol is a covalent step.

**Not signal-to-noise.** The landmark is the first crossing of half the *largest*
rolling slope, so on a curve with no signal the largest rolling slope is a noise
excursion and the landmark could be measuring the spectrophotometer. Regressed
on each curve's own `net/noise`, with the same per-experiment offsets:

| block | slope on signal-to-noise | curves |
|---|---|---|
| **4OMe catalysed** | **+0.003 ± 0.149** | 147 |
| BnOH, exps 135–151 | +0.619 ± 0.228 | 110 |
| 4OMe peroxide pair, exps 127–131 | +0.702 ± 0.241 | 15 |

**The block this analysis rests on passes; the other two fail.** That is why §3
is written from the catalysed 4OMe archive and why the next two paragraphs
report negatives rather than results.

### 4a. The peroxide question is open, and the archive cannot close it

If what the induction times is the catalyst binding the oxidant —
`K + H₂O₂ ⇌ KP`, `MECHANISM.md` step 4 — then the relaxation rate is
`1/τ = k_f[H₂O₂] + k_r`, the induction **time** has an order in `[H₂O₂]` between
0 and −1, and it cannot have any other sign. That is a clean falsifiable
prediction and this archive cannot test it.

Every 4OMe run in the archive sits at 82.5 mM H₂O₂ **except exps 127–131**,
which step 3.879 against 195.882 mM at five pH values, two cuvettes per level,
19 curves of which 15 are live. There the induction time's order in `[H₂O₂]` is
**+0.302 ± 0.092** — the wrong sign, three standard errors from zero — while the
rate's is +0.804 ± 0.168. But that same block's induction time regresses on
signal-to-noise at +0.702 ± 0.241, and `[H₂O₂]` is what sets the signal, so the
two cannot be separated on 15 curves. **The sign is suggestive and it is not
evidence.**

The same objection sinks the BnOH comparison. Exps 135–151 move `[H₂O₂]`
thirty-fold inside every run, which is the design this question wants, and the
induction time there carries an order of +0.447 ± 0.155 in it — but so does
signal-to-noise, at +0.619 ± 0.228. **This archive holds no clean peroxide
lever on an induction period.** One pre-incubation run would settle it and there
will not be one, so it goes to `COMPUTATIONAL.md` instead.

## 5. What the activation parameters say

Both rows are `arrhenius.activation_parameters`; neither is refitted here.
`1/τ` is the one parameter in that table whose Eyring entropy rests on no
assumption about the rate law — it is already a first-order relaxation rate in
s⁻¹, needing neither an extinction coefficient nor an enzyme concentration.

| | E<sub>a</sub> kJ/mol | ΔH‡ kJ/mol | ΔS‡ J/mol/K | ΔG‡(298) kJ/mol | curves |
|---|---|---|---|---|---|
| **the induction**, `1/τ` | 95.0 ± 15.7 | 92.6 ± 15.7 | +3.7 ± 53.2 | **91.48 ± 0.62** | 16, 4 T |
| **the turnover**, `v_peak` | 88.8 ± 1.8 | 86.3 ± 1.8 | −57.7 ± 5.9 | **103.47 ± 0.10** | 24, 6 T |

**The free-energy gap is the solid number: −11.99 ± 0.62 kJ/mol.** As a rate
constant at 298 K the induction step is **126× faster** than turnover — which is
what a step that has to finish before turnover can begin is required to be, and
it is the one prediction in this document that was made before the number was
looked at.

**Its decomposition is not solid, and should not be quoted as though it were.**
The point estimates put the whole gap in the entropy: +61.4 J/mol/K, an
associative transition state against a loose one, with the two enthalpies 6.3
apart. But the entropy gap carries **± 53.6** and the enthalpy gap **± 15.8**, so
each is about one standard error from zero. The Eyring numbers are *consistent*
with a unimolecular step and prove nothing on their own.

**What does carry that conclusion is §3**, from data that share no assumption
with the Eyring fit: the induction time has no order in the substrate, its
amplitude has none either, and the same is true of every rung of the ladder in
every run. Whatever the induction is waiting for, it is not waiting for anything
whose concentration this archive varies.

## 6. What this settles, and what it does not

**Settled.**

- The induction is a property of the **catalysed** reaction. It is absent from
  all 49 enzyme-free 4OMe curves at matched composition, including a 5-hour run.
- It ends on a **clock**, not at a product threshold: the coefficient is
  −0.025 ± 0.109 where product control requires −1, on 147 curves, and the
  bias-free route agrees at −0.26 ± 0.32.
- Its **amplitude is a fraction** of the eventual rate and does not scale with
  the substrate.
- It is a **covalent step**: 95 ± 16 kJ/mol is four times too large for
  dissolution, diffusion or thermal equilibration.
- It is **faster than turnover** by 126× in free energy, as it must be.

**Not settled, and the reasons are different in each case.**

- **Which step it is.** "The catalyst becomes active" is a shape, not a
  mechanism. The candidates that survive §3 are all unimolecular in everything
  the cuvette holds: the ketone's gem-diol **hydrate** dehydrating to the free
  ketone, the perhydrate collapsing to the dioxirane, or a conformational change
  of the cyclodextrin. Absorbance at one wavelength cannot choose between them.
- **Whether the peroxide is involved at all.** §4a: the only lever is confounded
  with signal-to-noise. The observed sign is the *wrong* one for the adduct
  being the step timed, which would make the perhydrate a resting state the
  catalyst has to leave rather than the activation itself — an inversion worth
  taking seriously and not worth asserting on 15 curves.
- **Whether the same thing happens with BnOH.** Exps 135–151 have the peroxide
  design and not the signal-to-noise; their induction statistic tracks the noise
  (+0.619 ± 0.228). The substrate comparison `product_fate` could make at the
  end of the curve cannot be made at the start.

**The three things that would finish it.**

1. `COMPUTATIONAL.md` **C7** — the hydration equilibrium and the dehydration
   barrier of the chemzyme's ketone in water, against the barrier for adding
   H₂O₂ to it. If the hydrate dominates and dehydration costs ≈ 95 kJ/mol, the
   induction has a name.
2. `COMPUTATIONAL.md` **C8** — the free-energy profile from the perhydrate to
   the dioxirane, to say whether the peroxide adduct is on the activation path
   or off it.
3. The experiment that will not be run: **pre-incubate the catalyst with H₂O₂,
   then add substrate.** If the induction is catalyst activation it disappears;
   if it is anything that waits for product it does not. One cuvette.
