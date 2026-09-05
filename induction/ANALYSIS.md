# What happens before the 4OMe curves start

Every catalysed 4-methoxybenzyl alcohol run from 15 to 30 °C begins slowly and
takes **thousands of seconds** to reach the rate `../temperature_series/` puts
on an Arrhenius plot. That document names the induction and times it. This one
asks what it *is*.

The answer the archive supports: **the catalyst is not in its active form when
the run starts, and it converts into it on its own clock, over a covalent
barrier, slowed by base and hurried by the substrate.** The induction needs
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

**Figures**: [`index.html`](index.html) is the presentation — fourteen figures, A
to N, one per claim below. [`progress_curves.html`](progress_curves.html)
carries **both channels of the 4OMe block**, 147 catalysed and 49 enzyme-free,
each with the form it earned and its induction landmark — §1's claim is a
contrast, so a page showing only the catalysed half would show only the half
that agrees. Both are rebuilt by `build_figures.py`, which computes nothing.

**§7 is the survey of all seven experimental variables** — temperature, pH,
`[buf]`, buffer salt, `[S]`, substrate and `[H2O2]` — and of which of them this
archive can answer for at all. It needed a lag statistic with no window in it,
because five of the seven move only between runs; `scope.frame`'s `lag_depth`
and `lag_half_s` are that statistic.

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

*Route two has a buffer riding on it, and correcting for that costs it.* In
every 4OMe run the substrate was added by volume and displaced buffer, so
`[buf]` falls 80 → 50 mM as `[S]` rises: the two correlate at **−0.96 in logs**
across all 28 runs with a ladder, at a slope of `d log[buf]/d log[S] = −0.264`.
Route two's "order in [S]" is therefore an order in the *pair*, exactly as
`temperature_series` §3 found for the rate — and §6 now measures the induction's
own buffer order, **−0.433 ± 0.201**, so the correction can be made:

| | order in [S] | distance from the −0.471 a product threshold needs |
|---|---|---|
| as measured | −0.121 ± 0.148 | 2.4σ |
| **with the buffer taken out** | **−0.235 ± 0.158** | **1.5σ** |

**So route two no longer excludes product control.** It is reported that way
rather than dropped, because the correction is real and small-n both — the
buffer order behind it is eight curves — and because the reader is entitled to
see the route that weakens as well as the one that holds.

**Route three — ask it in the units the hypothesis is stated in.** Both routes
above ask whether a faster cuvette's induction is *shorter*. What the product
hypothesis actually claims is that it ends at the same *product*, and
`InductionPoint.made` — the absorbance built up by the landmark — has been
computed on every curve since this module was written and was read by nothing
until 2026-09-05.

> A clock predicts **+1**: the time is fixed, so the product is whatever the
> rate makes in it.
> A threshold predicts **0**: the product is fixed, so the time is whatever it
> takes.

| block | d log(product at the landmark) / d log(rate) | curves |
|---|---|---|
| **4OMe catalysed** | **+0.885 ± 0.125** | 112 of 147 |
| BnOH two-axis | +0.756 ± 0.115 | 81 of 110 |
| the temperature series | +1.908 ± 0.474 | 22 of 24 |

*It is not an independent witness and must not be quoted as one.* Over the
induction the product is about the rate times the time, so this slope is
`1 + `route one's by construction — route one's −0.025 ± 0.109 predicts +0.975.
What it adds is the units, the spread comparison below, and a calibration that
turns out to matter more than either.

***A planted product threshold does not read 0 through this landmark. It reads
+0.4.*** The rolling slope's first centre sits half a window into the run, so
part of the induction is over before the landmark's own clock starts — and how
much depends on τ, which is exactly what a threshold varies. Writing it out for
a planting with `v·τ` held fixed, `made = C(ln 2 − ½e^(−c₀/τ))`, which rises by
a factor of 3.6 across the band. `induction.product_recovery` plants both
hypotheses at four run/τ geometries: a clock reads back at **+1.000** exactly,
and a threshold at **+0.31 to +0.42**.

So the honest exclusion is **3.7σ from a threshold**, not the 7.1σ the nominal
0 would give, and 0.9σ from a clock. A correction that halves the result is
worth more than the result.

*And the spread comparison, which is not a regression at all.* Inside one run
the schedule, pH, temperature, buffer and catalyst are fixed and only the
substrate ladder moves, changing the rate about four-fold. So ask which
candidate constant is actually the more nearly constant across those cuvettes:

| pooled within-run sd | 4OMe catalysed | a planted clock | a planted threshold |
|---|---|---|---|
| log(t<sub>ind</sub>) | **0.667** | 0.000 | 0.44 |
| log(product at it) | 1.091 | 0.454 | 0.19 |

The archive holds the **time** more nearly constant than the product, and the
plantings say that is the clock's signature and not the threshold's. It is a
weaker separation than the raw numbers suggest — the same half-window bias
stops a planted threshold reaching zero here too — so it is a direction, and it
agrees with route one and route three.

*What route three costs.* It is the one statistic in this folder that **drops
rows**: a curve with no landmark has made no product by it, so 35 of the 147
live catalysed 4OMe curves do not appear. That is a selection towards curves
that have an induction — which is the population the question is about, but it
is a selection, and the count is reported rather than absorbed.

**Route one is not reached by any of this.** Its regressor is the curve's own
measured rate rather than a concentration: it asks whether a faster cuvette's
induction is shorter, and that question is well posed whatever is making the
cuvette faster — buffer included. It is the load-bearing result of this section
and it stands at −0.025 ± 0.109, nine standard errors from −1.

Exps 135–151, where `[buf]` is constant in all seventeen runs, are the block
where `[S]` moves alone.

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
water run at 15–25 kJ/mol, so this is a covalent step. §7f measures the same
barrier on the window-free clock over all six temperatures and gets 77 ± 12;
the two agree, and neither is within a factor of three of a physical process.

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
`K + H₂O₂ ⇌ KP`, `MECHANISM.md` step 4 — then with `h = [H₂O₂]` in 100–6000×
excess over the catalyst the forward leg is pseudo-first-order and

    approach       1/τ = k_f·h + k_r
    destination    [KP]/E₀ = Kh/(1 + Kh),     K = k_f/k_r

**`1/τ` increases with `h` whatever the constants are.** There is no K and no
concentration at which more peroxide makes the approach slower, so a *positive*
order on the induction time is not a statement about being in the wrong regime —
it falsifies the scheme outright. That is what makes this worth measuring even
on a bad block.

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
signal-to-noise, at +0.619 ± 0.228. **This archive holds no clean peroxide lever
on an induction period.** One pre-incubation run would settle it and there will
not be one, so it goes to `COMPUTATIONAL.md` instead.

### 4b. The same scheme constrains both orders at once, and they violate it

The two rows of the table above are not independent. Taking log-log slopes,

    d ln v   / d ln h  =  1/(1 + Kh)        +1 unsaturated, → 0 saturated
    d ln τ   / d ln h  = −Kh/(1 + Kh)        0 unsaturated, → −1 saturated

so their **difference is 1 identically, for every K and every h**. "Is the rate
first order in peroxide" and "does peroxide shorten the induction" are one
question asked twice, and the constraint can be tested without knowing where on
the saturation curve the design sits. It is tested as *one* regression on
`log(v/t_ind)`, which also disposes of the correlation between two coefficients
fitted to the same curves.

| block | order of `v / t_ind` in [H₂O₂] | curves | from the required +1 |
|---|---|---|---|
| 4OMe peroxide, exps 127–131 | **+0.502 ± 0.194** | 15 | 2.6σ |
| BnOH two-axis, exps 135–151 | **+0.304 ± 0.188** | 110 | 3.7σ |

Both blocks fall short, in the same direction, by about the same amount.

*The floor moves this one.* A curve whose induction is shorter than one reading
is placed at `INDUCTION_FLOOR`, and `test_induction` plants an adduct fast
enough to be clipped and shows the bias is then **downward**, towards this
section's own reading. In this archive it runs the other way — the short
inductions are at *low* peroxide, so clipping pushes the coefficient **up**
towards +1 — and the sweep says so:

| floor | 1 s | 30 s | 60 s | 120 s | 300 s |
|---|---|---|---|---|---|
| 4OMe peroxide | +0.328 | +0.473 | **+0.502** | +0.542 | +0.600 |
| BnOH two-axis | −0.038 | +0.246 | **+0.304** | +0.404 | +0.566 |

The deviation from +1 is 1.9–3.7σ at every floor and never changes sign. The
module's 60 s floor is the conservative end of that range, not the flattering
one.

### And the rate is not first order in peroxide either

Which is worth stating separately, because "first order in H₂O₂" is the natural
expectation and it is the *unsaturated* limit of the scheme above, not a
general consequence of it. The clean design is the two-axis ladder — **63
curves in 17 runs over 2.45–163 mM**, one free level per run:

| | |
|---|---|
| free power law | a = **0.654** (0.534 to 0.775) |
| strict first order, a = 1 | **rejected, F = 32.0** |
| the scheme's own form, `v ∝ Kh/(1+Kh)` with the exponent held at 1 | fits **worse** than a free power, 11.36 against 10.48 |

So the order is fractional, and it is fractional *all the way across* a 67-fold
range rather than falling from 1 towards 0 the way one binding equilibrium
saturating would make it. That is the signature of several H₂O₂-dependent steps
in a chain — which is what the mechanism has, since peroxide is also the
stoichiometric oxidant and HOO⁻ is the nucleophile in steps 1–2.

### What the positive sign would mean, if it survives

Make the perhydrate a **trap** rather than the activation:

    K + H₂O₂ ⇌ KP      fast, OFF the activation path
    K → K*             slow, unimolecular — this is the induction

Only free catalyst can activate, so `1/τ = k_act/(1 + Kh)` and

    d ln τ / d ln h = +Kh/(1 + Kh),   which lies in (0, +1)

positive and bounded by 1, which is where both measurements sit. It keeps
everything §3 establishes — `k_act` carries no concentration, so no substrate
order, unimolecular, a clock, an amplitude that is a fraction — and the rate
keeps its own peroxide order from the oxidant being consumed downstream, which
is a different place in the scheme entirely. **It inverts step 4**: the
perhydrate becomes the resting state the catalyst has to leave.

Inverting the observed orders for K, and reading the same constant off the
rates through the saturating fit:

| route | K | as a molar constant | ΔG° |
|---|---|---|---|
| 4OMe induction order, at 40.8 mM | 0.0106 ± 0.0046 /mM | 11 /M | **−5.85 kJ/mol** |
| BnOH induction order, at 28.3 mM | 0.0286 ± 0.0180 /mM | 29 /M | **−8.31 kJ/mol** |
| BnOH rates, profiled | 0.0291 /mM (0.0134–0.0543) | 29 /M (13–54) | −8.36 (−6.43 to −9.90) |

Three routes within a factor of three, on two substrates and two blocks. **That
agreement is not evidence and should not be read as any.** Two of the three come
from blocks whose landmark tracks its own signal-to-noise, so their agreement is
also exactly what one shared artefact looks like; the third comes from a form
that fits worse than a free power law. What the range is good for is a target:
**a computed ΔG° of perhydrate formation near −6 to −10 kJ/mol would corroborate
the trap from a direction with no spectrophotometer in it**, and a strongly
negative one — a perhydrate that dominates at millimolar peroxide — would
contradict all three at once. That is `COMPUTATIONAL.md` C8's gate.

## 5. What the activation parameters say

Both rows are `arrhenius.activation_parameters`; neither is refitted here.
`1/τ` is the one parameter in that table whose Eyring entropy rests on no
assumption about the rate law — it is already a first-order relaxation rate in
s⁻¹, needing neither an extinction coefficient nor an enzyme concentration.

| | E<sub>a</sub> kJ/mol | ΔH‡ kJ/mol | ΔS‡ J/mol/K | ΔG‡(298) kJ/mol | curves |
|---|---|---|---|---|---|
| **the induction**, `1/τ` | 95.0 ± 15.7 | 92.6 ± 15.7 | +3.7 ± 53.2 | **91.48 ± 0.62** | 16, 4 T |
| **the turnover**, `v_peak` | 88.8 ± 1.8 | 86.3 ± 1.8 | −57.7 ± 5.9 | **103.47 ± 0.10** | 24, 6 T |

*The first row is over FOUR temperatures, 15 to 30 °C, because `inverse_tau` is
the ONE-phase fit's τ and above 32 °C a decelerating curve gives that form no
lag to find (`arrhenius.BURST_TRUSTWORTHY_BELOW_C`). §7f measures the same
barrier on the two-phase form over all six and gets **77 ± 12**; the two agree
and the six-temperature number is the one to quote.*

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

## 6. Which way the induction points, and what moves it

Everything above measures how *long* the induction is on curves that have one.
The archive also holds curves that begin **fast** and slow down, and `depth`
cannot see them: it is `1 − start/peak` and cannot go below zero, so a curve
that is fastest in its first window and one that starts a hair below its peak
both read as "no induction".

**Nothing needed refitting.** The two-phase form has carried both signs since it
was written — `A = c + v_ss·t − B₁(1−e^(−t/τ₁)) − B₂(1−e^(−t/τ₂))` with `B > 0`
a lag and `B < 0` a burst, and `TwoPhaseFit.kind` naming all four combinations.
What was missing is that `scope.frame` carried the time constants, the peak rate
and no sign at all, so no analysis in the package could ask which way a curve
pointed. It now carries `progress_kind`, `B_fast` and `B_slow`.

| block | curves | lag-first | the shapes |
|---|---|---|---|
| **4OMe catalysed** | 147 | **98** | 50 lag then fall, 46 lag, 21 burst, 20 mixed, 8 burst then fall, 2 two lags |
| the temperature series | 24 | **22** | 12 lag, 10 lag then fall, 1 mixed, 1 burst |
| **BnOH two-axis (135–151)** | 110 | **46** | 34 burst, 26 lag then fall, 20 lag, 19 mixed, 11 burst then fall |
| 4OMe enzyme-free | 49 | **10** | 28 burst, 11 burst then fall, 6 lag then fall, 4 lag |

**So yes, the two-axis curves have an induction — on 46 of 110 of them.** The
block is not one population: it splits almost evenly between curves that begin
below their eventual rate and curves that begin above it, which is why a single
induction time averaged over it has never meant much. The 4OMe blocks are not
split at all (98 of 147, and 22 of 24 in the temperature series), and the
enzyme-free block is split the other way (10 of 49) — the same
catalysed/enzyme-free contrast §2 makes, seen through the shape rather than
through the depth.

The sign is a **binary** here on purpose. The obvious continuous statistic is
the fast phase's amplitude normalised by the run's own signal, and it is
unusable: when the two exponentials are nearly degenerate the linear solve
trades enormous opposite amplitudes between them — exp 135 sample 3 returns
`B_fast = −241` and `B_slow = +303` on a curve that moves 0.06 AU — so
`B_fast/net` has an interquartile range of 1.8. The sign of that trade is stable
where its size is not.

### What moves it, and the two blocks disagree

A linear probability model: least squares of the 0/1 on `log(axis)` with one
offset per run. Not a logit, because 63 rows carrying 17 run offsets separate
perfectly in several runs and a logit's coefficient runs to infinity there. Read
it as a change in probability per e-fold, and no more than that.

| | axis | alone | with signal-to-noise |
|---|---|---|---|
| two-axis, **substrate arm** | [S] | −0.054 ± 0.042 | **−0.112 ± 0.052** |
| two-axis, **peroxide arm** | [H₂O₂] | +0.082 ± 0.046 | **−0.011 ± 0.069** |
| **4OMe catalysed** | [S] | +0.252 ± 0.059 | **+0.182 ± 0.073** |

Three things follow.

**The L has to be split.** Pooled over all seven cuvettes the sign looks like it
tracks peroxide hard — 10% lag-first below 10 mM against 56% above 30. Inside
the peroxide arm, where `[S]` is held at the run's top, that effect is gone
(−0.011 ± 0.069). The pooled version was reading the substrate arm's low-`[S]`
cuvettes, which all sit at the run's top peroxide.

**Signal-to-noise leans burst, and the substrate result runs against it.** The
share of lag-first curves rises from 0.29 to 0.50 across the two-axis
signal-to-noise quartiles: a curve with little signal gives the two-phase fit
little to choose between the shapes. That is the same objection figure F raises
against the induction *time*, and here it works in the analysis's favour — more
substrate means more signal, so the confound pushes towards *lag*, and the
measured substrate effect is towards **burst** anyway.

**And the two blocks give opposite signs**, −0.112 ± 0.052 against
+0.182 ± 0.073. They differ in exactly one structural way.

### The buffer is the candidate, and the one direct lever agrees

`[S]` and `[buf]` correlate at **−0.96** inside every 4OMe run and at **0.00**
inside every two-axis run — `[buf]` is constant across all seven cuvettes of all
seventeen. So in the 4OMe block "more substrate" also means "less buffer", and
the two blocks can disagree about the substrate without either being wrong:

> two-axis, buffer fixed: **more substrate → burst**
> 4OMe, buffer falling with substrate: **more substrate → lag**, which is the
> same as **more buffer → burst**

Read as a buffer effect, that is what **general acid/base catalysis of E → E\***
predicts — and [`buffer/`](../buffer/ANALYSIS.md) is where that hypothesis is
taken apart. It survives, along with general acid and "spectator": the archive's
two-pH species test measures +1.06 ± 0.77 against 1.76 and 0.52, and excludes
nothing: more buffer base, faster activation, and the lag over before the run
is properly under way. It is the mechanism the temperature series' own numbers
already hint at from a different direction — the catalysed buffer order of the
*rate* is +0.400 ± 0.028 at 50–200 mM and +0.803 ± 0.173 below 25, so the buffer
is doing something to this chemistry that saturates.

**The one direct lever, and what it took to read it.** Exps 32 and 34 hold
substrate, peroxide, pH and temperature fixed and step `[buf]` 3.125 → 200 mM,
which is 64-fold and is the only such design in the archive. Eight curves in two
runs, and reading them together takes two corrections.

*The two runs earn different model forms.* Every curve of exp 34 earns the
two-phase form (**F = 71 to 819** against a threshold of 12) and every curve of
exp 32 earns the one-phase form (**F = 1.6 to 4.7**), so `tau_fast` is τ₁ of a
two-phase fit in one run and τ of a one-phase fit in the other. The break is at
the run boundary because the *schedule* is: exp 34 ran **5280 s** and exp 32
**1767 s**, so exp 34's runs are long enough to contain the slow fall and exp
32's end before it. Regressing `tau_fast` inside each run — which is what this
folder did first — reported +0.457 ± 0.097 and −1.052 ± 0.469 and called it a
disagreement. It was never a comparison of two measurements of one thing.

*And the two runs sit at different levels.* `v_peak` — the one quantity defined
identically on both forms — **falls 1.80×** across the join, from 7.41 × 10⁻⁵ at
25 mM to 4.13 × 10⁻⁵ at 50 mM, while the buffer doubles. Whatever separates the
two days is larger than the effect being measured, so a fit with a shared
intercept measures the day. This is the same reason `temperature_series` §3
quotes the buffer order of the rate as two range-specific numbers rather than
one pooled one.

Read instead through a landmark on the *readings*, with a window in **seconds**
common to both runs (450 s; these are 40 °C curves reaching 0.05–0.31 AU, so
they can afford a fixed window where the 15 °C curves cannot), and pooled with
one level per run:

| | d log t<sub>ind</sub> / d log[buf] | d log depth / d log[buf] |
|---|---|---|
| exp 34, 3.125–25 mM | −0.300 ± 0.269 | −1.426 ± 0.398 |
| exp 32, 50–200 mM | −0.728 ± 0.316 | −0.171 ± 0.415 |
| **pooled, one level each** | **−0.433 ± 0.201** | **−1.036 ± 0.364** |

**The two runs agree, and they agree in the direction general base catalysis
predicts**: more buffer, shorter and shallower induction. Across windows of
300 to 1200 s the pooled slope runs −0.33 to −0.72 and never changes sign.

It is eight curves. Two runs, one temperature, one substrate, one pH, and the
step from 25 to 50 mM is also a step between experiments — the pooling handles
the level and cannot handle a *slope* that differs between days. **Treat it as
one measurement with a direction and not as a buffer order.** What makes it
worth having at all is that it agrees with the 28-run indirect signal above,
which is a different block, a different substrate and a different statistic.

**This is a design, not a result.** The measurement that would settle it is one
run with `[buf]` stepped at fixed `[S]` on the *catalysed* 4OMe system at 25 °C,
where the induction is thousands of seconds long instead of nearly over. It
belongs beside the pre-incubation experiment in §6 as something the archive
should have and does not.

### And on this axis §4b's constraint is met

§4b's `+1` is not a fact about hydrogen peroxide. It is a fact about any species
X held in excess that draws the catalyst into its active form: `1/τ = k_f[X] +
k_r` and `[E*]/E₀ = K[X]/(1 + K[X])` give `d ln v/d ln[X] − d ln τ/d ln[X] = 1`
for every K and every `[X]`, with nothing left to fit. Put the buffer in that
role — as a general base, or through a peroxo adduct of the buffer itself, the
algebra does not distinguish them — and the constraint transfers unchanged.
`joint_buffer_order` asks it as one regression on the same eight curves, so the
covariance between the two orders is the fit's own:

| window | order(*v*/t<sub>ind</sub>) | from the required +1 | signal control |
|---|---|---|---|
| 300 s | **+0.989 ± 0.159** | 0.1σ | −0.145 ± 0.266 passes |
| 450 s | **+1.094 ± 0.150** | 0.6σ | −0.275 ± 0.268 passes |
| 600 s | +1.379 ± 0.180 | 2.1σ | −0.660 ± 0.230 **fails** |
| 900 s | +1.339 ± 0.207 | 1.6σ | −0.643 ± 0.225 **fails** |
| 1200 s | +1.232 ± 0.236 | 1.0σ | −0.565 ± 0.224 **fails** |

**The buffer axis meets the constraint that the peroxide axis misses** at 2.6σ
and 3.7σ (§4b). The windows that fail their own signal control are exactly the
windows that overshoot, and they overshoot in the direction the artefact
predicts — more buffer, more signal, an earlier landmark, a larger gap — so the
two windows that pass are the ones to read.

It is still eight curves, and a constraint being met is weaker evidence than one
being violated: `+1` is where a *rate* order of about +0.66 and an *induction*
order of about −0.43 happen to land, and neither was chosen. What it does say is
that if any species in this archive is in pre-equilibrium with the catalyst
before it turns over, the buffer behaves like that species and H₂O₂ does not.
[`../buffer/`](../buffer/ANALYSIS.md) §6 takes it from here, including the
scheme in which the buffer's role is to carry the peroxide rather than to be a
base — which this archive cannot separate from general base catalysis, because
**0 of its 88 runs step `[buf]` and `[H₂O₂]` at once**.

## 7. The seven variables, and where the archive can answer for each

Everything above measures the induction against **two** of the seven things the
experiment moves: `[S]` and `[H2O2]`. That is not a choice of interest, it is
the statistic. `induction_point` reads through a rolling window a tenth of the
run, so it is safe inside a run and not between runs (§1) — and **temperature,
pH, `[buf]`'s partner the buffer salt, and the substrate are one value per run
in this archive, always.** Five of the seven questions were unaskable.

`scope.frame` now carries `lag_depth` and `lag_half_s`, read off the **fitted
rate** rather than a rolling one: how far below its peak the rate starts, and
how long the fit takes to cover half that gap. On a one-phase lag the second is
exactly τ·ln2; on a two-phase curve it is the mixture a single time constant
cannot be. `summary_kinetics.ProgressFit.lag_profile` is the definition, it
recovers a planted τ to within 6% at 200, 800 and 3000 s, and it returns zero on
a curve that begins at its fastest. Both come off the gas-corrected readings,
for the reason `tau_corrected` does (§4b).

**It carries no window. It still carries the schedule and the signal**, and
this section is built around both.

- The **schedule**, because the fit's τ grid runs from span/300 to 2·span: a
  1260 s run cannot report a 4000 s clock. `induction.lag_window_frame` refits
  every run of a block on a window they all share.
- The **signal**, because a curve with little signal gives the fit little to
  choose between the shapes. Every result below is quoted twice, with and
  without the curve's own log(net/noise) held alongside, and `lag_orders`
  reports the partial correlation that says whether holding it could matter.

### 7a. Which axes the archive moves, and where

`induction.lag_identifiability`, over all 88 experiments:

| axis | runs moving it internally | widest span inside one run |
|---|---|---|
| `[S]` | **77 of 88** | 50.1× |
| `[buf]` | **53 of 88** | 8.0× |
| `[H2O2]` | **20 of 88** | 66.7× |
| pH | **0** | — |
| temperature | **0** | — |
| buffer salt | **0** | — |
| substrate | **0** | — |

Three axes are within-run and four are not. The four have to be read off
between-run designs, and the archive holds exactly these: four pH ladders
(§7e), one temperature series (§7f), and five matched **pairs** — one for the
buffer salt, two for the substrate, two for the catalyst loading (§7g).

### 7b. The lag is not one phenomenon, and the substrates split it

`induction.lag_channel_table`, over the whole archive:

| substrate | channel | curves | with a lag | median depth | median clock |
|---|---|---|---|---|---|
| 4OMe | enzyme-free | 49 | **10** | 0.000 | 0 s |
| 4OMe | **catalysed** | 147 | **109** | **0.257** | 304 s |
| BnOH | enzyme-free | 26 | **14** | **0.138** | 69 s |
| BnOH | catalysed | 164 | **77** | 0.000 | 0 s |

The 4OMe row is §2 again through a different statistic and it holds. **The BnOH
row does not, and that is a second phenomenon rather than a counter-example.**
The enzyme-free BnOH lags are exps 3 and 6 — pH 6.71 phosphate, 11025 and
17934 s, 8 of their 10 curves — and exp 65, whose four cuvettes share the
break `scope.synchronised_break` measures.

The acceleration test says the same thing the shape does, on the same curves:
**4 of 26 enzyme-free BnOH curves accelerate past 3σ and 0 of 49 enzyme-free
4OMe curves do**. Steps 1–3 of `MECHANISM.md` are autocatalytic in the product
and need no catalyst at all, and the substituent argument under `MECHANISM.md`
S2 says why they should run on one substrate and not the other: 4-methoxy makes
the aldehyde a better target for an electrophilic oxidant and a worse hydride
*donor*, and the hydride transfer is step 2's rate-determining step. An
accelerating enzyme-free BnOH curve is what that predicts.

So: on 4OMe a lag is the catalyst becoming active; on BnOH a lag can also be the
catalyst-free loop finding its own product. **Do not pool them**, and read every
number below with its substrate attached.

### 7c. The replicate floor — what two identical runs do anyway

`scope.REPLICATE_RUNS` is exps 2, 4, 5 and 7: four repeats of one composition,
each 17934 s, the only four-fold repeat in the archive. It is the bar every
between-run number in this section has to clear.

| | exp 2 | exp 4 | exp 5 | exp 7 | spread |
|---|---|---|---|---|---|
| clock, s | 3675 | 3682 | 3388 | 2951 | **1.25×** |
| depth | 0.388 | 0.734 | 1.000 | 1.000 | **0.612** |

**The clock reproduces between runs and the depth does not.** A factor of 1.25
on the clock is 0.22 in log units, so a between-run clock effect has to beat
that; a depth spread of 0.61 on a statistic bounded by [0, 1] means the depth
is a within-run statistic and nothing else. Every between-run row below is read
on the clock for that reason.

### 7d. Within runs: `[S]`, `[H2O2]`, `[buf]`

`induction.lag_orders`, which is `scope.orders` — the package's own
within-experiment log-log fit, one offset per experiment — with the signal
pairing added. **Every axis the block moves is fitted at once**: one axis at a
time is a different regression on an L, and a substrate-only fit of the
two-axis clock reads −0.453 ± 0.107 where the joint fit reads −0.225 ± 0.115.

| block | axis | r(axis, signal) | clock order | with the signal held |
|---|---|---|---|---|
| 4OMe catalysed, 147 curves | `[S]` | +0.60 | +0.003 ± 0.164 | +0.017 ± 0.222 |
| | `[buf]` | −0.09 | −0.203 ± 0.297 | −0.194 ± 0.313 |
| | `[H2O2]` | +0.33 | +0.194 ± 0.117 | +0.199 ± 0.130 |
| **BnOH two-axis, 110 curves** | **`[S]`** | +0.04 | **−0.225 ± 0.115** | **−0.336 ± 0.127** |
| | `[H2O2]` | +0.57 | +0.679 ± 0.171 | +0.377 ± 0.232 |
| 4OMe peroxide, 127–131, 15 | `[H2O2]` | +0.67 | +0.194 ± 0.165 | −0.134 ± 0.160 |
| buffer titrations, 20 | `[buf]` | +0.63 | −0.349 ± 0.215 | +0.076 ± 0.222 |

**Read the third column first.** The two-axis block fails its own signal control
outright (+0.916 ± 0.259) and can *still* be asked about substrate, because in
that block substrate buys no signal: `vmax_corrected` carries a substrate order
of only +0.09 ± 0.05 over the same 110 curves — the published +0.01 ± 0.04 is
`vmax` over the eleven strong runs, a different cut of the same flatness — so
log[S] and log(net/noise) correlate at **+0.04** once the run offsets are out.
Peroxide is a different matter everywhere — it *is* the signal, at +0.57 and
+0.67 — and holding it takes the peroxide coefficient from +0.68 to +0.38 on one
block and from +0.19 to −0.13 on the other. **This archive cannot measure a
peroxide dependence of the induction**, through the landmark (§4a) or through
the fit, and the reason is the same both times.

*Two of these rows are the same measurement.* With one offset per experiment,
the only `[H2O2]` variation the 4OMe block has is inside exps 127–131 — every
other catalysed 4OMe run sits at 82.5 mM — so the two peroxide rows return the
same +0.194 with different errors. They are one result reported twice, not two
that agree.

The same objection sinks the buffer titrations: `[buf]` correlates with the
signal at +0.63 there, and holding it moves −0.349 ± 0.215 to +0.076 ± 0.222.
**The direct buffer result of §6 is not reproduced by the window-free clock.**

The 4OMe block's own `[buf]` row is the only one whose coefficient is unmoved by
the control (−0.203 → −0.194), and it is worth saying exactly why rather than
calling it clean. The block holds **two** sources of buffer variation: the five
titrations, where `[buf]` moves alone and tracks the signal at +0.63, and every
other 4OMe run, where substrate volume displaced buffer volume so `[buf]` falls
as `[S]` rises at −0.96 in logs. Neither is usable on its own — drop the
titrations and the remaining collinearity puts the buffer order at
+2.20 ± 1.27 — and pooled, their opposite signal correlations very nearly
cancel, leaving `[buf]` almost orthogonal to the signal (r = −0.09). So the
−0.20 ± 0.30 is a real fit rather than a confound, and its error covers
everything from §6's −0.433 to zero.

*The floor moves the sizes and not the signs.* 48 of the two-axis block's 110
live curves begin at their fastest and sit on the floor honestly. Sweeping it
from 1 s to 300 s:

| floor | 1 s | 30 s | 60 s | 120 s | 300 s |
|---|---|---|---|---|---|
| two-axis `[S]`, held | −0.715 | −0.400 | **−0.336** | −0.272 | −0.178 |
| two-axis `[H2O2]`, held | +0.755 | +0.441 | **+0.377** | +0.313 | +0.243 |
| 4OMe `[S]`, held | +0.620 | +0.124 | **+0.017** | −0.058 | −0.097 |

The substrate coefficient on BnOH is negative at every floor and the 4OMe one
changes sign, which is what "a null" and "a weak effect" look like from the
inside. Quote the BnOH substrate order as a range, **−0.18 to −0.71**, not as
a number.

**And product control is excluded a second time, on both substrates.** §3's
route one, re-run on the window-free clock:

| block | d log(clock) / d log(rate) | curves |
|---|---|---|
| 4OMe catalysed | **+0.245 ± 0.127** | 147 |
| BnOH two-axis | **+0.897 ± 0.178** | 110 |
| the temperature series | −0.073 ± 0.468 | 24 |

A clock predicts 0 and product control −1. The 4OMe row excludes product control
at 9.8σ, agreeing with the landmark's −0.025 ± 0.109 from a statistic that
shares no window with it. The BnOH row excludes it at 10.6σ and misses the clock
too, in the direction its own peroxide confound predicts — faster cuvettes there
are the high-peroxide ones, and peroxide lengthens the apparent clock.

### 7e. Between runs: pH, on four ladders

pH is one value per run everywhere, so a pH ladder is a set of runs sharing
everything else. There are four, and until now the project had used none of
them for this question. **A per-experiment offset is not available here** — it
would *be* the pH axis; `lag_ladder` raises rather than returning the
minimum-norm split, which is `scope._moves` one level out and which reported a
forty-sigma pH effect for one afternoon before the guard existed.

Every run is refit on the window they all share, and each ladder is confounded
differently — which is the reason to have four:

| ladder | runs in it | kept | window | pH | r(pH, length) | r(pH, signal) | clock, held |
|---|---|---|---|---|---|---|---|
| 4OMe phosphate | 9 | **6** | 1470 s | 5.64–8.95 | −0.25 | **+0.87** | +0.094 ± 0.303 |
| 4OMe boric | 9 | 9 | 1260 s | 8.46–10.34 | **+0.71** | −0.65 | +0.286 ± 0.322 |
| BnOH pyrophosphate 136–142 | 7 | 7 | 3720 s | 6.95–9.43 | −0.53 | +0.77 | +0.472 ± 0.343 |
| BnOH pyrophosphate 143–151 | 9 | 9 | 3000 s | 5.47–9.73 | **−0.79** | +0.79 | +0.405 ± 0.187 |
| **pooled** | | | | | | | **+0.336 ± 0.132** |

The phosphate ladder is the one that pays for the common window: three of its
nine runs are shorter than 1470 s and leave fewer than eight readings when
truncated, so they drop out. The other three ladders keep every run.

The four agree (χ² = 0.95 on 3) across two substrates, three buffers, and
confounds of both signs. **The pooled value moves with the window, though never
in sign**: at 0.75 and 0.5 of the shared window it is +0.119 ± 0.115 and
+0.176 ± 0.113. So the honest quote is

> **d ln τ / d pH = +0.12 to +0.34, positive at every window and never past
> 2.6σ.**

Positive means **more alkaline, longer induction** — the same direction §4a's
peroxide sign points, on an axis that is not the peroxide. As a saturation
fraction (`saturation_fraction`, dividing by ln 10) it is **0.05 to 0.15**: if
something is holding the catalyst off its activation path, the archive's whole
pH range moves it across at most a seventh of its range. The depth carries
nothing at all — pooled +0.020 ± 0.039.

### 7f. Temperature: the barrier, and what the schedule was doing to it

This is the row the schedule bites hardest, and the one where controlling for it
and refitting for it give different answers.

Cold runs are the long runs here — 1/T against log(run length) at **+0.66** —
and a clock cannot exceed its run. On whole runs the barrier is
**83.7 ± 8.9 kJ/mol**; put log(run length) into that regression and it drops to
**55.3 ± 24.3** on a six-point collinearity, which is not a measurement. Refit
every run on a window all six share and the collinearity is gone by
construction:

| window | 0.3× | 0.5× | 0.75× | 1.0× | 2× | whole runs |
|---|---|---|---|---|---|---|
| seconds | 1470 | 2499 | 3773 | 5047 | 10094 | 17934 |
| temperatures | 4 | 5 | 5 | 6 | 6 | 6 |
| E<sub>a</sub>, kJ/mol | 8.8 ± 37.5 | 1.7 ± 35.7 | 81.6 ± 22.0 | **76.7 ± 12.2** | 76.6 ± 8.7 | 83.7 ± 8.9 |

Stable at **77–84 kJ/mol** from three quarters of the shared window upwards, and
collapsing below it — where the 15 and 20 °C clocks are longer than the window
and the fit cannot see them, which the temperature count reports rather than
hides. The half-rise runs 2212, 2000, 2204, 630, 446 and 334 s from 15 to 40 °C.

**This replaces §5's 95.0 ± 15.7 kJ/mol as the number to quote.** That figure
comes from `arrhenius.inverse_tau`, which is the *one-phase* fit's τ and
therefore reaches only 15–30 °C: above 32 °C a decelerating curve gives the
one-phase form no lag to find and τ lands on the top of its grid
(`arrhenius.BURST_TRUSTWORTHY_BELOW_C`). The two-phase form does find it. The two
agree inside their errors and the new one is over six temperatures instead of
four, so **quote 77 ± 12 with a window systematic of about ±7** — and the
conclusion §5 drew from it is unchanged: a barrier of that size is three to four
times what dissolution, diffusion or thermal equilibration cost, and it is
covalent.

### 7g. The three pairs: buffer salt, substrate, catalyst

None of these three was ever laddered, so each is a pair read at a common
window, and the bar is §7c's replicate floor — **1.25× on the clock**.

| contrast | runs | clock | depth | verdict |
|---|---|---|---|---|
| **buffer salt**, 4OMe pH 8.96 | 12 phosphate / 42 boric | 0 s / 56 s | 0.000 / 0.054 | neither has a lag; nothing measured |
| **substrate**, boric pH 9.0 | 42 4OMe / 51 BnOH | 57 s / 0 s | 0.054 / 0.000 | both inside the floor |
| **substrate**, boric pH 9.7 | 45 4OMe / 55 BnOH | 853 s / 1644 s | 0.217 / 0.250 | 1.9×, just past the floor, one pair |
| **catalyst**, BnOH boric pH 8.51 | 59 (0.028) / 60 (0.014) | 0 s / 0 s | 0.000 / 0.000 | neither has a lag |
| **catalyst**, two-axis | 140 (0.034) / 141 (0.014) | 913 s / 517 s | 0.618 / 0.162 | 1.8×, and the wrong way |

**The substrate pairs look like a null and are not.** §7b's table has 109 of
147 catalysed 4OMe curves carrying a lag against 77 of 164 catalysed BnOH ones,
and the obvious reading is that the induction is bigger on the 4-methoxy
substrate. The two matched pairs say it is not the substrate: at pH 9.01 the two
give **2 of 4 against 1 of 4** lag-first, and at pH 9.70 they give **4 of 4
against 3 of 4**, with depths of 0.054/0.000 and 0.217/0.250 — inside the
replicate floor on both. **What separates the two catalysed blocks is their
conditions**, and they are not close: median `[enz]` 0.241 against 0.025 mM,
median signal 0.129 against 0.029 AU, median pH 7.5 against 8.5. So the archive-wide contrast is a
conditions contrast, the pairs are what shows it, and **the activation is a
property of the catalyst rather than of the alcohol** — which is what an
activation step upstream of turnover has to be.

**The catalyst pair is the one that matters and the archive spends it badly.** A
unimolecular activation of the catalyst predicts the clock does *not* depend on
how much catalyst there is — the one prediction of the activation reading that a
concentration can falsify. Exps 59 and 60 hold everything but `[enz]` and
neither has a lag to time. Exps 140 and 141 do, and 2.4× more catalyst gives a
1.8× **longer** clock, which is past the replicate floor and points the wrong
way for a unimolecular step. It is two runs, 0.07 pH units apart, in the block
whose signal control fails. **It is a flag, not a result** — and it is the
cheapest experiment left undone in this folder: one substrate ladder at two
catalyst loadings on the catalysed 4OMe system at 25 °C.

### 7h. What the dependences add up to

Collecting the rows that survive their own controls:

| variable | the induction clock |
|---|---|
| temperature | **77 ± 12 kJ/mol**, six temperatures, common window |
| `[S]`, BnOH, buffer fixed | **−0.18 to −0.71** across floors, always negative |
| `[S]`, 4OMe, buffer falling with it | **0.00 ± 0.16**, sign changes across floors |
| pH, four ladders | **+0.12 to +0.34** per unit, positive at every window |
| `[buf]`, 4OMe, unconfounded axis | −0.20 ± 0.30 |
| `[H2O2]` | **not identifiable** — it is the signal, in every block that moves it |
| buffer salt | not measured — one pair, neither run has a lag |
| substrate | **no effect**, at two matched pairs — the archive-wide contrast is conditions |
| `[enz]` | not measured, and the one pair points the wrong way |

**Three of the seven are measured, three are not, and `[S]` is measured on one
substrate only.** That is the honest state, and it is worth being explicit that
the two which would decide the most — `[H2O2]` and `[enz]` — are exactly the two
the archive cannot answer.

What the three that survive say, taken together:

**The activation is covalent, does not wait for product, is slowed by base and
is hurried by the substrate.** The barrier is covalent-sized on both routes that
can measure it. Nothing about it waits for product: regressed on the curve's own
rate the clock gives +0.245 ± 0.127 on 4OMe and +0.897 ± 0.178 on BnOH, where a
product threshold needs −1, so product control is excluded on both substrates by
a statistic that shares no window with §3's. (Neither row is a clock either —
the BnOH one is far from 0, in the direction its own peroxide confound predicts,
since the faster cuvettes there are the high-peroxide ones.) It is slowed by
alkali on four independent ladders. And on the one block where the substrate
moves with the buffer held fixed, more substrate *shortens* it.

Two bounded schemes and no others fit that shape, and they are §4a's pair with
the species left open:

    a TRAP        K·X ⇌ K + X,  K → K*      d ln τ / d ln[X] ∈ (0, +1)
    an ACTIVATOR  K + Y ⇌ K·Y → K*          d ln τ / d ln[Y] ∈ (−1, 0)

Neither bound contains a rate constant, so a coefficient outside it falsifies the
scheme rather than fitting it. The pH result puts **X at 5–15% saturated** across
the archive's whole pH range; the BnOH substrate result puts **Y at 18–72%**,
with the caveat that it rests on one block. A species whose abundance rises with
pH and holds the catalyst back, and a substrate that pulls it forward, is exactly
what a cyclodextrin with a ketone at its mouth would give: **the resting form is
an anionic adduct at the carbonyl — the gem-diolate, or the hydroperoxide adduct
of §4a — and the alcohol binding in the cavity is what displaces it.**

That is a reading and not a result. Three things would decide it, and two of
them are already in `COMPUTATIONAL.md`: **C7**, the hydration equilibrium and the
pK<sub>a</sub> of the hydrate — a hydrate pK<sub>a</sub> near 10 would put the
gem-diolate exactly where the boric ladder sits, and one near 13 would rule it
out; **C8**, the free energy of the perhydrate, whose gate §4b already set at
−6 to −10 kJ/mol; and the `[enz]` ladder of §7g, which the archive does not
hold.

## 8. What this settles, and what it does not

**Settled.**

- The induction is a property of the **catalysed** reaction. It is absent from
  all 49 enzyme-free 4OMe curves at matched composition, including a 5-hour run.
- It ends on a **clock**, not at a product threshold: the coefficient is
  −0.025 ± 0.109 where product control requires −1, on 147 curves, and the
  bias-free route agrees at −0.26 ± 0.32. Asked in the units the product
  hypothesis is stated in — is the product made by the landmark fixed, or is it
  whatever the rate made in a fixed time — the answer is **+0.885 ± 0.125**
  against the +1.000 a clock reads back and the +0.42 a *planted* threshold
  reads back through this landmark, so **3.7σ**. §3, route three; and read the
  calibration before quoting it, because the nominal value would say 7.1σ.
- Its **amplitude is a fraction** of the eventual rate and does not scale with
  the substrate.
- It is a **covalent step**: **77 ± 12 kJ/mol** over six temperatures at a
  window all six runs share (§7f), with a window systematic of about ±7. That
  is three to four times what dissolution, diffusion or thermal equilibration
  cost. It replaces the 95 ± 16 of §5, which is the same barrier measured on
  the one-phase fit's τ and therefore over four temperatures rather than six;
  the two agree inside their errors.
- It is **faster than turnover** by 126× in free energy, as it must be.
- The **sign** of the early curve is a property of the block: 98 of 147
  catalysed 4OMe curves begin below their eventual rate, against 10 of 49
  enzyme-free ones. §6.
- **The lag is two phenomena and the substrate separates them.** On 4OMe it
  needs the catalyst (10 of 49 enzyme-free curves, median depth 0.000); on BnOH
  it does not (14 of 26, median 0.138), and those sit in the long low-pH runs
  and in exp 65's shared break. The acceleration test agrees with the shape:
  4 of 26 enzyme-free BnOH curves accelerate past 3σ and 0 of 49 enzyme-free
  4OMe curves do, which is the catalyst-free loop of `MECHANISM.md` steps 1–3
  running on one substrate and not the other. §7b. Do not pool the two
  substrates.
- **`[H2O2]` cannot be asked**, and now for a second reason. The landmark
  tracks signal-to-noise on both blocks that move peroxide (§4a) and so does
  the window-free clock, because in those blocks the peroxide *is* the signal:
  the partial correlation is +0.57 and +0.67, and holding it takes the clock's
  peroxide order from +0.68 to +0.38 and from +0.19 to −0.13. §7d.

**Not settled, and the reasons are different in each case.**

- **Which step it is.** "The catalyst becomes active" is a shape, not a
  mechanism. The candidates that survive §3 are all unimolecular in everything
  the cuvette holds: the ketone's gem-diol **hydrate** dehydrating to the free
  ketone, the perhydrate collapsing to the dioxirane, or a conformational change
  of the cyclodextrin. Absorbance at one wavelength cannot choose between them.
- **Whether the peroxide is involved at all.** §4a–4b. Three things point the
  same way and none of them is clean: the induction's peroxide order has the
  wrong *sign* for an adduct, the joint constraint the scheme puts on both
  orders at once falls short by 2.6σ and 3.7σ in the two blocks that can test
  it, and the rate is not first order in peroxide either (a = 0.654, first
  order rejected at F = 32). But both induction orders come from blocks whose
  landmark tracks its own signal-to-noise, so their agreement is also what one
  shared artefact looks like. What the sign would mean if real is an inversion
  of step 4 — the perhydrate as a resting state rather than the activation —
  and C8 now has a number to test it with.
- **Whether the same thing happens with BnOH.** Exps 135–151 have the peroxide
  design and not the signal-to-noise; their induction statistic tracks the noise
  (+0.619 ± 0.228). The substrate comparison `product_fate` could make at the
  end of the curve cannot be made at the start. What can be said is that the
  block is not one population: 46 of its 110 live curves begin below their
  eventual rate and 45 begin above it (§6).
- **Whether the catalyst loading moves the clock**, which a unimolecular
  activation says it must not. This is the one prediction of the whole
  activation reading that a concentration can falsify, and the archive spends
  its two chances badly: exps 59 and 60 hold everything but `[enz]` and neither
  run has a lag to time, and exps 140/141 give a 1.8× **longer** clock at 2.4×
  the catalyst — past the replicate floor, the wrong way, on two runs 0.07 pH
  units apart in the block whose signal control fails. §7g. **A flag, not a
  result**, and the cheapest experiment this folder still wants.
- **Whether the buffer is what carries E → E\*.** §6, and this is the most
  promising open lead in the folder. §7d weakens the direct half of it: the
  five titrations' `[buf]` axis correlates with the signal at +0.63, and
  holding the signal moves the window-free clock's buffer order from
  −0.349 ± 0.215 to +0.076 ± 0.222. What survives is the 4OMe archive's own
  `[buf]` axis, the one buffer axis in the archive that is not
  signal-confounded (r = −0.09), at −0.20 ± 0.30 — the right sign, no
  significance. Three things point the same way: the 4OMe
  and two-axis blocks give *opposite* substrate effects on the sign and differ
  in exactly one structural way (`[buf]` falls with `[S]` at −0.96 in one and is
  constant in the other); the one direct buffer ladder gives
  **−0.433 ± 0.201** on the induction time, more buffer meaning a shorter
  induction; and general acid/base catalysis is what a unimolecular activation
  on a carbonyl centre would be expected to show. It is eight curves in two
  runs at one temperature, so it is a direction and not a buffer order.

**The three things that would finish it.**

1. `COMPUTATIONAL.md` **C7** — the hydration equilibrium and the dehydration
   barrier of the chemzyme's ketone in water, against the barrier for adding
   H₂O₂ to it. If the hydrate dominates and dehydration costs ≈ 95 kJ/mol, the
   induction has a name.
2. `COMPUTATIONAL.md` **C8** — the free-energy profile from the perhydrate to
   the dioxirane, to say whether the peroxide adduct is on the activation path
   or off it. §4b gives it a gate it did not have: a computed ΔG° of perhydrate
   formation of **−6 to −10 kJ/mol** would corroborate the trap, and anything
   much more negative would contradict all three routes at once.
3. Two experiments that will not be run. **Pre-incubate the catalyst with H₂O₂,
   then add substrate.** If the induction is catalyst activation it disappears;
   if it is anything that waits for product it does not. One cuvette. And
   **step `[buf]` at fixed `[S]` on the catalysed 4OMe system at 25 °C**, where
   the induction is thousands of seconds long rather than nearly over — the
   archive's only buffer titration sits at 40 °C, where it is not.
