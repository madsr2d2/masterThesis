# What the 4OMe product does

Every catalysed 4-methoxybenzyl alcohol run in this archive rises to a maximum
rate and then falls, at **0.3–3% conversion** — far too little for the substrate
to be running out. This is what the fall is, and why benzyl alcohol's curves do
the opposite.

The evidence is drawn from the **whole 4OMe archive**, not from one block. The
temperature series (exps 14–19) is where the question came from and is 23 live
curves; the discrimination that answers it needs run length and product to vary
independently, which takes 84.

    data/slowdown.py               the models, the landmarks, the regressions
    python data/slowdown.py        the whole argument, printed
    python data/test_slowdown.py   parameter recovery and the published numbers
    python product_fate/build_figures.py
    python product_fate/check_numbers.py

**Figures**: [`index.html`](index.html) is the presentation — six figures, A to
F, one per claim below. It is rebuilt by `build_figures.py`, which computes
nothing.

Related: [`../temperature_series/`](../temperature_series/ANALYSIS.md) for the
activation parameters this bears on, [`../background_reaction/`](../background_reaction/ANALYSIS.md)
for the enzyme-free channel, and `COMPUTATIONAL.md` C5 and C6 for the two
calculations that would finish it.

## 1. The four candidates, and why two were already gone

| | | |
|---|---|---|
| the substrate runs out | **excluded** | conversion is 0.3–3% |
| the catalyst dies | **excluded** for BnOH | Selwyn on exps 59/60; the ratio never falls below 1 on 20 comparisons, and inactivation requires it to |
| a **clock** — something shared is consumed on its own schedule | `A′ = v(t)·e^(−kt)` | untested until now |
| the **product** sets the rate | `A′ = v(t) − kA` or `A′ = v(t)/(1 + A/K_i)` | untested until now |

**One progress curve cannot separate the last two.** Within a curve the product
only ever grows with time, so "the rate fell after 2000 s" and "the rate fell
after 0.1 AU" are the same sentence. The separation has to come from a design
where run length and product move independently — and the archive has one, by
accident: it holds 1470 s runs that reached 0.27 AU and 17934 s runs that
reached 0.045. Across the catalysed 4OMe phosphate set the two regressors
correlate only **−0.39**.

## 2. It is the product, and it is not the clock

Each curve's deceleration is `scope`'s own `late_over_early` — the last fifth's
slope over the first fifth's — regressed on the run's length and on the product
it made, in **mM** rather than absorbance, so the two substrates are compared as
concentrations across their factor-of-6.1 difference in ε.

| block | curves | fell with **run length** | fell with **product** |
|---|---|---|---|
| **the temperature series** | 23 | +0.118 ± 0.122 | **−0.525 ± 0.071** |
| 4OMe catalysed, phosphate | 84 | −0.124 ± 0.079 | **−0.598 ± 0.053** |
| …with pH and temperature in the model | 84 | −0.020 ± 0.107 | **−0.670 ± 0.107** |
| …with one dummy per experiment | 84 | — | **−0.919 ± 0.161** |
| **4OMe, no enzyme at all** | 49 | **−0.361 ± 0.047** | +0.025 ± 0.048 |
| BnOH, exps 135–151 | 84 | **−0.697 ± 0.136** | +0.283 ± 0.129 |
| BnOH catalysed, every buffer | 137 | −0.138 ± 0.075 | +0.119 ± 0.096 |
| …with one dummy per experiment | 137 | — | +0.245 ± 0.121 |

Read the first row and the fifth together. **The temperature series slows in
proportion to what it has made and not at all with how long it has run.** The
same 4OMe chemistry with **no enzyme in the cuvette** does the exact opposite:
it is a clock, at 7.7σ, with no product term at all. That block is not a
curiosity — it is the reference channel every catalysed curve is measured
against, and its decay is subtracted out per cuvette.

**It is not a range artefact.** The two blocks reach different product levels,
so restrict both to the window they share, 0.004–0.077 mM: the catalysed set
still gives **−0.560 ± 0.090** on 58 curves and the enzyme-free set
**+0.100 ± 0.092** on 37.

**It is not the photometer either.** This is the objection the design invites and
it took a while to see: ε is **7.53** for the 4OMe aldehyde against **1.23** for
benzaldehyde, so at the same product *concentration* a 4OMe run sits at six times
the *absorbance*. Detector compression at high absorbance would bend a curve over
exactly as a sink does — `dA/dt` falling as `A` rises — and would do it to the
4OMe runs first. Two controls kill it.

- **Compare the substrates at matched absorbance, not matched concentration.**
  Over the window they share, 0.015–0.475 AU, the catalysed 4OMe curves give
  **−0.621 ± 0.099** on product and the catalysed BnOH curves **+0.386 ± 0.136**
  — a gap of 1.01 ± 0.17, **6.0σ**. The photometer does not know which alcohol
  is in the cuvette.
- **The enzyme-free 4OMe curves cover the same absorbances.** They run
  0.012–0.581 AU, against the catalysed set's 0.012–1.252, and their product
  coefficient is **+0.025 ± 0.048**. Same instrument, same wavelength, same
  substrate, same absorbance range, no product dependence. Compression cannot
  be selective for the presence of enzyme.

The pH-and-temperature row matters because `net` is bigger where the reaction is
faster, and the reaction is faster at high pH. Putting pH in the model moves the
product coefficient the *wrong* way for that objection, and pH's own coefficient
is +0.003 ± 0.157.

The last row of each block is the strongest form of the test. One dummy per
experiment absorbs temperature, pH, buffer, enzyme, run length, day and cell, so
the only thing left to carry the product term is **the substrate ladder inside
each run** — four cuvettes that differ in nothing but composition and sit in the
same cell block for the same number of seconds. On that comparison alone the
coefficient is **−0.919 ± 0.161**, indistinguishable from the −1 a rate that
falls as `1/(product)` would give.

**One bias cannot be removed, and it runs the safe way.** `net` is the integral
of the rate, so a curve that decelerates *less* makes *more* product at the same
starting rate; that pushes the product coefficient upwards, towards zero. Every
negative coefficient above is a floor.

## 3. What the product does: consumed, not blocking

Two mechanisms fit "the product sets the rate", and they put a straight line
through different transforms of the same two columns:

    A' = v - k A          the rate is linear in A          a SINK
    A' = v / (1 + A/Ki)   1/rate is linear in A            INHIBITION

Both are two-parameter lines on the same points, read off the rolling rate after
each curve's maximum, so the comparison is fair. Of the **29** catalysed 4OMe
phosphate curves whose decline is deep enough to have a shape (R² > 0.95),
**24 favour the rate, 0 favour the reciprocal** and 5 tie within 1%; the median
R² is **0.989** against **0.971**. The discrimination is carried by the curves
that fall furthest — over a 30% decline both transforms are nearly straight and
neither can win.

So the rate is linear in the product. That is production minus a **first-order
loss of the very thing being measured**, and it has one immediate consequence: a
stationary level `A∞ = v(S)/k`, which the run is heading for.

### The level it is heading for was predicted before it was measured

If the signal is a species the oxidant makes from the substrate and then
destroys, `A∞ = v(S)/k` must carry whatever substrate order the **rate** has —
not first order, because the production is partly saturated. The temperature
series measured that order on its own rates, **+0.577**, with nothing in this
document existing yet.

**The plateau's own order is +0.610 ± 0.067**, over a 37× range in [S] and 29
curves. Nothing was fitted to make those agree.

Where the two meet gives the one number a calculation can be checked against:

    k_A/k_S = [S] / [A]inf     median 54, IQR 42-81

the oxidant attacks the **aldehyde about fifty times faster than the alcohol**,
which at 25 °C is a barrier difference of about **9.9 kJ/mol**. It is an
**upper** bound, and the reason is worth stating: the plateau is read off the
catalytic *increment*, while the aldehyde in the cuvette is the increment plus
whatever the enzyme-free background made, and the background is not measured at
these compositions.

### There is not enough product to inhibit the catalyst

Blocking a fraction of the catalyst takes that fraction of `[enz]` **bound**,
whatever the binding constant — tightening the binding moves an equilibrium, it
does not create inhibitor.

| exp | T | product made | `[enz]` | product as a share of it | late/early |
|---|---|---|---|---|---|
| 19 | 15 °C | 5.5 µM | 241 µM | 2.3% | 2.32 |
| 18 | 20 °C | 5.0 µM | 273 µM | 1.8% | 2.02 |
| 14 | 25 °C | 24.5 µM | 273 µM | 9.0% | 1.11 |
| 17 | 30 °C | 11.8 µM | 273 µM | 4.3% | 1.43 |
| 15 | 35 °C | 31.9 µM | 273 µM | 11.7% | 0.88 |
| 16 | 40 °C | 39.8 µM | 241 µM | 16.5% | 0.84 |

At 40 °C the numbers are within a factor of about one of each other, so this
argument **narrows the field rather than closing it** — and the loophole is that
the background's aldehyde is in the cuvette and not in the signal. It bites
harder elsewhere: **exp 21 loses 40% of its rate on 36 µM of product against
241 µM of catalyst**, where blocking that fraction would take 97 µM bound.

## 4. And it is specific to this substrate

The BnOH blocks reach **0.386 mM** of product against this block's **0.214 mM**,
and their product coefficient is **positive** — more product, *less*
deceleration, which is the autocatalysis the in-scope analysis reports. On the
one-dummy-per-experiment comparison the two substrates differ by

    -0.919 +/- 0.161   against   +0.245 +/- 0.121     a gap of 1.16 +/- 0.20

**5.8σ.** Note what the in-scope BnOH block does instead: it decelerates with
**run length** (−0.697 ± 0.136) while accelerating with product. A clock and an
autocatalysis at once, which is a different problem from this one and is not
solved here.

### The chemistry that turns both signs

The peroxide adduct at the aldehyde — `MECHANISM.md`'s `C1` — has two exits, and
they load charge on different atoms:

    Cannizzaro (step 2)   hydride leaves C1 and adds to a second ArCHO
                          the ACCEPTOR's carbonyl carbon takes delta-minus
                          -> favoured by electron WITHDRAWAL, rho > 0

    Dakin                 the aryl group migrates from C to the peroxide O
                          the MIGRATING carbon takes delta-plus
                          -> favoured by electron DONATION, rho < 0 on sigma+

A 4-methoxy group therefore **shuts the first exit and opens the second**, and
the Dakin reaction's textbook requirement is exactly a *para* or *ortho* donor.
Benzaldehyde has none, so its adduct expels hydride and gives benzoic acid;
4-methoxybenzaldehyde's migrates and gives the formate and then 4-methoxyphenol,
which absorbs far less at 300 nm. One substituent, two opposite signs, and the
two substrates' curves bend in opposite directions.

**This is inference from the electronics, not observation.** No product was
identified in these runs and the archive holds one wavelength. It is also
already in `MECHANISM.md` as side reaction **S2**, which that document called
"the sink most likely to bite quantitatively" before any of this was measured.

The archive's own bound on the BnOH side is weak, and should be quoted as such:
pushing the same plateau arithmetic onto the in-scope BnOH curves — generously
attributing *all* of their deceleration to a sink, when it is mostly clock —
gives `k_A/k_S` ≈ **185** on the 3 of 110 curves that yield a plateau at all,
against 54 for 4OMe. The evidence for the contrast is the regression, not that
ratio.

## 5. What it does to the activation parameters

Every rate in the temperature series is already net of this loss, so none of
them is the production rate. Recovering it needs care, because the obvious
recipe is wrong: adding `k·A` back at the moment the rate peaks silently swaps
`v_peak` — an asymptote, and truncation-free by construction — for a value read
at the last reading on the four coldest runs, which reintroduces exactly the
cold-end truncation `v_peak` exists to avoid. That mistake moves the activation
energy by **+2.4 kJ/mol of pure bookkeeping**.

What is done instead is to fit the sink model itself, `A′ = v(1−e^(−t/τ)) − kA`,
with **k pinned** at the value `slowdown.sink_activation` gives for that
temperature. `v` is then a parameter of the model rather than a reconstruction,
and pinning k removes a degeneracy that makes the free four-parameter fit
useless here — a curve that has not turned over is fitted equally well by any
(v, k) holding v/k at the plateau, and three of these 24 come back with k a
hundredfold out.

| estimator | E<sub>a</sub> kJ/mol | ΔH‡ kJ/mol | ΔS‡ J/mol/K | ΔG‡(298) kJ/mol | Arrhenius rms |
|---|---|---|---|---|---|
| `v_peak`, published | 88.77 ± 1.77 | 86.28 | −57.7 ± 5.9 | 103.47 | 0.099 |
| `v_prod`, sink model with k pinned | 85.82 ± 3.40 | 83.32 | −67.2 ± 11.3 | 103.34 | 0.189 |

The production rate is **6.7% above** the published estimator on the median
curve, and between −6.6% and +12.7% depending on temperature **with no order to
it**. A factor that does not order in temperature cancels out of a slope, and
the activation energy duly moves by **−2.96 ± 3.83 kJ/mol** — nothing. The level
shift is worth **+0.54 J/mol/K** on ΔS‡ against a ±5.9 error.

**So the temperature series' section 4 stands, and now with a bound on this
rather than a hope.** The corrected route is the more nearly correct quantity
and the noisier one — Arrhenius scatter 0.189 against 0.099, because the sink
shape is one the four coldest runs cannot test — so quoting it as the headline
would trade a bias smaller than the error for a variance larger than it.

*What that rests on*: the sink's own activation energy, **72.3 ± 10.0 kJ/mol**,
fitted to three temperatures (25, 35 and 40 °C — the only ones whose curves turn
over far enough to give a slope) and extrapolated down to 15 °C.

## 6. What this settles, and what it does not

**Settled.** The fall is set by the product; it belongs to the catalysed pathway
and not to the solution, because the enzyme-free channel does something else
over the same product range; the rate is linear in the product rather than
hyperbolic in it; the stationary level carries the substrate order it has to;
and the effect is specific to 4OMe-BnOH at product concentrations where BnOH
shows nothing.

**Not settled, and not settleable from absorbance.** Whether the aldehyde is
*consumed* by the oxidant or merely *scavenges* it. Those are the same reaction
seen from two ends — an oxidant that attacks the electron-rich aldehyde is both
destroying the chromophore and being diverted from the alcohol — and a
single-wavelength trace cannot say which half is being watched. Both give
`A′ = v − kA`; both give `A∞ ∝ v(S)`.

Also not settled: **what the enzyme-free clock is.** `background_reaction` §7
left it open for BnOH ("something else decays during these runs — peroxide, or
the cell"), and this adds only that it is rate-independent and
product-independent, which rules product effects out of it and leaves the rest.

**Three things would finish it**, in order of cost:

1. **4-methoxyphenol's ε at 300 nm.** If the absorbance loss is the Dakin
   product replacing the aldehyde, the size of the loss is predictable from the
   two coefficients. Folded into `COMPUTATIONAL.md` C1, and the cheapest test
   here by a wide margin.
2. **C6** — ΔG‡ for the oxidant attacking ArCH₂OH against ArCHO, for Ar = phenyl
   and 4-methoxyphenyl. Gate: ≤ 9.9 kJ/mol for 4-methoxy, and enough smaller for
   phenyl that 0.386 mM of benzaldehyde does nothing.
3. **C5** — the substituent effect on the hydride-acceptor step, which is the
   other half of the same electronic argument and the half the BnOH
   autocatalysis rests on.
