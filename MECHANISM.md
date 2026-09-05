# Proposed reaction mechanism

This document lays out the mechanistic hypothesis for the chemzyme-catalyzed oxidation
of benzylic alcohols (benzyl alcohol, 4-methoxybenzyl alcohol) by H2O2, the reasoning
and literature support behind each step, and the open questions that still need
resolving before/while fitting a kinetic (ODE) model to the progress-curve data. See
`DATA_VERIFICATION.md` for the data-cleaning side of this work, and
`COMPUTATIONAL.md` for quantum-chemistry tasks that would settle questions
the literature and the dataset cannot.

[Structural analysis and ODE reduction](#structural-analysis-and-ode-reduction)
reduces the mechanism to a 3-ODE / 4-parameter system and records the tests it
has been put to against the progress curves — read that section before fitting.

**That reduction has since been implemented and fitted, and it does not fit.**
See `FITTING.md`: the model is first order in `[S]` by construction while the
data is roughly half order, the `E0 = 0` limit as written below is a fixed point
that never starts, and `signal = [A] + r[BA]` can only produce an induction
period if `r > 1`. The chemistry below is not yet convicted — the deficiencies
found so far sit in the reduction and the observation equation — but nothing in
this document should be fitted without reading that one first.

The catalyst ("chemzyme") is a cyclodextrin scaffold bearing a ketone at its active
site — this is the artificial-enzyme platform developed by the Bols group (see
[References](#references), items 1–9), reported to accelerate benzylic alcohol
oxidation by H2O2 up to ~60,000-fold over background.

## Species

| Symbol | Species |
|---|---|
| K | active-site ketone (catalyst), the ACTIVE form |
| Kh | the ketone's gem-diol **hydrate** — proposed resting form, added 2026-09-05 |
| Kh⁻ | its conjugate base, the gem-**diolate** — proposed off-path sink |
| KP | ketone + H2O2 adduct ("perhydrate", a gem-diol hydroperoxide) |
| KD | active-site dioxirane |
| S | substrate benzylic alcohol |
| A | benzaldehyde (the measured product) |
| HOO⁻ | hydroperoxide anion (conjugate base of H2O2, pKa ≈ 11.6–11.75) |
| C1 | Cannizzaro-type tetrahedral adduct of A + HOO⁻ |
| PBA | perbenzoic acid (peroxybenzoic acid) |
| BA | benzoic acid |

## The mechanism

**Catalyst-independent loop** (operates with or without the chemzyme present —
proposed to explain why enzyme-free control experiments also show accelerating,
autocatalytic kinetics):

1. `A + HOO⁻ ⇌ C1` — fast addition of hydroperoxide anion to the aldehyde,
   forming a Cannizzaro-type tetrahedral intermediate.
2. `C1 + A → PBA + S` — rate-determining step. C1's aldehydic C–H transfers as a
   hydride to a second molecule of A (the classic Cannizzaro disproportionation
   mechanism, just with HOO⁻ instead of OH⁻ as the initiating nucleophile). C1 is
   oxidized to perbenzoic acid; the second A is reduced back to the alcohol S. Net:
   `2A + HOO⁻ → PBA + S`.
3. `PBA + S → A + BA (+ H2O)` — slow, uncatalyzed direct oxidation of the alcohol
   by perbenzoic acid (peracids are known to oxidize alcohols, just more slowly
   than a dioxirane does). This closes the catalyst-free loop (1 → 2 → 3): each
   turn nets one BA and regenerates the A needed to feed the next turn.

**Catalyst-dependent loop:**

0. `Kh ⇌ K + H2O` — **the induction.** The catalyst is not in its active form
   when the run starts and converts into it on its own clock, before any
   turnover. *(Added 2026-09-05. The STEP is required by the data; its IDENTITY
   is not. What is measured is a unimolecular activation of the catalyst over a
   covalent barrier; that it is a hydrate dehydrating is the candidate that fits
   the dependences. `induction/ANALYSIS.md` §7.)*

   What the archive measures:

   - it **needs the catalyst** — 10 of 49 enzyme-free 4OMe curves show any lag
     at all, median depth 0.000, against 109 of 147 catalysed;
   - it is **not waiting for product** — the clock regressed on the curve's own
     rate gives `+0.245 ± 0.127` over 147 curves where a product threshold
     needs `−1`, and the BnOH block agrees at `+0.897 ± 0.178`;
   - its barrier is **77 ± 12 kJ/mol** over six temperatures at a window all
     six runs share, three to four times what dissolution, diffusion or thermal
     equilibration cost;
   - it is **slowed by base**: `d ln τ/d pH = +0.12 to +0.34`, positive on all
     four of the archive's pH ladders, in three buffers and on both substrates;
   - it is **hurried by the substrate** where the buffer is held fixed:
     `d ln τ/d ln[S] = −0.18 to −0.71` on exps 135–151, the only block in which
     `[S]` moves without `[buf]` moving against it.

   The last two are BOUNDED, and the bounds are what make them mechanistic. For
   a species X that holds the catalyst OFF the activation path, `1/τ =
   k_act/(1 + K[X])` gives `d ln τ/d ln[X] ∈ (0, +1)`; for a species Y whose
   bound form activates, `d ln τ/d ln[Y] ∈ (−1, 0)`. Neither bound contains a
   rate constant, so a coefficient outside it falsifies the scheme rather than
   fitting it. The pH result puts X at **5–15 % saturated** across the archive's
   whole pH range; the substrate result puts Y at **18–72 %**.

   **A resting anionic adduct at the carbonyl, displaced by the alcohol binding
   in the cavity, is the shape that fits.** Ketone hydration is fast and the
   hydrate of an electron-poor ketone is acidic; a gem-diolate cannot dehydrate,
   so more base means more of it and a longer wait — the sign measured. A
   hydrophobic guest in the cyclodextrin cavity shifts hydration towards the
   free ketone — the substrate sign. **`COMPUTATIONAL.md` C7 is the test**: a
   hydrate pK<sub>a</sub> near 10 puts the diolate exactly where the boric
   ladder (pH 8.46–10.34) sits, and one near 13 rules it out and leaves the
   perhydrate trap of C8 as the only survivor.

   **What is NOT established.** That the clock is independent of `[enz]`, which
   is what "unimolecular" means and the one prediction here that a concentration
   could falsify. The archive holds two `[enz]` pairs: one has no lag to time,
   and the other gives a 1.8× *longer* clock at 2.4× the catalyst — past the
   replicate floor, in the wrong direction, on two runs in the block whose
   signal control fails. `induction/ANALYSIS.md` §7g.
4. `K + H2O2 ⇌ KP` — the ketone forms its gem-diol hydroperoxide ("perhydrate")
   with H2O2.
5. `KP + S → K + A` — the perhydrate directly oxidizes the alcohol substrate to
   the aldehyde. Hypothesized as the slow, non-autocatalytic "seed" step that
   produces the first trace of A needed to start the autocatalytic loop, before
   enough product has accumulated for steps 1–2 to run.

   *2026-09-02, extended 2026-09-05.* The induction period of the catalysed
   4OMe curves is **not** this step running to its threshold: it ends on a clock
   rather than at a product level, it has no substrate order, and it is absent
   from every enzyme-free 4OMe curve in the archive. **That last clause is a
   4OMe clause.** 14 of 26 enzyme-free *BnOH* curves do begin below their
   eventual rate, at a median depth of 0.138 — in exps 3 and 6, three and five
   hour runs at pH 6.71, and in exp 65, whose four cuvettes share the break
   `scope.synchronised_break` measures. That is steps 1–3 above, which are
   autocatalytic in the product and need no catalyst, and which S2 already finds
   switched on for benzaldehyde and off for the 4-methoxy aldehyde. **A lag is
   two different things on the two substrates and they must not be pooled.** Whatever v5 does, something slower and
   unimolecular has to happen to the catalyst first. Step 4's own status is
   affected too and in the opposite direction to the obvious one — see
   `COMPUTATIONAL.md` C8 and `induction/ANALYSIS.md` §4a.
6. `PBA + K → KD + BA` — the ketone catalyst reacts with the perbenzoic acid
   produced in step 2 (not with H2O2/Oxone directly), forming a Criegee-type
   adduct that collapses by expelling benzoic acid as the leaving group and
   leaving a dioxirane at the former ketone carbon.
7. `KD + S → K + A` — the dioxirane oxidizes the substrate's benzylic C–H bond,
   regenerating free ketone catalyst and producing a fresh molecule of A. Closes
   the catalytic cycle.

Total catalyst is conserved: `[K] + [KP] + [KD] = [enz]0` throughout (no rate
constant needed to enforce this — it falls out of steps 4–7 by construction).
**With step 0 the law is `[Kh] + [Kh⁻] + [K] + [KP] + [KD] = [enz]0`**, and the
reduction below has not been redone under it: every fit in `FITTING.md` starts
from a fully active catalyst at `t = 0`, which is the one thing the induction
says is false. A known gap, listed in the open questions.

## Structural analysis and ODE reduction

*(Added 2026-08-30.)* Before fitting, the mechanism above was analysed as a
stoichiometric system. The reduction is worth doing for its own sake — it cuts
the fitting problem from 9 states and 9 rate constants to 3 ODEs and 4 rate
constants — but it also forced a correction to how the measured signal is
interpreted. Both are recorded here.

### Conservation laws (exact, no approximation)

The stoichiometric matrix of the 7 steps over the 9 species has rank **6, not
7**, because step 3 is exactly the sum of steps 6 and 7:

```
R3  =  R6 + R7        (PBA + S -> A + BA)  =  (PBA + K -> KD + BA) + (KD + S -> K + A)
```

i.e. the catalyst-free peracid oxidation and the catalysed dioxirane route have
**identical net stoichiometry**; the catalyst changes the rate, not the
bookkeeping. Three independent conservation laws follow:

| | |
|---|---|
| catalyst | `[K] + [KP] + [KD] = E0` |
| aryl | `[S] + [A] + [C1] + [PBA] + [BA] = S0` |
| oxidizing equivalents | `[H2O2] + [KP] + [KD] + 2[C1] + [A] + 2[BA] + 3[PBA] = const` |

The last uses the assignment S = 0, A = +1, BA = +2, PBA = +3, H2O2 = KP = KD =
+1, C1 = +2 (one "equivalent" = one O atom of oxidizing power); every step
balances under it. These three laws remove three ODEs at zero cost in
assumptions: **9 -> 6 states.**

### Four-stage reduction

| stage | assumption | states | rate constants |
|---|---|---|---|
| 0 | full mass-action system | 9 | 9 |
| 1 | the three conservation laws above | 6 | 9 |
| 2 | C1 pre-equilibrium (step 1 fast, step 2 rate-determining) | 5 | 7 |
| 3 | catalyst QSSA (KP fast equilibrium, KD reactive intermediate) | 3 | 6 |
| 4 | both saturation terms negligible | **3** | **4** |

**Stage 2** is the mechanism's own claim and the textbook Cannizzaro treatment.
Setting `[C1] = K1[A][HOO-]` collapses `k1f, k1r, k2` into a single constant:

```
v_can = k_can * [A]^2 * p(pH) * [H2O2]      k_can = k2*k1f/k1r,  p = Ka/(Ka + 10^-pH)
```

which is literally the classical Cannizzaro rate law `k[ArCHO]^2[OH-]` with HOO-
substituted for OH- — so the lumped constant is the quantity the literature
reports anyway, and the three microscopic constants were never separately
identifiable from progress curves.

**Stage 3** makes all three catalyst states algebraic:

```
[KD] = k6[PBA][K] / (k7[S])
[K]  = E0 / (1 + K4[H2O2] + (k6/k7)*[PBA]/[S])
```

E0 runs 0.014–0.27 mM against tens of mM H2O2, so sequestration of peroxide into
KP is negligible and free H2O2 ~ total H2O2.

**Stage 4** drops `(k6/k7)[PBA]/[S]` (a trace intermediate over a mM substrate)
and, if `K4[H2O2] << 1`, sets `[K] = E0`, merging `k5` and `K4` into `k5'`.

### The reduced system

```
v_can = k_can [A]^2 p [H2O2]        v3 = k3 [PBA][S]
v5    = k5' E0 [H2O2][S]            v6 = k6 E0 [PBA]   (= v7)

dA/dt   = -2 v_can + v3 + v5 + v6
dPBA/dt =   +v_can - v3      - v6
dS/dt   =   +v_can - v3 - v5 - v6

[H2O2] ~ const,   [BA] = S0 - S - A - PBA
```

Three ODEs, four rate constants (`k_can, k3, k5', k6`), **exactly linear in
E0**, and non-stiff — so `solve_ivp` with RK45 rather than Radau.

The identifiability payoff: at **E0 = 0** this collapses to 2 ODEs and 2
parameters (`k_can, k3`). Fit those on the enzyme-free controls, freeze them,
and only `k5'` and `k6` remain for the catalysed runs. That sequential strategy
is the difference between a well-posed fit and a hopeless one.

### The exact identity, and what it says about the observable

Summing the A and PBA balances, everything cancels except two terms. With no
approximations whatsoever (keeping C1 and KD explicit):

```
d( [A] + [C1] + [PBA] + [KD] ) / dt  =  v5 - v_can
```

Two consequences, both structural rather than parameter-dependent:

**1. Step 5 is the only net source of aldehyde.** Summing steps 1+2+3 gives
`A + HOO- -> BA + H2O`; summing 1+2+6+7 gives the same thing. Step 2 destroys
two aldehydes per peracid and steps 3/7 return only one. So within this
mechanism A is an *intermediate* on the route `S -> A -> BA`, not an
autocatalyst.

**2. The mechanism IS autocatalytic — in benzoic acid, not in benzaldehyde.**
The terminal sink obeys `dBA/dt = (k3[S] + k6*E0)*[PBA]`, which starts at
exactly zero (PBA = 0 at t = 0) and accelerates as A accumulates: a genuine
sigmoid whose lag is set by how fast the seed step v5 makes the first aldehyde.
The aldehyde pool does the opposite — since `dPBA/dt >= 0` during any growth
phase, `dA/dt <= v5 - v_can <= v5(0)`, and v5 only ever decreases as S and H2O2
deplete. **[A] must therefore have its steepest point at t = 0**: concave
throughout, no lag, no upward inflection.

So the mechanism predicts an accelerating total-product curve and a
decelerating aldehyde curve. Which one the absorbance tracks decides whether it
fits the data.

### Empirical test: the data shows the accelerating curve

Position of maximum slope, all usable curves, known-inverted experiments (50,
58, 77–79, 84, 85) excluded, curves smoothed over ~5% of the run before
differentiating (n = 326):

| position of peak slope | curves | |
|---|---|---|
| < 5% into the run | 105 | 32% |
| 5–15% | 52 | 16% |
| 15–50% | 67 | 21% |
| past halfway | 102 | 31% |

**52% of curves reach peak slope more than 15% into the run**, 129 of them at
1.5–5x the initial slope and many with an initial slope of essentially zero.
That is unambiguous autocatalysis with a real induction period, and it is
incompatible with the bound `dA/dt <= v5(0)`.

A separate check for the predicted rise-and-fall of the aldehyde pool found
**381 curves that rise and plateau with no turnover** against only **30 with a
genuine interior maximum followed by decay** (clustered in experiments 128–151);
28 further apparent turnovers are the already-excluded inverted curves.

**Conclusion: the aldehyde reading of the absorbance is falsified, not the
mechanism.** The chemistry above survives intact; what does not survive is the
assumption that the observable is `[A]`.

### The extinction coefficients support the reinterpretation

*(Updated 2026-08-30 — the wavelength is no longer an inference. Each `.xls`
sheet declares it in its header: **285 nm for BnOH** with `e = 1.23`, **300 nm
for 4OMe-BnOH** with `e = 7.53`. Literature puts benzaldehyde in water at
**eps ~ 1400 M^-1 cm^-1 at 278-279 nm**, so `e = 1230 M^-1 cm^-1` at 285 nm sits
on the falling edge of that same n->pi* band: **`e` is benzaldehyde's own
extinction coefficient, not a differential one**, and the abs -> [P] conversion
does assume the aldehyde is the sole absorber. What is still missing is
eps(benzoate) at 285 nm; benzoate's weak band lies near 268 nm and no reliable
aqueous value at 285 nm could be retrieved (the two best sources are paywalled).
A band-shape bracket gives eps ~ 100-400 M^-1 cm^-1, i.e. **r ~ 0.08-0.33**,
lower than the ~0.5 guessed below - though the signal contribution is
`r x [BA]/[A]`, and A is an intermediate whose pool stays small while BA
accumulates, so a small r does not by itself restore the aldehyde reading.
One UV spectrum of benzoic acid at working pH settles it. See
`DATA_VERIFICATION.md` 2026-08-30.)*

`experiment_data.csv` carries exactly one `e` per substrate: **1.23 mM^-1 cm^-1
for BnOH** and **7.53 mM^-1 cm^-1 for 4OMe-BnOH**. **Confirmed: the BnOH
experiments are monitored at 285 nm.** Benzaldehyde's strong pi->pi* band at
250 nm has eps ~ 12,000–14,000 M^-1 cm^-1 — about ten times higher than 1230 —
so the assay is sitting on the weak n->pi* band, and 285 nm is consistent with
that band's tail. There, benzoic acid (eps ~ 800 M^-1 cm^-1 at 273 nm) absorbs
*comparably*, not negligibly. The abs -> [P] conversion in the notebook
assumes the aldehyde is the sole absorber; it very likely is not. *(Both
literature eps values are from memory and should be verified against a
spectrum at 285 nm specifically before citing — the 273 nm BA value in
particular should be re-checked at 285 nm, since eps can shift meaningfully
over a 12 nm span on a band tail.)*

### Consequence for fitting: one extra parameter in the observation equation

No change to the chemistry is required — only to how the model is compared to
the data:

```
signal(t) = [A] + r * [BA]          r = eps_BA / eps_A,  fitted or measured
```

`r -> 0` recovers the pure-aldehyde reading, so **the fit adjudicates the
question**. If `r` converges on a physically sensible value (~0.5 at 280 nm),
that is independent corroboration of the whole picture. Final model: **3 ODEs,
4 rate constants + 1 chromophore ratio.**

### The observation equation does not rescue the lag

*(Added 2026-08-31, when the reduced system was implemented as
`data/kinetic_model.py`. Two structural results came out of writing it down;
neither depends on any fitted number, and both are pinned as tests in
`data/test_kinetic_model.py`.)*

**1. The enzyme-free limit needs a third constant, not two.** The reduction
above says that at `E0 = 0` the system "collapses to 2 ODEs and 2 parameters
(`k_can, k3`)". That is true of the algebra and false of the trajectory: with
`A = PBA = 0` at `t = 0`, both `v_can` (which goes as `[A]^2`) and `v3` (which
goes as `[PBA]`) are zero, so the system sits at a fixed point and never starts.
Seeding it with a trace of aldehyde does not help either — steps 1–2 destroy two
aldehydes per peracid while step 3 returns one, so the catalyst-free loop is a
net aldehyde **sink** and the trace decays. This is the same fact stated above
as "step 5 is the only net source of aldehyde", followed to its conclusion.

Since the enzyme-free controls demonstrably do react, an **E0-independent
source** is required: the uncatalysed direct oxidation `S + H2O2 -> A`, i.e. the
`E0 -> 0` limit of step 5. Writing the seed as `(k0 + k5' E0)[H2O2][S]` keeps the
model exactly linear in E0 as the reduction promises, with an intercept rather
than through the origin. Stage 1 therefore fits **three** rate constants plus
`r`, not two.

**2. The observable can only accelerate if `r > 1`.** Differentiating
`signal = [A] + r[BA]` along the reduced system, with `W = v3 + v6`:

```
d(signal)/dt = v_seed + (1 + r) W - 2 v_can
dPBA/dt      = v_can - W  >=  0        while peracid accumulates
=>  d(signal)/dt  <=  v_seed + (r - 1) v_can  <=  v_seed(0)     for r <= 1
```

and `v_seed = (k0 + k5' E0)[H2O2][S]` only ever decreases, since `S` and `H2O2`
only deplete. So for any `r <= 1` the signal's steepest point is at `t = 0`:
**concave throughout, no lag, no upward inflection** — precisely the bound
derived above for `[A]` alone. The observation equation inherits the
falsification rather than escaping it.

A random search over 227 parameter sets spanning eight orders of magnitude in
every rate constant, at both `E0 = 0` and `E0 > 0`, produced **not one**
accelerating curve at `r <= 1`. At `r > 1` acceleration appears immediately, up
to 4.7x the initial slope.

*(The 52% figure in the table above was re-measured on 2026-08-31 over the 402
curves the fitting code selects, by the same smoothed method, and comes out at
**37.6%** — 151/402. The selections differ: the n = 326 above predates the
carbonate rule, the exclusions of exps 50, 64 and 85, and the cuvette exclusions
of 25,2 and 25,4. The figure was 136/402 until 2026-08-31, when the readings
moved to the instrument's own .rre files; the export's 0.001 AU rounding had
been flattening fifteen lags below the threshold.)*

**And the first fit shows the error runs the other way.** Fitted on
BnOH/25 C/phosphate, the model at its best-fit `r = 1.52` lags in 19 of 23
enzyme-free curves and 19 of 20 catalysed ones, where that block's data lags in
7 of 43. So `r <= 1` produces no lag anywhere and the fitted `r > 1` produces one
almost everywhere, with the data in between. See `DATA_VERIFICATION.md`
2026-08-31.

**What this costs the mechanism.** Against the measured 52% of curves that reach
peak slope more than 15% into the run, the model requires
`eps(benzoate) > eps(benzaldehyde)` at 285/300 nm. The band-shape bracket
recorded above puts `r ~ 0.08–0.33`, and benzaldehyde's n->pi* band is the
stronger of the two by every account — so `r > 1` is not a value the
spectroscopy will support. The three ways out, in order of how cheaply they can
be tested:

- **The chromophore is neither A nor BA.** Perbenzoic acid is the one species in
  the mechanism whose concentration genuinely accelerates, and no extinction
  coefficient for it at 285 nm has been sought. Adding `r_PBA [PBA]` to the
  observation equation is a one-line change to `kinetic_model.observable`.
- **A step is missing.** Any route in which the product catalyses its own
  formation — rather than merely accumulating — would lift the bound. The
  radical/O2 chain of S3 is the standing candidate and is not in the model.
- **The reduction is too aggressive.** Stage 2's C1 pre-equilibrium removes the
  one species whose build-up could itself produce a lag. Restoring C1 explicitly
  costs two rate constants and is worth trying before abandoning the mechanism.

This is a sharper statement of the open question below than the one it replaces:
it is no longer "what does the absorbance measure?" but "**the absorbance cannot
be a non-decreasing combination of A and BA alone**".

### Lag-scaling test of the seed step

Within this mechanism the induction period is ended by v5, so lag length should
shorten as E0, [H2O2] and [S] rise. Measuring lag as the time of maximum slope,
**within each experiment** (so pH, T, buffer identity and run duration are all
held fixed):

| condition | experiments varying it | mean rho(lag, x) | verdict |
|---|---|---|---|
| `[S]` | 58 | **-0.26** (40/58 negative, sign-test p = 0.005) | correct sign, significant |
| `[H2O2]` | 22 | +0.31 (6/22 negative, p = 0.053) | marginal, opposite sign |
| `E0` | **1** | — | untestable |

Reading these:

- **The `[S]` result supports the seed step**, but `[S]` is collinear with
  `[buf]` in every titration series (see `DATA_VERIFICATION.md`), so it could be
  a buffer effect wearing substrate's clothes. Consistent, not confirmatory.
- **The `[H2O2]` sign is not a falsification.** H2O2 accelerates both v5 (which
  seeds A) and v_can (which drains it), while v3 and v6 — the steps that
  actually generate the observable BA — carry no H2O2 dependence at all. The
  mechanism makes no clean sign prediction here; this is a quantitative
  constraint for the fit to satisfy, not a test it failed.
- **E0 is untestable from this dataset**: only *one* experiment varies enzyme
  concentration across its own samples. A cross-experiment correlation gives
  rho = +0.19 (p = 0.018, wrong sign), but pH (rho = -0.43, p = 5e-9) and T
  (rho = -0.31, p = 5e-5) also vary between experiments, so that number is pure
  confounding and should not be used.

**2026-09-02 — the sign survives and the magnitude does not.** `induction/`
re-ran this question with an induction statistic that is defined for every curve
rather than only where a landmark exists, and with the two channels separated.
The `[S]` correlation is still weakly negative — the induction time's substrate
order is **−0.121 ± 0.148** over the 147 live catalysed 4OMe curves — but that is
the wrong *size* for the seed step. If the induction ended when v5 had made
enough A, the induction time would carry the whole of the rate's substrate
order with the sign reversed, and that order is measured on the same curves at
**+0.471 ± 0.084**. Regressing the induction time directly on each curve's own
rate gives **−0.025 ± 0.109** where a product threshold requires −1.

So this row should no longer be read as support for the seed step. Two further
results decide it. **The induction does not exist without the catalyst at all**
— 0 of 49 enzyme-free 4OMe curves have one, at matched substrate, buffer,
peroxide, pH and temperature, including a 17934 s run — which is fatal to any
account in which the induction is the catalyst-free loop of steps 1–3 waiting
for its first aldehyde. And **its barrier is a covalent one**: **77 ± 12 kJ/mol** over six
temperatures at a window all six runs share, three to four times too large for
anything physical. What the induction times is a unimolecular change on the
catalyst, before turnover begins — **step 0 above**. See
`induction/ANALYSIS.md` §3 and §7f, and `COMPUTATIONAL.md` C7.

*(The barrier read 95 ± 16 until 2026-09-05. That figure is
`arrhenius.inverse_tau`, the ONE-phase fit's τ, which reaches only 15–30 °C
because above 32 °C a decelerating curve gives the one-phase form no lag to
find. The two agree inside their errors; the six-temperature number is the one
to quote, and it carries a window systematic of about ±7.)*

## Step-by-step reasoning and literature support

**Steps 1–2 (Cannizzaro-type disproportionation of the product aldehyde, HOO⁻ as
nucleophile).** *This is a working hypothesis of this thesis, not a
literature-grounded step — state it as such.* A dedicated deep literature search
found **no precedent anywhere** for an aromatic aldehyde + HOO⁻/H2O2
disproportionating via a bimolecular Cannizzaro-style hydride transfer to give
peracid + alcohol, nor any measured reaction order in [aldehyde] for peracid
formation from an aromatic aldehyde + H2O2. Every classical/industrial route to
perbenzoic acid found instead starts from an acyl derivative (benzoyl
peroxide/benzoyl chloride routes — items 25–26) or uses radical O2 autoxidation
of the aldehyde, never base-catalyzed aldehyde + H2O2 disproportionation.

What *is* solid is the analogy's baseline: the classical hydroxide-mediated
Cannizzaro reaction genuinely does run rate = k[ArCHO]²[OH⁻], and its
rate-determining step is confirmed by isotope labelling and solvent isotope
effects (k(D2O)/k(H2O) = 1.9) to be a real intermolecular hydride transfer from
the tetrahedral gem-diolate to a *second, separate* aldehyde molecule (Swain,
Powell, Sheppard & Morgan, item 27; earlier kinetics by Alexander, item 28).
Benzaldehyde is the textbook substrate. So bimolecularity itself is not the
issue — the open question is only whether the *HOO⁻* adduct does the same
chemistry as the OH⁻ adduct.

Two arguments for caution, both worth stating in the write-up:

- The classical Cannizzaro's driving force is formation of a fully
  resonance-stabilized carboxylate. The analogous product here, a
  percarboxylate (ArC(=O)OO⁻), does not delocalize symmetrically in the same
  way, so the thermodynamic push for the peroxide variant is weaker. (Reasoning,
  not literature-sourced.) The thesis author's position — that this makes the
  reaction *slower* than the OH⁻ analog but does not rule it out — is
  reasonable.
- The unimolecular alternative is well characterized for essentially this same
  adduct. Ogata & Sawaki (item 29) measured pH-dependent product profiles and
  migration isotope effects (kH/kD = 1.4–3.0 for H-migration, ≈1.0 for aryl
  migration) for benzaldehyde + perbenzoic acid, the signature of unimolecular
  Criegee-adduct collapse via migration to the peroxidic oxygen. For an
  aldehyde-derived adduct, H has the highest migratory aptitude in the
  Baeyer-Villiger series, so collapse straight to the ordinary carboxylic acid
  is the default expectation. This is a **competing branch at C1**, not a
  refutation — see [Competing side reactions](#competing-side-reactions).

Note: 4-methoxybenzaldehyde (the oxidation product of the other substrate,
4-methoxybenzyl alcohol) carries the para-alkoxy group that activates the
*aryl*-migration (Dakin) pathway, unlike plain benzaldehyde — so the two
substrates may not share identical background kinetics.

**The direct experimental test**, if bench time ever allows: measure the initial
rate of peracid appearance as a function of [ArCHO] at fixed [H2O2] and pH.
Second-order-in-aldehyde behaviour would be the first genuine direct evidence
for this step that the literature currently lacks entirely.

**Step 3 (PBA + S → A + BA, uncatalyzed peracid oxidation of the alcohol).**
Peracids are documented to oxidize alcohols directly, more slowly than a
dioxirane. This step is needed to close the catalyst-free version of the loop and
explain autocatalysis in enzyme-free controls — without it, steps 1–2 alone only
consume A to regenerate S, they never produce net new A.

**Step 4 (K + H2O2 ⇌ KP, ketone perhydrate formation).** Well precedented for
simple ketones (cyclohexanone forms the analogous 1-hydroperoxycyclohexanol
adduct reversibly — item 14). Cyclodextrin-cavity encapsulation is independently
documented to shift ketone/gem-diol-type equilibria substantially (item 15),
supporting that this equilibrium is real and cavity-perturbed for the chemzyme
specifically.

**Step 5 (KP + S → K + A, the "seed" step).** The weakest-precedented step,
mechanistically. The closest literature analogy found is metal-peroxo hydride
abstraction from a benzylic alcohol by a rhenium-diperoxo species (item 16,
k ~ 10⁻⁴–10⁻⁵ M⁻¹s⁻¹) — this establishes that a peroxide-bearing center
oxidizing a benzylic alcohol by formal hydride transfer is mechanistically real
and slow, giving plausibility and an order-of-magnitude rate anchor, but no
literature was found describing an organic ketone-perhydrate itself acting as
the hydride acceptor for a separate substrate molecule. Treat this rate constant
as free/unconstrained.

**Step 6 (PBA + K → KD + BA, peracid converts the ketone to its active
dioxirane).** Still the step needing the most support, but materially better
founded after a third targeted literature pass — *and* now with a specific
counter-experiment that has to be addressed head-on.

*The mechanistic picture.* Dioxirane formation from a Criegee adduct is not a
migration at all: it is a **3-exo-tet intramolecular displacement**, where the
adduct's own alkoxide oxygen attacks the distal peroxide oxygen and expels the
leaving group (bisulfate for Oxone, carboxylate here). Because that needs an
**anionic** peroxide species, the dioxirane channel is favoured at high pH while
the competing Baeyer-Villiger migration is favoured at low pH — so the branching
should be **pH-controlled**. This is the mechanism drawn in essentially every
secondary source and is consistent with all the primary evidence found, but see
the caveat below.

*Positive evidence:*

- **Item 32 (Schulz, Liebsch, Kluge & Adam, *J. Org. Chem.* 1997) is the key
  citation.** An organic **peroxy acid** (arenesulfonyl, generated in situ from
  an ArSO2-imidazole + H2O2) plus a ketone plus **NaOH** generates the
  corresponding dioxirane, established by **18-O labelling**, with the dioxirane
  pathway "virtually the exclusive one." This refutes the earlier conclusion
  that only peroxymonosulfate can convert a ketone to a dioxirane, and it needs
  base, exactly as the pH picture predicts. Caveat: the leaving group there is a
  **sulfonate**, not a carboxylate — that extrapolation remains unproven.
- **Item 33 (Porter, Yin & Pratt, *J. Am. Chem. Soc.* 2000)** is titled "The
  Peroxy Acid Dioxirane Equilibrium: Base-Promoted Exchange of Peroxy Acid
  Oxygens" — i.e. a peroxy acid and a dioxirane are a published, named,
  **base-promoted** equilibrium pair. Structurally the same elementary step, in
  its intramolecular form. Its SI reportedly contains computed transition-state
  structures.
- **Items 34–35 (Lange & Brauer 1996; Lange, Hild & Brauer 1999)** are the only
  kinetic studies of dioxirane *formation* (as opposed to dioxirane-mediated
  oxidation): rate = k[ketone][HSO5-]·F with F = Kw/([H+]+Ka2), for eleven
  ketones. The 1999 paper shows that for strongly electron-poor ketones the
  keto/gem-diol equilibrium **and the gem-diol's first ionization** must be
  included in the kinetic model — the closest literature precedent for treating
  the adduct's ionization quantitatively.
- **Item 14 (Rozhko et al. 2015)**: DFT finds dioxirane formation proceeds via
  **hydroperoxide-anion** addition to the ketone, and that the dioxirane sits in
  **equilibrium with the ketone** — notably with no good leaving group at all.

*The counter-experiment that must be addressed:*

- **Item 36 (Doering & Dorfman, *J. Am. Chem. Soc.* 1953)** ran 18-O labelling on
  **benzophenone + perbenzoic acid** — essentially this exact reaction — and
  **excluded both the dioxirane and the carbonyl-oxide routes**, establishing
  the Criegee/Baeyer-Villiger mechanism instead. Item 37 (Armstrong et al.
  1994/1996) reached the same conclusion by 18-O in a second peracid+ketone
  system.
- The available defence is that both were run in aprotic solvent **with no
  base** — precisely the regime in which the pH picture predicts BV should win.
  That defence is legitimate but **circular unless supported independently**,
  since the pH hypothesis is the thing being tested. State this openly in the
  write-up rather than omitting the 1953 result.

*Two things that weakened since the last pass:*

- The Kokotos precedent (item 11) is thinner than it looked. Their Criegee
  adduct does **not** collapse to the dioxirane on its own — it needs the Payne
  intermediate (from MeCN + H2O2 at pH 11), and the peracid is *regenerated*
  rather than expelled as carboxylate. Their steep pH dependence (98% at pH 11,
  41% at pH 10) is attributed to Payne-intermediate generation, not Criegee
  deprotonation. The same group later (item 38, *Green Chem.* 2025) explicitly
  **ruled out** a peracid intermediate in a closely related system.
- **The Bols group does not propose this step.** Item 4 turned out to be
  retrievable after all (hybrid open access), and its stated mechanism is only
  "covalent binding of H2O2 to the carbonyl group forming an activated
  cyclodextrin (CD*) followed by oxidation of the nearby amine or alcohol" — no
  peracid, no dioxirane, run at pH 7 phosphate. So steps 1–2 and 6 are this
  thesis's own addition to the group's published position, not an elaboration
  of it.

*Leaving-group penalty (reasoning, not literature).* Benzoate's conjugate acid
(pKa 4.2) is ~2 pKa units weaker-leaving than bisulfate (pKa2 ~ 2.0) and 6–10
units weaker than sulfonate. With a typical Bronsted beta_lg of -0.5 to -1.0 for
O-O displacement, the peracid case might run 10^2–10^4 slower than the Oxone case
— a real penalty, but potentially offset by high effective molarity in a
preorganized active site. Working the other way: the acyl group that makes
benzoate a poor leaving group also makes the peroxyester a *better* substrate
for anti-periplanar BV migration. No experimental or computational study
comparing dioxirane-forming barriers across leaving groups appears to exist.

**Step 7 (KD + S → K + A, dioxirane oxidizes the substrate).** Well precedented
by analogy to DMDO (dimethyldioxirane) chemistry: second-order kinetics
(rate = k[ROH][dioxirane]) for benzylic-alcohol oxidation by dioxiranes are
established, with KIE evidence (primary KIE ≈ 5.2 at the α-C–H) supporting a
concerted or oxygen-rebound-type C–H oxidation mechanism (items 18–20).

## Competing side reactions

These are not part of the proposed catalytic cycle, but the literature says they
are real and they compete for the same intermediates. They belong in the ODE
model as explicit sink terms rather than being argued about in the abstract — if
the fitted branching ratios come out physically absurd, *that* is the evidence
against the mechanism.

- **S1. `C1 → BA + OH⁻`** (unimolecular collapse of the Cannizzaro adduct
  straight to benzoic acid, via H-migration — item 29). Competes directly with
  step 2 for C1, and is the literature-default expectation. Consumes A without
  producing any PBA. For the 4-methoxy substrate, an aryl-migration (Dakin)
  variant giving 4-methoxyphenol competes as well.
- **S2. `PBA + A → 2 BA`** (Baeyer-Villiger oxidation of the aldehyde by the
  peracid — item 29 studies exactly this reaction). Competes with step 6 for
  PBA. **This is the sink most likely to bite quantitatively:** `[enz]` runs
  0.014–0.27 mM across the dataset while accumulated aldehyde can reach the
  low-mM range, so PBA sees 10–100× more aldehyde than active-site ketone. For
  the catalytic cycle to run at all, `k6` must beat this sink by roughly that
  factor. Worth watching closely in the fit.

  **Found, 2026-09-02, and it is the biggest thing in the 4OMe progress
  curves.** `product_fate/ANALYSIS.md` shows the catalysed 4OMe-BnOH
  rate falls **linearly in the accumulated product**, which is `A′ = v − kA` —
  production minus first-order consumption of the aldehyde — on 24 of 29 curves
  against 0 for the hyperbolic law that reversible product inhibition would
  give. The stationary level `A∞ = v(S)/k` carries substrate order
  **+0.610 ± 0.067** against the **+0.577** measured on the rates beforehand,
  and the implied selectivity is `k_A/k_S` ≈ **54**: the oxidant attacks the
  aldehyde about fifty times faster than the alcohol. The same statistic on the
  catalysed BnOH blocks is **+0.119 ± 0.096** — no product-driven deceleration
  at all, at 0.386 mM of benzaldehyde against 4OMe's 0.214. So S2 is switched on
  by the 4-methoxy group and effectively off without it, which is what an
  electron-rich ring should do to an electrophilic oxidant and to
  Baeyer–Villiger migration alike. Note that this cuts against steps 1–2 for
  4OMe in the same stroke: the substituent that makes the aldehyde a better
  target for the oxidant makes it a worse **hydride acceptor**, and the 4OMe
  curves show no autocatalysis while the BnOH curves do. `COMPUTATIONAL.md` C5
  and C6 are the two calculations that would confirm the pair of signs.
- **S4. `2 H2O2 -> 2 H2O + O2`, catalysed.** *(Found 2026-09-03; the last of
  these sinks to be identified and the only one whose product leaves the
  solution.)* Many BnOH pyrophosphate curves rise, fall by more in one 60 s
  reading than the reaction moves in five, and resume. Absorbance that goes
  away was never product: it is an O2 bubble in the SAMPLE beam, scattering
  light out of the aperture while it grows and releasing it when it detaches.
  `two_axis/ANALYSIS.md` section 5 carries the identification; three things make
  it this reaction rather than any other.

  **It needs the catalyst, and the double-beam layout is what says so.** The
  reference cuvette of a catalysed run omits only the enzyme -- same peroxide,
  substrate, buffer and pH, verified structurally across the archive by
  `verify_enzyme.py` -- so if the gas came from peroxide standing in solution
  the reference would bubble too and the difference would cancel. It does not:
  the archive's large steps run **122 falls beyond 20 sigma against 23 rises**
  (`bubble_step_asymmetry`), which puts the gas in the cuvette that holds the
  enzyme. Every run is therefore its own matched +/- enzyme control.

  Two weaker arguments were carrying this until 2026-09-04 and should not be
  quoted for it again. `bubble_turnover_control` -- the two weakest runs at
  73.4 mM shedding nothing -- is **confounded with pH**: those runs sit at pH
  6.95 and 7.53, and pH alone predicts their zero
  (`turnover_control_confound`). And the archive's four enzyme-free runs are
  too short to say anything: matched to their own buffers' catalysed rates they
  are worth about **one expected event** (`gas_enzyme_control`).

  **What that does NOT establish is that the KETONE is the catalyst.** Anything
  arriving with the enzyme stock would do the same job -- the cyclodextrin
  scaffold, or a trace transition metal, which is catalase-like and would show
  exactly this pH dependence. The archive holds no run with the cyclodextrin
  alone, the ketone alone, or a chelator added. The attribution to the ketone
  rests on the literature precedent (refs 34-35) and on the coincidence of
  conditions with the productive chemistry, not on a measurement here.
  `COMPUTATIONAL.md` C10.

  **It is not a property of one block, and pH is the trigger rather than
  peroxide.** The same detachment test over all **402 curves of 88
  experiments** (`scope.gas_curves`) finds the gas with **both substrates** --
  20 of 58 4OMe curves against 22 of 68 BnOH, catalysed, above 40 mM and pH 8,
  in three buffers (`gas_substrate_control`). That is the prediction S4 makes,
  since a catalyst decomposing peroxide involves no alcohol at all. Meanwhile
  the archive's median [H2O2] is 82.5 mM and 278 of 402 curves sit above 80,
  yet only 28 of those chop: what separates them is pH, and inside every buffer
  the rate climbs with it from a hard floor of **zero detachments in 270 hours
  below pH 7.5** (`gas_survey`). A reaction consuming HOO- is what that
  describes.

  **It is made from the peroxide and not from the alcohol.** The production rate
  is fitted from the timing and size of the detachments alone -- the fit never
  sees a concentration -- and comes out **+1.389 +/- 0.251 in [H2O2]** against
  **-0.307 +/- 0.103 in [S]** (`scope.gas_rate_drivers`). First order in
  peroxide, weakly negative in substrate: a catalase-like disproportionation
  competing with the productive cycle for the same oxidant, and one the alcohol
  slows rather than feeds.

  **The detachments track peroxide across the block**: 0 of 7 curves below 5 mM
  detach, 5 of 5 above 80 mM, monotone across six bands (`bubble_ladder`), and
  17 coincidences over 357 cuvette pairs against 16.0 expected, so it is the
  chemistry and not the lamp (`bubble_synchrony`).

  **THE GAS HAS NEVER BEEN MEASURED.** There is no headspace analysis, no
  manometry and no oxygen electrode anywhere in this project, so "O2" is an
  inference. What supports it: the budget below; first order in peroxide and
  negative in substrate. What excludes the obvious alternative is a pKa
  argument rather than a measurement -- CO2 generated at these pH values is not
  volatile, since H2CO3/HCO3- has pKa1 6.35 and above pH 8 essentially all
  dissolved inorganic carbon is bicarbonate. The gas also becomes MORE common
  as pH rises, which is backwards for carbonate, and that disposes of the one
  CO2 route the composition cannot otherwise exclude -- oxidative degradation
  of the cyclodextrin itself, which would carry the same three signatures but
  would still be captured as bicarbonate at pH 9. The defensible statement is
  a non-condensable gas, made in the enzyme-containing cuvette, first order in
  peroxide and rising steeply with pH.

  **Why it matters mechanistically, beyond the artefact.** It is an
  unproductive, catalyst-dependent drain on H2O2 that runs fastest at exactly
  the high-peroxide, high-pH conditions where the productive chemistry is
  strongest -- so it competes with step 4 for the oxidant and works against any
  enhancement being visible. It is not in the seven steps above and no rate
  constant here accounts for it. Whether it proceeds through the same perhydrate
  KP (a Criegee-type collapse expelling O2 rather than oxidising S) or through a
  separate path is not established by anything in this archive, and the
  distinction matters: the first would make the sink share step 4's
  pre-equilibrium and its saturation, the second would not. `solution_chemistry.
  oxygen_budget` disposes of the "too small to see" objection -- the solution
  saturates on 1.5% of the peroxide at the top of the ladder and cannot saturate
  at all at the bottom.

- **S3. Radical-chain autoxidation.** Benzoylperoxy radicals from aerobic
  benzaldehyde autoxidation abstract H from benzyl alcohol (Sankar et al., item
  21 — 2% benzyl alcohol suffices to suppress benzaldehyde autoxidation at room
  temperature), and photoexcited benzaldehyde can drive related radical
  chemistry (item 22, item 30). Either would produce autocatalytic-looking
  kinetics with no ionic mechanism involved at all. Not modelable from the
  existing data; needs dark/anaerobic/radical-scavenger controls to exclude.

## Buffer chemistry — the buffers are not innocent

This matters directly for fitting: the four buffer systems in the dataset
(phosphate, pyrophosphate, boric, carbonate) are **chemically different
reagents**, not just different ways of setting pH. Rates measured in different
buffers at the same nominal pH are not comparable without accounting for this.

- **Borate is the most compromised, in three separate ways.** (i) It forms
  **peroxoborate** with H2O2 (K = 2.0e-8, so significant above pH ≈ 7.7); the
  anionic peroxoborates are *much* faster oxidants than H2O2 itself and
  "deliver the hydroperoxide anion at a lower pH than when H2O2 is used"
  (item 39). (ii) It generates **dioxaborirane**, a highly reactive
  three-membered cyclic peroxide with a barrier of only 2.8 kcal/mol — a
  competing electrophilic oxidant with nothing to do with the catalyst
  (items 40–41). (iii) Most damaging for this mechanism: **boric acid catalyses
  peroxyacid hydrolysis**, ~12-fold for both peracetic acid and mCPBA, with a
  maximum at **pH 8.4–9** (item 42). If the mechanism runs through PBA, borate
  buffer is actively destroying the key intermediate, with its maximum right in
  the middle of the dataset's pH range. Treat all boric-buffer points as
  suspect and say why.
- **Carbonate forms a different oxidant.** Bicarbonate + H2O2 gives
  **peroxymonocarbonate (HCO4-)**, formed within minutes near neutral pH, a
  two-electron oxidant ~300× faster than H2O2 for sulfide oxidation (items
  43–44). In carbonate buffer the effective oxidant is partly HCO4-, not H2O2.
  Interestingly HCO4- is itself a peroxyacid-type species with a *carbonate*
  leaving group — arguably a closer analogue to the PBA case than Oxone is.
- **Phosphate catalyses the very first step.** Sander & Jencks (item 45) showed
  H2O2 addition to aldehydes is subject to **both general acid and general base
  catalysis** — so buffer species accelerate perhydrate/adduct formation
  independently of pH. This is exactly the general acid/base catalysis expected
  at steps 1, 4 and 6, and it means buffer concentration is a real kinetic
  variable, not a nuisance parameter.
- **WHICH step the buffer catalyses in the ENZYME-FREE background is NOT
  settled** (noted 2026-09-01; reasoning, not sourced). Sander & Jencks is
  cited above for the first-order buffer term, and it is the right physics --
  a buffer species in the transition state of the rate-determining step -- but
  it is catalysis of H2O2 addition to a CARBONYL, and the background has no
  carbonyl at t = 0. There are two ways out and the data prefers the second:
    * *Benzaldehyde, the product, supplies the carbonyl as it accumulates.*
      That route is autocatalytic, and the background does not accelerate:
      1 of 16 curves at the fixed-buffer anchor, against 50 of 110 two-axis
      catalysed increments, and exps 67, 69 and 70 actively DECELERATE
      (background_reaction/ANALYSIS.md section 7). Disfavoured.
    * *The buffer catalyses the oxidant's formation, before the substrate is
      involved at all.* The measured rate law fits that shape: total peroxide
      order ~ +2.3 at fixed pH (+1.49 in [H2O2] plus +0.82 in [HOO-]) with a
      substrate order of only +0.32. Two peroxide-derived species in the
      rate-determining step and the substrate barely in it -- making the
      oxidant is slow, oxidising BnOH is fast. A general base assisting the
      proton transfer that makes it is then first order in buffer, which is
      what is measured.
  The second is a hypothesis consistent with the orders, not a result. What is
  established is the order itself (+1.19 to +1.31, 5 sigma from zero and
  indistinguishable from exactly 1) and that it is catalysis rather than the
  buffer supplying the oxidant -- the boric probe, next.
- **Does phosphate form a peroxo species too?** Raised 2026-09-01, and it is
  the right question to ask of a first-order buffer dependence: if
  HPO4^2- + H2O2 <-> HPO5^2- + H2O ran at all, the buffer would be a REAGENT
  making an oxidant, not a catalyst, and the enzyme-free rate law would mean
  something different.

  Peroxomonophosphoric acid (H3PO5) is real and a competent oxidant, but the
  case is much weaker than for borate or carbonate. It is prepared from P4O10
  or from concentrated H2O2 with strong acid, and in dilute near-neutral
  solution it HYDROLYSES to phosphate + H2O2 -- that is the downhill
  direction. Borate is different because B(OH)3 is a Lewis acid that adds HOO-
  directly; carbonate because its carbon is electrophilic (K ~ 0.3 M^-1, formed
  in minutes). Phosphate at pH 7-8 is H2PO4-/HPO4^2-, an anion, so attack by
  HOO- at tetrahedral phosphorus is electrostatically disfavoured and slow.
  Expect a very small equilibrium constant. That is an argument, not a
  measurement, and it is not conclusive: a tiny concentration of a much faster
  oxidant can still carry a rate.

  **The phosphate kinetics cannot decide it.** Both this and general acid/base
  catalysis are first order in a buffer SPECIES, and the enzyme-free design
  cannot resolve a species from the total: within exps 3 and 6, the only sweep
  of [buf] at fixed pH, log[buf], log[H2PO4-] and log[HPO4^2-] are the same
  variable (correlation 1.000000). Across pH there are two phosphate levels and
  everything moves together -- substituting the basic form for the total drives
  its VIF to 30.2 against 2.8. See background_reaction/ANALYSIS.md section 6b.

  **The BORIC run can decide it, and the answer is no** (2026-09-01,
  `scope.peroxo_buffer_test`). Borate is where the peroxo route is not a
  hypothesis: item 39 above, K = 2.0e-8, significant above pH ~7.7, and the
  peroxoborates are much faster oxidants than H2O2 itself. Exp 65 is boric at
  pH 8.51. Exps 65 and 67 are matched cuvette for cuvette -- same substrate
  ladder, same 122.426 mM H2O2, same temperature, instrument and .rre source,
  87.5 against 85.0 mM buffer -- so predicting exp 65 from a law fitted on
  phosphate alone needs only the [buf] and [HOO-] orders. Excess over
  prediction: 0.66x (v0), 1.35x (vmax), 0.83x (v0_whole), 0.20x (v0_quad).
  NO EXCESS, where tens of mM of a much faster oxidant would show as a
  multiple. And the first-order buffer term is itself a phosphate number:
  dropping exp 65 from the anchor too leaves +1.31 +/- 0.23 (vmax).

  **WEAKENED 2026-09-01, same day, from the plots.** Exp 65's four cuvettes
  share a mid-run breakpoint at 504-560 s -- a span of 56 s, two sampling
  intervals -- across which every one of them STEEPENS, by 1.82, 2.04, 5.59
  and 15.94x, most at the LOWEST substrate. So a single rate number does not
  describe exp 65 -- vmax reads the post-break stretch and v0 the pre-break
  one -- and the excesses above compare two different processes.

  **CORRECTED the same day.** The first version of this note argued the break
  was not borate chemistry because exps 59-62 (boric, same and HIGHER pH) are
  smooth. That does not follow: those are with_E runs, whose reference channel
  omits the ENZYME, so the background is in both beams and CANCELS and the
  curve is a catalytic increment. They cannot show a background shape at all.
  Restricted to the 21 enzyme-free experiments -- the whole population in
  which a background feature is visible -- 19 are phosphate or pyrophosphate
  and none breaks (0.22-1.23), and exactly ONE is boric: exp 65, which breaks
  on all four cuvettes. One boric run of one, against 19 of 19. That is
  CONSISTENT WITH BORATE CHEMISTRY and rests on a single run. Exp 64, the only
  other boric background run, was aborted at 448 s -- before exp 65's break --
  and is dead besides.

  **The second probe was withdrawn as blind.** Exps 66 and 68 (same design
  with enzyme) were added to avoid exp 65 and cannot do the job: both are
  catalysed, so a buffer-made oxidant acts on the background, sits in both
  beams of each run, and cancels. Kept as the catalysed comparison it is --
  boric/phosphate = 0.66x (v0), 0.65x (vmax), 0.54x (v0_whole), 1.11x
  (v0_quad) against 2.19x the [HOO-], i.e. the CATALYSED reaction is 2-4x
  slower in borate, which is what item 42's peroxyacid hydrolysis predicts.

  NET: the peroxo question is OPEN, and the one piece of direct evidence now
  points mildly TOWARD a buffer-made oxidant rather than against it. The
  buffer term is still most likely catalysis on the chemistry above, but that
  is an argument, not a measurement. Two experiments would close it: the 31P
  NMR below, and A REPEAT OF EXP 65 -- enzyme-free, boric, run long.

  **Decisive test: 31P NMR of the buffer under run conditions.**
  Peroxomonophosphate is a distinct resonance; one spectrum settles presence or
  absence with no kinetics at all. Do this before writing the interpretation up
  either way. Supporting tests: saturation in [H2O2] (a pre-equilibrium adduct
  saturates, general catalysis does not -- our peroxide order is >= 1
  everywhere, which is weak evidence against); a Brønsted plot across buffers
  at one pH; and [buf] swept at a THIRD pH, which would break the
  species/total degeneracy outright.

- **The progress curves decay, and it is NOT catalyst inactivation** (tested
  2026-09-02, `scope.selwyn_test`). Fourteen of the temperature series' 24
  curves have a rate that RISES to a maximum and then FALLS -- exp 16 peaks at
  1666-2303 s and ends 30% below its peak. Four candidate causes, three of them
  now excluded:
    * *substrate depletion* -- no: conversion at the peak is 0.5-1.1%, and
      [H2O2] is 82.5 mM against ~0.01 mM of product.
    * *catalyst inactivation* -- no, by the classical Selwyn test. Exps 59 and
      60 are the archive's only pair differing in [enz] and nothing else that
      matters (0.028 against 0.014, exactly 2.000x). Inactivation requires
      P(low [E]) / P(high [E]) at matched [E]0 t to fall BELOW 1, because the
      slow run has twice the real time to decay. It is 1.10 to 2.20, median
      1.58, and never below 1 on any of 20 comparisons.
    * *a product threshold* -- not supported: the product at the rate maximum
      spans about fourfold across the runs where the peak is not truncated
      (0.024 to 0.098 AU), rather than sitting at one level.
    * *loss of oxidant to catalytic H2O2 decomposition* -- not tested, and the
      remaining candidate along with product inhibition proper.

  WHAT THIS MEANS FOR FITTING. A second phase belongs in the rate law
  empirically -- the one-exponential burst/lag form cannot hold a rate that
  rises then falls, and forcing it produces residuals of 3.8x noise at 40 C
  against 1.4x for a two-exponential -- but its time constant must NOT be read
  as a catalyst decay rate constant or given an activation energy, because the
  process behind it is unidentified. See temperature_series/ANALYSIS.md.

- **The rate may be markedly SUB-FIRST-ORDER in catalyst, and this is
  unresolved.** The same Selwyn pair implies an order of **+0.34**: at matched
  [E]0 t a rate going as [E]^n gives P proportional to [E]^(n-1), so a ratio of
  1.58 across a twofold pair gives n = 1 + ln(1.58)/ln(0.5).

  It matters because every turnover number and every Eyring ENTROPY in this
  project divides by [enz] to the first power. dH is untouched -- it is a slope
  -- but dS on the temperature series would move -53.5 -> -46.4 J/mol/K and
  dG(298) 103.6 -> 101.4 kJ/mol.

  AGAINST IT: the pair is BnOH in boric at pH 8.51, one pair, on days that are
  not controlled, with [buf] 77.0 against 75.0; and the 4OMe temperature series
  prefers the opposite -- normalising by [enz]^1 gives an Arrhenius rms of
  0.078 against 0.109 for [enz]^0.34, a 40% worse fit, though its own [enz]
  lever is only 11.7%. A dedicated [enz] ladder at one condition would settle
  it in an afternoon and is the single most valuable missing experiment for the
  activation parameters.

- **Pyrophosphate is probably acting as a metal chelator** (reasoning, not
  sourced): trace Fe/Cu catalyse H2O2 decomposition, which is why essentially
  every dioxirane paper adds EDTA. Pyrophosphate chelates; phosphate and
  carbonate do not, to the same degree. A fourth, non-mechanistic reason for
  buffer-to-buffer differences.

  There is weak indirect evidence, from two directions. Our enzyme-free runs —
  all phosphate and boric — run a median **34×** above the literature's
  uncatalysed rate once pH is matched (open questions, below), which is the
  size of excess an unsuppressed trace-metal path could produce. And matched on
  [HOO⁻] and corrected for [H2O2], catalysed **phosphate** (exp 68) runs
  **2.3× faster** than catalysed **pyrophosphate** (exp 138) despite carrying
  *less* enzyme, 0.028 against 0.034 mM. The chelating buffer is slower in both
  comparisons.

  Neither is proof — no metal was measured, buffer catalysis (Sander & Jencks,
  above) would produce the same signs, and an EDTA control would separate them
  in an afternoon. Note also that the excess is **not** what makes the paired
  controls come out flat: that is the catalyst loading, and it would be flat at
  a background 34× lower too.

Practical consequences for the fitting work: report rate constants against
**buffer concentration at fixed pH** within each system before interpreting the
pH profile at all; treat boric and carbonate points as suspect for the reasons
above; and note that a buffer-dilution series is the only clean way to separate
buffer catalysis from pH — which the existing titration data cannot do, since
`[buf]` and `[sub]` were varied together (see `DATA_VERIFICATION.md`).

**Measured 2026-09-01, and it is first order.** The dilution series does not
exist, but the buffer order is still recoverable, because the archive holds
*two* enzyme-free BnOH designs that disagree: exps 67/69/70 hold `[buf]`
constant along their substrate ladder and read a substrate order of
**+0.343 ± 0.071**, while exps 3/6 let `[buf]` fall 85 → 25 mM as `[sub]` rises
and read **−0.298 ± 0.087**. The order changes *sign* with the buffer design,
and the gap divided by the runs' own `dlog[buf]/dlog[sub]` gives an order in
buffer of **+1.31 ± 0.23** — confirmed independently at **+0.83 ± 0.29** on the
4OMe-BnOH / 40 °C block, whose buffer contrast sits between experiments at
fixed pH and fixed `[H2O2]`. First order in buffer is exactly what Sander &
Jencks predict. `python data/scope.py --buffer`; `scope.buffer_dependence`.

This does **not** reach exps 135–151: `[buf]` is 75.013 mM in all 119 of their
curves, so the two-axis block's substrate order carries none of this. It also does not
change what the catalysed curves *are* — the reference channel omits only the
enzyme, so the background is already subtracted at matched buffer.

**Noted 2026-09-02 — the same argument has never been made for phosphate, and
this section is the reason it should be.** Two of the four buffers above are
already granted a peroxo species that is the effective oxidant: borate forms
peroxoborate, and bicarbonate forms peroxymonocarbonate, "itself a peroxyacid-
type species with a *carbonate* leaving group — arguably a closer analogue to
the PBA case than Oxone is". Phosphate is treated in the same list as a general
acid/base only. But the leaving group is exactly what step 6 is short of:
closing a dioxirane out of a Criegee adduct made from plain H₂O₂ means expelling
hydroxide, which is why this chemistry is normally run on peroxymonosulfate or a
peracid. A **phosphate perhydrate** — `H₂O₂ + P ⇌ P–OOH`, then `P–OO⁻`
delivering the oxygen and leaving as phosphate — would play the peracid's part
in step 6 with the buffer in place of the product, and it would do it at 50–200
mM instead of at whatever the aldehyde has reached.

The kinetic consequence is a term in `[buf][H₂O₂]` where a general base gives a
term in `[buf]` alone, and the archive cannot see the difference: **0 of its 88
runs step both**, and all five buffer titrations sit at 82.5 mM peroxide. What
it does have is `induction/` §6's parameter-free pre-equilibrium constraint,
which the **buffer axis meets** (+1.094 ± 0.150 against a required +1) and the
**peroxide axis misses** (2.6σ and 3.7σ). `buffer/` §6 has the argument, and
`COMPUTATIONAL.md` **C9** the calculation — whose gate is that orthophosphate's
perhydrate equilibrium must come out unfavourable, since peroxomonophosphate is
made from P₂O₅ and not from phosphate in water. Pyrophosphate, which has a
leaving group where orthophosphate has none, is the sharp prediction: it should
be much the better catalyst at matched pH and matched peroxide, and `buffer/` §5
is why the archive cannot check.

## What the data says without a model

Everything below is measured, not fitted. It comes from `data/scope.py` and
`data/curve_metrics.py`; `DATA_VERIFICATION.md` (2026-08-31) carries the
working. None of it depends on the reduced system in this document being right,
so none of it moves if the fit changes.

### Reaction orders, measured within runs

The primary scope (exps 135–151, BnOH / 25 °C / pyrophosphate) is the only
block in the archive that varies both substrate and peroxide *inside* single
runs: 100.0% of its log[S] variance and 94.1% of its log[H2O2] variance is
within-experiment, so neither order can be absorbed by a per-experiment offset.
Orders are read off ladders — cuvettes of one run sharing a value of the
partner axis — over the eleven runs whose own cuvettes predict their own rates
(`scope.strong_runs()`: `concentration_agreement` ≥ 0.70, which separates the
block cleanly — the eleven score 0.72–0.97, the five excluded score 0.61 down
to 0.005, and exp 151 cannot be scored at all, its cuvettes scattering 234-fold
with two negative rates):

| | in [BnOH] | in [H2O2] |
|---|---|---|
| `v0`, the pre-catalytic rate | +0.43 | +0.69 |
| `vmax`, the developed rate | **+0.01 ± 0.04** | **+0.87 ± 0.06** |
| `gain` = `vmax`/`v0` | **−0.37 ± 0.07** | +0.18 |

Fitting `vmax` within experiments gives R² = 0.88 against 0.53 for the same
eleven runs pooled; across all seventeen it is 0.82 against 0.29. That gap is
the whole argument for scoping — pooled, most of the variance is between runs,
and a per-experiment offset would absorb the order being measured.
`python data/scope.py --orders` reproduces the table.

**The `vmax` substrate order of +0.01 is not saturation.** Saturation and
inhibition both flatten an average order, and only the order *as a function of
concentration* tells them apart — which is what `scope.local_orders` returns,
as the log-log slope between adjacent rungs:

| [BnOH] band | `v0` | `vmax` |
|---|---|---|
| below 1 mM | +1.04 | +0.30 |
| 1–3 mM | +0.60 | +0.35 |
| above 3 mM | +0.17 (6 of 15 negative) | **−0.386 (13 of 15 negative)** |

(The bands above pool all 17 runs, to keep the rung counts usable. On the 11
strong runs alone the top band is sharper still: **−0.457, negative in 11 of
11**.)

`v0` decays toward zero and stays positive — ordinary saturation, apparent Km
around 2–3 mM. `vmax` **turns over**. The developed catalyst is inhibited by
its own substrate at concentrations where the pre-catalytic rate is merely
saturating.

The obvious reading is substrate crowding the cyclodextrin cavity so the
peracid cannot bind — but see the paired controls below, which show the same
turnover with no cyclodextrin present. It should not be asserted as the
explanation.

### Autocatalysis

**51 of 110 live two-axis curves are steeper later than at the start by more
than 3σ.** It tracks the peroxide anion: 87% of the 30 curves above 0.1 mM
[HOO⁻] accelerate, against 31% of the 80 below. Exps 140–143 accelerate hardest (median z of +5 to
+14, 19 of their 28 live curves); exps 146, 149 and 150 decelerate outright
(median z near −20).

(It read 48 and 28% until 2026-09-01, when `line_fit`'s variance floor stopped
being hardcoded at the .txt export's quantisation on curves read from the
instrument's own `.rre`. The z-score divides by two standard errors that floor
sets, so the export's rounding had been suppressing it on data a thousand times
finer. `DATA_VERIFICATION.md` carries the working; three curves changed verdict
and nothing else in this section moved.) The block-slope statistic (`curve_metrics.acceleration`) is used
rather than a point-wise gradient because a point-wise gradient is too noisy to
locate anything here: its own scatter is a median 12% of the largest gradient
in the curve, and resampling a curve's noise moves the gradient's argmax by a
median 14% of the run. Where the rate peaks is partly a draw; whether the curve
is steeper later than at the start is not.

The induction period — the time to reach half the eventual rate,
`curve_metrics.lag_time` — has a **median of 23 minutes and an interquartile
range of 18–35**, measured over the 35 accelerating curves in the 8 strong runs
that carry at least three of them. Two things about it matter.

**It does not track substrate**: r = +0.03 against log[BnOH]. And **it does not
track conversion**: the fraction converted when the rate switches spans 332-fold
across those same curves, 0.11% to 37%. A product-derived autocatalyst —
benzaldehyde → perbenzoic acid, the obvious candidate — would have to fire at a
roughly fixed conversion, and does not.

**Nor is it a per-run clock.** If the induction were set by mixing or thermal
equilibration it would be a property of the run, but the spread *within* a run
(median 27 min) is as large as the spread of run medians (48 min, 16 to 64).

What it does track, weakly, is peroxide: r = −0.35 against log[HOO⁻] and +0.40
against log[H2O2] — one effect, not two, since the two axes move together, and
weak enough that it constrains little. What survives is that the active oxidant
must be *assembled* before it works, on a timescale set by something other than
how much substrate has been consumed.

### What draws the catalyst into its active form

*(Added 2026-09-04.)* Step 4 puts the catalyst into a **pre-equilibrium with a
species held in excess**, and that shape makes a prediction with no free
parameters. For `E + X ⇌ E*` with X in excess,

    1/τ = k_f[X] + k_r,     [E*]/E₀ = K[X]/(1 + K[X])

so `d ln v/d ln[X] = 1/(1 + K[X])` and `d ln τ/d ln[X] = −K[X]/(1 + K[X])`, and
their **difference is exactly 1 for every K and every [X]**. It is not a fact
about H₂O₂: it holds for *any* species that draws the catalyst into its active
form, so it can be asked of any axis the archive moves.
`induction.joint_order` asks it as one regression rather than two differenced by
hand, since the two orders come off the same curves on the same design.

**On the buffer axis it is met: +1.094 ± 0.150.** Read through a landmark with a
window in seconds common to both runs and one free level per run
(`induction.joint_buffer_order`, eight curves of exps 32 and 34). Eight curves
is a weak test — but it is a weak test of a parameter-free prediction, which is
not the same as no test, and its planted-recovery check reads back the scheme it
is built from. **This is the strongest statement the archive makes about what
E → E\* runs on**, and it points at the buffer rather than at the peroxide.

**On the peroxide axis it falls short, and how far depends on the clock.** The
induction landmark is a rolling window a tenth of the run wide, so on the
two-axis block — run length spanning 9.6× — it is not comparable between runs,
and `signal_control` fails there at +0.619 ± 0.228. `tau` and `tau_slow` come
from the progress fit instead, carry no window, and are subject to the same
identity. Both sides of the ratio are read off the O₂-corrected curves (S4
above), which matters here because the gas is made *from* peroxide and would
otherwise inflate the rate's order and shorten the clock, flattering the +1:

| clock | curves | order in [H₂O₂] | from +1 | control, [S] |
|---|---|---|---|---|
| `t_ind`, windowed | 110 | +0.304 ± 0.188 | 3.7σ | 6.4σ |
| `tau`, from the fit | 65 | +0.707 ± 0.158 | 1.9σ | 8.5σ |
| `tau_slow`, from the fit | 32 | +0.669 ± 0.242 | 1.4σ | 6.0σ |

**The substrate axis is the control and it must miss**, because the alcohol is
not the activating species and the clock carries no substrate order. It misses
by 4.2σ to 8.5σ in every cut. A reading where both axes met +1 would be a
regression that had stopped discriminating rather than a mechanism.

**Nothing is concluded from the peroxide column.** `tau_slow` is resolved on 32
of 110 live curves and 20 of the 77 strong ones, and across cuts the estimate
moves +0.67 to +0.85 — short of +1 everywhere and nowhere by enough to reject
it. What the two axes together say is that *something* in excess activates the
catalyst, that the buffer meets the constraint where the peroxide does not, and
that the archive cannot yet choose between them: `induction.peroxide_crossing`
finds that of 88 runs, 53 step `[buf]`, 20 step `[H2O2]` and **0 step both**.

**2026-09-05 — and the peroxide axis is not merely unresolved, it is
unmeasurable here.** Every block in the archive that moves `[H2O2]` has the
peroxide as its own signal: once the run offsets are projected out, log[H2O2]
correlates with log(net/noise) at **+0.57** on the two-axis block and **+0.67**
on exps 127–131. Holding the signal takes the induction clock's peroxide order
from +0.68 to +0.38 on one and from +0.19 to −0.13 on the other. That is true of
the windowed landmark (§4a of `induction/ANALYSIS.md`) *and* of the window-free
clock, so it is a property of the design and not of the statistic, and no
refinement of the statistic will fix it. **A peroxide dependence of the
induction cannot be measured in this archive at all.**

**What can be measured is pH, and it says the same thing the peroxide sign was
pointing at.** pH is one value per run everywhere here, so it needed a lag
statistic with no window in it (`scope.frame.lag_half_s`) and a refit of every
run on a window its ladder shares. There are four such ladders — 4OMe/phosphate
pH 5.64–8.95, 4OMe/boric pH 8.46–10.34, and the two-axis block's own two, BnOH
in pyrophosphate — and they are confounded in *opposite* directions, pH against
log(run length) running −0.25, +0.71, −0.53, −0.79 and against
log(signal/noise) +0.87, −0.65, +0.77, +0.79. All four coefficients on the clock
are positive and they agree (χ² = 0.95 on 3):

> **d ln τ / d pH = +0.12 to +0.34** across windows, pooled over four ladders,
> three buffers and both substrates. More alkaline, longer induction.

Positive is the sign a **trap** requires and the wrong sign for the species that
activates. Converted to a saturation fraction it is **5–15 %**: whatever holds
the catalyst back is only slightly engaged over the archive's whole pH range, so
this cannot yet name the species — but it removes the one objection that made
§4a's peroxide sign worthless, which was that both blocks carrying it were
signal-confounded. The pH ladders are not, in the same direction, and they
agree. `induction/ANALYSIS.md` §7e.

**And the substrate is on the other side of it.** On exps 135–151 — the only
block where `[S]` moves with `[buf]` held fixed — the clock **shortens** with
substrate at every floor, `d ln τ/d ln[S] = −0.18 to −0.71`. That block fails
its own signal control and can still carry this, because its rate order in `[S]`
is +0.09 ± 0.05: substrate buys no signal there (r = +0.04), which is the
structural reason the same block cannot be asked about peroxide and can be asked
about substrate. A negative order is the **activating** branch of the same
bounded algebra, 18–72 % engaged. Whatever the catalyst is waiting for, base
holds it back and the alcohol pulls it forward.

Two consequences for this document. The rate is **not** first order in H₂O₂ —
`peroxide_saturation` rejects `a = 1` at F = 41 on the two-axis ladder — so step
4's pre-equilibrium is saturating, and "first order in H₂O₂" is the
*unsaturated* limit of the scheme rather than a consequence of it. And every
buffer order here is an order in **total** buffer: at one pH the acid, the base
and the total are proportional, so no titration can name the species, and
`buffer/ANALYSIS.md` records that the two-pH test the archive does hold excludes
nothing (+1.06 ± 0.77 against 1.76 for general base and 0.52 for general acid).

### The +/− chemzyme controls

Exps 65–71 run the same substrate ladder, `[H2O2]`, buffer, pH and temperature
twice, once without the chemzyme and once with it at 0.028 mM
(`scope.PAIRED_CONTROLS`, `python data/scope.py --controls`). They sit outside
the primary scope — phosphate and boric buffer — so they inform interpretation
and cannot be pooled into a fit. Two results bear directly on this document:

**Autocatalysis needs the chemzyme.** At matched [HOO⁻] of 0.03–0.10 mM, 0 of
16 enzyme-free curves accelerate against 7 of 23 catalysed (Fisher p = 0.029).
No enzyme-free run reaches [HOO⁻] above 0.1 mM, where the catalysed block
reaches 87%, so the strongest regime is untested without enzyme. This is the
strongest evidence in the archive that the chemzyme is chemically active, and
it is the one result a fast background cannot manufacture: a background makes
curves faster, not sigmoidal.

**The rate enhancement, however, is invisible here** — 0.63× over 9 live rungs,
inside the 1.55× that separates exps 69 and 70. That is a statement about these
conditions and not about the catalyst: at 0.028 mM the literature's own `kcat`
predicts an enhancement of only 1.3× over the background these runs have, which
is smaller than their own reproducibility. See the open question below.

**The uncatalysed reaction is first order in the peroxide anion.** Fitted
across all six enzyme-free BnOH runs, which span pH 6.71 to 8.51
(`scope.background_orders()`, 24 live curves, R² 0.96): **+0.20 ± 0.06** in
[BnOH], +0.96 ± 0.31 in [H2O2], **+0.84 ± 0.04** in [HOO⁻]. So the background
climbs roughly tenfold per pH unit — measured directly, the rate rises 10× for
a 23× rise in [HOO⁻] between exps 3/6 and exps 67/69/70 at matched [H2O2].
Any comparison of backgrounds across pH that does not correct for this is
meaningless, and the near-zero substrate order is itself a puzzle: the
literature's `kuncat` is first order in substrate and ours is not, so the two
may not be the same reaction.

**The substrate turnover does not.** Above 3 mM the clean enzyme-free runs give
a local `vmax` order of −0.245, negative in both available pairs, against
−0.386 catalysed; below 3 mM, +0.182 against +0.338. Two pairs settle nothing
and this does not exclude cavity binding — but the turnover appears with no
cyclodextrin in the cuvette, so **substrate crowding the cavity is not the
required explanation** and should not be asserted as one. Cuvette position and
signal starvation at the top rung were both excluded rather than assumed
(`DATA_VERIFICATION.md`, 2026-08-31).

## Open questions

- **Does the catalyst loading move the induction clock?** A unimolecular
  activation says it must not, and that is the only prediction of step 0 a
  concentration can falsify. The archive's two `[enz]` pairs cannot answer:
  exps 59 and 60 hold everything but `[enz]` and neither run has a lag to time,
  and exps 140/141 give a **1.8× longer** clock at 2.4× the catalyst — past the
  replicate floor, in the wrong direction, on two runs 0.07 pH units apart in
  the block whose signal control fails. One substrate ladder at two catalyst
  loadings on the catalysed 4OMe system at 25 °C settles it, and it is the
  cheapest experiment this document still wants.
  `induction/ANALYSIS.md` §7g.
- **The reduction has not been redone with step 0 in it.** Every fit in
  `FITTING.md` starts from a fully active catalyst at `t = 0`, and the induction
  says that is false — on the catalysed 4OMe block the catalyst takes
  thousands of seconds to arrive. Adding `Kh ⇌ K` costs one rate constant and
  one state and would change the `E0 = 0` limit not at all, so the sequential
  strategy survives; what it changes is every catalysed trajectory's first
  hour. Until that is done, a fitted `k5'` or `k6` from this archive is a rate
  constant for a catalyst that was assumed present and was not.
- **Is the ketone the catalyst that destroys the peroxide, or something else
  in the stock?** S4 establishes that the gas forms in the cuvette that holds
  the enzyme -- the double-beam layout makes every run a matched control, and
  the large steps run 122 falls against 23 rises. It does not establish *what
  in that cuvette* is responsible. The cyclodextrin scaffold, or a trace
  transition metal carried in with the stock, would each give a catalase-like
  decomposition with the same first-order peroxide dependence and the same
  steep rise with pH. Three cheap experiments separate them, none of which the
  archive holds: **cyclodextrin without the ketone**, **the ketone without the
  cyclodextrin**, and **the full chemzyme with a chelator (EDTA)**, all at high
  peroxide and pH 9. Until one is run, "the chemzyme's ketone catalyses the
  disproportionation" is the best available reading of the data and not a
  measured fact -- and neither is the identity of the gas, which nothing in
  this project has measured. An oxygen electrode settles both halves at once.
- **Does the perhydrate shed the O2, or does something else?** S4 establishes
  that the catalyst destroys H2O2. It does not establish *where in the
  cycle*. If KP collapses to `K + O2 + H2O`,
  the sink shares step 4's pre-equilibrium and its saturation, competes with
  step 5 for the same intermediate, and the two should scale together across
  the peroxide axis. If instead a second peroxide attacks KP, or the dioxirane
  KD is reduced by H2O2, the sink has its own concentration dependence and can
  be tuned independently. **This archive cannot separate them**: the gas rate is
  measured from detachment timing alone, it is first order in peroxide either
  way, and there is no run in which the productive and unproductive routes are
  moved against each other. `COMPUTATIONAL.md` C10 is the calculation that
  would; the experiment would be a peroxide ladder with O2 evolution measured
  directly rather than inferred from scattering.
- **Which species draws the catalyst into its active form — buffer, peroxide,
  or neither?** The pre-equilibrium constraint above is met on the buffer axis
  (+1.094 +/- 0.150) and not met on the peroxide axis (1.4 to 1.9 sigma short
  through the unwindowed clocks, 3.7 through the landmark). Taken at face value
  that says the buffer, not H2O2, is what E -> E* runs on -- which would make
  step 4 the wrong entry point and put a buffer adduct or a general-base
  deprotonation ahead of it. Three things stop that being a conclusion. The
  buffer result is **eight curves**. Every buffer order in this project is an
  order in TOTAL buffer, so it cannot name the acid, the base, or a perhydrate
  of either. And `induction.peroxide_crossing` reports that of 88 runs, 53 step
  `[buf]`, 20 step `[H2O2]`, and **0 step both** -- so no run in the archive
  separates a buffer term from a `[buf][H2O2]` term, which is exactly the
  difference between a general base and a buffer perhydrate delivering the
  oxygen. **The missing experiment is a two-dimensional buffer x peroxide grid
  at fixed pH**, and it is cheap.
- **What does the absorbance actually measure?** This is now the single most
  consequential open question, ahead of any mechanistic detail. The structural
  analysis shows the mechanism predicts a decelerating aldehyde curve and an
  accelerating total-product curve, and the data shows the latter. Resolving it
  needs two cheap facts from the lab: the **monitoring wavelength**, and the
  **extinction coefficients of benzaldehyde, benzoic acid and the starting
  alcohol at that wavelength**. The wavelength is now settled (285 nm for BnOH,
  300 nm for 4OMe-BnOH -- each sheet declares it) and benzaldehyde's eps is known
  from the literature; **eps of benzoate at 285 nm is the one missing number**.
  Tracked as task C1 in `COMPUTATIONAL.md`, which also records why a naive
  TD-DFT calculation would answer it wrongly. Fitting `r` is possible but rests
  on varying the pH (established 2026-08-31): at a **single** pH and `[H2O2]`,
  `k_can` and `r` come out -0.999 correlated at a Jacobian condition number of
  6e7 — they enter the observable almost purely as the product `k_can*r`.
  Varying the pH breaks that, because pH enters only through `[HOO-]`, which
  multiplies `k_can` and not `r`. The real BnOH/25 C/phosphate block spans
  pH 6.71–8.01 and does separate them (`corr = +0.05`, condition number 523).
  So a designed pH series is what makes `r` estimable at all, and a measured
  `eps(benzoate)` would still be worth far more than a fitted one — the fitted
  value is only as good as the model producing it, and that model currently
  misfits. `data/fit_kinetics.py --profile-r` reports the cost profile over `r`
  alongside the fitted value. ORCA ground-state and TD-DFT calculations for
  benzaldehyde and benzoic acid (gas phase and solvated) are in `uvvis/`,
  started toward task C1.
- **No rate enhancement is visible in the only paired controls, and the
  arithmetic says the experiment could not have seen one.** Over the 9 substrate
  rungs where both sides carry a live signal, `vmax` with 0.028 mM chemzyme is
  **0.63× the enzyme-free value (range 0.31–1.41×)**. This is *not* a
  measurement of retardation — exps 69 and 70 are the same experiment run twice
  and disagree by up to 1.55× — and it is *not* evidence that the chemzyme is
  inactive.

  Take item 4's `kcat` (44×10⁻⁵ s⁻¹, Km 1.25 mM, at 0.4 mM catalyst) and ask
  what enhancement it predicts **at this archive's loading**, over the
  background these same runs have (`scope.predicted_enhancement()`,
  `python data/scope.py --literature`):

  | | ratio to background alone |
  |---|---|
  | predicted at 0.028 mM, our loading | **1.32× median (1.15–1.87×)** |
  | observed | 0.71× median |
  | predicted at the literature's 0.4 mM | 6× here, 13× at pH 7.5 |

  **The predicted effect is smaller than the run-to-run reproducibility.** No
  BnOH run in this archive is loaded above 0.069 mM, against the literature's
  0.4 mM — **6–29× below** — and the paired controls sit at pH 8.0–8.5 where
  the uncatalysed path is fastest. The enhancement is hidden by design, not
  absent.

  The background itself is not wildly anomalous either, once pH is handled.
  It runs a median **34× (12–88×)** above the literature's uncatalysed rate at
  matched pH and `[H2O2]` — a real but modest excess, plausibly buffer
  catalysis or trace metal. An earlier version of this passage said 876×; that
  compared a pH 8.01 background against a pH 7.0 literature value and was
  corrected on 2026-09-01, see `DATA_VERIFICATION.md`.

  So the missing experiment is **a catalysed series near 0.4 mM at neutral pH**,
  where the predicted enhancement is above tenfold and the background is
  lowest. An enzyme-free pyrophosphate run is still wanted for the fit's stage
  1 but would not by itself have made the enhancement visible. The primary
  scope cannot substitute: every run in it carries enzyme and spans only 4.9×
  in `[enz]`, all between-run, and regressing `log vmax` on `log[enz]` gives
  anything from +0.04 ± 0.54 to +1.21 ± 0.37 depending on how peroxide is
  controlled for.
- **The largest catalysed block has no background.** Rate constants may be
  pooled only within one (substrate, temperature, buffer) cell — Arrhenius,
  different molecules, and the buffer section above all forbid pooling across
  them. `data/fit_kinetics.py --list` shows that of eleven such cells, only two
  hold enzyme-free controls alongside catalysed runs: BnOH/25 C/phosphate
  (23 vs 20 curves) and 4OMe-BnOH/40 C/phosphate (37 vs 24 — it read 59 vs 4
  until exps 32 and 34-37 were ruled catalysed on 2026-08-31). The archive's
  biggest catalysed block — **127 BnOH pyrophosphate curves** — has no
  enzyme-free control at all, so `k_can, k3, k0` could only be imported from a
  different buffer, which is exactly what the buffer section says must not be
  done. An enzyme-free pyrophosphate control is a cheap experiment and would
  unlock the largest single block in the dataset.
- **No E0 titration exists.** Only one experiment varies enzyme concentration
  across its own samples, so the mechanism's cleanest prediction — that the
  reduced model is exactly linear in E0, and that the induction period shortens
  as E0 rises — cannot be tested against the current dataset at all. A single
  E0 series at otherwise fixed conditions (one pH, one buffer, one T, fixed
  [S] and [H2O2], 5–6 enzyme concentrations spanning a decade) would be the most
  informative new experiment available, because E0 is the only variable that
  cleanly separates the catalysed loop from the catalyst-free background.
- **Bols-group full text: item 4 retrieved, items 1–3 still paywalled.** Item 4
  (ChemCatChem 2025) turned out to be hybrid open access. Its kinetics for
  benzyl alcohol are now confirmed: diketone **8**, kcat = 44 ± 3.2 x10^-5 s^-1,
  Km = 1.25 ± 0.47 mM, kcat/kuncat = 28,000; monoketone **10**, kcat =
  12.7 ± 1.4 x10^-5 s^-1, Km = 2.42 ± 0.99 mM, kcat/kuncat = 8,000. Conditions:
  pH 7 phosphate, 25 C, 72 mM H2O2, 0.4 mM catalyst. It reports **no pH-rate
  profile and no pKa determination**, and proposes neither a peracid nor a
  dioxirane. Items 1–3 (2006 Angew., 2009 OBC, 2009 ChemBioChem) remain
  unretrieved — worth pulling via library access to see whether the earlier
  papers say more about mechanism than the 2025 one does.
- **Step 6's dioxirane-vs-Baeyer-Villiger branching is pH-controlled in
  principle, but the dataset may not be able to prove it.** The specific claim
  that the pH-sensitive event is deprotonation of a *neutral Criegee adduct
  before collapse* is not established anywhere: no pKa has ever been measured
  for the relevant adduct, and Lange & Brauer's rate law is algebraically
  equally consistent with the peroxide **anion/dianion adding to the ketone**,
  which would put the pH control one step earlier, in the addition. Worse for
  us: if HOO- (pKa 11.6) is the nucleophile in the addition step, rate rises
  across pH 5.5→11.8 with an inflection near 11.6 **from the addition alone** —
  so the pH data cannot by itself distinguish the two. The discriminator is
  *where the inflection sits*: near 11.6 points to peroxyanion nucleophilicity;
  well below 10 would support adduct ionization. An electron-poor sugar ketone
  might plausibly put the adduct pKa in the 8–10 range (reasoning, not sourced),
  which would be inside the dataset's window — worth looking for.
- **Steps 1–2 have no literature precedent at all** (see above). The strongest
  available move is to model them alongside the literature-favoured alternative
  (S1: unimolecular collapse straight to benzoic acid, no peracid at all) and
  let the progress-curve fits adjudicate. If the no-peracid variant fits the
  enzyme-free controls just as well, the Cannizzaro step is unfalsifiable from
  this dataset and should be presented as speculative.
- **A catalyst-independent radical/O2-chain pathway has not been ruled out** as
  an alternative (or additional) explanation for background autocatalysis in
  enzyme-free controls — benzaldehyde/O2 radical-chain autoxidation to
  perbenzoic acid is separately well documented (items 21–22) and would be
  killed by dark or anaerobic controls but not by the ionic mechanism proposed
  here. Worth checking whether any experiments were run under inert atmosphere,
  or whether dissolved O2 is a live confound.
- **Two papers worth retrieving via library access**, in priority order:
  item 33 (Porter/Yin/Pratt, JACS 2000) **and its Supporting Information**,
  which reportedly holds computed transition-state structures for a
  peroxyacid/dioxirane closure — the closest thing in existence to the missing
  DFT comparison; and item 32 (Schulz/Adam, JOC 1997), for the mechanism they
  propose for the sulfonic-peracid → dioxirane step and any leaving-group
  discussion.

## References

1. Marinescu, L. G.; Bols, M. "Very High Rate Enhancement of Benzyl Alcohol
   Oxidation by an Artificial Enzyme." *Angew. Chem. Int. Ed.* **2006**, *45*
   (28), 4590–4593. DOI: 10.1002/anie.200600812.
2. Fenger, T. H.; Marinescu, L. G.; Bols, M. "Cyclodextrin ketones as oxidation
   catalysts: Investigation of bridged derivatives." *Org. Biomol. Chem.*
   **2009**, *7* (5), 933–943. DOI: 10.1039/b814245a. PMID: 19225677.
3. Hauch Fenger, T.; Bjerre, J.; Bols, M. "Cyclodextrin Aldehydes are Oxidase
   Mimics." *ChemBioChem* **2009**, *10* (15), 2494–2503.
   DOI: 10.1002/cbic.200900448. PMID: 19739193.
4. Zimmermann, M. L.; Friis, V.; Zorck, W. F.; Bols, M. "Improvement of the
   Catalytic Effect of Cyclodextrin-Based Oxidase Mimics by Introduction of a
   Molecular Cap." *ChemCatChem* **2025**, *17* (8), e202500263.
   DOI: 10.1002/cctc.202500263.
5. Rousseau, C.; Christensen, B.; Petersen, T. E.; Bols, M. "Cyclodextrins
   containing an acetone bridge. Synthesis and study as epoxidation catalysts."
   *Org. Biomol. Chem.* **2004**, *2*, 3476–3482. DOI: 10.1039/b410098k.
6. Bjerre, J.; Fenger, T. H.; Marinescu, L. G.; Bols, M. "Synthesis of some
   trifluoromethylated cyclodextrin derivatives and analysis of their
   properties as artificial glycosidases and oxidases." *Eur. J. Org. Chem.*
   **2007**, *2007* (4), 704–710. DOI: 10.1002/ejoc.200600762.
7. Bjerre, J.; Rousseau, C.; Marinescu, L. G.; Bols, M. "Artificial enzymes,
   'Chemzymes': current state and perspectives." *Appl. Microbiol. Biotechnol.*
   **2008**, *81* (1), 1–11. DOI: 10.1007/s00253-008-1653-5. PMID: 18787819.
8. Marinescu, L. G.; Doyagüez, E. G.; Petrillo, M.; Fernández-Mayoralas, A.;
   Bols, M. "Amino–Acetone-Bridged Cyclodextrins — Artificial Alcohol
   Oxidases." *Eur. J. Org. Chem.* **2010**, *2010* (1), 157–167 (published
   online 16 Dec 2009). DOI: 10.1002/ejoc.200901099.
9. Fenger, T. H.; Marinescu, L. G.; Bols, M. "Cyclodextrin Ketones with the
   Catalytic Group at the Secondary Rim and Their Effectiveness in Enzyme-Like
   Epoxidation of Stilbenes." *Eur. J. Org. Chem.* **2011**, *2011* (12),
   2339–2345. DOI: 10.1002/ejoc.201001696.
10. Deary, M. E.; Davies, M. I. "Evidence for cyclodextrin dioxiranes."
    *Carbohydr. Res.* **1998**, *309* (1), 17–29.
    DOI: 10.1016/S0008-6215(98)00111-6. (Non-Bols precedent; Oxone-based, not
    peracid-based.)
11. Triandafillidi, I.; Kokotou, M. G.; Lotter, D.; Sparr, C.; Kokotos, C. G.
    "Aldehyde-catalyzed epoxidation of unactivated alkenes with aqueous
    hydrogen peroxide." *Chem. Sci.* **2021**, *12* (30), 10191–10196.
    DOI: 10.1039/D1SC02360H. Open access, PMC8336450. Best available precedent
    for step 6's peracid+carbonyl-catalyst→dioxirane bond reorganization, via
    an aldehyde catalyst rather than a ketone.
12. Poursaitidis, E. T.; Gkizis, P. L.; Triandafillidi, I.; Kokotos, C. G.
    "Organocatalytic activation of hydrogen peroxide: towards green and
    sustainable oxidations." *Chem. Sci.* **2024**, *15* (4), 1177–1203.
    DOI: 10.1039/D3SC05618J. Open access, PMC10806817. Review; confirms
    ketone-catalyzed epoxidations in this literature typically go via direct
    perhydrate formation (item 4/step), not a peracid intermediate.
13. "2,2,2-Trifluoroacetophenone: An Organocatalyst for an Environmentally
    Friendly Epoxidation of Alkenes." *J. Org. Chem.* **2014**.
    DOI: 10.1021/jo5003938. *Author list not independently confirmed — verify
    before citing.* Oxone/perhydrate-based mechanism, not peracid-based.
14. Rozhko, E.; Solmi, S.; Cavani, F.; Albini, A.; Righi, P.; Ravelli, D.
    "Revising the Role of a Dioxirane as an Intermediate in the Uncatalyzed
    Hydroperoxidation of Cyclohexanone in Water." *J. Org. Chem.* **2015**,
    *80*, 6425–6431. DOI: 10.1021/acs.joc.5b00861. (Full citation confirmed on
    the third research pass, superseding the partial one recorded earlier.)
    DFT: dioxirane formation proceeds via hydroperoxide-anion addition, and the
    dioxirane sits in equilibrium with the ketone.
15. Thangavel, A.; Elder, D.; Sotiriou-Leventis, C.; Dawes, J.; Leventis, N.
    "Control of the Ketone to gem-Diol Equilibrium by Host–Guest
    Interactions." *Org. Lett.* **2008**, *10*. *Exact issue/page not
    independently confirmed — verify before citing.*
16. Zauche, T. H.; Espenson, J. H. "Oxidation of Alcohols by Hydrogen
    Peroxide, Catalyzed by Methyltrioxorhenium (MTO): A Hydride Abstraction."
    *Inorg. Chem.* **1998**, *37* (26), 6827–6831.
17. Renz, M.; Meunier, B. "100 Years of Baeyer–Villiger Oxidations." *Eur. J.
    Org. Chem.* **1999**. *Exact issue/pages not independently confirmed —
    verify before citing.* General Baeyer-Villiger mechanism review.
18. Bach, R. D. et al. "Molecular Dynamics of Dimethyldioxirane C–H
    Oxidation." *J. Am. Chem. Soc.* **2016**, *138*. *Exact issue/pages not
    independently confirmed — verify before citing.*
19. "Oxidation of alcohols by dimethyldioxirane." *Russ. Chem. Bull.*
    *Full citation not retrieved — needs a proper literature search before
    citing.*
20. Toteva, M. M.; Richard, J. P. "Oxidation of benzyl alcohols by
    dimethyldioxirane." *Tetrahedron* **2001/2002**. *Exact year/issue/pages
    not independently confirmed — verify before citing.*
21. "The benzaldehyde oxidation paradox…" *Nat. Commun.* **2014**, *5*, 3350.
    *Author list not retrieved — needs a proper literature search before
    citing.*
22. Wang et al. "HOO• as the Chain Carrier for the Autocatalytic
    Photooxidation of Benzylic Alcohols." *Molecules* **2024**, *29* (14),
    3429. *Full author list not independently confirmed — verify before
    citing.*
23. "The Over-Riding Role of Autocatalysis in Allylic Oxidation." *Catal.
    Lett.* **2021**. *Full citation not retrieved — needs a proper literature
    search before citing.*
24. "Shedding Light on Synthetic Autocatalysis: From Conventional Closed-Shell
    Chemistries to Overlooked Open-Shell Occurrences." *Chem. Eur. J.*
    **2025**. *Full citation not retrieved — needs a proper literature search
    before citing.*

25. Braun, G. "Perbenzoic acid." *Org. Synth.* **1941**, Coll. Vol. 1, 431;
    *Org. Synth.* **1925**, *5*, 90. Classical benzoyl peroxide + alkoxide
    route to perbenzoic acid.
26. US Patent 3,321,512, "Manufacture of perbenzoic acids." Industrial
    benzoyl chloride → benzoyl peroxide → perbenzoic acid route (pH ≥ 10
    during reaction, acidification to pH ≤ 2 for isolation).
27. Swain, C. G.; Powell, A. L.; Sheppard, W. A.; Morgan, C. R. "Mechanism of
    the Cannizzaro Reaction." *J. Am. Chem. Soc.* **1979**, *101* (13),
    3576–3583. DOI: 10.1021/ja00507a023. Isotope-dilution and solvent-isotope
    evidence that the RDS is intermolecular hydride transfer from the
    gem-diolate to a second aldehyde; rules out an ester intermediate.
28. Alexander, E. R. "Studies on the Mechanism of the Cannizzaro Reaction."
    *J. Am. Chem. Soc.* **1947**, *69*, 289–294. *Cited via a later paper —
    verify directly before citing.* Source of the rate = k[ArCHO]²[OH⁻] rate
    law.
29. Ogata, Y.; Sawaki, Y. "Kinetics and mechanism of the Baeyer-Villiger
    reaction of benzaldehydes with perbenzoic acids." *J. Am. Chem. Soc.*
    **1972**, *94* (12), 4189–4196. The most directly relevant kinetics paper
    for side reactions S1/S2: pH–product profiles and migration isotope
    effects for benzaldehyde + perbenzoic acid.
30. "Benzaldehyde-Promoted (Auto)Photocatalysis under Visible Light."
    ChemRxiv preprint, **2022**. *Full citation not retrieved — verify before
    citing.* Photoexcited benzaldehyde oxidising benzyl alcohol while
    generating H2O2; a third route to autocatalytic-looking kinetics.
31. Kharasch, M. S.; Foy, M. "The Peroxide Effect in the Cannizzaro Reaction."
    *J. Am. Chem. Soc.* **1935**, *57* (8), 1510. Located during the search
    for peroxide-Cannizzaro precedent; on inspection almost certainly concerns
    trace peroxide acting as a radical initiator/inhibitor in the ordinary
    NaOH-Cannizzaro reaction, **not** HOO⁻ replacing OH⁻ as nucleophile.
    Full text not retrieved — do not cite as precedent for steps 1–2 without
    reading it.

32. Schulz, M.; Liebsch, S.; Kluge, R.; Adam, W. "Organo Sulfonic Peracids. 4.
    The Reaction of Arenesulfonylimidazoles with H2O2 in the Presence of
    Ketones. A New Entry to Dioxiranes." *J. Org. Chem.* **1997**, *62*,
    188–193. DOI: 10.1021/jo9610164. **The key citation for step 6** — an
    organic peroxy acid + ketone + NaOH gives a dioxirane, proven by 18-O
    labelling, dioxirane pathway "virtually the exclusive one." Leaving group
    is a sulfonate, not a carboxylate.
33. Porter, N. A.; Yin, H.; Pratt, D. A. "The Peroxy Acid Dioxirane
    Equilibrium: Base-Promoted Exchange of Peroxy Acid Oxygens." *J. Am. Chem.
    Soc.* **2000**, *122*, 11272–11273. DOI: 10.1021/ja005648q. SI reportedly
    contains computed transition-state structures. **Retrieval priority 1.**
34. Lange, A.; Brauer, H.-D. "On the formation of dioxiranes and of singlet
    oxygen by the ketone-catalysed decomposition of Caro's acid." *J. Chem.
    Soc., Perkin Trans. 2* **1996**, 805. DOI: 10.1039/P29960000805.
    Rate law for dioxirane *formation*, eleven ketones.
35. Lange, A.; Hild, M.; Brauer, H.-D. "New aspects concerning the mechanism
    of the ketone-catalysed decomposition of Caro's acid." *J. Chem. Soc.,
    Perkin Trans. 2* **1999**, 1343–1350. DOI: 10.1039/A901976F. Shows the
    gem-diol equilibrium and its first ionization must be modelled for
    electron-poor ketones.
36. Doering, W. von E.; Dorfman, E. "Mechanism of the Peracid Ketone–Ester
    Conversion. Analysis of Organic Compounds for Oxygen-18." *J. Am. Chem.
    Soc.* **1953**, *75*, 5595–5598. DOI: 10.1021/ja01118a035. **The
    counter-experiment**: benzophenone + perbenzoic acid, 18-O labelling
    excludes both dioxirane and carbonyl-oxide routes.
37. Armstrong, A.; Barsanti, P. A.; Clarke, P. A.; Wood, A. "Ketone-directed
    peracid epoxidation." *Tetrahedron Lett.* **1994**, *35*, 6155–6158.
    DOI: 10.1016/0040-4039(94)88103-0; and *J. Chem. Soc., Perkin Trans. 1*
    **1996**, 1373. DOI: 10.1039/P19960001373. 18-O labelling excludes a
    dioxirane in a second peracid+ketone system.
38. Poursaitidis, E. T.; Mantzourani, C.; Triandafillidi, I.; Kokotou, M. G.;
    Kokotos, C. G. "Green epoxidation of unactivated alkenes via the catalytic
    activation of hydrogen peroxide by 4-hydroxybenzaldehyde." *Green Chem.*
    **2025**, *27*, 11192–11202. DOI: 10.1039/D5GC02537K. Explicitly rules out
    a peracid intermediate in a closely related system.
39. Davies, D. M.; Deary, M. E.; Quill, K.; Smith, R. A. "Borate-Catalyzed
    Reactions of Hydrogen Peroxide: Kinetics and Mechanism of the Oxidation of
    Organic Sulfides by Peroxoborates." *Chem. Eur. J.* **2005**, *11*,
    3552–3558. DOI: 10.1002/chem.200401209.
40. Durrant, M. C.; Davies, D. M.; Deary, M. E. "Dioxaborirane: a highly
    reactive peroxide that is the likely intermediate in borate catalysed
    electrophilic reactions of hydrogen peroxide in alkaline aqueous
    solution." *Org. Biomol. Chem.* **2011**, *9*, 7249–7254.
    DOI: 10.1039/C1OB06142A.
41. Deary, M. E.; Durrant, M. C.; Davies, D. M. *Org. Biomol. Chem.* **2013**,
    *11*, 309–317. DOI: 10.1039/C2OB26842F. Extends item 40 across a range of
    substrates.
42. Deary, M. E. "Boric acid catalysed hydrolysis of peroxyacids." *RSC Adv.*
    **2023**, *13*, 11826–11837. DOI: 10.1039/D3RA01046E. ~12-fold rate
    enhancement of peracid hydrolysis, maximal at pH 9 (peracetic) / 8.4
    (mCPBA).
43. Richardson, D. E.; Yao, H.; Frank, K. M.; Bennett, D. A. "Equilibria,
    Kinetics, and Mechanism in the Bicarbonate Activation of Hydrogen
    Peroxide: Oxidation of Sulfides by Peroxymonocarbonate." *J. Am. Chem.
    Soc.* **2000**, *122*, 1729–1739. DOI: 10.1021/ja9927467.
44. Bakhmutova-Albert, E. V.; Yao, H.; Denevan, D. E.; Richardson, D. E.
    "Kinetics and Mechanism of Peroxymonocarbonate Formation." *Inorg. Chem.*
    **2010**, *49*, 11287–11296. DOI: 10.1021/ic1007389.
45. Sander, E. G.; Jencks, W. P. "General acid and base catalysis of the
    reversible addition of hydrogen peroxide to aldehydes." *J. Am. Chem.
    Soc.* **1968**, *90*, 4377–4386. DOI: 10.1021/ja01018a032. **The citation
    for buffer catalysis of the adduct-forming steps.**
46. Cremer, D.; Schmidt, T.; Gauss, J.; Radhakrishnan, T. P. "Formation of
    Dioxirane from Carbonyl Oxide." *Angew. Chem. Int. Ed. Engl.* **1988**,
    *27*, 427–428. DOI: 10.1002/anie.198804271. Carbonyl oxide → dioxirane
    cyclization is easy (~9–20 kcal/mol) and exothermic, so a stepwise variant
    of step 6 is not thermodynamically forbidden — though expelling benzoate to
    give a free zwitterionic carbonyl oxide in water looks less plausible than
    direct 3-exo-tet displacement.
47. Wang, Z.-X.; Shi, Y. "A pH Study on the Chiral Ketone Catalyzed Asymmetric
    Epoxidation of Hydroxyalkenes." *J. Org. Chem.* **1998**, *63*, 3099–3104.
    DOI: 10.1021/jo972106r. Shi's own published rationale for high pH is
    suppression of *uncatalysed background* epoxidation — not the deprotonation
    argument that secondary sources give.
48. Ball, D. L.; Edwards, J. O. "The Kinetics and Mechanism of the
    Decomposition of Caro's Acid. I." *J. Am. Chem. Soc.* **1956**, *78*,
    1125–1129. DOI: 10.1021/ja01587a011. pKa2(HSO5-) ≈ 9.3–9.4.
49. Denmark, S. E.; Forbes, D. C.; Hays, D. S.; DePue, J. S.; Wilde, R. G.
    "Catalytic Epoxidation of Alkenes with Oxone." *J. Org. Chem.* **1995**,
    *60*, 1391–1407. DOI: 10.1021/jo00110a049. Optimum pH 7.5–8.0.
50. Denmark, S. E.; Wu, Z. "Dioxiranes Are the Active Agents in
    Ketone-Catalyzed Epoxidations with Oxone." *J. Org. Chem.* **1997**, *62*,
    8964–8965. DOI: 10.1021/jo971500m.
51. Adam, W.; Curci, R.; Edwards, J. O. "Dioxiranes: a new class of powerful
    oxidants." *Acc. Chem. Res.* **1989**, *22*, 205–211.
    DOI: 10.1021/ar00162a002.

Items 14, 15, 17–24, 28, 30, 31 have incomplete author/page metadata from the research
passes that produced this document (search-snippet-only, full text not
retrieved) — verify each against the primary source before using it in the
thesis itself. Items 1–13 have DOIs and are believed complete and correct as
listed, but were themselves largely sourced from abstracts/metadata rather
than full text (see Open questions) — the mechanistic *claims* attributed to
items 1–4 in this document's step-by-step section are reasoning by the thesis
author and research assistant, not confirmed to be stated in those papers.
