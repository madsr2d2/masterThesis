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
| K | active-site ketone (catalyst) |
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

4. `K + H2O2 ⇌ KP` — the ketone forms its gem-diol hydroperoxide ("perhydrate")
   with H2O2.
5. `KP + S → K + A` — the perhydrate directly oxidizes the alcohol substrate to
   the aldehyde. Hypothesized as the slow, non-autocatalytic "seed" step that
   produces the first trace of A needed to start the autocatalytic loop, before
   enough product has accumulated for steps 1–2 to run.
6. `PBA + K → KD + BA` — the ketone catalyst reacts with the perbenzoic acid
   produced in step 2 (not with H2O2/Oxone directly), forming a Criegee-type
   adduct that collapses by expelling benzoic acid as the leaving group and
   leaving a dioxirane at the former ketone carbon.
7. `KD + S → K + A` — the dioxirane oxidizes the substrate's benzylic C–H bond,
   regenerating free ketone catalyst and producing a fresh molecule of A. Closes
   the catalytic cycle.

Total catalyst is conserved: `[K] + [KP] + [KD] = [enz]0` throughout (no rate
constant needed to enforce this — it falls out of steps 4–7 by construction).

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
for BnOH** and **7.53 mM^-1 cm^-1 for 4OMe-BnOH**. Benzaldehyde's strong pi->pi*
band at 250 nm has eps ~ 12,000–14,000 M^-1 cm^-1 — about ten times higher than
1230 — so the assay must be sitting on the weak n->pi* band near 280 nm. There,
benzoic acid (eps ~ 800 M^-1 cm^-1 at 273 nm) absorbs *comparably*, not
negligibly. The abs -> [P] conversion in the notebook assumes the aldehyde is
the sole absorber; it very likely is not. *(Both literature eps values are from
memory and should be verified against a spectrum before citing; the monitoring
wavelength should be confirmed from the lab notebooks.)*

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

**51 of 110 live in-scope curves are steeper later than at the start by more
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
  alongside the fitted value.
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
