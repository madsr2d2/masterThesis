# The early trough

A catalysed curve's absorbance is already reference-subtracted: the recorded
value is sample cuvette minus an enzyme-free reference cuvette holding the
same composition, so whatever reaction the two cuvettes share cancels
between them — that shared reaction is exactly what `background_reaction/`
measures directly. It cancels only while the two cuvettes see the same free
concentration of whatever it runs on. If the catalyst engages one reactant
fast enough to measurably deplete it in the **sample** cuvette alone, the
sample's own share of that shared background chemistry runs slower than the
reference's for a while, and the reported curve dips **below zero** before
the catalysed rate overtakes it and the curve turns positive.

Seventeen curves, across both substrates, do exactly this. `data/early_trough.py`
has the machinery; `progress_curves.html` shows every one of the seventeen
with its full fit.

## 1 · A trough that survives three separate ways to be spurious

`curve_metrics.early_trough` smooths the first half of a curve's run with a
9-reading rolling mean and reports its minimum, in units of the curve's own
noise. Over every live catalysed curve in the archive (`scope.archive()`,
**311 curves**), **17** clear a smoothed trough of −4.0 sigma or deeper.

A smoothed mean alone cannot tell a genuine, sustained decline from one deep
outlier reading dragging it down. Three independent screens separate them:

- **Sustained.** A LEAVE-ONE-OUT test: drop the trough window's single worst
  reading and require the remaining 8 to still average below −3.0 sigma. A
  real decline is robust to this; a single bad reading collapses without it.
  This design replaced an earlier one that required a fixed count of the
  window's 9 readings to individually clear a per-reading depth bar — and
  that version wrongly rejected two real curves, exps **4.1** (smoothed z
  **−7.31**) and **22.2** (z **−6.23**), whose declines are genuine (9 of 9
  window readings negative in both) but shallower per reading than the
  block's more dramatic examples. Leave-one-out passes both: dropping the
  single worst reading moves exp 4.1's mean only −7.31 → −6.12 sigma, and
  exp 22.2's only −6.23 → −5.10 — nothing like what a real single-outlier
  case does (a synthetic planted spike collapses from −10.8 to −0.7 under
  the same test). See `DATA_VERIFICATION.md` for the correction.
- **No bubble overlap.** The trough window must not overlap a detected O2
  detachment (`curve_metrics.detachments`). None of the archive's candidates
  do.
- **Survives debubble correction.** The trough recomputed on the
  `curve_metrics.debubble`-corrected curve must still clear −4.0 sigma.

**17 of 17** candidates pass all three — every candidate the scan found is
now confirmed genuine.

## 2 · The driver is the oxidant, not the substrate, in both substrates independently

Define `dominance = max([enz]/[S], [enz]/[HOO⁻])` — whichever reactant the
catalyst's own fixed concentration comes closest to overwhelming, taking the
**larger** of the two ratios because [enz] is negligible against most of what
this archive puts in a cuvette and the reactant that matters is whichever one
is scarce enough for [enz] to be a large fraction of it.

Spearman rank correlation of the trough depth against log-dominance, over
every scanned curve, by substrate:

| substrate | n | ρ against log[enz]/[HOO⁻] | p | ρ against log[enz]/[S] | p |
|---|---|---|---|---|---|
| 4OMe-BnOH | 147 | **−0.621** | **4.7 × 10⁻¹⁷** | +0.045 | 0.59 |
| BnOH | 164 | **−0.315** | **4.0 × 10⁻⁵** | −0.090 | 0.25 |

[enz]/[HOO⁻] is significant at p < 10⁻⁴ in **both** substrates scanned
independently; [enz]/[S] is not significant in either. The hydroperoxide
anion is nanomolar across most of the archive's pH range, so a catalyst held
at 0.014–0.273 mM is routinely a 100–4000-fold molar excess over it — the
same regime `MECHANISM.md`'s S4 and `BUBBLES.md`'s gas already put the
catalyst in, engaging the peroxide non-productively. This is that same
engagement's signature on the *other* cuvette's signal, not the sample's own
gas.

## 3 · The species test: HOO⁻, not H₂O₂ in general

If the effect tracked total peroxide regardless of protonation state,
`[enz]/[H₂O₂]` (un-weighted by pH) should predict the trough at least as
well as the anion-specific ratio:

| substrate | ρ against log[enz]/[HOO⁻] | ρ against log[enz]/[H₂O₂] |
|---|---|---|
| 4OMe-BnOH | **−0.621** | +0.160 (p = 0.05, **wrong sign**) |
| BnOH | **−0.315** | −0.234 (same sign, weaker) |

In 4OMe-BnOH the total-peroxide ratio is not even correctly signed, while
the pH-weighted anion ratio is the strongest correlation in the whole
analysis. `[HOO⁻]` and `[H₂O₂]` are not fully independent — `[HOO⁻] =
[H₂O₂]·f(pH)` — but weighting by pH should not matter if it were really
total peroxide being consumed, and it clearly does. That is consistent with
**HOO⁻'s much greater nucleophilicity toward a carbonyl** (the α-effect)
than the neutral molecule — the same chemistry a ketone catalyst forming a
Criegee-type peroxide adduct would be expected to run through.

## 4 · A second, smaller cluster: the substrate itself, at its scarcest

**14 of 17** genuine curves are oxidant-dominated. The remaining **3** —
exps **141.4**, **142.4**, **143.4**, all BnOH at the two-axis block's lowest
substrate rung, 0.216 mM — are substrate-dominated instead: the pH there
(9.15–9.73) is high enough that [HOO⁻] is abundant and `[enz]/[HOO⁻]` is
small, 0.022–0.041, smaller than `[enz]/[S]` (6.5–9.7%) for the first time in
the genuine set. Every other genuine curve has `[enz]/[HOO⁻]` **above 95** —
more than three orders of magnitude past these three's ceiling of **0.041** —
which is the one composition in the archive where [enz] is instead a large
enough fraction of [S] for simple stoichiometric substrate binding to
plausibly starve the shared background reaction. It is why [enz]/[S] carries
no *overall* signal in Section 2 — three curves cannot move a correlation
dominated by the other fourteen — without the cluster being any less real on
its own.

## 5 · Not the O2 artefact

A growing bubble raises absorbance (`BUBBLES.md` §1, light scattered out of
the detector), so it cannot produce a dip on its own, and none of the
seventeen troughs overlaps a detected detachment (Section 1). The sharper
test is what `debubble` correction does to the two curves that carry any O2
event at all:

| curve | events | raw trough | corrected trough |
|---|---|---|---|
| exp 141.4 | 3 | −13.55σ | **−30.65σ** |
| exp 142.4 | 7 | −25.58σ | **−36.89σ** |

Correction **deepens** both troughs; it does not weaken either. The O2
correction had been partially masking the dip with later bubble-driven rise,
not causing it.

## 6 · How fast — a rate constant from the fit already computed

A trough IS the "lag" shape of `summary_kinetics.fit_progress`'s one/two-phase
functional form taken far enough that `B/τ` exceeds `v_ss` — the same fit
already computed for every curve in the archive, so `τ_fast` (`scope.frame`'s
own column) is already the relaxation time of whatever produces the trough,
with no new fitting needed.

Pseudo-first-order: `k_obs = 1/τ_fast = k_on·[excess reagent]`, where the
excess reagent is `[enz]` for the oxidant cluster (ratio 95–4263, a clean
excess) and `[S]` for the substrate cluster (ratio only ~10–15, a rougher
approximation). Across all 17 genuine curves:

- **k_on ranges 0.51 to 32.7 M⁻¹s⁻¹ — a spread of only 1.8 orders of
  magnitude**, despite [enz] varying 20-fold, pH varying 5.47–9.73,
  temperature varying 15–40 °C, both clusters, both substrates and two
  buffers all pooled together.
- That is **8–9 orders of magnitude below the diffusion limit**
  (~10⁹–10¹⁰ M⁻¹s⁻¹) — the regime of a chemically-controlled, bond-forming
  step, not a barrierless encounter.

## 7 · Retracting a two-point activation energy

Exps 19.1 (15 °C) and 34.4 (40 °C) are an almost-matched pair — identical
pH (7.00), identical [enz] (0.241 mM), identical [H₂O₂] (82.5 mM), only
temperature and [S] differ. The obvious two-point Arrhenius estimate from
them: **Ea = 85.7 kJ/mol**.

That number does not survive its own self-consistency check. Solving for
the implied pre-exponential factor from the same two points gives
**2.1 × 10¹⁵ M⁻¹s⁻¹** — about 2 × 10⁵ times the diffusion limit — and the
equivalent Eyring ΔS‡ is **+40 J/mol/K**, where a genuine bimolecular
association (losing translational and rotational freedom to form one
ordered transition state) must be **negative**. A positive activation
entropy for a simple two-body association is not physically possible.

The honest version uses every available point, not the two extremes.
`arrhenius.arrhenius_fit` (the same regression and the same standard error
every temperature-series result in this project already uses) applied to
all 14 oxidant-cluster curves and their own actual temperatures gives:

| quantity | value |
|---|---|
| Ea | **78.9 ± 48.6 kJ/mol** |
| t-statistic on the slope | **1.6** (below the ~2 needed for significance) |

That is not the same estimate tightened — it is a demonstration that the
slope is not resolved. Among the **12** curves that happen to share the
single most common temperature (298.15 K), `k_on` alone spans **1.43**
orders of magnitude — almost as wide as the **1.74** spanned across the
*entire* 15–40 °C range. Most of what looked like a temperature trend in the
two-point comparison is ordinary curve-to-curve scatter a single point per
temperature cannot distinguish from a real one. **The 86 kJ/mol figure is
retracted**, and the archive does not support a replacement number — only
the (honest, wide) statement that it cannot resolve one.

## 8 · The buffer leaves a footprint the simple story doesn't predict

If the trough were a clean elementary step between the catalyst and free
HOO⁻, `k_on` should not depend on which buffer holds the pH, once `[enz]`
is accounted for. It does, among the oxidant-cluster curves:

| buffer | n | geometric mean k_on | range |
|---|---|---|---|
| Phosphate | 11 | 2.85 M⁻¹s⁻¹ | 0.59 – 18.0 |
| Pyrophosphate | 3 | **17.75 M⁻¹s⁻¹** | 7.88 – 32.7 |

Pyrophosphate runs roughly **6×** phosphate's geometric mean, even after
normalising by `[enz]` identically throughout. The likeliest reading
connects to something already established elsewhere in this project:
`induction.joint_buffer_order` finds the catalyst's own E → E* activation
step satisfies the pre-equilibrium "+1" rule specifically on the **buffer**
axis (+1.094 ± 0.150), not the peroxide axis. A two-step picture fits both
findings without contradiction: a fast, cuvette-symmetric buffer–HOO⁻
pre-equilibrium (`buffer + HOO⁻ ⇌ buffer–OOH`) sets how much reactive
material is available — a fast, symmetric equilibrium present identically
in both cuvettes cancels in the reference subtraction on its own, so it
does not by itself cause the trough — while the catalyst's own, slower,
genuinely asymmetric engagement with whatever sits in that pool is the
trough itself, and the already-established buffer-driven step converts the
loaded intermediate into the active catalyst afterward.

**Only the first and third pieces have real statistical power behind
them in this archive** (this section's `[enz]/[HOO⁻]` result, and
`induction/`'s own buffer-axis result). The middle piece — that a buffer
runs its own fast HOO⁻ equilibrium, and that different buffers hold
different-sized reservoirs — is inferred from the buffer-identity effect
above, not directly measured; an attempt to check it against the curves'
own slower relaxation (`tau_slow`) came back with only 6 of 14
oxidant-cluster curves carrying a resolved value, too few to test anything.

## 9 · The seventeen curves

| curve | substrate | pH | [S] mM | [H₂O₂] mM | [enz] mM | [enz]/[S] | [enz]/[HOO⁻] | trough, σ | cluster |
|---|---|---|---|---|---|---|---|---|---|
| exp 4.2 | 4OMe-BnOH | 6.71 | 3.228 | 82.5 | 0.175 | 0.054 | 141 | **−41.5** | oxidant |
| exp 5.2 | 4OMe-BnOH | 6.71 | 3.242 | 82.5 | 0.175 | 0.054 | 141 | −40.1 | oxidant |
| exp 151.7 | BnOH | 5.47 | 10.816 | 5.14 | 0.021 | 0.002 | 4263 | −28.7 | oxidant |
| exp 142.4 | BnOH | 9.43 | 0.216 | 73.424 | 0.014 | 0.065 | 0.022 | −25.6 | substrate |
| exp 143.4 | BnOH | 9.73 | 0.216 | 35.244 | 0.021 | 0.097 | 0.035 | −22.1 | substrate |
| exp 5.3 | 4OMe-BnOH | 6.71 | 4.864 | 82.5 | 0.175 | 0.036 | 144 | −15.5 | oxidant |
| exp 7.2 | 4OMe-BnOH | 6.71 | 3.011 | 82.5 | 0.175 | 0.058 | 141 | −14.4 | oxidant |
| exp 141.4 | BnOH | 9.15 | 0.216 | 73.424 | 0.014 | 0.065 | 0.041 | −13.6 | substrate |
| exp 5.1 | 4OMe-BnOH | 6.71 | 1.621 | 82.5 | 0.175 | 0.108 | 139 | −11.8 | oxidant |
| exp 34.4 | 4OMe-BnOH | 7.00 | 8.251 | 82.5 | 0.241 | 0.029 | 140 | −11.5 | oxidant |
| exp 151.5 | BnOH | 5.47 | 10.816 | 24.964 | 0.021 | 0.002 | 878 | −10.6 | oxidant |
| exp 19.1 | 4OMe-BnOH | 7.00 | 1.85 | 82.5 | 0.241 | 0.130 | 96 | −8.2 | oxidant |
| exp 10.1 | 4OMe-BnOH | 5.87 | 1.85 | 82.5 | 0.273 | 0.148 | 1564 | −8.2 | oxidant |
| exp 22.1 | 4OMe-BnOH | 6.50 | 1.85 | 82.5 | 0.241 | 0.130 | 314 | −8.2 | oxidant |
| exp 4.1 | 4OMe-BnOH | 6.71 | 1.614 | 82.5 | 0.175 | 0.108 | 139 | −7.3 | oxidant |
| exp 127.3 | 4OMe-BnOH | 6.94 | 9.47 | 3.879 | 0.032 | 0.003 | 331 | −6.5 | oxidant |
| exp 22.2 | 4OMe-BnOH | 6.50 | 3.7 | 82.5 | 0.241 | 0.065 | 321 | −6.2 | oxidant |

Exps 4, 5, 7 are three of the archive's four-fold `REPLICATE_RUNS` — the
same composition run four times — and carry the three deepest, cleanest
examples in the whole scan.

## What this does and does not establish

**Established.** The dip is a real, sustained, multi-reading decline in the
raw readings — first seen by eye through this session's derivative-of-fit
panel on exps 135.5 and 151.5–.7 — not a fitted-curve extrapolation artefact
and not the O2 artefact. Archive-wide it tracks `[enz]/[HOO⁻]` specifically
— not total `[H₂O₂]`, not `[S]` — far more strongly than either alternative,
in both substrates independently. Its own relaxation time gives a
self-consistent pseudo-first-order rate constant (0.5–33 M⁻¹s⁻¹, 1.8 orders
of magnitude, well below the diffusion limit) with no new fitting required,
and that rate depends on buffer identity in a way a clean elementary step
should not.

**Not established.** Whether the catalyst is engaging HOO⁻ directly or via
a buffer-perhydrate intermediate — this archive has no headspace or
manometric measurement of anything, the same limit `BUBBLES.md` states for
the gas itself. No activation energy: the two-point estimate (86 kJ/mol) is
retracted on its own Eyring self-consistency check (an unphysical
pre-exponential factor, a positive activation entropy), and the honest
14-point regression (Section 7) is not statistically significant. And exp
151 (three of the strongest oxidant-cluster curves) is one of the two-axis
block's own weakest, most drift-dominated runs (`scope.concentration_agreement`'s
screen does not even resolve it) — real chemistry there is hardest to pull
apart from the cell's own wander, which is exactly why the 4OMe
`REPLICATE_RUNS` curves matter: a different substrate, a different buffer,
none of the two-axis block's own caveats, and the three deepest troughs in
the archive.

## Reproducing

```
python data/test_early_trough.py
python early_trough/build_figures.py
python early_trough/check_numbers.py
```

`check_numbers.py` re-derives every number above from `early_trough.py` and
fails if the prose and the code disagree.
