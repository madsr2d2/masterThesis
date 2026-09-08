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

Fifteen curves, across both substrates, do exactly this. `data/early_trough.py`
has the machinery; `progress_curves.html` shows every one of the fifteen with
its full fit, plus the two curves that looked like more of the same and are
rejected on inspection.

## 1 · A trough that survives three separate ways to be spurious

`curve_metrics.early_trough` smooths the first half of a curve's run with a
9-reading rolling mean and reports its minimum, in units of the curve's own
noise. Over every live catalysed curve in the archive (`scope.archive()`,
**311 curves**), **17** clear a smoothed trough of −4.0 sigma or deeper.

A smoothed mean alone cannot tell a genuine, sustained decline from one deep
outlier reading dragging it down. Three independent screens separate them:

- **Sustained.** At least 5 of the 9 readings inside the trough window must
  themselves sit below 2 sigma, not just the smoothed mean. This alone
  removes two candidates: exp **4.1** (smoothed z **−7.31**, only 3 of 9
  readings individually qualify) and exp **22.2** (z **−6.23**, only 4 of 9).
- **No bubble overlap.** The trough window must not overlap a detected O2
  detachment (`curve_metrics.detachments`). None of the archive's candidates
  do.
- **Survives debubble correction.** The trough recomputed on the
  `curve_metrics.debubble`-corrected curve must still clear −4.0 sigma.

**15 of 17** candidates pass all three (Figure A shows the excluded two
beside the fifteen that survive).

## 2 · The driver is the oxidant, not the substrate, in both substrates independently

![figure A](index.html#figure-a)

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

## 3 · A second, smaller cluster: the substrate itself, at its scarcest

![figure B](index.html#figure-b)

**12 of 15** genuine curves are oxidant-dominated. The remaining **3** —
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
dominated by the other twelve — without the cluster being any less real on
its own.

## 4 · Not the O2 artefact

![figure C](index.html#figure-c)

A growing bubble raises absorbance (`BUBBLES.md` §1, light scattered out of
the detector), so it cannot produce a dip on its own, and none of the fifteen
troughs overlaps a detected detachment (Section 1). The sharper test is what
`debubble` correction does to the two curves that carry any O2 event at all:

| curve | events | raw trough | corrected trough |
|---|---|---|---|
| exp 141.4 | 3 | −13.55σ | **−30.65σ** |
| exp 142.4 | 7 | −25.58σ | **−36.89σ** |

Correction **deepens** both troughs; it does not weaken either. The O2
correction had been partially masking the dip with later bubble-driven rise,
not causing it.

## 5 · The fifteen curves

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
| exp 127.3 | 4OMe-BnOH | 6.94 | 9.47 | 3.879 | 0.032 | 0.003 | 331 | −6.5 | oxidant |

Exps 4, 5, 7 are three of the archive's four-fold `REPLICATE_RUNS` — the
same composition run four times — and carry the three deepest, cleanest
examples in the whole scan.

## What this does and does not establish

**Established.** The dip is a real, sustained, multi-reading decline in the
raw readings — first seen by eye through this session's derivative-of-fit
panel on exps 135.5 and 151.5–.7 — not a fitted-curve extrapolation artefact
and not the O2 artefact, and archive-wide it tracks `[enz]/[HOO⁻]` far more
strongly than `[enz]/[S]`, in both substrates independently.

**Not established.** Whether the catalyst is engaging the peroxide, the
hydroperoxide anion specifically, or changing its own resting state some
other way — this archive has no headspace or manometric measurement of
anything, the same limit `BUBBLES.md` states for the gas itself. Nor can it
separate that from a background reaction that runs on something else the
catalyst also touches. And exp 151 (three of the strongest oxidant-cluster
curves) is one of the two-axis block's own weakest, most drift-dominated
runs, per `scope.concentration_agreement`'s screen — real chemistry there is
hardest to pull apart from the cell's own wander, which is exactly why the
4OMe `REPLICATE_RUNS` curves matter: a different substrate, a different
buffer, none of the two-axis block's own caveats, and the three deepest
troughs in the archive.

## Reproducing

```
python data/test_early_trough.py
python early_trough/build_figures.py
python early_trough/check_numbers.py
```

`check_numbers.py` re-derives every number above from `early_trough.py` and
fails if the prose and the code disagree.
