---
name: analyse-scope
description: Use when analysing, measuring, plotting or answering any question about the in-scope kinetics experiments (exps 135-151) in this thesis repo — including rates, noise, lag, autocatalysis, substrate/peroxide/pH order, curve quality, or "how many curves do X". Load before writing any analysis script against data/.
---

# Analysing the in-scope experiments

Fitting and analysis are scoped to **exps 135-151** — 119 curves, BnOH / 25 °C /
pyrophosphate. See `FITTING.md` for why, and `DATA_VERIFICATION.md` (2026-08-31)
for the evidence.

## The rule

**Import from `data/scope.py` and `data/curve_metrics.py`. Do not re-derive.**

Every measurement that describes a curve — its noise, its initial rate, where
its slope peaks — has exactly one definition, in `curve_metrics`. Every
selection of the in-scope data has exactly one definition, in `scope`.

This is not a style preference. Six functions in this repo were once defined
twice and four had diverged; the lag statistic's two copies disagreed on 96 of
402 curves and one of them would have put 21% in a thesis where the evidence
said 34%. `data/test_curve_metrics.py::test_no_duplicate_definitions` now fails
if a duplicate reappears.

That archive figure is now **37.6%** (151/402), not because the statistic
changed but because the readings did: since 2026-08-31 they come from the
instrument's own `.rre` files rather than the 0.001 AU `.txt` exports. Use
`curve.noise`, never a fresh `curve_noise(values)` — the floor depends on
`curve.source` and `frame()` already applies the right one.

The same failure has already happened in prose: the within-experiment contrast
figures were computed once in a throwaway script over the wrong curve set and
written into four documents as 98.4% / 82.4%, when over the actual scope they
are 100.0% / 94.1%.

## Start here

```python
import sys; sys.path.insert(0, "data")
from scope import frame, design, curves, ladder, within_experiment_share

df = frame()      # one row per in-scope curve, every derived column filled
```

`frame()` columns: `experiment sample pH s0 h2o2 e0 hoo duration_s points
noise net live v0 v0_stderr v0_rms peak lags late_over_early`.

Most questions are a groupby on that frame, not a new script:

```python
df[df.live].groupby("experiment").v0.median()        # rate by run
df[(df.pH >= 9) & df.live].lags.mean()               # lag fraction at high pH
design()                                             # per-experiment design table
within_experiment_share("s0")                        # the number the scope rests on
```

From the shell: `python data/scope.py` and `python data/scope.py --design`.

## If the quantity you need is missing

Add it to `scope.py` (a selection or a derived column) or `curve_metrics.py`
(a measurement of curve shape). **Do not compute it inline in a script**, even
a throwaway one — throwaway numbers end up in documents. If it was worth
computing once it will be asked for again.

## The +/- chemzyme controls

The scope has no enzyme-free run — every one of exps 135-151 carries enzyme.
The nearest thing is `scope.PAIRED_CONTROLS` (exps 65/66, 67/68, 69+70/71):
the same ladder run with and without 0.028 mM chemzyme, in phosphate and boric
buffer at pH 8.0-8.5. `python data/scope.py --controls`.

Use them for interpretation, never for a fit or a pooled constant -- different
buffer, and one `[H2O2]` per run. And when reading a substrate order from
enzyme-free runs, use `scope.FREE_BNOH` (exps 65, 67, 69, 70) and **not**
exps 3 and 6, which are buffer titrations whose rate falls with substrate for
a reason that is not substrate. `scope.FREE_BNOH_BUFFER_TITRATIONS` records it.

## Conventions that are already settled

- **Sheet over filename.** A declared sheet value beats a filename; filenames
  get copied forward between runs. Inverting this caused the worst error in
  `DATA_VERIFICATION.md`.
- **The scope is `fit_dataset.PRIMARY_SCOPE`, never a folder.**
  `data/Mads/good data BnOH/` is *not* the scope — it also holds exp 50
  (excluded), exp 51 (a different cell) and exp 134 (no instrument data).
- **Curve shape is never a defect.** Lags, dips and bursts are reported, not
  excluded.
- **There is no background run in this block.** Exps 150 and 151 were called
  the in-cell background; they are not. Every scoped run carries enzyme, and
  cuvette 1 of exps 150-151 sits within 2x of the [HOO⁻] trend the other runs
  define, so the reaction is still on — they are the bottom rung of the pH
  ladder, not a blank. Do not subtract them from anything.
- **Six runs move with something that is not the reaction.** In exps 136, 137,
  147, 149, 150 and 151, `scope.concentration_agreement` is 0.61 or below and
  as low as 0.005: their cuvettes' rates bear almost no relation to their
  cuvettes' concentrations, while exps 135, 138, 139, 140 and 142 run 0.93 to
  0.97. Their rates sit at the cell's own drift, a few times 1e-7 AU/s.
  Every conclusion in FITTING.md is *stronger* with them dropped, so use
  `concentration_agreement` before letting a weak run carry an argument.
- **Concentrations are mM, time is s**, throughout.
- Quote a rate constant only where it is identifiable; two of the six are lower
  bounds (`FITTING.md` F4).

## Before reporting a number

Numbers that reach a document must come from `scope`/`curve_metrics`, not from
a scratch calculation. If a number contradicts one already written down, say so
and check which selection each was computed over — that has been the cause both
times it has happened.
