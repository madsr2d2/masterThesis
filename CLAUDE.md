# Working in this repository

Master's thesis pipeline: oxidation of benzyl alcohol and 4-methoxybenzyl
alcohol by H<sub>2</sub>O<sub>2</sub>, catalysed by a cyclodextrin–ketone
"chemzyme". `README.md` has the layout; this file has the rules.

## The scope

Fitting and analysis are scoped to **exps 135-151** (119 curves, BnOH / 25 °C /
pyrophosphate), defined by `fit_dataset.PRIMARY_SCOPE` and checked by
`test_fit_kinetics.test_scope`. The scope is **not** the hand-sorted
`data/Mads/good data BnOH/` folder, which also holds an excluded run, a run
from a different cell, and a sheet with no instrument data.

For any analysis of the in-scope data, invoke the **`analyse-scope`** skill.

## Do not re-derive measurements

Curve measurements live in `data/curve_metrics.py`; in-scope selection and
derived columns live in `data/scope.py`. Import them.

Six functions were once defined in two modules each and four had diverged. The
lag statistic's copies disagreed on 96 of 402 curves — one would have reported
21% where the evidence said 34%. `data/test_curve_metrics.py` fails if a name
is defined at top level in two modules.

Readings come from the instrument's own `.rre` where one exists and from the
`.txt` export otherwise; `Curve.source` says which, and the floor differs by a
factor of 1096 between them. Never floor a noise **or a residual variance** at
`QUANTISATION_SIGMA` without checking the source: call
`fit_dataset.source_floor(curve.source)` and pass the result. Everything in
`curve_metrics` that floors takes the value as an argument, and the default is
the `.txt` export's.

This is not only about the noise column. `line_fit` floors the variance behind
every standard error in the package, and `acceleration` divides by two of them,
so a floor left at the export's quantisation on `.rre` data suppresses the
z-score it is measured by — that is how the in-scope acceleration count read
48/110 until 2026-09-01 when the instrument's own readings say 51/110.

If you need a quantity that does not exist, add it to `scope.py` or
`curve_metrics.py`. Never compute it inline in a script, including a throwaway
one: throwaway numbers end up in documents, which is how 98.4% / 82.4% reached
four files when the scoped figures are 100.0% / 94.1%.

## Precedence and provenance

- **Sheet over filename, and sheet over the instrument's header.** A declared
  sheet value beats a filename. Filenames get copied forward between runs and
  only partly updated. Inverting this is how the worst error in
  `DATA_VERIFICATION.md` happened. The `.txt` exports' `Substrate Conc.` header
  copies forward the same way — exp 72's is a truncated copy of exp 71's — so
  it loses to the sheet too. It is still worth reading: it is the only
  concentration record that is not the workbook, it corroborates 60 of 75
  experiments, and it is what caught exps 69 and 70
  (`data/verify_instrument.py`).
- Every judgement call is recorded **twice**: as a `RULINGS` or
  `KNOWN_EXCLUSIONS` entry in `data/build_manifest.py` with its reasoning, and
  as a dated entry in `DATA_VERIFICATION.md` with the evidence. Nothing is
  corrected silently.
- A check is only worth its runtime if its source is **independent** of the
  extraction's. Reading a value back from the cell it was written from confirms
  nothing.
- Hand patches go in `kinetics_io.EXPERIMENT_CORRECTIONS`, never into
  `data/experiment_data.csv` — that file is rebuilt from `data/data`.

## Before trusting anything

```bash
python data/validate_dataset.py --deep    # 0 errors expected
python data/test_curve_metrics.py         # duplicate guard + the lag statistic
python data/test_fit_kinetics.py          # selection, scope, parameter recovery
python data/test_validator.py             # fault injection
python data/test_slowdown.py              # the slowdown models and their regressions
```

Units: concentrations mM, time s.

## The temperature series

`temperature_series/` holds the archive's only temperature block (exps 14-19,
4OMe-BnOH / pH 7.00 / phosphate, 15-40 °C) — the sole route to activation
parameters. `data/arrhenius.py` has the fits; `scope.TEMPERATURE_SERIES` the
scope. Use `vmax`, never `v0`: 20 of its 24 curves accelerate past 3σ, so an
initial rate there is the induction rate.

`data/verify_enzyme_stock.py` recomputes every experiment's `[enz]` from the
weighing recorded beside it — an independent source, the way
`verify_instrument.py` is for concentrations.

## What the product does

`product_fate/` holds why every catalysed 4OMe curve rises to a maximum rate and
then falls: the rate declines **linearly in the product it has made**, while the
same chemistry with no enzyme declines on a clock instead. It is its own folder
and not part of `temperature_series/` because the discrimination needs 84 curves
across the whole 4OMe archive; the temperature series' 23 cannot make it.

`data/slowdown.py` has the machinery. Use `deceleration_drivers` before
asserting that anything in this archive slowed down "over time" — one progress
curve cannot tell time from product, because inside a curve the product only
grows with time, and the separation exists only across curves.
