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
`.txt` export otherwise; `Curve.source` says which, and the noise floor differs
by a factor of 1000 between them. Never floor a noise at `QUANTISATION_SIGMA`
without checking the source — see `data/read_rre.py`.

If you need a quantity that does not exist, add it to `scope.py` or
`curve_metrics.py`. Never compute it inline in a script, including a throwaway
one: throwaway numbers end up in documents, which is how 98.4% / 82.4% reached
four files when the scoped figures are 100.0% / 94.1%.

## Precedence and provenance

- **Sheet over filename.** A declared sheet value beats a filename. Filenames
  get copied forward between runs and only partly updated. Inverting this is
  how the worst error in `DATA_VERIFICATION.md` happened.
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
```

Units: concentrations mM, time s.
