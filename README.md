# Chemzyme kinetics

Master's thesis data pipeline: oxidation of benzyl alcohol and
4-methoxybenzyl alcohol by H<sub>2</sub>O<sub>2</sub>, catalysed by a
cyclodextrin–ketone "chemzyme".

The experiments were run in 2010. This repository takes the surviving archive —
instrument exports and the recipe spreadsheets that go with them — and turns it
into a dataset that can carry a kinetic model, with every concentration
traceable to something that was weighed or measured.

## Layout

```
data/Mads/                  the archive as delivered, in its original folders
data/data/                  flat copy of what the pipeline reads (.txt + .xls)
Mads-...-001.zip            the pristine 2024-12-07 delivery, kept as an anchor

data/kinetics_io.py         sheet + instrument-export reader; builds the dataset
data/build_manifest.py      declared ground truth: rulings, exclusions, provenance
data/solution_chemistry.py  ionic strength, pKa activity correction, [HOO-]

data/validate_dataset.py    the gate: run this before trusting anything
data/verify_*.py            four independent cross-checks (see below)
data/test_*.py              three test suites, including fault injection
data/build_dossier.py       one HTML page per experiment, for review by eye

MECHANISM.md                the 7-step mechanism, its reduction, and the evidence
COMPUTATIONAL.md            quantum-chemistry task register (C1 pending)
DATA_VERIFICATION.md        dated log of every check and every ruling
computational/hellowater/   ORCA smoke test, proves the toolchain runs
```

## Running the checks

From the repository root:

```bash
python data/validate_dataset.py           # fast: metadata, optics, exclusions
python data/validate_dataset.py --deep    # adds the four independent chains
python data/test_validator.py             # fault injection: corrupt, expect a catch
python data/test_solution_chemistry.py
python data/test_curve_flags.py
python data/build_manifest.py --write     # rebuild data/manifest.csv
python data/build_dossier.py              # rebuild dossier.html
```

Current state:

```
compiled   454 rows / 100 experiments      data/experiment_data.csv
fittable   404 rows /  88 experiments      after clean_experiment_dataframe
manifest   89 use, 11 excluded, 15 rulings, 0 open questions
checks     0 errors; 17/17 fault injection; all suites pass
```

## Why there are four verification modules

A check is only worth its runtime if its source is independent of the
extraction's. Reading a value back from the cell it was written from confirms
nothing — a mistake made once here, and the reason this is stated as a rule.

| module | independent source |
|---|---|
| `verify_enzyme.py` | the weighed mass of catalyst, through the recorded stock and cuvette volumes |
| `verify_rate_workbook.py` | `Rate(pH).xls`, the experimenter's own analysis, written for another purpose |
| `verify_dilutions.py` | the recorded dilution series, traced back to weighed grams |
| `verify_buffer.py` | each buffer stock recovered from the sheet's own cuvette volumes |

`data/experiment_data.csv` is fully reproducible from `data/data`: rebuilding it
gives zero differing cells. Rebuilds are therefore safe, and any hand-patch has
to be encoded in `kinetics_io.EXPERIMENT_CORRECTIONS` rather than applied to the
file.

## Where the judgement lives

Metadata that could not be settled by rule was settled by hand, and every such
decision is recorded twice: as a `RULINGS` or `KNOWN_EXCLUSIONS` entry in
`data/build_manifest.py` with its reasoning, and as a dated entry in
`DATA_VERIFICATION.md` with the evidence. Nothing is corrected silently.

The standing precedence is **sheet over filename**. A filename gets copied
forward between runs and only partly updated; a declared sheet value is what was
measured. Three experiments were caught this way (75, 76, 78), and the rule is
enforced by `verify_buffer.py`.

## Known limitations

- **Activity model.** `[HOO-]` uses Davies. 70% of live rows exceed the 100 mM
  where extended Debye–Hückel holds, and Davies itself is bounded rather than
  physical above ~500 mM. Pitzer is the principled fix.
- **Mixed buffers.** Exps 135–151 are two-salt pyrophosphate/phosphate mixtures
  treated as pure pyrophosphate.
- **Exps 75, 76** name their buffer hexametaphosphate but are given
  pyrophosphate pK<sub>a</sub> values.
- **No ODE fitter yet.** `MECHANISM.md` reduces the 7-step system to 6 states
  with three conservation laws and works out the observation equation; the
  implementation does not exist.
