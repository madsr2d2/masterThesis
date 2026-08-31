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

data/kinetic_model.py       the reduced mechanism as ODEs; no I/O, no fitting
data/fit_dataset.py         which curves are fittable, and what each one is
data/fit_kinetics.py        the sequential fit: enzyme-free first, then catalysed
data/plot_fit.py            draws a saved fit against the curves it was fitted to
data/fits/*.json            saved fit results (a fit costs ~30 min; these do not)

data/validate_dataset.py    the gate: run this before trusting anything
data/verify_*.py            four independent cross-checks (see below)
data/test_*.py              five test suites, including fault injection
data/build_dossier.py       one HTML page per experiment, for review by eye

MECHANISM.md                the 7-step mechanism, its reduction, and the evidence
FITTING.md                  what has been fitted, and what the fits established
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
python data/test_kinetic_model.py
python data/test_fit_kinetics.py
python data/build_manifest.py --write     # rebuild data/manifest.csv
python data/build_dossier.py              # rebuild dossier.html
```

Fitting:

```bash
python data/fit_dataset.py                # what is fittable, by block
python data/fit_kinetics.py --list        # which blocks support both stages
python data/fit_kinetics.py --substrate BnOH --temperature 25 --buffer Phosphate
python data/fit_kinetics.py --substrate BnOH --profile-r    # profile r instead of fitting it
python data/plot_fit.py data/fits/BnOH_25C_Phosphate.json   # -> figures/
```

`plot_fit.py` refits nothing, so it is seconds. It writes three figures: the
enzyme-free curves against the model, the catalysed ones, and a six-panel
diagnostic showing *how* the fit fails — parity, residuals in units of each
curve's own noise, lag position, the species the model proposes, and the
reaction order in substrate.

Current state:

```
compiled   454 rows / 100 experiments      data/experiment_data.csv
fittable   404 rows /  88 experiments      after clean_experiment_dataframe
manifest   89 use, 11 excluded, 15 rulings, 0 open questions
checks     0 errors, 11 warnings, 9 notes; 17/17 fault injection
suites     test_kinetic_model 24/24; test_fit_kinetics 39/39; all others pass
```

A full sequential fit takes roughly 30 minutes: the model is integrated once per
curve per residual evaluation, and the optimiser is deliberately started from
many points (see `fit_kinetics._guaranteed_points`).

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
- **The reduced mechanism does not fit**, and the sharpest reason needs no fit
  at all: every rate term carries `[S]` linearly, so the model is first order in
  substrate by construction, while the data is roughly half order. On the fitted
  block it misfits at **20–24x the curves' own noise**. What this does and does
  not license saying — the reduction and the observation equation are deficient,
  the chemistry is not yet convicted — is set out in `FITTING.md`.
- **Two of the six rate constants are lower bounds, not values.** `k3` and `k5'`
  each stop affecting the observable once their step is no longer rate-limiting.
  Only `k_can`, `k0`, `r` and `k6` are determined.
- **Only two blocks support the sequential fit.** Rate constants may only be
  pooled within one (substrate, temperature, buffer) cell, and enzyme-free
  controls exist in the same cell as catalysed runs for just two of eleven:
  BnOH/25 C/phosphate and 4OMe-BnOH/40 C/phosphate. The largest catalysed
  block, 127 BnOH pyrophosphate curves, has no background data at all.
