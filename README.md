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
data/verify_*.py            five independent cross-checks (see below)
data/read_rre.py            reads the instrument binaries the .txt exports came from
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
python data/validate_dataset.py --deep    # adds the five independent chains
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
python data/fit_kinetics.py --list                    # blocks in scope (exps 135-151)
python data/fit_kinetics.py --buffer Pyrophosphate    # fit the scope
python data/fit_kinetics.py --scope all --list        # every block, ignoring the scope
python data/fit_kinetics.py --scope all --substrate BnOH --temperature 25 --buffer Phosphate
python data/plot_fit.py data/fits/BnOH_25C_Phosphate.json   # -> figures/
```

`plot_fit.py` refits nothing, so it is seconds. It writes three figures: the
enzyme-free curves against the model, the catalysed ones, and a six-panel
diagnostic showing *how* the fit fails — parity, residuals in units of each
curve's own noise, lag position, the species the model proposes, and the
reaction order in substrate.

Current state:

```
readings   402 rre /   0 txt curves        every cuvette from the instrument's own file
compiled   454 rows / 100 experiments      data/experiment_data.csv
fittable   402 rows /  88 experiments      clean_experiment_dataframe, less 25,2 and 25,4
fit scope  119 curves / 17 experiments     exps 135-151; fit_dataset.PRIMARY_SCOPE
manifest   89 use, 11 excluded, 19 ruled experiments, 0 open questions, 0 conflicts
checks     0 errors, 23 warnings, 15 notes; 20/20 fault injection
suites     test_kinetic_model 29/29; test_fit_kinetics 59/59;
           test_summary_kinetics 73/73; test_curve_screen 42/42;
           all others pass
```

## What the fitting is scoped to

Fitting is scoped to **exps 135-151** — 119 curves, BnOH / 25 °C /
pyrophosphate. They are the only runs in the archive that vary *both* the
substrate and the peroxide inside a single run (100.0% of the scope's log[S]
variance and 94.1% of its log[H₂O₂] variance is within-experiment), and they
span 19 pH values from 5.47 to 9.73 — 5.1 decades of [HOO⁻] — in one block,
with no exclusions and no open questions. `FITTING.md` sets out the full case,
the two conditions that must be met before a fit here is quotable, and the one
thing the scope costs: it holds no enzyme-free curves at all.

The scope lives in `fit_dataset.PRIMARY_SCOPE` and is re-derived from the
designs by `test_fit_kinetics.test_scope`, which fails if any run outside it
ever turns out to carry the same two-axis design. It is deliberately **not**
the hand-sorted `data/Mads/good data BnOH/` folder, which also contains an
already-excluded run, a run from a different cell, and a sheet with no
instrument data — see `DATA_VERIFICATION.md`, 2026-08-31.

A full sequential fit takes roughly 30 minutes: the model is integrated once per
curve per residual evaluation, and the optimiser is deliberately started from
many points (see `fit_kinetics._guaranteed_points`).

## Why there are four verification modules

A check is only worth its runtime if its source is independent of the
extraction's. Reading a value back from the cell it was written from confirms
nothing — a mistake made once here, and the reason this is stated as a rule.

| module | independent source |
|---|---|
| `verify_enzyme.py` | the weighed mass of catalyst, through the recorded stock and cuvette volumes — and, separately, the *shape* of the cuvette table, which says whether an experiment had catalyst at all without reading a single concentration |
| `verify_rate_workbook.py` | `Rate(pH).xls`, the experimenter's own analysis, written for another purpose |
| `verify_dilutions.py` | the recorded dilution series, traced back to weighed grams |
| `verify_buffer.py` | each buffer stock recovered from the sheet's own cuvette volumes |
| `verify_instrument.py` | the instrument's own `Substrate Conc.` export header — the only concentration record that is not the workbook, so the only one that can catch a workbook copied forward from another run |

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

Inverting that precedence is how the worst error in this log happened: exps 32
and 34–37 were forced to `[enz] = 0` on 2026-08-30 because their filenames say
`with_NO_E`, and were ruled catalysed on 2026-08-31 when the sheets were read
properly. Twenty catalysed curves had been sitting in the enzyme-free set. The
check that now prevents it reads the cuvette table's *layout* rather than its
numbers — rows 5–8 are the reference channel, and what the reference omits says
whether the experiment had catalyst.

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
