# The temperature series: exps 14-19

4OMe-BnOH + H<sub>2</sub>O<sub>2</sub>, chemzyme-catalysed, in 65 mM phosphate
at pH 7.00 with 82.5 mM H<sub>2</sub>O<sub>2</sub>, run at **15, 20, 25, 30, 35
and 40 °C**. The same four-rung substrate ladder (1.850, 3.700, 5.549,
7.399 mM) in every one, so the block is four independent Arrhenius fits over
six temperatures, not one fit over 24 points.

**It is the only temperature series in the archive.** Nothing in the BnOH set
varies temperature at all, so this is the sole route to activation parameters
and the only kinetic quantity that can be compared against the ORCA barriers in
`COMPUTATIONAL.md`.

    scope.TEMPERATURE_SERIES        the six experiments, in temperature order
    data/arrhenius.py               the fits
    python data/verify_enzyme_stock.py --sequence

| exp | 14 | 15 | 16 | 17 | 18 | 19 |
|---|---|---|---|---|---|---|
| temperature | 25 °C | 35 °C | **40 °C** | 30 °C | 20 °C | 15 °C |
| `[enz]` mM | 0.273 | 0.273 | **0.241** | 0.273 | 0.273 | **0.241** |
| curves, live | 4 | 4 | 4 | 4 | 4 | 4 |

All 24 curves are live. **18 come from the instrument's own `.rre` and 6 from
the `.txt` export** — and the six are not scattered: they are the
5.549 mM rung in every one of the six runs. That rung is therefore the coarse
one throughout, its noise pinned at the export's quantisation (identical to six
figures across all six runs, which is the 0.001 AU lattice, not a coincidence).
`scope.frame` already floors each curve by its own source
(`fit_dataset.source_floor`), so the rates are treated correctly; it is worth
knowing because a per-rung conclusion drawn from 5.549 mM alone rests entirely
on export-rounded data. It is currently the *best*-behaved rung
(Arrhenius rms 0.073), so nothing turns on it yet.

## 1. The enzyme mismatch is real, not a typing mistake

Two of the six runs record an enzyme concentration **11.7% lower** than the
other four, and one of them — exp 16, the 40 °C run — sits *between* two runs at
the higher value. That had to be settled before any fit, because the two
readings call for opposite treatments: a real restock is divided out, a
transcription error is corrected.

**Three independent lines agree that the values are real.**

### The weighing chain, which the pipeline never reads

`[enz]` reaches `experiment_data.csv` from the sheet's `kuv` cell. Beside it
each sheet also records the preparation that produced that cell:

    mass / molar mass / volume        ->  stock concentration
    stock x Enz[ml] / Vol[ml]         ->  kuv

Four numbers written at four different moments, none of them the one the
pipeline reads. `data/verify_enzyme_stock.py` recomputes `kuv` from the mass and
the cuvette volumes for every experiment in the archive:

| stock | weighed | made up to | cuvette mM | experiments |
|---|---|---|---|---|
| 1 | 0.0122 g | 3.3 mL | 0.17533 | 2-9 |
| 2 | 0.023 g | 4 mL | 0.272695 | 10-15, **17, 18** |
| 3 | 0.0203 g | 4 mL | 0.240683 | **16**, 19-22, 32, 34-36 |
| 4 | 0.0228 g | 4 mL | 0.270324 | 37, 41-49 |
| 5 | 0.0236 g | 4 mL | 0.167885 | 50-59, 85 |
| 6 | 0.0236 g | 40 mL | 0.0279809 | 60-84 |
| 7 | 0.0027 g | 4 mL | 0.0320121 | 127-131 |

**62 of the 63 chains are self-consistent** and **63 of 63 compiled values equal
the sheet's own `kuv`**. The single break is exp 58, which is already excluded
for running backwards, and its `kuv` misses its own preparation by 0.5%.

So both values in question descend from a real weighing — 0.0203 g against
0.023 g — each consistent to six figures, and stock 3 is used by **nine**
experiments. A transcription error into one cell would break the chain. Neither
chain is broken.

### But the ordering is genuinely odd, and the workbooks cannot settle it

Stocks are made once and used until exhausted, so in experiment order the mass
should change rarely and never change back. **Exp 16 is the only experiment in
the entire campaign whose stock interrupts a run** (`--sequence`): stock 3
between exp 15 and exp 17, both stock 2.

That is exactly the shape a copied workbook leaves, and the file evidence is
ambiguous. Exps **16, 17, 18 and 20 are byte-identical in size (92672) and
structure** — one lineage — and it straddles the boundary: 16 and 20 carry
stock 3, 17 and 18 carry stock 2. Whichever direction the copying went, someone
edited an enzyme block; the documents cannot say who was edited into what.

The instrument files add nothing here. The `.rre` binaries store an internal
path, and it does catch two copies elsewhere in the archive — exp 003's file is
internally named `mads_t006.rre` and exp 029's is `mads_t023.rre` — but exp
016's is named `mads_t016.rre`. The `161008` field is a firmware build code,
identical in all 143 files, not a date.

### The kinetics decide it

Rate divided by enzyme concentration should fall on a straight line against
1/T. A concentration that is wrong for one run displaces that run's four points
together, off the line the other five temperatures define. So: divide by each
candidate and ask which is most collinear (`arrhenius.enzyme_hypotheses`).

On `vmax`, residual rms in ln units — a relative scatter:

| hypothesis | mean rms | worst rung | rungs won | E<sub>a</sub> kJ/mol |
|---|---|---|---|---|
| **as recorded** | **0.078** | **0.085** | **4 of 4** | 90.2 |
| exp 16 restocked to 0.273 | 0.108 | 0.115 | 0 of 4 | 87.6 |
| exps 16 and 19 both 0.273 | 0.129 | 0.143 | 0 of 4 | 90.4 |

**The recorded values win every rung.** On `v0_quad` as a sensitivity they win
every rung again — 0.236 against 0.263 and 0.275 — so the verdict does not
depend on the estimator, and the four rungs are independent fits rather than
one result counted four times.

Per run it is clearer still. Exp 16's mean distance from the line:

- as recorded (0.241): **−0.060** — about 6% low, inside the others' scatter
  (exp 15 is +0.050, exp 17 +0.025)
- forced to 0.273: **−0.121** — twice as far out, and further out on all four
  rungs

Forcing the higher concentration makes exp 16 a **worse** outlier, in the
direction a too-large divisor produces.

### Verdict

**Real. Use the recorded values and divide by them.** Exp 16 and exp 19 were run
with a genuinely different enzyme stock, and the fits below normalise by
`[enz]` per run.

**What this does not establish**: *why* exp 16 sits out of sequence. Most likely
it was re-run later with the 19-22 session and kept its design number, but
nothing in the archive records that, and it does not change the treatment. Note
also that the method is weakest against an error common to *both* stock-3 runs —
displacing 15 °C and 40 °C the same way tilts the line hardly at all — so
"exps 16 and 19 both wrong" is rejected less firmly than "exp 16 alone", though
it loses every rung too.

## 2. Activation energy, first pass

From `arrhenius.rung_fits()` on `vmax`, per substrate rung:

| `[S]` mM | E<sub>a</sub> kJ/mol | rms |
|---|---|---|
| 1.850 | 94.5 | 0.078 |
| 3.700 | 86.0 | 0.085 |
| 5.549 | 89.0 | 0.073 |
| 7.399 | 91.5 | 0.076 |

Mean **90.2 kJ/mol**, spread 8.5 across four rungs that should agree. This is a
first pass and is not yet a result: it is an *apparent* activation energy of
`vmax`, which is a composite of every step, and the background is subtracted out
of these curves so it is the catalysed increment's temperature dependence, not
the reaction's. The Eyring form, the enthalpy/entropy split, and the comparison
against `COMPUTATIONAL.md`'s barriers are not done here.

**A caution for whatever comes next.** 20 of these 24 curves clear the 3sigma
acceleration gate and 8 exceed z = +15, so `v0` is an *induction* rate
here and must not be used — `vmax` is the estimator throughout. And none of
these curves has been through `curve_metrics.segmented_fit`; exp 65 in the BnOH
background taught us that the start-versus-end statistics step straight over a
mid-run break.
