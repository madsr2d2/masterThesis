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

## 2. The breakpoint screen: the breaks are the induction, not an artefact

Run first, because exp 65 in the BnOH background showed that the start-versus-end
shape statistics step straight over a mid-run break, and every conclusion here
rests on a single rate per curve. `scope.synchronised_break`:

| T | 15 °C | 20 °C | 25 °C | 30 °C | 35 °C | 40 °C |
|---|---|---|---|---|---|---|
| slope after / before | 1.50-18.08 | 1.32-2.27 | 0.99-1.81 | 1.27-1.36 | 0.89-1.00 | 0.70-0.86 |
| cuvettes steepening | 3 of 4 | 2 of 4 | 1 of 4 | 0 of 4 | 0 of 4 | 0 of 4 |

**Nothing here looks like exp 65.** That run's four cuvettes broke within 56 s of
each other regardless of a 20-fold substrate range, which is the signature of
something that is not the reaction. Here the ratio falls with temperature — cold
runs accelerate throughout, hot runs decelerate — and a second, unrelated
statistic says the same thing: every run from 15 to 30 °C classifies as a **lag**
curve, with τ falling **6489 → 3666 → 945 → 876 s** as it warms, and 35 and 40 °C
flip to `burst`. Two independent statistics agreeing on an induction period that
shortens with heating is chemistry, not an artefact.

**The fall is not quite monotone, and the exception is instructive.** Median
ratios run 1.75, 1.56, 1.19, **1.33**, 0.91, 0.82 — 30 °C sits above 25 °C. The
break ratio does not track temperature directly, it tracks **how much of the run
the induction occupies**, and the runs are not the same length: measured in units
of its own τ, the 25 °C run is **19τ** long and the 30 °C run only **5.8τ**. The
quantity that *is* monotone in temperature is τ itself.

**A limitation of the statistic, seen here for the first time.** Exp 19's
lowest-substrate cuvette reports a ratio of **18.08**, far outside everything
else. It is a near-zero denominator: that curve is flat at 5.6e-08 AU/s for the
first 6076 s — the coldest temperature at the lowest substrate, the slowest
condition in the block — and then rises at 1.0e-06. The reading is "the
pre-break slope is indistinguishable from flat", not "the post-break slope is
remarkable". `break_ratio` is unstable wherever `slope_before` approaches zero,
and a lag curve at the slowest condition is exactly where that happens.

### But it exposed a real problem with `vmax`

`vmax_where` — where the steepest block sits, as a fraction of the run:

| T | 40 °C | 35 °C | 30 °C | 25 °C | 20 °C | 15 °C |
|---|---|---|---|---|---|---|
| `vmax_where` | 0.10-0.30 | 0.29-0.49 | 0.48-0.87 | 0.30-0.50 | 0.70-0.90 | 0.70-0.90 |

At 15 and 20 °C the rate is **still rising when the run ends**, so `vmax` there
is not a maximum — it is wherever the measurement stopped. Under-reading the
cold end tilts the Arrhenius line and inflates E<sub>a</sub>.

**The two estimators fail at opposite ends**, which is the awkward part. On a lag
curve the burst form's `v_ss` is the rate `vmax` is trying to reach, and at
15-30 °C that fit is sound (tau resolved on 15 of 16 curves, residuals 0.94-1.64x
noise). At 35–40 °C it degenerates — τ unresolved on 7 of 8, the form flipping
to `burst`, tau unresolved on 7 of 8, because a decelerating curve has no lag to
measure — so `v_ss` there
is a meaningless late rate, and `vmax` is the good one.

**How big is the bias?** `arrhenius.truncation_sensitivity` substitutes `v_ss`
for `vmax` at 15 and 20 °C only, purely to size it:

| | E<sub>a</sub> kJ/mol |
|---|---|
| `vmax` throughout | 90.2 |
| `v_ss` at the cold end | 88.3 |
| **inflation from truncation** | **1.9** |

`v_ss`/`vmax` is 1.042 at 15 C and 1.065 at 20 C — a 4–6.5% under-read — and
0.98, 0.99 at 25 and 30 C where the two should agree, which is the check that
the substitution means anything. (It is 0.67 and 0.29 at 35 and 40 °C, the
degenerate end, and is not used there.)

**So truncation costs about 1.9 kJ/mol, against a spread of 8.5 across the four
substrate rungs that should agree.** It is real, it is in a known direction, and
it is *not* the limiting problem. Mixing estimators is normally the error rather
than the fix, and it is used here only to size the bias, never for a headline.

## 3. Activation energy, first pass

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
