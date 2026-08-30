# Data Verification Log

Record of checks performed on the raw kinetics dataset (`data/data/*.txt` + `*.xls`)
and its derived artifacts, and any corrections made as a result. New entries go at
the top. See `MECHANISM.md` for the chemistry and `COMPUTATIONAL.md` for pending
quantum-chemistry tasks.

---

## 2026-08-31 — Repository cleanup, and what the pristine zip caught

Housekeeping, with one finding worth recording.

### The delivered zip is not redundant, and it earns its place

`Mads-20241207T151327Z-001.zip` is the archive exactly as delivered on
2024-12-07. It looked like a 14 MB duplicate of `data/Mads`, so it was checked
before being considered for removal.

67 of its 716 files differ byte-wise from the working copies. A first comparison
put 35 of those down as real content differences — **wrong**, because comparing
`df.astype(str)` makes a float `nan` and an object-column `None` look different.
Compared null-aware, only **three** files differ in content, and all three are
stray keystrokes in cells nothing reads:

| file | difference |
|---|---|
| `mads_t002...xls` | R24CR holds a pasted list `1.5199,3.0397,4.5596,6.0795` |
| `mads_t014...xls` | R38CR holds `j` |
| `mads_t057...xls` | R1CA holds **`jjjj`** where the zip has the experiment number 57 |

The last one matters: exp 57's experiment number was typed over in the
`data/Mads` copy. The copy under `data/data` — the one the pipeline actually
reads — still holds 57, so nothing downstream was affected, and exp 57 is
excluded in any case.

**The zip stays.** It is the only way to detect edits like these, and deleting it
from the working tree would not shrink the repository anyway, since git history
holds it either way.

### Removed

Dead artefacts, all recoverable from history:

| | |
|---|---|
| `outputs/` | 60 directories, 42 MB of PySR symbolic-regression runs from 2024-12-15, referenced nowhere |
| `model.ipynb` | a Michaelis–Menten ODE prototype predating the 7-step mechanism. Both cells are broken — cell 0 unpacks 5 states from a 4-element `y0`, cell 1 returns 6 derivatives for 5 states |
| `test.ipynb` | one cell of spline scratch work |
| `Sample00*_time_series.csv` | four stray exports from 2024-12-14 |

### Moved and added

`hellowater.*` (the ORCA smoke test `COMPUTATIONAL.md` cites as proof the
toolchain runs) moved from the top level to `computational/hellowater/`, and the
reference updated.

`README.md` added — the repository had none. It carries the layout, the commands,
the current counts, the sheet-over-filename precedence rule, and the four known
limitations.

`data/movedata.py` gained a docstring. It is the one-off that flattened
`data/Mads` into `data/data`; it is imported by nothing and has no CLI, and was
kept because it records how `data/data` came to exist.

`.gitignore` gained `outputs/` and `.ipynb_checkpoints/`. Two LibreOffice lock
files were cleared from the working tree.

### After

The top level is now nine entries: three documents, the README, the notebook,
`data/`, `computational/`, the delivered zip, and the generated dossier. Every
module under `data/` is reachable, every check still passes, and
`experiment_data.csv` still reproduces from `data/data` with zero differing
cells.

---

## 2026-08-31 — The pH series is fixed, and [HOO⁻] now uses Davies

Two corrections before any fitting begins, both to quantities the model reads
directly.

### 1. Exps 127-131 now carry their own pH

Fixed at the root rather than as five per-experiment rulings.
`find_pH_value_in_range` now searches the more specific `buffer pH` label before
the bare `pH`, so the unused phosphate block can no longer win:

| exp | was | now | source |
|---|---|---|---|
| 127 | 7.29 | **6.94** | `buffer pH` |
| 128 | 7.29 | **7.41** | `buffer pH` |
| 129 | 7.29 | **8.11** | `buffer pH` |
| 130 | 7.29 | **8.56** | `buffer pH` |
| 131 | 7.29 | **8.98** | `buffer pH`, given as the range "8.88-9.07" |

Exp 131 states its buffer pH as a written range rather than a number, so
`_parse_pH_cell` takes a two-number cell at its midpoint. That is why the
earlier numeric scan found only four of the five.

Exactly **21 cells** changed in `experiment_data.csv`, all of them `pH` in these
five experiments. Every other value is byte-identical, and the manifest's
provenance for all five moved from `filename` to `sheet`. The open question is
closed; `KNOWN_OPEN_QUESTIONS` is empty again.

The choice was the buffer pH as prepared, not the post-run readings. The
post-run pair straddles it in every case, so the run-average is close either
way, and the prepared value is the one condition that is the same for all
cuvettes.

### 2. Davies replaces extended Debye-Hückel

`ACTIVITY_MODEL = "davies"`, resolved at call time so either treatment can be
selected: `effective_pka_h2o2(I, model="debye")` reproduces every earlier number
exactly.

The reason is coverage: 70% of live rows exceed the 100 mM where extended
Debye-Hückel is reliable, and the pyrophosphate runs reach 1069 mM. The effect
on the driving concentration is not small, and it is not random:

```
I (mM)        n    [HOO-] Davies / Debye
0-100       118          0.907
100-200     107          0.783
200-400      66          0.680
400-700      94          0.489
700-1200     13          0.248
```

Median 0.760 across 406 live rows; 62% change by more than 20%. **The error was
monotonic in ionic strength**, which is the shape that would have distorted the
buffer-concentration dependence specifically — the thing the buffer titrations
exist to measure.

Ionic strength itself is unchanged; only pKa_eff(H2O2) moves.

### Davies' own ceiling, stated rather than hidden

The -0.3I term eventually dominates, so the predicted shift turns around near
0.5 mol/L:

```
I (mM)      100     400     700    1069
pKa_eff  11.536  11.478  11.500  11.559
```

Above ~500 mM the curve is bounded rather than physical. What changed is that
the error stops growing without limit — it did not become correct.
`out_of_range_fraction` and the validator still report those rows, and a Pitzer
treatment with real ion-interaction parameters remains the principled fix.

`data/test_solution_chemistry.py` now pins **both** models against hand
calculations at I = 75 mM (Davies 11.554, Debye-Hückel 11.494), asserts they
agree at I = 0, asserts Davies shifts less wherever it matters, and requires an
unknown model to raise. The Δ(z²) = 2 regression is written against whichever
model is active, so switching cannot quietly retire it.

---

## 2026-08-31 — Exps 127-131 are a pH series recorded as a single pH

A re-audit of the metadata, done from the sheets rather than from the manifest,
found that **five live experiments carry the wrong pH**.

Exps 127–131 are a pyrophosphate pH series. Their sheets record:

| exp | buffer pH | pH after run (1,2,5,6) | pH after run (3,4,7,8) | manifest |
|---|---|---|---|---|
| 127 | 6.94 | 6.86 | 6.99 | **7.29** |
| 128 | 7.41 | 7.57 | 7.64 | **7.29** |
| 129 | 8.11 | 8.03 | 8.11 | **7.29** |
| 130 | 8.56 | 8.47 | 8.56 | **7.29** |
| 131 | — | 8.99 | 9.10 | **7.29** |

A series spanning roughly 6.9 to 9.1 is being modelled at a single point.

### Why it happened

Each of these sheets contains **two** buffer blocks:

```
I22 Buffer      J22 pyrophosphate      <- the one used
I23 mmol/l      J23 225.159
I24 buffer pH   J24 6.94

F24 Buffer      G24 Phosphate          <- an unused second block
F25 mol/l       G25 0.2
F27 pH          G27 7.29
```

`kinetics_io.find_pH_value_in_range` matches a **bare** `pH` label. `buffer pH`
is not a bare match, so the only cell it can see is F27 — which belongs to the
unused phosphate block and holds 7.29 in all five sheets.

The pyrophosphate block is demonstrably the buffer in use: 225.159 mM x the
0.86/1.0 volume ratio gives 193.637 mM, which is exactly the `[buf]` the dataset
carries for all five.

### The check that should have caught this was vacuous

The 2026-08-31 entry above reports "all 100 sheets declare pH, all agreeing with
the manifest to within 0.02". That comparison read the **same cell the pipeline
reads**, so it could only ever agree. A verification that consults the same
source as the thing it verifies confirms nothing.

The same trap applies elsewhere and is worth stating as a standing rule: a check
is only worth its runtime if its source is independent of the extraction's. The
`[enz]` chain (weighed catalyst, then `Rate(pH).xls`) and the buffer stock
recovery (cuvette volumes vs a declared column) are independent. A field read
back from the cell it was written from is not.

### Recorded, not fixed

Logged as an open question on all five experiments in
`build_manifest.KNOWN_OPEN_QUESTIONS`. It needs a ruling on which value to
adopt: the `buffer pH` as prepared, or the mean of the two post-run readings.
The post-run pair straddles the buffer value in every case, so the drift is
small — but choosing between them is a judgement, and the two differ by up to
0.16 units.

---

## 2026-08-31 — Correction: every sheet declares its pH, and the provenance column was hiding it

A status summary on 2026-08-30 said pH was "~92% filename-sourced" and called
that the dataset's largest remaining exposure. **That was wrong**, and the
manifest's own provenance column is what made it look true.

### All 100 sheets declare pH

Every sheet carries a bare `pH` label with the value in the cell to its right,
in one of four places:

| cell | sheets |
|---|---|
| K11 | 57 |
| N5 | 21 |
| I70 | 17 |
| F27 | 5 |

**All 100 agree with the manifest**, to within 0.02. The I70 group is the
135–151 series, which is why a first scan limited to the top 30 rows found only
83 — the sheets were not sparse, the scan was shallow.

### The pipeline was already reading them

`kinetics_io.find_pH_value_in_range` matches `\bpH\b` over a 100 x 100 window
and takes the neighbouring cell, falling back to the filename only if that
fails. The same is true of `find_temperature_value_in_range`, `find_buffer_type`
and `find_substrate_type`. So the sheet has always been the primary source for
all four fields.

What went wrong was the label. `build_manifest` stamped `provenance = filename`
whenever a filename agreed with the extraction, which read as "the filename is
where this came from" when it actually meant "the filename agreed with the
sheet". The strongest state the metadata reaches — two independent sources
saying the same thing — was recorded as if it were the weakest.

### The provenance column now distinguishes them

`read_sheet_claims()` calls each `find_*` helper with `filename=None`, isolating
what the sheet alone declares, and the resolution stamps `sheet+filename` where
both agree, `sheet` where only the sheet has it.

| field | sheet+filename | sheet only | filename only | folder | ruling |
|---|---|---|---|---|---|
| pH | 82 | 5 | **0** | 0 | 2 |
| substrate | 83 | 6 | **0** | 0 | 0 |
| T | 60 | 17 | 12 | 0 | 0 |
| buffer | 33 | 17 | 39 | 0 | 0 |
| has_enzyme | 0 | 0 | 42 | 24 | 23 |

(89 live experiments.)

**pH and substrate have no filename-only cases at all.** Every live experiment
has both from the sheet, and in 82 and 83 of them the filename independently
agrees. Buffer is the weakest of the four at 39 filename-only — and that is the
field where a stale filename has already been caught three times, in exps 75, 76
and 78.

Only the `provenance` column changed. `experiment_data.csv` still reproduces
from `data/data` with zero differing cells, and every check still passes.

---

## 2026-08-30 — A declared sheet value outranks a filename, and reading the sheets properly changed three things

**Ruling: where the sheet declares a value, that is ground truth; the filename
is not.** A filename is copied forward between runs and only partly updated,
which this archive does repeatedly and demonstrably — exp 85's sheet still says
experiment 83, exps 2–10 carry a stale 285 nm, exp 57's substrate block came
from t056. A declared value is what the experimenter measured and wrote down.

Acting on that ruling meant reading a part of the sheet the pipeline had never
looked at.

### Twenty-one sheets state the buffer stock outright

The header block carries a name and a molarity:

```
exp 75    N3 Buffer   O3 hexametaphosphate (pyrophosphate)
          N4 mol/l    O4 0.033

exp 79    N3 Buffer   O3 carbonate
          N4 mol/l    O4 0.1
```

The sheet's `[buffer]` column is computed from that cell. **In all 21 the two
agree** — the sheet is internally consistent, and only the filename dissents.

### So exps 75, 76 and 78 are settled

All three say `0.1M` in the filename over a sheet that says **0.033 M**, and the
`[buf]` column follows the sheet: 24.75 mM in the cuvette, not 75. The compiled
dataset already used the sheet, so no number changes; what changes is that this
is now adjudicated rather than open.

Exp 78 is the cleanest demonstration. Exp 79 is the same buffer, same pH 9.40,
same day, the same 1.55 ml of buffer in 2 ml — and declares 77.5 mM, exactly
what 0.1 M gives. Exp 78 declares 25.575, exactly what 0.033 M gives. Two
sheets, one condition, and each is internally consistent with its own stated
stock. Both filenames say `0.1M`.

`verify_buffer.py` now reports this as `bufstale` — a note, not a warning — and
gains `bufdeclared`, which errors if a sheet's stated stock and its own `[buf]`
column ever disagree. That check fires nowhere today.

### Exps 75 and 76 are hexametaphosphate, not pyrophosphate

Their sheets name the buffer `hexametaphosphate (pyrophosphate)`. A search of
every sheet in the archive finds that word in exactly these two files.

`solution_chemistry.py` classifies them as `Pyrophosphate` and applies
pK_a = [0.85, 1.96, 6.60, 9.41] with charges 0 to −4. Sodium hexametaphosphate,
(NaPO₃)₆, is a polymeric metaphosphate and not that species at all. Both
experiments are **live**, and I and [HOO⁻] depend on the choice.

**Open: which speciation to use for exps 75 and 76.** Not resolved here, since
it is a chemistry judgement rather than a data-provenance one.

### The later pyrophosphate series is genuine — and mixed

Exps 135–151 name their salts explicitly and carry gram-level preparations:

```
Buf 1.  Na4P2O7*10H2O   M.W. 446.06    30.1495 g in 500 ml  ->  135.18 mmol/l
Buf 2.  Na2HPO4*2H2O    M.W. 156.02    39.0052 g in 500 ml
```

(exps 143–151 use `NaH2PO4*2H2O` as Buf 2 instead.) So they are real sodium
pyrophosphate, correctly classified — but every one of them is a **two-salt
mixture of pyrophosphate and phosphate**, treated as pure pyrophosphate.

That puts a number on a limitation `solution_chemistry.py` has documented since
it was written: the mixed-buffer approximation applies to **17 experiments**,
the largest single BnOH series in the dataset.

**This corrects an earlier entry.** The 2026-08-30 note on the 0.1 M buffer
stock called exp 13's borate preparation "the only buffer recipe in the
archive". It is not: 17 more sit in the 135–151 sheets, and they are the better
ones — two named salts, weighed masses, final volumes, and a computed molarity
that checks out (30.1495 / 446.06 / 0.5 L = 135.18 mM).

### Standing

```
                                              before   after
sheets read for a declared buffer stock            0      21
experiments whose stock is confirmed twice         0      21
open filename-vs-sheet disagreements               3       0
fault-injection cases                             16      17
```

Nothing in `experiment_data.csv` changed. What changed is how much of it is
evidenced rather than inferred.

---

## 2026-08-30 — The 0.1 M buffer stock is evidenced, and [buf] now has its own verification

[buf] was the least-verified concentration in the dataset. [sub], [h2o2] and
[enz] each trace to a weighed mass through a recorded dilution chain; [buf] did
not, because the buffer was made up once and used across many runs. Thirty-three
experiments — the phosphate campaign of 2010-04-08 to 2010-05-06, which predates
the `0.1M` filename convention starting at t041 — rested on a fallback that
hardcodes both a 2 ml cuvette and a 0.1 M stock. I and [HOO⁻] scale linearly
with it.

### The experimenter's own arithmetic says 0.1 M

Exp 14's workbook carries the calculation in a scratch column beside the cuvette
table. In the `.ods` copy under `data/Mads/Variable Temperature/` it is still a
live formula:

```
P27  =(C20/H20)*0.1  = 0.08        C = Buf [ml] (1.6, 1.4, 1.2, 1.0)
P28  =(C21/H21)*0.1  = 0.07        H = Vol [ml] (2)
P29  =(C22/H22)*0.1  = 0.06
P30  =(C23/H23)*0.1  = 0.05
```

That is `(V_buf / V_total) × 0.1 mol/L` — the fallback exactly, in the
experimenter's hand. The `.xls` copies keep the values but drop the formula,
which is why it had never been seen.

Exp 14 ran 2010-04-18, in the middle of the campaign, and **21 of the 33 share
its cuvette recipe to the millilitre** (V_buf 1.6/1.4/1.2/1.0 in 2 ml); all 33
use a 2 ml cuvette. A scan of every sheet of all 33 workbooks found this
arithmetic in exp 14 alone, so it is one instance — but it is a direct one, and
it is what the extraction already assumed.

### Every filename that states a molarity agrees

Dividing a declared `[buf]` by the cuvette's volume ratio recovers the stock, and
it is identical across cuvettes in all 43 experiments that declare it — one
bottle per experiment, as it should be. Of the 21 whose filename says `0.1M`,
**18 recover exactly 100 mM.**

The stocks actually in use across the project:

| stock | experiments | how known |
|---|---|---|
| 33 mM | 3 | declared |
| 100 mM | 19 | declared, plus exp 14's own arithmetic |
| 150 mM | 17 | declared |
| 225 mM | 5 | declared |

So 0.1 M was never a universal convention — which is why assuming it needed
evidence rather than habit.

### `Buffer.xls` gives the recipe but cannot date-match

`data/Mads/Buffer.xls` is the preparation calculator: NaH₂PO₄·2H₂O (156.02) and
Na₂HPO₄·12H₂O (358.14) at pKa 7.20, Henderson–Hasselbalch, reproducing its own
gram figures to ten significant figures. It was **created 2010-04-15**, in the
middle of the campaign.

Its surviving values are 0.1 L at pH 8.00 and **0.4 mol/L** — but it was **last
saved 2010-08-30**, four months later, so those are August's numbers and April's
were overwritten. It confirms how the buffer was made, not how strong the April
stock was. (0.4 M is one of the stocks the late-April buffer titration used:
exps 32 and 34–37 step 0.1/0.2/0.3/0.4 M.)

### A new module: `data/verify_buffer.py`

Four checks, wired into `validate_dataset.py --deep`:

| check | what it catches |
|---|---|
| `bufstock` | a declared `[buf]` column implying a stock that changes between cuvettes |
| `buflabel` | a filename and its sheet disagreeing about the stock |
| `bufarith` | exp 14's arithmetic no longer readable — the campaign's only direct evidence |
| `bufvolume` | the fallback's hardcoded 2 ml applied to a sheet whose cuvettes are not 2 ml |
| `bufcompiled` | a compiled `[buf]` that is not what that stock gives at the sheet's own volumes |

All 100 experiments now recover a stock and pass `bufcompiled`. Two fault
injections in `data/test_validator.py` pin it: a wholesale ×4 rescaling (the
shape the 0.1-vs-0.4 M question would take) and a single wrong cuvette.

### It found a real defect on its first run: exp 58

Exp 58's cuvettes are **2.10–2.11 ml, not 2 ml**, so the fallback's division by 2
makes every `[buf]` about 5% high — 90.0 where the sheet's own volumes give
85.71. It is the only experiment in the archive where the hardcoded 2 ml is
wrong.

Exp 58 is excluded for running backwards, so no result depends on it. Recorded
as an accepted deviation rather than corrected, since correcting an excluded run
would change the compiled CSV for no benefit.

### And a filename/sheet disagreement awaiting a ruling

Exps **75, 76** (pyrophosphate) and **78** (carbonate) have `0.1M` in the
filename but a declared `[buf]` implying a **33 mM** stock — 24.75 mM in the
cuvette where 0.1 M would give 75 mM, a factor of three.

The compiled dataset follows the sheet, which is the standing precedence
everywhere else, so nothing is wrong with it as built. But exps 75 and 76 are
live and are the only two pyrophosphate runs of that era, and their neighbours
77, 79 and 80 — same week, same `0.1M` naming — recover 100 mM. Reported as a
warning rather than resolved.

### Two reader bugs fixed on the way

Both were in `build_dossier.py`, so they affected the dossier's recipe tables
as well as this module.

**The header merge ate a cuvette.** `column_labels` joins the header row with the
one below to catch labels split across two rows (`[sub]` over `mmol/l`), and
`recipe_rows` skipped that row if *any* cell in it was text. The 135–151 series
has a one-row header, so its first cuvette — whose name cell, `Kuv. 1 (a1,a2)`,
is text — was folded into the labels (`Buf [ml] 0.5`) and dropped. Seventeen
experiments were showing 6 of their 7 cuvettes. Both now use a majority rule:
the row is a units line only if it is mostly words.

**An unlabelled totals row passed as a cuvette.** Exp 9 ends with a row whose
`kuv` cell is blank and whose Buf and Vol cells hold 10.4 and 16 — the column
sums of the eight cuvettes above. It has a filled total-volume cell, so the
existing `Sum:` guard let it through. A cuvette must carry an identifier.

Neither bug reached `experiment_data.csv`, which `kinetics_io.py` builds with its
own reader: the dataset still reproduces from `data/data` with **zero differing
cells**, and the manifest is byte-identical apart from exp 58's new accepted
deviation.

---

## 2026-08-30 — The curve flag was measuring reaction rate, not data quality

The dossier's `curve_flags` marked 72 curves across 31 experiments. A survey of
all 421 curves in the live experiments found that most of those marks were
wrong, and that a whole class of real defect was invisible to the rule.

### It condemned 33 healthy curves

The rule was `max < 0` — "never rises above baseline". It fires on any curve
whose absorbance offset is negative, whatever the curve then does. It fired on
36 live curves and **32 of them have a genuine positive rise.** Exp 67 s2 climbs
+0.037 at 60σ; it merely started at −0.044. Baseline offset is an instrument
zero, not a measurement failure.

### The curves that really are flat are the informative ones

Ranking every cuvette by the rate its own recipe predicts, **35 of 39** flagged
curves sit in the slowest-expected half of their experiment. They are the lowest
rung of a deliberate titration — exp 136 s7 is 3.67 mM H₂O₂ against s1's 73.4 —
or the zero-enzyme control, exp 128 s5.

Those points carry the K_M information. Discarding a cuvette for being slow
would drop the low-concentration end and bias every fit toward saturation. Flat
is therefore reported as a **note**, never a defect.

### The rule could not see a curve that thrashes

It read only the endpoints, so a curve that bounces wildly and happens to land
in the right place passed clean. Measuring backtracking — absorbance travelled
against the curve's own direction, after a 5-point median filter — finds **37
curves in 20 live experiments moving more than 0.02 AU backwards, of which the
old rule flagged 2.**

```
exp 131 s2   net +0.106   backtracks 0.474 AU     4.5x more motion backwards than forwards
exp 131 s1   net +0.069   backtracks 0.391
exp 135 s2   net +0.113   backtracks 0.353
```

Exp 131 s1 reads `+0.00 +0.03 +0.09 +0.10 +0.10 +0.15 +0.04 +0.09 +0.13 +0.15
+0.07` — not a progress curve. In exps 135 and 138 samples 1–4 thrash while
samples 5–7 of the same run are smooth.

### The thrashing tracks peroxide, and it is not a cuvette-position effect

Within an experiment, the highest-[H₂O₂] cuvette backtracks more than the lowest
in **19 of 22** experiments where [H₂O₂] varies (binomial p = 0.0009).

The obvious confound is that the high-peroxide cuvettes are also the early
positions in every one of these designs. It is ruled out by the 68 experiments
where [H₂O₂] is constant across all cuvettes: there, position does nothing —
early beats late in 24, late beats early in 23, p = 1.00.

This is consistent with O₂ from peroxide decomposition crossing the light path.
**For fitting, the affected points are not biased on average but their residuals
are far from Gaussian**, so unweighted least squares will over-trust them.

### The replacement

`curve_findings()` in `data/build_dossier.py` returns `(level, message)`, where
`defect` means the measurement is wrong and `note` means something worth seeing
that may be the experiment working as designed.

| finding | level | test |
|---|---|---|
| runs backwards | defect | net < 0, beyond 5σ, and at least 0.01 AU |
| backtracks | defect | ≥ 0.02 AU and ≥ half the net |
| only N points | defect | fewer than 20 |
| flat | note | net within 0.01 AU |

Noise is the median absolute second difference (which annihilates any trend),
floored at the quantisation σ of the instrument's three-decimal output, since a
still curve has first differences of exactly zero and any spread-based estimate
would collapse to nothing.

Live experiments now carry **40 defective curves across 18 experiments** and 61
noted flat, against 54 flagged before. `data/test_curve_flags.py` pins this with
28 checks: synthetic curves whose defect is known by construction, a regression
test for each way the old rule was wrong, and the archive counts themselves.

### Two whole-run failures surfaced. Both are now excluded (ruled 2026-08-30)

**Exp 50** — all four curves descend 0.075–0.105 AU, and the descent has *no
ordering by substrate* (2.06 mM → −0.094, 4.12 → −0.075, 6.19 → −0.105, 8.25 →
−0.086). Chemistry would order them; an instrument fault would not.

```
exp 55  same buffer, pH, ladder and [enz]   +0.145 +0.281 +0.392 +0.476   monotone
exp 51  same day as 50                      +0.054 +0.072 +0.082 +0.100   monotone
exp 50                                      -0.094 -0.075 -0.105 -0.086   no ordering
```

A same-day sibling is normal and the repeat under identical conditions is
textbook — the criterion that excluded exps 72 and 82. It survived until now
because it is hand-sorted into `data/Mads/good data BnOH/`.

**Exp 85** — all seven curves crash, and the drop is *anti*-correlated with
substrate: 14.30 mM falls −0.298, 0.51 mM falls −0.604. Nearly all of it happens
in the first 15 of 60 minutes and then flattens, which is a decay rather than a
reaction. The sheet gives pH **11.84** in Na₂CO₃, a full 1.5 units above
carbonate's pK_a2, so the run is barely buffered. Exps 86–109 at pH 11 were
hand-sorted into `data/Mads/bad data pH ca. 11/`; exp 85 escaped only because it
sits in `data/Mads/carbonate buffer/`.

### Substrate inhibition was considered as an alternative and does not fit

The proposal: substrate crowding the cavity blocks H₂O₂ from reaching the
ketone, so less dioxirane forms and the rate falls at high substrate. It is the
right shape of idea for a negative substrate dependence, but it cannot produce
these curves.

Inhibition drives the rate **toward zero**; its limit is a flat line, never a
descending one. Exp 85 falls 0.30–0.63 AU. It also predicts the *low*-substrate
cuvette to be the fastest and therefore the **largest positive** — exp 85's
low-substrate cuvette is the most negative of the seven. And exp 50 shows no
substrate dependence at all (ρ = 0.00), which is the absence of the effect
rather than an instance of it.

Tested independently on the healthy data, pooling replicate series with each run
normalised by its own mean rate so no rung is favoured:

```
borate 4OMe, exps 41-49 (n=9)          BnOH 7-cuvette, exps 135-151 (n=16)
   1.53 mM  0.58 +/- 0.04                 0.22 mM  0.14 +/- 0.13
   3.06 mM  0.89 +/- 0.02                 0.86 mM  0.85 +/- 0.19
   4.59 mM  1.15 +/- 0.02                 3.03 mM  1.53 +/- 0.11
   6.11 mM  1.38 +/- 0.04                10.82 mM  1.50 +/- 0.17
```

The widest-range series (50×) saturates by 3 mM and stays flat to 10.8 mM,
p = 0.90 — ordinary Michaelis–Menten with no downturn. The tightest-signal
series is still climbing at the top of its range. **There is no substrate
inhibition to find below ~11 mM**, which is a useful negative result in itself.

Settling the competitive question properly needs a 2-D grid in [sub] × [H₂O₂].
The 7-cuvette design gives two 1-D slices meeting at a single point, so the
archive cannot answer it.

### The notebook had been sign-flipping both runs, and that manufactured the effect

`clean_experiment_dataframe` carried a block negating `[P]` for exps 50 and 85.
It has been removed.

Negation leaves the ordering of |net| against [sub] untouched, so it rescues
neither run — exp 50 stays at ρ = 0.00, exp 85 at −0.71, against +1.00 for every
healthy run. Flipped, exp 50 contributes four cuvettes rising by the *same*
amount whatever the substrate.

Flipped exp 85 is worse than useless: a rate that falls with substrate is
exactly the substrate-inhibition signature, produced here out of a decaying
baseline. Anyone fitting it would have found strong apparent inhibition that the
real ladders say is not there.

Exp 85 was already being dropped by the carbonate-buffer rule at the end of the
function, so only exp 50 reached the fits — **four sign-flipped curves with no
substrate dependence, in the dataset the catalysed constants are fitted on.**

### What changed

| | before | after |
|---|---|---|
| manifest exclusions | 9 | 11 |
| `clean_experiment_dataframe` output | 408 rows, 89 experiments | 404 rows, 88 experiments |
| `[P]` sign flips | 2 experiments | none |

`experiments_to_remove` in the notebook and `KNOWN_EXCLUSIONS` in
`data/build_manifest.py` now hold the same eleven numbers.

---

## 2026-08-30 — The 62 uncompiled files: only two were recoverable, and both are now in

Inventory of everything in `data/Mads` that had never reached the dataset.

| | n | |
|---|---|---|
| analysis workbooks, not runs | 5 | `Buffer`, `Hammett`, `Rate(T)`, `Rate(pH)`, `Rate_uncat` — two now serve as verification sources |
| duplicate copies of compiled experiments | 6 | #143 (+autosaved), t018(1), t019(1), t042 variant, t052, t053 |
| **never compiled, has raw curve data** | **2** | **exps 3 and 64** |
| never compiled, **no raw curve data** | 49 | t001, t054, t056, t063, t081, the pH-11 batch t086–t109, the solvent-isotope batch t110–t126 and t132, #133, #134 |

### The archive holds measurements for exactly 100 experiments

All 105 `.txt` files under `data/Mads` parse, and between them they cover 100
distinct experiment numbers: the 98 already compiled plus exps 3 and 64. Every
file in the last row above has a **recipe sheet and no measurement**. Those runs
are not missing from the dataset; they are missing from the archive.

**This closes the Hammett question as impossible rather than pending.** t054 and
t056 are genuine 4-bromobenzyl alcohol runs, but neither has an instrument
export, so the bromo arm cannot be built from this archive at all. It is no
longer a scope decision.

### Exp 3 added — the largest enzyme-free experiment in the dataset

`t003, BnOH, phosphate pH 6.71, 25 °C, no enzyme.` Seven cuvettes, a substrate
ladder 1.28 → 8.98 mM at fixed H₂O₂ = 82.5 mM, 227 points over 185 minutes at
dt = 49 s. All seven curves rise, none flat, none backwards; the enzyme block
reads `kuv = 0` explicitly.

Every other enzyme-free experiment in the dataset has **four** samples — 23 of
them, one with three. Exp 3 nearly doubles that, and the catalyst-free constants
`k_can` and `k3` are fitted on precisely this set before anything catalysed is
touched.

It uses M = 109.13 for benzyl alcohol, the early value that overstates the true
108.14, so its `[sub]` is 0.9% low — the same defect already recorded for exp 6.

### Exp 64 added and excluded

`t064, BnOH, boric pH 8.51, no enzyme, H₂O₂ = 244.9 mM` — the highest peroxide
in the dataset, which is why it was worth looking at.

```
Sample001  17 pts   7 min   net +0.006
Sample002  17 pts   7 min   net +0.000
Sample003  16 pts   7 min   net -0.001
Sample004  16 pts   7 min   net -0.007
```

Seven minutes, and three of four curves flat or backwards — the criterion exps
72 and 82 were excluded on. Its session was troubled throughout: the sibling
sheet `mads_t063_..._no_E_94` is named **NO_DATA_FILE**, and t063 has no export
at all.

Compiled anyway and marked excluded, so the archive-to-dataset mapping is
complete and the exclusion is recorded rather than silent.

### The rebuild was safe because the dataset is reproducible

Before touching anything, the compiled CSV was regenerated from `data/data` and
diffed against the stored one: 443 rows, 12 columns, **zero differing cells**.
Every hand-patch made today — exps 79/80's `[enz]`, exps 57/58's `[sub]` revert
— is encoded in `EXPERIMENT_CORRECTIONS` and reproduces from the sheets. So
adding experiments is a rebuild rather than a merge, and the check was repeated
afterwards: all 443 pre-existing rows are byte-identical, with 11 new ones
appended.

**Both new experiments pass every deep check with no accepted deviation.**

| | before | after |
|---|---|---|
| rows / experiments | 443 / 98 | **454 / 100** |
| traced `[enz]` to weighed catalyst | 98 | **100** |
| concentration columns from volume tables | 248 | **250** |
| cuvette values through dilution chains | 341 | **345** |

### Noted, not fixed

`data/data` holds two copies of exp 32's sheet, one with a truncated filename.
Their `Sheet1` contents are identical, and `find_and_parse_experiment_file`
takes the first match, so nothing depends on which is read.

Errors 0, warnings 9, notes 11. With `--deep`: 0 errors.

---

## 2026-08-30 — A third, independent record of [enz]: Rate(pH).xls confirms all nine points

`verify_enzyme.py` traces `[enz]` through the sheet that ran the experiment — its
weighed catalyst, its stock, its cuvette volumes. That is **one document**. If a
sheet's enzyme block were itself wrong — copied from another run, which is
exactly what happened to exps 57/58's substrate block — the whole chain would
agree with itself and still be wrong.

`Rate(pH).xls`, in `data/Mads`, is the experimenter's own pH-dependence analysis
of the 4OMe-BnOH series. It tabulates the enzyme concentration used at each pH,
written down for a different purpose, at a different time, from the analyst's
side rather than the bench's. **`data/verify_rate_workbook.py`** checks the
compiled dataset against it.

### All nine points confirm

| pH | workbook `[E]` mM | candidates | matched |
|---|---|---|---|
| 5.64 | 0.175330 | 9, 11 | **9** |
| 5.87 | 0.272695 | 10 | 10 |
| 6.50 | 0.240683 | 22 | 22 |
| 6.71 | 0.175330 | 2, 4, 5, 7 | 2, 4, 5, 7 |
| 7.00 | 0.272695 | 14, 15, 16, 17, 18, 19 | 14, 15, 17, 18 |
| 7.50 | 0.240683 | 20, 84 | 20 |
| 8.01 | 0.175330 | 8 | 8 |
| 8.50 | 0.240683 | 21, 41 | 21 |
| 8.95 | 0.272695 | 12, 42 | 12, 42 |

The workbook labels rows by pH rather than by experiment, so each row is matched
to every 4OMe enzyme experiment within 0.05 pH and passes if any carries the
tabulated value. **Nine of nine.** Where several experiments share a pH the
`[E]` value discriminates between them — at pH 5.64 it picks exp 9 (0.17533)
over exp 11 (0.272695), and at pH 7.00 it separates exps 14, 15, 17 and 18 from
exps 16 and 19.

### Two mapping notes, reported not raised

- **exp 9** — the workbook labels its row pH **5.64** where the dataset has
  **5.67**. That is the filename value, and we ruled to the sheet earlier today.
- **exp 42** — the workbook labels its row pH **8.95** where both the dataset
  *and the filename* say **8.98**.

Exp 42 is what makes exp 9 easy to read: the workbook's pH labels drift by up to
0.03 in **both** directions and disagree with the filename as readily as with the
sheet, so they are round numbers for plotting rather than a pH source. They do
not disturb the exp 9 ruling. Reported at NOTE level, since they record which
value the experimenter's own analysis used rather than any disagreement about a
measurement.

### What this catches that the sheet chain cannot

A fault-injection case pins it (14 total, all passing): setting exp 8's `[enz]`
to a wrong value now raises under `ratee`. The sheet chain would still pass such
a fault if the sheet itself were the source of the error — a copied enzyme block
would be internally consistent all the way down. Only a document written
independently of that sheet can catch it.

`[enz]` now has **three** sources: the cuvette table, the weighed-catalyst chain
in the same sheet, and this workbook for the nine pH points of the 4OMe series.

---

## 2026-08-30 — [enz] traced back to the weighed catalyst, in all 98

`[enz]` was the one concentration with no second source. `[buf]`, `[h2o2]` and
`[sub]` are checked against the cuvette volume tables and, where serially
diluted, traced through the recorded dilution chain back to weighed grams.
`[enz]` was read from a single cell of the cuvette table, and when that cell was
wrong nothing noticed — exps 79 and 80 held `0.000001` there, roughly 14,000×
too low, and the error survived until it was spotted by eye.

The sheets carry a second route. Every one prepares the catalyst in a header
block recording the molar mass, the mass weighed, the volume and the resulting
stock, then the concentration in the cuvette:

```
g/mol   1054.29          <- the catalyst
g       0.0236
l       0.004
mol/l   0.005596
mmol/l  5.596183         <- stock
kuv     0.139905 mmol/l  <- in the cuvette
```

**`data/verify_enzyme.py`** checks every link of that chain:

| link | check |
|---|---|
| `enzmw` | the block really is the catalyst — M = 1054.29 g/mol |
| `enzstock` | the block's own arithmetic: `g / (g/mol) / l` = the declared mM |
| `enzkuv` | stock × `V_enz / V_total` = the declared `kuv` |
| `enzuse` | the compiled dataset agrees with that cuvette value |

**All 98 experiments now trace `[enz]` back to a weighed mass of catalyst**, and
the first three links pass everywhere. `[enz]` has the same standing as the
other three concentrations.

### Anchoring

The first version keyed on the `kuv` label and found only 81 of 98. The later
series (exps 135–151) lays the block out differently — volume in ml rather than
l, and no `kuv` row at all, the cuvette concentration living in the table
instead. Anchoring on **the catalyst's molar mass, 1054.29**, finds all 98: that
number appears in no other block, the substrate blocks carrying 108.14, 138.17
or 187.03. Where there is no `kuv` row, the volume route supplies the cuvette
value.

### What fires, and why it is right

Only exps 32, 34, 35, 36 and 37 — the recovered buffer titrations. Their enzyme
block describes the four **planned** with-enzyme cuvettes, which were never run:
these are enzyme-free titrations, their filenames say `with_NO_E`, all five sit
in `data/Mads/'No enzyme'/`, and only cuvettes 5–8 were measured. `[enz] = 0` is
correct and the block's 0.24–0.27 mM is the plan. Declared as an accepted
deviation with that reason.

### Proof that it bites

Two fault-injection cases, run against the deep check directly (13 total, all
passing):

- **the bug it was built for** — restoring exps 79 and 80 to `[enz] = 0` now
  raises immediately, where before it survived unnoticed through every check;
- **a silently halved `[enz]`** on exp 14 — the failure a scanner cannot catch,
  because the value is neither missing nor absurd, only wrong.

Errors 0. With `--deep`: 0 errors, 10 accepted deviations, all explained.

---

## 2026-08-30 — The 0.1 M buffer stock: one recipe exists, and it computes to 100.5 mM

`[buf]` is inferred as `(V_buf/2)·100` in 56 experiments, and for 33 of them the
0.1 M stock that rule assumes was stated nowhere. Searched every sheet in the
archive — 98 compiled, 62 uncompiled, and the summary workbooks — for free text.
Six substantive notes exist in the whole collection. **One is a buffer recipe.**

### Exp 13

> *"The Buffer was prepeared by mixing 0.0699 g NaBH4 with 0.1964 g B(OH)3 in
> 0.05 l water and the pH was adjusted with a NaOH solution."*

```
B(OH)3   0.1964 g / 61.83 g/mol = 3.176 mmol
NaBH4    0.0699 g / 37.83 g/mol = 1.848 mmol
                          total = 5.024 mmol boron in 0.05 l
                                = 100.5 mM
```

**0.1 M to within half a percent.**

Two caveats belong with it:

- **The salt name is almost certainly mis-transcribed.** NaBH₄ is sodium
  borohydride — a reducing agent that would not survive contact with H₂O₂, and
  that hydrolyses in water. But every plausible substitute falls short:
  NaBO₂ → 84.8 mM, NaBO₂·4H₂O → 73.7, borax → 78.2, plain NaOH → 63.5. That the
  literal reading is the one landing on 0.1 M argues the masses are right and
  only the formula was written down carelessly.
- **It is a borate buffer**, and exp 13 is the only boric-acid experiment among
  the 33. The other 32 are phosphate, so this documents the borate stock and
  supports the phosphate one by house-standard analogy alone.

### The rest of the evidence

**The volume scheme is confirmed textually** in 29 of the 33 sheets:

> *"The buffer was added as 1\*1ml plus x\*0,2 ml, with x=3,2,1 respectivly."*

That gives 1.6, 1.4, 1.2, 1.0 ml — exactly the recorded `Buf [ml]` column. It
pins the ratios between cuvettes, which is what makes the buffer titration real.
It says nothing about the scale.

**65 filenames declare `0.1M`, from t041 onward**, across B(OH)₃, phosphate and
carbonate. The convention is documented — starting one experiment after the
assumed set ends. The 33 are t002–t040 plus t052, and not one declares a
molarity in its filename.

**`Buffer.xls` is a different, later preparation** — 0.4 M phosphate from
NaH₂PO₄·2H₂O and Na₂HPO₄·12H₂O in 0.1 l, a concentrate using the dodecahydrate
of the later era. Nothing in `Rate(pH).xls`, `Rate(T).xls` or `Hammett.xls`;
they are derived-rate summaries.

### Recorded

The manifest gains two columns, both computed once in
`build_manifest.classify_buffer`:

| `buf_provenance` | n | meaning |
|---|---|---|
| `declared` | 42 | the sheet states `[buf]` per cuvette |
| `sheet-recipe` | **1** | a weighed recipe in free text — exp 13 |
| `filename` | 18 | the stock molarity is in the filename |
| `sheet-text` | 5 | a label beside a cuvette, e.g. `1 (0.1M)` — exps 32, 34–37 |
| `assumed` | 32 | stated nowhere |

| `design` | n |
|---|---|
| `volume` | 51 — buffer volume traded against substrate |
| `stock` | 47 — fixed volumes, diluted stocks |

The recipe itself is stored in exp 13's `notes`, so the evidence travels with
the classification.

### A drift the columns exposed

`build_dossier` had its own copy of this classification, and the two had already
diverged: it counted a trailing sum row in exps 127–131 and called five
fixed-volume experiments volume titrations (56/42 against the correct 51/47).
The dossier's copy is deleted; it now reads the manifest's columns, like
everything else.

### Two side notes

- The recipe confirms **NaOH was used for pH adjustment**, which adds sodium and
  therefore ionic strength that `solution_chemistry` does not model. Previously
  a suspected limitation, now a documented one.
- **`Rate(pH).xls` carries an `[E] [mmol/l]` column per pH** for the early
  campaign (0.17533, 0.272695, 0.24068 …) — a second independent record of
  `[enz]` for exactly the experiments where the enzyme column has been trouble.

---

## 2026-08-30 — Exp 80 was an enzyme run posing as a control; last questions closed

Four rulings, one of which was a live defect in an experiment that is in use.
After them the dataset carries no unresolved per-experiment questions.

### Exps 79 and 80: enzyme runs, `[enz]` was zero

Both filenames say `with_E`; both compiled with `[enz] = 0`.

Their cuvette tables **do** carry an `[Enz] mmol/l` column — it holds
**0.000001** in every measured row. That is a broken formula, roughly 14,000×
too low. The extraction read it faithfully and rounding to three decimals turned
it into `0.0`.

The right value sits in the sheet's header block, on the `kuv` row of the enzyme
stock calculation, and it reconciles with the recorded volumes:

| | stock | Enz [ml] / Vol | header `kuv` | check |
|---|---|---|---|---|
| exp 79 | 0.559618 mM | 0.05 / 2 | 0.01399 | 0.559618 × 0.025 = 0.013990 ✓ |
| exp 80 | 0.559618 mM | 0.1 / 2 | 0.027981 | 0.559618 × 0.05 = 0.027981 ✓ |

That header row is the same one that agrees with the table exactly in every
healthy sheet — exp 2 declares `kuv` = 0.17533 and its table column reads
0.17533.

**Why this one mattered.** Exp 79 is already excluded for running backwards, so
it is consequence-free. **Exp 80 is in use**, and with `[enz] = 0` it was
sitting in the dataset indistinguishable from an enzyme-free control — which is
precisely the set the catalyst-independent constants `k_can` and `k3` are meant
to be fitted on before `k5'` and `k6` are touched. An enzyme run in that set
would have inflated the uncatalysed rate and then been subtracted from the
catalysed one twice over.

Corrected in `kinetics_io.EXPERIMENT_CORRECTIONS`; the deep check reads the
table, so it necessarily disagrees and the deviation is declared with that
reason.

**Scope, checked across all 98 sheets:** 63 declare a header `kuv`; **58 agree**
with the compiled `[enz]`; the only other disagreements are exps 32 and 34–37,
where the `kuv` belongs to the planned-but-unmeasured with-enzyme rows and is
already handled. So the broken column is confined to these two.

That the header `kuv` agrees with the table in 58 of 63 sheets makes it a usable
**independent cross-check on `[enz]`**, in the same way the volume tables are for
the other three concentrations — currently unexploited.

### Exps 9 and 38: pH ruled to the sheet

Both disagree with their filename by hundredths, and both are ruled to the
**sheet**: it carries the reading taken on the day, and the filename is typed
from it.

| | filename | sheet | ruled |
|---|---|---|---|
| exp 9 | 5.64 | 5.67 | **5.67** |
| exp 38 | 6.97 | 7.00 | **7.00** |

The dataset already held the sheet value in both cases, so **no number changes**
— what changes is that the provenance is stated rather than left as an
unadjudicated conflict.

**With these two, the dataset has no unresolved questions left.** Every
disagreement between the filenames, the hand-sorted folders, the sheets and the
extraction has been either ruled on or turned into an exclusion. Fourteen
experiments carry a recorded ruling (2, 4, 5, 7, 8, 9, 10, 38, 57, 58, 79, 80,
84, 85) and eight are excluded (57, 58, 72, 77, 78, 79, 82, 84).

What remains open is not about individual experiments: the buffer stock assumed
for 33 of them and the Debye–Hückel range. The design classification and the
uncompiled files have since been resolved; see the entries above.

### A `RULINGS` bug this exposed

Exp 9 needed a second ruling on top of its existing `abs_nm` one, and the
plain assignment `RULINGS[9] = {...}` silently discarded the first — visible
only as the note count dropping from 9 to 8. Every ruling now merges into the
experiment's entry rather than replacing it.

Errors 0, warnings 12 → **8**, notes 9 — the 8 remaining warnings are all
exclusions working as designed. With `--deep`: 0 errors.

---

## 2026-08-30 — Exps 57/58 ruled 4OMe-BnOH; exp 57 excluded, [sub] not recoverable

Ruled by the experimenter: **exps 57 and 58 are 4-methoxybenzyl alcohol runs.**
The workbook was copied from t056, a 4-bromobenzyl alcohol run, with the
substrate label in cell C4 and the filename updated and the stock block and
optics left behind.

### What the copy left behind

| | t056 (4-brom) | t057 | t058 |
|---|---|---|---|
| stock mass / volume | 0.3511 g / 0.1 L | **same** | **same** |
| molar mass used | 187.03 | **187.03** | **187.03** |
| stock molarity | 18.772 mM | **18.772 mM** | **18.772 mM** |
| wavelength / `e` | 285 nm / 1.59 | 285 nm / 1.59 | 285 nm / 1.59 |
| label in cell C4 | `4-brom-BnOH` | `4-MeO-BnOH` | `4-MeO-BnOH` |

The optics being stale is unremarkable — it is the same pattern as exps 2–10,
and the dataset's 300 nm and `e` = 7.53 are correct by the substrate convention,
so **no ruling was needed there**.

The **molar mass** is the consequential one. The sheet divides 0.3511 g by
187.03 to get 18.772 mM. For 4-methoxybenzyl alcohol, M = 138.17, so the stock
was **25.411 mM** and every cuvette is higher than recorded by

> 187.03 / 138.17 = **1.3536**

That rescaling was applied, then **withdrawn** — see below. `[sub]` in the
dataset remains the sheet's own value, and **exp 57 is now excluded** alongside
exp 58.

### Why the correction did not survive

Rescaling by 187.03/138.17 assumes the mass and volume are the real methoxy
weighing and only `M` stayed stale. Cataloguing every stock preparation in the
archive argues against exactly that assumption.

**Stock reuse is real and common** — six long-lived preparations run across the
dataset, and t055 reuses the BnOH stock first made at t050/t051 after skipping
over t052–t054, which used other substrates. So a bottle genuinely does get
carried forward and picked up later. **But whenever that happens the molar mass
travels with the mass and volume and is correct for the substrate.** From t001
to t062 there is no case where the label says one compound and `M` says another
— except t057 and t058.

**And no earlier 4OMe preparation matches.** Every methoxy stock ever recorded:

```
0.021 / 0.0223 / 0.0224 / 0.0208 / 0.0219 g in 0.01 L
0.15335 g / 0.06 L    0.1425 g / 0.05 L    0.0528 g / 0.05 L
0.0158 g / 0.06 L     1.0 g / 0.05 L       0.2112 g / 0.05 L
0.0158 g / 0.12 L
```

None is 0.3511 g / 0.1 L. The only preparation in the archive with those numbers
is the bromo stock at t054/t056. The most recent methoxy stock before t057 was
0.0158 g / 0.12 L = **0.953 mM**, a month older and 27× too dilute for these
runs, so that bottle cannot have been used either.

So a fresh methoxy stock must have been weighed, and the sheet records t056's
instead. Weighing a different compound and landing on 0.3511 g — matching t056
to a tenth of a milligram — is not credible, which means `g` and `L` are copied
too and the real stock is **unrecorded**. `[sub]` for these two is therefore
unknown, not merely wrong by a known factor.

### Outcome: exp 57 excluded

Exp 58 was already excluded for running backwards; **exp 57 is now excluded as
well**, with `[sub] is not recoverable` as its reason. `[sub]` in the dataset is
left at the sheet's own value rather than at an estimate nothing supports, and
`CONCENTRATION_RESCALINGS` — which existed only for this case — is removed
rather than left as machinery for inventing numbers. The accepted deviations it
required are removed with it, so `--deep` is clean on the strength of the data
rather than on a declaration.

The cost is one carbonate experiment. The notebook already discards every
carbonate run, so nothing downstream loses anything.

### Two side findings

- **`Hammett.xls` has its first block transposed.** It lists file 50 as p-OMe
  and 45 as p-H, while the filenames are t050 = BnOH and t045 = 4OMe — swapped.
  Its second block (42 = p-OMe, 51 = p-H) matches the filenames exactly. The
  file is an uncompiled analysis workbook, so nothing in the dataset depends on
  it, but it should not be trusted as a substrate source.
- **The notebook drops every carbonate experiment** (`buffer != "Carbonate"`),
  which is exps 57, 58, 77, 78, 79, 80, 85 — 31 rows. So this error never
  reached the notebook's analysis. It was live only in
  `data/experiment_data.csv`, which is the artifact everything is being moved
  onto, so it would have mattered from the next step onward.

Errors 0, warnings 11, notes 9. With `--deep`: 0 errors, all five deviations
accepted and explained.

---

## 2026-08-30 — The stock-solution block names the substrate; exps 57/58 are a third compound

Pointed at the `Stamopløsning` (stock solution) block, which every sheet uses to
weigh out its substrate. It carries a **molar mass**, and the sheet's own
concentrations are computed from it — so it is not a label but the arithmetic
the experiment actually ran on. That makes it the strongest substrate evidence
in the files, and it is machine-readable.

| substrate | M (g/mol) |
|---|---|
| benzyl alcohol | 108.14 |
| 4-methoxybenzyl alcohol | 138.17 |
| 4-bromobenzyl alcohol | 187.03 |

A molar mass was found in **97 of 98** sheets.

### Exps 84 and 85: closed, 4OMe-BnOH

Their filenames say BnOH. Everything inside the sheets says otherwise:

- the stock block is labelled **`Stamopløsning 4 / 4-MeO-BnOH [g]`**;
- it computes from **M = 138.17**, so every concentration already assumes 4OMe;
- the method string reads **`1h_kcat(4OMeBnOH)_7cuv`**;
- the optics are 300 nm and `e` = 7.53, the 4OMe convention.

Four independent sources against the filename alone. Both sheets also carry
**83** as their experiment number — the workbook was copied from t083, a genuine
BnOH run, and the header was never updated, which is almost certainly where the
filename's "BnOH" came from as well. The dataset already carries 4OMe-BnOH, so
**no number changes**; the question is simply closed.

### Exps 57 and 58: this morning's ruling was wrong

The same scan found the declared molar mass disagreeing with the manifest in
exactly two experiments — and it is **187.03**, 4-bromobenzyl alcohol.

Laying the sheets beside the two `4-brom-BnOH` runs that sit in `data/Mads` and
were never compiled:

| sheet | label in the sheet | nm | e | M | g | L | stock |
|---|---|---|---|---|---|---|---|
| t054 (uncompiled) | `4-brom-BnOH` | 285 | 1.23 | 187.03 | 0.3511 | 0.1 | 18.772 mM |
| t056 (uncompiled) | `4-brom-BnOH` | 285 | **1.59** | 187.03 | 0.3511 | 0.1 | 18.772 mM |
| **t057** (in use) | `4-MeO-BnOH` | 285 | **1.59** | **187.03** | **0.3511** | **0.1** | **18.772 mM** |
| **t058** (excluded) | `4-MeO-BnOH` | 285 | **1.59** | **187.03** | **0.3511** | **0.1** | **18.772 mM** |

Exps 57 and 58 use the **same stock solution as the confirmed 4-bromo runs** —
same compound, same weighed mass, same volume, same resulting molarity. The
`4-MeO-BnOH` text in cell C4, and the filenames' `4-MeOH-BnOH`, are the stale
parts; t058 likewise carries `57` as its experiment number.

**So `e` = 1.59 is not a stale workbook artifact. It is the 4-bromo
convention**, and this morning's ruling — that 57/58 were the same
copied-template case as exps 2–10 — was wrong. The pattern it followed is real
but belongs to t054: that sheet pairs M = 187.03 with the BnOH template's
`e` = 1.23, and the value is corrected to 1.59 by t056, exactly mirroring the
285 → 300 correction in the early 4OMe series.

### What this costs

**Exp 57 is in use**, recorded as 4OMe-BnOH at 300 nm with `e` = 7.53. All three
are wrong; on the substrate convention it should be 4-bromobenzyl alcohol at
285 nm with `e` = 1.59, so its `[P]` is understated by **7.53 / 1.59 = 4.74×**.

The `[sub]` concentrations are unaffected — they come from the sheet's own
volume table, computed with M = 187.03, which is correct for the bromo compound.

Two further consequences:

- `SUBSTRATE_PROPERTIES` has **no 4Br-BnOH entry**, so the substrate cannot
  currently be represented at all;
- t054 and t056, the two genuine 4-bromo runs, are among the 62 files in
  `data/Mads` that were never compiled. With `Hammett.xls` sitting in the same
  folder, the dataset evidently contains the beginnings of a **Hammett series**
  — BnOH, 4-MeO, 4-Br — of which only two mislabelled members are present.

Held for a ruling rather than patched: whether to relabel 57/58 and add the
substrate, or exclude them and leave the bromo series out of scope, is a
decision about what the thesis covers.

---

## 2026-08-30 — Optics settled for the whole dataset: wavelength follows the substrate, ε is a convention

Two general rules, supplied by the experimenter, close every remaining optics
question at once:

> **The monitoring wavelength is a property of the substrate, not of the run.**
> BnOH is read at 285 nm and 4OMe-BnOH at 300 nm, throughout.
>
> **`e` is a conversion factor applied uniformly at analysis time, not a
> per-experiment measurement.** A sheet's `e` cell is that workbook's own
> working note and is authoritative for nothing.

### The rules reproduce the dataset exactly

| | |
|---|---|
| dataset, by construction | 4OMe-BnOH → (300 nm, `e` = 7.53), 220 rows; BnOH → (285 nm, `e` = 1.23), 223 rows |
| sheets that deviate | **9 of 98**, every one a 4OMe workbook carrying the BnOH template's 285 nm |
| exps 2, 4, 5, 7, 8, 9, 10 | 285 nm with `e` = 7.53 — ruled earlier today |
| exps 57, 58 | 285 nm with `e` = 1.59 — the only two where the `e` cell was changed as well |

So exps 57/58 are the same stale-template artifact as the other seven. Both are
now ruled to the substrate convention.

### This was a wrong contract, not just a wrong value

`abs_nm` and `e_declared` had been treated as **ground truth read from the
sheet**, which the dataset must match — which is why exp 57 was being reported
as "every concentration scaled by 4.74×". If the wavelength is fixed by the
substrate and `e` is an analysis convention, that test is simply the wrong one.

What actually matters is **uniformity**: every experiment on a given substrate
must use the same `(abs, e)`. That was true by construction and asserted
nowhere. `validate_dataset.py` now checks it at dataset level, and the
per-experiment check compares against `SUBSTRATE_PROPERTIES` rather than
against the sheet. Two fault-injection cases pin it (11 total).

### The evidence is not erased

A ruling used to overwrite the extracted value, leaving the sheet's own reading
only in prose. The manifest now carries **`abs_nm_sheet`** and **`e_sheet`**
alongside `abs_nm` and `e_declared`: the first pair is the observation, the
second the adjudicated value, and a ruling never touches the observation. The
validator reads the `*_sheet` columns and reports all nine deviations at a new
**NOTE** level — visible, not defects, never gating a rebuild.

### What still propagates, and is now a standing assumption rather than a question

Within a substrate, an error in `e` is harmless: it is a single global scale
factor absorbed into the fitted rate constants, identical for every run.

**Between substrates it is not.** BnOH uses 1.23 and 4OMe-BnOH uses 7.53, a
factor of **6.1**. Any comparison of rate constants across the two substrates —
a substrate effect, a Hammett-type argument — inherits the ratio 7.53/1.23
directly. Neither number is a defect and neither needs changing, but their
*relative* accuracy is a live assumption for any cross-substrate conclusion.

These are not arbitrary numbers: 1.23 mM⁻¹cm⁻¹ = 1230 M⁻¹cm⁻¹ is consistent
with benzaldehyde's own ε at 285 nm (lit. ≈ 1400 at 278–279 nm in water), so
the convention was chosen to be physical. No source has been found for
4-methoxybenzaldehyde at 300 nm, which makes the ratio the natural thing for
`COMPUTATIONAL.md` task C1 to pin down alongside its main question.

And the deeper issue is untouched: if the signal is `[A] + r·[BA]` rather than
`[A]` alone, no single `e` is correct for either substrate, because the
conversion assumes the aldehyde is the sole absorber. That stays with C1.

Errors 0, warnings 21 → **13**, notes 9.

---

## 2026-08-30 — Wavelength question closed for exps 2–10; exps 57/58 stay open

The seven earliest 4OMe-BnOH runs declared **285 nm** while the dataset used
**300 nm**. Ruled: **300 nm. The sheet label is stale.**

### The evidence

Reading the declared optics in run order — including t001 and t003, which are
in `data/Mads` and were never compiled — makes the pattern plain:

| run | substrate | declares |
|---|---|---|
| t001 | BnOH | 285 nm, e = 1.23 |
| **t002** | **4OMe** | **285 nm, e = 7.53** |
| t003 | BnOH | 285 nm, e = 1.23 |
| **t004, t005** | **4OMe** | **285 nm, e = 7.53** |
| t006 | BnOH | 285 nm, e = 1.23 |
| **t007–t010** | **4OMe** | **285 nm, e = 7.53** |
| t011 … t022 | 4OMe | 300 nm, e = 7.53 |

The early series alternates BnOH and 4OMe runs. The 4OMe workbooks were copied
from the BnOH one with the substrate and `e` changed and the `abs [nm]` cell
left at 285. At t011 the cell is corrected, and it reads 300 for every
subsequent 4OMe run.

The decisive detail is that **`e` does not change across the t010 → t011
boundary.** ε is wavelength-dependent, so a genuine retune from 285 to 300
would have forced `e` to change with it — and this operator does exactly that
when they really change wavelength, which is what exps 57/58 show by pairing
285 nm with `e = 1.59`.

### What could not be used

Stated so the ruling's basis is not overstated later:

- the `.txt` instrument export records batch, instrument, sample names, data
  mode, smoothing, dates and times — **but not the wavelength**;
- conversion in these runs stays under 1% (exp 2 sample 4 ends at A = 0.093,
  i.e. 0.012 mM of 6.08 mM substrate), so no internal consistency test — a
  conversion ceiling, a rate comparison against the 300 nm runs — has any
  power. The early curves sit near background and scatter hugely.

The argument is from the workbook's own edit history, not from a measurement.

### Consequence: none, either way

The sheet and the dataset both carry `e = 7.53` for these seven, so **no
concentration changes under either reading.** This was a question about what
the signal physically is, not a numerical error — the opposite of exp 57, where
the sheet says 1.59, the dataset uses 7.53, and every concentration is off by
4.74×.

### How the ruling is recorded

`build_manifest.py` gains a `RULINGS` table: a decided value that overrides the
extracted one, marks the field's provenance as `ruling` so it can never be
mistaken for something read off a file, and writes the sheet's own value plus
the reasoning into `notes` so the evidence is not erased. Rebuilding the
manifest changed exactly 7 rows and only the columns `abs_nm`, `provenance`,
`open_questions` and `notes`.

### A second question surfaced by enforcing it

`abs_nm` had been declared in the manifest since the wavelength work but
**compared against nothing** — `validate_dataset.py` read it only to
interpolate into an error message. Adding the check (parallel to the existing
`e` check, and pinned by a tenth fault-injection case) immediately raised two
errors: **exps 57 and 58 disagree on wavelength as well as on `e`.**

That half of their conflict had been invisible. The sheet states the pair
(285 nm, `e` = 1.59); the dataset states the pair (300 nm, `e` = 7.53). Unlike
exps 2–10, the `e` differs too, so a stale template cell cannot explain it and
it needs a ruling of its own. Both fields are now listed in those experiments'
open questions, so the validator warns rather than errors. **Exp 57 is in use.**

Warnings: 24 → 21. Errors: 0.

---

## 2026-08-30 — Ionic strength and [HOO⁻] moved out of the notebook

`I` and `[HOO-]` are the only quantities in this dataset that are **computed**
rather than measured or read off a sheet, and until now both lived in notebook
cells 37 and 38. Every one of the four bugs previously recorded against that
block survived because a notebook cell cannot be imported, diffed against a
second implementation, or tested.

They now live in **`data/solution_chemistry.py`**, with
**`data/test_solution_chemistry.py`** pinning them.

### The physics is unchanged

The module reproduces the corrected notebook implementation **exactly** — max
absolute difference 2.3e-13 mM on `I` and 0 on `[HOO-]` across all 443 rows.
Nothing that has already been looked at changes.

### What the tests pin

71 checks, all passing. The load-bearing ones are the two worked by hand, which
give the module a validation gate rather than mere self-consistency:

| check | expected | why it is checkable by hand |
|---|---|---|
| 100 mM phosphate, pH 7.00 | **I = 177.37 mM** | pKa₂ = 7.20 ⇒ fractions 0.61313 / 0.38684; Na⁺ = 1.38682 per phosphate; I = ½·100·(0.61313 + 4·0.38684 + 1.38682) |
| 100 mM boric, pH = pKa | **I = 50 mM** | half dissociated ⇒ ½·100·(0.5·1 + 0.5·1) |
| pKa(H₂O₂) at I = 75 mM | **11.494** | 11.75 − 0.509·2·√0.075/(1 + 0.328·√0.075) |

Plus one regression per historical bug, each written to fail if the old
behaviour returns:

1. an unrecognised buffer name now **raises** instead of falling through to a
   single-species default — that default is what turned the `'Boric Acid'` /
   `'Boric'` key mismatch into a silent `I = 0` for 26 experiments;
2. `I` must exceed the anion-only value by more than 1.5×, catching a dropped
   counter-ion term;
3. the Debye–Hückel shift must be under a fifth of what feeding mM directly
   would give, catching the unit error;
4. the shift must equal exactly twice the `z² = 1` value, pinning Δ(z²) = 2.

### New finding: the Debye–Hückel equation is being extrapolated

Stamping the column across the whole dataset made the range visible for the
first time:

| | rows | experiments | buffers |
|---|---|---|---|
| I > 100 mM (usual validity limit) | **308 of 443 (70%)** | 71 | Carbonate, Phosphate, Pyrophosphate |
| I > 300 mM | 133 | 24 | Phosphate, Pyrophosphate |
| I > 500 mM (the usual stretch limit) | 63 | 11 | Pyrophosphate |
| I > 1000 mM | 21 | 5 | Pyrophosphate |

The maximum is **1069 mM** (exps 128, 131 — pyrophosphate, 194 mM, pH 7.29).
Pyrophosphate drives this by construction: a tetraprotic acid contributes z²
up to 16, so its `I` starts at 158 mM and never drops below it.

The consequence is bounded — the correction only ever moves pKa(H₂O₂) from
11.75 to 10.96 across the entire dataset — but it is a **systematic** error in
`[HOO-]`, not a random one, and it is largest exactly where `[HOO-]` is
largest. A Davies or Pitzer treatment would be the principled fix. Recorded as
an open item rather than patched, and exposed as
`DEBYE_HUCKEL_RELIABLE_mM` / `out_of_range_fraction()` so any result depending
on it can state its own exposure.

### Other limitations, documented rather than modelled

The module docstring records four more, none of which the dataset carries
enough information to fix:

- **mixed buffers** — several sheets prepare one "buffer" from two salts
  (exp 146 mixes Na₄P₂O₇·10H₂O with NaH₂PO₄·2H₂O), but the dataset carries a
  single name and a single `[buf]`; for that pyrophosphate/phosphate mixture
  near pH 8.7 the error in `I` is roughly 10–20%;
- **titrant** — HCl/NaOH used to reach the target pH is recorded nowhere;
- **thermodynamic pKa values** used at finite ionic strength for the buffers
  while H₂O₂'s pKa *is* activity-corrected — deliberate, so the module stays
  numerically identical to results already recorded;
- **temperature** — 25 °C constants applied across 15–40 °C.

### Notebook

Cells 37 and 38 no longer define the physics; they import the module and keep
thin wrappers so downstream cells still run. Verified: the notebook contains
**zero** copies of `buffer_pKa`, `species_charges` or `pKa_H2O2`, and both
cells execute and agree with `add_solution_columns` to floating-point
equality. This removes the last duplicated calculation between the notebook and
the pipeline.

---

## 2026-08-30 — Wavelength and extinction coefficient: the sheets declare them, the pipeline ignores them

Chasing `MECHANISM.md`'s open observable question turned up the answer in the
sheets themselves, plus a live defect.

### What the sheets declare

Every sheet states its monitoring wavelength and extinction coefficient in its
header (`abs [nm]`, `e [U/mM]`). Across the 98 experiments:

| substrate | nm | e (mM⁻¹cm⁻¹) | experiments |
|---|---|---|---|
| BnOH | 285 | 1.23 | 43 |
| 4OMe-BnOH | 300 | 7.53 | 46 |
| 4OMe-BnOH | 285 | 7.53 | 7 — exps 2, 4, 5, 7, 8, 9, 10 |
| 4OMe-BnOH | 285 | **1.59** | 2 — exps 57, 58 |

`kinetics_io` ignores all of it, hardcoding `e` per substrate in
`SUBSTRATE_PROPERTIES` (BnOH 1.23 / 285; 4OMe-BnOH 7.53 / 300).

### Live defect: exps 57 and 58

Their sheets declare `e = 1.59` at 285 nm; the compiled dataset uses **7.53**.
Every concentration in those experiments is therefore scaled by **4.74×**.
Exp 58 is already excluded, so it is consequence-free — **exp 57 is in use.**
Recorded as an open question rather than corrected, pending a ruling.

> **Reopened 2026-08-30, same day.** The ruling below and its successor were
> both wrong: these two are **4-bromobenzyl alcohol** runs, not 4OMe. Their
> sheets compute from M = 187.03 and use the same weighed stock as the confirmed
> 4-brom-BnOH runs t054/t056, so `e` = 1.59 is the 4-bromo convention, not a
> stale cell. Exp 57 is in use and its `[P]` is understated 4.74×. See the entry
> at the top of this log.

### Second question: the seven earliest 4OMe runs

> **Closed 2026-08-30 — ruled 300 nm, the sheet label is stale.** See the entry
> at the top of this log for the run-order evidence. The paragraph below records
> the question as it stood.

Exps 2, 4, 5, 7, 8, 9, 10 declare **285 nm** but carry `e = 7.53`, which is the
value every **300 nm** sheet uses, while exps 57/58 pair 285 nm with 1.59.
Either the wavelength cell is a stale template value or `e` was never updated
when the wavelength changed. An empirical test — comparing raw absorbance slope
per unit substrate and enzyme between the two groups, using the near-matched
pair exp 9 (pH 5.67, "285") against exp 11 (pH 5.64, 300) — was **inconclusive**:
the nominal replicates at pH 6.71 (exps 2, 4, 5, 7) scatter from −2.8e-6 to
+3.7e-6, so the comparison has no power. Recorded as an open question.

### The observable question is now partly answered

Literature: benzaldehyde in water has **ε ≈ 1400 M⁻¹cm⁻¹ at 278–279 nm**, the
weak n→π* band of the aldehyde carbonyl (the strong π→π* sits at 248 nm with
ε ≈ 12,000–14,000). The sheets' `e = 1.23 mM⁻¹cm⁻¹` = **1230 M⁻¹cm⁻¹ at 285 nm**
sits exactly on the falling edge of that band.

Two consequences:

1. **`e` is benzaldehyde's own ε, not a differential coefficient.** The
   abs → `[P]` conversion therefore assumes benzaldehyde is the *sole* absorber,
   which is what `MECHANISM.md` suspected. Confirmed, not merely inferred.
2. The earlier inference that 285 nm is the weak n→π* band is confirmed
   directly from the recorded wavelength, not deduced from the magnitude of ε.

**Still missing: ε of benzoate at 285 nm.** Benzoic acid's strong band is at
224–230 nm; the weak "C band" sits near 273 nm (ε ≈ 2000 M⁻¹cm⁻¹ in ethanol),
and at pH 5.5–11.8 the species present is benzoate, whose weak band is near
268 nm. At 285 nm we are on its tail. No reliable aqueous value at 285 nm was
found in the open literature — the two most promising sources (RSC
*Absorption Spectra of Benzoic Acid in Water at Different pH*, and the 1958
*Light Absorption Studies* tabulation) are both paywalled (HTTP 403).

Rough bracket from band shape: ε(benzoate, 285 nm) plausibly 100–400 M⁻¹cm⁻¹,
giving `r = ε_BA/ε_A` ≈ **0.08–0.33** — appreciably lower than the ~0.5 guessed
in `MECHANISM.md`. That does **not** settle the observable question, because the
signal contribution is `r × [BA]/[A]`, and within the mechanism benzaldehyde is
an intermediate whose pool stays small while benzoic acid accumulates. A small
`r` with a large `[BA]/[A]` still lets the acid dominate late in a run.

The decisive measurement is cheap and local: **run a UV spectrum of benzoic
acid at the working pH and read ε at 285 nm.** One cuvette settles what the
literature would not give up.

### Changes made

- `data/manifest.csv` gains `abs_nm` and `e_declared`, read from each sheet, so
  the optics become declared ground truth like every other field. Also gains
  `accepted_deviations` / `accepted_deviation_reason` seeded in the bootstrap so
  they survive a rebuild.
- `validate_dataset.py` gains an `optics` check comparing the compiled `e`
  against the sheet's declared value and reporting the exact scaling error.
- Nine new open questions recorded (exps 57, 58 on `e`; exps 2, 4, 5, 7, 8, 9,
  10 on wavelength). No data changed.

    python data/validate_dataset.py          ->  0 errors, 24 warnings
    python data/validate_dataset.py --deep   ->  0 errors, 27 warnings
    python data/test_validator.py            ->  9 passed, 0 failed

---

## 2026-08-30 — The 59 "unverifiable" concentrations are verifiable after all

The previous entry concluded that 59 concentration columns could not be checked
from the files, because they were made by serially diluting the stock rather
than by varying the volume, and the volume table does not record the stock.
**That conclusion was wrong.** The sheets do record it — in a dilution-series
table (`Fortyndingsrække` / `opløsning`) that the earlier check never looked at.
`data/verify_dilutions.py` reads it, and all 59 are now covered.

### The chain is traceable back to a weighed mass

Both sheet generations record, per dilution, the stock volume taken, the final
volume made up to, and the resulting concentration — and above it the master
stock traced to a mass and a molar mass. Every link reproduces exactly:

**Exp 65** (earlier layout, mol/L):

    0.1581 g / 108.14 g/mol / 0.01 L        = 0.14619937 M   sheet: 0.14619937118550025
    0.146199 x 0.001/0.001, /0.002, /0.004, 0.0005/0.01      all four match to 1e-9
    0.146199 M x 0.1/2 ml                   = 7.3100 mM      compiled [sub]: 7.31

**Exp 135** (later layout, mmol/L, with the dilutions labelled a1/b1/c1/d1 and
the cuvette labels indexing into them — `Kuv. 1 (a1,a2)`):

    0.148 g / 108.14 / 50 ml                = 27.3719253 mM  sheet: 27.371925282041794
    a1 = 27.3719 x 20/25 = 21.89754023      sheet: 21.89754023   (b1, c1, d1 likewise)
    Kuv.1 [sub]  = 21.8975 x 0.4/1          = 8.7590 mM      compiled: 8.759
    Kuv.1 [H2O2] = 3263.2949 x 0.05/1       = 163.1647 mM    compiled: 163.16

This is a **stronger** check than the volume route, which only confirms internal
proportionality: this one starts from grams on a balance.

### Coverage

The two routes partition the dataset exactly. The 56 experiments with no
dilution table (2–62) are the earlier volume-titration design already confirmed
by the volume route; the 42 with one are the stock-dilution design.

| | |
|---|---|
| dilution tables internally consistent | **42 / 42** |
| `[sub]` cuvettes traced to a recorded dilution | **202**, across 37 experiments |
| `[h2o2]` cuvettes traced | **139**, across 22 experiments |
| remaining findings | **2**, both exp 128 s5 |

37 + 22 = 59 — **exactly the gap the previous entry left open, now closed.** The
two residuals are the known reference-row case (`ref 5` has its own cuvette
volume as well as its own concentrations), recorded in the manifest's
`accepted_deviations`.

Two scoping errors were made and fixed along the way, both in the check rather
than the data: the chain test initially ran on species the volume route had
already confirmed (in several experiments H2O2 was taken straight from the 30%
stock, so it appears in no dilution table while being perfectly verified), which
produced 102 spurious findings; and the parser initially recognised only
tabular dilution series, missing the species-labelled standalone declarations
(`[H2O2]` / `mmol/l` / value) that exps 127–131 use for both peroxide stocks —
a further 20.

### A quantitative by-product

Exp 135's buffer traces end to end: Na4P2O7 30.1495 g → 135.181 mM and
Na2HPO4 39.0052 g → 500.003 mM, mixed 74 + 20 → 100 ml, then 9 → 12 ml, giving
150.026 mM, and 0.5/1 ml into the cuvette = **75.013 mM**, matching the compiled
`[buf]` exactly. So the two-component buffer of exps 135–151 is
**37.51 mM pyrophosphate + 37.51 mM phosphate**, not 75 mM of one species —
which is what the ionic-strength calculation needs, and what it currently
assumes wrongly.

Also visible in these sheets and worth noting for the mechanism work: the
monitoring wavelength and extinction coefficient are recorded directly
(`abs 285 nm`, `e 1.23 U/mM` for BnOH), and exps 135–151 log **pH before and
after each run** (exp 135: 7.63 → 7.77), which is a direct measurement of the
pH drift previously only estimated.

### Standing state

    python data/validate_dataset.py          ->  0 errors, 13 warnings
    python data/validate_dataset.py --deep   ->  0 errors, 16 warnings
    python data/test_validator.py            ->  9 passed, 0 failed

Every concentration in the compiled dataset is now independently corroborated —
by volume proportionality (248 columns) or by the recorded dilution chain (341
cuvette values) — except exp 128 sample 5, which is a known artefact.

---

## 2026-08-30 — Concentrations re-derived from the volume tables

The manifest closed the metadata gap but left the concentration columns
(`[enz]`, `[buf]`, `[h2o2]`, `[sub]`) with nothing independent checking them.
`data/recompute_concentrations.py` closes it by taking a different path through
the same sheets: instead of reading the pre-computed mmol/l columns that
`kinetics_io` reads, it reads the per-cuvette component **volumes** and the
total volume, back-calculates the implied stock, and re-derives each
concentration. Available as `validate_dataset.py --deep` (~2.5 s).

The volume table was locatable in **all 98 experiments** (76 sheets carry all
six volume headers, 22 lack only `H2O [ml]`, which the calculation does not
need).

### Result 1 — block selection was right in 97 of 98

This was the open worry. `find_numeric_values_below_header` takes the first
`sample_num` numeric rows below the header, and **97 of 98 sheets plan more
cuvettes than were measured** (typically 8 planned / 4 measured; 16/7 for
exps 135–151, 14/7 for exps 84–85). Selecting the measured block *explicitly* —
dropping rows labelled `ref` and honouring the manifest's `has_enzyme` — and
comparing against the compiled data produces **exactly one disagreement**:

    exp 128  s5 [sub]: sheet 8 vs compiled 9.47

which is the already-documented Round 4 finding (sample 5 is the reference row
`ref 5`, not a fifth titration condition). Recorded in the manifest's new
`accepted_deviations` column so the deep check stays green rather than
permanently red.

So the first-N-rows assumption, although structurally load-bearing in 97 of 98
experiments, produced a wrong answer only in the five cases already fixed
(32/34/35/36/37) plus this one. **The extraction was not silently damaged
elsewhere.** That is the main thing this check was built to establish.

Caveat on its strength: the block selector consults the manifest's
`has_enzyme`, which for the 23 experiments whose filename declares no enzyme
status was itself seeded from the extraction. For those the check is partly
circular and confirms less than it appears to.

### Result 2 — 248 concentration columns independently confirmed, 59 unverifiable

Counting species × experiment, the implied stock `c · V_total / V_component` is
**constant across the measured cuvettes in 248 cases** — the sheet's own numbers
reproduce exactly from the volumes, so those concentrations are independently
confirmed.

The other **59 are unverifiable by this route, not wrong**. Two dilution designs
exist in the dataset and only one is checkable:

- **Volume design** (e.g. exp 2): substrate volume steps 0.2 → 0.8 ml at a fixed
  15.2 mM stock. The implied stock is identical on every row. Verifiable.
- **Stock design** (e.g. exp 65): volume fixed at 0.1 ml while the *stock* is
  serially diluted (7.31 → 3.66 → 1.83 → 0.37 mM). The stock appears nowhere in
  the volume table, so the computed concentration is the only record of it.

Unverifiable by species: `[sub]` in 37 experiments (65–83 and the 135–151
block), `[h2o2]` in 22 (127–131, 135–151). **`[enz]` and `[buf]` are confirmed
wherever the sheet records them** — those were never serially diluted.

An earlier version of this check reported 64 "inconsistencies" that were all
this design difference, plus noise from planned-but-unmeasured rows belonging to
other designs recorded on the same sheet (exps 135–151 carry D2O variants with
a `[D2O]%` column). Both were errors in the check, not the data; restricting to
measured rows and distinguishing the two designs removes all 64.

### Result 3 — 56 experiments' `[buf]` still rests on a hardcoded assumption

`kinetics_io` reads `[buf]` from the sheet's own `[buf]`/`[buffer]` column where
one exists (**42 experiments**), and otherwise falls back to `(V_buf / 2) * 100`
— which hardcodes a **0.1 M buffer stock** (**56 experiments**: 2–62).

No contradiction was found: all 38 filenames that declare a stock molarity say
0.1 M. But that leaves the assumption unverified for the rest, and exps 32–37
are proof that non-0.1 M stocks were used in this series — they are simply
recorded as text labels (`1 (0.1M)`, `2 (0.2M)`) rather than in a column. The
five known cases are patched via `EXPERIMENT_CORRECTIONS`; whether any of the
remaining 51 used a non-0.1 M stock is **open and cannot be settled from the
files alone**.

### Standing state

    python data/validate_dataset.py          ->  0 errors, 13 warnings
    python data/validate_dataset.py --deep   ->  0 errors, 14 warnings
    python data/test_validator.py            ->  9 passed, 0 failed

### What this does and does not license

Confirmed: metadata, block selection, and 248 of 307 concentration columns.
Not established: the 59 serially-diluted concentrations, the 0.1 M stock for 51
experiments, and the six manifest open questions. No extraction code was
changed — the recomputation found nothing to fix, which is itself the result.

---

## 2026-08-30 — Validation layer: the dataset is now checked against declared ground truth

Every defect found in this log so far was found **by accident**, during
analysis aimed at something else. The extraction scans heterogeneous
spreadsheets for anything pH-shaped, buffer-name-shaped or table-shaped; when
it guesses wrong it returns a plausible number rather than raising. Nothing
ever contradicted it. That failure mode, not the parsing quality, is what this
adds a fix for.

The architecture is inverted: metadata is now **declared** in a checked-in
manifest, and the extraction is **validated** against it.

### `data/manifest.csv` — declared ground truth (98 rows, one per experiment)

Columns: `substrate, buffer, pH, T, has_enzyme, n_samples, status,
exclude_reason, provenance, open_questions, xls_file, notes`.

Bootstrapped by `data/build_manifest.py` from three sources independent of the
sheets' interiors:

1. **Filenames** — confirmed as ground truth by the experimentalist. Recover
   75–92% of fields (pH 92/98, T 81/98, substrate 92/98, buffer 81/98,
   enzyme-presence 75/98).
2. **`data/Mads/` hand-sorted folders** — `No enzyme` states enzyme status
   outright; `bad data` and `bad data pH ca. 11` record curation judgements.
3. **The extraction** — used only where the two above are silent.

Exclusions previously living as a bare list literal in the notebook
(`experiments_to_remove = [84, 58, 77, 78, 79, 72, 82]`) are now declared rows
with reasons attached.

**Where a declared source and the extraction disagree, the manifest keeps the
extracted value** and records the disagreement in `open_questions`. Adopting
the manifest therefore changed no data; adjudicating means editing the value
and clearing the note. Six such questions stand — see below.

### `data/validate_dataset.py` — the validator

Checks coverage (every compiled experiment declared and vice versa), metadata
agreement, structure (sample counts; pH/T/substrate/buffer constant within an
experiment) and invariants (extinction coefficient matches the substrate;
`[enz]` is zero exactly when the manifest says no enzyme; concentrations
finite, non-negative, non-null).

Findings are graded: **ERROR** for a disagreement nothing has accounted for,
**WARN** for one already recorded as an open question or an expected exclusion.
Exit status is non-zero on any ERROR, so it can gate a rebuild.

Current state: **0 errors, 13 warnings** — 7 expected (excluded experiments
retained by the raw compile by design) and 6 unresolved questions.

### `data/test_validator.py` — fault injection

A validator that has never failed is indistinguishable from one that cannot
fail. Nine corruptions are injected into a copy of the dataset and each must be
caught: wrong pH, wrong buffer, substrate swapped without updating `e`, enzyme
concentration lost, negative concentration, null concentration, a dropped
sample, pH varying within an experiment, an experiment appearing from nowhere.
**9/9 caught.**

Retrospective check: restoring the pre-fix `[enz]` values for experiments
32/34/35/36/37 raises **10 errors immediately**. The bug that took days to
surface would have been caught on the first run.

### Six open questions awaiting adjudication

| exp | field | filename says | extraction says | in use? |
|---|---|---|---|---|
| 9 | pH | 5.64 | 5.67 | yes |
| 38 | pH | 6.97 | 7.00 | yes |
| 79 | has_enzyme | True | False | no (excluded) |
| **80** | **has_enzyme** | **True** | **False** | **yes** |
| 84 | substrate | BnOH | 4OMe-BnOH | no (excluded) |
| **85** | **substrate** | **BnOH** | **4OMe-BnOH** | **yes** |

The two in bold are consequential:

- **exp 85** — if the filename is right, `e` should be 1.23 rather than 7.53
  and **every concentration in that experiment is off by 6.1×**. Seven samples,
  currently in use, and one of the two sign-flip-corrected runs.
- **exp 80** — if the filename is right, a run that had catalyst is currently
  being treated as an enzyme-free control, which would contaminate exactly the
  `E0 = 0` block the mechanism fit depends on.

The pH pairs (9, 38) are minor but should be settled for consistency. Note that
79 and 84 show the same two failure modes as 80 and 85 respectively, in
experiments already excluded — consistent with the extraction, not the
filenames, being the unreliable party in both cases.

### Not done

The root cause behind the worst extraction bugs — `find_numeric_values_below_header`
taking the **first** N rows of the concentration table rather than the rows
actually measured — is still unfixed. It caused both the exp 32–37 and exp 128
defects. Contained fix, not attempted here.

---

## 2026-08-30 — Buffer-titration recovery: five experiments had two wrong columns

Triggered by a claim that a buffer-concentration titration at constant `[sub]`
existed in the dataset. An exhaustive search of `data/experiment_data.csv` found
none — including across experiments (matched on substrate, buffer, T, pH within
0.15, `[sub]`/`[h2o2]` within 10%, `[enz]` within 20%, only one pair anywhere
differed in `[buf]` by more than 20%, and only 1.37-fold). The claim was correct
anyway: **the titration exists in the raw data and had been erased by the
extraction.** Three separate defects, all now fixed.

### 1. Experiments 32, 34, 35, 36, 37 — `[enz]` and `[buf]` both wrong

These five are phosphate / 4OMe-BnOH / 40 °C buffer-concentration titrations at
constant `[sub]`, `[h2o2]`, `[enz]`, pH and T. One bug corrupted two columns:

Their `.xls` sheets lay out **eight planned cuvettes** — rows 1–4 with
`Enz 0.1 ml`, rows 5–8 with `Enz 0` — while only four channels were ever
measured. These were no-enzyme days, so the four that ran were cuvettes **5–8**,
but `find_numeric_values_below_header` reads the *first four rows* of the table
and so picked up the with-enzyme plan rows.

- **`[enz]`** was extracted as 0.241 / 0.270 mM. The runs were enzyme-free. Two
  independent confirmations: the filenames carry `with_NO_E`, and all five sit
  in the hand-sorted `data/Mads/"No enzyme"/` folder. That folder holds 27
  experiments; **22 extract correctly as `[enz] = 0` and exactly these five did
  not.** As extracted they were the *highest*-enzyme runs in the whole dataset.
- **`[buf]`** was extracted as a flat 50 mM, because every cuvette receives the
  same buffer *volume* (1 ml into 2 ml total) and only the stock differs — and
  the stock exists solely as a text label in the `kuv` column (`1 (0.1M)`,
  `2 (0.2M)`, …), invisible to the volume-based extraction.

Corrected cuvette concentrations (= half the stock, from the confirmed 1 ml →
2 ml dilution):

| experiments | `[enz]` | `[buf]` per sample |
|---|---|---|
| 32, 35, 36, 37 | 0 | 50, 100, 150, 200 mM |
| 34 | 0 | 25, 12.5, 6.25, 3.125 mM |

A prior partial patch in `clean_experiment_dataframe` set `[100, 200, 300, 400]`
and `[50, 25, 12.5, 6.25]` — those are the **stock** concentrations, uniformly
2× the cuvette values every other row in the dataset uses — and never touched
`[enz]` at all. Replaced.

**Why this matters more than any other correction so far.** Corrected, these
five are enzyme-free buffer titrations spanning **3.125–200 mM (64-fold) at
constant substrate**. They are simultaneously the `E0 = 0` data needed to fit
the catalyst-free loop (`k_can`, `k3` in `MECHANISM.md`) and the only real
evidence on buffer catalysis anywhere in the project — with no catalyst present
to confound either. Scope caveat: 4OMe-BnOH at 40 °C, pH 7.00–7.53, so the
buffer constant they pin is for the methoxy substrate at that temperature.

Fixed in `data/kinetics_io.py` as `EXPERIMENT_CORRECTIONS` +
`apply_experiment_corrections()`, applied by both
`populate_experimental_data_from_directory` and `load_experiment`; mirrored
idempotently in the notebook's `clean_experiment_dataframe` (cell 16).
`data/experiment_data.csv` regenerated — a full recompile changed **36 cells and
nothing else**: 20 `[enz]` and 16 `[buf]`, all inside these five experiments.

### 2. Ionic strength was silently zero for 26 experiments

`add_ionic_strength_column` keyed its pKa table on `'Boric Acid'` while
`find_buffer_type` emits `'Boric'`, and omitted `'Carbonate'` altogether. Both
fell through to the `return [1]` / `charges = [0]` fallback, giving **`I = 0`
for 107 rows across 26 experiments** with no warning — precisely the high-pH
runs (Boric 8.46–10.34, Carbonate 9.40–11.84) where the correction to
pKa(H2O2), and hence `[HOO-]`, matters most. Since the notebook's feature list
drops `[buf]` and `buffer` and relies on `I` to carry buffer information, for
those 26 experiments the buffer effect was being loaded entirely onto `[sub]`.

Fixed: keys corrected to `'Boric'` (pKa 9.24) and `'Carbonate'` (pKa 6.35,
10.33) with matching charge lists. Two further corrections in the same chain,
without which fixing the keys would have made the numbers worse rather than
better:

- **Counter-ions.** `I = ½Σc·z²` was summed over the buffer anions only. The
  Na+ required for electroneutrality contributes too; omitting it understates
  `I` by roughly a factor of two for a 1:1 salt. Now included.
- **Units and charge factor in the Debye–Hückel step** (cell 38). `I` is
  returned in mM (because `[buf]` is in mM) but was fed straight into a formula
  whose constant `A = 0.509` is defined for mol/L — inflating √I by ~31× and the
  pKa shift by ~9× (1.15 vs 0.13 units at I = 75 mM). Separately, for
  `H2O2 ⇌ H+ + HOO-` the correct term is Δ(z²) = (+1)² + (−1)² − 0 = **2**, not
  z² = 1. Both fixed.

Net effect on the computed `[HOO-]`, which feeds `v_can` directly in the
mechanism: **median ratio new/old = 0.17** (i.e. the old values were ~6× too
high), ranging 0.14–2.02 — the values above 1 being the Boric and Carbonate rows
that previously got no correction at all.

**Known limitation, not fixed:** the corrected `I` reaches ~1.07 M for the
193.6 mM pyrophosphate experiments (127–131). The extended Debye–Hückel form is
only valid to ~0.1 M and the Davies equation to ~0.5 M, so those points are an
empirical extrapolation. Choosing an activity model for that range is a
modelling decision, deliberately left open.

### 3. `NaH2PO4*2H2O` was mapped to Pyrophosphate

The buffer map in `find_buffer_type` contained
`"NaH2PO4*2H2O": "Pyrophosphate"`. NaH2PO4·2H2O is monosodium **phosphate**.
Corrected to `"Phosphate"` in both `data/kinetics_io.py` and notebook cell 10.

No labels changed as a result, because the 17 experiments that hit this key
(135–151) have `Na4P2O7*10H2O` at sheet row 0 and `NaH2PO4*2H2O` at row 2, and
the scan is row-major with early return — so they still resolve to
`Pyrophosphate`. Verified explicitly after the change.

That co-occurrence is itself worth recording: **experiments 135–151 use a
two-component pyrophosphate/phosphate buffer**, so their single `[buf]` value of
75.013 mM stands for two species, and the ionic-strength calculation (which uses
pyrophosphate pKa's only) is incomplete for all 17. They remain the
best-designed block in the dataset — `[buf]` fixed, samples 1–4 titrating
`[sub]` at fixed `[h2o2]`, samples 5–7 titrating `[h2o2]` at fixed `[sub]`,
across 17 pH values from 5.47 to 9.73 — but the buffer composition is not what
the single label implies.

### Design classes, for the record

Recount of what the 98 compiled experiments actually vary:

| design | experiments | |
|---|---|---|
| `[buf]` and `[sub]` displaced together | 49 | confounded (mostly exps 2–62) |
| `[buf]` fixed, `[sub]` varied | 19 | clean (exps 65–85) |
| `[buf]` fixed, `[sub]` + `[h2o2]` varied | 17 | clean (exps 135–151) |
| `[buf]` varied at constant `[sub]` | 5 | **exps 32, 34, 35, 36, 37 — recovered here** |
| `[h2o2]` or `[enz]` only | 5 | |
| nothing varies | 7 | |

The `[buf]`/`[sub]` confound is therefore **quarantined to 49 experiments, not
pervasive**. Note also that only **one** experiment varies `[enz]` across its own
samples, so the E0 dependence remains untestable from this dataset.

### Raw material present but never compiled

49 experiment numbers have raw files under `data/Mads` but no row in the CSV.
Most are deliberate and should stay out — the whole `bad data pH ca. 11` folder
(86–126: the Na2HPO4 pH 11 series and every solvent-isotope run), and 54/56
which are **4-bromo-benzyl alcohol**, a third substrate outside this dataset's
scope. Two groups are worth a second look:

- **exp 134** sits in the hand-sorted `good data BnOH` folder alongside 135–151
  and has an `.xls`, but never reached `data/data` or the CSV. If it belongs
  with that block it is an 18th member of the best-designed series.
- **exps 3, 53, 64** are no-enzyme runs with raw `.xls` in `data/Mads` but no
  `dataN.txt` in `data/data`. (Exp 63 is explained — its own filename ends
  `NO_DATA_FILE`.)

---

## 2026-08-30 — Mechanism research: consequences for how this data can be fitted

Three literature research passes were run to pin down the reaction mechanism
before fitting (full write-up, reasoning and 51 references in `MECHANISM.md`).
Several findings bear directly on the data itself and on what can and cannot be
concluded from it, so they are recorded here too.

### 1. The four buffers are chemically different reagents, not just pH setters

This is the most important finding for data interpretation. It also gives a
mechanistic explanation for the buffer-dependent behaviour noticed by eye in
the dataset.

- **Boric buffer points should be treated as suspect.** Borate does three
  separate things to this chemistry: it forms **peroxoborate** with H2O2
  (significant above pH ~7.7), whose anions are much faster oxidants than H2O2
  and "deliver the hydroperoxide anion at a lower pH than when H2O2 is used";
  it generates **dioxaborirane**, a highly reactive cyclic peroxide that is a
  competing oxidant unrelated to the catalyst; and — most damaging — **boric
  acid catalyses peroxyacid hydrolysis ~12-fold, with a maximum at pH 8.4–9**.
  The dataset's Boric experiments span pH 8.46–10.34, i.e. straight through
  that maximum. If the proposed mechanism's peracid intermediate is real,
  borate buffer is actively destroying it.
- **Carbonate buffer points should also be treated as suspect.** Bicarbonate +
  H2O2 forms **peroxymonocarbonate (HCO4-)**, a two-electron oxidant ~300x
  faster than H2O2 for sulfide oxidation, formed within minutes near neutral
  pH. In the dataset's 7 Carbonate experiments (pH 9.40–11.84) the effective
  oxidant is partly HCO4-, not H2O2.
- **Phosphate catalyses the first step of the mechanism directly.** H2O2
  addition to a carbonyl is subject to **both general acid and general base
  catalysis** (Sander & Jencks 1968), so buffer concentration is a genuine
  kinetic variable, not a nuisance parameter.
- **Pyrophosphate is probably chelating trace metals** (reasoning, not
  sourced) — trace Fe/Cu catalyse H2O2 decomposition, which is why dioxirane
  papers routinely add EDTA. Pyrophosphate chelates; phosphate and carbonate
  do not, to the same degree.

**Consequence:** rate constants from different buffer systems at the same
nominal pH are **not directly comparable**. Any pH-rate profile built by
pooling across buffers is confounded. This should be stated explicitly
alongside any such plot.

### 2. `[buf]` and `[sub]` are collinear within every titration experiment — a buffer-concentration effect cannot be isolated from the existing data

Checked directly: in every titration experiment, `[buf]` decreases in lockstep
as `[sub]` increases across samples 1→4 (substrate stock dilutes the buffer).
Example, experiment 2: `[buf]` 80→70→60→50 mM while `[sub]` 1.52→3.04→4.56→6.08
mM. 50 of 98 experiments vary `[buf]` within the experiment, and **not one of
them holds `[sub]` fixed while doing so**.

So any within-experiment "buffer effect" is perfectly confounded with the
substrate effect. Separating them requires either a multi-variable regression
of v0 against `[sub]`, `[h2o2]`, `[enz]`, `[buf]` jointly across experiments,
or (cleanly, but this needs new bench work) a buffer-dilution series at fixed
pH, fixed `[sub]` and fixed ionic strength.

Buffer *type* is the cleaner comparison available now: the pH ranges overlap
across buffers (Phosphate 5.64–8.95, Pyrophosphate 5.47–9.73, Boric 8.46–10.34,
Carbonate 9.40–11.84), so different buffer species can be compared at matched
pH — the classic diagnostic for general acid/base catalysis. Bear §1 in mind
when interpreting the result, though: borate and carbonate bring their own
chemistry, so a buffer-type difference is not automatically general acid/base
catalysis.

### 3. The pH range may not be able to discriminate the mechanism's key branch

The mechanism predicts that dioxirane formation (its central catalytic step) is
pH-controlled, because it needs an anionic peroxide species. But the same pH
dependence arises trivially from HOO- (pKa 11.6) being the nucleophile in the
*addition* step. Both predict rate rising across pH 5.5→11.8. The discriminator
is **where the inflection sits**: near 11.6 points to peroxyanion
nucleophilicity; well below 10 would point to ionization of the tetrahedral
adduct. Worth fitting the pH profile carefully enough to locate the inflection
— but note that no pKa has ever been measured for the relevant adduct, so this
is an open question, not a calibrated test.

### 4. Reminder: the differential-measurement design still governs everything

Nothing in the mechanism work changes the round-2 §4 finding — with-enzyme
`[P]` is already background-subtracted at the source. But it is now clearer
*why* the background is substantial and accelerating: the proposed
catalyst-free loop (aldehyde + HOO- → peracid → oxidises more alcohol) is
autocatalytic on its own, with no catalyst required. Whether that specific
ionic mechanism is right is unresolved (it has no literature precedent — see
`MECHANISM.md`), and a radical/O2-chain alternative is equally consistent with
the data and cannot be excluded without dark/anaerobic controls.

## 2026-08-29 — Round 3: row-level block-structure audit (all 98 experiments)

Round 2 confirmed `.txt`↔`.xls` sample alignment indirectly (structurally on
two examples, statistically via rate-vs-concentration correlation on 12
experiments). This round checks it **directly, on every experiment**, now
that the paired-reference-cuvette design (round 2 §4) is understood: for
each of the 98 experiments, located the `[Enz]` concentration column in the
raw `.xls`, read its full numeric run (not just the truncated
`sample_num`-length slice the pipeline keeps), and split it into the
"leading" block (what the pipeline actually extracts) vs. whatever follows
(the reference block, if recorded).

Checked for two failure modes that would be serious silent labeling errors:
the leading block containing a mix of zero/nonzero `[Enz]` (would mean a
reference row leaked into the real samples), and an all-zero leading block
sitting in front of a nonzero trailing block (would mean the block order is
reversed and the pipeline is reading the reference cuvette instead of the
real one).

**Result: 97 of 98 experiments clean.** No reversed blocks, no reference
row ever extracted in place of a real sample. Two apparent flags turned out
to be benign and are noted for completeness:

- **Exp 6, exp 65** ("no_E" designs): these have no separate `[Enz]`
  concentration column at all (enzyme volume is 0 throughout, so no
  concentration is computed), or only a single shared reference row instead
  of one per sample. Both fall back correctly to `[enz] = 0` — already
  understood from round 2, not a new issue.

- **Exp 128 — real finding, not benign.** Its leading (extracted) block is
  `[0.032, 0.032, 0.032, 0.032, 0.0]` — 4 real enzyme-containing values
  followed by a **0** still inside the block the pipeline keeps, because
  this `.txt` file declares **5** samples while its sibling experiments in
  the same series (127, 129, 130, 131 — same `sub_H2O2` design, same day)
  all declare 4. Read the raw sheet directly (`mads_t128_..._sub H2O2_
  pyrophosphate_02.xls`, rows 32–39): it has the standard 4 "kuv" (real)
  rows followed by 4 "ref" (matched, no-enzyme) rows — but unlike other
  experiments, where the reference block is written to the sheet purely as
  a record and never separately monitored over time, **this run's 5th
  measured channel (`Sample005`) is the first reference row itself
  (`ref 5`)**, not a 5th titration condition. Its `[enz]=0` and its
  `[buf]`/`[h2o2]`/`[sub]` (193.6 / 195.9 / 9.47) exactly match `ref 5`,
  which is the matched blank for `kuv 1`/`kuv 2` (H2O2 = 195.9 mM).

  Current `data/experiment_data.csv` (rows 155–159) already reflects this:
  sample 5 has `[enz] = 0.000` while samples 1–4 have `[enz] = 0.032`.
  Samples 2 and 3 of this experiment were already excluded by the existing
  `clean_experiment_dataframe` (pre-dating this verification effort, for an
  unrelated reason — a backwards/negative time-series trend, see the
  reaction-direction check above). **Resolved in the follow-up investigation
  below** — see "Round 4" for the final disposition of all five samples.

---

## 2026-08-29 — Round 4: plotted experiment 128 directly, resolved samples 2/3/4/5

Built two small reusable modules in `data/` for this and future
investigations, rather than another one-off scratchpad script:

- `data/kinetics_io.py` — the notebook's extraction functions (cells 3,
  5–12 of `masterThesis.ipynb`), copied verbatim so there is one canonical
  implementation, plus a new `load_experiment(exp_num)` convenience wrapper
  that returns one experiment's metadata *and* raw time series together
  (the whole-dataset `populate_experimental_data_from_directory` returns
  metadata only).
- `data/plot_kinetics.py` — `plot_experiment(exp_num, mark_samples={...})`,
  raw absorbance and Beer-Lambert-converted `[P]` side by side, one fixed
  color per sample index, optional dashed styling + legend annotation for
  samples under investigation. Runnable as `python data/plot_kinetics.py
  <exp_num>` from the repo root, or imported directly.

Plotting experiment 128 (`python data/plot_kinetics.py 128`) showed all 5
samples together for the first time and changed the picture from round 3:

| sample | condition | shape | current status |
|---|---|---|---|
| 1 | [enz]=0.032, H2O2=196mM | clean rise, 0.030→0.047 | kept |
| 2 | [enz]=0.032, H2O2=196mM (same as 1) | flat/noisy, net trend −0.008 | **excluded** |
| 3 | [enz]=0.032, H2O2=3.88mM | flat/noisy, net trend −0.004 | was excluded, **reinstated** |
| 4 | [enz]=0.032, H2O2=3.88mM (same as 3) | flat/noisy, net trend +0.003 | kept |
| 5 | [enz]=0, H2O2=196mM (ref channel, see round 3 §4) | flat, net trend +0.002 | **excluded** |

The key observation: samples 3 and 4 are two *independent replicates* of
the same low-H2O2 condition, and both go flat — reproducible, not a
one-off. Meanwhile at the high-H2O2 condition, one replicate (1) works
cleanly and the other (2) doesn't, on the same enzyme aliquot — that's the
signature of a single failed measurement, not a real effect (the aliquot
demonstrably works, since sample 1 proves it).

**Decision:** treat 3.88 mM H2O2 as a genuine, reproducible near-zero-rate
condition rather than two failed measurements. Sample 3's original
exclusion (under the general "any net-negative time-series trend is
unphysical" rule from the earlier reaction-direction check) was conflating
a small, noise-scale negative trend with the large, clearly-wrong-sign
trends that rule was designed to catch (e.g. exp 50/85 pre-correction, or
58/77/78/79/84). Reinstating sample 3 alongside 4 avoids cherry-picking
whichever replicate's noise happened to land on the "nice" side of zero.
Sample 2 stays excluded — unlike 3/4, it disagrees with its own matched
replicate (1) at a condition proven capable of strong signal, which is a
single-measurement failure, not a low-rate result.

**Fix applied** to `clean_experiment_dataframe` (`masterThesis.ipynb`,
cell 16): `exp_num_sample_num` changed from `[[128, 2], [128, 3]]` to
`[[128, 2], [128, 5]]` — drops the sample-3 exclusion, adds the sample-5
(reference channel) exclusion from round 3. Net effect on experiment 128:
samples 1, 3, 4 feed the fit as real data points; 2 and 5 are excluded.

**Open question flagged, not yet chased:** the same "small negative trend
excluded under the blanket backwards-trend rule" pattern that applied to
128,3 may exist elsewhere in the dataset — worth a follow-up pass
distinguishing large-magnitude (real sign-flip) exclusions from
small/noise-scale ones before finalizing the excluded-sample list.

---

## 2026-08-29 — Round 5: broad pass on all "backwards trend" exclusions — closed, one mischaracterization fixed

Followed up on round 4's open question: checked every sample in
experiments 58, 77, 78, 79, 84 (the ones grouped together as pre-existing
"backwards trend" removals) plus 50/85 (sign-flip corrected, as a sanity
check) with a proper significance test instead of eyeballing the sign of
the net trend. For each sample, ran linear regression (`scipy.stats.
linregress`) of raw absorbance vs. time and compared the fitted slope to
its standard error (t-test for slope ≠ 0), plus amplitude vs. the residual
noise level. The question: is each "backwards" trend a real, high-SNR
inverted signal (genuine error, correctly excluded) or a small one
indistinguishable from noise (candidate to reinstate, per the 128,3
reasoning)?

**Result: no new reinstatement candidates.** 128,3 was the exception, not
part of a wider pattern:

- **77, 78:** every sample, large and highly significant negative slopes
  (t down to −95, amplitude 0.02–0.06 vs. noise ~0.001–0.003). Real,
  high-confidence inverted signal — correctly excluded. (Also Carbonate
  buffer, so excluded via that filter regardless of this list.)
- **58:** mixed — 3 samples strongly negative, 1 strongly positive, all
  large real signal (not noise), internally inconsistent between samples
  (possibly a cuvette-position mixup). Correctly excluded; also Carbonate.
- **79:** all 4 samples tiny amplitude (0.001–0.003, at the noise floor) —
  but unlike 128, there's no sibling sample in the group with clean signal
  to contrast against, so this reads as a group dead-run (same flavor as
  72/82) rather than "noise obscuring a real, localized effect." Left
  excluded; also Carbonate.
- **50, 85** (sign-flip corrected, not excluded): every sample in both,
  large amplitude (0.09–0.63), highly significant, consistently inverted.
  Confirms the existing sign-flip fix is well-founded, not a noise
  artifact.

**One mischaracterization found and fixed here:** round 1's original audit
grouped **84** in with the "backwards trend" removals, and round 4 repeated
that framing. It's wrong — checked directly, 84's one recorded sample has
a small but *correctly-signed, positive* trend (t = +22), not a backwards
one. Plotting it (`python data/plot_kinetics.py 84`) shows why it's still
correctly excluded, for an entirely different reason: its `.xls` shows a
planned 6-sample substrate titration (14.3→1.0 mM `[sub]`), but the `.txt`
only ever recorded **one** of those six samples, and that lone curve is a
coarse, quantized 5-step staircase over just 385 s (56 points, dt=7s —
every other experiment in this set uses dt=28–31s). This is an incomplete/
truncated run, not usable kinetic data, regardless of its trend direction.

No code changes from this round — every disposition above already matches
what `clean_experiment_dataframe` currently does. This closes the open
question from round 4.

---

## 2026-08-29 — Round 2: concentration/curve cross-checks

A second verification pass, this time checking things that need *both* data
sources together (round 1 checked the csv against the xls, and the raw
curves against themselves, but not whether a given `.txt` sample and its
`.xls` row are actually the same physical sample).

### 1. `.txt` sample ↔ `.xls` row alignment — verified, no mismatches found

The extraction pipeline assumes the *n*-th `Sample00n` in a `.txt` file is
the *n*-th data row under the concentration header in the matching `.xls`
sheet. This was never directly checked before. Two lines of evidence:

- **Structural**: many `.xls` sheets (e.g. exp 5, exp 135–151) list each
  condition twice — a block of "Kuv." (cuvette, the real monitored sample)
  rows immediately followed by an equal-sized block of "Ref." (no-enzyme
  reference/blank) rows in the *same* concentration columns. Hand-checked
  exp 5 and exp 135 directly against the raw sheet: the `.txt` file only
  declares as many samples as there are Kuv rows, and the pipeline's
  truncate-to-sample-count logic grabs exactly that leading Kuv block, in
  order — it does not spill into the Ref block or reorder anything.
- **Statistical**: for every experiment with exactly one concentration field
  varying across its samples (a clean titration), computed each sample's
  initial reaction rate (linear fit, first 30% of the curve) and checked
  its Spearman correlation with the varied concentration. If `.txt` order
  and `.xls` row order were ever out of sync for a given experiment, this
  would show up as a scrambled or inverted (rate falls as substrate rises)
  correlation. Result, 12 substrate-titration experiments (65, 66, 68, 69,
  70, 71, 73, 74, 75, 76, 83 — excluding the already-known no-enzyme
  controls) all show positive correlation (rho 0.6–1.0), consistent with
  expected Michaelis-Menten behavior. No inverted or scrambled cases.

  A separate "numeric run length" check was tried first (count how many
  numeric cells actually sit under each concentration header, independent
  of sample count, to catch a mismatch directly) but it produces mostly
  false positives, because of the legitimate Kuv+Ref double-block layout
  above — abandoned in favor of the two checks described here.

### 2. Full range/outlier sweep — clean

Checked all 405 non-removed sample rows: pH in [2, 12], T in [0, 80] °C, and
no negative `[enz]`/`[buf]`/`[h2o2]`/`[sub]` values. Zero violations.

The only rows flagged were 72 samples across experiments 6, 23–31, 38–40,
52, 65, 67, 69, 70, 128 (`[enz] = 0` for every sample in the experiment).
Spot-checked several filenames directly — e.g.
`mads_t065_..._BnOH_no_E_H2O2_122.xls` — confirming these are deliberate
no-enzyme blank/background-control runs, not an extraction failure.

### 3. Baseline check (`[P]` near zero at t=0) — clean

Converted each curve to `[P]` and checked the first point sits near the
curve's minimum. 3 borderline flags out of 405, all in samples with a tiny
absolute amplitude (near the noise floor at the lowest substrate
concentration in a titration) — not real calibration offsets, just a small
denominator inflating the fractional metric.

### 4. Explained: `[P]` is a *differential* (reference-subtracted) signal, not an absolute one — this is why no-enzyme "blank" runs can rate as high as or higher than their catalyzed pair

Initially flagged as an anomaly (no-enzyme blanks 65/67/69/70 showing rates
at or above their enzyme-containing pairs 66/68/71) and logged as an open
item. Root cause identified by reading the raw cuvette layout directly, not
just the extracted concentration columns — **every kinetics run in this
dataset is a two-channel differential measurement, not a single absolute
one.** Each `.xls` sheet lists twice as many cuvette rows as the matching
`.txt` file has samples: a leading block of "real" sample rows (labelled
`kuv`/`Kuv.`/`prøve`), immediately followed by an equal-sized block of
reference rows (labelled `ref`/`Ref.`). The instrument (or the operator,
manually) reads each sample cuvette *against* its paired reference cuvette,
and it is that difference — not the sample's raw absorbance — that ends up
in the `.txt` progress curve. The extraction pipeline only ever reads the
leading `sample_num` rows (the real samples); the trailing reference block
is present in the spreadsheet purely as a record of what the reference
cuvette contained, and is correctly never pulled into `[enz]`/`[buf]`/
`[h2o2]`/`[sub]`.

**What differs between experiment types is what the reference cuvette omits:**

- **"with_E" experiments (the bulk of the dataset, e.g. exp 10, exp 66):**
  the reference row for each sample has **identical `Sub [ml]` and
  `H2O2 [ml]` to that sample**, and only `Enz [ml]` set to 0. Example, exp 10
  (`mads_t010_..._with_E.xls`), rows 19–26:

  | kuv | Buf[ml] | H2O[ml] | Enz[ml] | Sub[ml] | H2O2[ml] | role |
  |---|---|---|---|---|---|---|
  | 1–4 | 1.6→1.0 | 0 | 0.1 | 0.2→0.8 | 0.1 | real sample (enzyme present) |
  | 5–8 | 1.6→1.0 | 0.1 | **0** | 0.2→0.8 | 0.1 | matched reference (no enzyme, same Sub/H2O2) |

  So the reported `[P]` for a with_E sample is **already net of the
  non-enzymatic background** — the enzyme-free reaction is running in the
  paired reference cuvette under identical Sub/H2O2 conditions and is
  subtracted out by construction. This is the design in exp 5, exp 10, exp
  66, exp 135–151, and appears to be the general pattern for the with-enzyme
  half of the dataset.

- **"no_E" experiments (65, 67, 69, 70 — standalone background-characterization
  runs, not part of a with/without-enzyme pair in the same sheet):** the
  reference cuvette is missing a *different* reagent, not the enzyme (there
  is no enzyme in either channel):
  - Exp 65 (Boric, pH 8.51): reference has `Sub [ml] = 0` (water makes up the
    volume), but **H2O2 is still present** (`H2O2 = 122.4 mM`, same as the
    sample). So the reference isolates "does H2O2 alone drift?" — it does
    **not** subtract the substrate+H2O2 background reaction, since the
    reference has no substrate for that reaction to run on.
  - Exp 67/69/70 (Phosphate, pH 8.01): reference has `H2O2 = 0`, but
    **substrate is still present** at the same concentration as the sample.
    So here the reference isolates "does substrate alone drift?" — it does
    **not** subtract the substrate+H2O2 background reaction either, since
    the reference has no H2O2 to drive it.

  Either way, the `.txt` curve for a no_E experiment is the **raw,
  undiminished non-enzymatic (substrate + H2O2) background reaction** —
  nothing chemically equivalent to it is subtracted. This background is
  expected to be substantial in its own right and likely autocatalytic
  (accelerating over the course of the reaction rather than a small linear
  drift), which further explains why its measured rate can look as large
  as, or larger than, the *net* (background-subtracted) rate reported for
  a with_E experiment.

**Practical implications for rate-expression fitting:**

- `[P]` for with_E samples throughout the main dataset should be treated as
  **already background-corrected** — do **not** additionally subtract a
  no_E curve from a with_E curve; that would double-subtract, since the
  with_E curve's own paired reference already performed that subtraction
  at the instrument level.
- The four no_E experiments (65, 67, 69, 70) are **not on the same
  reference basis** as the rest of the dataset (no matched enzyme-free
  reference of their own) and are not directly comparable, sample-for-sample,
  to any with_E experiment's raw rate. They exist to characterize the size
  of the background reaction on its own terms, not to be arithmetically
  combined with catalyzed runs.
- This was previously logged as an unresolved "open item" — it no longer
  is. No further action needed on this point.

---

## 2026-08-29 — Raw [P] time-series audit

Audited the raw progress-curve data (`data/data/*.txt`, parsed by
`parse_experiment_data`) that feeds the `[P]` column, independently of the
metadata audit above.

### 1. Structural integrity — clean

Reparsed all 98 `.txt` files (443 sample-series) standalone. No parse failures,
no empty series, no backwards or duplicate timestamps, uniform sampling interval
within every series. Point counts range 10–481, all plausible.

### 2. Reaction-direction check — confirms existing fixes, no new cases

Flagged every sample series with a net-negative absorbance trend (signal
decreasing over time, which is wrong for product formation). This reproduced
*exactly* the set already special-cased in `clean_experiment_dataframe`:
experiments 50 and 85 (already sign-flipped), 58/77/78/79/84 (already removed),
experiment 128 samples 2 and 3 (already removed), and experiment 80 (already
dropped by the "remove all Carbonate buffer" filter). No new sign-flip
candidates found.

### 3. Physical-plausibility check — clean

Converted each curve's amplitude to Δ[P] via its extinction coefficient (`e`)
and confirmed it never exceeds the sample's initial `[sub]` (product can't
exceed starting substrate) across all 443 samples. Zero violations.

### 4. Experiments 32 and 82 are dead runs — removed

Experiments 72 and 82 show essentially zero signal (amplitude 0.001–0.008
absorbance units, vs. a dataset median of 0.044) in *every* sample, despite
normal, nonzero `[enz]`/`[h2o2]`/`[sub]` per the (trusted) xls values. Compared
against same-day, same-instrument sibling runs:
- Exp 72 (6/8/2010, 7:13pm) vs exp 71 (same day, 4:51pm): exp 71 shows real
  signal (amplitude 0.017–0.031) at *lower* substrate concentrations than exp 72
  used. Exp 72 is flat despite higher substrate.
- Exp 82 (6/14/2010, 1:55pm) vs exp 83 (same day, 5:03pm): exp 83 shows strong,
  clean signal (amplitude 0.031–0.043) at comparable/lower substrate. Exp 82 is
  flat despite comparable/higher substrate.

Same-day neighbors with comparable reagents produced normal curves, so this
isn't a bad-reagent-stock day or a data-entry problem — the working hypothesis
is a failed individual run (e.g. an inactive enzyme aliquot), the same failure
mode presumed for the already-excluded 58/77/78/79/84.

**Fix applied:** added 72 and 82 to `experiments_to_remove` in
`clean_experiment_dataframe`, alongside 58/77/78/79/84. See
`masterThesis.ipynb`, cell 16.

### 5. Open item, not chased: `.txt` "Substrate Conc." line disagrees with xls `[sub]`

Each `.txt` file has its own `Substrate Conc.` annotation, independent of the
`.xls` parsing and **not read by the pipeline at all**. Cross-checked it against
the xls-derived `[sub]` anyway as a bonus sanity check: 51 of 313 comparable
rows disagree, with no single clean pattern (includes a clear raw transcription
typo — exp 30 sample 3 literally reads `2856.0000 mmol/l` against neighbors
`0.0953, 0.1906, _, 0.3812`, almost certainly meant `0.2856`).

**Disposition: the `.xls` files are the source of truth for concentrations —
this `.txt` field was not reliably filled in and should not be trusted.** Not
worth further investigation since the pipeline never reads it; logged here only
in case it resurfaces.

---

## 2026-08-29 — Metadata sweep + buffer-titration fix

### 1. `data/experiment_data.csv` verified against the raw `.xls` files

Re-implemented the notebook's own extraction pipeline (`parse_experiment_data`,
`find_and_parse_experiment_file`, `find_header_row`, `find_numeric_values_below_header`,
`find_pH_value_in_range`, `find_temperature_value_in_range`, `find_buffer_type`,
`find_substrate_type` — notebook cells 3, 5–11) as a standalone script and ran it
fresh against all 98 experiment `.txt`/`.xls` pairs in `data/data/`.

**Result: all 443 rows match `experiment_data.csv` exactly, field-by-field
(`substrate`, `abs`, `e`, `buffer`, `pH`, `T`, `[enz]`, `[buf]`, `[h2o2]`, `[sub]`).
Zero mismatches.** Spot-checked the raw cell layout of experiment 10 by hand to
confirm the extraction logic itself (not just internal consistency) is reading the
right cells.

No warnings or extraction failures were emitted across any of the 98 files.

### 2. `data/data.toml` is stale and unreliable — do not use

- Does not parse as TOML: 8 experiments have `pH = [6.71, 6.71, 6.71. 6.71]` (a
  stray period instead of a comma), and `experiment_022.references` has
  `pH = [.50, ...]` (missing leading digit).
- Not referenced anywhere in `masterThesis.ipynb` — the notebook builds its
  dataframe directly from `data/data/`, never from `data.toml` or
  `experiment_data.csv`.
- Where it disagrees with the raw `.xls`/filename ground truth, it's the one
  that's wrong:
  - exp 2: toml says pH 7.0; raw filename and csv both say 6.71 (correct).
  - exp 38: toml lists 4 samples; the raw `.txt` file only has 3
    (`Sample001–003`) — toml has a phantom 4th sample.
- Decision: **not fixing this file for now.** It's slated for a full rewrite, so
  effort goes into the source data instead. Do not treat it as authoritative in
  the meantime.

### 3. Buffer-titration experiments were silently flattened — fixed

Experiments 32, 34, 35, 36, 37 are buffer-concentration titrations: their `.xls`
files label each sample `1 (0.1M)`, `2 (0.2M)`, `3 (0.3M)`, `4 (0.4M)` (confirmed
with the user: these are the actual buffer **stock** concentrations used). But
that molarity only exists as text in the label column — there's no numeric `[buf]`
column like there is for `[sub]`/`[h2o2]`. The generic `[buf]` extraction fallback
(`find_numeric_values_below_header`, notebook cell 7) instead derives `[buf]` from
sample volume assuming a fixed 100 mM buffer:
```python
adjusted_value = round((value / 2) * 100, 3)  # assumes the buffer is 100 mM
```
Since the pipetted volume is constant across all 4 samples in these experiments
(only the stock concentration changes), this produces an identical, wrong
`[buf] = 50.0` for every sample, regardless of which stock was actually used.

`clean_experiment_dataframe` (notebook cell 16) already had a hardcoded override
for this exact bug on experiments 34 and 35/36 — but 32 and 37 were missed and
still carried the flat wrong value.

**Confirmed these were meant to vary, not a mislabeled flat-concentration run.**
This experiment type is a *buffer-concentration titration*, distinct from the
bulk of the dataset (which titrate substrate at a nominal fixed buffer level).
Evidence:
- Filenames encode a list of four distinct concentrations, not a single value —
  `Phosphat_0.1_0.2_0.3_0.4` (32, 35, 36, 37) and
  `Phosphat_0.05_0.025_0.0125_0.00625` (34, a 2-fold serial dilution) — unlike
  ordinary experiment filenames, which name just one buffer.
- The raw sheet labels each sample individually: `1 (0.1M)`, `2 (0.2M)`,
  `3 (0.3M)`, `4 (0.4M)`.
- `Sub [ml]` and `Buf [ml]` volumes are both held constant across all 4 samples
  (e.g. `Sub [ml] = 0.8`, `Buf [ml] = 1` for every sample) — the opposite of a
  normal substrate titration (e.g. exp 10), where `Sub [ml]` climbs
  (`0.2, 0.4, 0.6, 0.8`) and `Buf [ml]` shrinks to compensate
  (`1.6, 1.4, 1.2, 1.0`). Here substrate is pinned and only the buffer stock
  pipetted in differs sample-to-sample — the signature of buffer concentration
  being the deliberately titrated variable.
- User confirmed 0.1M/0.2M/0.3M/0.4M as the actual buffer stock concentrations
  used.

**Fix applied:** extended the existing override in `clean_experiment_dataframe`
to cover experiments 32 and 37 as well, using the same `[100, 200, 300, 400]` mM
convention already established for 35/36. See `masterThesis.ipynb`, cell 16.

### 4. Full sweep for the same class of bug across all 98 experiments

Checked every experiment for any concentration field (`[enz]`, `[buf]`, `[h2o2]`,
`[sub]`) that's flat across all samples, then checked whether the filename/raw
label implies it should vary (the tell-tale sign of this bug class).

- **Experiment 26**: all four fields flat. Checked the raw sheet directly —
  genuinely a 4-way replicate design at fixed conditions (no enzyme; comparing
  "with H2O2" vs "without H2O2" blanks), not a titration. Not a bug.
- **Experiments 32, 34, 35, 36, 37**: buffer titrations — fixed (see §3).
- **Experiments 127–131** ("sub H2O2" filenames): `[sub]` flat, `[h2o2]` varies —
  correct by design, these substitute an H2O2 concentration range for the usual
  substrate range.
- No other experiment showed an unexplained flat field. No further metadata bugs
  found in this pass.

### 5. Experiment 32 has a stray duplicate `.xls` file

`data/data/` contains two files for experiment 32:
`mads_t032_..._Phosphat_0.1_0.2_0.xls` (looks like a truncated/interrupted save)
and `mads_t032_..._0.1_0.2_0.3_0.4_..._with_NO_E.xls` (the complete one). Checked
both — they contain byte-identical concentration tables, so whichever one the
pipeline's glob picks doesn't affect results. Not fixed/removed, just noted here
in case the stray file causes confusion later.

---

## Open items

- ~~**Exps 57 and 58 are 4-bromobenzyl alcohol runs.**~~ **CLOSED 2026-08-30** —
  ruled 4OMe-BnOH, and both excluded: their substrate stock block was copied
  from t056 and the real preparation is unrecorded, so `[sub]` is not
  recoverable. See the entry at the top of this log.
- ~~**The header `kuv` row is an unused independent check on `[enz]`.**~~
  **CLOSED 2026-08-30** — built as `data/verify_enzyme.py` and folded into
  `--deep`. All 98 experiments now trace `[enz]` back to a weighed mass of
  catalyst.
- ~~**`Rate(pH).xls` carries an `[E] [mmol/l]` column per pH.**~~
  **CLOSED 2026-08-30** — built as `data/verify_rate_workbook.py` and folded
  into `--deep`; all nine tabulated points confirm the dataset.
- **Exp 6 uses M = 109.13 for benzyl alcohol** where every later BnOH run uses
  the correct 108.14, making its `[sub]` 0.9% low. The only compiled experiment
  affected (t001 and t003 share it but were never compiled). Recorded, not
  corrected — it is well inside the other uncertainties on that run.
- ~~**t054 and t056 are genuine 4-bromobenzyl alcohol runs and are not
  compiled.**~~ **CLOSED 2026-08-30 — impossible, not pending.** Neither has an
  instrument export, so the bromo arm cannot be built from this archive at all.
  The Hammett series stays at two substrates.
- **The cross-substrate ε ratio 7.53/1.23 is a standing assumption.** Within a
  substrate an error in `e` is one global scale factor and harmless; between
  substrates it enters any comparison of rate constants directly. Not a defect,
  but it should be stated wherever a substrate effect is claimed.
- **The extended Debye–Hückel equation is out of range for 70% of the dataset**
  (308 rows, 71 experiments above I = 100 mM; 21 rows above 1 M, all
  pyrophosphate). Systematic, not random, and largest where `[HOO-]` is
  largest. A Davies or Pitzer treatment is the fix; until then any
  `[HOO-]`-dependent result should state its exposure via
  `solution_chemistry.out_of_range_fraction()`.
- **The buffer stock molarity is an assumption for 33 experiments** (32 plus
  exp 3, added 2026-08-30, which is the same early phosphate era). `[buf]` is
  read from a declared `[buf] mmol/l` column in 42, and inferred as
  `(V_buf/2)·100` in the other 56 — hardcoding a 2 ml cuvette (true in 55 of 56;
  exp 58 is 2.1 ml and already excluded) and a 0.1 M stock. That stock is now
  confirmed by a weighed recipe in 1 (exp 13, computing to 100.5 mM), by the
  filename in 18, and by in-sheet text in 5 (exps 32, 34–37). It is stated
  **nowhere** for the remaining 33 — exps 2, 3, 4–12, 14–31, 38, 39, 40 and 52,
  which are all phosphate and all predate the filename convention that starts at
  t041. `I` and `[HOO⁻]` scale linearly with any error there. Tracked as the
  `buf_provenance` column.
- `data.toml` still doesn't parse and isn't reconciled with the verified csv
  (deprioritized — planned rewrite of the toml-generating script).
- Metadata (pH, T, buffer, substrate, `[enz]`/`[buf]`/`[h2o2]`/`[sub]`), the
  raw `[P]` time-series data, and the cross-linkage between the two (sample
  alignment, range sanity, baseline, rate-vs-concentration monotonicity)
  have now all been audited (see entries above). Remaining known-unreliable
  field: the `.txt` files' own `Substrate Conc.` line — not used by the
  pipeline, don't trust it if it resurfaces.
- **Read this before fitting rate expressions:** `[P]` for with-enzyme
  samples throughout the dataset is a *differential* signal — each sample
  cuvette is measured against a paired, matched, no-enzyme reference
  cuvette (same `[sub]`/`[h2o2]`), so the non-enzymatic background is
  already subtracted at the source. Do not subtract a no-enzyme experiment
  (65, 67, 69, 70) from a with-enzyme one — that double-subtracts. See
  round 2 §4 for the full explanation and raw-cell evidence.
- Row-level block structure (which rows in each `.xls` the pipeline
  actually extracts vs. the paired reference rows it correctly leaves
  behind) has now been checked directly for all 98 experiments — see
  round 3. 97/98 clean; experiment 128 was the one exception, resolved in
  round 4 (samples 2 and 5 excluded, 3 and 4 kept as a reproducible
  near-zero-rate replicate pair).
- ~~Round 4 open question about other "backwards trend" exclusions~~ —
  **closed in round 5.** Checked every sample in 58/77/78/79/84 (plus 50/85
  as a sanity check) with a proper significance test; 128,3 was the only
  noise-scale case. One mischaracterization fixed: 84 was never actually a
  backwards-trend case — it's an incomplete run (1 of a planned 6-sample
  titration was ever recorded, and that one is a coarse quantized-staircase
  artifact) — but it's still correctly excluded either way. No code changes.
- `data/kinetics_io.py` and `data/plot_kinetics.py` now exist as reusable,
  single-source-of-truth tools for loading/plotting one experiment at a
  time (`python data/plot_kinetics.py <exp_num>`) — use these instead of
  writing new one-off extraction scripts for future investigations.
- **Buffer systems are not interchangeable** — borate, carbonate, phosphate
  and pyrophosphate each bring their own chemistry to an H2O2 reaction (see
  the 2026-08-30 mechanism entry, §1). Boric and Carbonate points in
  particular should be flagged as suspect in any pooled analysis. Rate
  constants from different buffers at the same pH are not directly
  comparable.
- **A buffer-concentration effect cannot be isolated from the current data**
  — `[buf]` and `[sub]` are perfectly collinear within every titration
  experiment (50 of 98 experiments vary `[buf]`; none hold `[sub]` fixed
  while doing so). Buffer *type* at matched pH is the comparison that is
  available; buffer *concentration* needs either a joint multi-variable
  regression or new bench work.
- `MECHANISM.md` now holds the proposed 7-step mechanism, its per-step
  literature support (51 references), the competing side reactions that
  belong in the ODE model as sink terms, and the open questions. Read it
  before building the kinetic model.

- **Experiments 32/34/35/36/37 are now the highest-value block in the dataset**
  (enzyme-free, 3.125–200 mM buffer at constant substrate) and have never been
  used as such — the prior analyses saw them as top-of-range enzyme runs at a
  flat 50 mM. Any conclusion drawn from `[enz]` or `[buf]` before 2026-08-30
  should be re-derived.
- **All `[HOO-]`-derived results predate the ionic-strength fix** and used
  values ~6× too high (median), with 26 experiments getting no correction at
  all. The RandomForest / PySR feature work on `I` and `[HOO-]` needs re-running.
- **Activity model above 0.1 M is unresolved** — corrected `I` reaches ~1.07 M
  for exps 127–131, outside the validity of both Debye–Hückel and Davies.
- **Exps 135–151 use a two-component pyrophosphate/phosphate buffer**, so their
  single `[buf]` value and their ionic strength are both incomplete.
- **Exp 134, and exps 3/53/64**, have raw material but no compiled row — see the
  2026-08-30 buffer-titration entry for whether they should be recovered.
- **Six manifest open questions need your ruling** — most urgently exp 85's
  substrate (6.1x concentration error if the filename is right) and exp 80's
  enzyme status. Run `python data/validate_dataset.py` to see them.
- **`find_numeric_values_below_header` still takes the first N rows** of the
  concentration table rather than the measured ones — the root cause of both
  the exp 32-37 and exp 128 defects.
- **51 experiments' `[buf]` assumes a 0.1 M stock and cannot be checked from
  the files** (the 56 without a `[buf]` column, less the 5 patched). Exps 32-37
  prove non-0.1 M stocks were used in this series. Resolvable only from lab
  notebooks.
- ~~59 concentration columns unverifiable~~ **CLOSED 2026-08-30** via the
  recorded dilution series - see `data/verify_dilutions.py`.
- **Exps 135-151's buffer is 37.51 mM pyrophosphate + 37.51 mM phosphate**, not
  75 mM of a single species. The ionic-strength calculation still treats it as
  pure pyrophosphate.
- **The sheets record the monitoring wavelength and extinction coefficient**
  (285 nm, e = 1.23 U/mM for BnOH) and exps 135-151 log pH before/after each
  run. Both bear directly on `MECHANISM.md`'s open observable question.
