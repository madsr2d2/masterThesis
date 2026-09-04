# Working in this repository

Master's thesis pipeline: oxidation of benzyl alcohol and 4-methoxybenzyl
alcohol by H<sub>2</sub>O<sub>2</sub>, catalysed by a cyclodextrin–ketone
"chemzyme". `README.md` has the layout; this file has the rules.

## There is no privileged block

**Each question names the block it needs, and every block is defined in
`data/scope.py`** — never as a list of experiment numbers in an analysis
script. `TEMPERATURE_SERIES`, `FREE_BNOH_ALL`, `FREE_BNOH_PHOSPHATE`,
`BUFFER_FIXED`, `buffer_role.TITRATIONS`, `slowdown.substrate_blocks`,
`induction.induction_blocks` are the vocabulary; add to it rather than
selecting inline.

This changed on 2026-09-03. Until then `fit_dataset.TWO_AXIS_BLOCK` named exps
135-151 and the whole project was declared scoped to it — and that had stopped
being true. Four of the five analysis folders work outside it; the temperature
series is the **sole** route to activation parameters; the deceleration needed
84 curves the block does not hold; the buffer question needed five titrations,
none of them in it. Meanwhile the block called "primary" had no folder of its
own and the only saved mechanism fit is on phosphate runs sharing no experiment
with it. A name asserting a priority the work contradicts is a hazard, not a
label.

**Exps 135-151 are now `fit_dataset.TWO_AXIS_BLOCK`**, named for the design
that selects them and checked by `test_fit_kinetics.test_two_axis_block`: the
archive's only runs that move BOTH concentration axes inside a single run —
seven cuvettes each, `[S]` over ≥ 39.9x, `[H2O2]` over ≥ 6.8x, spanning three
pH units and three decades of [HOO⁻]. **Do not call it a chemistry cell.** It
is not "BnOH catalysed" (the archive holds 175 such curves across three
buffers) and not "BnOH catalysed pyrophosphate at 25 °C" either — exps 75 and
76 are exactly that and are not in it. The test's last check is that no run
outside the block carries the design; if one ever does, it fails and forces a
decision.

The block is **not** the hand-sorted `data/Mads/good data BnOH/` folder, which
also holds an excluded run, a run from a different cell, and a sheet with no
instrument data.

**`data/Mads` is the WORKED archive, not the delivered one.** It is what the
pipeline reads (`read_rre.ARCHIVE_DIR`), and 67 of its 717 files differ from
the delivery: mostly a LibreOffice re-save (last-bit floats, added columns,
every original cell intact), but `mads_t057…CO3…with_E.xls` has `'jjjj'` typed
over its experiment number and two computed cells blanked. Nothing published
depends on that -- exp 57 is already excluded -- but the delivered archive is
`Mads-20241207T151327Z-001.zip`, which is the ONLY pristine copy and is tracked
for that reason. Do not delete it as a duplicate; it is not one.

For any analysis of the kinetics data, invoke the **`analyse-kinetics`** skill.

## Do not re-derive measurements

Curve measurements live in `data/curve_metrics.py`; block selection and derived
columns live in `data/scope.py`. Import them.

Six functions were once defined in two modules each and four had diverged. The
lag statistic's copies disagreed on 96 of 402 curves — one would have reported
21% where the evidence said 34%. `data/test_curve_metrics.py` fails if a name
is defined at top level in two modules.

Readings come from the instrument's own `.rre` where one exists and from the
`.txt` export otherwise; `Curve.source` says which, and the floor differs by a
factor of 1096 between them. **All 402 curves are `.rre` today** — the last 32,
exps 2-32, arrived when `read_rre` stopped matching only `rate<n>.rre` — so
nothing in the archive currently sits on the export's floor. The rule below
still stands, because the DEFAULT is the export's and a returning `.txt` curve
would be silent. Never floor a noise **or a residual variance** at
`QUANTISATION_SIGMA` without checking the source: call
`fit_dataset.source_floor(curve.source)` and pass the result. Everything in
`curve_metrics` that floors takes the value as an argument, and the default is
the `.txt` export's.

This is not only about the noise column. `line_fit` floors the variance behind
every standard error in the package, and `acceleration` divides by two of them,
so a floor left at the export's quantisation on `.rre` data suppresses the
z-score it is measured by — that is how the two-axis block's acceleration count read
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
python data/test_scope.py                 # the order machinery and the pH ladders
python data/test_fit_kinetics.py          # selection, scope, parameter recovery (9 min)
python data/test_validator.py             # fault injection
python data/test_slowdown.py              # the slowdown models and their regressions
python data/test_induction.py             # the induction landmark and its controls
python data/test_buffer_role.py           # the species test, planted both ways
python test_doc_check.py                  # the contract every check_numbers runs on
```

Or all of it, which is what `python run_gates.py` is for -- **20 gates in about
40 seconds**, non-zero if any fails, `--all` to add the slow optimiser suite
(9 minutes, and it IS the wall time), `--only two_axis` to narrow and
`--jobs 1` when a failure needs reading in order. Gates run in PARALLEL because
they are independent processes: nothing here builds a page, and the only three
that write anything write into their own `tempfile` directories. **It DISCOVERS the gates rather than listing
them.** The list above named 9 and the repository has 20: `test_curve_flags`,
`test_curve_screen`, `test_kinetic_model`, `test_read_rre`,
`test_solution_chemistry` and `test_summary_kinetics` were in the tree and in
no documented suite. A hardcoded list is a list that drifts.

And each analysis folder's own `check_numbers.py`, which re-derives every number
in its `ANALYSIS.md` from the modules. About twenty seconds each.

`scope.frame` is MEMOISED and hands out a COPY (`scope._frame` is the cache).
Every row of it costs a progress fit and a debubble, so one build of the
two-axis block is about four seconds -- and the checkers ask for the same few
blocks over and over: `two_axis` called it 79 times over 2 distinct scopes,
`background_reaction` 123 times over 10. Unmemoised, those two runs were 290 s
and 465 s against 18 s and 25 s now. The copy is not optional: an lru_cache
handing out a shared DataFrame is one in-place edit from a silent wrong answer,
and this frame is passed into five folders. `_gas_curves` is the same pattern.

Units: concentrations mM, time s.

## One contract for the folder documents

`doc_check.py` is the comparison behind every `check_numbers.py`. There were
five copies of it until 2026-09-02 and no two were the same -- 7 substitutions
in `background_reaction`, 17 in `induction` -- so five documents were held to
five standards, and the folder with the most numbers ran the weakest one.

- **Typography folds; emphasis does not.** A hyphen against U+2212, `tau`
  against `τ`, `10-5` against `10⁻⁵`, backticks, a rewrapped line: all noise.
  But `**` survives, because in these documents bold marks the number a section
  argues for. Three folders stripped it while their claims BUILT it
  (`bold = "**" if ...`), so the emphasis was asserted and then discarded --
  live-looking assertions that could not fail. Build the markers and mean them.
  The contract is one-directional: a claim that carries `**` requires the
  document to carry it, a claim that omits it matches either way.
- **Read figure letters off the RENDERED page**, never off the builder source.
  `Checker.figures` does. Reading the source cannot see a letter used twice or
  out of order, which is how `induction/index.html` came to draw A B C D E G H
  F I and no check noticed.
- Add an assertion with `doc.claim` / `doc.check`; never write a private
  comparison in a folder. `test_doc_check` fails if a folder defines one.

**The duplicate guard covers all of it.** `test_curve_metrics.
test_no_duplicate_definitions` globbed `data/*.py` alone until 2026-09-02, and
the drift simply moved to where it could not see -- five copies of the document
contract, five palettes, and `_table` meaning a memoised frame in two folders
and an HTML renderer in a third. It now covers `data/`, the repository root and
every folder's `build_figures.py` and `check_numbers.py`: 58 modules. When it
fires, DELETE one copy and import the other.

**It reads CONSTANTS as well, since 2026-09-04.** It read only functions and
classes until then -- so no constant had ever been covered, which is the whole
category the rule exists for: the five palettes and the `TEMPERATURES` ramp are
all constants. 27 were defined in more than one of the 58 modules and five had
genuinely diverged. `temperature_series/build_figures.py` imported five names
from `figure_kit` and redeclared all five underneath, shadowing the import --
the 2026-09-02 palette drift, still live in one folder two days later.
`BURST_COLOUR` was `#12856a` in `curve_dossier` and `#7a4bb8` in
`background_reaction/build_figures.py`, which was also a colour declared in a
folder; `INDUCTION_FLOOR` was `60.0` SECONDS in `induction` and `1/300` of a
RUN in `slowdown`; `PALETTE` was three copies of which `svgplot`'s had eight
colours and the dossiers' nine.

Three shapes are skipped structurally rather than by name: a leading `_`
(private to its module), `X = X` (a re-export of an imported name, not a copy),
and lowercase module state. `PERMITTED_DUPLICATE_CONSTANTS` holds the rest --
each a module's own LOCATION or HARNESS STATE (`HERE`, `FAILURES`,
`DATASET_PATH`, `DOCUMENT`), never a measurement, a threshold or a colour.
**Where two constants are genuinely different things, the new name says which
one it is** -- `INDUCTION_GRID_FLOOR`, `DOSSIER_MINIMUM_POINTS`,
`DILUTION_TOLERANCE`, `BUFFER_FILENAME_MOLARITY`, `DOSSIER_PALETTE`,
`DOSSIER_BURST_COLOUR`, `REVIEW_STYLE`/`CURVE_DOSSIER_STYLE` -- and where one
was simply a copy, it is deleted and imported: `build_dossier` now takes
`ABSORBANCE_QUANTUM` and `QUANTISATION_SIGMA` from `curve_metrics` instead of
recomputing the noise floor, and `test_fit_kinetics` takes its ladder bar from
`scope.LADDER_MINIMUM` instead of restating it. Rename only when the two are
genuinely different things that shared a name, and then the new name has to say
which one it is.

`figure_kit.py` is the same idea for the drawing: `svgplot` has the primitives,
`figure_kit` the palettes, `fig`, `styled` and `write_pages`. The copies had
drifted here too -- `induction/` carried a different six-step `TEMPERATURES`
ramp from `temperature_series/`, and never drew with it, and
`background_reaction/` was the one folder whose figures did NOT sit on the fixed
light `SURFACE` its palette was validated against.

- **Ordered variables get sequential ramps** (`RUNGS`, `TEMPERATURES`,
  `PH_RAMP`); only genuinely categorical sets get `CATEGORY`, which is the
  validated trio. Never declare a colour in a folder.
- **Check a build's exit status.** `write_pages` returns non-zero on a clipped
  mark and the process exits non-zero on a crash. On 2026-09-03 a builder was
  broken for a commit because it was rebuilt with its output piped away, and
  the "page is byte-identical" check then passed for the wrong reason: the
  crash meant the file was never written at all.
- **A clip id names the REGION, not the object.** SVG ids are DOCUMENT-scoped,
  so `clip-path='url(#c)'` resolves to the FIRST clipPath of that name on the
  page -- not the one emitted beside the marks. `svgplot` named them
  `f"clip{id(self)}"` until 2026-09-04, and CPython reuses an address the
  moment a panel's `Axes` is freed, so `two_axis/progress_curves.html` drew 119
  clipPaths under 5 ids and FIVE OF THE SIX index pages bound one id to
  conflicting rectangles -- clipping a figure's marks to another figure's plot
  area and deleting whatever fell outside, silently, on the published page.
  `clipped_marks` cannot see this: it reads each figure's own first clipPath,
  which is the one the browser ignores. The id is now a hash of the rectangle,
  `svgplot.colliding_clips` is the check, and `write_pages` runs it. It also
  makes the build DETERMINISTIC -- with an address in the markup, the
  "page is byte-identical" check could never have passed.
- **`write_pages` reports clipped marks at build time**, which two of the five
  folders did and three did not. A mark outside the axis limits vanishes
  silently and the figure still looks finished, so the build has to say so.
- **Every folder has a `progress_curves.html`**, and it is the audit surface:
  `index.html` presents the argument, the curves page shows the fits the
  argument is read off. Build one with `progress_axes`, `progress_overlay` and
  `panel`; `progress_overlay` draws WHICHEVER FORM THE CURVE EARNED and asserts
  the fit is narrower than a mark. Draw the window a statistic was read
  through, and take it from the fit -- `sink_fit.tail_start` exists because a
  page that guessed the tail would draw a different window from the one that
  was fitted.
- **A curves page shows the whole block, including the control.** `induction/`
  draws both channels because its first claim is a contrast, and a page with
  only the catalysed half would show only the half that agrees. Each folder's
  `check_numbers.py` asserts one panel per live curve, so a page cannot lose
  one quietly.

## The two-axis block

`two_axis/` is exps 135-151. Two rules come out of building it.

- **A regressor is identified only where the fit's own offsets cannot absorb
  it.** `scope.orders` with `within=True` carries one indicator per experiment,
  so an axis that is constant inside every run and differs only between them is
  collinear with them -- and `scope._moves` now tests inside runs rather than
  over the whole column. Before 2026-09-03 it tested the column, and the block's
  peroxide arm (where `[S]` is fixed per run and differs only between runs) came
  back with a substrate order of **-6.094 +/- 0.078** on 63 curves. Nothing
  published moved, because every quoted order is measured where both axes move
  within runs. Read the arms with `scope.arm_orders`, which is also the check
  that the joint fit's additivity is free: worst gap 1.4 sigma of eight.
- **The block is a pH ladder as well as an L**, and that is the design nobody
  had used. Exps 136-142 carry one set of seven compositions and exps 143-151
  another -- same substrate ladder, different peroxide levels -- so inside
  either set `scope.ph_order` measures the pH axis with one offset per
  **cuvette** -- the mirror of the per-experiment offset the concentration
  orders use, and necessary because only 11.5% of the block's log[HOO-]
  variance is within-run against 100.0% and 94.1% for the concentration axes.
  A ladder is runs sharing a composition **and** an enzyme loading: `[enz]`
  steps 0.014 to 0.069 mM between runs, exactly where pH does, and a cuvette
  offset cannot absorb it.

Three things follow that are worth not re-deriving.

- **The pooled pH order is +0.554 +/- 0.040 on `vmax` -- a HALF order in
  [HOO-]**, over the strong runs. Over all 17 it is +0.400 +/- 0.034 and the two
  ladders disagree at 5 sigma instead of 2: the weak runs FLATTEN the ladder,
  because exps 149-151 sit at its bottom and measure the cell's own wander.
  Run `scope.strong_runs` before quoting anything from this block.
- **The two ladders control each other for free.** They walk pH in opposite
  directions against the schedule (+0.75 and -0.74), so a monotone drift over
  the twelve days enters their slopes with opposite signs; both are positive.
  `scope.ph_schedule_control`, and `scope.run_dates` for the instrument's own
  `Date Collected`, which is the only record of when a run happened.
- **[HOO-] is not the whole story.** `scope.hoo_consistency` moves it two ways
  and the orders part at **3.4 sigma** on `vmax`, peroxide the stronger lever.
  Either `[H2O2]` contributes beyond its hydroperoxide content -- the perhydrate
  of `buffer/` 6 would, since `[buf]` is fixed here -- or the rate saturates.
  The block cannot separate them: the L has NO INTERIOR POINT, so no cuvette
  crosses `[S]` with `[H2O2]` and no run crosses `[buf]` with `[H2O2]`.

- **The chop is O2, and the curves may NOT be stitched.** Many curves rise,
  fall by more in one 60 s reading than the reaction moves in five, and resume
  -- absorbance that goes away was never product. It is O2 from the catalysed
  decomposition of the peroxide: `bubble_ladder` (0 of 7 curves below 5 mM
  detach, 5 of 5 above 80 mM, monotone across six bands),
  `bubble_turnover_control` (exps 136 and 137 sit at 73.4 mM with NONE, being
  the two weakest runs there -- but read its pH confound below before quoting
  it) and `bubble_synchrony` (17 coincidences over 357 cuvette pairs against
  16.0 expected, so it is not the instrument).
  **The bubble is in the SAMPLE beam and it SCATTERS**, so absorbance climbs
  while it grows and drops when it goes. The reference omits only the enzyme,
  so it holds the same peroxide and would bubble too if the gas came from
  peroxide standing in solution. `bubble_step_asymmetry` is the measurement:
  122 falls beyond 20 sigma against 23 rises, so the gas is in the cuvette that
  holds the enzyme and EVERY RUN IS ITS OWN +/- ENZYME CONTROL. Do not use
  `bubble_turnover_control` for this -- it is CONFOUNDED WITH pH (exps 136 and
  137 sit at pH 6.95 and 7.53, and `turnover_control_confound` prices their
  silence at 0.68 and 1.48 expected events from pH alone), and do not use the
  archive's enzyme-free runs either: `gas_enzyme_control` puts the whole set at
  about ONE expected event. `solution_chemistry.oxygen_budget` kills the "too small to
  see" objection -- the solution saturates on 1.5% of the peroxide at the top
  of the ladder and cannot saturate at all at the bottom.
  **THE GAS IS A REACTION, NOT ONLY AN ARTEFACT.** It is the catalysed
  disproportionation of the peroxide and it is `MECHANISM.md` **S4**, a fourth
  competing sink alongside S1-S3: an unproductive, catalyst-dependent drain on
  the oxidant, fastest where the productive chemistry is strongest, and one no
  rate constant in the seven steps accounts for. WHERE it sits in the cycle is
  open -- KP shedding O2, KP + H2O2, or KD reduced before it can oxidise --
  and the archive cannot separate them, because no run moves the productive and
  unproductive routes against each other. `COMPUTATIONAL.md` C10.
  **IT IS NOT A PROPERTY OF THE BLOCK, AND pH IS THE TRIGGER RATHER THAN
  PEROXIDE.** `scope.gas_curves` runs the same test over all 402 curves of 88
  experiments, keeping the ones that did NOT bubble because those are the
  control. It appears with BOTH substrates -- 20 of 58 4OMe against 22 of 68
  BnOH, catalysed, above 40 mM and pH 8, in three buffers
  (`gas_substrate_control`), which is the prediction S4 makes since a catalyst
  decomposing peroxide involves no alcohol. And the archive's median [H2O2] is
  82.5 mM with 278 of 402 curves above 80, of which only 28 chop: inside every
  buffer the rate climbs with pH from a hard floor of ZERO detachments in 270
  hours below pH 7.5, over 23 experiments (`gas_survey`). Use `archive()` or
  `parse_scope("archive")` for the widest scope; nothing else in the project
  needs it.
  **THE KETONE IS NOT ESTABLISHED AS THE CATALYST, AND THE GAS IS NOT
  ESTABLISHED AS O2.** The beam asymmetry says the gas forms where the enzyme
  is; it does not say what in that cuvette does it, and a cyclodextrin or a
  trace transition metal would give the same peroxide order and the same pH
  rise. Nothing in this project has ever measured the gas -- no headspace, no
  manometry, no electrode. CO2 is close to excluded on a pKa argument (pKa1
  6.35, so above pH 8 it is bicarbonate and stays in solution, and the gas gets
  MORE common as pH rises, which is backwards for carbonate). Say "a
  non-condensable gas, made in the enzyme-containing cuvette, first order in
  peroxide and rising steeply with pH" and quote O2 as the reading, not the
  finding.
  **One curve is more than one bubble**: `bubble_record(141, 3)` sheds its
  LARGEST drop after its SHORTEST growth window (r = -0.91) and ends below
  every earlier post-release level. So THE GAS CARRIES OVER: a detachment
  empties one bubble of several, and dating each bubble from the previous drop
  cannot work.
  **Adding each step back is the one repair that must not be used**: a bubble
  that costs delta when it leaves contributed delta of rise while it grew, so
  stitching keeps the artefact's whole upward half -- it puts exp 135 cuvette 4
  at 1.26x the absorbance its own substrate could make, and against planted
  sawtooths it never beats DOING NOTHING (1.15/1.32/1.64/2.34 against
  1.12/1.24/1.58/2.26). Stitching IS `debubble` with its rate set to zero.
  **`curve_metrics.debubble` splits the readings into `f + b`** -- a
  non-decreasing chemistry and a non-negative gas made at a steady rate. THREE
  clauses: the gas may not outrun the curve it rides on (that IS `f' >= 0`), it
  must pay for every detachment out of gas already made, and it may hold no
  more than `unreleased_gas` -- the total of the detachments STILL TO COME, so
  a reading after the last fall carries no gas at all. One parameter,
  `bubble_rate`, pinned to the least rate that pays.
  **ONLY GAS THAT WAS WATCHED TO LEAVE IS SUBTRACTED**, and that third clause
  is where the recovery table splits in two. On a planting whose run stops
  making gas at its last release the repair is exact -- 0.99/0.98/1.00/0.99 and
  1.01/1.02/1.02/1.07 -- and on one still making gas at the last reading it
  KEEPS the bubble that never detached: 1.00/1.00/1.15/1.65 and
  1.03/1.04/1.17/1.79. `scope.bubble_recovery(ends_holding=...)` is both, the
  gap between them IS the systematic, and `curve_metrics.quiet_tail` says which
  curves could be paying it. On a curve with no detachment `debubble` returns
  the readings UNCHANGED.
  **THE BUBBLE THAT NEVER LEFT IS BRACKETED, NOT EXCLUDED.** That price is now
  per curve. `curve_metrics.terminal_gas` bounds what the beam may still hold
  at the last reading -- the fitted rate over the quiet tail, capped by what the
  tail rose -- and `scope.terminal_bubbles` is the table: 38 of 110 live curves
  carry one, 13 above a fifth of their rise. The bound CANNOT tell a run that
  ended mid-bubble from one that stopped making gas, because it asks the rate
  and not the readings; `curve_metrics.tail_excess` can, and it is the tail's
  slope minus the body's. Exp 140.4 runs +4.4e-05 AU/s faster, 0.83 of its own
  gas rate; exp 149.4 runs SLOWER, and 7 of the 8 curves past two shedding
  intervals do. It is ONE-SIDED -- an accelerating curve ends steeper with no
  gas in it and a decelerating one hides a bubble, both planted -- so a positive
  excess is evidence, a negative one is not. `vmax_terminal` is the far end of
  the bracket and NOTHING published lives inside it: +0.016 in substrate and
  -0.005 in peroxide against errors of 0.047 and 0.071. The two curves it bites
  hardest on, 140.4 and 142.4, are already `bubble_load` > 1.
  **A MONOTONE RECONSTRUCTION CAN STILL BE THE WRONG ONE** -- pulling a curve
  down by a smooth ramp leaves it smooth. Three faults hid behind that, and all
  three were found by eye, not by a test.
  Carrying the rate across a quiet tail invents gas nothing saw leave: exp 149
  cuvette 4 sheds 0.0031 AU once and ended holding 0.0273, rebuilding to
  -0.0209 against a raw rise of +0.0064. A ceiling at the most the beam had
  carried BOUNDED that without curing it -- exp 149 cuvette 3 still sat a flat
  0.0022 AU under its own readings for 82% of the run, and exp 150 cuvette 1
  under its by 99% of everything it rose. THE READINGS REFUTE THE
  EXTRAPOLATION: on 18 of 44 repairable curves the rate would have made more
  gas over the tail than the trace rose in total, and 11 of 44 ran more than a
  full shedding interval past their last detachment without one. Check
  `rebuild_smoothness`'s `gas_at_end` -- it is zero on every curve, so every
  reconstruction lands back ON the readings.
  And **A FALL THAT COMES STRAIGHT BACK IS NOT GAS**: gas that leaves the beam
  does not return, and a bubble cannot grow half its size in one 60 s reading,
  so `detachments` rejects a fall that a single adjacent reading undoes by more
  than `BUBBLE_RECOVERY_FRACTION`. Exp 149 cuvette 5's two "detachments" are
  9.3 and 8.2 sigma and neither is gas -- the first falls 0.00206 and the next
  reading climbs 0.00222 back -- and they were licensing the removal of 0.0097
  AU from a curve that rose 0.0262. 34 of 214 candidate falls are rejected;
  five curves lose all of theirs and are returned untouched.
  **`local_outlier_z` CANNOT be used for this**: its window spans the fall, so
  a genuine step flags itself (exp 135 cuvette 2's 0.1196 AU detachment scores
  +130). The test looks only at the two readings either side.
  So `rebuild_smoothness`'s guarantee is `worst_at_event`, NOT `rebuilt_worst`
  -- rejected excursions stay in the curve on purpose, and `isolated_outliers`
  is what nominates those. Quote the
  gap to `monotone_bound` as its systematic. Rewritten 2026-09-03: the segment
  ramp it replaced subtracted a ramp steeper than the curve rises (five steps
  past 8 sigma in 141.3) and skipped the second of two adjacent falls entirely
  (-0.0165 at 60 sigma in 144.2).
  **`rebuild_smoothness` is the test a repair has to pass**: the rebuilt curves'
  worst fall is -9.6 sigma against -260.4 as read and -5.8 on the curves that
  never bubbled. The ONE survivor is exp 135 cuvette 6, whose fall is in the
  first interval -- no rate explains a bubble grown before the run, and
  `debubble` returns that curve untouched.
  **`gas_rate_drivers`**: the fitted rate is +1.389 +/- 0.251 in peroxide and
  -0.307 +/- 0.103 in substrate, from a fit that never saw a concentration.
  **Read `bubble_load` before quoting a rate**: 13 of 110 live curves sit above
  1 and carry no measurable rate -- all four substrate rungs of exp 135, plus
  inner rungs of 138, 140, 141, 142 and 150. They are FLAGGED, NOT EXCLUDED.
  The SUBSTRATE order moves under no repair, but the PEROXIDE order does --
  +0.794 to +0.688 over all live and +0.871 to +0.760 over the strong runs,
  1.0 and 1.2 sigma -- which is what an artefact made from peroxide requires.
  Do not repeat the older claim that no order moves.
- **The early rise is counted by the CATALYST, not the substrate.**
  `curve_metrics.burst_amplitude` reads it off the FITTED CURVE, because the
  two-phase solve trades amplitude between its exponentials without moving the
  curve -- exp 135 cuvette 3 returns `B_fast` = -241 against `B_slow` = +303 on
  a curve that moves 0.06 AU, so neither `B_fast` nor their sum is the burst.
  `burst_bounded` says whether the rise finished inside the run: only 17 of 110
  live curves do, and the bound bites BETWEEN runs, not within them (inside a
  run the length is constant to one 60 s reading, between them it spans 9.6x).
  `scope.turnovers` divides by the catalyst; below pH 9 the median is 1.41 and
  above it 3.31.
- **`[enz]` is never identified within a run**, so `burst_drivers` gives its
  orders in `[S]` and `[H2O2]` (+0.136 +/- 0.045 and +0.793 +/- 0.063 over the
  strong runs) and refuses the catalyst. The block's one lever is
  `scope.enzyme_pair` -- exps 140 and 141, matched composition, 0.034 against
  0.014 mM, 0.07 pH units apart -- read at the largest window the two runs share
  and corrected for that pH gap with the block's own pH order. The rise is
  **proportional to the catalyst**: no dependence excluded at 3.6 sigma, first
  order admitted at 1.0 sigma. It rests on TWO RUNS, and
  `enzyme_pair_sensitivity` moves the apparent order from +0.77 to +1.35 across
  windows -- quote the range, not one window.

And what it cannot do: no enzyme-free curve in the whole pyrophosphate cell; no
windowed statistic across it (run length spans 9.6x); no LANDMARK induction
analysis (`signal_control` fails at +0.619 +/- 0.228) -- but that closes a
statistic and not the block, and `joint_clocks` asks the +1 rule through the
progress fit's own time constants instead. **No mechanism fit has ever been
run on it** -- `data/fits/` holds one, on BnOH/25 C/*phosphate*, sharing no
experiment with the block.

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

The sink's own activation energy carries a **window systematic larger than its
error** — 72 kJ/mol at the published window and 102 at 0.30, because a wider
window smooths the slow cold curves more than the fast warm ones. Quote it as
"72 with a systematic of about +30", never as 72.3 ± 10.0.
`slowdown.sink_window_sensitivity` is the sweep, and the null it feeds is
unaffected.

`data/slowdown.py` has the machinery. Use `deceleration_drivers` before
asserting that anything in this archive slowed down "over time" — one progress
curve cannot tell time from product, because inside a curve the product only
grows with time, and the separation exists only across curves.

## What the catalyst does first

`induction/` holds the other end of the same curves: every catalysed 4OMe run
begins slowly, and the induction is **the catalyst becoming active on its own
clock**, not the product seeding anything. It needs the catalyst (0 of 49
enzyme-free curves have one), it has no substrate order, and its barrier is
95 kJ/mol — four times too large to be physical.

`data/induction.py` has the machinery, and two rules come with it.

- **The induction landmark is safe within a run and not between runs.** Its
  rolling window is a tenth of the run, so a between-run comparison compares
  windows: at 25 °C the induction time regresses on run length at
  +0.437 ± 0.181 and on pH at −0.004 ± 0.123. Every concentration order is
  measured with one offset per experiment, and the temperature dependence comes
  from `arrhenius`'s fitted `inverse_tau`, which is not windowed.
- **A peroxide adduct constrains both orders at once.** `K + H2O2` in
  pre-equilibrium fixes `d ln v/d ln h - d ln tau/d ln h = 1` identically, for
  every K and every h, so the two questions are one and `joint_order` asks it as
  one regression -- `joint_peroxide_order` and `joint_buffer_order` are that
  function with the species filled in. Through the LANDMARK both blocks that
  can test it fall short (2.6σ and 3.7σ), and the rate is not first order in
  peroxide either -- `peroxide_saturation` rejects a = 1 at F = 32 on the
  two-axis ladder. Do not assume "first order in H2O2": that is the
  *unsaturated* limit of the scheme, not a consequence of it.
- **The +1 does not belong to the landmark, and the clock decides who can be
  asked.** `t_ind` is a rolling window a tenth of the run wide; `tau` and
  `tau_slow` come from the progress fit and carry none, so a block that fails
  `signal_control` or spans run lengths can still be asked through them.
  `joint_clocks` runs every clock on an axis BESIDE ITS CONTROL AXIS, and the
  control is the point: the +1 belongs to the activating species, so the
  substrate axis must MISS it, and it does by 4.2σ to 8.5σ. On the two-axis
  block the two routes disagree -- the peroxide axis falls 3.7σ short through
  the landmark, 1.9σ and 1.4σ through the fitted clocks. **Conclude nothing
  from that yet**: `tau_slow` is resolved on 32 of 110 live curves and the
  estimate moves +0.67 to +0.85 across cuts. Pass `gate=` and not `floor=` for
  a fitted clock -- a floor puts an unresolved constant ON the floor and calls
  it the fastest curve in the block.
- **CORRECT THE CLOCK, NOT JUST THE RATE.** `vmax_corrected` sat beside `vmax`
  from the start while `tau` and `tau_slow` were fitted to the readings alone,
  so until 2026-09-04 every question asked of a time constant was asked of a
  curve with the O2 still in it -- and 40% of the curves with a resolved
  `tau_slow` carry detachments. On a peroxide axis that is not a wash: the gas
  is MADE from peroxide, so it inflates the rate's order AND shortens the
  apparent clock, both pushing `d ln v - d ln tau` towards the +1 under test.
  Asked of the readings the block's `tau_slow` row sat 0.3σ from +1; asked of
  the rebuilt curves, 1.4σ. `frame` now carries `tau_corrected`,
  `tau_slow_corrected` and their resolved flags, `joint_clocks` DEFAULTS to
  them, and `JOINT_CLOCKS_RAW` is kept so the difference can be shown. The
  repair costs no resolution -- it buys some (62 to 65 and 25 to 32 curves),
  because the artefact was what those fits could not pin.
- **That `+1` is not about H2O2.** It holds for ANY species held in excess that
  draws the catalyst into its active form, so it transfers to any axis the
  archive moves. On the BUFFER axis it is MET -- `joint_buffer_order` gives
  +1.094 +/- 0.150 at the 450 s window -- which is the strongest thing this
  archive says about what E -> E* runs on. Read only the windows whose
  `signal_control` passes: 300 s and 450 s do, 600-1200 s do not, and the ones
  that fail overshoot +1 in exactly the direction the artefact predicts.
- **`scope.frame` carries the SIGN of the early curve** -- `progress_kind`,
  `B_fast`, `B_slow`, from the form the curve earned. `depth` is floored at
  zero and cannot see a curve that begins fast; the two-axis block splits
  almost evenly (46 lag-first against 45 burst-first of 110 live), so an
  induction time averaged over it means little.
- **Every buffer order in this project is an order in TOTAL buffer.** At one
  pH the acid, the base and the total are proportional, so a titration cannot
  name the species. `buffer/` has the two-pH test the archive does hold and it
  excludes nothing: +1.06 +/- 0.77 against 1.76 for general base and 0.52 for
  general acid. There are FIVE buffer titrations (exps 32, 34, 35, 36, 37), not
  the two that were being used.
- **And the archive cannot say whether the buffer acts on the catalyst or on
  the peroxide.** A general base is a term in `[buf]`; a buffer perhydrate
  delivering the oxygen is a term in `[buf][H2O2]`. They differ by that
  interaction and by nothing else at one pH, and `peroxide_crossing` is the
  answer: of 88 runs, 53 step `[buf]`, 20 step `[H2O2]`, and **0 step both**.
  The species test does not separate them either -- nucleophilic attack by the
  dianion has general base's pH signature -- so the perhydrate scheme is a
  FOURTH survivor of `buffer/` 3, not a resolution of it. `buffer/` 6 and
  `COMPUTATIONAL.md` C9.
- **Do not attribute the pH rise of the buffer-free route to [HOO-] alone.**
  `free_route_order` divides one by the other on the one matched-substrate pair
  the archive has: the level rises 7.90x over 0.50 pH units where [HOO-] gives
  3.23x, an apparent order of +1.79 rather than +1.02. It is two runs on two
  days against a measured between-day step of 1.80x, so it is a pointer -- but
  the earlier text asserted the two agreed, and they do not.
- **`[S]` and `[buf]` move together in every 4OMe run**, at -0.96 in logs,
  because substrate volume displaced buffer volume. A substrate order measured
  there is an order in the pair, and `induction.composition_collinearity` is
  the number. Exps 135-151 are the block where `[buf]` is constant and `[S]`
  moves alone -- and the two blocks give OPPOSITE substrate effects on the
  sign of the induction, which is what makes the buffer a live candidate for
  the E -> E* step.
- **Never compare a time constant across two runs of different length.** The
  archive's one buffer titration (exps 32 and 34) earns the two-phase form in
  one run and the one-phase form in the other, because the runs are 5280 s and
  1767 s, so `tau_fast` is tau1 of one fit and tau of another. Read that way it
  looked like the two runs disagreed in sign; read through a landmark with a
  window in SECONDS common to both, they agree at -0.433 +- 0.201 -- more
  buffer, shorter induction. `induction.buffer_lever` and `buffer_order`.
- **Run `signal_control` before believing any induction result.** The landmark
  is the first crossing of half the largest rolling slope, so on a curve with no
  signal it measures the spectrophotometer. The catalysed 4OMe block passes
  (+0.003 ± 0.149); both blocks that carry a peroxide ladder fail, which is why
  this archive cannot say whether the induction is the catalyst binding
  H<sub>2</sub>O<sub>2</sub>.
