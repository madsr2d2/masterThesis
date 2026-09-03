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
factor of 1096 between them. Never floor a noise **or a residual variance** at
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
python data/test_fit_kinetics.py          # selection, scope, parameter recovery
python data/test_validator.py             # fault injection
python data/test_slowdown.py              # the slowdown models and their regressions
python data/test_induction.py             # the induction landmark and its controls
python data/test_buffer_role.py           # the species test, planted both ways
python test_doc_check.py                  # the contract every check_numbers runs on
```

And each analysis folder's own `check_numbers.py`, which re-derives every number
in its `ANALYSIS.md` from the modules. About a minute each.

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
every folder's `build_figures.py` and `check_numbers.py`: 55 modules. When it
fires, DELETE one copy and import the other. Rename only when the two are
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
windowed statistic across it (run length spans 9.6x); no induction analysis
(`signal_control` fails at +0.619 +/- 0.228). **No mechanism fit has ever been
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
  every K and every h, so the two questions are one and
  `joint_peroxide_order` asks it as one regression. Both blocks that can test
  it fall short (2.6σ and 3.7σ), and the rate is not first order in peroxide
  either -- `peroxide_saturation` rejects a = 1 at F = 32 on the two-axis
  ladder. Do not assume "first order in H2O2": that is the *unsaturated* limit
  of the scheme, not a consequence of it.
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
