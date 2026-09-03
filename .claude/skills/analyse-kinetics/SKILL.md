---
name: analyse-kinetics
description: Use when analysing, measuring, plotting or answering any question about the kinetics experiments in this thesis repo — any block, in either substrate — including rates, noise, lag, induction, autocatalysis, substrate/peroxide/buffer/pH order, curve quality, or "how many curves do X". Load before writing any analysis script against data/.
---

# Analysing the kinetics data

**There is no privileged block.** Each question names the block it needs, and
every block is defined in `data/scope.py`: `TWO_AXIS_BLOCK`,
`TEMPERATURE_SERIES`, `FREE_BNOH_ALL`, `FREE_BNOH_PHOSPHATE`, `BUFFER_FIXED`,
plus `buffer_role.TITRATIONS`, `slowdown.substrate_blocks` and
`induction.induction_blocks`. Add to that vocabulary rather than selecting
experiment numbers inline.

The **two-axis block** (`fit_dataset.TWO_AXIS_BLOCK`, exps 135-151, 119 curves
of BnOH / 25 °C / pyrophosphate) is the block the mechanism fitting uses and
the one `frame()` defaults to. It is named for its design — the archive's only
runs that move both concentration axes inside a single run — not for its
chemistry, which would not select it: exps 75 and 76 are also BnOH, catalysed,
pyrophosphate at 25 °C and are not in it. It was called `PRIMARY_SCOPE` until
2026-09-03; see `DATA_VERIFICATION.md` that date, and `FITTING.md` for the case.

Most of this thesis's analysis happens **outside** that block — the temperature
series is the only route to activation parameters, the deceleration needs 84
curves across the 4OMe archive, the buffer question needs five titrations. The
rules below apply to all of it.

## The rule

**Import from `data/scope.py` and `data/curve_metrics.py`. Do not re-derive.**

Every measurement that describes a curve — its noise, its initial rate, where
its slope peaks — has exactly one definition, in `curve_metrics`. Every block
selection has exactly one definition, in `scope`.

This is not a style preference. Six functions in this repo were once defined
twice and four had diverged; the lag statistic's two copies disagreed on 96 of
402 curves and one of them would have put 21% in a thesis where the evidence
said 34%. `data/test_curve_metrics.py::test_no_duplicate_definitions` now fails
if a duplicate reappears.

That archive figure is now **37.6%** (151/402), not because the statistic
changed but because the readings did: since 2026-08-31 they come from the
instrument's own `.rre` files rather than the 0.001 AU `.txt` exports. All 119
curves of the two-axis block are `.rre`; only exps 2–32 are still on the export, and no
`.rre` survives for them.

Use `curve.noise`, never a fresh `curve_noise(values)` — the floor depends on
`curve.source` and `frame()` already applies the right one. **If you call
anything in `curve_metrics` directly, pass `floor=source_floor(curve.source)`**
(from `fit_dataset`). The default is the `.txt` export's quantisation, which is
1096× too coarse for a `.rre` curve, and it reaches every standard error, not
just the noise: until 2026-09-01 `line_fit` hardcoded it and the
two-axis block's acceleration count read 48/110 where the instrument's own readings say
**51/110**.

The same failure has already happened in prose: the within-experiment contrast
figures were computed once in a throwaway script over the wrong curve set and
written into four documents as 98.4% / 82.4%, when over the actual block they
are 100.0% / 94.1%.

## Start here

```python
import sys; sys.path.insert(0, "data")
from scope import frame, design, curves, ladder, within_experiment_share

df = frame()      # the two-axis block; pass a block for any other
```

`frame()` columns: `experiment sample pH s0 h2o2 e0 hoo duration_s points
noise net live v0 v0_stderr v0_rms peak lags late_over_early`.

Most questions are a groupby on that frame, not a new script:

```python
df[df.live].groupby("experiment").v0.median()        # rate by run
df[(df.pH >= 9) & df.live].lags.mean()               # lag fraction at high pH
design()                                             # per-experiment design table
within_experiment_share("s0")                        # the number the block rests on
```

From the shell: `python data/scope.py` and `python data/scope.py --design`.

### The two-axis block has a second design, and its own accessors

Exps 136-142 carry one set of seven compositions and exps 143-151 another, so
inside either set a cuvette is matched across runs and the block is a pH ladder
as well as an L. `two_axis/` is the folder; these are the functions, all in
`scope`:

```python
arm_orders()              # each order from the arm that moves only its axis
ph_ladders(strong_runs()) # runs sharing a composition AND an enzyme loading
ph_order("vmax", scope=strong_runs())   # one offset per CUVETTE, not per run
ph_schedule_control()     # the two ladders walk pH in opposite directions
hoo_consistency()         # move [HOO-] two ways; they disagree at 3.4 sigma
run_dates()               # the instrument's own Date Collected
acceleration_by_ph()      # the acceleration is high-pH, not long-run
burst_table()             # the early rise, and it in units of catalyst
burst_drivers()           # its orders in [S] and [H2O2], within runs
enzyme_pair()             # the block's one lever on [enz]: exps 140 vs 141
enzyme_pair_sensitivity() # and the window it is read at, swept
bubble_table()            # every curve's O2 load, and its rate under each repair
bubble_ladder()           # detachments against [H2O2]: the evidence
bubble_turnover_control() # and that peroxide alone does not do it
bubble_synchrony()        # not the lamp: 17 coincidences against 16.0 expected
bubble_mass_balance()     # what stitching would claim of the substrate
bubble_recovery()         # recovered vmax / true vmax, at a planted truth
rebuild_smoothness()      # does a repaired curve look like a clean one?
gas_rate_drivers()        # what the fitted gas rate depends on
bubble_sensitivity()      # every order under every repair
bubble_step_asymmetry()   # which beam: 122 falls beyond 20 sigma against 23 rises
bubble_record(141, 3)     # one curve's detachments, one row each
```

The bubble is in the **sample** beam and it **scatters** light out of the
aperture, so absorbance climbs while it grows and drops when it detaches. The
reference omits only the enzyme, so it holds the same peroxide and would bubble
too if the gas came from peroxide standing in solution; it needs turnover,
which only the sample has. `solution_chemistry.oxygen_budget` answers "could a
side reaction really make a bubble" -- the solution saturates on 1.5% of the
peroxide at the top of the ladder and cannot saturate at all at the bottom.

**The chop is gas, and the curves may NOT be stitched.** Many of the block's
curves rise, fall by more in one 60 s reading than the reaction moves in five,
and resume -- absorbance that goes away was never product. Adding each step
back is the repair that suggests itself and the one that must not be used: a
bubble that costs delta when it leaves contributed delta of rise while it grew,
so stitching keeps the artefact's whole upward half. It puts exp 135 cuvette 4
at 1.26x the absorbance its own substrate could make, and against planted
sawtooths it never beats DOING NOTHING at any severity -- stitching IS
`debubble` with its one parameter set to zero.

Use `curve_metrics.debubble`. It splits the readings into `f + b`: a
non-decreasing chemistry and a non-negative gas made at a steady rate, which
may never outrun the curve it rides on, must pay for each detachment out of gas
already made, and after the LAST detachment may not exceed `bubble_ceiling` --
the most the beam carried while detachments were still happening. That last
clause is not decoration: without it a rate fixed by early drops is
extrapolated across hours with none, and exp 149 cuvette 4 rebuilt to -0.0209
against a raw rise of +0.0064. A MONOTONE RECONSTRUCTION CAN STILL BE THE WRONG
ONE -- read `rebuild_smoothness`'s `gas_held` against `biggest_bubble`. `bubble_rate` pins the one parameter to the least rate that
pays, so the result is an UPPER bound on the chemistry -- quote the gap to
`monotone_bound` as its systematic. It recovers a planted `vmax` to within a
tenth at every severity up to 2x whether the bubbles empty or only partly
empty, and on a curve with NO detachment it returns the readings unchanged.

`rebuild_smoothness` is the check: rebuilt curves fall by at worst -9.6 sigma
against -260.4 as read and -5.8 on curves that never bubbled. The one survivor,
exp 135 cuvette 6, has its fall in the FIRST interval -- a bubble grown before
the run leaves no rise to date it from, and that curve is returned untouched.
`gas_rate_drivers` is the independent check on the diagnosis: the fitted rate
is +1.203 +/- 0.221 in peroxide, from a fit that never saw a concentration.

**Read `bubble_load` before quoting a rate off this block.** Thirteen of 110
live curves sit above 1 and carry no measurable rate -- all four substrate rungs
of exp 135, plus inner rungs of 138, 140, 141, 142 and 150. They are FLAGGED,
NOT EXCLUDED: they stay in the frame, the live counts and the curves page.
Nothing published moves under any repair (`bubble_sensitivity`).

`burst` is `curve_metrics.burst_amplitude` on the FITTED curve, never on
`B_fast` -- the two-phase solve trades amplitude between its exponentials
without moving the curve. Check `burst_bounded` before comparing two runs:
only 17 of 110 live curves finish their rise, and the rest are lower bounds.

Two cautions come with them. **Filter to `strong_runs()` before quoting a pH
order** -- the weak runs sit at the bottom of the ladder and flatten it, from
+0.554 +/- 0.040 to +0.400 +/- 0.034. And **a regressor is identified only
where the fit's offsets cannot absorb it**: `orders(within=True)` carries one
indicator per experiment, so an axis constant inside every run returns NaN
rather than a number (`scope._moves`, fixed 2026-09-03).

## If the quantity you need is missing

Add it to `scope.py` (a selection or a derived column) or `curve_metrics.py`
(a measurement of curve shape). **Do not compute it inline in a script**, even
a throwaway one — throwaway numbers end up in documents. If it was worth
computing once it will be asked for again.

## The +/- chemzyme controls

The two-axis block has no enzyme-free run — every one of exps 135-151
carries enzyme.
The nearest thing is `scope.PAIRED_CONTROLS` (exps 65/66, 67/68, 69+70/71):
the same ladder run with and without 0.028 mM chemzyme, in phosphate and boric
buffer at pH 8.0-8.5. `python data/scope.py --controls`.

Use them for interpretation, never for a fit or a pooled constant -- different
buffer, and one `[H2O2]` per run. They show no rate enhancement (0.63x over 9
live rungs). **Do not quote that as evidence the chemzyme is inactive**: at
0.028 mM the literature's kcat predicts only a 1.3x enhancement over these
runs' own background, which is smaller than the 1.55x by which exps 69 and 70
-- the same experiment twice -- disagree. `python data/scope.py --literature`.

For a background rate at neutral pH use `scope.FREE_BNOH_NEUTRAL` (exps 3, 6).
They are the buffer titrations, so they give a RATE but never an ORDER. And when reading a substrate order from
enzyme-free runs, use `scope.FREE_BNOH` (exps 65, 67, 69, 70) and **not**
exps 3 and 6, which are buffer titrations whose rate falls with substrate for
a reason that is not substrate. `scope.FREE_BNOH_BUFFER_TITRATIONS` records it.

## Conventions that are already settled

- **Sheet over filename.** A declared sheet value beats a filename; filenames
  get copied forward between runs. Inverting this caused the worst error in
  `DATA_VERIFICATION.md`.
- **A block is a constant in `scope.py`, never a folder.**
  `data/Mads/good data BnOH/` is *not* the two-axis block — it also holds exp 50
  (excluded), exp 51 (a different cell) and exp 134 (no instrument data).
- **Curve shape is never a defect.** Lags, dips and bursts are reported, not
  excluded.
- **There is no background run in this block.** Exps 150 and 151 were called
  the in-cell background; they are not. Every run in the block carries enzyme, and
  cuvette 1 of exps 150-151 sits within 2x of the [HOO⁻] trend the other runs
  define, so the reaction is still on — they are the bottom rung of the pH
  ladder, not a blank. Do not subtract them from anything.
- **Six runs move with something that is not the reaction.** In exps 136, 137,
  147, 149, 150 and 151, `scope.concentration_agreement` is 0.61 or below and
  as low as 0.005: their cuvettes' rates bear almost no relation to their
  cuvettes' concentrations, while exps 135, 138, 139, 140 and 142 run 0.93 to
  0.97. Their rates sit at the cell's own drift, a few times 1e-7 AU/s.
  Every conclusion in FITTING.md is *stronger* with them dropped, so use
  `concentration_agreement` before letting a weak run carry an argument.
- **Concentrations are mM, time is s**, throughout.
- Quote a rate constant only where it is identifiable; two of the six are lower
  bounds (`FITTING.md` F4).

## Before reporting a number

Numbers that reach a document must come from `scope`/`curve_metrics`, not from
a scratch calculation. If a number contradicts one already written down, say so
and check which selection each was computed over — that has been the cause both
times it has happened.
