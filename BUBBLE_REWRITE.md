# Rewriting the bubble machinery: segmentation instead of stitching

**Status: a plan, not a record.** Nothing here has been done. When it has,
this file is superseded by a `DATA_VERIFICATION.md` entry and should be
deleted rather than left to rot. Numbers quoted below were measured on
2026-09-10; the ones that survive the work belong in `BUBBLES.md`, gated,
not here.

---

## 1. What is actually wrong

Not tidiness. **The abstraction is wrong.**

The module detects events *pointwise* — a fall here, a rise there, each judged
by its own local veto — and then tries to reconstruct the physical picture by
stitching those points together afterwards. But the physics is **phases**: the
beam alternates between accumulating gas and releasing it, over runs of
readings. Nothing in the code represents a phase. `detachments` gets closest,
grouping strictly consecutive falls, and every failure found on 2026-09-10 is
a case where that grouping was not enough.

The symptom is accretion. Fixing the 2026-09-10 defects properly needs four
more rules on top of the existing ones:

| rule | what it patches |
|---|---|
| merge consecutive anomalous rises | a bubble growing over several readings is credited with one step of itself (exp 44.1: 0.0231 of a real 0.0667) |
| merge falls across a one-reading gap | a stuttering release fragments, and the uptick between fragments scores as an arrival (exps 43.1, 135.2) |
| merge rises across a one-reading gap | the mirror, and its absence splits exp 135.1's 274–277 in two |
| an excursion guard | …which then had to be **retracted**: it struck out exps 135.1 and 138.4, both real |

Four rules, three new constants, one already withdrawn — all doing the job
that one segmentation should do. That is the signature of a wrong
decomposition, not of missing care. The functions listed in §6 come to **810
lines**, and the reason they are hard to reason about is that the phase
structure is implicit in all of them and owned by none.

**And the vetoes are biased in opposite directions.** `_is_excursion` is
called on the curve for a fall and on the *negated* curve for a rise. Its
recovery test compares absolute absorbance while the chemistry climbs
underneath, so the reaction *adds* to apparent recovery for a fall and
*subtracts* for a rise. `_local_step_scale` compensates for one direction.
Nothing compensates the other. (`BUBBLES.md` §6.6 has the measurement, and why
the obvious fix was rejected.)

---

## 2. The design

**Segment once, classify segments, read events off the classification.**

### 2.1 One anomaly score

`step_anomaly(times, values, noise)` → a signed score per interval: how far
this step departs from what the curve does locally, in units of the local
step size (`_local_step_scale`'s estimator, which already exists). This is the
single place where "is this step unusual" is decided.

It replaces two thresholds that currently answer that question in different
currencies — `BUBBLE_DROP_SIGMA` in units of curve noise, `EXCURSION_LOCAL_SIGMA`
in units of local step. **Absolute σ is the wrong currency and exp 144.2
proves it**: every step of its real 14-reading acceleration clears 6σ (7–12σ),
and *not one* is anomalous against its own neighbourhood. That curve is the
whole reason arrivals are currently forbidden to merge; scored locally, it
needs no special rule at all.

### 2.2 Segments

`segments(times, values, noise)` → maximal runs of consistently-signed
anomalous intervals, each an `accumulate` or a `release`, tolerating an
interruption of up to `SEGMENT_GAP` readings that gives back less than
`SEGMENT_GAP_RESTORE` of the run so far.

One gap rule, applied symmetrically. It replaces the two merge rules above and
the strict-consecutive grouping in `detachments`.

### 2.3 Classification is a property of a segment PAIR

This is the part that earns the rewrite. A spike is not a single event that
"comes back" — it is **two adjacent segments that cancel**. Expressed that
way, three cases that currently need separate machinery fall out of one rule,
and the order of the pair carries the physics:

| pair | reading | verdict |
|---|---|---|
| `accumulate` → `release` | gas arrives, then leaves | both real — the ordinary bubble lifecycle |
| `release` → `accumulate`, level **returns** | the trace goes back to where it was | both artefact — an instrument spike |
| `release` → `accumulate`, level **overshoots and holds** | a bubble left and a new one grew | both real |

Row two is exp 149.5 — the curve that forced the recovery clause. Row three is
exp 138.4, which the current code gets wrong in both halves: it rejects the
real 2-reading detachment because a rise follows, then (with the guard) refuses
the real 4-reading arrival for following a rejected fall. One rule, correctly
scoped, decides both.

Gas must arrive before it can leave, so the asymmetry is physical rather than
fitted: `accumulate → release` is never a spike, whatever the sizes.

### 2.4 The public API does not move

`detachments`, `bubble_arrivals`, `split_arrivals`, `debubble`,
`bubble_profile`, `bubble_rate`, `apply_gains` all keep their signatures and
their meanings. `detachments` becomes "the release segments that classified as
real"; `bubble_arrivals` becomes "the accumulate segments that did". The
rewrite is **internal**, which is what bounds the blast radius: `scope.py`, the
eight folders and every existing test call the same names.

`_is_excursion` is deleted. So is the guard, which was never shipped.

---

## 3. Success is already defined

`data/bubble_cases.py` (2026-09-10) is the fixture, and it exists precisely so
this work can be judged instead of argued about. 25 rows, each a curve, an
event, a verdict and the evidence for it.

- **The 18 pinned rows must still pass.** These are the curves the current
  thresholds were pinned to. If segmentation cannot hold them, the design is
  wrong and the work stops — that is the abandon condition, and it is checkable
  on day one of phase 3.
- **The 7 open rows should flip.** They are asserted to *fail* today; `test_bubble_cases`
  reports when one starts agreeing. Five of seven is a good outcome; seven is
  the target.
- `worst_at_event` (−9.61) and `rebuilt_worst` (−61.1) must not degrade.
- Nothing published may move outside its own standard error.
- The constant count must go **down**.

---

## 4. Phases

Each phase ends green on `run_gates.py`, and nothing is deleted until the
thing replacing it is checked.

1. **Extract `step_anomaly`.** Pure function, nothing calls it yet. Assert it
   reproduces the existing sign/threshold decisions where they are currently
   made. No behaviour change; verify the archive rebuilds bit-identical.
2. **Build `segments()` beside the old code, unused.** Test it directly
   against `BUBBLE_CASES` — every pinned event must fall inside a segment of
   the right kind, before any classification exists.
3. **Reimplement `detachments` over segments.** Run the 18 pinned rows. This
   is the go/no-go gate: 135.1's two falls, 130.2's two, 149.5's four, and the
   four SNR-floor curves either hold or the design is refuted.
4. **Reimplement `bubble_arrivals` over segments.** Check the 7 open rows flip
   and that 144.2 still yields nothing in readings 29–43.
5. **Sweep the constants against the fixture.** Each new constant gets the
   `DETACHMENT_SNR_FLOOR` treatment: show where it sits in the distribution
   and whether that is a gap or a continuum. Record the sweep. Where it is a
   continuum — as `BUBBLES.md` §6.5 already found for the kink bar — say so in
   the constant's own comment rather than implying a calibration that does not
   exist.
6. **Cut over and clean up.** Delete `_is_excursion`, the strict-consecutive
   grouping and `arrival_candidates`' per-step scoring. Re-run every
   `check_numbers.py`, rebuild all eight folders' figures, update `BUBBLES.md`
   §4.4–§6.6 and the `CLAUDE.md`/skill passages, and write the
   `DATA_VERIFICATION.md` entry. Delete this file.

---

## 5. What will move, and what is exposed

Measured for the prototype on 2026-09-10; the rewrite should land near these
but will not match them exactly.

- Archive detachments 374 → ~356 (**same falls, fewer events** — grouping, not
  detection).
- Arrivals 89 → ~87, with 14 of them credited in full rather than by their
  largest step.
- Two-axis peroxide order +0.7042 → ~+0.712, gas rate +1.343 → ~+1.352, the
  `tau_slow` row +0.757/0.84σ → ~+0.752/0.86σ. All inside their own errors.

**The boric ladder is the exposed block.** Exps 43, 44, 45 and 49 carry 13 of
the 24 curves with arrivals, and `ph/ANALYSIS.md`'s turnover claim is read off
exps 43 (the peak) through 49. It survived the 2026-09-10 change and was
slightly strengthened; it must be re-checked here, not assumed.

**One case to watch.** Exp 139.2's 3720 s jump is currently unresolvable (
`BUBBLES.md` §6.5) and sits just under the kink bar. Segmentation may pair it
with the sub-threshold −4.1σ fall two readings later and admit it. That would
be a *decision the rewrite makes by accident*, so it needs a deliberate look
when phase 4 lands rather than being accepted silently.

---

## 6. The code as it stands

For scale, and so the cleanup can be measured:

| function | lines | fate |
|---|---|---|
| `local_outlier_z` | 40 | kept — the kink test still needs it |
| `_local_step_scale` | 23 | kept, promoted into `step_anomaly` |
| `bubble_drops` | 24 | absorbed into segmentation |
| `detachments` | 43 | reimplemented over segments |
| `_is_excursion` | 75 | **deleted** — replaced by the pair rule |
| `bubble_arrivals` + `arrival_candidates` | 167 | reimplemented over segments |
| `split_arrivals` | 27 | kept unchanged |
| `unreleased_gas`, `bubble_profile`, `bubble_shortfall`, `bubble_rate`, `bubble_onset`, `apply_gains`, `debubble` | 411 | **untouched** — the model layer is not what is broken |

The model layer is sound: `A_obs = f + b`, one fitted rate, the three bounding
clauses, and the 2026-09-10 arrival routing. Only the detection layer is being
replaced.
