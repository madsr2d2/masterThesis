# C7 — what the catalyst does during the induction period

Handover plan, written 2026-09-16. It executes `COMPUTATIONAL.md` C7 under that
register's conventions. Read this file from the top.

> **The induction period is the one step in the cycle with a measured barrier
> and no assignment.** The candidate is that the active-site ketone sits mostly
> as its gem-diol hydrate, and that the catalyst wakes up by dehydrating.
> This plan tests that candidate where it is cheapest to kill, and stops.

## The goal, in one paragraph

Decide whether the hydrate candidate survives its own pre-registered gate, and
— if it does — whether the buffer is what dehydrates it. The second question is
the one the kinetics could not answer: `MECHANISM_EVIDENCE.md` §4 records that
candidate C0 (activation unimolecular, held back by base) and C1 (the buffer
drives activation) tie, and that the design cannot separate them. A
water-assisted barrier against a phosphate-assisted barrier separates them
computationally.

## 0. Rules for the executor

**Repository rules** (`CLAUDE.md`, `COMPUTATIONAL.md` Conventions):
- Commit straight to master. Subject `area: what changed`, one paragraph, then
  exactly these two lines:
  ```
  Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01Vht919cHDPUdJWghusZFdS
  ```
- **Track inputs and small `.xyz` only.** Outputs and scratch stay
  homelab-local, as `.gitignore` already enforces. `uvvis/` is exempt by a
  prior decision; this task is not.
- **Layout.** `computational/C7_hydrate_clock/<reaction>/<stage>/`, with
  reusable species geometries in `computational/species/`. Reaction names use
  `MECHANISM.md`'s shorthand (`K+H2O_to_Kh`), never an ad hoc English name.
- **Never edit** `MECHANISM.md`, `FITTING.md`, `CLAUDE.md`, `BUBBLES.md`,
  `MECHANISM_EVIDENCE.md` or `COMPUTATIONAL.md`'s existing entries. Propose
  text in a stop message. You WILL append to `COMPUTATIONAL.md`'s Log, and
  only there, and only what the templates in section 7 say.
- **A number that reaches a document comes from a job that is in the
  repository**, named by its path, never from a recollection.

**Rules that exist because they were broken before** (three of them in the
kinetics work, one in `computational/` itself):
1. **A gate is stated before the run, and never widened afterwards.** If a gate
   fails, that is the result. Report it and stop.
2. **Never trust a geometry you have not checked.** After every optimisation,
   run the connectivity check of section 3.6: the formula, the bonds that were
   supposed to form or break, and the frequency count. A geometry that
   optimised to something else is the single most common way a computed number
   becomes fiction — `computational/C8_perhydrate_trap/` holds a whole reaction
   folder renamed `UNIDENTIFIED` for exactly this.
3. **Never trust one initial Hessian across a long optimisation.**
   `Calc_Hess true` with `Recalc_Hess 50` for a minimisation and `10` for a TS
   search, per `COMPUTATIONAL.md`.
4. **Validate every input before running it.** Use the ORCA input-syntax
   validator available in this environment, and keep the smallest possible
   test job as the first run of any keyword combination this project has not
   used before (`computational/hellowater/` is the precedent).
5. **Report wall time and `nprocs` for every job**, in the log entry.

## 1. What this has to reproduce

From `induction/ANALYSIS.md`, `COMPUTATIONAL.md` C7 and
`MECHANISM_EVIDENCE.md`. None of it is assumed here; each is a number the
calculation is measured against.

| measured | value |
|---|---|
| the induction's barrier, window-free clock, six temperatures | **77 ± 12 kJ/mol** |
| the same through the one-phase τ, four temperatures | 95.0 ± 15.7 kJ/mol (ΔH‡ 92.6) |
| the induction's depth at 15 °C, a LOWER bound | 0.79 |
| the same at 40 °C | 0.06 |
| the clock against pH, four ladders | `d ln τ/d pH` = +0.12 to +0.34 |
| the clock against substrate, buffer fixed | `d ln τ/d ln[S]` = −0.18 to −0.71 |
| what the kinetics cannot decide | buffer-driven against unimolecular activation |

**The two signs matter as much as the barrier.** Base slowing the clock is what
an unreactive conjugate base of the hydrate predicts; substrate hurrying it is
what a hydrophobic guest displacing the hydration equilibrium predicts.

## 2. The chemistry, and the model

**The active site.** The catalyst model this project already uses
(`computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/rc/job.xyz`,
140 atoms) carries its ketone as a **1,3-dialkoxy ketone**: the carbonyl carbon
(index 128, 0-based, in that file) is flanked by two CH₂ groups, each bearing a
ring ether oxygen. It is not fluorinated. Two β-alkoxy groups withdraw
inductively, so a raised hydrate fraction is plausible — and a very large one is
not. **The gate can fail, and it is supposed to be able to.**

**The species** (`MECHANISM.md`'s shorthand):

| symbol | species |
|---|---|
| K | the active-site ketone, free carbonyl |
| Kh | its gem-diol hydrate |
| Kh⁻ | the gem-diolate, Kh less one hydroxyl proton |
| KP | the perhydrate, K + H₂O₂ |

**The truncation, and why it is legitimate.** Hydration is a local change at one
carbon. Tasks 1 and 2 use **1,3-dimethoxypropan-2-one**,
`CH₃-O-CH₂-C(=O)-CH₂-O-CH₃`, as the truncated active site: the same carbonyl
with the same two β-alkoxy substituents, capped with methyl where the macrocycle
continues. Full thermochemistry is affordable there and is not affordable on the
macrocycle. **Task 3 then asks what the cavity does to the answer**, and its
result is a correction to Task 2, not a replacement for it. If Task 3 shows the
cavity moves the equilibrium by more than the gate's own width, say so and stop:
that finding outranks the truncated number.

## 3. Method, standard states and checks

**3.1 Tiers.** Isolated small species (Tasks 1, 2, 5): geometry and frequencies
with `! r2SCAN-3c CPCM(water) Opt Freq TightSCF`, then a single point
`! DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C CPCM(water) TightSCF`. The 140-atom
catalyst (Task 3, 4): the QM/XTB tier of `COMPUTATIONAL.md`
(`! QM/XTB r2SCAN-3c ddCOSMO(Water)`), with the DLPNO-CCSD(T) ONIOM-style
correction on the extracted `job.QMRegion.xyz`, exactly as C8 does.

**3.2 The free energy of every species** is the DLPNO-CCSD(T) electronic energy
plus the r2SCAN-3c thermal correction to G at the temperature in question.

**3.3 Temperatures.** 288.15, 298.15 and 313.15 K. The measured block runs
15–40 °C. Obtain the three from ONE Hessian per species rather than three
optimisations; confirm with the validator how this ORCA version is asked for
that, and record the keyword you used in the log.

**3.4 Standard states — get these right or nothing below means anything.**
- Every species' gas-phase G is on a 1 atm standard state; add
  `RT ln(24.46)` = 7.93 kJ/mol at 298.15 K (scale with T) per species to put it
  on 1 mol/L.
- **Water as a reagent is 55.34 mol/L, not 1 mol/L.** Subtract `RT ln(55.34)` =
  9.95 kJ/mol at 298.15 K per water consumed (9.94 until 2026-09-17 -- RT ln(55.34)
  is 9.9493, and §10.3 takes it from `orca_io.solvent_standard_state_correction`
  rather than from this line).
- State both corrections explicitly in every reported ΔG, as a line of
  arithmetic, so a reader can check them.

**3.5 Explicit water.** Run every hydration equilibrium BOTH ways:
`n = 0` (the bare cluster-continuum reaction `K + H₂O ⇌ Kh`) and `n = 2` (two
spectator waters hydrogen-bonded to the carbonyl/diol, the SAME two on both
sides). Task 1 decides which of the two reproduces experiment; Task 2 then uses
that one and reports the other beside it.

**3.6 The connectivity check, after every optimisation.** A small script in
`computational/monitor/` (add it there; that directory already exists) that
reads the final `.xyz` and reports: the molecular formula; every bond under the
covalent cutoffs; and specifically whether the C–O(H) bonds that were supposed
to form did. Plus: a minimum has **zero** imaginary frequencies, a TS has
**exactly one**, and its mode must be the motion claimed — visualise it before
believing it.

## 4. Assumptions and their stop triggers

| id | assumption | trigger → STOP and report |
|---|---|---|
| A1 | the method reproduces known hydration equilibria | any calibrant's computed log₁₀ K_hyd differs from its measured value by more than 1.5 log units, at both n = 0 and n = 2 |
| A2 | the truncated active site is hydrated enough to be the induction's amplitude | the hydrate fraction at 288.15 K is below 0.8 |
| A3 | the cavity does not overturn the truncated answer | Task 3 moves log₁₀ K_hyd by more than 1.5 log units |
| A4 | the dehydration barrier matches the measurement | ΔG‡ differs from 92.6 kJ/mol by more than 20 kJ/mol |
| A5 | adding peroxide is not the slow step | ΔG‡(K + H₂O₂ → KP) exceeds ΔG‡(dehydration) |

**A1 and A2 are the cheap ones and they come first.** A failure of either
eliminates the hydrate candidate rather than weakening it, and that is a result
worth having for two days of compute.

## 5. Tasks and stops

| task | content | ends with |
|---|---|---|
| 1 | calibration: K_hyd of four small carbonyls, n = 0 and n = 2 | **STOP 1** |
| 2 | the truncated active site's K_hyd at three temperatures | **STOP 2** |
| 3 | the cavity correction on the 140-atom catalyst | STOP 3 |
| 4 | dehydration barriers: water-assisted against phosphate-assisted; and K + H₂O₂ | STOP 4 |
| 5 | the hydrate's pK_a by isodesmic exchange | STOP 5 |

**Tasks 3-5 are specified here only to the level of their gates.** They are
written out properly in an amendment after STOP 2, because what they should be
depends on what Tasks 1 and 2 return. Do not start them.

---

## 6. Task 1 — does the method reproduce known hydration equilibria?

**Goal.** Fix the method against measurement before any unknown is computed.

**The calibrants**, chosen to span about six orders of magnitude and to bracket
the active site electronically:

| molecule | why |
|---|---|
| formaldehyde | almost completely hydrated; the top of the range |
| acetaldehyde | near unity |
| acetone | a plain dialkyl ketone; the bottom of the range |
| 1,3-dichloropropan-2-one | the electron-poor analogue of the active site |

**Measured values are NOT stated in this plan on purpose.** Take each
molecule's K_hyd from a primary source, record the value, the source and the
temperature in the log entry, and flag any you cannot source rather than
substituting a remembered number. The gate is applied to what you record.

**Build.** Write each geometry by hand at rough standard bond lengths, optimise,
then run the section 3.6 check — a hand-built geometry that optimised to the
wrong isomer is the failure mode here. For each molecule: the carbonyl, the
gem-diol, and both with n = 2 spectator waters.

**Run.** Section 3.1's isolated-species tier, at 298.15 K.

**Report** (template C, section 7), one table:

| molecule | measured log₁₀ K_hyd | computed, n = 0 | computed, n = 2 | difference |

**Gate.** A1 → STOP and report. Otherwise commit, then **STOP 1**, whose
question is: "The calibration is above. n = <0 or 2> reproduces measurement
best. Run Task 2 with it?"

## 7. Task 2 — is the active site hydrated enough to be the induction?

**Goal.** The one number that decides whether the hydrate candidate survives.

**Species.** 1,3-dimethoxypropan-2-one and its gem-diol, with the explicit-water
count Task 1 selected, at 288.15, 298.15 and 313.15 K.

**Report.** K_hyd and the hydrate fraction `K_hyd/(1 + K_hyd)` at each
temperature, with the standard-state arithmetic written out; the van 't Hoff
ΔH° across the three temperatures; and the same at the other n, beside it.

**Gate.** A2 → STOP and report: the hydrate candidate is eliminated, C8's
perhydrate trap is the next in line, and the plan ends there. Otherwise commit,
then **STOP 2**, whose question is: "The truncated active site's hydrate
fraction is above. Write the amendment for Tasks 3-5?"

## 8. Tasks 3-5 — specified only by their gates

- **Task 3, the cavity.** Build Kh from the tracked 140-atom geometry by adding
  water across C128's carbonyl, on both faces, and optimise at the QM/XTB tier
  with the DLPNO-CCSD(T) correction on the QM region. Compare with Task 2. A3
  is the trigger. No full Hessian on the macrocycle is attempted; say so, and
  treat the comparison as electronic plus solvation only.
- **Task 4, the barriers, and the question the kinetics could not answer.**
  ΔG‡ for `Kh → K + H₂O` water-assisted, and the same with **one phosphate
  dianion** in the model, at 50-80 mM-relevant geometry — the block is
  phosphate at pH 7.00. A4 is the trigger on the water-assisted barrier.
  **The phosphate-assisted barrier against the water-assisted one is the C0/C1
  discriminator**: a general base that lowers this barrier materially puts the
  buffer in the activation step, which is what candidate C1 asserts and what
  the curves could not test. Also `K + H₂O₂ → KP` for A5.
- **Task 5, the pK_a.** Isodesmic proton exchange against TWO references whose
  pK_a brackets the question (one near 9-10, one near 12-13), reporting both.
  Near 10 puts the diolate where the boric ladder sits and explains the pH sign;
  near 13 rules the diolate out. State which reference you used and its source.

## 9. Report templates

### Template C — the COMPUTATIONAL.md Log entry

Append newest-first under `## Log`, nothing else in that file:

```
### <date> — C7 Task <n>: <the gate's verdict in one clause>

What was asked: <the task's Goal, one sentence>.
Jobs: <paths under computational/C7_hydrate_clock/, one line each, with
nprocs and wall time>.
ORCA: 6.1.1.

<the task's table>

Standard states: <the arithmetic, explicitly>.
Checks: <formula, connectivity, imaginary-frequency counts>.
Gate: <A1/A2 with its number, and pass or fail>.
Could overturn this: <conformational sampling was <what you did>; the
truncation; the solvation model; any calibrant value you could not source>.

Nothing is adopted.
```

### Template S — the stop message (in chat; nothing else)

```
STOP <n> — C7 Task <n>

Done:
- <at most five bullets, each naming a job path or a file>

<the task's table>

Gate: <verbatim, pass or fail>
Checks: <formula and frequency counts, or the failure>
Could overturn this: <the task's list>

Question: <the stop's question>
```

---

## 10. Amendment 1 — every number comes through `computational/orca_io.py`

Added 2026-09-17, before Task 1 was started. It changes no chemistry and no
gate; it changes how the numbers are read, and one of the two parsers it
replaces would have silently answered this plan's own three-temperature jobs
for the wrong temperature.

**Invoke the `run-orca` skill before writing any input or reading any output.**
It carries the conventions in full; what follows is only what binds this plan.

**10.1 Do not parse an ORCA output file.** Not with a regex, not with `grep`,
not with a throwaway script. `computational/orca_io.py` reads every number,
through OPI (`orca-pi`), from ORCA's own structured property JSON.
`test_orca_io.py` is its gate and `run_gates.py` discovers it.

Two parsers were in use until 2026-09-17 and both were wrong. A regex on
`FINAL SINGLE POINT ENERGY` cannot match a **labelled** line, and ORCA closes a
QM/XTB run with four of them — on the C8 reaction complex it returned the
r2SCAN-3c QM1 region alone, `-418.878772917932`, against a true QM/QM2 total of
`-659.180810832202`, while reporting a free energy built on the total in the
same dict.

**10.2 The three-temperature hazard, which is this plan's specifically.**
§3.3 runs 288.15 / 298.15 / 313.15 K. On a job with three thermochemistry
blocks:

| | answers for |
|---|---|
| the superseded regex parser | 313.15 K — the last block |
| OPI's own `get_free_energy()` | 288.15 K — entry `[0]` |
| `orca_io.free_energy(job, T)` | the T you asked for, or it raises |

Neither of the first two says which it gave. **Never call an OPI convenience
getter in this task.** `orca_io.temperatures(job)` lists what a job computed.

**10.3 Standard states come from the functions, not from the arithmetic.**
§3.4 already says to scale with T; do it with
`orca_io.gas_to_molar_correction(T)` and
`orca_io.solvent_standard_state_correction(T)`. Over this plan's own range the
gas correction runs 7.578 / 7.926 / 8.452 kJ/mol and the water term
9.616 / 9.949 / 10.450, so the 298.15 K constants are more than a third of a
kJ/mol wrong at both ends — a systematic in the van 't Hoff slope, entering
where §1's ±20 kJ/mol gate cannot see it.

**10.4 The hydrate fraction comes from the functions too.**
`orca_io.equilibrium_constant(dG_kJ, T)` then
`orca_io.fraction_from_equilibrium_constant(K)`. Task 2's bar of 0.8 is K = 4;
the gate is checked against those two calls, not against a hand conversion.

**10.5 Report the stationary point through `orca_io.stationary_point(job)`,**
which applies `IMAGINARY_CUTOFF_CM1` (50 cm⁻¹) and **names every imaginary mode
it ignored**. Do not apply the older skills' rule that more than one imaginary
mode means the geometry is bad: the C8 reaction complex is a converged minimum
carrying imaginary modes at −19.76 and −3.55 cm⁻¹, which ORCA's own
thermochemistry discards. Quote the verdict string verbatim in the report,
including the ignored modes, so the judgement stays visible.

**10.6 Cost reporting.** §0 asks for wall time and `nprocs`. `job.cpu_seconds`
is ORCA's own summed module time; take `nprocs` from the input you wrote. Do
not read either off the log by hand.

**10.7 `load_job` refuses a job that did not converge**, including one whose
geometry optimisation stopped short, and refuses a path with no `.out` rather
than reporting it as a crash. If it raises, that is the answer — report it,
do not pass `require_converged=False` to get past it.
