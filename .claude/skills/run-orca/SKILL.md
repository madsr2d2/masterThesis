---
name: run-orca
description: Use when setting up, running, restarting or reading ANY ORCA quantum-chemistry job in this repo — the C1–C10 tasks in COMPUTATIONAL.md, anything under computational/, geometry or TS optimisation, frequencies, DLPNO single points, hydration or perhydrate equilibria, barriers, or extracting an energy, free energy, frequency or wall time from a finished job. Load before writing any ORCA input or reading any ORCA output.
---

# Running and reading ORCA jobs

`COMPUTATIONAL.md` is the task register and the conventions; this is how to
execute against it without re-deriving anything.

## The rule

**Every number that comes out of an ORCA job comes out through
`computational/orca_io.py`. Never parse an output file yourself.**

Not with a regex, not with `grep`, not with a throwaway script — the same rule
`analyse-kinetics` states for `data/curve_metrics.py`, for the same reason and
with a worse precedent behind it. Two parsers read this project's jobs until
2026-09-17 and both were wrong, silently, in opposite directions:

- A regex on `FINAL SINGLE POINT ENERGY` **cannot match a labelled line**, and
  ORCA closes a QM/XTB run with four of them. On the C8 reaction complex it
  returned `-418.878772917932` — the r2SCAN-3c QM1 region **alone** — while
  reporting a Gibbs free energy built on the QM/QM2 total in the same breath.
  The correct total is `-659.180810832202`. Two numbers, two different
  systems, 240 Eh apart, nothing saying so.
- On a job run at more than one temperature, that parser answered for the
  **last** block and OPI's own `get_free_energy()` answers for **entry [0]**.
  On a real three-temperature job those are 313.15 K and 288.15 K. Neither
  says which it gave.

`orca_io` is built on OPI (`orca-pi`), ORCA's own Python interface, which reads
the structured property JSON `orca_2json` writes rather than the human-readable
log. `test_orca_io.py` is the gate and `run_gates.py` discovers it.

**Do not use the `orca-mcp-server` MCP tools.** All three were measured against
this project's own working input on 2026-09-17: `validate_input_syntax`
returned three false errors and zero true findings on it; `generate_input_file`
injected `PBE` and `def2-SVP` into an r2SCAN-3c job and reported no warnings;
`suggest_keywords` recommends `def2-SVP` when `def2-TZVP` is already present.

## Reading a finished job

```python
import sys; sys.path.insert(0, "computational")
import orca_io

job = orca_io.load_job("computational/C8_perhydrate_trap/.../rc")
```

`load_job` raises unless the job terminated normally, its SCF converged and —
for an optimisation — the geometry converged. It also raises when the `.out` is
missing, because OPI's status checks return `False` rather than raising for an
absent file, so a mistyped path would otherwise be indistinguishable from a
crash. Pass `require_converged=False` only to inspect a job you already know
failed.

- `job.final_energy_eh` — the **QM/QM2 total** on a QM/XTB job.
- `job.version`, `job.natoms`, `job.charge`, `job.cpu_seconds` — the version
  and cost every `COMPUTATIONAL.md` log entry has to report, structurally
  rather than off the log's last line.
- `orca_io.free_energy(job, 298.15)` — **and there is no version of this that
  does not take a temperature.** `orca_io.temperatures(job)` lists what the job
  computed; asking for one it did not raises and names the ones it has.
  Likewise `enthalpy`, `entropy`, `zero_point_energy`.
- `orca_io.stationary_point(job)` — `"minimum"`, `"transition state"` or
  `"saddle point of order N"`, past `IMAGINARY_CUTOFF_CM1` (50 cm⁻¹).

**An imaginary frequency needs a cutoff, and this is why.** The C8 reaction
complex is a converged minimum that ORCA's own thermochemistry treats as one
(`The first frequency considered to be a vibration is 8`) and it carries two
imaginary modes at −19.76 and −3.55 cm⁻¹. On a floppy 140-atom macrocycle that
is numerical noise. A rule reading "more than one imaginary mode means this is
not a clean saddle point" — which both of the older skills state — condemns
every real TS this project will produce. `stationary_point` names what it
ignored, so a mode near the cutoff stays a visible judgement call.

## Standard states are functions of temperature

`+7.93 kJ/mol` (1 atm → 1 M) and `−9.95 kJ/mol` (water at 55.34 M) are the
**298.15 K** values. Over C7's own 288.15–313.15 K range they run 7.58–8.45 and
9.62–10.45, so a single constant is more than a third of a kJ/mol wrong at both
ends of a three-temperature van 't Hoff. Use
`orca_io.gas_to_molar_correction(T)` and
`orca_io.solvent_standard_state_correction(T)`, never a literal.

`orca_io.equilibrium_constant(dG_kJ, T)` and
`fraction_from_equilibrium_constant(K)` are the route from a free energy to the
kind of number C7's gate is stated in (a hydrate **fraction**, not a rate).

## The composite energy

Every final ΔG/ΔG‡ in the C7–C10 register is
`G(QM/XTB at T) + [E(DLPNO-CCSD(T)) − E(r2SCAN-3c)]` on the extracted
`job.QMRegion.xyz` — ORCA's own subtractive QM1/QM2 scheme applied post hoc,
because r2SCAN-3c carries ~20 kJ/mol RMSE against CCSD(T) barriers, at C7's own
gate width. That is `orca_io.composite_free_energy_eh(...)`, which raises if the
two fragment jobs are not the same system. It lived only as prose in `job.inp`
comment headers until 2026-09-17, i.e. it was re-derived by hand each time.

## Writing an input

Build it with OPI's `Calculator` and its **typed** keywords (`Opt.OPTTS`,
`Scf.TIGHTSCF`, `Dlpno`, `DEF2_TZVPP_C`), so a typo is an `AttributeError`
rather than a silently different calculation. `qmatoms` takes an
`IntGroupEnd`, not a string. OPI adds `%output jsonpropfile true end` itself,
which is what makes the job readable afterwards — keep it.

**Two keywords this project needs are not in OPI's enums** and have to be
`SimpleKeyword("...")` raw strings: `ddcosmo(water)` and plain `qm/xtb`
(OPI has `QM_XTB0/1/2` and `CPCM`/`SMD` only). Those two are therefore
unvalidated — check them by eye against a known-good `job.inp` before running.

The conventions, all settled and all in `COMPUTATIONAL.md`:

| | |
|---|---|
| geometry / TS / Hessian on the full catalyst | `! QM/XTB r2SCAN-3c ddCOSMO(Water) ...` |
| solvent, QM/XTB tier | **ddCOSMO** — ALPB leaves the QM1 region in vacuum; CPCM/SMD are rejected by ORCA outright |
| solvent, isolated fragment | **CPCM(water)** |
| final energies | `! DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C CPCM(Water) TightSCF` on the QM region |
| Hessian refresh | `Calc_Hess true` **plus** `Recalc_Hess 10` for a TS, `50` for a minimum |
| resources | `nprocs 8`, `maxcore 4000` — benchmarked; more ranks is *slower* on this chip |
| layout | `computational/C<n>_<slug>/<reaction>/<stage>/`, shared species in `computational/species/` |

Use the machine by running **different reactions concurrently**, not by raising
one job's `nprocs`.

## Before trusting a geometry

**Check the connectivity, every time.** The most expensive failure this project
has had is `K+peroxide+BnOH_bridged_UNIDENTIFIED/` — an entire reaction folder
invalidated because the optimised geometry was not the species it was labelled.
A converged optimisation with a clean Hessian tells you nothing about *which*
minimum you converged to. Read the bonds off `job.xyz` and say what the
structure is, in the `job.inp` header, before any energy from it is quoted.

Keep a failed or superseded attempt in its own subdirectory with a `README.md`
saying why (`attempt1_stalled/`, `attempt2_bad_hessian/`), rather than
overwriting `job.inp` in place.

## While a job runs

`computational/monitor/` is the live view — it tails the `.out` incrementally
for cycle count, gradient history, negative-eigenvalue count and stall
detection. `orca_io` reads **finished** jobs and needs the property JSON, which
is written at the end; the two do not overlap. Run
`.venv/bin/python -m computational.monitor.validate` to check the monitor's own
parsing.

Do not kill a running ORCA job on your own initiative.

## Not installed

xtb, CREST, Psi4, NWChem, Gaussian, pyscf, ASE, RDKit, cclib. ORCA 6.1.1 at
`~/orca_6_1_1/orca`, and OPI. Anything a conformer search would want, this
machine does not have — say so rather than routing around it.
