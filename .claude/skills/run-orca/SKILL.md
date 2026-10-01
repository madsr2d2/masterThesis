---
name: run-orca
description: This project's ORCA conventions — use when setting up, launching or restarting ANY ORCA quantum-chemistry job in this repo, or taking a number from a FINISHED one for the thesis — the C1–C10 tasks in COMPUTATIONAL.md, anything under computational/, geometry or TS optimisation, frequencies, DLPNO single points, hydration or perhydrate equilibria, barriers, or an energy, free energy, frequency or wall time to quote. Load before writing any ORCA input or quoting any ORCA number. Not for checking on, triaging or waiting for a job — that is the orcamon skill.
---

# Running and reading ORCA jobs

`COMPUTATIONAL.md` is the task register and the conventions; this is how to
execute against it without re-deriving anything.

## The rule

**Every number that comes out of an ORCA job comes out through
`computational/orca_io.py`. Never parse an output file yourself.** Watching a
job that has not finished is the one exception, and it goes through orcamon
(see "While a job runs"): its numbers are for monitoring and are never quoted.

Not with a regex, not with `grep`, not with a throwaway script — the same rule
`analyse-kinetics` states for `data/curve_metrics.py`, for the same reason and
with a worse precedent behind it. Two parsers read this project's jobs until
2026-09-17 and both were wrong, silently, in opposite directions:

- A regex on `FINAL SINGLE POINT ENERGY` **cannot match a labelled line**, and
  ORCA closes a QM/XTB run with four of them. On the C8 reaction complex it
  returned `-647.909099472771` — the r2SCAN-3c QM1 region **alone** — while
  reporting a Gibbs free energy built on the QM/QM2 total in the same breath.
  The correct total is `-873.7566473144415`. Two numbers, two different
  systems, 226 Eh apart, nothing saying so.
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
(`The first frequency considered to be a vibration is 8`) and it carries an
imaginary mode at −16.27 cm⁻¹. On a floppy 140-atom macrocycle that
is numerical noise. A rule reading "an imaginary mode means this is
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

Every final ΔG/ΔG‡ in the C7–C10 register refines the cheap geometry tier with a
correlated method, because the geometry functional (B97-3c; r2SCAN-3c before
2026-09-21) carries ~20 kJ/mol RMSE against CCSD(T) barriers, at C7's own gate
width. There is **ONE route** since 2026-09-21:

**The single point is ALWAYS DLPNO-CCSD(T) as the QM1 of a second QM/XTB job on
the whole system**, so it is electrostatically embedded in the host:

    ! QM/XTB DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C ddCOSMO(Water) TightSCF

ORCA prints `Embedding Scheme ... electrostatic` and `Point charges in QM calc.
from MM atoms    120`, so the DLPNO density is polarised by the 120 host atoms
inside the whole-system ddCOSMO cavity. The number is

    orca_io.composite_free_energy_eh(geom_job, T, dlpno_qm1_job, geom_job)

— the geometry job passed twice, because the DLPNO job is a single point with no
thermochemistry: G comes from the geometry tier's frequencies plus an electronic
level change. Do NOT add a `direct_free_energy_eh` or a `qm1_level_correction_eh`
— the first is a synonym of this call, the second of `oniom_correction_eh`; both
were withdrawn on 2026-09-20 as duplicates.

**A DLPNO single point on an extracted `job.QMRegion.xyz` in CPCM is NOT used.**
The embedding is the whole point: on C8 item 1 the isolated fragment gives a
6.4 kJ/mol different correction on the RC (and 5–7 kJ/mol on the barrier) and
cannot see the host at all. That route, and `energy_dlpno-ccsdt/`, are retired.

Open, and not to be assumed: a species with NO host (C1's small molecules, C9's
`H₂O₂ + HPO₄²⁻`) has no QM2 to embed in, so a QM/XTB "ONIOM" covering the whole
system is just plain DLPNO. Whether those tasks get a host model, or are the one
place plain DLPNO + CPCM survives, is undecided.

`oniom_correction_eh(high, low)` is the level-change half on its own, and it
lived only as prose in `job.inp` headers until 2026-09-17.

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

Set the QM2 topology builder with `AutoFF_QM2_Method GFNFF` — OPI does type it,
as `BlockQmmm(autoff_qm2_method="gfnff")`. The default `XTB` builder
misperceives non-bonded QM···QM2 contacts as bonds on this folded host and
ORCA's link-atom pre-optimisation then aborts; the switch changes the topology
only (QM2 level and its Hirshfeld embedding charges are unaffected). See
`COMPUTATIONAL.md` Conventions and the 2026-09-17 log.

**A charged QM region needs its charge in TWO places.** The `* xyzfile <charge>
<mult>` line is the **QM region's** charge/mult and is what the SCF uses;
`%QMMM Charge_total <q>` is the **whole system's**. Every earlier C8 job had a
neutral QM region, so both were `0` and the distinction was invisible. The anion
model (ketone + HOO⁻) had `Charge_total -1` but `* xyzfile 0 1`, so ORCA built a
neutral QM region with 93 electrons against a singlet and aborted:
`Error : multiplicity (1) is odd and number of electrons (93) is odd ->
impossible`. The fix is `* xyzfile -1 1` (94 electrons, even); the split is then
−1 (QM) + 0 (XTB) = −1 (total).

The conventions, all settled and all in `COMPUTATIONAL.md`:

| | |
|---|---|
| geometry / TS / Hessian on the full catalyst | `! QM/XTB B97-3c ddCOSMO(Water) ...` — B97-3c since 2026-09-21 (~21% cheaper per atom than r2SCAN-3c, and the only arm in the three-arm test that gave a clean minimum) |
| QM2 topology builder | **`AutoFF_QM2_Method GFNFF`** — the default `XTB` builder misread a non-bonded QM···QM2 contact as a bond and the link-atom step aborted; topology only, QM2 level and charges unchanged |
| solvent, QM/XTB tier | **ddCOSMO** — ALPB leaves the QM1 region in vacuum; CPCM/SMD are rejected by ORCA outright |
| final energies, full catalyst | `! QM/XTB DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C ddCOSMO(Water) TightSCF` on the whole system — the **ONLY** route; the QM1 is electrostatically embedded in the host, which is the point |
| final energies, host-free species | **OPEN** — a species with no QM2 reduces to plain DLPNO (no embedding), so whether C1/C9 get a host model or keep plain DLPNO + CPCM is undecided |
| optimisation + frequency | **one job**: `! Opt Freq` / `! OptTS Freq` **`TightOpt`** — the frequency runs at the converged geometry, ORCA aborts before it if the geometry did not converge, and `TightOpt` is required for anything a barrier is quoted from |
| Hessian refresh | `Calc_Hess true` **plus** `Recalc_Hess 10` for a TS, `50` for a minimum |
| resources | `nprocs 8`, `maxcore 4000` — benchmarked; more ranks is *slower* on this chip |
| layout | `computational/C<n>_<slug>/<reaction>/<stage>/`, shared species in `computational/species/` |

Use the machine by running **different reactions concurrently**, not by raising
one job's `nprocs`.

## One job for the geometry and the frequency

**`! OptTS Freq` (and `! Opt Freq` for a minimum) runs the frequency at the
converged geometry, in the same job and the same `.out`.** Every geometry that
needs a frequency is written this way — never as an optimisation followed by a
hand-copied `NumFreq`. Two reasons, one of which already bit us:

- **There is no geometry file to carry between jobs**, so there is no way to
  point the frequency at the wrong structure. On 2026-09-25 a hand-copied
  `NumFreq` was launched with the wrong working directory and silently became a
  re-run of the *previous* step; the geometry it was meant to check sat in a
  sibling directory, untouched, with no `.out` at all. `optts_freq/` has
  nothing to get wrong: run `orca job.inp` in its own directory and the `.xyz`
  is already beside it.
- **ORCA refuses to run the frequency if the optimisation did not converge**
  — `As a subsequent Frequencies calculation has been requested / ORCA will
  abort at this point of the run`. A two-job split has no such guard: the
  second job will compute a Hessian on whatever unconverged geometry it is
  handed and report it as a finished frequency.

`Freq` is the accepted keyword (confirmed against ORCA 6.1.1 on 2026-09-25);
`NumFreq` is the same calculation. `Calc_Hess true` still supplies the seed
Hessian `OptTS` needs to identify the imaginary mode, and the trailing `Freq`
recomputes the Hessian at the final geometry — **two Hessians, exactly what the
two-job split cost**, so the win is provenance and the convergence guard, not
time. `orca_io.stationary_point` reads the resulting `.out` as it does any
other.

**For anything a barrier is quoted from, add `TightOpt`.** The default
tolerance leaves the gradients near 1e-4 Eh/bohr, and the ZPE and G are read off
the Hessian at that point, so under-convergence shows up in the thermochemistry
rather than in `E`. `TightOpt` tightens the default about three-fold
(`TolRMSG 1.0e-4 → 3.0e-5`, `TolMAXG 3.0e-4 → 1.0e-4` Eh/bohr; confirmed in the
job's own output). On the C8 TS it is what separated two consistent runs from an
anomalous one: two independent `! OptTS Freq` runs (default and `TightOpt`)
agree to 0.09 kJ/mol on G, while the older two-job optimisation had converged to
a *different low-frequency well* and differed by 30 kJ/mol on G — all of it in
the ZPE (lowest real mode 49 cm⁻¹ against 130–137). Cost is negligible (9:15
against 9:26 on that system); the geometry it buys is the one the barrier
belongs to.

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

Checking on a job, triaging one or waiting for it belongs to the **orcamon**
skill, not this one. Here the command is `.venv/bin/orcamon`, because the
venv is not on PATH. `.venv/bin/orcamon tui computational/` is the live view
for a person.

orcamon reads the live log, so its numbers are for monitoring and triage.
A number this thesis quotes comes from `orca_io` on the finished job (above),
never from orcamon.

Running the monitor headless (Textual `run_test`, screenshots) needs
`HERDR_ENV=0`: the process inherits the pane's herdr socket and would
otherwise draw the molecule over the pane you are working in.

Do not kill a running ORCA job on your own initiative.

## Not installed

Psi4, NWChem, Gaussian, pyscf, ASE, RDKit, cclib. ORCA 6.1.1 at
`~/orca_6_1_1/orca`, and OPI. **Standalone xtb 6.7.1 and CREST 3.0.2 ARE
installed** (`~/xtb-6.7.1/bin/xtb`, `~/crest-3.0.2/crest`, on no `PATH`) --
this paragraph said they were not until 2026-09-30 while the C8 CREST runs
were already in the tree. The `run-crest-xtb` skill covers them; ORCA's
SOLVATOR seeds a cluster and CREST samples it.
