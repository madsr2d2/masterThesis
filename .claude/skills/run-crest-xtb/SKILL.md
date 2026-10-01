---
name: run-crest-xtb
description: Use when setting up, launching or reading a STANDALONE CREST or xtb run in this repo -- conformer/cluster searches (crest -nci, iMTD-GC), constrained or frozen-host sampling, xtb single points, optimisations, Hessians, or checking whether a solvation cluster kept its waters in a cavity. Not for ORCA jobs (run-orca) or for watching an ORCA job (orcamon). Load before writing a crest command line or a constraint file, and before taking any number from a crest.out / xtb.out.
---

# Running and reading standalone CREST and xtb

Everything below was checked against CREST 3.0.2 and xtb 6.7.1 on this machine,
on the C8 host + 2 waters cluster (136 atoms) and a one-water xtb job. What is
only from the documentation is marked *(docs)*.

## The rule

**Numbers come out through `computational/crest_io.py`; do not grep an output.**
`test_crest_io.py` is the gate (real-run fixtures in `computational/fixtures/`)
and `run_gates.py` discovers it.

- `crest_io.load_run(dir)` raises unless `crest.out` ends `CREST terminated
  normally.` A crashed run still leaves `crest_conformers.xyz` from an earlier
  stage, so the ensemble file is never evidence that a run finished.
- `crest_io.load_xtb(dir)` raises unless `xtb.out` has `normal termination of
  xtb` **and** `xtbout.json` exists, so **always pass `--json` to xtb**. A run
  that stops early still writes `xtbopt.xyz`.
- `host_rmsd`, `cavity_occupancy` answer "did the host move" and "is the water
  still in the cavity" -- the two questions a solvation search has to answer.

## Binaries

Neither is on `PATH`. `~/crest-3.0.2/crest`, `~/xtb-6.7.1/bin/xtb`. Use the full
path in scripts and READMEs. 32 cores here; the C8 jobs use `-T 8`.

## CREST

`crest <struc.xyz> [flags]`, atom order preserved throughout, energies Eh.

| Flag | Meaning |
|---|---|
| `--gfnff` / `--gfn2` | method. GFN-FF for large clusters, GFN2 for small ones |
| `--alpb water` | implicit solvent (`-g` = GBSA). Applies to every stage |
| `-chrg N`, `-uhf N` | only when charge != 0 / unpaired electrons |
| `-T N` | threads |
| `-nci` | NCI mode: ellipsoid wall around the input, MTD bias suited to clusters; the wall is off during optimisation. `-wscal x` scales the wall axes |
| `-quick` `-squick` `-mquick` | reduced-settings searches, crude to cruder |
| `-mdlen x0.3` or `-mdlen 10` | MTD length: multiplier of the default, or ps |
| `-cinp file` | constraint file, applies to **every** xtb call: MD, MTD, opt, singlepoint |
| `-cregen f -ewin 6 -rthr 0.125 -ethr 0.05` | resort an ensemble (kcal/mol, A, kcal/mol) |
| `-mdopt f` | reoptimise every frame of an ensemble/trajectory |
| `-rmsd a.xyz b.xyz`, `-testtopo f` | standalone checks |
| `-constrain 1-130` | writes `.xcontrol.sample` and `coord.ref`; see below |

`crest --help conf|other|thermo|qcg` lists the rest. A TOML input
(`crest input.toml`, `runtype = "nci-mtd" | "imtd-gc" | "ancopt" | "sp"`,
`[[calculation.level]] method = "gfnff"`) is the documented alternative
*(docs)*; this repo has used the command line.

### Constraint file: `-cinp`

xtb's xcontrol syntax. `crest struc.xyz -constrain <atoms>` writes it for you:

```
$constrain
  atoms: 1-130
  force constant=0.05
  reference=coord.ref
$metadyn
  atoms: 131-136
$end
```

- `$constrain atoms:` restrains the listed atoms to `reference`. CREST reports
  it as `# bond constraints : 8526` -- one distance restraint per pair, so it
  holds shape, not position. Other forms *(docs)*: `distance: 2, 12, 1.85`,
  `angle:`, `dihedral:`, `bond: 1,2,auto`, `force constant=` in Eh/Bohr^2.
- `$metadyn atoms:` restricts the RMSD bias to those atoms. The log confirms
  it: `# of atoms considered for RMSDs: 6`.
- **The default force constant of 0.5 is too stiff for a 130-atom restraint.**
  `Initial geometry optimization failed!` and `crestopt.log` with no reason.
  0.05 ran. Lower it before suspecting the input.
- Verified effect: host heavy atoms 0.003 A from the seed in all conformers
  (unfrozen: 1.4-2.2 A); both seed waters stayed in the cavity (unfrozen: 0 of 58).
- A constraint biases the energies. Do not compare a restrained run's
  `crest.energies` to an unrestrained run's.

### Reading a run

| File | Content |
|---|---|
| `crest.out` | log; last line `CREST terminated normally.`; per-stage `* wall-time:` (the last is the total); `population of lowest in %` |
| `crest_conformers.xyz` | ensemble, lowest first; comment line = total energy in Eh |
| `crest_rotamers.xyz` | conformers plus rotamers (degenerate by symmetry) |
| `crest_best.xyz` | the lowest, same format |
| `crest.energies` | `index  relative kcal/mol`, one per conformer |
| `cre_members` | conformer index, degeneracy |
| `crest_dynamics.trj` | MTD trajectory of the last stage |
| `crest.restart`, `crest_0.mdrestart` | restart files |
| `crestopt.log` | the initial optimisation; read it first when that stage fails |

Lines of the summary worth quoting: `E lowest`, `ensemble free energy
(kcal/mol)`, `population of lowest in %`, `number of unique conformers`.
**The population is a GFN-FF-level Boltzmann weight of a sampled ensemble, not
a measurement**; a value near 78% from a search that never held the water where
you need it says nothing about that water.

### What `-nci` does and does not do

It samples the aggregate, and with the whole cluster in the RMSD it also samples
the host: for the cyclodextrin the 58 unconstrained conformers were mostly host
variants and the waters left the cavity. To permute waters only, freeze the host
as above.

## xtb

Always: `xtb struc.xyz --json > xtb.out`, plus what the job needs.

| Job | Flags |
|---|---|
| single point | *(none)* |
| optimise | `--opt vtight` |
| opt + frequencies | `--ohess vtight` |
| Hessian only | `--hess` |
| method | `--gfn 2` (default), `--gfnff` |
| solvent | `--alpb water [reference\|bar1M\|gsolv]`, `--gbsa`, `--cosmo` |
| charge / spin | `-c N`, `-u N` |
| constraints | `--input xcontrol.inp` with `$fix` (frozen in opt), `$constrain`, `$wall` *(docs)* |
| MD | `--omd` with a `$md` block: `temp`, `time` (ps), `step` (fs), `dump` (fs), `shake`, `nvt` *(docs)* |

`xtb --help` lists the rest. Outputs: `xtb.out`, `xtbout.json`
(`total energy`, `HOMO-LUMO gap / eV`, `partial charges`, dipole), `xtbopt.xyz`
(comment line `energy: <Eh> gnorm: <Eh/A> xtb: <version>`), `xtbopt.log` (the
trajectory), `charges`, `wbo`, and for a Hessian `hessian`, `vibspectrum`
(frequencies in cm-1, six translations/rotations first, written as 0.00) and
`g98.out`.

- The summary box `| TOTAL ENERGY`, `| TOTAL FREE ENERGY`, `| GRADIENT NORM` is
  the number that matters; `:: -> Gsolv` lines above it are components.
  `TOTAL FREE ENERGY` exists only after a Hessian: `load_xtb` returns `None`
  otherwise, never a zero.
- `# imaginary freq.` with `imag. cutoff -20 cm-1` is xtb's own count. The
  ORCA skill's 50 cm-1 cutoff argument applies here too on a floppy macrocycle.
- On a large system xtb wants the usual stack settings
  (`ulimit -s unlimited`, `OMP_STACKSIZE`) *(standard xtb advice, not needed
  on the 3-atom test)*.

## Per-job directory

One directory per job, like the ORCA ones: the seed, the constraint file and its
`coord.ref`, a `run.sh` with the exact command, and a `README.md` saying what
the seed is, what the run is for, and what was measured. See
`computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/crest_nci_fixhost_apo/`.
Run the quickest version first (`-mquick -mdlen x0.3`, tens of seconds) to prove
the input works before spending hours.

## Traps met so far

- **Force constant 0.5 on a big `$constrain atoms:`** aborts the initial
  optimisation with no message (above).
- **A stale ensemble is not a finished run.** Check the termination line.
- **Every stage prints `* wall-time:`**; the total is the last one.
- **The seed is not necessarily a GFN-FF minimum**, so a restrained run's
  energies sit above an unrestrained run's and are not comparable.
- **CREST's `coord.ref` is written in the working directory** by
  `-constrain`; it must sit next to the constraint file at run time.
