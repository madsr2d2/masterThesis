---
name: orcamon
description: Monitor and inspect ORCA quantum-chemistry jobs (running or finished) from the terminal -- status, convergence, energies, geometries, frequencies, errors, and waiting for a job to finish. Use instead of grep/tail/cat on ORCA .out files.
---

# orcamon: ORCA jobs from the command line

`orcamon` reads ORCA inputs and outputs and answers in a few bounded lines.
Use it instead of `grep`, `tail` or `cat` on `.out` files. Those files run to
hundreds of MB and are easy to misread: a multilayer job prints several
energies per geometry, a relaxed scan restarts its cycle count at every step,
and an output ending "ORCA TERMINATED NORMALLY" can still hold an
unconverged optimization.

A JOB argument is a job directory, a `.inp`/`.out` file, or part of a job's
label as `orcamon ls` prints it. Several matches are refused, with the
candidates listed. Exit 2 means "name it more precisely", so do not pick one
yourself. Commands look under `--root` (default: the current directory).

## Start here

```
orcamon ls --attention        # every job that needs a look, and why
orcamon show JOB              # one job: what it is, where it got, what is wrong
```

In `ls`, PROGRESS `cyc 41/200` is cycle 41 of at most 200 (MaxIter), and
`step 4/15·7` is scan step 4 of 15, cycle 7 within it. IMAG counts the
imaginary modes of the job's final frequencies: `0` means none, `-` means
there are no final frequencies. `ls` exits 1 when any listed job is failed or
stopped. That is information about the jobs, not an error in the command.

Read the plain text yourself. Add `--json` when a program parses the output.
Each command then prints one JSON document with a `schema` key; units are in
the key names (`_eh`, `_kj_mol`, `_s`), and `attention` is a list of
`{code, message}`.

## Status

| status | means | look next |
|---|---|---|
| `running` | alive (or output still fresh when liveness is unknown) | `orcamon conv JOB` |
| `stalled?` | alive, but the RMS gradient has stopped improving | `orcamon conv JOB`, `orcamon energies JOB` |
| `quiet` | liveness unknown and no output for a while; could be a long silent step | `orcamon tail JOB` |
| `queued` | the scheduler has it pending | nothing yet |
| `finished` | ORCA terminated normally. It does NOT mean the result is good: read the flags | the flags |
| `failed` | ORCA printed an error | `orcamon errors JOB` |
| `stopped` | gone without a termination line (killed, time limit, lost node) | `orcamon tail JOB` |
| `not run` | an input with no output | `orcamon input JOB` |

## Attention flags

| code | means | look next |
|---|---|---|
| `failed` | ORCA error, first error line quoted | `orcamon errors JOB` |
| `stopped` | ended without a termination line | `orcamon tail JOB` |
| `stalled` | gradient not improving | `orcamon conv JOB --last 30` |
| `quiet` | no output for a long time | `orcamon tail JOB` |
| `opt_not_converged` | hit MaxIter (a scan names the step) | `orcamon conv JOB`, `orcamon energies JOB` |
| `scan_incomplete` | a relaxed scan terminated normally before its last step | `orcamon energies JOB`, `orcamon tail JOB` |
| `irc_not_converged` | an IRC direction hit MaxIter before reaching a minimum | `orcamon energies JOB` |
| `ts_hessian` | a TS search whose Hessian does not have exactly one negative eigenvalue | `orcamon conv JOB` (neg eig column) |
| `ts_imaginary` | final frequencies of a TS with other than one imaginary mode | `orcamon freqs JOB` |
| `minimum_imaginary` | final frequencies of a minimum search with an imaginary mode | `orcamon freqs JOB` |
| `qm2_errors` | the low-level (QM2) region's calculation reported errors | `orcamon errors JOB` |

Imaginary modes are counted and never excused as "small". Whether a
-15 cm⁻¹ mode matters is a chemistry judgement, so report it and let the
person decide. Frequencies from a Hessian computed mid-optimization are
shown as `(Hessian at cycle N, not the result)` and raise no flag.

## Looking closer

```
orcamon conv JOB --last 20     # convergence per cycle; * marks a criterion met
orcamon energies JOB           # energy per cycle, or the scan profile per step
orcamon freqs JOB              # imaginary modes, then the lowest real ones
orcamon input JOB              # run types, method, charge and multiplicity per layer
orcamon errors JOB             # error lines with context (--context N widens it), QM2 error count
orcamon geom JOB > latest.xyz  # the latest geometry as XYZ
```

On a multilayer (QM/MM, QM/XTB, ONIOM) job the energy is the combined total,
labelled `QM/QM2`, never the high-level region alone.

`orcamon snapshot JOB -o view.png` renders the geometry to a PNG for a person
to look at (needs the `images` extra). It is for showing someone a structure,
not for answering a question about it.

`orcamon freqs JOB --mode N` names the atoms that move most in mode N — ORCA's
mode number, imaginary or real — each with its displacement relative to the
largest (`--top K`, default 8). Use it to say in words what a mode does.
`orcamon snapshot JOB --mode N -o file.png` draws any of those modes for a
person.

`orcamon geom` takes a cycle, or a scan step and its cycle. ORCA does not
print coordinates for every cycle. When the one asked for was not printed,
the command exits 2 and names the nearest earlier cycle that was. Report
that; do not quietly use the other geometry in its place. A job that prints
no coordinates at all (ORCA's SOLVATOR writes its cluster to a file) gives the
geometry file the log names, or else its input geometry. The XYZ comment line
and the pane title say which, so report it as that file, not as a cycle.

On an IRC or NEB job, `orcamon energies JOB` is the path, one row per point:
the dE is from the TS (IRC) or from image 0 (NEB), and an IRC row carries its
monitored internals. `orcamon geom JOB --point N` gives one point's geometry:
N is the signed IRC point (0 is the TS, negative is backward) or the NEB
image. An NEB path is the CURRENT one: ORCA rewrites it every iteration.

## Waiting

```
orcamon wait JOB --until done
```

`done` means finished, failed or stopped; `finished` means ORCA terminated
normally. The exit code follows the status only, so a job can finish, exit
0 and still carry a flag: read the flags in the summary `wait` prints.
Waiting blocks and polls. Its default timeout sits below a 10-minute
tool-call ceiling. On timeout it exits 4 and prints where each job stands.
To wait longer, run the command in the background and let its exit wake you.
Never loop `sleep` and `tail` against the output. The conditions are listed
in the reference below. `--until change` returns on a new status or a new
flag.

## The log is for monitoring

Energies and frequencies here are read off the live log: they are for
watching and triage. A number a project will quote, compare or publish
should come from the project's own reader of finished jobs, if it has one
(for example ORCA's property JSON), not from this tool.

## When the commands do not answer

When `orcamon errors` shows only ORCA's error line and nothing before it
explains it, the log records no cause. Report that line as the failure and
say the log gives no more; do not go looking for a cause the output does
not contain.

`orcamon tail JOB --grep REGEX` searches the end of the output and is the
sanctioned way to look at raw lines. Read the `.out` directly only after that
and only a bounded slice of it, never the whole file.

## Other rules

- Another machine: `ssh HOST orcamon ls --json` works the same way.
- Where `squeue` is installed, liveness comes from SLURM, with this machine's
  own process table for anything SLURM does not list. `show` prints a
  `liveness:` line whenever the answer did not come from a local process,
  and says so when it could not decide. See `--liveness` in the reference.
- `orcamon tui` is the interactive view for a person. An agent should not
  start it.
- Never kill, signal, edit or resubmit a running job on your own initiative.
  orcamon only reads.
- `orcamon skill --check` says whether this copy of the skill matches the
  installed orcamon; if it is stale, `orcamon skill install` refreshes it.
