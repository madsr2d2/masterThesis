---
name: orcamon
description: Monitor and inspect ORCA quantum-chemistry jobs (running or finished) from the terminal -- status, convergence, energies, geometries, frequencies, errors, and waiting for a job to finish. Use instead of grep/tail/cat on ORCA .out files.
metadata:
  orcamon_version: 0.1.0
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

`orcamon geom` takes a cycle, or a scan step and its cycle. ORCA does not
print coordinates for every cycle. When the one asked for was not printed,
the command exits 2 and names the nearest earlier cycle that was. Report
that; do not quietly use the other geometry in its place. A job that prints
no coordinates at all (ORCA's SOLVATOR writes its cluster to a file) gives the
geometry file the log names, or else its input geometry. The XYZ comment line
and the pane title say which, so report it as that file, not as a cycle.

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

## Command reference

Generated from the installed orcamon's own parser; the flags below are the ones that exist.

Global options (before or after the command):

- `--root DIR`: where jobs are looked for (default: $ORCAMON_ROOT or .)
- `--liveness auto|process|slurm|mtime`: how a job is known to be running (default: auto)
- `--quiet-after S`: seconds without output before a job of unknown liveness is called quiet (default: 1800)
- `--no-cache`: parse every output from scratch; neither read nor write the state cache

### orcamon tui

the terminal UI: every job under ROOT, live

- `ROOT`: directory to watch (default: --root)
- `--exclude GLOB`: skip jobs and directories matching GLOB (repeatable)
- `--graphics auto|herdr|kitty|text`: geometry pane: herdr/kitty pixels or text; auto asks the terminal and picks text inside tmux or screen (default: auto)
- `--notify off|bell|osc|all`: how finished/failed/stopped/stalled jobs and new flags are announced: terminal bell, OSC 9/777 desktop notification, or both. Inside tmux the OSC needs `set -g allow-passthrough on` (default: all)
- `--on-event CMD`: run CMD through the shell per event, with ORCAMON_EVENT, ORCAMON_JOB, ORCAMON_STATUS and ORCAMON_MESSAGE set (e.g. to push to a phone)

### orcamon ls

list every job under the root with its status and flags

- `--status S[,S]`: only jobs with these statuses
- `--attention`: only jobs with attention flags
- `--since DUR`: only jobs whose output changed within DUR (30m, 2h, 1d)
- `--sort path|age|status`: row order (default: path)
- `--exclude GLOB`: skip jobs and directories matching GLOB (repeatable)
- `--json`: print one JSON document instead of text
- `--max-lines N`: cap on text output lines (default: 60)

### orcamon show

one job's summary: identity, status, criteria, energy, flags

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `--json`: print one JSON document instead of text
- `--max-lines N`: cap on text output lines (default: 60)

### orcamon conv

the geometry convergence table, one row per cycle

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `--last N`: cycles to show (default: 10)
- `--json`: print one JSON document instead of text
- `--max-lines N`: cap on text output lines (default: 60)

### orcamon energies

the energy of every geometry, or of every scan step

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `--last N`: rows to show (default: 50 for an optimization, every step for a scan)
- `--json`: print one JSON document instead of text
- `--max-lines N`: cap on text output lines (default: 60)

### orcamon geom

a geometry as XYZ: the latest, or a given cycle or scan step

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `--cycle C`: optimization cycle (within --step for a scan)
- `--step S`: scan step (its last cycle unless --cycle)
- `--region all|qm`: qm: only a multilayer job's high-level (QM1) atoms (default: all)
- `--json`: print one JSON document instead of text

### orcamon snapshot

render a job's geometry to a PNG (needs the images extra)

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `-o, --output FILE`: where to write the PNG
- `--size WxH`: image size in pixels (default: 900x750)
- `--elev DEG`: elevation angle (default: 20)
- `--azim DEG`: azimuth angle (default: -60)
- `--no-labels`: omit the atom index labels
- `--no-fog`: do not fade distant atoms
- `--distances`: label the QM-QM bond distances

### orcamon freqs

imaginary and lowest real vibrational frequencies

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `--lowest K`: real modes to list (default: 6)
- `--json`: print one JSON document instead of text
- `--max-lines N`: cap on text output lines (default: 60)

### orcamon input

what the input asks for: run types, method, charge and multiplicity

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `--json`: print one JSON document instead of text
- `--max-lines N`: cap on text output lines (default: 60)

### orcamon errors

ORCA error and QM2-error lines with context

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `--context N`: lines around each match (default: 3)
- `--json`: print one JSON document instead of text
- `--max-lines N`: cap on text output lines (default: 60)

### orcamon tail

the last lines of the output, optionally filtered

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `-n N`: lines (default: 40, at most 400)
- `--grep REGEX`: only lines matching REGEX among the last 2000
- `--json`: print one JSON document instead of text

### orcamon wait

block until a job meets a condition, then show it

- `JOB`: a job directory, a .inp/.out file, or (part of) a job's label
- `--until COND`: done | finished | change | attention | cycle>=N | step>=N
- `--all`: wait for every JOB, not the first
- `--interval S`: seconds between polls (default: 15.0)
- `--timeout S`: give up after S seconds, exit 4 (default: 540.0, under a 10-minute tool ceiling)
- `--json`: print one JSON document instead of text
- `--max-lines N`: cap on text output lines (default: 60)

### orcamon skill

print the agent skill, or install it for Claude Code

- `install`: write the skill to .claude/skills/orcamon/SKILL.md
- `--check`: exit 0 if the installed skill is current, 1 if stale or edited, 2 if absent
- `--user`: under ~/.claude (the default)
- `--project DIR`: under DIR/.claude instead
- `--force`: replace a SKILL.md orcamon did not write

### Exit codes

- `0`: success; for show and wait, the job is not failed or stopped
- `1`: the command worked and the job (for ls: any listed job) is failed or stopped
- `2`: usage error, or a job not found or ambiguous
- `4`: wait timed out
- `130`: interrupted (Ctrl-C)
