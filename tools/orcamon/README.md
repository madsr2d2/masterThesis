# orcamon

Watch and query [ORCA](https://www.faccts.de/orca/) quantum-chemistry jobs
from the terminal. There are two tiers, both built on one parser and one report model:

- **`orcamon tui`, for people**: every job under a directory, live, over
  ssh, with convergence charts, the geometry, the output tail, and a
  notification when something needs you.
- **`orcamon <command>`, for agents** (and scripts, and people in a shell):
  short, bounded answers about a job, with `--json` on every command and
  meaningful exit codes. An agent no longer greps an 85 MB `.out` file.

A person and an agent always see the same numbers: the TUI's summary pane
and `orcamon show` are two renderings of one `JobReport`.

orcamon reads the live log, so its numbers are for monitoring and triage.
For numbers you will quote, read the finished job through ORCA's own
structured output (the property JSON) instead.

## Install

```
uv tool install './tools/orcamon[tui,images]'    # or: pipx install './tools/orcamon[tui,images]'
pip install ./tools/orcamon                      # an agents-only server: standard library only
```

| extra | brings | for |
|---|---|---|
| (none) | nothing: the core and CLI are standard library only | every `orcamon <command>` |
| `tui` | textual, textual-plotext, plotext | `orcamon tui` |
| `images` | matplotlib, numpy, pillow | the pixel geometry pane (without it the pane is text) |
| `procs` | psutil | liveness on a platform without `/proc` |

Python 3.10 or later.

## For agents: the shipped skill

```
orcamon skill                        # print the agent skill (instructions + generated command reference)
orcamon skill install                # write it to ~/.claude/skills/orcamon/SKILL.md
orcamon skill install --project .    # or to ./.claude/skills/orcamon/SKILL.md
orcamon skill --check --project .    # 0 current, 1 stale or edited, 2 absent
```

The skill ships inside the package and is released with the commands it
describes. Its command reference is generated from the argparse parsers
each time, so it cannot list a flag that does not exist. The installed copy
carries `orcamon_version` in its frontmatter, and `--check` catches an
upgrade that was not followed by a reinstall.

## Status and attention flags

| status | meaning |
|---|---|
| `not run` | an input with no output |
| `queued` | the scheduler has it pending |
| `failed` | ORCA printed an error |
| `finished` | ORCA terminated normally; the flags say whether the result is good |
| `stalled?` | alive, and the RMS gradient has stopped improving |
| `running` | alive, or output still fresh when liveness is unknown |
| `quiet` | liveness unknown and no output for `--quiet-after` seconds |
| `stopped` | gone without a termination line: killed, a time limit, a lost node |

Flag codes are part of the JSON contract: `failed`, `stopped`, `stalled`,
`quiet`, `opt_not_converged`, `scan_incomplete`, `ts_hessian`, `ts_imaginary`,
`minimum_imaginary`, `qm2_errors`. Imaginary modes are counted without a
cutoff. Frequencies from a Hessian computed mid-optimization are shown with
their cycle and raise no flag.

## Exit codes

| code | meaning |
|---|---|
| 0 | success; for `show` and `wait`, the job is not failed or stopped |
| 1 | the command worked, and the job (for `ls`: any listed job) is failed or stopped |
| 2 | usage error, or a job not found or ambiguous |
| 4 | `wait` timed out |
| 130 | interrupted |

## The JSON contract

Every command's `--json` prints one document, `{"schema": 1, ...}`, and
nothing else on stdout. Units are in the key names (`_eh`, `_kj_mol`,
`_s`); an absent value is `null`. `show --json` carries the full report
under `"job"`, with these keys:

`schema`, `label`, `path`, `stem`, `run_types`, `method`, `basis`,
`charge`, `mult`, `layers`, `multilayer`, `n_atoms`, `n_qm_atoms`,
`nprocs`, `maxcore_mb`, `status`, `liveness_source`, `sched_id`, `cycle`,
`max_cycles`, `scan_step`, `scan_total`, `scan_coordinate`,
`wall_time_s`, `last_output_age_s`, `opt_converged`,
`opt_maxiter_reached`, `criteria`, `negative_eigenvalues`,
`imaginary_freqs`, `freq_cycle`, `energy_eh`, `energy_label`,
`delta_kj_mol`, `scan_max_kj_mol`, `scan_max_at`, `qm2_errors`,
`crash_lines`, `attention` (a list of `{code, message}`).

`ls --json` gives `{"schema", "root", "jobs": [...]}`, the same keys
without `crash_lines` and `criteria`. Changing the key set bumps
`schema`.

## Liveness: is the job running?

`--liveness auto|process|slurm|mtime` (default `auto`: `slurm` where
`squeue` is on PATH, `process` elsewhere).

- `process`: this machine's process table, read from `/proc` (psutil
  where there is none). A job directory with no `orca*` process is not
  running.
- `slurm`: `squeue --me`, matched to a job directory by working
  directory, exactly or as the ancestor of exactly one unfinished job. An
  ambiguous match is left unmatched and says so. A job squeue does not
  mention is looked up in the local process table. If no process is found,
  it is unknown while its output is fresh and not running after that.
- `mtime`: no process information (a copied or mounted tree). It never
  says a job is dead: a job is `running` while its output is fresh and
  `quiet` after `--quiet-after` (default 1800 s).

## The TUI over ssh

- **Geometry**: `--graphics auto|herdr|kitty|text`. `auto` uses herdr's
  graphics API if it answers, and text inside tmux or screen. Otherwise it
  asks the terminal whether it speaks the Kitty graphics protocol and uses
  kitty or text accordingly. The pixel panes are z-buffered: near atoms hide
  far ones, spheres are lit, distant atoms fade (fog) and overlapping
  silhouettes are outlined. `v` cycles the representation -- ball-and-stick,
  licorice, space-filling and wireframe -- and `h` hides hydrogens. The text
  pane draws bonds in braille and atoms as element symbols; it cannot shade
  or occlude, so it offers ball-and-stick and wireframe only, works in any
  terminal, and costs a few KB a frame.
- **Notifications**: `--notify off|bell|osc|all` (default `all`) announces a
  job that finishes, fails, stops, stalls or goes quiet, and any new flag.
  It uses the terminal bell and OSC 9 / OSC 777 desktop notifications.
  **Inside tmux the OSC sequences need `set -g allow-passthrough on`.**
  `--on-event CMD` runs a command per event with `ORCAMON_EVENT`,
  `ORCAMON_JOB`, `ORCAMON_STATUS` and `ORCAMON_MESSAGE` set: ntfy, mail,
  a webhook.
- New jobs appear within 30 s (`r` forces a rescan). A job whose
  directory disappears is shown as `gone`.
- Keys: `q` quit, `r` refresh, `m` maximize the geometry. In the chart,
  left/right scrub through geometries, `end` follows the newest, and `a`
  shows all cycles of a scan. In the geometry pane, arrows rotate,
  `[`/`]` zoom, ctrl+arrows pan, `d` shows distances, `l` hides the atom
  labels (the index numbers in pixel mode; in text mode atoms become dots),
  `v` cycles the representation, `f` toggles fog, `h` toggles hydrogens,
  `p` steps through face-on principal-axis views, and `0` resets the camera.
  `o` rocks the molecule gently (+/-15 degrees over 6 s) to show it in the
  round; it is off by default, because a steady 5 frames a second at up to
  35 KB each is noticeable over a slow ssh link.

## Saving a picture

`orcamon snapshot JOB -o view.png` renders a job's latest geometry (or the
file it wrote, when the log prints none) to a PNG -- the same renderer the
pixel panes use, with the same `--representation`, `--no-fog`, `--no-labels`
and `--distances` options and the `--elev`/`--azim` view. It needs the
`images` extra, and is how a person shares a structure or checks the pane
without a terminal.

## State cache

Parsed state is kept in `$XDG_CACHE_HOME/orcamon/` (default
`~/.cache/orcamon/`), so a finished job's output is read once. An entry is
dropped automatically when the parser's source changes or the output is
replaced or truncated. `--no-cache` bypasses it.

## Development

orcamon is a general ORCA tool, not a reader for any one project's jobs, so
its tests never read real jobs. Both are repository-root gates, discovered by
`run_gates.py` like any other:

- **`test_monitor.py`** feeds the parser snippets of ORCA's own output. It
  holds the multilayer energy (the QM/QM2 total, not the QM1 region), scans
  keyed by step and cycle, the chunk parser against the line-by-line one,
  status and attention flags, and the palette-PNG frame.
- **`test_orcamon.py`** runs every command, the report contract, the cache,
  the SLURM probe, a headless TUI and the shipped skill against a synthetic
  tree. It also checks that the repository's installed skill is current:
  after changing a command or a flag, run `orcamon skill install --project .`.

`python -m orcamon.validate ROOT` runs the parser over a real tree, for a
person to read against the outputs, and checks that every marker line
reaches the parser past its prefilters.

## Command reference

Pasted from `--help`.

```
usage: orcamon [-h] [--version] [--root DIR] [--liveness {auto,process,slurm,mtime}]
               [--quiet-after S] [--no-cache]
               COMMAND ...

`orcamon` -- the command line. Standard library only.

    orcamon [ROOT]                 the TUI (needs the `tui` extra)
    orcamon tui [ROOT]             the same
    orcamon <command> [args]       one bounded answer about a job, for agents,
                                   scripts and people in a shell

Every command takes `--json` and prints one JSON document; without it, plain
text of at most `--max-lines` lines. Exit status: 0 success, 1 the job (or,
for `ls`, any listed job) is failed or stopped, 2 usage error or no such
job, 4 `wait` timed out, 130 interrupted.

positional arguments:
  COMMAND
    tui                 the terminal UI: every job under ROOT, live
    ls                  list every job under the root with its status and flags
    show                one job's summary: identity, status, criteria, energy, flags
    conv                the geometry convergence table, one row per cycle
    energies            the energy of every geometry, or of every scan step
    geom                a geometry as XYZ: the latest, or a given cycle or scan step
    snapshot            render a job's geometry to a PNG (needs the images extra)
    freqs               imaginary and lowest real vibrational frequencies
    input               what the input asks for: run types, method, charge and
                        multiplicity
    errors              ORCA error and QM2-error lines with context
    tail                the last lines of the output, optionally filtered
    wait                block until a job meets a condition, then show it
    skill               print the agent skill, or install it for Claude Code

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
```

<details>
<summary>Every command's --help</summary>

### `orcamon tui`

```
usage: orcamon tui [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                   [--quiet-after S] [--no-cache] [--exclude GLOB]
                   [--graphics {auto,herdr,kitty,text}] [--notify {off,bell,osc,all}]
                   [--on-event CMD]
                   [ROOT]

the terminal UI: every job under ROOT, live

example: orcamon tui runs/

positional arguments:
  ROOT                  directory to watch (default: --root)

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --exclude GLOB        skip jobs and directories matching GLOB (repeatable)
  --graphics {auto,herdr,kitty,text}
                        geometry pane: herdr/kitty pixels or text; auto asks the
                        terminal and picks text inside tmux or screen (default: auto)
  --notify {off,bell,osc,all}
                        how finished/failed/stopped/stalled jobs and new flags are
                        announced: terminal bell, OSC 9/777 desktop notification, or
                        both. Inside tmux the OSC needs `set -g allow-passthrough on`
                        (default: all)
  --on-event CMD        run CMD through the shell per event, with ORCAMON_EVENT,
                        ORCAMON_JOB, ORCAMON_STATUS and ORCAMON_MESSAGE set (e.g. to
                        push to a phone)
```

### `orcamon ls`

```
usage: orcamon ls [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                  [--quiet-after S] [--no-cache] [--status S[,S]] [--attention]
                  [--since DUR] [--sort {path,age,status}] [--exclude GLOB] [--json]
                  [--max-lines N]

list every job under the root with its status and flags

example: orcamon ls --attention

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --status S[,S]        only jobs with these statuses
  --attention           only jobs with attention flags
  --since DUR           only jobs whose output changed within DUR (30m, 2h, 1d)
  --sort {path,age,status}
                        row order (default: path)
  --exclude GLOB        skip jobs and directories matching GLOB (repeatable)
  --json                print one JSON document instead of text
  --max-lines N         cap on text output lines (default: 60)
```

### `orcamon show`

```
usage: orcamon show [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                    [--quiet-after S] [--no-cache] [--json] [--max-lines N]
                    JOB

one job's summary: identity, status, criteria, energy, flags

example: orcamon show opt/ts

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --json                print one JSON document instead of text
  --max-lines N         cap on text output lines (default: 60)
```

### `orcamon conv`

```
usage: orcamon conv [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                    [--quiet-after S] [--no-cache] [--last N] [--json] [--max-lines N]
                    JOB

the geometry convergence table, one row per cycle

example: orcamon conv opt/ts --last 20

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --last N              cycles to show (default: 10)
  --json                print one JSON document instead of text
  --max-lines N         cap on text output lines (default: 60)
```

### `orcamon energies`

```
usage: orcamon energies [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                        [--quiet-after S] [--no-cache] [--last N] [--json]
                        [--max-lines N]
                        JOB

the energy of every geometry, or of every scan step

example: orcamon energies scan

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --last N              rows to show (default: 50 for an optimization, every step for a
                        scan)
  --json                print one JSON document instead of text
  --max-lines N         cap on text output lines (default: 60)
```

### `orcamon geom`

```
usage: orcamon geom [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                    [--quiet-after S] [--no-cache] [--cycle C] [--step S]
                    [--region {all,qm}] [--json]
                    JOB

a geometry as XYZ: the latest, or a given cycle or scan step

example: orcamon geom opt/ts --cycle 12 > ts12.xyz

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --cycle C             optimization cycle (within --step for a scan)
  --step S              scan step (its last cycle unless --cycle)
  --region {all,qm}     qm: only a multilayer job's high-level (QM1) atoms (default:
                        all)
  --json                print one JSON document instead of text
```

### `orcamon snapshot`

```
usage: orcamon snapshot [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                        [--quiet-after S] [--no-cache] -o FILE
                        [--representation {ball-and-stick,licorice,space-filling,wireframe}]
                        [--size WxH] [--elev DEG] [--azim DEG] [--no-labels] [--no-fog]
                        [--distances]
                        JOB

render a job's geometry to a PNG (needs the images extra)

example: orcamon snapshot opt/ts -o ts.png

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  -o FILE, --output FILE
                        where to write the PNG
  --representation {ball-and-stick,licorice,space-filling,wireframe}
                        how to draw the high-level region (default: ball-and-stick)
  --size WxH            image size in pixels (default: 900x750)
  --elev DEG            elevation angle (default: 20)
  --azim DEG            azimuth angle (default: -60)
  --no-labels           omit the atom index labels
  --no-fog              do not fade distant atoms
  --distances           label the QM-QM bond distances
```

### `orcamon freqs`

```
usage: orcamon freqs [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                     [--quiet-after S] [--no-cache] [--lowest K] [--json]
                     [--max-lines N]
                     JOB

imaginary and lowest real vibrational frequencies

example: orcamon freqs opt/ts

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --lowest K            real modes to list (default: 6)
  --json                print one JSON document instead of text
  --max-lines N         cap on text output lines (default: 60)
```

### `orcamon input`

```
usage: orcamon input [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                     [--quiet-after S] [--no-cache] [--json] [--max-lines N]
                     JOB

what the input asks for: run types, method, charge and multiplicity

example: orcamon input opt/ts

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --json                print one JSON document instead of text
  --max-lines N         cap on text output lines (default: 60)
```

### `orcamon errors`

```
usage: orcamon errors [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                      [--quiet-after S] [--no-cache] [--context N] [--json]
                      [--max-lines N]
                      JOB

ORCA error and QM2-error lines with context

example: orcamon errors opt/ts --context 5

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --context N           lines around each match (default: 3)
  --json                print one JSON document instead of text
  --max-lines N         cap on text output lines (default: 60)
```

### `orcamon tail`

```
usage: orcamon tail [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                    [--quiet-after S] [--no-cache] [-n N] [--grep REGEX] [--json]
                    JOB

the last lines of the output, optionally filtered

example: orcamon tail opt/ts --grep 'FINAL SINGLE'

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  -n N                  lines (default: 40, at most 400)
  --grep REGEX          only lines matching REGEX among the last 2000
  --json                print one JSON document instead of text
```

### `orcamon wait`

```
usage: orcamon wait [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                    [--quiet-after S] [--no-cache] --until COND [--all] [--interval S]
                    [--timeout S] [--json] [--max-lines N]
                    JOB [JOB ...]

block until a job meets a condition, then show it

example: orcamon wait opt/ts --until done

positional arguments:
  JOB                   a job directory, a .inp/.out file, or (part of) a job's label

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --until COND          done | finished | change | attention | cycle>=N | step>=N
  --all                 wait for every JOB, not the first
  --interval S          seconds between polls (default: 15.0)
  --timeout S           give up after S seconds, exit 4 (default: 540.0, under a
                        10-minute tool ceiling)
  --json                print one JSON document instead of text
  --max-lines N         cap on text output lines (default: 60)
```

### `orcamon skill`

```
usage: orcamon skill [-h] [--root DIR] [--liveness {auto,process,slurm,mtime}]
                     [--quiet-after S] [--no-cache] [--check] [--user | --project DIR]
                     [--force]
                     [install]

print the agent skill, or install it for Claude Code

example: orcamon skill install --project .

positional arguments:
  install               write the skill to .claude/skills/orcamon/SKILL.md

options:
  -h, --help            show this help message and exit
  --root DIR            where jobs are looked for (default: $ORCAMON_ROOT or .)
  --liveness {auto,process,slurm,mtime}
                        how a job is known to be running (default: auto)
  --quiet-after S       seconds without output before a job of unknown liveness is
                        called quiet (default: 1800)
  --no-cache            parse every output from scratch; neither read nor write the
                        state cache
  --check               exit 0 if the installed skill is current, 1 if stale or edited,
                        2 if absent
  --user                under ~/.claude (the default)
  --project DIR         under DIR/.claude instead
  --force               replace a SKILL.md orcamon did not write
```

</details>
