# orcamon — implementation plan

> **Status, 2026-09-28: implemented**, one commit per phase, `a8eca67`
> through `694cb28`. Read the commit messages for what changed against this
> plan. Three points matter most:
> - **§1's premise about `pc` and `ts_anion` was wrong.** ORCA aborted the
>   requested frequencies after MaxIter. The imaginary modes came from
>   Hessians recomputed mid-optimization, and those are now labelled rather
>   than flagged.
> - **`squeue` is installed on this machine,** so the SLURM probe also reads
>   the local process table.
> - **A `scan_incomplete` flag was added.**
>
> Not built: the optional `orcamon watch`, and PBS. The `computational/monitor`
> paths below describe the move and are kept as written.

Written 2026-09-28 for an implementation agent. It turns `computational/monitor/`
(the Textual ORCA job monitor) into **orcamon**, one tool with two tiers:

1. **The TUI, for people.** Watch ORCA jobs running on a remote machine, over
   ssh, without leaving the terminal.
2. **`orcamon <command>`, for agents** (and for scripts and humans in a shell).
   Short, bounded, tested answers about ORCA inputs and outputs, so an agent
   stops writing its own `grep`/`tail` commands against 85 MB `.out` files.

Both tiers sit on **one parser** and **one report model**, so a person and an
agent always see the same numbers.

The app **ships its own agent skill**. `orcamon skill` prints it, and
`orcamon skill install` puts it where Claude Code finds it. The instructions
are released with the commands they describe, so they cannot drift apart
(Phase 7).

Work through the phases in order. Each phase ends in a working, committed
state with its gates passing. Section 9 lists the decisions still open. Each
has a default, so proceed with the default unless the user has said
otherwise.

---

## 1. The user story, and why each piece exists

The user runs ORCA on a remote machine (today: the homelab, reached over ssh;
possibly a SLURM cluster later) and works from a terminal, often next to a
coding agent.

| Need | Today | After this plan |
|---|---|---|
| See every job at a glance, over ssh | TUI, works | Same, plus jobs submitted after launch appear, and the status is correct |
| Know when something needs me | Only by looking | Desktop notification or bell through ssh, and an optional hook command |
| See the molecule | Pixel images only (herdr or Kitty); blank in tmux or mosh | Pixel where the terminal supports it, a text rendering everywhere else |
| Agent checks a job | Ad-hoc `grep`/`tail`, easy to misread (four energies per QM/XTB geometry, scan cycles restarting) | `orcamon show/ls/conv/energies/...`, with `--json` |
| Agent waits for a job | `sleep` + `tail` loops | `orcamon wait <job> --until done` |
| Teach an agent to use it | Nothing; the agent improvises | `orcamon skill` / `orcamon skill install`, versioned with the tool and tested against it |
| Install on another server | Lives inside the thesis repo, needs matplotlib and friends | `uv tool install`: a stdlib-only core, with the TUI and images as extras |

### Status is wrong today, not just incomplete

Fix this in Phase 2:

- `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/pc/job.out`
  hit MaxIter ("The optimization did not converge but reached the maximum
  number of optimization cycles."). ORCA then ran the requested frequencies,
  found **9 imaginary modes**, and terminated normally. The table shows it as
  green **"converged"**, because `status.py` maps "ORCA TERMINATED NORMALLY"
  to `CONVERGED`. `geometry_b973c-xtb/ts_anion/job.out` is the same case.
- Any job whose process is gone and whose output has no termination line is
  shown as **"crashed"**. On a cluster login node, `procs.py` can't see jobs
  running on compute nodes, so every running job would read "crashed".

### Boundary with `computational/orca_io.py`

It must stay. `orca_io` reads **finished** jobs through ORCA's own property
JSON (OPI), and it is the only source for numbers quoted in the thesis. Its
docstring forbids regex parsing of the log for that purpose. orcamon parses
the **live log** for monitoring and triage. The `run-orca` skill already
states the split ("the two do not overlap"); keep it true:

- orcamon's output never claims to be a publishable number. `show` and
  `energies` are for monitoring.
- Do not import `orca_io` or `orca-pi` into orcamon.
- Do not add property-JSON reading to orcamon (listed as out of scope).

---

## 2. Ground rules for the implementer

- **Interpreter:** `.venv/bin/python`, never bare `python`. See CLAUDE.md.
- **Headless TUI runs need `HERDR_ENV=0`.** The session inherits the user's
  herdr pane, and without it the app draws the molecule over the user's
  pane. If it happens anyway, clear it with:
  `.venv/bin/python -c "from computational.monitor import herdr_graphics as h; h.clear('geometry')"`
  (the module path changes after Phase 1).
- **Never kill or signal a running ORCA process.** Jobs may be running while
  you work. `ps`/`/proc` reads are fine.
- **Leave the user's uncommitted work alone.** At the time of writing:
  `COMPUTATIONAL.md`, `PLAN_C7_INDUCTION.md`, `computational/orca_io.py`,
  `test_orca_io.py`, the C8 job files and untracked C8 directories, and part of
  `.claude/skills/run-orca/SKILL.md`. Stage files by path; never `git add -A`
  or `git commit -a`. When you edit SKILL.md, commit only your hunks
  (`git add -p` is interactive and unavailable, so build a patch with
  `git diff` and apply it with `git apply --cached`), or ask the user.
- **Commit straight to `master`**, one commit per phase or sub-phase, message
  prefixed `orcamon:`. There are no feature branches in this repo.
- **Gates:**
  - `.venv/bin/python test_monitor.py` after every change.
  - `.venv/bin/python run_gates.py --only monitor` and `--only orcamon` while
    iterating.
  - The full `.venv/bin/python run_gates.py --jobs 4` before each commit. It
    takes about 11 minutes, so run it in the background. Do not run `--all`
    (the slow optimiser suite is unrelated).
- **The duplicate guard** (`data/test_curve_metrics.py`) scans every root
  `*.py`. A new root test file's top-level names must not repeat a name in
  another root module. Only `check`, `FAILURES`, `HERE` and `main` are
  permitted; prefix helpers with `_`.
- **`test_root_documents.py` checks numbers in SKILL.md.** Don't add bare
  numbers to the skill; if one is unavoidable, add a claim for it.
- **The monitor is a general ORCA tool.** No thresholds, paths or names
  specific to this thesis in orcamon's code. Show all imaginary frequencies,
  never a cutoff.
- **Match the existing style:** comments explain *why*, with the failure that
  motivated them; dataclasses; no new dependencies in the core.

---

## 3. Target architecture

```
tools/orcamon/
  pyproject.toml
  README.md                  # install, the two tiers, the command reference
  src/orcamon/
    __init__.py              # __version__
    __main__.py              # python -m orcamon  ->  cli.main()
    SKILL.md                 # the agent skill, shipped as package data (Phase 7)
    core/                    # STANDARD LIBRARY ONLY
      parser.py              # from computational/monitor/parser.py
      orca_input.py          # from computational/monitor/orca_input.py
      discovery.py           # from computational/monitor/discovery.py
      liveness.py            # NEW: /proc (+ optional psutil), SLURM, mtime
      status.py              # REWRITTEN: new Status set + attention flags
      job.py                 # Job, moved out of app.py (no Textual)
      report.py              # NEW: JobReport + plain / markup / dict renderers
      resolve.py             # NEW: job argument -> Job(s)
      cache.py               # NEW (Phase 4): on-disk JobState cache
      units.py               # EH_TO_KJ_PER_MOL, format_age, format_wall_time
    cli/                     # STANDARD LIBRARY ONLY
      __init__.py            # main(): argparse, dispatch, exit codes
      commands.py            # one function per subcommand
      skill.py               # `orcamon skill`: print / install / --check (Phase 7)
    tui/                     # needs the [tui] extra; pixel images need [images]
      app.py                 # MonitorApp, ConvergencePlot, widgets
      geometry_text.py       # NEW (Phase 5): braille/half-block molecule
      geometry_render.py     # matplotlib renderer, from the monitor
      herdr_graphics.py      # from the monitor
      kitty.py               # the Kitty escape writer, split out of app.py
      notify.py              # NEW (Phase 5): bell / OSC 9 / OSC 777 / hook
      graphics_probe.py      # NEW (Phase 5): pick herdr / kitty / text
    validate.py              # the real-tree probe, from the monitor
```

### Dependency rule

This rule is enforced by a test. `orcamon.core` and `orcamon.cli` import only
the standard library. `orcamon.tui` may import Textual, textual-plotext and
plotext. Only `tui/geometry_render.py`, `tui/herdr_graphics.py` and the Kitty
pixel path may import numpy, matplotlib or PIL, and they are imported lazily,
only when pixel graphics are selected. psutil is optional everywhere (Linux
uses `/proc`).

### Data flow

```
 files ──► core.parser.JobState ─┐
 .inp  ──► core.orca_input       ├─► core.job.Job ─► core.report.JobReport ─┬─► plain text   (cli)
 procs/scheduler/mtime ──────────┘   (+ liveness,         (one model)          ├─► dict / JSON  (cli --json)
                                       status, flags)                          └─► Rich markup  (tui summary)
```

`summary_text` in today's `app.py` becomes `report.render_markup(report)`. Its
content is unchanged, but it is built from the same `JobReport` the CLI
prints.

---

## Phase 1 — Make it a package

**Goal:** the same behaviour, now installable and importable as `orcamon`.

### Tasks

1. `git mv computational/monitor tools/orcamon/src/orcamon` to keep the
   history, then arrange the files into `core/`, `tui/` and the top level as
   in section 3.
2. Split `app.py`:
   - `Job`, `format_age`, `format_wall_time`, `cycle_label`, `geometry_shown`,
     `describe_point`, `summary_text`, `steps_text`, `EH_TO_KJ_PER_MOL`,
     `TS_RUN_TYPES`, `_CRITERIA` and `STEP_ROWS` move to `core/` (`job.py`,
     `units.py`, and for now `report.py`).
   - `summary_text` and `steps_text` keep producing Rich markup for now;
     Phase 2 splits them.
   - `STATUS_STYLE` stays in `tui/`.
3. Swap `psutil` for a stdlib `/proc` reader in `core/liveness.py`
   (`running_orca_cwds()` keeps the same signature and return value).
   - For each numeric entry in `/proc`: read `comm` (first 15 characters, so
     match `startswith("orca")` as today), `os.readlink(f"/proc/{pid}/cwd")`,
     and the start time from field 22 of `stat` divided by
     `os.sysconf("SC_CLK_TCK")`, plus the boot time from `/proc/stat` `btime`.
   - Skip entries raising `FileNotFoundError`, `PermissionError` or
     `ProcessLookupError`.
   - Fall back to psutil only when `/proc` does not exist and psutil imports.
   - The oldest `orca` process per cwd wins, as today.
4. Add `pyproject.toml`:
   - `name = "orcamon"`, `requires-python = ">=3.10"`, no core dependencies.
   - Extras: `tui = ["textual>=8", "textual-plotext>=1", "plotext>=5"]`,
     `images = ["matplotlib>=3.8", "numpy>=1.26", "pillow>=10"]`,
     `procs = ["psutil>=5"]`.
   - Script `orcamon = "orcamon.cli:main"`, src layout (setuptools or
     hatchling; pick one and keep it minimal).
5. Add a minimal `cli.main()`. For now it only launches the TUI:
   `orcamon [ROOT]` and `orcamon tui [ROOT]`. The root defaults to the
   **current directory**; the old default, the package's parent, makes no
   sense for an installed tool. If Textual is missing, exit 2 with:
   "the TUI needs the tui extra: uv tool install 'orcamon[tui]'".
6. Hook it into the repo:
   - Add `-e ./tools/orcamon[tui,images]` to `requirements.txt`, and remove
     the monitor's own pins there if they would conflict, keeping the
     textual/plotext pins as constraints.
   - Move the monitor block's comment in `requirements.txt` over.
   - Check that a fresh venv built with the README's two lines imports
     `orcamon`.
7. Update imports in `test_monitor.py` (root) to `orcamon.core.*` and
   `orcamon.tui.*`. Keep its name, since `run_gates.py` discovers it.
8. Update every reference to the old path:
   - SKILL.md "While a job runs": `orcamon tui computational/` and
     `python -m orcamon.validate`.
   - The `requirements.txt` comment.
   - The memory file `~/.claude/projects/-home-madsr2d2-masterThesis/memory/monitor-headless-no-herdr.md`
     (the clear command's module path).
   - `grep -rn "computational.monitor\|computational/monitor"` must come back
     empty outside `git log`.

### Acceptance

- `.venv/bin/python -m orcamon tui computational/` behaves exactly as
  `python -m computational.monitor.app` did.
- `test_monitor.py` passes, and so does the full gate run.
- New gate test `test_the_core_needs_only_the_standard_library` (goes in
  `test_orcamon.py`, section 8): in a subprocess, install a `sys.meta_path`
  finder that raises on `textual`, `rich`, `plotext`, `textual_plotext`,
  `matplotlib`, `numpy`, `PIL` and `psutil`. Then import `orcamon.core.*` and
  `orcamon.cli` and run `orcamon ls` on a temporary tree.

---

## Phase 2 — Correct status, attention flags, and one report model

**Goal:** status means what it says; a job's problems are listed explicitly;
the TUI and the CLI render the same `JobReport`.

### 2a. Parser additions (`core/parser.py`)

Each addition gets a `MARKER_CASES` row in `validate.py` (both `feed_line` and
`feed_text` paths) and, where it has behaviour, a `test_monitor.py` case.

- `opt_maxiter_reached: bool`: set by a line containing `The optimization did
  not converge but reached the maximum`. Real ORCA wraps the sentence onto a
  second line, so match only this prefix.
- `frequencies: list[float] | None`: **every** frequency of the last
  `VIBRATIONAL FREQUENCIES` block, in printed order, including the zero
  translations and rotations. `imaginary_freqs` stays and is derived from the
  same block. A new block replaces the old one.
- `crash_context`: when a crash marker line is fed, store it with the 5
  lines before it.
  - The `feed_text` fast path computes `tail` in bulk, so the tail deque is not
    in sync with marker processing there. Either take the context from the
    chunk text around the match, or have `orcamon errors` re-scan the file on
    demand (section 4, `errors`).
  - Pick one. On-demand scanning is simpler and is the default.
- **SCF non-convergence:** do NOT add a marker unless you have the exact text
  from a real ORCA output or the ORCA manual. None exists in this tree (`grep`
  found none). Leave a `TODO` naming the missing sample.

### 2b. Liveness (`core/liveness.py`)

```python
@dataclass
class Liveness:
    alive: bool | None      # None = this source cannot say
    source: str             # "process" | "slurm" | "mtime"
    since: float | None     # process/job start (epoch s), for wall time
    queued: bool = False    # scheduler has it pending
    sched_id: str | None = None
```

- `LivenessProbe` has one method, `snapshot() -> dict[Path, Liveness]`, taken
  once per refresh and keyed by resolved job directory.
- `ProcessProbe` wraps `running_orca_cwds()`. A job directory absent from it
  is `alive=False`.
- `MtimeProbe` never says False: `alive=None`, and the status step decides.
- The SLURM probe is Phase 6.
- Mode selection: `--liveness auto|process|slurm|mtime`, default `auto`.
  `auto` means `slurm` when `squeue` is on `PATH` (after Phase 6), otherwise
  `process`.

### 2c. Status (`core/status.py`)

Replace the enum.

| Status | value (shown) | Condition, first match wins |
|---|---|---|
| `NOT_RUN` | `not run` | no `.out` and not alive/queued |
| `QUEUED` | `queued` | no `.out`, scheduler says pending |
| `FAILED` | `failed` | an ORCA error marker was read (`crashed_marker`) |
| `FINISHED` | `finished` | `ORCA TERMINATED NORMALLY` |
| `STALLED` | `stalled?` | alive and `possibly_stalled()` |
| `RUNNING` | `running` | alive, or (`alive is None` and output younger than `QUIET_AFTER_S`) |
| `QUIET` | `quiet` | `alive is None` and output older than `QUIET_AFTER_S` |
| `STOPPED` | `stopped` | not alive, no termination line, no error marker |

- `QUIET_AFTER_S = 1800` by default, with a CLI/TUI option. Justify it in a
  comment: a DLPNO step can be silent for a long time, so the value is
  deliberately generous, and QUIET is a flag, not a verdict.
- `compute_status(state, liveness) -> Status`.
- Keep the old "crashed without marker" case as `STOPPED`: the job may have
  been killed, hit a time limit or lost its node, and none of those is an ORCA
  error.

### 2d. Attention flags (`core/status.py`)

`attention(state, inp, status) -> list[Flag]`, with `Flag(code: str,
message: str)`. The codes are part of the JSON contract.

| code | when | message (example) |
|---|---|---|
| `failed` | status FAILED | `ORCA error: <first crash line>` |
| `stopped` | status STOPPED | `ended without a termination line` |
| `stalled` | status STALLED | `RMS gradient has not improved over the last N cycles` |
| `quiet` | status QUIET | `no output for 42 min` |
| `opt_not_converged` | `opt_maxiter_reached` and not `opt_converged` | `optimization hit MaxIter (200) without converging` |
| `ts_hessian` | TS run type, not finished, latest negative-eigenvalue count ≠ 1 | `Hessian has 2 negative eigenvalues (a TS search wants 1)` |
| `ts_imaginary` | TS run type, frequencies read, count ≠ 1 | `2 imaginary frequencies (a TS wants 1)` |
| `minimum_imaginary` | an optimization that is not a TS search, frequencies read, count ≥ 1 | `9 imaginary frequencies after a minimum search` |
| `qm2_errors` | `qm2_error_count > 0` | `3 QM2 errors` |

- The TS run types are `TS_RUN_TYPES`; optimization run types are `OPT`,
  `COPT`, `ZOPT`, `GDIIS-OPT`, `LOOSEOPT`, `NORMALOPT`, `TIGHTOPT` and
  `VERYTIGHTOPT` (from `orca_input.RUN_TYPE_KEYWORDS`).
- Take run types from the input. A job with no readable input gets only the
  flags that don't need them.
- There are no numeric thresholds beyond the ones listed; the tool stays
  general.

**Acceptance:** `.venv/bin/python -m orcamon.validate` shows `pc/job.out` and
`ts_anion/job.out` as `finished` with `opt_not_converged` and
`minimum_imaginary` (or `ts_imaginary`, per run type). These are homelab-local
outputs, so this is a manual check reported in the commit message; the gate
case uses a synthetic snippet.

### 2e. Report model (`core/report.py`)

```python
@dataclass
class JobReport:
    schema: int                      # REPORT_SCHEMA = 1
    label: str; path: str; stem: str
    # identity (input first, so it exists before any output)
    run_types: list[str]; method: str; basis: str | None
    charge: int | None; mult: int | None
    layers: dict[str, dict]          # {"total": {"charge": 0, "mult": 1}, ...}
    multilayer: bool
    n_atoms: int | None; n_qm_atoms: int | None
    nprocs: int | None; maxcore_mb: int | None
    # progress
    status: str; liveness_source: str | None; sched_id: str | None
    cycle: int | None; max_cycles: int | None
    scan_step: int | None; scan_total: int | None; scan_coordinate: str | None
    wall_time_s: float | None; last_output_age_s: float | None
    opt_converged: bool; opt_maxiter_reached: bool
    # results so far
    criteria: list[dict] | None      # latest cycle: name, value, tol, ok
    negative_eigenvalues: int | None
    imaginary_freqs: list[float] | None
    energy_eh: float | None; energy_label: str | None
    delta_kj_mol: float | None       # since first geometry (non-scan)
    scan_max_kj_mol: float | None; scan_max_at: float | int | None
    qm2_errors: int
    crash_lines: list[str]
    attention: list[dict]            # [{"code":..., "message":...}]
```

- `build_report(job, now=None) -> JobReport` is the only place any of these
  values is computed. `now` is injectable for tests.
- `report.to_dict()` gives the JSON form. Keys are stable, `None` is `null`,
  and units are in the key name (`_eh`, `_kj_mol`, `_s`).
- `render_plain(report) -> str` for the CLI and `render_markup(report) -> str`
  for the TUI. Both walk the same list of `(key, text, style)` lines; the
  styles are dropped for plain.
- `steps_rows(state, n) -> list[dict]` and `render_steps_plain` /
  `render_steps_markup` do the same for the convergence table.
- The TUI's summary and steps panes call the markup renderers. Show the
  attention flags in the summary under the rule, one per line, in yellow, or
  red for `failed` and `stopped`.

### Tests (added to `test_monitor.py`, or `test_orcamon.py` if it grows too large)

- `test_status_names_the_outcome`: normal termination plus a MaxIter line
  gives `FINISHED` with an `opt_not_converged` flag; a process gone without
  a termination line gives `STOPPED`, not `FAILED`; an error marker gives
  `FAILED`; `alive=None` with old output gives `QUIET`.
- `test_attention_flags`: one synthetic case per flag code, and one negative
  case each (for example a TS with exactly 1 imaginary has no flag).
- `test_plain_and_markup_say_the_same_thing`: strip Rich markup from
  `render_markup` and compare line-for-line with `render_plain`.
- `test_report_keys_are_stable`: the set of `to_dict()` keys equals a literal
  set in the test. Changing the contract has to be deliberate: bump
  `REPORT_SCHEMA` and edit the test.

---

## Phase 3 — The agent tier: `orcamon <command>`

**Goal:** the commands below, stdlib-only, compact by default, `--json` on
every one, bounded output, meaningful exit codes. Build `ls`, `show` and
`wait` first and commit them, then the rest.

### 3a. Common behaviour

- `orcamon [--root DIR] <command> [args]`. `--root` defaults to the current
  directory. `ORCAMON_ROOT` in the environment overrides that default.
- **Job arguments** (`core/resolve.py`). A `JOB` argument is tried in this
  order:
  1. An existing path. A directory is the job in that directory: with
     several stems, it must be unique, otherwise exit 2 listing them. A
     `.inp` or `.out` file names its stem directly.
  2. Otherwise, a case-sensitive substring of a job label under `--root`.
     Exactly one match resolves; several exit 2 with up to 20 candidates;
     none exits 2.
  - Paths with spaces exist (`orca_stuff/cat/cat+H2O/ONIOM/B97-3c_XTB /ts1`);
    test one.
- **Discovery changes** (`core/discovery.py`):
  - Every `<stem>.inp` is a job, not one stem per directory.
  - Exclude ORCA's generated inputs: `_D\d+\.scfgrad\.inp$` as today, plus
    `\.scfgrad\.inp$`. `job.scfgrad.inp` exists in this tree.
  - Skip hidden directories.
  - A job's label is its directory relative to the root when the directory
    holds one job, and `dir/stem` when it holds several.
  - `--exclude GLOB` (repeatable) on `ls` and the TUI.
- **Output:** plain text, no colour, no Rich. At most `--max-lines`
  (default 60) lines; when cut, the last line says so and how to see more.
  Numbers are formatted as in the TUI.
- **`--json`:** a single JSON document on stdout, `{"schema": 1, ...}`, no
  other stdout output. Errors go to stderr as text.
- **Exit codes:**

  | code | meaning |
  |---|---|
  | 0 | success; for `show` and `wait`, the job is not failed or stopped |
  | 1 | the command worked and the job (or, for `ls`, any listed job) has status `failed` or `stopped` |
  | 2 | usage error, or a job not found or ambiguous |
  | 4 | `wait` timed out |

- **Liveness** comes from `--liveness` (Phase 2b). One probe snapshot per
  invocation (per poll for `wait`).

### 3b. Commands

In each example the job path is real and the numbers are illustrative.

**`orcamon ls [--status S[,S]] [--attention] [--since DUR] [--sort path|age|status] [--json]`**

One line per job:

```
STATUS    PROGRESS     ENERGY/Eh          IMAG  LAST   JOB
running   cyc 41/200   -873.756647 QM/QM2 -     2 min  C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts
finished  cyc 200/200  -873.801230 QM/QM2 9     3.1 d  C8_perhydrate_trap/.../geometry_r2scan3c-xtb/pc  ! opt_not_converged minimum_imaginary
stopped   step 4/15·7  -873.790001 QM/QM2 -     5.2 h  C8_perhydrate_trap/.../geometry_b973c-xtb/ts/scan_qm  ! stopped
```

- `--attention` lists only jobs with flags.
- `--since 2h` lists only jobs whose output changed in that window. Durations
  take `s`, `m`, `h` and `d`.
- Labels are shortened with `short_label` to fit `COLUMNS` (default 120).
  `--json` always carries the full label and path.
- JSON: `{"schema":1,"root":...,"jobs":[<the JobReport dict without
  crash_lines and criteria>]}`.

**`orcamon show JOB [--json]`**: the summary pane as plain text. The identity
block, the status line, the criteria, the Hessian, the imaginary frequencies,
the energy, the attention flags and the last 3 crash lines. JSON is the full
`JobReport`.

**`orcamon conv JOB [--last N] [--json]`**: the convergence table:

- Default `--last 10`.
- Columns `cycle` (`step·cycle` for scans), `dE`, `RMS grad`, `MAX grad`,
  `RMS step` and `MAX step`, each value with a `*` suffix when converged.
- A tolerance row, and a `neg eig` column when Hessian readings exist.
- JSON is a list of per-cycle dicts with value, tolerance and ok for each
  criterion.

**`orcamon energies JOB [--last N] [--json]`**: energy per geometry.

- For an optimization: `cycle`, `energy_eh`, `label`, and `dE_kj_mol` from
  the first point, default `--last 50`.
- For a scan: one row per step (from `scan_points`): `step`, the coordinate
  value when it is a one-parameter scan, `energy_eh`, `dE_kj_mol` from step 1
  and from the lowest step, with the maximum marked `<- max`. It shows all
  steps by default.
- When `points` has hit `HISTORY_LEN`, say `earliest points not kept` (JSON:
  `"truncated": true`).
- Energies are always the QM/QM2 total on multilayer jobs; print the label.

**`orcamon geom JOB [--cycle C | --step S [--cycle C]] [--region all|qm] [--json]`**

- XYZ on stdout: count line, a comment line (`label · scan step 3 · cycle 7 ·
  E = -873.75 Eh (QM/QM2)`), then atoms.
- The default is the latest point with coordinates.
- If the requested point's coordinates were never printed, exit 2 with
  `cycle 60: coordinates not printed; nearest earlier with coordinates:
  cycle 51`. Never substitute silently.
- `--region qm` writes only the QM1 atoms (`qm_atom_indices`) of a
  multilayer job.
- This is the one command exempt from `--max-lines`, because an XYZ file is
  its product. JSON: `{"atoms":[["C",x,y,z],...], "point":{...}}`.

**`orcamon freqs JOB [--lowest K] [--json]`**: every imaginary frequency, then
the `K` lowest real non-zero ones (default 6), with the count of modes. Say
`no frequency block yet` when there is none.

**`orcamon input JOB [--json]`**: the parsed input:

- Run types, method, `charge`/`mult` per layer, whether it is multilayer,
  nprocs, maxcore.
- The coordinate source: inline `* xyz`, or `*xyzfile` with its filename.
- The basis ORCA reported, if an output exists.
- It works on a job with no output.

**`orcamon errors JOB [--context N] [--json]`**:

- Scan the `.out` on demand for `_CRASH_RE` and QM2-error lines.
- Print up to 20 matches, each with `N` lines of context before and after
  (default 3), plus the QM2 error count.
- The scan reads in `READ_CHUNK_BYTES` chunks and keeps a ring buffer for
  context, so memory stays bounded on an 85 MB file.

**`orcamon tail JOB [-n N] [--grep REGEX]`**: the last `N` lines (default 40,
hard cap 400). With `--grep`, only the matching lines among the last 2000,
capped at `N`. This is the sanctioned way for an agent to look at raw output.

**`orcamon wait JOB [JOB...] --until COND [--interval S] [--timeout S] [--json]`**

- Conditions:
  - `done` (FINISHED, FAILED or STOPPED)
  - `finished`
  - `change` (any status change, or a new attention flag)
  - `cycle>=N`
  - `step>=N`
  - `attention` (any flag appears)
- With several jobs, it returns when **any** meets the condition (`--all` to
  require every one).
- Polling:
  - One `Job` per argument is kept across polls, so each poll reads only new
    bytes.
  - `--interval` defaults to 15 s.
  - The liveness snapshot is taken once per poll.
  - A SLURM query is made at most every 30 s (Phase 6).
- On return it prints which job met which condition, then `show` for that
  job. Exit codes as in 3a (4 on timeout).
- `--timeout` defaults to 540 s, below the agent Bash tool's 10-minute
  ceiling, so a foreground call can't be killed mid-wait.
- Document in the skill (Phase 7) that longer waits run in the background.
- It must handle SIGINT cleanly (exit 130, no traceback).

### Tests (`test_orcamon.py`, a new root gate)

- Build a temporary tree with synthetic outputs (reuse
  `test_monitor.SYNTHETIC_OUTPUT` by importing it, not by copying): a
  finished opt, a MaxIter opt with frequencies, a running scan (the "alive"
  state comes from an injected probe), a failed job, a job with no output,
  and a directory with a space and two stems.
- Run every command through `cli.main(argv)` in-process, capturing stdout,
  and check:
  - exit codes;
  - line bounds;
  - that `--json` parses and carries `schema`;
  - that `show`'s energy equals `ls`'s for the same job;
  - that `geom` of a point without coordinates exits 2 with the
    nearest-earlier message;
  - job resolution (unique, ambiguous, not found, path with a space).
- `wait`:
  - A background thread appends a termination line to a synthetic `.out`
    after about 1 s; `wait --until done --interval 0.2` returns 0 within
    5 s.
  - With nothing appended, `--timeout 1` returns 4.
- Rerun the stdlib-only test from Phase 1 with `ls`, `show` and `wait` in it.

---

## Phase 4 — State cache

**Goal:** `orcamon ls` over a large tree and the TUI's startup stop reparsing
every finished job's whole output.

For scale: this tree's 64 jobs are 667 MB of `.out`, and a full parse is
about 4.4 s. Every `ls` would pay that. Measure before and after, and report
both in the commit message.

- Location: `$XDG_CACHE_HOME/orcamon/` (default `~/.cache/orcamon/`), one
  file per job, named by the sha1 of the absolute `.out` path.
- Content: `pickle` of `(CACHE_KEY, JobState)`.
  - `CACHE_KEY` is the sha1 of the source of `core/parser.py`, computed at
    import. Any parser change then invalidates every entry automatically, and
    nobody has to remember to bump a version.
  - Store the `.out`'s `st_ino`, and the `offset` already in the state.
- Load: use the entry only if `CACHE_KEY` matches, the inode matches, and the
  current size ≥ `offset`. Then `update_job` reads only the appended bytes,
  as it already does. Otherwise discard the entry and parse from scratch.
- Save: after a refresh that read new bytes. Write to a temporary file, then
  `os.replace`.
  - Never cache a state that was mid-block (`_in_block()` true): the partial
    block state would be restored incorrectly. Skip the save; the next call
    reparses from the last good cache or from zero.
- `--no-cache` on every command and in the TUI. Any exception reading a
  cache entry means discard it and continue, never fail.
- The TUI uses the cache at startup and saves on exit and every 60 s.
- Tests:
  - `test_a_cached_state_equals_a_fresh_parse`: parse, save, append more
    output, load and update, then compare field by field (reuse
    `_state_fields` from `test_monitor.py`) with a fresh full parse.
  - A changed `CACHE_KEY` invalidates the entry.
  - A replaced file (new inode) invalidates it.

---

## Phase 5 — The human tier, fit for ssh

### 5a. New jobs appear

- The TUI rediscovers every `REDISCOVER_SECONDS = 30` (the `r` key forces
  it), inside the existing scan worker.
  - Add rows for new jobs; do not reorder existing rows.
  - Mark rows whose directory has gone as `gone` (dim) rather than removing
    them under the cursor.
  - Discovery costs about 0.03 s on this tree (measured), but on a
    Lustre/NFS tree `rglob` is metadata-heavy, which is why it runs every
    30 s rather than every 3 s.
- Only `stat` `.out` files of jobs not yet in a terminal status every tick.
  Recheck terminal jobs every `REDISCOVER_SECONDS`: a re-run replaces the
  file, and `update_job` already handles that.

### 5b. Notifications (`tui/notify.py`)

- Events: a job's status changes to `finished`, `failed`, `stopped`,
  `stalled?` or `quiet`, or a new attention flag appears. **No events from
  the first full scan** (it would announce every existing job).
- Delivery (`--notify off|bell|osc|all`, default `all`):
  - `bell`: `\a`.
  - `osc`: OSC 9 `\x1b]9;<text>\x07` and OSC 777
    `\x1b]777;notify;orcamon;<text>\x07`. Terminals ignore the ones they
    don't know.
  - Inside tmux (`$TMUX` set), wrap each in tmux passthrough
    `\x1bPtmux;` + sequence with every `\x1b` doubled + `\x1b\\`. Say in
    `--help` that tmux needs `set -g allow-passthrough on`.
  - Write through the same raw-output path the Kitty widget uses
    (`self.app._driver.write` or whatever `_write` does now).
- Hook: `--on-event CMD` runs `CMD` through the shell with environment
  variables `ORCAMON_EVENT`, `ORCAMON_JOB`, `ORCAMON_STATUS` and
  `ORCAMON_MESSAGE`, non-blocking (a `subprocess.Popen`, reaped later),
  with output discarded. This is how a user pushes to a phone (ntfy, email)
  from a remote server.
- The event text reads like `ts_anion: finished · 9 imaginary after a minimum
  search`.
- Event detection lives in `core/` (`events(prev_report, report) ->
  list[Event]`), so a future `orcamon watch --on-event` CLI can reuse it.
  Build the CLI `watch` only if time allows; it is optional.
- Tests:
  - `events()` on report pairs, including no events on the first scan.
  - The tmux wrapping of an OSC string (a byte-exact expected value).

### 5c. Graphics that degrade (`tui/graphics_probe.py`, `tui/geometry_text.py`)

- `--graphics auto|herdr|kitty|text`, default `auto`:
  1. `herdr` if `herdr_graphics.available()`.
  2. Otherwise `text` if `$TMUX` or `$STY` is set. Kitty graphics through a
     multiplexer is fragile; `--graphics kitty` forces it.
  3. Otherwise probe for Kitty before `App.run()`:
     - With stdin/stdout a TTY, put the terminal in raw mode (`termios`/`tty`)
       and write `\x1b_Gi=31,s=1,v=1,a=q,t=d,f=24;AAAA\x1b\\` then
       `\x1b[c`.
     - Read until the DA1 reply (`\x1b[?...c`) or 500 ms.
     - If a `\x1b_Gi=31;OK` reply arrived before DA1, use `kitty`, otherwise
       `text`.
     - Restore the terminal in a `finally`.
  4. If `[images]` is not installed, use `text`, with a one-line note in the
     geometry border title.
- **Text renderer** (`geometry_text.py`, stdlib plus Rich, no numpy):
  - Project with the same camera as `geometry_render.camera_basis`
    (reimplement the 3×3 rotation in pure Python; atom counts are in the
    hundreds).
  - Bonds are braille lines: a 2×4 dot grid per cell, Bresenham, with the
    bond cutoff from `geometry_render.BOND_CUTOFF` (move the constant to
    `core/` so both import it).
  - Atoms are their element symbol, coloured with `ELEMENT_COLORS` (also
    moved to `core/`), drawn over bonds and depth-sorted so nearer atoms win
    a cell.
  - QM-region atoms are bold; the rest dim.
  - `h` toggles hydrogens (default: shown when the atom count ≤ 60).
  - It reuses `RotatableGeometryImage`'s view state (rotate, zoom, pan and
    distances keys work the same).
  - It is a normal Textual widget that returns a `rich.text.Text`, so it
    costs ordinary cell updates, a few KB per frame.
- Kitty pixel path over ssh: give `KittyGeometryImage` the preview-then-settle
  behaviour `HerdrGeometryImage` already has (`PREVIEW_SCALE` during a key
  repeat, full quality after `_push_settled`). Frames are about 9.5 KB
  (preview) and 35 KB (full) today.
- The geometry border title names the mode: `geometry (text)`,
  `geometry (kitty)`, `geometry (herdr)`.
- Tests (`HERDR_ENV=0` set inside the test process before import):
  - `auto` picks `text` with `$TMUX` set.
  - The probe's reply parser, on recorded byte strings (a Kitty OK then DA1;
    DA1 alone; nothing).
  - The text renderer draws a 3-atom water in a 20×10 region with the O and
    both H present and at least one braille character between them.
  - A render of 300 atoms takes under 50 ms.

---

## Phase 6 — Schedulers (for cluster use)

**Goal:** jobs on compute nodes read `running`/`queued` from the login node.

- `SlurmProbe` (`core/liveness.py`):
  - Run `squeue --me --noheader --format=%i|%T|%Z|%S` with a 10 s timeout.
    `%Z` is the working directory, `%T` the state, `%S` the start time.
  - `PENDING` gives `queued=True`, `alive=False`. `RUNNING`, `COMPLETING` and
    `CONFIGURING` give `alive=True`. Anything else is ignored.
  - Match a job directory to a SLURM job:
    1. an exact working directory match;
    2. otherwise, the working directory is an ancestor of exactly one job
       directory that is not in a terminal status.
  - Several jobs under one working directory with no exact match means
    unmatched, so the job falls back to the mtime rule. Say so in `show`
    (`liveness: mtime (squeue workdir ambiguous)`).
  - Cache the result for 30 s across calls in one process.
  - Failure (non-zero exit, timeout, not found) gives an empty snapshot, with
    a stderr note once per process; never crash.
- `auto` becomes `slurm` when `squeue` is on `PATH`, else `process`.
- In `slurm` mode a job with no match is `alive=False` if it has not written
  within `QUIET_AFTER_S`, otherwise `alive=None`. A job started outside SLURM
  on the login node then isn't declared stopped.
- PBS/Torque (`qstat -f -F json`, `PBS_O_WORKDIR`): only if the user asks.
  Leave a clear extension point (`PROBES = {"process": ..., "slurm": ...}`).
- `show` and `ls --json` carry `sched_id` and `liveness_source`.
- Tests: put a fake `squeue` shell script on `PATH` in a temporary directory,
  returning canned lines (pending, running, an ambiguous workdir, garbage,
  exit 1). Check the statuses, the ambiguity fallback, and that
  `--liveness process` ignores it.

---

## Phase 7 — The shipped agent skill, documentation, and agent integration

### 7a. The skill ships inside the package

**Why it's in the package:** the skill describes commands, flags, exit codes
and JSON keys. Kept anywhere else, it goes stale the first time a command
changes. In the package it is released with the commands and tested against
them.

- **Source:** `src/orcamon/SKILL.md`, declared as package data in
  `pyproject.toml` and read at run time with
  `importlib.resources.files("orcamon") / "SKILL.md"`. It is standard library
  only, like the rest of the CLI.
  - Claude Code skill format: YAML frontmatter with `name: orcamon` and a
    `description` that says when to load it. For example: "Monitor and inspect
    ORCA quantum-chemistry jobs (running or finished) from the terminal:
    status, convergence, energies, geometries, frequencies, errors, and
    waiting for a job to finish. Use instead of grep/tail on ORCA .out
    files." Then the markdown body.
  - The body is agent-agnostic markdown, so a non-Claude agent can use the
    printed text as-is.
- **What the hand-written body says:** the judgement calls, not a flag list.
  - Start with `orcamon ls --attention` or `orcamon show JOB`, not the `.out`.
  - Pass `--json` when parsing output programmatically; read the plain text
    when reading it yourself.
  - What each status and attention flag means, and what to check next for
    each. For example, `opt_not_converged` → `conv` and `energies`;
    `ts_imaginary` → `freqs`; `failed` → `errors`.
  - Waiting: `orcamon wait JOB --until done`. For anything longer than the
    default timeout, run it in the background and let its exit wake you;
    never loop `sleep` + `tail`.
  - `orcamon geom` exits 2 rather than substituting a geometry; report that,
    don't work around it.
  - Energies are monitoring values from the log. For numbers a project will
    quote, use the project's own finished-job reader (in this repo:
    `orca_io`). Keep this generic: "the project's own reader, if it has
    one".
  - `orcamon tail --grep` is the fallback for questions the commands don't
    answer. Read the `.out` directly only after that, and bounded.
  - Remote use: `ssh HOST orcamon ...` works the same way.
  - Never kill or modify a running job.
  - The skill must stay general: no thesis-specific paths or thresholds.
- **Generated command reference:** `orcamon skill` appends a
  `## Command reference` section generated from the argparse parsers: every
  subcommand, its one-line purpose, its flags and the exit-code table. The
  hand-written body never lists flags itself, so they cannot drift.
  `skill.py` builds the section by walking the parser's subparsers; a
  `--help` dump is not good enough, so format it compactly.
- **`orcamon skill`:** prints the full skill (frontmatter, body and the
  generated reference) to stdout. It exits 0 and is exempt from
  `--max-lines`.
- **`orcamon skill install [--user | --project DIR] [--force]`:**
  - Writes the same text to `~/.claude/skills/orcamon/SKILL.md` (`--user`,
    the default) or `DIR/.claude/skills/orcamon/SKILL.md`.
  - The installed copy's frontmatter carries `metadata: {orcamon_version:
    X.Y.Z}`.
  - It refuses to overwrite a file that lacks that key (someone's own skill)
    unless `--force`.
  - It prints the path written.
- **`orcamon skill --check [--user | --project DIR]`:**
  - Exit 0 when the installed copy's `orcamon_version` equals the running
    package's `__version__` and its text is identical to what `orcamon skill`
    prints.
  - Exit 1 with a one-line reason when stale or modified ("installed skill is
    0.3.0, orcamon is 0.4.0: run `orcamon skill install`").
  - Exit 2 when not installed.
  - A full install, not a stub that says "run `orcamon skill`": the agent gets
    the instructions without an extra round trip, and `--check` covers
    staleness.
- **Tests** (in `test_orcamon.py`):
  - `test_the_skill_names_every_command`: every subcommand registered in the
    parser appears in the hand-written body at least once, and every
    `orcamon <word>` in the body is a registered subcommand.
  - `test_the_skill_examples_run`: every fenced `orcamon ...` line in the
    body runs through `cli.main(argv)` against the synthetic tree. It must not
    raise, and it returns the exit code written beside it, as a trailing
    `# exit N` comment when it isn't 0. Placeholders such as `JOB` are
    substituted by a fixed mapping in the test.
  - `test_the_skill_reference_matches_the_parser`: every flag of every
    subcommand appears in the generated reference.
  - `test_skill_install_and_check`: install into a temporary `--project`
    directory, then `--check` gives 0. Edit the file and `--check` gives 1;
    remove it and `--check` gives 2. A pre-existing file without the version
    key is refused without `--force`.
  - `test_the_skill_is_general`: the body contains no string from a small
    deny-list taken from this repo (`masterThesis`, `computational/C`,
    `C8_`, `orca_io`, `homelab`). The `orca_io` mention belongs in the repo's
    own skill.
- **Agent evaluation.** This is manual, costs a model run, and is not a gate.
  Do it once before calling Phase 7 done, and again after any substantial
  edit to the body.
  - Give a fresh subagent (no conversation context, the installed skill
    only) the synthetic tree and five questions:
    1. Which jobs need attention, and why?
    2. Has the TS search finished, and how many imaginary modes does it
       have?
    3. Where is the scan's energy maximum, in kJ/mol above step 1?
    4. Why did the failed job fail?
    5. Wait for the running job to finish, then report its final energy.
  - Pass when every answer is correct and the transcript shows `orcamon`
    commands, not `grep` or `tail` on a `.out`.
  - Record the result, and any wording change it prompted, in the commit
    message.

### 7b. Documentation and this repo's skill


- `tools/orcamon/README.md`:
  - Install: `uv tool install ./tools/orcamon[tui,images]`, or `pipx`;
    plain `pip install ./tools/orcamon` for agents-only servers.
  - The two tiers.
  - A command reference generated from `--help` (paste it in; don't invent).
  - The JSON contract with the key list.
  - The exit codes.
  - The graphics modes and the tmux passthrough note.
  - The liveness modes.
- `--help` for every command: one sentence of purpose and one example each.
- Install the shipped skill for this project:
  `orcamon skill install --project .`, and commit
  `.claude/skills/orcamon/SKILL.md`.
  - Have `test_orcamon.py` run `orcamon skill --check --project` on the repo,
    so an upgrade without a reinstall fails the gate.
  - Check that `test_root_documents.py` does not start reading the new
    file's numbers; it reads a fixed list today.
- `run-orca` SKILL.md, "While a job runs", becomes project-specific only:
  - "Use the **orcamon** skill to check or wait on a job" (no command docs
    duplicated here).
  - The `orca_io` boundary sentence stays: orcamon is for monitoring; quoted
    numbers come from `orca_io`.
  - The `HERDR_ENV=0` note stays, and the TUI command becomes
    `orcamon tui computational/`.
  - Keep bare numbers out (`test_root_documents.py`).
- Memory files:
  - Update `monitor-headless-no-herdr.md` (module paths).
  - Add one `reference` memory: "orcamon CLI: agent-facing ORCA job queries;
    prefer it to grep on .out files".
  - Add a pointer line in `MEMORY.md`.
- CLAUDE.md: no change needed, unless the gate count text is being corrected
  anyway (it says 25 gates in about 80 s; the runner finds 37). That text
  belongs to the user; mention it, don't edit it.

---

## 8. Test plan summary

| Gate file | New or changed | Covers |
|---|---|---|
| `test_monitor.py` (root, existing) | imports updated; parser cases added | parser markers, fast path equals line path, multilayer energy, scan keying, frames, input parsing, status and attention, plain and markup equivalence, report keys |
| `test_orcamon.py` (root, new) | new | stdlib-only import rule, job resolution, every CLI command, exit codes, `--json` contracts, `wait`, cache, events, tmux OSC wrapping, graphics selection, text renderer, SLURM probe |
| `test_orcamon.py`, skill section | new | the skill names every command and nothing else; its examples run with their stated exit codes; the reference matches the parser; install / `--check` / `--force`; no repo-specific strings; this repo's installed copy is current |
| Agent evaluation (manual, not a gate) | new | a fresh subagent answers five questions on the synthetic tree using `orcamon`, not grep (Phase 7a) |
| `orcamon.validate` (manual) | extended | the real tree: status and flags on `pc/job.out` and `ts_anion/job.out`, discovery counts, marker coverage |

Both root gates are discovered by `run_gates.py` automatically; confirm with
`run_gates.py --only orcamon`. The duplicate guard also sees root test files,
so keep their top-level names unique (section 2).

A headless TUI smoke test: in `test_orcamon.py`, with `HERDR_ENV=0` and
`--graphics text`, run `MonitorApp.run_test()` on the synthetic tree. Press
down, `m` and `r`, and check that the summary pane contains the selected
job's label and status. Keep it under 10 s.

---

## 9. Open decisions (defaults the implementer should use)

| Decision | Default | Why |
|---|---|---|
| Where the package lives | `tools/orcamon/` in this repo | Extracting it later to its own repo is a `git subtree split`; the user decides when. |
| Package and command name | `orcamon` | Short. Check `pip index versions orcamon` before any public release; it is local-only until then. |
| Python floor | 3.10 | `uv` supplies a modern Python on old clusters; the core stays simple enough to lower it later. |
| Default `QUIET_AFTER_S` | 1800 s | Generous, because long quiet correlated steps are normal; it is a flag, not a verdict. |
| `wait` default timeout | 540 s | Under the agent Bash tool's 10-minute ceiling. |
| Notifications default | `all` (bell + OSC) | Remote use is the point; `--notify off` exists. |
| How the skill is distributed | Inside the package: `orcamon skill` prints it, `skill install` writes the full text with a version stamp, `--check` detects staleness | Released with the commands it describes; no extra round trip for the agent. A Claude Code plugin can wrap it later. |
| SLURM | yes (Phase 6); PBS no | The homelab needs neither today; SLURM is the likely cluster. |

## 10. Out of scope

- Reading ORCA's property JSON or using OPI (that is `orca_io`'s job).
- Orbitals, densities, vibration animation, spectra: hand off to Avogadro or
  Chemcraft.
- Editing or submitting inputs, killing jobs.
- A web UI, an MCP server, several servers in one view. `ssh host orcamon ls
  --json` covers the last one well enough for now.
- Any threshold specific to this thesis's chemistry.
