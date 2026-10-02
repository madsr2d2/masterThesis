# orcamon — the agent starts and stops ORCA when the person says so (AR1–AR3)

**Written:** 2026-10-02. **Verified against:** `master` at `d232ec5` (the Act-geometry run, `PLAN_ORCAMON_ACT_GEOMETRY.md`, is committed and complete). The runner below was prototyped on a scratch copy against this commit with the fake ORCA of § Facts — started, found, stopped by SIGTERM and by SIGKILL, and failing at once — and the ORCA, MPI and SLURM facts were read off this machine.
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_ACT_RUN_PROGRESS.md` (template at the end).

The last step of the loop is the one that spends compute: the person says "run it" and the agent starts ORCA, or "stop that one" and the agent stops it. Today the agent does both by hand — `orca job.inp > job.out &` in some shell, `kill` on some pid — with nothing between it and the mistakes the project has already paid for: an input with an error, a job started twice, an output overwritten, a run killed by the agent's own shell going away, a process tree left half alive. This plan adds `orcamon run`, which checks the input first, refuses a job that is running or already has an output (unless `--rerun`, which keeps the old output), starts ORCA by its full path, detached from the agent's shell, with the output beside the input, watches the first seconds for an immediate failure, and shows the new job in the person's TUI; and `orcamon stop`, which ends a job's whole process tree, gently first. The skill's rule becomes: run or stop a job only when the person has asked for it. **The goal a deviation must still meet:** an agent can start and stop ORCA on this machine with one command each, and nothing it starts or stops is anything other than the job it named; no test or manual command in this plan starts the real ORCA or signals a process it did not start; and every existing command, key, output and JSON document behaves as before.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 266 `pass` lines, last line `0 failure(s)` | ~3 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_run/xdg XDG_STATE_HOME=/tmp/orcamon_run/state .venv/bin/python test_orcamon.py` | 285 `pass` lines, last line `0 failure(s)` | ~35 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~45 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_run/xdg XDG_STATE_HOME=/tmp/orcamon_run/state .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (the same for `test_orcamon.py`, with its environment variables). Each task states the counts it must end at.

Every task runs the first three commands. Only AR3 runs `run_gates.py`: every task touches only `tools/orcamon/`, `test_orcamon.py` and `.claude/skills/orcamon/SKILL.md`, which no other gate reads except the duplicate guard.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's; a new gate file would change the count and break that gate. New tests go in `test_orcamon.py`.

`test_orcamon.py` runs its tests from an explicit list in `if __name__ == "__main__":` at the bottom of the file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- a hand-derived expected value in this plan does not hold after the change is implemented as written (report the measured value; do not change the expected value or the design to make it pass);
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance;
- a test leaves a fake-ORCA process running after it ends (report the pids; do not leave it);
- a real-tree command differs from this plan's (report it);
- in AR3, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- **Never start the real ORCA.** Every `orcamon run` against the real tree in this plan carries `--dry-run`; tests start only the fake ORCA of § Facts, through `--orca` or `runner.start_orca`, and every test stops what it started in a `finally`.
- **Never run `orcamon stop` against the real tree**, and never kill or signal a process the task did not start.
- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's. (`tools/orcamon/README.md` and `tools/orcamon/src/orcamon/SKILL.md` are orcamon's; `.claude/skills/orcamon/SKILL.md` changes ONLY through `orcamon skill install --project .`.)
- Never `git add -A`, `git add .` or `git commit -a`; stage the files the task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file.
- Never define the same top-level name in both `test_monitor.py` and `test_orcamon.py`, private or not.
- Never let a test or a manual command read or write the user's own state or cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_run/xdg` and `XDG_STATE_HOME=/tmp/orcamon_run/state` (so `run` never queues a focus into the user's own TUI). Never write into `computational/`.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `core` and `cli` import only the standard library at module level (`test_the_core_needs_only_the_standard_library` imports every `orcamon.core` and `orcamon.cli` module with the TUI's dependencies refused).
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill.
- Never change the key set of an existing `--json` document, and never change `REPORT_SCHEMA` or `steer.STEER_SCHEMA`. `run --json` and `stop --json` are new documents.
- Never add a fenced code block containing an `orcamon …` line to `tools/orcamon/src/orcamon/SKILL.md`: `test_orcamon.test_the_skill_examples_run` executes every such line — a fenced `orcamon run` would start a job. Mention commands in prose, in backticks.

## Scope

**In scope:** `core/runner.py` (finding ORCA, telling a SLURM login node, finding a job's ORCA processes and their descendants, starting ORCA detached, stopping a process tree); `orcamon run`; `orcamon stop`; the skill's run/stop rule; README; tests with a fake ORCA.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | `run` starts ORCA on THIS machine, never through a scheduler (SLURM `submit` is later). | user, 2026-10-01 |
| D2 | The skill's rule: start a job with `orcamon run` and stop one with `orcamon stop` only when the person has asked for that job; never kill, signal, edit or resubmit a running job any other way. orcamon cannot enforce it; the skill states it. | user, 2026-10-01 |
| D3 | `run` refuses (exit 2, nothing started, nothing renamed), in this order: a SLURM login node (D6) unless `--local`; no ORCA found (D5); the job already running; an input with `check` errors; `nprocs` above 1 with no `mpirun` on PATH; an existing `<stem>.out` without `--rerun`. With `--rerun` the old output is renamed to the first free `<stem>.out.N` (N = 1, 2, …), never deleted. `--dry-run` does every check, renames nothing, starts nothing, and says what it would do. | planner, 2026-10-02 |
| D4 | ORCA is started as `[<absolute ORCA path>, "<stem>.inp"]` in the job's directory, with stdout and stderr into a NEW `<stem>.out` (opened with mode `"x"`), stdin from `/dev/null`, `start_new_session=True` (its own session, so the agent's shell exiting or timing out does not take it down) and the caller's environment. `run` then watches it for `--settle` seconds (default 3): an ORCA that exits in that time is reported with its exit status and the output's last lines, and `run` exits 1. | planner, 2026-10-02 |
| D5 | The ORCA binary is `--orca PATH`, else `$ORCAMON_ORCA`, else `orca` on PATH; resolved with `os.path.realpath`; refused unless it is an executable file whose directory holds at least 10 other entries named `orca_*` — the real ORCA's modules (92 beside this machine's), which no other program called `orca` (the GNOME screen reader is one) has. | planner, 2026-10-02 |
| D6 | A SLURM login node is a machine where `squeue` and `sinfo` are on PATH and `sinfo -h -N -o %N` (10 s timeout) lists nodes that do not include this machine's hostname or its short name. Running ORCA there would run it on the login node, so `run` refuses unless `--local`. Where `sinfo` fails or times out, nothing is refused. This machine (`homelab`) is its own only SLURM node, so it is not a login node. | planner, 2026-10-02 |
| D7 | `stop` signals only the job's own processes: every process named `orca` (its `/proc/<pid>/comm`) whose working directory is the job's directory and whose command line names `<stem>.inp` (`liveness.input_stem`), plus all their descendants. SIGTERM to all of them, then, after `--grace` seconds (default 10), SIGKILL to any still alive. A zombie counts as gone. A job SLURM runs is refused with the `scancel` command to use instead. | planner, 2026-10-02 |
| D8 | After starting a job, `run` queues an `orcamon focus` of it for the person's TUI (no waiting) when one is running on this machine and watches a root holding the job, unless `--no-focus`. | planner, 2026-10-02 |
| D9 | Commits go straight to `master`, one per task. | user, 2026-10-02 |

## Facts (verified 2026-10-02)

**Environment.** `/home/madsr2d2/masterThesis/.venv/bin/python` is Python 3.12.3; orcamon is installed editable and on PATH as `orcamon`. `git status --porcelain` was empty at `d232ec5`. Baseline: `test_monitor.py` 266 pass, `test_orcamon.py` 285 pass (both re-run for this plan), `run_gates.py` 38 gates, 0 failed. README `--help` blocks are pasted at `COLUMNS=90`.

**This machine.** ORCA 6.1.1 is `/home/madsr2d2/orca_6_1_1/orca` (an x86-64 ELF executable; `which orca` finds it); its directory holds 92 entries named `orca_*` (among them `orca_startup`, `orca_leanscf`, `orca_leanscf_mpi`, `orca_mdci`, `orca_plot`, `orca_2mkl`; there is no `orca_scf` in 6.1.1). `mpirun` is `/usr/bin/mpirun` (Open MPI 4.1.6). SLURM is installed (`/usr/bin/squeue`, `/usr/bin/sbatch`): one partition `debug`, one node `homelab`, idle; `sinfo -h -N -o %N` prints `homelab`, and `socket.gethostname()` is `homelab`. The project runs ORCA directly, `orca job.inp` in the job's own directory (`.claude/skills/run-orca/SKILL.md`), not through `sbatch`.

**The code** (at `d232ec5`). `core/liveness.py`: `is_orca_process(name)`, `input_stem(args) -> str | None` (the first argument after `argv[0]` that ends in `.inp`, as a stem), `running_orca_stems()`; it reads `/proc/<pid>/comm`, `cwd` and `cmdline`. A job's `Liveness` has `alive` (True/False/None), `source` (`process`, `slurm`, `mtime`), `queued` and `sched_id`; `_with_job` loads the job with the `--liveness` probe the command was given. `core/inputcheck.py`: `check_input(job) -> InputCheck` with `errors` and `warnings` lists of `Finding(code, message)`; the `output_exists` warning reads `f"{stem}.out exists: running this job again overwrites it"`. `cli/commands.py`: `cmd_focus` (line ~840) checks the TUI's root and writes the request with `steer.write_json(steer.focus_path(pid), {"schema", "id", "file", "point", "step", "cycle", "mode", "atoms", "requested_unix_s"})` (lines 869–880); `steer.current_view()`; `_job_exit`, `_error`, `UsageError`, `_emit_json`. `cli/__init__.py`: `EXIT_CODES` (exit 1 means `the command worked and the job (for ls: any listed job) is failed or stopped, or (check) the input has an error`) and the module docstring's exit sentence (`… 1 the job (or, for ls, any listed job) is failed or stopped, or (check) the input has an error, …`).

**The skill now** (`tools/orcamon/src/orcamon/SKILL.md`). "Changing a geometry" ends with the paragraph that starts `When the person asks for a change` and contains the words `` `orcamon focus` the draft so the person sees the new structure; and run it only when they say so. After the run, `orcamon connectivity` says what it did to the bonds.``. "Other rules" has the bullet `- Never kill, signal, edit or resubmit a running job on your own initiative.` followed by `  orcamon never changes a job: it reads, and `orcamon focus` only moves what` / `  the person's TUI shows.`

**The fake ORCA** (prototyped). A directory holding this file as `orca` (mode 755) and ten empty files `orca_0` … `orca_9`:
```python
#!/usr/bin/env python3
"""A stand-in for ORCA in orcamon's tests: it calls itself `orca`, starts one
child, and does what words in its input tell it to."""
import ctypes, signal, subprocess, sys, time
ctypes.CDLL(None).prctl(15, b"orca", 0, 0, 0)
text = open(sys.argv[1]).read()
print("fake ORCA running", sys.argv[1], flush=True)
if "FAKE_FAIL" in text:
    print("ORCA finished by error termination in SCF", flush=True)
    sys.exit(1)
if "FAKE_STUBBORN" in text:
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
subprocess.Popen(["sleep", "60"])
time.sleep(60)
```
`prctl(15, …)` is `PR_SET_NAME`: after it the process's `comm` is `orca` and its command line `python3 /…/orca job.inp`, so `input_stem` reads `job`. Measured with a prototype runner that starts it as D4 says, in-process: a plain fake is running after 0.5 s with `comm` `orca`, one driver pid for (its directory, `job`), one descendant (`sleep`), its output `fake ORCA running job.inp`; stopping it sends SIGTERM to 2 processes and both are gone within 0.1 s, none needing SIGKILL. A `FAKE_STUBBORN` fake ignores SIGTERM and so does its `sleep` child (an ignored signal stays ignored across exec), so a 0.5 s grace ends in SIGKILL to both 2. A `FAKE_FAIL` fake exits with status 1 inside 0.5 s and its output is the two lines `fake ORCA running job.inp` and `ORCA finished by error termination in SCF`. A process the test itself started becomes a zombie when killed until it is reaped: `state` `Z` in `/proc/<pid>/stat` (the field after the last `)`); `os.waitpid(pid, os.WNOHANG)` reaps it, and raises `ChildProcessError` for a pid that is not the caller's child.

**Real-tree targets for dry runs** (read-only). `C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/ts5` is "not run": only `job.inp`, `job.QMRegion.xyz`, `job.xtbhess.xyz`, `job.xyz`; `orcamon check` on it prints `0 error(s), 0 warning(s) · 140 atoms, 594 electrons, charge 0, multiplicity 1`; its input has `nprocs 12` and `%maxcore 2200`. `C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts_freq_tight` is finished and has `job.out`.

**Test helpers** (`test_orcamon.py`): `check`, `_orcamon(argv)` (in-process), `_json(argv)`, `_water_job(root, name="water")` (a finished 3-atom water job with `job.inp` and `job.out`). The last tests in the `__main__` list before the skill tests are `test_check_finds_what_would_fail` and `test_connectivity_names_bonds_made_and_broken`. Count: 285.

## Conventions

- Match the surrounding code: one-line comments where the reason is not obvious; docstrings say WHY.
- New tests go in `test_orcamon.py`, use `check(...)` for every assertion, build trees under `tempfile.TemporaryDirectory()`, and are added to its `__main__` list right after the test the task names. CLI calls in new tests pass `["--root", str(root), "--liveness", "process", "--no-cache", …]` and, for `run`, `--orca <fake>`, `--settle 0.5` and `--local` (so the tests never call this machine's `sinfo`), unless the check says otherwise.
- A draft in a test is a directory with only `job.inp` = `"! B97-3c Opt\n* xyz 0 1\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n*\n"`, with `# FAKE_FAIL` or `# FAKE_STUBBORN` as its first line where a check says so.
- Every test that starts a fake keeps the pids it saw and, in a `finally`, calls `runner.stop_tree(pids, 1.0)` on any still alive.
- A new command's README entry is a `### \`orcamon <name>\`` section in `## Command reference` holding `COLUMNS=90 .venv/bin/orcamon <name> --help` in a fenced block; whenever a command is added, the reference's first block is replaced with `COLUMNS=90 .venv/bin/orcamon --help`. A task that adds a command or edits `SKILL.md` mentions the command in `SKILL.md` in prose as `` `orcamon <name>` `` and ends with `cd /home/madsr2d2/masterThesis && .venv/bin/orcamon skill install --project /home/madsr2d2/masterThesis` (with `--force` only if it refuses).
- No network access.

---

# Phase 1 — the runner

### AR1 — `core/runner.py`: find ORCA, tell a login node, start ORCA detached, stop a process tree

**Why:** starting and stopping are the two places a mistake costs the most, so they live in one module with one way to do each (D4–D7), tested against a fake ORCA because the real one must never be started by a test.
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/runner.py` (new), `test_orcamon.py`.

**Do:**
1. Create `core/runner.py` (standard library and `orcamon.core` only; Linux `/proc`) with a module docstring giving D4–D7 in prose, `class RunError(Exception)` (the message says what to do), `MIN_ORCA_MODULES = 10` and `SINFO_TIMEOUT_S = 10.0`, each commented with why.
2. `find_orca(explicit: str | None = None) -> Path`: the candidate is `explicit`, else `os.environ.get("ORCAMON_ORCA")`, else `shutil.which("orca")`; none → `RunError("ORCA not found: pass --orca PATH or set ORCAMON_ORCA")`; resolve with `os.path.realpath`; not an executable file → `RunError(f"{path} is not an executable file")`; fewer than `MIN_ORCA_MODULES` entries named `orca_*` in its directory → `RunError(f"{path} is not ORCA: no orca_* programs beside it (pass --orca PATH or set ORCAMON_ORCA)")`. Return the resolved `Path`.
3. `slurm_login_node() -> str | None` (D6): None unless `shutil.which("squeue")` and `shutil.which("sinfo")`; run `["sinfo", "-h", "-N", "-o", "%N"]` with `capture_output=True, text=True, timeout=SINFO_TIMEOUT_S`; None on `OSError`, `subprocess.TimeoutExpired` or a non-zero exit; `nodes` = the non-empty stripped lines; `host = socket.gethostname()`, `short = host.split(".")[0]`; None when either is in `nodes`; else `f"this machine ({short}) is not one of SLURM's nodes ({', '.join(sorted(nodes)[:5])}{', ...' if len(nodes) > 5 else ''}): ORCA started here would run on the login node. Submit it to the queue instead (orcamon cannot yet), or pass --local if this machine is meant to run it"`.
4. Process helpers: `_state(pid) -> str | None` (the field after the last `)` of `/proc/<pid>/stat`, None when unreadable); `alive(pid) -> bool` (`_state` is not None and not `"Z"`); `driver_pids(job_dir, stem) -> list[int]` (D7: every numeric `/proc` entry whose `comm` is exactly `orca`, whose resolved `cwd` equals `Path(job_dir).resolve()`, whose `cmdline` split on NUL gives `liveness.input_stem(args) == stem`, and which is `alive`; unreadable entries skipped; sorted); `descendants(pid) -> list[int]` (build parent → children from every `/proc/<n>/stat`'s second field after the last `)`, then walk from `pid`; sorted, `pid` itself excluded).
5. `start_orca(orca: Path, job_dir: Path, stem: str, settle: float) -> tuple[int, int | None, list[str]]` (D4): open `job_dir / f"{stem}.out"` with mode `"x"`; `subprocess.Popen([str(orca), f"{stem}.inp"], cwd=job_dir, stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT, start_new_session=True, close_fds=True)`; close the file in the parent; poll every 0.05 s for `settle` seconds; return `(pid, returncode, last 10 lines of the output)` when it exited, else `(pid, None, [])`.
6. `stop_tree(roots: list[int], grace: float) -> tuple[list[int], list[int]]` (D7): `pids` = sorted alive members of `roots` plus all their `descendants`, collected BEFORE signalling; SIGTERM each (`ProcessLookupError` ignored; `PermissionError` re-raised as `RunError(f"cannot signal pid {pid}: not yours")`); poll every 0.1 s up to `grace`, reaping with `os.waitpid(p, os.WNOHANG)` (`ChildProcessError` ignored), until none is `alive`; SIGKILL those still alive and wait up to 2 s more, reaping the same way. Return `(pids, the ones that needed SIGKILL)`.
7. `test_orcamon.py`: `_fake_orca(directory) -> Path` (writes § Facts' fake as `directory / "orca"`, mode 755, and the ten `orca_N` files; returns the path) and `test_the_runner_finds_orca_and_stops_what_it_started()`, registered right after `test_connectivity_names_bonds_made_and_broken`. Exactly 5 checks (`from orcamon.core import runner`):
   - `runner.find_orca(str(fake)) == fake.resolve()`; with `ORCAMON_ORCA` set to the fake and no argument, the same; and a lone executable file `orca` in an otherwise empty directory raises `RunError` whose message contains `is not ORCA: no orca_* programs beside it` (restore `ORCAMON_ORCA` in a `finally`);
   - `pid, code, tail = runner.start_orca(fake, draft, "job", 0.5)` on a plain draft: `code is None`, `runner.alive(pid)`, `/proc/<pid>/comm` reads `orca`, `runner.driver_pids(draft, "job") == [pid]`, `len(runner.descendants(pid)) == 1`, and `job.out` contains `fake ORCA running job.inp`;
   - `runner.stop_tree([pid], 2.0)` returns 2 pids and `[]`, and afterwards `runner.alive(p)` is False for both and `runner.driver_pids(draft, "job") == []`;
   - a `# FAKE_STUBBORN` draft: started the same way, `runner.stop_tree([pid], 0.5)` returns 2 pids and both of them as needing SIGKILL, and neither is alive afterwards; a `# FAKE_FAIL` draft: `code == 1` and `tail[-1] == "ORCA finished by error termination in SCF"`;
   - `slurm_login_node()` with `PATH` set (and restored in a `finally`) to a temporary directory holding executable `squeue` (`#!/bin/sh\nexit 0\n`) and `sinfo` (`#!/bin/sh\necho node01\necho node02\n`): the message contains `is not one of SLURM's nodes (node01, node02)`; with that `sinfo` printing `socket.gethostname()` instead: None; with `PATH` set to an empty temporary directory: None.
8. Run the first three verification commands, and `ps -eo pid,comm,args | grep -E "fake|sleep 60" | grep -v grep` afterwards.

**Acceptance:**
- `test_the_runner_finds_orca_and_stops_what_it_started` passes its 5 checks, and the `ps` line shows no fake ORCA or `sleep 60` left behind.
- `.venv/bin/python -c "from orcamon.core import runner; print(runner.find_orca(), runner.slurm_login_node())"` prints `/home/madsr2d2/orca_6_1_1/orca None`.
- `test_orcamon.py`: 290 `pass`; `test_monitor.py`: 266 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

# Phase 2 — the commands

### AR2 — `orcamon run`

**Why:** "run it" should be one command that refuses what would go wrong and leaves ORCA running safely after the agent's shell is gone (D3, D4, D8).
**Depends on:** AR1.
**Files:** `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. `cli/__init__.py`, after `connectivity`: `add("run", "start ORCA on a job on this machine, after checking it", "orcamon run opt/ts", commands.cmd_run)`, `_job(p)`, `--orca` (`metavar="PATH"`, help `"the ORCA binary (default: $ORCAMON_ORCA, else orca on PATH)"`), `--rerun` (`store_true`, help `"run a job that has an output: the old one is kept as <stem>.out.N"`), `--local` (`store_true`, help `"run here even where SLURM says this is not one of its nodes"`), `--dry-run` (`store_true`, help `"check everything and say what would run; start nothing"`), `--settle` (`type=float`, `default=3.0`, `metavar="S"`, help `"seconds to watch for ORCA exiting at once (default: %(default)s)"`), `--no-focus` (`store_true`, help `"do not show the started job in the person's TUI"`), `_output(p, lines=False)`. Exit 1's meaning in `EXIT_CODES` gains `, or (run) ORCA exited at once`, and the docstring's exit sentence the same words.
2. `commands.py`, a helper `_queue_focus(view, file, atoms=(), point=None, step=None, cycle=None, mode=None) -> str` holding the lines that build and write the request in `cmd_focus` (the root check stays in `cmd_focus`); `cmd_focus` calls it and behaves exactly as before.
3. `@_with_job def cmd_run(args, job)`, refusing with `UsageError` in D3's order: `runner.slurm_login_node()` (skipped with `--local`); `runner.find_orca(args.orca)` (`RunError` → `UsageError`); `job.liveness` with `alive is True` or `queued`, or `runner.driver_pids(job.state.path, job.state.stem)` non-empty → `f"{job.label} is already running: orcamon stop it first, if the person asks"`; `check_input(job)` with errors → print each `f"error: {message}"` and refuse `f"the input has {n} error{'s' if n != 1 else ''}: fix them before running (orcamon check {job.label})"`; `job.input.nprocs > 1` and `shutil.which("mpirun") is None` → `f"nprocs {n} needs mpirun on PATH: ORCA runs its parallel modules through it"`; `job.out_path.exists()` and not `--rerun` → `f"{stem}.out exists: --rerun moves it to {stem}.out.{k} and runs again"` (`k` the first free number).
4. Then: the warnings of `check_input` other than `output_exists` are printed as `f"warning: {message}"`. With `--dry-run`: print `f"would move {stem}.out to {stem}.out.{k}"` when it applies, then `f"would start: {orca} {stem}.inp in {job.state.path}"` and return 0, having written nothing. Otherwise rename the old output when `--rerun` applies (`os.replace`, to the first free `.N`), call `runner.start_orca(orca, job.state.path, stem, args.settle)`; when it exited, print its tail lines and `return _error(…)`-style exit 1: stderr `f"orcamon: ORCA exited with status {code} within {args.settle:g} s; the output's last lines are above"` and `return EXIT_JOB_BAD`. When it is running: print `f"started {job.label}: pid {pid}, {orca} {stem}.inp, output {job.out_path}"`; unless `--no-focus`, when `steer.current_view()` is not None and `job.inp_path.resolve()` is inside its `root`, call `_queue_focus(view, job.inp_path.resolve())` and print `f"queued for the TUI (pid {view['pid']})"`. Return 0.
5. JSON instead of text: `{"job", "dry_run": bool, "orca": str, "output": str, "previous_output": str or None, "pid": int or None, "exited": int or None, "tail": [...], "warnings": [{"code", "message"}], "focus_queued": bool}`; a refusal prints nothing on stdout but the `error:` lines, and its message on stderr.
6. `SKILL.md`: in "Changing a geometry", replace `and run it only when they say so.` with `and when they say so, start it with `orcamon run` and follow it with `orcamon wait`.`; add a paragraph after it: "`orcamon run JOB` starts ORCA on this machine, by its full path, in the job's directory, with the output beside the input and detached from your shell. It checks the input first and refuses a job that is already running, has errors, needs `mpirun` that is missing, or already has an output — `--rerun` keeps the old output as `job.out.1` and runs again. `--dry-run` says what it would do. An ORCA that dies in its first seconds is reported at once (exit 1) with the output's last lines. On a cluster's login node it refuses: `run` does not submit to a queue." README: a section `## Running and stopping` after "Changing a geometry" with the same content in the README's voice, the exit-code table's row 1 as in step 1, the `run` reference section and the refreshed first block. Reinstall the skill.
7. `test_orcamon.py`: `test_run_starts_orca_only_when_it_should()`, registered right after `test_the_runner_finds_orca_and_stops_what_it_started`. A temp root with drafts `water`, `fresh`, `bad` (its `* xyz 0 2`), `fails` (`# FAKE_FAIL`) and `wide` (`%pal nprocs 4 end` before the geometry), the fake in `tmp / "orca_bin"`, and `PATH` and `XDG_STATE_HOME` restored in a `finally`. Exactly 9 checks:
   - `run water`: exit 0; line 1 starts `started water: pid `; that pid's `/proc/<pid>/comm` is `orca`; `job.out` contains `fake ORCA running job.inp`;
   - `show water --json`: `doc["status"] == "running"`;
   - `run water` again: exit 2, stderr contains `water is already running`;
   - after `runner.stop_tree(<that pid>, 2.0)`, `run water`: exit 2, stderr contains `job.out exists: --rerun moves it to job.out.1 and runs again`; then `run water --rerun`: exit 0, `job.out.1` contains `fake ORCA running job.inp` and a new `job.out` exists (stop the new pid);
   - `run fresh --dry-run`: exit 0, the last line is `f"would start: {fake.resolve()} job.inp in {root.resolve() / 'fresh'}"`, no `job.out` exists in `fresh`, and `runner.driver_pids` there is `[]`;
   - `run bad`: exit 2, stdout contains `error: 10 electrons cannot have multiplicity 2: it needs an odd number`, no `job.out`;
   - `run fails`: exit 1, stdout contains `ORCA finished by error termination in SCF`, stderr contains `ORCA exited with status 1 within 0.5 s`;
   - `run wide` with `PATH` set to an empty temporary directory: exit 2, stderr contains `nprocs 4 needs mpirun on PATH`, and no `job.out` in `wide`;
   - with `PATH` set to a directory holding the fake `squeue` and the `node01`/`node02` `sinfo` of AR1: `run fresh` WITHOUT `--local`: exit 2, stderr contains `is not one of SLURM's nodes`; `run fresh --local --dry-run`: exit 0.
8. Run the first three verification commands, the `ps` check of AR1, and on the real tree (`cd /home/madsr2d2/masterThesis`, the two environment variables of § Never): `orcamon --root computational run C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/ts5 --dry-run` and `orcamon --root computational run C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts_freq_tight --dry-run`.

**Acceptance:**
- The `ts5` dry run exits 0 and its last line is `would start: /home/madsr2d2/orca_6_1_1/orca job.inp in /home/madsr2d2/masterThesis/computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/ts5`, and `ts5` still holds only its four files; the `optts_freq_tight` dry run exits 2 with `job.out exists: --rerun moves it to job.out.1 and runs again`. Neither started a process (`ps` shows no new `orca`).
- `test_run_starts_orca_only_when_it_should` passes its 9 checks; `test_focus_waits_for_the_tui_to_take_it` and `test_the_tui_takes_the_agents_focus` pass unchanged; no fake process is left.
- `test_orcamon.py`: 299 `pass`; `test_monitor.py`: 266 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### AR3 — `orcamon stop`, the skill's rule, then every gate

**Why:** "stop that one" must end the job's whole process tree and nothing else (D7), and the skill must say when the agent may run or stop at all (D2).
**Depends on:** AR2.
**Files:** `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. `cli/__init__.py`, after `run`: `add("stop", "stop a job's ORCA run on this machine: its whole process tree", "orcamon stop opt/ts", commands.cmd_stop)`, `_job(p)`, `--grace` (`type=float`, `default=10.0`, `metavar="S"`, help `"seconds between SIGTERM and SIGKILL (default: %(default)s)"`), `_output(p, lines=False)`.
2. `@_with_job def cmd_stop(args, job)`: when `job.liveness.source == "slurm"` and `job.liveness.sched_id` → `UsageError(f"{job.label} runs under SLURM as job {id}: stop it with scancel {id}")`; `roots = runner.driver_pids(job.state.path, job.state.stem)`; none → `UsageError(f"{job.label} is not running: no ORCA in {job.state.path} names {job.state.stem}.inp")`; `terminated, killed = runner.stop_tree(roots, args.grace)` (`RunError` → `UsageError`). Text: `f"stopped {job.label}: SIGTERM to {n} process{'es' if n != 1 else ''} ({', '.join(map(str, terminated))})"`, then, when `killed`, `f"SIGKILL to {m} that outlived {args.grace:g} s ({', '.join(map(str, killed))})"`. JSON `{"job", "terminated": [...], "killed": [...]}`. Return 0.
3. `SKILL.md`, "Other rules": replace the bullet `- Never kill, signal, edit or resubmit a running job on your own initiative.` and its two continuation lines with `- Start a job with `orcamon run` and stop one with `orcamon stop` only when the person has asked for that job. Never kill, signal, edit or resubmit a running job any other way. `orcamon stop JOB` ends the job's whole ORCA process tree on this machine — SIGTERM, then SIGKILL after `--grace` seconds — and nothing else; a job SLURM runs is stopped with `scancel`.` (wrapped like its neighbours). README: `stop` in "Running and stopping", its reference section and the refreshed first block. Reinstall the skill.
4. `test_orcamon.py`: `test_stop_ends_the_whole_run()`, registered right after `test_run_starts_orca_only_when_it_should`. Drafts `water` and `stubborn` (`# FAKE_STUBBORN`), the fake, `--orca`, `--settle 0.5`, `--local`, a `finally` that stops anything left. Exactly 4 checks:
   - `run water`, then `stop water --grace 2`: exit 0, line 1 starts `stopped water: SIGTERM to 2 processes (`, and neither the started pid nor its `sleep` child is `runner.alive` afterwards;
   - `stop water` again: exit 2, stderr contains `water is not running: no ORCA in` and `names job.inp`;
   - `run stubborn`, then `stop stubborn --grace 0.5 --json`: exit 0, `len(doc["terminated"]) == 2`, `len(doc["killed"]) == 2`, and nothing in `doc["terminated"]` is alive afterwards;
   - `test_the_skill_names_every_command` holds for `run` and `stop`: `` `orcamon run` `` and `` `orcamon stop` `` both appear in `skill.body()`, and the body no longer contains `Never kill, signal, edit or resubmit a running job on your own initiative.`
5. Run all four verification commands and the `ps` check of AR1.

**Acceptance:**
- `test_stop_ends_the_whole_run` passes its 4 checks; no fake process is left.
- `test_orcamon.py`: 303 `pass`; `test_monitor.py`: 266 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`. `run_gates.py`: `38 gates`, `0 failed`.
- `git diff d232ec5 -- '*test*'` across the whole plan shows only added tests, fixtures, imports and `__main__` entries; no existing assertion or expected value changed.

---

## Out of scope

1. Submitting to SLURM (`orcamon submit`, a job-script template), stopping a SLURM job (`scancel` stays the person's or the agent's own command), and SLURM liveness changes.
2. Restarts that rewrite the input (`MORead`, `%moinp`, reading a `.gbw`): the agent writes the input; `run` starts what is there.
3. Running on another machine (`ssh HOST orcamon run …` works only as far as the remote orcamon does).
4. Queueing several jobs, priorities, or limiting how many run at once; `nice`, CPU pinning, or memory limits beyond `check`'s warnings.
5. Any change to the TUI beyond the focus `run` queues.
6. Platforms without `/proc`: `run` still starts ORCA there, but `stop` and the running check need `/proc` and are not made portable here.

## Progress template

```markdown
# orcamon — the agent starts and stops ORCA when the person says so — Progress

Plan: `PLAN_ORCAMON_ACT_RUN.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| AR1 the runner | TODO | | | |
| AR2 orcamon run | TODO | | | |
| AR3 orcamon stop, the skill's rule, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 266 | 285 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
