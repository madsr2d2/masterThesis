# orcamon — fix the defects the end-to-end review found (FX1–FX5)

**Written:** 2026-10-01. **Verified against:** `master` at `090bae7` (every defect below was reproduced or read on this commit; suite counts and times were measured on it).
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_DEFECTS_PROGRESS.md` (template at the end).

An end-to-end review of orcamon (2026-10-01) found five defects in what the tool already does, before any new feature is built. Liveness is decided per DIRECTORY, so a job that never ran reads `running` when a sibling job in the same directory runs, and orcamon counts its OWN process as a running ORCA. `ls` shortens labels with `…`, so an agent cannot pass them back. The pixel TUI asks for matplotlib, which the `images` extra no longer installs. Any `.inp` is taken for an ORCA job, CREST's included. And a converged optimization shows unmet gradient criteria without saying that ORCA's relaxed rule converged it. **The goal a deviation must still meet:** every status, label and line orcamon prints names the right job and tells the truth about it, and nothing that works today changes behaviour except where a task says so.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 155 `pass` lines, last line `0 failure(s)` | ~1 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_defects/xdg .venv/bin/python test_orcamon.py` | 164 `pass` lines, last line `0 failure(s)` | ~9 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~45 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_defects/xdg .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (same for `test_orcamon.py`, with its environment variables). Each task states the counts it must end at.

Every task runs the first three commands. Only FX5 runs `run_gates.py`: FX1–FX4 touch only `tools/orcamon/` and `test_orcamon.py`, which no other gate reads except the duplicate guard.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's; a new gate file would change the count and break that gate. New tests go in `test_monitor.py` or `test_orcamon.py`.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- a hand-derived expected value in this plan does not hold after the change is implemented as written (report the measured value; do not change the expected value or the design to make it pass);
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance;
- in FX1, `/proc/<pid>/comm` of the fake ORCA process is not `orca`, or its `/proc/<pid>/cmdline` does not contain `opt.inp` (the test's premise is then wrong);
- in FX4, `orcamon --root computational ls --json` lists a number of jobs other than 80, and the difference is not explained by job directories added to `computational/` since `090bae7` (report the labels that differ);
- in FX5, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's. (`tools/orcamon/README.md` is orcamon's and may be edited where a task says so.)
- Never `git add -A`, `git add .` or `git commit -a`; stage the files the task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file.
- Never run `orcamon` against the user's cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_defects/xdg`. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never write into `computational/`, never start a real ORCA, and never kill or signal a process the task did not start itself. Tests use synthetic job trees in temporary directories; the real tree appears only in read-only acceptance commands.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `core` and `cli` import only the standard library at module level.
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill.
- Never change the JSON contract: no key added to or removed from any `--json` document, and `REPORT_SCHEMA` stays `1`.

## Scope

**In scope:** per-job liveness and the self-exclusion (`core/liveness.py` and its four `lookup` call sites); full labels in `ls`; the `images` availability check; content-based input discovery; the convergence reason in the text report; the tests for each; `tools/orcamon/README.md` where a task names it.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | This plan is the defect fixes only. IRC/NEB, all normal modes, and the Act and Steer layers (geometry edits, `check`, `run`/`submit`, shared selection) are later plans. | user, 2026-10-01 |
| D2 | A process is ORCA when its `comm` is exactly `orca` (the driver) or starts with `orca_` (a module such as `orca_scf_mpi`). `orcamon` itself, and any other `orca…` name, is not. | planner, 2026-10-01 |
| D3 | Which job in a directory is running is read from the `orca` driver's command line: its first argument that is non-empty, does not start with `-` and ends in `.inp`, taken as a basename without `.inp`. When at least one driver in a directory names its input, the stems it names are running and the directory's other jobs are NOT running (`Liveness(False, <the entry's source>)`). When no driver in the directory names an input (only modules visible), the old directory-level answer stands. SLURM entries (matched by working directory) stay directory-level. | planner, 2026-10-01 |
| D4 | `ls` prints every label in full, however long the line gets: a label is the handle an agent passes back. The TUI's job table keeps shortening (it has a fixed column). | planner, 2026-10-01 |
| D5 | The pixel graphics need exactly the `images` extra's modules: `numpy` and `PIL`. matplotlib is not checked. | planner, 2026-10-01 |
| D6 | A `.inp` is an ORCA job only when one of its lines, within the first 262144 bytes, has `!`, `%` or `*` as its first non-blank character. An unreadable file is still listed (unchanged behaviour). A `.inp` FILE named explicitly on the command line is obeyed whatever it holds. ORCA's generated inputs (`*.scfgrad.inp`) stay excluded by name, as now. | planner, 2026-10-01 |
| D7 | The convergence reason goes into the TEXT report only (`show`, `wait`, the TUI summary), held on `JobReport` as the private field `_converged_reason`, which `to_dict` already drops. The JSON contract does not change. | planner, 2026-10-01 |
| D8 | Commits go straight to `master`; one commit per task. The workspace was clean at `090bae7` (`git status --porcelain` empty). | user, 2026-10-01 |

## Facts (verified 2026-10-01)

**Environment.** `/home/madsr2d2/masterThesis/.venv/bin/python` is Python 3.12.3 with numpy, PIL, matplotlib, psutil and textual installed. orcamon is installed editable; its console script is `/home/madsr2d2/masterThesis/.venv/bin/orcamon` (not on PATH). `/bin/tail` exists. No ORCA job is running on this machine. `squeue` is installed and the user has no SLURM jobs.

**Liveness** (`tools/orcamon/src/orcamon/core/liveness.py`):
- `running_orca_cwds() -> dict[Path, float]` (line 22) maps each working directory holding an ORCA process to the run's start time; it calls `_proc_entries()` (line 53) where `/proc` exists, else `_psutil_entries()` (line 82). Both return `list[tuple[str, Path, float]]` = (comm name, resolved cwd, start epoch seconds).
- The name filter is `if not name.startswith("orca"):` at line 66 (`/proc`) and `if not name.startswith("orca") or cwd is None:` at line 92 (psutil, which asks `process_iter(["name", "cwd", "create_time"])` at line 88).
- `Liveness` (dataclass, line 108) fields: `alive`, `source`, `since=None`, `queued=False`, `sched_id=None`, `note=None` (line 121 is the last).
- `ProcessProbe.snapshot` (line 143) returns `{cwd: Liveness(True, "process", since=t)}` from `running_orca_cwds()`.
- `SlurmProbe.snapshot` (line 219) adds local processes in a loop at line 249: `for cwd, started in running_orca_cwds().items(): if cwd not in live: live[cwd] = Liveness(True, "process", since=started)`.
- `lookup(snapshot, probe, job_dir)` (line 344) is `snapshot.get(_resolve(job_dir)) or probe.default(job_dir)`. Its callers: `cli/commands.py:99` (`_load_ref`, has `ref.stem`), `cli/commands.py:809` (`cmd_wait`'s `poll`, has `job.state.stem`), `tui/app.py:1283` (`_scan`, has `job.state.stem`), `validate.py:29` (has `ref.stem`).
- `test_orcamon.py:1048` monkeypatches `liveness.running_orca_cwds` with a lambda returning `{path: start}`; the SLURM probe must keep calling `running_orca_cwds()` for that existing test to keep working.

**The two liveness defects, reproduced on `090bae7`:**
- A directory holding `opt.inp`, `opt.out` and `freq.inp`, with a process named `orca` running in it, listed `j/freq` as `running` (it has no output and never ran).
- A directory holding `job.inp` and `job.out` (no termination line), with no ORCA anywhere, listed as `running` when `.venv/bin/orcamon --liveness process --no-cache ls` was run FROM that directory: the console script's own `comm` is `orcamon`, which `startswith("orca")`.
- A copy of `/bin/tail` named `orca`, started as `./orca -f opt.inp` in a directory, shows `comm` = `orca` and `cmdline` = `./orca\0-f\0opt.inp\0` and its `cwd` is that directory.

**ls** (`tools/orcamon/src/orcamon/cli/commands.py::cmd_ls`): line 43 `LS_DEFAULT_COLUMNS = 120`; lines 271–275 read `COLUMNS` and compute `label_width`; line 279 prints `short_label(r.label, label_width)`. `short_label` is imported from `core.discovery` at line 22 and used nowhere else in `commands.py`. On the real tree, `ls` printed `K+H2O2_water-relay_to_KP/…/rc_minus/dlpno-ccsdt` and `show` of that string exited 2; the job's full label is `C8_perhydrate_trap/K+H2O2_water-relay_to_KP/energy_dlpno-ccsdt/rc_minus/dlpno-ccsdt`. The TUI uses `short_label` itself (`tui/app.py`) and is not touched.

**images** (`tools/orcamon/src/orcamon/tui/graphics_probe.py::images_available`, line 45) returns `all(importlib.util.find_spec(m) is not None for m in ("matplotlib", "numpy", "PIL"))`, importing `importlib.util` inside the function. `tools/orcamon/pyproject.toml`'s `images` extra is `["numpy>=1.26", "pillow>=10.1"]`; nothing in the package imports matplotlib (`tui/geometry_render.py` says so in its docstring). `tools/orcamon/README.md` line 31 is the table row `| \`images\` | matplotlib, numpy, pillow | the pixel geometry pane (without it the pane is text) |`. `tui/app.py` line 35 is the comment `# \`geometry_render\` (matplotlib, numpy, PIL -- the \`images\` extra) is imported`.

**Discovery** (`tools/orcamon/src/orcamon/core/discovery.py`): `stems_in(job_dir)` (line 32) returns `sorted(n[:-4] for n in names if n.endswith(".inp") and not is_generated_input(n))` (line 38); `discover(root, exclude)` (line 41) builds stems at line 58 with the same test. The module imports `fnmatch, os, re`, `dataclass` and `Path`. `resolve._from_path` (`core/resolve.py`) accepts any explicit `.inp`/`.out` file and calls `stems_in` for a directory argument. On the real tree, `orcamon --root computational ls --json` lists **85** jobs, of which these **5** have no line starting `!`, `%` or `*` (CREST constraint files beginning `$constrain`): `C8_perhydrate_trap/K+H2O2_water-relay_to_KP/crest_nci_fixhost_apo`, `…/crest_nci_fixhost_cavrim`, `…/crest_nci_fixhost_droplet`, `…/crest_nci_fixhost_droplet/quick`, `fixtures/crest_frozen_host`. So **80** after the fix.

**Convergence reason.** ORCA prints, between a cycle's convergence table and `THE OPTIMIZATION HAS CONVERGED`, one of three reasons when it converges without all five criteria; counted over every `.out` under `computational/`, the line before `Convergence will therefore be signaled now` is reached from:
- `       The gradient convergence is overachieved with ` (35 times; next line `reasonable convergence on the displacements`)
- `       The step convergence is overachieved with ` (4 times; next line `reasonable convergence on the gradient`)
- `       Everything but the energy has converged. However, the energy` (38 times)
When all five criteria are met it prints none of them. `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts_freq_tight/job.out` lines 54953–54973: last table `RMS gradient … NO`, `MAX gradient … NO`, then `The step convergence is overachieved with`, then `THE OPTIMIZATION HAS CONVERGED`. `orcamon show geometry_b973c-xtb/ts/optts_freq_tight` currently prints `finished · cycle 41/400 · 0:09:15 · last output … ago · optimization converged` and `criteria   dE ✓  RMS grad ✗  MAX grad ✗  RMS step ✓  MAX step ✓`.

**Parser** (`tools/orcamon/src/orcamon/core/parser.py`): `_CYCLE_RE` (line 19), `_OPT_DONE_RE` (line 26, matches `OPTIMIZATION HAS CONVERGED|OPTIMIZATION RUN DONE`). `_RARE_MARKERS_RE` (lines 110–118) is the prefilter: a marker whose literal is missing from it is never read; `_LINE_OF_INTEREST_RE` includes `_RARE_MARKERS_RE.pattern`, so adding the literal there covers both paths. `JobState._feed_marker` (line 422) handles each marker; the cycle banner sets `self.cycle` and `self._pending_eig = None` (lines 425–428); `_OPT_DONE_RE` sets `self.opt_converged = True` (lines 445–446). `reset_for_restart` copies every field from a fresh instance, so new fields need no reset code. `core/cache.py` keys entries on a hash of `parser.py`'s source, so a parser edit invalidates every cache entry by itself. `tools/orcamon/src/orcamon/validate.py::MARKER_CASES` lists `(line, predicate)` pairs and `test_monitor.test_every_marker_reaches_its_parser` checks them as ONE check (`validate.MARKER_CASES all read`).

**Report** (`tools/orcamon/src/orcamon/core/report.py`): `JobReport` ends with two private fields, `_input_read` and `_liveness_note`, declared with `field(default=…, repr=False)`; `to_dict` drops every key starting with `_`. `build_report` sets them by keyword. `report_lines` lines 279–280: `if report.opt_converged: segs += [(" · ", None), ("optimization converged", "green")]`. `test_monitor.test_report_keys_are_stable` pins the JSON key set.

**Tests.** `test_monitor.py`: 155 `pass`; `test_orcamon.py`: 164 `pass`; `data/test_curve_metrics.py`: `0 failure(s)`; `run_gates.py`: 38 gates, 0 failed. Both files define `check(label, ok, detail="")`. `test_orcamon.py` has `_orcamon(argv) -> (exit code, stdout, stderr)`, which runs the CLI in-process; `test_monitor.py` has `_feed(state, text)` (feeds each line through `feed_line`) and `_report_job(state, inp_text, liveness, now=None)`. `test_orcamon.py`'s `__main__` list contains, in order, `test_a_job_argument_names_one_job`, `test_ls_lists_every_job_boundedly`, …, `test_the_geometry_pane_degrades_to_text`, …, `test_squeue_lines_parse`, …; `test_monitor.py`'s contains `test_attention_flags`.

## Conventions

- Match the surrounding code: one-line comments where the reason is not obvious; docstrings say WHY.
- New tests go in the file the task names, use `check(...)` for every assertion, build their job trees under `tempfile.TemporaryDirectory()`, and are added to that file's `__main__` list right after the test named in the task.
- Expected values are the hand-derived ones written in this plan, never values copied from a run.
- Any environment variable or module attribute a test changes is restored in a `finally:`.
- No network access.

---

# Phase 1 — name the right job

### FX1 — Liveness belongs to the job, and orcamon is not ORCA

**Why:** Liveness is keyed by directory, so a never-run `freq.inp` beside a running `opt.inp` reads `running`; and `orcamon` run from inside a job directory counts itself as an ORCA process, so a stopped job there reads `running`. Both reproduced on `090bae7` (§ Facts, "The two liveness defects").
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/liveness.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/src/orcamon/validate.py`, `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. In `core/liveness.py` add, above `running_orca_cwds`:
   ```python
   def is_orca_process(name: str) -> bool:
       """The ORCA driver (`orca`) or one of its modules (`orca_scf_mpi`, ...).
       Not `orcamon`, whose console script would otherwise count itself."""
       return name == "orca" or name.startswith("orca_")


   def input_stem(args: list[str]) -> str | None:
       """The job an `orca` driver was started on: its first argument that is
       an input file, as a stem (`/a/b/job.inp` -> `job`); None when none is."""
       for arg in args[1:]:
           if arg and not arg.startswith("-") and arg.endswith(".inp"):
               return Path(arg).name[:-4]
       return None
   ```
2. Make `_proc_entries` and `_psutil_entries` return `list[tuple[str, Path, float, list[str]]]` — the fourth item is the argument vector. In `_proc_entries`, replace the name test at line 66 with `if not is_orca_process(name):` and read `args = open(f"{base}/cmdline", "rb").read().decode("utf-8", "replace").split("\0")` inside the same `try` (a failure to read it gives `args = []`, not a skipped process). In `_psutil_entries`, ask `process_iter(["name", "cwd", "create_time", "cmdline"])`, test `if not is_orca_process(name) or cwd is None:`, and use `proc.info["cmdline"] or []`. Add `def _orca_entries()` returning `_proc_entries()` when `_PROC.is_dir()` else `_psutil_entries()`, and make `running_orca_cwds` iterate `_orca_entries()` (unpacking four items; its return value is unchanged).
3. Add `def running_orca_stems() -> dict[Path, frozenset[str]]`: for each cwd, the set of `input_stem(args)` over its entries whose name is exactly `"orca"` and whose `input_stem` is not None; a cwd with no such entry is left out of the dict. Docstring: a directory left out is one where no driver names its input, so its liveness stays directory-level.
4. Add the field `stems: frozenset | None = None` as the LAST field of `Liveness`, with the comment `# the jobs (stems) known to run in this directory; None = cannot tell which`. In `ProcessProbe.snapshot`, call `stems = running_orca_stems()` and build `Liveness(True, self.name, since=t, stems=stems.get(cwd))`. In `SlurmProbe.snapshot`, call `stems = running_orca_stems()` before the loop at line 249 and build `Liveness(True, "process", since=started, stems=stems.get(cwd))` there; keep that loop calling `running_orca_cwds()`.
5. Change `lookup` to:
   ```python
   def lookup(snapshot, probe, job_dir, stem=None) -> Liveness:
       entry = snapshot.get(_resolve(job_dir))
       if entry is None:
           return probe.default(job_dir)
       # A driver named the job it runs: this directory's other jobs are not it.
       if stem is not None and entry.stems is not None and stem not in entry.stems:
           return Liveness(False, entry.source)
       return entry
   ```
   (keep the existing type annotations). Pass the stem at every call site: `cli/commands.py:99` `lookup(snapshot, probe, ref.path, ref.stem)`; `cli/commands.py:809` `lookup(snapshot, probe, job.state.path, job.state.stem)`; `tui/app.py:1283` `lookup(snapshot, self.probe, job.state.path, job.state.stem)`; `validate.py:29` `lookup(snapshot, probe, ref.path, ref.stem)`.
6. In `tools/orcamon/README.md`, section "Liveness: is the job running?", replace the `process` bullet's second sentence ("A job directory with no `orca*` process is not running.") with: "A job is running when an `orca` driver or an `orca_*` module runs in its directory; where the driver names its input (`orca job.inp`), only that job is, and the directory's other jobs are not."
7. Add `test_liveness_names_the_job_not_only_the_directory()` to `test_orcamon.py`, registered right after `test_squeue_lines_parse`. Exactly 6 `check` calls:
   - `is_orca_process("orca") and is_orca_process("orca_scf_mpi") and not is_orca_process("orcamon") and not is_orca_process("orcabox")`;
   - `input_stem(["./orca", "-f", "opt.inp", ""]) == "opt" and input_stem(["/x/orca", "/a/b/job.inp"]) == "job" and input_stem(["orca_scf_mpi", "job.scfinp", "job"]) is None`;
   - a temporary root holding directory `j` with `opt.inp`, `freq.inp` and `old.inp` (each `"! B97-3c Opt\n* xyz 0 1\nH 0 0 0\n*\n"`), `opt.out` and `old.out` (each `"partial output\n"`). Copy `/bin/tail` to `j/orca` (`shutil.copy`, then `os.chmod(…, 0o755)`) and start it with `subprocess.Popen([str(j / "orca"), "-f", "opt.inp"], cwd=j, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)`; kill and `wait()` it in a `finally:`. While it runs, `_orcamon(["--root", str(root), "--liveness", "process", "--no-cache", "ls", "--json"])`, parse stdout, and map label -> status. Check `j/opt` is `running`; check `j/freq` is `not run`; check `j/old` is `stopped`;
   - a second temporary directory `d` holding `job.inp` (same text) and `job.out` (`"partial output\n"`), no fake ORCA. Run `subprocess.run([str(Path(sys.executable).with_name("orcamon")), "--liveness", "process", "--no-cache", "ls", "--json"], cwd=d, capture_output=True, text=True, timeout=60, env={**os.environ, "XDG_CACHE_HOME": <a temp dir>, "HERDR_ENV": "0"})`, parse stdout; check its one job's status is `stopped`.
8. Run the first three verification commands, and the manual probe: `cd /tmp/orcamon_defects && mkdir -p self && cd self && printf '! SP\n* xyz 0 1\nH 0 0 0\n*\n' > job.inp && printf 'partial output\n' > job.out && XDG_CACHE_HOME=/tmp/orcamon_defects/xdg /home/madsr2d2/masterThesis/.venv/bin/orcamon --liveness process --no-cache ls`.

**If unsure:** If `/proc/<pid>/cmdline` of the fake process does not contain `opt.inp` (for example the copy is a symlink that changes `argv`), stop: the test premise is the same as the real `orca job.inp` and must hold. If `Path(sys.executable).with_name("orcamon")` does not exist, use `/home/madsr2d2/masterThesis/.venv/bin/orcamon` and report it as a deviation.

**Acceptance:**
- The manual probe in step 8 prints one row whose STATUS is `stopped` (it printed `running` on `090bae7`).
- `test_liveness_names_the_job_not_only_the_directory` passes all 6 checks: `j/opt` running, `j/freq` not run, `j/old` stopped, the self-run job stopped.
- `test_orcamon.py`: 170 `pass`, `0 failure(s)` (including the existing `in slurm mode a local ORCA process still counts, quiet or not`). `test_monitor.py`: 155 `pass`, `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.
- `git diff HEAD -- test_orcamon.py` shows only the added test and its `__main__` line.

### FX2 — `ls` prints labels an agent can pass back

**Why:** `ls` shortens labels with `…` (`cli/commands.py:279`), and the shortened text is not a job name: `show "K+H2O2_water-relay_to_KP/…/rc_minus/dlpno-ccsdt"` exits 2 on the real tree. The shipped skill tells agents a JOB argument is a label "as `orcamon ls` prints it".
**Depends on:** FX1.
**Files:** `tools/orcamon/src/orcamon/cli/commands.py`, `test_orcamon.py`.

**Do:**
1. In `cmd_ls`, delete the `COLUMNS` read and the `label_width` computation (lines 271–275) and print `r.label` in place of `short_label(r.label, label_width)`. Add a one-line comment: the label is printed whole because it is the handle the other commands take.
2. Delete `LS_DEFAULT_COLUMNS` (line 43) and remove `short_label` from the `core.discovery` import at line 22 (nothing else in `commands.py` uses either).
3. Add `test_ls_prints_labels_an_agent_can_pass_back()` to `test_orcamon.py`, registered right after `test_ls_lists_every_job_boundedly`. In a temporary root, create `reaction_with_a_long_descriptive_name/geometry_b973c-xtb_with_implicit_solvent/ts/optts_freq_tight_attempt2/job.inp` containing `"! B97-3c Opt\n* xyz 0 1\nH 0 0 0\n*\n"`; its label is that directory path, 107 characters. Set `os.environ["COLUMNS"] = "120"` for the test (restore it in `finally:`). Run `_orcamon(["--root", str(root), "--liveness", "mtime", "--no-cache", "ls"])`. Exactly 3 `check` calls:
   - the full 107-character label appears in the `ls` stdout;
   - `"…"` does not appear in the `ls` stdout;
   - `_orcamon(["--root", str(root), "--liveness", "mtime", "--no-cache", "show", <the full label>])` exits 0.
4. Run the first three verification commands, and `cd /home/madsr2d2/masterThesis && XDG_CACHE_HOME=/tmp/orcamon_defects/xdg .venv/bin/orcamon --root computational show "C8_perhydrate_trap/K+H2O2_water-relay_to_KP/energy_dlpno-ccsdt/rc_minus/dlpno-ccsdt"; echo "exit $?"`.

**Acceptance:**
- `XDG_CACHE_HOME=/tmp/orcamon_defects/xdg .venv/bin/orcamon --root computational ls --max-lines 200 | grep -c '…'` prints `0`.
- The `show` in step 4 prints the job's summary and `exit 0`.
- `test_ls_prints_labels_an_agent_can_pass_back` passes its 3 checks; the existing `test_ls_lists_every_job_boundedly` still passes unchanged.
- `test_orcamon.py`: 173 `pass`, `0 failure(s)`. `test_monitor.py`: 155 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- `tui/app.py` is not in the diff (the TUI keeps `short_label`).

# Phase 2 — the install and the tree

### FX3 — The pixel panes need numpy and Pillow, not matplotlib

**Why:** `images_available()` (`tui/graphics_probe.py:45`) requires matplotlib, which the `images` extra no longer installs (`pyproject.toml`: numpy, pillow), so on a clean `[tui,images]` install the TUI would always fall back to text and say "images extra not installed". It went unseen because this venv has matplotlib for the thesis.
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/tui/graphics_probe.py`, `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. In `images_available`, check `("numpy", "PIL")` only. Rewrite its comment: the modules are the `images` extra's (numpy, Pillow); they are found, not imported, so startup does not pay for importing them.
2. In `tui/app.py` line 35, change `(matplotlib, numpy, PIL -- the \`images\` extra)` to `(numpy, PIL -- the \`images\` extra)`.
3. In `tools/orcamon/README.md` line 31, change the `images` row's middle cell from `matplotlib, numpy, pillow` to `numpy, pillow`.
4. Add `test_the_images_extra_needs_no_matplotlib()` to `test_orcamon.py`, registered right after `test_the_geometry_pane_degrades_to_text`. Replace `importlib.util.find_spec` with a wrapper that returns `None` for the names in a set `hidden` and calls the real function otherwise; restore the real one in `finally:`. Exactly 2 `check` calls:
   - with `hidden = {"matplotlib"}`, `graphics_probe.images_available()` is `True`;
   - with `hidden = {"PIL"}`, it is `False`.
5. Run the first three verification commands.

**Acceptance:**
- `grep -rn matplotlib tools/orcamon/src/orcamon/tui/graphics_probe.py tools/orcamon/src/orcamon/tui/app.py tools/orcamon/README.md` prints nothing.
- `test_the_images_extra_needs_no_matplotlib` passes both checks; `test_the_geometry_pane_degrades_to_text` still passes unchanged.
- `test_orcamon.py`: 175 `pass`, `0 failure(s)`. `test_monitor.py`: 155 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.

### FX4 — Only ORCA inputs are jobs

**Why:** `discover` and `stems_in` take every `.inp` for an ORCA job, so CREST's constraint files (`$constrain … $end`) appear as 5 `not run` jobs on the real tree (§ Facts, "Discovery").
**Depends on:** FX2.
**Files:** `tools/orcamon/src/orcamon/core/discovery.py`, `test_orcamon.py`.

**Do:**
1. In `core/discovery.py` add, after `_GENERATED_INP_RE`:
   ```python
   # How much of an input is read to decide it is ORCA's. The `!` line or the
   # first `%` block is near the top of any real input; this bounds the read.
   ORCA_INPUT_PROBE_BYTES = 256 * 1024
   # An ORCA input has a `!` keyword line, a `%` block or a `*` coordinate line.
   # CREST and xtb inputs (`$constrain ... $end`) have none of them.
   _ORCA_INPUT_LINE_RE = re.compile(rb"^[ \t]*[!%*]", re.M)


   def looks_like_orca_input(path: Path) -> bool:
       """Whether `path` reads as an ORCA input. A file that cannot be read is
       kept: hiding a job is worse than listing one that is not."""
       try:
           with open(path, "rb") as f:
               head = f.read(ORCA_INPUT_PROBE_BYTES)
       except OSError:
           return True
       return bool(_ORCA_INPUT_LINE_RE.search(head))
   ```
2. Apply it after the name tests in both places: in `stems_in` (line 38) add `and looks_like_orca_input(job_dir / n)`; in `discover` (line 58) add `and looks_like_orca_input(here / f)`. Change the module docstring's first paragraph to say a job is a `<stem>.inp` that reads as an ORCA input. Do not change `core/resolve.py`: an explicit `.inp` file argument is still obeyed (D6).
3. Add `test_only_orca_inputs_are_jobs()` to `test_orcamon.py`, registered right after `test_a_job_argument_names_one_job`. In a temporary root create these directories, each with one `job.inp` except `crest`:
   - `orca/job.inp`: `"! B97-3c SP\n* xyz 0 1\nH 0 0 0\n*\n"`
   - `indented/job.inp`: `"   ! B97-3c SP\n"`
   - `blockonly/job.inp`: `"%pal nprocs 4 end\n"`
   - `commented/job.inp`: `"# a comment\n$metadyn\n  atoms: 1-3\n$end\n"`
   - `empty/job.inp`: `""`
   - `crest/fixhost.inp`: `"$constrain\n  atoms: 1-130\n  force constant=0.05\n$end\n"`

   Exactly 3 `check` calls:
   - `[r.label for r in discover(root)] == ["blockonly", "indented", "orca"]`;
   - `stems_in(root / "crest") == []`;
   - `_orcamon(["--root", str(root), "--liveness", "mtime", "--no-cache", "ls", "--json"])` lists exactly 3 jobs.
4. Run the first three verification commands, and `cd /home/madsr2d2/masterThesis && XDG_CACHE_HOME=/tmp/orcamon_defects/xdg .venv/bin/orcamon --root computational ls --json | .venv/bin/python -c "import json,sys; j=json.load(sys.stdin)['jobs']; print(len(j), [x['label'] for x in j if 'crest' in x['label']])"`.

**If unsure:** If an existing test builds a job whose `.inp` has no `!`, `%` or `*` line and now loses it, stop: that is a test whose expected value would change, which § Stop conditions forbids. (None was found on `090bae7`: every fixture input in `test_orcamon.py` and `test_monitor.py` starts with `!`.)

**Acceptance:**
- The real-tree command in step 4 prints `80 []`.
- `test_only_orca_inputs_are_jobs` passes its 3 checks.
- `test_orcamon.py`: 178 `pass`, `0 failure(s)`. `test_monitor.py`: 155 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- `core/resolve.py` is not in the diff.

# Phase 3 — say why it converged

### FX5 — Name ORCA's relaxed convergence rule, then run every gate

**Why:** `show` on `geometry_b973c-xtb/ts/optts_freq_tight` prints `optimization converged` beside `RMS grad ✗  MAX grad ✗`. Both are true: ORCA signalled convergence because the step criteria were overachieved (`job.out` lines 54953–54973). The report should say so, or a person and an agent read the two lines as a contradiction. ORCA prints one of three reasons (§ Facts, "Convergence reason").
**Depends on:** FX4.
**Files:** `tools/orcamon/src/orcamon/core/parser.py`, `tools/orcamon/src/orcamon/core/report.py`, `tools/orcamon/src/orcamon/validate.py`, `test_monitor.py`.

**Do:**
1. In `core/parser.py`, after `_OPT_MAXITER_RE`, add:
   ```python
   # Why ORCA called an optimization converged without all five criteria met,
   # printed between the cycle's convergence table and the HURRAY banner:
   #     The gradient convergence is overachieved with
   #     The step convergence is overachieved with
   #     Everything but the energy has converged. However, the energy
   # All five met prints none of them.
   _CONVERGED_REASON_RE = re.compile(
       r"The (gradient|step) convergence is overachieved|Everything but the energy has converged"
   )
   ```
   and append `r"|convergence is overachieved|Everything but the energy has converged"` to `_RARE_MARKERS_RE`.
2. Add two `JobState` fields: public `opt_converged_reason: str | None = None` (after `opt_maxiter_reached`, with the comment `# ORCA's reason when it converged without all five criteria; None when all five were met`) and private `_pending_converged_reason: str | None = None` (beside `_pending_eig`).
3. In `_feed_marker`: where the cycle banner matches, also set `self._pending_converged_reason = None`. Add a block: when `_CONVERGED_REASON_RE` matches, set `self._pending_converged_reason` to `"gradient overachieved"` if group 1 is `gradient`, `"step overachieved"` if group 1 is `step`, and `"energy nearly converged"` if group 1 is None. Where `_OPT_DONE_RE` matches, also set `self.opt_converged_reason = self._pending_converged_reason`.
4. In `core/report.py`: add `_converged_reason: str | None = field(default=None, repr=False)` after `_liveness_note`; set `_converged_reason=state.opt_converged_reason` in `build_report`; in `report_lines`, after the `("optimization converged", "green")` segment, append `(f" (on ORCA's relaxed rule: {report._converged_reason})", "dim")` when `report._converged_reason` is not None.
5. In `validate.py`, append three rows to `MARKER_CASES`:
   - `("       The gradient convergence is overachieved with ", lambda s: s._pending_converged_reason == "gradient overachieved")`
   - `("       The step convergence is overachieved with ", lambda s: s._pending_converged_reason == "step overachieved")`
   - `("       Everything but the energy has converged. However, the energy", lambda s: s._pending_converged_reason == "energy nearly converged")`
6. Add `test_orca_says_why_it_converged()` to `test_monitor.py`, registered right after `test_attention_flags`. Use the banner lines `B = "         *              GEOMETRY OPTIMIZATION CYCLE   5            *"`, `B6` (the same with `6`) and `H = "      ***        THE OPTIMIZATION HAS CONVERGED     ***"`, each fed with `_feed` into a fresh `JobState(path=Path("/nonexistent"), stem="job")`. Exactly 6 `check` calls:
   - `B`, `"       The step convergence is overachieved with "`, `H` -> `opt_converged_reason == "step overachieved"`;
   - `B`, `"       Everything but the energy has converged. However, the energy"`, `H` -> `"energy nearly converged"`;
   - `B`, `"       The gradient convergence is overachieved with "`, `H` -> `"gradient overachieved"`;
   - `B`, `"       The step convergence is overachieved with "`, `B6`, `H` -> `opt_converged is True and opt_converged_reason is None` (a reason from an earlier cycle is not carried over);
   - for the first state, `render_plain(build_report(_report_job(state, "! B97-3c Opt\n* xyz 0 1\n*\n", Liveness(False, "process"))))` contains `optimization converged (on ORCA's relaxed rule: step overachieved)`;
   - the same report's `to_dict()` has no key containing `reason`.
7. Run all four verification commands, and `cd /home/madsr2d2/masterThesis && XDG_CACHE_HOME=/tmp/orcamon_defects/xdg .venv/bin/orcamon --root computational show geometry_b973c-xtb/ts/optts_freq_tight`.

**Acceptance:**
- The real-tree `show` prints a status line ending `optimization converged (on ORCA's relaxed rule: step overachieved)`, and the `criteria` line is unchanged (`dE ✓  RMS grad ✗  MAX grad ✗  RMS step ✓  MAX step ✓`).
- `test_orca_says_why_it_converged` passes its 6 checks; `validate.MARKER_CASES all read`, `test_report_keys_are_stable` and `test_plain_and_markup_say_the_same_thing` still pass.
- `test_monitor.py`: 161 `pass`, `0 failure(s)`. `test_orcamon.py`: 178 `pass`, `0 failure(s)` (including `test_this_repository_has_the_current_skill`: no command or flag changed, so the installed skill stays current). `data/test_curve_metrics.py`: `0 failure(s)`.
- `run_gates.py`: `38 gates`, `0 failed`.
- `git diff 090bae7 -- '*test*'` across the whole plan shows only added tests, added `__main__` lines and the three `MARKER_CASES` rows; no existing assertion or expected value changed.

---

## Out of scope

1. IRC and NEB support (the IRC path summary, `*_IRC_Full_trj.xyz`, NEB images) — the next plan.
2. Keeping real normal modes as well as imaginary ones.
3. Index labels on environment (non-QM) atoms, atom picking, and any shared selection or focus between the TUI and an agent (the Steer layer).
4. Geometry edits, input checks, connectivity diffs, `run`/`submit`/`stop` (the Act layer), and any change to the shipped skill's rules.
5. Resolving a job by its bare stem when its label is a directory (for example `show job` on a root-level job labelled `.`).
6. A JSON key for the convergence reason (it would bump `REPORT_SCHEMA`).
7. PBS or any scheduler other than SLURM; per-stem liveness for SLURM jobs.
8. Moving orcamon to its own repository, or its tests out of the thesis root.

## Progress template

```markdown
# orcamon — fix the defects the end-to-end review found — Progress

Plan: `PLAN_ORCAMON_DEFECTS.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| FX1 Liveness per job, orcamon is not ORCA | TODO | | | |
| FX2 Full labels in ls | TODO | | | |
| FX3 No matplotlib for the pixel panes | TODO | | | |
| FX4 Only ORCA inputs are jobs | TODO | | | |
| FX5 Convergence reason, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 155 | 164 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
