# orcamon — what the first agent session in the Steer loop got wrong (FS1–FS4)

**Written:** 2026-10-01. **Verified against:** `master` at `a69969f` (the Steer run, `PLAN_ORCAMON_STEER.md`, is committed and complete; the skill description fix is in). The system-size rule was prototyped on a scratch script against this commit and run over all 80 jobs under `computational/`; the `b` bug was reproduced headless on this commit; the fixture values below were measured on those prototypes.
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_STEER_FIXES_PROGRESS.md` (template at the end).

The first real session of the two-pane loop (the person in the TUI, an agent beside it) ran `orcamon view` correctly and still told the person two false things: that the job has "23 atoms, 21 of them in the QM region" (the structure on screen has 137), and that the TUI "will jump to whatever job changes next" (it never changes job by itself). The first is an orcamon bug: for a QM/XTB job whose log prints only the QM block, `show` takes the printed block as the system size. The second is wording: `view` says "following the newest". The agent's `measure` call also failed, because without `--point` it measured the 23-atom printed block, not the 137-atom structure the person picked atoms in. This plan fixes the system size, makes `view` say which numbering its atom numbers belong to and hand over a ready `measure` command, and clears three small defects the Steer run left in its backlog. **The goal a deviation must still meet:** an agent that reads `orcamon show` and `orcamon view` is told the true size of the system and which geometry the picked atom numbers belong to; every existing command, key, JSON key set and test assertion behaves as before; and nothing here starts, edits or stops an ORCA job.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 237 `pass` lines, last line `0 failure(s)` | ~3 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_fix/xdg XDG_STATE_HOME=/tmp/orcamon_fix/state .venv/bin/python test_orcamon.py` | 248 `pass` lines, last line `0 failure(s)` | ~25 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~45 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_fix/xdg XDG_STATE_HOME=/tmp/orcamon_fix/state .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (the same for `test_orcamon.py`, with its environment variables). Each task states the counts it must end at.

Every task runs the first three commands. Only FS4 runs `run_gates.py`: every task touches only `tools/orcamon/`, `test_monitor.py`, `test_orcamon.py` and `.claude/skills/orcamon/SKILL.md`, which no other gate reads except the duplicate guard.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's; a new gate file would change the count and break that gate. New tests go in `test_monitor.py` or `test_orcamon.py`.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- a hand-derived expected value in this plan does not hold after the change is implemented as written (report the measured value; do not change the expected value or the design to make it pass);
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance — in particular `test_the_tui_tells_the_agent_what_it_shows`, `test_the_tui_takes_the_agents_focus`, `test_the_steer_files_name_one_live_tui`, `test_report_keys_are_stable` and `test_a_job_that_prints_no_geometry_shows_the_one_it_wrote` pin today's behaviour and must pass unchanged;
- a real-tree value in FS1 differs from this plan's (report it, and say whether that job's files changed since `a69969f`);
- in FS4, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's. (`tools/orcamon/README.md` and `tools/orcamon/src/orcamon/SKILL.md` are orcamon's; `.claude/skills/orcamon/SKILL.md` changes ONLY through `orcamon skill install --project .`.)
- Never `git add -A`, `git add .` or `git commit -a`; stage the files the task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file.
- Never define the same top-level name in both `test_monitor.py` and `test_orcamon.py`, private or not: the duplicate guard reads private names too. A fixture both files need lives in `test_monitor.py` and `test_orcamon.py` imports it.
- Never let a test or a manual command read or write the user's own state or cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_fix/xdg` and `XDG_STATE_HOME=/tmp/orcamon_fix/state`; a test that touches the shared state files sets `XDG_STATE_HOME` to its own temporary directory and restores it, as `test_the_tui_tells_the_agent_what_it_shows` does. The user's own TUI is running on this machine: a test or a manual command that wrote a `focus` request into `~/.local/state/orcamon` would steer it. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never start the TUI against `computational/`, never write into `computational/`, never start ORCA, and never kill or signal a process the task did not start. The real tree appears only in FS1's read-only acceptance commands.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `core` and `cli` import only the standard library at module level (`test_the_core_needs_only_the_standard_library` enforces it).
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill.
- Never change the key set of an existing `--json` document, and never change `REPORT_SCHEMA` or `steer.STEER_SCHEMA` (both stay `1`). `n_atoms`' VALUE changes (D1); no key is added to or removed from the report. The view document gains three top-level keys (D2) and changes none of its existing ones.
- Never add a fenced code block containing an `orcamon …` line to `tools/orcamon/src/orcamon/SKILL.md`: `test_orcamon.test_the_skill_examples_run` executes every such line, and `view` exits 2 when no TUI runs. Mention commands in prose, in backticks.
- Never change what `i`, `I`, `s`, the chart keys, `focus`, `geom`, `snapshot` or `freqs` do.

## Scope

**In scope:** the system size in `show` and the report; the view document's three new keys, `orcamon view`'s new lines and wording, and the skill and README text for them; `b` clearing a mode the agent started; the shadowed `events` import; a test for a followed path point's address.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | `JobReport.n_atoms` is the size of the WHOLE system: the largest of every count the job offers, each of which is a lower bound — the atoms the log printed, the highest QM atom index + 1, the input's inline atoms, the first line of the input's `*xyzfile`, and the atoms of a file geometry. The key set does not change; `n_atoms`' value does (it was the printed block). | planner, 2026-10-01 |
| D2 | The view document gains three top-level keys and changes none: `atoms_shown` (the atom count of the geometry the pane draws, or None without a job), `selection_atoms` (`[{"atom": i, "element": "O" or None}, …]` for the picked atoms, `[]` when none), and `as_measure` (the `orcamon measure` command for the picked atoms on that geometry, or None unless 2–4 atoms are picked and all exist in it). `orcamon view` prints, in this order, `geometry`, `atoms`, `mode` (when set), `selection`, `elements`, `as measure`, `as XYZ`; the new lines appear only when the document has the keys, and `as XYZ` stays the last line. The `geometry` line's text for a following pane becomes `(following the job)`; a scrubbed one stays `(scrubbed)`. | planner, 2026-10-01 |
| D3 | `b` restores exactly what was on screen before the agent's last `focus`: when no mode was running then, the mode the agent started is stopped. | planner, 2026-10-01 |
| D4 | Commits go straight to `master`, one per task. | user, 2026-10-01 |

## Facts (verified 2026-10-01)

**Environment.** `/home/madsr2d2/masterThesis/.venv/bin/python` is Python 3.12.3 with Textual 8.2.8; orcamon is installed editable and its script is on PATH as `orcamon` (a symlink in `~/.local/bin` to the venv). `git status --porcelain` was empty at `a69969f`. Baseline: `test_monitor.py` 237 pass, `test_orcamon.py` 248 pass (both re-run for this plan), `run_gates.py` 38 gates, 0 failed (the Steer run's progress file). README `--help` blocks are pasted at `COLUMNS=90`.

**The system-size bug** (`tools/orcamon/src/orcamon/core/report.py:171`): `n_atoms=(len(state.atoms) if state.atoms else len(job.file_geometry.atoms) if getattr(job, "file_geometry", None) else None)`. `state.atoms` is the last geometry the LOG printed; a QM/XTB job's log prints only the QM block. `render_plain`/`render_markup` print `QM region … (N_qm of N atoms)` (line 267) for a multilayer job and `· N atoms` otherwise (line 273); `n_qm_atoms` is `len(state.qm_atom_indices)`. `core/job.py`: `Job.__init__` holds `state`, `input` (a `JobInput` with `coords_file`, `coords_atoms`, `multilayer`), `file_geometry` (set only while nothing is printed) and `inp_path`; `refresh` sets them. `core/geometry.py` has `read_xyz(path)` (the whole first frame). `test_monitor.py`'s `_report_job(state, inp_text, liveness)` builds a `Job` WITHOUT calling `refresh` (it sets `job.state`, `liveness`, `parsed`, `input` by hand), so the new count must be a method that works from the current state and input, not a field `refresh` fills. `state.path` is the job directory; `state.qm_atom_indices` comes from the log's `QM1 Subsystem ... 0 1 2` line (regex `^QM1 Subsystem\s*\.\.\.\s*(.*)$` on the stripped line).

**The fixture** (measured on the prototype), in a temp dir `d`: `job.inp` = `"! QM/XTB r2SCAN-3c Opt PAL8\n%maxcore 4000\n%qmmm\n  QMAtoms {0:2} end\n  Charge_Total 0\n  Mult_Total 1\nend\n* xyzfile 0 1 start.xyz\n"`; `job.out` = `"QM1 Subsystem                    ...  0 1 2\nCARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------\n  O      0.000000    0.000000    0.000000\n  H      0.960000    0.000000    0.000000\n  H     -0.240000    0.930000    0.000000\n\n                             ****ORCA TERMINATED NORMALLY****\n"`; `start.xyz` = `"8\nc\n"` + eight lines `C {i * 1.5} 0.0 0.0` for `i` in 0–7. Parsed with `Job(d, "job", d.parent)` and `refresh(Liveness(False, "process"))`: `len(state.atoms) == 3`, `qm_atom_indices == {0, 1, 2}`, `input.coords_file == "start.xyz"`, `input.multilayer` True, `file_geometry` None. **Today** `build_report(job)` gives `n_atoms == 3`, `n_qm_atoms == 3` and the summary line `system charge 0 · singlet   QM region charge 0 · singlet (3 of 3 atoms)   8 procs, 4000 MB/core`. **After FS1** it must give `n_atoms == 8` and `(3 of 8 atoms)` in `render_plain`.

**The real tree** (`computational/`, 80 jobs, read-only; the prototype took the largest of the five counts of D1): 14 jobs change, every one to a LARGER number, all of them QM/XTB jobs whose log prints only the QM block: `…/geometry_b973c-xtb/irc_bridge` 23 → 137 (its `n_qm_atoms` stays 21); `…/geometry_r2scan3c-xtb/ts/verify_exact_hessian` 22 → 140; `…/geometry_b973c-xtb/neb_anion` None → 133; the three `energy_qm1_dlpno` jobs 22 → 140; the five `energy_qm1_dlpno_b973c` jobs 19 → 133 or 137; two `numfreq` jobs 19 → 137. The other 66 do not change, among them `…/geometry_b973c-xtb/ts/optts_freq_tight` (137, `n_qm_atoms` 17). Full labels start `C8_perhydrate_trap/K+H2O2_water-relay_to_KP/`.

**The `b` bug**, reproduced headless on `a69969f` (job `modes_all` of the modes plan: `job.inp` = `_OPT_FREQ`, `job.out` = `_opt(2, [-80.0, -80.1]) + _CONVERGED + "\n" + _freqs(_MODE_FREQS) + _normal_modes(_MODE_VECTORS) + _DONE + "\n"`): with no mode running, `orcamon focus modes_all --atoms 1,0 --mode 5 --timeout 0` makes the pane show mode 5 and selection `(1, 0)`; after `b` the selection is `()` but `mode_index` is STILL 5 and the title still says `mode 5 50.0 cm-1 (real 1/5)`. With the person's own mode 3 running first (`i`), a focus on mode 5 and then `b` correctly returns to `mode 3 -150.0 cm-1 (1/1)`. In `tui/app.py`, `MonitorApp._check_focus` records `doc["_back"] = (self.selected_label, geometry.selection, geometry.mode_index)`; `action_agent_back` (`label, selection, mode = self._back`) calls `geometry.set_selection(selection)` and `geometry.show_mode(mode)` only `if mode is not None`. `RotatableGeometryImage._stop_modes()` stops the timer and clears `mode_index`, `mode_label` and `_mode_s`; the off branch of `_cycle_through` does `self._stop_modes()`, `self._mode_alternates = []`, then `self._input.request()` and `self.app.update_detail()`.

**The view document and command.** `MonitorApp._publish_view(job, point)` (`tui/app.py`) builds `{"schema", "pid", "host", "root", "job", "geometry", "mode", "selection", "measurement", "as_xyz", "focus"}`; for a job it sets `atoms, shown = geometry_shown(job, point)`, `following = chart.following`, `args = steer.geometry_address(shown, following)`, `file = str(job.inp_path)`, `"job": {"label", "file"}`, `"geometry": {"shown", "following", "args"}`, `"as_xyz": shlex.join(["orcamon", "geom", file, *args])`, `"selection": list(geometry.selection)` and, for 2–4 picked atoms, `"measurement"` from `measure(atoms, geometry.selection)`. `cli/commands.py::cmd_view` prints `TUI pid … `, then `job`, `file`, `geometry   {shown} ({'following the newest' if following else 'scrubbed'})`, `mode` when set, `selection` (the measurement text, else the numbers, else `none`) and `as XYZ` — each label left-justified to 11 columns. **Today's tests pin** (and must keep passing unchanged): `geometry   cycle 1 (scrubbed)` and `selection  angle 1-0-2: 104.5°` are lines of `orcamon view`'s output, its LAST line is the `as XYZ     …` line, and the document's `geometry` value is exactly `{"shown": "cycle 2", "following": True, "args": []}` — which is why D2 adds top-level keys and does not touch `geometry`.

**The shadowed import** (`tui/app.py`): line 12 `from textual import events, work` and line 24 `from ..core.events import events` — the second rebinds `events` to orcamon's own function; line 846 `def on_click(self, event: events.Click)` works only because `from __future__ import annotations` never evaluates the hint. `events(…)` (the function) is called once, in `_scan`.

**`geometry_address`** (`core/steer.py`): `["--point", str(shown.index)]` for a `paths.PathPoint` (tested before the `following` guard), `[]` when `following` is true for anything else, when `shown` is None or a `geometry.FileGeometry`, else `--step`/`--cycle`. `test_the_steer_files_name_one_live_tui` covers a PathPoint with `following=False` and a GeometryPoint with both; no test covers a PathPoint with `following=True` (the Steer run's backlog). `PathPoint(kind="irc", index=0, label="IRC TS", energy=None, de_kj_mol=0.0, atoms=[])`; `GeometryPoint(scan_step=None, cycle=4)`; `FileGeometry([], "input geometry")`.

**Test helpers.** `test_orcamon.py`: `check`, `_orcamon(argv)`, `_json(argv)`, `_water_job(root, name="water")` (3-atom water, cycles 1 and 2; cycle 2 is `O (0.2, 0, 0)`, `H (1.16, 0, 0)`, `H (-0.04, 0.93, 0)`), `_opt`, `_CONVERGED`, `_DONE`, `_OPT_FREQ`, `_freqs`, `_normal_modes`, `from test_monitor import _MODE_FREQS, _MODE_VECTORS` inside tests; the headless pattern with a private `XDG_STATE_HOME` is `test_the_tui_tells_the_agent_what_it_shows` (water job, `s`, `*"1 0 2"`, `enter`, reading `steer.read_json(steer.view_path())`) and `test_the_tui_takes_the_agents_focus` (`cli = lambda *a: _orcamon(["--root", str(root), "--liveness", "mtime", "--no-cache", *a])`, an `until(cond, seconds)` poller). `test_monitor.py`: `check`, `Job`, `Liveness`, `build_report`, `render_plain` are imported at the top. Counts: `test_monitor.py` 237, `test_orcamon.py` 248.

## Conventions

- Match the surrounding code: one-line comments where the reason is not obvious; docstrings say WHY.
- New tests go in the file the task names, use `check(...)` for every assertion, build trees under `tempfile.TemporaryDirectory()`, and are added to that file's `__main__` list right after the test named in the task.
- Expected values are the hand-derived ones in this plan, never values copied from a run.
- A README `--help` block changes only if a command's flags change; none do here.
- No network access.

---

# Phase 1 — the system size

### FS1 — `show` counts the whole system, not the block the log printed

**Why:** `orcamon show` printed `QM region … (21 of 23 atoms)` for a 137-atom job and an agent repeated it to the person (§ Facts, "The system-size bug"); 14 of the 80 jobs in the real tree are misdescribed the same way.
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/geometry.py`, `tools/orcamon/src/orcamon/core/job.py`, `tools/orcamon/src/orcamon/core/report.py`, `test_monitor.py`.

**Do:**
1. `core/geometry.py`: add `xyz_atom_count(path) -> int | None` beside `read_xyz`: open the file with `errors="replace"`, read ONLY its first line, return `int(line.split()[0])` when that is at least 1; `None` for a missing or unreadable file, a first line that is not an integer, or a count below 1 (`OSError`, `ValueError`, `IndexError`). Docstring: the atom count an XYZ file announces, without reading its atoms — a 137-atom trajectory's count costs one line.
2. `core/job.py`: `from .geometry import FileGeometry, read_xyz, xyz_atom_count`. In `Job.__init__`, `self._xyz_count: tuple | None = None` with a comment: `((path, mtime_ns, size), count)` of the input's `*xyzfile`, so a re-scan does not reopen an unchanged file. Add the method `system_atoms(self) -> int | None` after `_update_file_geometry`, with a docstring giving D1 in prose (why: ORCA's log may print only the QM block of a multilayer job; every source below is a lower bound on the system's size, so the largest is taken). It collects: `len(self.state.atoms)`; `max(self.state.qm_atom_indices) + 1` when that set is non-empty; `len(self.file_geometry.atoms)` when `file_geometry` is not None; and, when `self.input` is not None, `len(self.input.coords_atoms)` and — when `self.input.coords_file` is set — the input xyz's count: `path = self.state.path / self.input.coords_file`, `st = path.stat()` (an `OSError` gives no count), `stamp = (str(path), st.st_mtime_ns, st.st_size)`, reuse `self._xyz_count[1]` when `self._xyz_count` is not None and its stamp equals `stamp`, else `count = xyz_atom_count(path)` and `self._xyz_count = (stamp, count)`; a None count is left out. Return the largest collected count, or `None` when it is 0 or nothing was collected. Reading nothing but the first line and one `stat` per job keeps `ls` cheap.
3. `core/report.py`: replace the `n_atoms=(…)` expression (line 171) with `n_atoms=job.system_atoms()`. `n_qm_atoms` and every other field are unchanged.
4. `test_monitor.py`: add `test_the_summary_counts_the_whole_system()`, registered right after `test_a_job_that_prints_no_geometry_shows_the_one_it_wrote`. Build the § Facts fixture in a `tempfile.TemporaryDirectory()` (`from orcamon.core.geometry import xyz_atom_count` inside the test; a local `xyz(path, n)` helper writing `f"{n}\nc\n"` + `n` lines `C {i * 1.5} 0.0 0.0`), `gone = Liveness(False, "process")`, `job = Job(d, "job", d.parent)`, `job.refresh(gone)`, `report = build_report(job)`. Exactly 6 checks:
   - `report.n_atoms == 8`;
   - `report.n_qm_atoms == 3`;
   - `"(3 of 8 atoms)" in render_plain(report)`;
   - after `xyz(d / "start.xyz", 10)` and `job.refresh(gone)`: `build_report(job).n_atoms == 10` (the unchanged-file cache does not serve a rewritten file);
   - `xyz_atom_count(d / "start.xyz") == 10`, `xyz_atom_count(d / "missing.xyz") is None`, and a file containing `hello\n` and a file containing `0\n` both give None;
   - a second job (`d2`, `job.inp` = `"! B97-3c Opt\n* xyz 0 1\nO 0.0 0.0 0.0\n*\n"`, the same `job.out` without its first line): `build_report(job2).n_atoms == 3` (a log that prints the whole system, larger than the input's one inline atom, still wins).
5. Run the first three verification commands, and on the real tree (`cd /home/madsr2d2/masterThesis`, `XDG_CACHE_HOME=/tmp/orcamon_fix/xdg XDG_STATE_HOME=/tmp/orcamon_fix/state`) for each label `L` in `geometry_b973c-xtb/irc_bridge`, `geometry_r2scan3c-xtb/ts/verify_exact_hessian`, `geometry_b973c-xtb/neb_anion`, `geometry_b973c-xtb/ts/optts_freq_tight`: `.venv/bin/orcamon --root computational show C8_perhydrate_trap/K+H2O2_water-relay_to_KP/L --json | .venv/bin/python -c "import json,sys; d=json.load(sys.stdin); print(d['n_atoms'], d['n_qm_atoms'])"`, and `.venv/bin/orcamon --root computational show C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/irc_bridge | grep "QM region"`.

**Acceptance:**
- The four real-tree commands print `137 21`, `140 …`, `133 …` and `137 17` (the `…` are `n_qm_atoms`, unchanged by this task; report them), and the `grep` line contains `(21 of 137 atoms)`.
- `test_the_summary_counts_the_whole_system` passes its 6 checks; `test_report_keys_are_stable`, `test_a_job_that_prints_no_geometry_shows_the_one_it_wrote` and `test_plain_and_markup_say_the_same_thing` pass unchanged.
- `test_monitor.py`: 243 `pass`; `test_orcamon.py`: 248 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

# Phase 2 — `view` says what its numbers mean

### FS2 — `view` gives the atom count, the elements and a ready `measure` command

**Why:** the agent measured the picked atoms on the wrong geometry (23 printed atoms, not the 137 the numbers belong to), had to run `measure` just to learn the elements, and read "following the newest" as "will change job" (§ Facts, "The view document and command"). D2.
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. `tui/app.py::_publish_view`: add `"atoms_shown": None`, `"selection_atoms": []` and `"as_measure": None` to the initial document, right after `"as_xyz"` (so `"focus"` stays last). In the `job is not None` branch: `doc["atoms_shown"] = len(atoms)`; `doc["selection_atoms"] = [{"atom": i, "element": atoms[i][0] if 0 <= i < len(atoms) else None} for i in geometry.selection]`; and when `2 <= len(geometry.selection) <= 4` and every picked index is inside `range(len(atoms))`, `doc["as_measure"] = shlex.join(["orcamon", "measure", file, *(str(i) for i in geometry.selection), *args])`. Nothing else in the function changes (the existing `geometry`, `selection`, `measurement` and `as_xyz` values are exactly as today).
2. `cli/commands.py::cmd_view`: change the `geometry` line's text to `'following the job'` where it said `'following the newest'` (the `scrubbed` text stays). After that line, when `doc.get("atoms_shown")` is truthy, print `f"{'atoms':<11}{doc['atoms_shown']} in this geometry"`. After the `selection` line, when `doc.get("selection_atoms")` is non-empty, print `f"{'elements':<11}" + "-".join(a.get("element") or "?" for a in doc["selection_atoms"])`; then, when `doc.get("as_measure")`, print `f"{'as measure':<11}{doc['as_measure']}"`; the `as XYZ` line stays last. A document without the new keys (an older TUI still running) prints exactly as before. Update `cmd_view`'s docstring to mention the atom count, the elements and the measure command.
3. `tools/orcamon/src/orcamon/SKILL.md`, "Working beside a person": replace the first paragraph's second sentence onward (from "It says which job…" to "…do not guess.") with: "It says which job and which geometry their TUI shows and how many atoms that geometry has, the atoms they picked — ORCA's 0-based numbers, as the TUI labels them — with their elements and the distance, angle or dihedral between them, and two ready commands: `as measure` for the picked atoms and `as XYZ` for the geometry. An atom number means something only in the geometry it was picked in: a QM/MM job's printed geometry can hold fewer atoms than the system, so keep the `--point`, `--step` or `--cycle` those commands carry when you run `orcamon measure`, `orcamon geom` or `orcamon snapshot` yourself. "Following the job" means their pane follows the job's own newest geometry (an IRC's TS); it never switches jobs by itself. `orcamon view` exits 2 when no TUI is running for this user on this machine: then ask which job and atoms they mean; do not guess." In "Other rules", after the `orcamon tui` bullet, add a bullet: "Only `orcamon tui` needs a terminal; every other command runs as it is from an agent's shell, with no special environment." README, "Working beside an agent": after the sentence ending "…removes its own file.", add "The view carries the number of atoms in the geometry shown, the picked atoms' elements and ready `orcamon measure` and `orcamon geom` commands for exactly that geometry." Then `cd /home/madsr2d2/masterThesis && .venv/bin/orcamon skill install --project /home/madsr2d2/masterThesis` (with `--force` only if it refuses).
4. `test_orcamon.py`: `test_view_names_the_geometry_its_numbers_belong_to()`, registered right after `test_the_tui_tells_the_agent_what_it_shows`. No Textual. Own temporary `XDG_STATE_HOME` set and restored in a `finally`; `from orcamon.core import steer`; `file = "/r/irc/job.inp"`; `base` = `{"schema": 1, "pid": os.getpid(), "host": steer.host_name(), "root": "/r", "updated_unix_s": time.time(), "job": {"label": "irc", "file": file}, "geometry": {"shown": "IRC TS", "following": True, "args": ["--point", "0"]}, "atoms_shown": 137, "mode": None, "selection": [133, 134, 135], "selection_atoms": [{"atom": 133, "element": "O"}, {"atom": 134, "element": "O"}, {"atom": 135, "element": "H"}], "measurement": {"kind": "angle", "value": 98.1764, "unit": "degree"}, "as_measure": f"orcamon measure {file} 133 134 135 --point 0", "as_xyz": f"orcamon geom {file} --point 0", "focus": None}`; `view()` = `_orcamon(["view"])` with `lines = out.strip().split("\n")` after `steer.write_json(steer.view_path(), doc)`. Exactly 6 checks:
   - with `base`: exit 0 and `"geometry   IRC TS (following the job)" in lines`;
   - `"atoms      137 in this geometry" in lines`;
   - `"selection  angle 133-134-135: 98.2°"` is in `lines` and the very next line is `"elements   O-O-H"`;
   - `lines[-2] == "as measure orcamon measure /r/irc/job.inp 133 134 135 --point 0"` and `lines[-1] == "as XYZ     orcamon geom /r/irc/job.inp --point 0"`;
   - with `following` False: `"geometry   IRC TS (scrubbed)" in lines`;
   - with `atoms_shown`, `selection_atoms` and `as_measure` removed from the document: exit 0, no line starts with `atoms `, `elements ` or `as measure`, and `lines[-1]` is still the `as XYZ` line.
5. `test_orcamon.py`: `test_the_tui_publishes_the_numbering_it_shows()`, registered right after `test_view_names_the_geometry_its_numbers_belong_to`. Own temporary `XDG_STATE_HOME` and tree holding `_water_job(root, "water")`; `file = str((root / "water" / "job.inp").resolve())`; headless as `test_the_tui_tells_the_agent_what_it_shows` (same wait until the job is parsed and shown; `geometry.focus()`; `s`, `*"1 0 2"`, `enter`, `await pilot.pause(0.2)`; `doc = lambda: steer.read_json(steer.view_path())`). Exactly 5 checks:
   - `doc()["atoms_shown"] == 3`;
   - `doc()["selection_atoms"] == [{"atom": 1, "element": "H"}, {"atom": 0, "element": "O"}, {"atom": 2, "element": "H"}]`;
   - `doc()["as_measure"] == shlex.join(["orcamon", "measure", file, "1", "0", "2"])`;
   - after `app.query_one("#chart").focus()`, `left`, `await pilot.pause(0.2)`: `doc()["as_measure"] == shlex.join(["orcamon", "measure", file, "1", "0", "2", "--cycle", "1"])`;
   - after `geometry.focus()`, `s`, five `backspace`, `enter` (clearing it): `doc()["as_measure"] is None and doc()["selection_atoms"] == []`.
6. Run the first three verification commands.

**Acceptance:**
- `test_view_names_the_geometry_its_numbers_belong_to` passes its 6 checks and `test_the_tui_publishes_the_numbering_it_shows` its 5; `test_the_tui_tells_the_agent_what_it_shows`, `test_the_skill_names_every_command`, `test_the_skill_examples_run`, `test_the_skill_is_general` and `test_this_repository_has_the_current_skill` pass unchanged.
- `test_orcamon.py`: 259 `pass`; `test_monitor.py`: 243 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

# Phase 3 — small defects from the Steer run

### FS3 — `b` stops a mode the agent started

**Why:** `b` is promised to return the person to what was on screen before the agent's last `focus` (D3); today it leaves an agent-started mode animating (§ Facts, "The `b` bug").
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `test_orcamon.py`.

**Do:**
1. `RotatableGeometryImage`: add `clear_mode(self) -> None` right after `show_mode`: docstring "Stop animating a mode, as the `i` and `I` cycles do after their last one"; body `self._stop_modes()`, `self._mode_alternates = []`, `self._input.request()`, `self.app.update_detail()`.
2. `MonitorApp.action_agent_back`: its last statement is `if mode is not None: geometry.show_mode(mode)`; give that `if` an `else: geometry.clear_mode()`. Nothing else changes (one level, the `job is being read` retry, `set_selection`).
3. `test_orcamon.py`: `test_b_does_not_leave_the_agents_mode_running()`, registered right after `test_the_tui_takes_the_agents_focus`. Own temporary `XDG_STATE_HOME` as in `test_the_tui_takes_the_agents_focus`; a tree holding only `modes_all` (`job.inp` = `_OPT_FREQ`, `job.out` as in § Facts "The `b` bug"; import `_MODE_FREQS`, `_MODE_VECTORS` from `test_monitor`); headless with `MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")`, the same wait until `geometry._job` is `modes_all`, and that test's `cli` and `until` helpers. Exactly 4 checks:
   - `cli("focus", "modes_all", "--atoms", "1,0", "--mode", "5", "--timeout", "0")` exits 0 and within 3 s `geometry.mode_index == 5 and geometry.selection == (1, 0)`;
   - after `b` and `await pilot.pause(0.4)`: `geometry.mode_index is None`, `geometry.selection == ()` and `"cm-1" not in geometry.border_title`;
   - after `geometry.focus()`, `i`, `await pilot.pause(0.2)` (the person's own mode: `geometry.mode_index == 3`), then `cli("focus", "modes_all", "--mode", "5", "--timeout", "0")` exits 0 and within 3 s `geometry.mode_index == 5`;
   - after `b` and `await pilot.pause(0.4)`: `geometry.mode_index == 3` and `"mode 3 -150.0 cm-1 (1/1)" in geometry.border_title`.
4. Run the first three verification commands.

**Acceptance:**
- `test_b_does_not_leave_the_agents_mode_running` passes its 4 checks; `test_the_tui_takes_the_agents_focus` and `test_the_tui_cycles_real_modes` pass unchanged.
- `test_orcamon.py`: 263 `pass`; `test_monitor.py`: 243 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### FS4 — Stop shadowing Textual's `events`, cover a followed path point's address, then every gate

**Why:** two backlog items from the Steer run: `events` is bound to Textual's module and then to orcamon's function in the same file, which works only because a type hint is never evaluated; and `geometry_address` has no test for a path point while the pane follows the path.
**Depends on:** FS3 (same file).
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `test_monitor.py`.

**Do:**
1. `tui/app.py`: line 12 becomes `from textual import work`; add `from textual.events import Click` beside the other `textual` imports (alphabetical order among them: after `textual.css.query`); line 846's annotation becomes `event: Click`. `from ..core.events import events` and its one call stay exactly as they are. No other change.
2. `test_monitor.py`: add `test_a_followed_path_point_is_addressed()`, registered right after `test_the_steer_files_name_one_live_tui`. Imports inside the test as that test's: `steer`, `FileGeometry`, `GeometryPoint`, `PathPoint`. Exactly 3 checks:
   - `steer.geometry_address(PathPoint(kind="irc", index=0, label="IRC TS", energy=None, de_kj_mol=0.0, atoms=[]), True) == ["--point", "0"]`;
   - `steer.geometry_address(GeometryPoint(scan_step=None, cycle=4), True) == []`;
   - `steer.geometry_address(FileGeometry([], "input geometry"), False) == []`.
3. Run all four verification commands.

**Acceptance:**
- `grep -n "events" tools/orcamon/src/orcamon/tui/app.py` shows no `events.Click` and no `from textual import events`; `from textual import work` is the only `from textual import` line.
- `test_a_followed_path_point_is_addressed` passes its 3 checks.
- `test_monitor.py`: 246 `pass`; `test_orcamon.py`: 263 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`. `run_gates.py`: `38 gates`, `0 failed`.
- `git diff a69969f -- '*test*'` across the whole plan shows only added tests, imports and `__main__` entries; no existing assertion or expected value changed.

---

## Out of scope

1. Restoring the chart's scrubbed point on `b`: `b` returns the job, the atoms and the mode, but a person who had scrubbed an optimization's chart to cycle 5 is back on the job's newest geometry. A real gap in "back to where you were", for the user to decide on.
2. An IRC job whose first forward point is 21.8 kJ/mol above its TS point (`irc_bridge`), or any chemistry judgement about a job; orcamon reports, and flags only what the existing attention flags say.
3. Making `measure`, `geom` or `snapshot` default to the path's geometry on an IRC or NEB job; the `as measure` and `as XYZ` lines carry the right flags (D2).
4. Everything that acts on ORCA (the Act plan): geometry edits, `check`, `connectivity`, `run`, `stop`.
5. Deleting an old view file's keys, bumping `STEER_SCHEMA`, or changing `REPORT_SCHEMA`.

## Progress template

```markdown
# orcamon — what the first agent session in the Steer loop got wrong — Progress

Plan: `PLAN_ORCAMON_STEER_FIXES.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| FS1 show counts the whole system | TODO | | | |
| FS2 view gives atoms, elements, measure | TODO | | | |
| FS3 b stops the agent's mode | TODO | | | |
| FS4 events import, path address test, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 237 | 248 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
