# orcamon — the person and the agent share one view (ST1–ST8)

**Written:** 2026-10-01. **Verified against:** `master` at `91cb1d9` (the normal-modes run, `PLAN_ORCAMON_MODES.md`, is committed and complete). The measurement, projection, picking and selection-ring rules below were prototyped on a scratch copy against this commit's renderer; every hand-derived value in § Facts was measured on that prototype, and the Textual behaviours were probed on the installed Textual 8.2.8.
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_STEER_PROGRESS.md` (template at the end).

orcamon is meant to be used in two panes over ssh: the TUI for the person, an agent (Claude Code, opencode) beside it driving ORCA through the shipped skill. Today the two cannot see each other. When the person says "rotate *this* bond", the agent cannot tell which job, which geometry or which atoms they mean; when the agent has something to show, it cannot put it in front of them. This plan adds the shared view and nothing that touches ORCA: the person picks atoms in the TUI (by number with `s`, or by clicking) and sees them ringed and measured; the TUI publishes what it shows to a small state file; `orcamon view` reads it for the agent; `orcamon focus JOB` makes the TUI jump to a job, a geometry, atoms and a mode, and `b` takes the person back; `orcamon measure` gives the agent the same distance, angle or dihedral the person sees. **The goal a deviation must still meet:** an agent beside a running TUI can learn exactly which geometry and atoms the person means and can show them a job, a geometry, atoms and a mode, through files under `$XDG_STATE_HOME/orcamon` with no daemon, socket or multiplexer coupling; atom numbers everywhere are ORCA's 0-based numbers, the ones the TUI labels; nothing in this plan starts, edits or stops an ORCA job; and every existing key, command, output and JSON document behaves as before.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 208 `pass` lines, last line `0 failure(s)` | ~2 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_steer/xdg XDG_STATE_HOME=/tmp/orcamon_steer/state .venv/bin/python test_orcamon.py` | 212 `pass` lines, last line `0 failure(s)` | ~15 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~45 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_steer/xdg XDG_STATE_HOME=/tmp/orcamon_steer/state .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (the same for `test_orcamon.py`, with its environment variables). Each task states the counts it must end at.

Every task runs the first three commands. Only ST8 runs `run_gates.py`: every task touches only `tools/orcamon/`, `test_monitor.py`, `test_orcamon.py` and `.claude/skills/orcamon/SKILL.md`, which no other gate reads except the duplicate guard.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's; a new gate file would change the count and break that gate. New tests go in `test_monitor.py` or `test_orcamon.py`.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- a hand-derived expected value in this plan does not hold after the change is implemented as written (report the measured value; do not change the expected value or the design to make it pass);
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance (`test_geom_never_substitutes_a_geometry`, `test_snapshot_writes_a_png`, `test_the_tui_runs_headless`, `test_the_tui_scrubs_a_path`, `test_the_tui_cycles_real_modes` and the renderer tests pin today's behaviour and must pass unchanged);
- a real-tree value in ST1 differs from this plan's (report it, and say whether that job's files changed since `91cb1d9`);
- a Textual behaviour in § Facts does not hold (report what Textual did instead; do not pick another key, event or widget on your own);
- in ST8, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's. (`tools/orcamon/README.md` and `tools/orcamon/src/orcamon/SKILL.md` are orcamon's; `.claude/skills/orcamon/SKILL.md` changes ONLY through `orcamon skill install --project .`.)
- Never `git add -A`, `git add .` or `git commit -a`; stage the files the task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file.
- Never define the same top-level name in both `test_monitor.py` and `test_orcamon.py`, private or not: the duplicate guard reads private names too. A fixture both files need lives in `test_monitor.py` and `test_orcamon.py` imports it.
- Never let a test or a manual command read or write the user's own state or cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_steer/xdg` and `XDG_STATE_HOME=/tmp/orcamon_steer/state`; `test_orcamon.py` points `XDG_STATE_HOME` at a temporary directory at import (ST2), and the `test_monitor.py` steer test sets and restores it itself. The user's own TUI must never receive a `focus` request from this run. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never start the TUI against `computational/`, never write into `computational/`, never start ORCA, and never kill or signal a process the task did not start (`steer.alive` sends signal 0, which delivers nothing; nothing else sends a signal). The real tree appears only in ST1's read-only acceptance commands.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `core` and `cli` import only the standard library at module level; `test_the_core_needs_only_the_standard_library` imports every `orcamon.core` and `orcamon.cli` module with the TUI's dependencies refused, so `core/steer.py` must import only the standard library and other `orcamon.core` modules.
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill.
- Never change the key set of an existing `--json` document, and never change `REPORT_SCHEMA` (it stays `1`). `measure --json`, `view --json` and `focus --json` are new documents.
- Never add a fenced code block containing an `orcamon …` line to `tools/orcamon/src/orcamon/SKILL.md`: `test_orcamon.test_the_skill_examples_run` executes every such line, and `view` and `focus` exit 2 when no TUI runs. Mention the new commands in prose, in backticks.
- Never change what `i`, `I`, the chart keys, `geom`, `snapshot` or `freqs` do.

## Scope

**In scope:** `measure` and the atom-list rules in `core/geometry.py`; the `orcamon measure` command; `core/steer.py` (the state files); the selection drawn by both renderers; `s` and click picking in the geometry pane, with the measurement in its title; the TUI publishing its view; `orcamon view`; `orcamon focus` and the TUI taking it; `b`; tests, README and the shipped skill.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | `orcamon focus` makes the TUI jump straight to what the agent names — no confirmation — with a toast saying the agent did it, and `b` returns the person to where they were before the agent's last focus (one level). | user, 2026-10-01 |
| D2 | The agent→TUI command is `orcamon focus`, because `orcamon show` already names the one-job summary. The TUI→agent command is `orcamon view`. | planner, 2026-10-01 |
| D3 | Shared state is files under `$XDG_STATE_HOME/orcamon`, or `~/.local/state/orcamon` when `XDG_STATE_HOME` is unset or empty: `view-<host>-<pid>.json`, written only by the TUI with that pid, and `focus-<host>-<pid>.json`, written only by the CLI and taken (renamed, read, deleted) by that TUI. Writes are atomic (a temporary file in the same directory, then `os.replace`). No daemon, no socket, no multiplexer coupling: it works the same under tmux, herdr, screen or plain ssh, and across two ssh sessions on one machine. | planner, 2026-10-01 |
| D4 | A TUI is live when its view file is on THIS host (`<host>` is `socket.gethostname()` with every character outside `A-Za-z0-9._-` replaced by `_`) and `os.kill(pid, 0)` does not raise `ProcessLookupError` (a `PermissionError` counts as live). The CLI talks to the live TUI with the largest `updated_unix_s`: the one most recently changed. Another host's files are never read. A TUI deletes its own two files when it exits and, when it starts, the files of dead pids on this host. | planner, 2026-10-01 |
| D5 | Atom numbers are ORCA's 0-based numbers in the geometry the pane shows: the labels the pixel pane draws (`str(i)`), `geom`'s line order, and what ORCA's `%geom` blocks take. A selection holds 1–4 distinct atoms; 2 measure a distance, 3 an angle, 4 a dihedral. It is cleared when the pane changes job and kept when the geometry within a job changes (scrubbing, following a running job), and is then re-measured on the new geometry. | planner, 2026-10-01 |
| D6 | Measurements: distance in Å, written with 3 decimals (`0.960 Å`); angle and dihedral in degrees with 1 decimal (`104.5°`, `-120.0°`). The dihedral uses the IUPAC sign (positive when, looking from the 2nd atom to the 3rd, the 1st must turn clockwise to eclipse the 4th), range (-180, 180]. JSON carries `value` rounded to 4 decimals and `unit` `"angstrom"` or `"degree"`. | planner, 2026-10-01 |
| D7 | Picked atoms are always ringed (cyan, `(0, 229, 255)`) and labelled with their number in the pixel panes, whatever `l` says and including environment atoms and hidden hydrogens; in the text pane they are drawn as their number in reverse video. | planner, 2026-10-01 |
| D8 | Picking: `s` in the geometry pane opens a prompt holding the current selection; Enter applies it, an empty entry clears, Escape cancels, a bad entry is refused with a notification and changes nothing. A click toggles the atom under it; a fifth atom drops the oldest; a click on nothing does nothing. "Under it" is the atom whose projected centre is nearest the clicked cell's centre within 1.5 cell widths (vertical distance doubled, as a cell is twice as tall as wide), ties going to the atom nearest the viewer; hidden hydrogens cannot be clicked. | planner, 2026-10-01 |
| D9 | `focus` waits for the TUI to take the request, up to `--timeout` seconds (default 5): exit 0 when shown, 2 when the TUI could not show it, 4 when it did not take it in time (the request stays queued), and `--timeout 0` queues and returns 0 at once. A TUI that cannot find the job keeps the request for 10 s, asking for one rescan, so a draft input the agent has just written can still be shown. | planner, 2026-10-01 |
| D10 | The skill's rule becomes: "Never kill, signal, edit or resubmit a running job on your own initiative. orcamon never changes a job: it reads, and `orcamon focus` only moves what the person's TUI shows." The rewrite for running jobs belongs to the Act plan. | planner, 2026-10-01 |
| D11 | Commits go straight to `master`, one per task. | user, 2026-10-01 |

## Facts (verified 2026-10-01)

**Environment.** `/home/madsr2d2/masterThesis/.venv/bin/python` is Python 3.12.3 with Textual 8.2.8 and Pillow; orcamon is installed editable, console script `/home/madsr2d2/masterThesis/.venv/bin/orcamon`. `git status --porcelain` was empty at `91cb1d9`. Baseline: `test_monitor.py` 208 pass, `test_orcamon.py` 212 pass, `run_gates.py` 38 gates, 0 failed (from the modes run's progress file at `91cb1d9`; the first two re-run for this plan). README `--help` blocks are pasted at `COLUMNS=90` (`COLUMNS=90 orcamon geom --help` reproduces the README's block byte for byte).

**Textual 8.2.8, probed:**
- `App.on_unmount` runs when an `async with app.run_test()` block exits.
- `await pilot.click("#id", offset=(x, y))` delivers a `Click` to that widget with `event.x, event.y == (x, y)` relative to the widget's outer region, border included; `event.get_content_offset(widget)` returns the offset inside the border (`(x - 1, y - 1)` for a one-cell border) and `None` outside the content. Clicking a `can_focus` widget focuses it.
- A `ModalScreen` holding `Input(value="1 2")`: the cursor starts at the end; `await pilot.press("backspace", …)` deletes; `await pilot.press(*"7, 8", "enter")` types `7, 8` (space included) and `on_input_submitted(self, event)` sees `event.value == "7, 8"`; `def key_escape(self): self.dismiss(None)` on the screen cancels. `self.app.push_screen(screen, callback)` calls `callback(value)` with what `dismiss` was given, and focus returns to the widget that had it.
- `table.move_cursor(row=table.get_row_index(row_key))` moves a `DataTable`'s cursor and posts `RowHighlighted` with that key.

**The TUI** (`tools/orcamon/src/orcamon/tui/app.py`):
- `GEOMETRY_BINDINGS` (line 360) ends `("i", "cycle_modes", "Modes")`, `("I", "cycle_real_modes", "Real modes")`, then `v h p o 0 [ ]` and ctrl+arrows; `s`, `b` and mouse handlers are unused anywhere in the app. `MonitorApp.BINDINGS` (line 1188) is `q`, `r`, `m`.
- `RotatableGeometryImage.__init__` (line 417) holds the view state; `view()` (line 460) builds a `core.geometry.View`; `_atoms_for(job, point)` (static) returns `geometry_shown(job, point)[0]`; `self._job`, `self._point` are what the pane shows; `_input.request()` redraws (coalesced); `self.app.update_detail()` rebuilds the title and panes.
- `_cycle_through(self, modes, kind)` (line 650) steps the mode cycle: its `else` branch sets `self.mode_index`, calls `self._start_modes()`, sets `self._mode_alternates = vibrations.alternate_geometries(self._job)` and `self.mode_label = f"mode {i} {cm1:.1f} cm-1 ({kind}{position + 1}/{len(indices)})"`; both branches end with `self._input.request()` and `self.app.update_detail()`.
- Each backend's `show_job(job, point)` — `KittyGeometryImage` (line 809), `TextGeometry` (line 931), `HerdrGeometryImage` (line 1007) — has an `if job is not self._job:` branch that calls `self._stop_modes()`.
- `KittyGeometryImage.FULL_PX = (900, 750)` is the frame it renders and the terminal scales into the pane; `HerdrGeometryImage._push` renders `(region.width * cells.width_px, region.height * cells.height_px)` with `cells = herdr_graphics.cell_size()` (None when unknown); the text pane draws `2 * width` by `4 * height` braille dots.
- `ConvergencePlot` (line 122): `_points`, `_mode` (`"irc"`, `"neb"`, `"scan"`, `"all"`, `"opt"` or `"scf"`), `_following_latest`, `selected_point()`, `_select(index)` (sets `_following_latest = index == len(self._points) - 1` except on a path, where it sets False, then calls `self.app.update_detail()`); a new job resets `_following_latest = True` in `show_job`.
- `MonitorApp`: `self.root` (resolved), `self.jobs`, `self.selected_label`, `self._gone`, `self._force_discover`, `refresh_all()`, `selected_job()`, `update_detail()` (skips the render when the scan holds the job's lock), `_render_detail(job, placeholder=…)` (line 1541; calls `chart.show_job(job)`, `point = chart.selected_point()`, `geometry.show_job(job, point)`, `self._title_geometry(geometry, job, point)`), `_title_geometry` (line 1563; `_atoms, shown = geometry_shown(job, point)`, then appends `· {describe_point(shown)}`). `on_mount` (line 1260) sets the 3 s refresh interval. A new job found by `_rediscover` is appended to `self.jobs` on the main thread (`call_from_thread` blocks) before the same scan reads it, so it is parsed in that scan.
- `Job.inp_path` is `state.path / f"{stem}.inp"`; `state.path` is absolute and resolved.

**The renderers.** `core/geometry.py`: `View` dataclass (fields `elev, azim, zoom, pan, representation, fog, show_labels, show_distances, show_hydrogens, see_through`), `camera_basis`, `camera_forward`, `VDW_RADII`, `DEFAULT_VDW_RADIUS`, `FileGeometry(atoms, source)`. `tui/raster.py`: `render(atoms, view, *, qm_atom_indices, size_px, bond_list)` (line 462) draws labels last via `_draw_labels(image, atoms, view, sx, sy, scale, is_qm, visible, has_bond, owner, size_px)` (line 355; skips `not visible[i] or not is_qm[i]`, then occluded atoms; `str(i)`, `LABEL_RGB`, black stroke) and then `_draw_distances`; `project(atoms, view, size_px)` returns `sx, sy, depth, scale`; `_font(size_px)`; `BALL_SCALE = 0.28`. `tui/geometry_text.py`: `render(atoms, width, height, view, qm_atom_indices, bond_list)` builds `glyphs[k] = (ch, style)` and `cell_depth`, dims cells with `cell_depth[k] <= far_depth`, and appends runs of identical style objects to a `rich.text.Text`.

**The CLI** (`cli/commands.py`, `cli/__init__.py`): `add(name, help_, example, func)`, `_job(p)`, `_output(p, lines=True)` (`--json`, `--max-lines`), `_with_job` (resolve, load, `UsageError` → exit 2 with `orcamon: <message>` on stderr), `_emit_json(doc)` (prepends `"schema": 1`), `Out(max_lines)`, `_find_point(state, step, cycle)` (`cycle 9: not in the kept history (cycles 1-2)` for a missing cycle), `_path_point(job, n)` (`--point needs an IRC or NEB job; this job has no reaction path`), `describe_point`, `format_age` (`"<1 min"` under 60 s). `cmd_geom` (line 470) resolves its point (`--point` / `_find_point` / `job.file_geometry`, refusing `no coordinates printed yet, and no geometry file to read` and `<point>: coordinates not printed; nearest earlier with coordinates: …`), then builds the JSON `point` dict `{"scan_step", "cycle", "energy_eh", "energy_label", "source"}`. `cmd_snapshot`'s mode refusals: `this job has no frequency block (no Freq or NumFreq step)` and `no mode {N} in this job; {vibrations.available_text(state)}`. `EXIT_CODES` and the module docstring in `cli/__init__.py` describe exit 2 as `usage error, or a job not found or ambiguous` and 4 as `wait timed out`; `test_the_skill_reference_matches_the_parser` requires the reference to contain `` `4`: wait timed out``.

**The skill tests.** `test_the_skill_names_every_command` requires every registered command to appear in `SKILL.md`'s body as `` `orcamon <name>`` (backticked) or at a line start inside a fenced block, and every such mention to be a registered command. `test_this_repository_has_the_current_skill` requires `.claude/skills/orcamon/SKILL.md` to match, so every task that adds a command or edits `SKILL.md` ends with `cd /home/madsr2d2/masterThesis && .venv/bin/orcamon skill install --project /home/madsr2d2/masterThesis` (with `--force` only if it refuses).

**Test helpers.** `test_orcamon.py`: `check`, `_orcamon(argv)` (in-process; `(exit, stdout, stderr)`), `_json(argv)`, `_opt(cycles, energies)`, `_coords(x)` (`O (x, 0, 0)`, `H (x + 0.96, 0, 0)`, `H (x - 0.24, 0.93, 0)`, printed at `x = n / 10` for cycle n), `_CONVERGED`, `_DONE`, `_OPT_FREQ`, `_freqs`, `_normal_modes`; `from test_monitor import _irc_tree, _MODE_FREQS, _MODE_VECTORS` inside tests; the headless pattern is `test_the_tui_cycles_real_modes` (line 1298): `MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")`, `app.run_test(size=(180, 50))`, a 5 s wait until every job is parsed and the pane shows the first one, `await pilot.pause(0.2)` after keys. The module sets `os.environ["HERDR_ENV"] = "0"` at import. `test_monitor.py`: `check`, renderer tests import `raster` from `orcamon.tui`. Counts: `test_monitor.py` 208, `test_orcamon.py` 212.

**The water fixture** (`job.inp` = `_OPT_FREQ`, `job.out` = `_opt(2, [-80.00, -80.10]) + _CONVERGED + "\n" + _DONE + "\n"`): two cycles; cycle 2 is `O (0.2, 0, 0)`, `H (1.16, 0, 0)`, `H (-0.04, 0.93, 0)`; cycle 1 is the same shifted by -0.1 in x. On both: distance 0–1 is 0.96 → `distance 0-1: 0.960 Å`; angle 1–0–2 is 104.4702941° → JSON 104.4703, text `angle 1-0-2: 104.5°`. `modes_all` (the modes plan's fixture) is this geometry plus `_freqs(_MODE_FREQS) + _normal_modes(_MODE_VECTORS)`; its real modes lowest first are `[5, 4, 7, 6, 8]`, so mode 5 is `mode 5 50.0 cm-1 (real 1/5)`. `_irc_tree`'s path points are `-2 IRC backward 1`, `-1 IRC backward 0`, `0 IRC TS`, `1 IRC forward 0`, `2 IRC forward 1`, `3 IRC forward 2`; following, the pane shows `IRC TS`.

**Hand-derived measurement values** (prototype agrees to 1e-12): `(0,0,0)`–`(3,4,0)` distance 5.0; `(1,0,0)`, `(0,0,0)`, `(0,1,0)` angle 90.0; with `a=(1,0,0)`, `b=(0,0,0)`, `c=(0,0,1)`, `d=(cos φ, sin φ, 1)` the dihedral a-b-c-d is φ: `d=(0.5, 0.8660254037844386, 1.0)` gives 60.0 and `d=(-0.5, -0.8660254037844386, 1.0)` gives -120.0 (after `round(…, 4)`), and the reversed order d-c-b-a gives the same 60.0. The formula: `b1 = p1 - p0`, `b2 = p2 - p1`, `b3 = p3 - p2`, `y = |b2| * (b1 · (b2 × b3))`, `x = (b1 × b2) · (b2 × b3)`, dihedral `= degrees(atan2(y, x))`; the angle at the middle atom is `degrees(acos(clamp(u·v / (|u||v|), -1, 1)))` with `u = p0 - p1`, `v = p2 - p1`.

**Real-tree measurements** (`computational/`, read-only; `geometry_b973c-xtb/ts/optts_freq_tight`, latest = cycle 41, 137 atoms, label `C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts_freq_tight` under `--root computational`): 128–132 `C-H` 2.2988 Å; 128–132–133 `C-H-O` 46.3461°; 125–128–132–133 `H-C-H-O` 114.9543°.

**Projection and picking** (prototype): a pure-Python `screen_positions(atoms, view, frame)` written as ST8 specifies reproduces `raster.project`'s `sx, sy, depth` to within 2e-13 px for all four representations, at the default view and at `View(elev=35, azim=110, zoom=1.7, pan=(0.3, -0.2, 0.1))`, on a 5-atom C/O/H/N/K set at `(900, 750)`. The pick fixture `five = [("C", 0, 0, 0), ("O", 0, 2, 0), ("H", 0, -2, 0), ("N", 1, 0, 0), ("S", -1, 0, 0)]` under `View(elev=0, azim=0)` (camera on +x, screen right = +y, screen up = +z) in a `(100, 100)` frame: mean `(0, 0, 0)`, extent 2 + 0.5, scale 20 px/Å, positions `(50, 50, 0)`, `(90, 50, 0)`, `(10, 50, 0)`, `(50, 50, 1)`, `(50, 50, -1)`. With `cells = (20, 10)`: click `(9, 4)` → 3 (atoms 0, 3 and 4 tie at d = 1.118; 3 is nearest the viewer), `(17, 5)` → 1, `(14, 5)` → None, `(1, 4)` → 2, and `(1, 4)` with `show_hydrogens=False` → None.

**Selection rendering** (prototype, water cycle-2 geometry, default `View()`, `(300, 240)`): ring radius `max(8, round(1.4 * 0.28 * vdW * scale))` with scale 93.45 px/Å is 56 px for O and 40 px for H; selecting atom 1 gives 716 pixels exactly `(0, 229, 255)`, and with nothing selected no render has a pixel of that colour (checked in all four representations, with and without `qm_atom_indices={0}`); with `qm_atom_indices={0}` and atom 2 selected, the 9×9 box around atom 2's centre differs from the unselected render (its label); with `show_labels=False` and atom 1 selected, the box around atom 1 differs from the unselected render; with `show_hydrogens=False` and atom 1 selected, 716 ring pixels; a selection of `(7,)` draws nothing.

## Conventions

- Match the surrounding code: one-line comments where the reason is not obvious; docstrings say WHY.
- New tests go in the file the task names, use `check(...)` for every assertion, build trees under `tempfile.TemporaryDirectory()`, and are added to that file's `__main__` list right after the test the task names. CLI calls in new tests pass `["--root", str(root), "--liveness", "mtime", "--no-cache", …]` unless the task says otherwise.
- Expected values are the hand-derived ones in § Facts, never values copied from a run. Floats are compared after `round(…, 4)` unless § Facts gives them exactly.
- A new `orcamon` command's README entry is a `### \`orcamon <name>\`` section in `## Command reference`, holding the `COLUMNS=90 .venv/bin/orcamon <name> --help` output in a fenced block like its neighbours; whenever a command is added, the reference's first block is replaced with `COLUMNS=90 .venv/bin/orcamon --help`.
- No network access.

---

# Phase 1 — measure, and the state files

### ST1 — `measure`: one rule for distances, angles and dihedrals, and the command

**Why:** the person and the agent must get the same number for "the distance between 129 and 55", and the agent today would compute it by hand from `geom` output. The same task adds the one parser for atom lists that `s`, `focus` and `measure` share (D5, D6).
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/geometry.py`, `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_monitor.py`, `test_orcamon.py`.

**Do:**
1. `core/geometry.py` (pure Python): add `MAX_PICKED = 4` and
   - `measure(atoms, indices) -> tuple[str, float]`: raise `ValueError(f"measure takes 2 to 4 atoms, got {len(indices)}")` unless `2 <= len(indices) <= 4`; then `_check_in_range(indices, len(atoms))`; then `ValueError("an atom is named twice")` on a repeat. Return `("distance", Å)`, `("angle", degrees)` or `("dihedral", degrees)` by the formulas in § Facts "Hand-derived measurement values".
   - `_check_in_range(indices, n)`: `ValueError("there is no geometry to pick atoms from")` when `n == 0`; else for the first `i` outside `0 <= i < n`, `ValueError(f"atom {i} is out of range: the geometry has {n} atoms (0-{n - 1})")`.
   - `measurement_text(kind, indices, value) -> str`: `f"{kind} {'-'.join(str(i) for i in indices)}: "` + (`f"{value:.3f} Å"` for a distance, `f"{value:.1f}°"` otherwise).
   - `describe_measurement(atoms, indices) -> str`: `""` when `indices` is empty or any index is outside `range(len(atoms))`; `f"atom {i} {atoms[i][0]}"` for one atom; else `kind, value = measure(atoms, indices)` and `measurement_text(kind, indices, value)`.
   - `parse_atom_list(text, n) -> tuple[int, ...]`: split `text.strip()` on `re.split(r"[\s,]+", …)`, dropping empty tokens; `()` when none; each token through `int()`, else `ValueError(f"not an atom number: {token!r}")`; more than `MAX_PICKED` → `ValueError("at most 4 atoms: 2 measure a distance, 3 an angle, 4 a dihedral")`; a repeat → `ValueError("an atom is named twice")`; then, when `n is not None`, `_check_in_range`.
   Docstrings say the numbers are ORCA's 0-based atom numbers as the pane labels them, and give the dihedral sign rule (D6).
2. `cli/commands.py`: move the first part of `cmd_geom` — from `state = job.state` through the `coordinates not printed` refusal — into `_requested_point(args, job)`, which returns the point, and the JSON `point` dict into `_point_doc(point)`; `cmd_geom` keeps its own `state = job.state`, calls both, and is otherwise unchanged (same messages, same output). Add `@_with_job def cmd_measure(args, job)`: `point = _requested_point(args, job)`; `kind, value = measure(point.atoms, args.atoms)` with `ValueError` re-raised as `UsageError(str(exc))`; JSON `{"job": job.label, "point": _point_doc(point), "atoms": [{"atom": i, "element": point.atoms[i][0]}, …], "kind": kind, "value": round(value, 4), "unit": "angstrom" if kind == "distance" else "degree"}`; text, one line: `f"{measurement_text(kind, args.atoms, value)} ({'-'.join(elements)}) · {job.label}"` + `f" · {where}"` when `where = describe_point(point)` is non-empty. Return `_job_exit(job)`.
3. `cli/__init__.py`: register `measure` right after `geom`: `add("measure", "the distance, angle or dihedral between 2 to 4 atoms of a geometry", "orcamon measure opt/ts 12 15", commands.cmd_measure)`, `_job(p)`, `p.add_argument("atoms", nargs="+", type=int, metavar="ATOM", help="ORCA's atom numbers, 0-based as the TUI labels them: 2 give a distance, 3 an angle, 4 a dihedral")`, then `--cycle`, `--step` and `--point` exactly as `geom` declares them, then `_output(p, lines=False)`.
4. `SKILL.md`, "Looking closer", one prose paragraph after the `freqs --mode` one: "`orcamon measure JOB I J [K [L]]` gives the distance (2 atoms, Å), angle (3) or dihedral (4, degrees, IUPAC sign) between atoms by ORCA's 0-based numbers — the ones the TUI labels — on the latest geometry or on the one `--cycle`, `--step` or `--point` names. Use it rather than computing from `orcamon geom` output yourself." README: the `measure` reference section after `geom`'s and the refreshed first block (§ Conventions). Then reinstall the skill (§ Facts "The skill tests").
5. `test_monitor.py`: `test_measurements_follow_the_iupac_convention()`, registered right after `test_modes_are_found_and_described`. Exactly 8 checks, values from § Facts:
   - distance `(0,0,0)`–`(3,4,0)`: `("distance", 5.0)`;
   - angle `(1,0,0)`, `(0,0,0)`, `(0,1,0)`: kind `"angle"`, `round(value, 4) == 90.0`;
   - dihedral with `d=(0.5, 0.8660254037844386, 1.0)`: `round(…, 4) == 60.0`;
   - with `d=(-0.5, -0.8660254037844386, 1.0)`: `-120.0`;
   - the 60° set measured as `[3, 2, 1, 0]`: `60.0`;
   - on the water cycle-2 atoms `[("O", 0.2, 0.0, 0.0), ("H", 1.16, 0.0, 0.0), ("H", -0.04, 0.93, 0.0)]`: `measure(w, [0])`, `[0, 5]`, `[0, 0]` and `[0, 1, 2, 0, 1]` raise `ValueError` with exactly `measure takes 2 to 4 atoms, got 1`, `atom 5 is out of range: the geometry has 3 atoms (0-2)`, `an atom is named twice`, `measure takes 2 to 4 atoms, got 5`;
   - `describe_measurement(w, [1, 0, 2]) == "angle 1-0-2: 104.5°"`, `[0, 1]` → `"distance 0-1: 0.960 Å"`, `[2]` → `"atom 2 H"`, `[0, 9]` → `""`;
   - `parse_atom_list("129, 55", 200) == (129, 55)`, `parse_atom_list("  ", 3) == ()`, and `"1 2 3 4 5"`, `"a"`, `"1 1"`, `"7"` (n = 3) raise `at most 4 atoms: 2 measure a distance, 3 an angle, 4 a dihedral`, `not an atom number: 'a'`, `an atom is named twice`, `atom 7 is out of range: the geometry has 3 atoms (0-2)`.
6. `test_orcamon.py`: add `_water_job(root, name="water")` (writes `root/name/job.inp` and `job.out` as § Facts "The water fixture") and `test_measure_names_the_geometry_it_measured()`, registered right after `test_any_mode_is_drawn_and_described`, on a temp root holding `water`. Exactly 6 checks:
   - `measure water 0 1`: exit 0, stdout `distance 0-1: 0.960 Å (O-H) · water · cycle 2\n`;
   - `measure water 1 0 2 --json`: `doc["kind"] == "angle"`, `doc["value"] == 104.4703`, `doc["unit"] == "degree"`, `doc["atoms"] == [{"atom": 1, "element": "H"}, {"atom": 0, "element": "O"}, {"atom": 2, "element": "H"}]`, `doc["point"]["cycle"] == 2`, `doc["job"] == "water"`;
   - `measure water 0 1 --cycle 1`: stdout `distance 0-1: 0.960 Å (O-H) · water · cycle 1\n`;
   - `measure water 0 5`: exit 2, stderr contains `atom 5 is out of range: the geometry has 3 atoms (0-2)`;
   - `measure water 0`: exit 2, stderr contains `measure takes 2 to 4 atoms, got 1`;
   - `measure water 0 1 --cycle 9`: exit 2, stderr contains `cycle 9: not in the kept history`.
7. Run the first three verification commands, and on the real tree (`cd /home/madsr2d2/masterThesis`, the two environment variables of § Never): `.venv/bin/orcamon --root computational measure geometry_b973c-xtb/ts/optts_freq_tight 128 132`, then `… 128 132 133`, then `… 125 128 132 133`.

**Acceptance:**
- The real tree prints `distance 128-132: 2.299 Å (C-H) · C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts_freq_tight · cycle 41`, `angle 128-132-133: 46.3° (C-H-O) · … · cycle 41` and `dihedral 125-128-132-133: 115.0° (H-C-H-O) · … · cycle 41`, each exiting 0.
- `test_measurements_follow_the_iupac_convention` passes its 8 checks and `test_measure_names_the_geometry_it_measured` its 6; `test_geom_never_substitutes_a_geometry` passes unchanged.
- `test_monitor.py`: 216 `pass`; `test_orcamon.py`: 218 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### ST2 — The state files: `core/steer.py`

**Why:** the TUI and the CLI need one module that agrees on where the files are, how they are written, which TUI is live and how a request is taken, so neither side can drift from the other (D3, D4).
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/steer.py` (new), `test_monitor.py`, `test_orcamon.py`.

**Do:**
1. Create `core/steer.py` (standard library plus `.geometry` and `.paths` only) with a module docstring giving D3 and D4 in prose, `STEER_SCHEMA = 1`, and:
   - `state_dir() -> Path`: `Path(os.environ["XDG_STATE_HOME"]) / "orcamon"` when that variable is set and non-empty, else `Path.home() / ".local" / "state" / "orcamon"`.
   - `host_name() -> str`: `re.sub(r"[^A-Za-z0-9._-]", "_", socket.gethostname()) or "localhost"`.
   - `view_path(pid: int | None = None) -> Path` → `state_dir() / f"view-{host_name()}-{pid or os.getpid()}.json"`; `focus_path(pid: int) -> Path` → `state_dir() / f"focus-{host_name()}-{pid}.json"`.
   - `write_json(path, doc)`: `path.parent.mkdir(parents=True, exist_ok=True)`, write `json.dumps(doc)` to `path.with_name(f".{path.name}.{os.getpid()}.tmp")`, then `os.replace` it onto `path`.
   - `read_json(path) -> dict | None`: the parsed document when the file exists, parses and is a dict; else None (never raises for a missing, unreadable or malformed file).
   - `alive(pid: int) -> bool`: `os.kill(pid, 0)`; False on `ProcessLookupError`, True on `PermissionError`, True otherwise. Comment: signal 0 delivers nothing; it only asks whether the pid exists.
   - `live_views() -> list[dict]`: for every `state_dir().glob(f"view-{host_name()}-*.json")` whose suffix between the prefix and `.json` is an integer pid, keep `read_json(path)` when it is not None, its `"pid"` equals that pid and `alive(pid)`; sort by `doc.get("updated_unix_s", 0.0)`, largest first. `current_view() -> dict | None`: the first, or None.
   - `take_focus(pid: int) -> dict | None`: `os.replace(focus_path(pid), claimed)` with `claimed = focus_path(pid).with_name(focus_path(pid).name + ".taken")` (None on `FileNotFoundError`), `doc = read_json(claimed)`, `claimed.unlink(missing_ok=True)`, return `doc`. Comment: the rename claims the request atomically, so one the CLI writes meanwhile is kept for the next check.
   - `remove_stale() -> int`: unlink every `view-{host_name()}-<pid>.json` and `focus-{host_name()}-<pid>.json` (integer pid) whose pid is not `alive`; return how many were removed. Another host's files are never touched.
   - `geometry_address(shown, following: bool) -> list[str]`: `["--point", str(shown.index)]` for a `paths.PathPoint`; `[]` when `following`, when `shown` is None or a `geometry.FileGeometry`; else `["--step", str(shown.scan_step)]` when `scan_step is not None`, followed by `["--cycle", str(shown.cycle)]` when `shown.cycle` is truthy. Docstring: these are `orcamon geom`'s own flags for the geometry the pane shows.
2. `test_orcamon.py`: at module level, right after the `HERDR_ENV` line, `import tempfile` (at the top imports) and:
   ```python
   # The TUI publishes what it shows under $XDG_STATE_HOME/orcamon, and `view`
   # and `focus` read it there: a run must never see -- or steer -- the
   # person's own TUI.
   _STATE_HOME = tempfile.TemporaryDirectory(prefix="orcamon_state_")
   os.environ["XDG_STATE_HOME"] = _STATE_HOME.name
   ```
3. `test_monitor.py`: `test_the_steer_files_name_one_live_tui()`, registered right after `test_a_job_that_prints_no_geometry_shows_the_one_it_wrote` (the last in the list). It saves `XDG_STATE_HOME` and `HOME`, works in a `tempfile.TemporaryDirectory()` as `XDG_STATE_HOME`, and restores both in a `finally`. `dead` is the pid of a finished `subprocess.Popen([sys.executable, "-c", "pass"])` after `.wait()`; `me = os.getpid()`, `parent = os.getppid()`. A view document is `{"pid": p, "updated_unix_s": t}`. Exactly 8 checks:
   - `state_dir() == Path(tmp) / "orcamon"`; and with `XDG_STATE_HOME` set to `""` and `HOME` set to a second temp dir, `state_dir() == Path(home) / ".local" / "state" / "orcamon"`;
   - `write_json(p, {"a": 1})` then `read_json(p) == {"a": 1}`, and `sorted(os.listdir(p.parent)) == [p.name]`;
   - `read_json` of a file holding `not json`, of a file holding `[1]` and of a missing path are all None;
   - `alive(me) is True and alive(dead) is False`;
   - with view files for `me` (updated 100.0), `parent` (200.0) and `dead` (300.0) written at `view_path(pid)`, and one for `me` at `state_dir() / f"view-otherhost-{me}.json"` (400.0): `[d["pid"] for d in live_views()] == [parent, me]` and `current_view()["pid"] == parent`;
   - `write_json(focus_path(me), {"id": "x"})`: `take_focus(me) == {"id": "x"}`, then `focus_path(me).exists() is False` and `take_focus(me) is None`;
   - with `write_json(focus_path(dead), {})` added: `remove_stale() == 2`, `view_path(dead)` and `focus_path(dead)` are gone, and the `me`, `parent` and `otherhost` view files remain;
   - `geometry_address(PathPoint(kind="irc", index=-2, label="IRC backward 1", energy=None, de_kj_mol=0.0, atoms=[]), False) == ["--point", "-2"]`, `GeometryPoint(scan_step=3, cycle=7)` with `following=False` → `["--step", "3", "--cycle", "7"]` and with `True` → `[]`, `GeometryPoint(scan_step=None, cycle=4)` → `["--cycle", "4"]`, `FileGeometry([], "input geometry")` → `[]`, `None` → `[]`.
4. Run the first three verification commands.

**If unsure:** if the test runs with no parent process the user owns (`alive(parent)` False), report it; do not replace `parent` with another pid on your own.

**Acceptance:**
- `test_the_steer_files_name_one_live_tui` passes its 8 checks; `test_the_core_needs_only_the_standard_library` passes (it imports `orcamon.core.steer` with the TUI's dependencies refused).
- `test_monitor.py`: 224 `pass`; `test_orcamon.py`: 218 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

# Phase 2 — the person picks atoms

### ST3 — Both renderers draw the selection

**Why:** a picked atom must be findable in a 137-atom structure, including an environment atom, which today gets no label (`raster._draw_labels` skips `not is_qm[i]`), and a hidden hydrogen (D7).
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/geometry.py`, `tools/orcamon/src/orcamon/tui/raster.py`, `tools/orcamon/src/orcamon/tui/geometry_text.py`, `test_monitor.py`.

**Do:**
1. `View`: add `selected: tuple = ()` as the last field, commented: ORCA's numbers of the atoms the person picked; ringed and labelled whatever the other options say.
2. `raster.py`: constants `SELECTION_RGB = (0, 229, 255)`, `SELECTION_RING_SCALE = 1.4`, `SELECTION_RING_MIN_PX = 8`, `SELECTION_RING_WIDTH = 3`. Add `_draw_selection(image, atoms, view, sx, sy, scale, size_px)`: for each `i` in `view.selected` with `0 <= i < len(atoms)`, `r = max(SELECTION_RING_MIN_PX, int(round(SELECTION_RING_SCALE * BALL_SCALE * VDW_RADII.get(element, DEFAULT_VDW_RADIUS) * scale)))`, `draw.ellipse((sx[i] - r, sy[i] - r, sx[i] + r, sy[i] + r), outline=SELECTION_RGB, width=SELECTION_RING_WIDTH)`, then the label exactly as `_draw_labels` draws one (`str(i)`, `_font(size_px)`, `LABEL_RGB`, `anchor="mm"`, black stroke 1) — with no visibility, layer or occlusion test.
3. `raster.render`: after the distances, `if atoms and view.selected: _draw_selection(image, atoms, view, sx, sy, scale, (width, height))`. `_draw_labels` skips any `i` in `set(view.selected)` (it is labelled by `_draw_selection`).
4. `geometry_text.py`: `SELECTED_STYLE = Style(reverse=True, bold=True)`. In `render`, after the atom-glyph loop and before fog: for each `i` in `view.selected` with `0 <= i < n`, take `col, row` as the glyph loop does; when `0 <= row < height`, write each character of `str(i)` at `col + off` (inside `0 <= col + off < width`) as `glyphs[k] = (ch, SELECTED_STYLE)` with `cell_depth[k] = math.inf`, so fog never dims it — whatever `visible[i]` says.
5. `test_monitor.py`: `test_the_selection_is_ringed_and_labelled()`, registered right after `test_text_mode_draws_the_right_enantiomer`. `w = [("O", 0.2, 0.0, 0.0), ("H", 1.16, 0.0, 0.0), ("H", -0.04, 0.93, 0.0)]`, size `(300, 240)`, `ring(img)` = the count of pixels equal to `(0, 229, 255)` (via `numpy.asarray(img)`), `box(a, b, i)` = whether any pixel in the 9×9 box around `raster.project(w, view, (300, 240))`'s `(round(sx[i]), round(sy[i]))` differs between images `a` and `b`. Exactly 7 checks, values from § Facts "Selection rendering":
   - `View(selected=(1,))`: `ring > 0`; `View()`: `ring == 0`;
   - `qm_atom_indices={0}`: `box(render with selected=(2,), render with selected=(), 2)` is True;
   - `View(show_labels=False, selected=(1,))` against `View(show_labels=False)`: `box(…, 1)` is True;
   - `View(show_hydrogens=False, selected=(1,))`: `ring > 0`;
   - `View(selected=(7,))`: renders without raising and `ring == 0`;
   - `geometry_text.render(w, 40, 12, View(selected=(2,)))`: some span of the returned `Text` has a style with `reverse` True and covers text containing `2`;
   - the same with `View(show_hydrogens=False, selected=(2,))`.
6. Run the first three verification commands.

**Acceptance:**
- `test_the_selection_is_ringed_and_labelled` passes its 7 checks; `test_labels_follow_occlusion`, `test_wireframe_keeps_labels_and_distances` and every other renderer test pass unchanged.
- `test_monitor.py`: 231 `pass`; `test_orcamon.py`: 218 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### ST4 — `s` picks atoms in the geometry pane, and the title measures them

**Why:** the person needs to point at atoms by the numbers they read off the labels, and see the distance, angle or dihedral at once (D5, D8).
**Depends on:** ST3.
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. Imports: `from textual.screen import ModalScreen`, `Input` from `textual.widgets`, and `describe_measurement`, `parse_atom_list` from `..core.geometry`. Add `class AtomPrompt(ModalScreen)` before `RotatableGeometryImage`: docstring "`s` in the geometry pane: up to four atom numbers"; `DEFAULT_CSS` centring a `Vertical` 64 wide with a round `$accent` border; `__init__(self, value: str)`; `compose` yields, inside a `Vertical`, `Static("atom numbers as the labels show them: 2 give a distance, 3 an angle, 4 a dihedral; empty clears")` and `Input(value=value, id="atoms")`; `on_input_submitted(self, event)` → `self.dismiss(event.value)`; `key_escape(self)` → `self.dismiss(None)`.
2. `GEOMETRY_BINDINGS`: add `("s", "select_atoms", "Select atoms")` right after `("I", "cycle_real_modes", "Real modes")`. `RotatableGeometryImage.__init__`: `self.selection: tuple = ()`, commented (D5). `view()` passes `selected=self.selection`.
3. `RotatableGeometryImage`: `action_select_atoms` pushes `AtomPrompt(" ".join(str(i) for i in self.selection))` with callback `self._atoms_entered`; `_atoms_entered(value)`: return on None; `n = len(self._atoms_for(self._job, self._point))`; `parse_atom_list(value, n)`, on `ValueError` `self.notify(str(exc), severity="error")` and return; else `self.set_selection(chosen)`. `set_selection(atoms)`: `self.selection = tuple(atoms)`, `self._input.request()`, `self.app.update_detail()`.
4. In each of the three `show_job` methods, inside `if job is not self._job:`, add `self.selection = ()` (D5: a new job clears it; a new geometry of the same job keeps it).
5. `MonitorApp._title_geometry`: rename `_atoms` to `atoms`; after the `where` part, `picked = describe_measurement(atoms, geometry.selection)` and, when non-empty, `title += f" · {picked}"`.
6. README, the `Keys:` bullet of "The TUI over ssh": after the `I` clause add "`s` picks atoms by number (up to four, as the labels show them; the pane title gives their distance, angle or dihedral, and picked atoms are ringed and labelled, environment atoms and hidden hydrogens included)".
7. `test_orcamon.py`: `test_the_tui_selects_atoms_by_number()`, registered right after `test_the_tui_cycles_real_modes`. Temp root holding `_water_job(root, "water")` and `_water_job(root, "water2")`; headless as § Facts "Test helpers"; wait until both are parsed and the pane shows `water`; `geometry.focus()`; `await pilot.pause(0.2)` after each key sequence; `title = lambda: geometry.border_title`. Exactly 6 checks:
   - `s`, `*"1 0 2"`, `enter`: `geometry.selection == (1, 0, 2)` and `"angle 1-0-2: 104.5°" in title()`;
   - `chart = app.query_one("#chart")`, `chart.focus()`, `left`: `"cycle 1" in title()` and `"angle 1-0-2: 104.5°" in title()` and the selection is unchanged;
   - `geometry.focus()`, `s`, five `backspace`, `7`, `enter`: `geometry.selection == (1, 0, 2)` (refused);
   - `s`, five `backspace`, `enter`: `geometry.selection == ()` and `"angle" not in title()`;
   - `s`, `*"0 1"`, `enter`: `"distance 0-1: 0.960 Å" in title()`;
   - `app.query_one("#job_table").focus()`, `down`, `await pilot.pause(0.4)`: `app.selected_label == "water2"` and `geometry.selection == ()`.
8. Run the first three verification commands.

**If unsure:** if `pilot.press(*"1 0 2")` does not type the spaces (§ Facts says it does), stop and report what the `Input` held.

**Acceptance:**
- `test_the_tui_selects_atoms_by_number` passes its 6 checks; `test_the_tui_runs_headless`, `test_the_tui_scrubs_a_path` and `test_the_tui_cycles_real_modes` pass unchanged.
- `test_orcamon.py`: 224 `pass`; `test_monitor.py`: 231 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

# Phase 3 — the agent sees and shows

### ST5 — The TUI publishes what it shows; `orcamon view` reads it

**Why:** "this bond" resolves only if the agent can read which job, which geometry and which atoms the person's TUI shows (D3, D4, D5).
**Depends on:** ST2, ST4.
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. `ConvergencePlot`: a read-only property `following` returning `self._following_latest`.
2. `MonitorApp` (`import os, shlex`; `from ..core import steer`; `measure` from `..core.geometry`): `self._published: dict | None = None` and `self._last_focus: dict | None = None` in `__init__`. `on_mount` calls `steer.remove_stale()` (OSError ignored) before the first refresh. Add `on_unmount(self)`: unlink `steer.view_path()` and `steer.focus_path(os.getpid())` with `missing_ok=True`, `OSError` ignored.
3. `_publish_view(self, job, point)`, called as the last statement of `_render_detail` on both of its paths (with `job` None on the placeholder path). It builds `{"schema": steer.STEER_SCHEMA, "pid": os.getpid(), "host": steer.host_name(), "root": str(self.root), "job": None, "geometry": None, "mode": None, "selection": [], "measurement": None, "as_xyz": None, "focus": self._last_focus}`; when `job` is not None: `atoms, shown = geometry_shown(job, point)`, `following = chart.following`, `args = steer.geometry_address(shown, following)`, `file = str(job.inp_path)`, `"job": {"label": job.label, "file": file}`, `"geometry": {"shown": describe_point(shown) or "latest", "following": following, "args": args}`, `"as_xyz": shlex.join(["orcamon", "geom", file, *args])`, `"mode": {"index": m.index, "cm1": m.cm1}` when `geometry.mode_index` names a mode (`vibrations.find_mode`), `"selection": list(geometry.selection)`, and, for 2–4 picked atoms, `"measurement": {"kind": kind, "value": round(value, 4), "unit": "angstrom" if kind == "distance" else "degree"}` from `measure(atoms, selection)` (None on `ValueError`). When the document differs from `self._published`, set `"updated_unix_s": time.time()`, `steer.write_json(steer.view_path(), doc)` (`OSError` ignored: a read-only state directory must not stop the TUI), and keep the document without `updated_unix_s` in `self._published`.
4. `cli/commands.py`: `_no_tui()` returns `f"no orcamon TUI is running for this user on {steer.host_name()}; ask the person to start `orcamon tui`"` (the backticks are part of the message). `cmd_view(args)` (not `_with_job`): `doc = steer.current_view()`; None → `return _error(_no_tui())`; `age = max(0.0, time.time() - doc.get("updated_unix_s", time.time()))`; JSON: `_emit_json({**doc, "age_s": round(age, 1)})`; text through `Out(args.max_lines)`, each label left-justified to 11 characters: `f"TUI pid {pid} on {host}, watching {root}, last change {format_age(age)} ago"`; `job` (the label, or `none` and stop there); `file`; `geometry` = `f"{shown} ({'following the newest' if following else 'scrubbed'})"`; `mode` = `f"{index} ({cm1:.1f} cm-1)"` only when set; `selection` = `measurement_text(kind, selection, value)` when there is a measurement, else the numbers space-separated, else `none`; `as XYZ` = `as_xyz`. Return 0.
5. `cli/__init__.py`: register `add("view", "what the person's orcamon TUI shows: job, geometry, picked atoms", "orcamon view --json", commands.cmd_view)` after `measure`, with `_output(p)`. Exit 2's meaning becomes `usage error, a job not found or ambiguous, or (view, focus) no TUI running` in `EXIT_CODES`, and the module docstring's exit sentence says `2 usage error, no such job or no running TUI`.
6. `SKILL.md`: a new section `## Working beside a person`, after "When the commands do not answer": "A person may be watching the same jobs in `orcamon tui`, usually in the pane beside yours. When they say "this", "these atoms", "the bond I picked" or "what I'm looking at", run `orcamon view` first. It says which job and which geometry their TUI shows, the atoms they picked — ORCA's 0-based numbers, as the TUI labels them — with the distance, angle or dihedral between them, and the `orcamon geom` command that prints exactly that geometry. It exits 2 when no TUI is running for this user on this machine: then ask which job and atoms they mean; do not guess." README: a section `## Working beside an agent` after "The TUI over ssh" saying where the TUI writes its view file (D3), that `orcamon view` reads the newest live one, and that nothing there touches a job; the exit-code table's 2 row as in step 5; the `view` reference section and refreshed first block (§ Conventions). Reinstall the skill.
7. `test_orcamon.py`: `test_the_tui_tells_the_agent_what_it_shows()`, registered right after `test_the_tui_selects_atoms_by_number`. It runs in its own `tempfile.TemporaryDirectory()` set as `XDG_STATE_HOME` (restored to the previous value in a `finally`). Temp root holding `_water_job(root, "water")`; `file = str((root / "water" / "job.inp").resolve())`; `doc = lambda: steer.read_json(steer.view_path())`; before starting the app, write `{"pid": dead, "updated_unix_s": 1.0}` to `steer.view_path(dead)` and `{}` to `steer.focus_path(dead)` (`dead` as in ST2). Exactly 7 checks:
   - after the app has parsed `water`: neither dead-pid file exists;
   - `doc()` has `pid == os.getpid()`, `root == str(root.resolve())`, `job == {"label": "water", "file": file}`, `geometry == {"shown": "cycle 2", "following": True, "args": []}`, `selection == []`, `measurement is None`, `mode is None`, `focus is None`, `as_xyz == shlex.join(["orcamon", "geom", file])`;
   - after `s`, `*"1 0 2"`, `enter` in the geometry pane: `selection == [1, 0, 2]` and `measurement == {"kind": "angle", "value": 104.4703, "unit": "degree"}`;
   - after `left` in the chart: `geometry == {"shown": "cycle 1", "following": False, "args": ["--cycle", "1"]}` and `as_xyz == shlex.join(["orcamon", "geom", file, "--cycle", "1"])`;
   - `_orcamon(["view"])`: exit 0; the first line starts `f"TUI pid {os.getpid()} on "`; the lines include `job        water`, `geometry   cycle 1 (scrubbed)`, `selection  angle 1-0-2: 104.5°`, and the last is `"as XYZ     " + shlex.join(["orcamon", "geom", file, "--cycle", "1"])`;
   - `_json(["view"])`: `doc["schema"] == 1`, `doc["job"]["label"] == "water"`, and `"age_s"` is a key;
   - after the `async with` block has exited: `steer.view_path().exists() is False`, and `_orcamon(["view"])` exits 2 with `no orcamon TUI is running` in stderr.
8. Run the first three verification commands.

**Acceptance:**
- `test_the_tui_tells_the_agent_what_it_shows` passes its 7 checks; every existing TUI and skill test passes unchanged.
- `test_orcamon.py`: 231 `pass`; `test_monitor.py`: 231 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### ST6 — `orcamon focus` writes the request and waits for the TUI

**Why:** the agent's half of "show me": it names a job, a geometry, atoms and a mode, checked against the job before anything reaches the person's screen (D1, D9).
**Depends on:** ST5.
**Files:** `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. `cli/__init__.py`: register after `view`: `add("focus", "show a job in the person's orcamon TUI: a geometry, picked atoms, a mode", "orcamon focus opt/ts --atoms 12,15", commands.cmd_focus)`, `_job(p)`, `--atoms` (`metavar="I,J,…"`, help `"atoms to ring and measure: ORCA's 0-based numbers, comma-separated, up to 4"`), `--point`, `--step`, `--cycle` exactly as `geom` declares them, `--mode` (`type=int`, `metavar="N"`, help `"animate this normal mode (ORCA's mode number)"`), `--timeout` (`type=float`, `default=5.0`, `metavar="S"`, help `"seconds to wait for the TUI to take it; 0 queues it and returns (default: %(default)s)"`), `_output(p, lines=False)`. Exit 4's meaning becomes `wait timed out, or focus was not taken by the TUI in time` in `EXIT_CODES`, and the docstring's exit sentence says ``4 `wait` timed out or `focus` not taken``.
2. `commands.py`: `FOCUS_POLL_S = 0.1`. `@_with_job def cmd_focus(args, job)`, in this order: `view = steer.current_view()`, None → `UsageError(_no_tui())`; more than one of `--point/--step/--cycle` → `UsageError("focus takes one of --point, --step or --cycle")`; `point = _requested_point(args, job)`; `atoms = parse_atom_list(args.atoms or "", len(point.atoms))` (`ValueError` → `UsageError`); with `--mode`, refuse exactly as `cmd_snapshot` does (no frequency block; `find_mode` None with `available_text`); `file = job.inp_path.resolve()`, `tui_root = Path(view["root"])`, and when `not file.is_relative_to(tui_root)` → `UsageError(f"the TUI (pid {pid}) watches {tui_root}; {file} is outside it")`.
3. Then write `{"schema": steer.STEER_SCHEMA, "id": f"{time.time_ns()}-{os.getpid()}", "file": str(file), "point": args.point, "step": args.step, "cycle": args.cycle, "mode": args.mode, "atoms": list(atoms), "requested_unix_s": time.time()}` to `steer.focus_path(pid)`. With `args.timeout <= 0`: text `f"queued for the TUI (pid {pid}); `orcamon view` shows what it took"`, exit 0. Else poll `steer.read_json(steer.view_path(pid))` every `FOCUS_POLL_S` until its `"focus"` has this `"id"` or `args.timeout` passes: not taken → stderr `f"orcamon: the TUI (pid {pid}) has not taken the request within {args.timeout:g} s; it stays queued"`, return `EXIT_TIMEOUT`; taken with `ok` → stdout `f"{message} in the TUI (pid {pid})"`, exit 0; taken not ok → `return _error(f"the TUI (pid {pid}) could not show it: {message}")`. JSON instead of text: `{"tui_pid": pid, "request": id, "taken": bool, "ok": bool | None, "message": str | None}`, same exit codes.
4. `SKILL.md`, in "Working beside a person", a second paragraph: "`orcamon focus JOB` puts a job in front of the person: their TUI jumps to it, with `--atoms 129,55` ringed and measured, `--point`, `--step` or `--cycle` for one geometry, and `--mode N` animated. Use it to show what you are talking about, and to show a proposed structure before anything runs: write the draft input beside the others and focus it — the TUI draws a job's input geometry before there is any output. It waits up to `--timeout` seconds (default 5) for the TUI to take the request; exit 4 means it did not, exit 2 that it could not show it." In "Other rules", replace the two sentences "Never kill, signal, edit or resubmit a running job on your own initiative. orcamon only reads." with D10's wording. README: the "Working beside an agent" section gains a paragraph on `focus` (what it shows, the 10 s a TUI keeps a request it cannot match yet, `--timeout`); the exit-code table's 4 row as in step 1; the `focus` reference section and refreshed first block. Reinstall the skill.
5. `test_orcamon.py`: `test_focus_waits_for_the_tui_to_take_it()`, registered right after `test_the_tui_tells_the_agent_what_it_shows`. No Textual: a fake TUI. Own `XDG_STATE_HOME` as in ST5; temp root holding `_water_job(root, "water")`; `fake = {"schema": 1, "pid": os.getpid(), "host": steer.host_name(), "root": str(root.resolve()), "updated_unix_s": time.time(), "job": None, "geometry": None, "mode": None, "selection": [], "measurement": None, "as_xyz": None, "focus": None}`; `ack(ok, message)` runs in a `threading.Thread`: waits up to 3 s for `steer.focus_path(os.getpid())`, `req = steer.take_focus(os.getpid())`, appends `req` to a list, and writes `fake` with `"focus": {"id": req["id"], "ok": ok, "message": message}` to `steer.view_path()`. Exactly 7 checks:
   - with no view file: `focus water` exits 2, stderr contains `no orcamon TUI is running`;
   - `fake` written, `ack(True, "showing water (test)")` started, `focus water --atoms 1,0 --timeout 3`: exit 0, stdout contains `showing water (test) in the TUI (pid {os.getpid()})`, and the taken request had `file == str((root / "water" / "job.inp").resolve())`, `atoms == [1, 0]`, `point is None`, `mode is None`;
   - `ack(False, "nope")`, `focus water --timeout 3`: exit 2, stderr contains `could not show it: nope`;
   - no ack, `focus water --timeout 0.3`: exit 4, stderr contains `has not taken the request within 0.3 s`, and `steer.focus_path(os.getpid()).exists()`;
   - `focus water --timeout 0`: exit 0, stdout contains `queued for the TUI`;
   - `focus water --atoms 0,9`, `--mode 3`, `--step 1 --cycle 1` and `--point 0`: each exits 2, with stderr containing respectively `atom 9 is out of range: the geometry has 3 atoms (0-2)`, `this job has no frequency block`, `focus takes one of --point, --step or --cycle`, `--point needs an IRC or NEB job`;
   - `fake` rewritten with `"root"` = a second temp directory: `focus water --timeout 0` exits 2 with `outside it` in stderr.
6. Run the first three verification commands.

**Acceptance:**
- `test_focus_waits_for_the_tui_to_take_it` passes its 7 checks; `test_the_skill_names_every_command`, `test_the_skill_reference_matches_the_parser` and `test_this_repository_has_the_current_skill` pass.
- `test_orcamon.py`: 238 `pass`; `test_monitor.py`: 231 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### ST7 — The TUI takes the agent's focus, and `b` goes back

**Why:** the person's half of "show me": the TUI jumps where the agent pointed, says the agent did it, and lets the person return (D1, D9).
**Depends on:** ST6.
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. `ConvergencePlot.select_point(self, point=None, step=None, cycle=None) -> bool`: the first `i` with (`point is not None` and `self._mode in ("irc", "neb")` and `p.index == point`) or (`step is not None` and `self._mode == "scan"` and `p.scan_step == step`) or (`cycle is not None` and `self._mode == "opt"` and `p.cycle == cycle`); none → False; else `self._following_latest = False`, `self._select(i)`, True.
2. `RotatableGeometryImage`: move the `else` branch of `_cycle_through` (the four statements of § Facts) into `_start_mode_at(self, modes, kind, position)` and call it there, so `i` and `I` behave byte for byte as before. Add `show_mode(self, index) -> bool`: for `kind, modes` in `("", imaginary_modes(state))`, `("real ", real_modes(state))`, when `index` is among the mode indices, `_start_mode_at(modes, kind, that position)`, `self._input.request()`, `self.app.update_detail()`, return True; False otherwise (and when `self._job` is None).
3. `MonitorApp`: `FOCUS_CHECK_S = 0.5` and `FOCUS_WAIT_S = 10.0` (module constants); `self._pending_focus = None` and `self._back = None` in `__init__`; `on_mount` adds `self.set_interval(FOCUS_CHECK_S, self._check_focus)`; binding `("b", "agent_back", "Back")` after `m`. `_select_label(label)`: `table.move_cursor(row=table.get_row_index(label))` (ignore `RowDoesNotExist` from `textual.widgets.data_table`), `self.selected_label = label`, `self.update_detail()`.
4. `_check_focus`: when nothing is pending, `doc = steer.take_focus(os.getpid())`; return on None; set `doc["_since"] = time.monotonic()` and keep it in `self._pending_focus`. `path = Path(str(doc.get("file", ""))).resolve()`. Not `path.is_relative_to(self.root)` → `_finish_focus(doc, False, f"{path} is outside the TUI's root {self.root}")`. `job` = the job in `self.jobs`, not in `self._gone`, whose `inp_path.resolve() == path`. When `job` is None or not `job.parsed`: after `FOCUS_WAIT_S` → `_finish_focus(doc, False, f"{path} not found under {self.root} within 10 s")`; else, once per request (`doc["_rescan"]`) when `job` is None, `self._force_discover = True` and `self.refresh_all()`; return. Then: record `doc["_back"] = (self.selected_label, geometry.selection, geometry.mode_index)` if not yet recorded; `self._select_label(job.label)`; if `geometry._job is not job`, return (the scan holds the job; the next check retries). `notes = []`; for `point`, `step`, `cycle` (named `point`, `scan step`, `cycle`) given in `doc`, when `chart.select_point(**{key: value})` is False append `f"{name} {value} is not on the chart, showing the newest"`; when `doc.get("atoms")`, `geometry.set_selection(parse_atom_list(" ".join(map(str, doc["atoms"])), len(geometry_shown(job, chart.selected_point())[0])))`, a `ValueError` appending `f"atoms not picked: {exc}"`; when `doc.get("mode") is not None` and `geometry.show_mode(doc["mode"])` is False, append `f"mode {doc['mode']} not shown: no such mode"`. `self._back = doc["_back"]`; `self.notify(f"agent: showing {job.label} (b goes back)", timeout=8)`; `_finish_focus(doc, True, "; ".join([f"showing {job.label}", *notes]))`.
5. `_finish_focus(doc, ok, message)`: `self._pending_focus = None`; `self._last_focus = {"id": doc.get("id"), "ok": ok, "message": message}`; when not ok, `self.notify(f"agent's focus not shown: {message}", severity="warning", timeout=8)`; `self.update_detail()` (which publishes it). `action_agent_back`: with `self._back` None, `self.notify("nothing to go back to")` and return; else `label, selection, mode = self._back`; `self._select_label(label)`; when the pane does not now show that label's job, `self.notify("the job is being read; press b again")` and return, keeping `_back`; else `self._back = None`, `geometry.set_selection(selection)`, and `geometry.show_mode(mode)` when `mode is not None`.
6. README `Keys:` bullet: "`b` returns to where you were before the agent's last `orcamon focus`". `SKILL.md`, at the end of the `focus` paragraph: "The person can press `b` to go back." Reinstall the skill.
7. `test_orcamon.py`: `test_the_tui_takes_the_agents_focus()`, registered right after `test_focus_waits_for_the_tui_to_take_it`. Own `XDG_STATE_HOME` as in ST5. Temp root holding `_irc_tree(root / "irc")` and `modes_all` (built as in `test_the_tui_cycles_real_modes`); wait until both are parsed and the pane shows `irc`. `cli = lambda *a: _orcamon(["--root", str(root), "--liveness", "mtime", "--no-cache", *a])`; `until(cond, s)` polls `await pilot.pause(0.1)` up to `s` seconds; `title = lambda: geometry.border_title`. Exactly 6 checks:
   - `cli("focus", "modes_all", "--atoms", "1,0,2", "--mode", "5", "--timeout", "0")` exits 0; within 3 s `app.selected_label == "modes_all"`, `geometry.selection == (1, 0, 2)`, `geometry.mode_index == 5`, and the title contains `mode 5 50.0 cm-1 (real 1/5)` and `angle 1-0-2: 104.5°`;
   - `steer.read_json(steer.view_path())["focus"]` has `ok is True` and `message == "showing modes_all"`;
   - `cli("focus", "irc", "--point", "1", "--timeout", "0")`; within 3 s `app.selected_label == "irc"`, the title ends with `· IRC forward 0`, and the view document's `geometry["args"] == ["--point", "1"]`;
   - `b`, `await pilot.pause(0.4)`: `selected_label == "modes_all"`, `geometry.selection == (1, 0, 2)`, `geometry.mode_index == 5`;
   - `b` again: `selected_label == "modes_all"` (nothing to go back to);
   - write `root / "draft" / "job.inp"` = `"! B97-3c Opt\n* xyz 0 1\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n*\n"`, then `cli("focus", "draft", "--atoms", "0,1", "--timeout", "0")` exits 0; within 8 s `selected_label == "draft"` and the title contains `input geometry` and `distance 0-1: 0.960 Å`.
8. Run the first three verification commands.

**If unsure:** if the draft is not shown within 8 s, report whether it was discovered (`[j.label for j in app.jobs]`) and parsed; do not lengthen `FOCUS_WAIT_S` or the test's wait on your own.

**Acceptance:**
- `test_the_tui_takes_the_agents_focus` passes its 6 checks; `test_the_tui_cycles_real_modes` passes unchanged (`i`, `I` and their labels).
- `test_orcamon.py`: 244 `pass`; `test_monitor.py`: 231 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### ST8 — A click picks the atom under it, then every gate

**Why:** finding atom 129 by its label in a crowded structure is slow; pointing at it is how a person means "this one" (D8).
**Depends on:** ST7.
**Files:** `tools/orcamon/src/orcamon/core/geometry.py`, `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/README.md`, `test_monitor.py`, `test_orcamon.py`.

**Do:**
1. `core/geometry.py` (pure Python): `PICK_RADIUS_CELLS = 1.5` and `screen_positions(atoms, view, frame) -> list[tuple[float, float, float]]`, the pixel renderer's projection without numpy: `mean` of the coordinates; `radius` = the largest distance from the mean; `pad` = `max(d_i + vdW_i) - radius` for `space-filling` (vdW from `VDW_RADII`, default `DEFAULT_VDW_RADIUS`), else 0.5; `center = mean + pan`; `right, up = camera_basis(elev, azim)`, `forward = camera_forward(elev, azim)`; `span = (radius + pad) / zoom`; `scale = min(W, H) / (2 * span)` (1.0 when `span` is 0); for each atom, `p = xyz - center`, `(W / 2 + p·right * scale, H / 2 - p·up * scale, p·forward)`. Docstring: it must stay equal to `tui.raster.project`, which `test_a_click_picks_the_atom_under_it` holds.
2. `pick_atom(atoms, view, frame, cells, click) -> int | None` (D8): None for no atoms or any non-positive size; `px, py = click[0] + 0.5, click[1] + 0.5`; for each atom (skipping hydrogens when `not view.show_hydrogens`), `d = hypot(sx / W * w - px, 2 * (sy / H * h - py))` with `(W, H) = frame`, `(w, h) = cells`; among `d <= PICK_RADIUS_CELLS` the smallest `(d, -depth)`. And `toggle_selection(selection, atom) -> tuple`: without `atom` when present; else `selection[1:] + (atom,)` when it already holds `MAX_PICKED`; else `selection + (atom,)`.
3. `RotatableGeometryImage._pick_frame(self, width, height)`: `(2 * width, 4 * height)` (the text pane's braille dots); `KittyGeometryImage` overrides with `self.FULL_PX`; `HerdrGeometryImage` with `(width * cells.width_px, height * cells.height_px)` from `herdr_graphics.cell_size()`, or `(width * 10, height * 20)` when that is None. `on_click(self, event)`: `offset = event.get_content_offset(self)`, return on None; `w, h = self.content_size.width, self.content_size.height`; `atom = pick_atom(self._atoms_for(self._job, self._point), self.view(), self._pick_frame(w, h), (w, h), (offset.x, offset.y))`; return on None; `self.set_selection(toggle_selection(self.selection, atom))`.
4. README `Keys:` bullet: "a click on an atom picks it and a second click drops it (a fifth drops the oldest)".
5. `test_monitor.py`: `test_a_click_picks_the_atom_under_it()`, registered right after `test_the_selection_is_ringed_and_labelled`. Exactly 6 checks, values from § Facts "Projection and picking":
   - `screen_positions` equals `raster.project`'s `(sx, sy, depth)` within 1e-9 for the 5-atom C/O/H/N/K set of § Facts in all four representations, at `View()` and at `View(elev=35, azim=110, zoom=1.7, pan=(0.3, -0.2, 0.1))`, `(900, 750)` — the atoms `[("C", 0.1, 0.2, 0.3), ("O", 1.3, -0.4, 0.2), ("H", -0.9, 0.8, -0.5), ("N", 0.4, 1.5, 1.1), ("K", -1.2, -1.0, 0.6)]`;
   - `five` under `View(elev=0, azim=0)` at `(100, 100)`: positions, each rounded to 6, `[(50.0, 50.0, 0.0), (90.0, 50.0, 0.0), (10.0, 50.0, 0.0), (50.0, 50.0, 1.0), (50.0, 50.0, -1.0)]`;
   - with `cells=(20, 10)`: clicks `(9, 4)`, `(17, 5)`, `(14, 5)`, `(1, 4)` give `3`, `1`, `None`, `2`;
   - `(1, 4)` under `View(elev=0, azim=0, show_hydrogens=False)` gives None;
   - `toggle_selection((), 3) == (3,)`, `((3,), 3) == ()`, `((0, 1, 2, 3), 4) == (1, 2, 3, 4)`, `((0, 1), 2) == (0, 1, 2)`;
   - `pick_atom([], …)` and `pick_atom(five, View(), (100, 100), (0, 10), (1, 1))` are None.
6. `test_orcamon.py`: `test_the_tui_picks_atoms_by_click()`, registered right after `test_the_tui_takes_the_agents_focus`. Temp root holding `_water_job(root, "water")`; headless; wait until parsed and shown. `w, h` = the pane's `content_size`; `atoms = geometry._atoms_for(geometry._job, geometry._point)`; `frame = geometry._pick_frame(w, h)`; `pos = screen_positions(atoms, geometry.view(), frame)`; `cell(i) = (int(pos[i][0] / frame[0] * w), int(pos[i][1] / frame[1] * h))`; `click(c)` = `await pilot.click("#geometry", offset=(c[0] + 1, c[1] + 1))` then `await pilot.pause(0.2)` (the `+ 1` is the border). Exactly 4 checks:
   - click `cell(2)`: `geometry.selection == (2,)`;
   - click `cell(2)` again: `()`;
   - click `cell(0)`, then `cell(1)`: `(0, 1)`, and the title contains `distance 0-1: 0.960 Å`;
   - `pick_atom(atoms, geometry.view(), frame, (w, h), (0, 0)) is None`, and a click at content `(0, 0)` leaves the selection `(0, 1)`.
7. Run all four verification commands.

**If unsure:** if `event.get_content_offset(self)` is None for a click inside the pane (§ Facts says it is not), stop and report the event's `x`, `y` and the pane's `region`.

**Acceptance:**
- `test_a_click_picks_the_atom_under_it` passes its 6 checks and `test_the_tui_picks_atoms_by_click` its 4.
- `test_monitor.py`: 237 `pass`; `test_orcamon.py`: 248 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`. `run_gates.py`: `38 gates`, `0 failed`.
- `git diff 91cb1d9 -- '*test*'` across the whole plan shows only added tests, fixtures, imports, the `XDG_STATE_HOME` lines and `__main__` entries; no existing assertion or expected value changed.

---

## Out of scope

1. Everything that acts on ORCA — the Act plan: `geom --rotate-bond` / `--set-dihedral` / `--set-distance` (decided 2026-10-01 by the user: a bond rotation moves the SMALLER fragment), `check`, `connectivity`, `run` and `stop` (decided 2026-10-01 by the user: local runs only at first; SLURM `submit` later, with its own liveness), and the rewrite of the skill's rule for running jobs.
2. Showing a bare `.xyz` file, side-by-side or overlaid proposals, and a diff of two geometries; a proposal is shown as a draft input (D9).
3. A TUI on another host, several TUIs addressed by name, and any network transport; the newest live TUI on this host is the one (D4).
4. A confirm-before-jump mode, more than one level of `b`, and an activity strip of the agent's commands in the TUI.
5. Dragging to rotate, rubber-band selection, and picking bonds rather than atoms.
6. Any change to `snapshot` (it does not draw a selection), `geom`, `freqs`, the chart keys, `i` and `I`.

## Progress template

```markdown
# orcamon — the person and the agent share one view — Progress

Plan: `PLAN_ORCAMON_STEER.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| ST1 measure and the command | TODO | | | |
| ST2 the state files | TODO | | | |
| ST3 renderers draw the selection | TODO | | | |
| ST4 s picks atoms | TODO | | | |
| ST5 the TUI publishes; orcamon view | TODO | | | |
| ST6 orcamon focus | TODO | | | |
| ST7 the TUI takes focus; b | TODO | | | |
| ST8 click picking, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 208 | 212 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
