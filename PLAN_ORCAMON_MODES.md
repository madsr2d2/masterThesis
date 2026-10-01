# orcamon — every normal mode, drawn and described (NM1–NM4)

**Written:** 2026-10-01. **Verified against:** `master` at `1c60096` (the reaction-path run, `PLAN_ORCAMON_PATHS.md`, is committed). The parser change in NM1 was prototyped on a scratch copy of this commit's package: the full `test_monitor.py` passed against it (193/193), and the real-tree sizes, counts and atom rankings below were measured on that prototype.
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_MODES_PROGRESS.md` (template at the end).

orcamon keeps only the imaginary columns of a `NORMAL MODES` block, so neither a person nor an agent can look at a real mode: not the 19 cm⁻¹ wobble that decides a ZPE, nor the 130 cm⁻¹ mode next to a TS's reaction coordinate. This plan keeps every mode of a job's FINAL Hessian, stored compactly. `snapshot --mode N` then draws any of them; a new `freqs --mode N` names the atoms that move most, so an agent can say in words what a mode does; and `I` in the TUI cycles the real modes lowest-frequency first, as `i` cycles the imaginary ones. On the way, the modes parser stops rebuilding the set of imaginary modes once per printed row, which makes the largest output in the tree parse about 2.8 times faster. **The goal a deviation must still meet:** every mode of a final Hessian can be drawn and described, using ORCA's own mode number; everything an imaginary mode does today — `state.modes`, `i`, the existing `snapshot --mode` behaviour and messages — still holds; and a Hessian from the middle of an optimization stays unanimated, as it is now.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 193 `pass` lines, last line `0 failure(s)` | ~2 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_modes/xdg .venv/bin/python test_orcamon.py` | 197 `pass` lines, last line `0 failure(s)` | ~12 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~45 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_modes/xdg .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (same for `test_orcamon.py`, with its environment variables). Each task states the counts it must end at.

Every task runs the first three commands. Only NM4 runs `run_gates.py`: every task touches only `tools/orcamon/`, `test_monitor.py`, `test_orcamon.py` and `.claude/skills/orcamon/SKILL.md`, which no other gate reads except the duplicate guard.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's; a new gate file would change the count and break that gate. New tests go in `test_monitor.py` or `test_orcamon.py`.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- a hand-derived expected value in this plan does not hold after the change is implemented as written (report the measured value; do not change the expected value or the design to make it pass);
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance (`test_normal_modes_*`, the `snapshot --mode` checks in `test_snapshot_writes_a_png`, and the `i` checks in `test_the_tui_runs_headless` pin today's imaginary-mode behaviour and must pass unchanged);
- a real-tree value differs from this plan's (report it, and say whether that job's files changed since `1c60096`);
- in NM4, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's. (`tools/orcamon/README.md` and `tools/orcamon/src/orcamon/SKILL.md` are orcamon's; `.claude/skills/orcamon/SKILL.md` changes ONLY through `orcamon skill install --project .`.)
- Never `git add -A`, `git add .` or `git commit -a`; stage the files the task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file.
- Never define the same top-level name in both `test_monitor.py` and `test_orcamon.py`, private or not: the duplicate guard reads private names too (found by the previous run). A fixture both files need lives in `test_monitor.py` and `test_orcamon.py` imports it (`from test_monitor import …`), as `_irc_tree` already is.
- Never run `orcamon` against the user's cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_modes/xdg`. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never write into `computational/`, never start ORCA, and never kill or signal a process the task did not start. Tests use the synthetic fixtures in § Facts; the real tree appears only in read-only acceptance commands.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `core` and `cli` import only the standard library at module level (`array` is standard library).
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill.
- Never change the key set of an existing `--json` document, and never change `REPORT_SCHEMA` (it stays `1`). `freqs --mode N --json` is a new document (D6).
- Never add a fenced code block containing an `orcamon …` line to `tools/orcamon/src/orcamon/SKILL.md`: `test_orcamon.test_the_skill_examples_run` executes every such line.

## Scope

**In scope:** keeping every mode of a final Hessian (`JobState.normal_modes`); hoisting the imaginary set out of the modes row loop; `vibrations.real_modes`, `find_mode`, `participation`, `available_text`; `snapshot --mode` for any mode; `freqs --mode N [--top K]`; the TUI's `I`; tests, README and the shipped skill.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | A person picks a real mode in the TUI with `I` (Shift+I): off → the lowest-frequency real mode → the next → … → the highest → off. `i` keeps cycling the imaginary modes exactly as today. Any mode by number stays available through `snapshot --mode N`, and later through the agent's focus command (Steer plan). No number prompt. | user, 2026-10-01 |
| D2 | `orcamon freqs JOB --mode N [--top K]` lists the atoms that move most in mode N — ORCA's atom index (the one the pixel pane labels), its element, and its displacement relative to the largest — so an agent can describe a mode in words. | user, 2026-10-01 |
| D3 | All modes are kept in a NEW field, `JobState.normal_modes`: every mode with a non-zero frequency of the last `NORMAL MODES` block, and ONLY when that block belongs to a FINAL frequency block (`freqs_final` true when the block opens); otherwise `None`. `JobState.modes` (the imaginary modes, list vectors, also kept for intermediate Hessians) is unchanged, because existing tests pin it and the `ts_hessian` story depends on it. | planner, 2026-10-01 |
| D4 | The vectors in `normal_modes` are `array('d')`, not lists: a 420-mode Hessian is 174k numbers, 1.4 MB as `array` against about 5.6 MB as a list of float objects, held in memory for every job the TUI watches. `NormalMode.vector` may therefore be a list (imaginary, `modes`) or an array (`normal_modes`); everything that reads it indexes it and takes its `len`, which both support. | planner, 2026-10-01 |
| D5 | Modes are always named by ORCA's own mode number (the `freqs` and `NORMAL MODES` column number), never by a rank. "Real" means `cm1 > 0`; zero-frequency modes (translations and rotations) are never kept, listed or animated. Real modes are ordered by frequency (ties by index) wherever they are cycled. | planner, 2026-10-01 |
| D6 | `freqs --mode N --json` prints `{"schema": 1, "job": …, "mode": {"index": N, "cm1": …, "atoms": [{"atom": i, "element": "H", "relative": 1.0}, …]}}`, with `relative` rounded to 4 decimals. Plain `freqs --json` is unchanged. | planner, 2026-10-01 |
| D7 | The "no such mode" message keeps the two substrings existing tests pin (`available: [0]`, `available: none`) by starting with the imaginary list: `available: [3] imaginary; real: 5 modes, 4-8`. | planner, 2026-10-01 |
| D8 | Commits go straight to `master`, one per task. | user, 2026-10-01 |

## Facts (verified 2026-10-01)

**Environment.** `/home/madsr2d2/masterThesis/.venv/bin/python` is Python 3.12.3 with Textual 8.2.8; orcamon is installed editable, console script `/home/madsr2d2/masterThesis/.venv/bin/orcamon`. `git status --porcelain` was empty at `1c60096`. In a Textual pilot, `await pilot.press("I")` fires a binding declared as `("I", …)` (probed; `"shift+i"` does not).

**The parser** (`tools/orcamon/src/orcamon/core/parser.py`):
- `NormalMode(index, cm1, vector)` dataclass; `JobState.frequencies` (every frequency of the last block, printed order, so a mode's position is ORCA's index), `JobState.modes` (imaginary modes of the last `NORMAL MODES` block, `vector` a list), `JobState.freqs_final`.
- In `_feed_blocks`, the frequency block closes on the first non-row, non-blank line after its rows — which, in ORCA's output, is the `NORMAL MODES` header itself (or the `------------` line before it) — and sets `frequencies`, `imaginary_freqs`, `freqs_final` and `freq_cycle` there. The `NORMAL MODES` header test comes AFTER that in the same call, so `freqs_final` is already set when the modes block opens. A new `VIBRATIONAL FREQUENCIES` header sets `self.modes = None`.
- The modes block: the header sets `_in_modes_block`, `_modes_columns = []`, `_modes_values = {}`, `_modes_seen_row = False`. Each row matched by `_MODE_ROW_RE` currently runs `imaginary = {i for i, f in enumerate(self.frequencies or []) if f < 0}` — once PER ROW, over every frequency — then stores the imaginary columns in `_modes_values[mode][coord]`. `_finish_modes` (called on the first non-row line after a row) keeps a mode only if every one of its `3N` coordinates was printed, where `3N == len(self.frequencies)`.
- The prototype of NM1 (exactly the steps NM1 lists) changed nothing else and passed all 193 `test_monitor.py` checks.

**Real-tree measurements** (`computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/`, read-only; parsed with `JobState` + `read_appended`, the pickle being `len(pickle.dumps(state))`):

| job | modes (zero) | final | imaginary | `normal_modes` after NM1 | pickle now → after | parse now → after |
|---|---|---|---|---|---|---|
| `geometry_r2scan3c-xtb/ts/verify_exact_hessian` | 420 (6) | yes | [6] | 414 modes, array of 420 | 24,911 → 1,435,618 | 0.39 → 0.15 s |
| `geometry_b973c-xtb/ts/optts_bridge` | 411 (6) | yes | [6, 7] | 405 modes, array of 411 | 469,000 → 1,828,247 | 0.80 → 0.33 s |
| `geometry_r2scan3c-xtb/pc` (67 MB, 9 Hessians) | 420 (6) | NO (intermediate) | [6 … 14] | `None` | 394,517 → 394,625 | 3.52 → 1.24 s |

Times are one machine's and vary with load (the same `pc` parse also took 4.7–7.1 s through the CLI); the pickle sizes and counts are deterministic.

**Real mode rankings** (top atoms by relative displacement, rounded to 2; the mode is drawn on the 137-atom structure through `sorted(qm_atom_indices)`, as `snapshot --mode` already does): `geometry_b973c-xtb/ts/optts_freq_tight` has 51 modes, 44 real, lowest real `7: 130.43`, `8: 147.17`, `9: 173.90`; mode 6 (`-369.03`) → `132 H 1.00`, `128 C 0.70`, `133 O 0.50`, `125 H 0.43`, `135 H 0.42`, `122 H 0.40`, `131 O 0.31`, `126 H 0.27`; mode 7 (`130.43`) → `130 H 1.00`, `131 O 0.45`, `134 O 0.41`, `132 H 0.26`, `135 H 0.24`, `121 H 0.21`, `129 O 0.19`, `133 O 0.19`.

**Vibrations** (`core/vibrations.py`): `imaginary_modes(state)` (`[]` unless `state.modes` and `state.freqs_final`; most negative first), `target_indices`, `alternate_geometries(job)`, `mode_geometry(candidates, qm_atom_indices, mode) -> (atoms, qm) | None`, `offsets(atoms, qm, mode, amplitude=0.20) -> list | None` (None for a zero vector), `displaced`, `phase_sine`.

**The CLI** (`cli/commands.py`, `cli/__init__.py`):
- `cmd_snapshot` with `--mode`: raises `UsageError("this job has no frequency block (no Freq or NumFreq step)")` when `state.frequencies is None`; then looks the mode up in `vibrations.imaginary_modes(state)` and otherwise raises `UsageError(f"no imaginary mode {args.mode} in this job; available: {available or 'none'}")`; then `mode_geometry` / `offsets` (refusing `f"no geometry in this job holds mode {N} ({k} atoms)"`); it prints `f"{output}: {W}x{H}, {n} atoms, {source}, mode {N} at {phase:g} deg"`.
- `cmd_freqs` prints `"{n} modes ({zero} zero, {k} imaginary)"`, the imaginary modes, then `lowest {k} real:`; its JSON keys are `job, n_modes, imaginary, lowest_real`. `build_parser` declares `freqs` with `--lowest K` only.
- Existing checks that must keep passing: `snapshot mode_ts --mode 999` exits 2 with `available: [0]` in stderr; `snapshot opt_done --mode 0` exits 2 with `available: none`; `snapshot opt_maxiter --mode 0` exits 2 with `no frequency block`; `snapshot mode_ts --mode 0` prints `… 3 atoms, cycle 2, mode 0 at 90 deg`. (`mode_ts` has frequencies `[-100, 10, 20, …, 80]` with only mode 0 non-zero; `opt_done` has frequencies `[0, 0, 1600, 3700, 3800]` and no `NORMAL MODES` block.)

**The TUI** (`tui/app.py`): `GEOMETRY_BINDINGS` holds `("i", "cycle_modes", "Modes")`; `RotatableGeometryImage.action_cycle_modes` builds `modes = vibrations.imaginary_modes(...)`, notifies `"this job computes no frequencies (no Freq or NumFreq step)"` or `"no imaginary mode in this job"` when empty, steps `position` through `indices` (off after the last), and sets `self.mode_index`, `self._mode_alternates = vibrations.alternate_geometries(self._job)` and `self.mode_label = f"mode {i} {cm1:.1f} cm-1 ({position + 1}/{len(indices)})"`, then `self._input.request()` and `self.app.update_detail()`. `_apply_mode(atoms, job)` looks the mode up with `next((m for m in vibrations.imaginary_modes(job.state) if m.index == self.mode_index), None)`. The title shows `mode_label` (`_geometry_title`). `test_the_tui_runs_headless` presses `i`, `i` on `mode_ts` and `i` on `opt_maxiter`.

**Test helpers.** `test_orcamon.py`: `check`, `_orcamon(argv)`, `_freqs(values)`, `_normal_modes(vectors)` (`{mode index: 3N values}`), `_opt(cycles, energies, …)` (prints coordinates `_coords(n / 10)` per cycle: at cycle 2, `O (0.2, 0, 0)`, `H (1.16, 0, 0)`, `H (-0.04, 0.93, 0)`), `_CONVERGED`, `_DONE`, `_OPT_FREQ`; `test_the_tui_scrubs_a_path` shows the headless pattern on a private temp root. `test_monitor.py`: `check`, `_feed(state, text)` (line by line). Counts: `test_monitor.py` 193, `test_orcamon.py` 197.

**The synthetic modes** (`_MODE_FREQS`, `_MODE_VECTORS` and `_modes_block` in `test_monitor.py`; `test_orcamon.py` imports `_MODE_FREQS` and `_MODE_VECTORS`). Three atoms (O, H, H), so 9 modes:
```python
_MODE_FREQS = [0.0, 0.0, 0.0, -150.0, 100.0, 50.0, 300.0, 200.0, 400.0]
_MODE_VECTORS = {
    0: [0.0] * 9, 1: [0.0] * 9, 2: [0.0] * 9,
    3: [0.1, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.3],
    4: [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0],
    5: [0.0, 0.0, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0],
    6: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0],
    7: [0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    8: [0.0, 0.1, 0.0, 0.0, 0.1, 0.0, 0.0, 0.1, 0.0],
}


def _modes_block(freqs, vectors):
    """A VIBRATIONAL FREQUENCIES block and its NORMAL MODES block as ORCA
    prints them: six mode columns per sub-block, every coordinate repeated."""
    lines = ["VIBRATIONAL FREQUENCIES", "-----------------------", "",
             "Scaling factor for frequencies =  1.000000000  (already applied!)", ""]
    for i, f in enumerate(freqs):
        lines.append(f"{i:5d}:  {f:8.2f} cm**-1" + ("  ***imaginary mode***" if f < 0 else ""))
    lines += ["", "NORMAL MODES", "------------", "",
              "These modes are the Cartesian displacements weighted by the diagonal matrix",
              "M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom",
              "Thus, these vectors are normalized but *not* orthogonal", ""]
    n = len(freqs)
    for start in range(0, n, 6):
        cols = list(range(start, min(start + 6, n)))
        lines.append("            " + "".join(f"{c:11d}" for c in cols))
        for coord in range(n):
            lines.append(f"{coord:6d}     " + "".join(f"{vectors[c][coord]:11.6f}" for c in cols))
    lines += ["", ""]
    return "\n".join(lines) + "\n"
```
`_modes_block(_MODE_FREQS, _MODE_VECTORS)` was fed through the prototype: `normal_modes` indices `[3, 4, 5, 6, 7, 8]`, `modes` `[3]` with a list vector, `feed_text` identical to the line path; preceded by a `GEOMETRY OPTIMIZATION CYCLE   3` banner, `freqs_final` is False, `normal_modes` None and `modes` still `[3]`.

**Hand-derived values for the synthetic modes** (atoms `O, H, H` at indices 0, 1, 2; per-atom norm of each mode's (x, y, z) triple):
- real modes lowest first: `[5, 4, 7, 6, 8]` (50, 100, 200, 300, 400 cm⁻¹); by index `4-8`; imaginary `[3]`;
- `available_text` → `available: [3] imaginary; real: 5 modes, 4-8`;
- participation of mode 3: norms 0.1, 0.2, 0.3 → `[(2, "H", 1.0), (1, "H", 0.6667), (0, "O", 0.3333)]` (rounded 4);
- mode 5: norms 0.4, 0.3, 0 → `[(0, "O", 1.0), (1, "H", 0.75)]` (a zero-amplitude atom is left out);
- mode 8: norms 0.1, 0.1, 0.1 → `[(0, "O", 1.0), (1, "H", 1.0), (2, "H", 1.0)]` (ties by atom index);
- mode 5 drawn at full amplitude (`sine = 1.0`) on the cycle-2 coordinates: the scale is `0.20 / 0.4 = 0.5`, so O moves `(0, 0, +0.2)` to `(0.2, 0.0, 0.2)` and the first H moves `(+0.15, 0, 0)` to `(1.31, 0.0, 0.0)`.

## Conventions

- Match the surrounding code: one-line comments where the reason is not obvious; docstrings say WHY.
- New tests go in the file the task names, use `check(...)` for every assertion, build trees under `tempfile.TemporaryDirectory()`, and are added to that file's `__main__` list right after the test named in the task.
- Expected values are the hand-derived ones in § Facts, never values copied from a run. Floats are compared after `round(…, 4)` unless § Facts gives them exactly.
- No network access.

---

# Phase 1 — keep every mode

### NM1 — Keep every mode of a final Hessian, compactly, and stop rebuilding the imaginary set per row

**Why:** the parser drops every real column of `NORMAL MODES`, and rebuilds the imaginary set once per printed row: a 420-mode block has 28,980 rows, so `pc/job.out`'s nine blocks cost about a hundred million set operations (§ Facts, real-tree table).
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/parser.py`, `test_monitor.py`.

**Do:**
1. `core/parser.py`: `from array import array`. Add the `JobState` field `normal_modes: list | None = None` after `modes`, commented: every mode with a non-zero frequency of the last `NORMAL MODES` block, kept only when that block is the job's final Hessian (`freqs_final`); vectors are `array('d')` because a 420-mode Hessian is 174k numbers (D3, D4). Update `NormalMode`'s docstring to say `vector` is a list for `modes` and an `array('d')` for `normal_modes`.
2. Private fields beside `_modes_seen_row`: `_modes_imaginary: frozenset = frozenset()`, `_modes_keep_all: bool = False`, `_modes_all: dict = field(default_factory=dict)`, `_modes_filled: dict = field(default_factory=dict)`.
3. Where the `VIBRATIONAL FREQUENCIES` header sets `self.modes = None`, also set `self.normal_modes = None`.
4. When the `NORMAL MODES` header opens the block, also set `self._modes_imaginary = frozenset(i for i, f in enumerate(self.frequencies or []) if f < 0)`, `self._modes_keep_all = bool(self.freqs_final)`, `self._modes_all = {}`, `self._modes_filled = {}`.
5. In the row branch, delete the per-row `imaginary = {…}` line and test `mode_index in self._modes_imaginary` instead. In the same column loop, when `self._modes_keep_all and coord < n3 and mode_index < n3 and freqs[mode_index] != 0.0` (with `freqs = self.frequencies or []`, `n3 = len(freqs)`, taken once per row before the loop): create `self._modes_all[mode_index] = array("d", bytes(8 * n3))` and `self._modes_filled[mode_index] = 0` on first sight, then store the value at `[coord]` and add 1 to `_modes_filled[mode_index]`.
6. In `_finish_modes`, after `self.modes = modes or None`: `normal_modes` = the `NormalMode(index=i, cm1=self.frequencies[i], vector=self._modes_all[i])` for `i` in `sorted(self._modes_all)` whose `_modes_filled[i] == 3 * coords` (and `coords` non-zero), or `None` when there are none; then reset `_modes_all` and `_modes_filled` to `{}`. Update its docstring: imaginary modes go to `modes`; every non-zero mode of a final block to `normal_modes`; a partly printed mode to neither.
7. `test_monitor.py`: add `_MODE_FREQS`, `_MODE_VECTORS` and `_modes_block` exactly as § Facts gives them, and `test_every_mode_of_a_final_block_is_kept()`, registered right after `test_normal_modes_take_the_last_block`. Exactly 7 `check` calls, each on a fresh `JobState(path=Path("/nonexistent"))`:
   - fed (`_feed`) `_modes_block(_MODE_FREQS, _MODE_VECTORS)`: `[m.index for m in state.normal_modes] == [3, 4, 5, 6, 7, 8]`;
   - same state: `state.normal_modes[2].vector == array("d", [0.0, 0.0, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0])` and its `typecode == "d"`;
   - same state: `[m.index for m in state.modes] == [3]` and `isinstance(state.modes[0].vector, list)`;
   - fed `"         *                GEOMETRY OPTIMIZATION CYCLE   3            *"` and then the block: `state.freqs_final is False and state.normal_modes is None and [m.index for m in state.modes] == [3]`;
   - fed the block, then a second `_modes_block` whose mode 4 vector is `[0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0]` (the rest unchanged): `state.normal_modes[1].vector[4] == 0.7`;
   - `state.feed_text(<the block>)` gives the same indices and equal vectors as the line-fed state;
   - fed the block's text cut right after its first sub-block's row `     4 …` and then one blank line: `state.normal_modes is None`.
8. Run the first three verification commands, and the real-tree measurement: `cd /home/madsr2d2/masterThesis && .venv/bin/python -c "import pickle, time; from pathlib import Path; from orcamon.core.parser import JobState, read_appended
base = Path('computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP')
for rel in ['geometry_r2scan3c-xtb/ts/verify_exact_hessian', 'geometry_b973c-xtb/ts/optts_bridge', 'geometry_r2scan3c-xtb/pc']:
    p = base / rel; t = time.perf_counter(); s = JobState(path=p, stem='job'); read_appended(s, p / 'job.out'); nm = s.normal_modes
    print(rel.split('/')[-1], round(time.perf_counter() - t, 2), None if nm is None else (len(nm), nm[0].index, nm[0].vector.typecode, len(nm[0].vector)), len(pickle.dumps(s)), [m.index for m in s.modes or []])"`.

**Acceptance:**
- The real-tree measurement prints, apart from the times: `verify_exact_hessian … (414, 6, 'd', 420)` with a pickle between 1,400,000 and 1,480,000 and modes `[6]`; `optts_bridge … (405, 6, 'd', 411)`, pickle 1,780,000–1,880,000, modes `[6, 7]`; `pc … None`, pickle 394,000–395,500, modes `[6, 7, 8, 9, 10, 11, 12, 13, 14]`. Report the three times (expected about 0.15, 0.33 and 1.2 s on an idle machine; not a pass condition).
- The diff has no set comprehension over `self.frequencies` inside the row branch.
- `test_every_mode_of_a_final_block_is_kept` passes its 7 checks; `test_normal_modes_are_parsed`, `test_normal_modes_keep_only_imaginary`, `test_normal_modes_take_the_last_block` and `test_a_cached_state_equals_a_fresh_parse` pass unchanged.
- `test_monitor.py`: 200 `pass`; `test_orcamon.py`: 197 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### NM2 — Find any mode, and rank the atoms it moves

**Why:** the CLI and the TUI need one way to find a mode by ORCA's number — imaginary or real — one message listing what exists, and one ranking of the atoms a mode moves. Written once here, so `snapshot`, `freqs` and the TUI cannot disagree.
**Depends on:** NM1.
**Files:** `tools/orcamon/src/orcamon/core/vibrations.py`, `test_monitor.py`.

**Do:**
1. Add to `core/vibrations.py`, after `imaginary_modes`:
   ```python
   def real_modes(state) -> list:
       """The real modes of the FINAL Hessian, lowest frequency first; [] until one."""
       if state.normal_modes is None or not state.freqs_final:
           return []
       return sorted((m for m in state.normal_modes if m.cm1 > 0), key=lambda m: (m.cm1, m.index))


   def find_mode(state, index: int):
       """ORCA's mode `index`, imaginary or real, or None."""
       return next((m for m in imaginary_modes(state) + real_modes(state) if m.index == index), None)


   def available_text(state) -> str:
       """What a 'no such mode' message offers instead."""
       imaginary = [m.index for m in imaginary_modes(state)]
       real = sorted(m.index for m in real_modes(state))
       text = f"available: {imaginary or 'none'} imaginary"
       if real:
           text += f"; real: {len(real)} modes, {real[0]}-{real[-1]}"
       return text
   ```
2. Add `participation(atoms, qm_atom_indices, mode, top=None) -> list | None`: `None` when `target_indices(len(atoms), len(mode.vector) // 3, qm_atom_indices)` is None or every row's norm is 0; else, for each Hessian row `j` mapped to atom `indices[j]`, the norm `sqrt(x² + y² + z²)` of its three values divided by the largest norm; drop entries whose relative is exactly 0.0; sort by `(-relative, atom)`; keep the first `top` when `top` is given; return `[(atom, element, relative), …]` with `element = atoms[atom][0]` and `relative` unrounded. Docstring: what it ranks (the printed displacement pattern, the same one `offsets` scales) and that the atom number is ORCA's index in the drawn structure.
3. Update the module docstring's first paragraph to say every mode of a final Hessian can be drawn, not only the imaginary ones.
4. `test_monitor.py`: add `test_modes_are_found_and_described()`, registered right after `test_every_mode_of_a_final_block_is_kept`, on a state fed `_modes_block(_MODE_FREQS, _MODE_VECTORS)`, with `atoms3 = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0)]` and `rnd = lambda rows: [(a, e, round(r, 4)) for a, e, r in rows]`. Exactly 8 `check` calls, against § Facts "Hand-derived values for the synthetic modes":
   - `[m.index for m in real_modes(state)] == [5, 4, 7, 6, 8]`;
   - `find_mode(state, 3).cm1 == -150.0 and find_mode(state, 5).cm1 == 50.0 and find_mode(state, 0) is None and find_mode(state, 99) is None`;
   - `rnd(participation(atoms3, None, find_mode(state, 3))) == [(2, "H", 1.0), (1, "H", 0.6667), (0, "O", 0.3333)]`;
   - mode 5 → `[(0, "O", 1.0), (1, "H", 0.75)]`;
   - mode 8 → `[(0, "O", 1.0), (1, "H", 1.0), (2, "H", 1.0)]`;
   - mode 3 with `top=1` → `[(2, "H", 1.0)]`;
   - `available_text(state) == "available: [3] imaginary; real: 5 modes, 4-8"`;
   - a state fed the cycle-3 banner then the block: `real_modes(state) == [] and find_mode(state, 5) is None`.
5. Run the first three verification commands.

**Acceptance:**
- `test_modes_are_found_and_described` passes its 8 checks; every existing `test_mode_offsets_*` and `test_the_qm_map_needs_indices_that_fit` passes unchanged.
- `test_monitor.py`: 208 `pass`; `test_orcamon.py`: 197 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

# Phase 2 — use them

### NM3 — `snapshot --mode` for any mode, and `freqs --mode N` for words

**Why:** today `snapshot --mode 7` on a TS refuses a real mode, and an agent asked "what does mode 6 do" can only guess from an animation it cannot see. `freqs --mode N` gives it the atoms (D2).
**Depends on:** NM2.
**Files:** `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/README.md`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `test_orcamon.py`.

**Do:**
1. `cli/__init__.py`: change `snapshot`'s `--mode` help to `"animate a normal mode, imaginary or real: ORCA's mode number (the one `orcamon freqs` prints)"`; change `freqs`'s description to `"imaginary and lowest real vibrational frequencies, or the atoms one mode moves"`, and add to `freqs` `--mode N` (`type=int`, help `"describe one mode, imaginary or real: the atoms it moves most, by ORCA's atom number"`) and `--top K` (`type=int`, `default=None`, help `"atoms to list with --mode (default: 8)"` written literally, since the parser default is None so step 3 can tell an explicit `--top` apart).
2. `cmd_snapshot`: keep the `frequencies is None` refusal first and unchanged; replace the imaginary-only lookup with `mode = vibrations.find_mode(state, args.mode)` and, when it is None, `raise UsageError(f"no mode {args.mode} in this job; {vibrations.available_text(state)}")`. The rest (geometry choice, offsets, output line) is unchanged.
3. `cmd_freqs`: when `args.mode is None`, raise `UsageError("--top needs --mode")` if `args.top is not None`, else behave exactly as today. In the mode branch, `top = 8 if args.top is None else args.top`. When `args.mode` is set: refuse a missing frequency block with the same message `cmd_snapshot` uses; look the mode up with `find_mode` (refuse as in step 2); choose the structure as `cmd_snapshot` does (`[latest printed atoms or the file geometry] + vibrations.alternate_geometries(job)`, `mode_geometry`, refusing `f"no geometry in this job holds mode {N} ({k} atoms)"`); rank with `participation(atoms, qm, mode, top)`. Text: `f"mode {mode.index}  {mode.cm1:.2f} cm-1 ({'imaginary' if mode.cm1 < 0 else 'real'}), drawn on {len(atoms)} atoms"`, then `f"{'atom':>6}  {'el':<2}  {'relative':>8}"`, then one `f"{atom:>6}  {element:<2}  {relative:>8.2f}"` per atom, then `relative: each atom's displacement in this mode over the largest one's`. JSON: D6, `relative` rounded to 4. Exit through `_job_exit(job)`.
4. `tools/orcamon/src/orcamon/SKILL.md`, "Looking closer", one prose paragraph (no fenced block): "`orcamon freqs JOB --mode N` names the atoms that move most in mode N — ORCA's mode number, imaginary or real — each with its displacement relative to the largest (`--top K`, default 8). Use it to say in words what a mode does. `orcamon snapshot JOB --mode N -o file.png` draws any of those modes for a person." `tools/orcamon/README.md`: in "Saving a picture" change "an imaginary normal mode" to "any normal mode, imaginary or real"; replace the pasted `orcamon freqs` and `orcamon snapshot` `--help` blocks with the new `--help` output (at the width the other blocks use); add one sentence after "Saving a picture": "`orcamon freqs JOB --mode N` lists the atoms mode N moves most, by ORCA's atom number, with `--json`." Then run `cd /home/madsr2d2/masterThesis && .venv/bin/orcamon skill install --project /home/madsr2d2/masterThesis` (with `--force` only if it refuses).
5. `test_orcamon.py`: `from test_monitor import _MODE_FREQS, _MODE_VECTORS` inside the new test, and add `test_any_mode_is_drawn_and_described()`, registered right after `test_an_neb_is_listed_and_drawn`. Temp root with job `modes_all`: `job.inp` = `_OPT_FREQ`, `job.out` = `_opt(2, [-80.00, -80.10]) + _CONVERGED + "\n" + _freqs(_MODE_FREQS) + _normal_modes(_MODE_VECTORS) + _DONE + "\n"`. Every call passes `["--root", str(root), "--liveness", "mtime", "--no-cache", …]`. Exactly 8 `check` calls:
   - `snapshot modes_all --mode 5 -o <tmp>/p.png --size 120x100`: exit 0, the file starts `b"\x89PNG"`, stdout ends with `3 atoms, cycle 2, mode 5 at 90 deg`;
   - `snapshot modes_all --mode 0 -o <tmp>/q.png`: exit 2, stderr contains `available: [3] imaginary; real: 5 modes, 4-8`;
   - `freqs modes_all --mode 3`: exit 0; the first line is `mode 3  -150.00 cm-1 (imaginary), drawn on 3 atoms`, and the three lines after the `  atom  el  relative` header are `     2  H       1.00`, `     1  H       0.67`, `     0  O       0.33`;
   - `freqs modes_all --mode 3 --json`: `doc["mode"] == {"index": 3, "cm1": -150.0, "atoms": [{"atom": 2, "element": "H", "relative": 1.0}, {"atom": 1, "element": "H", "relative": 0.6667}, {"atom": 0, "element": "O", "relative": 0.3333}]}`;
   - `freqs modes_all --mode 5 --top 1`: exactly one atom line, `     0  O       1.00`;
   - `freqs modes_all --mode 0`: exit 2, stderr contains `available: [3] imaginary`;
   - `freqs modes_all --top 3`: exit 2, stderr contains `--top needs --mode`;
   - `freqs modes_all` (no flags): stdout contains `9 modes (3 zero, 1 imaginary)` and `lowest 5 real:` (unchanged behaviour).
6. Run the first three verification commands, and on the real tree (`cd /home/madsr2d2/masterThesis`, `XDG_CACHE_HOME=/tmp/orcamon_modes/xdg`):
   - `.venv/bin/orcamon --root computational freqs geometry_b973c-xtb/ts/optts_freq_tight --mode 6 --json | .venv/bin/python -c "import json,sys; m=json.load(sys.stdin)['mode']; print(m['index'], m['cm1'], [(a['atom'], a['element'], round(a['relative'], 2)) for a in m['atoms']])"` and the same with `--mode 7`;
   - `.venv/bin/orcamon --root computational snapshot geometry_b973c-xtb/ts/optts_freq_tight --mode 7 -o /tmp/orcamon_modes/mode7.png --size 300x250`.

**Acceptance:**
- Mode 6 prints `6 -369.03 [(132, 'H', 1.0), (128, 'C', 0.7), (133, 'O', 0.5), (125, 'H', 0.43), (135, 'H', 0.42), (122, 'H', 0.4), (131, 'O', 0.31), (126, 'H', 0.27)]`; mode 7 prints `7 130.43 [(130, 'H', 1.0), (131, 'O', 0.45), (134, 'O', 0.41), (132, 'H', 0.26), (135, 'H', 0.24), (121, 'H', 0.21), (129, 'O', 0.19), (133, 'O', 0.19)]`; the snapshot exits 0 and reports `137 atoms` and `mode 7 at 90 deg`.
- `test_any_mode_is_drawn_and_described` passes its 8 checks; `test_snapshot_writes_a_png` passes unchanged (`available: [0]`, `available: none`, `no frequency block`, `mode 0 at 90 deg`).
- `test_orcamon.py`: 205 `pass` (including the skill tests); `test_monitor.py`: 208 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### NM4 — `I` cycles the real modes in the TUI, then every gate

**Why:** the person following the agent should be able to look at a real mode without leaving the terminal (D1); `i` already does this for the imaginary ones.
**Depends on:** NM3.
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. `GEOMETRY_BINDINGS`: add `("I", "cycle_real_modes", "Real modes")` right after `("i", "cycle_modes", "Modes")`.
2. Move the body of `action_cycle_modes` after its empty-list check into `_cycle_through(self, modes, kind)` — `kind` is `""` for imaginary and `"real "` for real — which steps `position` exactly as today and sets `mode_label = f"mode {i} {cm1:.1f} cm-1 ({kind}{position + 1}/{len(indices)})"`. `action_cycle_modes` keeps its docstring, its two notifications and calls `_cycle_through(modes, "")`, so its label stays `mode 0 -100.0 cm-1 (1/1)` byte for byte. Add `action_cycle_real_modes` (docstring: `I`: off -> the lowest real mode -> ... -> the highest -> off): `modes = vibrations.real_modes(self._job.state) if self._job is not None else []`; when empty notify `"this job computes no frequencies (no Freq or NumFreq step)"` if `self._job is not None and self._job.state.frequencies is None`, else `"no real mode with a displacement pattern in this job"`, and return; else `_cycle_through(modes, "real ")`. Pressing `I` while an imaginary mode runs (or `i` while a real one does) starts that kind's cycle at its first mode, because the running index is not in the other list.
3. `_apply_mode`: look the mode up with `vibrations.find_mode(job.state, self.mode_index)` instead of searching `imaginary_modes`.
4. `tools/orcamon/README.md`: where the TUI keys say `i` animates the job's imaginary normal modes, add that `I` steps through its real modes, lowest frequency first.
5. `test_orcamon.py`: add `test_the_tui_cycles_real_modes()`, registered right after `test_the_tui_scrubs_a_path`. Temp root holding only `modes_all` (built as in NM3; import `_MODE_FREQS`, `_MODE_VECTORS` from `test_monitor`); `MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")`, headless at `(180, 50)` as `test_the_tui_scrubs_a_path` does; wait until the job is parsed and the geometry pane shows it; `geometry = app.query_one("#geometry")`; `geometry.focus()`; after each key `await pilot.pause(0.2)`. Exactly 7 `check` calls:
   - `I`: `geometry.mode_index == 5` and the title contains `mode 5 50.0 cm-1 (real 1/5)`;
   - `I`: `mode_index == 4` and the title contains `(real 2/5)`;
   - `i`: `mode_index == 3` and the title contains `mode 3 -150.0 cm-1 (1/1)`;
   - `I`: `mode_index == 5` (the real cycle starts again at its lowest);
   - `I` four more times: `mode_index == 8` and the title contains `mode 8 400.0 cm-1 (real 5/5)`;
   - `I`: `mode_index is None` and `"cm-1" not in` the title;
   - with `mode_index` set back to 5 by pressing `I` once, set `geometry._mode_s = 1.0` and, with no `await` in between (so the animation timer cannot overwrite it), call `geometry._apply_mode(list(app.selected_job().state.atoms), app.selected_job())`: its atom 0 is `("O", 0.2, 0.0, 0.2)` and atom 1 `("H", 1.31, 0.0, 0.0)`, each coordinate compared after `round(…, 6)`.
6. Run all four verification commands.

**If unsure:** if `pilot.press("I")` does not reach the pane's binding (§ Facts says it does on Textual 8.2.8), stop and report the key event Textual delivered instead; do not change the binding to another key on your own.

**Acceptance:**
- `test_the_tui_cycles_real_modes` passes its 7 checks; `test_the_tui_runs_headless` passes unchanged (its `i` checks included).
- `test_orcamon.py`: 212 `pass`; `test_monitor.py`: 208 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`. `run_gates.py`: `38 gates`, `0 failed`.
- `git diff 1c60096 -- '*test*'` across the whole plan shows only added tests, fixtures and `__main__` lines; no existing assertion or expected value changed.

---

## Out of scope

1. A typed mode-number prompt in the TUI (D1), and the agent-driven `focus --mode N` (the Steer plan).
2. Modes of an intermediate Hessian (mid-optimization): still not animated, listed by `freqs --mode`, or kept in `normal_modes` (D3).
3. Arrows, amplitude or speed controls, side-by-side modes, or IR intensities and thermochemistry from the output.
4. Internal-coordinate descriptions of a mode (which bonds stretch, which angles bend) — the atom ranking is the description this plan ships.
5. Any change to `JobState.modes`, `imaginary_modes`, the `i` key's behaviour, or the existing `freqs` output without `--mode`.
6. Reading `.hess` files or `orca_pltvib` output.

## Progress template

```markdown
# orcamon — every normal mode, drawn and described — Progress

Plan: `PLAN_ORCAMON_MODES.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| NM1 Keep every mode, compactly | TODO | | | |
| NM2 Find, list and rank modes | TODO | | | |
| NM3 snapshot/freqs --mode for any mode | TODO | | | |
| NM4 I in the TUI, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 193 | 197 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
