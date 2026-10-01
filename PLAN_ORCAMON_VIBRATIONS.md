# orcamon — animate imaginary normal modes (MO1–MO5)

**Written:** 2026-10-01. **Verified against:** `master` at `a77e316` (the last commit; `tools/orcamon` is at `112397c`).
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_VIBRATIONS_PROGRESS.md` (template at the end).

orcamon reads ORCA jobs but throws away the normal-mode displacement vectors, so a transition state's reaction coordinate cannot be seen. This plan parses the `NORMAL MODES` block an ORCA frequency job already prints, maps a mode onto the geometry the pane is showing, and animates it: `i` in the TUI cycles through the imaginary modes (most negative first, then off), and `snapshot --mode N` renders one phase of one mode to a PNG. **The goal a deviation must still meet:** with no mode selected nothing about the current rendering changes; with a mode selected the geometry oscillates along that mode's own displacement pattern at a fixed peak amplitude, and a mode that cannot be mapped to the displayed structure is refused with a message rather than drawn wrongly.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 136 `pass` lines, last line `0 failure(s)` | ~1 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python test_orcamon.py` | 157 `pass` lines, last line `0 failure(s)` | ~6 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~17 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (same for `test_orcamon.py`).

Every task runs the first three commands. Only MO5 runs `run_gates.py`. MO1–MO4 touch only `tools/orcamon/`, `test_monitor.py` and `test_orcamon.py`.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's and must not be edited; a new gate file would change the count and break that gate. All new tests go in `test_monitor.py` or `test_orcamon.py`, which are already gates.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- a hand-derived value written in this plan does not hold after the algorithm is implemented as written (report the measured value; do not change the expected value, the tolerance, or the algorithm to make it pass);
- the MO2 ground-truth probe does not reproduce `verify_exact_hessian/job.hess.v006.xyz` (report the ratio spread; the mass-weighting or the vector convention is then wrong);
- the MO2 mapping probe shows that `sorted(qm_atom_indices)` does not reproduce `optts_freq_tight/job.activeRegion.xyz`;
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance;
- in MO5, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's.
- Never `git add -A`, `git add .` or `git commit -a`; stage the files a task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file (see the note above).
- Never run `orcamon` against the user's cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg`. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never run `orca_pltvib` or any other ORCA tool, and never write into `computational/`: the mode data is read from the `.out` only. Real jobs appear only in the MO2 probes and the MO5 pictures, all outside the repository.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `core` and `cli` import only the standard library at module level (a test enforces this).
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill. Tests use synthetic atoms only.

## Scope

**In scope:** parsing ORCA's `NORMAL MODES` block; a `core/vibrations.py` that maps a mode onto the displayed geometry; `snapshot --mode N [--phase DEG]`; the TUI `i` key, its animation timer and its pane-title mark; tests; README and installed-skill updates.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | The mode data comes from the `NORMAL MODES` block already in the `.out`, not from `orca_pltvib` or the `.hess` (orcamon must work on copied trees, live jobs and hosts with no ORCA install, and parsing gives every mode). | user (chose the source), 2026-10-01 |
| D2 | Only modes with `frequency < 0` are kept — the same rule `cmd_freqs` uses, with no cutoff. | user, 2026-10-01 |
| D3 | The displacement is the printed vector scaled so its largest atom moves `MODE_AMPLITUDE_ANGSTROM = 0.20` Å. **No mass factor**: `orca_pltvib` writes `2.0 ×` the printed vector for both C and H (measured, MO2 probe), i.e. it scales the printed pattern, it does not un-weight it. | planner, 2026-10-01 (measured) |
| D4 | A mode maps to the displayed atoms directly when `len(vector)//3 == len(atoms)`, and through `sorted(qm_atom_indices)` when it equals `len(qm_atom_indices)`. When the pane's geometry matches neither, the mode is drawn on the job's own structure instead — `job.file_geometry` if set, else the coordinate file the input names, read with `core.geometry.read_xyz` — with the same `qm_atom_indices` split; otherwise it is refused (`None`), never guessed. | user (chose the fallback), 2026-10-01 (verified on two real jobs; MO2 found the pane's geometry is not the Hessian's) |
| D5 | Animation is `sin(2π (t − t0) / MODE_PERIOD_S)` with `MODE_PERIOD_S = 2.5`, at the pane's own `ROCK_FPS`; it runs in the text/braille pane as well, because that renderer is cheap and a mode that "does nothing" in text would repeat the see-through complaint. `0` (reset view) leaves the mode running. | planner, 2026-10-01 |
| D6 | `i` cycles: off → most-negative imaginary mode → next → … → off. The pane title carries `mode <index> <cm-1> cm⁻¹ (<j>/<n>)`. With no imaginary mode, `i` shows a Textual notification and changes nothing. | user (key) and planner (cycle), 2026-10-01 |
| D7 | `snapshot --mode N [--phase DEG]` renders one frame; `--phase` is degrees, default 90 (the `+` extreme). `N` is ORCA's own mode index, the one `orcamon freqs` prints. | planner, 2026-10-01 |
| D8 | The three fast verification commands per task; `run_gates.py` only in MO5. | planner, 2026-10-01 |
| D9 | Commit messages start with `orcamon:`; one commit per task. | planner, 2026-10-01 |
| D10 | Because the fallback can draw `job.xyz` instead of the parsed block, the pane's structure may change while a mode is on and the chart's scrub selection is ignored until the mode is turned off. | planner, 2026-10-01 |

## Facts (verified 2026-10-01)

**Environment.** Repository root `/home/madsr2d2/masterThesis`. Interpreter `.venv/bin/python` (never bare `python`); the CLI is `.venv/bin/orcamon`. orcamon is installed editable from `tools/orcamon/src/orcamon/`, so edits take effect without reinstalling.

**ORCA's `NORMAL MODES` block.** Every job that diagonalized a Hessian prints it into the `.out`, immediately after its `VIBRATIONAL FREQUENCIES` block. Format: an all-integer header line giving the mode indices of the block (6 per block, fewer in the last), then one row per Cartesian coordinate (`0 .. 3N-1`) with that many float values, then the next header, until the block ends at a blank line. Real excerpt (`verify_exact_hessian/job.out`, N=140, 420 coordinates, the first two blocks):
```
------------
NORMAL MODES
------------

These modes are the Cartesian displacements weighted by the diagonal matrix
M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
Thus, these vectors are normalized but *not* orthogonal

                  0          1          2          3          4          5    
      0       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
...
                   6          7          8          9         10         11    
      0      -0.004219   0.010924   0.018078   0.029563   0.000175  -0.007542
```
The same file's `VIBRATIONAL FREQUENCIES` block (`job.out` line 2851) has mode `6: -576.33 cm**-1 ***imaginary mode***`, and modes 0-5 are the zero translations/rotations. `optts_freq_tight/job.out` line 2183 has `6: -362.52 cm**-1 ***imaginary mode***` and 51 modes (N=17).

**The convention, measured.** `verify_exact_hessian/job.hess.v006.xyz` is `orca_pltvib`'s 20-frame trajectory of mode 6: 2840 lines = 20 frames of (2 header + 140 atom) lines, each atom line `element x y z dx dy dz`. Frame 1's `(dx,dy,dz)` equals **2.000 ×** the printed mode-6 vector for every atom checked, **for both C and H** (C coordinate 0: printed `-0.004219`, frame `-0.008437`; H atom 66 coordinate 198: printed `0.002789`, frame `0.005578`). So the printed vector is the Cartesian pattern and `orca_pltvib` scales it globally; there is no `sqrt(m)` factor to apply. The equilibrium in that file equals `job.xyz` (140 atoms) to 6 decimals.

**The mapping, measured.** The Hessian covers the ACTIVE region, which may differ from the displayed structure:
- `geometry_r2scan3c-xtb/ts/verify_exact_hessian`: `job.out` line 325 `optimized atoms = activeRegion ... 140`, line 328 `QM1 Subsystem ... 61 64 120 … 139` (20 atoms), Hessian N = 140. The `.out` prints exactly ONE `CARTESIAN COORDINATES (ANGSTROEM)` block (line 359) and it is 22 atoms, so `state.atoms` is 22 (measured); the 140-atom structure the Hessian belongs to is in `job.xyz`, which the input's `* xyzfile` names. The mode therefore maps through the fallback, and `qm_atom_indices` (20) still gives the host/guest split on the 140-atom render.
- `geometry_b973c-xtb/ts/optts_freq_tight`: `job.out` line 340 `optimized atoms = activeRegion ... 17`, line 343 `QM1 Subsystem ... 61 64 120 … 136` (17 atoms), Hessian N = 17 = `len(state.qm_atom_indices)`, while `state.atoms` and `job.xyz` are both 137. `job.activeRegion.xyz` is the 17-atom geometry the Hessian belongs to: `sorted(qm_atom_indices)` applied to `job.xyz` reproduces it exactly, and applied to `state.atoms` it is 1.79e-4 A off (the pane is one optimiser step behind).

**Parser architecture** (`tools/orcamon/src/orcamon/core/parser.py`):
- `feed_line` is the reference path; `feed_text` is the fast one and bounces whole lines past `_LINE_OF_INTEREST_RE` unless `_in_block()`. So a new block needs its start marker added to `_LINE_OF_INTEREST_PARTS` (line 126) and its open flag added to `_in_block()` (line 354).
- The block state machines are `_feed_blocks` (line 495), called from `_process` line 390 for every line. `_FREQ_HEADER_RE = re.compile(r"^VIBRATIONAL FREQUENCIES\s*$")`, `_FREQ_LINE_RE` matches `^\s*\d+:\s+(-?[\d.]+)\s+cm\*\*-1(...)`, `_GEOM_HEADER_RE = re.compile(r"^CARTESIAN COORDINATES \(ANGSTROEM\)\s*$")`.
- The frequency block closes on the first non-empty line that is not a table row after at least one row; it then sets `imaginary_freqs`, `frequencies`, `freqs_final` and `freq_cycle`. `frequencies[i]` is the value of ORCA mode `i`.
- `JobState` fields include `frequencies: list | None`, `imaginary_freqs`, `freqs_final`, `qm_atom_indices: set | None`, `points: deque` of `GeometryPoint`, `atoms: list`; atom entries are `(element, x, y, z)` tuples. `_IDENTITY_FIELDS` is `{"path", "stem", "has_out", "offset", "inode", "mtime"}`.
- `cache.py` keys an entry on the sha1 of `parser.py`'s source, so adding a field invalidates old entries by itself; it refuses to save while `_in_block()` is true.
- `validate.py` `MARKER_CASES` (line 63) drives each marker line through both `feed_line` and `feed_text` and asserts a field; `("----- Orbital basis set information -----", lambda s: s._in_orbital_basis)` is the shape to copy.

**`tools/orcamon/src/orcamon/tui/app.py`**:
- `GEOMETRY_BINDINGS` (line 302) ends `… ("f", "toggle_fog", "Fog"), ("x", "toggle_see_through", "See-through"), ("v", "next_representation", "View"), ("h", "toggle_hydrogens", "Hydrogens"), ("p", "next_axis", "Principal"), ("o", "toggle_rock", "Rock"), ("0", "reset_view", "Reset"), …`. `i` is free and `i` is not bound anywhere in the app.
- `RotatableGeometryImage.__init__` (line 357) sets `self.elev/azim/show_distances/show_labels/representation/fog/see_through/zoom/pan` and `ROCK_FPS = 5`. `view()` (line 390) builds the `View`. `action_toggle_rock` / `_rock_tick` / `_stop_rock` / `_rock_redraw` (lines 525-548) are the timer pattern to copy; `_rock_redraw` = `self._push(self._next_seq(), quality="full")`, and `TextGeometry` overrides it as `self.refresh()` with `ROCK_FPS = 10`.
- Kitty `_push` (line 697) and Herdr `_push` (line ~900) both do `atoms = self._atoms_for(job, point)` and then `geometry_render.render(atoms, self.view(), qm_atom_indices=job.state.qm_atom_indices, ...)`. `TextGeometry.render` (line 775) uses the cached `self._atoms`.
- `MonitorApp._geometry_title()` (line 1124) reads `geometry = self.query_one("#geometry")` and appends `" · fog"` and `" · see-through"`.
- `on_unmount` stops rock; `action_reset_view` stops rock but must not touch the mode.

**`tools/orcamon/src/orcamon/cli/__init__.py`**: the `snapshot` parser (line 126) has `-o/--output`, `--representation`, `--size`, `--elev`, `--azim`, `--no-labels`, `--no-fog`, `--see-through`, `--distances`.
**`cli/commands.py`**: `cmd_snapshot` (line 471) builds `View(...)`, calls `geometry_render.render(point.atoms, view, qm_atom_indices=state.qm_atom_indices, size_px=size)` and prints `"{output}: {w}x{h}, {n} atoms, {source}"`. `cmd_freqs` (line 513) already indexes modes as `enumerate(freqs)`.

**Skill.** `test_orcamon.py` runs `orcamon skill --check --project <repo>`; adding a CLI flag without `.venv/bin/orcamon skill install --project .` makes it fail. `COLUMNS=90 .venv/bin/orcamon snapshot --help` reproduces the fenced block under `### \`orcamon snapshot\`` in `tools/orcamon/README.md`.

**Tests.** `test_monitor.py` 136 `pass`, `test_orcamon.py` 157 `pass`, both `0 failure(s)`; `data/test_curve_metrics.py` `0 failure(s)`; `run_gates.py` 38 gates, 0 failed. Both root files define `check(label, ok, detail="")` printing `  pass  <label>` or `  FAIL  <label>`. Test names must be unique across all root `.py` files and `data/*.py`. `test_orcamon.py` imports `SYNTHETIC_OUTPUT` from `test_monitor` (line 181) and builds a synthetic job tree; its skill test runs every `orcamon …` line in a fenced block of the skill body.

**Keep the gate count at 38.** `test_root_documents.py` line 435 claims `run_gates.gate_paths()`'s lengths against `README.md` and `CLAUDE.md`; `CLAUDE.md` is the user's. No new root gate file.

## Conventions

- Match the surrounding code: module-level constants each with a one-line comment; private helpers start with `_`; docstrings and comments say WHY, naming the failure that motivated the code.
- numpy arrays in `tui/raster.py` stay float64; `core/` and `cli/` stay standard-library only.
- New tests go in the file the task names, use `check(...)` for every assertion, and are added to that file's `__main__` list right after the test named in the task.
- Expected values in tests are the hand-derived ones written in this plan, never values copied from a run.
- No network access in any task.

---

# Phase 1 — the mode data

### MO1 — Parse ORCA's `NORMAL MODES` block into the job state

**Why:** The displacement vectors are the one piece of a frequency job orcamon does not read; every mode cue needs them. The block is already in the `.out`, so no external tool is required.
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/parser.py`, `tools/orcamon/src/orcamon/validate.py`, `test_monitor.py`.

**Do:**
1. Add a `NormalMode` dataclass near `JobState` (same file): fields `index: int` (ORCA's mode number), `cm1: float` (its signed frequency) and `vector: list` (the `3N` printed values, coordinate-major). Add a docstring saying one normal mode of one Hessian.
2. Add `modes: list | None = None` to `JobState` (a list of `NormalMode`, the imaginary modes of the LAST `NORMAL MODES` block), with a comment mirroring `imaginary_freqs`'s ("... of the LAST frequency block"). Add the private block-state fields `_in_modes_block: bool = False`, `_modes_columns: list = field(default_factory=list)`, `_modes_values: dict = field(default_factory=dict)` (mode index -> `{coordinate: value}`), `_modes_seen_row: bool = False`.
3. Add regexes beside the frequency ones:
   ```python
   _NORMAL_MODES_HEADER_RE = re.compile(r"^NORMAL MODES\s*$")
   _MODE_COLUMNS_RE = re.compile(r"^\s*(\d+(?:\s+\d+)*)\s*$")
   _MODE_ROW_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+(?:\s+-?\d+\.\d+)*)\s*$")
   ```
   A column header is all integers; a data row has a coordinate index then floats with decimal points, so the two never collide.
4. In `_feed_blocks`, after the geometry branch, add:
   ```python
   if _NORMAL_MODES_HEADER_RE.match(stripped):
       self._in_modes_block = True
       self._modes_columns = []
       self._modes_values = {}
       self._modes_seen_row = False
   elif self._in_modes_block:
       m = _MODE_COLUMNS_RE.match(stripped)
       if m:
           self._modes_columns = [int(tok) for tok in m.group(1).split()]
           self._modes_seen_row = False
       else:
           m = _MODE_ROW_RE.match(line)
           if m:
               coord = int(m.group(1))
               values = [float(tok) for tok in m.group(2).split()]
               imaginary = {i for i, f in enumerate(self.frequencies or []) if f < 0}
               for column, mode_index in enumerate(self._modes_columns):
                   if mode_index in imaginary:
                       self._modes_values.setdefault(mode_index, {})[coord] = values[column]
               self._modes_seen_row = True
           elif not self._modes_seen_row:
               pass                     # the block's own preamble
           else:
               self._finish_modes()
   ```
5. Add `_finish_modes(self) -> None`: set `_in_modes_block = False`; let `coords = len(self.frequencies) // 3` if `self.frequencies` else `0`; for each index in `sorted(self._modes_values)`, keep it only when `coords` and every `c` in `range(3 * coords)` is present, building `NormalMode(index=index, cm1=self.frequencies[index], vector=[row[c] for c in range(3 * coords)])`; assign `self.modes = modes or None` (only if `modes`); clear the three private fields. Its docstring says the block's rows are incomplete until the blank line closes it, so nothing partial is stored.
6. On `_FREQ_HEADER_RE.match(stripped)` (the existing branch), also set `self.modes = None` — a later Hessian replaces the earlier one's modes, and a job whose second Hessian has none must not keep the first's.
7. Add `"NORMAL MODES"` to `_LINE_OF_INTEREST_PARTS` and `self._in_modes_block` to `_in_block()`.
8. In `validate.py`, add to `MARKER_CASES`: `("NORMAL MODES", lambda s: s._in_modes_block),`.
9. Add three tests to `test_monitor.py`, each added to the `__main__` list right after the last parser test there:
   - `test_normal_modes_are_parsed()` — feed a two-atom synthetic output (6 modes, one imaginary at index 0, one real at index 1, the rest any sign) followed by a `NORMAL MODES` block whose mode-0 column is `0.100000 0.200000 -0.300000 0.400000 -0.500000 0.600000` and whose other columns are zero. Exactly 3 checks: `len(state.modes) == 1`; `state.modes[0].index == 0 and state.modes[0].cm1 == -100.0`; `state.modes[0].vector == [0.1, 0.2, -0.3, 0.4, -0.5, 0.6]`.
   - `test_normal_modes_keep_only_imaginary()` — the same output with two imaginary modes at indices 0 and 1; exactly 1 check: `[m.index for m in state.modes] == [0, 1]`.
   - `test_normal_modes_take_the_last_block()` — two frequency+modes pairs, the second's imaginary mode at index 2 with value `-200.0`; exactly 1 check: `state.modes[0].index == 2 and state.modes[0].cm1 == -200.0`.
   Feed with `new_state(...)`/`feed_line` or `feed_text` as the surrounding tests do; give the exact snippet inline in each test.

**If unsure:** If `_MODE_COLUMNS_RE` also matches something the real output prints inside the block that is not a column header, widen it to require the line's whole token list to be integers and add the real line under **If unsure** in the report.

**Acceptance:**
- The three new tests pass, with 5 new `check` calls in total.
- `test_monitor.py`: 141 `pass`, `0 failure(s)`. `test_orcamon.py`: 157 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- The diff touches `parser.py`, `validate.py` and `test_monitor.py` only.

### MO2 — Map a mode onto the displayed geometry and displace it

**Why:** The printed vector is indexed against the Hessian's atoms, which may be the whole structure or just the QM layer; the pane's geometry can be either. The mapping, the scale and the phase belong in one tested place.
**Depends on:** MO1.
**Files:** `tools/orcamon/src/orcamon/core/vibrations.py` (new), `test_monitor.py`.

**Do:**
1. Create `core/vibrations.py`, standard library only, with a module docstring naming the `orca_pltvib` convention (it writes `2.0 *` the printed vector for C and H alike, so the printed vector is the Cartesian pattern and no mass factor is applied) and why the mapping is needed. Import `math` at module level and `from .geometry import read_xyz`. Contents:
   ```python
   MODE_AMPLITUDE_ANGSTROM = 0.20   # how far the largest-moving atom travels
   MODE_PERIOD_S = 2.5              # one full oscillation, seconds

   def imaginary_modes(state) -> list:
       """The imaginary modes, most negative first; [] until a FINAL block."""
       if state.modes is None or not state.freqs_final:
           return []
       return sorted(state.modes, key=lambda mode: mode.cm1)

   def target_indices(n_atoms: int, n_mode_atoms: int, qm_atom_indices) -> list | None:
       """Which displayed atom each mode row belongs to, or None."""
       if n_mode_atoms == n_atoms:
           return list(range(n_atoms))
       if qm_atom_indices and len(qm_atom_indices) == n_mode_atoms:
           return sorted(qm_atom_indices)
       return None

   def alternate_geometries(job) -> list:
       """Structures to try when the pane's own does not hold the mode.

       A multilayer `.out` can print only the QM block while the Hessian covers
       the whole model, so the mode's own structure is read from the file the
       input names. `job.file_geometry` is only ever set when the log printed
       nothing, so it is tried first but is usually None here."""
       candidates = []
       if getattr(job, "file_geometry", None) is not None:
           candidates.append(job.file_geometry.atoms)
       coords_file = getattr(getattr(job, "input", None), "coords_file", None)
       if coords_file:
           atoms = read_xyz(job.state.path / coords_file)
           if atoms:
               candidates.append(atoms)
       return candidates

   def mode_geometry(candidates: list, qm_atom_indices, mode) -> tuple | None:
       """The `(atoms, qm_atom_indices)` to draw `mode` on, from the candidates
       in order, or None when none of them holds the mode."""
       n_mode_atoms = len(mode.vector) // 3
       for atoms in candidates:
           if atoms and target_indices(len(atoms), n_mode_atoms, qm_atom_indices) is not None:
               return list(atoms), qm_atom_indices
       return None

   def offsets(atoms: list, qm_atom_indices, mode, amplitude=MODE_AMPLITUDE_ANGSTROM) -> list | None:
       """Per-atom (dx, dy, dz): `mode` scaled so its largest atom moves `amplitude`."""
       n_mode_atoms = len(mode.vector) // 3
       indices = target_indices(len(atoms), n_mode_atoms, qm_atom_indices)
       if indices is None:
           return None
       raw = [(mode.vector[3 * j], mode.vector[3 * j + 1], mode.vector[3 * j + 2])
              for j in range(n_mode_atoms)]
       biggest = max((math.sqrt(x * x + y * y + z * z) for x, y, z in raw), default=0.0)
       if biggest == 0.0:
           return None
       scale = amplitude / biggest
       result = [(0.0, 0.0, 0.0)] * len(atoms)
       for j, atom_index in enumerate(indices):
           x, y, z = raw[j]
           result[atom_index] = (x * scale, y * scale, z * scale)
       return result

   def displaced(atoms: list, mode_offsets: list, sine: float) -> list:
       """`atoms` moved `sine` of the way along the offsets."""
       return [(el, x + sine * dx, y + sine * dy, z + sine * dz)
               for (el, x, y, z), (dx, dy, dz) in zip(atoms, mode_offsets)]

   def phase_sine(elapsed_s: float, period_s: float = MODE_PERIOD_S) -> float:
       return math.sin(2.0 * math.pi * elapsed_s / period_s)
   ```
2. Add five tests to `test_monitor.py` (after MO1's tests), each added to the `__main__` list:
   - `test_mode_offsets_fill_the_whole_structure()` — a `NormalMode(index=0, cm1=-100.0, vector=[0.0,0.0,0.0, 3.0,4.0,0.0])` (two atoms, the second's vector length 5) mapped onto two atoms with `qm_atom_indices=None`. Exactly 3 checks: `offsets(...)` is not None; the second atom's offset has length `MODE_AMPLITUDE_ANGSTROM`; the first atom's offset is `(0.0, 0.0, 0.0)`.
   - `test_mode_offsets_map_onto_the_qm_subset()` — a 4-atom structure, `qm_atom_indices={1, 3}`, a `NormalMode` whose `vector` has 6 entries (2 mode atoms) with a nonzero first atom. Exactly 3 checks: the offsets are not None; atoms 1 and 3 carry the vector (atom 1's length is `MODE_AMPLITUDE_ANGSTROM`); atoms 0 and 2 are exactly `(0.0, 0.0, 0.0)`.
   - `test_mode_offsets_refuse_a_mismatched_mode()` — 3 displayed atoms, `qm_atom_indices=None`, a mode with 2 mode atoms. Exactly 1 check: `offsets(...) is None`.
   - `test_the_sine_phase_displaces_atoms()` — two atoms and known offsets. Exactly 2 checks: `phase_sine(0.0) == 0.0` and `displaced(atoms, offs, 1.0)` equals `atoms + offs`.
   - `test_the_mode_geometry_falls_back_to_an_alternate()` — a 2-atom pane, a 4-atom alternate, a 4-atom `NormalMode` and a 2-atom `NormalMode`. Exactly 3 checks: `mode_geometry([pane, alternate], None, mode4) == (alternate, None)`; `mode_geometry([pane], None, mode4) is None`; `mode_geometry([pane, alternate], None, mode2) == (pane, None)`.
3. Write the ground-truth probe below to `/tmp/orcamon_mode_parity.py` (outside the repository) and run it from the repository root with `XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python /tmp/orcamon_mode_parity.py`. It must print four lines: `verify_exact_hessian: N=<n>, imaginary modes <list>`; `parity: max |d - 2*printed| <dev> over <n> components`; `optts_freq_tight: N=<n>, qm=<n>, atoms=<n>`; `qm map: coordinates match <True|False>`.
   ```python
   """Ground-truth check of the parsed normal modes against orca_pltvib's own xyz."""
   from pathlib import Path
   from orcamon.cli import build_parser, commands
   from orcamon.core import vibrations

   REPO = Path("/home/madsr2d2/masterThesis")
   C = "computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/"

   def load(job):
       args = build_parser().parse_args(["snapshot", C + job, "-o", "/dev/null", "--root", str(REPO)])
       return commands._load(args, C + job)

   # 1. The 140-atom Hessian: the `.out` prints only a 22-atom QM block, so the
   #    mode is drawn on job.xyz; orca_pltvib's own frame 1 writes 2.0 * printed.
   job = load("geometry_r2scan3c-xtb/ts/verify_exact_hessian")
   state = job.state
   modes = vibrations.imaginary_modes(state)
   mode = modes[0]
   chosen = vibrations.mode_geometry([state.atoms] + vibrations.alternate_geometries(job),
                                     state.qm_atom_indices, mode)
   atoms, qm = chosen
   frame = (REPO / C / "geometry_r2scan3c-xtb/ts/verify_exact_hessian/job.hess.v006.xyz").read_text().splitlines()
   diffs = []
   for i, line in enumerate(frame[2:2 + len(atoms)]):
       dx, dy, dz = (float(t) for t in line.split()[-3:])
       printed = mode.vector[3 * i:3 * i + 3]
       diffs.extend(d - 2.0 * p for d, p in zip((dx, dy, dz), printed))
   print(f"verify_exact_hessian: N={len(mode.vector) // 3}, imaginary modes {[m.index for m in modes]}")
   print(f"parity: max |d - 2*printed| {max(abs(x) for x in diffs):.3g} over {len(diffs)} components")

   # 2. The 17-atom Hessian: sorted(qm) applied to the FULL structure must be the
   #    activeRegion xyz ORCA wrote.
   job = load("geometry_b973c-xtb/ts/optts_freq_tight")
   state = job.state
   mode = vibrations.imaginary_modes(state)[0]
   full = next((a for a in vibrations.alternate_geometries(job) if len(a) == 137), state.atoms)
   qm = sorted(state.qm_atom_indices)
   picked = [full[i] for i in qm]
   ref = (REPO / C / "geometry_b973c-xtb/ts/optts_freq_tight/job.activeRegion.xyz").read_text().splitlines()
   same = len(ref) - 2 == len(picked) and all(
       line.split()[0] == el and max(abs(float(a) - b) for a, b in zip(line.split()[1:4], (x, y, z))) < 1e-4
       for line, (el, x, y, z) in zip(ref[2:], picked))
   print(f"optts_freq_tight: N={len(mode.vector) // 3}, qm={len(qm)}, atoms={len(full)}")
   print(f"qm map: coordinates match {same}")
   ```
   Adjust only the `commands._load`/`build_parser` call names if the real API differs (read `commands.py` first); the four printed quantities are the check.

**If unsure:** If `commands._load(args, path)` is not the loader the CLI uses, read `commands.py` and use the real one — that is not a deviation, but say which you used. The check is the convention itself, `d == 2 * printed`; both files are written to six decimals, so the only tolerance is `|d - 2*printed| <= 1.5e-6`. If `max |d - 2*printed|` exceeds `2e-6`, STOP and report NEEDS_USER: the convention in D3 is then wrong, and it is the user's call.

**Acceptance:**
- The five new tests pass, with 12 new `check` calls in total.
- The probe prints `verify_exact_hessian: N=140, imaginary modes [6]`, a `parity: max |d - 2*printed|` at or below `2e-6` (six-decimal printing is the only tolerance), `optts_freq_tight: N=17, qm=17, atoms=137`, and `qm map: coordinates match True`.
- `test_monitor.py`: 153 `pass`, `0 failure(s)`. `test_orcamon.py`: 157 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- The diff touches `core/vibrations.py` and `test_monitor.py` only.

# Phase 2 — the surfaces

### MO3 — `snapshot --mode N --phase DEG`

**Why:** A still of one extreme of one mode is how the result is shared and checked without a terminal; it is also the only way to test the whole path headlessly.
**Depends on:** MO2.
**Files:** `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/README.md`, `.claude/skills/orcamon/SKILL.md`, `test_orcamon.py`.

**Do:**
1. In the `snapshot` parser, after `--see-through`, add:
   ```python
   p.add_argument("--mode", type=int, metavar="N",
                  help="animate an imaginary normal mode: ORCA's mode number (the one `orcamon freqs` prints)")
   p.add_argument("--phase", type=float, default=90.0, metavar="DEG",
                  help="where in the oscillation to draw, 0-360 degrees (default: %(default)s)")
   ```
2. In `cmd_snapshot`, after `point` is chosen and before `geometry_render.render`, when `args.mode is not None`: import `from ..core import vibrations`; take `modes = vibrations.imaginary_modes(state)`; if `args.mode` is not in `[m.index for m in modes]`, raise `UsageError(f"no imaginary mode {args.mode} in this job; available: {[m.index for m in modes] or 'none'}")`; else choose the drawn structure first: `candidates = [point.atoms] + vibrations.alternate_geometries(job)`, `chosen = vibrations.mode_geometry(candidates, state.qm_atom_indices, mode)`; if `chosen is None` raise `UsageError(f"no geometry in this job holds mode {args.mode} ({len(mode.vector) // 3} atoms)")`; then `atoms, qm = chosen`, `offs = vibrations.offsets(atoms, qm, mode)`, and if `offs is None` raise `UsageError(...)` as above, else `atoms = vibrations.displaced(atoms, offs, math.sin(math.radians(args.phase)))`. Import `math` at module level if it is not already imported. Pass `atoms` and `qm_atom_indices=qm` to `geometry_render.render` instead of `point.atoms` and `state.qm_atom_indices`, and append `f", mode {args.mode} at {args.phase:g} deg"` to the printed line.
3. Run `.venv/bin/orcamon skill install --project .` from the repository root. Then regenerate the fenced block under `### \`orcamon snapshot\`` in `tools/orcamon/README.md` with the exact output of `COLUMNS=90 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/orcamon snapshot --help`. Change no other help block.
4. Add to `test_orcamon.py`, in the snapshot test (the one that already writes a PNG), 3 checks using a synthetic job whose output carries a frequency block and a `NORMAL MODES` block (add such a job to the synthetic tree; reuse `test_monitor.SYNTHETIC_OUTPUT` plus the mode block, or a small helper in `test_orcamon.py`):
   - `snapshot <job> -o <tmp>/mode.png --mode <imaginary index> --size 120x100` exits 0 and the file starts with `b"\x89PNG"`;
   - `snapshot <job> --mode 999 …` exits 2 and the message names the available mode;
   - `snapshot <job-without-frequencies> --mode 0 …` exits 2.

**If unsure:** If the synthetic tree cannot host a frequency job without disturbing the existing tests, add the smallest standalone synthetic job directory the tree builder allows, and say so in the report.

**Acceptance:**
- `test_orcamon.py`: 160 `pass`, `0 failure(s)`; the skill check passes against the reinstalled skill.
- `COLUMNS=90 … orcamon snapshot --help` equals the README's fenced snapshot block byte for byte, and the README diff touches no other block.
- `test_monitor.py`: 153 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.

### MO4 — The TUI `i` key: cycle the imaginary modes and animate

**Why:** The mode is for a person watching a job; the key and the oscillation are the feature.
**Depends on:** MO3.
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `test_orcamon.py`.

**Do:**
1. In `GEOMETRY_BINDINGS`, add `("i", "cycle_modes", "Modes")` immediately after the `"x"` entry.
2. In `RotatableGeometryImage.__init__`, add `self.mode_index = None`, `self._mode_alternates = []`, `self.mode_label = ""`, `self._mode_s = 0.0`, `self._mode_t0 = 0.0`, `self._mode_timer = None`.
3. Add the methods, mirroring rock:
   - `action_cycle_modes`: `modes = vibrations.imaginary_modes(self._job.state) if self._job is not None else []`; if empty, `self.notify("no imaginary mode in this job")` and return. Otherwise let `indices = [m.index for m in modes]`; if `self.mode_index is None` take `indices[0]`, else take the next one or `None` past the end. When it becomes `None` call `_stop_modes()` and set `self._mode_alternates = []`; otherwise call `_start_modes()`, set `self._mode_alternates = vibrations.alternate_geometries(self._job)`, and set `self.mode_label` to `f"mode {index} {cm1:.1f} cm-1 ({k + 1}/{len(indices)})"`. Finish with `self._input.request()` and `self.app.update_detail()`.
   - `_start_modes` resets `_mode_t0 = time.monotonic()` and starts `self.set_interval(1.0 / self.ROCK_FPS, self._mode_tick)`; `_mode_tick` sets `self._mode_s = vibrations.phase_sine(time.monotonic() - self._mode_t0)` and calls `_rock_redraw()`; `_stop_modes` stops the timer, sets `_mode_timer = None`, `mode_index = None`, `mode_label = ""` and `_mode_s = 0.0`. Add `vibrations` to the module imports (`from ..core import vibrations`) — `app.py` is TUI code and may import it.
   - `_apply_mode(self, atoms, job)`: return `atoms` when `self.mode_index is None` or `job is None`; find the `NormalMode` with that index in `vibrations.imaginary_modes(job.state)`; return `atoms` if it is None; `chosen = vibrations.mode_geometry([atoms] + self._mode_alternates, job.state.qm_atom_indices, mode)`; return `atoms` if `chosen is None`; `drawn, qm = chosen`; `offs = vibrations.offsets(drawn, qm, mode)`; return `atoms` if `offs is None`; else `vibrations.displaced(drawn, offs, self._mode_s)`. The renderer keeps `job.state.qm_atom_indices` for the host/guest split, which is the same global index set on either structure.
4. Call `_apply_mode` on the atoms right before each render: in Kitty `_push` and Herdr `_push`, wrap the `self._atoms_for(job, point)` result; in `TextGeometry.render`, wrap `self._atoms`. In `show_job`, when the job actually changed (`job is not self._job`), call `_stop_modes()` first. `on_unmount` calls `_stop_modes()` as well as `_stop_rock()`. `action_reset_view` must NOT touch the mode.
5. In `MonitorApp._geometry_title()`, after the see-through mark, add `mode = f" · {geometry.mode_label}" if geometry.mode_label else ""` and include it in the returned string.
6. In `test_orcamon.py`'s headless TUI test (`test_the_tui_runs_headless`), after the hydrogen check, add 3 checks on the synthetic frequency job: press `i`, `await pilot.pause(0.2)`, and check `geometry.mode_index == <imaginary index> and "cm-1" in geometry.border_title`; press `i` again and check `geometry.mode_index is None and "cm-1" not in geometry.border_title`; select a job with no frequencies, press `i`, and check `geometry.mode_index is None`.

**If unsure:** If `self.notify(...)` is not available on the widget or disrupts headless tests, use `self.app.notify(...)`; if neither is safe, drop the notification and record it under NOTES — the mode staying off is the asserted behaviour.

**Acceptance:**
- The 3 new checks pass and `test_orcamon.py` is 163 `pass`, `0 failure(s)`.
- The TUI test proves the title carries the mode when on and not when off.
- `test_monitor.py`: 153 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- The diff touches `tui/app.py` and `test_orcamon.py` only.

# Phase 3 — documentation, gates, pictures

### MO5 — Document the modes, run every gate, render the comparison pictures

**Why:** The README says nothing about modes, and the user needs pictures to judge the result.
**Depends on:** MO4.
**Files:** `tools/orcamon/README.md`.

**Do:**
1. In `tools/orcamon/README.md`, in the pixel-pane paragraph, add that `i` animates the job's imaginary normal modes (the most negative first) and that the text pane animates them too.
2. In the key list paragraph starting "`[`/`]` zoom, ctrl+arrows pan", add "`i` animates the job's imaginary normal modes" after the see-through key.
3. In `## Saving a picture`, add `--mode` / `--phase` to the list of options.
4. Run `run_gates.py` (§ Verification commands; at least 900000 ms timeout).
5. Render pictures into `/tmp/orcamon_depth_cues/modes/`, outside the repository, from the repository root with `XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg`, for the 140-atom job `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/ts/verify_exact_hessian`:
   - `orcamon snapshot <job> -o …/mode6_plus.png --mode 6 --phase 90` and `--phase 270` (`mode6_minus.png`);
   - measure the two PNGs' differing pixel count and report it.
6. Report in the commit message the two picture paths, their differing-pixel count, the run_gates count and time, and that the TUI animation could not be checked headless (Kitty and herdr pixel panes are unavailable; the braille pane needs a terminal).

**Acceptance:**
- `run_gates.py`: `38 gates`, `0 failed`.
- `test_monitor.py`: 153 `pass`. `test_orcamon.py`: 163 `pass`.
- The two PNGs exist in `/tmp/orcamon_depth_cues/modes/` and differ; give their sizes.
- The README diff touches only the three places named in **Do**.

## Out of scope

1. Animating real (non-imaginary) modes, or a mode browser for them (D2).
2. Reading the `.hess` file or running `orca_pltvib`/`orca_vib` (D1).
3. Mass-un-weighting the printed vector; it would not match `orca_pltvib` (D3).
4. An amplitude or speed control, or a per-mode `--amplitude` flag.
5. Drawing the mode as arrows or a static `+`/`−` overlay; only the oscillation and `snapshot` frames are built.
6. A GIF export, or a mode animation in `snapshot` beyond one phase.
7. Persisting the chosen mode across job selections.
8. Any change to the braille renderer's own drawing code; it animates through the shared atom path or not at all.

## Progress template

```markdown
# orcamon — animate imaginary normal modes — Progress

Plan: `PLAN_ORCAMON_VIBRATIONS.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| MO1 Parse NORMAL MODES | TODO | | | |
| MO2 Map and displace a mode | TODO | | | |
| MO3 `snapshot --mode` | TODO | | | |
| MO4 TUI `i` animation | TODO | | | |
| MO5 Docs, gates, pictures | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 136 | 157 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
