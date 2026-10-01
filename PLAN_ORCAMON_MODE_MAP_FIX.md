# orcamon — draw a mode on a structure it fits (MF1–MF1)

**Written:** 2026-10-01. **Verified against:** `master` at `abcc25b` (the `PLAN_ORCAMON_VIBRATIONS.md` run is committed at `86f6fba`).
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_MODE_MAP_FIX_PROGRESS.md` (template at the end).

A user selected `K+H2O2_water-relay_to_KP/…/ts/optts` and `i` correctly reported no imaginary mode (that job is an `OptTS`, not a `Freq` job). Checking the sibling `ts/optts_numfreq` showed a real defect in the shipped feature: `core/vibrations.target_indices` maps a mode through `qm_atom_indices` whenever the **counts** match, without checking that those indices fit the structure being drawn. `optts_numfreq`'s `.out` prints only a 19-atom QM+link block, while its 17-atom Hessian's ORCA indices run to 136, so `offsets` indexes a 19-element list with 136 and raises `IndexError`. This plan adds the missing fit check and makes the "no modes" message name the real cause. **The goal a deviation must still meet:** a mode is drawn only on a structure whose atom count and index range it actually fits — the pane, else the job's own `.xyz`, else refused — and a job with no frequency block says so rather than "no imaginary mode".

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 153 `pass` lines, last line `0 failure(s)` | ~1 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python test_orcamon.py` | 163 `pass` lines, last line `0 failure(s)` | ~6 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~17 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (same for `test_orcamon.py`). This plan is one task and runs `run_gates.py` in it.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's; a new gate file would change the count and break that gate. New tests go in `test_monitor.py` or `test_orcamon.py`.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- the probe in MF1 step 6 does not raise before the fix or does not print the expected values after it (report the real output);
- a hand-derived value written in this plan does not hold after the change (report the measured value; do not change the expected value or the fix to make it pass);
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance;
- a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this task.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's.
- Never `git add -A`, `git add .` or `git commit -a`; stage the files this task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file.
- Never run `orcamon` against the user's cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg`. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never write into `computational/`, and never run `orca_pltvib` or any other ORCA tool. The real job appears only in the MF1 probe, outside the repository; tests use synthetic atoms only.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `core` and `cli` import only the standard library at module level.
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill.

## Scope

**In scope:** the fit check in `core/vibrations.py::target_indices`; the "no frequency block" wording in `tui/app.py::action_cycle_modes` and `cli/commands.py::cmd_snapshot`; tests; the probe.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | The QM-subset branch of `target_indices` requires `max(qm_atom_indices) < n_atoms` as well as the count match: the QM indices are indices into the FULL structure, so they only address a list that is the full system. A structure they overrun is skipped, and the file fallback is tried instead. | user (reported the defect), 2026-10-01 |
| D2 | The message distinguishes "this job has no frequency block" from "no imaginary mode": a job like `ts/optts` computes no frequencies at all and should say so. | user, 2026-10-01 |
| D3 | One task; `run_gates.py` runs in it. | planner, 2026-10-01 |
| D4 | The direct map (`n_mode_atoms == n_atoms`) and the `alternate_geometries` file fallback are unchanged. | planner, 2026-10-01 |

## Facts (verified 2026-10-01)

**The shipped code.** `tools/orcamon/src/orcamon/core/vibrations.py` at `f9c1e48`:
- `target_indices(n_atoms, n_mode_atoms, qm_atom_indices)` (lines 31-37) is the only place the mapping is decided: `n_mode_atoms == n_atoms` -> `range(n_atoms)`; `qm_atom_indices and len(qm_atom_indices) == n_mode_atoms` -> `sorted(qm_atom_indices)`; else `None`.
- `alternate_geometries(job)` (line 40) returns `job.file_geometry.atoms` when set, then the atoms of the coordinate file `job.input.coords_file` (read with `read_xyz`).
- `mode_geometry(candidates, qm_atom_indices, mode)` (line 58) returns the first candidate for which `target_indices` is not None.
- `offsets(atoms, qm_atom_indices, mode, ...)` (line 68) re-derives the indices and writes `result[atom_index]` for each; a `sorted(qm)` whose values exceed `len(result)` raises `IndexError`.
- `imaginary_modes(state)` (line 24) returns `[]` when `state.modes is None` OR `not state.freqs_final`.

**The TUI hook.** `tui/app.py`: `action_cycle_modes` at line 561; when `imaginary_modes` is empty it calls `self.notify("no imaginary mode in this job")` at line 569 and returns. `_apply_mode` at line 609 calls `mode_geometry([atoms] + self._mode_alternates, job.state.qm_atom_indices, mode)` and then `offsets`.

**The CLI hook.** `cli/commands.py::cmd_snapshot`: the `args.mode` block added by MO3 raises `UsageError(f"no imaginary mode {args.mode} in this job; available: {[m.index for m in modes] or 'none'}")` when the requested mode is not in `imaginary_modes(state)`.

**The failing job.** `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts_numfreq/job.out`:
- line 318 `Size of QMMM System ... 137`, line 320 `Size of QM1 Subsystem ... 17`, line 322 `Size of QM1 Subsystem plus link atoms ... 19`, line 324 `optimized atoms = activeRegion ... 17`, lines 327-330 `QM1 Subsystem ... 61 64 120 121 122 124 125 126 128 129 130 131 132 133 134 135 136`.
- exactly ONE `CARTESIAN COORDINATES (ANGSTROEM)` block, at line 359, and it is the 19-atom `QM - SMALL SYSTEM`, so `state.atoms` is 19.
- its one `VIBRATIONAL FREQUENCIES` block (line 1850) has 51 modes; line 1861 `6: -379.28 cm**-1 ***imaginary mode***`.
- `job.xyz` is 137 atoms (10309 bytes) and the input names it (`* xyzfile 0 1 job.xyz`), so `alternate_geometries` finds it.
- So `target_indices(19, 17, <17 indices to 136>)` takes the QM branch and `offsets` writes `result[136]` on a 19-element list.

**The two jobs that were used as ground truth and are unaffected:**
- `geometry_r2scan3c-xtb/ts/verify_exact_hessian`: 22-atom printed block, `n_mode_atoms = 140`, 20 QM indices -> the direct map fails and the QM count differs, so the file fallback (140) is used.
- `geometry_b973c-xtb/ts/optts_freq_tight`: `state.atoms` 137, `n_mode_atoms = 17`, 17 QM indices with max 136 < 137 -> the QM subset applies.

**Tests.** `test_monitor.py` 153 `pass`; `test_orcamon.py` 163 `pass`; `data/test_curve_metrics.py` `0 failure(s)`; `run_gates.py` 38 gates, 0 failed. Both root files define `check(label, ok, detail="")`. Existing tests this task must not break: `test_monitor.test_mode_offsets_map_onto_the_qm_subset` (a 4-atom structure with `qm_atom_indices={1, 3}`, indices fit) and `test_the_mode_geometry_falls_back_to_an_alternate` (`qm_atom_indices=None`).

## Conventions

- Match the surrounding code: one-line comments where the reason is not obvious; docstrings say WHY.
- New tests go in the file the task names, use `check(...)` for every assertion, and are added to that file's `__main__` list right after the test named in the task.
- Expected values are the hand-derived ones written in this plan, never values copied from a run.
- No network access.

---

# Phase 1 — fit the mapping, name the cause

### MF1 — Require the QM indices to fit, and say when there is no frequency block

**Why:** `optts_numfreq` prints a 19-atom QM+link block for a 17-atom Hessian whose indices run to 136; `target_indices` accepts the QM branch on the count alone, so `offsets` writes outside the drawn list. And `ts/optts` (no `Freq` step at all) is told "no imaginary mode", which names the wrong cause.
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/vibrations.py`, `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `test_monitor.py`, `test_orcamon.py`.

**Do:**
1. In `core/vibrations.py::target_indices`, require the QM indices to fit the list they address:
   ```python
   if (qm_atom_indices and len(qm_atom_indices) == n_mode_atoms
           and max(qm_atom_indices) < n_atoms):
       return sorted(qm_atom_indices)
   ```
   Add one comment line saying the QM indices are ORCA's indices into the FULL system, so they only address a list that is the full system. Change nothing else in the module.
2. In `tui/app.py::action_cycle_modes`, replace the single notification with the two causes:
   ```python
   if self._job.state.frequencies is None:
       self.notify("this job computes no frequencies (no Freq or NumFreq step)")
   else:
       self.notify("no imaginary mode in this job")
   ```
3. In `cli/commands.py::cmd_snapshot`, when `args.mode is not None` and `state.frequencies is None`, raise `UsageError("this job has no frequency block (no Freq or NumFreq step)")` before the existing "no imaginary mode {N} … available …" check.
4. Add `test_the_qm_map_needs_indices_that_fit()` to `test_monitor.py`, registered in `__main__` right after `test_the_mode_geometry_falls_back_to_an_alternate`. It has exactly 2 `check` calls:
   - a 19-atom `pane`, a 137-atom `alternate`, a `NormalMode(index=6, cm1=-379.28, vector=[0.01] * 51)` (17 mode atoms), and `qm = {61, 64, 120, 121, 122, 124, 125, 126, 128, 129, 130, 131, 132, 133, 134, 135, 136}`;
   - `target_indices(19, 17, qm) is None`;
   - `mode_geometry([pane, alternate], qm, mode) == (alternate, qm)`.
5. Add 1 check to the snapshot test in `test_orcamon.py`: `snapshot <the job with no frequency block> --mode 0` exits 2 and its output names the missing frequency block (contains `no frequency block`).
6. Write the probe below to `/tmp/orcamon_qm_map_probe.py` (outside the repository). Run it BEFORE applying step 1 and record the exception, then after all edits and record the values.
   ```python
   """The QM mapping on a job whose .out prints only the QM+link block."""
   from pathlib import Path
   from orcamon.cli import build_parser, commands
   from orcamon.core import vibrations

   REPO = Path("/home/madsr2d2/masterThesis")
   C = "computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/"
   job_path = "geometry_b973c-xtb/ts/optts_numfreq"
   args = build_parser().parse_args(["snapshot", C + job_path, "-o", "/dev/null", "--root", str(REPO)])
   job = commands._load(args, C + job_path)
   state = job.state
   mode = vibrations.imaginary_modes(state)[0]
   chosen = vibrations.mode_geometry([state.atoms] + vibrations.alternate_geometries(job),
                                     state.qm_atom_indices, mode)
   print(f"pane atoms {len(state.atoms)}, qm {len(state.qm_atom_indices)}, mode {mode.index} {mode.cm1:.2f}")
   atoms, qm = chosen
   offs = vibrations.offsets(atoms, qm, mode)
   biggest = max((dx * dx + dy * dy + dz * dz) ** 0.5 for dx, dy, dz in offs)
   print(f"drawn atoms {len(atoms)}, max |offset| {biggest:.4f}")
   ```
   Before the fix it must raise `IndexError`; after it must print `drawn atoms 137` and `max |offset| 0.2000` (the `MODE_AMPLITUDE_ANGSTROM` peak) for `mode 6`.
7. Run all four verification commands.

**If unsure:** If the synthetic "job with no frequency block" used by step 5 is awkward to reach in `test_orcamon.py` without disturbing an existing test's counts, reuse the one the MO3 test already uses for the no-frequencies case and say so in the report.

**Acceptance:**
- The probe raises `IndexError` before step 1 and prints `drawn atoms 137` and `max |offset| 0.2000` after; paste both.
- `test_the_qm_map_needs_indices_that_fit` passes both checks; `test_monitor.py`: 155 `pass`, `0 failure(s)`.
- `test_orcamon.py`: 164 `pass`, `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.
- `run_gates.py`: `38 gates`, `0 failed`.
- The diff touches the five named files only, and `git diff HEAD -- '*test*'` shows the two added tests and no change to an existing assertion or expected value.

## Out of scope

1. The rest of the normal-mode feature (parsing, mapping rules, animation, `snapshot --mode`) beyond the fit check.
2. Any change to the direct map or the `alternate_geometries` fallback.
3. Animating real modes, amplitude/speed controls, or arrows (the parent plan's out-of-scope list still stands).
4. Wording changes beyond the two messages named in D2.
5. Any new CLI flag or command.

## Progress template

```markdown
# orcamon — draw a mode on a structure it fits — Progress

Plan: `PLAN_ORCAMON_MODE_MAP_FIX.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| MF1 Fit the QM map, name the cause | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 153 | 163 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
