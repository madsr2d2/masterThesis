# orcamon — see-through ghosts the environment layer (SG1–SG1)

**Written:** 2026-10-01. **Verified against:** `master` at `bbab83d` (the parent plan's DC6 progress commit; the renderer is at `50826d7`).
**Branch:** `master` (the agents commit straight to `master`, as the parent plan did).
**Progress file:** `PLAN_ORCAMON_SEE_THROUGH_GHOST_PROGRESS.md` (template at the end).

The `x` / `snapshot --see-through` toggle added by `PLAN_ORCAMON_DEPTH_CUES.md` (DC5) is a no-op in practice: it blends only the pixels where the environment layer *strictly wins* over a guest surface, measured at **688 of 675,000 pixels (0.10%)** for the default view on `geometry_b973c-xtb/ts/optts` at 900x750, so toggling it looks like nothing happens. This plan makes see-through ghost the **whole environment layer**: every environment pixel is composited at its own opacity against whatever lies behind it, so the host fades to a faint skeleton and the guest reads through it. **The goal a deviation must still meet:** with see-through OFF nothing changes at all; with it ON an environment-only pixel is faded toward the background, an environment pixel covering a guest keeps the existing 0.3/0.7 guest blend, and no pixel the guest wins changes.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 136 `pass` lines, last line `0 failure(s)` | ~1 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python test_orcamon.py` | 157 `pass` lines, last line `0 failure(s)` | ~6 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~17 s |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (same for `test_orcamon.py`).

`run_gates.py` is **not** required for SG1: the change touches only `tools/orcamon/` and `test_monitor.py`, which no other gate reads except the duplicate guard, and the gate count (38) is unaffected because `run_gates.py` discovers gates, not tests. This is the parent plan's D10 (only its docs/gates task ran `run_gates.py`).

`test_monitor.py` runs its tests from an explicit list in `if __name__ == "__main__":` at the bottom. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- the hand-derived value written in this plan does not hold after the algorithm is implemented as written (report the measured value; do not change the expected value, the tolerance, or the algorithm to make it pass);
- `test_see_through_shows_the_guest_behind_the_host` fails after the change;
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance;
- the working tree has changes not made by this task.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's.
- Never `git add -A`, `git add .` or `git commit -a`; stage the files this task lists, by path.
- Never push.
- Never change the CLI parser, the `--see-through` help string, `tools/orcamon/README.md` or `.claude/skills/orcamon/SKILL.md`: the flag keeps its name and help text, so the installed skill stays current and no regeneration is needed.
- Never run `orcamon` against the user's cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg`. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `raster.py` already uses numpy; no new import is needed.
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill. Tests use synthetic atoms only.

## Scope

**In scope:** the see-through branch of `render` in `tui/raster.py`; the `SEE_THROUGH_ALPHA` comment; one new test in `test_monitor.py` and its `__main__` entry.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | See-through ghosts the **whole environment layer**: every environment-winning pixel is composited at `SEE_THROUGH_ALPHA` against the guest colour where a guest surface is behind it, and against the background elsewhere. This replaces the strict-win-only blend. | user (chose the ghost over removal or the projection blend), 2026-10-01 |
| D2 | `SEE_THROUGH_ALPHA` stays **0.3**. That is the value already wired; lowering it would require editing the DC5 test's expected value, which the Never list forbids, and 0.3 already reads as a faint skeleton against the dark background. | planner, 2026-10-01 |
| D3 | With see-through OFF the render is byte-identical to `bbab83d`; the change is confined to pixels the environment layer wins. | planner, 2026-10-01 |
| D4 | Docs, README and CLI help are unchanged: the one-line `--see-through` help still describes showing the guest through the environment, and changing it would force a skill reinstall and a README help-block regeneration for no behaviour change. | planner, 2026-10-01 |
| D5 | The three fast verification commands only; no `run_gates.py` (parent plan D10). | planner, 2026-10-01 |

## Facts (verified 2026-10-01)

**Environment.** Repository root `/home/madsr2d2/masterThesis`. Interpreter `.venv/bin/python` (never bare `python`); the CLI is `.venv/bin/orcamon`. orcamon is installed editable from `tools/orcamon/src/orcamon/`, so edits take effect without reinstalling.

**`tools/orcamon/src/orcamon/tui/raster.py` at `50826d7`** (the DC5 commit):
- Constants include `BACKGROUND = (30, 30, 30)`; `FOG = 0.6`; `HOST_FOG = 0.85`; `HOST_FOG_POWER = 2.0`; `OUTLINE_JUMP = 0.35`; `HALO_MAX_PX = 5`; `HALO_FULL_GAP = 2.5`; `BOUNDARY_BALL_SCALE = 0.6`; `SEE_THROUGH_ALPHA = 0.3` (line 61, comment "how much of the host shows where it covers the QM region in see-through").
- `render(atoms, view=None, *, qm_atom_indices=None, size_px=(900, 750), bond_list=None) -> PIL.Image` order (lines 462-511): build `env`/`qm` `_Layer`s → `_draw_geometry(..., env, qm)` → fog each layer when `view.fog` → `frame, env_wins = _composite(env, qm)` → the see-through branch (lines 489-498) → `foreground = frame.zbuf != -np.inf` → background fill → `_apply_halos` → `_draw_labels` → `_draw_distances`.
- Current see-through branch, verbatim:
  ```python
      if view.see_through:
          # A host stick in front of a QM ball hides most of it from some
          # angles, so where the two overlap the pixel becomes a blend towards
          # the QM layer's colour. The owner changes with the colour, so a
          # covered QM atom keeps its label; the depth stays the HOST's, so the
          # halo still rings the host that is in front.
          blend = env_wins & (qm.zbuf != -np.inf)
          frame.color[blend] = (SEE_THROUGH_ALPHA * env.color[blend]
                                + (1.0 - SEE_THROUGH_ALPHA) * qm.color[blend])
          frame.owner[blend] = qm.owner[blend]
  ```
- `_composite(env: _Layer, qm: _Layer) -> tuple[_Layer, np.ndarray]` returns `env_wins = env.zbuf > qm.zbuf` (strict; the QM layer wins a tie) and `frame` with the winning depth, colour and owner. `_Layer.zbuf` is `-inf` where nothing was drawn.
- `project(atoms, view, size_px) -> (sx, sy, depth, scale)`; tests use `raster.render`, `raster.project`, `raster.BACKGROUND`.

**`tools/orcamon/src/orcamon/core/geometry.py`**: `View` is a plain (not frozen) dataclass, last field `see_through: bool = False`. So two `View(...)` instances can be constructed for the same geometry with the flag on and off.

**Measured on the working tree (read-only diagnostic, 2026-10-01)**, job `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts` (137 atoms, 17 QM), 900x750, `View()` vs `View(see_through=True)`: `env_wins & (qm.zbuf != -inf)` is 688 pixels (0.10%); the two renders differ in 602 pixels. At a space-filling view (`elev=-30, azim=200`) it is 3,052. This is why the toggle reads as inert.

**Tests.** `test_monitor.py` now has 133 `pass`; `test_orcamon.py` 157 `pass`; `data/test_curve_metrics.py` `0 failure(s)`; `run_gates.py` 38 gates, 0 failed (not run by this task). Both root test files define `check(label, ok, detail="")`, which prints `  pass  <label>` or `  FAIL  <label>`. Test names must be unique across all root `.py` files and `data/*.py`; `data/test_curve_metrics.py` fails otherwise. `test_see_through_shows_the_guest_behind_the_host` is at `test_monitor.py:769`; the `__main__` list calls it at line 1197. `test_host_fog_keeps_near_host_bright` (line 647) uses the atoms this task's test reuses.

## Conventions

- Match the surrounding code: module-level constants each with a one-line comment; private helpers start with `_`; docstrings and comments say WHY, naming the failure that motivated the code.
- numpy arrays in `raster.py` stay float64.
- New tests go in the file the task names, use `check(...)` for every assertion, and are added to that file's `__main__` list right after the test named in the task.
- Expected values in tests are the hand-derived ones written in this plan, never values copied from a run.
- No network access in any task.

---

# Phase 1 — ghost the environment layer

### SG1 — See-through composites the whole environment layer at its own opacity

**Why:** DC5's see-through blends only `env_wins & (qm.zbuf != -inf)` — a measured 688 of 675,000 pixels (0.10%) at the default view on `ts/optts` — so pressing `x` looks like nothing happens. The environment layer is the host, and the intent is that it becomes a ghost so the guest reads through it.
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/tui/raster.py`, `test_monitor.py`.

**Do:**
1. In `render` (`raster.py`), replace the see-through branch (currently lines 489-498) with one that composites the whole environment layer at `SEE_THROUGH_ALPHA` against whatever lies behind it:
   ```python
       frame, env_wins = _composite(env, qm)
       if view.see_through:
           # See-through ghosts the ENTIRE environment layer, not only the
           # pixels it wins over a guest surface: at the default view that
           # strict-win set is 0.10% of the frame, so blending only it looked
           # like a no-op. Each environment pixel is composited against what
           # lies behind it -- the guest's colour where a guest surface is
           # behind, the background elsewhere -- so the host fades to a faint
           # skeleton and a guest atom behind a host stick reads through it.
           # The depth is left as the environment's, so the halo still rings
           # the host that is in front; the owner changes only where a guest
           # really is behind, so a covered guest atom keeps its label.
           qm_behind = qm.zbuf != -np.inf
           behind = np.where(qm_behind[..., None], qm.color,
                             np.asarray(BACKGROUND, dtype=float))
           frame.color[env_wins] = (SEE_THROUGH_ALPHA * env.color[env_wins]
                                    + (1.0 - SEE_THROUGH_ALPHA) * behind[env_wins])
           frame.owner[env_wins & qm_behind] = qm.owner[env_wins & qm_behind]
   ```
   Nothing else in `render` changes: the background fill, halos, labels and distances keep their order. Where `env_wins` is false the pixel is the guest's, untouched.
2. Change the `SEE_THROUGH_ALPHA` comment (line 61) to: `# the environment layer's opacity when see-through ghosts it: against the guest where the guest is behind, the background elsewhere`.
3. Add `test_see_through_ghosts_the_environment()` to `test_monitor.py` and add it to the `__main__` list right after `test_see_through_shows_the_guest_behind_the_host()`. It has exactly 3 `check` calls:
   - Atoms, `qm_atom_indices`, `size_px` and view are exactly `test_host_fog_keeps_near_host_bright`'s: `[("C", 3.0, -3.0, 0.0), ("C", 3.0, -1.6, 0.0), ("C", -3.0, 1.6, 0.0), ("C", -3.0, 3.0, 0.0), ("O", 0.0, 0.0, 2.5)]`, `qm_atom_indices={4}`, `size_px=(400, 300)`, `View(elev=0, azim=0, fog=False, show_labels=False)`. The near host stick (atoms 0-1) and the far host stick (atoms 2-3) have no guest behind their midpoints; atom 4 is the guest oxygen, at its own pixel. `fog=False` because the ghost blends against the raw background.
   - Render `off` with that view and `on` with `View(elev=0, azim=0, fog=False, show_labels=False, see_through=True)`, both as `np.asarray(raster.render(...), float)`.
   - Get the two stick midpoints from `raster.project(atoms, view, (400, 300))`, as `test_host_fog_keeps_near_host_bright` does.
   - Check the near stick's midpoint, the far stick's midpoint and their ghost: `abs(on[y, x, c] - (0.3 * off[y, x, c] + 0.7 * 30.0)) <= 2` for every channel `c` (one `check` per stick).
   - Check the guest oxygen's own pixel is unchanged: `max(abs(on[y, x] - off[y, x])) <= 2` (the guest layer is not the environment, so no ghost and no blend where the guest wins).
   - The comment above the checks must give the derivation: with see-through off an environment pixel is its full colour; with it on the environment layer is composited at `SEE_THROUGH_ALPHA = 0.3` against `BACKGROUND = (30, 30, 30)`, so each channel is `0.3 * off + 0.7 * 30`; the guest oxygen is in the guest layer, not the environment, so its colour is unchanged.

**Acceptance:**
- `test_see_through_ghosts_the_environment` passes all 3 checks.
- `test_see_through_shows_the_guest_behind_the_host` passes unchanged (2 checks): the guest-behind-a-host blend at `0.3 * host + 0.7 * guest` is unchanged.
- `test_host_fog_keeps_near_host_bright`, `test_halos_widen_with_the_depth_gap` and `test_boundary_host_atoms_get_a_ball` pass unchanged.
- `test_monitor.py`: 136 `pass`, `0 failure(s)`. `test_orcamon.py`: 157 `pass`, `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.
- The diff touches only `raster.py` and `test_monitor.py`; `git diff HEAD -- '*test*'` shows one added test and its `__main__` entry and no change to any existing assertion or expected value.

## Out of scope

1. Changing `SEE_THROUGH_ALPHA` from 0.3, or exposing its value as a CLI option (D2).
2. Any change to `tools/orcamon/README.md`, the root `README.md`, the `--see-through` help string or `.claude/skills/orcamon/SKILL.md` (D4).
3. The four depth cues of `PLAN_ORCAMON_DEPTH_CUES.md` (fog, halos, boundary balls, the per-layer compositing).
4. The braille text renderer (`tui/geometry_text.py`); `x` still just toggles the flag and title there.
5. Running or maintaining `run_gates.py`.
6. Making the see-through state persist across job selections, or binding it to anything other than `x`.

## Progress template

```markdown
# orcamon — see-through ghosts the environment layer — Progress

Plan: `PLAN_ORCAMON_SEE_THROUGH_GHOST.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| SG1 Ghost the environment layer | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics |
|---|---|---|---|
| baseline | 133 | 157 | 0 failure(s) |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
