# orcamon — depth cues for multilayer structures (DC1–DC6)

**Written:** 2026-10-01. **Verified against:** `master` at `f055638` (every expected value below was measured on a prototype of exactly these algorithms built on `f055638`'s `raster.py`; the parity probe in DC1 was run; suite counts and times were measured).
**Branch:** `master` (the agents commit straight to `master`; all earlier work was committed and pushed before the run, in `8961e1b` and `f85609a`).
**Progress file:** `PLAN_ORCAMON_DEPTH_CUES_PROGRESS.md` (template at the end).

orcamon's pixel renderer (`tools/orcamon/src/orcamon/tui/raster.py`) draws occlusion correctly, but in QM/XTB jobs the low-level layer ("host", drawn as thin sticks) often reads as wrongly in front of the QM layer ("guest", drawn as balls). It was checked pixel by pixel on `ts/optts`: of 423 host-bond points crossing a QM ball, the 57 drawn over it are 2–4 Å nearer and the 366 behind it are all hidden. The depth is right; the picture does not say which object is in front. This plan adds four cues that make it say so: per-layer depth fog that fades the far host hard and keeps the near host bright, dark halos whose width grows with the depth gap, small balls on host atoms bonded to the QM region, and a see-through toggle. The goal a deviation must still meet: a host stick that is in front of a guest atom must be visibly in front of it, and nothing about which pixel wins may change except where this plan says so.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 119 `pass` lines, last line `0 failure(s)` | ~1 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python test_orcamon.py` | 154 `pass` lines, last line `0 failure(s)` | ~6 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~17 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (same for `test_orcamon.py`). Each task states how many `pass` lines it adds.

Every task runs the first three commands. Only DC6 runs `run_gates.py`: tasks DC1–DC5 touch only `tools/orcamon/`, `test_monitor.py`, `test_orcamon.py` and `.claude/skills/orcamon/SKILL.md`, which no other gate reads except the duplicate guard.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- the DC1 parity probe prints `DIFFERS` for any case;
- an expected value written in this plan does not hold after the algorithm is implemented as written (report the measured value; do not change the expected value, the tolerance, or the algorithm to make it pass);
- `test_monitor.test_the_renderer_is_fast_enough` fails;
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance;
- in DC6, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's.
- Never `git add -A`, `git add .` or `git commit -a`; stage the files a task lists, by path.
- Never push.
- Never kill or signal a running ORCA process (reading `/proc` is fine).
- Never put thesis paths, job names or chemistry names into orcamon's package code (`tools/orcamon/src/`) or its shipped `tools/orcamon/src/orcamon/SKILL.md`. Tests use synthetic atoms only. Real structures appear only in the DC1 parity probe and the DC6 pictures, both outside the repository.
- Never run `orcamon` against the user's cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg`. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. numpy and PIL are allowed only in `tools/orcamon/src/orcamon/tui/`; `orcamon.core` and `orcamon.cli` import only the standard library at module level (a test enforces this).

## Scope

**In scope:** the four cues in `tui/raster.py`; the `see_through` view flag, its `x` key, its pane-title mark and its `snapshot --see-through` flag; tests; README and installed-skill updates.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | Build all four cues: per-layer depth fog, depth-gap halos, boundary host balls, see-through toggle. | user, 2026-10-01 |
| D2 | Run on `master` from a clean tree. All earlier uncommitted work was committed before the run (the CREST tooling in `8961e1b`; the C8 job directories in `f85609a`, with CREST scratch now in `.gitignore`). | user, 2026-10-01 |
| D3 | "Host" = the layer drawn as thin sticks (any bond with at least one non-QM end, plus DC4's boundary balls). "Guest" = QM–QM bonds and QM atom spheres. In a job that is not multilayer every atom is guest and the host layer is empty. In code the two layers are named `env` and `qm`, matching the existing names `ENV_BOND_RADIUS` and `qm_atom_indices`. | planner, 2026-10-01 |
| D4 | The frame is drawn into two layers, each with its own depth, colour and owner buffer, then composited per pixel. Where depths tie the QM layer wins (`env.zbuf > qm.zbuf`, strict). The old code's comment claimed ties went to QM, but its strict `>` gave them to the environment; parity was measured identical on 12 cases (DC1) regardless. | planner, 2026-10-01 |
| D5 | Fog is applied to each layer BEFORE compositing. QM: `mix = FOG * t` (unchanged, `FOG = 0.6`). Host: `mix = HOST_FOG * t ** HOST_FOG_POWER` with `HOST_FOG = 0.85`, `HOST_FOG_POWER = 2.0`. `t` is the existing fog parameter: `(near - z) / (near - far)` clipped to [0, 1], with `near, far` from `_fog_depths`. The `f` key still turns all fog off. | planner, 2026-10-01 (values prototyped) |
| D6 | Halos replace the 1-px depth outline. Background silhouettes stay 1 px. A foreground pixel is darkened when, within L1 distance `d` pixels, there is a pixel nearer by more than `T_d`. `T_d` runs linearly from `OUTLINE_JUMP` (0.35 Å) at `d = 1` to `HALO_FULL_GAP` (2.5 Å) at `d = levels`, where `levels = max(1, round(HALO_MAX_PX * min(width, height) / 750))` and `HALO_MAX_PX = 5`. Darkening reuses `OUTLINE_DARKEN = 0.25`. Halos run on the composited depth buffer. | planner, 2026-10-01 (values prototyped) |
| D7 | Boundary balls: in a multilayer job drawn as ball-and-stick, every visible non-QM atom bonded to a visible QM atom gets a sphere in the host layer, radius `BOUNDARY_BALL_SCALE * BALL_SCALE * vdW` with `BOUNDARY_BALL_SCALE = 0.6`, coloured `_desaturate(_rgb(element))`, owner = its atom index. Not in licorice, space-filling or wireframe. Boundary atoms get no label (labels stay QM-only). | planner, 2026-10-01 |
| D8 | See-through: `View.see_through: bool = False`. Where the host wins a pixel and a QM surface exists behind it, colour = `SEE_THROUGH_ALPHA * host + (1 - SEE_THROUGH_ALPHA) * qm` with `SEE_THROUGH_ALPHA = 0.3`, and the pixel's owner becomes the QM owner (so a covered QM atom keeps its label). The depth buffer keeps the host's depth, so halos still ring the host. Off by default; key `x`; `0` (reset) leaves it alone, as it does fog. | user (scope) and planner (design), 2026-10-01 |
| D9 | The braille text renderer (`tui/geometry_text.py`) gets none of these. `x` is a shared binding and still toggles the flag and the title in text mode; the text renderer ignores the flag. | planner, 2026-10-01 |
| D10 | Per-task gates: the three fast commands in § Verification commands. Full `run_gates.py` only in DC6. | planner, 2026-10-01 |
| D11 | Commit messages start with `orcamon:`; one commit per task. | planner, 2026-10-01 |

## Facts (verified 2026-10-01)

**Environment.** Repository root `/home/madsr2d2/masterThesis`. Interpreter `.venv/bin/python` (never bare `python`); the CLI is `.venv/bin/orcamon`. orcamon is installed editable from `tools/orcamon/src/orcamon/`, so edits take effect without reinstalling.

**`tools/orcamon/src/orcamon/tui/raster.py` at `f055638`** (last commit to touch it):
- Constants: `BACKGROUND = (30, 30, 30)`; `BALL_SCALE = 0.28`; `BOND_RADIUS = {"ball-and-stick": 0.10, "licorice": 0.16, "wireframe": 0.09}`; `LICORICE_RADIUS = 0.16`; `WIREFRAME_ATOM_RADIUS = 0.12`; `ENV_BOND_RADIUS = 0.06`; `ENV_DESATURATE = 0.5`; `BALL_BOND_RGB = (184, 184, 184)`; `SAMPLE_SPACING = 1.2`; `FOG = 0.6`; `OUTLINE_JUMP = 0.35`; `OUTLINE_DARKEN = 0.25`.
- `render(atoms, view=None, *, qm_atom_indices=None, size_px=(900, 750), bond_list=None) -> PIL.Image` (RGB). Order today: allocate `zbuf = full((h, w), -inf)`, `color = zeros((h, w, 3))`, `owner = full((h, w), -1, int)` → `project` → `_draw_geometry` → `foreground = zbuf != -inf` → fill background with `BACKGROUND` → if `view.fog`: `_apply_fog(zbuf, color, foreground, *_fog_depths(atoms, view))` → `_apply_outlines(zbuf, foreground, color)` → `Image.fromarray(clip(color, 0, 255).astype(uint8))` → `_draw_labels` → `_draw_distances`.
- `project(atoms, view, size_px) -> (sx, sy, depth, scale)` (numpy arrays; depth in Å, larger = nearer). Public; tests use it.
- `_draw_geometry(atoms, view, pairs, has_bond, sx, sy, depth, scale, is_qm, visible, zbuf, color, owner)`: draws every bond with a non-QM end at `ENV_BOND_RADIUS` in desaturated element colours, then QM–QM bonds (skipped in space-filling), then a sphere per visible QM atom with radius `_atom_radius(element, representation, has_bond[i])` (None = no sphere). Its only caller is `render`.
- `_draw_bond(zbuf, color, owner, a, b, length, radius, scale, rgb_a, rgb_b, owner_a, owner_b)` and `_draw_sphere(zbuf, color, owner, cx, cy, z0, radius, scale, rgb, atom_index)` write through slice views into whichever three buffers they are given.
- `_bond_owner(i, n) = n + i` is the owner of a bond half touching atom `i`.
- `_apply_fog(zbuf, color, foreground, near, far)`, `_fog_depths(atoms, view) -> (near, far)` (the framing sphere's front and back), `_apply_outlines(zbuf, foreground, color)`.
- `_rgb(element) -> np.ndarray` of floats, `_desaturate(rgb)`.
- Tests use only `raster.render`, `raster.project` and `raster.BACKGROUND`.

**`tools/orcamon/src/orcamon/core/geometry.py`**: `View` dataclass with fields `elev=20.0, azim=-60.0, zoom=1.0, pan=(0.0, 0.0, 0.0), representation="ball-and-stick", fog=True, show_labels=True, show_distances=False, show_hydrogens=True`. `VDW_RADII` (K 2.75, C 1.70, O 1.52). `ELEMENT_COLORS`: K `#8F40D4`, O `#FF0D0D`, C `#909090`.

**Camera** (`core.geometry.camera_basis`). With `View(elev=0, azim=0)` the eye is on +x: screen right is +y, screen up is +z, and depth is x minus the mean x. With `View(elev=90, azim=-90)` the eye is on +z: screen right is +x, screen up is +y.

**`tools/orcamon/src/orcamon/tui/app.py`**:
- `GEOMETRY_BINDINGS` is a list of `(key, action, label)` tuples and includes `("f", "toggle_fog", "Fog")`. `x` is not bound anywhere in the app.
- `RotatableGeometryImage.__init__` sets `self.fog = True` among the view fields.
- `RotatableGeometryImage.view()` builds the `View`.
- `action_toggle_fog` flips `self.fog`, then calls `self._input.request()` and `self.app.update_detail()`.
- `action_reset_view` resets only the camera and helpers.
- `MonitorApp._geometry_title()` builds `geometry ({graphics}{note} · {representation}{' · fog' if fog})`.

**`tools/orcamon/src/orcamon/cli/__init__.py`**: `build_parser()` registers `snapshot` with `-o/--output`, `--representation`, `--size`, `--elev`, `--azim`, `--no-labels`, `--no-fog`, `--distances`. **`cli/commands.py`**: `cmd_snapshot` builds `View(elev=args.elev, azim=args.azim, representation=args.representation, fog=not args.no_fog, show_distances=args.distances, show_labels=not args.no_labels)`.

**Skill.** `test_orcamon.py` runs `orcamon skill --check --project <repo>`, which fails when the parser has a flag that the installed `.claude/skills/orcamon/SKILL.md` lacks. `.venv/bin/orcamon skill install --project .` (from the repository root) regenerates that file.

**README.** `tools/orcamon/README.md` has:
- the TUI description (the paragraph starting "kitty or text accordingly. The pixel panes are z-buffered");
- the key list (starting "`[`/`]` zoom, ctrl+arrows pan");
- `## Saving a picture`;
- a fenced block under `### \`orcamon snapshot\``.

`COLUMNS=90 .venv/bin/orcamon snapshot --help` reproduces that fenced block byte for byte today.

**Tests.** Both root test files define `check(label, ok, detail="")`, which prints `  pass  <label>` or `  FAIL  <label>`. Test names must be unique across all root `.py` files and `data/*.py`; `data/test_curve_metrics.py` fails otherwise. None of the names this plan introduces exists today. `test_monitor.test_the_renderer_is_fast_enough` requires a synthetic 140-atom, 20-QM frame at 940x900 to render in under 150 ms.

**Prototype measurements** (`f055638` plus these algorithms, 900x750, `ts/optts`, 137 atoms of which 17 QM): committed renderer 46 ms per frame; all four cues 71 ms; with see-through 72 ms.

## Conventions

- Match the surrounding code: module-level constants each with a one-line comment; private helpers start with `_`; docstrings and comments say WHY, naming the failure that motivated the code (here: a host stick 2–4 Å in front of a QM ball read as a rendering error because nothing marked it as in front).
- numpy arrays in `raster.py` stay float64 except where a task says float32.
- New tests go in the file the task names, use `check(...)` for every assertion, and are added to that file's `__main__` list right after the test named in the task.
- Expected values in tests are the hand-derived ones written in this plan, never values copied from a run.
- No network access in any task.

---

# Phase 1 — layers and depth shading

### DC1 — Draw the host and guest into separate layers, then composite

**Why:** Every cue in this plan needs to know, per pixel, which layer won and what lies behind it. Today the frame has a single depth buffer, so the guest surface under a host stick is overwritten and lost.
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/tui/raster.py`.

**Do:**
1. Add a private dataclass in `raster.py`, with a docstring saying it holds one layer of the frame:
   ```python
   @dataclass
   class _Layer:
       zbuf: np.ndarray    # (h, w) float, -inf where nothing is drawn
       color: np.ndarray   # (h, w, 3) float
       owner: np.ndarray   # (h, w) int, -1 where nothing is drawn

       @classmethod
       def empty(cls, height: int, width: int) -> "_Layer":
           return cls(np.full((height, width), -np.inf), np.zeros((height, width, 3)),
                      np.full((height, width), -1, dtype=int))
   ```
   (`from dataclasses import dataclass`.)
2. Change `_draw_geometry`'s last three parameters `zbuf, color, owner` to `env: _Layer, qm: _Layer`. Draw every bond with at least one non-QM end into `env` (`_draw_bond(env.zbuf, env.color, env.owner, ...)`). Draw QM–QM bonds and the QM atom spheres into `qm`. Make no other change to what is drawn or how. Delete the comment "Environment first, so when a QM and an environment sample tie on depth the QM region wins the pixel"; the composite now decides ties.
3. Add `_composite(env: _Layer, qm: _Layer) -> tuple[_Layer, np.ndarray]`. It returns the composited layer and the boolean mask `env_wins = env.zbuf > qm.zbuf` (strict, so the QM layer wins a tie; § Decisions D4). Composite with `np.where(env_wins, env.zbuf, qm.zbuf)`, `np.where(env_wins[..., None], env.color, qm.color)` and `np.where(env_wins, env.owner, qm.owner)`.
4. In `render`: create `env = _Layer.empty(height, width)` and `qm = _Layer.empty(height, width)`, and call `_draw_geometry(..., env, qm)` when there are atoms. Then `frame, _env_wins = _composite(env, qm)`. From there use `frame.zbuf`, `frame.color` and `frame.owner` where the old code used `zbuf`, `color` and `owner`. The background fill, fog, outlines, labels and distances stay exactly as they are, in the same order.
5. Write the parity probe below to `/tmp/orcamon_layer_parity.py` (outside the repository) and run it from the repository root with `XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/python /tmp/orcamon_layer_parity.py`. It compares the working-tree renderer against `raster.py` at `f055638`, loaded from git:
   ```python
   """Pixel parity of the working-tree raster against raster.py at f055638."""
   import subprocess, types
   import numpy as np
   from orcamon.cli import build_parser, commands
   from orcamon.core.geometry import View
   from orcamon.tui import raster

   REPO = "/home/madsr2d2/masterThesis"
   source = subprocess.run(["git", "-C", REPO, "show", "f055638:tools/orcamon/src/orcamon/tui/raster.py"],
                           capture_output=True, text=True, check=True).stdout
   source = source.replace("from ..core.geometry import", "from orcamon.core.geometry import")
   base = types.ModuleType("raster_f055638")
   exec(compile(source, "raster_f055638.py", "exec"), base.__dict__)

   C = "computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/"
   views = [View(), View(elev=60, azim=30, representation="licorice", fog=False),
            View(elev=-30, azim=200, representation="space-filling", show_distances=True),
            View(elev=10, azim=100, representation="wireframe", show_distances=True)]
   for job in ("geometry_b973c-xtb/ts/optts", "geometry_r2scan3c-xtb/pc", "solv_apo50_b"):
       args = build_parser().parse_args(["snapshot", C + job, "-o", "/dev/null", "--root", REPO])
       loaded = commands._load(args, C + job)
       point = commands._find_point(loaded.state, None, None)[0]
       atoms = point.atoms if point is not None else loaded.file_geometry.atoms
       qm = loaded.state.qm_atom_indices
       for view in views:
           a = np.asarray(base.render(atoms, view, qm_atom_indices=qm, size_px=(600, 500)))
           b = np.asarray(raster.render(atoms, view, qm_atom_indices=qm, size_px=(600, 500)))
           differ = int((a != b).any(-1).sum())
           print(f"{job} {view.representation}: {'identical' if differ == 0 else f'DIFFERS on {differ} px'}")
   ```

**If unsure:** If the probe cannot load a job (the paths are the user's thesis data and should exist), report it as NEEDS_USER rather than substituting another structure.

**Acceptance:**
- The parity probe prints 12 lines, every one ending in `identical` (3 jobs × 4 views: `ball-and-stick`, `licorice`, `space-filling`, `wireframe` for each of `geometry_b973c-xtb/ts/optts`, `geometry_r2scan3c-xtb/pc`, `solv_apo50_b`).
- `test_monitor.py`: 119 `pass`, `0 failure(s)`. `test_orcamon.py`: 154 `pass`, `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.
- The diff touches only `raster.py`, and `_draw_geometry` has exactly one caller (`render`).

### DC2 — Fog each layer separately: the near host stays bright, the far host fades hard

**Why:** The host is drawn uniformly desaturated, which the eye reads as "far" everywhere. That is exactly wrong for a host stick in front of the guest. With one linear fog for both layers, a host stick near the front is faded about 11% and one at the back about 48%, so nearness is barely visible.
**Depends on:** DC1.
**Files:** `tools/orcamon/src/orcamon/tui/raster.py`, `test_monitor.py`.

**Do:**
1. Add constants after `FOG`, each with a one-line comment: `HOST_FOG = 0.85` (how far the farthest host pixel fades to BACKGROUND) and `HOST_FOG_POWER = 2.0` (the host fade grows as t**2, so the near host stays bright).
2. Change `_apply_fog` to `_apply_fog(layer: _Layer, near: float, far: float, amount: float, power: float) -> None`. It computes `foreground = layer.zbuf != -np.inf`, returns if that is empty, sets `t = np.clip((near - layer.zbuf[foreground]) / max(near - far, 1e-6), 0.0, 1.0) ** power` and `mix = (amount * t)[:, None]`, then `layer.color[foreground] = layer.color[foreground] * (1 - mix) + np.asarray(BACKGROUND, float) * mix`. Rewrite the docstring to say why the layers are fogged separately (the reason in **Why**).
3. In `render`, between `_draw_geometry` and `_composite`, when `view.fog` and there are atoms: `near, far = _fog_depths(atoms, view)`, then `_apply_fog(qm, near, far, FOG, 1.0)` and `_apply_fog(env, near, far, HOST_FOG, HOST_FOG_POWER)`. Remove the fog call that ran on the composited frame. Background fill, outlines, labels and distances are unchanged.
4. Add `test_host_fog_keeps_near_host_bright()` to `test_monitor.py` and add it to the `__main__` list right after `test_fog_holds_still_under_rotation()`. It has exactly 3 `check` calls:
   - Atoms: `[("C", 3.0, -3.0, 0.0), ("C", 3.0, -1.6, 0.0), ("C", -3.0, 1.6, 0.0), ("C", -3.0, 3.0, 0.0), ("O", 0.0, 0.0, 2.5)]` with `qm_atom_indices={4}` and `size_px=(400, 300)`. That gives a near host bond at depth +3, a far host bond at depth −3, and a lone QM oxygen; nothing else is bonded.
   - Render twice, with `View(elev=0, azim=0, show_labels=False)` and with the same view plus `fog=False`, both as numpy float arrays. Get pixel positions from `raster.project(atoms, View(elev=0, azim=0), (400, 300))`.
   - Define `fade(x, y) = (off[y, x, 0] - on[y, x, 0]) / (off[y, x, 0] - 30.0)`, using the red channel at `(int(round(x)), int(round(y)))`.
   - Check the near host at the midpoint of atoms 0 and 1: `abs(fade - 0.027) < 0.02`.
   - Check the far host at the midpoint of atoms 2 and 3: `abs(fade - 0.555) < 0.02`.
   - Check the QM oxygen at atom 4's position: `abs(fade - 0.273) < 0.02`.

   The comment above the checks must give the derivation: extent = |(3, −3, −0.5)| + 0.5 = 4.772, so near = 4.772 and far = −4.772. The near stick's surface at depth ≈ 3.06 gives t = 0.179 and 0.85·t² = 0.027 (the old linear fog gives 0.108). The far stick at −2.94 gives t = 0.808 and 0.85·t² = 0.555. The O ball at 0.426 gives t = 0.455 and 0.6·t = 0.273.

**Acceptance:**
- `test_host_fog_keeps_near_host_bright` passes all 3 checks. On the prototype the measured fades were 0.036, 0.564 and 0.275.
- `test_monitor.py`: 122 `pass`, `0 failure(s)`. `test_orcamon.py`: 154 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- `test_fog_dims_the_far_side` and `test_fog_holds_still_under_rotation` still pass unchanged.

# Phase 2 — halos and boundary balls

### DC3 — Depth-gap halos: a dark rim on the far object, wider the larger the gap

**Why:** The only per-pixel occlusion cue today is a 1-px darkened line wherever a neighbour is nearer by more than 0.35 Å. A host stick 3 Å in front of a ball gets the same faint line as two atoms touching, so the stick does not read as floating in front. Measured on the committed renderer, a stick 0.6, 1.2 and 3.0 Å in front of a ball darkens 1, 1 and 1 rows of the ball on each side.
**Depends on:** DC2.
**Files:** `tools/orcamon/src/orcamon/tui/raster.py`, `test_monitor.py`.

**Do:**
1. Add constants after `OUTLINE_DARKEN`, each with a one-line comment: `HALO_MAX_PX = 5` (the widest halo, in pixels, on a frame whose short side is 750 px) and `HALO_FULL_GAP = 2.5` (the depth gap, Å, that earns the widest halo).
2. Replace `_apply_outlines` with `_apply_halos(zbuf, foreground, color)` and call it where `_apply_outlines` was called, on the composited frame (`frame.zbuf`, `foreground`, `frame.color`). Implement exactly this:
   ```python
   height, width = zbuf.shape
   edge = np.zeros_like(foreground)
   padded_fg = np.pad(foreground, 1, constant_values=False)
   for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):            # 1 px against the background
       edge |= foreground & ~padded_fg[1 + dy:1 + dy + height, 1 + dx:1 + dx + width]
   levels = max(1, round(HALO_MAX_PX * min(width, height) / 750))
   z = np.where(foreground, zbuf, -1e30).astype(np.float32)
   nearest = z.copy()          # after step d: the nearest depth within L1 distance d
   for d in range(1, levels + 1):
       padded = np.pad(nearest, 1, constant_values=-1e30)
       nearest = np.maximum.reduce([nearest, padded[0:height, 1:width + 1], padded[2:height + 2, 1:width + 1],
                                    padded[1:height + 1, 0:width], padded[1:height + 1, 2:width + 2]])
       threshold = (OUTLINE_JUMP if levels == 1
                    else OUTLINE_JUMP + (d - 1) * (HALO_FULL_GAP - OUTLINE_JUMP) / (levels - 1))
       edge |= foreground & (nearest - z > threshold)
   color[edge] *= OUTLINE_DARKEN
   ```
   The docstring must say: a pixel is darkened when something nearer by more than `T_d` lies within d pixels; `T_d` grows from `OUTLINE_JUMP` to `HALO_FULL_GAP` across the levels, so the halo's width tells the eye how far in front the near object is; the width scales with the frame's short side, so a 329x315 preview frame gets 2 levels, not 5. It must also give the reason in **Why**. The `-1e30` sentinel (not `-inf`) avoids `inf - inf = nan`.
3. Add `test_halos_widen_with_the_depth_gap()` to `test_monitor.py` and add it to the `__main__` list right after `test_host_fog_keeps_near_host_bright()`. It has exactly 6 `check` calls, two per gap:
   - For each `gap` in `(0.6, 1.2, 3.0)`: atoms `[("K", 0.0, 0.0, 0.0), ("C", 2.75 + gap, -0.7, 0.0), ("C", 2.75 + gap, 0.7, 0.0)]`, `qm_atom_indices={0}`, `View(elev=0, azim=0, representation="space-filling", fog=False, show_labels=False)`, `size_px=(750, 750)`. The K ball (vdW 2.75 Å) faces the eye. A host C–C stick runs across its centre `gap` Å in front of the ball's surface; the C atoms are ≥ 3.35 Å from K, so they are not bonded to it.
   - `sx, sy, _, _ = raster.project(atoms, view, (750, 750))`, `cx, cy = int(round(sx[0])), int(round(sy[0]))`, `im = np.asarray(raster.render(...), float)`, `ref = im[cy - 60, cx].sum()`.
   - For each row `y` in `cy - 14 .. cy + 14`, classify pixel `(r, g, b) = im[y, cx]`:
     - **K-hue** if `g < 0.6 * r` (K is `#8F40D4`, so g/r = 0.45, and darkening keeps that ratio);
     - **halo** if K-hue and `sum < 0.5 * ref`;
     - **stick** if not K-hue (the grey stick has g/r = 1).
   - Find the topmost and bottommost stick rows. Count consecutive halo rows directly above the topmost and directly below the bottommost.
   - Check that both counts equal the hand-derived width: 1 for gap 0.6, 2 for gap 1.2, 5 for gap 3.0. The thresholds are 0.35, 0.8875, 1.425, 1.9625 and 2.5 Å, and the width is the number of thresholds below the gap.

**Acceptance:**
- `test_halos_widen_with_the_depth_gap` passes all 6 checks with exactly 1, 2 and 5 rows each side. The prototype gave exactly these; the committed renderer gives 1, 1 and 1.
- `test_outlines_separate_overlapping_atoms`, `test_near_atoms_hide_far_ones` and `test_the_renderer_is_fast_enough` pass unchanged.
- `test_monitor.py`: 128 `pass`, `0 failure(s)`. `test_orcamon.py`: 154 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- `_apply_outlines` no longer exists (`grep -n "_apply_outlines" tools/orcamon/src/orcamon/tui/raster.py` prints nothing).

### DC4 — Small balls on host atoms bonded to the QM region

**Why:** A host atom in front of a QM ball is drawn only as the junction of thin sticks, so the eye sees a stick slicing through a ball rather than an atom in front of it. The case that prompted this is host C58 in front of O129 in a QM/XTB transition state. Most host atoms in the way are the ones bonded to the QM region.
**Depends on:** DC3.
**Files:** `tools/orcamon/src/orcamon/tui/raster.py`, `test_monitor.py`.

**Do:**
1. Add the constant `BOUNDARY_BALL_SCALE = 0.6` with the comment "a host atom bonded to the QM region: this fraction of its ball-and-stick ball".
2. In `_draw_geometry`, after the bonds and before the QM spheres, draw the boundary balls when all of these hold: `representation == "ball-and-stick"`, and some atom is not QM (`not all(is_qm)`). The boundary set is every `i` with `not is_qm[i]` and `visible[i]` that has a pair `(i, j)` or `(j, i)` in `pairs` with `is_qm[j]` and `visible[j]`. Build the set with one pass over `pairs`. Draw each one with `_draw_sphere(env.zbuf, env.color, env.owner, sx[i], sy[i], depth[i], BOUNDARY_BALL_SCALE * BALL_SCALE * VDW_RADII.get(element, DEFAULT_VDW_RADIUS), scale, _desaturate(_rgb(element)), i)`.
3. Update the comment "The environment never gets a ball" in `_draw_geometry`. The environment gets no full ball, for the reasons it already gives, except the atoms bonded to the QM region in ball-and-stick, which get a small desaturated one so a host atom in front of the guest reads as an atom.
4. Add `test_boundary_host_atoms_get_a_ball()` to `test_monitor.py` and add it to the `__main__` list right after `test_halos_widen_with_the_depth_gap()`. It has exactly 3 `check` calls:
   - Atoms `[("O", 0.0, 0.0, 0.0), ("C", 1.4, 0.0, 0.0), ("C", 2.8, 0.0, 0.0)]`, `qm_atom_indices={0}`, `size_px=(400, 400)`, `View(elev=90, azim=-90, fog=False, show_labels=False)`. The eye is on +z with x to the right. Atom 1 is a boundary atom; atom 2 is not.
   - Measure the column extent: the count of rows in image column `int(round(sx[i]))` whose pixel differs from `(30, 30, 30)` by more than 6 summed over channels.
   - Hand derivation: mean x = 1.4, extent = 1.4 + 0.5 = 1.9, scale = 400 / 3.8 = 105.26 px/Å. The boundary ball radius is 0.6 × 0.28 × 1.70 = 0.2856 Å = 30 px, so 61 rows. A host stick end of radius 0.06 Å = 6 px gives 13 rows.
   - Check that atom 1's column in ball-and-stick is 57–65 rows.
   - Check that atom 2's column in ball-and-stick is 11–15 rows.
   - Check that atom 1's column with `representation="licorice"` is 11–15 rows (no boundary ball outside ball-and-stick).

**Acceptance:**
- `test_boundary_host_atoms_get_a_ball` passes all 3 checks. The prototype measured 61, 13 and 13.
- `test_monitor.py`: 131 `pass`, `0 failure(s)`. `test_orcamon.py`: 154 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- In the diff, the boundary-ball drawing sits inside a condition containing both `representation == "ball-and-stick"` and `not all(is_qm)`. A non-multilayer job has every atom QM, so it is unaffected.

# Phase 3 — see-through, documentation

### DC5 — A see-through toggle: `x` in the TUI, `--see-through` on `snapshot`

**Why:** Even with good cues, a QM region sitting inside its host is partly covered from most angles. On the reported transition state, host atoms lie in front of the QM centre from 56% of viewing directions. A see-through mode shows the covered QM atoms while still showing the host as in front.
**Depends on:** DC4.
**Files:** `tools/orcamon/src/orcamon/core/geometry.py`, `tools/orcamon/src/orcamon/tui/raster.py`, `tools/orcamon/src/orcamon/tui/app.py`, `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/README.md`, `.claude/skills/orcamon/SKILL.md`, `test_monitor.py`, `test_orcamon.py`.

**Do:**
1. `core/geometry.py`: add `see_through: bool = False` as the last field of `View`.
2. `raster.py`: add the constant `SEE_THROUGH_ALPHA = 0.3` with the comment "how much of the host shows where it covers the QM region in see-through". In `render`, right after `frame, env_wins = _composite(env, qm)` (rename `_env_wins` to `env_wins`), when `view.see_through`:
   - `blend = env_wins & (qm.zbuf != -np.inf)`;
   - `frame.color[blend] = SEE_THROUGH_ALPHA * env.color[blend] + (1 - SEE_THROUGH_ALPHA) * qm.color[blend]`;
   - `frame.owner[blend] = qm.owner[blend]`.

   `frame.zbuf` keeps the host depth, so halos still ring the host. Add a comment giving both reasons: the owner change keeps a covered QM atom's label, and the kept depth keeps the halo.
3. `tui/app.py`:
   - In `GEOMETRY_BINDINGS`, add `("x", "toggle_see_through", "See-through")` right after the `"f"` entry.
   - In `RotatableGeometryImage.__init__`, add `self.see_through = False` next to `self.fog = True`.
   - In `view()`, pass `see_through=self.see_through`.
   - Add `action_toggle_see_through`, written like `action_toggle_fog`: flip the flag, call `self._input.request()`, then `self.app.update_detail()`.
   - In `MonitorApp._geometry_title()`, after the fog mark, add `" · see-through"` when `geometry.see_through`. The title then reads `geometry (… · ball-and-stick · fog · see-through)`.
   - `action_reset_view` must not touch `see_through`.
4. `cli/__init__.py`: in the `snapshot` parser, after `--no-fog`, add `p.add_argument("--see-through", action="store_true", help="show the QM region through the environment where the environment covers it")`. In `cli/commands.py` `cmd_snapshot`, pass `see_through=args.see_through` to `View(...)`.
5. Run `.venv/bin/orcamon skill install --project .` from the repository root. Then regenerate the fenced block under `### \`orcamon snapshot\`` in `tools/orcamon/README.md` with the exact output of `COLUMNS=90 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/orcamon snapshot --help`. Change no other help block.
6. Add `test_see_through_shows_the_guest_behind_the_host()` to `test_monitor.py` and add it to the `__main__` list right after `test_boundary_host_atoms_get_a_ball()`. It has exactly 2 `check` calls:
   - `front = [("O", 0.0, 0.0, 0.0), ("C", 2.0, -0.7, 0.0), ("C", 2.0, 0.7, 0.0)]` is a host stick 2 Å in front of a QM oxygen.
   - `back` is the same with the C atoms at x = −2.0, the stick behind the oxygen. It has the same framing: same scale, and the O at the same pixel.
   - Use `qm_atom_indices={0}`, `size_px=(400, 300)` and `View(elev=0, azim=0, fog=False, show_labels=False)`. The pixel is O's projected centre.
   - Take `c_env` from rendering `front`, `c_qm` from rendering `back`, and `c_x` from rendering `front` with `see_through=True`, all at that pixel.
   - Check that every channel of `c_x` is within 2 of `0.3 * c_env + 0.7 * c_qm`. The prototype measured (76, 76, 76), (212, 16, 16) and (171, 34, 34).
   - Check that `c_x[0] > c_env[0] + 50` (visibly red).
7. In `test_orcamon.py`:
   - In `test_the_tui_runs_headless`, right after the check `"h flips the hydrogen flag"`, press `x`, `await pilot.pause(0.2)`, then check `geometry.see_through is True and "see-through" in geometry.border_title`. Press `x` again, pause, and check `geometry.see_through is False and "see-through" not in geometry.border_title`. These are 2 checks.
   - In `test_snapshot_writes_a_png`, after the `--representation space-filling` check, run `snapshot opt_done -o <tmp>/see.png --size 120x100 --see-through` and check exit 0 and that the file starts with `b"\x89PNG"`. That is 1 check.

**Acceptance:**
- `test_monitor.py`: 133 `pass`, `0 failure(s)`. `test_orcamon.py`: 157 `pass`, `0 failure(s)`; this includes `orcamon skill --check` passing against the reinstalled skill. `data/test_curve_metrics.py`: `0 failure(s)`.
- `COLUMNS=90 XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg .venv/bin/orcamon snapshot --help` equals the README's fenced snapshot block byte for byte, and the diff of `tools/orcamon/README.md` touches no other block.
- `git diff` for this task shows `.claude/skills/orcamon/SKILL.md` gaining `--see-through`.

### DC6 — Document the cues, run every gate, render the comparison pictures

**Why:** The README still describes the 1-px outline and uniform fog, and the user needs pictures to judge the result.
**Depends on:** DC5.
**Files:** `tools/orcamon/README.md`.

**Do:**
1. In `tools/orcamon/README.md`, in the paragraph starting "kitty or text accordingly. The pixel panes are z-buffered", replace "distant atoms fade (fog) and overlapping silhouettes are outlined" with a sentence saying four things:
   - distant atoms fade, and the low-level (environment) layer fades harder with depth, so its near parts stay bright;
   - a nearer object casts a dark halo on what is behind it, wider the larger the depth gap;
   - environment atoms bonded to the QM region are drawn as small balls in ball-and-stick;
   - `x` shows the QM region through the environment.
2. In the key list paragraph starting "`[`/`]` zoom, ctrl+arrows pan", add "`x` shows the QM region through the environment (pixel panes)" after the fog key.
3. In `## Saving a picture`, add `--see-through` to the list of options.
4. Run `run_gates.py` (§ Verification commands; at least 900000 ms timeout).
5. Render pictures into `/tmp/orcamon_depth_cues/pictures/`, outside the repository, from the repository root with `XDG_CACHE_HOME=/tmp/orcamon_depth_cues/xdg`:
   - for `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts`: `orcamon snapshot <job> -o …/optts_default.png`, then the same with `--see-through` (`optts_see_through.png`), with `--no-fog` (`optts_no_fog.png`), and with `--representation licorice` (`optts_licorice.png`);
   - for `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/pc`: `pc_default.png`.
6. Measure render time for `ts/optts` at 900x750 with `geometry_render.render(atoms, View(), qm_atom_indices=qm, size_px=(900, 750))`, as the minimum of 3 runs. Measure it again at 329x315. Put both in the commit message beside the prototype's 71 ms.

**Acceptance:**
- `run_gates.py`: `38 gates`, `0 failed`.
- `test_monitor.py`: 133 `pass`. `test_orcamon.py`: 157 `pass`.
- The five PNGs exist in `/tmp/orcamon_depth_cues/pictures/`; list them with sizes.
- The commit message lists the picture paths and the two timings, and says that the Kitty and herdr pixel panes could not be checked headless.
- The README diff touches only the three places named in **Do**.

---

## Out of scope

1. Real cylinder shading for bonds (today a bond is a row of sphere samples, which looks beaded up close).
2. Performance work: per-frame fixed cost (about 28 ms at 900x750 in full-image passes) and batching the sphere draws.
3. A rotation-matrix camera (roll, rocking about the screen's vertical, principal-axis views with the long axis horizontal).
4. Label and distance-label collision avoidance.
5. Any of these cues in the braille text renderer.
6. Perspective projection and ambient occlusion.
7. Boundary balls in licorice, space-filling or wireframe.
8. `snapshot --zoom` and `snapshot --no-hydrogens`.
9. Changing `FOG` for the QM layer, or the colours of either layer.

## Progress template

```markdown
# orcamon — depth cues for multilayer structures — Progress

Plan: `PLAN_ORCAMON_DEPTH_CUES.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| DC1 Two layers, composited | TODO | | | |
| DC2 Per-layer depth fog | TODO | | | |
| DC3 Depth-gap halos | TODO | | | |
| DC4 Boundary host balls | TODO | | | |
| DC5 See-through toggle | TODO | | | |
| DC6 Docs, gates, pictures | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 119 | 154 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
