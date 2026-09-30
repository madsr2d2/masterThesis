# orcamon geometry view: implementation plan

Written 2026-09-30 for an implementation agent. It makes orcamon's
geometry pane readable as a 3D structure, the way a molecular viewer is:
- correct occlusion, shaded spheres and outlines;
- depth fog;
- switchable representations (ball-and-stick, licorice, space-filling,
  wireframe);
- view helpers: reset, principal-axis views, gentle rocking;
- `orcamon snapshot`, to render any job's geometry to a PNG.

**Read sections 1–3 before touching any code.** Then do the phases in order.
Each phase ends committed, with its gates passing. Do not start a phase until
the previous one is committed.

---

## 1. Why, in one paragraph each

**The text view is a mirror image.** `core/geometry.camera_basis()` returns
a "right" vector that points to screen LEFT. This was checked against
matplotlib at three view angles, and the result was always the opposite
sign. The text renderer uses that vector as screen +x, so text mode draws
the molecule's enantiomer. The pixel views are drawn by matplotlib and are
correct, but their pan keys use the same vector. **This is a correctness
bug and is fixed first (Phase 0).**

**The pixel view hides depth.** `tui/geometry_render.py` is a matplotlib
3D scatter. It has three defects:
- It draws all bonds first and all spheres after. A bond passing IN FRONT
  of an atom is therefore painted behind it, and the eye reads contradictory
  depth.
- Spheres are flat discs with no lighting.
- Nothing distinguishes near atoms from far ones.

Real viewers use a depth buffer, shaded spheres, silhouette outlines and
depth fog. Those four things are most of what makes a VMD, PyMOL or Jmol
picture readable.

**Baseline to beat.** Measured 2026-09-30 on this machine with the current
renderer:

| structure | atoms | frame | render | palette PNG |
|---|---|---|---|---|
| C8 RC, QM/XTB (`geometry_r2scan3c-xtb/rc/displace_minus`) | 140 (20 QM) | 940x900 | 28 ms | 25 KB |
| same | | 329x315 (preview) | 24 ms | 7 KB |
| B97-3c TS (`geometry_b973c-xtb/ts`) | 137 (17 QM) | 940x900 | 26 ms | 25 KB |
| SOLVATOR cluster (`solv_apo50_b/job.solvator.xyz`) | 280 (no layers, all balls) | 940x900 | 204 ms | 75 KB |

---

## 2. Ground rules (read all of them)

- **Interpreter:** `.venv/bin/python`, never bare `python`. The command is
  `.venv/bin/orcamon`.
- **Headless TUI runs need `HERDR_ENV=0`.** Without it the app draws the
  molecule over the user's terminal pane. If it happens anyway:
  `.venv/bin/python -c "from orcamon.tui import herdr_graphics as h; h.clear('geometry')"`.
- **Never kill or signal a running ORCA process.** Reading `/proc` is fine.
- **Leave the user's uncommitted work alone.** Run `git status --short`
  first and write down which files are already modified. Never
  `git add -A`, `git add .` or `git commit -a`. Stage files by path, only
  the ones this plan changes. If a file you must change already has
  uncommitted edits, stop and ask the user.
- **Back up before any scripted edit:**
  `cp FILE /tmp/<scratch>/FILE.bak`. Prefer the Edit tool.
  **Never write `open(p, "w").write(f(open(p).read()))`**: the write-mode
  open runs first and empties the file. Read into a variable, then write.
- **Use a scratch cache:** set `XDG_CACHE_HOME=/tmp/<scratch>/xdg` on every
  manual `orcamon` run, so the user's `~/.cache/orcamon` is not touched.
- **Commit straight to `master`,** one commit per phase, message prefixed
  `orcamon:`. **Do not push.** The user decides when.
- **Gates:**
  - After every change: `.venv/bin/python test_monitor.py` and
    `.venv/bin/python test_orcamon.py`.
  - Before every commit: `.venv/bin/python run_gates.py`, the full suite
    (about 8 minutes; run it in the background and wait for it). All gates
    must pass. Do not use `--all`.
  - Put the gate count and time in the commit message.
- **The dependency rule is enforced by a test.** `orcamon.core` and
  `orcamon.cli` import ONLY the standard library. numpy and PIL may be used
  in `orcamon.tui` only, and only in modules the pixel paths import lazily.
  `test_orcamon.test_the_core_needs_only_the_standard_library` fails
  otherwise. Anything the text renderer needs (camera, colours, radii,
  principal axes) goes in `core/geometry.py`, in pure Python.
- **The duplicate guard.** `data/test_curve_metrics.py` fails if a
  top-level name in a ROOT `.py` file (including `test_monitor.py` and
  `test_orcamon.py`) is defined in another root or `data/` module. That
  includes private `_helpers` that are functions. Give every new test
  function and helper a unique, descriptive name.
- **orcamon is a general ORCA tool.** Put no thesis paths or names in its
  code or its shipped `SKILL.md`. Use real thesis structures only for manual
  checks and screenshots.
- **The shipped skill has gates.** Every subcommand must be named in
  `tools/orcamon/src/orcamon/SKILL.md`. Every fenced `orcamon ...` example
  there is RUN by a test, so do not put a command that writes files in a
  fenced block; mention it inline in backticks. After changing commands or
  flags, run `.venv/bin/orcamon skill install --project .` and commit
  `.claude/skills/orcamon/SKILL.md` with the change.
- **Style:** match the existing code. Comments say WHY, with the failure
  that motivated them. Use dataclasses. No new dependencies.

---

## 3. Facts about the current code

| what | where |
|---|---|
| camera | `core/geometry.py: camera_basis(elev, azim) -> (right, up)`. WRONG handedness, see Phase 0 |
| element colours | `core/geometry.py: ELEMENT_COLORS` (H, O, C, N only), `DEFAULT_ELEMENT_COLOR`, `BOND_CUTOFF = 1.7` |
| pixel renderer | `tui/geometry_render.py: render(atoms, elev, azim, dpi, show_distances, show_labels, zoom, pan, qm_atom_indices, size_px) -> PIL.Image` (RGB) and `frame_png(image) -> bytes` (256-colour palette PNG) |
| text renderer | `tui/geometry_text.py: View` dataclass, `bonds(atoms)`, `render(atoms, width, height, view, qm_atom_indices, bond_list) -> rich.text.Text` |
| widgets | `tui/app.py: RotatableGeometryImage` (shared state and keys), `KittyGeometryImage`, `HerdrGeometryImage`, `TextGeometry`. Keys are in `GEOMETRY_BINDINGS` plus `TextGeometry.BINDINGS` (`h`) |
| widget state | `elev=20, azim=-60, zoom=1.0, pan=(0,0,0), show_distances=False, show_labels=True`; text: `_hydrogens` |
| render calls | `KittyGeometryImage._push` (size `FULL_PX = (900, 750)` or `PREVIEW_SCALE` of it) and `HerdrGeometryImage._push` (pane pixel size) |
| pane title | `MonitorApp._geometry_title()` and `MonitorApp._title_geometry()` |
| existing keys (do not reuse) | app: `q r m`; chart: `left right end a`; geometry: arrows, ctrl+arrows, `d [ ] l`, `h` (text) |
| framing | the circumscribed sphere: `center = mean(coords) + pan`, `span = (max distance from mean + 0.5) / zoom`, and neither changes with rotation |
| multilayer | `qm_atom_indices` (0-based set) is the high-level region. Other atoms are the "environment" and are drawn as thin lines today |
| tests | `test_monitor.py` (parser and renderers), `test_orcamon.py` (CLI, TUI and skill). Both are root gates |

---

## 4. Target design

### 4.1 One view description, shared by both renderers

`core/geometry.py` gets a `View` dataclass (pure Python). It REPLACES
`geometry_text.View`: delete that class and import this one.

```python
REPRESENTATIONS = ("ball-and-stick", "licorice", "space-filling", "wireframe")
TEXT_REPRESENTATIONS = ("ball-and-stick", "wireframe")

@dataclass
class View:
    elev: float = 20.0
    azim: float = -60.0
    zoom: float = 1.0
    pan: tuple[float, float, float] = (0.0, 0.0, 0.0)
    representation: str = "ball-and-stick"
    fog: bool = True
    show_labels: bool = True
    show_distances: bool = False
    show_hydrogens: bool = True
```

The widget keeps these values as attributes. It builds a `View` with a new
method, `RotatableGeometryImage.view() -> View`, which adds the rocking
offset to `azim` (Phase 5). Every render call takes `view=self.view()`.

### 4.2 The pixel renderer

`tui/raster.py` (NEW; numpy and PIL) is a z-buffered sprite renderer.
`geometry_render.render(atoms, view, *, qm_atom_indices=None, size_px=(900, 750))`
becomes a thin wrapper around it. Keep `frame_png` unchanged. When
parity is reached, the matplotlib code is deleted (Phase 2).

### 4.3 New keys (geometry pane focused)

| key | action | modes |
|---|---|---|
| `v` | next representation | pixel: all four; text: ball-and-stick and wireframe |
| `f` | depth fog on/off | all |
| `h` | hydrogens on/off | all (today: text only) |
| `0` | reset the view | all |
| `p` | view along the next principal axis | all |
| `o` | rocking on/off | all |

The pane title names the view, e.g.
`geometry (herdr · licorice · fog) · cycle 12`.

### 4.4 A new command

`orcamon snapshot JOB -o FILE.png [--representation R] [--no-fog]
[--no-labels] [--size WxH] [--elev E] [--azim A]` renders the job's latest
geometry, or its file geometry when nothing was printed, to a PNG. It needs
the `images` extra and exits 2 with a clear message without it. It is how
every phase's result is checked by eye, and how the user shares a picture.

---

## Phase 0: fix the camera's handedness

**Goal:** screen right really is right, so no view is mirrored.

### Background (do not skip)

The eye sits at `+f` from the centre, where
`f = (cos e cos a, cos e sin a, sin e)` for elevation `e` and azimuth `a`.
This is matplotlib's `view_init` convention. The correct basis is:
- `R = normalize(cross(z_hat, f))`: screen right. When `f` is parallel to
  z, use `R = (1, 0, 0)`.
- `U = cross(f, R)`: screen up.

Check: `cross(R, U) == f` (the screen's right-handed normal points at the
viewer). The current code returns `cross(f, z_hat)`, which is `-R`.

### Tasks

1. In `core/geometry.py`, fix `camera_basis` to return `(R, U)` as above.
   Update its docstring to state the convention and the check.
2. Add `camera_forward(elev, azim) -> f` beside it; `geometry_text.render`
   currently recomputes it inline.
3. Find every caller: `grep -rn "camera_basis" tools/`. That is
   `geometry_text.render` and `RotatableGeometryImage._pan_by`.
4. Fix panning so the keys move the MOLECULE in the arrow's direction.
   `self.pan` shifts the camera target, and the picture moves the opposite
   way. So in `_pan_by`, `pan_right` must SUBTRACT `R * step` and `pan_up`
   must subtract `U * step`. Adjust the signs in `action_pan_*` or in
   `_pan_by` (only one of them) and add a comment saying why.

### Tests (in `test_monitor.py`)

- `test_the_view_is_not_mirrored`:
  - For `(elev, azim)` in `[(0, 0), (20, -60), (35, 110), (-40, 200)]`,
    assert `cross(R, U)` equals `camera_forward` within 1e-9.
  - Add a matplotlib cross-check, because matplotlib is still installed in
    this phase. Project `R` and the origin with `mpl_toolkits.mplot3d.proj3d`
    after `view_init`; the projected x of `R` must be greater than the
    origin's. Mark this sub-check with a comment saying it is deleted in
    Phase 2 along with matplotlib.
- `test_text_mode_draws_the_right_enantiomer`:
  - Place atoms `A=("O", 0,0,0)`, `B=("N", 0,1,0)` and `C=("S", 0,0,1)`.
  - Render with `View(elev=0, azim=0)` in 21x11 cells. The eye is on +x, so
    +y is screen right and +z is screen up.
  - Assert that the `N` glyph's column is greater than the `O` glyph's
    column, and that the `S` glyph's row is smaller than the `O` glyph's row.

**Acceptance:** tests pass. Also check by eye in text mode on a chiral
structure: its handedness matches the pixel pane in the same terminal.

**Commit:** `orcamon: the camera's right is right -- text mode drew the mirror image`.

---

## Phase 1: `orcamon snapshot`, and baseline pictures

**Goal:** render any job to a PNG from the command line, using the current
renderer, and keep "before" pictures to compare against later.

### Tasks

1. In `cli/__init__.py: build_parser()`, register `snapshot`, using the
   `add(...)` helper like the other commands:
   - help: `"render a job's geometry to a PNG (needs the images extra)"`
   - example: `"orcamon snapshot opt/ts -o ts.png"`
   - arguments:
     - `job` (`_job(p)`);
     - `-o/--output` (required, metavar `FILE`);
     - `--size` (default `"900x750"`, metavar `WxH`);
     - `--elev` (float, default 20);
     - `--azim` (float, default -60);
     - `--no-labels`;
     - `--distances`.
   - Representation and fog flags are added in Phases 3 and 4.
2. In `cli/commands.py`, add `cmd_snapshot(args, job)` with the
   `@_with_job` decorator.
   - Pick the geometry the same way `cmd_geom` does with no `--cycle`: the
     latest printed point with atoms, else `job.file_geometry`, else exit 2
     with `"no geometry to draw"`.
   - Import the renderer INSIDE the function:
     `from ..tui import geometry_render`. Catch `ImportError`, print
     `"snapshot needs the images extra: pip install 'orcamon[images]'"` and
     exit 2.
   - Parse `--size` with `re.fullmatch(r"(\d+)x(\d+)", ...)`; a bad value is
     a `UsageError`.
   - Render with the options, then write `geometry_render.frame_png(image)`
     to the file.
   - Print one line: `FILE: WxH, N atoms, <source>`, where `<source>` is
     `describe_point(point)` or `"latest geometry"`.
   - Exit with `_job_exit(job)`.
3. In `SKILL.md` (package), name it INLINE, not in a fenced block. For
   example, in "Looking closer": "`orcamon snapshot JOB -o view.png` renders
   the geometry to a PNG for a person to look at." Reinstall the project
   copy.
4. Render baseline pictures with the CURRENT renderer into
   `/tmp/<scratch>/view_baseline/`, NOT into the repository. Use these three
   structures:
   - `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/rc/displace_minus`
   - `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts`
   - `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/solv_apo50_b`

   Render each at `--size 900x750`, at the default view and at
   `--elev 60 --azim 30`.

### Tests (in `test_orcamon.py`)

`test_snapshot_writes_a_png`, using the `_Tree` fixture:
- `snapshot opt_done -o <tmp>/x.png` exits 0. The file starts with
  `b"\x89PNG"`, and PIL opens it at the requested size.
- A job with NO geometry anywhere exits 2 with "no geometry to draw". The
  fixture's jobs all have one: every input has an inline atom, which the
  input-geometry fallback draws. So create one inside the test:
  `t.root / "nothing" / "job.inp"` containing
  `"! XTB SP\n* xyzfile 0 1 missing.xyz\n"`, with no `.out` and no
  `missing.xyz`.
- `--size 12` exits 2.

The stdlib-only test must still pass, because the renderer import is lazy.

**Commit:** `orcamon: snapshot -- render a job's geometry to a PNG`.

---

## Phase 2: the z-buffered renderer (ball-and-stick parity)

**Goal:** replace matplotlib with a renderer that draws correct
occlusion, shaded spheres and cylinder bonds. It keeps today's features:
index labels, distances, and the thin environment layer.

### 2.1 Shared data in `core/geometry.py` (pure Python)

- Replace `ELEMENT_COLORS` with the Jmol/CPK colours for at least these
  elements. The text renderer's `TEXT_ELEMENT_COLORS` override for C stays.
  ```
  H #FFFFFF  C #909090  N #3050F8  O #FF0D0D  F #90E050  P #FF8000
  S #FFFF30  Cl #1FF01F Br #A62929 I #940094  B #FFB5B5 Si #F0C8A0
  Na #AB5CF2 K #8F40D4  Mg #8AFF00 Ca #3DFF00 Fe #E06633 Cu #C88033
  Zn #7D80B0
  ```
  `DEFAULT_ELEMENT_COLOR` stays `#c060c0`.
- Add `VDW_RADII` (Å, Bondi):
  ```
  H 1.10  C 1.70  N 1.55  O 1.52  F 1.47  P 1.80  S 1.80  Cl 1.75
  Br 1.85 I 1.98  B 1.92  Si 2.10 Na 2.27  K 2.75  Mg 1.73 Ca 2.31
  Fe 2.00 Cu 1.40 Zn 1.39
  ```
  and `DEFAULT_VDW_RADIUS = 1.70`.
- Add `bonds(atoms)`: move `geometry_text.bonds` here unchanged (grid
  binning, H–H excluded) and import it back into `geometry_text`. Both
  renderers must bond the same pairs.
- Add `REPRESENTATIONS`, `TEXT_REPRESENTATIONS` and `View` (section 4.1).

### 2.2 `tui/raster.py` (NEW)

Module-level constants, each with a one-line comment:

```python
BACKGROUND = (30, 30, 30)            # "#1e1e1e", as today
LIGHT = normalize((-0.4, 0.55, 0.73)) # camera space: from upper left, in front
AMBIENT, DIFFUSE, SPECULAR, SHININESS = 0.28, 0.72, 0.45, 40.0
BALL_SCALE = 0.28                    # ball radius = BALL_SCALE * vdW radius
BOND_RADIUS = {"ball-and-stick": 0.10, "licorice": 0.16, "wireframe": 0.04}
ENV_BOND_RADIUS = 0.06               # the environment layer's thin sticks
ENV_DESATURATE = 0.5                 # mix environment colours halfway to grey
BALL_BOND_RGB = (184, 184, 184)      # ball-and-stick bonds: light grey
SAMPLE_SPACING = 0.8                 # bond sample spacing, in bond radii
MAX_SAMPLES_PER_BOND = 64
```

**Projection** (`project(atoms, view, size_px) -> (sx, sy, depth, scale)`,
numpy arrays). This function is public because the tests use it.

1. `R, U = camera_basis(view.elev, view.azim)`,
   `F = camera_forward(view.elev, view.azim)`.
2. `center = mean(coords) + pan`. `radius = max |coords - mean(coords)|`,
   ignoring pan, as today.
3. `pad = the largest drawn radius` (for space-filling, the largest vdW
   radius among the atoms drawn; otherwise 0.5). `span = (radius + pad) / zoom`.
4. `W, H = size_px`, `scale = min(W, H) / (2 * span)`, in pixels per Å.
5. `p = coords - center`. `sx = W/2 + (p·R) * scale`,
   `sy = H/2 - (p·U) * scale` (screen y points down),
   `depth = p·F` (in Å; LARGER is NEARER).

**Sprites.** Every drawn primitive is a shaded sphere. A bond is a row of
spheres along its axis, which is how a cylinder with round caps looks and
is automatically correct under the z-buffer. Cache one "stamp" per integer
pixel radius:

```python
@functools.lru_cache(maxsize=512)
def _stamp(r: int):
    """mask, nz (the unit-sphere height 0..1) and shade (ambient + diffuse)
    and spec (the specular term), each (2r+1, 2r+1), for a sphere of pixel
    radius r lit by LIGHT and seen along +z."""
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1].astype(float)
    rr = max(r, 0.5)
    d2 = (xx * xx + yy * yy) / (rr * rr)
    mask = d2 <= 1.0
    nz = np.sqrt(np.clip(1.0 - d2, 0.0, 1.0))
    nx, ny = xx / rr, -yy / rr                    # screen y is down
    lx, ly, lz = LIGHT
    diffuse = np.clip(nx * lx + ny * ly + nz * lz, 0.0, 1.0)
    hx, hy, hz = normalize((lx, ly, lz + 1.0))    # half-vector, viewer on +z
    spec = np.clip(nx * hx + ny * hy + nz * hz, 0.0, 1.0) ** SHININESS
    return mask, nz, AMBIENT + DIFFUSE * diffuse, SPECULAR * spec
```

**Drawing a sphere** at float screen position `(cx, cy)`, depth `z0` (Å),
radius `ra` (Å), colour `rgb` (0..255 floats), and owner id `k` (the atom
index, used for labels and outlines):

```python
def _draw_sphere(zbuf, color, owner, cx, cy, z0, ra, scale, rgb, k):
    r = max(1, int(round(ra * scale)))
    mask, nz, shade, spec = _stamp(r)
    H, W = zbuf.shape
    x0, y0 = int(round(cx)) - r, int(round(cy)) - r      # top-left of the box
    x1, y1 = x0 + 2 * r + 1, y0 + 2 * r + 1               # exclusive
    # Clip the box to the image, and cut the stamp the same way.
    bx0, by0 = max(x0, 0), max(y0, 0)
    bx1, by1 = min(x1, W), min(y1, H)
    if bx0 >= bx1 or by0 >= by1:
        return
    sx0, sy0 = bx0 - x0, by0 - y0
    sx1, sy1 = sx0 + (bx1 - bx0), sy0 + (by1 - by0)
    z_view = zbuf[by0:by1, bx0:bx1]          # basic slices are VIEWS: writing
    c_view = color[by0:by1, bx0:bx1]         # into them writes the buffers
    o_view = owner[by0:by1, bx0:bx1]
    zs = z0 + nz[sy0:sy1, sx0:sx1] * ra
    m = mask[sy0:sy1, sx0:sx1] & (zs > z_view)
    z_view[m] = zs[m]
    c_view[m] = np.asarray(rgb) * shade[sy0:sy1, sx0:sx1][m][:, None] \
        + 255.0 * spec[sy0:sy1, sx0:sx1][m][:, None]
    o_view[m] = k
```

Write through the named slice views, as above. `zbuf[by0:by1, bx0:bx1][m] =
...` also works, because a basic slice is a view. But indexing with index
ARRAYS (`zbuf[ys, xs]` where ys and xs are arrays) makes a COPY, and writing
into that copy silently changes nothing. The occlusion test below catches
this.

**Buffers:** `zbuf = full((H, W), -inf)`, `color = zeros((H, W, 3), float)`,
`owner = full((H, W), -1, int)`.

**What is drawn.** Let `is_qm[i]` be `i in qm_atom_indices` when
`0 < len(qm) < n`, else True for all atoms.
- `view.show_hydrogens` False removes H atoms, and bonds to them, from
  everything.
- QM atoms, ball-and-stick:
  - A sphere per atom, radius `BALL_SCALE * vdW`, element colour.
  - A bond between two QM atoms is a row of samples with radius
    `BOND_RADIUS["ball-and-stick"]`, colour `BALL_BOND_RGB`, running from
    atom i's centre to atom j's.
  - Spacing: `SAMPLE_SPACING * bond_radius` in Å, at most
    `MAX_SAMPLES_PER_BOND`, always including both ends.
  - Sample depth is interpolated linearly between the two ends.
- Environment atoms, meaning any bond with at least one non-QM end: samples
  with radius `ENV_BOND_RADIUS`. The half nearer each atom takes that
  atom's colour, desaturated: `rgb = rgb*(1-d) + gray*d`, where `gray` is
  the mean of rgb and `d = ENV_DESATURATE`. Environment atoms get no ball.
- Draw order does not matter; the z-buffer decides. Draw environment first
  anyway, so ties go to the QM region.
- Draw the bond samples of a bond as one loop over positions, reusing one
  stamp (they share a radius).

**Labels** (`view.show_labels`, QM atoms only, as today):
- Draw after shading, with `PIL.ImageDraw.text` at `(sx, sy)`, anchor `"mm"`.
- Font: `ImageFont.load_default(size=max(9, round(H / 70)))`. This needs
  Pillow ≥ 10.1; update the extra's pin.
- White text with a 1-px black stroke (`stroke_width=1, stroke_fill="black"`).
- Draw a label only if its atom is at least partly visible: at least one
  pixel of `owner` equals `i` inside that atom's own sprite box (the same
  box `_draw_sphere` used, clipped to the image). A label must not float
  over an atom hidden behind another. That is the occlusion rule again, and
  the reason for the `owner` buffer.

**Distances** (`view.show_distances`): for each bond with both ends in the
QM region, draw `f"{d:.2f}"` at the projected midpoint, in yellow
`(255, 224, 102)` with a black stroke, using the same visibility rule
applied to either end's owner.

**Output:** fill background pixels (`zbuf == -inf`) with `BACKGROUND`, clip
to 0..255, convert to uint8, and return `Image.fromarray(..., "RGB")`.

**Performance rules:**
- No Python loop over PIXELS, ever. Loops over spheres and samples are fine.
- The bond list comes from `core.geometry.bonds`, which is linear time.
- Pass `bond_list=` in from the widget when it has one cached (the text
  widget already caches by `id(atoms)`); the pixel widgets should cache the
  same way.

### 2.3 Wire it in

1. `geometry_render.render(atoms, view=None, *, qm_atom_indices=None, size_px=(900, 750), bond_list=None)`
   calls `raster.render(...)`. An empty atom list returns a
   background-filled image of `size_px`.
2. Delete the matplotlib code from `geometry_render.py`:
   - `_new_figure`, `_to_image` and the old body;
   - `ELEMENT_RADII`, `DEFAULT_ELEMENT_RADIUS`, `BOX_ZOOM` and
     `_OVERLAY_ZORDER`.

   Keep `frame_png`. Delete the Phase 0 matplotlib cross-check in
   `test_monitor.py`.
3. Update both call sites in `tui/app.py` (`KittyGeometryImage._push`,
   `HerdrGeometryImage._push`) to `render(atoms, self.view(),
   qm_atom_indices=..., size_px=..., bond_list=...)`.
4. Add `RotatableGeometryImage.view()` (section 4.1) and make the text
   widget use it too (`geometry_text.render(atoms, w, h, self.view(), ...)`).
   `TextGeometry._show_hydrogens()` feeds `view.show_hydrogens`.
5. The old test `test_a_frame_is_a_small_faithful_palette_png` checks that
   no pyplot figure is left open. Replace that sub-check with "`render`
   imports no matplotlib": run it in a subprocess with a `sys.meta_path`
   finder refusing `matplotlib`, the same technique as
   `test_orcamon._STDLIB_ONLY_PROBE`. Keep its other checks: RGB, size,
   palette, and within about 1 level of 255.

### Tests (in `test_monitor.py`)

Use unique names.

- `test_near_atoms_hide_far_ones`: `View(elev=0, azim=0)` puts the eye on
  +x.
  - An `O` at (0,0,0) and a `C` at (2,0,0) overlap on screen. The centre
    pixel must be carbon-coloured (grey, R≈G≈B), not red.
  - Add an `N` at (0,0,0), and two `C` atoms at (2,−1,0) and (2,1,0) with
    the bond between them passing IN FRONT of the N. The centre pixel must
    be `BALL_BOND_RGB`-ish (grey), not blue.
- `test_spheres_are_shaded`: one atom; the pixel at the upper-left of its
  disc (toward `LIGHT`) is brighter than the lower-right one.
- `test_labels_follow_occlusion`: an atom hidden completely behind a bigger
  near one gets no label. Render with and without labels and compare the
  pixels around the hidden atom's projected position: they must be equal.
- `test_the_renderer_is_fast_enough`, with generous bounds so it does not
  flake:
  - a synthetic 140-atom structure with 20 QM atoms at 940x900 renders in
    under 150 ms;
  - its `frame_png` is under 80 KB.

  Put the measured numbers in the commit message.

### Acceptance

- Render the three structures from Phase 1 into
  `/tmp/<scratch>/view_phase2/` with `orcamon snapshot` and look at them.
  Bonds in front of atoms must be drawn in front, and spheres must look
  round.
- Measure and report the same table as section 1: time and PNG size, full
  and preview, for the three structures.
  - Targets: 140-atom full frame ≤ 40 ms, preview ≤ 12 ms, 280-atom cluster
    ≤ 150 ms.
  - If a target is missed, report it honestly. Do NOT add threads or C
    extensions; try a larger `SAMPLE_SPACING` (up to 1.2) first.
- Headless TUI with `--graphics text` still passes `test_orcamon`. Pixel
  modes cannot be checked headless here; say so in the commit.

**Commit:** `orcamon: a z-buffered renderer -- bonds in front are in front, spheres are lit`.

---

## Phase 3: depth fog and outlines

**Goal:** far atoms recede and overlapping atoms separate.

### Pixel (`tui/raster.py`)

- **Fog**, when `view.fog`:
  - Over foreground pixels, `zmin, zmax = zbuf.min(), zbuf.max()`.
  - `t = (zmax - zbuf) / max(zmax - zmin, 1e-6)`; 0 is nearest, 1 farthest.
  - `color = color * (1 - FOG * t) + BACKGROUND * (FOG * t)` with
    `FOG = 0.6`.
  - Apply fog BEFORE labels, so labels stay legible.
- **Outlines**, always on:
  - For each of the 4 neighbour shifts (use `np.roll`, or slicing with an
    edge pad), a pixel is an edge if the neighbour is background, OR the
    neighbour is nearer by more than `OUTLINE_JUMP = 0.35` Å
    (`zbuf_nb - zbuf > OUTLINE_JUMP`).
  - Edge pixels (foreground only) are multiplied by `OUTLINE_DARKEN = 0.25`.
  - This darkens the FAR side of every silhouette, so a near atom is ringed
    where it overlaps a far one.
  - Do the outline pass after fog and before labels.

### Text (`tui/geometry_text.py`)

- Track depth per cell. For every plotted bond dot, compute its depth by
  linear interpolation between the two atoms' depths, and keep the NEAREST
  depth per cell (`cell_depth[k] = max(...)`). Atom glyph cells take the
  atom's depth.
- When `view.fog`: split the scene's depth range into thirds. Cells in the
  far third get `dim=True` added to their style (Rich `Style(..., dim=True)`
  or `style + Style(dim=True)`). The near two thirds are unchanged. A
  terminal's background colour is unknown, so dimming is the only safe
  "fade".

### Keys and flags

- `("f", "toggle_fog", "Fog")` in `GEOMETRY_BINDINGS`. Its action flips
  `self.fog` and calls `self._input.request()`, like `toggle_labels`.
- `orcamon snapshot --no-fog`.
- The pane title shows `· fog` when fog is on (see Phase 4 for how the
  title is built).

### Tests

- `test_fog_dims_the_far_side` (`test_monitor.py`):
  - Two identical `C` atoms side by side on screen, one at depth +2 and one
    at −2.
  - With fog, the near one's brightest pixel is brighter than the far
    one's.
  - Without fog, they are equal within 2 levels.
- `test_outlines_separate_overlapping_atoms`: two same-colour atoms
  overlapping at different depths. Along the line between their centres
  there is at least one pixel darker than both centres' colours times 0.5.
- Text: `test_text_fog_dims_only_the_far_third` (`test_orcamon.py`): three
  atoms at depths −3, 0, +3. With fog, only the −3 atom's glyph cell has a
  dim style. Read styles from `Text.spans`, or render with
  `Console(record=True)`.

**Commit:** `orcamon: depth fog and outlines, so near and far read apart`.

---

## Phase 4: representations, and hydrogens everywhere

**Goal:** `v` cycles ball-and-stick → licorice → space-filling →
wireframe (pixel), or ball-and-stick → wireframe (text).

### Pixel (`tui/raster.py`), per QM atom and QM–QM bond

| representation | atom sphere radius | QM–QM bond | bond colour |
|---|---|---|---|
| ball-and-stick | `BALL_SCALE * vdW` | samples, r = 0.10 | `BALL_BOND_RGB` |
| licorice | 0.16 (cap) | samples, r = 0.16 | each half its atom's colour |
| space-filling | `vdW` | none | – |
| wireframe | 0.12, only atoms with NO bond (otherwise none) | samples, r = 0.04 | each half its atom's colour |

The environment layer is the same thin desaturated sticks in every
representation. Space-filling the host would bury the QM region that the
pane exists to show; say this in a code comment.

Framing `pad` (Phase 2.2 step 3) must use the largest radius of the
CURRENT representation, so space-filling never clips at the frame edge.

### Text (`tui/geometry_text.py`)

- ball-and-stick: today's behaviour.
- wireframe: bonds only. Atoms with no bond are drawn as `ATOM_DOT`;
  bonded atoms get no glyph.
- If `view.representation` is licorice or space-filling (the widget never
  sends these, but a caller might), render as ball-and-stick.

### Widget and keys (`tui/app.py`)

- `RotatableGeometryImage` gains `self.representation = "ball-and-stick"`,
  and a class attribute `SUPPORTED = REPRESENTATIONS`. `TextGeometry` sets
  `SUPPORTED = TEXT_REPRESENTATIONS`.
- `("v", "next_representation", "View")`: step to the next entry of
  `SUPPORTED`, wrapping around, then `self._input.request()` and
  `self.app.update_detail()`, so the title changes.
- Move `h` from `TextGeometry.BINDINGS` into `GEOMETRY_BINDINGS`, so it
  works in pixel mode too. State: `self._hydrogens: bool | None = None`
  (None means "by atom count": shown up to
  `geometry_text.SHOW_HYDROGENS_UP_TO` atoms in TEXT mode, always shown in
  pixel mode). Move `_show_hydrogens()` up to the base class, with the
  pixel/text default as a class attribute `HYDROGENS_BY_DEFAULT_UP_TO`
  (`None` for pixel, meaning always; 60 for text).
- **Title:** `MonitorApp._geometry_title()` becomes
  `f"geometry ({self.graphics}{note} · {g.representation}{' · fog' if g.fog else ''})"`,
  where `g = self.query_one("#geometry")`. Guard for the widget not being
  mounted yet (`on_mount` calls it) with `try/except NoMatches` or by
  building the title after `compose`.
- `orcamon snapshot --representation R` (choices `REPRESENTATIONS`,
  default ball-and-stick).

### Tests

- `test_representations_change_what_is_drawn` (`test_monitor.py`): water at
  a fixed view; count foreground pixels. Space-filling must cover more than
  ball-and-stick, which must cover more than wireframe, and licorice and
  ball-and-stick must differ.
- `test_space_filling_is_framed` (`test_monitor.py`): space-filling
  foreground pixels do not touch the image border at `zoom=1`.
- `test_orcamon.test_the_tui_runs_headless`: add two checks.
  - Focus `#geometry`, press `v`: the title contains `wireframe` (text mode
    has two representations).
  - Press `v` again: `ball-and-stick`.
  - Press `h`: the widget's hydrogen flag flips.
- `test_snapshot_writes_a_png`: add `--representation space-filling`.

**Commit:** `orcamon: v cycles representations; h works in every mode`.

---

## Phase 5: view helpers (reset, principal axes, rocking)

### 5.1 Reset (`0`)

`("0", "reset_view", "Reset")` sets `elev=20, azim=-60, zoom=1,
pan=(0,0,0)`, turns rocking off and clears the principal-axis cycle. It
leaves representation, fog, labels and hydrogens alone.

### 5.2 Principal axes (`p`), pure Python in `core/geometry.py`

```python
def principal_axes(coords) -> list[tuple[tuple[float, float, float], float]]:
    """(unit axis, variance) for the three principal axes of `coords`
    (unweighted), largest variance first. A 3x3 symmetric eigenproblem,
    solved by cyclic Jacobi rotations -- the core may not import numpy."""

def view_along(axis) -> tuple[float, float]:
    """(elev, azim) in degrees for an eye on +axis:
    elev = degrees(asin(z)), azim = degrees(atan2(y, x))."""
```

- **Jacobi:**
  - Start with `A` = the covariance matrix of the coordinates and `V = I`.
  - Repeat up to 50 sweeps: for each pair `(p, q)` in `(0,1), (0,2), (1,2)`
    with `|A[p][q]| > 1e-12`, compute
    `theta = 0.5 * atan2(2*A[p][q], A[q][q] - A[p][p])`, apply the rotation
    `J(p, q, theta)` as `A = Jᵀ A J` and `V = V J`.
  - Stop when the off-diagonal sum is below 1e-12.
  - The eigenvalues are `diag(A)`; the eigenvectors are `V`'s COLUMNS.
  - If you are unsure of the rotation's sign convention, test it: the result
    must satisfy `A_original · v ≈ λ v` for each pair. Write that test first.
- **The key:** the first press looks along the SMALLEST-variance axis, so
  the molecule lies face-on. The next presses step to the middle axis, then
  the largest, then wrap around.
  - Use the atoms currently shown: the QM region when there is one,
    otherwise all atoms.
  - Set `elev, azim = view_along(axis)` and `pan = (0, 0, 0)`.
  - Keep a counter `self._axis_step`, reset by `0` and by a new job
    selection.

### 5.3 Rocking (`o`)

- `("o", "toggle_rock", "Rock")`. When on, start
  `self._rock_timer = self.set_interval(1 / ROCK_FPS, self._rock_tick)`.
  When off, stop it and redraw once.
- Constants: `ROCK_FPS = 5` for the pixel widgets and 10 for `TextGeometry`;
  `ROCK_AMPLITUDE_DEG = 15`; `ROCK_PERIOD_S = 6`.
- `_rock_tick` sets
  `self._rock_offset = ROCK_AMPLITUDE_DEG * sin(2π * (monotonic() - self._rock_t0) / ROCK_PERIOD_S)`
  and triggers a redraw:
  - pixel widgets: `self._push(self._next_seq(), quality="full")`, NOT the
    preview-then-settle path, which would never settle;
  - text: `self.refresh()`.
- `view()` adds `self._rock_offset` to `azim`. The stored `self.azim` is not
  changed, so turning rocking off returns to where you were.
- Stop the timer in `on_unmount`. Rocking off is the default. Write a
  comment explaining why: a steady stream of frames, roughly 5 per second
  at up to 35 KB, is noticeable on a slow ssh link.

### Tests

- `test_principal_axes_face_a_plane_on` (`test_monitor.py`):
  - Six points on a ring of radius 1.4 Å in a tilted plane (normal
    `(1, 1, 1)/√3`).
  - `principal_axes`: the smallest-variance axis is parallel to the normal
    (`|dot| > 0.999`), and `A v = λ v` holds for all three axes.
  - Projecting with `view_along(smallest)` gives screen extents in x and y
    within 20% of each other: face-on.
- `test_orcamon.test_the_tui_runs_headless`: `0` resets `elev` and `azim`
  after an arrow press. `o` starts the rock timer, and `azim` in `view()`
  changes within 0.5 s while `self.azim` does not. `o` again stops it.

**Commit:** `orcamon: reset, principal-axis views and rocking`.

---

## Phase 6: documentation, extras, and the user's review

1. In `tools/orcamon/pyproject.toml`, set the `images` extra to
   `["numpy>=1.26", "pillow>=10.1"]` (matplotlib is no longer needed by
   orcamon). Leave `requirements.txt`'s matplotlib pin alone: other code in
   the repo uses it. Run `.venv/bin/pip install -e "./tools/orcamon[tui,images]"`.
2. `tools/orcamon/README.md`:
   - the key list gains `v f h 0 p o`;
   - the TUI section describes the representations, fog, outlines and
     rocking (with its bandwidth note);
   - document `orcamon snapshot`;
   - regenerate the pasted `--help` block for the changed commands.
3. In `SKILL.md` (package), keep it short: agents should not start the TUI.
   One inline sentence on `orcamon snapshot` for showing a person a
   structure. Reinstall the project copy.
4. Final pictures: render the three structures at the default view in all
   four representations, plus one with fog off, into
   `/tmp/<scratch>/view_final/`. Put the paths of the baseline and the final
   PNGs, and the performance table, in the commit message AND the final
   report to the user, so they can compare.
5. Tell the user what could NOT be checked here: the Kitty and herdr pixel
   paths in a real terminal. Ask them to try `v`, `f`, `p`, `o` and `0` in
   `orcamon tui computational/`.

**Commit:** `orcamon: document the new view; the images extra drops matplotlib`.

---

## 7. Out of scope (do not build)

- A QM-focus mode (hide or fog everything beyond N Å of the QM region).
- An axis triad.
- Perspective projection. Orthographic stays: bonds stay parallel under
  rotation.
- Supersampling or anti-aliasing beyond what outlines give.
- Stereo, ambient occlusion, cartoon or ribbon representations.
- Any change to the parser, the report model, or commands other than
  `snapshot`.

## 8. If something does not work

- **Occlusion test fails:** the boolean assignment went into a copy. Assign
  through a named slice view (`zv = zbuf[a:b, c:d]; zv[m] = ...`) and do
  the same for `color` and `owner`.
- **The picture is mirrored or upside down:** re-read Phase 0. Screen
  `y = H/2 - p·U`, with a MINUS sign, because image rows grow downward.
- **Space-filling is clipped at the edge:** the framing `pad` does not use
  the current representation's radii.
- **The stdlib-only test fails:** something in `core/` or `cli/` imports
  numpy or PIL at module level. Move the import into the function that
  needs it, or into `tui/`.
- **The duplicate-guard gate fails:** a new root-test helper shares a name
  with one in another root or `data/` module. Rename it.
- **A gate unrelated to orcamon fails:** do not "fix" it. Stop and report
  it; it is not caused by this work.
