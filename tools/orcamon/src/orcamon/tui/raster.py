"""The geometry pane's pixels as a z-buffered sprite renderer.

Every drawn thing is a shaded sphere. A bond is a row of spheres along its
axis, which is how a cylinder with round caps looks and is automatically
correct under the z-buffer -- where the old matplotlib scatter was not: it
drew all bonds first and all spheres after, so a bond passing IN FRONT of an
atom was painted behind it and the eye read contradictory depth. Now near
spheres hide far ones, spheres are lit, and the same depth buffer tells a
label whether its own atom is at all visible.

numpy and PIL only. No Python loop ever walks pixels: the loops are over
atoms and over a bond's samples, and each sample writes a small cached stamp
into the buffers through a slice VIEW (basic slicing yields a view; indexing
with index arrays yields a copy and would silently change nothing).
"""
from __future__ import annotations

import functools
import math
from dataclasses import dataclass

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from ..core.geometry import (
    DEFAULT_ELEMENT_COLOR, DEFAULT_VDW_RADIUS, ELEMENT_COLORS, VDW_RADII, View, bonds,
    camera_basis, camera_forward,
)


def _normalize(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v) if n else v


BACKGROUND = (30, 30, 30)             # "#1e1e1e", as today
LIGHT = _normalize((-0.4, 0.55, 0.73))  # camera space: upper left, in front
AMBIENT, DIFFUSE, SPECULAR, SHININESS = 0.28, 0.72, 0.45, 40.0
BALL_SCALE = 0.28                     # ball radius = BALL_SCALE * vdW radius
# Wireframe's QM sticks must stay THICKER than ENV_BOND_RADIUS. They were
# 0.04 against the environment's 0.06, so wireframe drew the reactive centre
# as the faintest thing in the pane -- the reverse of what the layer is for.
BOND_RADIUS = {"ball-and-stick": 0.10, "licorice": 0.16, "wireframe": 0.09}
LICORICE_RADIUS = 0.16                # licorice atom cap
WIREFRAME_ATOM_RADIUS = 0.12          # the dot an unbonded atom keeps
ENV_BOND_RADIUS = 0.06                # the environment layer's thin sticks
ENV_DESATURATE = 0.5                  # mix environment colours halfway to grey
BOUNDARY_BALL_SCALE = 0.6             # a host atom bonded to the QM region: this fraction of its ball-and-stick ball
BALL_BOND_RGB = (184, 184, 184)       # ball-and-stick bonds: light grey
SAMPLE_SPACING = 1.2                  # bond sample spacing, in bond radii
MAX_SAMPLES_PER_BOND = 64
LABEL_RGB = (255, 255, 255)
DISTANCE_RGB = (255, 224, 102)
FOG = 0.6                             # how far the farthest pixel fades to BACKGROUND
HOST_FOG = 0.85                       # how far the farthest host pixel fades to BACKGROUND
HOST_FOG_POWER = 2.0                  # the host fade grows as t**2, so the near host stays bright
OUTLINE_JUMP = 0.35                   # angstrom a nearer neighbour must be to ring
OUTLINE_DARKEN = 0.25                 # silhouette pixels are multiplied by this
HALO_MAX_PX = 5                       # the widest halo, in pixels, on a frame whose short side is 750 px
HALO_FULL_GAP = 2.5                   # the depth gap, angstrom, that earns the widest halo
SEE_THROUGH_ALPHA = 0.3               # how much of the host shows where it covers the QM region in see-through


@dataclass
class _Layer:
    """One layer of the frame, drawn on its own before compositing.

    The two layers shared a single depth buffer, so the guest surface under a
    host stick was overwritten and lost: nothing could look at a pixel and say
    which layer won it and what lies behind. Each layer now keeps its own
    depth, colour and owner; `_composite` decides which one owns a pixel.
    """

    zbuf: np.ndarray    # (h, w) float, -inf where nothing is drawn
    color: np.ndarray   # (h, w, 3) float
    owner: np.ndarray   # (h, w) int, -1 where nothing is drawn

    @classmethod
    def empty(cls, height: int, width: int) -> "_Layer":
        return cls(np.full((height, width), -np.inf), np.zeros((height, width, 3)),
                   np.full((height, width), -1, dtype=int))


def _rgb(element: str) -> np.ndarray:
    colour = ELEMENT_COLORS.get(element, DEFAULT_ELEMENT_COLOR)
    return np.array([int(colour[i:i + 2], 16) for i in (1, 3, 5)], dtype=float)


def _desaturate(rgb: np.ndarray) -> np.ndarray:
    grey = rgb.mean()
    return rgb * (1.0 - ENV_DESATURATE) + grey * ENV_DESATURATE


def _framing_pad(atoms: list, centered: np.ndarray, view: View) -> float:
    """How much to add to the molecule's own radius so nothing clips.

    Space-filling draws each atom at its full vdW radius, so the pad is
    `max(|p_i| + r_i) - max(|p_i|)` over the atoms -- the ACTUAL extent, not
    merely the largest radius: an atom at the frame edge is what clips, and
    the largest atom is not always the one there. (The plan's step 3 said
    "the largest radius", which for a structure whose extreme atom is small
    leaves the pad too small and clips.) Every other representation needs
    only half an angstrom of slack."""
    distances = np.linalg.norm(centered, axis=1)
    radius = float(distances.max())
    if view.representation == "space-filling":
        radii = np.array([VDW_RADII.get(el, DEFAULT_VDW_RADIUS) for el, *_ in atoms])
        return float((distances + radii).max() - radius)
    return 0.5


def _frame_sphere(atoms: list, view: View):
    """`(coords, center, extent)`: the atoms as an array, the point the camera
    looks at (`mean + pan`) and the radius, in angstrom and before zoom, of
    the sphere around the MEAN that holds everything drawn. Neither the
    sphere nor its radius changes with rotation, which is what both the
    framing and the fog need from it."""
    coords = np.array([(x, y, z) for _el, x, y, z in atoms], dtype=float)
    mean = coords.mean(axis=0)
    centered = coords - mean
    center = mean + np.asarray(view.pan, dtype=float)
    extent = float(np.linalg.norm(centered, axis=1).max()) + _framing_pad(atoms, centered, view)
    return coords, center, extent


def project(atoms: list, view: View, size_px: tuple[int, int]):
    """Screen x, screen y, depth (angstrom, larger is NEARER) and pixels per
    angstrom for every atom. Public because the tests read it.

    The frame is the CIRCUMSCRIBED sphere -- centre `mean + pan`, radius the
    farthest atom -- so no rotation can push an atom out of the pane, and
    neither changes with rotation."""
    coords, center, extent = _frame_sphere(atoms, view)
    right, up = camera_basis(view.elev, view.azim)
    forward = camera_forward(view.elev, view.azim)
    right, up, forward = np.array(right), np.array(up), np.array(forward)
    span = extent / view.zoom
    width, height = size_px
    scale = min(width, height) / (2 * span) if span > 0 else 1.0
    p = coords - center
    sx = width / 2 + (p @ right) * scale
    sy = height / 2 - (p @ up) * scale  # screen y points down
    depth = p @ forward
    return sx, sy, depth, scale


@functools.lru_cache(maxsize=512)
def _stamp(radius: int):
    """`mask`, `nz` (the unit-sphere height 0..1) and the `shade` and `spec`
    terms, each (2r+1, 2r+1), for a sphere of pixel radius r lit by LIGHT and
    seen along +z. Cached per integer radius: a bond reuses one stamp for all
    of its samples."""
    yy, xx = np.mgrid[-radius:radius + 1, -radius:radius + 1].astype(float)
    rr = max(radius, 0.5)
    d2 = (xx * xx + yy * yy) / (rr * rr)
    mask = d2 <= 1.0
    nz = np.sqrt(np.clip(1.0 - d2, 0.0, 1.0))
    nx, ny = xx / rr, -yy / rr                     # screen y is down
    lx, ly, lz = LIGHT
    diffuse = np.clip(nx * lx + ny * ly + nz * lz, 0.0, 1.0)
    hx, hy, hz = _normalize((lx, ly, lz + 1.0))    # half-vector, viewer on +z
    spec = np.clip(nx * hx + ny * hy + nz * hz, 0.0, 1.0) ** SHININESS
    return mask, nz, AMBIENT + DIFFUSE * diffuse, SPECULAR * spec


def _draw_sphere(zbuf, color, owner, cx, cy, z0, radius, scale, rgb, atom_index):
    r = max(1, int(round(radius * scale)))
    mask, nz, shade, spec = _stamp(r)
    height, width = zbuf.shape
    x0, y0 = int(round(cx)) - r, int(round(cy)) - r      # top-left of the box
    x1, y1 = x0 + 2 * r + 1, y0 + 2 * r + 1              # exclusive
    bx0, by0 = max(x0, 0), max(y0, 0)
    bx1, by1 = min(x1, width), min(y1, height)
    if bx0 >= bx1 or by0 >= by1:
        return
    sx0, sy0 = bx0 - x0, by0 - y0
    sx1, sy1 = sx0 + (bx1 - bx0), sy0 + (by1 - by0)
    z_view = zbuf[by0:by1, bx0:bx1]          # basic slices are VIEWS: writing
    c_view = color[by0:by1, bx0:bx1]         # into them writes the buffers
    o_view = owner[by0:by1, bx0:bx1]
    zs = z0 + nz[sy0:sy1, sx0:sx1] * radius
    near = mask[sy0:sy1, sx0:sx1] & (zs > z_view)
    z_view[near] = zs[near]
    c_view[near] = (np.asarray(rgb) * shade[sy0:sy1, sx0:sx1][near][:, None]
                    + 255.0 * spec[sy0:sy1, sx0:sx1][near][:, None])
    o_view[near] = atom_index


def _draw_bond(zbuf, color, owner, a, b, length, radius, scale, rgb_a, rgb_b,
               owner_a, owner_b):
    """A bond is a row of spheres of `radius` from screen point `a` to `b`
    (`(x, y, depth)`), spaced SAMPLE_SPACING * radius apart along its true 3D
    length, with the depth interpolated between the two ends. One loop, one
    stamp, radius shared.

    Each half is owned by the atom it touches (`owner_a`, `owner_b`, see
    `_bond_owner`), so an atom drawn with no sphere of its own -- a bonded
    atom in wireframe -- can still be asked whether it is visible."""
    count = int(length / (SAMPLE_SPACING * radius)) + 1
    count = min(max(count, 2), MAX_SAMPLES_PER_BOND)
    ax, ay, az = a
    bx, by, bz = b
    for s in range(count):
        t = s / (count - 1)
        rgb, k = (rgb_a, owner_a) if t <= 0.5 else (rgb_b, owner_b)
        _draw_sphere(zbuf, color, owner,
                     ax + (bx - ax) * t, ay + (by - ay) * t, az + (bz - az) * t,
                     radius, scale, rgb, k)


def _bond_owner(atom_index: int, n: int) -> int:
    """The `owner` value of a bond half touching `atom_index`: offset by the
    atom count so it never reads as the atom's own sphere. A ball whose own
    sprite is hidden must not be labelled because a bond stub poking out from
    behind its occluder happens to lie inside its box."""
    return n + atom_index


def _atom_radius(element: str, representation: str, bonded: bool):
    if representation == "wireframe":
        return None if bonded else WIREFRAME_ATOM_RADIUS
    if representation == "licorice":
        return LICORICE_RADIUS
    if representation == "space-filling":
        return VDW_RADII.get(element, DEFAULT_VDW_RADIUS)
    return BALL_SCALE * VDW_RADII.get(element, DEFAULT_VDW_RADIUS)


def _bonded_mask(pairs, visible: list, n: int) -> list:
    """Which atoms have any bond to a VISIBLE partner -- wireframe keeps a dot
    only where there is none."""
    has_bond = [False] * n
    for i, j in pairs:
        if visible[i] and visible[j]:
            has_bond[i] = True
            has_bond[j] = True
    return has_bond


def _draw_geometry(atoms, view, pairs, has_bond, sx, sy, depth, scale,
                   is_qm, visible, env: _Layer, qm: _Layer):
    coords = np.array([(x, y, z) for _el, x, y, z in atoms], dtype=float)
    representation = view.representation
    for qm_bond in (False, True):
        for i, j in pairs:
            if not (visible[i] and visible[j]):
                continue
            if (is_qm[i] and is_qm[j]) is not qm_bond:
                continue
            layer = qm if qm_bond else env
            if qm_bond:
                radius = BOND_RADIUS.get(representation)
                if representation == "space-filling" or radius is None:
                    continue
                if representation == "ball-and-stick":
                    rgb_a = rgb_b = np.asarray(BALL_BOND_RGB, dtype=float)
                else:
                    rgb_a, rgb_b = _rgb(atoms[i][0]), _rgb(atoms[j][0])
            else:
                radius = ENV_BOND_RADIUS
                rgb_a = _desaturate(_rgb(atoms[i][0]))
                rgb_b = _desaturate(_rgb(atoms[j][0]))
            length = float(np.linalg.norm(coords[j] - coords[i]))
            _draw_bond(layer.zbuf, layer.color, layer.owner,
                       (sx[i], sy[i], depth[i]), (sx[j], sy[j], depth[j]),
                       length, radius, scale, rgb_a, rgb_b,
                       _bond_owner(i, len(atoms)), _bond_owner(j, len(atoms)))
    # The environment never gets a full ball: on space-filling the host would
    # bury the QM region the pane exists to show, and on the other
    # representations its thin desaturated sticks are what marks it as the
    # lower layer. The exception is a host atom bonded to the QM region in
    # ball-and-stick: it gets a small desaturated ball, because a host atom in
    # front of a QM ball was otherwise drawn only as the junction of thin
    # sticks, and the eye read a stick slicing through a ball rather than an
    # atom in front of it.
    if representation == "ball-and-stick" and not all(is_qm):
        boundary = set()
        for i, j in pairs:
            if not (visible[i] and visible[j]):
                continue
            if not is_qm[i] and is_qm[j]:
                boundary.add(i)
            elif is_qm[i] and not is_qm[j]:
                boundary.add(j)
        for i in sorted(boundary):
            element = atoms[i][0]
            _draw_sphere(env.zbuf, env.color, env.owner, sx[i], sy[i], depth[i],
                         BOUNDARY_BALL_SCALE * BALL_SCALE
                         * VDW_RADII.get(element, DEFAULT_VDW_RADIUS),
                         scale, _desaturate(_rgb(element)), i)
    for i, (element, *_rest) in enumerate(atoms):
        if not visible[i] or not is_qm[i]:
            continue
        radius = _atom_radius(element, representation, has_bond[i])
        if radius is None:
            continue
        _draw_sphere(qm.zbuf, qm.color, qm.owner, sx[i], sy[i], depth[i], radius, scale,
                     _rgb(element), i)


def _composite(env: _Layer, qm: _Layer) -> tuple[_Layer, np.ndarray]:
    """Merge the two layers into one frame, and say which layer won each pixel.

    `env_wins` is strict (`>`), so a depth tie goes to the QM layer: the guest
    is the subject the pane exists to show, and a host sample exactly on its
    surface must not paint over it. Where nothing was drawn both depths are
    -inf, the comparison is False, and the background fill takes the pixel.
    """
    env_wins = env.zbuf > qm.zbuf
    frame = _Layer(
        zbuf=np.where(env_wins, env.zbuf, qm.zbuf),
        color=np.where(env_wins[..., None], env.color, qm.color),
        owner=np.where(env_wins, env.owner, qm.owner),
    )
    return frame, env_wins


def _sprite_box(cx, cy, radius, scale, size_px):
    r = max(1, int(round(radius * scale)))
    width, height = size_px
    x0, y0 = int(round(cx)) - r, int(round(cy)) - r
    bx0, by0 = max(x0, 0), max(y0, 0)
    bx1, by1 = min(x0 + 2 * r + 1, width), min(y0 + 2 * r + 1, height)
    if bx0 >= bx1 or by0 >= by1:
        return None
    return bx0, by0, bx1, by1


def _is_visible(owner, atoms, atom_index, view, has_bond, sx, sy, scale, size_px) -> bool:
    """At least one pixel of this atom's own sprite survived occlusion. A
    label must not float over an atom hidden behind another.

    An atom drawn with no sphere -- a bonded atom in wireframe -- has no
    sprite to ask, and this returned False for every one of them, so `l` and
    `d` silently drew nothing in wireframe. There the question goes to the
    bond halves that meet at the atom, inside a box one stick wide around
    its centre: the junction the label sits on."""
    element = atoms[atom_index][0]
    radius = _atom_radius(element, view.representation, has_bond[atom_index])
    target = atom_index
    if radius is None:
        radius = max(BOND_RADIUS.get(view.representation, 0.0), ENV_BOND_RADIUS)
        target = _bond_owner(atom_index, len(atoms))
    box = _sprite_box(sx[atom_index], sy[atom_index], radius, scale, size_px)
    if box is None:
        return False
    bx0, by0, bx1, by1 = box
    return bool((owner[by0:by1, bx0:bx1] == target).any())


def _font(size_px):
    return ImageFont.load_default(size=max(9, round(size_px[1] / 70)))


def _draw_labels(image, atoms, view, sx, sy, scale, is_qm, visible, has_bond, owner, size_px):
    draw = ImageDraw.Draw(image)
    font = _font(size_px)
    for i, (element, *_rest) in enumerate(atoms):
        if not visible[i] or not is_qm[i]:
            continue
        if not _is_visible(owner, atoms, i, view, has_bond, sx, sy, scale, size_px):
            continue
        draw.text((sx[i], sy[i]), str(i), font=font, fill=LABEL_RGB, anchor="mm",
                  stroke_width=1, stroke_fill="black")


def _draw_distances(image, atoms, view, pairs, sx, sy, scale, is_qm, visible, has_bond, owner, size_px):
    draw = ImageDraw.Draw(image)
    font = _font(size_px)
    coords = np.array([(x, y, z) for _el, x, y, z in atoms], dtype=float)
    for i, j in pairs:
        if not (visible[i] and visible[j] and is_qm[i] and is_qm[j]):
            continue
        if not (_is_visible(owner, atoms, i, view, has_bond, sx, sy, scale, size_px)
                or _is_visible(owner, atoms, j, view, has_bond, sx, sy, scale, size_px)):
            continue
        distance = float(np.linalg.norm(coords[j] - coords[i]))
        draw.text(((sx[i] + sx[j]) / 2, (sy[i] + sy[j]) / 2), f"{distance:.2f}",
                  font=font, fill=DISTANCE_RGB, anchor="mm",
                  stroke_width=1, stroke_fill="black")


def _fog_depths(atoms: list, view: View) -> tuple[float, float]:
    """`(near, far)`: the depths fog runs between -- the front and back of
    the framing sphere (`_frame_sphere`), whose centre sits at depth
    `-pan . forward`.

    This was the frame's own nearest and farthest pixel, which moves every
    time the view turns: the same atom at the same depth faded by a
    different amount at each angle, so the whole picture pulsed in
    brightness under rotation and visibly under `o`. The sphere does not
    turn, so an atom's fade now changes only when its depth does. Zoom is
    left out on purpose: it changes the picture's scale, not the depths."""
    _coords, _center, extent = _frame_sphere(atoms, view)
    forward = np.array(camera_forward(view.elev, view.azim))
    mid = -float(np.asarray(view.pan, dtype=float) @ forward)
    return mid + extent, mid - extent


def _apply_fog(layer: _Layer, near: float, far: float, amount: float, power: float) -> None:
    """Fade a layer's foreground toward BACKGROUND between `near` and `far`:
    `near` untouched, `far` at `amount`, the fade growing as `t ** power`.

    The layers are fogged separately, not the composited frame, because the
    host is drawn in uniform desaturated colours the eye already reads as
    "far" everywhere -- and one linear fog for both layers left a host stick
    in front of the guest faded by about a tenth and one behind it by about a
    half, so a host in front did not read as nearer. The host therefore gets
    its own harder amount (HOST_FOG) on a `t ** power` curve
    (HOST_FOG_POWER), which keeps its near end bright; the guest keeps FOG,
    linear. Applied before compositing, so the outline pass and the labels
    see the final colours."""
    foreground = layer.zbuf != -np.inf
    if not foreground.any():
        return
    span = max(near - far, 1e-6)
    t = np.clip((near - layer.zbuf[foreground]) / span, 0.0, 1.0) ** power
    mix = (amount * t)[:, None]
    layer.color[foreground] = (layer.color[foreground] * (1.0 - mix)
                               + np.asarray(BACKGROUND, dtype=float) * mix)


def _apply_halos(zbuf, foreground, color):
    """Darken the far side of every silhouette, as wide as its depth gap.

    A pixel is darkened when something nearer by more than `T_d` lies within
    d pixels. At d = 1 that is the old silhouette: any 4-neighbour is
    background, or is nearer by more than OUTLINE_JUMP. `T_d` then grows
    linearly from OUTLINE_JUMP to HALO_FULL_GAP across `levels`, so the
    halo's width tells the eye HOW FAR in front the near object is -- the
    larger the gap, the more levels clear their threshold and the wider the
    ring. The width scales with the frame's short side (a 329x315 preview
    frame gets 2 levels, not 5), so it stays a cue at pane size instead of
    swallowing the picture.

    This replaced a single 1-px line that gave a host stick 3 angstrom in
    front of a QM ball the same faint edge as two atoms touching, so the
    stick did not read as floating in front and the picture looked like a
    rendering error. Halos are always on -- a readability fix, not an option.

    The `-1e30` sentinel, not `-inf`, keeps the max-propagation finite:
    `inf - inf` is `nan`, and a nan test never darkens.
    """
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


def render(atoms: list, view: View | None = None, *,
           qm_atom_indices: set | None = None,
           size_px: tuple[int, int] = (900, 750),
           bond_list: list[tuple[int, int]] | None = None) -> Image.Image:
    """A shaded, z-buffered frame of `atoms` under `view`. An empty atom list
    gives a background-filled image of `size_px`."""
    view = view or View()
    width, height = max(int(size_px[0]), 1), max(int(size_px[1]), 1)
    env = _Layer.empty(height, width)
    qm = _Layer.empty(height, width)

    if atoms:
        sx, sy, depth, scale = project(atoms, view, (width, height))
        n = len(atoms)
        multilayer = qm_atom_indices is not None and 0 < len(qm_atom_indices) < n
        is_qm = [(i in qm_atom_indices) if multilayer else True for i in range(n)]
        visible = [view.show_hydrogens or element != "H" for element, *_ in atoms]
        pairs = bond_list if bond_list is not None else bonds(atoms)
        has_bond = _bonded_mask(pairs, visible, n)
        _draw_geometry(atoms, view, pairs, has_bond, sx, sy, depth, scale,
                       is_qm, visible, env, qm)
        if view.fog:
            near, far = _fog_depths(atoms, view)
            _apply_fog(qm, near, far, FOG, 1.0)
            _apply_fog(env, near, far, HOST_FOG, HOST_FOG_POWER)

    frame, env_wins = _composite(env, qm)
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
    foreground = frame.zbuf != -np.inf
    frame.color[~foreground] = np.asarray(BACKGROUND, dtype=float)
    if foreground.any():
        _apply_halos(frame.zbuf, foreground, frame.color)
    image = Image.fromarray(np.clip(frame.color, 0, 255).astype(np.uint8), "RGB")

    if atoms and view.show_labels:
        _draw_labels(image, atoms, view, sx, sy, scale, is_qm, visible,
                     has_bond, frame.owner, (width, height))
    if atoms and view.show_distances:
        _draw_distances(image, atoms, view, pairs, sx, sy, scale, is_qm, visible,
                        has_bond, frame.owner, (width, height))
    return image
