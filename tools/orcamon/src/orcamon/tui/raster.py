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
BOND_RADIUS = {"ball-and-stick": 0.10, "licorice": 0.16, "wireframe": 0.04}
LICORICE_RADIUS = 0.16                # licorice atom cap
WIREFRAME_ATOM_RADIUS = 0.12          # the dot an unbonded atom keeps
ENV_BOND_RADIUS = 0.06                # the environment layer's thin sticks
ENV_DESATURATE = 0.5                  # mix environment colours halfway to grey
BALL_BOND_RGB = (184, 184, 184)       # ball-and-stick bonds: light grey
SAMPLE_SPACING = 1.2                  # bond sample spacing, in bond radii
MAX_SAMPLES_PER_BOND = 64
LABEL_RGB = (255, 255, 255)
DISTANCE_RGB = (255, 224, 102)


def _rgb(element: str) -> np.ndarray:
    colour = ELEMENT_COLORS.get(element, DEFAULT_ELEMENT_COLOR)
    return np.array([int(colour[i:i + 2], 16) for i in (1, 3, 5)], dtype=float)


def _desaturate(rgb: np.ndarray) -> np.ndarray:
    grey = rgb.mean()
    return rgb * (1.0 - ENV_DESATURATE) + grey * ENV_DESATURATE


def _framing_pad(atoms: list, view: View) -> float:
    """The largest radius the CURRENT representation draws, so space-filling
    never clips at the frame edge; the smaller representations need only half
    an angstrom of slack."""
    if view.representation == "space-filling":
        return max((VDW_RADII.get(el, DEFAULT_VDW_RADIUS) for el, *_ in atoms),
                   default=DEFAULT_VDW_RADIUS)
    return 0.5


def project(atoms: list, view: View, size_px: tuple[int, int]):
    """Screen x, screen y, depth (angstrom, larger is NEARER) and pixels per
    angstrom for every atom. Public because the tests read it.

    The frame is the CIRCUMSCRIBED sphere -- centre `mean + pan`, radius the
    farthest atom -- so no rotation can push an atom out of the pane, and
    neither changes with rotation."""
    coords = np.array([(x, y, z) for _el, x, y, z in atoms], dtype=float)
    right, up = camera_basis(view.elev, view.azim)
    forward = camera_forward(view.elev, view.azim)
    right, up, forward = np.array(right), np.array(up), np.array(forward)
    mean = coords.mean(axis=0)
    radius = float(np.linalg.norm(coords - mean, axis=1).max())
    center = mean + np.asarray(view.pan, dtype=float)
    span = (radius + _framing_pad(atoms, view)) / view.zoom
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


def _draw_bond(zbuf, color, owner, a, b, length, radius, scale, rgb_a, rgb_b):
    """A bond is a row of spheres of `radius` from screen point `a` to `b`
    (`(x, y, depth)`), spaced SAMPLE_SPACING * radius apart along its true 3D
    length, with the depth interpolated between the two ends. One loop, one
    stamp, radius shared."""
    count = int(length / (SAMPLE_SPACING * radius)) + 1
    count = min(max(count, 2), MAX_SAMPLES_PER_BOND)
    ax, ay, az = a
    bx, by, bz = b
    for s in range(count):
        t = s / (count - 1)
        rgb = rgb_a if t <= 0.5 else rgb_b
        _draw_sphere(zbuf, color, owner,
                     ax + (bx - ax) * t, ay + (by - ay) * t, az + (bz - az) * t,
                     radius, scale, rgb, -1)


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
                   is_qm, visible, zbuf, color, owner):
    coords = np.array([(x, y, z) for _el, x, y, z in atoms], dtype=float)
    representation = view.representation
    # Environment first, so when a QM and an environment sample tie on depth
    # the QM region wins the pixel.
    for qm_bond in (False, True):
        for i, j in pairs:
            if not (visible[i] and visible[j]):
                continue
            if (is_qm[i] and is_qm[j]) is not qm_bond:
                continue
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
            _draw_bond(zbuf, color, owner,
                       (sx[i], sy[i], depth[i]), (sx[j], sy[j], depth[j]),
                       length, radius, scale, rgb_a, rgb_b)
    # The environment never gets a ball: on space-filling the host would bury
    # the QM region the pane exists to show, and on the other representations
    # its thin desaturated sticks are what marks it as the lower layer.
    for i, (element, *_rest) in enumerate(atoms):
        if not visible[i] or not is_qm[i]:
            continue
        radius = _atom_radius(element, representation, has_bond[i])
        if radius is None:
            continue
        _draw_sphere(zbuf, color, owner, sx[i], sy[i], depth[i], radius, scale,
                     _rgb(element), i)


def _sprite_box(cx, cy, radius, scale, size_px):
    r = max(1, int(round(radius * scale)))
    width, height = size_px
    x0, y0 = int(round(cx)) - r, int(round(cy)) - r
    bx0, by0 = max(x0, 0), max(y0, 0)
    bx1, by1 = min(x0 + 2 * r + 1, width), min(y0 + 2 * r + 1, height)
    if bx0 >= bx1 or by0 >= by1:
        return None
    return bx0, by0, bx1, by1


def _is_visible(owner, atom_index, cx, cy, radius, scale, size_px) -> bool:
    """At least one pixel of this atom's own sprite survived occlusion. A
    label must not float over an atom hidden behind another."""
    box = _sprite_box(cx, cy, radius, scale, size_px)
    if box is None:
        return False
    bx0, by0, bx1, by1 = box
    return bool((owner[by0:by1, bx0:bx1] == atom_index).any())


def _font(size_px):
    return ImageFont.load_default(size=max(9, round(size_px[1] / 70)))


def _draw_labels(image, atoms, view, sx, sy, scale, is_qm, visible, has_bond, owner, size_px):
    draw = ImageDraw.Draw(image)
    font = _font(size_px)
    for i, (element, *_rest) in enumerate(atoms):
        if not visible[i] or not is_qm[i]:
            continue
        radius = _atom_radius(element, view.representation, has_bond[i])
        if radius is None or not _is_visible(owner, i, sx[i], sy[i], radius, scale, size_px):
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
        radius_i = _atom_radius(atoms[i][0], view.representation, has_bond[i])
        radius_j = _atom_radius(atoms[j][0], view.representation, has_bond[j])
        shown = (radius_i is not None
                 and _is_visible(owner, i, sx[i], sy[i], radius_i, scale, size_px)) or (
                 radius_j is not None
                 and _is_visible(owner, j, sx[j], sy[j], radius_j, scale, size_px))
        if not shown:
            continue
        distance = float(np.linalg.norm(coords[j] - coords[i]))
        draw.text(((sx[i] + sx[j]) / 2, (sy[i] + sy[j]) / 2), f"{distance:.2f}",
                  font=font, fill=DISTANCE_RGB, anchor="mm",
                  stroke_width=1, stroke_fill="black")


def render(atoms: list, view: View | None = None, *,
           qm_atom_indices: set | None = None,
           size_px: tuple[int, int] = (900, 750),
           bond_list: list[tuple[int, int]] | None = None) -> Image.Image:
    """A shaded, z-buffered frame of `atoms` under `view`. An empty atom list
    gives a background-filled image of `size_px`."""
    view = view or View()
    width, height = max(int(size_px[0]), 1), max(int(size_px[1]), 1)
    zbuf = np.full((height, width), -np.inf)
    color = np.zeros((height, width, 3), dtype=float)
    owner = np.full((height, width), -1, dtype=int)

    if atoms:
        sx, sy, depth, scale = project(atoms, view, (width, height))
        n = len(atoms)
        multilayer = qm_atom_indices is not None and 0 < len(qm_atom_indices) < n
        is_qm = [(i in qm_atom_indices) if multilayer else True for i in range(n)]
        visible = [view.show_hydrogens or element != "H" for element, *_ in atoms]
        pairs = bond_list if bond_list is not None else bonds(atoms)
        has_bond = _bonded_mask(pairs, visible, n)
        _draw_geometry(atoms, view, pairs, has_bond, sx, sy, depth, scale,
                       is_qm, visible, zbuf, color, owner)

    foreground = zbuf != -np.inf
    color[~foreground] = np.asarray(BACKGROUND, dtype=float)
    image = Image.fromarray(np.clip(color, 0, 255).astype(np.uint8), "RGB")

    if atoms and view.show_labels:
        _draw_labels(image, atoms, view, sx, sy, scale, is_qm, visible,
                     has_bond, owner, (width, height))
    if atoms and view.show_distances:
        _draw_distances(image, atoms, view, pairs, sx, sy, scale, is_qm, visible,
                        has_bond, owner, (width, height))
    return image
