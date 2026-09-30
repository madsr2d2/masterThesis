"""The molecule as text: bonds in braille, atoms as their element symbols.

For every terminal the pixel paths cannot reach -- tmux, screen, mosh, a
plain xterm, an ssh session through anything that strips escapes -- this is
the geometry pane. It is ordinary cell updates, a few KB a frame, so it
rotates as fast as the link carries text.

Same camera as the pixel renderer (`core.geometry.camera_basis`, the one
matplotlib's `view_init` uses) and the same bond rule, so switching mode
does not change what is bonded or which way is up. Orthographic, framed on
the circumscribed sphere like the pixel frame, so no rotation pushes an atom
out of the pane.

Pure Python plus Rich: no numpy. The bond search bins atoms into a grid of
BOND_CUTOFF-sized cells, so it is linear in the atom count rather than
quadratic -- 300 atoms is ~2,000 distance checks, not 45,000.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from rich.style import Style
from rich.text import Text

from ..core.geometry import BOND_CUTOFF, DEFAULT_ELEMENT_COLOR, ELEMENT_COLORS, camera_basis

# The pixel renderer draws carbon dark grey on its own dark background, where
# it reads through the sphere's shading and black edge. A bare letter in that
# colour all but vanishes on a dark terminal, so text mode lifts it.
TEXT_ELEMENT_COLORS = {**ELEMENT_COLORS, "C": "#a8a8a8"}
QM_BOND_STYLE = Style(color="#cccccc")
ENV_BOND_STYLE = Style(color="#6a6a6a")
SHOW_HYDROGENS_UP_TO = 60  # atoms; above this, hydrogens start hidden
ATOM_DOT = "●"  # an atom with its label off

# Braille: each cell is a 2x4 dot grid, dot (x, y) -> bit.
_BRAILLE_BIT = {
    (0, 0): 0x01, (0, 1): 0x02, (0, 2): 0x04, (1, 0): 0x08,
    (1, 1): 0x10, (1, 2): 0x20, (0, 3): 0x40, (1, 3): 0x80,
}


@dataclass
class View:
    elev: float = 20.0
    azim: float = -60.0
    zoom: float = 1.0
    pan: tuple[float, float, float] = (0.0, 0.0, 0.0)
    show_hydrogens: bool = True
    # Off: each atom is a dot in its element colour instead of its symbol,
    # so a crowded region reads as structure rather than as letters.
    show_labels: bool = True


def bonds(atoms: list) -> list[tuple[int, int]]:
    """Index pairs closer than BOND_CUTOFF, H-H excluded, by grid binning."""
    cell = BOND_CUTOFF
    grid: dict[tuple[int, int, int], list[int]] = {}
    for i, (_el, x, y, z) in enumerate(atoms):
        grid.setdefault((int(math.floor(x / cell)), int(math.floor(y / cell)), int(math.floor(z / cell))), []).append(i)
    cutoff2 = cell * cell
    found = []
    for (cx, cy, cz), members in grid.items():
        near = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    near.extend(grid.get((cx + dx, cy + dy, cz + dz), ()))
        for i in members:
            el_i, xi, yi, zi = atoms[i]
            for j in near:
                if j <= i:
                    continue
                el_j, xj, yj, zj = atoms[j]
                if el_i == "H" and el_j == "H":
                    continue
                if (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2 < cutoff2:
                    found.append((i, j))
    return found


def render(
    atoms: list,
    width: int,
    height: int,
    view: View | None = None,
    qm_atom_indices: set | None = None,
    bond_list: list[tuple[int, int]] | None = None,
) -> Text:
    """`atoms` drawn into `width` x `height` cells. `bond_list` may be passed
    in when the caller caches it -- the bonds of a geometry do not change as
    it is rotated."""
    view = view or View()
    if not atoms or width <= 0 or height <= 0:
        return Text("\n".join(" " * max(0, width) for _ in range(max(0, height))))
    n = len(atoms)
    multilayer = qm_atom_indices is not None and 0 < len(qm_atom_indices) < n
    is_qm = [(i in qm_atom_indices) if multilayer else True for i in range(n)]
    visible = [view.show_hydrogens or el != "H" for el, *_ in atoms]

    right, up = camera_basis(view.elev, view.azim)
    forward = (up[1] * right[2] - up[2] * right[1],
               up[2] * right[0] - up[0] * right[2],
               up[0] * right[1] - up[1] * right[0])
    cx = sum(a[1] for a in atoms) / n + view.pan[0]
    cy = sum(a[2] for a in atoms) / n + view.pan[1]
    cz = sum(a[3] for a in atoms) / n + view.pan[2]
    radius = max(math.dist((a[1], a[2], a[3]), (cx - view.pan[0], cy - view.pan[1], cz - view.pan[2]))
                 for a in atoms)
    span = (radius + 0.5) / view.zoom

    # Dots are square when a cell is twice as tall as it is wide, which is
    # the usual terminal font: 2 dots across a cell, 4 down.
    dots_w, dots_h = 2 * width, 4 * height
    scale = (min(dots_w, dots_h) - 1) / (2 * span)
    ox, oy = dots_w / 2, dots_h / 2

    projected = []
    for _el, x, y, z in atoms:
        px, py, pz = x - cx, y - cy, z - cz
        sx = px * right[0] + py * right[1] + pz * right[2]
        sy = px * up[0] + py * up[1] + pz * up[2]
        depth = px * forward[0] + py * forward[1] + pz * forward[2]
        projected.append((ox + sx * scale, oy - sy * scale, depth))

    cells = [0] * (width * height)
    bond_style: list[Style | None] = [None] * (width * height)

    def plot(dx: int, dy: int, style: Style) -> None:
        if 0 <= dx < dots_w and 0 <= dy < dots_h:
            k = (dy // 4) * width + (dx // 2)
            cells[k] |= _BRAILLE_BIT[(dx % 2, dy % 4)]
            if bond_style[k] is not QM_BOND_STYLE:
                bond_style[k] = style

    for i, j in (bond_list if bond_list is not None else bonds(atoms)):
        if not (visible[i] and visible[j]):
            continue
        style = QM_BOND_STYLE if (is_qm[i] and is_qm[j]) else ENV_BOND_STYLE
        _line(int(projected[i][0]), int(projected[i][1]), int(projected[j][0]), int(projected[j][1]), style, plot)

    # Atoms over bonds, nearest last so it wins a shared cell.
    glyphs: dict[int, tuple[str, Style]] = {}
    for i in sorted(range(n), key=lambda k: projected[k][2]):
        if not visible[i]:
            continue
        el = atoms[i][0]
        col, row = int(projected[i][0]) // 2, int(projected[i][1]) // 4
        if not (0 <= row < height):
            continue
        color = TEXT_ELEMENT_COLORS.get(el, DEFAULT_ELEMENT_COLOR)
        style = Style(color=color, bold=is_qm[i], dim=not is_qm[i])
        for off, ch in enumerate(el[:2] if view.show_labels else ATOM_DOT):
            if 0 <= col + off < width:
                glyphs[row * width + col + off] = (ch, style)

    text = Text(no_wrap=True, overflow="crop")
    for row in range(height):
        run_chars: list[str] = []
        run_style: Style | None = None
        for col in range(width):
            k = row * width + col
            if k in glyphs:
                ch, style = glyphs[k]
            elif cells[k]:
                ch, style = chr(0x2800 + cells[k]), bond_style[k]
            else:
                ch, style = " ", None
            if style is not run_style and run_chars:
                text.append("".join(run_chars), run_style)
                run_chars = []
            run_style = style
            run_chars.append(ch)
        if run_chars:
            text.append("".join(run_chars), run_style)
        if row != height - 1:
            text.append("\n")
    return text


def _line(x0: int, y0: int, x1: int, y1: int, style, plot) -> None:
    """Bresenham, all octants."""
    dx, dy = abs(x1 - x0), -abs(y1 - y0)
    sx, sy = (1 if x0 < x1 else -1), (1 if y0 < y1 else -1)
    err = dx + dy
    while True:
        plot(x0, y0, style)
        if x0 == x1 and y0 == y1:
            return
        e2 = 2 * err
        if e2 >= dy:
            err += dy
            x0 += sx
        if e2 <= dx:
            err += dx
            y0 += sy
