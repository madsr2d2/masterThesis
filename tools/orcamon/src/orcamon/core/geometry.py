"""Molecule-drawing facts both geometry renderers share: which atoms are
bonded, what colour an element is, how big its ball is, and where the camera
points.

The pixel renderer (`tui/raster.py`) and the text one (`tui/geometry_text.py`)
must agree on all of it, or rotating the same molecule in the two modes would
show different bonds, colours or angles. Pure Python so the text renderer
needs no numpy -- anything the text mode uses (camera, colours, radii, the
`View`) lives here; the z-buffer and the sprites do not.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

# Jmol/CPK colours, so a structure drawn here reads the same as in a real
# viewer. Carbon is mid grey rather than the near-black it used to be: on the
# fixed dark background a shaded sphere needs its own tone to read through,
# and the text renderer lifts carbon further still (TEXT_ELEMENT_COLORS).
ELEMENT_COLORS = {
    "H": "#FFFFFF", "C": "#909090", "N": "#3050F8", "O": "#FF0D0D",
    "F": "#90E050", "P": "#FF8000", "S": "#FFFF30", "Cl": "#1FF01F",
    "Br": "#A62929", "I": "#940094", "B": "#FFB5B5", "Si": "#F0C8A0",
    "Na": "#AB5CF2", "K": "#8F40D4", "Mg": "#8AFF00", "Ca": "#3DFF00",
    "Fe": "#E06633", "Cu": "#C88033", "Zn": "#7D80B0",
}
DEFAULT_ELEMENT_COLOR = "#c060c0"
BOND_CUTOFF = 1.7  # angstrom, generous single-bond distance cutoff

# Bondi van der Waals radii (angstrom), the space-filling and ball scale.
VDW_RADII = {
    "H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "B": 1.92, "Si": 2.10,
    "Na": 2.27, "K": 2.75, "Mg": 1.73, "Ca": 2.31, "Fe": 2.00, "Cu": 1.40,
    "Zn": 1.39,
}
DEFAULT_VDW_RADIUS = 1.70

REPRESENTATIONS = ("ball-and-stick", "licorice", "space-filling", "wireframe")
# Text cells cannot shade or occlude, so only the two the braille renderer can
# actually tell apart are offered there.
TEXT_REPRESENTATIONS = ("ball-and-stick", "wireframe")


@dataclass
class View:
    """One view description, shared by both renderers.

    Kept as plain fields (no methods) so the text widget, the two pixel
    widgets and the `snapshot` command all describe a view the same way, and
    so a new option is added in one place.
    """

    elev: float = 20.0
    azim: float = -60.0
    zoom: float = 1.0
    pan: tuple[float, float, float] = (0.0, 0.0, 0.0)
    representation: str = "ball-and-stick"
    fog: bool = True
    show_labels: bool = True
    show_distances: bool = False
    show_hydrogens: bool = True


def camera_forward(elev_deg: float, azim_deg: float) -> tuple[float, float, float]:
    """The unit vector from the molecule's centre TO the eye, for
    `(elev, azim)` in degrees (matplotlib's `view_init` convention):
    `f = (cos e cos a, cos e sin a, sin e)` for elevation e, azimuth a.
    A larger `p . f` is nearer the viewer."""
    elev, azim = math.radians(elev_deg), math.radians(azim_deg)
    return (math.cos(elev) * math.cos(azim), math.cos(elev) * math.sin(azim), math.sin(elev))


def camera_basis(elev_deg: float, azim_deg: float) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Screen-space (right, up) unit vectors for matplotlib's
    `view_init(elev_deg, azim_deg)`.

    `right = (-sin a, cos a, 0)` and `up = cross(forward, right)`, so
    `cross(right, up) == forward`: the screen's right-handed normal points AT
    the viewer. This returned `cross(forward, z_hat)`, which is `-right`,
    until 2026-09-30 -- and the text renderer uses it as screen +x, so text
    mode drew the molecule's mirror image. The six-pair check in
    `test_monitor.test_the_view_is_not_mirrored` pins the sign.

    `right` is NOT `normalize(cross(z_hat, forward))`, though it equals it
    for |elev| < 90. That cross product is `cos(e) * (-sin a, cos a, 0)`, so
    normalising it flips the sign the moment the elevation passes a pole --
    and the elevation keys run through 360 degrees, so holding one tumbled the
    molecule smoothly to 89 degrees and then snapped it upside down at 91. It
    also vanished AT the pole, where the old fallback `(1, 0, 0)` was a jump
    of its own for every azimuth but -90. The un-normalised direction is
    continuous everywhere, and is already a unit vector perpendicular to
    `forward`. `test_the_camera_tumbles_over_the_pole` holds it.

    Panning applies a step along these at the moment a pan key is pressed,
    then stores the result as a fixed WORLD-space offset (see
    `RotatableGeometryImage._pan_by`) -- exactly like translating a camera
    before turning it, so a later rotation swings the view around the
    panned-to point rather than re-deriving "screen right" from the new angle
    and drifting the pan with it."""
    forward = camera_forward(elev_deg, azim_deg)
    azim = math.radians(azim_deg)
    right = (-math.sin(azim), math.cos(azim), 0.0)
    up = _cross(forward, right)
    return right, up


def _cross(a, b) -> tuple[float, float, float]:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def bonds(atoms: list) -> list[tuple[int, int]]:
    """Index pairs closer than BOND_CUTOFF, H-H excluded, by grid binning.

    Both renderers and `snapshot` must bond the same pairs, so this is the one
    implementation; the search bins atoms into BOND_CUTOFF-sized cells, making
    it linear in the atom count rather than quadratic (300 atoms is ~2,000
    distance checks, not 45,000)."""
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


def principal_axes(coords: list) -> list[tuple[tuple[float, float, float], float]]:
    """(unit axis, variance) for the three principal axes of `coords`
    (unweighted), largest variance first.

    A 3x3 symmetric eigenproblem, solved by cyclic Jacobi rotations: the core
    may not import numpy. The covariance matrix is `A`, the rotations are
    accumulated into `V`, and the eigenvectors are V's COLUMNS. Each sweep
    visits the three off-diagonal pairs in turn and zeroes each with
    `theta = 0.5 * atan2(2 A[p][q], A[q][q] - A[p][p])`, applied as
    `A <- J^T A J` and `V <- V J`. `test_principal_axes_face_a_plane_on`
    checks `A v = lambda v` for the result, which pins the sign convention."""
    n = len(coords)
    if n == 0:
        return []
    mean = [sum(p[i] for p in coords) / n for i in range(3)]
    a = [[0.0, 0.0, 0.0] for _ in range(3)]
    for p in coords:
        d = (p[0] - mean[0], p[1] - mean[1], p[2] - mean[2])
        for i in range(3):
            for j in range(3):
                a[i][j] += d[i] * d[j] / n
    v = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    for _sweep in range(50):
        if abs(a[0][1]) + abs(a[0][2]) + abs(a[1][2]) < 1e-12:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            if abs(a[p][q]) < 1e-12:
                continue
            theta = 0.5 * math.atan2(2 * a[p][q], a[q][q] - a[p][p])
            c, s = math.cos(theta), math.sin(theta)
            for k in range(3):
                if k != p and k != q:
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = a[p][k] = c * akp - s * akq
                    a[k][q] = a[q][k] = s * akp + c * akq
            app, aqq, apq = a[p][p], a[q][q], a[p][q]
            a[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq
            a[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq
            a[p][q] = a[q][p] = 0.0
            for k in range(3):
                vkp, vkq = v[k][p], v[k][q]
                v[k][p] = c * vkp - s * vkq
                v[k][q] = s * vkp + c * vkq
    axes = [(tuple(v[k][i] for k in range(3)), a[i][i]) for i in range(3)]
    axes.sort(key=lambda axis: axis[1], reverse=True)
    return axes


def view_along(axis) -> tuple[float, float]:
    """(elev, azim) in degrees for an eye on +axis:
    `elev = degrees(asin(z))`, `azim = degrees(atan2(y, x))` -- the inverse
    of `camera_forward`."""
    x, y, z = axis
    return (math.degrees(math.asin(max(-1.0, min(1.0, z)))), math.degrees(math.atan2(y, x)))


@dataclass
class FileGeometry:
    """A geometry read from a FILE rather than printed in the log, and where
    from. It stands in for a GeometryPoint wherever the panes and `geom` take
    one, with no cycle and no energy, so nothing presents it as a printed
    optimization step: `source` is what the pane title and the XYZ comment
    say it is."""

    atoms: list
    source: str
    scan_step: None = None
    cycle: int = 0
    energy: None = None
    energy_label: None = None
    key: None = None


def read_xyz(path: Path) -> list | None:
    """The first frame of an XYZ file as (element, x, y, z) tuples, or None
    when the file is missing or is not XYZ. Tolerant of what ORCA writes:
    an element may carry a suffix (`C1`, `H:`), and a trajectory's later
    frames are ignored."""
    try:
        with open(path, errors="replace") as f:
            n = int(f.readline().split()[0])
            f.readline()
            atoms = []
            for _ in range(n):
                parts = f.readline().split()
                element = "".join(ch for ch in parts[0] if ch.isalpha())
                atoms.append((element.capitalize(), float(parts[1]), float(parts[2]), float(parts[3])))
    except (OSError, ValueError, IndexError):
        return None
    return atoms or None
