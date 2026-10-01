"""Molecule-drawing facts both geometry renderers share: which atoms are
bonded, what colour an element is, and where the camera points.

The pixel renderer (`tui/geometry_render.py`, matplotlib) and the text one
(`tui/geometry_text.py`, braille) must agree on all three, or rotating the
same molecule in the two modes would show different bonds from different
angles. Pure Python so the text renderer needs no numpy."""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

ELEMENT_COLORS = {"H": "#f2f2f2", "O": "#e04040", "C": "#4a4a4a", "N": "#4060e0"}
DEFAULT_ELEMENT_COLOR = "#c060c0"
BOND_CUTOFF = 1.7  # angstrom, generous single-bond distance cutoff


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

    `right = normalize(cross(z_hat, forward))` and `up = cross(forward,
    right)`, so `cross(right, up) == forward`: the screen's right-handed
    normal points AT the viewer. This returned `cross(forward, z_hat)`,
    which is `-right`, until 2026-09-30 -- and the text renderer uses it as
    screen +x, so text mode drew the molecule's mirror image. The six-pair
    check in `test_monitor.test_the_view_is_not_mirrored` pins the sign.

    Panning applies a step along these at the moment a pan key is pressed,
    then stores the result as a fixed WORLD-space offset (see
    `RotatableGeometryImage._pan_by`) -- exactly like translating a camera
    before turning it, so a later rotation swings the view around the
    panned-to point rather than re-deriving "screen right" from the new angle
    and drifting the pan with it."""
    forward = camera_forward(elev_deg, azim_deg)
    right = _cross((0.0, 0.0, 1.0), forward)
    norm = math.sqrt(sum(c * c for c in right))
    right = tuple(c / norm for c in right) if norm > 1e-6 else (1.0, 0.0, 0.0)
    up = _cross(forward, right)
    return right, up


def _cross(a, b) -> tuple[float, float, float]:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


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
