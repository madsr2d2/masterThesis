from __future__ import annotations

import io
import itertools

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import proj3d
from PIL import Image

ELEMENT_COLORS = {"H": "#f2f2f2", "O": "#e04040", "C": "#4a4a4a", "N": "#4060e0"}
ELEMENT_RADII = {"H": 0.32, "O": 0.66, "C": 0.70, "N": 0.68}
DEFAULT_ELEMENT_COLOR = "#c060c0"
DEFAULT_ELEMENT_RADIUS = 0.6
BOND_CUTOFF = 1.7  # angstrom, generous single-bond distance cutoff

# Atom indices are 0-based to match ORCA's own convention for referencing
# atoms (internal coordinate definitions, QM region selections, etc.), so a
# number read off the image can be pasted straight into an ORCA input file.
# Comfortably above any zorder Axes3D.computed_zorder could ever assign a
# per-atom scatter collection (zorder_offset + atom count -- see render()).
_OVERLAY_ZORDER = 100_000


def camera_basis(elev_deg: float, azim_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Screen-space right/up unit vectors for `view_init(elev_deg, azim_deg)`.

    Panning applies a step along these at the moment a pan key is pressed,
    then stores the result as a fixed WORLD-space offset (see
    `RotatableGeometryImage._pan_by` in app.py) -- exactly like translating a
    camera before turning it, so a later rotation swings the view around the
    panned-to point rather than re-deriving "screen right" from the new
    angle and drifting the pan with it."""
    elev = np.radians(elev_deg)
    azim = np.radians(azim_deg)
    forward = np.array(
        [np.cos(elev) * np.cos(azim), np.cos(elev) * np.sin(azim), np.sin(elev)]
    )
    world_up = np.array([0.0, 0.0, 1.0])
    right = np.cross(forward, world_up)
    right_norm = np.linalg.norm(right)
    right = right / right_norm if right_norm > 1e-6 else np.array([1.0, 0.0, 0.0])
    up = np.cross(right, forward)
    up = up / np.linalg.norm(up)
    return right, up


def render(
    atoms: list[tuple[str, float, float, float]],
    elev: float = 20,
    azim: float = -60,
    dpi: int = 150,
    show_distances: bool = False,
    zoom: float = 1.0,
    pan: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> Image.Image:
    """Ball-and-stick render of a geometry. Not to scale of any published
    figure -- this feeds a terminal halfcell widget whose real resolution is
    the panel's own cell count, so the source only needs to out-resolve that."""
    fig = plt.figure(figsize=(6, 5), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor("#1e1e1e")
    ax.set_facecolor("#1e1e1e")

    coords = np.array([(x, y, z) for _, x, y, z in atoms])
    bonds = []
    for (i1, a1), (i2, a2) in itertools.combinations(enumerate(atoms), 2):
        _, x1, y1, z1 = a1
        _, x2, y2, z2 = a2
        d = np.linalg.norm([x1 - x2, y1 - y2, z1 - z2])
        if d < BOND_CUTOFF:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color="#cccccc", linewidth=2, zorder=1)
            bonds.append((x1, y1, z1, x2, y2, z2, d))

    for element, x, y, z in atoms:
        color = ELEMENT_COLORS.get(element, DEFAULT_ELEMENT_COLOR)
        radius = ELEMENT_RADII.get(element, DEFAULT_ELEMENT_RADIUS)
        ax.scatter(
            [x], [y], [z], s=(radius * 400) ** 1.1, color=color,
            edgecolor="black", linewidth=0.5, depthshade=True, zorder=2,
        )

    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    span = (np.ptp(coords, axis=0).max() / 2 + 0.5) / zoom
    center = coords.mean(axis=0) + np.array(pan)
    ax.set_xlim(center[0] - span, center[0] + span)
    ax.set_ylim(center[1] - span, center[1] + span)
    ax.set_zlim(center[2] - span, center[2] + span)

    # Axes3D only finalizes its internal 2D data transform (the one
    # `xycoords="data"` below resolves against) DURING a draw call -- calling
    # `ax.get_proj()` beforehand can read a projection that hasn't been
    # reconciled with that transform yet, so the same atom drifts out from
    # under its label depending on view angle. Forcing one draw here, before
    # projecting anything, makes get_proj() and the annotations agree with
    # what savefig's own draw will actually place the spheres at.
    fig.canvas.draw()

    # Atom index labels are placed with ax.annotate at each atom's PROJECTED
    # 2D position, computed once here with the view/limits already final --
    # NOT ax.text at the 3D point. mplot3d re-sorts every 3D artist (bonds,
    # atom spheres, 3D text) by camera depth on each draw, and a label
    # sitting exactly at its atom's own point ties that sort, so whichever
    # one wins the pixel flipped with the view angle -- that was the "labels
    # disappear" bug. A plain 2D annotation isn't part of that resort, so it
    # always draws on top of every 3D artist already added, regardless of
    # view. Must come after set_xlim/ylim/zlim and view_init above: both
    # feed ax.get_proj(), which the projection below depends on.
    #
    # Axes3D.computed_zorder (on by default) ALSO reassigns the zorder of
    # every Collection -- and each atom's ax.scatter call is its own
    # Path3DCollection -- to zorder_offset..zorder_offset+n_atoms, sorted by
    # camera depth, EVERY draw. A fixed zorder on the label (10, originally)
    # only beats that for small molecules; past ~10 atoms some sphere's
    # reassigned zorder exceeds it, and WHICH one does depends on the
    # depth order, i.e. the view angle -- so a label would flicker in and
    # out from behind its own atom as you rotated. `_OVERLAY_ZORDER` is
    # picked far above anything that reassignment could ever produce.
    for index, (_, x, y, z) in enumerate(atoms):
        lx, ly, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
        ax.annotate(
            str(index), xy=(lx, ly), xycoords="data", ha="center", va="center",
            color="white", fontsize=13, fontweight="bold", zorder=_OVERLAY_ZORDER,
            path_effects=[pe.withStroke(linewidth=3, foreground="black")],
            bbox=dict(boxstyle="round,pad=0.15,rounding_size=0.4", fc="black", ec="none", alpha=0.65),
        )

    if show_distances:
        for x1, y1, z1, x2, y2, z2, d in bonds:
            mx, my, mz = (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2
            lx, ly, _ = proj3d.proj_transform(mx, my, mz, ax.get_proj())
            ax.annotate(
                f"{d:.2f}", xy=(lx, ly), xycoords="data", ha="center", va="center",
                color="#ffe066", fontsize=10, fontweight="bold", zorder=_OVERLAY_ZORDER - 1,
                path_effects=[pe.withStroke(linewidth=3, foreground="black")],
                bbox=dict(boxstyle="round,pad=0.15,rounding_size=0.4", fc="black", ec="none", alpha=0.65),
            )

    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf)
