from __future__ import annotations

import io

import matplotlib.patheffects as pe
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from PIL import Image

from ..core.geometry import BOND_CUTOFF, DEFAULT_ELEMENT_COLOR, ELEMENT_COLORS

ELEMENT_RADII = {"H": 0.32, "O": 0.66, "C": 0.70, "N": 0.68}
DEFAULT_ELEMENT_RADIUS = 0.6

# How much of the frame the data cube fills. mplot3d sizes the cube so its
# DIAGONAL fits at any view angle, so sqrt(3) is exactly the factor that wastes
# -- and exactly the safe ceiling, PROVIDED the axis limits come from the
# circumscribed sphere rather than the bounding box (see render()).
BOX_ZOOM = 3 ** 0.5

# Atom indices are 0-based to match ORCA's own convention for referencing
# atoms (internal coordinate definitions, QM region selections, etc.), so a
# number read off the image can be pasted straight into an ORCA input file.
# Comfortably above any zorder Axes3D.computed_zorder could ever assign a
# per-atom scatter collection (zorder_offset + atom count -- see render()).
_OVERLAY_ZORDER = 100_000


def _new_figure(figsize: tuple[float, float], dpi: float) -> Figure:
    """A figure on its own Agg canvas, outside pyplot. Renders run in worker
    threads -- two at once when a rotation overlaps the periodic refresh --
    and pyplot keeps a process-global registry of open figures that is not
    thread-safe; a bare Figure touches no shared state."""
    fig = Figure(figsize=figsize, dpi=dpi)
    FigureCanvasAgg(fig)
    return fig


def _to_image(fig: Figure) -> Image.Image:
    """The figure's pixels, straight off the canvas.

    This used to `savefig` a PNG and `Image.open` it -- encoding and decoding
    a frame that the caller then encoded AGAIN for transmission. The pixels
    are identical (checked: max difference 0 on a 940x900 frame), and a
    full-pane render went from 56 ms to 38 ms."""
    fig.canvas.draw()
    size = fig.canvas.get_width_height()
    return Image.frombuffer("RGBA", size, fig.canvas.buffer_rgba(), "raw", "RGBA", 0, 1).convert("RGB")


def frame_png(image: Image.Image) -> bytes:
    """A rendered frame as the PNG sent to the terminal: 256-colour palette,
    undithered.

    A frame is a flat background, flat element colours and their antialiased
    edges -- about 1,700 distinct colours on a 140-atom QM/XTB frame -- and a
    palette PNG carries it in 35 KB of base64 against 129 KB as truecolour,
    at half the encode time. Mean error is ~1 level in 255 and side by side
    the two are indistinguishable. On an ssh-relayed pane, bytes ARE the
    latency (see herdr_graphics), so this is the frame rate of a rotation."""
    buf = io.BytesIO()
    image.convert("RGB").quantize(
        colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE
    ).save(buf, format="PNG")
    return buf.getvalue()


def render(
    atoms: list[tuple[str, float, float, float]],
    elev: float = 20,
    azim: float = -60,
    dpi: int = 150,
    show_distances: bool = False,
    show_labels: bool = True,
    zoom: float = 1.0,
    pan: tuple[float, float, float] = (0.0, 0.0, 0.0),
    qm_atom_indices: set[int] | None = None,
    size_px: tuple[int, int] | None = None,
) -> Image.Image:
    """Ball-and-stick render of a geometry. Not to scale of any published
    figure -- this feeds a terminal halfcell widget whose real resolution is
    the panel's own cell count, so the source only needs to out-resolve that.

    `qm_atom_indices` (0-based, matching ORCA's own numbering) marks a
    multilayer (QM/MM, QM/XTB, ONIOM) job's high-level layer: those atoms get
    the full ball-and-stick treatment (sphere, index label, bond distances);
    every other atom -- the lower-level layer -- gets drawn wireframe-only
    (thin bond lines, no sphere, no label), so a 100+-atom MM environment
    doesn't drown out the QM region the panel exists to show. None (an
    ordinary, non-multilayer job) treats every atom as the QM layer, i.e.
    the plain ball-and-stick rendering this always did.

    `size_px` renders straight to a target pixel size instead of the default
    900x750. A caller that is going to downscale the result anyway -- the
    herdr backend's low-latency preview frame most of all, which lands around
    100px on its side -- should ask for the size it wants rather than pay to
    draw pixels it will throw away."""
    figsize = (6.0, 5.0)
    if size_px is not None:
        # Scale by DPI, keeping the figure's size in INCHES fixed. Fonts,
        # line widths and marker sizes are all in points, i.e. inches, so this
        # shrinks them along with the frame -- exactly as downscaling the big
        # render did. Deriving figsize from a fixed dpi instead would hold
        # them at their absolute size, and a 118px preview would be almost
        # entirely atom-number label.
        width_px = max(int(size_px[0]), 1)
        height_px = max(int(size_px[1]), 1)
        figsize = (6.0, 6.0 * height_px / width_px)
        dpi = width_px / 6.0

    if not atoms:
        # Both call sites guard this, but the signature is all-defaults and
        # the framing maths below (ptp/mean over the coordinates) raises on an
        # empty array rather than producing an empty picture.
        fig = _new_figure(figsize, dpi)
        fig.patch.set_facecolor("#1e1e1e")
        return _to_image(fig)

    coords = np.array([(x, y, z) for _, x, y, z in atoms], dtype=float)
    elements = np.array([element for element, *_ in atoms])
    n_atoms = len(atoms)

    is_multilayer = qm_atom_indices is not None and 0 < len(qm_atom_indices) < n_atoms
    if is_multilayer:
        qm_mask = np.zeros(n_atoms, dtype=bool)
        inside = [i for i in qm_atom_indices if 0 <= i < n_atoms]
        qm_mask[inside] = True
    else:
        qm_mask = np.ones(n_atoms, dtype=bool)

    fig = _new_figure(figsize, dpi)
    ax = fig.add_subplot(111, projection="3d")
    # Orthographic, not the mplot3d default perspective: parallel bonds stay
    # parallel and a farther atom doesn't shrink relative to a nearer one,
    # which makes rotation easier to read as a rigid tumble rather than a
    # constantly-reshaping shape.
    ax.set_proj_type("ortho")
    fig.patch.set_facecolor("#1e1e1e")
    ax.set_facecolor("#1e1e1e")
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)

    # Framing. Two things buy back the ~90% of the canvas this used to spend
    # on empty background -- which over ssh is 90% of the bytes:
    #
    #   * `set_position` takes the axes to the whole figure. The default
    #     subplot margins exist to leave room for the tick labels and title an
    #     axis-off 3D plot does not have.
    #   * `set_box_aspect(zoom=)` fills the frame with the data cube. mplot3d
    #     sizes the cube so its DIAGONAL fits at every view angle, which
    #     wastes a factor of sqrt(3) on a shape that never reaches the
    #     corners.
    #
    # BOX_ZOOM is safe only because `span` below is the radius of the
    # CIRCUMSCRIBED SPHERE rather than half the bounding box: a sphere
    # projects to the same circle from every direction, so no rotation can
    # push an atom out of frame. Measured over 48 elev/azim combinations,
    # sphere span at zoom sqrt(3) clips nothing while bounding-box span at the
    # same zoom clips 2 of 48 -- and mean ink area goes 17.5% -> 52.3%.
    center = coords.mean(axis=0)
    radius = float(np.linalg.norm(coords - center, axis=1).max())
    span = (radius + 0.5) / zoom
    center = center + np.array(pan)
    ax.set_box_aspect([1, 1, 1], zoom=BOX_ZOOM)
    ax.set_position([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(center[0] - span, center[0] + span)
    ax.set_ylim(center[1] - span, center[1] + span)
    ax.set_zlim(center[2] - span, center[2] + span)
    # Every artist added to an Axes3D otherwise re-runs auto_scale_xyz over
    # all three axes. With the limits already final that is pure waste, and it
    # was the single largest cost in the old per-bond ax.plot loop.
    ax.set_autoscale_on(False)
    # One Line3DCollection and one scatter cannot be depth-sorted against each
    # other per-bond the way 280 individual artists were, so pin the order
    # explicitly instead of letting computed_zorder reassign it every draw:
    # bonds behind, spheres in front. This is also what lets the overlay
    # labels below use a plain fixed zorder.
    ax.computed_zorder = False

    # Bond search, vectorised. The pairwise Python loop this replaces called
    # np.linalg.norm once per pair -- 5.8 ms at 140 atoms against 0.34 ms
    # here, and O(N^2) *in Python*, so a 1000-atom lower layer spent about
    # 0.3 s in that loop alone before anything was drawn.
    deltas = coords[:, None, :] - coords[None, :, :]
    distances = np.sqrt(np.einsum("ijk,ijk->ij", deltas, deltas))
    is_h = elements == "H"
    # H-H is never a real bond in this project's chemistry -- just clutter.
    bonded = (distances < BOND_CUTOFF) & ~(is_h[:, None] & is_h[None, :])
    i_idx, j_idx = np.where(np.triu(bonded, 1))

    qm_bond = qm_mask[i_idx] & qm_mask[j_idx]
    segments: list = []
    seg_colors: list[str] = []
    seg_widths: list[float] = []
    seg_alphas: list[float] = []

    qi, qj = i_idx[qm_bond], j_idx[qm_bond]
    if len(qi):
        segments.append(np.stack([coords[qi], coords[qj]], axis=1))
        seg_colors += ["#cccccc"] * len(qi)
        seg_widths += [2.0] * len(qi)
        seg_alphas += [1.0] * len(qi)
    # QM-QM bonds only -- these are the ones show_distances labels.
    bonds = [
        (*coords[a], *coords[b], float(distances[a, b])) for a, b in zip(qi, qj)
    ]

    wi, wj = i_idx[~qm_bond], j_idx[~qm_bond]
    if len(wi):
        # Wireframe atoms get no sphere to carry their element colour, so each
        # half of the bond is coloured from its own endpoint -- the same
        # ELEMENT_COLORS the ball-and-stick spheres use -- to keep elements
        # distinguishable without one.
        start, end = coords[wi], coords[wj]
        middle = (start + end) / 2
        segments.append(np.stack([start, middle], axis=1))
        segments.append(np.stack([middle, end], axis=1))
        seg_colors += [ELEMENT_COLORS.get(e, DEFAULT_ELEMENT_COLOR) for e in elements[wi]]
        seg_colors += [ELEMENT_COLORS.get(e, DEFAULT_ELEMENT_COLOR) for e in elements[wj]]
        seg_widths += [1.0] * (2 * len(wi))
        seg_alphas += [0.7] * (2 * len(wi))

    if segments:
        rgba = np.array([to_rgba(c, a) for c, a in zip(seg_colors, seg_alphas)])
        ax.add_collection3d(
            Line3DCollection(
                np.concatenate(segments), colors=rgba, linewidths=seg_widths, zorder=1,
            )
        )

    qm_indices = np.where(qm_mask)[0]
    if len(qm_indices):
        sphere_coords = coords[qm_indices]
        ax.scatter(
            sphere_coords[:, 0], sphere_coords[:, 1], sphere_coords[:, 2],
            s=[
                (ELEMENT_RADII.get(e, DEFAULT_ELEMENT_RADIUS) * 180) ** 1.1
                for e in elements[qm_indices]
            ],
            color=[ELEMENT_COLORS.get(e, DEFAULT_ELEMENT_COLOR) for e in elements[qm_indices]],
            edgecolor="black", linewidth=0.5, depthshade=True, zorder=2,
        )

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
    # Axes3D.computed_zorder used to reassign the zorder of every Collection
    # -- and each atom's own ax.scatter call was its own Path3DCollection --
    # to zorder_offset..zorder_offset+n_atoms, sorted by camera depth, EVERY
    # draw. A fixed zorder on the label (10, originally) only beat that for
    # small molecules; past ~10 atoms some sphere's reassigned zorder exceeded
    # it, and WHICH one did depended on the view angle, so a label flickered
    # in and out from behind its own atom as you rotated. computed_zorder is
    # off now and the drawing is two artists rather than n_atoms + n_bonds, so
    # nothing reassigns anything -- but `_OVERLAY_ZORDER` stays deliberately
    # far above both of them rather than becoming a third fragile "3".
    projection = ax.get_proj()
    # `show_labels` off leaves the spheres alone: on a crowded QM region the
    # index labels cover the structure they are there to identify.
    for index in (qm_indices if show_labels else ()):
        x, y, z = coords[index]
        lx, ly, _ = proj3d.proj_transform(x, y, z, projection)
        ax.annotate(
            str(index), xy=(lx, ly), xycoords="data", ha="center", va="center",
            color="white", fontsize=7, fontweight="bold", zorder=_OVERLAY_ZORDER,
            path_effects=[pe.withStroke(linewidth=2, foreground="black")],
            bbox=dict(boxstyle="round,pad=0.1,rounding_size=0.3", fc="black", ec="none", alpha=0.65),
        )

    if show_distances:
        for x1, y1, z1, x2, y2, z2, d in bonds:
            mx, my, mz = (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2
            lx, ly, _ = proj3d.proj_transform(mx, my, mz, projection)
            ax.annotate(
                f"{d:.2f}", xy=(lx, ly), xycoords="data", ha="center", va="center",
                color="#ffe066", fontsize=7, fontweight="bold", zorder=_OVERLAY_ZORDER - 1,
                path_effects=[pe.withStroke(linewidth=2, foreground="black")],
                bbox=dict(boxstyle="round,pad=0.1,rounding_size=0.3", fc="black", ec="none", alpha=0.65),
            )

    return _to_image(fig)
