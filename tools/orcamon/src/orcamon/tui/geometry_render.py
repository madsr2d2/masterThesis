"""The geometry pane's pixels: a thin wrapper over `raster`, the z-buffered
sprite renderer.

Imported LAZILY by the two pixel widgets and by `orcamon snapshot`, so the
TUI and the CLI start without the `images` extra; without it the import
raises and the caller says what to install. There is no matplotlib any more
-- the frame is numpy and PIL (`raster`), and `frame_png` turns it into the
palette PNG the pane transmits, exactly as it did before.
"""
from __future__ import annotations

import io

from PIL import Image

from ..core.geometry import View
from . import raster


def frame_png(image: Image.Image) -> bytes:
    """A rendered frame as the PNG sent to the terminal: 256-colour palette,
    undithered.

    A frame is a flat background, flat element colours and their shaded
    spheres -- far fewer distinct colours than a truecolour render of the
    same size -- and a palette PNG carries it in a few tens of KB against
    several times that as truecolour, at half the encode time. Mean error is
    ~1 level in 255 and side by side the two are indistinguishable. On an
    ssh-relayed pane, bytes ARE the latency (see herdr_graphics), so this is
    the frame rate of a rotation."""
    buf = io.BytesIO()
    image.convert("RGB").quantize(
        colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.NONE
    ).save(buf, format="PNG")
    return buf.getvalue()


def render(
    atoms: list[tuple[str, float, float, float]],
    view: View | None = None,
    *,
    qm_atom_indices: set[int] | None = None,
    size_px: tuple[int, int] = (900, 750),
    bond_list: list[tuple[int, int]] | None = None,
) -> Image.Image:
    """A frame of `atoms` under `view`.

    `qm_atom_indices` (0-based, matching ORCA's numbering) marks a multilayer
    job's high-level layer: those atoms get the full representation, every
    other atom only the thin desaturated environment bonds. None (an ordinary
    job) treats every atom as the high-level layer.

    `bond_list` may be passed in when the caller caches it -- a geometry's
    bonds do not change as it is rotated."""
    return raster.render(atoms, view, qm_atom_indices=qm_atom_indices,
                         size_px=size_px, bond_list=bond_list)
