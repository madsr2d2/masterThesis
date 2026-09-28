"""The Kitty graphics protocol, as the strings the geometry pane writes.

Absolute placement: the cursor is saved, moved to the pane's top-left cell,
the image is drawn scaled into `cols` x `rows` cells, and the cursor is
restored. `q=2` suppresses the terminal's per-command acknowledgement, which
otherwise comes back up the wire and lands in Textual's own input parser.
Pure string building, so it is testable without a terminal."""
from __future__ import annotations

import base64

CHUNK = 4096  # the protocol's maximum base64 payload per escape


def _at(col: int, row: int, body: str) -> str:
    return f"\x1b[s\x1b[{row + 1};{col + 1}H{body}\x1b[u"


def transmit(png: bytes, image_id: int, col: int, row: int, cols: int, rows: int) -> str:
    """Send `png` and display it (`a=T`), storing it under `image_id` so later
    re-placements need not resend the pixels."""
    data = base64.b64encode(png).decode("ascii")
    chunks = [data[i : i + CHUNK] for i in range(0, len(data), CHUNK)] or [""]
    parts = []
    for i, chunk in enumerate(chunks):
        controls = ["q=2"]
        if i == 0:
            controls += ["a=T", "f=100", "t=d", f"i={image_id}", f"c={cols}", f"r={rows}"]
        controls.append("m=1" if i != len(chunks) - 1 else "m=0")
        parts.append(f"\x1b_G{','.join(controls)};{chunk}\x1b\\")
    return _at(col, row, "".join(parts))


def place(image_id: int, col: int, row: int, cols: int, rows: int) -> str:
    """Re-place an image the terminal already holds: ~50 bytes."""
    return _at(col, row, f"\x1b_Ga=p,i={image_id},q=2,c={cols},r={rows}\x1b\\")


def delete(image_id: int) -> str:
    return f"\x1b_Ga=d,d=i,i={image_id}\x1b\\"
