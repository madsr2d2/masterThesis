"""Client for herdr's own pane.graphics.* socket API (herdr 0.9.0+). Real
pixel images, composited by herdr itself over a screen region -- not a
terminal image protocol (Sixel/Kitty), which do not survive this project's
herdr+ssh transport. herdr 0.8.2's graphics API existed in schema but its
cell-size negotiation never completed (`cell_size_unavailable`); 0.9.0 fixed
it.

Frames go as `png`. This module used to send raw `rgba` on the finding that
`png` "renders solid black" -- that was retested against herdr 0.9.0 on
2026-09-09 and PNG now displays correctly, so the note it rested on is gone.
The format is the whole ballgame on a slow link: the same 411x291 frame is
30,956 bytes of base64 as PNG against 637,872 as RGBA, a factor of 20.6, and
round-trips in 7.4 ms against 193.9 ms. That 193.9 ms is the "~100-200 ms"
this module used to attribute to herdr's relay of a full-quality frame; it
was the payload all along.

Supported formats are `png`, `rgb`, `rgba` and `bgra` -- herdr's own error
text enumerates them."""
from __future__ import annotations

import base64
import io
import json
import os
import socket
import time
from dataclasses import dataclass

from PIL import Image

_TIMEOUT_S = 2.0

# The socket API rejects a frame with "image_too_large" between 679,828 bytes
# (accepted) and 719,996 bytes (rejected) of base64 payload, and drops the
# connection outright above roughly 1.5 MB instead of answering. Measured by
# direct sweep against the running server; not documented anywhere.
#
# The limit is on the PAYLOAD, not the picture: a 2400x1700 frame -- 16.3 MB
# once decoded -- is accepted happily at 281,740 bytes of PNG. So there is no
# resolution ceiling worth designing around any more, only a byte one, and a
# geometry frame at full pane resolution is about 95,000 bytes. This cap
# exists to keep a pathological frame from tripping the limit, not because
# anything normal approaches it.
MAX_PNG_B64_BYTES = 600_000


def _encode_png(image: Image.Image) -> str:
    buf = io.BytesIO()
    image.convert("RGB").save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _encode_within(image: Image.Image, max_b64_bytes: int) -> tuple[str, Image.Image]:
    """PNG the image, shrinking it until the payload fits.

    Compression makes the encoded size unpredictable from the pixel count --
    a molecule on a flat background compresses about 100x, noise not at all --
    so the only honest way to respect a byte cap is to encode and look."""
    data = _encode_png(image)
    while len(data) > max_b64_bytes and min(image.width, image.height) > 16:
        image = image.resize(
            (max(1, image.width * 3 // 4), max(1, image.height * 3 // 4)),
            Image.LANCZOS,
        )
        data = _encode_png(image)
    return data, image


def _endpoint() -> tuple[str, str] | None:
    if os.environ.get("HERDR_ENV") != "1":
        return None
    socket_path = os.environ.get("HERDR_SOCKET_PATH")
    pane_id = os.environ.get("HERDR_PANE_ID")
    if not socket_path or not pane_id:
        return None
    return socket_path, pane_id


def _call(method: str, params: dict) -> dict | None:
    endpoint = _endpoint()
    if endpoint is None:
        return None
    socket_path, _ = endpoint
    request = {"id": f"monitor:{time.time_ns()}", "method": method, "params": params}
    try:
        client = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        client.settimeout(_TIMEOUT_S)
        client.connect(socket_path)
        client.sendall((json.dumps(request) + "\n").encode())
        raw = client.recv(1 << 20)
        client.close()
    except OSError:
        return None
    try:
        return json.loads(raw.decode())
    except json.JSONDecodeError:
        return None


@dataclass
class CellSize:
    width_px: int
    height_px: int


_cell_size: CellSize | None = None
_cell_size_probed = False


def cell_size(*, refresh: bool = False) -> CellSize | None:
    """Pixel size of one terminal cell, cached.

    This was probed on EVERY frame, which put a full request/response over the
    herdr socket in front of every rotation before the image was even
    rendered -- a round trip of pure latency on a link where latency is the
    thing being economised. It only changes when the terminal's font size
    does, which also resizes the pane, so `refresh=True` from a resize handler
    is the invalidation this needs."""
    global _cell_size, _cell_size_probed
    if _cell_size_probed and not refresh:
        return _cell_size
    _cell_size, _cell_size_probed = _probe_cell_size(), True
    return _cell_size


def _probe_cell_size() -> CellSize | None:
    endpoint = _endpoint()
    if endpoint is None:
        return None
    _, pane_id = endpoint
    response = _call("pane.graphics.info", {"pane_id": pane_id})
    if response is None or "error" in response:
        return None
    result = response.get("result", {})
    width, height = result.get("cell_width_px"), result.get("cell_height_px")
    if not width or not height:
        return None
    return CellSize(width, height)


def available() -> bool:
    return cell_size() is not None


def set_image(
    image: Image.Image,
    *,
    grid_cols: int,
    grid_rows: int,
    viewport_col: int,
    viewport_row: int,
    layer_id: str,
    max_b64_bytes: int = MAX_PNG_B64_BYTES,
) -> bool:
    endpoint = _endpoint()
    if endpoint is None or grid_cols <= 0 or grid_rows <= 0:
        return False
    _, pane_id = endpoint
    data_b64, sent = _encode_within(image, max_b64_bytes)
    response = _call(
        "pane.graphics.set",
        {
            "pane_id": pane_id,
            "format": "png",
            "image_width": sent.width,
            "image_height": sent.height,
            "data_base64": data_b64,
            "layer_id": layer_id,
            "placement": {
                "grid_cols": grid_cols,
                "grid_rows": grid_rows,
                "viewport_col": viewport_col,
                "viewport_row": viewport_row,
            },
        },
    )
    return bool(response and response.get("result", {}).get("type") == "ok")


def clear(layer_id: str) -> None:
    endpoint = _endpoint()
    if endpoint is None:
        return
    _, pane_id = endpoint
    _call("pane.graphics.clear", {"pane_id": pane_id, "layer_id": layer_id})
