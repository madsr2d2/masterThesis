"""Client for herdr's own pane.graphics.* socket API (herdr 0.9.0+). Real
pixel images, composited by herdr itself over a screen region -- not a
terminal image protocol (Sixel/Kitty), which do not survive this project's
herdr+ssh transport. See computational/CLAUDE.md-adjacent conversation
history for why: herdr 0.8.2's graphics API existed in schema but its
cell-size negotiation never completed (`cell_size_unavailable`); 0.9.0 fixed
it, and raw `rgba` data works where `png` format still renders solid black
(PNG decoding appears unfinished as of 0.9.0)."""
from __future__ import annotations

import base64
import json
import os
import socket
import time
from dataclasses import dataclass

from PIL import Image

_TIMEOUT_S = 2.0

# The socket API rejects a frame with "image_too_large" somewhere between
# 498,436 bytes (works) and 547,600 bytes (rejected) of raw RGBA -- measured
# by a direct sweep against the running herdr server, not documented anywhere.
# Above roughly 850,000 bytes the connection is dropped outright mid-write
# instead of getting a clean error.
MAX_RGBA_BYTES = 480_000

# A frame this size (or the ~50,000-byte range generally) round-trips in
# single-digit milliseconds -- also measured directly. The real cost at
# MAX_RGBA_BYTES is the ~100-200ms herdr spends relaying the frame to the
# outer terminal over ssh, which is data-transfer time, not fixed overhead,
# so a small interim frame is genuinely fast rather than merely "smaller."
PREVIEW_MAX_RGBA_BYTES = 40_000


def _clamp_to_budget(image: Image.Image, max_bytes: int) -> Image.Image:
    raw_bytes = image.width * image.height * 4
    if raw_bytes <= max_bytes:
        return image
    scale = (max_bytes / raw_bytes) ** 0.5
    new_size = (max(1, int(image.width * scale)), max(1, int(image.height * scale)))
    return image.resize(new_size)


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


def cell_size() -> CellSize | None:
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
    max_bytes: int = MAX_RGBA_BYTES,
) -> bool:
    endpoint = _endpoint()
    if endpoint is None or grid_cols <= 0 or grid_rows <= 0:
        return False
    _, pane_id = endpoint
    rgba = _clamp_to_budget(image.convert("RGBA"), max_bytes)
    data_b64 = base64.b64encode(rgba.tobytes()).decode("ascii")
    response = _call(
        "pane.graphics.set",
        {
            "pane_id": pane_id,
            "format": "rgba",
            "image_width": rgba.width,
            "image_height": rgba.height,
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
