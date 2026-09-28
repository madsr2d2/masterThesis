"""Which geometry pane this terminal can show: herdr, kitty or text.

    auto  herdr when herdr's graphics API answers; text inside tmux or
          screen (Kitty graphics through a multiplexer is fragile, and
          `--graphics kitty` forces it); otherwise ASK the terminal, and use
          kitty only if it says it supports the protocol.
    text  always works, everywhere.

The pixel modes also need the `images` extra; without it the answer is
text, with a note saying why.

Asking is a Kitty graphics query followed by a primary device-attributes
request (DA1). Every terminal answers DA1; one that speaks the Kitty
protocol answers the query first. So: OK before DA1 means kitty, DA1 alone
means no, silence past the timeout means no. This runs before Textual
starts, with the terminal briefly in raw mode and restored in a `finally`.
"""
from __future__ import annotations

import os
import re
import select
import sys
import time

GRAPHICS_MODES = ("auto", "herdr", "kitty", "text")
PROBE_TIMEOUT_S = 0.5
KITTY_QUERY = "\x1b_Gi=31,s=1,v=1,a=q,t=d,f=24;AAAA\x1b\\"
DA1_QUERY = "\x1b[c"

_KITTY_OK_RE = re.compile(rb"\x1b_Gi=31;OK\x1b\\")
_DA1_RE = re.compile(rb"\x1b\[\?[\d;]*c")


def parse_probe_reply(reply: bytes) -> bool | None:
    """True: the Kitty OK arrived before DA1. False: DA1 without it. None:
    no DA1 yet (keep reading, or give up at the timeout)."""
    da1 = _DA1_RE.search(reply)
    if da1 is None:
        return None
    ok = _KITTY_OK_RE.search(reply)
    return ok is not None and ok.start() < da1.start()


def images_available() -> bool:
    # Found, not imported: importing matplotlib costs half a second of
    # startup, and the pixel widget imports it lazily anyway.
    import importlib.util
    return all(importlib.util.find_spec(m) is not None for m in ("matplotlib", "numpy", "PIL"))


def probe_kitty(timeout: float = PROBE_TIMEOUT_S) -> bool:
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        return False
    try:
        import termios
        import tty
    except ImportError:
        return False
    fd = sys.stdin.fileno()
    try:
        saved = termios.tcgetattr(fd)
    except termios.error:
        return False
    reply = b""
    try:
        tty.setraw(fd)
        sys.stdout.write(KITTY_QUERY + DA1_QUERY)
        sys.stdout.flush()
        deadline = time.monotonic() + timeout
        while True:
            left = deadline - time.monotonic()
            if left <= 0:
                return False
            ready, _, _ = select.select([fd], [], [], left)
            if not ready:
                return False
            reply += os.read(fd, 1024)
            verdict = parse_probe_reply(reply)
            if verdict is not None:
                return verdict
    except OSError:
        return False
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, saved)


def choose(requested: str = "auto", environ=None, herdr_available=None, kitty_probe=None,
           have_images=None) -> tuple[str, str | None]:
    """(mode, note). The hooks are for tests; each defaults to asking."""
    environ = os.environ if environ is None else environ
    if have_images is None:
        have_images = images_available()
    if requested == "text":
        return "text", None
    if requested in ("herdr", "kitty"):
        if not have_images:
            return "text", "images extra not installed"
        return requested, None
    if herdr_available is None:
        from . import herdr_graphics
        herdr_available = herdr_graphics.available()
    if herdr_available:
        return ("herdr", None) if have_images else ("text", "images extra not installed")
    if environ.get("TMUX") or environ.get("STY"):
        return "text", None
    if not have_images:
        return "text", "images extra not installed"
    kitty = (kitty_probe or probe_kitty)()
    return ("kitty", None) if kitty else ("text", None)
