"""Where the TUI and the CLI agree the other one is: the shared state files.

The person's TUI and the agent's CLI sit over the same machine -- under tmux,
herdr, screen or plain ssh, in one session or in two -- and are wired together
by files and nothing else: no daemon, no socket, no multiplexer coupling, so
the arrangement does not depend on which terminal holds which side (D3).

D3: the state lives under `$XDG_STATE_HOME/orcamon`, or
`~/.local/state/orcamon` when `XDG_STATE_HOME` is unset or empty.
`view-<host>-<pid>.json` is written only by the TUI with that pid and says
what it shows; `focus-<host>-<pid>.json` is written only by the CLI and read
(renamed, then deleted) by that TUI. Every write is atomic -- a temporary
file in the same directory, then `os.replace` -- so a reader never sees half
a document.

D4: a TUI is live when its view file is on THIS host -- `<host>` is
`socket.gethostname()` with every character outside `A-Za-z0-9._-` replaced
by `_` -- and its pid still exists (`os.kill(pid, 0)` does not raise
`ProcessLookupError`; a `PermissionError` counts as live). The CLI talks to
the live TUI with the largest `updated_unix_s`, the one changed most
recently; another host's files are never read. A TUI removes its own two
files when it exits and, when it starts, the files of dead pids on this
host.
"""
from __future__ import annotations

import json
import os
import re
import socket
from pathlib import Path

from .geometry import FileGeometry
from .paths import PathPoint

STEER_SCHEMA = 1


def state_dir() -> Path:
    """`$XDG_STATE_HOME/orcamon`, or `~/.local/state/orcamon` when the
    variable is unset or empty (D3). Not `mkdir`ed here: a reader must not
    create the directory just by looking."""
    root = os.environ.get("XDG_STATE_HOME")
    if root:
        return Path(root) / "orcamon"
    return Path.home() / ".local" / "state" / "orcamon"


def host_name() -> str:
    """The host component of every state file name (D4): the machine's name
    with anything that could complicate a file name -- or a glob -- replaced
    by `_`, and a fallback for a hostname that is empty."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", socket.gethostname()) or "localhost"


def view_path(pid: int | None = None) -> Path:
    """This TUI's view file, or another pid's when one is given (D3)."""
    return state_dir() / f"view-{host_name()}-{pid or os.getpid()}.json"


def focus_path(pid: int) -> Path:
    """The focus file the CLI writes for the TUI with `pid` (D3)."""
    return state_dir() / f"focus-{host_name()}-{pid}.json"


def write_json(path: Path, doc) -> None:
    """Write `doc` at `path` atomically: a temporary file beside it, then
    rename (D3). A reader never sees a half-written document, and a crash
    leaves the previous one in place."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(doc))
    os.replace(tmp, path)


def read_json(path: Path) -> dict | None:
    """The document at `path`, or None for a missing, unreadable, malformed
    or non-object file: a reader treats all of them as nothing to read."""
    try:
        doc = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    return doc if isinstance(doc, dict) else None


def alive(pid: int) -> bool:
    """Whether the process `pid` exists. Signal 0 delivers nothing; it only
    asks the kernel whether the pid does (D4). A process owned by another
    user is live but unsignallable, hence the PermissionError."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def live_views() -> list[dict]:
    """Every live TUI's view document on this host, most recently updated
    first (D4). A file whose name does not carry an integer pid, whose
    document is not a dict or whose pid does not match the name is not a
    view; another host's files are never read."""
    prefix = f"view-{host_name()}-"
    views = []
    for path in state_dir().glob(f"{prefix}*.json"):
        stem = path.name[len(prefix):-len(".json")]
        if not stem.isdigit():
            continue
        pid = int(stem)
        doc = read_json(path)
        if doc is None or doc.get("pid") != pid or not alive(pid):
            continue
        views.append(doc)
    views.sort(key=lambda doc: doc.get("updated_unix_s", 0.0), reverse=True)
    return views


def current_view() -> dict | None:
    """The live TUI the CLI talks to: the one updated most recently, or None
    when no TUI is running on this host (D4)."""
    views = live_views()
    return views[0] if views else None


def take_focus(pid: int) -> dict | None:
    """Claim the focus request waiting for TUI `pid`, or None when there is
    none. The rename claims it atomically, so a request the CLI writes
    meanwhile is kept for the next check."""
    path = focus_path(pid)
    claimed = path.with_name(path.name + ".taken")
    try:
        os.replace(path, claimed)
    except FileNotFoundError:
        return None
    doc = read_json(claimed)
    claimed.unlink(missing_ok=True)
    return doc


def remove_stale() -> int:
    """Remove this host's view and focus files whose pid no longer exists,
    and say how many went (D4). A dead TUI cannot clean up after itself, and
    another host's files are never touched."""
    removed = 0
    for kind in ("view", "focus"):
        prefix = f"{kind}-{host_name()}-"
        for path in state_dir().glob(f"{prefix}*.json"):
            stem = path.name[len(prefix):-len(".json")]
            if stem.isdigit() and not alive(int(stem)):
                path.unlink(missing_ok=True)
                removed += 1
    return removed


def geometry_address(shown, following: bool) -> list[str]:
    """`orcamon geom`'s own flags for the geometry the pane shows: the path
    point it has scrubbed to, or the step and cycle of the printed one.
    A path point is addressed even while following: on an IRC or NEB the
    pane draws the path's focus, not the job's latest printed geometry, so
    `geom`'s default would return a different structure. Empty while it
    follows the newest printed geometry, where that default is already
    what the pane draws."""
    if isinstance(shown, PathPoint):
        return ["--point", str(shown.index)]
    if following or shown is None or isinstance(shown, FileGeometry):
        return []
    args = []
    if shown.scan_step is not None:
        args += ["--step", str(shown.scan_step)]
    if shown.cycle:
        args += ["--cycle", str(shown.cycle)]
    return args
