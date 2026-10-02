"""Start and stop ORCA: find the binary, tell a login node, end a tree.

The two places a mistake costs the most are starting a job and stopping one,
so both live here, each with one way to do it (D4-D7); a wrong binary or a
login node is refused before anything runs.

D4 -- start. `start_orca` runs `[<absolute ORCA path>, "<stem>.inp"]` in the
job's directory, with stdout and stderr into a NEW `<stem>.out` (mode "x", so
an output that already exists is never overwritten), stdin from /dev/null and
`start_new_session=True` so the agent's shell exiting or timing out does not
take the job down. It watches for `settle` seconds: an ORCA that exits in that
time is reported with its status and the output's last lines.

D5 -- which binary. `--orca PATH`, else `$ORCAMON_ORCA`, else `orca` on PATH,
resolved with `os.path.realpath`. It must be an executable file whose
directory holds at least `MIN_ORCA_MODULES` entries named `orca_*` -- the real
ORCA's own modules, which no other program called `orca` (the GNOME screen
reader is one) has.

D6 -- where. A SLURM login node is a machine where `squeue` and `sinfo` are
on PATH and `sinfo -h -N -o %N` lists nodes that do not include this machine's
hostname or its short name. ORCA started there would run on the login node,
so `run` refuses unless `--local`; where `sinfo` cannot answer, nothing is
refused.

D7 -- stop. `stop_tree` signals only the job's own processes: every process
named `orca` (its `/proc/<pid>/comm`) whose working directory is the job's and
whose command line names `<stem>.inp`, plus all their descendants. SIGTERM to
all of them, then SIGKILL to any still alive after `grace` seconds. A zombie
counts as gone, and `os.waitpid` reaps one where it is ours to reap.

Linux `/proc`, standard library only: this is what an agents-only server has.
"""
from __future__ import annotations

import os
import shutil
import signal
import socket
import subprocess
import time
from pathlib import Path

from . import liveness

# The real ORCA ships ~90 `orca_*` modules beside the driver. Ten is the bar
# that admits a minimal installation while refusing another program called
# `orca`: the GNOME screen reader has no such neighbours.
MIN_ORCA_MODULES = 10

# `sinfo` can hang behind a busy scheduler; a login-node check must not hold
# every `run` up for longer than this.
SINFO_TIMEOUT_S = 10.0


class RunError(Exception):
    """A job cannot be run or stopped as asked. The message says what to do."""


def find_orca(explicit: str | None = None) -> Path:
    """The ORCA binary D5: `explicit`, else `$ORCAMON_ORCA`, else `orca`."""
    candidate = explicit or os.environ.get("ORCAMON_ORCA") or shutil.which("orca")
    if not candidate:
        raise RunError("ORCA not found: pass --orca PATH or set ORCAMON_ORCA")
    path = Path(os.path.realpath(candidate))
    if not path.is_file() or not os.access(path, os.X_OK):
        raise RunError(f"{path} is not an executable file")
    try:
        modules = [entry for entry in path.parent.iterdir()
                   if entry.name.startswith("orca_")]
    except OSError:
        modules = []  # a directory we cannot read has no modules we can see
    if len(modules) < MIN_ORCA_MODULES:
        raise RunError(f"{path} is not ORCA: no orca_* programs beside it "
                       "(pass --orca PATH or set ORCAMON_ORCA)")
    return path


def slurm_login_node() -> str | None:
    """D6: why this machine is a SLURM login node, or None when it is not one
    -- or when `sinfo` cannot say, which refuses nothing."""
    if not (shutil.which("squeue") and shutil.which("sinfo")):
        return None
    try:
        done = subprocess.run(["sinfo", "-h", "-N", "-o", "%N"],
                              capture_output=True, text=True,
                              timeout=SINFO_TIMEOUT_S)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if done.returncode != 0:
        return None
    nodes = [line.strip() for line in done.stdout.splitlines() if line.strip()]
    host = socket.gethostname()
    short = host.split(".")[0]
    if host in nodes or short in nodes:
        return None
    listed = ", ".join(sorted(nodes)[:5]) + (", ..." if len(nodes) > 5 else "")
    return (f"this machine ({short}) is not one of SLURM's nodes ({listed}): "
            "ORCA started here would run on the login node. Submit it to the queue "
            "instead (orcamon cannot yet), or pass --local if this machine is meant "
            "to run it")


def _state(pid: int) -> str | None:
    """The single-letter state in `/proc/<pid>/stat` (`R`, `S`, `Z`, ...) or
    None when it cannot be read. The command name in field 2 is parenthesised
    and may hold spaces and `)`, so the state is the first field after the
    LAST `)`."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return None
    end = stat.rfind(")")
    if end < 0:
        return None
    rest = stat[end + 1:].split()
    return rest[0] if rest else None


def alive(pid: int) -> bool:
    """Whether `pid` is a running process. A zombie (`Z`) is not alive: it has
    already exited, and only its parent's `wait` is owed."""
    state = _state(pid)
    return state is not None and state != "Z"


def driver_pids(job_dir: Path, stem: str) -> list[int]:
    """D7: the pids of the `orca` drivers running `stem.inp` in `job_dir`.

    A module process (`orca_scf_mpi`) names no input, and one directory can
    hold several jobs, so the driver's own command line decides. Another
    user's process files are unreadable and skipped -- those are not ours to
    signal."""
    wanted = Path(job_dir).resolve()
    found = []
    try:
        entries = [entry for entry in os.scandir("/proc") if entry.name.isdigit()]
    except OSError:
        return []
    for entry in entries:
        pid = int(entry.name)
        base = f"/proc/{pid}"
        try:
            with open(f"{base}/comm") as f:
                name = f.read().strip()
            if name != "orca":
                continue
            cwd = Path(os.readlink(f"{base}/cwd")).resolve()
            if cwd != wanted:
                continue
            with open(f"{base}/cmdline", "rb") as f:
                args = f.read().decode("utf-8", "replace").split("\0")
        except OSError:
            continue  # gone or not ours since the scan began
        if liveness.input_stem(args) == stem and alive(pid):
            found.append(pid)
    return sorted(found)


def descendants(pid: int) -> list[int]:
    """Every process below `pid` in the parent chain, `pid` itself excluded.

    Built from every `/proc/<n>/stat` in one pass: asking each live process
    for its children would race the tree changing under the walk."""
    children: dict[int, list[int]] = {}
    try:
        entries = [entry for entry in os.scandir("/proc") if entry.name.isdigit()]
    except OSError:
        return []
    for entry in entries:
        try:
            stat = Path(f"/proc/{entry.name}/stat").read_text()
        except OSError:
            continue  # exited between the scan and the read
        end = stat.rfind(")")
        rest = stat[end + 1:].split() if end >= 0 else []
        if len(rest) < 2:
            continue  # the second field after ')' is the parent pid
        children.setdefault(int(rest[1]), []).append(int(entry.name))
    found: list[int] = []
    seen: set[int] = set()
    stack = [pid]
    while stack:
        for child in children.get(stack.pop(), []):
            if child not in seen:
                seen.add(child)
                found.append(child)
                stack.append(child)
    return sorted(found)


def start_orca(orca: Path, job_dir: Path, stem: str,
               settle: float) -> tuple[int, int | None, list[str]]:
    """D4: start ORCA on `<stem>.inp` in `job_dir`, detached, and watch it for
    `settle` seconds. Returns `(pid, returncode, tail)`: the code is None and
    the tail empty while it is still running."""
    job_dir = Path(job_dir)
    out_path = job_dir / f"{stem}.out"
    # Mode "x": a second run must never overwrite the first one's output.
    with open(out_path, "x") as out:
        proc = subprocess.Popen(
            [str(orca), f"{stem}.inp"],
            cwd=job_dir,
            stdin=subprocess.DEVNULL,
            stdout=out,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            close_fds=True,
        )
    deadline = time.monotonic() + settle
    while proc.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    if proc.poll() is None:
        return proc.pid, None, []
    lines = out_path.read_text(errors="replace").splitlines()
    return proc.pid, proc.returncode, lines[-10:]


def stop_tree(roots: list[int], grace: float) -> tuple[list[int], list[int]]:
    """D7: SIGTERM the roots and their descendants, then SIGKILL whatever is
    still alive after `grace` seconds. Returns `(pids, killed)`; only the
    pids passed in -- found by `driver_pids` -- are ever signalled."""
    pids: set[int] = set()
    for root in roots:
        if alive(root):
            pids.add(root)
        pids.update(descendants(root))
    pids = sorted(pids)
    for pid in pids:
        _signal(pid, signal.SIGTERM)
    _wait_gone(pids, grace)
    killed = [pid for pid in pids if alive(pid)]
    for pid in killed:
        _signal(pid, signal.SIGKILL)
    _wait_gone(pids, 2.0)
    return pids, killed


def _signal(pid: int, sig: int) -> None:
    """One signal, with the two ways it can fail named: a pid that is already
    gone is fine, a pid that is not ours is not."""
    try:
        os.kill(pid, sig)
    except ProcessLookupError:
        pass
    except PermissionError:
        raise RunError(f"cannot signal pid {pid}: not yours") from None


def _wait_gone(pids: list[int], seconds: float) -> None:
    """Poll `alive` every 0.1 s until none of `pids` is, or `seconds` pass.
    `os.waitpid` is what turns our own dead child from a zombie into a gone
    pid; for a pid that is not ours it raises and there is nothing to do."""
    deadline = time.monotonic() + seconds
    while True:
        for pid in pids:
            try:
                os.waitpid(pid, os.WNOHANG)
            except ChildProcessError:
                pass
        if not any(alive(pid) for pid in pids):
            return
        if time.monotonic() >= deadline:
            return
        time.sleep(0.1)
