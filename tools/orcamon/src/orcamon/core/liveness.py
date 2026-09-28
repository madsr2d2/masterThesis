"""Is a job's ORCA process alive? Read from /proc, standard library only.

This used psutil, which made a third-party package a hard requirement of the
one question every refresh asks. On Linux the answer is three small files per
process; psutil stays as a fallback for platforms without /proc, and only
when it happens to be installed.
"""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

_PROC = Path("/proc")


def running_orca_cwds() -> dict[Path, float]:
    """Working directory -> start time (epoch s) of the ORCA run there.

    ORCA is a driver (`orca`) that spawns per-module children (`orca_scf_mpi`,
    `orca_leanscf`, ...), all in the job's directory. The driver's start time
    is the run's; when it cannot be seen (another user's, or already reaped),
    the oldest `orca*` process in the directory stands in for it."""
    if _PROC.is_dir():
        entries = _proc_entries()
    else:
        entries = _psutil_entries()

    by_cwd: dict[Path, list[tuple[str, float]]] = {}
    for name, cwd, started in entries:
        by_cwd.setdefault(cwd, []).append((name, started))

    result: dict[Path, float] = {}
    for cwd, found in by_cwd.items():
        main = [t for name, t in found if name == "orca"]
        result[cwd] = min(main) if main else min(t for _, t in found)
    return result


def _boot_time() -> float:
    with open(_PROC / "stat") as f:
        for line in f:
            if line.startswith("btime "):
                return float(line.split()[1])
    raise OSError("no btime in /proc/stat")


def _proc_entries() -> list[tuple[str, Path, float]]:
    ticks = os.sysconf("SC_CLK_TCK")
    boot = _boot_time()
    found = []
    for entry in os.scandir(_PROC):
        if not entry.name.isdigit():
            continue
        base = f"/proc/{entry.name}"
        try:
            # `comm` is truncated to 15 characters by the kernel, so a module
            # like `orca_scf_mpi` survives but only the prefix is reliable.
            with open(f"{base}/comm") as f:
                name = f.read().strip()
            if not name.startswith("orca"):
                continue
            cwd = Path(os.readlink(f"{base}/cwd"))
            with open(f"{base}/stat") as f:
                stat = f.read()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        # Field 22 is the start time in clock ticks since boot. The command
        # name (field 2) is parenthesised and may itself hold spaces or ')',
        # so count fields from the LAST ')'.
        rest = stat[stat.rfind(")") + 2 :].split()
        started = boot + int(rest[19]) / ticks
        found.append((name, _resolve(cwd), started))
    return found


def _psutil_entries() -> list[tuple[str, Path, float]]:
    try:
        import psutil
    except ImportError:
        return []
    found = []
    for proc in psutil.process_iter(["name", "cwd", "create_time"]):
        try:
            name = proc.info["name"] or ""
            cwd = proc.info["cwd"]
            if not name.startswith("orca") or cwd is None:
                continue
            found.append((name, _resolve(Path(cwd)), proc.info["create_time"]))
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            continue
    return found


def _resolve(path: Path) -> Path:
    try:
        return path.resolve()
    except OSError:
        return path


@dataclass
class Liveness:
    """What one source can say about whether a job is running.

    `alive` is three-valued on purpose. A process table can say "not here",
    but a login node cannot see a compute node's processes and a file's age
    cannot see anything -- and a source that cannot tell must say None, so
    that the status step never turns "I don't know" into "stopped"."""

    alive: bool | None
    source: str                 # "process" | "slurm" | "mtime"
    since: float | None = None  # process/job start (epoch s), for wall time
    queued: bool = False        # the scheduler has it pending
    sched_id: str | None = None
    note: str | None = None     # why the source could not decide, if it could not


class LivenessProbe:
    """One snapshot per refresh, keyed by resolved job directory; `default`
    answers for a directory the snapshot does not mention."""

    name = "?"

    def snapshot(self, job_dirs: list[Path] | None = None) -> dict[Path, Liveness]:
        raise NotImplementedError

    def default(self, job_dir: Path) -> Liveness:
        raise NotImplementedError


class ProcessProbe(LivenessProbe):
    """This machine's process table. A job directory with no `orca*` process
    in it is not running -- on this machine, which is the only claim made."""

    name = "process"

    def snapshot(self, job_dirs=None):
        return {cwd: Liveness(True, self.name, since=t) for cwd, t in running_orca_cwds().items()}

    def default(self, job_dir):
        return Liveness(False, self.name)


class MtimeProbe(LivenessProbe):
    """No process information at all -- a copied tree, a mounted remote
    filesystem. It never says a job is dead; the status step reads the
    output's age instead."""

    name = "mtime"

    def snapshot(self, job_dirs=None):
        return {}

    def default(self, job_dir):
        return Liveness(None, self.name)


# The extension point for schedulers: a name on the command line -> a probe.
PROBES: dict[str, type[LivenessProbe]] = {
    "process": ProcessProbe,
    "mtime": MtimeProbe,
}

LIVENESS_MODES = ("auto", *PROBES)


def make_probe(mode: str = "auto") -> LivenessProbe:
    if mode == "auto":
        mode = "slurm" if "slurm" in PROBES and shutil.which("squeue") else "process"
    try:
        return PROBES[mode]()
    except KeyError:
        raise ValueError(f"unknown liveness mode {mode!r}; one of {', '.join(LIVENESS_MODES)}") from None


def lookup(snapshot: dict[Path, Liveness], probe: LivenessProbe, job_dir: Path) -> Liveness:
    return snapshot.get(_resolve(job_dir)) or probe.default(job_dir)
