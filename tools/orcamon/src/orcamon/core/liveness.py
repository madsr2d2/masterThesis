"""Is a job's ORCA process alive? Read from /proc, standard library only.

This used psutil, which made a third-party package a hard requirement of the
one question every refresh asks. On Linux the answer is three small files per
process; psutil stays as a fallback for platforms without /proc, and only
when it happens to be installed.
"""
from __future__ import annotations

import os
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
