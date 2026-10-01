"""Is a job's ORCA process alive? Read from /proc, standard library only.

This used psutil, which made a third-party package a hard requirement of the
one question every refresh asks. On Linux the answer is three small files per
process; psutil stays as a fallback for platforms without /proc, and only
when it happens to be installed.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

_PROC = Path("/proc")


def is_orca_process(name: str) -> bool:
    """The ORCA driver (`orca`) or one of its modules (`orca_scf_mpi`, ...).
    Not `orcamon`, whose console script would otherwise count itself."""
    return name == "orca" or name.startswith("orca_")


def input_stem(args: list[str]) -> str | None:
    """The job an `orca` driver was started on: its first argument that is
    an input file, as a stem (`/a/b/job.inp` -> `job`); None when none is."""
    for arg in args[1:]:
        if arg and not arg.startswith("-") and arg.endswith(".inp"):
            return Path(arg).name[:-4]
    return None


def running_orca_cwds() -> dict[Path, float]:
    """Working directory -> start time (epoch s) of the ORCA run there.

    ORCA is a driver (`orca`) that spawns per-module children (`orca_scf_mpi`,
    `orca_leanscf`, ...), all in the job's directory. The driver's start time
    is the run's; when it cannot be seen (another user's, or already reaped),
    the oldest `orca*` process in the directory stands in for it."""
    by_cwd: dict[Path, list[tuple[str, float]]] = {}
    for name, cwd, started, _args in _orca_entries():
        by_cwd.setdefault(cwd, []).append((name, started))

    result: dict[Path, float] = {}
    for cwd, found in by_cwd.items():
        main = [t for name, t in found if name == "orca"]
        result[cwd] = min(main) if main else min(t for _, t in found)
    return result


def running_orca_stems() -> dict[Path, frozenset[str]]:
    """Working directory -> the job stems an `orca` driver there names.

    Only the driver's own command line says which of a directory's several
    jobs it runs; a module (`orca_scf_mpi`) does not name one. A directory
    left out of the dict is one where no driver names its input, so its
    liveness stays directory-level."""
    by_cwd: dict[Path, set[str]] = {}
    for name, cwd, _started, args in _orca_entries():
        if name != "orca":
            continue
        stem = input_stem(args)
        if stem is not None:
            by_cwd.setdefault(cwd, set()).add(stem)
    return {cwd: frozenset(stems) for cwd, stems in by_cwd.items()}


def _boot_time() -> float:
    with open(_PROC / "stat") as f:
        for line in f:
            if line.startswith("btime "):
                return float(line.split()[1])
    raise OSError("no btime in /proc/stat")


def _orca_entries() -> list[tuple[str, Path, float, list[str]]]:
    return _proc_entries() if _PROC.is_dir() else _psutil_entries()


def _proc_entries() -> list[tuple[str, Path, float, list[str]]]:
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
            if not is_orca_process(name):
                continue
            cwd = Path(os.readlink(f"{base}/cwd"))
            with open(f"{base}/stat") as f:
                stat = f.read()
            try:
                args = open(f"{base}/cmdline", "rb").read().decode("utf-8", "replace").split("\0")
            except OSError:
                args = []
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        # Field 22 is the start time in clock ticks since boot. The command
        # name (field 2) is parenthesised and may itself hold spaces or ')',
        # so count fields from the LAST ')'.
        rest = stat[stat.rfind(")") + 2 :].split()
        started = boot + int(rest[19]) / ticks
        found.append((name, _resolve(cwd), started, args))
    return found


def _psutil_entries() -> list[tuple[str, Path, float, list[str]]]:
    try:
        import psutil
    except ImportError:
        return []
    found = []
    for proc in psutil.process_iter(["name", "cwd", "create_time", "cmdline"]):
        try:
            name = proc.info["name"] or ""
            cwd = proc.info["cwd"]
            if not is_orca_process(name) or cwd is None:
                continue
            found.append((name, _resolve(Path(cwd)), proc.info["create_time"],
                          proc.info["cmdline"] or []))
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
    stems: frozenset | None = None  # the jobs (stems) known to run in this directory; None = cannot tell which


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
        stems = running_orca_stems()
        return {cwd: Liveness(True, self.name, since=t, stems=stems.get(cwd))
                for cwd, t in running_orca_cwds().items()}

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


# How long one `squeue` answer is reused within a process. A TUI refresh is
# every few seconds and `wait` polls; the scheduler is not asked more often
# than this whatever they do.
SLURM_CACHE_S = 30.0
SQUEUE_TIMEOUT_S = 10.0
SQUEUE_COMMAND = ("squeue", "--me", "--noheader", "--format=%i|%T|%Z|%S")
_SLURM_ALIVE = frozenset({"RUNNING", "COMPLETING", "CONFIGURING"})


class SlurmProbe(LivenessProbe):
    """SLURM's own queue, for a login node that cannot see the compute
    nodes' processes.

    A SLURM job is matched to a job directory by its working directory:
    exactly, or -- when a batch script `cd`s into a subdirectory -- as the
    ancestor of exactly ONE job directory that has not finished. Several
    candidates under one working directory are left unmatched rather than
    guessed at, and say so (`note`).

    A job directory SLURM does not mention is looked up in this machine's
    process table, since it may be running here outside the scheduler. If
    no process is found either, its liveness is unknown (None) while its
    output is fresher than QUIET_AFTER_S, and False after that."""

    name = "slurm"
    _cache: tuple[float, list] | None = None
    _dirs_cache: dict = {}
    _warned = False

    def __init__(self, quiet_after: float | None = None):
        from .status import QUIET_AFTER_S
        self.quiet_after = QUIET_AFTER_S if quiet_after is None else quiet_after
        self._ambiguous: dict[Path, str] = {}

    @classmethod
    def entries(cls) -> list[tuple[str, str, Path, float | None]]:
        """(job id, state, working directory, start) per job, cached."""
        now = time.monotonic()
        if cls._cache is not None and now - cls._cache[0] < SLURM_CACHE_S:
            return cls._cache[1]
        try:
            done = subprocess.run(SQUEUE_COMMAND, capture_output=True, text=True,
                                  timeout=SQUEUE_TIMEOUT_S)
            if done.returncode != 0:
                raise OSError(done.stderr.strip() or f"exit {done.returncode}")
            found = parse_squeue(done.stdout)
        except (OSError, subprocess.SubprocessError) as exc:
            if not cls._warned:
                print(f"orcamon: squeue failed ({exc}); liveness falls back to output age",
                      file=sys.stderr)
                cls._warned = True
            found = []
        cls._cache = (now, found)
        return found

    def snapshot(self, job_dirs=None):
        """Every SLURM job matched to a job directory. The candidates under a
        working directory are found by walking IT, not taken from the
        caller's list: `orcamon show` passes one directory, and a match that
        looked unique only because its siblings were not passed would be a
        guess."""
        live: dict[Path, Liveness] = {}
        self._ambiguous = {}
        for job_id, state, workdir, started in self.entries():
            if state == "PENDING":
                result = Liveness(False, self.name, queued=True, sched_id=job_id)
            elif state in _SLURM_ALIVE:
                result = Liveness(True, self.name, since=started, sched_id=job_id)
            else:
                continue
            candidates = self._job_dirs_under(workdir)
            if workdir in candidates or not candidates:
                live[workdir] = result
                continue
            under = [d for d in candidates if not _looks_finished(d)]
            if len(under) == 1:
                live[under[0]] = result
            elif len(under) > 1:
                for d in under:
                    self._ambiguous[d] = f"squeue workdir ambiguous: job {job_id}"
        # A machine with SLURM installed can still run ORCA directly -- the
        # login node, or a workstation whose queue is empty. `auto` picks
        # this probe wherever squeue exists, and without this a job started
        # by hand read as stopped once it had been quiet for --quiet-after,
        # with its process alive.
        stems = running_orca_stems()
        for cwd, started in running_orca_cwds().items():
            if cwd not in live:
                live[cwd] = Liveness(True, "process", since=started, stems=stems.get(cwd))
                self._ambiguous.pop(cwd, None)
        return live

    @classmethod
    def _job_dirs_under(cls, workdir: Path) -> list[Path]:
        """Job directories at or below a SLURM working directory, cached with
        the `squeue` answer they were asked for."""
        stamp = cls._cache[0] if cls._cache else None
        key = (workdir, stamp)
        if key not in cls._dirs_cache:
            from .discovery import discover
            cls._dirs_cache = {k: v for k, v in cls._dirs_cache.items() if k[1] == stamp}
            cls._dirs_cache[key] = sorted({ref.path for ref in discover(workdir)}) if workdir.is_dir() else []
        return cls._dirs_cache[key]

    def default(self, job_dir):
        note = self._ambiguous.get(_resolve(job_dir))
        age = _output_age(job_dir)
        if note is None and age is not None and age >= self.quiet_after:
            return Liveness(False, self.name)
        return Liveness(None, "mtime" if note else self.name, note=note)


def parse_squeue(text: str) -> list[tuple[str, str, Path, float | None]]:
    """`%i|%T|%Z|%S` lines -> entries; lines that do not parse are skipped."""
    found = []
    for line in text.splitlines():
        parts = line.strip().split("|")
        if len(parts) != 4 or not parts[0] or not parts[2].startswith("/"):
            continue
        job_id, state, workdir, start = (p.strip() for p in parts)
        found.append((job_id, state.upper(), _resolve(Path(workdir)), _parse_start(start)))
    return found


def _parse_start(text: str) -> float | None:
    try:
        return datetime.fromisoformat(text).timestamp()
    except ValueError:
        return None  # "N/A", "Unknown" for a pending job


def _looks_finished(job_dir: Path) -> bool:
    """Whether every output in the directory ends with ORCA's termination
    line -- cheap enough to ask of a handful of candidate directories."""
    outs = list(job_dir.glob("*.out"))
    if not outs:
        return False
    for out in outs:
        try:
            with open(out, "rb") as f:
                f.seek(max(0, out.stat().st_size - 4096))
                if b"ORCA TERMINATED NORMALLY" not in f.read():
                    return False
        except OSError:
            return False
    return True


def _output_age(job_dir: Path) -> float | None:
    ages = []
    for out in job_dir.glob("*.out"):
        try:
            ages.append(time.time() - out.stat().st_mtime)
        except OSError:
            pass
    return min(ages) if ages else None


# The extension point for schedulers: a name on the command line -> a probe.
# PBS/Torque (`qstat -f -F json`, matching on PBS_O_WORKDIR) would be one more
# entry here.
PROBES: dict[str, type[LivenessProbe]] = {
    "process": ProcessProbe,
    "slurm": SlurmProbe,
    "mtime": MtimeProbe,
}

LIVENESS_MODES = ("auto", *PROBES)


def make_probe(mode: str = "auto", quiet_after: float | None = None) -> LivenessProbe:
    """`auto` is SLURM where `squeue` is on PATH, the process table elsewhere."""
    if mode == "auto":
        mode = "slurm" if shutil.which("squeue") else "process"
    try:
        cls = PROBES[mode]
    except KeyError:
        raise ValueError(f"unknown liveness mode {mode!r}; one of {', '.join(('auto', *PROBES))}") from None
    return cls(quiet_after=quiet_after) if cls is SlurmProbe else cls()


def lookup(snapshot: dict[Path, Liveness], probe: LivenessProbe, job_dir: Path,
           stem: str | None = None) -> Liveness:
    entry = snapshot.get(_resolve(job_dir))
    if entry is None:
        return probe.default(job_dir)
    # A driver named the job it runs: this directory's other jobs are not it.
    if stem is not None and entry.stems is not None and stem not in entry.stems:
        return Liveness(False, entry.source)
    return entry
