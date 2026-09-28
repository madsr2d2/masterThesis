"""One ORCA job as the monitor tracks it: its parsed output, its parsed
input, how it stands and what about it needs attention. No Textual here --
the TUI and the CLI hold the same objects."""
from __future__ import annotations

import threading
import time
from pathlib import Path

from .discovery import relative_label
from .liveness import Liveness
from .orca_input import JobInput, read_input
from .parser import JobState, new_state, update_job
from .status import QUIET_AFTER_S, Flag, Status, attention, compute_status


class Job:
    def __init__(self, job_dir: Path, stem: str, root: Path, label: str | None = None):
        self.state: JobState = new_state(job_dir, stem)
        self.label = label if label is not None else relative_label(root, job_dir)
        self.status: Status = Status.NOT_RUN
        self.flags: list[Flag] = []
        self.liveness: Liveness | None = None
        self.wall_time_s: float | None = None
        self.input: JobInput | None = None
        self._input_mtime: float | None = None
        # False until the first read of job.out has finished: until then the
        # state is partial, and "not run" would be a lie about the job.
        self.parsed = False
        # Held by the scan thread while it mutates `state`. The UI never waits
        # on it -- it skips a render instead (see MonitorApp.update_detail),
        # because iterating a deque another thread is appending to raises.
        self.lock = threading.Lock()

    @property
    def out_path(self) -> Path:
        return self.state.path / f"{self.state.stem}.out"

    @property
    def inp_path(self) -> Path:
        return self.state.path / f"{self.state.stem}.inp"

    def refresh(
        self, liveness: Liveness, now: float | None = None, quiet_after: float = QUIET_AFTER_S,
    ) -> None:
        """Read what the output gained since last time and re-derive the
        status and flags. `liveness` is this job's entry from the refresh's
        one probe snapshot (`liveness.lookup`)."""
        now = time.time() if now is None else now
        with self.lock:
            self._read_input()
            update_job(self.state)
            self.liveness = liveness
            self.status = compute_status(self.state, liveness, now, quiet_after)
            if self.state.wall_time_s is not None:
                self.wall_time_s = self.state.wall_time_s
            elif liveness.alive and liveness.since is not None:
                self.wall_time_s = now - liveness.since
            else:
                self.wall_time_s = None
            self.flags = attention(self.state, self.input, self.status, now)
            self.parsed = True

    def pending_bytes(self) -> int:
        """How much of job.out this job has still to read -- the order key
        for a scan. 0 when there is no output yet."""
        try:
            return max(0, self.out_path.stat().st_size - self.state.offset)
        except OSError:
            return 0

    def _read_input(self) -> None:
        """Re-read the input whenever it changes -- it is edited between
        re-runs far more often than the output is replaced."""
        path = self.inp_path
        try:
            mtime = path.stat().st_mtime
        except OSError:
            self.input, self._input_mtime = None, None
            return
        if mtime != self._input_mtime:
            self.input, self._input_mtime = read_input(path), mtime
