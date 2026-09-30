"""One ORCA job as the monitor tracks it: its parsed output, its parsed
input, how it stands and what about it needs attention. No Textual here --
the TUI and the CLI hold the same objects."""
from __future__ import annotations

import threading
import time
from pathlib import Path

from .discovery import relative_label
from .geometry import FileGeometry, read_xyz
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
        self._input_mtime: tuple | None = None
        # The geometry to show when the log prints none (`file_geometry`),
        # and the file stamp it was read at.
        self.file_geometry: FileGeometry | None = None
        self._file_stamp: tuple | None = None
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
            self._update_file_geometry()
            self.parsed = True

    def _update_file_geometry(self) -> None:
        """A geometry for a job whose log prints no coordinates.

        ORCA's SOLVATOR, for one, only writes its cluster to a file. In
        order: the file the log announced it wrote, then the input's own
        geometry (its `*xyzfile`, or an inline `* xyz` block). Only asked
        while nothing has been printed -- a printed geometry always wins --
        and re-read only when the file changes. `source` says which it is."""
        if self.state.atoms or any(p.atoms for p in self.state.points):
            self.file_geometry, self._file_stamp = None, None
            return
        candidates = []
        if self.state.result_geometry_file:
            candidates.append((self.state.result_geometry_file, self.state.result_geometry_file))
        if self.input is not None and self.input.coords_file:
            candidates.append((self.input.coords_file, f"input geometry, {self.input.coords_file}"))
        for name, source in candidates:
            path = self.state.path / name
            try:
                st = path.stat()
            except OSError:
                continue
            stamp = (str(path), st.st_mtime, st.st_size)
            if stamp == self._file_stamp and self.file_geometry is not None:
                return
            atoms = read_xyz(path)
            if atoms:
                self.file_geometry, self._file_stamp = FileGeometry(atoms, source), stamp
                return
        if self.input is not None and self.input.coords_atoms:
            # Kept as the same object until the input changes: the pixel
            # panes decide whether to redraw by identity.
            stamp = ("inline", self._input_mtime)
            if stamp != self._file_stamp or self.file_geometry is None:
                self.file_geometry = FileGeometry(list(self.input.coords_atoms), "input geometry")
                self._file_stamp = stamp
            return
        self.file_geometry, self._file_stamp = None, None

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
            # mtime AND size: two writes inside the filesystem's timestamp
            # granularity share an mtime, and the second edit went unseen.
            st = path.stat()
            mtime = (st.st_mtime_ns, st.st_size)
        except OSError:
            self.input, self._input_mtime = None, None
            return
        if mtime != self._input_mtime:
            self.input, self._input_mtime = read_input(path), mtime
