from __future__ import annotations

from enum import Enum

from .parser import JobState


class Status(Enum):
    NOT_RUN = "not run"
    RUNNING = "running"
    POSSIBLY_STALLED = "possibly stalled"
    CONVERGED = "converged"
    CRASHED = "crashed"


def compute_status(state: JobState, is_running: bool) -> Status:
    if not state.has_out:
        return Status.NOT_RUN
    if state.crashed_marker:
        return Status.CRASHED
    if state.normal_completion:
        return Status.CONVERGED
    if is_running:
        return Status.POSSIBLY_STALLED if state.possibly_stalled() else Status.RUNNING
    return Status.CRASHED
