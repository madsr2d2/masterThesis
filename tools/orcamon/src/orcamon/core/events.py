"""What changed about a job between two reports, as events worth telling a
person about. In the core, not the TUI, so a command-line watcher could
raise the same events from the same reports."""
from __future__ import annotations

from dataclasses import dataclass

from .report import JobReport
from .status import Status

# Statuses a person wants to hear a job has reached.
NOTIFY_STATUSES = frozenset({
    Status.FINISHED.value, Status.FAILED.value, Status.STOPPED.value,
    Status.STALLED.value, Status.QUIET.value,
})


@dataclass(frozen=True)
class Event:
    kind: str     # "status" | "flag"
    job: str      # the job's label
    status: str   # its status now
    text: str     # one line for a notification


def _short(label: str) -> str:
    return label.rstrip("/").rsplit("/", 1)[-1] or label


def events(prev: JobReport | None, report: JobReport) -> list[Event]:
    """Events between `prev` and `report` for the same job. `prev` None is
    the first time the job was seen, which raises nothing: the first full
    scan would otherwise announce every job that already exists."""
    if prev is None:
        return []
    out: list[Event] = []
    name = _short(report.label)
    new_flags = [f for f in report.attention if f["code"] not in set(prev.flag_codes)]
    if report.status != prev.status and report.status in NOTIFY_STATUSES:
        detail = ""
        # Name the most important thing that came with it, if anything did.
        if report.attention:
            detail = " · " + report.attention[0]["message"]
        out.append(Event("status", report.label, report.status, f"{name}: {report.status}{detail}"))
        # Flags that arrived with the status change are in its text already.
        new_flags = new_flags[1:] if report.attention and new_flags[:1] == report.attention[:1] else new_flags
    for flag in new_flags:
        out.append(Event("flag", report.label, report.status, f"{name}: {flag['message']}"))
    return out
