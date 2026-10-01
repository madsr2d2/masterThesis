"""What a job's status IS, and what about it needs a person.

Status answers one question -- where is the job in its life -- and answers it
from evidence, first match wins. It used to answer a second question as well,
and got it wrong: "ORCA TERMINATED NORMALLY" was shown as a green
"converged", so an optimization that hit MaxIter, then ran the frequencies it
was asked for and found nine imaginary modes, read as a success. ORCA
terminating normally says the PROGRAM finished; whether the CHEMISTRY did is
a separate list of findings, `attention()`, and a finished job can carry any
number of them.

It also called every job "crashed" whose process was gone and whose output
had no termination line. A job killed by a person, a time limit or a lost
node did not crash -- ORCA reported no error -- and on a cluster login node,
where the compute node's processes cannot be seen, every running job would
have read "crashed". That case is STOPPED now, and a liveness source that
cannot tell (`alive is None`) never produces it.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from enum import Enum

from .liveness import Liveness
from .orca_input import JobInput
from .parser import STALL_WINDOW, JobState
from .units import format_age

# How long a job whose liveness is unknown may go without writing before it is
# called QUIET rather than RUNNING. Deliberately generous: a DLPNO-CCSD(T)
# step or a large numerical Hessian can be silent for a long time and be
# perfectly healthy, and QUIET is a prompt to look, not a verdict -- it is a
# flag, never FAILED or STOPPED.
QUIET_AFTER_S = 1800.0


class Status(Enum):
    NOT_RUN = "not run"
    QUEUED = "queued"
    FAILED = "failed"
    FINISHED = "finished"
    STALLED = "stalled?"
    RUNNING = "running"
    QUIET = "quiet"
    STOPPED = "stopped"


# Statuses a job does not leave without a re-run.
TERMINAL = frozenset({Status.FAILED, Status.FINISHED, Status.STOPPED})

TS_RUN_TYPES = frozenset({"optts", "scants", "neb-ts", "zoom-neb-ts"})
OPT_RUN_TYPES = frozenset({
    "opt", "copt", "zopt", "gdiis-opt", "looseopt", "normalopt", "tightopt", "verytightopt",
})


def output_age(state: JobState, now: float | None = None) -> float | None:
    if state.mtime is None:
        return None
    return max(0.0, (time.time() if now is None else now) - state.mtime)


def compute_status(
    state: JobState,
    liveness: Liveness | None,
    now: float | None = None,
    quiet_after: float = QUIET_AFTER_S,
) -> Status:
    alive = liveness.alive if liveness is not None else False
    queued = liveness is not None and liveness.queued
    if not state.has_out:
        if queued:
            return Status.QUEUED
        return Status.RUNNING if alive else Status.NOT_RUN
    if state.crashed_marker:
        return Status.FAILED
    if state.normal_completion:
        return Status.FINISHED
    if alive:
        return Status.STALLED if state.possibly_stalled() else Status.RUNNING
    if queued:
        return Status.QUEUED  # re-submitted over an old output
    if alive is None:
        age = output_age(state, now)
        return Status.QUIET if age is not None and age >= quiet_after else Status.RUNNING
    return Status.STOPPED


@dataclass(frozen=True)
class Flag:
    code: str
    message: str

    def to_dict(self) -> dict:
        return {"code": self.code, "message": self.message}


# Flag codes are part of the JSON contract (`orcamon ... --json`), in the order
# they are reported.
FLAG_CODES = (
    "failed", "stopped", "stalled", "quiet", "opt_not_converged", "scan_incomplete",
    "irc_not_converged", "ts_hessian", "ts_imaginary", "minimum_imaginary", "qm2_errors",
)


def _plural(n: int, word: str) -> str:
    return f"{n} {word}{'' if n == 1 else 's'}"


def attention(
    state: JobState, inp: JobInput | None, status: Status, now: float | None = None,
) -> list[Flag]:
    """Everything about a job a person should look at, one flag each.

    Run types come from the INPUT; a job whose input could not be read gets
    only the flags that do not need them. There are no thresholds here beyond
    ORCA's own verdicts and the counts a TS or a minimum requires -- the tool
    is general, and what counts as a worrying frequency is the chemist's call,
    so every imaginary mode is counted and none is excused as "small"."""
    flags: list[Flag] = []
    runs = {r.lower() for r in inp.run_types} if inp is not None else set()
    is_ts = bool(runs & TS_RUN_TYPES)
    is_opt = bool(runs & OPT_RUN_TYPES)

    if status is Status.FAILED:
        first = state.crash_lines[0] if state.crash_lines else "error marker in output"
        flags.append(Flag("failed", f"ORCA error: {first}"))
    if status is Status.STOPPED:
        flags.append(Flag("stopped", "ended without a termination line"))
    if status is Status.STALLED:
        flags.append(Flag("stalled", f"RMS gradient has not improved over the last {STALL_WINDOW} cycles"))
    if status is Status.QUIET:
        age = output_age(state, now)
        flags.append(Flag("quiet", f"no output for {format_age(age or 0)}"))

    maxiter = f" ({state.max_cycles})" if state.max_cycles else ""
    if state.maxiter_scan_steps:
        steps = ", ".join(str(s) for s in state.maxiter_scan_steps)
        word = "step" if len(state.maxiter_scan_steps) == 1 else "steps"
        flags.append(Flag("opt_not_converged", f"scan {word} {steps} hit MaxIter{maxiter} without converging"))
    elif state.opt_maxiter_reached and not state.opt_converged:
        flags.append(Flag("opt_not_converged", f"optimization hit MaxIter{maxiter} without converging"))

    # While a TS search runs, its Hessian is the only evidence; once it has
    # ended, a final frequency block supersedes it -- but a search that ended
    # without one (MaxIter aborts the requested Freq) is still judged by it.
    hessian_speaks = status is not Status.FINISHED or not state.freqs_final
    # A relaxed scan that terminated normally short of its last step.
    # Failed and stopped jobs already say they ended early.
    if (status is Status.FINISHED and state.scan_total
            and (state.scan_step or 0) < state.scan_total):
        flags.append(Flag("scan_incomplete",
                          f"scan ended at step {state.scan_step or 0} of {state.scan_total}"))

    # ORCA stops an IRC direction at MaxIter without an error; the run says
    # FINISHED and only this banner tells a person the path stopped short.
    if state.irc_maxiter:
        flags.append(Flag("irc_not_converged",
                          f"IRC {' and '.join(state.irc_maxiter)} hit MaxIter{maxiter} before a minimum"))

    if is_ts and hessian_speaks and state.eigen_history:
        n = state.eigen_history[-1][1]
        if n != 1:
            flags.append(Flag("ts_hessian", f"Hessian has {_plural(n, 'negative eigenvalue')} (a TS search wants 1)"))

    # Only a FINAL frequency block speaks for the result; one computed at an
    # intermediate geometry is shown, labelled, but not flagged.
    if state.imaginary_freqs is not None and state.freqs_final:
        n = len(state.imaginary_freqs)
        if is_ts and n != 1:
            flags.append(Flag("ts_imaginary", f"{_imag(n)} (a TS wants 1)"))
        elif is_opt and not is_ts and n >= 1:
            flags.append(Flag("minimum_imaginary", f"{_imag(n)} after a minimum search"))

    if state.qm2_error_count:
        flags.append(Flag("qm2_errors", _plural(state.qm2_error_count, "QM2 error")))
    return flags


def _imag(n: int) -> str:
    return f"{n} imaginary frequenc{'y' if n == 1 else 'ies'}"
