"""What the monitor SAYS about a job, built once and rendered for each reader.

`build_report` is the only place any reported value is computed. The TUI's
summary pane, `orcamon show`, `orcamon ls` and every `--json` document are
renderings of the one `JobReport` it returns, so a person and an agent cannot
be shown different numbers for the same job.

These are MONITORING values read off the live log. They are not a record of
the job's results: a project quoting a number should read the finished job
through ORCA's own structured output instead.
"""
from __future__ import annotations

import re
import time
from dataclasses import asdict, dataclass, field

from .geometry import FileGeometry
from .job import Job
from .orca_input import describe_spin
from .parser import JobState
from .status import TS_RUN_TYPES, Status
from .units import EH_TO_KJ_PER_MOL, format_age, format_wall_time

# The JSON contract's version. `test_report_keys_are_stable` pins the key set;
# changing it means bumping this and editing that test, deliberately.
REPORT_SCHEMA = 1

# Rich markup style per status. Strings, not Rich objects: the core imports
# nothing outside the standard library, and the TUI hands these to Rich.
STATUS_STYLE = {
    Status.NOT_RUN: "dim",
    Status.QUEUED: "cyan",
    Status.FAILED: "bold red",
    Status.FINISHED: "bold green",
    Status.STALLED: "bold yellow",
    Status.RUNNING: "bold cyan",
    Status.QUIET: "yellow",
    Status.STOPPED: "bold red",
}
# A finished job with something to look at is not shown in success green.
FINISHED_WITH_FLAGS_STYLE = "bold yellow"
RED_FLAGS = frozenset({"failed", "stopped"})

CRASH_LINES_SHOWN = 3

_CRITERIA = [
    ("dE", "energy_change", "energy_tol", "energy_conv"),
    ("RMS grad", "rms_grad", "rms_grad_tol", "rms_grad_conv"),
    ("MAX grad", "max_grad", "max_grad_tol", "max_grad_conv"),
    ("RMS step", "rms_step", "rms_step_tol", "rms_step_conv"),
    ("MAX step", "max_step", "max_step_tol", "max_step_conv"),
]
CRITERIA_NAMES = tuple(name for name, *_ in _CRITERIA)
STEP_ROWS = 6


@dataclass
class JobReport:
    schema: int
    label: str
    path: str
    stem: str
    # identity -- from the input, so it exists before any output
    run_types: list[str]
    method: str
    basis: str | None
    charge: int | None
    mult: int | None
    layers: dict[str, dict]
    multilayer: bool
    n_atoms: int | None
    n_qm_atoms: int | None
    nprocs: int | None
    maxcore_mb: int | None
    # progress
    status: str
    liveness_source: str | None
    sched_id: str | None
    cycle: int | None
    max_cycles: int | None
    scan_step: int | None
    scan_total: int | None
    scan_coordinate: str | None
    wall_time_s: float | None
    last_output_age_s: float | None
    opt_converged: bool
    opt_maxiter_reached: bool
    # results so far
    criteria: list[dict] | None
    negative_eigenvalues: int | None
    imaginary_freqs: list[float] | None
    freq_cycle: int | None           # the cycle an intermediate Hessian's
                                     # frequencies were computed at; None when
                                     # they are the job's final ones
    energy_eh: float | None
    energy_label: str | None
    delta_kj_mol: float | None
    scan_max_kj_mol: float | None
    scan_max_at: float | int | None
    qm2_errors: int
    crash_lines: list[str]
    attention: list[dict]
    # Not part of the JSON: what the text renderers need that the contract
    # does not carry.
    _input_read: bool = field(default=True, repr=False)
    _liveness_note: str | None = field(default=None, repr=False)

    def to_dict(self) -> dict:
        return {k: v for k, v in asdict(self).items() if not k.startswith("_")}

    @property
    def flag_codes(self) -> list[str]:
        return [f["code"] for f in self.attention]


def scan_profile(state: JobState) -> list:
    """The latest point of every scan step that has an energy, by step."""
    return [p for _, p in sorted(state.scan_points.items()) if p.energy is not None]


def build_report(job: Job, now: float | None = None) -> JobReport:
    now = time.time() if now is None else now
    state, inp = job.state, job.input
    liveness = job.liveness

    criteria = None
    if state.convergence_history:
        step = state.convergence_history[-1]
        criteria = [
            {"name": name, "value": getattr(step, v), "tol": getattr(step, t), "ok": getattr(step, c)}
            for name, v, t, c in _CRITERIA if getattr(step, v) is not None
        ]

    delta = scan_max = scan_at = None
    if state.scan_points:
        profile = scan_profile(state)
        if len(profile) > 1:
            top = max(profile, key=lambda p: p.energy)
            scan_max = (top.energy - profile[0].energy) * EH_TO_KJ_PER_MOL
            scan_at = state.scan_values.get(top.scan_step, top.scan_step)
    else:
        energies = [p.energy for p in state.points if p.energy is not None]
        if len(energies) > 1:
            delta = (energies[-1] - energies[0]) * EH_TO_KJ_PER_MOL

    layers = {}
    if inp is not None:
        layers = {name: {"charge": c, "mult": m} for name, (c, m) in inp.layers.items()}

    return JobReport(
        schema=REPORT_SCHEMA,
        label=job.label,
        path=str(state.path),
        stem=state.stem,
        run_types=list(inp.run_types) if inp else [],
        method=inp.method if inp else "",
        basis=state.basis,
        charge=inp.charge if inp else None,
        mult=inp.mult if inp else None,
        layers=layers,
        multilayer=bool(inp and inp.multilayer),
        n_atoms=(len(state.atoms) if state.atoms
                 else len(job.file_geometry.atoms) if getattr(job, "file_geometry", None) else None),
        n_qm_atoms=len(state.qm_atom_indices) if state.qm_atom_indices else None,
        nprocs=inp.nprocs if inp else None,
        maxcore_mb=inp.maxcore_mb if inp else None,
        status=job.status.value,
        liveness_source=liveness.source if liveness else None,
        sched_id=liveness.sched_id if liveness else None,
        cycle=state.cycle or None,
        max_cycles=state.max_cycles,
        scan_step=state.scan_step,
        scan_total=state.scan_total,
        scan_coordinate=state.scan_label,
        wall_time_s=job.wall_time_s,
        last_output_age_s=None if state.mtime is None else max(0.0, now - state.mtime),
        opt_converged=state.opt_converged,
        opt_maxiter_reached=state.opt_maxiter_reached,
        criteria=criteria,
        negative_eigenvalues=state.eigen_history[-1][1] if state.eigen_history else None,
        imaginary_freqs=None if state.imaginary_freqs is None else list(state.imaginary_freqs),
        freq_cycle=state.freq_cycle,
        energy_eh=state.final_energy,
        energy_label=state.final_energy_label,
        delta_kj_mol=delta,
        scan_max_kj_mol=scan_max,
        scan_max_at=scan_at,
        qm2_errors=state.qm2_error_count,
        crash_lines=list(state.crash_lines),
        attention=[f.to_dict() for f in job.flags],
        _input_read=inp is not None,
        _liveness_note=liveness.note if liveness else None,
    )


# --- rendering -------------------------------------------------------------
#
# Each renderer walks the same list of lines; a line is a key and a run of
# (text, style) segments. Plain text drops the styles, markup wraps each
# styled segment in a Rich tag -- so the two cannot say different things,
# which `test_plain_and_markup_say_the_same_thing` holds them to.

Segment = tuple[str, "str | None"]
Line = tuple[str, list]

# Rich's own escaping rule (rich.markup.escape), restated so the core needs
# no Rich: a '[' that could open a tag gets a backslash.
_MARKUP_TAG_RE = re.compile(r"(\\*)(\[[a-z#/@][^[]*?])")


def escape_markup(text: str) -> str:
    return _MARKUP_TAG_RE.sub(lambda m: f"{m.group(1)}{m.group(1)}\\{m.group(2)}", text)


def status_style(report: JobReport) -> str:
    status = Status(report.status)
    if status is Status.FINISHED and report.attention:
        return FINISHED_WITH_FLAGS_STYLE
    return STATUS_STYLE[status]


def report_lines(report: JobReport) -> list[Line]:
    """What the job is, how far it has got, and whether anything is wrong.

    Identity comes from the INPUT, so it is there before the job has written
    a line; everything after the rule comes from the output, and a line
    appears only when it has something to say."""
    lines: list[Line] = [("label", [(report.label, "b")])]

    if not report._input_read:
        lines.append(("run", [(f"no {report.stem}.inp read", "dim")]))
    else:
        run = " ".join(report.run_types) or "SP"
        segs = [(run, "b")]
        if report.method:
            segs.append((f"  {report.method}", None))
        lines.append(("run", segs))

    if report.basis:
        # From the OUTPUT, not the input: a composite method's basis is never
        # written on the `!` line, and ORCA's own report is what it used.
        lines.append(("basis", [(f"basis  {report.basis}", None)]))

    if report._input_read:
        spin = describe_spin(report.charge, report.mult)
        segs: list[Segment] = []
        if report.multilayer:
            total = report.layers.get("total")
            if total:
                segs.append((f"system {describe_spin(total['charge'], total['mult'])}", None))
            else:
                segs += [("system ", None), ("total charge/mult not set", "dim")]
            qm = f"   QM region {spin}"
            if report.n_qm_atoms and report.n_atoms:
                qm += f" ({report.n_qm_atoms} of {report.n_atoms} atoms)"
            segs.append((qm, None))
            for name, cm in report.layers.items():
                if name != "total":
                    segs.append((f"   {name} {describe_spin(cm['charge'], cm['mult'])}", None))
        else:
            segs.append((spin + (f" · {report.n_atoms} atoms" if report.n_atoms else ""), None))
        if report.nprocs:
            mem = f", {report.maxcore_mb} MB/core" if report.maxcore_mb else ""
            segs += [("   ", None), (f"{report.nprocs} procs{mem}", "dim")]
        lines.append(("identity", segs))

    lines.append(("rule", [("─" * 40, "dim")]))

    status = Status(report.status)
    segs = [(report.status, status_style(report))]
    if report.cycle:
        segs.append((f" · cycle {report.cycle}" + (f"/{report.max_cycles}" if report.max_cycles else ""), None))
    if report.scan_step is not None:
        segs.append((f" · scan step {report.scan_step}" + (f"/{report.scan_total}" if report.scan_total else ""), None))
    segs.append((f" · {format_wall_time(report.wall_time_s)}", None))
    if report.last_output_age_s is not None and status is not Status.NOT_RUN:
        segs.append((f" · last output {format_age(report.last_output_age_s)} ago", None))
    if report.opt_converged:
        segs += [(" · ", None), ("optimization converged", "green")]
    lines.append(("status", segs))

    # Once a job has ended, how it was known to be running no longer matters,
    # and "liveness: slurm" on a finished job only raises a question.
    ended = status in (Status.FINISHED, Status.FAILED, Status.NOT_RUN)
    if not ended and report.liveness_source and (report.liveness_source != "process" or report._liveness_note):
        note = f" ({report._liveness_note})" if report._liveness_note else ""
        sched = f" · job {report.sched_id}" if report.sched_id else ""
        lines.append(("liveness", [(f"liveness: {report.liveness_source}{note}{sched}", "dim")]))

    for flag in report.attention:
        style = "bold red" if flag["code"] in RED_FLAGS else "yellow"
        lines.append(("flag", [(f"! {flag['message']}", style)]))

    if report.criteria:
        segs = [("criteria   ", None)]
        for i, c in enumerate(report.criteria):
            if i:
                segs.append(("  ", None))
            segs.append((f"{c['name']} ✓", "green") if c["ok"] else (f"{c['name']} ✗", "red"))
        lines.append(("criteria", segs))

    is_ts = any(r.lower() in TS_RUN_TYPES for r in report.run_types)
    if report.negative_eigenvalues is not None:
        n = report.negative_eigenvalues
        text = f"{n} negative eigenvalue{'s' if n != 1 else ''}"
        if is_ts:
            segs = [("Hessian    ", None), (text, "green" if n == 1 else "yellow"), (" (a TS search wants 1)", None)]
        else:
            segs = [(f"Hessian    {text}", None)]
        lines.append(("hessian", segs))

    if report.imaginary_freqs is not None:
        freqs = report.imaginary_freqs
        if freqs:
            values = ", ".join(f"{f:.1f}" for f in freqs)
            segs = [("imaginary  ", None), (str(len(freqs)), "b"), (f": {values} cm⁻¹", None)]
        else:
            segs = [("imaginary  none", None)]
        if report.freq_cycle is not None:
            segs.append((f" (Hessian at cycle {report.freq_cycle}, not the result)", "dim"))
        lines.append(("imaginary", segs))

    if report.energy_eh is not None:
        label = f" ({report.energy_label})" if report.energy_label else ""
        text = f"energy     {report.energy_eh:.6f} Eh{label}"
        if report.scan_max_kj_mol is not None:
            at = report.scan_max_at
            where = f"{at:g}" if isinstance(at, float) else f"step {at}"
            text += f" · scan max {report.scan_max_kj_mol:+.1f} kJ/mol at {where}"
        elif report.delta_kj_mol is not None:
            text += f" · {report.delta_kj_mol:+.1f} kJ/mol since first geometry"
        lines.append(("energy", [(text, None)]))

    if report.crash_lines:
        lines.append(("crash", [("crash markers:", "bold red")]))
        for cl in report.crash_lines[-CRASH_LINES_SHOWN:]:
            lines.append(("crash", [(f"  {cl}", None)]))
    return lines


def _plain(lines: list[Line]) -> str:
    return "\n".join("".join(text for text, _ in segs) for _, segs in lines)


def _markup(lines: list[Line]) -> str:
    def seg(text, style):
        text = escape_markup(text)
        return f"[{style}]{text}[/{style}]" if style else text

    return "\n".join("".join(seg(t, s) for t, s in segs) for _, segs in lines)


def render_plain(report: JobReport) -> str:
    return _plain(report_lines(report))


def render_markup(report: JobReport) -> str:
    return _markup(report_lines(report))


# --- the convergence table --------------------------------------------------


def step_label(step) -> str:
    return f"{step.scan_step}·{step.cycle}" if step.scan_step is not None else str(step.cycle)


def steps_rows(state: JobState, n: int | None = STEP_ROWS) -> list[dict]:
    """The last `n` optimization cycles (all when `n` is None), oldest first:
    per criterion its value, tolerance and verdict, and the Hessian's
    negative-eigenvalue count where one was printed for that cycle."""
    history = state.convergence_history
    start = 0 if n is None else max(0, len(history) - n)
    rows = []
    for i in range(start, len(history)):
        step = history[i]
        rows.append({
            "cycle": step.cycle,
            "scan_step": step.scan_step,
            "label": step_label(step),
            "criteria": {
                name: {"value": getattr(step, v), "tol": getattr(step, t), "ok": getattr(step, c)}
                for name, v, t, c in _CRITERIA
            },
            "neg_eig": step.neg_eig,
        })
    return rows


def _steps_lines(rows: list[dict]) -> list[Line]:
    """One row per cycle, newest last, under a header and the tolerances.

    This was five rows per cycle for eight cycles -- 41 lines in a pane with
    room for 8, which clipped from the bottom and so only ever showed the
    OLDEST cycles of the eight, never the one running."""
    if not rows:
        return [("empty", [("no geometry convergence data yet", "dim")])]
    width = 10
    with_eig = any(r["neg_eig"] is not None for r in rows)
    head = f"{'cycle':>7} " + "".join(f"{name:>{width}}" for name in CRITERIA_NAMES)
    if with_eig:
        head += f"{'neg eig':>{width}}"
    last = rows[-1]["criteria"]
    tol = f"{'tol':>7} " + "".join(
        f"{last[name]['tol']:>{width}.1e}" if last[name]["tol"] is not None else " " * width
        for name in CRITERIA_NAMES
    )
    lines: list[Line] = [("head", [(head, "b")]), ("tol", [(tol.rstrip(), "dim")])]
    for row in rows:
        segs: list[Segment] = [(f"{row['label']:>7} ", None)]
        for name in CRITERIA_NAMES:
            c = row["criteria"][name]
            if c["value"] is None:
                segs.append((" " * width, None))
                continue
            # The mark is in the TEXT, not only the colour, so plain output
            # and a colour-blind reader see which criteria are met.
            mark = "*" if c["ok"] else " "
            segs.append((f"{c['value']:>{width - 1}.1e}{mark}", "green" if c["ok"] else "red"))
        if with_eig:
            segs.append((f"{'-' if row['neg_eig'] is None else row['neg_eig']:>{width}}", None))
        lines.append(("row", segs))
    return lines


def render_steps_plain(rows: list[dict]) -> str:
    return "\n".join(line.rstrip() for line in _plain(_steps_lines(rows)).split("\n"))


def render_steps_markup(rows: list[dict]) -> str:
    return _markup(_steps_lines(rows))


# --- the geometry a point stands for ---------------------------------------


def cycle_label(state: JobState) -> str:
    if not state.cycle:
        return "-"
    if state.scan_step is not None:
        return f"{state.scan_step}·{state.cycle}"
    return str(state.cycle)


def geometry_shown(job: Job | None, point) -> tuple[list, object]:
    """The atoms to draw for `point` (a scrub selection from the chart) and
    the point they actually belong to.

    `point` None means "follow the job": its newest geometry. A point whose
    coordinates were never printed -- ORCA stops printing them on some long
    optimizations -- borrows the nearest EARLIER point's, and the second value
    says so, so the pane can label what it is really showing instead of
    passing one cycle's structure off as another's."""
    if job is None:
        return [], None
    state = job.state
    if point is None:
        for p in reversed(state.points):
            if p.atoms:
                return p.atoms, p
        if state.atoms:
            return state.atoms, None
        # Nothing printed: the file the job wrote, or its input geometry,
        # carrying its source so the pane can say what it is showing.
        fallback = getattr(job, "file_geometry", None)
        if fallback is not None:
            return fallback.atoms, fallback
        return [], None
    if point.atoms:
        return point.atoms, point
    history = list(state.points)
    try:
        start = next(i for i, p in enumerate(history) if p is point)
    except StopIteration:
        start = len(history)
    for p in reversed(history[:start]):
        if p.atoms:
            return p.atoms, p
    return [], None


def describe_point(point) -> str:
    if point is None:
        return ""
    if isinstance(point, FileGeometry):
        return point.source
    if point.scan_step is not None:
        return f"scan step {point.scan_step} · cycle {point.cycle}"
    return f"cycle {point.cycle}" if point.cycle else ""
