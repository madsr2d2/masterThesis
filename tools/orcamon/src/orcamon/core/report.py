"""What the monitor SAYS about a job, built once and rendered for each
reader. The TUI's summary and steps panes and the CLI print from here, so a
person and an agent always see the same numbers."""
from __future__ import annotations

import time

from .job import Job
from .orca_input import describe_spin
from .parser import JobState
from .status import Status
from .units import EH_TO_KJ_PER_MOL, format_age, format_wall_time

# Rich markup style per status. Strings, not Rich objects: the core imports
# nothing outside the standard library, and the TUI hands these to Rich.
STATUS_STYLE = {
    Status.NOT_RUN: "dim",
    Status.RUNNING: "bold cyan",
    Status.POSSIBLY_STALLED: "bold yellow",
    Status.CONVERGED: "bold green",
    Status.CRASHED: "bold red",
}


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
        return state.atoms, None
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
    if point.scan_step is not None:
        return f"scan step {point.scan_step} · cycle {point.cycle}"
    return f"cycle {point.cycle}" if point.cycle else ""


TS_RUN_TYPES = {"optts", "scants", "neb-ts", "zoom-neb-ts"}

_CRITERIA = [
    ("dE", "energy_change", "energy_tol", "energy_conv"),
    ("RMS grad", "rms_grad", "rms_grad_tol", "rms_grad_conv"),
    ("MAX grad", "max_grad", "max_grad_tol", "max_grad_conv"),
    ("RMS step", "rms_step", "rms_step_tol", "rms_step_conv"),
    ("MAX step", "max_step", "max_step_tol", "max_step_conv"),
]
STEP_ROWS = 6


def summary_text(job: Job) -> str:
    """What the job is, how far it has got, and whether anything is wrong.

    Identity comes from the INPUT, so it is there before the job has written
    a line; everything after the rule comes from the output, and a line
    appears only when it has something to say."""
    state, inp = job.state, job.input
    lines = [f"[b]{job.label}[/b]"]

    if inp is None:
        lines.append(f"[dim]no {state.stem}.inp read[/dim]")
    else:
        run = " ".join(inp.run_types) or "SP"
        lines.append(f"[b]{run}[/b]  {inp.method}".rstrip())

    if state.basis:
        # From the OUTPUT, not the input: a composite method's basis is never
        # written on the `!` line, and ORCA's own report is what it used.
        lines.append(f"basis  {state.basis}")

    if inp is not None:
        n_atoms = len(state.atoms) if state.atoms else None
        spin = describe_spin(inp.charge, inp.mult)
        if inp.multilayer:
            qm = f"QM region {spin}"
            if state.qm_atom_indices and n_atoms:
                qm += f" ({len(state.qm_atom_indices)} of {n_atoms} atoms)"
            total = inp.layers.get("total")
            system = f"system {describe_spin(*total)}" if total else "system [dim]total charge/mult not set[/dim]"
            parts = [system, qm]
            parts += [f"{name} {describe_spin(*cm)}" for name, cm in inp.layers.items() if name != "total"]
        else:
            parts = [spin + (f" · {n_atoms} atoms" if n_atoms else "")]
        if inp.nprocs:
            mem = f", {inp.maxcore_mb} MB/core" if inp.maxcore_mb else ""
            parts.append(f"[dim]{inp.nprocs} procs{mem}[/dim]")
        lines.append("   ".join(parts))

    lines.append("[dim]" + "─" * 40 + "[/dim]")

    style = STATUS_STYLE[job.status]
    status = [f"[{style}]{job.status.value}[/{style}]"]
    if state.cycle:
        status.append(f"cycle {state.cycle}" + (f"/{state.max_cycles}" if state.max_cycles else ""))
    if state.scan_step is not None:
        status.append(f"scan step {state.scan_step}" + (f"/{state.scan_total}" if state.scan_total else ""))
    status.append(format_wall_time(job.wall_time_s))
    if state.mtime is not None and job.status is not Status.NOT_RUN:
        status.append(f"last output {format_age(time.time() - state.mtime)} ago")
    if state.opt_converged:
        status.append("[green]optimization converged[/green]")
    lines.append(" · ".join(status))

    if state.convergence_history:
        step = state.convergence_history[-1]
        chips = []
        for name, value_attr, _tol, conv_attr in _CRITERIA:
            if getattr(step, value_attr) is None:
                continue
            ok = getattr(step, conv_attr)
            chips.append(f"[green]{name} ✓[/green]" if ok else f"[red]{name} ✗[/red]")
        lines.append("criteria   " + "  ".join(chips))

    is_ts = inp is not None and any(r.lower() in TS_RUN_TYPES for r in inp.run_types)
    if state.eigen_history:
        n = state.eigen_history[-1][1]
        text = f"{n} negative eigenvalue{'s' if n != 1 else ''}"
        if is_ts:
            colour = "green" if n == 1 else "yellow"
            text = f"[{colour}]{text}[/{colour}] (a TS search wants 1)"
        lines.append(f"Hessian    {text}")

    if state.imaginary_freqs is not None:
        freqs = state.imaginary_freqs
        if freqs:
            values = ", ".join(f"{f:.1f}" for f in freqs)
            lines.append(f"imaginary  [b]{len(freqs)}[/b]: {values} cm⁻¹")
        else:
            lines.append("imaginary  none")

    if state.final_energy is not None:
        label = f" ({state.final_energy_label})" if state.final_energy_label else ""
        text = f"energy     {state.final_energy:.6f} Eh{label}"
        energies = [p.energy for p in state.points if p.energy is not None]
        if state.scan_points:
            profile = [p for _, p in sorted(state.scan_points.items()) if p.energy is not None]
            if len(profile) > 1:
                top = max(profile, key=lambda p: p.energy)
                rise = (top.energy - profile[0].energy) * EH_TO_KJ_PER_MOL
                where = state.scan_values.get(top.scan_step)
                at = f"{where:g}" if where is not None else f"step {top.scan_step}"
                text += f" · scan max {rise:+.1f} kJ/mol at {at}"
        elif len(energies) > 1:
            change = (energies[-1] - energies[0]) * EH_TO_KJ_PER_MOL
            text += f" · {change:+.1f} kJ/mol since first geometry"
        lines.append(text)

    if state.qm2_error_count:
        lines.append(f"[yellow]QM2 errors: {state.qm2_error_count}[/yellow]")
    if state.crash_lines:
        lines.append("[bold red]crash markers:[/bold red]")
        lines.extend(f"  {cl}" for cl in state.crash_lines)
    return "\n".join(lines)


def steps_text(state: JobState) -> str:
    """The last STEP_ROWS optimization cycles, one row each, newest last.

    This was five rows per cycle for eight cycles -- 41 lines in a pane with
    room for 8, which clipped from the bottom and so only ever showed the
    OLDEST cycles of the eight, never the one running."""
    history = state.convergence_history
    if not history:
        return "[dim]no geometry convergence data yet[/dim]"
    width = 10
    head = f"{'cycle':>7} " + "".join(f"{name:>{width}}" for name, *_ in _CRITERIA)
    last = history[-1]
    tol = f"{'tol':>7} " + "".join(
        f"{getattr(last, tol_attr):>{width}.1e}" if getattr(last, tol_attr) is not None else " " * width
        for _name, _v, tol_attr, _c in _CRITERIA
    )
    rows = [f"[b]{head}[/b]", f"[dim]{tol}[/dim]"]
    for i in range(max(0, len(history) - STEP_ROWS), len(history)):
        step = history[i]
        label = f"{step.scan_step}·{step.cycle}" if step.scan_step is not None else str(step.cycle)
        cells = []
        for _name, value_attr, _tol, conv_attr in _CRITERIA:
            value = getattr(step, value_attr)
            if value is None:
                cells.append(" " * width)
                continue
            colour = "green" if getattr(step, conv_attr) else "red"
            cells.append(f"[{colour}]{value:>{width}.1e}[/{colour}]")
        rows.append(f"{label:>7} " + "".join(cells))
    return "\n".join(rows)
