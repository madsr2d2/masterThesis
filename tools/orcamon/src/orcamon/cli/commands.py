"""One function per `orcamon` subcommand. Standard library only.

Every command reads jobs the same way (`_load`): resolve the argument, take
ONE liveness snapshot, read the output (resuming from the state cache when
its entry is still valid), and build the `JobReport` the TUI would show. What each prints is
a rendering of that report or of the parsed state behind it -- never a
second computation of the same number.
"""
from __future__ import annotations

import json
import re
import signal
import sys
import time
from collections import deque
from pathlib import Path

from .. import __version__  # noqa: F401 -- `orcamon --version` reads it here
from ..core import cache
from ..core.discovery import JobRef, discover, short_label
from ..core.job import Job
from ..core.liveness import LivenessProbe, lookup, make_probe
from ..core.parser import _CRASH_RE, _QM2_ERROR_RE, HISTORY_LEN, READ_CHUNK_BYTES
from ..core.report import (
    JobReport, build_report, describe_point, render_plain, render_steps_plain, scan_profile,
    steps_rows,
)
from ..core.resolve import ResolveError, resolve
from ..core.status import TERMINAL, Status
from ..core.units import EH_TO_KJ_PER_MOL, format_age

SCHEMA = 1
EXIT_OK, EXIT_JOB_BAD, EXIT_USAGE, EXIT_TIMEOUT = 0, 1, 2, 4
BAD = frozenset({Status.FAILED.value, Status.STOPPED.value})

TAIL_CAP = 400
TAIL_GREP_WINDOW = 2000
ERRORS_SHOWN = 20
ENERGIES_DEFAULT_LAST = 50
LS_DEFAULT_COLUMNS = 120


class UsageError(Exception):
    """Exit 2 with this message on stderr."""


# --- plumbing ---------------------------------------------------------------


class Out:
    """Collects text lines and prints at most `max_lines` of them; when it
    has to cut, the last line says how many were dropped and how to see
    them."""

    def __init__(self, max_lines: int | None):
        self.max_lines = max_lines
        self.lines: list[str] = []

    def __call__(self, text: str = "") -> None:
        self.lines.extend(text.split("\n"))

    def flush(self) -> None:
        lines = self.lines
        if self.max_lines is not None and len(lines) > self.max_lines:
            dropped = len(lines) - (self.max_lines - 1)
            lines = lines[: self.max_lines - 1] + [
                f"... {dropped} more lines (raise --max-lines, or use --json)"
            ]
        if lines:
            sys.stdout.write("\n".join(lines) + "\n")
        self.lines = []


def _emit_json(doc: dict) -> None:
    sys.stdout.write(json.dumps({"schema": SCHEMA, **doc}, indent=None, separators=(",", ":")) + "\n")


def _error(message: str) -> int:
    print(f"orcamon: {message}", file=sys.stderr)
    return EXIT_USAGE


def _probe(args) -> LivenessProbe:
    return make_probe(args.liveness)


def _root(args) -> Path:
    return Path(args.root).expanduser()


def _load_ref(ref: JobRef, root: Path, probe: LivenessProbe, snapshot, args, now=None) -> Job:
    job = Job(ref.path, ref.stem, root, label=ref.label)
    use_cache = not args.no_cache
    restored = use_cache and cache.restore(job.state)
    offset = job.state.offset
    job.refresh(lookup(snapshot, probe, ref.path), now=now, quiet_after=args.quiet_after)
    # Saved only when this call read something: an unchanged finished job
    # costs one stat and one small unpickle, and writes nothing.
    if use_cache and (not restored or job.state.offset != offset):
        cache.save(job.state)
    return job


def _load(args, arg: str) -> Job:
    root = _root(args)
    ref = resolve(arg, root)
    probe = _probe(args)
    return _load_ref(ref, root, probe, probe.snapshot([ref.path]), args)


def _job_exit(job: Job) -> int:
    return EXIT_JOB_BAD if job.status.value in BAD else EXIT_OK


def _with_job(func):
    """Resolve the JOB argument (exit 2 on no match or several) and hand the
    loaded job to `func`."""
    def run(args):
        try:
            job = _load(args, args.job)
        except ResolveError as exc:
            return _error(str(exc))
        try:
            return func(args, job)
        except UsageError as exc:
            return _error(str(exc))
    run.__name__ = func.__name__
    run.__doc__ = func.__doc__
    return run


def parse_duration(text: str) -> float:
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([smhd]?)\s*", text)
    if not m:
        raise UsageError(f"bad duration {text!r}: use a number with s, m, h or d (e.g. 30m)")
    return float(m.group(1)) * {"": 1, "s": 1, "m": 60, "h": 3600, "d": 86400}[m.group(2)]


# --- tui --------------------------------------------------------------------

TUI_MISSING = "the TUI needs the tui extra: uv tool install 'orcamon[tui]'"


def add_tui_options(p) -> None:
    """The TUI's own options. Declared here so the parser (and so the
    generated reference) lists them without importing Textual."""
    p.add_argument("--graphics", choices=("auto", "herdr", "kitty", "text"), default="auto",
                   help="geometry pane: herdr/kitty pixels or text; auto asks the terminal "
                        "and picks text inside tmux or screen (default: auto)")
    p.add_argument("--notify", choices=("off", "bell", "osc", "all"), default="all",
                   help="how finished/failed/stopped/stalled jobs and new flags are announced: "
                        "terminal bell, OSC 9/777 desktop notification, or both. Inside tmux "
                        "the OSC needs `set -g allow-passthrough on` (default: all)")
    p.add_argument("--on-event", metavar="CMD",
                   help="run CMD through the shell per event, with ORCAMON_EVENT, ORCAMON_JOB, "
                        "ORCAMON_STATUS and ORCAMON_MESSAGE set (e.g. to push to a phone)")


def cmd_tui(args) -> int:
    try:
        from ..tui.app import run
    except ImportError as exc:
        if exc.name and exc.name.split(".")[0] in {"textual", "textual_plotext", "plotext", "rich"}:
            print(TUI_MISSING, file=sys.stderr)
            return EXIT_USAGE
        raise
    root = Path(args.tui_root if args.tui_root else args.root).expanduser()
    run(root, args)
    return EXIT_OK


# --- ls ---------------------------------------------------------------------


def _progress(r: JobReport) -> str:
    if r.scan_step is not None:
        total = f"/{r.scan_total}" if r.scan_total else ""
        return f"step {r.scan_step}{total}·{r.cycle or 0}"
    if r.cycle:
        return f"cyc {r.cycle}" + (f"/{r.max_cycles}" if r.max_cycles else "")
    return "-"


def _energy_cell(r: JobReport) -> str:
    if r.energy_eh is None:
        return "-"
    return f"{r.energy_eh:.6f}" + (f" {r.energy_label}" if r.energy_label else "")


_STATUS_ORDER = {s.value: i for i, s in enumerate(
    (Status.FAILED, Status.STOPPED, Status.STALLED, Status.QUIET, Status.RUNNING,
     Status.QUEUED, Status.FINISHED, Status.NOT_RUN))}


def _status_filter(text: str | None) -> set[str] | None:
    if not text:
        return None
    wanted = set()
    names = {s.value.rstrip("?"): s.value for s in Status}
    for part in text.split(","):
        part = part.strip().rstrip("?")
        if part not in names:
            raise UsageError(f"unknown status {part!r}; one of {', '.join(names)}")
        wanted.add(names[part])
    return wanted


def cmd_ls(args) -> int:
    """One line per job: status, progress, energy, imaginary count, output
    age, label, and its attention flags after a `!`."""
    try:
        statuses = _status_filter(args.status)
        since = parse_duration(args.since) if args.since else None
    except UsageError as exc:
        return _error(str(exc))
    root = _root(args)
    if not root.is_dir():
        return _error(f"no such directory: {root}")
    probe = _probe(args)
    refs = discover(root, args.exclude)
    snapshot = probe.snapshot([ref.path for ref in refs])
    now = time.time()
    reports = []
    for ref in refs:
        job = _load_ref(ref, root, probe, snapshot, args, now=now)
        r = build_report(job, now=now)
        if statuses is not None and r.status not in statuses:
            continue
        if args.attention and not r.attention:
            continue
        if since is not None and (r.last_output_age_s is None or r.last_output_age_s > since):
            continue
        reports.append(r)

    if args.sort == "age":
        reports.sort(key=lambda r: (r.last_output_age_s is None, r.last_output_age_s or 0))
    elif args.sort == "status":
        reports.sort(key=lambda r: (_STATUS_ORDER[r.status], r.label))

    code = EXIT_JOB_BAD if any(r.status in BAD for r in reports) else EXIT_OK
    if args.json:
        jobs = []
        for r in reports:
            d = r.to_dict()
            del d["crash_lines"], d["criteria"]
            jobs.append(d)
        _emit_json({"root": str(root.resolve()), "jobs": jobs})
        return code

    out = Out(args.max_lines)
    if not reports:
        out("no jobs" + (" match" if (statuses or args.attention or since) else f" under {root}"))
        out.flush()
        return code
    # IMAG counts a job's FINAL frequencies only; an intermediate Hessian's
    # modes are explained by `show`, not tallied as the result.
    rows = [(r.status, _progress(r), _energy_cell(r),
             "-" if r.imaginary_freqs is None or r.freq_cycle is not None else str(len(r.imaginary_freqs)),
             "-" if r.last_output_age_s is None else format_age(r.last_output_age_s),
             r) for r in reports]
    widths = [max(len(h), *(len(row[i]) for row in rows))
              for i, h in enumerate(("STATUS", "PROGRESS", "ENERGY/Eh", "IMAG", "LAST"))]
    try:
        columns = int(__import__("os").environ.get("COLUMNS", LS_DEFAULT_COLUMNS))
    except ValueError:
        columns = LS_DEFAULT_COLUMNS
    label_width = max(20, columns - sum(widths) - 2 * len(widths))
    head = "  ".join(h.ljust(w) for h, w in zip(("STATUS", "PROGRESS", "ENERGY/Eh", "IMAG", "LAST"), widths))
    out(f"{head}  JOB")
    for *cells, r in rows:
        line = "  ".join(c.ljust(w) for c, w in zip(cells, widths)) + "  " + short_label(r.label, label_width)
        if r.attention:
            line += "  ! " + " ".join(r.flag_codes)
        out(line)
    out.flush()
    return code


# --- show -------------------------------------------------------------------


@_with_job
def cmd_show(args, job: Job) -> int:
    report = build_report(job)
    if args.json:
        _emit_json({"job": report.to_dict()})
    else:
        out = Out(args.max_lines)
        out(render_plain(report))
        out.flush()
    return _job_exit(job)


# --- conv -------------------------------------------------------------------


@_with_job
def cmd_conv(args, job: Job) -> int:
    rows = steps_rows(job.state, max(1, args.last))
    if args.json:
        _emit_json({"job": job.label, "cycles": rows})
    else:
        out = Out(args.max_lines)
        out(render_steps_plain(rows))
        if rows:
            out("(* = criterion met)")
        out.flush()
    return _job_exit(job)


# --- energies ---------------------------------------------------------------


@_with_job
def cmd_energies(args, job: Job) -> int:
    state = job.state
    label = state.final_energy_label
    rows: list[dict] = []
    truncated = False
    if state.scan_points:
        mode = "scan"
        profile = scan_profile(state)
        if profile:
            first = profile[0].energy
            lowest = min(p.energy for p in profile)
            top = max(profile, key=lambda p: p.energy)
            for p in profile:
                rows.append({
                    "step": p.scan_step,
                    "coordinate": state.scan_values.get(p.scan_step) if state.scan_params == 1 else None,
                    "energy_eh": p.energy,
                    "label": p.energy_label,
                    "dE_first_kj_mol": (p.energy - first) * EH_TO_KJ_PER_MOL,
                    "dE_lowest_kj_mol": (p.energy - lowest) * EH_TO_KJ_PER_MOL,
                    "max": p is top and len(profile) > 1,
                })
        if args.last:
            rows = rows[-args.last:]
    else:
        mode = "opt"
        points = [p for p in state.points if p.energy is not None]
        truncated = len(state.points) >= HISTORY_LEN
        if points:
            first = points[0].energy
            for p in points:
                rows.append({
                    "cycle": p.cycle,
                    "energy_eh": p.energy,
                    "label": p.energy_label,
                    "dE_kj_mol": (p.energy - first) * EH_TO_KJ_PER_MOL,
                })
        rows = rows[-(args.last or ENERGIES_DEFAULT_LAST):]

    if args.json:
        _emit_json({"job": job.label, "mode": mode, "truncated": truncated, "rows": rows})
        return _job_exit(job)

    out = Out(args.max_lines)
    if not rows:
        out("no energies yet")
    elif mode == "scan":
        coord = state.scan_label if state.scan_params == 1 and state.scan_label else None
        head = f"{'step':>5}  " + (f"{coord[:18]:>18}  " if coord else "") + \
            f"{'energy/Eh':>18}  {'dE step1':>10}  {'dE lowest':>10}"
        out(head + "   (kJ/mol)")
        for row in rows:
            c = f"{row['coordinate']:>18.4f}  " if coord else ""
            mark = "  <- max" if row["max"] else ""
            out(f"{row['step']:>5}  {c}{row['energy_eh']:>18.9f}  {row['dE_first_kj_mol']:>+10.2f}"
                f"  {row['dE_lowest_kj_mol']:>+10.2f}{mark}")
    else:
        if truncated:
            out(f"earliest points not kept (the last {HISTORY_LEN} geometries are); "
                "dE is from the first point shown")
        out(f"{'cycle':>6}  {'energy/Eh':>18}  {'dE (kJ/mol)':>12}")
        for row in rows:
            out(f"{row['cycle']:>6}  {row['energy_eh']:>18.9f}  {row['dE_kj_mol']:>+12.2f}")
    if rows and label:
        out(f"energies are the {label} total")
    out.flush()
    return _job_exit(job)


# --- geom -------------------------------------------------------------------


def _find_point(state, step: int | None, cycle: int | None):
    """The point asked for, and the name it was asked by. Raises UsageError
    when the point was never reached or has fallen out of the history."""
    if step is None and cycle is None:
        for p in reversed(state.points):
            if p.atoms:
                return p, "latest"
        return None, "latest"
    if step is not None and cycle is None:
        p = state.scan_points.get(step)
        if p is None:
            raise UsageError(f"scan step {step}: not reached (steps {_range(sorted(state.scan_points))})")
        return p, f"scan step {step}"
    key = (step, cycle)
    for p in state.points:
        if p.key == key:
            return p, describe_point(p)
    name = f"scan step {step} · cycle {cycle}" if step is not None else f"cycle {cycle}"
    kept = [p.cycle for p in state.points if p.scan_step == step]
    raise UsageError(f"{name}: not in the kept history (cycles {_range(kept)})")


def _range(values) -> str:
    return f"{values[0]}-{values[-1]}" if values else "none"


@_with_job
def cmd_geom(args, job: Job) -> int:
    state = job.state
    point, asked = _find_point(state, args.step, args.cycle)
    if point is None:
        raise UsageError("no coordinates printed yet")
    if not point.atoms:
        history = list(state.points)
        i = next((k for k, p in enumerate(history) if p is point), len(history))
        earlier = next((p for p in reversed(history[:i]) if p.atoms), None)
        hint = f"; nearest earlier with coordinates: {describe_point(earlier)}" if earlier else ""
        raise UsageError(f"{describe_point(point) or asked}: coordinates not printed{hint}")

    atoms = point.atoms
    region = "all"
    if args.region == "qm" and state.qm_atom_indices:
        atoms = [a for i, a in enumerate(atoms) if i in state.qm_atom_indices]
        region = "qm"
    energy = ""
    if point.energy is not None:
        label = f" ({point.energy_label})" if point.energy_label else ""
        energy = f" · E = {point.energy:.9f} Eh{label}"
    where = describe_point(point)
    if args.json:
        _emit_json({
            "job": job.label, "region": region,
            "point": {"scan_step": point.scan_step, "cycle": point.cycle,
                      "energy_eh": point.energy, "energy_label": point.energy_label},
            "atoms": [[el, x, y, z] for el, x, y, z in atoms],
        })
        return _job_exit(job)
    lines = [str(len(atoms)), f"{job.label}{' · ' + where if where else ''}{energy}"
             + (" · QM1 region" if region == "qm" else "")]
    lines += [f"{el:<2} {x:15.8f} {y:15.8f} {z:15.8f}" for el, x, y, z in atoms]
    sys.stdout.write("\n".join(lines) + "\n")
    return _job_exit(job)


# --- freqs ------------------------------------------------------------------


@_with_job
def cmd_freqs(args, job: Job) -> int:
    freqs = job.state.frequencies
    if freqs is None:
        if args.json:
            _emit_json({"job": job.label, "n_modes": None, "imaginary": None, "lowest_real": None})
        else:
            print("no frequency block yet")
        return _job_exit(job)
    imaginary = [(i, f) for i, f in enumerate(freqs) if f < 0]
    real = sorted(((i, f) for i, f in enumerate(freqs) if f > 0), key=lambda m: m[1])[: max(0, args.lowest)]
    if args.json:
        _emit_json({
            "job": job.label, "n_modes": len(freqs),
            "imaginary": [{"mode": i, "cm1": f} for i, f in imaginary],
            "lowest_real": [{"mode": i, "cm1": f} for i, f in real],
        })
        return _job_exit(job)
    out = Out(args.max_lines)
    zero = sum(1 for f in freqs if f == 0)
    out(f"{len(freqs)} modes ({zero} zero, {len(imaginary)} imaginary)")
    if imaginary:
        out("imaginary:")
        for i, f in imaginary:
            out(f"  mode {i:>4}  {f:>10.2f} cm-1")
    else:
        out("imaginary: none")
    if real:
        out(f"lowest {len(real)} real:")
        for i, f in real:
            out(f"  mode {i:>4}  {f:>10.2f} cm-1")
    out.flush()
    return _job_exit(job)


# --- input ------------------------------------------------------------------


@_with_job
def cmd_input(args, job: Job) -> int:
    inp = job.input
    if inp is None:
        raise UsageError(f"no readable input: {job.inp_path}")
    layers = {name: {"charge": c, "mult": m} for name, (c, m) in inp.layers.items()}
    doc = {
        "run_types": inp.run_types, "method": inp.method, "keywords": inp.keywords,
        "charge": inp.charge, "mult": inp.mult, "multilayer": inp.multilayer, "layers": layers,
        "nprocs": inp.nprocs, "maxcore_mb": inp.maxcore_mb,
        "coords": inp.coords_kind, "coords_file": inp.coords_file,
        "basis_reported": job.state.basis,
    }
    if args.json:
        _emit_json({"job": job.label, "input": doc})
        return _job_exit(job)
    from ..core.orca_input import describe_spin
    out = Out(args.max_lines)
    out(f"{job.label}  ({job.inp_path.name})")
    out(f"run types  {' '.join(inp.run_types) or 'SP (none named)'}")
    out(f"method     {inp.method or '-'}")
    out(f"spin       {describe_spin(inp.charge, inp.mult)}" + ("   (the QM region)" if inp.multilayer else ""))
    if inp.multilayer:
        for name, (c, m) in inp.layers.items():
            out(f"  {name:<8} {describe_spin(c, m)}")
        if "total" not in inp.layers:
            out("  total    not set")
    coords = inp.coords_kind or "-"
    if inp.coords_file:
        coords += f" {inp.coords_file}"
    out(f"coords     {coords}")
    if inp.nprocs or inp.maxcore_mb:
        out(f"resources  {inp.nprocs or '?'} procs, {inp.maxcore_mb or '?'} MB/core")
    if job.state.basis:
        out(f"basis      {job.state.basis}   (as ORCA reported it)")
    out.flush()
    return _job_exit(job)


# --- errors -----------------------------------------------------------------


_ERRORS_RE = re.compile(f"{_CRASH_RE.pattern}|{_QM2_ERROR_RE.pattern}")


def scan_errors(path: Path, context: int, keep: int = ERRORS_SHOWN) -> tuple[list[dict], int, int]:
    """Every crash-marker and QM2-error line in the file, read in bounded
    chunks with a ring buffer for the lines before each match. Returns the
    first `keep` crash matches followed by QM2 matches up to `keep` in all,
    the crash count and the QM2-error count. Memory stays bounded on an
    85 MB output: one chunk, `context` lines behind, the kept matches."""
    before: deque = deque(maxlen=max(0, context))
    crashes: list[dict] = []
    qm2: list[dict] = []
    open_after: list[dict] = []
    n_crash = n_qm2 = 0
    line_no = 0

    def take(line: str) -> None:
        nonlocal n_crash, n_qm2
        for m in open_after:
            m["after"].append(line)
        open_after[:] = [m for m in open_after if len(m["after"]) < context]
        if _ERRORS_RE.search(line):
            kind = "crash" if _CRASH_RE.search(line) else "qm2"
            match = {"line": line_no, "kind": kind, "text": line, "before": list(before), "after": []}
            if kind == "crash":
                n_crash += 1
                target = crashes
            else:
                n_qm2 += 1
                target = qm2
            if len(target) < keep:
                target.append(match)
                if context:
                    open_after.append(match)
        before.append(line)

    with open(path, "rb") as f:
        carry = b""
        while True:
            block = f.read(READ_CHUNK_BYTES)
            if not block:
                break
            data = carry + block
            cut = data.rfind(b"\n")
            if cut == -1:
                carry = data
                continue
            text = data[: cut + 1].decode("utf-8", errors="replace")
            carry = data[cut + 1 :]
            lines = text.split("\n")[:-1]
            if not open_after and not _ERRORS_RE.search(text):
                # Nothing here: only the ring buffer needs the chunk's end.
                line_no += len(lines)
                before.extend(lines[-context:] if context else [])
                continue
            for line in lines:
                line_no += 1
                take(line)
        if carry:
            line_no += 1
            take(carry.decode("utf-8", errors="replace"))
    shown = crashes + qm2[: max(0, keep - len(crashes))]
    return shown, n_crash, n_qm2


@_with_job
def cmd_errors(args, job: Job) -> int:
    if not job.out_path.exists():
        raise UsageError(f"no output yet: {job.out_path}")
    matches, n_crash, n_qm2 = scan_errors(job.out_path, max(0, args.context))
    if args.json:
        _emit_json({"job": job.label, "crash_count": n_crash, "qm2_error_count": n_qm2, "matches": matches})
        return _job_exit(job)
    out = Out(args.max_lines)
    out(f"{n_crash} error line{'s' if n_crash != 1 else ''}, {n_qm2} QM2 error{'s' if n_qm2 != 1 else ''}"
        f" in {job.out_path.name}")
    for m in matches:
        out(f"--- line {m['line']} ({m['kind']})")
        for line in m["before"]:
            out(f"  {line}")
        out(f"> {m['text']}")
        for line in m["after"]:
            out(f"  {line}")
    hidden = n_crash + n_qm2 - len(matches)
    if hidden > 0:
        out(f"({hidden} further matches not shown)")
    out.flush()
    return _job_exit(job)


# --- tail -------------------------------------------------------------------


def last_lines(path: Path, n: int, block: int = 65536) -> list[str]:
    """The last `n` complete-or-not lines of a file, read backwards in blocks
    so an 85 MB output costs a few reads."""
    with open(path, "rb") as f:
        f.seek(0, 2)
        end = f.tell()
        data = b""
        pos = end
        while pos > 0 and data.count(b"\n") <= n:
            step = min(block, pos)
            pos -= step
            f.seek(pos)
            data = f.read(step) + data
    lines = data.decode("utf-8", errors="replace").split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines[-n:] if n else []


@_with_job
def cmd_tail(args, job: Job) -> int:
    if not job.out_path.exists():
        raise UsageError(f"no output yet: {job.out_path}")
    n = max(0, min(args.n, TAIL_CAP))
    if args.grep:
        try:
            rx = re.compile(args.grep)
        except re.error as exc:
            raise UsageError(f"bad --grep pattern: {exc}") from None
        lines = [ln for ln in last_lines(job.out_path, TAIL_GREP_WINDOW) if rx.search(ln)][-n:] if n else []
    else:
        lines = last_lines(job.out_path, n)
    if args.json:
        _emit_json({"job": job.label, "lines": lines})
    elif lines:
        sys.stdout.write("\n".join(lines) + "\n")
    return _job_exit(job)


# --- wait -------------------------------------------------------------------


_COND_RE = re.compile(r"^(done|finished|change|attention|(cycle|step)>=(\d+))$")


def _condition(text: str):
    """COND -> a predicate over (job, report, first report seen)."""
    m = _COND_RE.match(text.replace(" ", ""))
    if not m:
        raise UsageError(f"bad --until {text!r}: done, finished, change, attention, cycle>=N or step>=N")
    name = m.group(1)
    done = {Status.FINISHED.value, Status.FAILED.value, Status.STOPPED.value}
    if name == "done":
        return lambda r, first: r.status in done
    if name == "finished":
        return lambda r, first: r.status == Status.FINISHED.value
    if name == "change":
        return lambda r, first: (r.status != first.status
                                 or bool(set(r.flag_codes) - set(first.flag_codes)))
    if name == "attention":
        return lambda r, first: bool(r.attention)
    n = int(m.group(3))
    if m.group(2) == "cycle":
        return lambda r, first: (r.cycle or 0) >= n
    return lambda r, first: (r.scan_step or 0) >= n


def cmd_wait(args) -> int:
    try:
        met = _condition(args.until)
    except UsageError as exc:
        return _error(str(exc))
    root = _root(args)
    try:
        refs = [resolve(a, root) for a in args.job]
    except ResolveError as exc:
        return _error(str(exc))

    probe = _probe(args)
    jobs = [Job(ref.path, ref.stem, root, label=ref.label) for ref in refs]
    use_cache = not args.no_cache
    if use_cache:
        for job in jobs:
            cache.restore(job.state)
    first: dict[int, JobReport] = {}
    deadline = time.monotonic() + max(0.0, args.timeout)

    def poll() -> tuple[list[JobReport], list[int]]:
        # One liveness snapshot per poll; each Job is kept across polls, so a
        # poll reads only the bytes appended since the last.
        snapshot = probe.snapshot([job.state.path for job in jobs])
        now = time.time()
        reports = []
        for i, job in enumerate(jobs):
            job.refresh(lookup(snapshot, probe, job.state.path), now=now, quiet_after=args.quiet_after)
            r = build_report(job, now=now)
            first.setdefault(i, r)
            reports.append(r)
        return reports, [i for i, r in enumerate(reports) if met(r, first[i])]

    previous = _on_sigterm_interrupt()
    try:
        while True:
            reports, hits = poll()
            if (args.all_jobs and len(hits) == len(jobs)) or (not args.all_jobs and hits):
                return _wait_done(args, jobs, reports, hits, use_cache)
            left = deadline - time.monotonic()
            if left <= 0:
                return _wait_timeout(args, jobs, reports, use_cache)
            # The last sleep is cut short so the final poll lands on the
            # deadline rather than an interval past it.
            time.sleep(min(max(0.01, args.interval), left))
    finally:
        _restore_sigterm(previous)


def _raise_interrupt(signum, frame):
    raise KeyboardInterrupt


def _on_sigterm_interrupt():
    """SIGTERM (a tool's own timeout, `timeout(1)`) ends a wait the way
    Ctrl-C does -- exit 130, no traceback. Only the main thread may set a
    handler; elsewhere (a test's thread) this is a no-op."""
    try:
        return signal.signal(signal.SIGTERM, _raise_interrupt)
    except ValueError:
        return None


def _restore_sigterm(previous) -> None:
    if previous is not None:
        try:
            signal.signal(signal.SIGTERM, previous)
        except ValueError:
            pass


def _save_all(jobs, use_cache) -> None:
    if use_cache:
        for job in jobs:
            cache.save(job.state)


def _wait_done(args, jobs, reports, hits, use_cache) -> int:
    _save_all(jobs, use_cache)
    code = EXIT_JOB_BAD if any(reports[i].status in BAD for i in hits) else EXIT_OK
    if args.json:
        _emit_json({"until": args.until, "met": [jobs[i].label for i in hits],
                    "jobs": [reports[i].to_dict() for i in hits]})
        return code
    out = Out(args.max_lines)
    for i in hits:
        out(f"{jobs[i].label}: met --until {args.until} ({reports[i].status})")
    for i in hits:
        out("")
        out(render_plain(reports[i]))
    out.flush()
    return code


def _wait_timeout(args, jobs, reports, use_cache) -> int:
    _save_all(jobs, use_cache)
    if args.json:
        _emit_json({"until": args.until, "met": [], "timed_out": True,
                    "jobs": [r.to_dict() for r in reports]})
        return EXIT_TIMEOUT
    out = Out(args.max_lines)
    out(f"timed out after {args.timeout:g} s waiting for --until {args.until}")
    for r in reports:
        out(f"  {r.status:<9} {_progress(r):<14} {r.label}")
    out.flush()
    return EXIT_TIMEOUT


# --- skill ------------------------------------------------------------------


def add_skill_parser(add) -> None:
    """Registered by `skill.py` once it exists (Phase 7)."""
