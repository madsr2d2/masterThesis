"""An IRC reaction path: its points in path order, built once.

A reaction path here is the sequence of geometries ORCA computed along an
IRC, each with the energy it printed for it -- the transition state between
the two directions. The energies come from the iteration ROWS in the log and
the geometries from the trajectory files the log names; D2-D6 of
`PLAN_ORCAMON_PATHS.md` hold. `reaction_path` is the one builder, cached on
the job, so the CLI and the TUI cannot be shown different lists.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .geometry import read_xyz, read_xyz_frames
from .units import KCAL_TO_KJ


@dataclass(eq=False)
class PathPoint:
    """One point of an IRC or NEB path. Stands in for a GeometryPoint in the
    panes and in `geom`, with no cycle: `label` is what it is called."""

    kind: str                 # "irc" | "neb"
    index: int                # IRC: signed (backward < 0, TS 0, forward > 0); NEB: the image
    label: str
    energy: float | None
    de_kj_mol: float          # IRC: from the TS; NEB: from image 0
    atoms: list
    monitors: dict = field(default_factory=dict)
    energy_label: str | None = None
    source: str = ""          # the file the atoms came from, "" when none
    scan_step: None = None
    cycle: int = 0
    key: None = None


@dataclass
class PathView:
    kind: str
    points: list              # PathPoint, in path order
    focus: int                # index into `points` the chart rings while following (D6)


# Trajectory file name -> its frames, keyed on (path, mtime_ns, size). The
# TUI asks for the same trajectory on every repaint, so re-reading 100 frames
# of 137 atoms per paint is what this avoids. At most 8 entries, oldest first.
_FRAME_CACHE: dict[tuple, list] = {}
_FRAME_CACHE_MAX = 8


def _frames(path: Path) -> list:
    """`read_xyz_frames` memoised on the file's stamp. A missing file gives
    [] and is cached as such, so a trajectory that is not written yet is not
    stat-ed twice for nothing."""
    try:
        st = path.stat()
        key: tuple = (str(path), st.st_mtime_ns, st.st_size)
    except OSError:
        key = (str(path), None, None)
    if key not in _FRAME_CACHE:
        _FRAME_CACHE[key] = read_xyz_frames(path)
        if len(_FRAME_CACHE) > _FRAME_CACHE_MAX:
            del _FRAME_CACHE[next(iter(_FRAME_CACHE))]
    return _FRAME_CACHE[key]


def _stamp(path: Path) -> tuple | None:
    """What a file's contents can be identified by, or None when it is gone."""
    try:
        st = path.stat()
    except OSError:
        return None
    return (st.st_mtime_ns, st.st_size)


def reaction_path(job) -> PathView | None:
    """The job's IRC path, or None when it has no IRC rows (NEB is later).

    Built once and cached on `job.path_view`, keyed by the signature of what
    the path is made from: the row counts, whether the job ended, and the
    stamp of every file an atom can come from. While that signature is
    unchanged the SAME view -- and the same PathPoint objects -- is returned,
    so a pane that decides to redraw by identity cannot redraw the same path
    twice and a rebuild cannot answer differently."""
    state = job.state
    backward = state.irc_rows.get("backward", [])
    forward = state.irc_rows.get("forward", [])
    if not backward and not forward:
        return None
    names = (
        state.irc_files.get("full"),
        state.irc_files.get("forward"),
        state.irc_files.get("backward"),
        job.input.coords_file if job.input is not None else None,
    )
    stamps = tuple(_stamp(state.path / name) if name else None for name in names)
    signature = (len(forward), len(backward), state.normal_completion,
                 state.crashed_marker, *stamps)
    cached = getattr(job, "path_view", None)
    if cached is not None and cached[0] == signature:
        return cached[1]
    points = _irc_points(job, backward, forward)
    view = PathView("irc", points, _irc_focus(state, backward, len(points)))
    job.path_view = (signature, view)
    return view


def _irc_focus(state, backward: list, n_points: int) -> int:
    """D6: the TS once the job has ended -- normal termination or an error
    line -- otherwise the newest point. ORCA runs forward first, so a
    backward row means the backward end; a job that has read only forward
    rows is following its forward end."""
    if state.normal_completion or state.crashed_marker:
        return len(backward)
    return 0 if backward else n_points - 1


def _irc_point(state, row, index: int, label: str, atoms: list, source: str) -> PathPoint:
    return PathPoint(
        kind="irc", index=index, label=label, energy=row.energy,
        de_kj_mol=row.de_kcal * KCAL_TO_KJ, atoms=atoms,
        monitors=dict(zip(state.irc_monitors, row.monitors)),
        energy_label=state.final_energy_label, source=source,
    )


def _irc_points(job, backward: list, forward: list) -> list:
    """The IRC's points in path order: the backward rows reversed, the TS,
    then the forward rows (D3), each with ORCA's own dE and monitors (D4) and
    the geometry D5 says it has. A point whose frame is not written yet has
    no atoms and is never borrowed from another point."""
    state = job.state
    full_name = state.irc_files.get("full")
    full_frames = _frames(state.path / full_name) if full_name else []

    if len(full_frames) == len(backward) + 1 + len(forward):
        # D5: the full trajectory, when it holds every point, is in path order.
        def frame(j):
            return full_frames[j][1], full_name
    else:
        fwd_name = state.irc_files.get("forward")
        bwd_name = state.irc_files.get("backward")
        fwd_frames = _frames(state.path / fwd_name) if fwd_name else []
        bwd_frames = _frames(state.path / bwd_name) if bwd_name else []
        first = fwd_frames[0] if fwd_frames else (bwd_frames[0] if bwd_frames else None)
        coords = job.input.coords_file if job.input is not None else None
        ts_atoms, ts_source = [], ""
        if coords and first is not None:
            ts = read_xyz(state.path / coords)
            # The input geometry counts only when it has the trajectories'
            # atoms: a geometry of another size is another structure.
            if ts is not None and len(ts) == len(first[1]):
                ts_atoms, ts_source = ts, coords

        def frame(j):
            """The (atoms, source) of the point at path position j."""
            if j < len(backward):
                k, frames, name = len(backward) - 1 - j, bwd_frames, bwd_name
            elif j == len(backward):
                return ts_atoms, ts_source
            else:
                k, frames, name = j - len(backward) - 1, fwd_frames, fwd_name
            if k < len(frames):
                return frames[k][1], name
            return [], ""

    points = []
    for k in reversed(range(len(backward))):
        atoms, source = frame(len(backward) - 1 - k)
        points.append(_irc_point(state, backward[k], -(k + 1), f"IRC backward {k}", atoms, source))
    ts_energy = next((p.energy for p in state.points if p.energy is not None), None)
    atoms, source = frame(len(backward))
    points.append(PathPoint("irc", 0, "IRC TS", ts_energy, 0.0, atoms,
                            {}, state.final_energy_label, source))
    for k, row in enumerate(forward):
        atoms, source = frame(len(backward) + 1 + k)
        points.append(_irc_point(state, row, k + 1, f"IRC forward {k}", atoms, source))
    return points


def path_summary(state) -> tuple[str, str] | None:
    """The one-line `show`/`ls` description of an IRC and its progress
    fraction, read from the parsed state alone -- `ls` calls it for every job,
    so it may not touch a file. None when the job has no IRC rows."""
    backward = state.irc_rows.get("backward", [])
    forward = state.irc_rows.get("forward", [])
    if not backward and not forward:
        return None
    parts = []
    for direction, rows in (("backward", backward), ("forward", forward)):
        if not rows:
            continue
        maxiter = " (MaxIter)" if direction in state.irc_maxiter else ""
        plural = "s" if len(rows) != 1 else ""
        parts.append(f"{direction} {len(rows)} point{plural}{maxiter}, "
                     f"end {rows[-1].de_kcal * KCAL_TO_KJ:+.1f} kJ/mol")
    return ("IRC        " + " · ".join(parts) + " (from the TS)",
            f"irc B{len(backward)}/F{len(forward)}")
