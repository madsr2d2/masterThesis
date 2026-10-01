from __future__ import annotations

import math
import threading
import time
from itertools import islice
from pathlib import Path

from rich.text import Text
from textual import work
from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical
from textual.css.query import NoMatches
from textual.widget import Widget
from textual.widgets import DataTable, Footer, Header, RichLog, Static
from textual_plotext import PlotextPlot

from ..core import cache, vibrations
from ..core.discovery import discover, short_label
from ..core.events import events
from ..core.geometry import (
    REPRESENTATIONS, TEXT_REPRESENTATIONS, View, bonds, camera_basis, principal_axes, view_along,
)
from ..core.job import Job
from ..core.liveness import lookup, make_probe
from ..core.parser import TAIL_LINES, JobState
from ..core.report import (
    FINISHED_WITH_FLAGS_STYLE, STATUS_STYLE, build_report, cycle_label, describe_point,
    geometry_shown, render_markup, render_steps_markup, steps_rows,
)
from ..core.status import QUIET_AFTER_S, TERMINAL, Status
from ..core.units import EH_TO_KJ_PER_MOL, format_wall_time
from . import geometry_text, graphics_probe, herdr_graphics, kitty, notify

# `geometry_render` (numpy, PIL -- the `images` extra) is imported
# inside the two pixel widgets' render paths, never here, so the TUI starts
# without it.

REFRESH_SECONDS = 3.0

# How often the tree is walked again for jobs submitted since launch, and
# finished jobs re-checked for a re-run. Discovery is ~0.03 s on a local tree,
# but on Lustre or NFS an rglob is metadata-heavy, which is why this is not
# every tick. `r` forces it.
REDISCOVER_SECONDS = 30.0

# How often a scan in progress pushes what it has so far to the table. The
# first scan reads every job.out end to end -- 667 MB and ~5 s on this tree --
# and the table used to be written once, at the end, so every row read
# "not run" until the last job was parsed. An idle scan finishes well inside
# this and still writes once.
PROGRESS_APPLY_S = 0.25

# How often the scan thread writes changed jobs to the state cache. It also
# writes on exit; this bounds what a killed TUI loses.
CACHE_SAVE_S = 60.0

READING = "reading…"
READING_STYLE = "dim italic"

JOB_COL_WIDTH = 40
STATUS_COL_WIDTH = 17
CYCLE_COL_WIDTH = 7
NEG_EIG_COL_WIDTH = 9
WALL_TIME_COL_WIDTH = 10

ROTATE_STEP_DEG = 5.0

# Window over which repeated input is folded into one action. Long enough to
# swallow a held key's repeat rate, short enough that a single press still
# feels immediate -- which it is, because the first press of a burst is acted
# on at once (see Coalescer).
COALESCE_S = 0.12

# Linear size of the herdr backend's interim preview frame, as a fraction of
# the pane. At a typical pane this is ~330x315 and ~10 KB of palette PNG.
PREVIEW_SCALE = 0.35


class Coalescer:
    """Leading-edge-then-trailing debounce for repeated input.

    Acts on the first event immediately, then folds everything arriving within
    `delay` into a single trailing action. Holding an arrow key otherwise
    queued one full geometry render -- and one image transmission -- per key
    repeat, which over ssh arrives faster than it can drain; Textual's output
    queue is bounded, so a backlog there stalls the writer thread and with it
    the UI. `owner` is any Textual MessagePump (a widget or the app), whose
    set_timer keeps the trailing call on the main thread."""

    def __init__(self, owner, delay: float, action) -> None:
        self._owner = owner
        self._delay = delay
        self._action = action
        self._timer = None
        self._last_fired = 0.0

    def request(self) -> None:
        self.cancel()
        if time.monotonic() - self._last_fired >= self._delay:
            self._fire()
        else:
            self._timer = self._owner.set_timer(self._delay, self._fire)

    def _fire(self) -> None:
        self._timer = None
        self._last_fired = time.monotonic()
        self._action()

    def cancel(self) -> None:
        if self._timer is not None:
            self._timer.stop()
            self._timer = None


SELECTION_MARKER = "○"
SELECTION_COLOR = "red"
# One braille sub-dot per point: a whole-cell marker quantises a 16-row chart
# to about a tenth of its range per row, and near-equal energies merge.
POINT_MARKER = "braille"

class ConvergencePlot(PlotextPlot):
    """Energy per geometry -- one dot each -- or, for a job with only one
    geometry, its SCF iterations.

    Energies are drawn as dE in kJ/mol from the first point plotted: absolute
    energies differ in the fifth significant figure, which a terminal axis
    cannot label. A relaxed scan is drawn as its PROFILE by default, one dot
    per scan step (the step's latest geometry, against the scanned
    coordinate's value for a one-parameter scan); `a` switches to every cycle
    of every step.

    For a multi-geometry job, left/right when this widget is focused scrubs a
    ring through the plotted points and the geometry pane follows it via
    `MonitorApp.update_detail`, which reads `selected_point()`. `end` resumes
    following the newest point."""

    can_focus = True
    BINDINGS = [
        ("left", "select_prev", "Prev geom"),
        ("right", "select_next", "Next geom"),
        ("end", "select_latest", "Latest"),
        ("a", "toggle_all_cycles", "All cycles"),
    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._job: Job | None = None
        self._points: list = []
        self.selected_index: int | None = None
        # The selected GeometryPoint itself, so a selection survives the
        # bounded history shifting underneath it -- an index would silently
        # slide onto a different geometry once the deque is full.
        self._selected = None
        # True until the user deliberately scrubs backward: keeps the ring
        # (and the geometry pane) tracking the newest point as a running
        # optimization writes more cycles, rather than freezing wherever it
        # happened to be the first time a multi-point series appeared.
        self._following_latest = True
        self.all_cycles = False
        self._signature: tuple | None = ()

    def _series(self, state) -> tuple[str, list]:
        """Which points to plot, and in what mode."""
        if state.scan_points and not self.all_cycles:
            return "scan", [p for _, p in sorted(state.scan_points.items()) if p.energy is not None]
        points = [p for p in state.points if p.energy is not None]
        if len(points) > 1 or state.cycle > 0:
            return ("all" if state.scan_points else "opt"), points
        return "scf", []

    def _signature_of(self, job) -> tuple | None:
        if job is None:
            return None
        state = job.state
        last = state.points[-1] if state.points else None
        return (
            self.all_cycles,
            len(state.points),
            len(state.scan_points),
            id(last),
            None if last is None else last.energy,
            len(state.scf_iterations),
            self.selected_index,
        )

    def show_job(self, job: Job | None) -> None:
        is_new_selection = job is not self._job
        # Rebuilding tears down and replots the whole plotext figure and then
        # refreshes the widget -- a repaint of the chart pane. On a tick where
        # no new point arrived and the ring has not moved there is nothing to
        # redraw, and this ran unconditionally every 3 seconds.
        if not is_new_selection and self._signature_of(job) == self._signature:
            return
        self._job = job
        self.plt.clear_data()
        self.plt.clear_figure()
        if is_new_selection:
            self._following_latest = True
            self._selected = None
        if job is None:
            self._points = []
            self.selected_index = None
            self._signature = self._signature_of(job)
            self.refresh()
            return

        state = job.state
        mode, points = self._series(state)
        self._points = points
        if mode == "scf":
            self.selected_index = None
            if state.scf_iterations:
                xs = [i for i, _ in state.scf_iterations]
                ys = [e for _, e in state.scf_iterations]
                self.plt.title("SCF: energy per iteration")
                self.plt.xlabel("iteration")
                self.plt.ylabel("E (Eh)")
                self.plt.scatter(xs, ys, marker=POINT_MARKER)
            else:
                self.plt.title("no energies yet")
        else:
            ref = points[0].energy if points else 0.0
            ys = [(p.energy - ref) * EH_TO_KJ_PER_MOL for p in points]
            if mode == "scan":
                if state.scan_values and all(p.scan_step in state.scan_values for p in points):
                    xs = [state.scan_values[p.scan_step] for p in points]
                    self.plt.xlabel(state.scan_label or "scanned coordinate")
                else:
                    xs = [p.scan_step for p in points]
                    self.plt.xlabel("scan step")
                total = f"/{state.scan_total}" if state.scan_total else ""
                self.plt.title(f"relaxed scan: {len(points)}{total} steps  [a: all cycles]")
            elif mode == "all":
                xs = list(range(1, len(points) + 1))
                self.plt.xlabel("geometry")
                self.plt.title(f"all cycles of {len(state.scan_points)} scan steps  [a: profile]")
            else:
                xs = [p.cycle for p in points]
                self.plt.xlabel("cycle")
                self.plt.title("optimization: energy per cycle")
            self.plt.ylabel("dE (kJ/mol)")
            if points:
                self.plt.scatter(xs, ys, marker=POINT_MARKER)
            self._place_selection()
            if self.selected_index is not None and len(points) > 1:
                i = self.selected_index
                self.plt.scatter([xs[i]], [ys[i]], marker=SELECTION_MARKER, color=SELECTION_COLOR)
        self._signature = self._signature_of(job)
        self.refresh()

    def _place_selection(self) -> None:
        if not self._points:
            self.selected_index = None
            return
        if not self._following_latest and self._selected is not None:
            for i, p in enumerate(self._points):
                if p is self._selected:
                    self.selected_index = i
                    return
        # Following, or the selected point fell out of the history.
        self.selected_index = len(self._points) - 1
        self._selected = self._points[-1]

    def selected_point(self):
        """The scrubbed-to point, or None while following the newest geometry
        -- which may be one whose energy has not been printed yet, and so has
        no dot to put a ring on."""
        if self._following_latest or self.selected_index is None:
            return None
        return self._points[self.selected_index]

    def _select(self, index: int) -> None:
        self.selected_index = index
        self._selected = self._points[index]
        self._following_latest = index == len(self._points) - 1
        self.app.update_detail()

    def action_select_prev(self) -> None:
        if self.selected_index is None or len(self._points) <= 1:
            return
        self._following_latest = False
        self._select(max(0, self.selected_index - 1))

    def action_select_next(self) -> None:
        if self.selected_index is None or len(self._points) <= 1:
            return
        self._select(min(len(self._points) - 1, self.selected_index + 1))

    def action_select_latest(self) -> None:
        if self._points:
            self._select(len(self._points) - 1)

    def action_toggle_all_cycles(self) -> None:
        if self._job is None or not self._job.state.scan_points:
            return
        self.all_cycles = not self.all_cycles
        self._following_latest = True
        self.app.update_detail()


GEOMETRY_BINDINGS = [
    ("left", "rotate_left", "Rotate -"),
    ("right", "rotate_right", "Rotate +"),
    ("up", "rotate_up", "Tilt +"),
    ("down", "rotate_down", "Tilt -"),
    ("d", "toggle_distances", "Distances"),
    ("l", "toggle_labels", "Labels"),
    ("f", "toggle_fog", "Fog"),
    ("x", "toggle_see_through", "See-through"),
    ("i", "cycle_modes", "Modes"),
    ("v", "next_representation", "View"),
    ("h", "toggle_hydrogens", "Hydrogens"),
    ("p", "next_axis", "Principal"),
    ("o", "toggle_rock", "Rock"),
    ("0", "reset_view", "Reset"),
    ("[", "zoom_out", "Zoom -"),
    ("]", "zoom_in", "Zoom +"),
    ("ctrl+left", "pan_left", "Pan -"),
    ("ctrl+right", "pan_right", "Pan +"),
    ("ctrl+up", "pan_up", "Pan up"),
    ("ctrl+down", "pan_down", "Pan down"),
]
ZOOM_STEP = 1.2
ZOOM_MIN = 0.2
ZOOM_MAX = 5.0
PAN_STEP = 0.6  # angstrom, per keypress -- about half a bond length
ROCK_AMPLITUDE_DEG = 15.0
ROCK_PERIOD_S = 6.0


class RotatableGeometryImage(Widget):
    """Shared elev/azim state, key bindings, and stale-completion guarding
    for both geometry-image backends.

    @work(thread=True, exclusive=True) only stops a push that has not yet
    STARTED; a push already running in its own thread completes regardless,
    so two pushes (a rotation and the 3s periodic refresh, or two rotations
    in quick succession) can run concurrently and finish in either order --
    an older, slower one finishing after a newer one visibly reverts the
    display to a stale angle. `_next_seq`/`_is_stale`, checked immediately
    before the actual write and while holding `_write_lock` (so the check
    and the write are one atomic step relative to any other push for this
    widget), fix that without needing true thread cancellation."""

    can_focus = True
    BINDINGS = GEOMETRY_BINDINGS
    # Which representations this backend can tell apart -- text has no shading
    # or occlusion, so it offers only two (see `core.geometry`).
    SUPPORTED = REPRESENTATIONS
    # None means hydrogens are always shown; an integer means shown up to that
    # many atoms until the person says otherwise.
    HYDROGENS_BY_DEFAULT_UP_TO: int | None = None
    # Rocking frames per second. Text is cheap cells and can afford more than
    # the pixel panes, whose frames are tens of KB each.
    ROCK_FPS = 5

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.elev = 20.0
        self.azim = -60.0
        self.show_distances = False
        # Atom labels: the index numbers on the pixel render, the element
        # symbols in text mode. Off shows the bare structure.
        self.show_labels = True
        self.representation = "ball-and-stick"
        self.fog = True
        self.see_through = False
        self.zoom = 1.0
        self.pan: tuple[float, float, float] = (0.0, 0.0, 0.0)
        self._hydrogens: bool | None = None
        # Added to `azim` by `view()` only, so turning rocking off returns to
        # where the person left the molecule.
        self._rock_offset = 0.0
        self._rock_timer = None
        self._rock_t0 = 0.0
        # The imaginary mode `i` is animating: ORCA's own mode index (None
        # when off), the structures tried when the pane's own geometry does
        # not hold it, and the title mark. `_mode_s` is the phase, a sine of
        # the elapsed time.
        self.mode_index = None
        self._mode_alternates: list = []
        self.mode_label = ""
        self._mode_s = 0.0
        self._mode_t0 = 0.0
        self._mode_timer = None
        # Which principal axis `p` looks along next, and the job it belongs to
        # so a new selection restarts the cycle.
        self._axis_step = 0
        self._axis_job: Job | None = None
        self._job: Job | None = None
        self._point = None
        self._request_seq = 0
        self._write_lock = threading.Lock()
        self._bond_cache: tuple[list, list] | None = None  # (atoms, bonds)
        self._input = Coalescer(self, COALESCE_S, self._on_change)

    def render(self) -> str:
        return ""

    def view(self) -> View:
        """The shared view description every render call takes."""
        return View(
            elev=self.elev, azim=self.azim + self._rock_offset, zoom=self.zoom,
            pan=self.pan, representation=self.representation, fog=self.fog,
            show_labels=self.show_labels, show_distances=self.show_distances,
            show_hydrogens=self._show_hydrogens(), see_through=self.see_through,
        )

    def _show_hydrogens(self) -> bool:
        if self._hydrogens is not None:
            return self._hydrogens
        if self.HYDROGENS_BY_DEFAULT_UP_TO is None:
            return True
        return len(self._atoms) <= self.HYDROGENS_BY_DEFAULT_UP_TO

    def _bonds_for(self, atoms: list) -> list:
        """The bond list, cached by atom-list identity -- a geometry's bonds
        do not change as it is rotated.

        The cache holds the LIST, not `id(atoms)`. An id is an address, and
        CPython hands a freed list's address to the next list it makes, so a
        superseded geometry's id can come back on a new one and be served the
        old bonds -- out of range, if the new structure is smaller, which
        kills the render worker. Holding the list keeps its address taken.
        The tuple is read and replaced whole because the pixel panes call
        this from worker threads."""
        cached = self._bond_cache
        if cached is None or cached[0] is not atoms:
            cached = (atoms, bonds(atoms))
            self._bond_cache = cached
        return cached[1]

    def _next_seq(self) -> int:
        self._request_seq += 1
        return self._request_seq

    def _is_stale(self, seq: int) -> bool:
        return seq != self._request_seq

    @staticmethod
    def _atoms_for(job: Job | None, point) -> list:
        return geometry_shown(job, point)[0]

    def action_rotate_left(self) -> None:
        self.azim = (self.azim - ROTATE_STEP_DEG) % 360
        self._input.request()

    def action_rotate_right(self) -> None:
        self.azim = (self.azim + ROTATE_STEP_DEG) % 360
        self._input.request()

    def action_rotate_up(self) -> None:
        # No clamp, same as azim: elev just keeps turning past the pole
        # rather than stopping there -- a full tumble, not a fixed
        # "top/bottom" limit. `camera_basis` keeps it continuous through the
        # pole; it snapped the picture upside down there until 2026-10-01.
        self.elev = (self.elev + ROTATE_STEP_DEG) % 360
        self._input.request()

    def action_rotate_down(self) -> None:
        self.elev = (self.elev - ROTATE_STEP_DEG) % 360
        self._input.request()

    def action_toggle_distances(self) -> None:
        self.show_distances = not self.show_distances
        self._input.request()

    def action_toggle_labels(self) -> None:
        self.show_labels = not self.show_labels
        self._input.request()

    def action_toggle_fog(self) -> None:
        self.fog = not self.fog
        self._input.request()
        # The pane title names the fog state, so it has to be rebuilt too.
        self.app.update_detail()

    def action_toggle_see_through(self) -> None:
        self.see_through = not self.see_through
        self._input.request()
        # The pane title names the see-through state, so it has to be rebuilt
        # too.
        self.app.update_detail()

    def action_next_representation(self) -> None:
        options = self.SUPPORTED
        try:
            index = options.index(self.representation)
        except ValueError:
            index = -1
        self.representation = options[(index + 1) % len(options)]
        self._input.request()
        self.app.update_detail()

    def action_toggle_hydrogens(self) -> None:
        self._hydrogens = not self._show_hydrogens()
        self._input.request()

    def action_next_axis(self) -> None:
        """Look along the next principal axis: smallest variance first, so the
        molecule lies face-on, then the middle, then the largest, then wrap.

        Read off the atoms CURRENTLY shown -- the QM region when there is one,
        else every atom -- so `p` frames what the pane is actually drawing.
        A new job selection restarts the cycle."""
        atoms = self._atoms_for(self._job, self._point)
        qm = self._job.state.qm_atom_indices if self._job is not None else None
        if self._job is not self._axis_job:
            self._axis_step = 0
            self._axis_job = self._job
        n = len(atoms)
        multilayer = qm is not None and 0 < len(qm) < n
        region = [(x, y, z) for i, (_el, x, y, z) in enumerate(atoms) if not multilayer or i in qm]
        if len(region) < 2:
            return
        axes = principal_axes(region)      # largest variance first
        axis = axes[(2 - self._axis_step % 3) % 3][0]  # 2, 1, 0, wrap
        self.elev, self.azim = view_along(axis)
        self.pan = (0.0, 0.0, 0.0)
        self._axis_step += 1
        self._input.request()
        self.app.update_detail()

    def action_reset_view(self) -> None:
        """Back to the opening view. Representation, fog, labels and hydrogens
        are left as the person set them; only the camera and the helpers are
        reset."""
        self.elev, self.azim, self.zoom = 20.0, -60.0, 1.0
        self.pan = (0.0, 0.0, 0.0)
        self._stop_rock()
        self._axis_step = 0
        self._input.request()
        self.app.update_detail()

    def action_toggle_rock(self) -> None:
        if self._rock_timer is not None:
            self._stop_rock()
            self._rock_redraw()
            return
        self._rock_t0 = time.monotonic()
        self._rock_timer = self.set_interval(1.0 / self.ROCK_FPS, self._rock_tick)

    def _rock_tick(self) -> None:
        self._rock_offset = ROCK_AMPLITUDE_DEG * math.sin(
            2.0 * math.pi * (time.monotonic() - self._rock_t0) / ROCK_PERIOD_S)
        self._rock_redraw()

    def _stop_rock(self) -> None:
        if self._rock_timer is not None:
            self._rock_timer.stop()
            self._rock_timer = None
        self._rock_offset = 0.0

    def _rock_redraw(self) -> None:
        """The pixel backends push a full frame directly. The preview-then-
        settle path would never settle here -- every frame differs, so the
        settle timer would be reset forever and only previews would show."""
        self._push(self._next_seq(), quality="full")

    def action_cycle_modes(self) -> None:
        """`i`: off -> the most negative imaginary mode -> ... -> off.

        The whole set comes from the job's own `NORMAL MODES` block, so a job
        with no Hessian (or one whose only Hessian printed no modes) shows
        nothing and says so rather than animating something else."""
        modes = vibrations.imaginary_modes(self._job.state) if self._job is not None else []
        if not modes:
            if self._job is not None and self._job.state.frequencies is None:
                self.notify("this job computes no frequencies (no Freq or NumFreq step)")
            else:
                self.notify("no imaginary mode in this job")
            return
        indices = [mode.index for mode in modes]
        if self.mode_index is None:
            position = 0
        else:
            position = indices.index(self.mode_index) + 1 if self.mode_index in indices else 0
        if position >= len(indices):
            self._stop_modes()
            self._mode_alternates = []
        else:
            self.mode_index = indices[position]
            self._start_modes()
            self._mode_alternates = vibrations.alternate_geometries(self._job)
            self.mode_label = (f"mode {indices[position]} {modes[position].cm1:.1f} cm-1 "
                               f"({position + 1}/{len(indices)})")
        self._input.request()
        self.app.update_detail()

    def _start_modes(self) -> None:
        # A cycle from one mode straight to the next must not leave the
        # previous mode's timer running: two ticks per frame, and only the
        # newest timer is the one `_stop_modes` can stop.
        if self._mode_timer is not None:
            self._mode_timer.stop()
        self._mode_t0 = time.monotonic()
        self._mode_timer = self.set_interval(1.0 / self.ROCK_FPS, self._mode_tick)

    def _mode_tick(self) -> None:
        self._mode_s = vibrations.phase_sine(time.monotonic() - self._mode_t0)
        self._rock_redraw()

    def _stop_modes(self) -> None:
        if self._mode_timer is not None:
            self._mode_timer.stop()
            self._mode_timer = None
        self.mode_index = None
        self.mode_label = ""
        self._mode_s = 0.0

    def _apply_mode(self, atoms: list, job: Job | None) -> list:
        """The atoms to draw: displaced along the current mode's own pattern
        at the current phase, or unchanged when there is no mode, the job has
        none, or none maps onto these atoms.

        A mode drawn on a fallback structure replaces `atoms` whole; the
        renderer keeps `job.state.qm_atom_indices` for the host/guest split,
        which is the same global index set on either structure."""
        if self.mode_index is None or job is None:
            return atoms
        mode = next((m for m in vibrations.imaginary_modes(job.state)
                     if m.index == self.mode_index), None)
        if mode is None:
            return atoms
        chosen = vibrations.mode_geometry([atoms] + self._mode_alternates,
                                          job.state.qm_atom_indices, mode)
        if chosen is None:
            return atoms
        drawn, qm = chosen
        offs = vibrations.offsets(drawn, qm, mode)
        if offs is None:
            return atoms
        return vibrations.displaced(drawn, offs, self._mode_s)

    def action_zoom_in(self) -> None:
        self.zoom = min(ZOOM_MAX, self.zoom * ZOOM_STEP)
        self._input.request()

    def action_zoom_out(self) -> None:
        self.zoom = max(ZOOM_MIN, self.zoom / ZOOM_STEP)
        self._input.request()

    def action_pan_left(self) -> None:
        self._pan_by(-PAN_STEP, 0.0)

    def action_pan_right(self) -> None:
        self._pan_by(PAN_STEP, 0.0)

    def action_pan_up(self) -> None:
        self._pan_by(0.0, PAN_STEP)

    def action_pan_down(self) -> None:
        self._pan_by(0.0, -PAN_STEP)

    def _pan_by(self, dx: float, dy: float) -> None:
        # `pan` shifts the camera TARGET, and the picture moves the opposite
        # way -- so to move the molecule in the arrow's direction the step has
        # to be SUBTRACTED. Adding it (the old form, on the old -right basis)
        # sent every pan key backwards.
        right, up = camera_basis(self.elev, self.azim)
        px, py, pz = self.pan
        self.pan = (
            px - right[0] * dx - up[0] * dy,
            py - right[1] * dx - up[1] * dy,
            pz - right[2] * dx - up[2] * dy,
        )
        self._input.request()

    def _on_change(self) -> None:
        raise NotImplementedError


class KittyGeometryImage(RotatableGeometryImage):
    """Real pixel image via the raw Kitty graphics protocol, ABSOLUTE
    placement -- transmitted straight to the terminal driver's own output
    queue, positioned by moving the cursor first. Confirmed to render
    correctly over plain ssh with no herdr in the way; textual-image's
    TGPImage only implements the placeholder/diacritic placement method,
    which is broken on WezTerm (see its own source comment) -- that is why
    this exists instead of just using the library.

    Used only when NOT running under herdr; HerdrGeometryImage handles that
    case, because herdr's own graphics API draws in a layer immune to
    Textual's own repaints, and raw escape sequences injected into the
    ordinary character grid are not -- Textual can and will overwrite this
    region on any repaint nearby (focus change, row selection, the periodic
    refresh), so the placement has to be re-asserted on every tick as a
    self-heal.

    Re-asserting it does NOT mean re-sending the pixels, which is what this
    used to do: ~58 KB of base64 every 3 seconds whether or not the picture
    had changed, about 155 kbit/s of ssh traffic for a screen sitting still.
    The Kitty protocol separates transmission from placement, so the pixels go
    once and each self-heal is a ~50-byte `a=p` re-placement of an image the
    terminal already holds. RETRANSMIT_EVERY bounds the one risk in that --
    if the terminal ever drops the stored image, place-only would leave the
    pane blank -- to one minute, at 1/20th of the old traffic."""

    IMAGE_ID = 1
    RETRANSMIT_EVERY = 20  # self-heal ticks between full re-transmissions
    # The frame the pixel renderer draws by default, which the terminal
    # scales into the pane. A preview is this times PREVIEW_SCALE: ~9.5 KB
    # against ~35 KB, which over ssh is the difference a held key feels.
    FULL_PX = (900, 750)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._settle_timer = None
        self._last_atoms: list | None = None
        # Region the currently-stored image was rendered for. Placing it into
        # a different one would scale the wrong pixels into the wrong box.
        self._placed_region: tuple[int, int, int, int] | None = None
        self._heals_since_transmit = 0

    def on_resize(self) -> None:
        self._input.request()

    def on_unmount(self) -> None:
        self._stop_modes()
        self._stop_rock()
        self._input.cancel()
        if self._settle_timer is not None:
            self._settle_timer.stop()
        # Taking a sequence number under the write lock makes every push
        # still in flight stale, so none can place an image after this delete.
        with self._write_lock:
            self._next_seq()
            self._write(kitty.delete(self.IMAGE_ID))

    def show_job(self, job: Job | None, point=None) -> None:
        is_new_selection = job is not self._job or point is not self._point
        if job is not self._job:
            self._stop_modes()
        self._job = job
        self._point = point
        if is_new_selection:
            self._last_atoms = None
            self._input.request()
            return
        # Same job, same cycle: the periodic refresh. The parser installs a
        # new list object for each geometry block it accepts, so identity
        # answers "did the picture change?" -- and when it did not, all that
        # is owed is the self-heal placement.
        atoms = self._atoms_for(job, point)
        if atoms is not self._last_atoms:
            self._last_atoms = atoms
            self._push(self._next_seq(), quality="full")
            return
        self._heal()

    def _heal(self) -> None:
        region = self.content_region
        if region.width <= 0 or region.height <= 0:
            return
        self._heals_since_transmit += 1
        if (
            self._placed_region != (region.x, region.y, region.width, region.height)
            or self._heals_since_transmit >= self.RETRANSMIT_EVERY
        ):
            self._push(self._next_seq(), quality="full")
            return
        self._write(self._placement(region))

    def _placement(self, region) -> str:
        return kitty.place(self.IMAGE_ID, region.x, region.y, region.width, region.height)

    def _on_change(self) -> None:
        """Interaction: a small frame now, the full one once input settles --
        the behaviour HerdrGeometryImage has, for the same reason."""
        if self._settle_timer is not None:
            self._settle_timer.stop()
        self._push(self._next_seq(), quality="preview")
        self._settle_timer = self.set_timer(HerdrGeometryImage.SETTLE_DELAY_S, self._push_settled)

    def _push_settled(self) -> None:
        self._settle_timer = None
        self._push(self._next_seq(), quality="full")

    def _write(self, data: str) -> None:
        driver = self.app._driver
        if driver is not None:
            driver.write(data)

    @work(thread=True, exclusive=True)
    def _push(self, seq: int, quality: str = "full") -> None:
        # Bind the selection ONCE. This body runs in a worker thread while the
        # main thread can still be reassigning `_job`/`_point` underneath it,
        # so reading them more than once could pair one job's atoms with
        # another job's QM indices -- or, if the selection cleared in between,
        # dereference None and kill the worker.
        job, point = self._job, self._point
        atoms = self._atoms_for(job, point)
        if not atoms:
            with self._write_lock:
                if not self._is_stale(seq):
                    self._write(kitty.delete(self.IMAGE_ID))
            return
        region = self.content_region
        if region.width <= 0 or region.height <= 0:
            return
        atoms = self._apply_mode(atoms, job)

        from . import geometry_render

        size = self.FULL_PX
        if quality == "preview":
            size = (max(1, int(size[0] * PREVIEW_SCALE)), max(1, int(size[1] * PREVIEW_SCALE)))
        image = geometry_render.render(
            atoms, self.view(), qm_atom_indices=job.state.qm_atom_indices,
            size_px=size, bond_list=self._bonds_for(atoms),
        )
        sequence = kitty.transmit(
            geometry_render.frame_png(image), self.IMAGE_ID,
            region.x, region.y, region.width, region.height,
        )

        with self._write_lock:
            if self._is_stale(seq):
                return
            self._write(sequence)
            # `a=T` stores under IMAGE_ID as well as displaying, so later
            # self-heals can re-place these pixels instead of resending them.
            self._placed_region = (region.x, region.y, region.width, region.height)
            self._heals_since_transmit = 0


class TextGeometry(RotatableGeometryImage):
    """The molecule in braille and element symbols (`geometry_text`), for
    terminals the pixel paths cannot reach -- tmux, screen, mosh, anything
    that drops image escapes. Ordinary cells: a few KB a frame.

    The atoms are bound in `show_job`, which runs under the job's lock,
    rather than read in `render`, which Textual calls whenever it likes --
    iterating a job's history while the scan thread appends to it raises."""

    SUPPORTED = TEXT_REPRESENTATIONS
    HYDROGENS_BY_DEFAULT_UP_TO = geometry_text.SHOW_HYDROGENS_UP_TO
    ROCK_FPS = 10  # a text frame is a few KB, not tens

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._atoms: list = []
        self._qm: set | None = None

    def on_unmount(self) -> None:
        self._stop_modes()
        self._stop_rock()

    def _rock_redraw(self) -> None:
        self.refresh()

    def show_job(self, job: Job | None, point=None) -> None:
        if job is not self._job:
            self._hydrogens = None
            self._stop_modes()
        self._job, self._point = job, point
        atoms = self._atoms_for(job, point)
        qm = job.state.qm_atom_indices if job is not None else None
        if atoms is not self._atoms or qm is not self._qm:
            self._atoms, self._qm = atoms, qm
            self.refresh()

    def _on_change(self) -> None:
        self.refresh()

    def render(self):
        size = self.content_size
        atoms = self._atoms
        if not atoms:
            return ""
        atoms = self._apply_mode(atoms, self._job)
        return geometry_text.render(atoms, size.width, size.height, self.view(),
                                    self._qm, self._bonds_for(atoms))


class HerdrGeometryImage(RotatableGeometryImage):
    """Real pixel image of the latest geometry, drawn by herdr itself as a
    compositor layer over this widget's own screen region (herdr 0.9.0+'s
    pane.graphics API) -- not through Textual's character grid at all.
    MonitorApp only mounts this when herdr_graphics.available() at startup;
    otherwise KittyGeometryImage (raw Kitty protocol, no herdr) is used."""

    LAYER_ID = "geometry"

    # A full-quality push costs ~100-200ms over this project's ssh-relayed
    # herdr link -- measured as data-transfer time, not fixed overhead (a
    # preview-sized frame round-trips in single-digit ms). Rather than pick
    # one size that has to compromise between "responsive to rotate" and
    # "sharp at rest", every change fires an immediate small preview and
    # schedules one full-quality frame once input goes quiet.
    #
    # PNG cut a full-quality frame from ~638 KB to ~31 KB, which on the local
    # socket is 194 ms down to 7 ms -- enough that the preview looks
    # unnecessary. It is kept because that measurement was taken WITHOUT ssh
    # in the path, and the link this exists for cannot be measured from here.
    # A preview costs one cheap render and about 10 KB; since frames went to
    # a 256-colour palette (`geometry_render.frame_png`) a full one is ~35 KB
    # at 940x900, against 129 KB as truecolour PNG.
    SETTLE_DELAY_S = 0.35

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._settle_timer = None
        self._last_atoms: list | None = None

    def on_resize(self) -> None:
        # A font-size change is the only thing that moves the cell size, and it
        # resizes the pane too -- so this is where the cached value is due for
        # re-probing.
        herdr_graphics.cell_size(refresh=True)
        self._input.request()

    def on_unmount(self) -> None:
        self._stop_modes()
        self._stop_rock()
        self._input.cancel()
        if self._settle_timer is not None:
            self._settle_timer.stop()
        # herdr draws the layer over the PANE, not the app, so one set after
        # the app has gone stays on screen over whatever the pane shows next.
        # A push already rendering in a worker thread would do exactly that:
        # it checks staleness under this lock, so taking a fresh sequence
        # number here, under it, makes the clear the last word.
        with self._write_lock:
            self._next_seq()
            herdr_graphics.clear(self.LAYER_ID)

    def show_job(self, job: Job | None, point=None) -> None:
        # `Job` instances persist for the app's lifetime and are mutated in
        # place (see Job.refresh), so an `is` check tells a genuine selection
        # change (worth the preview flash for snappy feedback) apart from the
        # periodic refresh re-showing the SAME job -- and a `cycle` change
        # (the chart scrubbing to a different geometry) is exactly as much a
        # "new selection" as a different job, interaction-wise.
        is_new_selection = job is not self._job or point is not self._point
        if job is not self._job:
            self._stop_modes()
        self._job = job
        self._point = point
        if is_new_selection:
            self._last_atoms = None
            self._input.request()
            return
        # A running optimization writes new coordinates between ticks, so a
        # periodic refresh is still worth a redraw -- but only when the
        # geometry actually shown advanced, not unconditionally every 3s. A
        # scrubbed-to historical geometry is frozen (its cycle is fixed), so
        # this also makes the periodic tick a no-op while reviewing one.
        #
        # The parser installs a NEW list object for every geometry block it
        # accepts, so object identity is exactly the question being asked --
        # "was the geometry replaced?" -- where the old
        # (len(atoms), cycle) fingerprint was blind to fresh coordinates
        # arriving at an unchanged atom count and cycle number.
        atoms = self._atoms_for(job, point)
        if atoms is self._last_atoms:
            return
        self._last_atoms = atoms
        self._push(self._next_seq(), quality="full")

    def _on_change(self) -> None:
        if self._settle_timer is not None:
            self._settle_timer.stop()
        self._push(self._next_seq(), quality="preview")
        self._settle_timer = self.set_timer(self.SETTLE_DELAY_S, self._push_settled)

    def _push_settled(self) -> None:
        """The deferred full-quality frame, on a FRESH sequence number.

        Reusing the seq `_on_change` captured meant anything that bumped
        `_request_seq` inside the settle window -- the periodic refresh's own
        full-quality push, most easily -- made this frame stale before it was
        ever sent, so it was dropped and the pane sat at preview resolution
        until the next interaction. When this fires it IS the newest intent:
        anything the user did since would have stopped the timer."""
        self._settle_timer = None
        self._push(self._next_seq(), quality="full")

    @work(thread=True, exclusive=True)
    def _push(self, seq: int, quality: str = "full") -> None:
        # Bound once, for the reason given on KittyGeometryImage._push.
        job, point = self._job, self._point
        atoms = self._atoms_for(job, point)
        if not atoms:
            with self._write_lock:
                if not self._is_stale(seq):
                    herdr_graphics.clear(self.LAYER_ID)
            return
        region = self.content_region
        if region.width <= 0 or region.height <= 0:
            return
        cells = herdr_graphics.cell_size()
        if cells is None:
            return
        atoms = self._apply_mode(atoms, job)
        # Render straight to the size that will actually be sent, which for a
        # full-quality frame is now simply the pane's own pixel extent. It
        # used to be whatever a 480 KB raw-RGBA budget allowed -- about
        # 411x291, well under the pane -- and that budget went away with the
        # switch to PNG, which carries a full-resolution frame in ~35 KB.
        width = max(1, region.width * cells.width_px)
        height = max(1, region.height * cells.height_px)
        if quality == "preview":
            width = max(1, int(width * PREVIEW_SCALE))
            height = max(1, int(height * PREVIEW_SCALE))
        target = (width, height)
        from . import geometry_render

        image = geometry_render.render(
            atoms, self.view(), qm_atom_indices=job.state.qm_atom_indices,
            size_px=target, bond_list=self._bonds_for(atoms),
        )
        with self._write_lock:
            if self._is_stale(seq):
                return
            herdr_graphics.set_image(
                image,
                grid_cols=region.width,
                grid_rows=region.height,
                viewport_col=region.x,
                viewport_row=region.y,
                layer_id=self.LAYER_ID,
            )


def _status_style(job: Job) -> str:
    if job.status is Status.FINISHED and job.flags:
        return FINISHED_WITH_FLAGS_STYLE
    return STATUS_STYLE[job.status]


class MonitorApp(App):
    CSS = """
    Screen {
        background: $surface;
    }

    #left { width: 96; }
    #job_table {
        height: 16;
        border: round $panel-darken-2;
        border-title-color: $text-muted;
    }
    #job_table:focus {
        border: round $accent;
        border-title-color: $accent;
        border-title-style: bold;
    }
    #job_table > .datatable--header {
        text-style: bold;
        background: $panel;
    }
    #geometry {
        height: 1fr;
        border: round $panel-darken-2;
        border-title-color: $text-muted;
    }
    #geometry:focus {
        border: round $accent;
        border-title-color: $accent;
        border-title-style: bold;
    }

    #detail { width: 1fr; }

    Screen.maximized #job_table { display: none; }
    Screen.maximized #detail { display: none; }
    Screen.maximized #left { width: 1fr; }
    Screen.maximized Header { display: none; }
    Screen.maximized Footer { display: none; }

    #summary {
        height: auto;
        max-height: 16;
        padding: 0 1;
        border: round $panel-darken-2;
        border-title-color: $text-muted;
    }
    #chart {
        height: 16;
        border: round $panel-darken-2;
        border-title-color: $text-muted;
    }
    #chart:focus {
        border: round $accent;
        border-title-color: $accent;
        border-title-style: bold;
    }
    #convergence {
        height: auto;
        max-height: 10;
        padding: 0 1;
        border: round $panel-darken-2;
        border-title-color: $text-muted;
    }
    #tail {
        height: 1fr;
        border: round $panel-darken-2;
        border-title-color: $text-muted;
        scrollbar-size: 1 1;
    }
    #tail:focus {
        border: round $accent;
        border-title-color: $accent;
        border-title-style: bold;
    }
    """

    BINDINGS = [
        ("q", "quit", "Quit"),
        ("r", "refresh_now", "Refresh"),
        ("m", "toggle_maximize", "Maximize geometry"),
    ]

    def __init__(
        self, root: Path, liveness: str = "auto", quiet_after: float = QUIET_AFTER_S,
        exclude: tuple[str, ...] = (), use_cache: bool = True,
        graphics: str = "text", graphics_note: str | None = None,
        notify_mode: str = "all", on_event: str | None = None,
    ):
        super().__init__()
        self.root = root.resolve()
        self.exclude = tuple(exclude)
        self.use_cache = use_cache
        # Chosen BEFORE the app starts (see `run`): asking the terminal what
        # it supports needs it in raw mode, which Textual then owns.
        self.graphics = graphics
        self.graphics_note = graphics_note
        self.notify_mode = notify_mode
        self._hook = notify.Hook(on_event)
        # The last report per job, to diff for events. A job's first report
        # has nothing to diff against, so the first scan announces nothing.
        self._reports: dict[int, object] = {}
        self._gone: set[str] = set()
        self._last_discover = time.monotonic()
        self._force_discover = False
        # The offset each job's cache entry was written at, so a save writes
        # only the jobs that read something since.
        self._cached_offsets: dict[int, int] = {}
        self._last_cache_save = time.monotonic()
        self.probe = make_probe(liveness, quiet_after=quiet_after)
        self.quiet_after = quiet_after
        self.jobs: list[Job] = []
        self.selected_label: str | None = None
        self.maximized = False
        # What each pane is currently SHOWING. Every widget update dirties the
        # widget and a dirty widget is a repaint, which over ssh is the whole
        # cost of a tick on which nothing actually changed -- so each of these
        # exists to let an unchanged pane be left alone entirely.
        self._rendered_rows: dict[str, tuple] = {}
        self._rendered_text: dict[str, str] = {}
        self._tail_state: JobState | None = None
        self._tail_seen = 0
        self._detail = Coalescer(self, COALESCE_S, self.update_detail)
        # Held for the duration of a scan. Ticks are 3s apart but a scan can
        # outlast one on a big tree, and two of them would be mutating the
        # same JobState objects from two threads.
        self._scan_lock = threading.Lock()

    def compose(self) -> ComposeResult:
        yield Header()
        with Horizontal():
            with Vertical(id="left"):
                yield DataTable(id="job_table")
                widget = {"herdr": HerdrGeometryImage, "kitty": KittyGeometryImage}.get(
                    self.graphics, TextGeometry)
                yield widget(id="geometry")
            with Vertical(id="detail"):
                yield Static(id="summary")
                yield ConvergencePlot(id="chart")
                yield Static(id="convergence")
                # max_lines matches the parser's own tail window: the log is
                # appended to rather than cleared and rewritten, so it needs
                # its own bound to stay the same length.
                yield RichLog(
                    id="tail", wrap=False, highlight=False, markup=False,
                    max_lines=TAIL_LINES,
                )
        yield Footer()

    def on_mount(self) -> None:
        refs = discover(self.root, self.exclude)
        self.jobs = [Job(ref.path, ref.stem, self.root, label=ref.label) for ref in refs]

        table = self.query_one("#job_table", DataTable)
        table.border_title = "jobs"
        table.cursor_type = "row"
        table.zebra_stripes = True
        table.add_column("job", width=JOB_COL_WIDTH, key="job")
        table.add_column("status", width=STATUS_COL_WIDTH, key="status")
        table.add_column("cycle", width=CYCLE_COL_WIDTH, key="cycle")
        table.add_column("neg eig", width=NEG_EIG_COL_WIDTH, key="neg_eig")
        table.add_column("wall time", width=WALL_TIME_COL_WIDTH, key="wall_time")
        for job in self.jobs:
            table.add_row(
                short_label(job.label, JOB_COL_WIDTH),
                Text(READING, style=READING_STYLE),
                "-",
                "-",
                "-",
                key=job.label,
            )
        if self.jobs:
            self.selected_label = self.jobs[0].label

        self.query_one("#geometry").border_title = self._geometry_title()
        self.query_one("#summary", Static).border_title = "summary"
        self.query_one("#chart", ConvergencePlot).border_title = "convergence"
        self.query_one("#convergence", Static).border_title = "geometry steps"
        self.query_one("#tail", RichLog).border_title = "output tail"

        table.focus()
        self.refresh_all()
        self.set_interval(REFRESH_SECONDS, self.refresh_all)

    def action_refresh_now(self) -> None:
        self._force_discover = True
        self.refresh_all()

    def _geometry_title(self) -> str:
        note = f", {self.graphics_note}" if self.graphics_note else ""
        try:
            geometry = self.query_one("#geometry")
        except NoMatches:
            return f"geometry ({self.graphics}{note})"
        fog = " · fog" if geometry.fog else ""
        see_through = " · see-through" if geometry.see_through else ""
        mode = f" · {geometry.mode_label}" if geometry.mode_label else ""
        return f"geometry ({self.graphics}{note} · {geometry.representation}{fog}{see_through}{mode})"

    def action_toggle_maximize(self) -> None:
        self.maximized = not self.maximized
        self.screen.set_class(self.maximized, "maximized")

    def refresh_all(self) -> None:
        self._scan()

    @work(thread=True, exclusive=True)
    def _scan(self) -> None:
        """Read every job's new output, off the main thread.

        The first pass reads each job.out end to end -- 667 MB across this
        tree -- and it used to run on the main thread inside on_mount, so
        nothing appeared until it finished. Later passes are cheap, but only
        because nothing has usually been appended; a job writing hard, or a
        tree with big outputs, would stall the UI for as long as the parse
        took.

        The selected job is read first, so the job being looked at is ready
        before the rest, then the others by how much each has left to read --
        the tree's outputs run from 1 KB to 85 MB, so most rows fill almost at
        once and the few big ones last. Results are pushed to the table every
        PROGRESS_APPLY_S rather than only at the end."""
        if not self._scan_lock.acquire(blocking=False):
            return  # a previous scan is still going; this tick can be skipped
        try:
            if self.use_cache:
                # Before the first read of each job, resume from its cache
                # entry: a finished job then costs a stat, not a parse. Done
                # before ordering, so the order sees what is really left.
                for job in self.jobs:
                    if not job.parsed and id(job) not in self._cached_offsets:
                        restored = cache.restore(job.state)
                        self._cached_offsets[id(job)] = job.state.offset if restored else -1
            full = self._force_discover or time.monotonic() - self._last_discover >= REDISCOVER_SECONDS
            if full:
                self._force_discover = False
                self._rediscover()
            snapshot = self.probe.snapshot([job.state.path for job in self.jobs])
            selected = self.selected_job()
            order = sorted(self.jobs, key=lambda j: (j is not selected, j.pending_bytes()))
            last_apply = time.monotonic()
            raised = []
            for job in order:
                if job.label in self._gone:
                    continue
                # A finished, failed or stopped job only changes if it is
                # re-run, and `update_job` notices a replaced file whenever it
                # is next read -- so between rediscoveries it is not stat'ed.
                if not full and job.parsed and job.status in TERMINAL:
                    continue
                job.refresh(lookup(snapshot, self.probe, job.state.path, job.state.stem),
                            quiet_after=self.quiet_after)
                report = build_report(job)
                raised += events(self._reports.get(id(job)), report)
                self._reports[id(job)] = report
                if time.monotonic() - last_apply >= PROGRESS_APPLY_S:
                    self.call_from_thread(self._apply_scan)
                    last_apply = time.monotonic()
            if self.use_cache and time.monotonic() - self._last_cache_save >= CACHE_SAVE_S:
                self.save_cache()
        finally:
            self._scan_lock.release()
        self.call_from_thread(self._apply_scan)
        if raised:
            self.call_from_thread(self._announce, raised)

    def _rediscover(self) -> None:
        """Add jobs that appeared since launch; mark those whose directory
        went as gone rather than pulling a row out from under the cursor.
        Runs in the scan thread; the job list and the table change on the
        main thread."""
        self._last_discover = time.monotonic()
        refs = discover(self.root, self.exclude)
        known = {job.label for job in self.jobs}
        current = {ref.label for ref in refs}
        new = [Job(ref.path, ref.stem, self.root, label=ref.label) for ref in refs if ref.label not in known]
        gone = known - current
        if new or gone != self._gone:
            self.call_from_thread(self._apply_discovery, new, gone)

    def _apply_discovery(self, new: list, gone: set) -> None:
        table = self.query_one("#job_table", DataTable)
        for job in new:
            self.jobs.append(job)
            table.add_row(short_label(job.label, JOB_COL_WIDTH), Text(READING, style=READING_STYLE),
                          "-", "-", "-", key=job.label)
        for label in gone - self._gone:
            self._rendered_rows.pop(label, None)
            table.update_cell(label, "status", Text("gone", style="dim"), update_width=False)
        self._gone = set(gone)
        for label in {job.label for job in self.jobs} - gone:
            if label in self._rendered_rows and self._rendered_rows[label][0] == "gone":
                self._rendered_rows.pop(label)

    def _announce(self, raised: list) -> None:
        """Tell the person: the terminal's own notification (bell, OSC 9/777
        -- through ssh, and through tmux with passthrough on), the
        `--on-event` hook, and a toast in the app itself."""
        driver = self._driver
        for event in raised:
            seq = notify.sequences(event, self.notify_mode)
            if seq and driver is not None:
                driver.write(seq)
            self._hook(event)
            severity = "error" if event.status in ("failed", "stopped") else "warning" if event.kind == "flag" else "information"
            self.notify(event.text, severity=severity, timeout=8)
    def save_cache(self, lock_timeout: float | None = None) -> None:
        """Write every job whose output was read since its last save. From
        the scan thread (which holds no job lock between jobs), or at exit
        with `lock_timeout` so a scan still unwinding cannot hang the quit."""
        for job in self.jobs:
            if not job.parsed or job.state.offset == self._cached_offsets.get(id(job)):
                continue
            if not job.lock.acquire(timeout=-1 if lock_timeout is None else lock_timeout):
                continue
            try:
                if cache.save(job.state):
                    self._cached_offsets[id(job)] = job.state.offset
            finally:
                job.lock.release()
        self._last_cache_save = time.monotonic()

    def _apply_scan(self) -> None:
        table = self.query_one("#job_table", DataTable)
        for job in self.jobs:
            if not job.parsed or job.label in self._gone:
                continue  # still showing READING from add_row, or gone
            neg_eig = "-"
            if job.state.eigen_history:
                neg_eig = str(job.state.eigen_history[-1][1])
            # A finished job with flags is marked, not shown as a plain
            # success: the whole reason status was rewritten.
            cells = (
                job.status.value + (" !" if job.flags else ""),
                cycle_label(job.state),
                neg_eig,
                format_wall_time(job.wall_time_s),
            )
            # Every one of these was rewritten every tick regardless of whether
            # anything had moved -- on this tree that is 15 jobs x 4 cells with,
            # typically, zero jobs whose data had actually changed, and each
            # update_cell repaints the whole table.
            if self._rendered_rows.get(job.label) == cells:
                continue
            self._rendered_rows[job.label] = cells
            status, cycle, neg, wall = cells
            table.update_cell(
                job.label,
                "status",
                Text(status, style=_status_style(job)),
                update_width=False,
            )
            table.update_cell(job.label, "cycle", cycle)
            table.update_cell(job.label, "neg_eig", neg)
            table.update_cell(job.label, "wall_time", wall)
        self.update_detail()

    def _show_text(self, widget: Static, text: str) -> None:
        """Update a Static only when its content actually differs."""
        if self._rendered_text.get(widget.id) == text:
            return
        self._rendered_text[widget.id] = text
        widget.update(text)

    def _show_tail(self, tail: RichLog, job: Job | None) -> None:
        """Bring the output tail up to date with the fewest writes.

        This used to clear the log and rewrite all 300 lines on every tick and
        every cursor move. Now a genuine selection change rewrites it and
        everything else appends only the lines the parser has read since --
        usually none at all, and on a running job a handful."""
        state = job.state if job is not None else None
        if state is not self._tail_state:
            self._tail_state = state
            self._tail_seen = state.lines_seen if state is not None else 0
            tail.clear()
            if state is not None:
                for line in state.tail:
                    tail.write(line)
            return
        if state is None:
            return
        new_lines = state.lines_seen - self._tail_seen
        if new_lines == 0:
            return
        self._tail_seen = state.lines_seen
        if new_lines < 0 or new_lines >= len(state.tail):
            # A re-run reset the counter, or more arrived than the deque holds:
            # nothing to append onto, so redraw the window we have.
            tail.clear()
            for line in state.tail:
                tail.write(line)
            return
        for line in islice(state.tail, len(state.tail) - new_lines, None):
            tail.write(line)

    def selected_job(self) -> Job | None:
        for job in self.jobs:
            if job.label == self.selected_label:
                return job
        return None

    def on_data_table_row_highlighted(self, event: DataTable.RowHighlighted) -> None:
        if event.row_key.value is not None:
            self.selected_label = event.row_key.value
        # Coalesced: scrolling the job list with a held arrow key would
        # otherwise rebuild the summary, the convergence table, the whole
        # output tail and the chart once per key repeat.
        self._detail.request()

    def update_detail(self) -> None:
        job = self.selected_job()
        if job is None:
            self._render_detail(None)
            return
        if not job.parsed:
            # First read still going: show that, and leave the half-built
            # state alone. The scan's next apply comes back here.
            self._render_detail(None, placeholder=f"[b]{job.label}[/b]\n[dim]{READING}[/dim]")
            return
        # A tick's re-read of this job is in progress: skip rather than wait
        # (the UI thread must never block on the scan). The scan re-renders
        # when it finishes.
        if not job.lock.acquire(blocking=False):
            return
        try:
            self._render_detail(job)
        finally:
            job.lock.release()

    def _render_detail(self, job: Job | None, placeholder: str = "No job selected") -> None:
        summary = self.query_one("#summary", Static)
        chart = self.query_one("#chart", ConvergencePlot)
        geometry = self.query_one("#geometry")
        convergence = self.query_one("#convergence", Static)
        tail = self.query_one("#tail", RichLog)

        chart.show_job(job)
        point = chart.selected_point()
        geometry.show_job(job, point)
        self._title_geometry(geometry, job, point)

        if job is None:
            self._show_text(summary, placeholder)
            self._show_text(convergence, "")
            self._show_tail(tail, None)
            return

        self._show_text(summary, render_markup(build_report(job)))
        self._show_text(convergence, render_steps_markup(steps_rows(job.state)))
        self._show_tail(tail, job)

    def _title_geometry(self, geometry, job: Job | None, point) -> None:
        """Say WHICH geometry the pane is drawing -- and, when the requested
        one had no coordinates printed, that it is drawing another."""
        _atoms, shown = geometry_shown(job, point)
        title = self._geometry_title()
        where = describe_point(shown)
        if point is not None and shown is not point:
            asked = describe_point(point)
            where = f"{asked}: not printed, showing {where or 'none'}"
        if where:
            title += f" · {where}"
        if geometry.border_title != title:
            geometry.border_title = title


def run(root: Path, args=None) -> None:
    """Watch every ORCA job (every `<stem>.inp`) under `root`. `args` is the
    parsed `orcamon tui` command line, when there is one."""
    graphics, note = graphics_probe.choose(getattr(args, "graphics", "auto"))
    app = MonitorApp(
        root,
        liveness=getattr(args, "liveness", "auto"),
        quiet_after=getattr(args, "quiet_after", QUIET_AFTER_S),
        exclude=tuple(getattr(args, "exclude", ()) or ()),
        use_cache=not getattr(args, "no_cache", False),
        graphics=graphics, graphics_note=note,
        notify_mode=getattr(args, "notify", "all"),
        on_event=getattr(args, "on_event", None),
    )
    try:
        app.run()
    finally:
        if app.use_cache:
            app.save_cache(lock_timeout=1.0)
