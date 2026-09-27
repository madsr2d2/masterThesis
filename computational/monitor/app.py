from __future__ import annotations

import base64
import io
import threading
import time
from itertools import islice
from pathlib import Path

from rich.text import Text
from textual import work
from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widget import Widget
from textual.widgets import DataTable, Footer, Header, RichLog, Static
from textual_plotext import PlotextPlot

from . import geometry_render, herdr_graphics
from .discovery import discover_jobs, relative_label, short_label
from .orca_input import JobInput, describe_spin, read_input
from .parser import TAIL_LINES, JobState, new_state, update_job
from .procs import running_orca_cwds
from .status import Status, compute_status

REFRESH_SECONDS = 3.0

JOB_COL_WIDTH = 40
STATUS_COL_WIDTH = 17
CYCLE_COL_WIDTH = 7
NEG_EIG_COL_WIDTH = 9
WALL_TIME_COL_WIDTH = 10

STATUS_STYLE = {
    Status.NOT_RUN: "dim",
    Status.RUNNING: "bold cyan",
    Status.POSSIBLY_STALLED: "bold yellow",
    Status.CONVERGED: "bold green",
    Status.CRASHED: "bold red",
}

ROTATE_STEP_DEG = 5.0

# Window over which repeated input is folded into one action. Long enough to
# swallow a held key's repeat rate, short enough that a single press still
# feels immediate -- which it is, because the first press of a burst is acted
# on at once (see Coalescer).
COALESCE_S = 0.12

# Linear size of the herdr backend's interim preview frame, as a fraction of
# the pane. At a typical pane this is ~336x238 and ~20 KB of PNG.
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


def format_wall_time(seconds: float | None) -> str:
    if seconds is None:
        return "-"
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{s:02d}"


class Job:
    def __init__(self, job_dir: Path, stem: str, root: Path):
        self.state: JobState = new_state(job_dir, stem)
        self.label = relative_label(root, job_dir)
        self.status: Status = Status.NOT_RUN
        self.wall_time_s: float | None = None
        self.input: JobInput | None = None
        self._input_mtime: float | None = None

    def refresh(self, running_cwds: dict[Path, float]) -> None:
        self._read_input()
        update_job(self.state)
        create_time = running_cwds.get(self.state.path.resolve())
        is_running = create_time is not None
        self.status = compute_status(self.state, is_running)
        if self.state.wall_time_s is not None:
            self.wall_time_s = self.state.wall_time_s
        elif is_running:
            self.wall_time_s = time.time() - create_time
        else:
            self.wall_time_s = None

    def _read_input(self) -> None:
        """Re-read the input whenever it changes -- it is edited between
        re-runs far more often than the output is replaced."""
        path = self.state.path / f"{self.state.stem}.inp"
        try:
            mtime = path.stat().st_mtime
        except OSError:
            self.input, self._input_mtime = None, None
            return
        if mtime != self._input_mtime:
            self.input, self._input_mtime = read_input(path), mtime


def format_age(seconds: float) -> str:
    """Minute resolution on purpose: the summary is only repainted when its
    text changes, and a seconds count would change it on every tick."""
    seconds = max(0, int(seconds))
    if seconds < 60:
        return "<1 min"
    if seconds < 90 * 60:
        return f"{seconds // 60} min"
    if seconds < 36 * 3600:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds / 86400:.1f} d"


def cycle_label(state: JobState) -> str:
    if not state.cycle:
        return "-"
    if state.scan_step is not None:
        return f"{state.scan_step}·{state.cycle}"
    return str(state.cycle)


SELECTION_MARKER = "○"
SELECTION_COLOR = "red"
# One braille sub-dot per point: a whole-cell marker quantises a 16-row chart
# to about a tenth of its range per row, and near-equal energies merge.
POINT_MARKER = "braille"

EH_TO_KJ_PER_MOL = 2625.4996394799


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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.elev = 20.0
        self.azim = -60.0
        self.show_distances = False
        self.zoom = 1.0
        self.pan: tuple[float, float, float] = (0.0, 0.0, 0.0)
        self._job: Job | None = None
        self._point = None
        self._request_seq = 0
        self._write_lock = threading.Lock()
        self._input = Coalescer(self, COALESCE_S, self._on_change)

    def render(self) -> str:
        return ""

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
        # rather than stopping there (matplotlib's view_init renders that
        # fine -- it's a full tumble, not a fixed "top/bottom" limit).
        self.elev = (self.elev + ROTATE_STEP_DEG) % 360
        self._input.request()

    def action_rotate_down(self) -> None:
        self.elev = (self.elev - ROTATE_STEP_DEG) % 360
        self._input.request()

    def action_toggle_distances(self) -> None:
        self.show_distances = not self.show_distances
        self._input.request()

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
        right, up = geometry_render.camera_basis(self.elev, self.azim)
        px, py, pz = self.pan
        self.pan = (
            px + right[0] * dx + up[0] * dy,
            py + right[1] * dx + up[1] * dy,
            pz + right[2] * dx + up[2] * dy,
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
    KITTY_CHUNK = 4096
    RETRANSMIT_EVERY = 20  # self-heal ticks between full re-transmissions

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._last_atoms: list | None = None
        # Region the currently-stored image was rendered for. Placing it into
        # a different one would scale the wrong pixels into the wrong box.
        self._placed_region: tuple[int, int, int, int] | None = None
        self._heals_since_transmit = 0

    def on_resize(self) -> None:
        self._input.request()

    def on_unmount(self) -> None:
        self._input.cancel()
        self._write(f"\x1b_Ga=d,d=i,i={self.IMAGE_ID}\x1b\\")

    def show_job(self, job: Job | None, point=None) -> None:
        is_new_selection = job is not self._job or point is not self._point
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
            self._on_change()
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
            self._on_change()
            return
        self._write(self._placement(region))

    def _placement(self, region) -> str:
        """Re-place the already-transmitted image. `q=2` suppresses the
        terminal's per-command acknowledgement, which otherwise comes back up
        the wire and lands in Textual's own input parser."""
        return (
            f"\x1b[s\x1b[{region.y + 1};{region.x + 1}H"
            f"\x1b_Ga=p,i={self.IMAGE_ID},q=2,c={region.width},r={region.height}\x1b\\"
            "\x1b[u"
        )

    def _on_change(self) -> None:
        self._push(self._next_seq())

    def _write(self, data: str) -> None:
        driver = self.app._driver
        if driver is not None:
            driver.write(data)

    @work(thread=True, exclusive=True)
    def _push(self, seq: int) -> None:
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
                    self._write(f"\x1b_Ga=d,d=i,i={self.IMAGE_ID}\x1b\\")
            return
        region = self.content_region
        if region.width <= 0 or region.height <= 0:
            return

        image = geometry_render.render(
            atoms, elev=self.elev, azim=self.azim, show_distances=self.show_distances,
            zoom=self.zoom, pan=self.pan, qm_atom_indices=job.state.qm_atom_indices,
        )
        buf = io.BytesIO()
        image.convert("RGB").save(buf, format="PNG")
        data_b64 = base64.b64encode(buf.getvalue()).decode("ascii")

        chunks = [data_b64[i : i + self.KITTY_CHUNK] for i in range(0, len(data_b64), self.KITTY_CHUNK)] or [""]
        sequence = [f"\x1b[s\x1b[{region.y + 1};{region.x + 1}H"]
        for i, chunk in enumerate(chunks):
            controls = ["q=2"]
            if i == 0:
                controls += ["a=T", "f=100", "t=d", f"i={self.IMAGE_ID}", f"c={region.width}", f"r={region.height}"]
            controls.append("m=1" if i != len(chunks) - 1 else "m=0")
            sequence.append(f"\x1b_G{','.join(controls)};{chunk}\x1b\\")
        sequence.append("\x1b[u")

        with self._write_lock:
            if self._is_stale(seq):
                return
            self._write("".join(sequence))
            # `a=T` stores under IMAGE_ID as well as displaying, so later
            # self-heals can re-place these pixels instead of resending them.
            self._placed_region = (region.x, region.y, region.width, region.height)
            self._heals_since_transmit = 0


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
    # A preview costs one cheap render and about 20 KB.
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
        self._input.cancel()
        if self._settle_timer is not None:
            self._settle_timer.stop()
        herdr_graphics.clear(self.LAYER_ID)

    def show_job(self, job: Job | None, point=None) -> None:
        # `Job` instances persist for the app's lifetime and are mutated in
        # place (see Job.refresh), so an `is` check tells a genuine selection
        # change (worth the preview flash for snappy feedback) apart from the
        # periodic refresh re-showing the SAME job -- and a `cycle` change
        # (the chart scrubbing to a different geometry) is exactly as much a
        # "new selection" as a different job, interaction-wise.
        is_new_selection = job is not self._job or point is not self._point
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
        # Render straight to the size that will actually be sent, which for a
        # full-quality frame is now simply the pane's own pixel extent. It
        # used to be whatever a 480 KB raw-RGBA budget allowed -- about
        # 411x291, well under the pane -- and that budget went away with the
        # switch to PNG, which carries a full-resolution frame in ~95 KB.
        width = max(1, region.width * cells.width_px)
        height = max(1, region.height * cells.height_px)
        if quality == "preview":
            width = max(1, int(width * PREVIEW_SCALE))
            height = max(1, int(height * PREVIEW_SCALE))
        target = (width, height)
        image = geometry_render.render(
            atoms, elev=self.elev, azim=self.azim, show_distances=self.show_distances,
            zoom=self.zoom, pan=self.pan, qm_atom_indices=job.state.qm_atom_indices,
            size_px=target,
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

    def __init__(self, root: Path):
        super().__init__()
        self.root = root.resolve()
        self.jobs: list[Job] = []
        self.selected_label: str | None = None
        self.maximized = False
        self.use_herdr_graphics = herdr_graphics.available()
        self._geometry_mode = "pixel, herdr graphics" if self.use_herdr_graphics else "pixel, raw kitty"
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
                if self.use_herdr_graphics:
                    yield HerdrGeometryImage(id="geometry")
                else:
                    yield KittyGeometryImage(id="geometry")
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
        job_dirs = discover_jobs(self.root)
        self.jobs = [Job(d, stem, self.root) for d, stem in job_dirs]

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
                job.status.value,
                "-",
                "-",
                "-",
                key=job.label,
            )
        if self.jobs:
            self.selected_label = self.jobs[0].label

        self.query_one("#geometry").border_title = f"geometry ({self._geometry_mode})"
        self.query_one("#summary", Static).border_title = "summary"
        self.query_one("#chart", ConvergencePlot).border_title = "convergence"
        self.query_one("#convergence", Static).border_title = "geometry steps"
        self.query_one("#tail", RichLog).border_title = "output tail"

        table.focus()
        self.refresh_all()
        self.set_interval(REFRESH_SECONDS, self.refresh_all)

    def action_refresh_now(self) -> None:
        self.refresh_all()

    def action_toggle_maximize(self) -> None:
        self.maximized = not self.maximized
        self.screen.set_class(self.maximized, "maximized")

    def refresh_all(self) -> None:
        self._scan()

    @work(thread=True, exclusive=True)
    def _scan(self) -> None:
        """Read every job's new output, off the main thread.

        The first pass reads each job.out end to end -- 17 MB across this tree
        -- and it used to run on the main thread inside on_mount, so nothing
        appeared until it finished. Later passes are cheap, but only because
        nothing has usually been appended; a job writing hard, or a tree with
        big outputs, would stall the UI for as long as the parse took."""
        if not self._scan_lock.acquire(blocking=False):
            return  # a previous scan is still going; this tick can be skipped
        try:
            running_cwds = running_orca_cwds()
            for job in self.jobs:
                job.refresh(running_cwds)
        finally:
            self._scan_lock.release()
        self.call_from_thread(self._apply_scan)

    def _apply_scan(self) -> None:
        table = self.query_one("#job_table", DataTable)
        for job in self.jobs:
            neg_eig = "-"
            if job.state.eigen_history:
                neg_eig = str(job.state.eigen_history[-1][1])
            cells = (
                job.status.value,
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
                Text(status, style=STATUS_STYLE[job.status]),
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
            self._show_text(summary, "No job selected")
            self._show_text(convergence, "")
            self._show_tail(tail, None)
            return

        self._show_text(summary, summary_text(job))
        self._show_text(convergence, steps_text(job.state))
        self._show_tail(tail, job)

    def _title_geometry(self, geometry, job: Job | None, point) -> None:
        """Say WHICH geometry the pane is drawing -- and, when the requested
        one had no coordinates printed, that it is drawing another."""
        _atoms, shown = geometry_shown(job, point)
        title = f"geometry ({self._geometry_mode})"
        where = describe_point(shown)
        if point is not None and shown is not point:
            asked = describe_point(point)
            where = f"{asked}: not printed, showing {where or 'none'}"
        if where:
            title += f" · {where}"
        if geometry.border_title != title:
            geometry.border_title = title


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


def main() -> None:
    """`python -m computational.monitor.app [ROOT]` -- watch every ORCA job
    (every directory holding a .inp) under ROOT; by default the directory
    this package sits in."""
    import sys

    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
    MonitorApp(root).run()


if __name__ == "__main__":
    main()
