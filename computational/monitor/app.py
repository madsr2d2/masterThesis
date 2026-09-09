from __future__ import annotations

import base64
import io
import threading
import time
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
from .parser import JobState, new_state, update_job
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

ROTATE_STEP_DEG = 15.0


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

    def refresh(self, running_cwds: dict[Path, float]) -> None:
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


SELECTION_MARKER = "○"
SELECTION_COLOR = "red"


class ConvergencePlot(PlotextPlot):
    """Cycle energy for a geometry optimization, or SCF iteration energy
    for a single-point job -- whichever series the selected job has.

    For a multi-geometry job (an optimization or scan), left/right when this
    widget is focused scrubs `selected_index` back and forth through the
    plotted points, drawn with a ring marker -- and the geometry pane follows
    it via `MonitorApp.update_detail`, which reads `selected_cycle()`. A
    single-geometry job (bare SCF iterations, or only one cycle so far) has
    nothing to scrub through, so navigation is a no-op there."""

    can_focus = True
    BINDINGS = [
        ("left", "select_prev", "Prev geom"),
        ("right", "select_next", "Next geom"),
    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._job: Job | None = None
        self._xs: list[int] = []
        self.selected_index: int | None = None
        # True until the user deliberately scrubs backward: keeps the ring
        # (and the geometry pane) tracking the newest point as a running
        # optimization writes more cycles, rather than freezing wherever it
        # happened to be the first time a multi-point series appeared.
        self._following_latest = True

    def show_job(self, job: Job | None) -> None:
        is_new_selection = job is not self._job
        self._job = job
        self.plt.clear_data()
        self.plt.clear_figure()
        if job is None:
            self._xs = []
            self.selected_index = None
            self._following_latest = True
            self.refresh()
            return

        state = job.state
        if state.cycle_energies:
            xs = [c for c, _ in state.cycle_energies]
            ys = [e for _, e in state.cycle_energies]
            self.plt.title("geometry optimization: energy vs cycle")
            self.plt.xlabel("cycle")
            self.plt.ylabel("E (Eh)")
            self.plt.plot(xs, ys, marker="braille")
            self._xs = xs
            if is_new_selection:
                self._following_latest = True
            if self._following_latest or self.selected_index is None or self.selected_index >= len(xs):
                self.selected_index = len(xs) - 1
            if len(xs) > 1:
                self.plt.scatter(
                    [xs[self.selected_index]], [ys[self.selected_index]],
                    marker=SELECTION_MARKER, color=SELECTION_COLOR,
                )
        elif state.scf_iterations:
            xs = [i for i, _ in state.scf_iterations]
            ys = [e for _, e in state.scf_iterations]
            self.plt.title("SCF convergence: energy vs iteration")
            self.plt.xlabel("iteration")
            self.plt.ylabel("E (Eh)")
            self.plt.plot(xs, ys, marker="braille")
            self._xs = []
            self.selected_index = None
            self._following_latest = True
        else:
            self.plt.title("no convergence data yet")
            self._xs = []
            self.selected_index = None
            self._following_latest = True
        self.refresh()

    def selected_cycle(self) -> int | None:
        if self.selected_index is None or not self._xs:
            return None
        return self._xs[self.selected_index]

    def action_select_prev(self) -> None:
        if self.selected_index is None or len(self._xs) <= 1:
            return
        self._following_latest = False
        self.selected_index = max(0, self.selected_index - 1)
        self.app.update_detail()

    def action_select_next(self) -> None:
        if self.selected_index is None or len(self._xs) <= 1:
            return
        self.selected_index = min(len(self._xs) - 1, self.selected_index + 1)
        self._following_latest = self.selected_index == len(self._xs) - 1
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
        self._cycle: int | None = None
        self._request_seq = 0
        self._write_lock = threading.Lock()

    def render(self) -> str:
        return ""

    def _next_seq(self) -> int:
        self._request_seq += 1
        return self._request_seq

    def _is_stale(self, seq: int) -> bool:
        return seq != self._request_seq

    @staticmethod
    def _atoms_for(job: Job | None, cycle: int | None) -> list:
        """The atoms for `cycle` (a scrub selection from the chart), or the
        job's latest geometry when `cycle` is None (nothing selected there,
        or the job has only one geometry to begin with)."""
        if job is None:
            return []
        if cycle is None:
            return job.state.atoms
        for c, atoms in job.state.geometry_history:
            if c == cycle:
                return atoms
        return job.state.atoms

    def action_rotate_left(self) -> None:
        self.azim = (self.azim - ROTATE_STEP_DEG) % 360
        self._on_change()

    def action_rotate_right(self) -> None:
        self.azim = (self.azim + ROTATE_STEP_DEG) % 360
        self._on_change()

    def action_rotate_up(self) -> None:
        # No clamp, same as azim: elev just keeps turning past the pole
        # rather than stopping there (matplotlib's view_init renders that
        # fine -- it's a full tumble, not a fixed "top/bottom" limit).
        self.elev = (self.elev + ROTATE_STEP_DEG) % 360
        self._on_change()

    def action_rotate_down(self) -> None:
        self.elev = (self.elev - ROTATE_STEP_DEG) % 360
        self._on_change()

    def action_toggle_distances(self) -> None:
        self.show_distances = not self.show_distances
        self._on_change()

    def action_zoom_in(self) -> None:
        self.zoom = min(ZOOM_MAX, self.zoom * ZOOM_STEP)
        self._on_change()

    def action_zoom_out(self) -> None:
        self.zoom = max(ZOOM_MIN, self.zoom / ZOOM_STEP)
        self._on_change()

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
        self._on_change()

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
    refresh), which is why every trigger below, including that periodic
    refresh, unconditionally re-sends the image rather than only sending it
    on change: it is a self-heal, not just a redraw."""

    IMAGE_ID = 1
    KITTY_CHUNK = 4096

    def on_resize(self) -> None:
        self._on_change()

    def on_unmount(self) -> None:
        self._write(f"\x1b_Ga=d,d=i,i={self.IMAGE_ID}\x1b\\")

    def show_job(self, job: Job | None, cycle: int | None = None) -> None:
        self._job = job
        self._cycle = cycle
        self._on_change()

    def _on_change(self) -> None:
        self._push(self._next_seq())

    def _write(self, data: str) -> None:
        driver = self.app._driver
        if driver is not None:
            driver.write(data)

    @work(thread=True, exclusive=True)
    def _push(self, seq: int) -> None:
        atoms = self._atoms_for(self._job, self._cycle)
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
            zoom=self.zoom, pan=self.pan,
        )
        buf = io.BytesIO()
        image.convert("RGB").save(buf, format="PNG")
        data_b64 = base64.b64encode(buf.getvalue()).decode("ascii")

        chunks = [data_b64[i : i + self.KITTY_CHUNK] for i in range(0, len(data_b64), self.KITTY_CHUNK)] or [""]
        sequence = [f"\x1b[s\x1b[{region.y + 1};{region.x + 1}H"]
        for i, chunk in enumerate(chunks):
            controls = []
            if i == 0:
                controls += ["a=T", "f=100", "t=d", f"i={self.IMAGE_ID}", f"c={region.width}", f"r={region.height}"]
            controls.append("m=1" if i != len(chunks) - 1 else "m=0")
            sequence.append(f"\x1b_G{','.join(controls)};{chunk}\x1b\\")
        sequence.append("\x1b[u")

        with self._write_lock:
            if self._is_stale(seq):
                return
            self._write("".join(sequence))


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
    SETTLE_DELAY_S = 0.35

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._settle_timer = None
        self._last_fingerprint: tuple[int, int] | None = None

    def on_resize(self) -> None:
        self._on_change()

    def on_unmount(self) -> None:
        if self._settle_timer is not None:
            self._settle_timer.stop()
        herdr_graphics.clear(self.LAYER_ID)

    def show_job(self, job: Job | None, cycle: int | None = None) -> None:
        # `Job` instances persist for the app's lifetime and are mutated in
        # place (see Job.refresh), so an `is` check tells a genuine selection
        # change (worth the preview flash for snappy feedback) apart from the
        # periodic refresh re-showing the SAME job -- and a `cycle` change
        # (the chart scrubbing to a different geometry) is exactly as much a
        # "new selection" as a different job, interaction-wise.
        is_new_selection = job is not self._job or cycle != self._cycle
        self._job = job
        self._cycle = cycle
        if is_new_selection:
            self._last_fingerprint = None
            self._on_change()
            return
        # A running optimization writes new coordinates between ticks, so a
        # periodic refresh is still worth a redraw -- but only when the
        # geometry actually shown advanced, not unconditionally every 3s. A
        # scrubbed-to historical geometry is frozen (its cycle is fixed), so
        # this also makes the periodic tick a no-op while reviewing one.
        atoms = self._atoms_for(job, cycle)
        fingerprint = (len(atoms), cycle if cycle is not None else job.state.cycle) if job is not None else None
        if fingerprint == self._last_fingerprint:
            return
        self._last_fingerprint = fingerprint
        self._push(self._next_seq(), quality="full")

    def _on_change(self) -> None:
        if self._settle_timer is not None:
            self._settle_timer.stop()
        seq = self._next_seq()
        self._push(seq, quality="preview")
        self._settle_timer = self.set_timer(self.SETTLE_DELAY_S, lambda: self._push(seq, quality="full"))

    @work(thread=True, exclusive=True)
    def _push(self, seq: int, quality: str = "full") -> None:
        atoms = self._atoms_for(self._job, self._cycle)
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
        image = geometry_render.render(
            atoms, elev=self.elev, azim=self.azim, show_distances=self.show_distances,
            zoom=self.zoom, pan=self.pan,
        )
        image.thumbnail((max(1, region.width * cells.width_px), max(1, region.height * cells.height_px)))
        max_bytes = (
            herdr_graphics.PREVIEW_MAX_RGBA_BYTES if quality == "preview" else herdr_graphics.MAX_RGBA_BYTES
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
                max_bytes=max_bytes,
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
        max-height: 12;
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
                yield RichLog(id="tail", wrap=False, highlight=False, markup=False)
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

        geometry_mode = "pixel, herdr graphics" if self.use_herdr_graphics else "pixel, raw kitty"
        self.query_one("#geometry").border_title = f"geometry ({geometry_mode})"
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
        running_cwds = running_orca_cwds()
        table = self.query_one("#job_table", DataTable)
        for job in self.jobs:
            job.refresh(running_cwds)
            neg_eig = "-"
            if job.state.eigen_history:
                neg_eig = str(job.state.eigen_history[-1][1])
            style = STATUS_STYLE[job.status]
            table.update_cell(
                job.label,
                "status",
                Text(job.status.value, style=style),
                update_width=False,
            )
            table.update_cell(job.label, "cycle", str(job.state.cycle))
            table.update_cell(job.label, "neg_eig", neg_eig)
            table.update_cell(
                job.label, "wall_time", format_wall_time(job.wall_time_s)
            )
        self.update_detail()

    def selected_job(self) -> Job | None:
        for job in self.jobs:
            if job.label == self.selected_label:
                return job
        return None

    def on_data_table_row_highlighted(self, event: DataTable.RowHighlighted) -> None:
        if event.row_key.value is not None:
            self.selected_label = event.row_key.value
        self.update_detail()

    def update_detail(self) -> None:
        job = self.selected_job()
        summary = self.query_one("#summary", Static)
        chart = self.query_one("#chart", ConvergencePlot)
        geometry = self.query_one("#geometry")
        convergence = self.query_one("#convergence", Static)
        tail = self.query_one("#tail", RichLog)

        chart.show_job(job)
        geometry.show_job(job, cycle=chart.selected_cycle())

        if job is None:
            summary.update("No job selected")
            convergence.update("")
            return

        state = job.state
        style = STATUS_STYLE[job.status]
        lines = [
            f"[b]{job.label}[/b]",
            f"status: [{style}]{job.status.value}[/{style}]",
            f"cycle: {state.cycle}",
            f"wall time: {format_wall_time(job.wall_time_s)}",
            f"QM2 errors: {state.qm2_error_count}",
        ]
        if state.final_energy is not None:
            lines.append(f"final single point energy: {state.final_energy:.10f} Eh")
        if state.n_imaginary is not None:
            lines.append(f"imaginary frequencies (latest Hessian): {state.n_imaginary}")
        if state.eigen_history:
            recent = list(state.eigen_history)[-10:]
            trail = ", ".join(f"c{c}:{n}" for c, n in recent)
            lines.append(f"Hessian neg. eigenvalues (recent): {trail}")
        if state.crash_lines:
            lines.append("[bold red]crash markers:[/bold red]")
            for cl in state.crash_lines:
                lines.append(f"  {cl}")
        summary.update("\n".join(lines))

        if state.convergence_history:
            rows = ["[b]cycle   item             value          tolerance      conv[/b]"]
            for step in list(state.convergence_history)[-8:]:
                for name, val, tol, conv in [
                    ("Energy change", step.energy_change, step.energy_tol, step.energy_conv),
                    ("RMS gradient", step.rms_grad, step.rms_grad_tol, step.rms_grad_conv),
                    ("MAX gradient", step.max_grad, step.max_grad_tol, step.max_grad_conv),
                    ("RMS step", step.rms_step, step.rms_step_tol, step.rms_step_conv),
                    ("MAX step", step.max_step, step.max_step_tol, step.max_step_conv),
                ]:
                    if val is None:
                        continue
                    mark = "YES" if conv else "no"
                    color = "green" if conv else "red"
                    rows.append(
                        f"{step.cycle:>5}   {name:<15} {val:>13.7f}  {tol:>13.7f}  "
                        f"[{color}]{mark}[/{color}]"
                    )
            convergence.update("\n".join(rows))
        else:
            convergence.update("[dim]no geometry convergence data yet[/dim]")

        tail.clear()
        for line in state.tail:
            tail.write(line)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    MonitorApp(root).run()


if __name__ == "__main__":
    main()
