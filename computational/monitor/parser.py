from __future__ import annotations

import re
from collections import deque
from dataclasses import dataclass, field, fields
from pathlib import Path

STALL_WINDOW = 15
STALL_IMPROVEMENT_FLOOR = 0.10
TAIL_LINES = 300
HISTORY_LEN = 500
# Bytes read per pass over a job.out. Reading a whole file at once held the
# bytes, the decoded text and a list of every line at the same time: parsing
# this tree's 64 jobs (667 MB) peaked at 334 MB resident. In 4 MB chunks it
# peaks at 121 MB -- most of it the parsed state the app keeps anyway -- and
# is no slower; 16 MB chunks peaked at 183 MB for nothing.
READ_CHUNK_BYTES = 4 * 1024 * 1024

_CYCLE_RE = re.compile(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)")
_EIGEN_RE = re.compile(r"Hessian has\s+(\d+)\s+negative eigenvalue")
_QM2_ERROR_RE = re.compile(r"error in the QM2 calculation")
_CRASH_RE = re.compile(
    r"Aborting|TERMINATING THE RUN|ORCA finished by error termination"
)
_NORMAL_DONE_RE = re.compile(r"ORCA TERMINATED NORMALLY")
_OPT_DONE_RE = re.compile(r"OPTIMIZATION HAS CONVERGED|OPTIMIZATION RUN DONE")
# The label is optional. A multilayer run prints one line per sub-calculation
# as well as the combined total, e.g. for QM/XTB, in this order:
#     FINAL SINGLE POINT ENERGY (L-QM2)     -266.681990679470
#     FINAL SINGLE POINT ENERGY (S-QM2)      -40.834442837800
#     FINAL SINGLE POINT ENERGY      -647.909099472771
#     FINAL SINGLE POINT ENERGY (QM/QM2)     -873.756647314441
# and a pattern that required the bare form took the third -- the high-level
# region ALONE -- as the job's energy. `_SUBSYSTEM_LABEL_RE` names the
# extrapolation's sub-calculations (Large/Small model); of what remains, the
# LAST line printed for a geometry is the combined total.
_FINAL_ENERGY_RE = re.compile(
    r"FINAL SINGLE POINT ENERGY(?:\s+\(([^)]*)\))?\s+(-?\d+\.\d+)"
)
_SUBSYSTEM_LABEL_RE = re.compile(r"^[A-Z]-")
_SCAN_STEP_RE = re.compile(r"RELAXED SURFACE SCAN STEP\s+(\d+)")
_SCAN_TOTAL_RE = re.compile(r"There will be\s+(\d+)\s+constrained geometry optimizations")
_SCAN_PARAMS_RE = re.compile(r"There (?:is|are)\s+(\d+)\s+parameters? to be scanned")
# Inside the step banner: "*   Bond (130, 128)  :   1.40000000   *"
_SCAN_VALUE_RE = re.compile(r"^\s*\*\s+(\S.*?)\s*:\s+(-?\d+\.\d+)\s+\*\s*$")
# The orbital basis ORCA used -- including the one a composite method
# (r2SCAN-3c, B97-3c, ...) brings with it, which the input never names.
# Auxiliary (RI) bases get sections of their own and are not recorded:
#     ----- Orbital basis set information -----
#     Your calculation utilizes the basis: def2-mTZVPP
_BASIS_SECTION_RE = re.compile(r"^-+\s*(\S+) basis set information\s*-+\s*$")
_BASIS_NAME_RE = re.compile(r"Your calculation utilizes the basis:\s*(\S.*?)\s*$")
_MAXITER_RE = re.compile(r"Max\. no of cycles\s+MaxIter\s+\.+\s+(\d+)")
_RUNTIME_RE = re.compile(
    r"TOTAL RUN TIME:\s*(\d+)\s*days\s*(\d+)\s*hours\s*(\d+)\s*minutes"
    r"\s*(\d+)\s*seconds\s*(\d+)\s*msec"
)
_FREQ_HEADER_RE = re.compile(r"^VIBRATIONAL FREQUENCIES\s*$")
_FREQ_LINE_RE = re.compile(
    r"^\s*\d+:\s+(-?[\d.]+)\s+cm\*\*-1(\s+\*\*\*imaginary mode\*\*\*)?"
)
_GEOM_HEADER_RE = re.compile(r"^CARTESIAN COORDINATES \(ANGSTROEM\)\s*$")
_GEOM_ATOM_RE = re.compile(
    r"^\s*([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)
_SCF_ITER_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d{4,})\s+-?\d+\.\d+e[+-]\d+")
# A multilayer (QM/MM, QM/XTB, ONIOM) job prints this ONCE at startup, under
# "Composition of different systems (atoms start counting at 0):" -- the
# high-level ("QM1") layer's 0-based atom indices within the FULL system's
# coordinate list, wrapped across as many indented continuation lines as
# needed. Everything not in this set is the lower-level layer.
_QM1_HEADER_RE = re.compile(r"^QM1 Subsystem\s*\.\.\.\s*(.*)$")
_INDEX_LINE_RE = re.compile(r"^\d+(?:\s+\d+)*$")

# Prefilters. Each is the union of the LITERAL substrings its group's patterns
# require, so a line that fails it cannot match any of them -- they narrow the
# work without deciding anything. Keep them in step with the patterns above:
# adding a marker pattern without adding its literal here would silently stop
# that marker being seen. `parser` has no gate of its own, so
# `computational/monitor/validate.py` is where that is checked.
_RARE_MARKERS_RE = re.compile(
    r"GEOMETRY OPTIMIZATION CYCLE|Hessian has|error in the QM2 calculation"
    r"|Aborting|TERMINATING THE RUN|ORCA finished by error termination"
    r"|ORCA TERMINATED NORMALLY|OPTIMIZATION HAS CONVERGED|OPTIMIZATION RUN DONE"
    r"|FINAL SINGLE POINT ENERGY|TOTAL RUN TIME|RELAXED SURFACE SCAN STEP"
    r"|constrained geometry optimizations|to be scanned|Max\. no of cycles"
    r"|basis set information|utilizes the basis:"
)
_CONV_HINT_RE = re.compile(r"gradient|step|Energy change")

_NUM = r"(-?[\d.]+(?:[eE][-+]?\d+)?)"

# Every line that can change parser state OUTSIDE a block, as ONE pattern run
# over a whole chunk at C speed. On an 85 MB output that is 11,745 lines of
# 1,269,463 found in 0.16 s, where handing every line to `_process` took
# 1.55 s -- the regexes were cheap, the Python call per line was not. Inside
# a block (coordinates, frequencies, the QM1 index list, a scan banner) every
# line matters and `feed_text` walks them one by one; see `_in_block`.
#
# Each alternative is a literal some per-line pattern REQUIRES, so a line
# this skips could not have matched anything -- the same contract as the
# prefilters above, one level out, and held the same way: `validate.py`
# drives each marker through `feed_text`, not only `feed_line`.
_LINE_OF_INTEREST_PARTS = (
    _RARE_MARKERS_RE.pattern,
    "QM1 Subsystem",
    "Energy change|RMS gradient|MAX gradient|RMS step|MAX step",
    "VIBRATIONAL FREQUENCIES",
    r"CARTESIAN COORDINATES \(ANGSTROEM\)",
)
_LINE_OF_INTEREST_RE = re.compile("|".join(_LINE_OF_INTEREST_PARTS))
# SCF iteration rows carry no literal, and are only read before the first
# optimization cycle -- so they join the scan only while `cycle == 0`, and a
# long optimization does not pay for them.
_LINE_OF_INTEREST_SCF_RE = re.compile(
    "|".join((*_LINE_OF_INTEREST_PARTS, _SCF_ITER_RE.pattern)), re.M
)

# label -> (regex, GeometryStep value/tol/conv attribute names)
_CONV_ITEMS = {
    "energy_change": (re.compile(rf"Energy change\s+{_NUM}\s+{_NUM}\s+(YES|NO)"), "energy_change", "energy_tol", "energy_conv"),
    "rms_grad": (re.compile(rf"RMS gradient\s+{_NUM}\s+{_NUM}\s+(YES|NO)"), "rms_grad", "rms_grad_tol", "rms_grad_conv"),
    "max_grad": (re.compile(rf"MAX gradient\s+{_NUM}\s+{_NUM}\s+(YES|NO)"), "max_grad", "max_grad_tol", "max_grad_conv"),
    "rms_step": (re.compile(rf"RMS step\s+{_NUM}\s+{_NUM}\s+(YES|NO)"), "rms_step", "rms_step_tol", "rms_step_conv"),
    "max_step": (re.compile(rf"MAX step\s+{_NUM}\s+{_NUM}\s+(YES|NO)"), "max_step", "max_step_tol", "max_step_conv"),
}


@dataclass
class GeometryStep:
    cycle: int
    scan_step: int | None = None
    energy_change: float | None = None
    energy_tol: float | None = None
    energy_conv: bool | None = None
    rms_grad: float | None = None
    rms_grad_tol: float | None = None
    rms_grad_conv: bool | None = None
    max_grad: float | None = None
    max_grad_tol: float | None = None
    max_grad_conv: bool | None = None
    rms_step: float | None = None
    rms_step_tol: float | None = None
    rms_step_conv: bool | None = None
    max_step: float | None = None
    max_step_tol: float | None = None
    max_step_conv: bool | None = None


@dataclass
class GeometryPoint:
    """One geometry the job computed an energy for: an optimization cycle, a
    scan step's cycle, or a single point's one structure.

    Numbered by `(scan_step, cycle)` because ORCA restarts the cycle counter
    at every scan step -- keyed by cycle alone, a 15-step scan's points ran
    1..9, 1..16, ... and "cycle 3" named fifteen different geometries."""

    scan_step: int | None
    cycle: int
    atoms: list = field(default_factory=list)
    energy: float | None = None
    energy_label: str | None = None

    @property
    def key(self) -> tuple[int | None, int]:
        return (self.scan_step, self.cycle)


@dataclass
class JobState:
    path: Path
    stem: str = "job"
    has_out: bool = False
    offset: int = 0
    # st_ino of the job.out `offset` counts into. A re-run that replaces the
    # file by rename gets a new inode at the same-or-larger size, which the
    # size check alone cannot see -- see `update_job`.
    inode: int | None = None
    # st_mtime of job.out at the last scan -- "how long since it last wrote".
    mtime: float | None = None
    cycle: int = 0
    max_cycles: int | None = None
    scan_step: int | None = None
    scan_total: int | None = None
    scan_params: int | None = None
    # A ONE-parameter scan's coordinate: its name ("Bond (130, 128)") and its
    # value at each step, read off the step banners. Left empty for a 2D scan,
    # whose steps have no single value to plot against.
    scan_label: str | None = None
    scan_values: dict = field(default_factory=dict)
    basis: str | None = None
    convergence_history: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    eigen_history: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    qm2_error_count: int = 0
    crashed_marker: bool = False
    crash_lines: deque = field(default_factory=lambda: deque(maxlen=10))
    normal_completion: bool = False
    opt_converged: bool = False
    final_energy: float | None = None
    final_energy_label: str | None = None
    # Every imaginary frequency (cm**-1, negative) of the LAST frequency block
    # printed, in the order printed; None until a block has been read.
    imaginary_freqs: list | None = None
    wall_time_s: float | None = None
    tail: deque = field(default_factory=lambda: deque(maxlen=TAIL_LINES))
    # Monotonic count of lines ever fed. `tail` is a bounded deque, so it
    # cannot say how much of itself is new; this can, which is what lets the
    # UI append the new lines instead of clearing and rewriting all of them.
    lines_seen: int = 0
    scf_iterations: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    atoms: list = field(default_factory=list)
    # One GeometryPoint per geometry, in the order computed -- lets the UI
    # scrub back through an optimization or scan instead of only ever showing
    # the latest. Bounded, so a long scan keeps its most recent cycles...
    points: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    # ...and, separately and unbounded, the latest point of every scan step:
    # the scan's profile, which is what the chart shows by default, survives
    # however many cycles the steps take between them.
    scan_points: dict = field(default_factory=dict)
    # 0-based indices of the high-level ("QM1") layer within `atoms`, for a
    # multilayer (QM/MM, QM/XTB, ONIOM) job -- None for an ordinary job with
    # no layering at all, in which case every atom renders as the QM layer.
    qm_atom_indices: set | None = None

    _pending_step: dict = field(default_factory=dict)
    _in_freq_block: bool = False
    _freq_seen_line: bool = False
    _in_geom_block: bool = False
    _pending_atoms: list = field(default_factory=list)
    _in_qm1_composition: bool = False
    _in_scan_banner: bool = False
    _in_orbital_basis: bool = False
    _freq_imaginary: list = field(default_factory=list)

    _IDENTITY_FIELDS = frozenset({"path", "stem", "has_out", "offset", "inode", "mtime"})

    @property
    def n_imaginary(self) -> int | None:
        return None if self.imaginary_freqs is None else len(self.imaginary_freqs)

    def current_point(self) -> GeometryPoint:
        """The point for the geometry being computed now, opened on first use.

        Both a coordinate block and an energy can be the first thing seen for
        a new geometry, so either opens it; everything later for the same
        `(scan_step, cycle)` -- a multilayer job's subsystem coordinate
        blocks, the final re-evaluation after convergence -- lands on it."""
        key = (self.scan_step, self.cycle)
        if self.points and self.points[-1].key == key:
            return self.points[-1]
        point = GeometryPoint(scan_step=self.scan_step, cycle=self.cycle)
        self.points.append(point)
        if self.scan_step is not None:
            self.scan_points[self.scan_step] = point
        return point

    def reset_for_restart(self) -> None:
        """Drop everything derived from the file's CONTENTS, keeping only the
        job's identity and the caller's fresh file bookkeeping.

        A re-run recreates job.out from scratch, and without this every
        content-derived field -- the status above all -- stays frozen at the
        previous run's final values for the rest of the session. Reading the
        replacement values off a fresh instance rather than restating them
        here means a field added to this dataclass is reset automatically
        instead of being silently kept across a restart."""
        fresh = JobState(path=self.path, stem=self.stem)
        for f in fields(self):
            if f.name not in self._IDENTITY_FIELDS:
                setattr(self, f.name, getattr(fresh, f.name))

    def feed_line(self, line: str) -> None:
        """One line, the reference path. `feed_text` is the fast one and must
        leave the state exactly as feeding its lines here would."""
        self.tail.append(line.rstrip("\n"))
        self.lines_seen += 1
        self._process(line)

    def feed_text(self, text: str) -> None:
        """A run of COMPLETE lines (`text` ends with a newline).

        The tail and the line count are taken in bulk; then only the lines
        that can change state reach `_process` -- a hit of the pattern, or
        any line while a block is open."""
        if not text:
            return
        self.lines_seen += text.count("\n")
        start = len(text) - 1
        for _ in range(TAIL_LINES):
            start = text.rfind("\n", 0, start)
            if start == -1:
                break
        self.tail.extend(text[start + 1 : -1].split("\n"))

        pos, end = 0, len(text)
        while pos < end:
            if not self._in_block():
                pattern = _LINE_OF_INTEREST_SCF_RE if self.cycle == 0 else _LINE_OF_INTEREST_RE
                m = pattern.search(text, pos)
                if m is None:
                    return
                pos = text.rfind("\n", 0, m.start()) + 1
            nl = text.find("\n", pos)
            self._process(text[pos:nl])
            pos = nl + 1

    def _in_block(self) -> bool:
        return (
            self._in_geom_block or self._in_freq_block
            or self._in_qm1_composition or self._in_scan_banner
        )

    def _process(self, line: str) -> None:
        # Almost nothing in a job.out matches almost any of these patterns --
        # of 249,401 lines in this project's largest output, the cycle banner
        # hits 85 times and the crash, QM2-error and eigenvalue markers zero --
        # yet all of them used to run against all of it: ~13 re.search calls a
        # line, 3.2 million in total, to parse 16.6 MB. The two prefilters
        # below decide in one pass each whether the group is worth trying.
        if _RARE_MARKERS_RE.search(line):
            self._feed_marker(line)

        if self._in_scan_banner:
            self._feed_scan_banner(line)

        if self._in_qm1_composition or "QM1 Subsystem" in line:
            self._feed_qm1(line.strip())

        if self.cycle == 0:
            m = _SCF_ITER_RE.match(line)
            if m:
                self.scf_iterations.append((int(m.group(1)), float(m.group(2))))

        if _CONV_HINT_RE.search(line):
            for label, (rx, *_attrs) in _CONV_ITEMS.items():
                m = rx.search(line)
                if m:
                    value, tol, conv = float(m.group(1)), float(m.group(2)), m.group(3) == "YES"
                    self._pending_step[label] = (value, tol, conv)
                    if len(self._pending_step) == len(_CONV_ITEMS):
                        self._finish_step()

        self._feed_blocks(line, line.strip())

    def _feed_marker(self, line: str) -> None:
        """The one-off and per-cycle markers, reached only when the prefilter
        says one of their literals is present."""
        m = _CYCLE_RE.search(line)
        if m:
            self.cycle = int(m.group(1))

        m = _EIGEN_RE.search(line)
        if m:
            self.eigen_history.append((self.cycle, int(m.group(1))))

        if _QM2_ERROR_RE.search(line):
            self.qm2_error_count += 1

        if _CRASH_RE.search(line):
            self.crashed_marker = True
            self.crash_lines.append(line.strip())

        if _NORMAL_DONE_RE.search(line):
            self.normal_completion = True

        if _OPT_DONE_RE.search(line):
            self.opt_converged = True

        m = _FINAL_ENERGY_RE.search(line)
        if m:
            label, energy = m.group(1), float(m.group(2))
            if label is None or not _SUBSYSTEM_LABEL_RE.match(label):
                self.final_energy, self.final_energy_label = energy, label
                point = self.current_point()
                point.energy, point.energy_label = energy, label

        m = _SCAN_STEP_RE.search(line)
        if m:
            self.scan_step = int(m.group(1))
            self._in_scan_banner = True

        m = _SCAN_TOTAL_RE.search(line)
        if m:
            self.scan_total = int(m.group(1))

        m = _SCAN_PARAMS_RE.search(line)
        if m:
            self.scan_params = int(m.group(1))

        m = _BASIS_SECTION_RE.match(line)
        if m:
            self._in_orbital_basis = m.group(1) == "Orbital"

        m = _BASIS_NAME_RE.search(line)
        if m and self._in_orbital_basis:
            # First report wins: a job that re-prints its basis (a compound
            # job) keeps the one it started with rather than flickering.
            if self.basis is None:
                self.basis = m.group(1)
            self._in_orbital_basis = False

        m = _MAXITER_RE.search(line)
        if m:
            self.max_cycles = int(m.group(1))

        m = _RUNTIME_RE.search(line)
        if m:
            d, h, mi, s, ms = (int(x) for x in m.groups())
            self.wall_time_s = d * 86400 + h * 3600 + mi * 60 + s + ms / 1000

    def _feed_scan_banner(self, line: str) -> None:
        """The lines of a scan step's banner after its title: the scanned
        coordinate's value for this step, then the closing row of stars."""
        if set(line.strip()) == {"*"}:
            self._in_scan_banner = False
            return
        m = _SCAN_VALUE_RE.match(line)
        if m and self.scan_params == 1 and self.scan_step is not None:
            self.scan_label = m.group(1)
            self.scan_values[self.scan_step] = float(m.group(2))

    def _feed_qm1(self, stripped: str) -> None:
        """The QM1 header and its continuation lines. Kept out of the marker
        prefilter because the continuations are bare digits, which carry no
        literal to filter on -- the `_in_qm1_composition` flag is the cheap
        test that admits them."""
        m = _QM1_HEADER_RE.match(stripped)
        if m:
            self.qm_atom_indices = {int(tok) for tok in m.group(1).split()}
            self._in_qm1_composition = True
        elif self._in_qm1_composition:
            if _INDEX_LINE_RE.match(stripped):
                self.qm_atom_indices.update(int(tok) for tok in stripped.split())
            else:
                self._in_qm1_composition = False

    def _feed_blocks(self, line: str, stripped: str) -> None:
        """The frequency and geometry block state machines. Both need their
        own flag checked on every line, so neither can sit behind a prefilter
        -- but they can share the one strip() the three of them used to do
        separately."""
        if _FREQ_HEADER_RE.match(stripped):
            self._in_freq_block = True
            self._freq_imaginary = []
            self._freq_seen_line = False
        elif self._in_freq_block:
            m = _FREQ_LINE_RE.match(line)
            if m:
                self._freq_seen_line = True
                if m.group(2):
                    self._freq_imaginary.append(float(m.group(1)))
            elif stripped == "":
                pass
            elif not self._freq_seen_line:
                # preamble between the header and the first numbered mode --
                # the closing "---" divider and the "Scaling factor for
                # frequencies = ..." line both land here. Only real table rows
                # (matched above) or a line seen AFTER at least one real row
                # should end the block.
                pass
            else:
                self._in_freq_block = False
                self.imaginary_freqs = self._freq_imaginary

        if _GEOM_HEADER_RE.match(stripped):
            self._in_geom_block = True
            self._pending_atoms = []
        elif self._in_geom_block:
            m = _GEOM_ATOM_RE.match(line)
            if m:
                el, x, y, z = m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))
                self._pending_atoms.append((el, x, y, z))
            elif stripped and set(stripped) == {"-"}:
                pass  # the header's dashed underline
            elif stripped == "":
                if self._pending_atoms:
                    # A multilayer job (QM/MM, QM/XTB, ONIOM) prints THREE of
                    # these identically-headed blocks per cycle: the full
                    # system first, then the QM-region-only subset (twice)
                    # while ORCA sets up the embedded calculation -- e.g. 140
                    # atoms, then 14, then 14 again, all under the same
                    # generic "CARTESIAN COORDINATES (ANGSTROEM)" header with
                    # no other marker distinguishing them. Keeping "whichever
                    # printed last" (the old behaviour) kept the QM-only
                    # subset and silently dropped every other layer. Keep the
                    # LARGEST block seen for this cycle instead -- the full
                    # system is always a superset of any QM-region subset, so
                    # it's always the biggest.
                    point = self.current_point()
                    if len(self._pending_atoms) > len(point.atoms):
                        self.atoms = self._pending_atoms
                        point.atoms = self._pending_atoms
                self._in_geom_block = False
            else:
                self._in_geom_block = False

    def _finish_step(self) -> None:
        step = GeometryStep(cycle=self.cycle, scan_step=self.scan_step)
        for label, (value, tol, conv) in self._pending_step.items():
            _rx, value_attr, tol_attr, conv_attr = _CONV_ITEMS[label]
            setattr(step, value_attr, value)
            setattr(step, tol_attr, tol)
            setattr(step, conv_attr, conv)
        self.convergence_history.append(step)
        self._pending_step = {}

    def possibly_stalled(self) -> bool:
        history = self.convergence_history
        if len(history) < STALL_WINDOW:
            return False
        # Index the deque from its end rather than copying all HISTORY_LEN
        # entries out just to slice the last STALL_WINDOW of them.
        recent = [history[i] for i in range(len(history) - STALL_WINDOW, len(history))]
        grads = [s.rms_grad for s in recent if s.rms_grad is not None]
        if len(grads) < STALL_WINDOW:
            return False
        if any(s.rms_grad_conv and s.max_grad_conv for s in recent):
            return False
        half = STALL_WINDOW // 2
        best_early = min(grads[:half])
        best_late = min(grads[half:])
        if best_early <= 0:
            return False
        improvement = (best_early - best_late) / best_early
        return improvement < STALL_IMPROVEMENT_FLOOR


def update_job(state: JobState) -> None:
    out_path = state.path / f"{state.stem}.out"
    try:
        st = out_path.stat()
    except OSError:
        state.has_out = False
        return
    state.has_out = True

    # A re-run recreates job.out from scratch -- the single most common thing
    # this tool watches somebody do. `offset` only ever grew, so the seek
    # landed past EOF and returned nothing forever after: the pane kept
    # showing the PREVIOUS run's cycle, geometry and (worst) its "crashed" or
    # "converged" status for the rest of the session, and never recovered as
    # the new run wrote output. A shrunk file catches a truncate-in-place; the
    # inode catches a rename-and-replace, which can land at a larger size.
    if st.st_size < state.offset or (state.inode is not None and st.st_ino != state.inode):
        state.reset_for_restart()
        state.offset = 0
    state.inode = st.st_ino
    state.mtime = st.st_mtime

    if st.st_size == state.offset:
        return  # nothing appended -- the common case, so don't even open it

    read_appended(state, out_path)


def read_appended(state: JobState, out_path: Path, chunk_bytes: int = READ_CHUNK_BYTES) -> None:
    """Feed everything appended since `state.offset`, in bounded chunks.

    Only complete lines are fed: a chunk is cut at its last newline and the
    remainder carried into the next read, and a line still being written at
    EOF is left for the next tick. Cutting at a newline is also what makes
    decoding a chunk on its own safe -- a UTF-8 multibyte sequence never
    contains the newline byte."""
    with open(out_path, "rb") as f:
        f.seek(state.offset)
        carry = b""
        while True:
            block = f.read(chunk_bytes)
            if not block:
                return
            data = carry + block if carry else block
            last_nl = data.rfind(b"\n")
            if last_nl == -1:
                carry = data
                continue
            state.feed_text(data[: last_nl + 1].decode("utf-8", errors="replace"))
            state.offset += last_nl + 1
            carry = data[last_nl + 1 :]


def new_state(job_dir: Path, stem: str = "job") -> JobState:
    return JobState(path=job_dir, stem=stem)
