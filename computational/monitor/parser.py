from __future__ import annotations

import re
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

STALL_WINDOW = 15
STALL_IMPROVEMENT_FLOOR = 0.10
TAIL_LINES = 300
HISTORY_LEN = 500

_CYCLE_RE = re.compile(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)")
_EIGEN_RE = re.compile(r"Hessian has\s+(\d+)\s+negative eigenvalue")
_QM2_ERROR_RE = re.compile(r"error in the QM2 calculation")
_CRASH_RE = re.compile(
    r"Aborting|TERMINATING THE RUN|ORCA finished by error termination"
)
_NORMAL_DONE_RE = re.compile(r"ORCA TERMINATED NORMALLY")
_OPT_DONE_RE = re.compile(r"OPTIMIZATION HAS CONVERGED|OPTIMIZATION RUN DONE")
_FINAL_ENERGY_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
_RUNTIME_RE = re.compile(
    r"TOTAL RUN TIME:\s*(\d+)\s*days\s*(\d+)\s*hours\s*(\d+)\s*minutes"
    r"\s*(\d+)\s*seconds\s*(\d+)\s*msec"
)
_FREQ_HEADER_RE = re.compile(r"^VIBRATIONAL FREQUENCIES\s*$")
_FREQ_LINE_RE = re.compile(
    r"^\s*\d+:\s+-?[\d.]+\s+cm\*\*-1(\s+\*\*\*imaginary mode\*\*\*)?"
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

_NUM = r"(-?[\d.]+(?:[eE][-+]?\d+)?)"

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
class JobState:
    path: Path
    stem: str = "job"
    has_out: bool = False
    offset: int = 0
    cycle: int = 0
    convergence_history: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    eigen_history: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    qm2_error_count: int = 0
    crashed_marker: bool = False
    crash_lines: deque = field(default_factory=lambda: deque(maxlen=10))
    normal_completion: bool = False
    opt_converged: bool = False
    final_energy: float | None = None
    n_imaginary: int | None = None
    wall_time_s: float | None = None
    tail: deque = field(default_factory=lambda: deque(maxlen=TAIL_LINES))
    cycle_energies: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    scf_iterations: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    atoms: list = field(default_factory=list)
    # (cycle, atoms) snapshot per printed geometry -- lets the UI scrub back
    # through an optimization or scan instead of only ever showing the latest
    # geometry. Keyed by the same `cycle` value as cycle_energies, updated the
    # same "replace if same cycle repeats, else append" way.
    geometry_history: deque = field(default_factory=lambda: deque(maxlen=HISTORY_LEN))
    # 0-based indices of the high-level ("QM1") layer within `atoms`, for a
    # multilayer (QM/MM, QM/XTB, ONIOM) job -- None for an ordinary job with
    # no layering at all, in which case every atom renders as the QM layer.
    qm_atom_indices: set | None = None

    _pending_step: dict = field(default_factory=dict)
    _in_freq_block: bool = False
    _freq_imaginary_count: int = 0
    _freq_seen_line: bool = False
    _in_geom_block: bool = False
    _pending_atoms: list = field(default_factory=list)
    _in_qm1_composition: bool = False

    def feed_line(self, line: str) -> None:
        self.tail.append(line.rstrip("\n"))

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

        m = _QM1_HEADER_RE.match(line.strip())
        if m:
            self.qm_atom_indices = {int(tok) for tok in m.group(1).split()}
            self._in_qm1_composition = True
        elif self._in_qm1_composition:
            stripped = line.strip()
            if _INDEX_LINE_RE.match(stripped):
                self.qm_atom_indices.update(int(tok) for tok in stripped.split())
            else:
                self._in_qm1_composition = False

        if _NORMAL_DONE_RE.search(line):
            self.normal_completion = True

        if _OPT_DONE_RE.search(line):
            self.opt_converged = True

        m = _FINAL_ENERGY_RE.search(line)
        if m:
            self.final_energy = float(m.group(1))
            if self.cycle:
                if self.cycle_energies and self.cycle_energies[-1][0] == self.cycle:
                    self.cycle_energies[-1] = (self.cycle, self.final_energy)
                else:
                    self.cycle_energies.append((self.cycle, self.final_energy))

        if self.cycle == 0:
            m = _SCF_ITER_RE.match(line)
            if m:
                self.scf_iterations.append((int(m.group(1)), float(m.group(2))))

        m = _RUNTIME_RE.search(line)
        if m:
            d, h, mi, s, ms = (int(x) for x in m.groups())
            self.wall_time_s = d * 86400 + h * 3600 + mi * 60 + s + ms / 1000

        if _FREQ_HEADER_RE.match(line.strip()):
            self._in_freq_block = True
            self._freq_imaginary_count = 0
            self._freq_seen_line = False
        elif self._in_freq_block:
            if _FREQ_LINE_RE.match(line):
                self._freq_seen_line = True
                if "imaginary mode" in line:
                    self._freq_imaginary_count += 1
            elif line.strip() == "":
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
                self.n_imaginary = self._freq_imaginary_count

        if _GEOM_HEADER_RE.match(line.strip()):
            self._in_geom_block = True
            self._pending_atoms = []
        elif self._in_geom_block:
            stripped = line.strip()
            if stripped and set(stripped) == {"-"}:
                pass  # the header's dashed underline
            elif _GEOM_ATOM_RE.match(line):
                m = _GEOM_ATOM_RE.match(line)
                el, x, y, z = m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))
                self._pending_atoms.append((el, x, y, z))
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
                    same_cycle = self.geometry_history and self.geometry_history[-1][0] == self.cycle
                    if not same_cycle:
                        self.atoms = self._pending_atoms
                        self.geometry_history.append((self.cycle, self._pending_atoms))
                    elif len(self._pending_atoms) > len(self.geometry_history[-1][1]):
                        self.atoms = self._pending_atoms
                        self.geometry_history[-1] = (self.cycle, self._pending_atoms)
                self._in_geom_block = False
            else:
                self._in_geom_block = False

        for label, (rx, *_attrs) in _CONV_ITEMS.items():
            m = rx.search(line)
            if m:
                value, tol, conv = float(m.group(1)), float(m.group(2)), m.group(3) == "YES"
                self._pending_step[label] = (value, tol, conv)
                if len(self._pending_step) == len(_CONV_ITEMS):
                    self._finish_step()

    def _finish_step(self) -> None:
        step = GeometryStep(cycle=self.cycle)
        for label, (value, tol, conv) in self._pending_step.items():
            _rx, value_attr, tol_attr, conv_attr = _CONV_ITEMS[label]
            setattr(step, value_attr, value)
            setattr(step, tol_attr, tol)
            setattr(step, conv_attr, conv)
        self.convergence_history.append(step)
        self._pending_step = {}

    def possibly_stalled(self) -> bool:
        history = list(self.convergence_history)
        if len(history) < STALL_WINDOW:
            return False
        recent = history[-STALL_WINDOW:]
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


def read_new_lines(path: Path, offset: int) -> tuple[list[str], int]:
    with open(path, "rb") as f:
        f.seek(offset)
        data = f.read()
    if not data:
        return [], offset
    last_nl = data.rfind(b"\n")
    if last_nl == -1:
        return [], offset
    complete = data[: last_nl + 1]
    new_offset = offset + last_nl + 1
    text = complete.decode("utf-8", errors="replace")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines, new_offset


def update_job(state: JobState) -> None:
    out_path = state.path / f"{state.stem}.out"
    state.has_out = out_path.exists()
    if not state.has_out:
        return
    lines, new_offset = read_new_lines(out_path, state.offset)
    state.offset = new_offset
    for line in lines:
        state.feed_line(line)


def new_state(job_dir: Path, stem: str = "job") -> JobState:
    return JobState(path=job_dir, stem=stem)
