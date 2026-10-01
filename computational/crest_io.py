"""
Every number this project reads out of a standalone CREST or xtb run.

The counterpart of `orca_io` for the two programs that have no structured
output of their own worth trusting beyond what is used here:

- CREST writes ensembles as pasted XYZ frames whose comment line is the total
  energy in Eh, `crest.energies` as `index  relative-kcal/mol`, and reports
  success only in the LAST line of `crest.out`. A run that dies still leaves
  `crest_conformers.xyz` from an earlier stage behind, so the ensemble file is
  never evidence that a run finished. `load_run` checks the termination line.
- xtb writes `xtbout.json` under `--json` (the energy, gap and charges as
  numbers) and a summary box in `xtb.out` for the rest. `normal termination of
  xtb` is the success line; a run that stops early still writes `xtbopt.xyz`.

`--json` is not optional in this project: the total energy is read from it and
never from the log. The two labelled lines the log carries (`TOTAL FREE ENERGY`
and the imaginary-mode count) exist only in the summary box, and are read from
that box by exact label.
"""
import json
import os
import re

import numpy as np

HARTREE_TO_KCAL = 627.5094740631
CREST_OK = "CREST terminated normally."
XTB_OK = "normal termination of xtb"


def read_frames(path):
    """(elements, [coordinates (n, 3) in Angstrom], [comment line energies in
    Eh or None]) of an XYZ file with any number of pasted frames. The atom
    order is the file's; CREST keeps it across an ensemble."""
    with open(path) as handle:
        lines = handle.read().split("\n")
    elements, frames, energies, i = None, [], [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        block = lines[i + 2:i + 2 + n]
        if len(block) < n:
            raise ValueError(f"{path}: frame at line {i + 1} is truncated")
        elements = elements or [row.split()[0] for row in block]
        frames.append(np.array([[float(x) for x in row.split()[1:4]] for row in block]))
        energies.append(_comment_energy(lines[i + 1]))
        i += n + 2
    return elements, frames, energies


def _comment_energy(comment):
    """CREST writes the bare energy; xtb writes `energy: <E> gnorm: ...`."""
    match = re.search(r"(?:energy[:=]\s*)?(-?\d+\.\d+)", comment)
    return float(match.group(1)) if match else None


class CrestRun:
    """A finished CREST run. `energies_eh` and `relative_kcal` are per
    conformer in the file's order, which CREST sorts lowest first."""

    def __init__(self, directory, elements, frames, energies_eh, relative_kcal,
                 wall_seconds, population_of_lowest):
        self.directory = directory
        self.elements = elements
        self.frames = frames
        self.energies_eh = energies_eh
        self.relative_kcal = relative_kcal
        self.wall_seconds = wall_seconds
        self.population_of_lowest = population_of_lowest

    @property
    def natoms(self):
        return len(self.elements)


def load_run(directory):
    """Raises unless `crest.out` ends `CREST terminated normally.`, so a
    crashed run, whose `crest_conformers.xyz` may be a stale earlier stage, and
    a mistyped path are both errors and not an empty answer."""
    out = os.path.join(directory, "crest.out")
    if not os.path.exists(out):
        raise FileNotFoundError(out)
    with open(out, errors="replace") as handle:
        text = handle.read()
    if CREST_OK not in text or "terminated with failures" in text:
        raise RuntimeError(f"{out}: CREST did not terminate normally")
    elements, frames, energies = read_frames(os.path.join(directory, "crest_conformers.xyz"))
    relative = _read_relative(os.path.join(directory, "crest.energies"))
    if len(relative) != len(frames):
        raise ValueError(f"{directory}: {len(frames)} conformers but "
                         f"{len(relative)} lines in crest.energies")
    # every stage prints a wall-time summary; the run's total is the last one
    walls = re.findall(r"\* wall-time:\s+(\d+) d,\s+(\d+) h,\s+(\d+) min,\s+([\d.]+) sec", text)
    wall = walls[-1] if walls else None
    seconds = (86400 * int(wall[0]) + 3600 * int(wall[1]) + 60 * int(wall[2])
               + float(wall[3])) if wall else None
    population = re.search(r"population of lowest in %\s*:\s*([\d.]+)", text)
    return CrestRun(directory, elements, frames, energies, relative, seconds,
                    float(population[1]) if population else None)


def _read_relative(path):
    with open(path) as handle:
        return [float(line.split()[1]) for line in handle if line.strip()]


def host_rmsd(frame_a, frame_b, atoms):
    """Kabsch RMSD in Angstrom over `atoms` (0-based indices) after optimal
    superposition, proper rotations only. Comparing a frozen host to its seed
    reads ~0.003; an unfrozen NCI run on a cyclodextrin reads 1.4-2.2."""
    a = frame_a[atoms] - frame_a[atoms].mean(0)
    b = frame_b[atoms] - frame_b[atoms].mean(0)
    u, _, vt = np.linalg.svd(a.T @ b)
    sign = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1.0, 1.0, sign]) @ vt
    return float(np.sqrt(((a @ rotation - b) ** 2).sum() / len(atoms)))


def cavity_occupancy(frame, host, oxygens, axial=3.0, radial=3.5):
    """How many of the water oxygens (0-based atom indices in `oxygens`) sit
    inside a ring host's cavity: within `axial` A of the host's centroid along
    its ring axis and `radial` A of the axis. The axis is the direction of
    least spread of the host's `host` atoms, which is a ring's normal. The
    cut-offs are a judgement about a cyclodextrin (heavy atoms 0.7-7 A from
    the axis, spanning about -4..5 A along it), not a fact about every host."""
    points = frame[host]
    centre = points.mean(0)
    axis = np.linalg.svd(points - centre)[2][2]
    inside = 0
    for i in oxygens:
        d = frame[i] - centre
        z = float(d @ axis)
        r = float(np.linalg.norm(d - z * axis))
        inside += abs(z) < axial and r < radial
    return int(inside)


class XtbRun:
    def __init__(self, directory, total_energy_eh, gap_ev, free_energy_eh,
                 gradient_norm, imaginary_modes, geometry):
        self.directory = directory
        self.total_energy_eh = total_energy_eh
        self.gap_ev = gap_ev
        self.free_energy_eh = free_energy_eh
        self.gradient_norm = gradient_norm
        self.imaginary_modes = imaginary_modes
        self.geometry = geometry


def load_xtb(directory, optimised=True):
    """Raises unless `xtb.out` carries the success line and `xtbout.json`
    exists. `free_energy_eh` and `imaginary_modes` are None unless the run
    included a Hessian (`--ohess`, `--hess`); the summary box only has them
    then. `geometry` is `xtbopt.xyz` as (elements, coordinates) when
    `optimised`."""
    out = os.path.join(directory, "xtb.out")
    if not os.path.exists(out):
        raise FileNotFoundError(out)
    with open(out, errors="replace") as handle:
        text = handle.read()
    if XTB_OK not in text:
        raise RuntimeError(f"{out}: xtb did not terminate normally")
    with open(os.path.join(directory, "xtbout.json")) as handle:
        data = json.load(handle)
    free = re.search(r"\| TOTAL FREE ENERGY\s+(-?\d+\.\d+) Eh", text)
    imaginary = re.search(r"# imaginary freq\.\s+(\d+)", text)
    gradient = re.search(r"\| GRADIENT NORM\s+(\d+\.\d+) Eh/", text)
    geometry = None
    if optimised:
        elements, frames, _ = read_frames(os.path.join(directory, "xtbopt.xyz"))
        geometry = (elements, frames[-1])
    return XtbRun(directory, data["total energy"], data["HOMO-LUMO gap / eV"],
                  float(free[1]) if free else None,
                  float(gradient[1]) if gradient else None,
                  int(imaginary[1]) if imaginary else None, geometry)


def read_vibspectrum(path):
    """Frequencies in cm-1 from xtb's `vibspectrum` (TM format), the six
    translations and rotations included, as xtb writes them (0.00)."""
    values = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(("$", "#")) or not line.strip():
                continue
            fields = line.split()
            values.append(float(fields[-3] if fields[1].isalpha() else fields[1]))
    return values
