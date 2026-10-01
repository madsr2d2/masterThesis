"""How a normal mode's printed vector becomes a displacement of the geometry
a pane is showing.

`orca_pltvib` writes 2.0 * the printed vector, for carbon and hydrogen alike
(measured against its own trajectory), so the printed vector is the Cartesian
pattern and no mass factor is applied. It is indexed against the Hessian's
atoms, which may be the whole structure (the direct map) or just the QM layer
the Hessian was computed for (the map through `qm_atom_indices`), while the
pane may be showing either -- and a multilayer `.out` can print only its QM
block (a measured 22 atoms for a 140-atom Hessian). A mode that fits neither
the pane nor its QM subset is drawn on the job's own structure, and one that
fits none of them is refused rather than drawn against the wrong atoms.
"""
from __future__ import annotations

import math

from .geometry import read_xyz

MODE_AMPLITUDE_ANGSTROM = 0.20   # how far the largest-moving atom travels
MODE_PERIOD_S = 2.5              # one full oscillation, seconds


def imaginary_modes(state) -> list:
    """The imaginary modes, most negative first; [] until a FINAL block."""
    if state.modes is None or not state.freqs_final:
        return []
    return sorted(state.modes, key=lambda mode: mode.cm1)


def target_indices(n_atoms: int, n_mode_atoms: int, qm_atom_indices) -> list | None:
    """Which displayed atom each mode row belongs to, or None."""
    if n_mode_atoms == n_atoms:
        return list(range(n_atoms))
    # ORCA's QM indices are into the FULL system, so only a full-system list fits.
    if (qm_atom_indices and len(qm_atom_indices) == n_mode_atoms
            and max(qm_atom_indices) < n_atoms):
        return sorted(qm_atom_indices)
    return None


def alternate_geometries(job) -> list:
    """Structures to try when the pane's own does not hold the mode.

    A multilayer `.out` can print only the QM block while the Hessian covers
    the whole model, so the mode's own structure is read from the file the
    input names. `job.file_geometry` is only ever set when the log printed
    nothing, so it is tried first but is usually None here."""
    candidates = []
    if getattr(job, "file_geometry", None) is not None:
        candidates.append(job.file_geometry.atoms)
    coords_file = getattr(getattr(job, "input", None), "coords_file", None)
    if coords_file:
        atoms = read_xyz(job.state.path / coords_file)
        if atoms:
            candidates.append(atoms)
    return candidates


def mode_geometry(candidates: list, qm_atom_indices, mode) -> tuple | None:
    """The `(atoms, qm_atom_indices)` to draw `mode` on, from the candidates
    in order, or None when none of them holds the mode."""
    n_mode_atoms = len(mode.vector) // 3
    for atoms in candidates:
        if atoms and target_indices(len(atoms), n_mode_atoms, qm_atom_indices) is not None:
            return list(atoms), qm_atom_indices
    return None


def offsets(atoms: list, qm_atom_indices, mode, amplitude=MODE_AMPLITUDE_ANGSTROM) -> list | None:
    """Per-atom (dx, dy, dz): `mode` scaled so its largest atom moves `amplitude`."""
    n_mode_atoms = len(mode.vector) // 3
    indices = target_indices(len(atoms), n_mode_atoms, qm_atom_indices)
    if indices is None:
        return None
    raw = [(mode.vector[3 * j], mode.vector[3 * j + 1], mode.vector[3 * j + 2])
           for j in range(n_mode_atoms)]
    biggest = max((math.sqrt(x * x + y * y + z * z) for x, y, z in raw), default=0.0)
    if biggest == 0.0:
        return None
    scale = amplitude / biggest
    result = [(0.0, 0.0, 0.0)] * len(atoms)
    for j, atom_index in enumerate(indices):
        x, y, z = raw[j]
        result[atom_index] = (x * scale, y * scale, z * scale)
    return result


def displaced(atoms: list, mode_offsets: list, sine: float) -> list:
    """`atoms` moved `sine` of the way along the offsets."""
    return [(el, x + sine * dx, y + sine * dy, z + sine * dz)
            for (el, x, y, z), (dx, dy, dz) in zip(atoms, mode_offsets)]


def phase_sine(elapsed_s: float, period_s: float = MODE_PERIOD_S) -> float:
    return math.sin(2.0 * math.pi * elapsed_s / period_s)
