"""What in an ORCA input would fail or mislead before the job runs.

`check` exists because an input is cheapest to fix while the job is still
"not run": an atom number past the end of the geometry, an electron count
the multiplicity cannot have, a geometry file that is not there, more
processes or memory than this machine has, or an output a re-run would
overwrite. Everything here is read-only.

The electron count checked is the WHOLE system's (D6): on a multilayer job
the QM region's count depends on the link atoms ORCA inserts, which cannot
be reproduced from the input, so `Charge_Total`/`Mult_Total` are what is
compared. A multilayer job that sets neither is a warning, not an error.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass

from .edit import bonded_pairs
from .geometry import ATOMIC_NUMBERS, read_xyz
from .orca_input import atom_references, qm_atoms


@dataclass
class Finding:
    code: str
    message: str


@dataclass
class InputCheck:
    atoms: int | None
    electrons: int | None
    charge: int | None
    mult: int | None
    errors: list
    warnings: list


def _machine_memory_mb() -> int | None:
    """This machine's physical memory in MB, or None where `sysconf` cannot
    say (a platform without it, or a quantity it does not support)."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // 2 ** 20
    except (ValueError, OSError):
        return None


def check_input(job) -> InputCheck:
    """Every finding `check` reads off a loaded job's input and machine.

    The geometry is the inline `* xyz` block or the file a `* xyzfile`
    names; any other form is skipped with a warning rather than guessed at.
    With a geometry, the elements are counted against the whole system's
    charge and multiplicity and every atom number the input names is checked
    against its size. Independently of the geometry, the process and memory
    requests are compared with this machine and an existing output is
    reported: a re-run would overwrite it."""
    inp = job.input
    text = job.inp_path.read_text(errors="replace")
    errors: list[Finding] = []
    warnings: list[Finding] = []

    atoms = None
    if inp.coords_kind == "xyz":
        atoms = list(inp.coords_atoms) or None
    elif inp.coords_kind == "xyzfile":
        if inp.coords_file:
            atoms = read_xyz(job.state.path / inp.coords_file)
        if atoms is None:
            errors.append(Finding("geometry_file",
                                  f"{inp.coords_file}: missing or not an XYZ file"))
    elif inp.coords_kind is None:
        errors.append(Finding("no_geometry", "the input has no * xyz or * xyzfile line"))
    else:
        warnings.append(Finding("geometry_not_read",
                                f"a {inp.coords_kind} geometry is not read: "
                                "the atom checks are skipped"))

    charge = mult = electrons = None
    numbers = None
    if atoms:
        numbers = []
        for i, atom in enumerate(atoms):
            z = ATOMIC_NUMBERS.get(atom[0])
            if z is None:
                errors.append(Finding("unknown_element", f"atom {i}: unknown element {atom[0]!r}"))
            else:
                numbers.append(z)
        if len(numbers) != len(atoms):
            numbers = None  # a partial sum would be a wrong electron count

    if numbers is not None:
        if inp.multilayer:
            charge, mult = inp.layers.get("total", (None, None))
            if charge is None or mult is None:
                warnings.append(Finding("total_unset",
                                        "Charge_Total and Mult_Total are not both set: "
                                        "the whole system's electrons are not checked"))
        else:
            charge, mult = inp.charge, inp.mult
        if charge is not None and mult is not None:
            electrons = sum(numbers) - charge
            if mult < 1 or mult - 1 > electrons:
                errors.append(Finding("parity", f"multiplicity {mult} is impossible "
                                                f"with {electrons} electrons"))
            elif (electrons - (mult - 1)) % 2 != 0:
                errors.append(Finding(
                    "parity", f"{electrons} electrons cannot have multiplicity {mult}: "
                              f"it needs an {'odd' if mult % 2 == 0 else 'even'} number"))

    if atoms:
        n = len(atoms)
        for name, indices in atom_references(text):
            first = next((i for i in indices if i >= n), None)
            if first is not None:
                errors.append(Finding(
                    "index_range", f"atom {first} in `{name}` is out of range: "
                                   f"the geometry has {n} atoms (0-{n - 1})"))

    qm = qm_atoms(text) if inp.multilayer and atoms else None
    if qm:
        for a, b in bonded_pairs(atoms):
            if (a in qm) == (b in qm):
                continue
            q, e = (a, b) if a in qm else (b, a)
            if atoms[q][0] == "H" or atoms[e][0] == "H":
                d = math.dist(atoms[q][1:], atoms[e][1:])
                warnings.append(Finding(
                    "qm_boundary_h", f"QM atom {q} {atoms[q][0]} is {d:.2f} Å from "
                                     f"environment atom {e} {atoms[e][0]}: "
                                     "a bond to hydrogen crosses the QM boundary"))

    cpus = os.cpu_count()
    if inp.nprocs is not None and cpus is not None and inp.nprocs > cpus:
        warnings.append(Finding("nprocs", f"nprocs {inp.nprocs} is more than "
                                          f"this machine's {cpus} CPUs"))
    if inp.maxcore_mb is not None:
        p = inp.nprocs or 1
        memory = _machine_memory_mb()
        if memory is not None and p * inp.maxcore_mb > memory:
            warnings.append(Finding(
                "memory", f"{p} process{'es' if p != 1 else ''} x {inp.maxcore_mb} MB maxcore = "
                          f"{p * inp.maxcore_mb} MB is more than this machine's {memory} MB"))
    if job.out_path.exists():
        warnings.append(Finding("output_exists", f"{job.state.stem}.out exists: "
                                                  "running this job again overwrites it"))

    return InputCheck(atoms=len(atoms) if atoms else None, electrons=electrons,
                      charge=charge, mult=mult, errors=errors, warnings=warnings)
