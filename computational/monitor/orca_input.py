"""What a job IS, read off its ORCA input rather than its output.

The summary leads with this because it is known before the job runs: a wrong
charge or multiplicity is cheapest to catch while the job is still "not run",
and nothing in the output can be trusted to have been written yet. The
reading is deliberately tolerant -- it pulls out a handful of facts by pattern
and never tries to validate the input, which ORCA itself does on startup.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# Simple-input keywords that say what the job DOES rather than how it models
# the system. Compared case-insensitively; everything else on the `!` lines is
# shown as the method.
RUN_TYPE_KEYWORDS = {
    k.lower()
    for k in (
        "SP", "ENGRAD", "NUMGRAD",
        "OPT", "COPT", "ZOPT", "GDIIS-OPT",
        "LOOSEOPT", "NORMALOPT", "TIGHTOPT", "VERYTIGHTOPT",
        "OPTTS", "SCANTS",
        "FREQ", "NUMFREQ", "ANFREQ",
        "IRC", "NEB", "NEB-CI", "NEB-TS", "ZOOM-NEB", "ZOOM-NEB-CI", "ZOOM-NEB-TS",
        "MD", "GOAT",
    )
}

MULTIPLICITY_NAMES = {
    1: "singlet", 2: "doublet", 3: "triplet", 4: "quartet", 5: "quintet",
    6: "sextet", 7: "septet",
}

_COMMENT_RE = re.compile(r"#.*$", re.M)
# `* xyz 0 1`, `*xyzfile -1 2 geom.xyz`, `* int 0 1`, `* pdbfile 0 1 x.pdb` ...
_COORDS_RE = re.compile(
    r"^\s*\*\s*(xyzfile|xyz|internal|int|gzmtfile|gzmt|pdbfile)\s+(-?\d+)\s+(\d+)",
    re.I | re.M,
)
# The %coords block form: `Charge 0` / `Mult 1` on lines of their own.
_COORDS_BLOCK_RE = re.compile(r"%coords\b", re.I)
_BLOCK_CHARGE_RE = re.compile(r"^\s*charge\s+(-?\d+)\s*$", re.I | re.M)
_BLOCK_MULT_RE = re.compile(r"^\s*mult\s+(\d+)\s*$", re.I | re.M)
# Multilayer (%qmmm: QM/MM, QM/QM2, ONIOM) layer settings: Charge_Total,
# Mult_Total, Charge_Medium, ...
_LAYER_CHARGE_RE = re.compile(r"\bcharge_(\w+)\s+(-?\d+)", re.I)
_LAYER_MULT_RE = re.compile(r"\bmult_(\w+)\s+(\d+)", re.I)
_QMMM_RE = re.compile(r"%qmmm\b", re.I)
_SCAN_RE = re.compile(r"%geom\b.*?^\s*scan\b", re.I | re.M | re.S)
_NPROCS_RE = re.compile(r"\bnprocs\s+(\d+)", re.I)
_PAL_RE = re.compile(r"^pal(\d+)$", re.I)
_MAXCORE_RE = re.compile(r"%maxcore\s+(\d+)", re.I)


@dataclass
class JobInput:
    keywords: list[str] = field(default_factory=list)
    run_types: list[str] = field(default_factory=list)
    # The `*` line's (or %coords') charge and multiplicity. On a multilayer
    # job this is the HIGH-LEVEL region's, and the whole system's is in
    # `layers["total"]`.
    charge: int | None = None
    mult: int | None = None
    multilayer: bool = False
    # layer name (lowercase, e.g. "total", "medium") -> (charge, mult); either
    # half may be None when only one of the pair is set.
    layers: dict[str, tuple[int | None, int | None]] = field(default_factory=dict)
    nprocs: int | None = None
    maxcore_mb: int | None = None

    @property
    def method(self) -> str:
        return " ".join(k for k in self.keywords if k.lower() not in RUN_TYPE_KEYWORDS
                        and not _PAL_RE.match(k))


def read_input(path: Path) -> JobInput | None:
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return None
    return parse_input(text)


def parse_input(text: str) -> JobInput:
    text = _COMMENT_RE.sub("", text)
    job = JobInput()

    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("!"):
            job.keywords.extend(stripped[1:].split())
    for k in job.keywords:
        if k.lower() in RUN_TYPE_KEYWORDS:
            job.run_types.append(k)
        m = _PAL_RE.match(k)
        if m:
            job.nprocs = int(m.group(1))
    if _SCAN_RE.search(text):
        job.run_types.append("Scan")

    m = _COORDS_RE.search(text)
    if m:
        job.charge, job.mult = int(m.group(2)), int(m.group(3))
    elif _COORDS_BLOCK_RE.search(text):
        c, mu = _BLOCK_CHARGE_RE.search(text), _BLOCK_MULT_RE.search(text)
        job.charge = int(c.group(1)) if c else None
        job.mult = int(mu.group(1)) if mu else None

    job.multilayer = bool(_QMMM_RE.search(text))
    if job.multilayer:
        for m in _LAYER_CHARGE_RE.finditer(text):
            name = m.group(1).lower()
            job.layers[name] = (int(m.group(2)), job.layers.get(name, (None, None))[1])
        for m in _LAYER_MULT_RE.finditer(text):
            name = m.group(1).lower()
            job.layers[name] = (job.layers.get(name, (None, None))[0], int(m.group(2)))

    m = _NPROCS_RE.search(text)
    if m:
        job.nprocs = int(m.group(1))
    m = _MAXCORE_RE.search(text)
    if m:
        job.maxcore_mb = int(m.group(1))
    return job


def describe_spin(charge: int | None, mult: int | None) -> str:
    """`charge 0 · singlet`, tolerating either half being unknown."""
    c = "charge ?" if charge is None else f"charge {charge:+d}" if charge else "charge 0"
    if mult is None:
        return f"{c} · mult ?"
    return f"{c} · {MULTIPLICITY_NAMES.get(mult, f'mult {mult}')}"
