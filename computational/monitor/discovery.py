from __future__ import annotations

import re
from pathlib import Path

# Numerical-Hessian displacement scratch inputs (job_D00159.scfgrad.inp, ...)
# live inside a real job's own directory and are not jobs of their own.
_SCRATCH_INP_RE = re.compile(r"_D\d+\.scfgrad\.inp$")


def discover_jobs(root: Path) -> list[tuple[Path, str]]:
    root = root.resolve()
    by_dir: dict[Path, list[str]] = {}
    for inp in root.rglob("*.inp"):
        if _SCRATCH_INP_RE.search(inp.name):
            continue
        by_dir.setdefault(inp.parent, []).append(inp.stem)

    jobs = []
    for job_dir, stems in by_dir.items():
        stem = "job" if "job" in stems else sorted(stems)[0]
        jobs.append((job_dir, stem))
    return sorted(jobs)


def relative_label(root: Path, job_dir: Path) -> str:
    return str(job_dir.resolve().relative_to(root.resolve()))


def short_label(label: str, width: int) -> str:
    """Keep the distinguishing END of a job path (plus a head segment for
    reaction context), not naive left-to-right truncation."""
    if len(label) <= width:
        return label
    parts = label.split("/")
    # the reaction/molecule folder, one level in -- distinguishes two
    # reactions whose stage names (rc/ts/pc/attempt*) otherwise match
    head = parts[1] if len(parts) > 2 else parts[0]

    def fit_tail(budget: int) -> list[str]:
        tail_parts: list[str] = []
        for part in reversed(parts[1:] if len(parts) > 2 else parts):
            candidate = "/".join([part, *tail_parts])
            if len(candidate) > budget:
                break
            tail_parts.insert(0, part)
        return tail_parts

    tail_parts = fit_tail(width - len(head) - 2)  # "/…" separator
    if not tail_parts:
        tail_parts = fit_tail(width - 2)  # drop the head, tail alone
        if not tail_parts:
            return "…" + parts[-1][-(width - 1):]
        return "…/" + "/".join(tail_parts)
    return f"{head}/…/{'/'.join(tail_parts)}"
