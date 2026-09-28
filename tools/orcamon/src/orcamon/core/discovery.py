"""Which ORCA jobs live under a directory, and what each is called.

A job is a `<stem>.inp`. It used to be a DIRECTORY -- one stem per directory,
`job` if present, otherwise the first alphabetically -- so a directory holding
`opt.inp` and `freq.inp` showed one of them and silently hid the other."""
from __future__ import annotations

import fnmatch
import os
import re
from dataclasses import dataclass
from pathlib import Path

# Inputs ORCA writes for itself inside a real job's directory, which are not
# jobs of their own: the numerical Hessian's displaced-geometry gradients
# (job_D00159.scfgrad.inp) and the plain gradient helper (job.scfgrad.inp).
_GENERATED_INP_RE = re.compile(r"(?:_D\d+)?\.scfgrad\.inp$")


@dataclass(frozen=True)
class JobRef:
    path: Path   # the job's directory, resolved
    stem: str
    label: str   # what the job is called: its directory relative to the root,
                 # plus `/stem` when the directory holds more than one job


def is_generated_input(name: str) -> bool:
    return bool(_GENERATED_INP_RE.search(name))


def stems_in(job_dir: Path) -> list[str]:
    """Every job stem in one directory, sorted."""
    try:
        names = os.listdir(job_dir)
    except OSError:
        return []
    return sorted(n[:-4] for n in names if n.endswith(".inp") and not is_generated_input(n))


def discover(root: Path, exclude: tuple[str, ...] | list[str] = ()) -> list[JobRef]:
    """Every job under `root`, sorted by label.

    Hidden directories are skipped (a `.git` or a `.venv` holds no job
    worth watching and can be large). `exclude` takes shell globs, matched
    against each job's label and against every directory's path relative to
    the root -- a matched directory is not descended into."""
    root = root.resolve()
    by_dir: dict[Path, list[str]] = {}
    for dirpath, dirnames, filenames in os.walk(root):
        here = Path(dirpath)
        rel = _rel(root, here)
        dirnames[:] = sorted(
            d for d in dirnames
            if not d.startswith(".") and not _excluded(_join(rel, d), exclude)
        )
        stems = sorted(
            f[:-4] for f in filenames if f.endswith(".inp") and not is_generated_input(f)
        )
        if stems:
            by_dir[here] = stems

    refs = []
    for job_dir, stems in by_dir.items():
        rel = _rel(root, job_dir)
        for stem in stems:
            label = rel if len(stems) == 1 else _join(rel, stem)
            if not _excluded(label, exclude):
                refs.append(JobRef(job_dir, stem, label))
    return sorted(refs, key=lambda r: r.label)


def label_for(root: Path, job_dir: Path, stem: str) -> str:
    """The label `discover` would give this job -- or, for a directory
    outside the root, its path as it stands."""
    job_dir = job_dir.resolve()
    try:
        rel = str(job_dir.relative_to(root.resolve())) or "."
    except ValueError:
        rel = str(job_dir)
    return rel if len(stems_in(job_dir)) <= 1 else _join(rel, stem)


def _rel(root: Path, path: Path) -> str:
    rel = str(path.relative_to(root))
    return "." if rel in ("", ".") else rel


def _join(rel: str, name: str) -> str:
    return name if rel == "." else f"{rel}/{name}"


def _excluded(path: str, patterns) -> bool:
    return any(fnmatch.fnmatchcase(path, p) or fnmatch.fnmatchcase(path.rsplit("/", 1)[-1], p)
               for p in patterns)


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
