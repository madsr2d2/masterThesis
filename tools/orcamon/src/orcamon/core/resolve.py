"""A JOB argument on the command line -> the job it names.

Tried in order: an existing path (a job directory, or a `.inp`/`.out` file),
then the job's label under the root -- exactly, then as a case-sensitive
substring. Anything ambiguous is refused with the candidates rather than
guessed at: an agent that asked about `ts` and was answered about `ts5`
would report the wrong job's numbers with confidence."""
from __future__ import annotations

from pathlib import Path

from .discovery import JobRef, discover, label_for, stems_in

MAX_CANDIDATES = 20


class ResolveError(Exception):
    """The argument names no job, or more than one. The message says which
    and lists up to MAX_CANDIDATES of them."""


class _NoJobHere(ResolveError):
    """An existing directory with no job in it -- not yet a failure."""


def _candidates(labels: list[str]) -> str:
    shown = "\n".join(f"  {label}" for label in labels[:MAX_CANDIDATES])
    more = len(labels) - MAX_CANDIDATES
    return shown + (f"\n  ... and {more} more" if more > 0 else "")


def resolve(arg: str, root: Path, exclude=()) -> JobRef:
    root = root.resolve()
    path_error = None
    for base in (Path.cwd(), root):
        path = (base / arg) if not Path(arg).is_absolute() else Path(arg)
        if path.exists():
            try:
                return _from_path(path, root)
            except _NoJobHere as exc:
                # `ts` names a directory, but the job is `ts/final`: a label
                # match can still find it, so try that before refusing.
                path_error = str(exc)
                break

    refs = discover(root, exclude)
    exact = [r for r in refs if r.label == arg]
    if exact:
        return exact[0]
    matches = [r for r in refs if arg in r.label]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ResolveError(path_error or f"no job matches {arg!r} under {root}")
    raise ResolveError(
        f"{arg!r} matches {len(matches)} jobs under {root}; name one:\n"
        + _candidates([r.label for r in matches])
    )


def _from_path(path: Path, root: Path) -> JobRef:
    path = path.resolve()
    if path.is_file():
        if path.suffix not in (".inp", ".out"):
            raise ResolveError(f"{path} is not a .inp or .out file")
        return JobRef(path.parent, path.stem, label_for(root, path.parent, path.stem))
    stems = stems_in(path)
    if not stems:
        # A directory with an output but no input -- a copied result, say.
        stems = sorted(p.stem for p in path.glob("*.out"))
    if not stems:
        raise _NoJobHere(f"no ORCA job (.inp or .out) in {path}")
    if len(stems) > 1:
        raise ResolveError(
            f"{path} holds {len(stems)} jobs; name the file instead:\n"
            + _candidates([f"{path / s}.inp" for s in stems])
        )
    return JobRef(path, stems[0], label_for(root, path, stems[0]))
