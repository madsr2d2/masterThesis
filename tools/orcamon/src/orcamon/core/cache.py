"""Parsed job state kept on disk, so a finished job's output is read once.

This tree's 64 jobs are 667 MB of output and a full parse takes seconds;
without a cache every `orcamon ls` and every TUI start paid that again, for
outputs that had not changed in days. An entry is the pickled `JobState` at
some offset into the file, and `update_job` then reads only what was appended
since -- the same incremental read the TUI does between ticks.

An entry is trusted only when all of these hold, and is otherwise discarded
and the output parsed from zero:

- `CACHE_KEY` matches. It is the hash of `parser.py`'s own SOURCE, so any
  change to the parser invalidates every entry by itself; there is no version
  number for anyone to forget to bump.
- the `.out` is the same file (inode) and at least as long as the entry's
  offset -- a re-run replaces or truncates it.

Any failure to read, unpickle or check an entry means "no entry", never an
error: the cache can only make orcamon faster, not make it fail.
"""
from __future__ import annotations

import hashlib
import os
import pickle
import tempfile
from pathlib import Path

from . import parser
from .parser import JobState


def _parser_key() -> str:
    try:
        source = Path(parser.__file__).read_bytes()
    except OSError:
        source = b""
    return hashlib.sha1(source).hexdigest()


CACHE_KEY = _parser_key()


def cache_dir() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return Path(base) / "orcamon"


def entry_path(out_path: Path) -> Path:
    return cache_dir() / (hashlib.sha1(str(out_path.resolve()).encode()).hexdigest() + ".pickle")


def restore(state: JobState) -> bool:
    """Replace `state`'s contents with the cached entry for its output, when
    one exists and is still valid. Returns whether it did."""
    out = state.path / f"{state.stem}.out"
    try:
        st = out.stat()
        with open(entry_path(out), "rb") as f:
            key, cached = pickle.load(f)
        if key != CACHE_KEY or not isinstance(cached, JobState):
            return False
        if cached.inode != st.st_ino or st.st_size < cached.offset:
            return False
    except Exception:  # noqa: BLE001 -- any bad entry is simply no entry
        return False
    # Keep the caller's identity (path as it resolved THIS time), take the
    # rest from the entry.
    cached.path, cached.stem = state.path, state.stem
    state.__dict__.update(cached.__dict__)
    return True


def save(state: JobState) -> bool:
    """Write `state` for its output. Skipped when there is nothing to save,
    and when the parser is inside a block: the partial block would be
    restored as if it were complete-so-far state, which it is not -- the next
    call reparses from the last good entry or from zero instead."""
    if not state.has_out or state.offset == 0 or state._in_block():
        return False
    target = entry_path(state.path / f"{state.stem}.out")
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=target.parent, prefix=".tmp-")
        try:
            with os.fdopen(fd, "wb") as f:
                pickle.dump((CACHE_KEY, state), f, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(tmp, target)
        except BaseException:
            os.unlink(tmp)
            raise
    except OSError:
        return False
    return True
