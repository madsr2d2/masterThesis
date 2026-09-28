"""Units and the human formatting of times, shared by every renderer."""
from __future__ import annotations

EH_TO_KJ_PER_MOL = 2625.4996394799


def format_wall_time(seconds: float | None) -> str:
    if seconds is None:
        return "-"
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{s:02d}"


def format_age(seconds: float) -> str:
    """Minute resolution on purpose: the summary is only repainted when its
    text changes, and a seconds count would change it on every tick."""
    seconds = max(0, int(seconds))
    if seconds < 60:
        return "<1 min"
    if seconds < 90 * 60:
        return f"{seconds // 60} min"
    if seconds < 36 * 3600:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds / 86400:.1f} d"
