from __future__ import annotations

from pathlib import Path

import psutil


def running_orca_cwds() -> dict[Path, float]:
    by_cwd: dict[Path, list[tuple[str, float]]] = {}
    for proc in psutil.process_iter(["name", "cwd", "create_time"]):
        try:
            name = proc.info["name"] or ""
            if not name.startswith("orca"):
                continue
            cwd = proc.info["cwd"]
            if cwd is None:
                continue
            resolved = Path(cwd).resolve()
            by_cwd.setdefault(resolved, []).append((name, proc.info["create_time"]))
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            continue

    result: dict[Path, float] = {}
    for cwd, entries in by_cwd.items():
        main = [t for name, t in entries if name == "orca"]
        result[cwd] = min(main) if main else min(t for _, t in entries)
    return result
