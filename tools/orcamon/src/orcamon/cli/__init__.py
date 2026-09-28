"""`orcamon` -- the command line. Standard library only.

    orcamon [ROOT]            the TUI (needs the `tui` extra)
    orcamon tui [ROOT]        the same
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

TUI_MISSING = "the TUI needs the tui extra: uv tool install 'orcamon[tui]'"


def _tui(root: str) -> int:
    try:
        from ..tui.app import run
    except ImportError as exc:
        if exc.name and exc.name.split(".")[0] in {"textual", "textual_plotext", "plotext", "rich"}:
            print(TUI_MISSING, file=sys.stderr)
            return 2
        raise
    run(Path(root))
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(prog="orcamon", description="Watch and query ORCA jobs.")
    sub = parser.add_subparsers(dest="command")
    tui = sub.add_parser("tui", help="the terminal UI")
    tui.add_argument("root", nargs="?", default=".", help="directory to watch (default: .)")

    # `orcamon [ROOT]` is the TUI: anything that is not a command is a root.
    if not argv or (argv[0] not in sub.choices and not argv[0].startswith("-")):
        argv = ["tui", *argv]
    args = parser.parse_args(argv)
    if args.command == "tui":
        return _tui(args.root)
    parser.print_help()
    return 2
