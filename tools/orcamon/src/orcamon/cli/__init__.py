"""`orcamon` -- the command line. Standard library only.

    orcamon [ROOT]                 the TUI (needs the `tui` extra)
    orcamon tui [ROOT]             the same
    orcamon <command> [args]       one bounded answer about a job, for agents,
                                   scripts and people in a shell

Every command takes `--json` and prints one JSON document; without it, plain
text of at most `--max-lines` lines. Exit status: 0 success, 1 the job (or,
for `ls`, any listed job) is failed or stopped, 2 usage error or no such
job, 4 `wait` timed out, 130 interrupted.
"""
from __future__ import annotations

import argparse
import os
import sys

from ..core import liveness
from ..core.geometry import REPRESENTATIONS
from ..core.status import QUIET_AFTER_S
from . import commands, skill

TUI_MISSING = "the TUI needs the tui extra: uv tool install 'orcamon[tui]'"

EXIT_OK, EXIT_JOB_BAD, EXIT_USAGE, EXIT_TIMEOUT, EXIT_INTERRUPTED = 0, 1, 2, 4, 130
EXIT_CODES = (
    (EXIT_OK, "success; for show and wait, the job is not failed or stopped"),
    (EXIT_JOB_BAD, "the command worked and the job (for ls: any listed job) is failed or stopped"),
    (EXIT_USAGE, "usage error, or a job not found or ambiguous"),
    (EXIT_TIMEOUT, "wait timed out"),
    (EXIT_INTERRUPTED, "interrupted (Ctrl-C)"),
)

DEFAULT_MAX_LINES = 60


_Formatter = argparse.RawDescriptionHelpFormatter

# Options (on the main parser and every subcommand) that take a value, so
# that value is not mistaken for the first positional word.
_VALUED_OPTIONS = {"--root", "--liveness", "--quiet-after"}


def _common(parser: argparse.ArgumentParser, suppress: bool) -> None:
    """Options every command accepts, before or after its name. The copy on
    each subcommand defaults to SUPPRESS, so an option given before the
    command is not overwritten by the subcommand's default."""
    d = (lambda v: argparse.SUPPRESS) if suppress else (lambda v: v)
    parser.add_argument("--root", default=d(os.environ.get("ORCAMON_ROOT", ".")), metavar="DIR",
                        help="where jobs are looked for (default: $ORCAMON_ROOT or .)")
    parser.add_argument("--liveness", choices=("auto", *liveness.PROBES), default=d("auto"),
                        help="how a job is known to be running (default: auto)")
    parser.add_argument("--quiet-after", type=float, default=d(QUIET_AFTER_S), metavar="S",
                        help=f"seconds without output before a job of unknown liveness "
                             f"is called quiet (default: {QUIET_AFTER_S:g})")
    parser.add_argument("--no-cache", action="store_true", default=d(False),
                        help="parse every output from scratch; neither read nor write the state cache")


def _output(parser: argparse.ArgumentParser, lines: bool = True) -> None:
    parser.add_argument("--json", action="store_true", help="print one JSON document instead of text")
    if lines:
        parser.add_argument("--max-lines", type=int, default=DEFAULT_MAX_LINES, metavar="N",
                            help="cap on text output lines (default: %(default)s)")


def _job(parser: argparse.ArgumentParser, many: bool = False) -> None:
    parser.add_argument("job", nargs="+" if many else None, metavar="JOB",
                        help="a job directory, a .inp/.out file, or (part of) a job's label")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="orcamon", description=__doc__, formatter_class=_Formatter,
    )
    parser.add_argument("--version", action="version", version=f"orcamon {commands.__version__}")
    _common(parser, suppress=False)
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    def add(name, help_, example, func):
        p = sub.add_parser(name, help=help_, description=f"{help_}\n\nexample: {example}",
                           formatter_class=_Formatter)
        _common(p, suppress=True)
        p.set_defaults(func=func, example=example)
        return p

    p = add("tui", "the terminal UI: every job under ROOT, live", "orcamon tui runs/", commands.cmd_tui)
    p.add_argument("tui_root", nargs="?", metavar="ROOT", help="directory to watch (default: --root)")
    p.add_argument("--exclude", action="append", default=[], metavar="GLOB",
                   help="skip jobs and directories matching GLOB (repeatable)")
    commands.add_tui_options(p)

    p = add("ls", "list every job under the root with its status and flags", "orcamon ls --attention", commands.cmd_ls)
    p.add_argument("--status", metavar="S[,S]", help="only jobs with these statuses")
    p.add_argument("--attention", action="store_true", help="only jobs with attention flags")
    p.add_argument("--since", metavar="DUR", help="only jobs whose output changed within DUR (30m, 2h, 1d)")
    p.add_argument("--sort", choices=("path", "age", "status"), default="path", help="row order (default: path)")
    p.add_argument("--exclude", action="append", default=[], metavar="GLOB",
                   help="skip jobs and directories matching GLOB (repeatable)")
    _output(p)

    p = add("show", "one job's summary: identity, status, criteria, energy, flags", "orcamon show opt/ts", commands.cmd_show)
    _job(p)
    _output(p)

    p = add("conv", "the geometry convergence table, one row per cycle", "orcamon conv opt/ts --last 20", commands.cmd_conv)
    _job(p)
    p.add_argument("--last", type=int, default=10, metavar="N", help="cycles to show (default: %(default)s)")
    _output(p)

    p = add("energies", "the energy of every geometry, or of every scan step", "orcamon energies scan", commands.cmd_energies)
    _job(p)
    p.add_argument("--last", type=int, metavar="N",
                   help="rows to show (default: 50 for an optimization, every step for a scan)")
    _output(p)

    p = add("geom", "a geometry as XYZ: the latest, or a given cycle or scan step", "orcamon geom opt/ts --cycle 12 > ts12.xyz", commands.cmd_geom)
    _job(p)
    p.add_argument("--cycle", type=int, metavar="C", help="optimization cycle (within --step for a scan)")
    p.add_argument("--step", type=int, metavar="S", help="scan step (its last cycle unless --cycle)")
    p.add_argument("--region", choices=("all", "qm"), default="all",
                   help="qm: only a multilayer job's high-level (QM1) atoms (default: all)")
    _output(p, lines=False)

    p = add("snapshot", "render a job's geometry to a PNG (needs the images extra)",
            "orcamon snapshot opt/ts -o ts.png", commands.cmd_snapshot)
    _job(p)
    p.add_argument("-o", "--output", required=True, metavar="FILE", help="where to write the PNG")
    p.add_argument("--representation", choices=REPRESENTATIONS, default="ball-and-stick",
                   help="how to draw the high-level region (default: %(default)s)")
    p.add_argument("--size", default="900x750", metavar="WxH",
                   help="image size in pixels (default: %(default)s)")
    p.add_argument("--elev", type=float, default=20, metavar="DEG",
                   help="elevation angle (default: %(default)s)")
    p.add_argument("--azim", type=float, default=-60, metavar="DEG",
                   help="azimuth angle (default: %(default)s)")
    p.add_argument("--no-labels", action="store_true", help="omit the atom index labels")
    p.add_argument("--no-fog", action="store_true", help="do not fade distant atoms")
    p.add_argument("--distances", action="store_true", help="label the QM-QM bond distances")

    p = add("freqs", "imaginary and lowest real vibrational frequencies", "orcamon freqs opt/ts", commands.cmd_freqs)
    _job(p)
    p.add_argument("--lowest", type=int, default=6, metavar="K", help="real modes to list (default: %(default)s)")
    _output(p)

    p = add("input", "what the input asks for: run types, method, charge and multiplicity", "orcamon input opt/ts", commands.cmd_input)
    _job(p)
    _output(p)

    p = add("errors", "ORCA error and QM2-error lines with context", "orcamon errors opt/ts --context 5", commands.cmd_errors)
    _job(p)
    p.add_argument("--context", type=int, default=3, metavar="N", help="lines around each match (default: %(default)s)")
    _output(p)

    p = add("tail", "the last lines of the output, optionally filtered", "orcamon tail opt/ts --grep 'FINAL SINGLE'", commands.cmd_tail)
    _job(p)
    p.add_argument("-n", type=int, default=40, metavar="N", dest="n",
                   help=f"lines (default: %(default)s, at most {commands.TAIL_CAP})")
    p.add_argument("--grep", metavar="REGEX",
                   help=f"only lines matching REGEX among the last {commands.TAIL_GREP_WINDOW}")
    _output(p, lines=False)

    p = add("wait", "block until a job meets a condition, then show it", "orcamon wait opt/ts --until done", commands.cmd_wait)
    _job(p, many=True)
    p.add_argument("--until", required=True, metavar="COND",
                   help="done | finished | change | attention | cycle>=N | step>=N")
    p.add_argument("--all", action="store_true", dest="all_jobs", help="wait for every JOB, not the first")
    p.add_argument("--interval", type=float, default=15.0, metavar="S", help="seconds between polls (default: %(default)s)")
    p.add_argument("--timeout", type=float, default=540.0, metavar="S",
                   help="give up after S seconds, exit 4 (default: %(default)s, under a 10-minute tool ceiling)")
    _output(p)

    skill.add_parser(add)
    return parser


def subcommands(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    """Subcommand name -> its parser, in registration order."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return dict(action.choices)
    return {}


def _first_positional(argv: list[str]) -> int | None:
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--":
            return i + 1 if i + 1 < len(argv) else None
        if a in _VALUED_OPTIONS:
            i += 2
            continue
        if not a.startswith("-"):
            return i
        i += 1
    return None


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    known = set(subcommands(parser))
    # `orcamon [ROOT]` is the TUI: a first positional word that is not a
    # command is a root.
    first = _first_positional(argv)
    if first is None and not any(a in ("-h", "--help", "--version") for a in argv):
        argv = [*argv, "tui"]
    elif first is not None and argv[first] not in known:
        argv = [*argv[:first], "tui", *argv[first:]]
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code or 0)
    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return EXIT_INTERRUPTED
    except BrokenPipeError:
        # `orcamon tail ... | head`: the reader left, which is not an error.
        try:
            sys.stdout = open(os.devnull, "w")
        except OSError:
            pass
        return EXIT_OK
