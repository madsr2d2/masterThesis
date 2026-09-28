"""`orcamon skill`: the agent skill that ships inside the package.

The skill describes commands, flags, exit codes and JSON keys. Kept anywhere
else it goes stale the first time a command changes, so it lives in the
package (`SKILL.md`, package data) and is released with the commands. Only
the judgement calls are written by hand; the command reference is generated
from the argparse parsers every time it is printed, so a flag cannot be
documented that does not exist, or exist undocumented.

    orcamon skill                          print it
    orcamon skill install [--user | --project DIR] [--force]
    orcamon skill --check [--user | --project DIR]
"""
from __future__ import annotations

import argparse
import re
import sys
from importlib import resources
from pathlib import Path

from .. import __version__

SKILL_NAME = "orcamon"
VERSION_KEY = "orcamon_version"
_FRONTMATTER_RE = re.compile(r"\A---\n(.*?)\n---\n", re.S)
_VERSION_RE = re.compile(rf"^\s*{VERSION_KEY}:\s*(\S+)\s*$", re.M)


def body() -> str:
    """The hand-written skill: frontmatter and body."""
    return resources.files("orcamon").joinpath("SKILL.md").read_text(encoding="utf-8")


def _usage(action: argparse.Action) -> str:
    if not action.option_strings:
        return action.metavar or action.dest.upper()
    flag = ", ".join(action.option_strings)
    if action.nargs == 0:
        return flag
    if action.choices and not action.metavar:
        return f"{flag} {'|'.join(str(c) for c in action.choices)}"
    return f"{flag} {action.metavar or action.dest.upper()}"


def _help(action: argparse.Action) -> str:
    return (action.help or "").replace("%(default)s", str(action.default))


def reference(parser: argparse.ArgumentParser) -> str:
    """The command reference, from the parser: every subcommand, its purpose,
    its arguments and flags, then the exit codes. Compact on purpose -- one
    line per flag -- because it is read by a model with a context budget."""
    from . import EXIT_CODES, subcommands

    lines = ["## Command reference", "",
             "Generated from the installed orcamon's own parser; the flags below are the ones that exist.",
             "", "Global options (before or after the command):", ""]
    for action in parser._actions:
        if action.option_strings and action.dest not in ("help", "version"):
            lines.append(f"- `{_usage(action)}`: {_help(action)}")
    for name, sub in subcommands(parser).items():
        lines += ["", f"### orcamon {name}", "", (sub.description or "").split("\n")[0], ""]
        for action in sub._actions:
            if action.dest in ("help", "root", "liveness", "quiet_after", "no_cache"):
                continue
            if isinstance(action, argparse._SubParsersAction):
                for sub_name, subsub in action.choices.items():
                    lines.append(f"- `{name} {sub_name}`: {(subsub.description or '').splitlines()[0]}")
                continue
            lines.append(f"- `{_usage(action)}`: {_help(action)}")
    lines += ["", "### Exit codes", ""]
    lines += [f"- `{code}`: {meaning}" for code, meaning in EXIT_CODES]
    return "\n".join(lines) + "\n"


def full_text(parser: argparse.ArgumentParser) -> str:
    return body().rstrip("\n") + "\n\n" + reference(parser)


def stamped(text: str) -> str:
    """The text with `metadata: {orcamon_version: X}` in its frontmatter --
    what `install` writes and `--check` compares against."""
    m = _FRONTMATTER_RE.match(text)
    if not m:
        raise ValueError("SKILL.md has no frontmatter")
    front = m.group(1).rstrip("\n") + f"\nmetadata:\n  {VERSION_KEY}: {__version__}"
    return f"---\n{front}\n---\n" + text[m.end():]


def target(args) -> Path:
    base = Path(args.project).expanduser() if args.project else Path.home()
    return base / ".claude" / "skills" / SKILL_NAME / "SKILL.md"


def cmd_skill(args) -> int:
    from . import build_parser

    parser = build_parser()
    if args.check:
        return _check(args, parser)
    if args.action is None:
        sys.stdout.write(full_text(parser))
        return 0
    if args.action == "install":
        return _install(args, parser)
    print(f"orcamon: unknown skill action {args.action!r}", file=sys.stderr)
    return 2


def _install(args, parser) -> int:
    path = target(args)
    if path.exists() and not args.force:
        try:
            current = path.read_text(encoding="utf-8")
        except OSError as exc:
            print(f"orcamon: cannot read {path}: {exc}", file=sys.stderr)
            return 2
        if not _VERSION_RE.search(_front(current)):
            print(f"orcamon: {path} exists and was not written by orcamon skill install; "
                  "--force to replace it", file=sys.stderr)
            return 2
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(stamped(full_text(parser)), encoding="utf-8")
    print(path)
    return 0


def _front(text: str) -> str:
    m = _FRONTMATTER_RE.match(text)
    return m.group(1) if m else ""


def _check(args, parser) -> int:
    path = target(args)
    try:
        installed = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        print(f"no orcamon skill installed at {path}: run `orcamon skill install`")
        return 2
    except OSError as exc:
        print(f"cannot read {path}: {exc}")
        return 2
    m = _VERSION_RE.search(_front(installed))
    if not m:
        print(f"{path} was not written by orcamon skill install")
        return 1
    if m.group(1) != __version__:
        print(f"installed skill is {m.group(1)}, orcamon is {__version__}: run `orcamon skill install`")
        return 1
    if installed != stamped(full_text(parser)):
        print(f"{path} differs from what orcamon {__version__} prints (edited, or a stale build): "
              "run `orcamon skill install`")
        return 1
    print(f"{path} is current (orcamon {__version__})")
    return 0


def add_parser(add) -> None:
    p = add("skill", "print the agent skill, or install it for Claude Code",
            "orcamon skill install --project .", cmd_skill)
    p.add_argument("action", nargs="?", choices=("install",), metavar="install",
                   help="write the skill to .claude/skills/orcamon/SKILL.md")
    p.add_argument("--check", action="store_true",
                   help="exit 0 if the installed skill is current, 1 if stale or edited, 2 if absent")
    where = p.add_mutually_exclusive_group()
    where.add_argument("--user", action="store_true", help="under ~/.claude (the default)")
    where.add_argument("--project", metavar="DIR", help="under DIR/.claude instead")
    p.add_argument("--force", action="store_true", help="replace a SKILL.md orcamon did not write")
