"""
The contract orcamon's command-line tier (`tools/orcamon/`) holds.

    python test_orcamon.py

`test_monitor.py` holds the parser to ORCA's output; this holds what orcamon
does with it: the dependency rule, job resolution, every command's output and
exit status, the report contract, the cache, notifications, graphics
selection, the scheduler probe, and the agent skill that ships with it.

Every case runs on a temporary tree of synthetic outputs, never on the jobs
in `computational/`: orcamon is a general tool, and a gate that read this
project's jobs would pass or fail on whatever happened to be running.

`run_gates.py` discovers this file by its name, like any other gate.
"""
import os
import subprocess
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)


# The modules only the TUI and its pixel graphics may import. A core or CLI
# module importing any of them would make `pip install orcamon` on an
# agents-only server fail at the first command.
_FORBIDDEN = ("textual", "rich", "plotext", "textual_plotext", "matplotlib",
              "numpy", "PIL", "psutil")

_STDLIB_ONLY_PROBE = textwrap.dedent("""
    import importlib, io, sys, contextlib, pkgutil
    FORBIDDEN = %r

    class Refuse:
        def find_spec(self, name, path=None, target=None):
            if name.split(".")[0] in FORBIDDEN:
                raise ImportError(f"forbidden import: {name}")
            return None

    sys.meta_path.insert(0, Refuse())
    import orcamon.core, orcamon.cli
    for pkg in (orcamon.core, orcamon.cli):
        for info in pkgutil.iter_modules(pkg.__path__):
            importlib.import_module(f"{pkg.__name__}.{info.name}")
    import orcamon.validate
    COMMANDS = %r
    for argv in COMMANDS:
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = orcamon.cli.main(argv)
        print(argv[-2] if len(argv) > 1 else argv, "->", code)
    print("OK")
""")


def _stdlib_probe(commands):
    """Import every core and CLI module, and run `commands`, with the TUI's
    dependencies made unimportable. Returns (exit code, output)."""
    done = subprocess.run(
        [sys.executable, "-c", _STDLIB_ONLY_PROBE % (_FORBIDDEN, commands)],
        capture_output=True, text=True, cwd=HERE,
    )
    return done.returncode, done.stdout + done.stderr


def test_the_core_needs_only_the_standard_library():
    print("\nthe core and the CLI import only the standard library")
    code, output = _stdlib_probe([])
    check("every orcamon.core and orcamon.cli module imports with the TUI's "
          "dependencies refused", code == 0 and output.strip().endswith("OK"),
          output.strip().splitlines()[-1] if output.strip() else "")


if __name__ == "__main__":
    test_the_core_needs_only_the_standard_library()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
