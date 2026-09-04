"""
Every gate in the repository, in one command.

    python run_gates.py              # the routine suite
    python run_gates.py --all        # including the slow optimiser suite
    python run_gates.py --only two_axis induction

CLAUDE.md listed nine commands and then said "and each analysis folder's own
check_numbers.py", which is six more. Fifteen commands run by hand is a suite
that gets run in part, and the part that gets skipped is the slow one.

THE GATES ARE DISCOVERED, NOT LISTED. A hardcoded list is a list that drifts:
a new `test_*.py` would be written, committed, and never run by anybody, which
is the same failure `test_curve_metrics.test_every_test_is_run` exists to catch
one level down. So this globs, and `test_the_runner_finds_every_gate` in
`data/test_scope.py` fails if the glob comes back with fewer gates than the
repository has files matching -- a widened glob that silently matched nothing
would otherwise pass by finding nothing to complain about.

Exit status is the point: non-zero if ANY gate fails, so this is usable as one
line in a pre-commit hook or a CI step.
"""
import argparse
import glob
import os
import subprocess
import sys
import time

REPOSITORY = os.path.dirname(os.path.abspath(__file__))

# Runs the optimiser rather than reading a fit back, so it is minutes where
# every other gate is seconds. Excluded from the routine suite and reachable
# with --all; `gate_paths` still asserts it EXISTS, so it cannot quietly go
# missing the way an unlisted gate would.
SLOW_GATES = ("data/test_fit_kinetics.py",)


def gate_paths(include_slow=False):
    """
    Every gate in the repository: the dataset validator, every `test_*.py`,
    and every folder's `check_numbers.py`. Relative to the repository root,
    in the order they are cheapest to fail in.
    """
    found = (["data/validate_dataset.py"]
             + sorted(_relative(glob.glob(os.path.join(REPOSITORY, "test_*.py"))))
             + sorted(_relative(glob.glob(os.path.join(REPOSITORY, "data", "test_*.py"))))
             + sorted(_relative(glob.glob(os.path.join(REPOSITORY, "*", "check_numbers.py")))))
    missing = [path for path in SLOW_GATES if path not in found]
    if missing:
        raise SystemExit(f"a declared slow gate does not exist: {missing}")
    if not include_slow:
        found = [path for path in found if path not in SLOW_GATES]
    return found


def _relative(paths):
    return [os.path.relpath(path, REPOSITORY) for path in paths]


def run_gate(path, argv=()):
    """One gate, its wall time and its exit status. Output is captured so a
    passing run stays one line; a failing one gets its tail printed."""
    started = time.time()
    done = subprocess.run([sys.executable, path, *argv], cwd=REPOSITORY,
                          capture_output=True, text=True)
    return time.time() - started, done.returncode, done.stdout + done.stderr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true",
                        help="include the slow optimiser suite")
    parser.add_argument("--only", nargs="+", metavar="NAME",
                        help="run only gates whose path contains one of these")
    args = parser.parse_args()

    gates = gate_paths(include_slow=args.all)
    if args.only:
        gates = [g for g in gates if any(name in g for name in args.only)]
        if not gates:
            raise SystemExit(f"--only {args.only} matched no gate")

    failed = []
    started = time.time()
    for path in gates:
        argv = ("--deep",) if path.endswith("validate_dataset.py") else ()
        elapsed, status, output = run_gate(path, argv)
        last = [line for line in output.strip().splitlines() if line.strip()]
        print(f"  {'pass' if not status else 'FAIL'}  {elapsed:6.1f}s  "
              f"{path:<38} {last[-1].strip() if last else ''}")
        if status:
            failed.append(path)
            print("\n".join("        " + line for line in last[-25:]) + "\n")

    print(f"\n{len(gates)} gates in {time.time() - started:.0f}s, "
          f"{len(failed)} failed" + (f": {', '.join(failed)}" if failed else ""))
    if not args.all:
        print(f"(--all adds {', '.join(SLOW_GATES)})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
