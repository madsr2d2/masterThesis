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
import tempfile
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# The TUI tests run headless. Inheriting a herdr pane would draw the
# molecule over the pane this runs in.
os.environ["HERDR_ENV"] = "0"

# The TUI publishes what it shows under $XDG_STATE_HOME/orcamon, and `view`
# and `focus` read it there: a run must never see -- or steer -- the
# person's own TUI.
_STATE_HOME = tempfile.TemporaryDirectory(prefix="orcamon_state_")
os.environ["XDG_STATE_HOME"] = _STATE_HOME.name

FAILURES = []


def check(label, ok, detail=""):
    # The detail is for reading a FAILURE; on a pass it would bury the list.
    print(f"  {'pass' if ok else 'FAIL'}  {label}" + (f": {detail}" if detail and not ok else ""))
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
    with _Tree() as t:
        root = ["--root", str(t.root), "--liveness", "process"]
        commands = [
            [*root, "ls"], [*root, "ls", "--json"],
            [*root, "show", "opt_done"], [*root, "show", "opt_done", "--json"],
            [*root, "wait", "opt_done", "--until", "done", "--timeout", "5"],
        ]
        code, output = _stdlib_probe(commands)
        ran = [line for line in output.splitlines() if "->" in line]
        check("and ls, show and wait run with them refused",
              code == 0 and output.strip().endswith("OK") and len(ran) == len(commands),
              output.strip()[-400:])




# --- the synthetic tree ------------------------------------------------------
#
# One of each kind of job an agent asks about. Built fresh per run in a
# temporary directory; the state cache is pointed there too, so a run never
# reads or writes the user's own ~/.cache.

_CYCLE = "         *                GEOMETRY OPTIMIZATION CYCLE  {n:>2}            *"
_TABLE = """          Energy change      -0.0000038757            0.0000050000      {c}
          RMS gradient        {g:.10f}            0.0001000000      {c}
          MAX gradient        0.0025000000            0.0003000000      {c}
          RMS step            0.0020825435            0.0020000000      {c}
          MAX step            0.0192646101            0.0040000000      YES"""
_DONE = "                             ****ORCA TERMINATED NORMALLY****"
_CONVERGED = "                  ***        THE OPTIMIZATION HAS CONVERGED     ***"
_MAXITER_TEXT = ("       The optimization did not converge but reached the maximum \n"
                 "       number of optimization cycles.")


def _coords(x):
    return ("CARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------\n"
            f"  O      {x:.6f}    0.000000    0.000000\n"
            f"  H      {x + 0.96:.6f}    0.000000    0.000000\n"
            f"  H      {x - 0.24:.6f}    0.930000    0.000000\n")


def _freqs(values):
    rows = "\n".join(f"   {i:>3}:   {f:>9.2f} cm**-1" + ("  ***imaginary mode***" if f < 0 else "")
                     for i, f in enumerate(values))
    return ("VIBRATIONAL FREQUENCIES\n-----------------------\n\n"
            "Scaling factor for frequencies =  1.000000000  (already applied!)\n\n"
            f"{rows}\n\n------------\n")


def _normal_modes(vectors):
    """A `NORMAL MODES` block: {mode index: its printed coordinate-major
    values}. ORCA prints six mode columns at a time and repeats every
    coordinate row under each block's own integer header line."""
    indices = sorted(vectors)
    n_coords = len(vectors[indices[0]])
    lines = ["NORMAL MODES", "------------", ""]
    for start in range(0, len(indices), 6):
        columns = indices[start:start + 6]
        lines.append("".join(f"{c:>11}" for c in columns))
        for row in range(n_coords):
            lines.append(f"{row:>7}" + "".join(f"{vectors[c][row]:12.6f}" for c in columns))
    return "\n".join(lines) + "\n\n"


def _opt(cycles, energies, coords_at=None, converged=True, hessian=None):
    out = ["Max. no of cycles        MaxIter  .... 3"]
    for n in range(1, cycles + 1):
        out.append(_CYCLE.format(n=n))
        if coords_at is None or n in coords_at:
            out.append(_coords(n / 10))
        out.append(f"FINAL SINGLE POINT ENERGY     {energies[n - 1]:.9f}")
        if hessian is not None:
            out.append(f"        Hessian has     {hessian} negative eigenvalue")
        last = n == cycles and converged
        out.append(_TABLE.format(g=1e-2 / n, c="YES" if last else "NO"))
    return "\n".join(out) + "\n"


def _scan(steps):
    """A finished one-parameter relaxed scan: (coordinate, energy) per step,
    two cycles each (coordinates printed on both)."""
    out = ["There is 1 parameter to be scanned.",
           f"There will be   {len(steps)} constrained geometry optimizations."]
    for k, (value, energy) in enumerate(steps, start=1):
        out += ["         *************************************************************",
                f"         *               RELAXED SURFACE SCAN STEP   {k}               *",
                f"         *                 Bond (0, 1)  :   {value:.8f}           *",
                "         *************************************************************"]
        for c in (1, 2):
            out += [_CYCLE.format(n=c), _coords(value), f"FINAL SINGLE POINT ENERGY     {energy + (2 - c) * 1e-3:.9f}"]
        out.append(_CONVERGED)
    return "\n".join(out) + "\n" + _DONE + "\n"


def _water_job(root, name="water"):
    """The water fixture: two cycles of one three-atom structure, the second
    shifted 0.1 A in x, so `--cycle 1` is a different geometry from the
    latest. Distance 0-1 is 0.960 A, angle 1-0-2 is 104.4702941 degrees."""
    d = root / name
    d.mkdir(parents=True)
    (d / "job.inp").write_text(_OPT_FREQ)
    (d / "job.out").write_text(_opt(2, [-80.00, -80.10]) + _CONVERGED + "\n" + _DONE + "\n")
    return d


_OPT_FREQ = "! B97-3c Opt Freq PAL4\n%maxcore 2000\n* xyz 0 1\nO 0 0 0\n*\n"
_OPTTS_FREQ = "! B97-3c OptTS Freq\n* xyzfile -1 2 start.xyz\n"
_SCAN_INPUT = "! B97-3c Opt\n%geom\n  Scan\n    B 0 1 = 1.0, 2.0, 3\n  end\nend\n* xyz 0 1\n*\n"
_ERROR_OUT = "\n".join([
    "some setup line",
    "     SCF ITERATIONS",
    "  line before the error 1",
    "  line before the error 2",
    "ORCA finished by error termination in SCF",
    "Calling Command: mpirun ...",
    "[file orca_tools/qcmsg.cpp, line 394]:",
]) + "\n"


def _tree(base):
    """Write the synthetic tree under `base`; returns {name: directory}."""
    from test_monitor import SYNTHETIC_OUTPUT

    jobs = {
        "opt_done": (_OPT_FREQ, _opt(3, [-76.10, -76.20, -76.25]) + _CONVERGED + "\n"
                     + _freqs([0.0, 0.0, 1600.0, 3700.0, 3800.0]) + _DONE + "\n"),
        "opt_maxiter": (_OPT_FREQ, _opt(3, [-76.10, -76.15, -76.18], coords_at={1}, converged=False)
                        + _MAXITER_TEXT + "\n" + _DONE + "\n"),
        "ts/final": (_OPTTS_FREQ, _opt(2, [-80.00, -79.99], hessian=1) + _CONVERGED + "\n"
                     + _freqs([0.0, -512.3, 45.0, 120.5]) + _DONE + "\n"),
        "scan_done": (_SCAN_INPUT.replace("! B97-3c Opt", "! QM/XTB B97-3c Opt"), SYNTHETIC_OUTPUT),
        "scan_running": (_SCAN_INPUT, SYNTHETIC_OUTPUT.split("VIBRATIONAL FREQUENCIES")[0]),
        "scan_profile": (_SCAN_INPUT, _scan([(1.0, -100.000), (1.5, -99.990), (2.0, -99.995)])),
        "broken": (_OPT_FREQ, _ERROR_OUT),
        "not_yet": (_OPT_FREQ, None),
    }
    dirs = {}
    for name, (inp, out) in jobs.items():
        d = base / name
        d.mkdir(parents=True)
        (d / "job.inp").write_text(inp)
        if out is not None:
            (d / "job.out").write_text(out)
        dirs[name] = d
    spaced = base / "B97-3c_XTB " / "pair"
    spaced.mkdir(parents=True)
    (spaced / "opt.inp").write_text(_OPT_FREQ)
    (spaced / "opt.out").write_text(_opt(2, [-10.0, -10.5]) + _CONVERGED + "\n" + _DONE + "\n")
    (spaced / "freq.inp").write_text("! B97-3c Freq\n* xyz 0 1\n*\n")
    # ORCA's own generated inputs are not jobs.
    (dirs["opt_done"] / "job.scfgrad.inp").write_text("! engrad\n")
    (dirs["opt_done"] / "job_D00001.scfgrad.inp").write_text("! engrad\n")
    dirs["spaced"] = spaced
    return dirs


class _AliveProbe:
    """Test liveness: every directory in ALIVE has a running ORCA."""
    name = "test"
    ALIVE: set = set()

    def snapshot(self, job_dirs=None):
        from orcamon.core.liveness import Liveness
        import time as _time
        return {d.resolve(): Liveness(True, self.name, since=_time.time() - 60) for d in self.ALIVE}

    def default(self, job_dir):
        from orcamon.core.liveness import Liveness
        return Liveness(False, self.name)


def _orcamon(argv):
    """`orcamon ARGV` in-process: (exit code, stdout, stderr)."""
    import contextlib
    import io
    from orcamon import cli
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = cli.main(list(argv))
    return code, out.getvalue(), err.getvalue()


def _json(argv):
    import json
    code, out, err = _orcamon([*argv, "--json"])
    try:
        return code, json.loads(out), err
    except ValueError:
        return code, None, out + err


class _Tree:
    """The synthetic tree, the test liveness registered, the cache isolated."""

    def __enter__(self):
        import tempfile
        from pathlib import Path
        from orcamon.core import liveness
        self._tmp = tempfile.TemporaryDirectory()
        base = Path(self._tmp.name)
        self._env = {k: os.environ.get(k) for k in ("XDG_CACHE_HOME", "ORCAMON_ROOT", "COLUMNS")}
        os.environ["XDG_CACHE_HOME"] = str(base / "cache")
        os.environ.pop("ORCAMON_ROOT", None)
        os.environ["COLUMNS"] = "120"
        self.root = base / "tree"
        self.dirs = _tree(self.root)
        liveness.PROBES["test"] = _AliveProbe
        _AliveProbe.ALIVE = {self.dirs["scan_running"]}
        self.args = ["--root", str(self.root), "--liveness", "test"]
        return self

    def run(self, *argv):
        return _orcamon([*self.args, *argv])

    def json(self, *argv):
        return _json([*self.args, *argv])

    def __exit__(self, *exc):
        from orcamon.core import liveness
        liveness.PROBES.pop("test", None)
        for k, v in self._env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        self._tmp.cleanup()


def test_a_job_argument_names_one_job():
    print("\na JOB argument resolves to one job, or is refused with the candidates")
    with _Tree() as t:
        code, out, err = t.run("show", "opt_done")
        check("a unique substring of a label", code == 0 and "opt_done" in out, err)
        code, out, err = t.run("show", "opt")
        check("an ambiguous one exits 2 and lists the candidates",
              code == 2 and "opt_done" in err and "opt_maxiter" in err, err)
        code, out, err = t.run("show", "no_such_job")
        check("no match exits 2", code == 2 and "no job matches" in err, err)
        code, out, err = t.run("show", "ts/final")
        check("an exact label wins over the labels it is a substring of", code == 0, err)
        code, out, err = t.run("show", "ts")
        check("a directory holding no job of its own falls through to the label it is part of",
              code == 0 and out.startswith("ts/final"), out + err)
        code, out, err = t.run("show", str(t.dirs["spaced"]))
        check("a directory with two jobs exits 2 naming both files",
              code == 2 and "opt.inp" in err and "freq.inp" in err, err)
        code, out, err = t.run("show", str(t.dirs["spaced"] / "opt.out"))
        check("a path with a space in it resolves, by its file", code == 0 and "pair/opt" in out, out + err)
        code, doc, err = t.json("ls")
        labels = [j["label"] for j in doc["jobs"]]
        check("two jobs in one directory are two jobs, labelled dir/stem",
              "B97-3c_XTB /pair/opt" in labels and "B97-3c_XTB /pair/freq" in labels, f"{labels}")
        check("and ORCA's generated inputs are not jobs",
              not any("scfgrad" in label for label in labels), f"{labels}")


def test_only_orca_inputs_are_jobs():
    print("\nonly a .inp that reads as ORCA's is a job")
    import tempfile
    from pathlib import Path
    from orcamon.core.discovery import discover, stems_in

    # A `!`/`%`/`*` line is ORCA's; CREST and xtb use `$constrain ... $end`.
    inputs = {
        "orca/job.inp": "! B97-3c SP\n* xyz 0 1\nH 0 0 0\n*\n",
        "indented/job.inp": "   ! B97-3c SP\n",
        "blockonly/job.inp": "%pal nprocs 4 end\n",
        "commented/job.inp": "# a comment\n$metadyn\n  atoms: 1-3\n$end\n",
        "empty/job.inp": "",
        "crest/fixhost.inp": "$constrain\n  atoms: 1-130\n  force constant=0.05\n$end\n",
    }
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for rel, text in inputs.items():
            path = root / rel
            path.parent.mkdir(parents=True)
            path.write_text(text)
        labels = [r.label for r in discover(root)]
        check("discover keeps only the three inputs that read as ORCA's",
              labels == ["blockonly", "indented", "orca"], f"{labels}")
        stems = stems_in(root / "crest")
        check("a CREST constraint file is not a stem", stems == [], f"{stems}")
        code, doc, err = _json(["--root", str(root), "--liveness", "mtime", "--no-cache", "ls"])
        listed = doc["jobs"] if doc else None
        check("ls lists exactly those three jobs",
              listed is not None and [j["label"] for j in listed] == ["blockonly", "indented", "orca"],
              f"{code} {listed if listed is not None else err}")


def test_ls_lists_every_job_boundedly():
    print("\nls: one line per job, bounded, with a JSON twin")
    with _Tree() as t:
        code, out, err = t.run("ls")
        lines = out.strip().split("\n")
        check("exit 1, because one listed job failed", code == 1, f"{code} {err}")
        check("a header and one line per job", len(lines) == 1 + 10, f"{len(lines)}\n{out}")
        check("the failed job carries its flag", any("broken" in ln and "! failed" in ln for ln in lines), out)
        check("the MaxIter job reads finished with its flag, not converged",
              any("opt_maxiter" in ln and ln.startswith("finished") and "opt_not_converged" in ln
                  for ln in lines), out)
        code, out, err = t.run("ls", "--max-lines", "4")
        lines = out.strip().split("\n")
        check("--max-lines holds, and the cut says so", len(lines) == 4 and "more lines" in lines[-1], out)
        code, doc, err = t.json("ls")
        check("--json parses and carries the schema", doc is not None and doc.get("schema") == 1, err)
        check("without crash lines or criteria per job",
              doc is not None and all("crash_lines" not in j and "criteria" not in j for j in doc["jobs"]))
        code, doc, err = t.json("ls", "--attention")
        check("--attention keeps only flagged jobs",
              doc is not None and doc["jobs"] and all(j["attention"] for j in doc["jobs"]))
        code, doc, err = t.json("ls", "--status", "running")
        check("--status filters", doc is not None and [j["label"] for j in doc["jobs"]] == ["scan_running"],
              f"{doc and [j['label'] for j in doc['jobs']]}")
        code, out, err = t.run("ls", "--status", "bogus")
        check("an unknown status is a usage error", code == 2, err)


def test_ls_prints_labels_an_agent_can_pass_back():
    print("\nls prints labels whole: the label is the handle the other commands take")
    import tempfile
    from pathlib import Path

    label = ("reaction_with_a_long_descriptive_name/geometry_b973c-xtb_with_implicit_solvent"
             "/ts/optts_freq_tight_attempt2")
    saved = os.environ.get("COLUMNS")
    try:
        os.environ["COLUMNS"] = "120"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            job_dir = root / label
            job_dir.mkdir(parents=True)
            (job_dir / "job.inp").write_text("! B97-3c Opt\n* xyz 0 1\nH 0 0 0\n*\n")
            argv = ["--root", str(root), "--liveness", "mtime", "--no-cache"]
            code, out, err = _orcamon([*argv, "ls"])
            check("the whole 107-character label is printed", label in out,
                  f"{code} {out[:300]} {err[:200]}")
            check("ls shortens nothing", "…" not in out, out[:300])
            code, out, err = _orcamon([*argv, "show", label])
            check("show takes the label ls printed", code == 0, f"{code} {out[:200]} {err[:200]}")
    finally:
        if saved is None:
            os.environ.pop("COLUMNS", None)
        else:
            os.environ["COLUMNS"] = saved


def test_ls_says_where_an_irc_is():
    print("\nls says where an IRC is and that a direction hit MaxIter")
    import tempfile
    from pathlib import Path
    # One IRC fixture, not two: the duplicate guard reads functions even when
    # they are private, and its rule is to import the other copy, not clone it.
    from test_monitor import _irc_tree

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _irc_tree(root / "irc")
        code, out, err = _orcamon(["--root", str(root), "--liveness", "mtime", "--no-cache", "ls"])
        row = next((line for line in out.split("\n") if "irc" in line), "")
        check("the IRC row carries its progress and its flag",
              "irc B2/F3" in row and "irc_not_converged" in row, f"{code} {out} {err}")


def test_show_and_ls_agree():
    print("\nshow and ls report the same numbers for the same job")
    with _Tree() as t:
        _, listed, _ = t.json("ls")
        by_label = {j["label"]: j for j in listed["jobs"]}
        for name in ("opt_done", "scan_done", "ts/final"):
            code, shown, err = t.json("show", name)
            job = shown["job"] if shown else {}
            check(f"{name}: the same energy and status",
                  job.get("energy_eh") == by_label[name]["energy_eh"]
                  and job.get("status") == by_label[name]["status"], err)
        code, out, err = t.run("show", "broken")
        check("show exits 1 on a failed job and prints its crash line",
              code == 1 and "error termination in SCF" in out, out)
        code, out, err = t.run("show", "scan_running")
        check("and 0 on a running one", code == 0 and out.split("\n")[5].startswith("running"), out)


def test_geom_never_substitutes_a_geometry():
    print("\ngeom: XYZ of a real point, or exit 2 naming the nearest earlier one")
    with _Tree() as t:
        code, out, err = t.run("geom", "opt_done")
        lines = out.split("\n")
        check("the latest geometry as XYZ", code == 0 and lines[0] == "3" and "cycle 3" in lines[1]
              and lines[2].startswith("O "), out)
        code, out, err = t.run("geom", "opt_maxiter", "--cycle", "3")
        check("a cycle whose coordinates were never printed exits 2",
              code == 2 and "cycle 3: coordinates not printed; nearest earlier with coordinates: cycle 1" in err,
              err)
        code, out, err = t.run("geom", "opt_done", "--cycle", "9")
        check("a cycle never reached exits 2", code == 2 and "not in the kept history" in err, err)

        # SOLVATOR: no coordinate block in the log, the cluster only in the
        # file the log names.
        solv = t.root / "solvated"
        solv.mkdir()
        (solv / "job.inp").write_text("! XTB ALPB(Water) SOLVATOR\n* xyzfile 0 1 start.xyz\n")
        (solv / "job.solvator.xyz").write_text("2\ncluster\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\n")
        (solv / "job.out").write_text("Final structured saved to        :             job.solvator.xyz\n"
                                      + _DONE + "\n")
        code, out, err = t.run("geom", "solvated")
        lines = out.split("\n")
        check("a job that printed no coordinates gives the file it wrote, and names it",
              code == 0 and lines[0] == "2" and lines[1] == "solvated · job.solvator.xyz", out + err)
        code, doc, err = t.json("geom", "solvated")
        check("--json says where it came from", doc is not None and doc["point"]["source"] == "job.solvator.xyz"
              and doc["point"]["cycle"] == 0 and len(doc["atoms"]) == 2, f"{doc}")
        code, out, err = t.run("geom", "solvated", "--cycle", "1")
        check("but a cycle that was never printed is still refused", code == 2, out + err)


def test_snapshot_writes_a_png():
    print("\nsnapshot: renders a job's geometry to a PNG, or exits 2 with no geometry")
    from PIL import Image

    with _Tree() as t:
        target = t.root.parent / "x.png"
        code, out, err = t.run("snapshot", "opt_done", "-o", str(target), "--size", "120x100")
        check("exit 0, and the file is a PNG at the size asked for",
              code == 0 and target.read_bytes()[:4] == b"\x89PNG"
              and Image.open(target).size == (120, 100), f"{code} {out} {err}")
        check("and the one line names the file, size, atom count and source",
              out.strip() == f"{target}: 120x100, 3 atoms, cycle 3", out)

        nothing = t.root / "nothing"
        nothing.mkdir()
        (nothing / "job.inp").write_text("! XTB SP\n* xyzfile 0 1 missing.xyz\n")
        empty = t.root.parent / "empty.png"
        code, out, err = t.run("snapshot", "nothing", "-o", str(empty))
        check("a job with no geometry anywhere exits 2 and says so",
              code == 2 and "no geometry to draw" in err and not empty.exists(), f"{code} {out} {err}")

        code, out, err = t.run("snapshot", "opt_done", "-o", str(target), "--size", "12")
        check("a bad --size is a usage error", code == 2 and "--size" in err, err)
        for too_big in ("100000x100000", "0x100"):
            code, out, err = t.run("snapshot", "opt_done", "-o", str(target), "--size", too_big)
            check(f"--size {too_big} is a usage error, not a MemoryError",
                  code == 2 and "--size" in err and "Traceback" not in err, err)

        space = t.root.parent / "space.png"
        code, out, err = t.run("snapshot", "opt_done", "-o", str(space), "--size", "120x100",
                               "--representation", "space-filling")
        check("--representation space-filling renders too",
              code == 0 and space.read_bytes()[:4] == b"\x89PNG", f"{code} {out} {err}")

        see = t.root.parent / "see.png"
        code, out, err = t.run("snapshot", "opt_done", "-o", str(see), "--size", "120x100",
                               "--see-through")
        check("--see-through renders too", code == 0 and see.read_bytes()[:4] == b"\x89PNG",
              f"{code} {out} {err}")

        # A frequency job with one imaginary mode (ORCA mode 0) whose 9
        # printed values cover the three drawn atoms, so --mode can draw a
        # phase displaced along that mode's own pattern.
        modes_job = t.root / "mode_ts"
        modes_job.mkdir()
        (modes_job / "job.inp").write_text(_OPT_FREQ)
        vectors = {0: [0.1, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.3]}
        vectors.update({i: [0.0] * 9 for i in range(1, 9)})
        (modes_job / "job.out").write_text(
            _opt(2, [-80.00, -80.10]) + _CONVERGED + "\n"
            + _freqs([-100.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0])
            + _normal_modes(vectors) + _DONE + "\n")
        mode_png = t.root.parent / "mode.png"
        code, out, err = t.run("snapshot", "mode_ts", "-o", str(mode_png),
                               "--mode", "0", "--size", "120x100")
        check("--mode renders a phase of the imaginary mode and says which",
              code == 0 and mode_png.read_bytes()[:4] == b"\x89PNG"
              and out.strip().endswith("3 atoms, cycle 2, mode 0 at 90 deg"), f"{code} {out} {err}")
        code, out, err = t.run("snapshot", "mode_ts", "-o", str(mode_png), "--mode", "999")
        check("a mode the job does not have exits 2 naming the ones it does",
              code == 2 and "available: [0]" in err, err)
        code, out, err = t.run("snapshot", "opt_done", "-o", str(mode_png), "--mode", "0")
        check("a job with a frequency block but no imaginary mode exits 2",
              code == 2 and "available: none" in err, err)
        code, out, err = t.run("snapshot", "opt_maxiter", "-o", str(mode_png), "--mode", "0")
        check("a job with no frequency block exits 2 naming the cause",
              code == 2 and "no frequency block" in err, err)


def test_a_path_is_listed_and_drawn():
    print("\nan IRC path is listed by energies and drawn by geom/snapshot --point")
    import tempfile
    from pathlib import Path
    # One IRC fixture, not two: the duplicate guard reads functions even when
    # they are private, and its rule is to import the other copy, not clone it.
    from test_monitor import _irc_tree

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _irc_tree(root / "irc")
        (root / "plain").mkdir()
        (root / "plain" / "job.inp").write_text(
            "! B97-3c SP\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")
        argv = ["--root", str(root), "--liveness", "mtime", "--no-cache"]

        code, doc, err = _json([*argv, "energies", "irc"])
        rows = doc["rows"] if doc else []
        first = rows[0] if rows else {}
        check("energies --json lists the path's points with their dE and monitors",
              code == 0 and doc is not None and doc["mode"] == "irc" and len(rows) == 6
              and {**first, "dE_kj_mol": round(first.get("dE_kj_mol", 0.0), 3)} == {
                  "point": -2, "label": "IRC backward 1", "energy_eh": -100.00239,
                  "dE_kj_mol": -6.276,
                  "monitors": {"B(O 0,H 1)": 0.94, "B(O 0,H 2)": 0.98}},
              f"{code} {doc} {err}")

        code, out, err = _orcamon([*argv, "energies", "irc"])
        check("energies prints the path with ORCA's own dE, and says where dE is from",
              code == 0 and "IRC TS" in out and "dE from the TS" in out, out + err)

        code, out, err = _orcamon([*argv, "geom", "irc", "--point", "1"])
        lines = out.split("\n")
        check("geom --point gives that point's geometry",
              code == 0 and len(lines) > 2 and "IRC forward 0" in lines[1]
              and lines[2].split()[-1] == "3.00000000", f"{code} {out} {err}")

        code, out, err = _orcamon([*argv, "geom", "irc", "--point", "9"])
        check("a point off the path exits 2 with the range it does have",
              code == 2 and "not on the path (points -2 to 3)" in err, err)

        code, out, err = _orcamon([*argv, "geom", "irc", "--point", "0", "--cycle", "1"])
        check("--point refuses to be combined with a cycle or a step",
              code == 2, f"{code} {out} {err}")

        code, out, err = _orcamon([*argv, "geom", "plain", "--point", "0"])
        check("a job with no path refuses --point",
              code == 2 and "no reaction path" in err, err)

        png = root / "p.png"
        code, out, err = _orcamon([*argv, "snapshot", "irc", "--point", "-2", "-o", str(png)])
        check("snapshot --point draws that point and names it",
              code == 0 and png.exists() and "IRC backward 1" in out, f"{code} {out} {err}")

        code, out, err = _orcamon([*argv, "snapshot", "irc", "--point", "0", "--mode", "6",
                                   "-o", str(root / "q.png")])
        check("snapshot --point refuses --mode", code == 2, f"{code} {out} {err}")


def test_an_neb_is_listed_and_drawn():
    print("\nan NEB path is listed by energies and drawn by geom, show and the TUI")
    import asyncio
    import tempfile
    import time as _time
    from pathlib import Path
    # One NEB fixture, not two: the duplicate guard reads functions even when
    # they are private, and its rule is to import the other copy, not clone it.
    from test_monitor import _neb_tree

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _neb_tree(root / "neb")
        argv = ["--root", str(root), "--liveness", "mtime", "--no-cache"]

        code, doc, err = _json([*argv, "energies", "neb"])
        rows = doc["rows"] if doc else []
        check("energies --json lists one row per image, in path order, "
              "dE from image 0 and CI on the latest row's image",
              doc is not None and doc["mode"] == "neb" and len(rows) == 4
              and [round(r["dE_kj_mol"], 4) for r in rows] == [0.0, 13.1275, 21.004, 5.251]
              and rows[2]["label"] == "NEB image 2 (CI)", f"{doc} {err}")

        code, out, err = _orcamon([*argv, "geom", "neb", "--point", "2"])
        lines = out.split("\n")
        check("geom --point gives the image's own geometry",
              code == 0 and len(lines) > 2 and lines[2].split()[-1] == "2.00000000",
              f"{code} {out} {err}")

        code, out, err = _orcamon([*argv, "show", "neb"])
        check("`show` prints where the NEB is and what barrier it has",
              "NEB        iteration 2 · climbing image 2 · E(CI)-E(0) +21.0 kJ/mol" in out,
              f"{code} {out} {err}")

        code, out, err = _orcamon([*argv, "ls"])
        row = next((line for line in out.split("\n") if "neb" in line), "")
        check("the ls row carries the NEB's progress", "neb 2 CI2" in row, f"{code} {out} {err}")

        async def drive():
            from orcamon.tui.app import MonitorApp
            app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
            async with app.run_test(size=(180, 50)) as pilot:
                chart = app.query_one("#chart")
                deadline = _time.monotonic() + 5
                while not (app.jobs and all(j.parsed for j in app.jobs)
                           and chart._points) and _time.monotonic() < deadline:
                    await pilot.pause(0.05)
                title = app.query_one("#geometry").border_title
                check("the chart is the NEB path and the pane follows the climbing image",
                      chart._mode == "neb" and title.endswith("· NEB image 2 (CI)"),
                      f"{chart._mode} {title!r}")

        asyncio.run(drive())


def test_any_mode_is_drawn_and_described():
    print("\nsnapshot --mode draws any mode; freqs --mode names the atoms it moves")
    import tempfile
    from pathlib import Path
    from test_monitor import _MODE_FREQS, _MODE_VECTORS

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        job_dir = root / "modes_all"
        job_dir.mkdir()
        (job_dir / "job.inp").write_text(_OPT_FREQ)
        (job_dir / "job.out").write_text(
            _opt(2, [-80.00, -80.10]) + _CONVERGED + "\n"
            + _freqs(_MODE_FREQS) + _normal_modes(_MODE_VECTORS) + _DONE + "\n")
        argv = ["--root", str(root), "--liveness", "mtime", "--no-cache"]

        png = root / "p.png"
        code, out, err = _orcamon([*argv, "snapshot", "modes_all", "--mode", "5",
                                   "-o", str(png), "--size", "120x100"])
        check("a real mode is drawn like an imaginary one, and named",
              code == 0 and png.read_bytes()[:4] == b"\x89PNG"
              and out.strip().endswith("3 atoms, cycle 2, mode 5 at 90 deg"),
              f"{code} {out} {err}")
        code, out, err = _orcamon([*argv, "snapshot", "modes_all", "--mode", "0",
                                   "-o", str(root / "q.png")])
        check("a zero-frequency mode is refused, with what the job does have",
              code == 2 and "available: [3] imaginary; real: 5 modes, 4-8" in err, err)

        code, out, err = _orcamon([*argv, "freqs", "modes_all", "--mode", "3"])
        lines = out.strip().split("\n")
        check("freqs --mode heads with the mode, then its atoms largest first",
              code == 0 and lines[0] == "mode 3  -150.00 cm-1 (imaginary), drawn on 3 atoms"
              and lines[1] == "  atom  el  relative"
              and lines[2:5] == ["     2  H       1.00", "     1  H       0.67", "     0  O       0.33"],
              f"{code} {out} {err}")
        code, doc, err = _json([*argv, "freqs", "modes_all", "--mode", "3"])
        check("--json gives the same ranking to four decimals",
              doc is not None and doc["mode"] == {
                  "index": 3, "cm1": -150.0,
                  "atoms": [{"atom": 2, "element": "H", "relative": 1.0},
                            {"atom": 1, "element": "H", "relative": 0.6667},
                            {"atom": 0, "element": "O", "relative": 0.3333}]},
              f"{code} {doc} {err}")
        code, out, err = _orcamon([*argv, "freqs", "modes_all", "--mode", "5", "--top", "1"])
        lines = out.strip().split("\n")
        check("--top bounds the ranking",
              code == 0 and lines[2] == "     0  O       1.00" and lines[3].startswith("relative:"),
              f"{code} {out} {err}")
        code, out, err = _orcamon([*argv, "freqs", "modes_all", "--mode", "0"])
        check("a mode the job does not have is refused naming the ones it does",
              code == 2 and "available: [3] imaginary" in err, err)
        code, out, err = _orcamon([*argv, "freqs", "modes_all", "--top", "3"])
        check("--top without --mode is a usage error", code == 2 and "--top needs --mode" in err, err)
        code, out, err = _orcamon([*argv, "freqs", "modes_all"])
        check("plain freqs is unchanged",
              code == 0 and "9 modes (3 zero, 1 imaginary)" in out and "lowest 5 real:" in out, out)


def test_measure_names_the_geometry_it_measured():
    print("\nmeasure: the distance, angle or dihedral on the geometry a flag names")
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _water_job(root)
        argv = ["--root", str(root), "--liveness", "mtime", "--no-cache"]

        code, out, err = _orcamon([*argv, "measure", "water", "0", "1"])
        check("a distance prints the atoms, their elements and the cycle",
              code == 0 and out == "distance 0-1: 0.960 Å (O-H) · water · cycle 2\n",
              f"{code} {(out + err)!r}")

        code, doc, err = _json([*argv, "measure", "water", "1", "0", "2"])
        check("--json gives the kind, the value to four decimals and the atoms",
              doc is not None and doc["kind"] == "angle" and doc["value"] == 104.4703
              and doc["unit"] == "degree"
              and doc["atoms"] == [{"atom": 1, "element": "H"}, {"atom": 0, "element": "O"},
                                   {"atom": 2, "element": "H"}]
              and doc["point"]["cycle"] == 2 and doc["job"] == "water",
              f"{code} {doc} {err}")

        code, out, err = _orcamon([*argv, "measure", "water", "0", "1", "--cycle", "1"])
        check("--cycle measures the geometry that cycle printed",
              code == 0 and out == "distance 0-1: 0.960 Å (O-H) · water · cycle 1\n",
              f"{code} {(out + err)!r}")

        code, out, err = _orcamon([*argv, "measure", "water", "0", "5"])
        check("an atom the geometry does not have is refused with its size",
              code == 2 and "atom 5 is out of range: the geometry has 3 atoms (0-2)" in err, err)

        code, out, err = _orcamon([*argv, "measure", "water", "0"])
        check("fewer than two atoms is a usage error",
              code == 2 and "measure takes 2 to 4 atoms, got 1" in err, err)

        code, out, err = _orcamon([*argv, "measure", "water", "0", "1", "--cycle", "9"])
        check("a cycle the kept history does not hold is refused",
              code == 2 and "cycle 9: not in the kept history" in err, err)


def test_the_bond_cache_holds_the_geometry_it_answers_for():
    print("\nthe geometry pane's bond cache is keyed by the atom list itself, not its id")
    from types import SimpleNamespace

    from orcamon.tui.app import RotatableGeometryImage

    # An id is an address CPython reuses as soon as a list is freed, so a
    # cache keyed by `id(atoms)` could serve a superseded geometry's bonds to
    # a new one. Holding the list keeps its address from being reused.
    pane = SimpleNamespace(_bond_cache=None)
    water = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0)]
    first = RotatableGeometryImage._bonds_for(pane, water)
    check("the cache holds the list it was computed for",
          pane._bond_cache[0] is water and sorted(first) == [(0, 1), (0, 2)], f"{pane._bond_cache}")
    check("the same list is answered from the cache",
          RotatableGeometryImage._bonds_for(pane, water) is first)
    pair = [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0)]
    check("a different list is bonded afresh",
          RotatableGeometryImage._bonds_for(pane, pair) == [(0, 1)] and pane._bond_cache[0] is pair,
          f"{pane._bond_cache}")


def test_wait_returns_when_the_job_is_done():
    print("\nwait: returns on the condition, times out with 4")
    import threading
    with _Tree() as t:
        out_file = t.dirs["scan_running"] / "job.out"

        def finish():
            import time as _time
            _time.sleep(1.0)
            with open(out_file, "a") as f:
                f.write(_DONE + "\n")
            _AliveProbe.ALIVE = set()

        writer = threading.Thread(target=finish)
        import time as _time
        started = _time.monotonic()
        writer.start()
        code, out, err = t.run("wait", "scan_running", "--until", "done", "--interval", "0.2", "--timeout", "10")
        elapsed = _time.monotonic() - started
        writer.join()
        check("exit 0 once the termination line lands", code == 0 and "met --until done" in out,
              f"{code} {out} {err}")
        check("within 5 s", elapsed < 5, f"{elapsed:.1f} s")

        _AliveProbe.ALIVE = {t.dirs["scan_running"]}
        code, out, err = t.run("wait", "scan_done", "scan_running", "--until", "done")
        check("with several jobs, any one done returns at once", code == 0 and "scan_done" in out, out + err)
        code, out, err = t.run("wait", "opt_done", "--until", "change", "--interval", "0.2", "--timeout", "1")
        check("nothing changing times out with 4", code == 4 and "timed out" in out, f"{code} {out}")
        code, out, err = t.run("wait", "opt_done", "--until", "sometime")
        check("a bad condition is a usage error", code == 2, err)


def test_the_other_commands_answer_boundedly():
    print("\nconv, energies, freqs, input, errors, tail")
    with _Tree() as t:
        code, out, err = t.run("conv", "opt_done")
        lines = out.strip().split("\n")
        check("conv: header, tolerances, one row per cycle, the legend",
              code == 0 and len(lines) == 2 + 3 + 1 and lines[-2].lstrip().startswith("3 "), out)
        check("with the met criteria starred in the text", lines[-2].count("*") == 5, lines[-2])
        code, doc, err = t.json("conv", "opt_done", "--last", "2")
        check("--last and --json: the latest cycles, each criterion with value, tol and ok",
              doc is not None and [r["cycle"] for r in doc["cycles"]] == [2, 3]
              and doc["cycles"][-1]["criteria"]["RMS grad"]["ok"] is True, f"{doc}")
        code, out, err = t.run("conv", "opt_done", "--max-lines", "3")
        check("and --max-lines holds there too", len(out.strip().split("\n")) == 3, out)

        code, doc, err = t.json("energies", "opt_done")
        check("energies of an optimization: one row per cycle, dE from the first",
              doc is not None and doc["mode"] == "opt" and len(doc["rows"]) == 3
              and abs(doc["rows"][-1]["dE_kj_mol"] - (-0.15 * 2625.4996394799)) < 1e-6, f"{doc}")
        code, out, err = t.run("energies", "scan_profile")
        lines = out.strip().split("\n")
        check("energies of a scan: a row per step with its coordinate, the maximum marked",
              code == 0 and len(lines) == 4 and "<- max" in lines[2] and "1.5" in lines[2], out)
        code, doc, err = t.json("energies", "scan_profile")
        top = max(doc["rows"], key=lambda r: r["energy_eh"]) if doc else {}
        check("and --json agrees with show's scan maximum",
              doc is not None and top.get("max") is True
              and abs(top["dE_first_kj_mol"] - t.json("show", "scan_profile")[1]["job"]["scan_max_kj_mol"]) < 1e-9,
              f"{doc}")

        code, out, err = t.run("freqs", "ts/final")
        check("freqs: the count, the imaginary mode with its index, then the lowest real",
              code == 0 and "4 modes (1 zero, 1 imaginary)" in out and "mode    1     -512.30" in out
              and "lowest 2 real" in out, out)
        code, doc, err = t.json("freqs", "ts/final", "--lowest", "1")
        check("--lowest bounds the real modes", doc is not None and len(doc["lowest_real"]) == 1
              and doc["lowest_real"][0]["cm1"] == 45.0, f"{doc}")
        code, out, err = t.run("freqs", "not_yet")
        check("no block yet says so", code == 0 and out.strip() == "no frequency block yet", out)

        code, out, err = t.run("input", "ts/final")
        check("input: run types, spin and the coordinate file",
              code == 0 and "OptTS Freq" in out and "charge -1 · doublet" in out and "xyzfile start.xyz" in out, out)
        code, doc, err = t.json("input", "not_yet")
        check("and it works on a job with no output",
              code == 0 and doc is not None and doc["input"]["nprocs"] == 4 and doc["input"]["maxcore_mb"] == 2000,
              f"{doc} {err}")

        code, out, err = t.run("errors", "broken", "--context", "2")
        check("errors: the error line with its context, exit 1 for a failed job",
              code == 1 and "> ORCA finished by error termination in SCF" in out
              and "  line before the error 2" in out and "  Calling Command: mpirun ..." in out
              and "line before the error 1" in out and "some setup line" not in out, out)
        code, doc, err = t.json("errors", "broken", "--context", "0")
        check("--json counts them", doc is not None and doc["crash_count"] == 1 and doc["qm2_error_count"] == 0
              and doc["matches"][0]["before"] == [], f"{doc}")
        code, out, err = t.run("errors", "not_yet")
        check("no output is a usage error, not a crash", code == 2 and "no output yet" in err, err)

        code, out, err = t.run("tail", "opt_done", "-n", "2")
        check("tail -n: the last lines", out.split("\n")[:2] == ["------------", _DONE], out)
        code, out, err = t.run("tail", "opt_done", "--grep", "FINAL SINGLE")
        check("tail --grep: the matching lines only", code == 0 and len(out.strip().split("\n")) == 3
              and all("FINAL SINGLE" in line for line in out.strip().split("\n")), out)
        code, out, err = t.run("tail", "opt_done", "--grep", "(")
        check("a bad pattern is a usage error", code == 2, err)

        code, out, err = t.run()
        check("bare `orcamon` without a terminal points to the commands instead of hanging in the TUI",
              code == 2 and "orcamon ls" in err, err)


def test_a_cached_state_equals_a_fresh_parse():
    print("\nthe state cache: a resumed parse equals a fresh one, and a stale entry is refused")
    import tempfile
    from pathlib import Path
    from test_monitor import SYNTHETIC_OUTPUT, _state_fields
    from orcamon.core import cache
    from orcamon.core.parser import new_state, update_job

    saved = os.environ.get("XDG_CACHE_HOME")
    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_CACHE_HOME"] = str(Path(tmp) / "cache")
        try:
            job_dir = Path(tmp) / "job"
            job_dir.mkdir()
            out = job_dir / "job.out"
            cut = SYNTHETIC_OUTPUT.index("VIBRATIONAL FREQUENCIES")
            head, rest = SYNTHETIC_OUTPUT[:cut], SYNTHETIC_OUTPUT[cut:]
            out.write_text(head)

            first = new_state(job_dir)
            update_job(first)
            check("an entry is written for a parsed output", cache.save(first))

            with open(out, "a") as f:
                f.write(rest)
            resumed = new_state(job_dir)
            check("and read back", cache.restore(resumed) and resumed.offset == len(head.encode()))
            update_job(resumed)
            fresh = new_state(job_dir)
            update_job(fresh)
            a, b = _state_fields(resumed), _state_fields(fresh)
            differ = [k for k in a if k not in ("mtime",) and a[k] != b[k]]
            check("restored, then updated with the appended bytes, it equals a fresh parse",
                  not differ, f"{differ}")

            original_key = cache.CACHE_KEY
            try:
                cache.CACHE_KEY = "a different parser"
                check("an entry written by another parser is refused",
                      not cache.restore(new_state(job_dir)))
            finally:
                cache.CACHE_KEY = original_key

            cache.save(fresh)
            replacement = job_dir / "job.out.new"
            replacement.write_text(SYNTHETIC_OUTPUT + "\n" * 10)
            os.replace(replacement, out)
            check("a replaced file (a new inode) is refused",
                  not cache.restore(new_state(job_dir)))

            long = new_state(job_dir)
            update_job(long)
            cache.save(long)
            out.write_text("short\n")  # truncated in place: same inode
            check("so is a file shorter than the entry's offset",
                  not cache.restore(new_state(job_dir)))

            cache.entry_path(out).write_bytes(b"not a pickle")
            check("and a corrupt entry is no entry, not an error", not cache.restore(new_state(job_dir)))

            mid = new_state(job_dir)
            mid.feed_text("CARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------\n")
            mid.has_out, mid.offset = True, 10
            check("a state caught inside a block is not saved", not cache.save(mid))
        finally:
            if saved is None:
                os.environ.pop("XDG_CACHE_HOME", None)
            else:
                os.environ["XDG_CACHE_HOME"] = saved


def test_the_cli_uses_the_cache_and_can_refuse_it():
    print("\nthe CLI reads through the cache, and --no-cache bypasses it")
    from orcamon.core import cache
    with _Tree() as t:
        entry = cache.entry_path(t.dirs["opt_done"] / "job.out")
        code, _, _ = t.run("--no-cache", "show", "opt_done")
        check("--no-cache writes nothing", code == 0 and not entry.exists())
        code, first, _ = t.json("show", "opt_done")
        check("a command writes an entry", code == 0 and entry.exists())
        code, second, _ = t.json("show", "opt_done")
        check("and a second run from it reports the same", first == second or
              {k: v for k, v in first["job"].items() if k != "last_output_age_s"}
              == {k: v for k, v in second["job"].items() if k != "last_output_age_s"})


def test_events_announce_changes_not_the_first_scan():
    print("\nevents: a change is announced, the first sighting is not")
    from dataclasses import replace
    from orcamon.core.events import events
    from orcamon.core.report import JobReport

    base = JobReport(
        schema=1, label="rxn/ts_anion", path="/x", stem="job", run_types=["OptTS"], method="", basis=None,
        charge=0, mult=1, layers={}, multilayer=False, n_atoms=None, n_qm_atoms=None, nprocs=None,
        maxcore_mb=None, status="running", liveness_source="process", sched_id=None, cycle=3,
        max_cycles=None, scan_step=None, scan_total=None, scan_coordinate=None, wall_time_s=None,
        last_output_age_s=None, opt_converged=False, opt_maxiter_reached=False, criteria=None,
        negative_eigenvalues=None, imaginary_freqs=None, freq_cycle=None, energy_eh=None,
        energy_label=None, delta_kj_mol=None, scan_max_kj_mol=None, scan_max_at=None, qm2_errors=0,
        crash_lines=[], attention=[],
    )
    flag = {"code": "ts_imaginary", "message": "2 imaginary frequencies (a TS wants 1)"}
    done = replace(base, status="finished", attention=[flag])
    check("the first report of a job raises nothing", events(None, done) == [])
    check("an unchanged report raises nothing", events(base, base) == [])
    raised = events(base, done)
    check("finishing raises one event that carries its flag",
          len(raised) == 1 and raised[0].text == "ts_anion: finished · 2 imaginary frequencies (a TS wants 1)",
          f"{raised}")
    later = replace(done, attention=[flag, {"code": "qm2_errors", "message": "1 QM2 error"}])
    raised = events(done, later)
    check("a new flag on its own raises a flag event",
          [(e.kind, e.text) for e in raised] == [("flag", "ts_anion: 1 QM2 error")], f"{raised}")
    check("reaching running is not news", events(replace(base, status="queued"), base) == [])


def test_notifications_reach_the_terminal_through_tmux():
    print("\nnotifications: bell and OSC, wrapped for tmux byte-exactly")
    from orcamon.core.events import Event
    from orcamon.tui import notify

    check("tmux passthrough doubles every ESC inside the envelope",
          notify.tmux_wrap("\x1b]9;hi\x07") == "\x1bPtmux;\x1b\x1b]9;hi\x07\x1b\\")
    event = Event("status", "a/b", "finished", "b: finished")
    check("all = bell, OSC 9, OSC 777",
          notify.sequences(event, "all", in_tmux=False)
          == "\a\x1b]9;b: finished\x07\x1b]777;notify;orcamon;b: finished\x07")
    check("inside tmux the OSCs are wrapped, the bell is not",
          notify.sequences(event, "all", in_tmux=True)
          == "\a" + notify.tmux_wrap("\x1b]9;b: finished\x07")
          + notify.tmux_wrap("\x1b]777;notify;orcamon;b: finished\x07"))
    check("off is silent", notify.sequences(event, "off") == "")
    check("a BEL or ESC in the text cannot end the OSC early",
          notify.sequences(Event("status", "j", "failed", "x\x07y\x1bz"), "osc", in_tmux=False)
          == "\x1b]9;xyz\x07\x1b]777;notify;orcamon;xyz\x07")


def test_the_geometry_pane_degrades_to_text():
    print("\ngraphics: auto picks text in a multiplexer; the probe and the text renderer")
    import random
    import time as _time
    from orcamon.tui import geometry_text, graphics_probe as gp

    def auto(env, herdr=False, images=True, kitty=True):
        return gp.choose("auto", environ=env, herdr_available=herdr, have_images=images,
                         kitty_probe=lambda: kitty)
    check("auto inside tmux is text", auto({"TMUX": "/tmp/x"}) == ("text", None))
    check("and inside screen", auto({"STY": "1.pts"}) == ("text", None))
    check("herdr wins when it answers", auto({"TMUX": "/tmp/x"}, herdr=True) == ("herdr", None))
    check("outside a multiplexer the terminal is asked", auto({}) == ("kitty", None)
          and auto({}, kitty=False) == ("text", None))
    check("no images extra means text, with the reason",
          auto({}, images=False) == ("text", "images extra not installed"))
    check("--graphics kitty forces it", gp.choose("kitty", environ={"TMUX": "x"}, have_images=True)
          == ("kitty", None))

    check("probe: a Kitty OK before DA1 is kitty",
          gp.parse_probe_reply(b"\x1b_Gi=31;OK\x1b\\\x1b[?62;22c") is True)
    check("probe: DA1 alone is not", gp.parse_probe_reply(b"\x1b[?62;22c") is False)
    check("probe: nothing yet is undecided", gp.parse_probe_reply(b"") is None
          and gp.parse_probe_reply(b"\x1b_Gi=31;OK\x1b\\") is None)

    water = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0)]
    text = geometry_text.render(water, 20, 10).plain
    rows = text.split("\n")
    check("water in 20x10: O and both H drawn, braille between them",
          len(rows) == 10 and all(len(r) == 20 for r in rows) and text.count("O") == 1
          and text.count("H") == 2 and any(0x2800 < ord(c) <= 0x28FF for c in text), repr(text))
    bare = geometry_text.render(water, 20, 10, geometry_text.View(show_labels=False)).plain
    check("labels off: every atom a dot, no symbols, the bonds unchanged",
          bare.count(geometry_text.ATOM_DOT) == 3 and not any(c.isalpha() for c in bare)
          and any(0x2800 < ord(c) <= 0x28FF for c in bare), repr(bare))
    check("H-H is never a bond", sorted(geometry_text.bonds(water)) == [(0, 1), (0, 2)],
          f"{geometry_text.bonds(water)}")
    random.seed(7)
    big = [(random.choice("CCCNOH"), random.uniform(0, 12), random.uniform(0, 12), random.uniform(0, 12))
           for _ in range(300)]
    started = _time.perf_counter()
    geometry_text.render(big, 80, 30)
    elapsed = (_time.perf_counter() - started) * 1000
    check("300 atoms in under 50 ms", elapsed < 50, f"{elapsed:.1f} ms")


def test_the_images_extra_needs_no_matplotlib():
    print("\nthe pixel panes need numpy and Pillow, not matplotlib")
    import importlib.util
    from orcamon.tui import graphics_probe

    real_find_spec = importlib.util.find_spec
    hidden = {"matplotlib"}

    def hiding_find_spec(name, *args, **kwargs):
        if name in hidden:
            return None
        return real_find_spec(name, *args, **kwargs)

    importlib.util.find_spec = hiding_find_spec
    try:
        check("with matplotlib hidden the images extra is still available",
              graphics_probe.images_available() is True)
        hidden = {"PIL"}
        check("with PIL hidden it is not", graphics_probe.images_available() is False)
    finally:
        importlib.util.find_spec = real_find_spec


def test_text_fog_dims_only_the_far_third():
    print("\ntext fog dims only the far third of the depth range")
    from orcamon.core.geometry import View
    from orcamon.tui import geometry_text

    # Depths (eye on +x) are the x coordinates: -3 far, 0 middle, +3 near.
    atoms = [("N", -3.0, 0.0, 0.0), ("O", 0.0, 2.0, 0.0), ("S", 3.0, -2.0, 0.0)]
    text = geometry_text.render(atoms, 30, 15, View(elev=0, azim=0))
    plain = text.plain

    def style_of(symbol):
        index = plain.index(symbol)
        for start, end, style in text.spans:
            if start <= index < end:
                return style
        return None

    far, middle, near = style_of("N"), style_of("O"), style_of("S")
    check("the farthest atom's glyph is dim, the other two are not",
          far is not None and far.dim and not (middle and middle.dim) and not (near and near.dim),
          f"{far} / {middle} / {near}")


def test_text_fog_holds_still_under_rotation():
    print("\ntext fog dims a cell by its depth, not by the frame's own depth range")
    from orcamon.core.geometry import View
    from orcamon.tui import geometry_text

    # The sulfur sits at depth 0 from both views (eye on +x, then on -x).
    # The carbons put the frame's own depth range at [-1, 3] from the first
    # and [-3, 1] from the second, whose far thirds start at +0.33 and -1.67:
    # read off those, the sulfur was dimmed from one side and not the other.
    atoms = [("S", 0.0, 0.0, 2.0), ("C", 3.0, 2.5, 0.0), ("C", -1.0, -2.5, 0.0),
             ("C", -1.0, 2.5, -1.0), ("C", -1.0, -2.5, -1.0)]

    def sulfur_is_dim(azim):
        text = geometry_text.render(atoms, 40, 20, View(elev=0, azim=azim))
        index = text.plain.index("S")
        return any(start <= index < end and style is not None and bool(style.dim)
                   for start, end, style in text.spans)

    front, back = sulfur_is_dim(0), sulfur_is_dim(180)
    check("the sulfur is dimmed the same from both sides", front == back, f"{front} vs {back}")


def test_the_tui_runs_headless():
    print("\nthe TUI, headless: selection, maximize, refresh, a new job, an event")
    import asyncio
    import time as _time
    from orcamon.tui.app import MonitorApp

    async def drive(t, hook_file):
        # The `i` checks need a job whose `NORMAL MODES` block holds an
        # imaginary mode; the tree's own frequency jobs carry none (and
        # adding one to the tree would change ls's row count).
        mode_dir = t.root / "mode_ts"
        mode_dir.mkdir()
        (mode_dir / "job.inp").write_text(_OPT_FREQ)
        vectors = {0: [0.1, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.3]}
        vectors.update({i: [0.0] * 9 for i in range(1, 9)})
        (mode_dir / "job.out").write_text(
            _opt(2, [-80.00, -80.10]) + _CONVERGED + "\n"
            + _freqs([-100.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0])
            + _normal_modes(vectors) + _DONE + "\n")
        app = MonitorApp(t.root, liveness="test", graphics="text", notify_mode="off",
                         on_event=f'echo "$ORCAMON_JOB $ORCAMON_STATUS" >> "{hook_file}"')
        async with app.run_test(size=(180, 50)) as pilot:
            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)) and _time.monotonic() < deadline:
                await pilot.pause(0.05)
            await pilot.press("down")
            await pilot.pause(0.3)
            label = app.selected_label
            job = app.selected_job()
            summary = str(app.query_one("#summary").render())
            check("the summary shows the selected job's label and status",
                  label in summary and job.status.value in summary, summary[:200])
            await pilot.press("m")
            await pilot.pause(0.1)
            check("m maximizes the geometry", app.maximized)
            await pilot.press("m")
            geometry = app.query_one("#geometry")
            geometry.focus()
            await pilot.press("l")
            await pilot.pause(0.2)
            check("l on the geometry pane turns the atom labels off", geometry.show_labels is False)
            await pilot.press("l")
            await pilot.press("v")
            await pilot.pause(0.2)
            check("v cycles the text pane to wireframe, and the title says so",
                  geometry.representation == "wireframe" and "wireframe" in geometry.border_title,
                  f"{geometry.representation} {geometry.border_title!r}")
            await pilot.press("v")
            await pilot.pause(0.2)
            check("and again wraps back to ball-and-stick",
                  geometry.representation == "ball-and-stick", geometry.representation)
            hydrogens = geometry._show_hydrogens()
            await pilot.press("h")
            await pilot.pause(0.2)
            check("h flips the hydrogen flag", geometry._show_hydrogens() is not hydrogens)

            # `i` cycles the imaginary modes: on with the title naming the
            # mode, off again, and nothing at all on a job with no Hessian.
            table = app.query_one("#job_table")
            mode_job = next(j for j in app.jobs if j.label == "mode_ts")
            table.move_cursor(row=app.jobs.index(mode_job))
            deadline = _time.monotonic() + 2
            while geometry._job is not mode_job and _time.monotonic() < deadline:
                await pilot.pause(0.05)
            await pilot.press("i")
            await pilot.pause(0.2)
            check("i starts the imaginary mode and the title names it",
                  geometry.mode_index == 0 and "cm-1" in geometry.border_title,
                  f"{geometry.mode_index} {geometry.border_title!r}")
            await pilot.press("i")
            await pilot.pause(0.2)
            check("and i again turns it off, title and all",
                  geometry.mode_index is None and "cm-1" not in geometry.border_title,
                  f"{geometry.mode_index} {geometry.border_title!r}")
            no_freq = next(j for j in app.jobs if j.label == "opt_maxiter")
            table.move_cursor(row=app.jobs.index(no_freq))
            deadline = _time.monotonic() + 2
            while geometry._job is not no_freq and _time.monotonic() < deadline:
                await pilot.pause(0.05)
            await pilot.press("i")
            await pilot.pause(0.2)
            check("i on a job with no frequencies changes nothing",
                  geometry.mode_index is None, f"{geometry.mode_index}")

            await pilot.press("x")
            await pilot.pause(0.2)
            check("x turns see-through on, and the title says so",
                  geometry.see_through is True and "see-through" in geometry.border_title,
                  f"{geometry.see_through} {geometry.border_title!r}")
            await pilot.press("x")
            await pilot.pause(0.2)
            check("and x again turns it off, title and all",
                  geometry.see_through is False and "see-through" not in geometry.border_title,
                  f"{geometry.see_through} {geometry.border_title!r}")

            await pilot.press("up")
            await pilot.pause(0.1)
            moved = (geometry.elev, geometry.azim)
            await pilot.press("0")
            await pilot.pause(0.1)
            check("0 resets the camera after a rotation",
                  moved != (20.0, -60.0) and geometry.elev == 20.0 and geometry.azim == -60.0
                  and geometry.zoom == 1.0 and geometry.pan == (0.0, 0.0, 0.0),
                  f"{moved} -> {(geometry.elev, geometry.azim, geometry.zoom, geometry.pan)}")
            await pilot.press("o")
            await pilot.pause(0.4)
            check("o rocks: view()'s azim moves while the stored azim does not",
                  geometry._rock_timer is not None and geometry.view().azim != geometry.azim
                  and geometry.azim == -60.0,
                  f"{geometry.view().azim} vs {geometry.azim}")
            offset = geometry.view().azim
            await pilot.pause(0.2)
            check("and the offset keeps changing", geometry.view().azim != offset)
            await pilot.press("o")
            await pilot.pause(0.1)
            check("o again stops rocking and clears the offset",
                  geometry._rock_timer is None and geometry._rock_offset == 0.0)
            app.query_one("#job_table").focus()

            new = t.root / "submitted_later"
            new.mkdir()
            (new / "job.inp").write_text("! B97-3c Opt\n* xyz 0 1\n*\n")
            with open(t.dirs["scan_running"] / "job.out", "a") as f:
                f.write(_DONE + "\n")
            _AliveProbe.ALIVE = set()
            await pilot.press("r")
            deadline = _time.monotonic() + 5
            while _time.monotonic() < deadline and not (
                    any(j.label == "submitted_later" for j in app.jobs) and hook_file.exists()):
                await pilot.pause(0.1)
            check("r finds a job submitted after launch",
                  any(j.label == "submitted_later" for j in app.jobs))
            fired = hook_file.read_text() if hook_file.exists() else ""
            check("and a job finishing runs the --on-event hook with its environment",
                  "scan_running finished" in fired, repr(fired))

    started = _time.monotonic()
    with _Tree() as t:
        asyncio.run(drive(t, t.root.parent / "hook.log"))
    elapsed = _time.monotonic() - started
    check("in under 10 s", elapsed < 10, f"{elapsed:.1f} s")


def test_the_tui_scrubs_a_path():
    print("\nthe TUI charts the IRC path, and left/right/end scrub the geometry pane")
    import asyncio
    import tempfile
    import time as _time
    from pathlib import Path

    from orcamon.tui.app import MonitorApp
    from test_monitor import _irc_tree

    async def drive(root):
        app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
        async with app.run_test(size=(180, 50)) as pilot:
            chart = app.query_one("#chart")
            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)
                       and chart._points) and _time.monotonic() < deadline:
                await pilot.pause(0.05)

            title = lambda: app.query_one("#geometry").border_title
            check("the chart is the IRC path, one dot per point",
                  chart._mode == "irc" and len(chart._points) == 6,
                  f"{chart._mode} {len(chart._points)}")
            check("while following, the geometry pane is on the TS",
                  title().endswith("· IRC TS"), title())
            chart.focus()
            await pilot.press("right")
            await pilot.pause(0.2)
            check("right walks the pane forward along the path",
                  title().endswith("· IRC forward 0"), title())
            await pilot.press("left")
            await pilot.press("left")
            await pilot.pause(0.2)
            check("left walks it back",
                  title().endswith("· IRC backward 0"), title())
            await pilot.press("end")
            await pilot.pause(0.2)
            check("end follows the focus again, which is the TS",
                  title().endswith("· IRC TS"), title())

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _irc_tree(root / "irc")
        asyncio.run(drive(root))


def test_the_tui_cycles_real_modes():
    print("\nthe TUI's I cycles the real modes, lowest frequency first; i keeps the imaginary ones")
    import asyncio
    import tempfile
    import time as _time
    from pathlib import Path

    from orcamon.tui.app import MonitorApp
    from test_monitor import _MODE_FREQS, _MODE_VECTORS

    async def drive(root):
        app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
        async with app.run_test(size=(180, 50)) as pilot:
            geometry = app.query_one("#geometry")
            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)
                       and geometry._job is app.jobs[0]) and _time.monotonic() < deadline:
                await pilot.pause(0.05)
            geometry.focus()
            title = lambda: geometry.border_title

            await pilot.press("I")
            await pilot.pause(0.2)
            check("I starts the real cycle at its lowest mode and the title names it",
                  geometry.mode_index == 5 and "mode 5 50.0 cm-1 (real 1/5)" in title(),
                  f"{geometry.mode_index} {title()!r}")
            await pilot.press("I")
            await pilot.pause(0.2)
            check("I steps on to the next real mode and counts the position",
                  geometry.mode_index == 4 and "(real 2/5)" in title(),
                  f"{geometry.mode_index} {title()!r}")
            await pilot.press("i")
            await pilot.pause(0.2)
            check("i on a real mode starts the imaginary cycle at its first",
                  geometry.mode_index == 3 and "mode 3 -150.0 cm-1 (1/1)" in title(),
                  f"{geometry.mode_index} {title()!r}")
            await pilot.press("I")
            await pilot.pause(0.2)
            check("I on an imaginary mode starts the real cycle at its lowest again",
                  geometry.mode_index == 5, f"{geometry.mode_index} {title()!r}")
            for _ in range(4):
                await pilot.press("I")
                await pilot.pause(0.2)
            check("the real cycle reaches the highest mode and says so",
                  geometry.mode_index == 8 and "mode 8 400.0 cm-1 (real 5/5)" in title(),
                  f"{geometry.mode_index} {title()!r}")
            await pilot.press("I")
            await pilot.pause(0.2)
            check("and I past the last turns the animation off, title and all",
                  geometry.mode_index is None and "cm-1" not in title(),
                  f"{geometry.mode_index} {title()!r}")

            await pilot.press("I")
            geometry._mode_s = 1.0
            drawn = geometry._apply_mode(list(app.selected_job().state.atoms),
                                         app.selected_job())
            rounded = lambda atom: (atom[0], round(atom[1], 6), round(atom[2], 6),
                                    round(atom[3], 6))
            check("the real mode displaces the cycle-2 geometry at full amplitude",
                  rounded(drawn[0]) == ("O", 0.2, 0.0, 0.2)
                  and rounded(drawn[1]) == ("H", 1.31, 0.0, 0.0),
                  f"{drawn[:2]}")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        job_dir = root / "modes_all"
        job_dir.mkdir()
        (job_dir / "job.inp").write_text(_OPT_FREQ)
        (job_dir / "job.out").write_text(
            _opt(2, [-80.00, -80.10]) + _CONVERGED + "\n"
            + _freqs(_MODE_FREQS) + _normal_modes(_MODE_VECTORS) + _DONE + "\n")
        asyncio.run(drive(root))


def test_the_tui_selects_atoms_by_number():
    print("\nthe TUI's s picks atoms by number and the geometry title measures them")
    import asyncio
    import tempfile
    import time as _time
    from pathlib import Path

    from orcamon.tui.app import MonitorApp

    async def drive(root):
        app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
        async with app.run_test(size=(180, 50)) as pilot:
            geometry = app.query_one("#geometry")
            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)
                       and geometry._job is not None
                       and geometry._job.label == "water") and _time.monotonic() < deadline:
                await pilot.pause(0.05)
            geometry.focus()
            title = lambda: geometry.border_title

            await pilot.press("s")
            await pilot.press(*"1 0 2")
            await pilot.press("enter")
            await pilot.pause(0.2)
            check("s and the numbers pick the three atoms and the title measures them",
                  geometry.selection == (1, 0, 2) and "angle 1-0-2: 104.5°" in title(),
                  f"{geometry.selection} {title()!r}")

            chart = app.query_one("#chart")
            chart.focus()
            await pilot.press("left")
            await pilot.pause(0.2)
            check("scrubbing the chart keeps the selection and measures it there",
                  "cycle 1" in title() and "angle 1-0-2: 104.5°" in title()
                  and geometry.selection == (1, 0, 2),
                  f"{geometry.selection} {title()!r}")

            geometry.focus()
            await pilot.press("s")
            for _ in range(5):
                await pilot.press("backspace")
            await pilot.press("7")
            await pilot.press("enter")
            await pilot.pause(0.2)
            check("an atom number the geometry does not have is refused, changing nothing",
                  geometry.selection == (1, 0, 2),
                  f"{geometry.selection} {title()!r}")

            await pilot.press("s")
            for _ in range(5):
                await pilot.press("backspace")
            await pilot.press("enter")
            await pilot.pause(0.2)
            check("an empty entry clears the selection and the title",
                  geometry.selection == () and "angle" not in title(),
                  f"{geometry.selection} {title()!r}")

            await pilot.press("s")
            await pilot.press(*"0 1")
            await pilot.press("enter")
            await pilot.pause(0.2)
            check("two atoms give the distance in the title",
                  "distance 0-1: 0.960 Å" in title(), f"{title()!r}")

            app.query_one("#job_table").focus()
            await pilot.press("down")
            await pilot.pause(0.4)
            check("moving to another job clears the selection",
                  app.selected_label == "water2" and geometry.selection == (),
                  f"{app.selected_label} {geometry.selection}")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _water_job(root, "water")
        _water_job(root, "water2")
        asyncio.run(drive(root))


def test_the_tui_tells_the_agent_what_it_shows():
    print("\nthe TUI publishes what it shows, and orcamon view reads it")
    import asyncio
    import shlex
    import tempfile
    import time as _time
    from pathlib import Path

    from orcamon.core import steer
    from orcamon.tui.app import MonitorApp

    def view():
        return steer.read_json(steer.view_path())

    async def drive(root, file):
        app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
        async with app.run_test(size=(180, 50)) as pilot:
            geometry = app.query_one("#geometry")
            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)
                       and geometry._job is not None
                       and geometry._job.label == "water") and _time.monotonic() < deadline:
                await pilot.pause(0.05)

            check("the TUI clears a dead pid's files when it starts",
                  not steer.view_path(dead).exists() and not steer.focus_path(dead).exists(),
                  f"{steer.view_path(dead)} {steer.focus_path(dead)}")

            document = view()
            check("and publishes the pid, root, job, geometry and as-XYZ command",
                  document is not None
                  and document.get("pid") == os.getpid()
                  and document.get("root") == str(root.resolve())
                  and document.get("job") == {"label": "water", "file": file}
                  and document.get("geometry") == {"shown": "cycle 2", "following": True,
                                                   "args": []}
                  and document.get("selection") == []
                  and document.get("measurement") is None
                  and document.get("mode") is None
                  and document.get("focus") is None
                  and document.get("as_xyz") == shlex.join(["orcamon", "geom", file]),
                  f"{document}")

            geometry.focus()
            await pilot.press("s")
            await pilot.press(*"1 0 2")
            await pilot.press("enter")
            await pilot.pause(0.2)
            document = view()
            check("picking atoms publishes the selection and the measurement",
                  document is not None and document.get("selection") == [1, 0, 2]
                  and document.get("measurement") == {"kind": "angle", "value": 104.4703,
                                                       "unit": "degree"},
                  f"{document}")

            chart = app.query_one("#chart")
            chart.focus()
            await pilot.press("left")
            await pilot.pause(0.2)
            document = view()
            check("scrubbing the chart publishes the geometry and its address",
                  document is not None
                  and document.get("geometry") == {"shown": "cycle 1", "following": False,
                                                   "args": ["--cycle", "1"]}
                  and document.get("as_xyz") == shlex.join(["orcamon", "geom", file,
                                                            "--cycle", "1"]),
                  f"{document}")

            code, out, err = _orcamon(["view"])
            lines = out.strip().split("\n")
            check("`orcamon view` prints the job, geometry, picked atoms and the geom command",
                  code == 0 and lines and lines[0].startswith(f"TUI pid {os.getpid()} on ")
                  and "job        water" in lines
                  and "geometry   cycle 1 (scrubbed)" in lines
                  and "selection  angle 1-0-2: 104.5°" in lines
                  and lines[-1] == "as XYZ     " + shlex.join(["orcamon", "geom", file,
                                                               "--cycle", "1"]),
                  f"{code} {out!r} {err!r}")

            code, doc, err = _json(["view"])
            check("`orcamon view --json` is a new document with the view and its age",
                  code == 0 and doc is not None and doc.get("schema") == 1
                  and doc.get("job", {}).get("label") == "water" and "age_s" in doc,
                  f"{code} {doc} {err}")

    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    dead = proc.pid
    saved = os.environ.get("XDG_STATE_HOME")
    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_STATE_HOME"] = tmp
        try:
            root = Path(tmp) / "tree"
            root.mkdir()
            _water_job(root, "water")
            file = str((root / "water" / "job.inp").resolve())
            steer.write_json(steer.view_path(dead), {"pid": dead, "updated_unix_s": 1.0})
            steer.write_json(steer.focus_path(dead), {})
            asyncio.run(drive(root, file))
            code, out, err = _orcamon(["view"])
            check("a TUI that exits removes its file, and view then exits 2",
                  steer.view_path().exists() is False and code == 2
                  and "no orcamon TUI is running" in err,
                  f"{steer.view_path().exists()} {code} {out!r} {err!r}")
        finally:
            if saved is None:
                os.environ.pop("XDG_STATE_HOME", None)
            else:
                os.environ["XDG_STATE_HOME"] = saved


def test_view_names_the_geometry_its_numbers_belong_to():
    print("\nview names the geometry its atom numbers belong to, in full")
    import tempfile
    import time as _time

    from orcamon.core import steer

    file = "/r/irc/job.inp"
    saved = os.environ.get("XDG_STATE_HOME")
    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_STATE_HOME"] = tmp
        try:
            base = {
                "schema": 1, "pid": os.getpid(), "host": steer.host_name(),
                "root": "/r", "updated_unix_s": _time.time(),
                "job": {"label": "irc", "file": file},
                "geometry": {"shown": "IRC TS", "following": True,
                             "args": ["--point", "0"]},
                "atoms_shown": 137, "mode": None,
                "selection": [133, 134, 135],
                "selection_atoms": [{"atom": 133, "element": "O"},
                                    {"atom": 134, "element": "O"},
                                    {"atom": 135, "element": "H"}],
                "measurement": {"kind": "angle", "value": 98.1764, "unit": "degree"},
                "as_measure": f"orcamon measure {file} 133 134 135 --point 0",
                "as_xyz": f"orcamon geom {file} --point 0",
                "focus": None,
            }

            def view():
                code, out, err = _orcamon(["view"])
                return code, out.strip().split("\n"), err

            steer.write_json(steer.view_path(), base)
            code, lines, err = view()
            check("view names the geometry and says the pane follows the job",
                  code == 0 and "geometry   IRC TS (following the job)" in lines,
                  f"{code} {lines} {err!r}")

            check("view carries the atom count of the geometry it names",
                  "atoms      137 in this geometry" in lines, f"{lines}")

            picked = "selection  angle 133-134-135: 98.2°"
            check("view prints the picked atoms' elements right under the measurement",
                  picked in lines and lines[lines.index(picked) + 1] == "elements   O-O-H",
                  f"{lines}")

            check("view hands over both ready commands, measure above XYZ",
                  lines[-2] == f"as measure orcamon measure {file} 133 134 135 --point 0"
                  and lines[-1] == f"as XYZ     orcamon geom {file} --point 0",
                  f"{lines}")

            base["geometry"] = {**base["geometry"], "following": False}
            steer.write_json(steer.view_path(), base)
            code, lines, err = view()
            check("a scrubbed pane says scrubbed, not following",
                  code == 0 and "geometry   IRC TS (scrubbed)" in lines,
                  f"{code} {lines} {err!r}")

            older = {k: v for k, v in base.items()
                     if k not in ("atoms_shown", "selection_atoms", "as_measure")}
            steer.write_json(steer.view_path(), older)
            code, lines, err = view()
            check("an older TUI's document prints exactly as before",
                  code == 0
                  and not any(line.startswith(("atoms ", "elements ", "as measure"))
                              for line in lines)
                  and lines[-1] == f"as XYZ     orcamon geom {file} --point 0",
                  f"{code} {lines} {err!r}")
        finally:
            if saved is None:
                os.environ.pop("XDG_STATE_HOME", None)
            else:
                os.environ["XDG_STATE_HOME"] = saved


def test_the_tui_publishes_the_numbering_it_shows():
    print("\nthe TUI publishes the numbering it shows, beside the geometry")
    import asyncio
    import shlex
    import tempfile
    import time as _time
    from pathlib import Path

    from orcamon.core import steer
    from orcamon.tui.app import MonitorApp

    async def drive(root, file):
        app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
        async with app.run_test(size=(180, 50)) as pilot:
            geometry = app.query_one("#geometry")
            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)
                       and geometry._job is not None
                       and geometry._job.label == "water") and _time.monotonic() < deadline:
                await pilot.pause(0.05)

            doc = lambda: steer.read_json(steer.view_path())

            geometry.focus()
            await pilot.press("s")
            await pilot.press(*"1 0 2")
            await pilot.press("enter")
            await pilot.pause(0.2)
            check("the TUI publishes the atom count of the geometry it draws",
                  (doc() or {}).get("atoms_shown") == 3, f"{doc()}")

            check("and the element of every picked atom at the number it shows",
                  (doc() or {}).get("selection_atoms") == [
                      {"atom": 1, "element": "H"}, {"atom": 0, "element": "O"},
                      {"atom": 2, "element": "H"}],
                  f"{doc()}")

            check("and a measure command for exactly those atoms on that geometry",
                  (doc() or {}).get("as_measure")
                  == shlex.join(["orcamon", "measure", file, "1", "0", "2"]),
                  f"{doc()}")

            chart = app.query_one("#chart")
            chart.focus()
            await pilot.press("left")
            await pilot.pause(0.2)
            check("scrubbing carries the geometry's own --cycle into the measure command",
                  (doc() or {}).get("as_measure")
                  == shlex.join(["orcamon", "measure", file, "1", "0", "2", "--cycle", "1"]),
                  f"{doc()}")

            geometry.focus()
            await pilot.press("s")
            for _ in range(5):
                await pilot.press("backspace")
            await pilot.press("enter")
            await pilot.pause(0.2)
            check("clearing the selection clears the elements and the command",
                  (doc() or {}).get("as_measure") is None
                  and (doc() or {}).get("selection_atoms") == [],
                  f"{doc()}")

    saved = os.environ.get("XDG_STATE_HOME")
    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_STATE_HOME"] = tmp
        try:
            root = Path(tmp) / "tree"
            root.mkdir()
            _water_job(root, "water")
            file = str((root / "water" / "job.inp").resolve())
            asyncio.run(drive(root, file))
        finally:
            if saved is None:
                os.environ.pop("XDG_STATE_HOME", None)
            else:
                os.environ["XDG_STATE_HOME"] = saved


def test_focus_waits_for_the_tui_to_take_it():
    print("\nfocus sends a checked request and waits for the TUI to take it")
    import tempfile
    import threading
    import time as _time
    from pathlib import Path

    from orcamon.core import steer

    saved = os.environ.get("XDG_STATE_HOME")
    with tempfile.TemporaryDirectory() as tmp, tempfile.TemporaryDirectory() as other:
        root = Path(tmp) / "tree"
        root.mkdir()
        _water_job(root, "water")
        cli = lambda *a: _orcamon(["--root", str(root), "--liveness", "mtime", "--no-cache", *a])
        os.environ["XDG_STATE_HOME"] = tmp
        try:
            code, out, err = cli("focus", "water")
            check("with no TUI running focus refuses and says so",
                  code == 2 and "no orcamon TUI is running" in err, f"{code} {out!r} {err!r}")

            fake = {"schema": 1, "pid": os.getpid(), "host": steer.host_name(),
                    "root": str(root.resolve()), "updated_unix_s": _time.time(),
                    "job": None, "geometry": None, "mode": None, "selection": [],
                    "measurement": None, "as_xyz": None, "focus": None}
            taken = []

            def ack(ok, message):
                # A fake TUI: take the request when it appears and publish the
                # answer through the view file, as MonitorApp will.
                def run():
                    deadline = _time.monotonic() + 3
                    while (not steer.focus_path(os.getpid()).exists()
                           and _time.monotonic() < deadline):
                        _time.sleep(0.02)
                    req = steer.take_focus(os.getpid())
                    taken.append(req)
                    fake["focus"] = {"id": req["id"], "ok": ok, "message": message}
                    fake["updated_unix_s"] = _time.time()
                    steer.write_json(steer.view_path(), fake)
                threading.Thread(target=run).start()

            steer.write_json(steer.view_path(), fake)
            ack(True, "showing water (test)")
            code, out, err = cli("focus", "water", "--atoms", "1,0", "--timeout", "3")
            req = taken[0] if taken else None
            check("the TUI takes the checked request and focus reports what it showed",
                  code == 0
                  and f"showing water (test) in the TUI (pid {os.getpid()})" in out
                  and req is not None
                  and req.get("file") == str((root / "water" / "job.inp").resolve())
                  and req.get("atoms") == [1, 0]
                  and req.get("point") is None and req.get("mode") is None,
                  f"{code} {out!r} {err!r} {req!r}")

            ack(False, "nope")
            code, out, err = cli("focus", "water", "--timeout", "3")
            check("a TUI that could not show it exits 2 with its message",
                  code == 2 and "could not show it: nope" in err, f"{code} {out!r} {err!r}")

            code, out, err = cli("focus", "water", "--timeout", "0.3")
            check("a TUI that does not take it times out, exit 4, and it stays queued",
                  code == 4 and "has not taken the request within 0.3 s" in err
                  and steer.focus_path(os.getpid()).exists(),
                  f"{code} {out!r} {err!r}")

            code, out, err = cli("focus", "water", "--timeout", "0")
            check("--timeout 0 queues the request and returns",
                  code == 0 and "queued for the TUI" in out, f"{code} {out!r} {err!r}")

            code, out, err = cli("focus", "water", "--atoms", "0,9")
            bad_atoms = code == 2 and "atom 9 is out of range: the geometry has 3 atoms (0-2)" in err
            code, out, err = cli("focus", "water", "--mode", "3")
            no_freq = code == 2 and "this job has no frequency block" in err
            code, out, err = cli("focus", "water", "--step", "1", "--cycle", "1")
            two_axes = code == 2 and "focus takes one of --point, --step or --cycle" in err
            code, out, err = cli("focus", "water", "--point", "0")
            no_path = code == 2 and "--point needs an IRC or NEB job" in err
            check("the geometry, the atoms and the mode are checked before anything is sent",
                  bad_atoms and no_freq and two_axes and no_path,
                  f"{bad_atoms} {no_freq} {two_axes} {no_path}")

            fake["root"] = str(Path(other).resolve())
            fake["updated_unix_s"] = _time.time()
            steer.write_json(steer.view_path(), fake)
            code, out, err = cli("focus", "water", "--timeout", "0")
            check("a job outside what the TUI watches is refused, not queued",
                  code == 2 and "outside it" in err, f"{code} {out!r} {err!r}")
        finally:
            if saved is None:
                os.environ.pop("XDG_STATE_HOME", None)
            else:
                os.environ["XDG_STATE_HOME"] = saved


def test_the_tui_takes_the_agents_focus():
    print("\nthe TUI takes the agent's orcamon focus, and b goes back")
    import asyncio
    import tempfile
    import time as _time
    from pathlib import Path

    from orcamon.core import steer
    from orcamon.tui.app import MonitorApp
    from test_monitor import _MODE_FREQS, _MODE_VECTORS, _irc_tree

    async def drive(root):
        app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
        async with app.run_test(size=(180, 50)) as pilot:
            geometry = app.query_one("#geometry")
            cli = lambda *a: _orcamon(["--root", str(root), "--liveness", "mtime",
                                       "--no-cache", *a])
            title = lambda: geometry.border_title

            async def until(cond, seconds):
                deadline = _time.monotonic() + seconds
                while not cond() and _time.monotonic() < deadline:
                    await pilot.pause(0.1)
                return cond()

            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)
                       and geometry._job is not None
                       and geometry._job.label == "irc") and _time.monotonic() < deadline:
                await pilot.pause(0.05)

            code, out, err = cli("focus", "modes_all", "--atoms", "1,0,2",
                                 "--mode", "5", "--timeout", "0")
            check("focus shows the job, its atoms and its mode, and the title measures them",
                  code == 0 and await until(
                      lambda: app.selected_label == "modes_all"
                      and geometry.selection == (1, 0, 2)
                      and geometry.mode_index == 5
                      and "mode 5 50.0 cm-1 (real 1/5)" in title()
                      and "angle 1-0-2: 104.5°" in title(), 3),
                  f"{code} {app.selected_label} {geometry.selection} "
                  f"{geometry.mode_index} {title()!r} {out!r} {err!r}")

            focus = (steer.read_json(steer.view_path()) or {}).get("focus")
            check("the TUI records the focus it took for the agent to read back",
                  isinstance(focus, dict) and focus.get("ok") is True
                  and focus.get("message") == "showing modes_all",
                  f"{focus}")

            code, out, err = cli("focus", "irc", "--point", "1", "--timeout", "0")
            check("focus on a path point scrubs the chart and the pane follows it",
                  code == 0 and await until(
                      lambda: app.selected_label == "irc"
                      and title().endswith("· IRC forward 0")
                      and (steer.read_json(steer.view_path()) or {}).get(
                          "geometry", {}).get("args") == ["--point", "1"], 3),
                  f"{code} {app.selected_label} {title()!r} {out!r} {err!r}")

            await pilot.press("b")
            await pilot.pause(0.4)
            check("b returns to the job, atoms and mode the agent interrupted",
                  app.selected_label == "modes_all" and geometry.selection == (1, 0, 2)
                  and geometry.mode_index == 5,
                  f"{app.selected_label} {geometry.selection} {geometry.mode_index} "
                  f"{title()!r}")

            await pilot.press("b")
            await pilot.pause(0.2)
            check("b with nothing left to go back to changes nothing",
                  app.selected_label == "modes_all",
                  f"{app.selected_label} {title()!r}")

            (root / "draft").mkdir()
            (root / "draft" / "job.inp").write_text(
                "! B97-3c Opt\n* xyz 0 1\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\n"
                "H -0.24 0.93 0.0\n*\n")
            code, out, err = cli("focus", "draft", "--atoms", "0,1", "--timeout", "0")
            check("focus finds a draft input written after the TUI started",
                  code == 0 and await until(
                      lambda: app.selected_label == "draft"
                      and "input geometry" in title()
                      and "distance 0-1: 0.960 Å" in title(), 8),
                  f"{code} {app.selected_label} {title()!r} "
                  f"{[j.label for j in app.jobs]} {out!r} {err!r}")

    saved = os.environ.get("XDG_STATE_HOME")
    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_STATE_HOME"] = tmp
        try:
            root = Path(tmp) / "tree"
            root.mkdir()
            _irc_tree(root / "irc")
            job_dir = root / "modes_all"
            job_dir.mkdir()
            (job_dir / "job.inp").write_text(_OPT_FREQ)
            (job_dir / "job.out").write_text(
                _opt(2, [-80.00, -80.10]) + _CONVERGED + "\n"
                + _freqs(_MODE_FREQS) + _normal_modes(_MODE_VECTORS) + _DONE + "\n")
            asyncio.run(drive(root))
        finally:
            if saved is None:
                os.environ.pop("XDG_STATE_HOME", None)
            else:
                os.environ["XDG_STATE_HOME"] = saved


def test_b_does_not_leave_the_agents_mode_running():
    print("\nb stops a mode the agent started when no mode ran before the focus")
    import asyncio
    import tempfile
    import time as _time
    from pathlib import Path

    from orcamon.tui.app import MonitorApp
    from test_monitor import _MODE_FREQS, _MODE_VECTORS

    async def drive(root):
        app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
        async with app.run_test(size=(180, 50)) as pilot:
            geometry = app.query_one("#geometry")
            cli = lambda *a: _orcamon(["--root", str(root), "--liveness", "mtime",
                                       "--no-cache", *a])
            title = lambda: geometry.border_title

            async def until(cond, seconds):
                deadline = _time.monotonic() + seconds
                while not cond() and _time.monotonic() < deadline:
                    await pilot.pause(0.1)
                return cond()

            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)
                       and geometry._job is not None
                       and geometry._job.label == "modes_all") and _time.monotonic() < deadline:
                await pilot.pause(0.05)

            code, out, err = cli("focus", "modes_all", "--atoms", "1,0", "--mode", "5",
                                 "--timeout", "0")
            check("focus starts the agent's mode 5 on the picked atoms",
                  code == 0 and await until(
                      lambda: geometry.mode_index == 5 and geometry.selection == (1, 0), 3),
                  f"{code} {geometry.mode_index} {geometry.selection} {out!r} {err!r}")

            await pilot.press("b")
            await pilot.pause(0.4)
            check("b stops the agent's mode: no mode, no selection, no cm-1 in the title",
                  geometry.mode_index is None and geometry.selection == ()
                  and "cm-1" not in title(),
                  f"{geometry.mode_index} {geometry.selection} {title()!r}")

            geometry.focus()
            await pilot.press("i")
            await pilot.pause(0.2)
            person = geometry.mode_index
            code, out, err = cli("focus", "modes_all", "--mode", "5", "--timeout", "0")
            check("with the person's own mode 3 running, focus starts mode 5",
                  person == 3 and code == 0 and await until(
                      lambda: geometry.mode_index == 5, 3),
                  f"{person} {code} {geometry.mode_index} {out!r} {err!r}")

            await pilot.press("b")
            await pilot.pause(0.4)
            check("b returns to the person's own mode 3, title and all",
                  geometry.mode_index == 3
                  and "mode 3 -150.0 cm-1 (1/1)" in title(),
                  f"{geometry.mode_index} {title()!r}")

    saved = os.environ.get("XDG_STATE_HOME")
    with tempfile.TemporaryDirectory() as tmp:
        os.environ["XDG_STATE_HOME"] = tmp
        try:
            root = Path(tmp) / "tree"
            root.mkdir()
            job_dir = root / "modes_all"
            job_dir.mkdir()
            (job_dir / "job.inp").write_text(_OPT_FREQ)
            (job_dir / "job.out").write_text(
                _opt(2, [-80.00, -80.10]) + _CONVERGED + "\n"
                + _freqs(_MODE_FREQS) + _normal_modes(_MODE_VECTORS) + _DONE + "\n")
            asyncio.run(drive(root))
        finally:
            if saved is None:
                os.environ.pop("XDG_STATE_HOME", None)
            else:
                os.environ["XDG_STATE_HOME"] = saved


def test_the_tui_picks_atoms_by_click():
    print("\na click in the geometry pane toggles the atom under it")
    import asyncio
    import tempfile
    import time as _time
    from pathlib import Path

    from orcamon.core.geometry import pick_atom, screen_positions
    from orcamon.tui.app import MonitorApp

    async def drive(root):
        app = MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")
        async with app.run_test(size=(180, 50)) as pilot:
            geometry = app.query_one("#geometry")
            deadline = _time.monotonic() + 5
            while not (app.jobs and all(j.parsed for j in app.jobs)
                       and geometry._job is not None
                       and geometry._job.label == "water") and _time.monotonic() < deadline:
                await pilot.pause(0.05)
            title = lambda: geometry.border_title
            width, height = geometry.content_size.width, geometry.content_size.height
            atoms = geometry._atoms_for(geometry._job, geometry._point)
            frame = geometry._pick_frame(width, height)
            positions = screen_positions(atoms, geometry.view(), frame)

            def cell(i):
                return (int(positions[i][0] / frame[0] * width),
                        int(positions[i][1] / frame[1] * height))

            async def click(c):
                # + 1: the click offset is to the pane's outer region, and the
                # first content cell sits one cell inside the round border.
                await pilot.click("#geometry", offset=(c[0] + 1, c[1] + 1))
                await pilot.pause(0.2)

            await click(cell(2))
            check("a click on an atom picks it",
                  geometry.selection == (2,),
                  f"{geometry.selection} {cell(2)} {width}x{height}")

            await click(cell(2))
            check("a second click on the same atom drops it",
                  geometry.selection == (), f"{geometry.selection}")

            await click(cell(0))
            await click(cell(1))
            check("clicking a second atom adds it, and the title measures the pair",
                  geometry.selection == (0, 1) and "distance 0-1: 0.960 Å" in title(),
                  f"{geometry.selection} {title()!r}")

            empty = pick_atom(atoms, geometry.view(), frame, (width, height), (0, 0))
            await click((0, 0))
            check("a click on nothing leaves the selection alone",
                  empty is None and geometry.selection == (0, 1),
                  f"{empty} {geometry.selection}")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _water_job(root, "water")
        asyncio.run(drive(root))


def test_edit_writes_a_checked_geometry():
    print("\nedit: a checked copy of a geometry, with the refusals that matter")
    import tempfile
    from pathlib import Path

    from orcamon.core.geometry import read_xyz

    F = [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0), ("H", -0.5, 1.0, 0.0),
         ("O", 2.0, 1.0, 0.0), ("H", 2.0, 1.0, 1.0)]
    G = [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0), ("O", 2.1, 1.2, 0.0),
         ("O", 2.1, -1.9, 0.0), ("H", -0.5, 0.9, 0.0), ("H", -0.5, -0.9, 0.0)]

    def block(atoms):
        return ("CARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------\n"
                + "".join(f"  {el}      {x:.6f}    {y:.6f}    {z:.6f}\n"
                          for el, x, y, z in atoms))

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        out = root / "out"
        out.mkdir()
        for name, atoms in (("f", F), ("g", G)):
            d = root / name
            d.mkdir()
            (d / "job.inp").write_text(
                "! B97-3c Opt\n* xyz 0 1\n"
                + "".join(f"{el} {x} {y} {z}\n" for el, x, y, z in atoms) + "*\n")
        ml = root / "ml"
        ml.mkdir()
        (ml / "job.inp").write_text(
            "! QM/XTB B97-3c Opt\n%qmmm\n  QMAtoms {0 1} end\n"
            "  Charge_Total 0\n  Mult_Total 1\nend\n* xyz 0 1\n*\n")
        (ml / "job.out").write_text("QM1 Subsystem                    ...  0 1\n"
                                    + block(G) + "\n" + _DONE + "\n")
        argv = ["--root", str(root), "--liveness", "mtime", "--no-cache"]

        code, out_text, err = _orcamon([*argv, "edit", "f", "--rotate-bond", "0", "1", "90",
                                        "-o", str(out / "r.xyz")])
        check("a rotation writes the file and names the motion, the atoms and the contact",
              code == 0 and out_text.splitlines() == [
                  "rotate bond 0-1 by +90.0°",
                  "moved 2 of 5 atoms, the 0 side; --move 1 moves the other 3",
                  "closest new contact: 2 H-3 O 2.872 Å",
                  f"wrote {out / 'r.xyz'} (5 atoms)"],
              f"{code} {(out_text + err)!r}")

        atoms = read_xyz(out / "r.xyz")
        check("the new XYZ holds the turned geometry and says where it came from",
              atoms is not None and len(atoms) == 5
              and (atoms[2][0],) + tuple(round(v, 6) for v in atoms[2][1:]) == ("H", -0.5, 0.0, -1.0)
              and (out / "r.xyz").read_text().splitlines()[1]
              == "f · input geometry · rotate bond 0-1 by +90.0°",
              f"{atoms}")

        code, doc, err = _json([*argv, "edit", "f", "--set-dihedral", "2", "0", "1", "3", "-60",
                                "-o", str(out / "d.xyz")])
        check("--json reports the value reached, the moved atoms and the file written",
              doc is not None and code == 0 and doc["after"] == -60.0 and doc["moved"] == [0, 2]
              and doc["written"] == str(out / "d.xyz") and doc["new_bonds"] == [],
              f"{code} {doc} {err}")

        code, out_text, err = _orcamon([*argv, "edit", "g", "--rotate-bond", "0", "1", "180",
                                        "-o", str(out / "x.xyz")])
        check("a new bond is named and refused, and nothing is written",
              code == 2 and "new bond: 2 O-3 O 0.700 Å" in out_text
              and "refused: the edit makes 1 new bond; nothing written "
                  "(--allow-new-bonds writes it anyway)" in err
              and not (out / "x.xyz").exists(),
              f"{code} {(out_text + err)!r}")

        code, out_text, err = _orcamon([*argv, "edit", "g", "--rotate-bond", "0", "1", "180",
                                        "--allow-new-bonds", "-o", str(out / "x.xyz")])
        check("--allow-new-bonds writes it anyway",
              code == 0 and (out / "x.xyz").exists(), f"{code} {(out_text + err)!r}")

        before = (out / "r.xyz").read_text()
        code, out_text, err = _orcamon([*argv, "edit", "f", "--rotate-bond", "0", "1", "10",
                                        "-o", str(out / "r.xyz")])
        check("an existing file is refused and left untouched",
              code == 2 and "exists: orcamon edit never overwrites a file" in err
              and (out / "r.xyz").read_text() == before,
              f"{code} {(out_text + err)!r}")

        code_a, out_a, err_a = _orcamon([*argv, "edit", "f", "--rotate-bond", "0", "3", "10",
                                         "-o", str(out / "n.xyz")])
        code_b, out_b, err_b = _orcamon([*argv, "edit", "f", "--rotate-bond", "0", "1", "abc",
                                         "-o", str(out / "n.xyz")])
        check("atoms that are not bonded, and a value that is not a number, are refused",
              code_a == 2 and "atoms 0 and 3 are not bonded: 2.24 Å apart" in err_a
              and code_b == 2 and "not a number: 'abc'" in err_b,
              f"{code_a} {err_a!r} {code_b} {err_b!r}")

        code, out_text, err = _orcamon([*argv, "edit", "ml", "--rotate-bond", "0", "1", "90",
                                        "-o", str(out / "m.xyz")])
        check("a multilayer job splits the moved atoms into QM and environment",
              code == 0 and out_text.splitlines()[1]
              == "moved 2 of 6 atoms, the 1 side (1 QM, 1 environment); "
                 "--move 0 moves the other 3",
              f"{code} {(out_text + err)!r}")


def test_check_finds_what_would_fail():
    print("\ncheck: the input faults that would waste a run, before it runs")
    import tempfile
    from pathlib import Path

    water = "O 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        argv = ["--root", str(root), "--liveness", "mtime", "--no-cache"]

        def draft(name, body, out=False):
            d = root / name
            d.mkdir()
            (d / "job.inp").write_text(body)
            if out:
                (d / "job.out").write_text("")

        draft("water", "! B97-3c Opt\n* xyz 0 1\n" + water + "*\n")
        code, out_text, err = _orcamon([*argv, "check", "water"])
        check("a clean input exits 0 and counts the whole system",
              code == 0 and out_text.splitlines()[-1]
              == "0 error(s), 0 warning(s) · 3 atoms, 10 electrons, charge 0, multiplicity 1",
              f"{code} {(out_text + err)!r}")

        draft("doublet", "! B97-3c Opt\n* xyz 0 2\n" + water + "*\n")
        code, out_text, err = _orcamon([*argv, "check", "doublet"])
        check("an even electron count cannot have an odd multiplicity",
              code == 1 and "error: 10 electrons cannot have multiplicity 2: "
                            "it needs an odd number" in out_text,
              f"{code} {(out_text + err)!r}")

        draft("scan", "! B97-3c Opt\n%geom\n  Scan\n    B 0 5 = 1.0, 2.0, 3\n  end\nend\n"
                      "* xyz 0 1\n" + water + "*\n")
        code, out_text, err = _orcamon([*argv, "check", "scan"])
        check("an atom number past the geometry names the group it came from",
              code == 1 and "error: atom 5 in `Scan B 0 5` is out of range: "
                            "the geometry has 3 atoms (0-2)" in out_text,
              f"{code} {(out_text + err)!r}")

        draft("missing", "! B97-3c Opt\n* xyzfile 0 1 missing.xyz\n")
        code, out_text, err = _orcamon([*argv, "check", "missing"])
        check("a geometry file that is not there is an error, not a silent skip",
              code == 1 and "error: missing.xyz: missing or not an XYZ file" in out_text,
              f"{code} {(out_text + err)!r}")

        draft("pal", "! B97-3c Opt\n%pal nprocs 100000 end\n* xyz 0 1\n" + water + "*\n")
        code, out_text, err = _orcamon([*argv, "check", "pal"])
        check("more processes than this machine has is a warning, not an error",
              code == 0 and any("warning:" in line
                                and "nprocs 100000 is more than this machine's" in line
                                for line in out_text.splitlines()),
              f"{code} {(out_text + err)!r}")

        draft("again", "! B97-3c Opt\n* xyz 0 1\n" + water + "*\n", out=True)
        code, doc, err = _json([*argv, "check", "again"])
        check("an output a run would overwrite is reported in the JSON",
              code == 0 and doc is not None and {"code": "output_exists",
                  "message": "job.out exists: running this job again overwrites it"}
              in doc["warnings"], f"{code} {doc} {err}")

        draft("boundary", "! QM/XTB B97-3c Opt\n%qmmm\n  QMAtoms {0 1} end\n"
                          "  Charge_Total 0\n  Mult_Total 1\nend\n* xyz 0 1\n" + water + "*\n")
        code, out_text, err = _orcamon([*argv, "check", "boundary"])
        check("a bond to a hydrogen across the QM boundary is named",
              code == 0 and "warning: QM atom 0 O is 0.96 Å from environment atom 2 H: "
                            "a bond to hydrogen crosses the QM boundary" in out_text,
              f"{code} {(out_text + err)!r}")

        code, doc, err = _json([*argv, "check", "doublet"])
        check("--json reports the parity code and the electron count",
              doc is not None and code == 1 and doc["errors"][0]["code"] == "parity"
              and doc["electrons"] == 10, f"{code} {doc} {err}")


def test_connectivity_names_bonds_made_and_broken():
    print("\nconnectivity: the bonds a run formed and broke, from what it started with")
    import tempfile
    from pathlib import Path

    from test_monitor import _irc_tree

    G = [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0), ("O", 2.1, 1.2, 0.0),
         ("O", 2.1, -1.9, 0.0), ("H", -0.5, 0.9, 0.0), ("H", -0.5, -0.9, 0.0)]
    G2 = [G[0], G[1], ("O", 2.1, 2.0, 0.0), ("O", 2.9, 2.6, 0.0), G[4], G[5]]
    water = "O 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"

    def block(atoms):
        return ("CARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------\n"
                + "".join(f"  {el}      {x:.6f}    {y:.6f}    {z:.6f}\n"
                          for el, x, y, z in atoms))

    def output(atoms):
        return "\n".join([_CYCLE.format(n=1), block(G),
                          "FINAL SINGLE POINT ENERGY     -80.000000000",
                          _CYCLE.format(n=2), block(atoms),
                          "FINAL SINGLE POINT ENERGY     -80.100000000", _DONE]) + "\n"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for name, atoms in (("c", G2), ("short", G[:2])):
            d = root / name
            d.mkdir()
            (d / "job.inp").write_text("! B97-3c Opt\n* xyz 0 1\n*\n")
            (d / "job.out").write_text(output(atoms))
        draft = root / "draft"
        draft.mkdir()
        (draft / "job.inp").write_text("! B97-3c Opt\n* xyz 0 1\n" + water + "*\n")
        _irc_tree(root / "irc")
        argv = ["--root", str(root), "--liveness", "mtime", "--no-cache"]

        code, text, err = _orcamon([*argv, "connectivity", "c"])
        lines = text.splitlines()
        check("the changes are named from the first geometry to the latest, with both distances",
              code == 0 and lines[0] == "first printed geometry -> cycle 2: 1 bond formed, 1 broken"
              and "formed  2 O-3 O  3.10 -> 1.00 Å" in lines
              and "broken  1 C-2 O  1.34 -> 2.09 Å" in lines
              and lines[-1] == "bonds: under 1.7 Å between heavy atoms, 1.4 Å to a hydrogen, "
                               "1.0 Å between hydrogens",
              f"{code} {text!r}")

        code, text, err = _orcamon([*argv, "connectivity", "c", "--cycle", "1"])
        check("a cycle with no change reports only the count of bonds",
              code == 0 and text.splitlines()[0]
              == "first printed geometry -> cycle 1: no bond formed or broken (4 bonds)",
              f"{code} {text!r}")

        code, doc, err = _json([*argv, "connectivity", "c"])
        check("--json gives the same changes as atom numbers and rounded distances",
              doc is not None
              and doc["formed"] == [{"atoms": [2, 3], "elements": ["O", "O"],
                                     "before": 3.1, "after": 1.0}]
              and doc["broken"] == [{"atoms": [1, 2], "elements": ["C", "O"],
                                     "before": 1.342, "after": 2.088}]
              and doc["bonds"] == 4, f"{code} {doc} {err}")

        code, text, err = _orcamon([*argv, "connectivity", "irc", "--point", "1"])
        check("on a path, point 0 (the TS) is the baseline",
              code == 0 and text.splitlines()[0]
              == "IRC TS -> IRC forward 0: no bond formed or broken (2 bonds)",
              f"{code} {text!r}")

        code, text, err = _orcamon([*argv, "connectivity", "draft"])
        check("an output that printed no geometry has nothing to compare",
              code == 2 and "this output has printed no geometry yet" in err,
              f"{code} {text!r} {err!r}")

        code, text, err = _orcamon([*argv, "connectivity", "short"])
        check("two printed parts of the system are refused, not compared atom for atom",
              code == 2 and "the first printed geometry has 6 atoms and cycle 2 has 2" in err,
              f"{code} {text!r} {err!r}")


def test_slurm_names_what_the_login_node_cannot_see():
    print("\nSLURM: squeue's view of a job, and the fallbacks when it has none")
    import stat
    import time as _time
    from orcamon.core import liveness

    def reset():
        liveness.SlurmProbe._cache = None
        liveness.SlurmProbe._dirs_cache = {}
        liveness.SlurmProbe._warned = False

    with _Tree() as t:
        bin_dir = t.root.parent / "bin"
        bin_dir.mkdir()
        fake = bin_dir / "squeue"
        fake.write_text('#!/bin/sh\ncat "$FAKE_SQUEUE_OUT"\nexit "${FAKE_SQUEUE_EXIT:-0}"\n')
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
        unfinished = _opt(2, [-5.0, -5.1], converged=False)
        for rel in ("batch/sub", "multi/a", "multi/b", "lonely"):
            d = t.root / rel
            d.mkdir(parents=True)
            (d / "job.inp").write_text(_OPT_FREQ)
            (d / "job.out").write_text(unfinished)
        canned = t.root.parent / "squeue.out"
        canned.write_text("\n".join([
            f"101|PENDING|{t.dirs['not_yet']}|N/A",
            f"102|RUNNING|{t.dirs['scan_running']}|2026-09-28T10:00:00",
            f"103|RUNNING|{t.root / 'batch'}|2026-09-28T10:00:00",
            f"104|RUNNING|{t.root / 'multi'}|2026-09-28T10:00:00",
            "this line is garbage",
            f"105|COMPLETED|{t.dirs['opt_done']}|2026-09-28T09:00:00",
        ]) + "\n")
        saved = {k: os.environ.get(k) for k in ("PATH", "FAKE_SQUEUE_OUT", "FAKE_SQUEUE_EXIT")}
        os.environ["PATH"] = f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
        os.environ["FAKE_SQUEUE_OUT"] = str(canned)
        slurm = ["--root", str(t.root), "--liveness", "slurm"]
        try:
            reset()
            check("auto picks slurm when squeue is on PATH", liveness.make_probe("auto").name == "slurm")

            def status(name, argv=slurm):
                code, doc, err = _json([*argv, "show", name])
                return (doc or {}).get("job", {}), err

            job, _ = status("not_yet")
            check("PENDING is queued, with the SLURM job id", job.get("status") == "queued"
                  and job.get("sched_id") == "101", f"{job.get('status')} {job.get('sched_id')}")
            job, _ = status("scan_running")
            check("RUNNING at the job's own directory is running", job.get("status") == "running"
                  and job.get("sched_id") == "102" and job.get("liveness_source") == "slurm")
            job, _ = status("batch/sub")
            check("a working directory above exactly one unfinished job matches it",
                  job.get("status") == "running" and job.get("sched_id") == "103", f"{job}")
            job, _ = status("multi/a")
            code, out, err = _orcamon([*slurm, "show", "multi/a"])
            check("above two unfinished jobs it matches neither, and says why",
                  job.get("sched_id") is None and job.get("liveness_source") == "mtime"
                  and "squeue workdir ambiguous" in out, out)
            check("which, output being fresh, is not called stopped", job.get("status") == "running")
            job, _ = status("lonely")
            check("a job squeue does not mention, fresh, is not declared dead", job.get("status") == "running")
            old = _time.time() - 7200
            os.utime(t.root / "lonely" / "job.out", (old, old))
            reset()
            job, _ = status("lonely")
            check("and, silent past --quiet-after, is stopped", job.get("status") == "stopped", f"{job}")
            job, _ = status("scan_running", ["--root", str(t.root), "--liveness", "process"])
            check("--liveness process ignores squeue", job.get("status") == "stopped" and job.get("sched_id") is None)
            # A machine with SLURM installed can still run ORCA by hand: the
            # local process table answers for what squeue does not mention.
            real = liveness.running_orca_cwds
            liveness.running_orca_cwds = lambda: {(t.root / "lonely").resolve(): _time.time() - 60}
            try:
                reset()
                job, _ = status("lonely")
                check("in slurm mode a local ORCA process still counts, quiet or not",
                      job.get("status") == "running" and job.get("liveness_source") == "process", f"{job}")
            finally:
                liveness.running_orca_cwds = real

            reset()
            os.environ["FAKE_SQUEUE_EXIT"] = "1"
            code, out, err = _orcamon([*slurm, "ls"])
            code2, out2, err2 = _orcamon([*slurm, "show", "scan_running"])
            check("squeue failing does not crash, and says so once",
                  code in (0, 1) and err.count("squeue failed") == 1 and "squeue failed" not in err2, err + err2)
            check("the fallback reads the output's age", "running" in out2.split("\n")[5], out2)
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
            reset()


def test_squeue_lines_parse():
    print("\nsqueue's lines parse, and garbage is skipped")
    from pathlib import Path
    from orcamon.core.liveness import parse_squeue
    rows = parse_squeue("1|RUNNING|/w|2026-09-28T10:00:00\n2|PENDING|/v|N/A\nnope\n3|x|relative|N/A\n")
    check("two good lines, two skipped", [(r[0], r[1], r[2]) for r in rows]
          == [("1", "RUNNING", Path("/w").resolve()), ("2", "PENDING", Path("/v").resolve())], f"{rows}")
    check("a start time is read, N/A is None", rows[0][3] is not None and rows[1][3] is None)


def test_liveness_names_the_job_not_only_the_directory():
    print("\nliveness names the job, not only the directory, and orcamon is not ORCA")
    import json
    import shutil
    import tempfile
    from pathlib import Path
    from orcamon.core.liveness import input_stem, is_orca_process

    check("only the driver and its modules count as ORCA",
          is_orca_process("orca") and is_orca_process("orca_scf_mpi")
          and not is_orca_process("orcamon") and not is_orca_process("orcabox"))
    check("the driver's first input file names the job it runs",
          input_stem(["./orca", "-f", "opt.inp", ""]) == "opt"
          and input_stem(["/x/orca", "/a/b/job.inp"]) == "job"
          and input_stem(["orca_scf_mpi", "job.scfinp", "job"]) is None)

    inp = "! B97-3c Opt\n* xyz 0 1\nH 0 0 0\n*\n"
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        j = root / "j"
        j.mkdir()
        for stem in ("opt", "freq", "old"):
            (j / f"{stem}.inp").write_text(inp)
        for stem in ("opt", "old"):
            (j / f"{stem}.out").write_text("partial output\n")
        # A copy of tail is an `orca` driver that names opt.inp and stays
        # up: the same argument vector `orca opt.inp` has.
        fake = j / "orca"
        shutil.copy("/bin/tail", fake)
        os.chmod(fake, 0o755)
        proc = subprocess.Popen([str(fake), "-f", "opt.inp"], cwd=j,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            _code, doc, err = _json(["--root", str(root), "--liveness", "process",
                                     "--no-cache", "ls"])
            status = {r["label"]: r["status"] for r in (doc or {}).get("jobs", [])}
            check("the driver's own job is running", status.get("j/opt") == "running", f"{status} {err}")
            check("a job the driver does not name is not running",
                  status.get("j/freq") == "not run", f"{status}")
            check("a stopped sibling stays stopped", status.get("j/old") == "stopped", f"{status}")
        finally:
            proc.kill()
            proc.wait()

    with tempfile.TemporaryDirectory() as tmp, tempfile.TemporaryDirectory() as cache:
        d = Path(tmp)
        (d / "job.inp").write_text(inp)
        (d / "job.out").write_text("partial output\n")
        # Run the console script FROM the job directory: its own `comm` is
        # `orcamon`, which is not an ORCA process.
        done = subprocess.run(
            [str(Path(sys.executable).with_name("orcamon")), "--liveness", "process",
             "--no-cache", "ls", "--json"],
            cwd=d, capture_output=True, text=True, timeout=60,
            env={**os.environ, "XDG_CACHE_HOME": cache, "HERDR_ENV": "0"})
        jobs = json.loads(done.stdout)["jobs"] if done.stdout else []
        check("orcamon does not count itself as the job's ORCA process",
              len(jobs) == 1 and jobs[0]["status"] == "stopped",
              f"{done.returncode} {done.stdout[:200]} {done.stderr[:200]}")


# --- the shipped skill -------------------------------------------------------

_SKILL_PLACEHOLDERS = {"JOB": "opt_done", "REGEX": "FINAL", "COND": "done"}
_SHELL_OPERATORS = {">", ">>", "|", "&&", ";", "&", "2>&1"}


def _skill_mentions(text):
    """Every `orcamon <word>` the body uses as a command: inline in
    backticks, or at the start of a line inside a fenced block."""
    import re
    inline = re.findall(r"`orcamon ([a-z][a-z-]*)", text)
    fenced = []
    in_block = False
    for line in text.splitlines():
        if line.startswith("```"):
            in_block = not in_block
            continue
        m = re.match(r"\s*orcamon ([a-z][a-z-]*)", line) if in_block else None
        if m:
            fenced.append(m.group(1))
    return set(inline) | set(fenced)


def _skill_examples(text):
    """(argv, expected exit) for every `orcamon ...` line in a fenced block.
    A trailing `# exit N` states a non-zero exit; any other comment, and
    anything after a shell operator, is not part of the command."""
    import re
    import shlex
    found, in_block = [], False
    for line in text.splitlines():
        if line.startswith("```"):
            in_block = not in_block
            continue
        if not in_block or not line.strip().startswith("orcamon "):
            continue
        command, _, comment = line.partition("#")
        m = re.search(r"exit (\d+)", comment)
        words = shlex.split(command)
        cut = next((i for i, w in enumerate(words) if w in _SHELL_OPERATORS), len(words))
        argv = [_SKILL_PLACEHOLDERS.get(w, w) for w in words[1:cut]]
        found.append((line.strip(), argv, int(m.group(1)) if m else 0))
    return found


def test_the_skill_names_every_command():
    print("\nthe skill names every command, and only commands that exist")
    from orcamon.cli import build_parser, skill, subcommands
    body = skill.body()
    registered = set(subcommands(build_parser()))
    mentioned = _skill_mentions(body)
    check("every registered command appears in the hand-written body",
          registered <= mentioned, f"missing {sorted(registered - mentioned)}")
    check("every `orcamon <word>` in it is a registered command",
          mentioned <= registered, f"unknown {sorted(mentioned - registered)}")


def test_the_skill_examples_run():
    print("\nevery example in the skill runs, with the exit code it states")
    from orcamon.cli import skill
    examples = _skill_examples(skill.body())
    check("the skill has runnable examples", len(examples) >= 5, f"{len(examples)}")
    with _Tree() as t:
        # Only jobs with nothing wrong, so an example's exit states the
        # command's behaviour rather than the tree's contents.
        clean = t.root.parent / "clean"
        clean.mkdir()
        for name in ("opt_done", "scan_profile", "ts"):
            os.rename(t.root / name, clean / name)
        for line, argv, expected in examples:
            try:
                code, out, err = _orcamon(["--root", str(clean), "--liveness", "process", *argv])
            except Exception as exc:  # noqa: BLE001 -- a raise is the failure being tested
                check(f"`{line}` runs", False, repr(exc))
                continue
            check(f"`{line}` exits {expected}", code == expected, f"exit {code}: {(out + err)[:300]}")


def test_the_skill_reference_matches_the_parser():
    print("\nthe generated reference lists every command and every flag")
    import argparse
    from orcamon.cli import build_parser, skill, subcommands
    parser = build_parser()
    reference = skill.reference(parser)
    missing = []
    for name, sub in subcommands(parser).items():
        if f"### orcamon {name}" not in reference:
            missing.append(name)
        for action in sub._actions:
            if isinstance(action, argparse._HelpAction):
                continue
            for flag in action.option_strings:
                if flag not in reference:
                    missing.append(f"{name} {flag}")
    for action in parser._actions:
        for flag in action.option_strings:
            if flag not in ("-h", "--help", "--version") and flag not in reference:
                missing.append(flag)
    check("nothing registered is missing from it", not missing, f"{missing}")
    check("and it ends with the exit codes", "### Exit codes" in reference and "`4`: wait timed out" in reference)


def test_skill_install_and_check():
    print("\nskill install writes a stamped copy; --check says current, stale or absent")
    import tempfile
    from pathlib import Path
    from orcamon import __version__
    with tempfile.TemporaryDirectory() as tmp:
        project = Path(tmp)
        where = ["--project", str(project)]
        path = project / ".claude" / "skills" / "orcamon" / "SKILL.md"
        code, out, err = _orcamon(["skill", "--check", *where])
        check("absent: --check exits 2", code == 2, out + err)
        code, out, err = _orcamon(["skill", "install", *where])
        check("install writes it and prints the path", code == 0 and out.strip() == str(path) and path.exists(),
              out + err)
        text = path.read_text()
        check("with the version in its frontmatter",
              text.startswith("---\nname: orcamon\n") and f"orcamon_version: {__version__}" in text.split("---")[1])
        code, out, err = _orcamon(["skill"])
        check("and otherwise exactly what `orcamon skill` prints",
              text.replace(f"metadata:\n  orcamon_version: {__version__}\n", "") == out, "")
        code, out, err = _orcamon(["skill", "--check", *where])
        check("current: --check exits 0", code == 0, out + err)
        path.write_text(text.replace("## Waiting", "## Waiting (edited)"))
        code, out, err = _orcamon(["skill", "--check", *where])
        check("edited: --check exits 1", code == 1 and "differs" in out, out)
        path.write_text(text.replace(f"orcamon_version: {__version__}", "orcamon_version: 0.0.1"))
        code, out, err = _orcamon(["skill", "--check", *where])
        check("stale: --check exits 1 naming both versions",
              code == 1 and f"installed skill is 0.0.1, orcamon is {__version__}" in out, out)
        path.write_text("---\nname: my-own\n---\nmine\n")
        code, out, err = _orcamon(["skill", "install", *where])
        check("someone else's SKILL.md is not overwritten", code == 2 and path.read_text().endswith("mine\n"), err)
        code, out, err = _orcamon(["skill", "install", "--force", *where])
        check("unless --force", code == 0 and "orcamon_version" in path.read_text())


def test_the_skill_is_general():
    print("\nthe shipped skill says nothing about this repository")
    from orcamon.cli import build_parser, skill
    text = skill.full_text(build_parser())
    deny = ("masterThesis", "computational/C", "C8_", "orca_io", "homelab")
    found = [word for word in deny if word in text]
    check("none of this repository's names or paths", not found, f"{found}")


def test_this_repository_has_the_current_skill():
    print("\nthis repository's installed copy is the one this orcamon prints")
    code, out, err = _orcamon(["skill", "--check", "--project", HERE])
    check("orcamon skill --check --project <repo> (run `orcamon skill install --project .` after an upgrade)",
          code == 0, out + err)

if __name__ == "__main__":
    test_the_core_needs_only_the_standard_library()
    test_a_job_argument_names_one_job()
    test_only_orca_inputs_are_jobs()
    test_ls_lists_every_job_boundedly()
    test_ls_prints_labels_an_agent_can_pass_back()
    test_ls_says_where_an_irc_is()
    test_show_and_ls_agree()
    test_geom_never_substitutes_a_geometry()
    test_snapshot_writes_a_png()
    test_a_path_is_listed_and_drawn()
    test_an_neb_is_listed_and_drawn()
    test_any_mode_is_drawn_and_described()
    test_measure_names_the_geometry_it_measured()
    test_the_bond_cache_holds_the_geometry_it_answers_for()
    test_wait_returns_when_the_job_is_done()
    test_the_other_commands_answer_boundedly()
    test_a_cached_state_equals_a_fresh_parse()
    test_the_cli_uses_the_cache_and_can_refuse_it()
    test_events_announce_changes_not_the_first_scan()
    test_notifications_reach_the_terminal_through_tmux()
    test_the_geometry_pane_degrades_to_text()
    test_the_images_extra_needs_no_matplotlib()
    test_text_fog_dims_only_the_far_third()
    test_text_fog_holds_still_under_rotation()
    test_the_tui_runs_headless()
    test_the_tui_scrubs_a_path()
    test_the_tui_cycles_real_modes()
    test_the_tui_selects_atoms_by_number()
    test_the_tui_tells_the_agent_what_it_shows()
    test_view_names_the_geometry_its_numbers_belong_to()
    test_the_tui_publishes_the_numbering_it_shows()
    test_focus_waits_for_the_tui_to_take_it()
    test_the_tui_takes_the_agents_focus()
    test_b_does_not_leave_the_agents_mode_running()
    test_the_tui_picks_atoms_by_click()
    test_edit_writes_a_checked_geometry()
    test_check_finds_what_would_fail()
    test_connectivity_names_bonds_made_and_broken()
    test_slurm_names_what_the_login_node_cannot_see()
    test_squeue_lines_parse()
    test_liveness_names_the_job_not_only_the_directory()
    test_the_skill_names_every_command()
    test_the_skill_examples_run()
    test_the_skill_reference_matches_the_parser()
    test_skill_install_and_check()
    test_the_skill_is_general()
    test_this_repository_has_the_current_skill()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
