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
# The TUI tests run headless. Inheriting a herdr pane would draw the
# molecule over the pane this runs in.
os.environ["HERDR_ENV"] = "0"

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
    test_ls_lists_every_job_boundedly()
    test_ls_prints_labels_an_agent_can_pass_back()
    test_show_and_ls_agree()
    test_geom_never_substitutes_a_geometry()
    test_snapshot_writes_a_png()
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
