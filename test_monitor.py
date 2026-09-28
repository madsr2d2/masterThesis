"""
The contract orcamon (`tools/orcamon/`) -- the terminal ORCA viewer -- holds
ORCA output to.

    python test_monitor.py

The monitor is a GENERAL viewer, not a reader for this project's jobs, so
every case here is a few lines of ORCA's own output fed through the parser,
not a job in the tree. Each is a way the viewer showed the wrong thing:

- the ENERGY of a multilayer job was the high-level region alone. QM/XTB
  prints four `FINAL SINGLE POINT ENERGY` lines per geometry and the pattern
  could only match the unlabelled one, which is the QM1 region's;
- a relaxed SCAN was keyed by cycle, which ORCA restarts at every step, so
  the chart's x-axis ran 1..9, 1..16, ... and scrubbing to "cycle 3" of a late
  step drew step 1's geometry;
- the GEOMETRY STEPS pane rendered 41 lines into 8 and so showed only the
  oldest of the cycles it listed, never the one running;
- a point whose coordinates ORCA never printed silently borrowed another
  geometry.

`run_gates.py` discovers this file by its name, like any other gate.
"""
import os
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from orcamon import validate  # noqa: E402
from orcamon.core.orca_input import describe_spin, parse_input  # noqa: E402
from orcamon.core.parser import JobState  # noqa: E402
from orcamon.core.report import (  # noqa: E402
    REPORT_SCHEMA, STEP_ROWS, build_report, geometry_shown, render_markup, render_plain,
    render_steps_markup, render_steps_plain, steps_rows,
)
from orcamon.core.job import Job  # noqa: E402
from orcamon.core.liveness import Liveness  # noqa: E402
from orcamon.core.status import QUIET_AFTER_S, Status, attention, compute_status  # noqa: E402

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)


def _feed(state, text):
    for line in text.strip("\n").split("\n"):
        state.feed_line(line)


def _geometry(state, n_atoms, x=0.0):
    _feed(state, "CARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------")
    for i in range(n_atoms):
        state.feed_line(f"  H      {x + i:.6f}    0.000000    0.000000")
    state.feed_line("")


def _cycle(state, n):
    state.feed_line(f"         *                GEOMETRY OPTIMIZATION CYCLE  {n:>2}            *")


class _Job:
    def __init__(self, state):
        self.state = state


def _report_job(state, inp_text, liveness, now=None):
    """A Job around a fed state and an input text, with status and flags
    derived exactly as `Job.refresh` derives them -- without the disk."""
    state.has_out = True
    job = Job(state.path, state.stem, Path("/"), label="x")
    job.state, job.liveness, job.parsed = state, liveness, True
    job.input = parse_input(inp_text) if inp_text is not None else None
    job.status = compute_status(state, liveness, now)
    job.flags = attention(state, job.input, job.status, now)
    return job


def test_a_multilayer_energy_is_the_total():
    print("\na multilayer energy is the total, not the high-level region")
    state = JobState(path=Path("/nonexistent"))
    _cycle(state, 1)
    _geometry(state, 5)
    _feed(state, """
FINAL SINGLE POINT ENERGY (L-QM2)     -266.681990679470
FINAL SINGLE POINT ENERGY (S-QM2)      -40.834442837800
FINAL SINGLE POINT ENERGY      -647.909099472771
FINAL SINGLE POINT ENERGY (QM/QM2)     -873.756647314441
""")
    check("the job's energy is the QM/QM2 total",
          state.final_energy == -873.756647314441, f"{state.final_energy}")
    check("and says which it is", state.final_energy_label == "QM/QM2")
    check("the plotted point carries the same", state.points[-1].energy == -873.756647314441)

    plain = JobState(path=Path("/nonexistent"))
    plain.feed_line("FINAL SINGLE POINT ENERGY      -75.959340")
    check("an ordinary job's unlabelled energy is still read",
          plain.final_energy == -75.959340 and plain.final_energy_label is None)


def test_a_scan_is_keyed_by_step_and_cycle():
    print("\na scan is keyed by (step, cycle), not by the cycle ORCA restarts")
    state = JobState(path=Path("/nonexistent"))
    _feed(state, """
There is 1 parameter to be scanned.
There will be   2 constrained geometry optimizations.
""")
    for step, value, cycles in ((1, 1.40, 2), (2, 1.50, 1)):
        _feed(state, f"""
         *************************************************************
         *               RELAXED SURFACE SCAN STEP   {step}               *
         *                                                           *
         *                 Bond (130, 128)  :   {value:.8f}           *
         *************************************************************
""")
        for c in range(1, cycles + 1):
            _cycle(state, c)
            _geometry(state, 3, x=10 * step + c)
            state.feed_line(f"FINAL SINGLE POINT ENERGY     -{100 + step + c / 10:.6f}")

    keys = [p.key for p in state.points]
    check("three distinct geometries, not two", keys == [(1, 1), (1, 2), (2, 1)], f"{keys}")
    check("one profile point per step", sorted(state.scan_points) == [1, 2])
    check("the profile point is the step's LAST geometry",
          state.scan_points[1].key == (1, 2))
    check("the scanned coordinate is read per step",
          state.scan_values == {1: 1.40, 2: 1.50} and state.scan_label == "Bond (130, 128)",
          f"{state.scan_label} {state.scan_values}")
    check("and the total step count", state.scan_total == 2)
    step2 = state.points[-1]
    atoms, shown = geometry_shown(_Job(state), step2)
    check("step 2's cycle 1 draws step 2's geometry, not step 1's",
          shown is step2 and atoms[0][1] == 21.0, f"x0 = {atoms[0][1]}")


def test_a_missing_geometry_is_labelled_not_borrowed_silently():
    print("\na geometry that was never printed is labelled")
    state = JobState(path=Path("/nonexistent"))
    _cycle(state, 1)
    _geometry(state, 3)
    state.feed_line("FINAL SINGLE POINT ENERGY     -100.0")
    _cycle(state, 2)
    state.feed_line("FINAL SINGLE POINT ENERGY     -100.1")
    first, second = state.points
    atoms, shown = geometry_shown(_Job(state), second)
    check("the point exists, with its energy", second.energy == -100.1 and not second.atoms)
    check("the pane is handed the geometry it really draws", shown is first and atoms is first.atoms)


def test_every_imaginary_frequency_is_kept():
    print("\nevery imaginary frequency is kept, with its value")
    state = JobState(path=Path("/nonexistent"))
    _feed(state, """
VIBRATIONAL FREQUENCIES
-----------------------

Scaling factor for frequencies =  1.000000000  (already applied!)

     0:       0.00 cm**-1
     6:    -229.18 cm**-1  ***imaginary mode***
     7:     -10.30 cm**-1  ***imaginary mode***
     8:       8.88 cm**-1

------------
""")
    check("both modes, the small one included", state.imaginary_freqs == [-229.18, -10.30],
          f"{state.imaginary_freqs}")
    check("and the count agrees", state.n_imaginary == 2)
    check("every mode is kept too, zero and real ones included, in ORCA's order",
          state.frequencies == [0.0, -229.18, -10.30, 8.88], f"{state.frequencies}")


def test_the_steps_pane_shows_the_latest_cycle():
    print("\nthe geometry-steps pane shows the latest cycle")
    state = JobState(path=Path("/nonexistent"))
    for c in range(1, 11):
        _cycle(state, c)
        _feed(state, f"""
          Energy change      -0.0000038757            0.0000050000      NO
          RMS gradient        {0.01 / c:.10f}            0.0001000000      NO
          MAX gradient        0.0025000000            0.0003000000      NO
          RMS step            0.0020825435            0.0020000000      NO
          MAX step            0.0192646101            0.0040000000      YES
""")
    lines = render_steps_markup(steps_rows(state)).split("\n")
    check("one row per cycle, header and tolerances above",
          len(lines) == STEP_ROWS + 2, f"{len(lines)} lines")
    check("and the last row is the newest cycle", lines[-1].lstrip().startswith("10 "),
          lines[-1])


def test_the_input_names_charge_and_multiplicity():
    print("\nthe input names charge and multiplicity, in each syntax")
    job = parse_input("""
# a comment with * xyz 9 9 in it
! B3LYP def2-SVP OptTS Freq PAL8
* xyz -1 2
H 0 0 0
*
""")
    check("the * xyz line", (job.charge, job.mult) == (-1, 2), f"{job.charge} {job.mult}")
    check("run types split from the method",
          job.run_types == ["OptTS", "Freq"] and job.method == "B3LYP def2-SVP",
          f"{job.run_types} / {job.method}")
    check("PAL8 is the processor count", job.nprocs == 8)

    job = parse_input("! r2SCAN-3c Opt\n*xyzfile 0 3 geom.xyz\n")
    check("*xyzfile with no space", (job.charge, job.mult) == (0, 3))

    job = parse_input("! HF\n%coords\n  CTyp xyz\n  Charge 1\n  Mult 2\n  coords\n  end\nend\n")
    check("the %coords block", (job.charge, job.mult) == (1, 2))

    job = parse_input("""
! QM/XTB r2SCAN-3c Opt
%qmmm
  QMAtoms {0:5} end
  Charge_Total -1
  Mult_Total 1
end
%geom
  Scan
    B 0 1 = 1.0, 2.0, 5
  end
end
* xyz 0 1
*
""")
    check("a multilayer job carries the whole system's charge and multiplicity too",
          job.multilayer and job.layers == {"total": (-1, 1)} and (job.charge, job.mult) == (0, 1),
          f"{job.layers}")
    check("and a %geom Scan is a run type", "Scan" in job.run_types)
    check("the spin is named", describe_spin(-1, 2) == "charge -1 · doublet")


def test_the_basis_is_the_one_orca_reports():
    print("\nthe basis is the orbital basis ORCA reports, a composite method's included")
    state = JobState(path=Path("/nonexistent"))
    _feed(state, """
----- Orbital basis set information -----
Your calculation utilizes the basis: def2-mTZVP
   Stefan Grimme, ...

----- AuxJ basis set information -----
Your calculation utilizes the auxiliary basis: def2-mTZVP/J
""")
    check("the orbital basis, and not the RI auxiliary one",
          state.basis == "def2-mTZVP", f"{state.basis}")

    job = _report_job(state, "! B97-3c Opt\n* xyz 0 1\n*\n", Liveness(True, "process"))
    text = render_markup(build_report(job))
    check("and the summary names it alone",
          "basis  def2-mTZVP\n" in text and "/J" not in text, text)


SYNTHETIC_OUTPUT = """\
----- Orbital basis set information -----
Your calculation utilizes the basis: def2-SVP
QM1 Subsystem      ...  0 1
                         2
Composition done
   1   -75.9000000000   0.000000e+00  0.1
   2   -75.9500000000  -5.000000e-02  0.01
There is 1 parameter to be scanned.
There will be   2 constrained geometry optimizations.
         *************************************************************
         *               RELAXED SURFACE SCAN STEP   1               *
         *                 Bond (0, 1)  :   1.00000000           *
         *************************************************************
some ordinary line that mentions a gradient and a step
         *                GEOMETRY OPTIMIZATION CYCLE   1            *
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.000000    0.000000    0.000000
  H      0.960000    0.000000    0.000000
  H     -0.240000    0.930000    0.000000

FINAL SINGLE POINT ENERGY (L-QM2)     -1.000000
FINAL SINGLE POINT ENERGY      -75.960000
FINAL SINGLE POINT ENERGY (QM/QM2)     -80.000000
          Energy change      -0.0000038757            0.0000050000      NO
          RMS gradient        0.0004500000            0.0001000000      NO
          MAX gradient        0.0025000000            0.0003000000      NO
          RMS step            0.0020825435            0.0020000000      NO
          MAX step            0.0192646101            0.0040000000      YES
        Hessian has     1 negative eigenvalue
   3   -80.1000000000   0.000000e+00  0.1
VIBRATIONAL FREQUENCIES
-----------------------

Scaling factor for frequencies =  1.000000000  (already applied!)

     0:       0.00 cm**-1
     6:    -229.18 cm**-1  ***imaginary mode***
     7:     155.00 cm**-1

------------
                 ****ORCA TERMINATED NORMALLY****
"""


def _state_fields(state):
    import dataclasses
    from collections import deque

    def norm(v):
        if dataclasses.is_dataclass(v):
            return tuple((f.name, norm(getattr(v, f.name))) for f in dataclasses.fields(v))
        if isinstance(v, (list, tuple, deque)):
            return tuple(norm(x) for x in v)
        if isinstance(v, dict):
            return tuple(sorted((k, norm(x)) for k, x in v.items()))
        if isinstance(v, set):
            return tuple(sorted(v))
        return v
    return {f.name: norm(getattr(state, f.name)) for f in dataclasses.fields(state)}


def test_the_fast_path_reads_what_the_line_path_reads():
    print("\nfeed_text and chunked reads leave the state feed_line would")
    import tempfile
    from orcamon.core.parser import read_appended

    reference = JobState(path=Path("/nonexistent"))
    for line in SYNTHETIC_OUTPUT.split("\n")[:-1]:
        reference.feed_line(line)
    expected = _state_fields(reference)
    check("the reference parse read what it should",
          reference.basis == "def2-SVP" and reference.qm_atom_indices == {0, 1, 2}
          and reference.final_energy == -80.0 and reference.imaginary_freqs == [-229.18]
          and len(reference.convergence_history) == 1 and len(reference.scf_iterations) == 2
          and reference.scan_values == {1: 1.0} and reference.normal_completion,
          f"{reference.basis} {reference.qm_atom_indices} {reference.final_energy} "
          f"{reference.imaginary_freqs} {len(reference.scf_iterations)}")

    whole = JobState(path=Path("/nonexistent"))
    whole.feed_text(SYNTHETIC_OUTPUT)
    differ = [k for k, v in _state_fields(whole).items() if v != expected[k]]
    check("one feed_text of the whole output", not differ, f"{differ}")

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "job.out"
        out.write_text(SYNTHETIC_OUTPUT)
        for chunk in (1, 7, 64, 1 << 20):
            state = JobState(path=Path("/nonexistent"))
            read_appended(state, out, chunk_bytes=chunk)
            got = _state_fields(state)
            differ = [k for k, v in got.items() if k != "offset" and v != expected[k]]
            check(f"read in {chunk}-byte chunks", not differ and state.offset == out.stat().st_size,
                  f"{differ}")


def test_a_frame_is_a_small_faithful_palette_png():
    print("\na frame goes out as a palette PNG, rendered outside pyplot")
    import io
    import matplotlib.pyplot as plt
    import numpy as np
    from PIL import Image
    from orcamon.tui import geometry_render

    water = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0)]
    before = plt.get_fignums()
    image = geometry_render.render(water, size_px=(300, 240))
    check("the render is RGB at the size asked for",
          image.mode == "RGB" and image.size == (300, 240), f"{image.mode} {image.size}")
    check("and left no figure in pyplot's global registry",
          plt.get_fignums() == before, f"{plt.get_fignums()}")

    sent = Image.open(io.BytesIO(geometry_render.frame_png(image)))
    error = np.abs(np.asarray(sent.convert("RGB"), int) - np.asarray(image, int)).mean()
    check("the frame is a 256-colour palette PNG", sent.mode == "P" and sent.size == image.size,
          f"{sent.mode} {sent.size}")
    check("within ~1 level in 255 of the render", error < 2.0, f"mean error {error:.2f}")


def test_every_marker_reaches_its_parser():
    print("\nevery marker line reaches its parser past the prefilter")
    check("validate.MARKER_CASES all read", validate.check_markers() == 0)



_TERMINATED = "                 ****ORCA TERMINATED NORMALLY****"
_MAXITER = ("       The optimization did not converge but reached the maximum \n"
            "       number of optimization cycles.")
_OLD = QUIET_AFTER_S + 60


def _freq_block(freqs):
    lines = ["VIBRATIONAL FREQUENCIES", "-----------------------", "",
             "Scaling factor for frequencies =  1.000000000  (already applied!)", ""]
    for i, f in enumerate(freqs):
        tag = "  ***imaginary mode***" if f < 0 else ""
        lines.append(f"   {i:>3}:   {f:>9.2f} cm**-1{tag}")
    return "\n".join(lines + ["", "------------"])


def _fed(text, mtime_age=0.0, now=1_000_000.0):
    state = JobState(path=Path("/nonexistent"))
    _feed(state, text)
    state.has_out, state.mtime = True, now - mtime_age
    return state


def test_status_names_the_outcome():
    print("\nstatus says where the job is; it does not call a bad result a success")
    now = 1_000_000.0
    alive, gone, unknown = Liveness(True, "process"), Liveness(False, "process"), Liveness(None, "mtime")

    state = _fed(f"{_MAXITER}\n{_freq_block([-40.0, 10.0])}\n{_TERMINATED}")
    job = _report_job(state, "! B97-3c Opt Freq\n* xyz 0 1\n*\n", gone, now)
    check("normal termination after MaxIter is FINISHED, not converged",
          job.status is Status.FINISHED, job.status.value)
    check("and carries opt_not_converged", "opt_not_converged" in [f.code for f in job.flags],
          f"{[f.code for f in job.flags]}")

    state = _fed("some output\n")
    check("a process gone without a termination line is STOPPED, not FAILED",
          compute_status(state, gone, now) is Status.STOPPED)
    state = _fed("ORCA finished by error termination in SCF\n")
    check("an ORCA error marker is FAILED", compute_status(state, gone, now) is Status.FAILED)
    check("even while a process is still there", compute_status(state, alive, now) is Status.FAILED)
    state = _fed("some output\n", mtime_age=_OLD)
    check("unknown liveness and old output is QUIET", compute_status(state, unknown, now) is Status.QUIET)
    state = _fed("some output\n", mtime_age=10)
    check("unknown liveness and fresh output is RUNNING", compute_status(state, unknown, now) is Status.RUNNING)
    check("unknown liveness never gives STOPPED",
          compute_status(_fed("x\n", mtime_age=10 * _OLD), unknown, now) is not Status.STOPPED)

    empty = JobState(path=Path("/nonexistent"))
    check("no output and no process is NOT_RUN", compute_status(empty, gone, now) is Status.NOT_RUN)
    check("no output and the scheduler has it pending is QUEUED",
          compute_status(empty, Liveness(False, "slurm", queued=True), now) is Status.QUEUED)
    check("no output yet but a live process is RUNNING", compute_status(empty, alive, now) is Status.RUNNING)


def _plateau(state, cycles=20, rms=4.5e-4):
    for c in range(1, cycles + 1):
        _cycle(state, c)
        _feed(state, f"""
          Energy change      -0.0000038757            0.0000050000      NO
          RMS gradient        {rms if rms else 0.01 / c:.10f}            0.0001000000      NO
          MAX gradient        0.0025000000            0.0003000000      NO
          RMS step            0.0020825435            0.0020000000      NO
          MAX step            0.0192646101            0.0040000000      NO
""")


def _codes(state, inp, liveness, now=1_000_000.0):
    return [f.code for f in _report_job(state, inp, liveness, now).flags]


def test_attention_flags():
    print("\nevery attention flag fires on its case and only on it")
    now = 1_000_000.0
    alive, gone, unknown = Liveness(True, "process"), Liveness(False, "process"), Liveness(None, "mtime")
    opt, optts = "! B97-3c Opt Freq\n* xyz 0 1\n*\n", "! B97-3c OptTS Freq\n* xyz 0 1\n*\n"

    codes = _codes(_fed("ORCA finished by error termination in SCF\n"), opt, gone)
    check("failed", codes == ["failed"], f"{codes}")
    flag = _report_job(_fed("ORCA finished by error termination in SCF\n"), opt, gone).flags[0]
    check("and it quotes the error line", "error termination in SCF" in flag.message, flag.message)
    check("stopped", _codes(_fed("x\n"), opt, gone) == ["stopped"])
    check("not stopped while alive", _codes(_fed("x\n"), opt, alive) == [])

    state = JobState(path=Path("/nonexistent"))
    _plateau(state)
    check("stalled", _codes(state, opt, alive) == ["stalled"], f"{_codes(state, opt, alive)}")
    state = JobState(path=Path("/nonexistent"))
    _plateau(state, rms=None)
    check("an improving gradient is not stalled", _codes(state, opt, alive) == [])

    check("quiet", _codes(_fed("x\n", mtime_age=_OLD), opt, unknown) == ["quiet"])
    check("not quiet while fresh", _codes(_fed("x\n", mtime_age=5), opt, unknown) == [])

    state = _fed(f"{_MAXITER}\n{_TERMINATED}")
    check("opt_not_converged", _codes(state, opt, gone) == ["opt_not_converged"], f"{_codes(state, opt, gone)}")
    state = _fed(f"      ***        THE OPTIMIZATION HAS CONVERGED     ***\n{_TERMINATED}")
    check("a converged optimization has none", _codes(state, opt, gone) == [])
    state = _fed(f"""
         *               RELAXED SURFACE SCAN STEP   4               *
         *************************************************************
{_MAXITER}
         *               RELAXED SURFACE SCAN STEP   5               *
         *************************************************************
      ***        THE OPTIMIZATION HAS CONVERGED     ***
{_TERMINATED}""")
    job = _report_job(state, opt, gone)
    check("a scan step that hit MaxIter is flagged though later steps converged",
          [f.code for f in job.flags] == ["opt_not_converged"] and "step 4" in job.flags[0].message,
          f"{[(f.code, f.message) for f in job.flags]}")

    running_ts = "        Hessian has     2 negative eigenvalues\n"
    check("ts_hessian", _codes(_fed(running_ts), optts, alive) == ["ts_hessian"])
    check("a TS Hessian with one negative eigenvalue has none",
          _codes(_fed("        Hessian has     1 negative eigenvalue\n"), optts, alive) == [])
    check("a minimum search's Hessian is not a TS's business", _codes(_fed(running_ts), opt, alive) == [])
    ended = _fed(f"{running_ts}{_MAXITER}\n{_TERMINATED}")
    check("a TS search that ended at MaxIter with no final frequencies is judged by its Hessian",
          _codes(ended, optts, gone) == ["opt_not_converged", "ts_hessian"], f"{_codes(ended, optts, gone)}")
    judged = _fed(f"{running_ts}      ***        THE OPTIMIZATION HAS CONVERGED     ***\n"
                  f"{_freq_block([0.0, -300.0, 50.0])}\n{_TERMINATED}")
    check("but a final frequency block supersedes the Hessian", _codes(judged, optts, gone) == [])

    two = _fed(f"{_freq_block([0.0, -300.0, -20.0, 50.0])}\n{_TERMINATED}")
    one = _fed(f"{_freq_block([0.0, -300.0, 50.0])}\n{_TERMINATED}")
    none = _fed(f"{_freq_block([0.0, 50.0])}\n{_TERMINATED}")
    check("ts_imaginary", _codes(two, optts, gone) == ["ts_imaginary"], f"{_codes(two, optts, gone)}")
    check("a TS with exactly one imaginary mode has none", _codes(one, optts, gone) == [])
    check("a TS with none is flagged", _codes(none, optts, gone) == ["ts_imaginary"])
    check("minimum_imaginary", _codes(one, opt, gone) == ["minimum_imaginary"])
    check("a minimum with none has none", _codes(none, opt, gone) == [])
    check("a small imaginary mode is not excused",
          _codes(_fed(f"{_freq_block([-5.0, 50.0])}\n{_TERMINATED}"), opt, gone) == ["minimum_imaginary"])
    check("without a readable input, no run-type flags", _codes(two, None, gone) == [])
    state = JobState(path=Path("/nonexistent"))
    _cycle(state, 1)
    _feed(state, _freq_block([-300.0, -40.0, -20.0, 50.0]))
    _cycle(state, 2)
    state.has_out = True
    check("an intermediate Hessian's frequencies are not flagged", _codes(state, optts, gone) == ["stopped"],
          f"{_codes(state, optts, gone)}")
    check("and are labelled with their cycle", state.freq_cycle == 1 and not state.freqs_final
          and "Hessian at cycle 1" in render_plain(build_report(_report_job(state, optts, gone))))

    check("qm2_errors", _codes(_fed(f"There was an error in the QM2 calculation\n{_TERMINATED}"), opt, gone)
          == ["qm2_errors"])
    scan_head = "There will be   3 constrained geometry optimizations.\n"
    step = "         *               RELAXED SURFACE SCAN STEP   {}               *\n"
    short = _fed(scan_head + step.format(1) + _TERMINATED)
    whole = _fed(scan_head + step.format(1) + step.format(2) + step.format(3) + _TERMINATED)
    check("scan_incomplete: a scan that terminated normally short of its last step",
          _codes(short, opt, gone) == ["scan_incomplete"], f"{_codes(short, opt, gone)}")
    check("a scan that reached its last step has none", _codes(whole, opt, gone) == [])
    check("a clean finished optimization has no flags at all",
          _codes(_fed(f"      ***        THE OPTIMIZATION HAS CONVERGED     ***\n{_freq_block([0.0, 40.0])}\n{_TERMINATED}"),
                 opt, gone) == [])


_MULTILAYER_INPUT = """
! QM/XTB r2SCAN-3c OptTS Freq PAL8
%maxcore 4000
%qmmm
  QMAtoms {0:2} end
  Charge_Total -1
  Mult_Total 1
end
* xyz 0 1
*
"""


def test_plain_and_markup_say_the_same_thing():
    print("\nthe plain and the markup renderings say the same thing")
    from rich.text import Text
    state = JobState(path=Path("/nonexistent"))
    state.feed_text(SYNTHETIC_OUTPUT.replace("                 ****ORCA TERMINATED NORMALLY****\n", ""))
    _feed(state, "There was an error in the QM2 calculation\nAborting the run [because of x]")
    job = _report_job(state, _MULTILAYER_INPUT, Liveness(False, "process"))
    job.label = "rxn/[odd] label/ts"
    report = build_report(job, now=1_000_000.0)
    plain = render_plain(report).split("\n")
    markup = render_markup(report).split("\n")
    stripped = [Text.from_markup(line).plain for line in markup]
    check("line for line, once the markup is stripped", stripped == plain,
          f"{[(a, b) for a, b in zip(stripped, plain) if a != b][:3]}")
    check("with the attention flags in both", any("QM2 error" in line for line in plain)
          and any(line.startswith("! ") for line in plain))
    check("and a literal '[' survives the markup", "rxn/[odd] label/ts" in stripped[0], stripped[0])

    rows = steps_rows(state)
    steps_plain = render_steps_plain(rows).split("\n")
    steps_markup = [Text.from_markup(line).plain.rstrip() for line in render_steps_markup(rows).split("\n")]
    check("and the steps table too", steps_markup == steps_plain, f"{steps_markup} / {steps_plain}")
    check("which marks the met criterion in the text, not only in colour",
          any("*" in line for line in steps_plain[2:]))


def test_report_keys_are_stable():
    print("\nthe report's JSON keys are a contract")
    keys = set(build_report(_report_job(JobState(path=Path("/nonexistent")), None,
                                        Liveness(False, "process"))).to_dict())
    expected = {
        "schema", "label", "path", "stem",
        "run_types", "method", "basis", "charge", "mult", "layers", "multilayer",
        "n_atoms", "n_qm_atoms", "nprocs", "maxcore_mb",
        "status", "liveness_source", "sched_id", "cycle", "max_cycles",
        "scan_step", "scan_total", "scan_coordinate", "wall_time_s", "last_output_age_s",
        "opt_converged", "opt_maxiter_reached",
        "criteria", "negative_eigenvalues", "imaginary_freqs", "freq_cycle", "energy_eh", "energy_label",
        "delta_kj_mol", "scan_max_kj_mol", "scan_max_at", "qm2_errors", "crash_lines", "attention",
    }
    check(f"the key set of schema {REPORT_SCHEMA} is unchanged (change it deliberately: "
          "bump REPORT_SCHEMA and edit this set)", keys == expected,
          f"added {sorted(keys - expected)}, removed {sorted(expected - keys)}")
    check("and the schema is in the document", REPORT_SCHEMA == 1)


if __name__ == "__main__":
    test_a_multilayer_energy_is_the_total()
    test_a_scan_is_keyed_by_step_and_cycle()
    test_a_missing_geometry_is_labelled_not_borrowed_silently()
    test_every_imaginary_frequency_is_kept()
    test_the_steps_pane_shows_the_latest_cycle()
    test_the_input_names_charge_and_multiplicity()
    test_the_basis_is_the_one_orca_reports()
    test_the_fast_path_reads_what_the_line_path_reads()
    test_a_frame_is_a_small_faithful_palette_png()
    test_every_marker_reaches_its_parser()
    test_status_names_the_outcome()
    test_attention_flags()
    test_plain_and_markup_say_the_same_thing()
    test_report_keys_are_stable()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
