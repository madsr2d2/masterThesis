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
from array import array
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from orcamon import validate  # noqa: E402
from orcamon.core.orca_input import describe_spin, parse_input  # noqa: E402
from orcamon.core.parser import IrcRow, JobState, NebRow  # noqa: E402
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


def test_normal_modes_are_parsed():
    print("\nan imaginary mode's printed displacement vector is read from NORMAL MODES")
    state = JobState(path=Path("/nonexistent"))
    state.feed_text("""\
VIBRATIONAL FREQUENCIES
-----------------------

Scaling factor for frequencies =  1.000000000  (already applied!)

    0:   -100.00 cm**-1  ***imaginary mode***
    1:     10.00 cm**-1
    2:     20.00 cm**-1
    3:     30.00 cm**-1
    4:     40.00 cm**-1
    5:     50.00 cm**-1

NORMAL MODES
------------

These modes are the Cartesian displacements weighted by the diagonal matrix
M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
Thus, these vectors are normalized but *not* orthogonal

                  0          1          2          3          4          5
      0       0.100000   0.000000   0.000000   0.000000   0.000000   0.000000
      1       0.200000   0.000000   0.000000   0.000000   0.000000   0.000000
      2      -0.300000   0.000000   0.000000   0.000000   0.000000   0.000000
      3       0.400000   0.000000   0.000000   0.000000   0.000000   0.000000
      4      -0.500000   0.000000   0.000000   0.000000   0.000000   0.000000
      5       0.600000   0.000000   0.000000   0.000000   0.000000   0.000000

""")
    check("one imaginary mode is stored",
          state.modes is not None and len(state.modes) == 1, f"{state.modes}")
    check("with ORCA's own mode number and its signed frequency",
          state.modes is not None and state.modes[0].index == 0 and state.modes[0].cm1 == -100.0,
          f"{state.modes}")
    check("and the 3N printed values, coordinate-major",
          state.modes is not None
          and state.modes[0].vector == [0.1, 0.2, -0.3, 0.4, -0.5, 0.6],
          f"{state.modes[0].vector if state.modes else None}")


def test_normal_modes_keep_only_imaginary():
    print("\nonly the imaginary modes' columns of a NORMAL MODES block are kept")
    state = JobState(path=Path("/nonexistent"))
    state.feed_text("""\
VIBRATIONAL FREQUENCIES
-----------------------

Scaling factor for frequencies =  1.000000000  (already applied!)

    0:   -100.00 cm**-1  ***imaginary mode***
    1:    -50.00 cm**-1  ***imaginary mode***
    2:     20.00 cm**-1
    3:     30.00 cm**-1
    4:     40.00 cm**-1
    5:     50.00 cm**-1

NORMAL MODES
------------

These modes are the Cartesian displacements weighted by the diagonal matrix
M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
Thus, these vectors are normalized but *not* orthogonal

                  0          1          2          3          4          5
      0       0.100000   0.000000   0.000000   0.000000   0.000000   0.000000
      1       0.200000   0.000000   0.000000   0.000000   0.000000   0.000000
      2      -0.300000   0.000000   0.000000   0.000000   0.000000   0.000000
      3       0.400000   0.000000   0.000000   0.000000   0.000000   0.000000
      4      -0.500000   0.000000   0.000000   0.000000   0.000000   0.000000
      5       0.600000   0.000000   0.000000   0.000000   0.000000   0.000000

""")
    check("both imaginary modes, and neither real one",
          state.modes is not None and [m.index for m in state.modes] == [0, 1], f"{state.modes}")


def test_normal_modes_take_the_last_block():
    print("\na later Hessian's NORMAL MODES block replaces the earlier one's")
    state = JobState(path=Path("/nonexistent"))
    state.feed_text("""\
VIBRATIONAL FREQUENCIES
-----------------------

Scaling factor for frequencies =  1.000000000  (already applied!)

    0:   -100.00 cm**-1  ***imaginary mode***
    1:     10.00 cm**-1
    2:     20.00 cm**-1
    3:     30.00 cm**-1
    4:     40.00 cm**-1
    5:     50.00 cm**-1

NORMAL MODES
------------

These modes are the Cartesian displacements weighted by the diagonal matrix
M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
Thus, these vectors are normalized but *not* orthogonal

                  0          1          2          3          4          5
      0       1.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      1       1.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      2       1.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      3       1.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      4       1.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      5       1.000000   0.000000   0.000000   0.000000   0.000000   0.000000

VIBRATIONAL FREQUENCIES
-----------------------

Scaling factor for frequencies =  1.000000000  (already applied!)

    0:     10.00 cm**-1
    1:     20.00 cm**-1
    2:   -200.00 cm**-1  ***imaginary mode***
    3:     30.00 cm**-1
    4:     40.00 cm**-1
    5:     50.00 cm**-1

NORMAL MODES
------------

These modes are the Cartesian displacements weighted by the diagonal matrix
M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
Thus, these vectors are normalized but *not* orthogonal

                  0          1          2          3          4          5
      0       0.000000   0.000000   0.100000   0.000000   0.000000   0.000000
      1       0.000000   0.000000   0.200000   0.000000   0.000000   0.000000
      2       0.000000   0.000000   0.300000   0.000000   0.000000   0.000000
      3       0.000000   0.000000   0.400000   0.000000   0.000000   0.000000
      4       0.000000   0.000000   0.500000   0.000000   0.000000   0.000000
      5       0.000000   0.000000   0.600000   0.000000   0.000000   0.000000

""")
    check("only the second Hessian's mode 2 survives",
          state.modes is not None and state.modes[0].index == 2 and state.modes[0].cm1 == -200.0,
          f"{state.modes}")


_MODE_FREQS = [0.0, 0.0, 0.0, -150.0, 100.0, 50.0, 300.0, 200.0, 400.0]
_MODE_VECTORS = {
    0: [0.0] * 9, 1: [0.0] * 9, 2: [0.0] * 9,
    3: [0.1, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.3],
    4: [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0],
    5: [0.0, 0.0, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0],
    6: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0],
    7: [0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    8: [0.0, 0.1, 0.0, 0.0, 0.1, 0.0, 0.0, 0.1, 0.0],
}


def _modes_block(freqs, vectors):
    """A VIBRATIONAL FREQUENCIES block and its NORMAL MODES block as ORCA
    prints them: six mode columns per sub-block, every coordinate repeated."""
    lines = ["VIBRATIONAL FREQUENCIES", "-----------------------", "",
             "Scaling factor for frequencies =  1.000000000  (already applied!)", ""]
    for i, f in enumerate(freqs):
        lines.append(f"{i:5d}:  {f:8.2f} cm**-1" + ("  ***imaginary mode***" if f < 0 else ""))
    lines += ["", "NORMAL MODES", "------------", "",
              "These modes are the Cartesian displacements weighted by the diagonal matrix",
              "M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom",
              "Thus, these vectors are normalized but *not* orthogonal", ""]
    n = len(freqs)
    for start in range(0, n, 6):
        cols = list(range(start, min(start + 6, n)))
        lines.append("            " + "".join(f"{c:11d}" for c in cols))
        for coord in range(n):
            lines.append(f"{coord:6d}     " + "".join(f"{vectors[c][coord]:11.6f}" for c in cols))
    lines += ["", ""]
    return "\n".join(lines) + "\n"


def test_every_mode_of_a_final_block_is_kept():
    print("\nevery non-zero mode of a final Hessian is kept, packed")
    state = JobState(path=Path("/nonexistent"))
    _feed(state, _modes_block(_MODE_FREQS, _MODE_VECTORS))
    state.feed_line("")  # `_feed` strips the fixture's trailing blank line
    check("every non-zero mode of the final block, by ORCA's index",
          state.normal_modes is not None
          and [m.index for m in state.normal_modes] == [3, 4, 5, 6, 7, 8],
          f"{state.normal_modes}")
    check("mode 5's printed vector, in a packed array",
          state.normal_modes is not None
          and state.normal_modes[2].vector
          == array("d", [0.0, 0.0, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0])
          and state.normal_modes[2].vector.typecode == "d")
    check("while the imaginary mode stays a list in `modes`",
          state.modes is not None and [m.index for m in state.modes] == [3]
          and isinstance(state.modes[0].vector, list))

    interim = JobState(path=Path("/nonexistent"))
    _feed(interim, "         *                GEOMETRY OPTIMIZATION CYCLE   3            *")
    _feed(interim, _modes_block(_MODE_FREQS, _MODE_VECTORS))
    interim.feed_line("")
    check("an intermediate Hessian's block keeps only the imaginary modes",
          interim.freqs_final is False and interim.normal_modes is None
          and interim.modes is not None and [m.index for m in interim.modes] == [3])

    second = dict(_MODE_VECTORS)
    second[4] = [0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0]
    replaced = JobState(path=Path("/nonexistent"))
    _feed(replaced, _modes_block(_MODE_FREQS, _MODE_VECTORS))
    _feed(replaced, _modes_block(_MODE_FREQS, second))
    replaced.feed_line("")
    check("a later final block replaces the earlier one's modes",
          replaced.normal_modes is not None and replaced.normal_modes[1].vector[4] == 0.7)

    whole = JobState(path=Path("/nonexistent"))
    whole.feed_text(_modes_block(_MODE_FREQS, _MODE_VECTORS))
    check("feed_text keeps the same modes as feeding the lines",
          whole.normal_modes is not None and state.normal_modes is not None
          and [m.index for m in whole.normal_modes] == [m.index for m in state.normal_modes]
          and all(a.vector == b.vector for a, b in zip(whole.normal_modes, state.normal_modes)))

    lines = _modes_block(_MODE_FREQS, _MODE_VECTORS).split("\n")
    after = lines.index("NORMAL MODES")
    row4 = next(i for i, ln in enumerate(lines) if i > after and ln[:6].strip() == "4")
    partial = JobState(path=Path("/nonexistent"))
    _feed(partial, "\n".join(lines[: row4 + 1]))
    partial.feed_line("")
    check("a block cut off mid-print stores no modes",
          partial.normal_modes is None)


def test_modes_are_found_and_described():
    print("\nany mode of a final Hessian is found by ORCA's number and its atoms ranked")
    from orcamon.core import vibrations

    atoms3 = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0)]
    rnd = lambda rows: [(a, e, round(r, 4)) for a, e, r in rows]

    state = JobState(path=Path("/nonexistent"))
    _feed(state, _modes_block(_MODE_FREQS, _MODE_VECTORS))
    state.feed_line("")  # `_feed` strips the fixture's trailing blank line

    check("real_modes lists the final Hessian's real modes lowest frequency first",
          [m.index for m in vibrations.real_modes(state)] == [5, 4, 7, 6, 8],
          f"{[m.index for m in vibrations.real_modes(state)]}")
    found = [vibrations.find_mode(state, i) for i in (3, 5)]
    check("find_mode returns ORCA's own mode, imaginary or real, and None otherwise",
          all(m is not None for m in found)
          and found[0].cm1 == -150.0 and found[1].cm1 == 50.0
          and vibrations.find_mode(state, 0) is None
          and vibrations.find_mode(state, 99) is None,
          f"{found}")
    check("participation ranks the imaginary mode's atoms largest first",
          rnd(vibrations.participation(atoms3, None, vibrations.find_mode(state, 3)))
          == [(2, "H", 1.0), (1, "H", 0.6667), (0, "O", 0.3333)])
    check("an atom the mode leaves still is dropped, not ranked last",
          rnd(vibrations.participation(atoms3, None, vibrations.find_mode(state, 5)))
          == [(0, "O", 1.0), (1, "H", 0.75)])
    check("equal displacements tie by atom index",
          rnd(vibrations.participation(atoms3, None, vibrations.find_mode(state, 8)))
          == [(0, "O", 1.0), (1, "H", 1.0), (2, "H", 1.0)])
    check("top keeps only the largest atom",
          rnd(vibrations.participation(atoms3, None, vibrations.find_mode(state, 3), top=1))
          == [(2, "H", 1.0)])
    check("available_text starts with the imaginary list, then the real range",
          vibrations.available_text(state) == "available: [3] imaginary; real: 5 modes, 4-8",
          f"{vibrations.available_text(state)!r}")

    interim = JobState(path=Path("/nonexistent"))
    _feed(interim, "         *                GEOMETRY OPTIMIZATION CYCLE   3            *")
    _feed(interim, _modes_block(_MODE_FREQS, _MODE_VECTORS))
    interim.feed_line("")
    check("an intermediate Hessian offers no real mode",
          vibrations.real_modes(interim) == []
          and vibrations.find_mode(interim, 5) is None)


def test_mode_offsets_fill_the_whole_structure():
    print("\na whole-structure mode moves its largest atom by MODE_AMPLITUDE_ANGSTROM")
    import math
    from orcamon.core import vibrations
    from orcamon.core.parser import NormalMode

    mode = NormalMode(index=0, cm1=-100.0, vector=[0.0, 0.0, 0.0, 3.0, 4.0, 0.0])
    atoms = [("H", 0.0, 0.0, 0.0), ("H", 1.0, 0.0, 0.0)]
    offs = vibrations.offsets(atoms, None, mode)
    check("a mode over two atoms maps directly onto the two displayed atoms",
          offs is not None)
    length = math.sqrt(sum(c * c for c in offs[1])) if offs else 0.0
    check("the 3-4-0 atom travels MODE_AMPLITUDE_ANGSTROM",
          abs(length - vibrations.MODE_AMPLITUDE_ANGSTROM) < 1e-12, f"{length}")
    check("and the atom the mode does not move stays at zero",
          offs is not None and offs[0] == (0.0, 0.0, 0.0),
          f"{offs[0] if offs else None}")


def test_mode_offsets_map_onto_the_qm_subset():
    print("\na QM-layer mode lands on the displayed QM atoms, not the environment")
    import math
    from orcamon.core import vibrations
    from orcamon.core.parser import NormalMode

    atoms = [("H", 0.0, 0.0, 0.0), ("H", 1.0, 0.0, 0.0),
             ("H", 2.0, 0.0, 0.0), ("H", 3.0, 0.0, 0.0)]
    mode = NormalMode(index=0, cm1=-100.0, vector=[3.0, 4.0, 0.0, 0.0, 0.0, 2.5])
    offs = vibrations.offsets(atoms, {1, 3}, mode)
    check("a two-atom mode maps through the two QM atoms of four",
          offs is not None)
    first = math.sqrt(sum(c * c for c in offs[1])) if offs else 0.0
    second = math.sqrt(sum(c * c for c in offs[3])) if offs else 0.0
    check("atom 1 carries the largest row at the amplitude and atom 3 the smaller one",
          abs(first - vibrations.MODE_AMPLITUDE_ANGSTROM) < 1e-12
          and 0.0 < second < vibrations.MODE_AMPLITUDE_ANGSTROM,
          f"{first} {second}")
    check("and the two environment atoms are unmoved",
          offs is not None and offs[0] == (0.0, 0.0, 0.0) and offs[2] == (0.0, 0.0, 0.0),
          f"{offs}")


def test_mode_offsets_refuse_a_mismatched_mode():
    print("\na mode matching neither the structure nor its QM subset is refused")
    from orcamon.core import vibrations
    from orcamon.core.parser import NormalMode

    mode = NormalMode(index=0, cm1=-100.0, vector=[1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
    atoms = [("H", 0.0, 0.0, 0.0), ("H", 1.0, 0.0, 0.0), ("H", 2.0, 0.0, 0.0)]
    check("three displayed atoms, no QM set, a two-atom mode: None",
          vibrations.offsets(atoms, None, mode) is None)


def test_the_sine_phase_displaces_atoms():
    print("\nthe sine phase starts at zero and displaces the atoms by the offsets")
    from orcamon.core import vibrations

    atoms = [("O", 0.0, 0.0, 0.0), ("H", 0.5, 0.0, 0.0)]
    offs = [(0.0, 0.0, 0.0), (0.25, -0.5, 0.0)]
    check("sin at t=0 is zero", vibrations.phase_sine(0.0) == 0.0)
    check("at sine 1 the atoms land on atoms + offsets",
          vibrations.displaced(atoms, offs, 1.0)
          == [("O", 0.0, 0.0, 0.0), ("H", 0.75, -0.5, 0.0)],
          f"{vibrations.displaced(atoms, offs, 1.0)}")


def test_the_mode_geometry_falls_back_to_an_alternate():
    print("\na mode that misses the pane's geometry is drawn on the alternate")
    from orcamon.core import vibrations
    from orcamon.core.parser import NormalMode

    pane = [("H", 0.0, 0.0, 0.0), ("H", 1.0, 0.0, 0.0)]
    alternate = [("H", 0.0, 0.0, 0.0), ("H", 1.0, 0.0, 0.0),
                 ("H", 2.0, 0.0, 0.0), ("H", 3.0, 0.0, 0.0)]
    mode4 = NormalMode(index=0, cm1=-100.0, vector=[3.0, 4.0, 0.0, 0.0, 0.0, 2.5,
                                                     0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
    mode2 = NormalMode(index=0, cm1=-100.0, vector=[3.0, 4.0, 0.0, 0.0, 0.0, 2.5])
    check("a four-atom mode takes the four-atom alternate over the two-atom pane",
          vibrations.mode_geometry([pane, alternate], None, mode4) == (alternate, None),
          f"{vibrations.mode_geometry([pane, alternate], None, mode4)}")
    check("with only the pane it is refused",
          vibrations.mode_geometry([pane], None, mode4) is None)
    check("and a two-atom mode still takes the pane",
          vibrations.mode_geometry([pane, alternate], None, mode2) == (pane, None),
          f"{vibrations.mode_geometry([pane, alternate], None, mode2)}")


def test_the_qm_map_needs_indices_that_fit():
    print("\nQM indices that run past the drawn structure do not map onto it")
    from orcamon.core import vibrations
    from orcamon.core.parser import NormalMode

    pane = [("H", float(i), 0.0, 0.0) for i in range(19)]
    alternate = [("H", float(i), 0.0, 0.0) for i in range(137)]
    mode = NormalMode(index=6, cm1=-379.28, vector=[0.01] * 51)
    qm = {61, 64, 120, 121, 122, 124, 125, 126, 128, 129, 130, 131, 132, 133, 134, 135, 136}
    check("a 17-index QM map running to 136 does not fit a 19-atom pane",
          vibrations.target_indices(19, 17, qm) is None,
          f"{vibrations.target_indices(19, 17, qm)}")
    check("so the mode is drawn on the 137-atom alternate instead",
          vibrations.mode_geometry([pane, alternate], qm, mode) == (alternate, qm),
          f"{vibrations.mode_geometry([pane, alternate], qm, mode)}")


def test_a_frame_is_a_small_faithful_palette_png():
    print("\na frame goes out as a palette PNG, rendered without matplotlib")
    import io
    import subprocess
    import textwrap
    import numpy as np
    from PIL import Image
    from orcamon.core.geometry import View
    from orcamon.tui import geometry_render

    water = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0)]
    image = geometry_render.render(water, size_px=(300, 240))
    check("the render is RGB at the size asked for",
          image.mode == "RGB" and image.size == (300, 240), f"{image.mode} {image.size}")

    sent = Image.open(io.BytesIO(geometry_render.frame_png(image)))
    error = np.abs(np.asarray(sent.convert("RGB"), int) - np.asarray(image, int)).mean()
    check("the frame is a 256-colour palette PNG", sent.mode == "P" and sent.size == image.size,
          f"{sent.mode} {sent.size}")
    check("within ~1 level in 255 of the render", error < 2.0, f"mean error {error:.2f}")

    bare = np.asarray(geometry_render.render(water, View(show_labels=False), size_px=(300, 240)), int)
    labelled = np.asarray(image, int)
    changed = (np.abs(bare - labelled).sum(axis=2) > 0).mean()
    check("show_labels=False drops the index labels and nothing else of the frame",
          0 < changed < 0.05, f"{changed:.3%} of pixels differ")

    # matplotlib is no longer a dependency: the render must work with it made
    # unimportable, the same technique as test_orcamon._STDLIB_ONLY_PROBE. A
    # stray top-level import would otherwise only surface on a machine that
    # installed the images extra without matplotlib.
    probe = textwrap.dedent("""
        import sys
        class Refuse:
            def find_spec(self, name, path=None, target=None):
                if name.split(".")[0] == "matplotlib":
                    raise ImportError("no matplotlib")
                return None
        sys.meta_path.insert(0, Refuse())
        from orcamon.tui import geometry_render
        img = geometry_render.render([("O", 0., 0., 0.), ("H", .96, 0., 0.)], size_px=(60, 50))
        print(img.mode, img.size[0], img.size[1])
    """)
    done = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True, cwd=HERE)
    check("render imports no matplotlib", done.returncode == 0 and done.stdout.strip() == "RGB 60 50",
          done.stdout + done.stderr)


def test_near_atoms_hide_far_ones():
    print("\nnear atoms hide far ones, and a bond in front wins its pixel")
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # The eye is on +x at (elev, azim) = (0, 0), so screen x carries +y and
    # screen y carries +z; both atoms below sit ON the view axis, so both
    # project to the image centre and only depth decides the pixel.
    image = raster.render([("O", 0.0, 0.0, 0.0), ("C", 2.0, 0.0, 0.0)],
                          View(elev=0, azim=0), size_px=(300, 240))
    pixel = image.getpixel((150, 120))
    check("the nearer carbon hides the farther oxygen: the centre pixel is grey",
          abs(int(pixel[0]) - int(pixel[2])) < 20 and abs(int(pixel[0]) - int(pixel[1])) < 20,
          f"{pixel}")

    # The two carbons are 1.6 A apart so they BOND (the cutoff is 1.7); their
    # midpoint is (2, 0, 0), on the same sight line as the nitrogen at the
    # origin but 2 A nearer, so the bond's samples must win the centre pixel.
    image = raster.render([("N", 0.0, 0.0, 0.0), ("C", 2.0, -0.8, 0.0), ("C", 2.0, 0.8, 0.0)],
                          View(elev=0, azim=0), size_px=(300, 240))
    pixel = image.getpixel((150, 120))
    check("a bond passing in front of an atom wins the pixel: grey, not blue",
          abs(int(pixel[0]) - int(pixel[2])) < 20 and int(pixel[2]) < 220, f"{pixel}")


def test_spheres_are_shaded():
    print("\nspheres are lit: the side toward LIGHT is brighter")
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    image = raster.render([("C", 0.0, 0.0, 0.0)], View(elev=0, azim=0), size_px=(300, 240))
    upper_left = image.getpixel((120, 90))
    lower_right = image.getpixel((180, 150))
    check("the upper-left of the disc (toward LIGHT) is brighter than the lower-right",
          sum(upper_left) > sum(lower_right), f"{upper_left} vs {lower_right}")


def test_fog_dims_the_far_side():
    print("\nfog dims the farther atom, and does nothing when off")
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # Eye on +x. The two atoms are side by side on screen (+/-y) and differ
    # only in depth (+/-x): the left is near, the right far.
    atoms = [("C", 2.0, -1.0, 0.0), ("C", -2.0, 1.0, 0.0)]
    view = View(elev=0, azim=0, show_labels=False)
    sx, sy, _depth, _scale = raster.project(atoms, view, (300, 240))
    near = (int(round(sx[0])), int(round(sy[0])))
    far = (int(round(sx[1])), int(round(sy[1])))
    foggy = raster.render(atoms, view, size_px=(300, 240))
    check("with fog, the near atom's centre is brighter than the far one",
          sum(foggy.getpixel(near)) > sum(foggy.getpixel(far)),
          f"{foggy.getpixel(near)} vs {foggy.getpixel(far)}")
    clear = raster.render(atoms, View(elev=0, azim=0, fog=False, show_labels=False), size_px=(300, 240))
    check("without fog they are equal within 2 levels",
          max(abs(a - b) for a, b in zip(clear.getpixel(near), clear.getpixel(far))) <= 2,
          f"{clear.getpixel(near)} vs {clear.getpixel(far)}")


def test_outlines_separate_overlapping_atoms():
    print("\nan outline darkens where a near atom overlaps a far one")
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    atoms = [("C", 2.0, 0.3, 0.0), ("C", -2.0, -0.3, 0.0)]
    view = View(elev=0, azim=0, fog=False, show_labels=False)
    image = raster.render(atoms, view, size_px=(300, 240))
    sx, sy, _depth, _scale = raster.project(atoms, view, (300, 240))
    centres = min(sum(image.getpixel((int(round(sx[i])), int(round(sy[i]))))) for i in (0, 1))
    row = int(round((sy[0] + sy[1]) / 2))
    lo, hi = sorted((int(round(sx[0])), int(round(sx[1]))))
    darkest = min(sum(image.getpixel((x, row))) for x in range(lo, hi + 1))
    check("between the centres there is a pixel under half their brightness",
          darkest < 0.5 * centres, f"darkest {darkest} vs centres {centres}")


def test_representations_change_what_is_drawn():
    print("\nthe four representations draw different amounts of ink")
    import numpy as np
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    water = [("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0)]
    background = np.asarray(raster.BACKGROUND)

    def ink(representation):
        image = np.asarray(raster.render(
            water, View(elev=20, azim=-60, representation=representation, show_labels=False),
            size_px=(300, 240)))
        return int((np.abs(image - background).sum(axis=2) > 0).sum())

    space, ball, wire, licorice = (ink(r) for r in
                                   ("space-filling", "ball-and-stick", "wireframe", "licorice"))
    check("space-filling covers more than ball-and-stick, which covers more than wireframe",
          space > ball > wire, f"space {space} ball {ball} wire {wire}")
    check("licorice differs from ball-and-stick", licorice != ball, f"{licorice} vs {ball}")


def test_space_filling_is_framed():
    print("\nspace-filling is framed so its spheres do not touch the border")
    import numpy as np
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    atoms = [("O", 0.0, 0.0, 0.0), ("H", 1.4, 0.0, 0.0), ("H", -1.4, 0.0, 0.0),
             ("H", 0.0, 1.4, 0.0), ("H", 0.0, -1.4, 0.0)]
    image = np.asarray(raster.render(
        atoms, View(representation="space-filling", show_labels=False), size_px=(200, 160)))
    background = np.asarray(raster.BACKGROUND)
    border = np.concatenate([image[0].reshape(-1, 3), image[-1].reshape(-1, 3),
                             image[:, 0].reshape(-1, 3), image[:, -1].reshape(-1, 3)])
    check("no foreground pixel on any edge",
          bool((np.abs(border - background).sum(axis=1) == 0).all()))


def test_principal_axes_face_a_plane_on():
    print("\nprincipal axes: the small one is the plane normal, and views face-on")
    import math
    from orcamon.core.geometry import View, principal_axes, view_along
    from orcamon.tui import raster

    normal = (1 / math.sqrt(3),) * 3
    u = (1 / math.sqrt(2), -1 / math.sqrt(2), 0.0)
    w = (1 / math.sqrt(6), 1 / math.sqrt(6), -2 / math.sqrt(6))
    coords = [tuple(1.4 * (math.cos(t) * u[i] + math.sin(t) * w[i]) for i in range(3))
              for t in (k * math.pi / 3 for k in range(6))]
    axes = principal_axes(coords)  # largest variance first
    smallest = axes[-1][0]
    dot = abs(sum(smallest[i] * normal[i] for i in range(3)))
    check("the smallest-variance axis is the plane normal", dot > 0.999, f"dot {dot}")

    n = len(coords)
    mean = [sum(p[i] for p in coords) / n for i in range(3)]
    covariance = [[sum((p[i] - mean[i]) * (p[j] - mean[j]) for p in coords) / n
                   for j in range(3)] for i in range(3)]
    ok = True
    for axis, variance in axes:
        for i in range(3):
            av = sum(covariance[i][j] * axis[j] for j in range(3))
            if abs(av - variance * axis[i]) > 1e-9:
                ok = False
    check("every axis satisfies A v = lambda v", ok)

    elev, azim = view_along(smallest)
    atoms = [("C",) + c for c in coords]
    sx, sy, _depth, _scale = raster.project(atoms, View(elev=elev, azim=azim), (300, 240))
    extent_x, extent_y = float(max(sx) - min(sx)), float(max(sy) - min(sy))
    check("viewing along it is face-on: screen extents within 20%",
          abs(extent_x - extent_y) / max(extent_x, extent_y) < 0.2,
          f"{extent_x:.1f} vs {extent_y:.1f}")


def test_labels_follow_occlusion():
    print("\na label does not float over an atom hidden behind another")
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # The big near potassium at x=3 completely covers the small far hydrogen
    # at the origin (3 A apart, so no bond): the hydrogen's projected centre
    # is well inside the potassium's disc.
    atoms = [("K", 3.0, 0.3, 0.0), ("H", 0.0, 0.0, 0.0)]
    view = View(elev=0, azim=0)
    sx, sy, _depth, _scale = raster.project(atoms, view, (300, 240))
    labelled = raster.render(atoms, view, size_px=(300, 240))
    bare = raster.render(atoms, View(elev=0, azim=0, show_labels=False), size_px=(300, 240))
    cx, cy = int(round(sx[1])), int(round(sy[1]))
    same = all(labelled.getpixel((x, y)) == bare.getpixel((x, y))
               for x in range(cx - 3, cx + 4) for y in range(cy - 3, cy + 4))
    check("the fully hidden atom gets no label", same, f"around {(cx, cy)}")


def test_wireframe_keeps_labels_and_distances():
    print("\nwireframe still draws the labels and distances it was asked for")
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # Water: every atom is bonded, so in wireframe NONE has a sphere -- which
    # is what made `_is_visible` refuse all of them and `l`/`d` draw nothing.
    atoms = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.96, 0.0), ("H", 0.0, -0.24, 0.93)]
    size = (300, 240)
    sx, sy, _depth, _scale = raster.project(atoms, View(elev=0, azim=0, representation="wireframe"), size)

    def differs(a, b, cx, cy):
        cx, cy = int(round(cx)), int(round(cy))
        return any(a.getpixel((x, y)) != b.getpixel((x, y))
                   for x in range(cx - 4, cx + 5) for y in range(cy - 4, cy + 5))

    def wire(**kwargs):
        return raster.render(atoms, View(elev=0, azim=0, representation="wireframe", fog=False,
                                         **kwargs), size_px=size)

    bare = wire(show_labels=False)
    check("the oxygen's label is drawn", differs(wire(), bare, sx[0], sy[0]))
    measured = wire(show_labels=False, show_distances=True)
    check("the O-H distance is drawn",
          differs(measured, bare, (sx[0] + sx[1]) / 2, (sy[0] + sy[1]) / 2))


def test_wireframe_draws_the_qm_layer_thicker():
    print("\nwireframe draws the QM layer's sticks thicker than the environment's")
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # Eye on +x: two equal-length bonds side by side on screen, the left one
    # QM, the right one environment, 3 A apart so they do not bond each other.
    atoms = [("C", 0.0, -3.0, 0.0), ("C", 0.0, -1.5, 0.0),
             ("C", 0.0, 1.5, 0.0), ("C", 0.0, 3.0, 0.0)]
    size = (400, 200)
    image = raster.render(atoms, View(elev=0, azim=0, representation="wireframe", fog=False,
                                      show_labels=False),
                          qm_atom_indices={0, 1}, size_px=size)

    def ink(x0, x1):
        return sum(1 for x in range(x0, x1) for y in range(size[1])
                   if sum(abs(c - b) for c, b in zip(image.getpixel((x, y)), raster.BACKGROUND)) > 6)

    qm, env = ink(0, size[0] // 2), ink(size[0] // 2, size[0])
    check("the QM bond covers more pixels than the environment bond", qm > env, f"{qm} vs {env}")


def test_fog_holds_still_under_rotation():
    print("\nfog fades an atom by its depth, not by the frame's own depth range")
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # An oxygen at the centroid of a carbon triangle (3 A out, no bonds): its
    # depth is zero from every angle, but the frame's nearest and farthest
    # pixels are not -- flat seen from above, +/-2.6 A seen edge-on. Fog read
    # off those dimmed the oxygen differently at the two angles.
    atoms = [("O", 0.0, 0.0, 0.0), ("C", 3.0, 0.0, 0.0),
             ("C", -1.5, 2.6, 0.0), ("C", -1.5, -2.6, 0.0)]
    size = (300, 240)
    colours = []
    for elev, azim in ((90, -60), (0, 90)):
        view = View(elev=elev, azim=azim, show_labels=False)
        sx, sy, _depth, _scale = raster.project(atoms, view, size)
        image = raster.render(atoms, view, size_px=size)
        colours.append(image.getpixel((int(round(sx[0])), int(round(sy[0])))))
    check("the oxygen's centre is the same colour from above and edge-on",
          max(abs(a - b) for a, b in zip(*colours)) <= 2, f"{colours[0]} vs {colours[1]}")


def test_host_fog_keeps_near_host_bright():
    print("\nthe host is fogged on its own: near host bright, far host faded hard")
    import numpy as np
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # Two host C-C sticks, one 3 A in front of the guest's plane and one 3 A
    # behind it, and a lone QM oxygen at the centroid. The framing sphere of
    # `_frame_sphere` has extent |(3, -3, -0.5)| + 0.5 = 4.772, so near =
    # 4.772 and far = -4.772.
    #
    # The near stick's surface sits at depth about 3.06 (its ENV_BOND_RADIUS
    # toward the eye): t = (4.772 - 3.06) / 9.544 = 0.179 and the host fog
    # 0.85 * t**2 = 0.027. The old single linear fog gave 0.6 * t = 0.108,
    # which is why a host in front read no nearer than one behind.
    # The far stick at -2.94: t = 0.808 and 0.85 * t**2 = 0.555.
    # The QM oxygen's ball, radius 0.28 * 1.52 = 0.426, gives t = 0.455 and
    # the guest's unchanged 0.6 * t = 0.273.
    atoms = [("C", 3.0, -3.0, 0.0), ("C", 3.0, -1.6, 0.0),
             ("C", -3.0, 1.6, 0.0), ("C", -3.0, 3.0, 0.0), ("O", 0.0, 0.0, 2.5)]
    view = View(elev=0, azim=0, show_labels=False)
    sx, sy, _depth, _scale = raster.project(atoms, view, (400, 300))
    on = np.asarray(raster.render(atoms, view, qm_atom_indices={4}, size_px=(400, 300)), float)
    off = np.asarray(raster.render(atoms, View(elev=0, azim=0, fog=False, show_labels=False),
                                   qm_atom_indices={4}, size_px=(400, 300)), float)

    def fade(x, y):
        x, y = int(round(x)), int(round(y))
        return (off[y, x, 0] - on[y, x, 0]) / (off[y, x, 0] - 30.0)

    near_mid = ((sx[0] + sx[1]) / 2, (sy[0] + sy[1]) / 2)
    far_mid = ((sx[2] + sx[3]) / 2, (sy[2] + sy[3]) / 2)
    check("the near host stick stays bright: 0.85 * 0.179**2 = 0.027",
          abs(fade(*near_mid) - 0.027) < 0.02, f"{fade(*near_mid):.3f}")
    check("the far host stick fades hard: 0.85 * 0.808**2 = 0.555",
          abs(fade(*far_mid) - 0.555) < 0.02, f"{fade(*far_mid):.3f}")
    check("the guest keeps the old linear fog: 0.6 * 0.455 = 0.273",
          abs(fade(sx[4], sy[4]) - 0.273) < 0.02, f"{fade(sx[4], sy[4]):.3f}")


def test_halos_widen_with_the_depth_gap():
    print("\nthe halo on a far atom widens with the depth gap in front of it")
    import numpy as np
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # A potassium ball (vdW 2.75 A) faces the eye at (elev, azim) = (0, 0)
    # with a host C-C stick across its centre `gap` A in front of the ball's
    # surface. On the committed renderer the single 1-px outline darkened one
    # row on each side for every gap, so a stick 3 A in front of the ball read
    # no nearer than one touching it.
    #
    # The halo has five levels on a 750-px frame, with thresholds T_d of
    # 0.35, 0.8875, 1.425, 1.9625 and 2.5 A at d = 1..5, so the ring is as
    # wide as the number of thresholds below the gap: 1, 2 and 5 here.
    for gap, width in ((0.6, 1), (1.2, 2), (3.0, 5)):
        atoms = [("K", 0.0, 0.0, 0.0), ("C", 2.75 + gap, -0.7, 0.0),
                 ("C", 2.75 + gap, 0.7, 0.0)]
        view = View(elev=0, azim=0, representation="space-filling", fog=False,
                    show_labels=False)
        sx, sy, _depth, _scale = raster.project(atoms, view, (750, 750))
        cx, cy = int(round(sx[0])), int(round(sy[0]))
        im = np.asarray(raster.render(atoms, view, qm_atom_indices={0},
                                      size_px=(750, 750)), float)
        ref = im[cy - 60, cx].sum()

        def kind(y):
            r, g, _b = im[y, cx]
            if g < 0.6 * r:       # K is #8F40D4 (g/r = 0.45), and darkened
                return "halo" if im[y, cx].sum() < 0.5 * ref else "ball"
            return "stick"        # K keeps that ratio; the grey stick has g/r = 1

        rows = [kind(y) for y in range(cy - 14, cy + 15)]
        top = rows.index("stick")
        bottom = len(rows) - 1 - rows[::-1].index("stick")
        above = 0
        while top - 1 - above >= 0 and rows[top - 1 - above] == "halo":
            above += 1
        below = 0
        while bottom + 1 + below < len(rows) and rows[bottom + 1 + below] == "halo":
            below += 1
        check(f"a {gap} A gap darkens {width} row(s) above the stick",
              above == width, f"{above}")
        check(f"a {gap} A gap darkens {width} row(s) below the stick",
              below == width, f"{below}")


def test_boundary_host_atoms_get_a_ball():
    print("\na host atom bonded to the QM region gets a small ball in ball-and-stick")
    import numpy as np
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # A QM oxygen with two host carbons along +x at 1.4 A spacing. The eye is
    # on +z (elev=90, azim=-90), so x runs to the screen's right. Atom 1 is
    # bonded to the QM atom and is a boundary atom; atom 2 is one bond further
    # out and keeps only its stick end.
    #
    # The framing sphere's mean x is 1.4 and its extent is 1.4 + 0.5 = 1.9,
    # so scale = 400 / 3.8 = 105.26 px/A. The boundary ball's radius is
    # 0.6 * 0.28 * 1.70 = 0.2856 A = 30 px, so its column spans 61 rows; a
    # host stick end of radius ENV_BOND_RADIUS = 0.06 A = 6 px gives 13 rows.
    atoms = [("O", 0.0, 0.0, 0.0), ("C", 1.4, 0.0, 0.0), ("C", 2.8, 0.0, 0.0)]
    ball_view = View(elev=90, azim=-90, fog=False, show_labels=False)
    sx, _sy, _depth, _scale = raster.project(atoms, ball_view, (400, 400))

    def column_rows(image, i):
        column = np.asarray(image, float)[:, int(round(sx[i]))]
        return int((abs(column - 30.0).sum(axis=1) > 6).sum())

    ball = raster.render(atoms, ball_view, qm_atom_indices={0}, size_px=(400, 400))
    check("the boundary carbon is drawn as a ball: 0.6 * 0.28 * 1.70 = 61 rows",
          57 <= column_rows(ball, 1) <= 65, f"{column_rows(ball, 1)}")
    check("the carbon one bond further out keeps only the stick end: 13 rows",
          11 <= column_rows(ball, 2) <= 15, f"{column_rows(ball, 2)}")
    licorice = raster.render(atoms, View(elev=90, azim=-90, representation="licorice",
                                         fog=False, show_labels=False),
                             qm_atom_indices={0}, size_px=(400, 400))
    check("licorice gets no boundary ball: the same 13 rows",
          11 <= column_rows(licorice, 1) <= 15, f"{column_rows(licorice, 1)}")


def test_see_through_shows_the_guest_behind_the_host():
    print("\nsee-through blends the covered guest through the host in front of it")
    import numpy as np
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # A guest oxygen with a host C-C stick 2 A in front of it, and the same
    # stick 2 A behind it. The two structures have the same mean distance, so
    # the same framing: same scale, and the oxygen at the same pixel. Without
    # see-through the first shows only the grey stick there and the second
    # only the red oxygen; with see-through the first shows the oxygen
    # through the stick. The prototype measured (76, 76, 76) for the stick,
    # (212, 16, 16) for the oxygen and (171, 34, 34) blended.
    front = [("O", 0.0, 0.0, 0.0), ("C", 2.0, -0.7, 0.0), ("C", 2.0, 0.7, 0.0)]
    back = [("O", 0.0, 0.0, 0.0), ("C", -2.0, -0.7, 0.0), ("C", -2.0, 0.7, 0.0)]
    view = View(elev=0, azim=0, fog=False, show_labels=False)
    size = (400, 300)
    sx, sy, _depth, _scale = raster.project(front, view, size)
    px, py = int(round(sx[0])), int(round(sy[0]))

    def pixel(atoms, see_through=False):
        image = raster.render(atoms, View(elev=0, azim=0, fog=False, show_labels=False,
                                          see_through=see_through),
                              qm_atom_indices={0}, size_px=size)
        return np.asarray(image, float)[py, px]

    c_env, c_qm = pixel(front), pixel(back)
    c_x = pixel(front, see_through=True)
    expected = 0.3 * c_env + 0.7 * c_qm
    check("the covered QM oxygen shows through at 0.3 host + 0.7 guest",
          all(abs(a - b) <= 2.0 for a, b in zip(c_x, expected)), f"{c_x} vs {expected}")
    check("and the blend is visibly red, not the grey stick",
          c_x[0] > c_env[0] + 50, f"{c_x[0]} vs {c_env[0]}")


def test_see_through_ghosts_the_environment():
    print("\nsee-through ghosts the whole environment layer, not only where it covers a guest")
    import numpy as np
    from orcamon.core.geometry import View
    from orcamon.tui import raster

    # `test_host_fog_keeps_near_host_bright`'s atoms: a near host C-C stick
    # (atoms 0-1) and a far one (atoms 2-3), with a lone guest oxygen (atom 4)
    # off to the side of both. The eye is on +x (elev=0, azim=0), so the near
    # stick is at depth +3 and the far one at -3; neither midpoint has a guest
    # surface behind it. fog=False because the ghost blends against the raw
    # background, and the host fog would otherwise fade `off` by an amount the
    # derivation below does not carry.
    atoms = [("C", 3.0, -3.0, 0.0), ("C", 3.0, -1.6, 0.0),
             ("C", -3.0, 1.6, 0.0), ("C", -3.0, 3.0, 0.0), ("O", 0.0, 0.0, 2.5)]
    view = View(elev=0, azim=0, fog=False, show_labels=False)
    ghost = View(elev=0, azim=0, fog=False, show_labels=False, see_through=True)
    off = np.asarray(raster.render(atoms, view, qm_atom_indices={4}, size_px=(400, 300)), float)
    on = np.asarray(raster.render(atoms, ghost, qm_atom_indices={4}, size_px=(400, 300)), float)
    sx, sy, _depth, _scale = raster.project(atoms, view, (400, 300))

    # With see-through off an environment pixel is its full colour. With it on
    # the environment layer is composited at SEE_THROUGH_ALPHA = 0.3 against
    # whatever lies behind: BACKGROUND = (30, 30, 30) at a stick midpoint, so
    # each channel is 0.3 * off + 0.7 * 30. The guest oxygen is in the QM
    # layer, not the environment, so where it wins its colour is unchanged.
    near_mid = ((sx[0] + sx[1]) / 2, (sy[0] + sy[1]) / 2)
    far_mid = ((sx[2] + sx[3]) / 2, (sy[2] + sy[3]) / 2)
    near_x, near_y = int(round(near_mid[0])), int(round(near_mid[1]))
    far_x, far_y = int(round(far_mid[0])), int(round(far_mid[1]))
    near_expected = 0.3 * off[near_y, near_x] + 0.7 * 30.0
    far_expected = 0.3 * off[far_y, far_x] + 0.7 * 30.0
    check("the near host stick ghosts to 0.3 * off + 0.7 * 30",
          max(abs(on[near_y, near_x] - near_expected)) <= 2.0,
          f"{on[near_y, near_x]} vs {near_expected}")
    check("the far host stick ghosts to 0.3 * off + 0.7 * 30",
          max(abs(on[far_y, far_x] - far_expected)) <= 2.0,
          f"{on[far_y, far_x]} vs {far_expected}")
    gx, gy = int(round(sx[4])), int(round(sy[4]))
    check("the guest oxygen's own pixel is unchanged",
          max(abs(on[gy, gx] - off[gy, gx])) <= 2.0, f"{on[gy, gx]} vs {off[gy, gx]}")


def test_the_camera_tumbles_over_the_pole():
    print("\nthe camera turns smoothly over the pole instead of snapping upside down")
    from orcamon.core.geometry import camera_basis, camera_forward

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

    worst_step, worst_hand = 1.0, 0.0
    for azim in (-60, 30, 110):
        for pole in (90, 270):
            previous = camera_basis(pole - 10, azim)
            for elev in range(pole - 9, pole + 11):
                right, up = camera_basis(elev, azim)
                worst_step = min(worst_step, dot(right, previous[0]), dot(up, previous[1]))
                normal = cross(right, up)
                worst_hand = max(worst_hand, max(abs(n - f) for n, f in
                                                 zip(normal, camera_forward(elev, azim))))
                previous = (right, up)
    # One degree of elevation turns `up` by one degree: cos(1 deg) = 0.99985.
    check("right and up move by about a degree per degree, through both poles",
          worst_step > 0.999, f"worst dot {worst_step:.4f}")
    check("and stay right-handed there", worst_hand < 1e-9, f"{worst_hand:.2e}")


def test_the_renderer_is_fast_enough():
    print("\nthe renderer keeps a rotation interactive")
    import time
    from orcamon.core.geometry import View
    from orcamon.tui import geometry_render

    # A compact 5x4x7 cubic lattice at 1.4 A spacings -- 140 atoms, of which
    # 20 are the QM layer -- so there are many bonds to draw, the slow case.
    atoms = [("C" if (ix + iy + iz) % 5 else "O", ix * 1.4, iy * 1.4, iz * 1.4)
             for ix in range(5) for iy in range(4) for iz in range(7)]
    qm = set(range(20))
    view = View()
    image = None
    runs = []
    for _ in range(3):
        started = time.perf_counter()
        image = geometry_render.render(atoms, view, qm_atom_indices=qm, size_px=(940, 900))
        runs.append((time.perf_counter() - started) * 1000)
    # The FASTEST of three: the gate suite runs eight gates in parallel, so an
    # average measures the machine's load, not the renderer. The bound is
    # deliberately loose (a real regression is a multiple, not a few percent).
    elapsed = min(runs)
    check("140 atoms at 940x900 in under 250 ms (generous; measured below)",
          elapsed < 250, f"{elapsed:.1f} ms of {[round(r) for r in runs]}")
    size = len(geometry_render.frame_png(image))
    check("and its frame_png is under 80 KB", size < 80 * 1024, f"{size} bytes")


def test_the_view_is_not_mirrored():
    print("\nthe camera basis is right-handed: screen right really is right")
    from orcamon.core.geometry import camera_basis, camera_forward

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    for elev, azim in [(0, 0), (20, -60), (35, 110), (-40, 200)]:
        right, up = camera_basis(elev, azim)
        forward = camera_forward(elev, azim)
        normal = cross(right, up)
        check(f"cross(right, up) == forward at ({elev}, {azim})",
              max(abs(normal[i] - forward[i]) for i in range(3)) < 1e-9,
              f"{normal} vs {forward}")


def test_text_mode_draws_the_right_enantiomer():
    print("\ntext mode draws the molecule's own enantiomer, not its mirror")
    from orcamon.tui import geometry_text

    atoms = [("O", 0.0, 0.0, 0.0), ("N", 0.0, 1.0, 0.0), ("S", 0.0, 0.0, 1.0)]
    text = geometry_text.render(atoms, 21, 11, geometry_text.View(elev=0, azim=0)).plain
    rows = text.split("\n")

    def cell(symbol):
        for row, line in enumerate(rows):
            col = line.find(symbol)
            if col >= 0:
                return row, col
        return None

    o, n, s = cell("O"), cell("N"), cell("S")
    check("with the eye on +x, +y is screen right: the N glyph is right of the O",
          o is not None and n is not None and n[1] > o[1], f"O {o} N {n}")
    check("and +z is screen up: the S glyph is above the O",
          o is not None and s is not None and s[0] < o[0], f"O {o} S {s}")


def test_every_marker_reaches_its_parser():
    print("\nevery marker line reaches its parser past the prefilter")
    check("validate.MARKER_CASES all read", validate.check_markers() == 0)


_IRC_OUT = """
Max. no of cycles        MaxIter    .... 3
Storing full IRC trajectory in      .... job_IRC_Full_trj.xyz
Storing forward trajectory in       .... job_IRC_F_trj.xyz
Storing backward trajectory in      .... job_IRC_B_trj.xyz
FINAL SINGLE POINT ENERGY      -100.000000000000

         *************************************************************
         *                          FORWARD IRC                      *
         *************************************************************

Iteration    E(Eh)      dE(kcal/mol)  max(|G|)   RMS(G)  B(O 0,H 1) B(O 0,H 2)
Convergence thresholds                0.002000  0.000500
    0     -100.001594   -1.000000    0.010000  0.001000      0.97         0.96
    1     -100.003187   -2.000000    0.008000  0.000800      0.98         0.95
    2     -100.004781   -3.000000    0.006000  0.000600      0.99         0.94

         *************************************************************
         *  MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC RUN  *
         *************************************************************


         *************************************************************
         *                          BACKWARD IRC                     *
         *************************************************************

Iteration    E(Eh)      dE(kcal/mol)  max(|G|)   RMS(G)  B(O 0,H 1) B(O 0,H 2)
Convergence thresholds                0.002000  0.000500
    0     -100.000797   -0.500000    0.010000  0.001000      0.95         0.97
    1     -100.002390   -1.500000    0.008000  0.000800      0.94         0.98

                             ****ORCA TERMINATED NORMALLY****
"""


def _irc_frame(z):
    """One frame of the synthetic tree, identified by the O atom's z."""
    lines = ["3", "Coordinates from ORCA-job job E -100.000000"]
    for el, x, y in (("O", 0.0, 0.0), ("H", 0.96, 0.0), ("H", -0.24, 0.93)):
        lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines) + "\n"


def _irc_tree(d: Path):
    """§ Facts' synthetic IRC tree: the finished output and four trajectories
    whose frame order is read off the O atom's z."""
    d.mkdir(parents=True, exist_ok=True)
    (d / "job.inp").write_text("! B97-3c IRC\n%irc MaxIter 3 end\n* xyzfile 0 1 start.xyz\n")
    (d / "job.out").write_text(_IRC_OUT.strip("\n") + "\n")
    (d / "start.xyz").write_text(_irc_frame(2.0))
    (d / "job_IRC_Full_trj.xyz").write_text("".join(_irc_frame(z) for z in range(6)))
    (d / "job_IRC_F_trj.xyz").write_text("".join(_irc_frame(z) for z in (3, 4, 5)))
    (d / "job_IRC_B_trj.xyz").write_text("".join(_irc_frame(z) for z in (1, 0)))


def test_irc_rows_are_parsed():
    print("\nan IRC's rows, monitors and trajectory files are read")
    state = JobState(path=Path("/nonexistent"), stem="job")
    _feed(state, _IRC_OUT)

    check("the three trajectory files the log names",
          state.irc_files == {"full": "job_IRC_Full_trj.xyz",
                              "forward": "job_IRC_F_trj.xyz",
                              "backward": "job_IRC_B_trj.xyz"})
    check("the monitored internals the header names",
          state.irc_monitors == ["B(O 0,H 1)", "B(O 0,H 2)"])
    check("the per-direction rows, with ORCA's own dE and monitors",
          len(state.irc_rows["forward"]) == 3 and len(state.irc_rows["backward"]) == 2
          and state.irc_rows["forward"][2] == IrcRow(2, -100.004781, -3.0, [0.99, 0.94]))
    check("and the backward block keeps its own rows in printed order",
          state.irc_rows["backward"][0].monitors == [0.95, 0.97])
    check("only the direction whose MaxIter banner was printed is listed",
          state.irc_maxiter == ["forward"])

    fast = JobState(path=Path("/nonexistent"), stem="job")
    fast.feed_text(_IRC_OUT.strip("\n") + "\n")
    check("the whole-chunk path reads the same rows",
          fast.irc_rows == state.irc_rows and fast.irc_maxiter == state.irc_maxiter
          and fast.irc_files == state.irc_files)

    lines = _IRC_OUT.strip("\n").split("\n")
    cut = next(i for i, line in enumerate(lines)
               if line.startswith("    0     -100.000797")) + 1
    partial = JobState(path=Path("/nonexistent"), stem="job")
    _feed(partial, "\n".join(lines[:cut]))
    check("a direction still being written holds only what has been read",
          len(partial.irc_rows["backward"]) == 1 and partial._in_block())


def test_the_irc_path_is_built():
    print("\nan IRC's path is built from its rows and trajectories, once")
    import tempfile

    from orcamon.core.geometry import read_xyz_frames
    from orcamon.core.paths import PathPoint, reaction_path
    from orcamon.core.report import describe_point

    gone = Liveness(False, "process")
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        _irc_tree(d)
        job = Job(d, "job", d.parent, label="irc")
        job.refresh(gone)
        view = reaction_path(job)

        check("the path is an IRC, in signed point order",
              view is not None and view.kind == "irc"
              and [p.index for p in view.points] == [-2, -1, 0, 1, 2, 3])
        check("labelled as ORCA numbers the iterations",
              [p.label for p in view.points] == [
                  "IRC backward 1", "IRC backward 0", "IRC TS",
                  "IRC forward 0", "IRC forward 1", "IRC forward 2"])
        check("with ORCA's own dE, converted from kcal/mol",
              [round(p.de_kj_mol, 3) for p in view.points]
              == [-6.276, -2.092, 0.0, -4.184, -8.368, -12.552])
        check("and each point's geometry from the full trajectory, in path order",
              [p.atoms[0][3] for p in view.points] == [0, 1, 2, 3, 4, 5])
        check("the TS carries the job's energy and is the focus of a finished job",
              view.focus == 2 and view.points[2].energy == -100.0)
        check("the last point's monitors are the row's, named by the header",
              view.points[-1].monitors == {"B(O 0,H 1)": 0.99, "B(O 0,H 2)": 0.94})
        check("a second build with nothing changed is the same object",
              reaction_path(job) is view)

        (d / "job_IRC_Full_trj.xyz").unlink()
        per_direction = reaction_path(job)
        check("without the full trajectory, the per-direction files and the input's TS",
              per_direction is not view and per_direction is not None
              and [p.atoms[0][3] for p in per_direction.points] == [0, 1, 2, 3, 4, 5])

        partial = d / "partial.xyz"
        partial.write_text(_irc_frame(0.0) + _irc_frame(1.0) + "3\ncomment\nO 0 0\n")
        check("a half-written frame is dropped, not read as atoms",
              len(read_xyz_frames(partial)) == 2)

        p = view.points[2]
        check("the panes take a path point's own geometry and label",
              geometry_shown(job, p) == (p.atoms, p) and describe_point(p) == "IRC TS"
              and geometry_shown(job, PathPoint("irc", 9, "IRC forward 8", None, 0.0, []))
              == ([], None))

    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        lines = _IRC_OUT.strip("\n").split("\n")
        end = next(i for i, line in enumerate(lines)
                   if line.startswith("    2     -100.004781")) + 1
        (d / "job.inp").write_text("! B97-3c IRC\n%irc MaxIter 3 end\n* xyzfile 0 1 start.xyz\n")
        (d / "job.out").write_text("\n".join(lines[:end]) + "\n")
        job = Job(d, "job", d.parent, label="irc")
        job.refresh(gone)
        running = reaction_path(job)
        check("an unfinished IRC follows its newest point, not the TS",
              running is not None and len(running.points) == 4 and running.focus == 3)


def test_the_irc_is_summarised_and_flagged():
    print("\nan IRC says where it is, and a direction that hit MaxIter is flagged")
    from orcamon.core.paths import path_summary
    from orcamon.core.status import FLAG_CODES

    state = JobState(path=Path("/nonexistent"), stem="job")
    _feed(state, _IRC_OUT)

    summary = ("IRC        backward 2 points, end -6.3 kJ/mol · "
               "forward 3 points (MaxIter), end -12.6 kJ/mol (from the TS)", "irc B2/F3")
    check("the summary names each direction, its size and its end",
          path_summary(state) == summary)

    inp = "! B97-3c IRC\n%irc MaxIter 3 end\n* xyzfile 0 1 start.xyz\n"
    report = build_report(_report_job(state, inp, Liveness(False, "process")))
    check("`show` prints the summary line", summary[0] in render_plain(report))
    check("a direction that hit MaxIter raises irc_not_converged",
          {"code": "irc_not_converged", "message": "IRC forward hit MaxIter (3) before a minimum"}
          in report.attention, f"{report.attention}")

    # ORCA's own MaxIter banner is the only evidence a direction stopped
    # short; without it a finished IRC raises nothing.
    lines = _IRC_OUT.strip("\n").split("\n")
    start = next(i for i, line in enumerate(lines)
                 if "MAXIMUM NUMBER OF ITERATIONS REACHED" in line) - 1
    clean = JobState(path=Path("/nonexistent"), stem="job")
    _feed(clean, "\n".join(lines[:start] + lines[start + 3:]))
    clean_report = build_report(_report_job(clean, inp, Liveness(False, "process")))
    check("without the banner the direction is not called unconverged",
          "irc_not_converged" not in [f["code"] for f in clean_report.attention],
          f"{clean_report.attention}")

    check("and the code sits beside scan_incomplete in the contract",
          FLAG_CODES.index("irc_not_converged") == FLAG_CODES.index("scan_incomplete") + 1)


_NEB_OUT = """
Current trajectory will be written to    ....  job_MEP_trj.xyz

Starting iterations:

Optim.  Iteration  HEI  E(HEI)-E(0)  max(|Fp|)   RMS(Fp)    dS
Switch-on CI threshold               0.020000
   LBFGS     0      2    0.010000    0.050000   0.005000  2.0000
   LBFGS     1      2    0.009000    0.030000   0.004000  1.9000

Image  2 will be converted to a climbing image in the next iteration (max(|Fp|) < 0.0200)

Optim.  Iteration  CI   E(CI)-E(0)   max(|Fp|)   RMS(Fp)    dS     max(|FCI|)   RMS(FCI)
Convergence thresholds               0.020000   0.010000            0.002000    0.001000
   LBFGS     2      2    0.008000    0.015000   0.003000  1.8000    0.010000    0.002000
"""

_NEB_ENERGIES = ("-100.000000000000", "-99.995000000000",
                 "-99.992000000000", "-99.998000000000")


def _neb_frame(z, energy):
    """One frame of the synthetic NEB tree, identified by the O atom's z and
    carrying the image's energy in the comment."""
    lines = ["3", f"Coordinates from ORCA-job job_MEP E {energy}"]
    for el, x, y in (("O", 0.0, 0.0), ("H", 0.96, 0.0), ("H", -0.24, 0.93)):
        lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines) + "\n"


def _neb_tree(d: Path):
    """§ Facts' synthetic NEB tree: a running NEB-CI whose current MEP
    trajectory holds four images, the O atom's z identifying each."""
    d.mkdir(parents=True, exist_ok=True)
    (d / "job.inp").write_text(
        '! B97-3c NEB-CI\n%neb NEB_End_XYZFile "product.xyz" end\n* xyzfile 0 1 start.xyz\n')
    (d / "job.out").write_text(_NEB_OUT.strip("\n") + "\n")
    (d / "start.xyz").write_text(_neb_frame(0.0, _NEB_ENERGIES[0]))
    (d / "job_MEP_trj.xyz").write_text(
        "".join(_neb_frame(z, e) for z, e in enumerate(_NEB_ENERGIES)))


def test_the_neb_is_parsed_and_built():
    print("\nan NEB's iteration rows and its current path are read")
    import tempfile

    from orcamon.core.paths import path_summary, reaction_path

    state = JobState(path=Path("/nonexistent"), stem="job")
    _feed(state, _NEB_OUT)

    check("the three iteration rows, the last with its image, barrier and forces",
          len(state.neb_rows) == 3
          and state.neb_rows[-1] == NebRow(2, "CI", 2, 0.008, 0.015, 0.003),
          f"{state.neb_rows}")
    check("each row carries the phase its block's header named",
          [r.phase for r in state.neb_rows] == ["HEI", "HEI", "CI"])
    check("and the MEP trajectory the log names",
          state.neb_file == "job_MEP_trj.xyz", f"{state.neb_file}")

    summary = ("NEB        iteration 2 · climbing image 2 · E(CI)-E(0) +21.0 kJ/mol",
               "neb 2 CI2")

    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        _neb_tree(d)
        job = Job(d, "job", d.parent, label="neb")
        job.refresh(Liveness(False, "process"))
        view = reaction_path(job)

        check("the current path is an NEB, one point per image, CI on the latest row's",
              view is not None and view.kind == "neb"
              and [p.index for p in view.points] == [0, 1, 2, 3]
              and [p.label for p in view.points] == [
                  "NEB image 0", "NEB image 1", "NEB image 2 (CI)", "NEB image 3"])
        check("its dE is measured from image 0",
              view is not None
              and [round(p.de_kj_mol, 4) for p in view.points]
              == [0.0, 13.1275, 21.004, 5.251])
        check("the image the latest row names is the focus",
              view is not None and view.focus == 2)
        check("the summary says which image and what barrier it has",
              path_summary(state) == summary)

        trj = d / "job_MEP_trj.xyz"
        trj.write_text("".join(_neb_frame(z, e) for z, e in enumerate(_NEB_ENERGIES[:3]))
                       + "3\nCoordinates from ORCA-job job_MEP E -99.998000000000\n"
                         "O 0.000000 0.000000 3.000000\n")
        cut = reaction_path(job)
        check("a half-written frame ends the path there",
              cut is not None and [p.index for p in cut.points] == [0, 1, 2])

        trj.unlink()
        check("without the MEP file there is no path, but the summary still says where it is",
              reaction_path(job) is None and path_summary(state) == summary)


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


def test_orca_says_why_it_converged():
    print("\nwhen ORCA converges on its relaxed rule, the report names which one")
    B = "         *              GEOMETRY OPTIMIZATION CYCLE   5            *"
    B6 = "         *              GEOMETRY OPTIMIZATION CYCLE   6            *"
    H = "      ***        THE OPTIMIZATION HAS CONVERGED     ***"

    first = JobState(path=Path("/nonexistent"), stem="job")
    _feed(first, B)
    _feed(first, "       The step convergence is overachieved with ")
    _feed(first, H)
    check("the step overachievement is named",
          first.opt_converged_reason == "step overachieved", f"{first.opt_converged_reason}")

    state = JobState(path=Path("/nonexistent"), stem="job")
    _feed(state, B)
    _feed(state, "       Everything but the energy has converged. However, the energy")
    _feed(state, H)
    check("the energy-only convergence is named",
          state.opt_converged_reason == "energy nearly converged", f"{state.opt_converged_reason}")

    state = JobState(path=Path("/nonexistent"), stem="job")
    _feed(state, B)
    _feed(state, "       The gradient convergence is overachieved with ")
    _feed(state, H)
    check("the gradient overachievement is named",
          state.opt_converged_reason == "gradient overachieved", f"{state.opt_converged_reason}")

    state = JobState(path=Path("/nonexistent"), stem="job")
    _feed(state, B)
    _feed(state, "       The step convergence is overachieved with ")
    _feed(state, B6)
    _feed(state, H)
    check("a reason from an earlier cycle is not carried over",
          state.opt_converged is True and state.opt_converged_reason is None,
          f"{state.opt_converged} {state.opt_converged_reason}")

    report = build_report(_report_job(first, "! B97-3c Opt\n* xyz 0 1\n*\n",
                                      Liveness(False, "process")))
    check("the text report names the rule",
          "optimization converged (on ORCA's relaxed rule: step overachieved)" in render_plain(report))
    check("and the JSON contract does not change", not any("reason" in k for k in report.to_dict()))


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


def test_a_job_that_prints_no_geometry_shows_the_one_it_wrote():
    print("\na job that prints no coordinates shows the file it wrote, or its input, and says which")
    import os
    import tempfile
    from orcamon.core.report import describe_point

    def xyz(path, atoms):
        path.write_text(f"{len(atoms)}\ncomment\n" + "".join(f"{e} {x} 0.0 0.0\n" for e, x in atoms))

    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        (d / "job.inp").write_text("! XTB ALPB(Water) SOLVATOR\n* xyzfile 0 1 start.xyz\n")
        xyz(d / "start.xyz", [("O", 0.0), ("H", 0.96), ("H", -0.24)])
        xyz(d / "job.solvator.xyz", [("O", 0.0), ("H", 0.96), ("H", -0.24), ("O", 3.0), ("H", 3.9), ("H", 2.7)])
        # What SOLVATOR prints: no coordinate block, one line naming the file.
        (d / "job.out").write_text("                       * ORCA Solvator *\n"
                                   "Final structured saved to        :             job.solvator.xyz\n"
                                   "		****ORCA-SOLVATOR TERMINATED NORMALLY****\n"
                                   f"{_TERMINATED}\n")
        job = Job(d, "job", d.parent)
        gone = Liveness(False, "process")
        job.refresh(gone)
        atoms, shown = geometry_shown(job, None)
        check("the file the log names is drawn", len(atoms) == 6 and describe_point(shown) == "job.solvator.xyz",
              f"{len(atoms)} {describe_point(shown)!r}")
        check("and counted in the summary", build_report(job).n_atoms == 6)
        first = job.file_geometry
        job.refresh(gone)
        check("an unchanged file is the same object, so the pixel panes do not redraw",
              job.file_geometry is first)

        os.remove(d / "job.solvator.xyz")
        job.refresh(gone)
        atoms, shown = geometry_shown(job, None)
        check("without it, the input's *xyzfile, named as the input",
              len(atoms) == 3 and describe_point(shown) == "input geometry, start.xyz", describe_point(shown))

        (d / "job.inp").write_text("! XTB SP\n* xyz 0 1\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\n*\n")
        job.refresh(gone)
        atoms, shown = geometry_shown(job, None)
        check("an inline * xyz block is the input geometry too",
              len(atoms) == 2 and describe_point(shown) == "input geometry", describe_point(shown))

        with open(d / "job.out", "a") as f:
            f.write("CARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------\n"
                    "  O      0.000000    0.000000    0.000000\n\n")
        job.refresh(gone)
        atoms, shown = geometry_shown(job, None)
        check("a printed geometry always wins", len(atoms) == 1 and job.file_geometry is None
              and not describe_point(shown).startswith("input"), f"{len(atoms)} {describe_point(shown)!r}")


if __name__ == "__main__":
    test_a_multilayer_energy_is_the_total()
    test_a_scan_is_keyed_by_step_and_cycle()
    test_a_missing_geometry_is_labelled_not_borrowed_silently()
    test_every_imaginary_frequency_is_kept()
    test_the_steps_pane_shows_the_latest_cycle()
    test_the_input_names_charge_and_multiplicity()
    test_the_basis_is_the_one_orca_reports()
    test_the_fast_path_reads_what_the_line_path_reads()
    test_normal_modes_are_parsed()
    test_normal_modes_keep_only_imaginary()
    test_normal_modes_take_the_last_block()
    test_every_mode_of_a_final_block_is_kept()
    test_modes_are_found_and_described()
    test_mode_offsets_fill_the_whole_structure()
    test_mode_offsets_map_onto_the_qm_subset()
    test_mode_offsets_refuse_a_mismatched_mode()
    test_the_sine_phase_displaces_atoms()
    test_the_mode_geometry_falls_back_to_an_alternate()
    test_the_qm_map_needs_indices_that_fit()
    test_a_frame_is_a_small_faithful_palette_png()
    test_near_atoms_hide_far_ones()
    test_spheres_are_shaded()
    test_fog_dims_the_far_side()
    test_outlines_separate_overlapping_atoms()
    test_representations_change_what_is_drawn()
    test_space_filling_is_framed()
    test_principal_axes_face_a_plane_on()
    test_labels_follow_occlusion()
    test_wireframe_keeps_labels_and_distances()
    test_wireframe_draws_the_qm_layer_thicker()
    test_fog_holds_still_under_rotation()
    test_host_fog_keeps_near_host_bright()
    test_halos_widen_with_the_depth_gap()
    test_boundary_host_atoms_get_a_ball()
    test_see_through_shows_the_guest_behind_the_host()
    test_see_through_ghosts_the_environment()
    test_the_camera_tumbles_over_the_pole()
    test_the_renderer_is_fast_enough()
    test_the_view_is_not_mirrored()
    test_text_mode_draws_the_right_enantiomer()
    test_every_marker_reaches_its_parser()
    test_irc_rows_are_parsed()
    test_the_irc_path_is_built()
    test_the_irc_is_summarised_and_flagged()
    test_the_neb_is_parsed_and_built()
    test_status_names_the_outcome()
    test_attention_flags()
    test_orca_says_why_it_converged()
    test_plain_and_markup_say_the_same_thing()
    test_report_keys_are_stable()
    test_a_job_that_prints_no_geometry_shows_the_one_it_wrote()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
