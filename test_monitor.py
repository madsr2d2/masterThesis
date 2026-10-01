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
    test_the_camera_tumbles_over_the_pole()
    test_the_renderer_is_fast_enough()
    test_the_view_is_not_mirrored()
    test_text_mode_draws_the_right_enantiomer()
    test_every_marker_reaches_its_parser()
    test_status_names_the_outcome()
    test_attention_flags()
    test_plain_and_markup_say_the_same_thing()
    test_report_keys_are_stable()
    test_a_job_that_prints_no_geometry_shows_the_one_it_wrote()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
