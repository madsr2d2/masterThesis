"""The parser on a real tree, and the check that every marker reaches it.

    python -m orcamon.validate [ROOT]

`MARKER_CASES` is the part with a verdict (exit 1 if any marker line is not
read); the rest prints what the parser made of every job under ROOT, for a
person to compare against the outputs themselves."""
from __future__ import annotations

import sys
from pathlib import Path

from .core.discovery import discover
from .core.job import Job
from .core.liveness import lookup, make_probe, running_orca_cwds
from .core.parser import JobState


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    root = Path(argv[0] if argv else ".").resolve()
    running_cwds = running_orca_cwds()
    probe = make_probe("auto")
    jobs = discover(root)
    snapshot = probe.snapshot([ref.path for ref in jobs])
    print(f"discover found {len(jobs)} jobs under {root}\n")
    for ref in jobs:
        job = Job(ref.path, ref.stem, root, label=ref.label)
        job.refresh(lookup(snapshot, probe, ref.path, ref.stem))
        state = job.state
        print(f"--- {ref.label} ({ref.stem}) ---")
        print(f"  status: {job.status.value}  (liveness: {job.liveness.source})")
        for flag in job.flags:
            print(f"  ! {flag.code}: {flag.message}")
        print(f"  bytes read: {state.offset}  cycle: {state.cycle}  "
              f"convergence steps: {len(state.convergence_history)}  "
              f"eigenvalue readings: {len(state.eigen_history)}")
        print(f"  final_energy: {state.final_energy} ({state.final_energy_label})  "
              f"n_imaginary: {state.n_imaginary}  wall_time_s: {state.wall_time_s}")
        print(f"  normal_completion: {state.normal_completion}  opt_converged: {state.opt_converged}  "
              f"opt_maxiter_reached: {state.opt_maxiter_reached}  "
              f"qm2_errors: {state.qm2_error_count}  possibly_stalled: {state.possibly_stalled()}")
        print()

    print("running orca* cwds detected on this machine:")
    for cwd, started in running_cwds.items():
        print(f"  {cwd}  (started {started:.0f})")

    print("\n--- synthetic stall-heuristic check ---")
    print("plateaued RMS gradient over 20 cycles:", feed_synthetic_cycles(plateaued=True).possibly_stalled())
    print("steadily improving RMS gradient over 20 cycles:", feed_synthetic_cycles(plateaued=False).possibly_stalled())

    print("\n--- prefilter coverage ---")
    return 1 if check_markers() else 0


# One line per marker `feed_line` reads, and the field it must land in. The
# parser runs its patterns behind literal prefilters (`_RARE_MARKERS_RE`,
# `_CONV_HINT_RE`), so a pattern whose literal is missing from its prefilter
# would simply stop matching -- silently, and only on real output. This is the
# check that says so: it drives the whole path, prefilter included, rather
# than inspecting the regexes.
MARKER_CASES = [
    ("         *              GEOMETRY OPTIMIZATION CYCLE   7            *",
     lambda s: s.cycle == 7),
    ("        Hessian has     1 negative eigenvalue",
     lambda s: list(s.eigen_history)[-1][1] == 1),
    ("There was an error in the QM2 calculation",
     lambda s: s.qm2_error_count == 1),
    ("ORCA finished by error termination in SCF",
     lambda s: s.crashed_marker),
    ("Aborting the run", lambda s: s.crashed_marker),
    ("TERMINATING THE RUN", lambda s: s.crashed_marker),
    ("                 ****ORCA TERMINATED NORMALLY****",
     lambda s: s.normal_completion),
    ("      ***        THE OPTIMIZATION HAS CONVERGED     ***",
     lambda s: s.opt_converged),
    ("Final structured saved to        :             job.solvator.xyz",
     lambda s: s.result_geometry_file == "job.solvator.xyz"),
    ("       The optimization did not converge but reached the maximum ",
     lambda s: s.opt_maxiter_reached),
    ("       The optimization did not converge but reached the maximum number of",
     lambda s: s.opt_maxiter_reached),
    ("                    *** OPTIMIZATION RUN DONE ***",
     lambda s: s.opt_converged),
    ("FINAL SINGLE POINT ENERGY      -647.909099472771",
     lambda s: s.final_energy == -647.909099472771),
    ("FINAL SINGLE POINT ENERGY (QM/QM2)     -873.756647314441",
     lambda s: s.final_energy == -873.756647314441 and s.final_energy_label == "QM/QM2"),
    ("         *               RELAXED SURFACE SCAN STEP   3               *",
     lambda s: s.scan_step == 3),
    ("There will be   15 constrained geometry optimizations.",
     lambda s: s.scan_total == 15),
    ("There is 1 parameter to be scanned.",
     lambda s: s.scan_params == 1),
    ("Max. no of cycles        MaxIter  .... 200",
     lambda s: s.max_cycles == 200),
    ("----- Orbital basis set information -----",
     lambda s: s._in_orbital_basis),
    ("TOTAL RUN TIME: 0 days 1 hours 2 minutes 3 seconds 456 msec",
     lambda s: s.wall_time_s == 3723.456),
    ("QM1 Subsystem      ...  128 129 130",
     lambda s: s.qm_atom_indices == {128, 129, 130}),
    ("          Energy change      -0.0000038757            0.0000050000      NO",
     lambda s: s._pending_step.get("energy_change") == (-3.8757e-06, 5e-06, False)),
    ("          RMS gradient        0.0004500000            0.0001000000      NO",
     lambda s: s._pending_step.get("rms_grad") == (0.00045, 0.0001, False)),
    ("          MAX gradient        0.0025000000            0.0003000000      NO",
     lambda s: s._pending_step.get("max_grad") == (0.0025, 0.0003, False)),
    ("          RMS step            0.0020825435            0.0020000000      NO",
     lambda s: s._pending_step.get("rms_step") == (0.0020825435, 0.002, False)),
    ("          MAX step            0.0192646101            0.0040000000      YES",
     lambda s: s._pending_step.get("max_step") == (0.0192646101, 0.004, True)),
    ("NORMAL MODES", lambda s: s._in_modes_block),
    ("       The gradient convergence is overachieved with ",
     lambda s: s._pending_converged_reason == "gradient overachieved"),
    ("       The step convergence is overachieved with ",
     lambda s: s._pending_converged_reason == "step overachieved"),
    ("       Everything but the energy has converged. However, the energy",
     lambda s: s._pending_converged_reason == "energy nearly converged"),
    ("         *                          FORWARD IRC                      *",
     lambda s: s.irc_direction == "forward"),
    ("Iteration    E(Eh)      dE(kcal/mol)  max(|G|)   RMS(G)  B(O 0,H 1)",
     lambda s: s._in_irc_rows and s.irc_monitors == ["B(O 0,H 1)"]),
    ("         *  MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC RUN  *",
     lambda s: len(s.irc_maxiter) == 1),
    ("Storing forward trajectory in       .... job_IRC_F_trj.xyz",
     lambda s: s.irc_files == {"forward": "job_IRC_F_trj.xyz"}),
    ("Optim.  Iteration  CI   E(CI)-E(0)   max(|Fp|)   RMS(Fp)    dS",
     lambda s: s._in_neb_rows and s._neb_phase == "CI"),
    ("Current trajectory will be written to    ....  job_MEP_trj.xyz",
     lambda s: s.neb_file == "job_MEP_trj.xyz"),
]


def check_markers() -> int:
    # Both paths: `feed_line` behind the per-line prefilters, and `feed_text`,
    # whose whole-chunk scan is one more prefilter that must admit the line.
    failures = 0
    for line, holds in MARKER_CASES:
        for path, feed in (("feed_line", lambda s: s.feed_line(line)),
                           ("feed_text", lambda s: s.feed_text(line + "\n"))):
            state = JobState(path=Path("/nonexistent"), stem="job")
            feed(state)
            if not holds(state):
                failures += 1
                print(f"  FAIL not read past the prefilter by {path}: {line.strip()!r}")
    label = "all read" if not failures else f"{failures} NOT READ"
    print(f"marker lines reaching their parser ({len(MARKER_CASES)} cases): {label}")
    return failures


def feed_synthetic_cycles(plateaued: bool) -> JobState:
    state = JobState(path=Path("/nonexistent"), stem="job")
    for cycle in range(1, 21):
        rms_grad = 4.5e-4 if plateaued else 1.0e-2 / cycle
        state.feed_line(f"         *                GEOMETRY OPTIMIZATION CYCLE  {cycle}            *")
        state.feed_line(f"          Energy change      -0.0000038757            0.0000050000      NO")
        state.feed_line(f"          RMS gradient        {rms_grad:.10f}            0.0001000000      NO")
        state.feed_line(f"          MAX gradient        0.0025000000            0.0003000000      NO")
        state.feed_line(f"          RMS step            0.0020825435            0.0020000000      NO")
        state.feed_line(f"          MAX step            0.0192646101            0.0040000000      NO")
    return state


if __name__ == "__main__":
    sys.exit(main())
