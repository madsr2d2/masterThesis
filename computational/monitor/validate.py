from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from computational.monitor.discovery import discover_jobs
from computational.monitor.parser import JobState, new_state, update_job
from computational.monitor.procs import running_orca_cwds
from computational.monitor.status import compute_status

ROOT = Path(__file__).resolve().parents[2] / "computational"

CASES = [
    (ROOT / "hellowater", "hellowater"),
    (
        ROOT / "C8_perhydrate_trap/K+peroxide+BnOH_bridged_UNIDENTIFIED/energy_dlpno-ccsdt/r2scan3c",
        "job",
    ),
    (
        ROOT / "C8_perhydrate_trap/K+peroxide+BnOH_bridged_UNIDENTIFIED/energy_dlpno-ccsdt/dlpno-ccsdt",
        "job",
    ),
    (
        ROOT / "C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/rc",
        "job",
    ),
]


def main() -> None:
    running_cwds = running_orca_cwds()
    for job_dir, stem in CASES:
        state = new_state(job_dir, stem)
        update_job(state)
        is_running = job_dir.resolve() in running_cwds
        status = compute_status(state, is_running)
        print(f"--- {job_dir.relative_to(ROOT.parent)} ---")
        print(f"  status: {status.value}")
        print(f"  has_out: {state.has_out}")
        print(f"  bytes read (offset): {state.offset}")
        print(f"  cycle: {state.cycle}")
        print(f"  convergence steps recorded: {len(state.convergence_history)}")
        print(f"  eigenvalue readings: {len(state.eigen_history)}")
        print(f"  qm2_error_count: {state.qm2_error_count}")
        print(f"  final_energy: {state.final_energy}")
        print(f"  n_imaginary: {state.n_imaginary}")
        print(f"  wall_time_s: {state.wall_time_s}")
        print(f"  normal_completion: {state.normal_completion}")
        print(f"  opt_converged: {state.opt_converged}")
        if state.convergence_history:
            last = state.convergence_history[-1]
            print(f"  last convergence table (cycle {last.cycle}):")
            print(
                f"    energy_change={last.energy_change} (tol {last.energy_tol}, "
                f"conv={last.energy_conv})"
            )
            print(
                f"    rms_grad={last.rms_grad} (tol {last.rms_grad_tol}, "
                f"conv={last.rms_grad_conv})"
            )
            print(
                f"    max_grad={last.max_grad} (tol {last.max_grad_tol}, "
                f"conv={last.max_grad_conv})"
            )
        if state.eigen_history:
            print(f"  last eigenvalue readings: {list(state.eigen_history)[-5:]}")
        print(f"  possibly_stalled(): {state.possibly_stalled()}")
        print()

    print("running orca* cwds detected on this machine:")
    for cwd, ctime in running_orca_cwds().items():
        print(f"  {cwd}  (create_time={ctime})")

    jobs = discover_jobs(ROOT)
    print(f"\ndiscover_jobs found {len(jobs)} job directories under computational/")
    for job_dir, stem in jobs:
        print(f"  {job_dir.relative_to(ROOT)}  (stem={stem})")

    print("\n--- synthetic stall-heuristic check (no real plateaued job.out survives) ---")
    print("plateaued RMS gradient over 20 cycles:", feed_synthetic_cycles(plateaued=True).possibly_stalled())
    print("steadily improving RMS gradient over 20 cycles:", feed_synthetic_cycles(plateaued=False).possibly_stalled())

    print("\n--- prefilter coverage ---")
    if check_markers():
        sys.exit(1)


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
    ("                    *** OPTIMIZATION RUN DONE ***",
     lambda s: s.opt_converged),
    ("FINAL SINGLE POINT ENERGY      -418.878772917932",
     lambda s: s.final_energy == -418.878772917932),
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
]


def check_markers() -> int:
    failures = 0
    for line, holds in MARKER_CASES:
        state = JobState(path=Path("/nonexistent"), stem="job")
        state.feed_line(line)
        if not holds(state):
            failures += 1
            print(f"  FAIL not read past the prefilter: {line.strip()!r}")
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
    main()
