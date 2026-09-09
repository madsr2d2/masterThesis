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
