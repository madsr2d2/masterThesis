"""
The contract `computational/crest_io.py` holds CREST and xtb output to.

    python test_crest_io.py

THE FIXTURES ARE REAL RUNS (computational/fixtures/). `crest_frozen_host` is the
136-atom cyclodextrin-ketone + 2 waters seed of C8 sampled with the host frozen
(`-nci --gfnff --alpb water -mquick -mdlen x0.3` and `fixhost.inp`, force
constant 0.05); `xtb_water` is `xtb --ohess vtight --alpb water --json` on one
water.

Two things here are measured and would otherwise be folklore:

- A frozen host stays frozen: its heavy atoms sit 0.003 A from the seed in
  every conformer, where the same search unfrozen moved the host 1.4-2.2 A.
- Both seed waters stay inside the cavity in all 9 conformers, where 58 of 58
  conformers of the unfrozen run had none.
"""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "computational"))
import crest_io

HERE = os.path.dirname(os.path.abspath(__file__))
CREST = os.path.join(HERE, "computational", "fixtures", "crest_frozen_host")
XTB = os.path.join(HERE, "computational", "fixtures", "xtb_water")
FAILURES = []


def check(name, condition, detail=""):
    print(("ok   " if condition else "FAIL ") + name + (f"  {detail}" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


def raises(function, *args, **kwargs):
    try:
        function(*args, **kwargs)
    except Exception:
        return True
    return False


def truncated_copy(source, filename, drop):
    """A copy of a fixture directory whose `filename` has `drop` removed."""
    target = tempfile.mkdtemp()
    for name in os.listdir(source):
        shutil.copy(os.path.join(source, name), target)
    path = os.path.join(target, filename)
    with open(path) as handle:
        text = handle.read()
    with open(path, "w") as handle:
        handle.write(text.replace(drop, ""))
    return target


def test_crest_run():
    run = crest_io.load_run(CREST)
    check("crest: 9 conformers of 136 atoms", len(run.frames) == 9 and run.natoms == 136)
    check("crest: lowest first, relative energy 0", run.relative_kcal[0] == 0.0
          and run.relative_kcal == sorted(run.relative_kcal))
    check("crest: energies in Eh from the comment line",
          abs(run.energies_eh[0] - -21.19630059) < 1e-8)
    check("crest: population of the lowest", run.population_of_lowest == 40.285)
    check("crest: wall time read", run.wall_seconds is not None and 30 < run.wall_seconds < 60)


def test_crest_refuses_a_run_that_did_not_finish():
    check("crest: missing termination line raises",
          raises(crest_io.load_run, truncated_copy(CREST, "crest.out", crest_io.CREST_OK)))
    check("crest: missing directory raises", raises(crest_io.load_run, "no/such/run"))


def test_frozen_host_and_cavity():
    run = crest_io.load_run(CREST)
    _, seed, _ = crest_io.read_frames(os.path.join(CREST, "solvcluster.xyz"))
    host = [i for i in range(130) if run.elements[i] != "H"]
    worst = max(crest_io.host_rmsd(f, seed[0], host) for f in run.frames)
    check("host: frozen to 0.01 A in every conformer", worst < 0.01, f"worst {worst:.4f}")
    inside = [crest_io.cavity_occupancy(f, host, (130, 133)) for f in run.frames]
    check("cavity: both waters inside in all 9", inside == [2] * 9, str(inside))
    check("cavity: the seed itself holds both", crest_io.cavity_occupancy(seed[0], host, (130, 133)) == 2)
    rotated = seed[0] @ [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    check("host_rmsd is rotation-invariant", crest_io.host_rmsd(rotated + 3.0, seed[0], host) < 1e-9)
    mirrored = seed[0] * [1, 1, -1]
    check("host_rmsd does not superpose a mirror image",
          crest_io.host_rmsd(mirrored, seed[0], host) > 0.5)


def test_xtb_run():
    run = crest_io.load_xtb(XTB)
    check("xtb: total energy from the json", run.total_energy_eh == -5.08502128)
    check("xtb: free energy from the summary box", run.free_energy_eh == -5.082974170718)
    check("xtb: no imaginary modes", run.imaginary_modes == 0)
    check("xtb: gap", abs(run.gap_ev - 14.76747854) < 1e-8)
    check("xtb: optimised geometry", run.geometry[0] == ["O", "H", "H"] and run.geometry[1].shape == (3, 3))
    frequencies = crest_io.read_vibspectrum(os.path.join(XTB, "vibspectrum"))
    check("xtb: vibspectrum has 9 modes, 6 of them zero",
          len(frequencies) == 9 and sum(abs(f) < 1 for f in frequencies) == 6)
    check("xtb: bend and both stretches", [round(f) for f in frequencies[6:]] == [1444, 3595, 3597])
    check("xtb: missing success line raises",
          raises(crest_io.load_xtb, truncated_copy(XTB, "xtb.out", crest_io.XTB_OK)))
    check("xtb: a run with no json raises", raises(crest_io.load_xtb, tempfile.mkdtemp()))


if __name__ == "__main__":
    test_crest_run()
    test_crest_refuses_a_run_that_did_not_finish()
    test_frozen_host_and_cavity()
    test_xtb_run()
    if FAILURES:
        print(f"\n{len(FAILURES)} failed: {FAILURES}")
        sys.exit(1)
    print("\nall passed")
