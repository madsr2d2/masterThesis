"""
The contract `computational/orca_io.py` holds ORCA output to.

    python test_orca_io.py

TWO PARSERS READ THIS PROJECT'S ORCA JOBS UNTIL 2026-09-17 AND BOTH WERE
WRONG, in opposite directions, silently.

`~/.claude/skills/lib/orca_thermo.py` matched `FINAL SINGLE POINT ENERGY` with
a regex that cannot match a LABELLED line, so on every QM/XTB job -- the tier
the whole C7-C10 register is built on -- it returned the r2SCAN-3c QM1 region
alone (-418.878772917932 on the C8 reaction complex) while reporting, in the
same dict, a Gibbs free energy built on the QM/QM2 total (-658.19992068). OPI
returns -659.180810832202 for that file.

And on a job run at more than one temperature the two disagree about WHICH
temperature they answered for. Measured on the fixture below, a real ORCA
6.1.1 water job at 288.15/298.15/313.15 K:

    the regex parser          -75.95690128   (the LAST block, 313.15 K)
    OPI's get_free_energy()   -75.95510681   (entry [0], 288.15 K)

Neither says which it gave. `orca_io.free_energy` demands a temperature and
raises for one the job did not compute, and the first test below is that
regression: C7 runs exactly these three temperatures.

THE FIXTURE IS A REAL JOB, NOT A SYNTHETIC ONE, and it is 43 KB of ORCA's own
structured property JSON plus the log. OPI parses it with no ORCA binaries
present -- checked below -- so this gate runs on a fresh clone where
`orca_2json` does not exist. The heavy QM/XTB regressions need jobs that stay
homelab-local by policy (`.gitignore`), so they run WHEN PRESENT and say so
when not, rather than being left unwritten.

`run_gates.py` discovers this file by its name, like any other gate.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "computational"))
sys.path.insert(0, HERE)

import orca_io
from doc_check import Checker

FAILURES = []

# The tracked three-temperature job. `%freq Temp 288.15, 298.15, 313.15`.
FIXTURE = os.path.join(HERE, "computational", "fixtures", "water_freq_3T")
FIXTURE_TEMPERATURES = (288.15, 298.15, 313.15)

# The free energies ORCA itself printed for that job, one per temperature,
# read off its own `Final Gibbs free energy` lines.
FIXTURE_FREE_ENERGIES = {288.15: -75.95510681,
                         298.15: -75.95582142,
                         313.15: -75.95690128}

# The C8 reaction complex: QM/XTB, 140 atoms, a converged MINIMUM that carries
# two imaginary modes at -19.76 and -3.55 cm-1 -- the case the cutoff exists
# for. Homelab-local; skipped where absent.
C8_RC = os.path.join(HERE, "computational", "C8_perhydrate_trap",
                     "K+H2O2_water-relay_to_KP", "geometry_r2scan3c-xtb", "rc")
C8_RC_TOTAL_ENERGY_EH = -659.180810832202
C8_RC_QM1_ONLY_EH = -418.878772917932

# The skill that tells an agent to use `orca_io` at all. CLAUDE.md's rule --
# "a number that reaches CLAUDE.md, BUBBLES.md or the skill now needs a claim
# in it" -- applies to a second skill as much as the first, and
# `test_root_documents.py` cannot see this one: its SKILL path is hardcoded
# to `analyse-kinetics`. So this file owns run-orca's numbers, which is the
# right place anyway: the gate that DERIVES a number is the gate that should
# check the document quoting it.
RUN_ORCA_SKILL = os.path.join(HERE, ".claude", "skills", "run-orca", "SKILL.md")


def check(label, ok, detail=""):
    print(f"  {'pass' if ok else 'FAIL'}  {label}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label)


def skip(label, why):
    print(f"  skip  {label}: {why}")


def test_a_temperature_is_never_guessed():
    """Every temperature is readable, and the right one comes back."""
    print("\na temperature is never guessed")
    job = orca_io.load_job(FIXTURE)
    found = orca_io.temperatures(job)
    check("the job carries all three temperatures",
          [round(t, 2) for t in found] == list(FIXTURE_TEMPERATURES), f"{found}")
    for temperature, expected in FIXTURE_FREE_ENERGIES.items():
        got = orca_io.free_energy(job, temperature)
        check(f"G at {temperature} K is the block ORCA printed for it",
              abs(got - expected) < 1e-7, f"{got:.8f} against {expected:.8f}")

    # The trap: OPI's own accessor answers for entry [0] whatever you meant.
    convenience = job.output.get_free_energy()
    check("OPI's own get_free_energy() still answers for entry [0]",
          abs(convenience - FIXTURE_FREE_ENERGIES[288.15]) < 1e-7,
          f"{convenience:.8f} -- which is why nothing here calls it")
    check("and it differs from the 298.15 K value a caller would expect",
          abs(convenience - FIXTURE_FREE_ENERGIES[298.15]) > 1e-6,
          "the two temperatures would be indistinguishable otherwise")


def test_the_entropy_is_not_the_entropy_term():
    """ORCA's `entropys` is T*S, and reading it as S is out by a factor of T."""
    print("\nthe entropy is not the entropy term")
    job = orca_io.load_job(FIXTURE)
    for temperature in FIXTURE_TEMPERATURES:
        term = orca_io.entropy_term(job, temperature)
        gap = orca_io.enthalpy(job, temperature) - orca_io.free_energy(job, temperature)
        check(f"T*S closes G = H - T*S at {temperature} K",
              abs(term - gap) < 1e-9, f"{term:.9f} against {gap:.9f}")
        check(f"and S is that divided by T at {temperature} K",
              abs(orca_io.entropy(job, temperature) - term / temperature) < 1e-15)


def test_an_absent_temperature_raises():
    """Asking for a temperature the job did not compute is an error."""
    print("\nan absent temperature raises")
    job = orca_io.load_job(FIXTURE)
    try:
        orca_io.free_energy(job, 273.15)
    except orca_io.OrcaJobError as error:
        message = str(error)
        check("it raises rather than returning a neighbouring temperature", True)
        check("and the message names what the job actually has",
              "288.15" in message and "313.15" in message, message)
    else:
        check("it raises rather than returning a neighbouring temperature", False,
              "returned a number for a temperature the job never computed")


def test_a_missing_output_is_not_a_failed_job():
    """OPI reads a missing `.out` as False, not as absent."""
    print("\na missing output is not a failed job")
    try:
        orca_io.load_job(os.path.join(HERE, "computational", "fixtures"), "no_such_job")
    except orca_io.OrcaJobError as error:
        check("load_job raises for a path with no output file", True)
        check("and says so, rather than reporting a crash",
              "no output file" in str(error), str(error))
    else:
        check("load_job raises for a path with no output file", False,
              "a mistyped path read as a job that crashed")


def test_the_fixture_parses_without_orca():
    """The gate must run on a clone with no ORCA installed."""
    print("\nthe fixture parses without orca")
    present = os.path.exists(os.path.join(FIXTURE, "job.property.json"))
    check("the structured property JSON is tracked beside the log",
          present, "OPI would otherwise have to shell out to orca_2json")
    job = orca_io.load_job(FIXTURE)
    check("and ORCA's own version comes off it",
          job.version == orca_io.ORCA_VERSION, f"{job.version}")
    check("with the job's atom count",
          job.natoms == 3, f"{job.natoms}")
    check("and its cost, which every log entry has to report",
          job.cpu_seconds is not None and job.cpu_seconds > 0,
          f"{job.cpu_seconds}")


def test_the_standard_states_move_with_temperature():
    """+7.93 and -9.94 kJ/mol are the 298.15 K values, not constants."""
    print("\nthe standard states move with temperature")
    at298 = orca_io.gas_to_molar_correction(298.15)
    check("1 atm -> 1 M is +7.93 kJ/mol at 298.15 K",
          abs(at298 - 7.93) < 0.01, f"{at298:.3f}")
    low = orca_io.gas_to_molar_correction(288.15)
    high = orca_io.gas_to_molar_correction(313.15)
    # To two decimals, not to two-and-a-bit: a +/- 0.02 window passes a value
    # that rounds the other way, which is exactly the slip the skill's own
    # "9.61" was until the document claim below caught it.
    check("but 7.58 at 288.15 K and 8.45 at 313.15 K",
          f"{low:.2f}" == "7.58" and f"{high:.2f}" == "8.45",
          f"{low:.4f} and {high:.4f}")
    check("so a single constant is over a third of a kJ/mol wrong at C7's ends",
          (at298 - low) > 0.3 and (high - at298) > 0.3,
          f"{at298 - low:.3f} and {high - at298:.3f}")

    water = orca_io.solvent_standard_state_correction(298.15)
    check("water at 55.34 M is 9.95 kJ/mol at 298.15 K",
          f"{water:.2f}" == "9.95", f"{water:.4f}")
    check("and it moves too: 9.62 at 288.15 K, 10.45 at 313.15 K",
          f"{orca_io.solvent_standard_state_correction(288.15):.2f}" == "9.62"
          and f"{orca_io.solvent_standard_state_correction(313.15):.2f}" == "10.45")


def test_an_equilibrium_becomes_a_fraction():
    """C7's gate is a hydrate FRACTION, which needs both steps."""
    print("\nan equilibrium becomes a fraction")
    check("dG = 0 is K = 1",
          abs(orca_io.equilibrium_constant(0.0, 298.15) - 1.0) < 1e-12)
    check("and a half-and-half mixture",
          abs(orca_io.fraction_from_equilibrium_constant(1.0) - 0.5) < 1e-12)
    # C7's own bar: a hydrate fraction of 0.8 is K = 4.
    bound = orca_io.fraction_from_equilibrium_constant(4.0)
    check("K = 4 is C7's 0.8 bar", abs(bound - 0.8) < 1e-12, f"{bound}")
    downhill = orca_io.equilibrium_constant(-10.0, 298.15)
    check("and a favourable dG gives K > 1", downhill > 1.0, f"{downhill:.2f}")


def test_a_fragment_pair_must_be_one_fragment():
    """The ONIOM correction's two jobs have to be the same system."""
    print("\na fragment pair must be one fragment")
    job = orca_io.load_job(FIXTURE)
    same = orca_io.oniom_correction_eh(job, job)
    check("the same job against itself corrects by nothing",
          abs(same) < 1e-12, f"{same}")

    class Mismatched:
        label, natoms, charge = "synthetic", 11, -1
        final_energy_eh = -1.0
    try:
        orca_io.oniom_correction_eh(Mismatched(), job)
    except orca_io.OrcaJobError as error:
        check("and a pair of different sizes raises", True)
        check("naming both atom counts", "11" in str(error) and "3" in str(error),
              str(error))
    else:
        check("and a pair of different sizes raises", False,
              "two different systems were subtracted from each other")


def test_the_water_job_is_a_minimum():
    """A converged minimum with no imaginary mode reads as one."""
    print("\nthe water job is a minimum")
    job = orca_io.load_job(FIXTURE)
    modes = orca_io.frequencies(job, 298.15)
    check("every mode is present, not a truncated window",
          len(modes) == 9, f"{len(modes)} of 9 for three atoms")
    check("none of them is imaginary past the cutoff",
          orca_io.imaginary_modes(job) == [], f"{orca_io.imaginary_modes(job)}")
    check("so the Hessian says minimum",
          orca_io.stationary_point(job) == "minimum",
          orca_io.stationary_point(job))


def test_the_qm_xtb_total_is_the_energy_returned():
    """The regression the regex parser failed, where the job is present."""
    print("\nthe qm/xtb total is the energy returned")
    if not os.path.exists(os.path.join(C8_RC, "job.out")):
        skip("the C8 reaction complex",
             "homelab-local by .gitignore policy; run this on the machine that has it")
        return
    job = orca_io.load_job(C8_RC)
    energy = job.final_energy_eh
    check("the energy is the QM/QM2 total",
          abs(energy - C8_RC_TOTAL_ENERGY_EH) < 1e-6, f"{energy}")
    check("and NOT the r2SCAN-3c QM1 region the regex parser returned",
          abs(energy - C8_RC_QM1_ONLY_EH) > 100.0,
          f"{abs(energy - C8_RC_QM1_ONLY_EH):.3f} Eh apart")
    check("it has 140 atoms", job.natoms == 140, f"{job.natoms}")
    verdict = orca_io.stationary_point(job)
    check("and it is a MINIMUM despite two imaginary modes under the cutoff",
          verdict.startswith("minimum"), verdict)
    check("which the verdict names rather than hiding",
          "-19.76" in verdict and "-3.55" in verdict, verdict)


def test_the_skill_quotes_what_the_code_says():
    """Every number in the run-orca skill, re-derived.

    One number is quoted there and NOT re-derived from a live job:
    `-418.878772917932`, what the superseded regex parser returned for the C8
    reaction complex. It is the bare `FINAL SINGLE POINT ENERGY` line, which
    includes r2SCAN-3c's gCP and dispersion terms, and OPI's `get_energies()`
    exposes only the SCF component (-418.88865521) -- so re-deriving it would
    mean regexing the log, which is the thing this module exists to stop. It
    is pinned instead by `test_the_qm_xtb_total_is_the_energy_returned`, which
    asserts against the live job that the true total sits 240 Eh away from it.
    Document -> constant -> live job is a real chain; reading the number back
    out of the line it was written from would not be.
    """
    print("\nthe skill quotes what the code says")
    doc = Checker(RUN_ORCA_SKILL, label="the run-orca skill",
                  document_label="SKILL.md")
    doc.section("The run-orca skill")

    doc.claim("the QM/XTB total OPI returns", f"{C8_RC_TOTAL_ENERGY_EH}")
    doc.claim("against the QM1 region the regex parser returned",
              f"{C8_RC_QM1_ONLY_EH}")

    doc.claim("the imaginary-mode cutoff",
              f"{orca_io.IMAGINARY_CUTOFF_CM1:g} cm")
    job = orca_io.load_job(FIXTURE)
    if os.path.exists(os.path.join(C8_RC, "job.out")):
        small = orca_io.small_imaginary_modes(orca_io.load_job(C8_RC))
        for mode in small:
            doc.claim(f"the C8 minimum's ignored mode {mode:.2f}", f"{mode:.2f}")
    else:
        skip("the C8 reaction complex's ignored modes", "homelab-local")

    doc.claim("water's molarity", f"{orca_io.WATER_MOLARITY_M}")
    doc.claim("1 atm -> 1 M at 298.15 K",
              f"{orca_io.gas_to_molar_correction(298.15):.2f}")
    doc.claim("and the range over C7's temperatures",
              f"{orca_io.gas_to_molar_correction(288.15):.2f}"
              f"\u2013{orca_io.gas_to_molar_correction(313.15):.2f}")
    doc.claim("the water term's range",
              f"{orca_io.solvent_standard_state_correction(288.15):.2f}"
              f"\u2013{orca_io.solvent_standard_state_correction(313.15):.2f}")

    for temperature in orca_io.temperatures(job):
        doc.claim(f"the fixture's {temperature} K rung", f"{temperature}")

    doc.claim("and the version every log entry reports", orca_io.ORCA_VERSION)

    if doc.failures:
        FAILURES.extend(doc.failures)
    print(f"  ({doc.claims} claims)")


if __name__ == "__main__":
    test_a_temperature_is_never_guessed()
    test_the_entropy_is_not_the_entropy_term()
    test_an_absent_temperature_raises()
    test_a_missing_output_is_not_a_failed_job()
    test_the_fixture_parses_without_orca()
    test_the_standard_states_move_with_temperature()
    test_an_equilibrium_becomes_a_fraction()
    test_a_fragment_pair_must_be_one_fragment()
    test_the_water_job_is_a_minimum()
    test_the_qm_xtb_total_is_the_energy_returned()
    test_the_skill_quotes_what_the_code_says()
    print(f"\n{len(FAILURES)} failure(s)"
          + (": " + ", ".join(FAILURES) if FAILURES else ""))
    raise SystemExit(1 if FAILURES else 0)
