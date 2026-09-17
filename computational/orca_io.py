"""
Every number this project reads out of an ORCA job, read through OPI.

THE OUTPUT FILE IS NOT PARSED WITH REGULAR EXPRESSIONS HERE, AND MUST NOT BE.
`~/.claude/skills/lib/orca_thermo.py` did that until 2026-09-17 and returned
the WRONG SYSTEM's energy on every QM/XTB job in this project -- the tier the
whole C7-C10 register is built on. ORCA closes a QM/XTB run with four
differently-labelled energies:

    FINAL SINGLE POINT ENERGY (L-QM2)     -266.683980845360
    FINAL SINGLE POINT ENERGY (S-QM2)      -26.381942931090
    FINAL SINGLE POINT ENERGY      -418.878772917932
    FINAL SINGLE POINT ENERGY (QM/QM2)     -659.180810832202

and `FINAL SINGLE POINT ENERGY\\s+(-?\\d+\\.\\d+)` cannot match a labelled
line, so it took the bare one -- the r2SCAN-3c QM1 region ALONE -- while
returning, in the same dict, a Gibbs free energy built on the QM/QM2 total.
The two numbers described different systems and nothing said so. OPI returns
-659.180810832202 for the same file, because it reads ORCA's own structured
property JSON rather than the human-readable log.

OPI (`orca-pi`, ORCA's own Python interface) shells out to `orca_2json` to
build `<basename>.property.json` and validates it through pydantic models, so
a field that moves between ORCA versions fails loudly instead of matching
nothing. `%output jsonpropfile true end` is what makes that JSON appear;
`build_input` adds it, and OPI's own `Calculator` does too.

WHAT THIS MODULE ADDS ON TOP OF OPI, and why it is not a thin pass-through:

1. `free_energy` DEMANDS A TEMPERATURE. OPI's own `get_free_energy()` reads
   `thermochemistry_energies[0]` -- its source comments "there should always
   be only one index" -- and on a multi-temperature job there are three.
   Measured on a three-temperature water job (288.15/298.15/313.15 K):
   `get_free_energy()` returned the 288.15 K value with no warning, and the
   regex parser it replaces returned the 313.15 K one. The two failed in
   OPPOSITE directions and both failed silently. C7 runs three temperatures,
   so there is no function here that returns a free energy without being told
   which one, and asking for a temperature the job did not compute raises.
2. THE STANDARD-STATE CORRECTIONS ARE FUNCTIONS OF TEMPERATURE, not the
   single constants they are usually quoted as. +7.93 kJ/mol (1 atm -> 1 M)
   and -9.95 kJ/mol (water at 55.34 M) hold at 298.15 K ONLY; over C7's own
   288.15-313.15 K range they run 7.58-8.45 and 9.62-10.45. Hardcoding the
   298.15 K value would put 0.5 kJ/mol of pure bookkeeping error into two of
   C7's three temperatures.
3. A JOB WITH NO `.out` FILE READS AS FAILED RATHER THAN MISSING. OPI's
   `terminated_normally()`, `scf_converged()` and
   `geometry_optimization_converged()` all return False -- not raise -- when
   the output file is absent, so a mistyped path looks exactly like a crashed
   job. `load_job` requires the file to exist.
4. AN IMAGINARY FREQUENCY NEEDS A CUTOFF. The C8 reaction complex -- a
   genuine minimum, converged, which ORCA's own thermochemistry treats as one
   ("The first frequency considered to be a vibration is 8") -- carries two
   imaginary modes at -19.76 and -3.55 cm-1. On a floppy 140-atom macrocycle
   that is numerical noise, and a rule reading "more than one imaginary mode
   means this is not a clean saddle point" condemns every real TS this
   project will produce. `stationary_point` applies `IMAGINARY_CUTOFF_CM1`
   and says what it ignored.

Units: energies in Eh from ORCA, kJ/mol everywhere this project quotes them
(`MECHANISM.md` and `COMPUTATIONAL.md` are in kJ/mol throughout; the
superseded skill's sanity range was in kcal/mol, which is an invitation to a
transcription error). Temperatures in K, concentrations in M.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

from opi.output.core import Output

# Hartree -> kJ/mol, from the CODATA values OPI itself is built on.
EH_TO_KJ_PER_MOL = 4.3597447222060e-18 * 6.02214076e23 / 1000.0
R_KJ_PER_MOL_K = 8.314462618 / 1000.0

# Molar volume of an ideal gas, L/mol/K, for the standard-state correction.
_GAS_CONSTANT_L_ATM = 0.082057366080960

# Pure water at 298.15 K and 1 atm: 997.047 g/L / 18.01528 g/mol.
WATER_MOLARITY_M = 55.34

# Below this, an imaginary mode is numerical noise rather than a reaction
# coordinate. Set from the C8 reaction complex, a converged minimum carrying
# imaginary modes at -19.76 and -3.55 cm-1 that ORCA's own thermochemistry
# discards; 50 cm-1 clears both with room to spare and is far below any real
# heavy-atom reaction coordinate (the C8 relay TS's own mode is >1000 cm-1).
# A mode between the cutoff and a real coordinate is a JUDGEMENT CALL, and
# `stationary_point` reports what it ignored so the call can be made.
IMAGINARY_CUTOFF_CM1 = 50.0

# The version every entry in COMPUTATIONAL.md's log is run under.
ORCA_VERSION = "6.1.1"


class OrcaJobError(RuntimeError):
    """A job cannot be read, or was read and did not finish cleanly."""


@dataclass(frozen=True)
class OrcaJob:
    """One finished ORCA job, with its own directory and OPI output attached.

    `output` is the OPI object; everything in this module takes an `OrcaJob`
    rather than a bare path so a job is read and validated exactly once.
    """

    directory: Path
    basename: str
    output: Output

    @property
    def label(self) -> str:
        return f"{self.directory}/{self.basename}"

    @property
    def version(self) -> str | None:
        status = self.output.results_properties.calculation_status
        return getattr(status, "version", None) if status else None

    @property
    def natoms(self) -> int | None:
        info = self.output.results_properties.calculation_info
        return getattr(info, "numofatoms", None) if info else None

    @property
    def charge(self) -> int | None:
        info = self.output.results_properties.calculation_info
        return getattr(info, "charge", None) if info else None

    @property
    def multiplicity(self) -> int | None:
        info = self.output.results_properties.calculation_info
        return getattr(info, "mult", None) if info else None

    @property
    def cpu_seconds(self) -> float | None:
        """Summed module time from ORCA's own timing block.

        This is the number `PLAN_C7_INDUCTION.md` asks to be reported
        alongside `nprocs`; it is CPU time across modules, close to but not
        identical with the wall time the log's last line prints.
        """
        timings = self.output.results_properties.calculation_timings
        return getattr(timings, "sum", None) if timings else None

    @property
    def method(self) -> str | None:
        """The `!` keyword lines of the sibling `.inp`, joined.

        OPI and this project both write one directive per line, so a single
        line is not the whole method.
        """
        inp = self.directory / f"{self.basename}.inp"
        if not inp.exists():
            return None
        keywords = [line.strip().lstrip("!").strip()
                    for line in inp.read_text(errors="replace").splitlines()
                    if line.strip().startswith("!")]
        return " ".join(keywords) if keywords else None

    @property
    def final_energy_eh(self) -> float | None:
        """The job's final energy -- the QM/QM2 TOTAL on a QM/XTB job."""
        return self.output.get_final_energy()


def load_job(directory, basename: str = "job", *, require_converged: bool = True) -> OrcaJob:
    """Read a finished ORCA job.

    `require_converged` raises unless the job terminated normally, its SCF
    converged and -- where the job optimised a geometry -- the optimisation
    converged. Pass False only to inspect a job that is KNOWN to have failed;
    every number that reaches a document comes from a converged job.

    Requires the `.out` to exist: OPI's status checks return False rather
    than raising when it is missing, so a mistyped path would otherwise be
    indistinguishable from a crash.
    """
    directory = Path(directory)
    outfile = directory / f"{basename}.out"
    if not outfile.exists():
        raise OrcaJobError(
            f"no output file at {outfile} -- OPI's status checks return False "
            "for a missing file, so this would otherwise read as a crashed job"
        )

    output = Output(basename=basename, working_dir=directory)
    output.parse()
    job = OrcaJob(directory=directory, basename=basename, output=output)

    if require_converged:
        problems = []
        if not output.terminated_normally():
            problems.append("did not terminate normally")
        if not output.scf_converged():
            problems.append("SCF did not converge")
        if _is_optimisation(job) and not output.geometry_optimization_converged():
            problems.append("the geometry optimisation did not converge")
        if problems:
            raise OrcaJobError(f"{job.label}: " + "; ".join(problems))
    return job


def _is_optimisation(job: OrcaJob) -> bool:
    method = (job.method or "").lower()
    return bool(re.search(r"\bopt\b|\boptts\b|\bcopt\b", method))


# --------------------------------------------------------------------------
# Thermochemistry. Every one of these takes a temperature.
# --------------------------------------------------------------------------

def temperatures(job: OrcaJob) -> list[float]:
    """Every temperature the job's final geometry has thermochemistry for."""
    return [entry.temperature for entry in _thermo_entries(job)
            if entry.temperature is not None]


def _thermo_entries(job: OrcaJob, *, index: int = -1):
    geometries = job.output.results_properties.geometries
    if not geometries:
        raise OrcaJobError(f"{job.label}: no geometries in the property JSON")
    entries = geometries[index].thermochemistry_energies
    if not entries:
        raise OrcaJobError(
            f"{job.label}: no thermochemistry -- was this job run with Freq? "
            "A single point or a plain Opt has no Hessian and no free energy."
        )
    return entries


def thermochemistry(job: OrcaJob, temperature_K: float, *, tolerance: float = 0.01):
    """The thermochemistry block computed at `temperature_K`.

    Raises rather than falling back to another temperature. OPI's own
    accessors take entry [0] and the regex parser it replaces took the last
    one; on a three-temperature job those are different answers and neither
    tells you which it gave.
    """
    entries = _thermo_entries(job)
    for entry in entries:
        if entry.temperature is not None and abs(entry.temperature - temperature_K) <= tolerance:
            return entry
    available = ", ".join(f"{t:g}" for t in temperatures(job)) or "none"
    raise OrcaJobError(
        f"{job.label}: no thermochemistry at {temperature_K:g} K "
        f"(the job computed: {available})"
    )


def free_energy(job: OrcaJob, temperature_K: float) -> float:
    """Gibbs free energy in Eh at `temperature_K`, ZPE and thermal included."""
    return thermochemistry(job, temperature_K).freeenergyg


def enthalpy(job: OrcaJob, temperature_K: float) -> float:
    return thermochemistry(job, temperature_K).enthalpyh


def entropy_term(job: OrcaJob, temperature_K: float) -> float:
    """T*S in Eh -- the term in G = H - T*S, which is what ORCA reports.

    ORCA's `entropys` field is the PRODUCT, not the entropy: on the tracked
    water fixture it equals H - G exactly at all three temperatures. Reading
    it as S would be wrong by a factor of T, which at 298 K is a factor of
    300 and looks nothing like a unit slip.
    """
    return thermochemistry(job, temperature_K).entropys


def entropy(job: OrcaJob, temperature_K: float) -> float:
    """Entropy in Eh/K, i.e. ORCA's reported T*S divided back out by T."""
    return entropy_term(job, temperature_K) / temperature_K


def zero_point_energy(job: OrcaJob, temperature_K: float) -> float:
    return thermochemistry(job, temperature_K).zpe


# --------------------------------------------------------------------------
# Frequencies and what kind of stationary point they say this is.
# --------------------------------------------------------------------------

def frequencies(job: OrcaJob, temperature_K: float | None = None) -> list[float]:
    """Every vibrational frequency in cm-1, including the six near-zero modes.

    The whole list, not a truncated window: the superseded parser read only
    the first 4000 characters after `VIBRATIONAL FREQUENCIES`, which on this
    project's 140-atom systems is about 80 of 420 modes.
    """
    entries = _thermo_entries(job)
    entry = (thermochemistry(job, temperature_K) if temperature_K is not None
             else entries[0])
    if not entry.freq:
        return []
    return [row[0] for row in entry.freq]


def imaginary_modes(job: OrcaJob, *, cutoff_cm1: float = IMAGINARY_CUTOFF_CM1,
                    temperature_K: float | None = None) -> list[float]:
    """Imaginary frequencies past the noise cutoff, most negative first."""
    below = [f for f in frequencies(job, temperature_K) if f < -abs(cutoff_cm1)]
    return sorted(below)


def small_imaginary_modes(job: OrcaJob, *, cutoff_cm1: float = IMAGINARY_CUTOFF_CM1,
                          temperature_K: float | None = None) -> list[float]:
    """Imaginary modes the cutoff discards -- reported, never hidden."""
    return sorted(f for f in frequencies(job, temperature_K)
                  if -abs(cutoff_cm1) <= f < 0.0)


def stationary_point(job: OrcaJob, *, cutoff_cm1: float = IMAGINARY_CUTOFF_CM1,
                     temperature_K: float | None = None) -> str:
    """What the Hessian says this geometry is.

    Returns one of "minimum", "transition state", or
    "saddle point of order N" -- and appends what the cutoff ignored, because
    a mode just under it is a judgement call that has to be visible.
    """
    real = imaginary_modes(job, cutoff_cm1=cutoff_cm1, temperature_K=temperature_K)
    small = small_imaginary_modes(job, cutoff_cm1=cutoff_cm1, temperature_K=temperature_K)
    if len(real) == 0:
        verdict = "minimum"
    elif len(real) == 1:
        verdict = "transition state"
    else:
        verdict = f"saddle point of order {len(real)}"
    if small:
        ignored = ", ".join(f"{f:.2f}" for f in small)
        verdict += f" (ignoring {len(small)} imaginary mode(s) under {cutoff_cm1:g} cm-1: {ignored})"
    return verdict


# --------------------------------------------------------------------------
# Standard states. Functions of temperature, not constants.
# --------------------------------------------------------------------------

def gas_to_molar_correction(temperature_K: float, pressure_atm: float = 1.0) -> float:
    """kJ/mol to move a species from 1 atm to 1 M at `temperature_K`.

    RT ln(V_m / 1 L) with V_m = RT/P. +7.93 kJ/mol at 298.15 K -- the value
    usually quoted, and the one `PLAN_C7_INDUCTION.md` states -- but 7.58 at
    288.15 K and 8.45 at 313.15 K, so over C7's own range a single constant
    is half a kJ/mol wrong at both ends.
    """
    molar_volume_l = _GAS_CONSTANT_L_ATM * temperature_K / pressure_atm
    return R_KJ_PER_MOL_K * temperature_K * math.log(molar_volume_l)


def solvent_standard_state_correction(temperature_K: float,
                                      molarity_M: float = WATER_MOLARITY_M) -> float:
    """kJ/mol to move the SOLVENT from 1 M to its own neat concentration.

    RT ln(molarity), to be SUBTRACTED from a free energy that treats water as
    a 1 M solute. -9.95 kJ/mol at 298.15 K for water at 55.34 M; -9.62 and
    -10.45 at C7's other two temperatures.
    """
    return R_KJ_PER_MOL_K * temperature_K * math.log(molarity_M)


def equilibrium_constant(delta_g_kj_per_mol: float, temperature_K: float) -> float:
    """K = exp(-dG / RT), with dG already on the standard state you want."""
    return math.exp(-delta_g_kj_per_mol / (R_KJ_PER_MOL_K * temperature_K))


def fraction_from_equilibrium_constant(k: float) -> float:
    """K / (1 + K) -- the bound fraction of a two-state equilibrium.

    C7's gate is a hydrate FRACTION, not a rate constant; this is the step
    the superseded skills had no route to at all.
    """
    return k / (1.0 + k)


# --------------------------------------------------------------------------
# The composite energy this project's two-tier method is built on.
# --------------------------------------------------------------------------

def oniom_correction_eh(high_fragment: OrcaJob, low_fragment: OrcaJob) -> float:
    """E(high) - E(low) on the SAME extracted QM-region fragment, in Eh.

    COMPUTATIONAL.md's Conventions: a DLPNO-CCSD(T)/def2-TZVPP single point on
    `job.QMRegion.xyz` minus an r2SCAN-3c single point on the same geometry,
    added to the full embedded energy as ORCA's own subtractive QM1/QM2
    scheme applied post hoc. Raises if the two jobs are not the same
    fragment, which is the mistake this function exists to make impossible.
    """
    if high_fragment.natoms != low_fragment.natoms:
        raise OrcaJobError(
            f"the two fragment jobs are different systems: "
            f"{high_fragment.label} has {high_fragment.natoms} atoms, "
            f"{low_fragment.label} has {low_fragment.natoms}"
        )
    if high_fragment.charge != low_fragment.charge:
        raise OrcaJobError(
            f"the two fragment jobs carry different charges: "
            f"{high_fragment.charge} and {low_fragment.charge}"
        )
    high = high_fragment.final_energy_eh
    low = low_fragment.final_energy_eh
    if high is None or low is None:
        raise OrcaJobError("a fragment job has no final energy")
    return high - low


def composite_free_energy_eh(embedded: OrcaJob, temperature_K: float,
                             high_fragment: OrcaJob, low_fragment: OrcaJob) -> float:
    """G(QM/XTB, at T) + [E(DLPNO) - E(r2SCAN-3c)] on the QM region, in Eh.

    The formula every final number in the C7-C10 register rests on. It lived
    only as prose in `job.inp` comment headers until 2026-09-17, which means
    it was re-derived by hand each time it was used.
    """
    return (free_energy(embedded, temperature_K)
            + oniom_correction_eh(high_fragment, low_fragment))


def to_kj_per_mol(energy_eh: float) -> float:
    return energy_eh * EH_TO_KJ_PER_MOL
