"""
The reduced kinetic model of MECHANISM.md, as code.

MECHANISM.md reduces the 7-step mechanism to three ODEs by three exact
conservation laws plus a C1 pre-equilibrium and a catalyst QSSA. This module is
that reduced system and nothing else: no I/O, no fitting, no data. It exists
separately from `fit_kinetics.py` so the chemistry can be tested against
conservation laws and limiting cases without a dataset in the way.

    v_can = k_can [A]^2 [HOO-]          steps 1-2, Cannizzaro-type
    v3    = k3 [PBA][S]                 step 3, uncatalysed peracid oxidation
    v_seed= (k0 + k5 E0) [H2O2][S]      step 5, plus its uncatalysed analogue
    v6    = k6 E0 [PBA]                 steps 6-7, the catalysed loop

    dA/dt   = -2 v_can + v3 + v_seed + v6
    dPBA/dt =    v_can - v3            - v6
    dS/dt   =    v_can - v3 - v_seed   - v6

    [BA] = S0 + A0 - S - A - PBA        aryl conservation
    [H2O2] ~ const                      (median conversion in this dataset is
                                         0.9%, so this is comfortable)

Why `k0` exists, and why it is not optional
-------------------------------------------
MECHANISM.md writes the enzyme-free limit as "2 ODEs and 2 parameters
(k_can, k3)". Structurally that is right, but the trajectory it produces is
identically zero: at E0 = 0 with A = PBA = 0 both rates vanish, so the system
sits at a fixed point and never starts. Seeding it with a trace of aldehyde does
not rescue it either -- steps 1-2 destroy two aldehydes per peracid and step 3
returns only one, so the catalyst-free loop is a net aldehyde *sink* and the
trace decays. This is the same fact MECHANISM.md states as "step 5 is the only
net source of aldehyde"; the consequence is that enzyme-free controls, which
demonstrably do react, require an E0-independent source.

`k0` is that source: the uncatalysed direct oxidation `S + H2O2 -> A`, the
E0 -> 0 limit of step 5. Writing the seed as `(k0 + k5 E0)` keeps the model
exactly linear in E0, as the reduction promises, with an intercept rather than
through the origin. Both `test_kinetic_model.py` and the docstring above are
kept honest by `test_enzyme_free_without_seed_is_frozen`.

Units: concentrations mM, time s. Hence k_can [mM^-2 s^-1], k3 [mM^-1 s^-1],
k0 [mM^-1 s^-1], k5 [mM^-2 s^-1], k6 [mM^-1 s^-1]; r is dimensionless.
"""
import warnings
from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import solve_ivp

# The fitted parameters, in the order the packing helpers use.
PARAMETER_NAMES = ("k_can", "k3", "k0", "k5", "k6", "r")

# Parameters that are rate constants and so are fitted in log10 space; `r` is a
# ratio of extinction coefficients, is allowed to be exactly zero (that is the
# pure-aldehyde reading the fit is meant to adjudicate), and so is fitted
# linearly.
LOG_PARAMETERS = ("k_can", "k3", "k0", "k5", "k6")

# The state vector integrated by `solve_ivp`.
STATE_NAMES = ("A", "PBA", "S")


@dataclass(frozen=True)
class RateConstants:
    """One set of rate constants, shared across every curve in a fit group."""
    k_can: float = 0.0   # v_can = k_can [A]^2 [HOO-]
    k3: float = 0.0      # v3    = k3 [PBA][S]
    k0: float = 0.0      # uncatalysed seed:  v = k0 [H2O2][S]
    k5: float = 0.0      # catalysed seed:    v = k5 E0 [H2O2][S]
    k6: float = 0.0      # v6    = k6 E0 [PBA]
    r: float = 0.0       # eps(benzoate) / eps(benzaldehyde) at the assay wavelength

    def replace(self, **changes):
        return replace(self, **changes)


@dataclass(frozen=True)
class Conditions:
    """
    One cuvette. `hoo` is [HOO-] in mM and carries all of the pH and
    ionic-strength dependence, so this module never sees a pH -- that
    conversion is `solution_chemistry.hydroperoxide`'s job.
    """
    s0: float            # [S] at t = 0, mM
    h2o2: float          # [H2O2], mM, held constant
    e0: float            # total catalyst, mM (0 for the enzyme-free controls)
    hoo: float           # [HOO-], mM
    a0: float = 0.0      # trace aldehyde present at t = 0, mM


def rates(state, constants, conditions):
    """The four lumped rates, given a state vector. Clamped non-negative."""
    aldehyde, peracid, substrate = (max(x, 0.0) for x in state)
    v_can = constants.k_can * aldehyde * aldehyde * conditions.hoo
    v3 = constants.k3 * peracid * substrate
    v_seed = (constants.k0 + constants.k5 * conditions.e0) * conditions.h2o2 * substrate
    v6 = constants.k6 * conditions.e0 * peracid
    return v_can, v3, v_seed, v6


def rhs(_time, state, constants, conditions):
    """dy/dt for y = (A, PBA, S). Signature matches `solve_ivp`'s."""
    v_can, v3, v_seed, v6 = rates(state, constants, conditions)
    return (
        -2.0 * v_can + v3 + v_seed + v6,   # dA/dt
        v_can - v3 - v6,                   # dPBA/dt
        v_can - v3 - v_seed - v6,          # dS/dt
    )


class _EvaluationLimit(Exception):
    """Raised inside the rhs to abort an integration that will not finish."""


# An optimiser exploring log-space will propose rate constants many orders of
# magnitude from anything physical, and those make the system stiff enough that
# RK45 grinds for minutes on a single curve. Capping rhs evaluations turns that
# from a hang into a failed trial the optimiser simply walks away from. Chosen
# well above what a healthy curve needs (a few thousand) and far below what a
# pathological one would burn.
MAX_RHS_EVALUATIONS = 200_000


def simulate(constants, conditions, times, rtol=1e-8, atol=1e-14,
             max_evaluations=MAX_RHS_EVALUATIONS):
    """
    Integrates one cuvette over `times` (seconds, ascending, starting at 0).

    The reduced system is non-stiff at physical parameter values, so RK45 is the
    documented choice; LSODA is kept as a fallback because a badly-scaled trial
    parameter set during fitting can still stall RK45, and a fit that dies on
    one bad step is worse than one that costs a little more per call. Both are
    capped at `max_evaluations` rhs calls -- see MAX_RHS_EVALUATIONS.

    Returns a dict of arrays keyed by 'A', 'PBA', 'S', 'BA', or None if the
    integration failed or hit the cap -- callers must handle None rather than
    trusting a partial trajectory.
    """
    times = np.asarray(times, dtype=float)
    span = (0.0, float(times[-1])) if times[-1] > 0 else (0.0, 1.0)
    initial = (conditions.a0, 0.0, conditions.s0)

    for method in ("RK45", "LSODA"):
        budget = [max_evaluations]

        def counted(time, state, *args):
            budget[0] -= 1
            if budget[0] < 0:
                raise _EvaluationLimit
            return rhs(time, state, *args)

        try:
            # A trial parameter set the optimiser will reject anyway makes LSODA
            # warn about repeated convergence failures. That branch is expected
            # and already handled -- `success` is False and the caller gets None
            # -- so the warning is noise on stderr, not information.
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning,
                                        module="scipy.integrate.*")
                solution = solve_ivp(
                    counted, span, initial, t_eval=times, args=(constants, conditions),
                    method=method, rtol=rtol, atol=atol,
                )
        except _EvaluationLimit:
            continue
        if solution.success:
            aldehyde, peracid, substrate = solution.y
            aryl_total = conditions.s0 + conditions.a0
            return {
                "A": aldehyde,
                "PBA": peracid,
                "S": substrate,
                "BA": aryl_total - substrate - aldehyde - peracid,
            }
    return None


def observable(constants, conditions, times, **kwargs):
    """
    The baseline-subtracted signal the spectrophotometer sees, in mM of
    aldehyde-equivalent:

        signal(t) = ([A](t) - [A](0)) + r * [BA](t)

    The `- [A](0)` matters. If a trace of aldehyde is present at t = 0 its
    absorbance is inside the baseline that gets subtracted from the data, so the
    model must subtract it too; with a0 = 0 the term is zero and this is
    MECHANISM.md's `signal = [A] + r [BA]` unchanged.

    Multiply by the substrate's extinction coefficient to compare against raw
    absorbance -- which is what `fit_kinetics.py` does, so that residuals carry
    the instrument's own roughly-constant noise rather than a noise rescaled by
    a factor of six between the two substrates.
    """
    trajectory = simulate(constants, conditions, times, **kwargs)
    if trajectory is None:
        return None
    return (trajectory["A"] - conditions.a0) + constants.r * trajectory["BA"]


def aryl_residual(trajectory, conditions):
    """
    Aryl conservation, which the reduction guarantees exactly: S + A + PBA + BA
    must equal S0 + A0 at every time. Returned as the worst absolute deviation,
    for tests and for a cheap integration-quality check during fitting.
    """
    total = (trajectory["S"] + trajectory["A"]
             + trajectory["PBA"] + trajectory["BA"])
    return float(np.max(np.abs(total - (conditions.s0 + conditions.a0))))


# --- packing, for the optimiser ------------------------------------------

def pack(constants, names):
    """Free parameters -> the vector the optimiser works in (log10 for rates)."""
    values = []
    for name in names:
        value = getattr(constants, name)
        values.append(np.log10(value) if name in LOG_PARAMETERS else value)
    return np.asarray(values, dtype=float)


def unpack(vector, names, base):
    """The optimiser's vector -> RateConstants, filling the rest from `base`."""
    changes = {}
    for name, value in zip(names, vector):
        changes[name] = float(10.0 ** value) if name in LOG_PARAMETERS else float(value)
    return base.replace(**changes)
