"""
Tests for kinetic_model.py.

Three kinds of test. The first pins the reduction's exact guarantees -- aryl
conservation, linearity in E0 -- so the code is checked against the algebra in
MECHANISM.md rather than against itself. The second pins the structural results
this module was written to establish: that the enzyme-free limit is frozen
without a seed, and that the observable cannot accelerate unless r > 1. The
third is ordinary robustness: the solver must fail loudly rather than hang.

    python data/test_kinetic_model.py
"""
import sys

import numpy as np

from kinetic_model import (
    LOG_PARAMETERS, PARAMETER_NAMES, Conditions, RateConstants,
    aryl_residual, observable, pack, rates, rhs, simulate, unpack,
)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


# A parameter set in the regime the data occupies: ~1% conversion over a
# 9000 s run, aldehyde staying small while benzoate accumulates.
NOMINAL = RateConstants(k_can=6.0, k3=1e-2, k0=1e-9, k5=1.0, k6=1.0, r=0.3)
CUVETTE = Conditions(s0=8.25, h2o2=82.5, e0=0.0, hoo=3e-3)
TIMES = np.linspace(0.0, 9000.0, 200)


def _slope_ratio(signal, times):
    """max slope / initial slope, and where the maximum sits in the run."""
    slope = np.gradient(signal, times)
    if slope[0] <= 0:
        return 0.0, 0.0
    return slope.max() / slope[0], times[np.argmax(slope)] / times[-1]


# --- the reduction's exact guarantees -------------------------------------

def test_conservation():
    print("\naryl conservation (exact in the reduction)")
    for e0 in (0.0, 0.1):
        conditions = Conditions(s0=8.25, h2o2=82.5, e0=e0, hoo=3e-3)
        trajectory = simulate(NOMINAL, conditions, TIMES)
        check(f"S + A + PBA + BA = S0 at E0 = {e0}",
              trajectory is not None and aryl_residual(trajectory, conditions) < 1e-9,
              f"residual {aryl_residual(trajectory, conditions):.2e}" if trajectory else "no solution")

    conditions = Conditions(s0=8.25, h2o2=82.5, e0=0.0, hoo=3e-3, a0=1e-3)
    trajectory = simulate(NOMINAL, conditions, TIMES)
    check("the law absorbs a trace aldehyde into the total (S0 + A0)",
          trajectory is not None and aryl_residual(trajectory, conditions) < 1e-9)


def test_linear_in_e0():
    print("\nlinearity in E0 (the reduction's identifiability payoff)")
    state = (1e-3, 1e-4, 8.0)
    base = Conditions(s0=8.0, h2o2=82.5, e0=0.0, hoo=3e-3)
    one = Conditions(s0=8.0, h2o2=82.5, e0=0.1, hoo=3e-3)
    two = Conditions(s0=8.0, h2o2=82.5, e0=0.2, hoo=3e-3)
    at_zero = np.array(rhs(0.0, state, NOMINAL, base))
    at_one = np.array(rhs(0.0, state, NOMINAL, one))
    at_two = np.array(rhs(0.0, state, NOMINAL, two))
    check("the rhs is exactly affine in E0",
          bool(np.allclose(at_two - at_one, at_one - at_zero, rtol=1e-12)),
          f"{at_two - at_one} vs {at_one - at_zero}")

    catalysed_only = NOMINAL.replace(k0=0.0)
    check("k5 and k6 do nothing at E0 = 0",
          bool(np.allclose(rhs(0.0, state, catalysed_only, base),
                           rhs(0.0, state, catalysed_only.replace(k5=0.0, k6=0.0), base))))


def test_states_stay_physical():
    print("\nphysicality")
    trajectory = simulate(NOMINAL, CUVETTE, TIMES)
    check("no state goes negative",
          trajectory is not None
          and all(trajectory[name].min() > -1e-12 for name in ("A", "PBA", "S", "BA")))
    check("substrate only depletes",
          trajectory is not None and np.all(np.diff(trajectory["S"]) <= 1e-15))
    check("benzoate only accumulates",
          trajectory is not None and np.all(np.diff(trajectory["BA"]) >= -1e-15))

    negative = rates((-1.0, -1.0, -1.0), NOMINAL, CUVETTE)
    check("rates clamp a negative state to zero rather than propagating it",
          all(v == 0.0 for v in negative), str(negative))


# --- the structural results ------------------------------------------------

def test_enzyme_free_without_seed_is_frozen():
    """
    MECHANISM.md writes the E0 = 0 limit as "2 ODEs and 2 parameters
    (k_can, k3)". This is why that limit needs a third: with no seed the
    trajectory is identically zero.
    """
    print("\nthe enzyme-free limit needs a seed")
    unseeded = RateConstants(k_can=6.0, k3=1e-2, k0=0.0, r=0.3)
    trajectory = simulate(unseeded, CUVETTE, TIMES)
    check("with k0 = 0 and A(0) = 0 the system never leaves its fixed point",
          trajectory is not None
          and trajectory["A"].max() < 1e-15 and trajectory["BA"].max() < 1e-15,
          f"A reached {trajectory['A'].max():.2e}" if trajectory else "no solution")

    seeded_by_trace = Conditions(s0=8.25, h2o2=82.5, e0=0.0, hoo=3e-3, a0=1e-3)
    trajectory = simulate(unseeded, seeded_by_trace, TIMES)
    check("a trace of aldehyde decays rather than growing (the loop is a net sink)",
          trajectory is not None and trajectory["A"][-1] <= trajectory["A"][0] + 1e-15,
          f"A: {trajectory['A'][0]:.3e} -> {trajectory['A'][-1]:.3e}" if trajectory else "")

    seeded = simulate(unseeded.replace(k0=1e-9), CUVETTE, TIMES)
    check("with k0 > 0 the reaction runs",
          seeded is not None and seeded["A"][-1] > 1e-6)


def test_acceleration_requires_r_above_one():
    """
    MECHANISM.md derives dA/dt <= v5(0) -- the aldehyde reading cannot show a
    lag -- and proposes signal = A + r*BA as the fix. The same bound survives
    the fix for every r <= 1:

        d(signal)/dt = v_seed + (1 + r)(v3 + v6) - 2 v_can
        dPBA/dt      = v_can - (v3 + v6) >= 0  while peracid accumulates
        =>  d(signal)/dt <= v_seed + (r - 1) v_can <= v_seed(0)

    So only r > 1 -- benzoate absorbing MORE strongly than benzaldehyde --
    lets the model produce the lag the data shows.
    """
    print("\nthe observation equation only rescues the lag if r > 1")
    for r in (0.0, 0.3, 0.9, 1.0):
        signal = observable(NOMINAL.replace(r=r), CUVETTE, TIMES)
        ratio, _ = _slope_ratio(signal, TIMES)
        check(f"r = {r:.1f}: signal is concave throughout (steepest at t = 0)",
              signal is not None and ratio <= 1.0 + 1e-9, f"slope ratio {ratio:.4f}")

    for r in (1.5, 3.0):
        signal = observable(NOMINAL.replace(r=r), CUVETTE, TIMES)
        ratio, position = _slope_ratio(signal, TIMES)
        check(f"r = {r:.1f}: signal accelerates",
              signal is not None and ratio > 1.05,
              f"slope ratio {ratio:.4f}")

    print("        for reference, the data: 52% of curves reach peak slope more")
    print("        than 15% into the run, 129 of them at 1.5-5x the initial slope.")


def test_seed_alone_is_linear():
    print("\nthe seed alone")
    seed_only = RateConstants(k_can=0.0, k3=0.0, k0=1e-9, r=0.0)
    trajectory = simulate(seed_only, CUVETTE, TIMES)
    expected = seed_only.k0 * CUVETTE.h2o2 * CUVETTE.s0 * TIMES
    check("with k_can = k3 = 0 the aldehyde grows at k0*[H2O2]*[S]",
          trajectory is not None
          and np.allclose(trajectory["A"], expected, rtol=2e-3),
          f"max rel. error {np.max(np.abs(trajectory['A'] - expected) / expected[-1]):.2e}"
          if trajectory else "")
    check("no peracid or benzoate is made without steps 1-2",
          trajectory is not None
          and trajectory["PBA"].max() < 1e-15 and trajectory["BA"].max() < 1e-15)


def test_observable():
    print("\nthe observation equation")
    signal = observable(NOMINAL, CUVETTE, TIMES)
    check("the baseline-subtracted signal starts at zero",
          signal is not None and abs(signal[0]) < 1e-15)

    with_trace = Conditions(s0=8.25, h2o2=82.5, e0=0.0, hoo=3e-3, a0=5e-4)
    signal = observable(NOMINAL, with_trace, TIMES)
    check("it also starts at zero when a trace aldehyde is present at t = 0",
          signal is not None and abs(signal[0]) < 1e-15,
          f"signal(0) = {signal[0]:.2e}" if signal is not None else "")

    pure = observable(NOMINAL.replace(r=0.0), CUVETTE, TIMES)
    trajectory = simulate(NOMINAL, CUVETTE, TIMES)
    check("r = 0 recovers the pure-aldehyde reading",
          pure is not None and np.allclose(pure, trajectory["A"], rtol=1e-8))


# --- robustness ------------------------------------------------------------

def test_solver_gives_up_cleanly():
    print("\nsolver robustness")
    absurd = RateConstants(k_can=1e12, k3=1e12, k0=1e-3, k5=1e12, k6=1e12, r=0.3)
    result = simulate(absurd, CUVETTE, TIMES, max_evaluations=200)
    check("an unintegrable parameter set returns None, not a partial trajectory",
          result is None, f"got {type(result).__name__}")
    check("observable() propagates that None",
          observable(absurd, CUVETTE, TIMES, max_evaluations=200) is None)


def test_packing():
    print("\npacking for the optimiser")
    names = ("k_can", "k3", "r")
    vector = pack(NOMINAL, names)
    check("rate constants are packed in log10, r linearly",
          np.allclose(vector, [np.log10(NOMINAL.k_can), np.log10(NOMINAL.k3), NOMINAL.r]))
    restored = unpack(vector, names, RateConstants())
    check("unpack inverts pack",
          all(np.isclose(getattr(restored, n), getattr(NOMINAL, n)) for n in names))
    check("unpack leaves untouched parameters at their base values",
          restored.k5 == 0.0 and restored.k6 == 0.0)
    check("every fitted parameter is either logged or explicitly linear",
          set(LOG_PARAMETERS) <= set(PARAMETER_NAMES) and "r" not in LOG_PARAMETERS)


if __name__ == "__main__":
    test_conservation()
    test_linear_in_e0()
    test_states_stay_physical()
    test_enzyme_free_without_seed_is_frozen()
    test_acceleration_requires_r_above_one()
    test_seed_alone_is_linear()
    test_observable()
    test_solver_gives_up_cleanly()
    test_packing()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
