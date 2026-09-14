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

from curve_metrics import ACCELERATION_SIGMA, acceleration
from fit_dataset import QUANTISATION_SIGMA
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


# --- the extensions, and the old model as their exact limit ----------------

# The saved fit's constants (data/fits/BnOH_25C_Phosphate.json). The arrays
# beside each Conditions below were computed from the UNMODIFIED model BEFORE
# the extensions were added, so this is a check on the code, not on itself.
FITTED_ONE = RateConstants(k_can=10.861875693541291, k3=23.813310077486502,
                           k0=2.1373986699478295e-09, r=1.5232234147468324)
FITTED_TWO = FITTED_ONE.replace(k5=1.4222875249708403e-07,
                                k6=0.004362701086786607)
FROZEN = (
    ("enzyme-free a", Conditions(s0=8.25, h2o2=82.5, e0=0.0, hoo=3e-3),
     FITTED_ONE,
     [0.0, 0.0001818673367089, 0.00036387136924551, 0.00054615199750525,
      0.00072884745091609, 0.00091209349813371, 0.00109602267974743,
      0.00128076357710748, 0.00146644011708472, 0.00165317095237555,
      0.00184106890847486, 0.00203024046776032, 0.00222078534470215,
      0.0024127961329384, 0.00260635802049177, 0.00280154857549017,
      0.00299843761341052, 0.00319708713358467, 0.00339755132106601,
      0.00359987661222762, 0.00380410179761462, 0.00401025822457742,
      0.00421837001144157, 0.00442845430253408, 0.00464052159351484]),
    ("enzyme-free b", Conditions(s0=2.0, h2o2=82.5, e0=0.0, hoo=1e-3),
     FITTED_ONE,
     [0.0, 4.4083820794724076e-05, 8.816942832570035e-05,
      1.3225958171442153e-04, 1.7635703742943521e-04, 2.2046454775964764e-04,
      2.6458485998246669e-04, 3.0872071427663124e-04, 3.5287484224733821e-04,
      3.9704996643368927e-04, 4.4124879874789106e-04, 4.8547403941198947e-04,
      5.2972837542936656e-04, 5.7401447922958155e-04, 6.1833500802714809e-04,
      6.6269260177537215e-04, 7.0708988230839805e-04, 7.5152945175593888e-04,
      7.9601389364677684e-04, 8.4054576855048632e-04, 8.8512761460198279e-04,
      9.2976194619489119e-04, 9.7445125180033274e-04, 1.0191979952287225e-03,
      1.0640046121924637e-03]),
    ("catalysed a", Conditions(s0=8.25, h2o2=82.5, e0=0.24, hoo=3e-3),
     FITTED_TWO,
     [0.0, 0.0030921277197304, 0.00622266100488822, 0.00942732648520284,
      0.01273496326012171, 0.016165743306201, 0.01973078424120396,
      0.02343294151205132, 0.02726836568858948, 0.03122837605214235,
      0.03530128768476624, 0.03947396531628972, 0.04373300805612914,
      0.04806556525963593, 0.05245983031080333, 0.05690529141483466,
      0.06139280643235355, 0.06591456682745087, 0.07046399489555706,
      0.07503561084418177, 0.0796248868472038, 0.0842281080057138,
      0.08884224292097022, 0.09346482765313703, 0.09809386622156055]),
    ("catalysed b", Conditions(s0=4.0, h2o2=35.0, e0=0.05, hoo=1e-3),
     FITTED_TWO,
     [0.0, 0.00016185755565004, 0.00032374574760195, 0.00048570170885661,
      0.00064776243924092, 0.00080996473403998, 0.00097234512680013,
      0.00113493982755728, 0.00129778465668464, 0.00146091498571072,
      0.0016243656840975, 0.00178817105507933, 0.00195236477863052,
      0.00211697986537249, 0.00228204860864318, 0.00244760253370953,
      0.00261367235182363, 0.00278028791688684, 0.00294747818863066,
      0.00311527119186763, 0.00328369398950112, 0.00345277264980683,
      0.00362253221805589, 0.00379299669318003, 0.00396418900762511]),
)


def test_defaults_are_the_old_model():
    print("\nthe extensions switch off exactly at their defaults")
    times = np.linspace(0.0, 3000.0, 25)
    for name, conditions, constants, expected in FROZEN:
        signal = observable(constants, conditions, times)
        worst = float(np.max(np.abs(signal - np.asarray(expected))))
        check(f"{name}: extended model = old model to 1e-10",
              worst < 1e-10, f"worst deviation {worst:.2e}")


def test_conservation_survives_the_extension():
    print("\nconservation with every extension on")
    constants = RateConstants(k_can=6.0, k3=1e-2, k0=1e-9, k5=1.0, k6=1.0,
                              r=0.3, k_sink=1e-3, K4=0.05, km_s=3.0,
                              k_act_r=5e-3, K_act=0.05, km_s_background=3.0)
    conditions = Conditions(s0=8.25, h2o2=82.5, e0=0.1, hoo=3e-3,
                            buf=50.0, species=3e-3)
    trajectory = simulate(constants, conditions, TIMES)
    residual = aryl_residual(trajectory, conditions) if trajectory else np.nan
    check("aryl conservation still holds with the sink and all factors on",
          trajectory is not None and residual < 1e-8, f"residual {residual:.2e}")


def test_activation_makes_a_lag_with_r_below_one():
    """
    The old model can only lag if r > 1 (test_acceleration_requires_r_above_one).
    A finite `k_act_r` gives the catalyst its own clock and a lag at any r.
    """
    print("\nactivation supplies a lag the r <= 1 reading could not")
    # k5 at the real fitted scale, so substrate is not burnt in the first
    # seconds -- a fast seed makes the curve saturate, which is concave and
    # would hide the very lag this test is about.
    constants = RateConstants(k5=1e-6, k_act_r=2e-3, r=0.2)
    conditions = Conditions(s0=8.25, h2o2=82.5, e0=0.1, hoo=3e-3)
    times = np.linspace(0.0, 3000.0, 200)
    signal = observable(constants, conditions, times)
    z, where = acceleration(times, signal, floor=QUANTISATION_SIGMA)
    check("r = 0.2 with a finite k_act_r accelerates",
          signal is not None and z > ACCELERATION_SIGMA,
          f"z = {z:.1f} at {where:.2f} of the run")
    flat = observable(constants.replace(k_act_r=float("inf")), conditions, times)
    z_flat, _ = acceleration(times, flat, floor=QUANTISATION_SIGMA)
    check("...and switching the activation off returns the concave reading",
          z_flat <= 1.0, f"z = {z_flat:.2f}")


def test_substrate_saturation_order():
    print("\nsubstrate saturation")
    km_s = 3.0
    constants = RateConstants(k5=1.0, km_s=km_s)
    conditions = Conditions(s0=5.0, h2o2=50.0, e0=0.05, hoo=1e-3)
    s0 = 5.0
    step = 1e-4
    low = rates((0.0, 0.0, s0), constants, conditions, time=0.0)[2]
    high = rates((0.0, 0.0, s0 * (1.0 + step)), constants, conditions, time=0.0)[2]
    order = float(np.log(high / low) / np.log(1.0 + step))
    expected = km_s / (km_s + s0)
    check("the initial rate's substrate order is km_s/(km_s + S)",
          abs(order - expected) < 0.01 * expected,
          f"{order:.4f} against {expected:.4f}")


def test_peroxide_binding_saturates():
    print("\nperoxide binding")
    K4 = 0.05
    constants = RateConstants(k5=1.0, K4=K4)
    substrate = 5.0
    for x in (1.0, 5.0, 20.0, 100.0):
        conditions = Conditions(s0=substrate, h2o2=x, e0=0.05, hoo=1e-3,
                                species=x)
        seed = rates((0.0, 0.0, substrate), constants, conditions, time=0.0)[2]
        prefactor = constants.k5 * conditions.e0 * substrate
        shape = seed / prefactor
        check(f"the catalysed seed at [H2O2] = {x:g} follows x/(1 + K4 x)",
              abs(shape - x / (1.0 + K4 * x)) < 0.01 * (x / (1.0 + K4 * x)),
              f"{shape:.5f} against {x / (1.0 + K4 * x):.5f}")


def test_sink_is_first_order():
    print("\nthe product sink")
    k_sink = 1e-2
    constants = RateConstants(k_sink=k_sink)
    conditions = Conditions(s0=8.0, h2o2=82.5, e0=0.0, hoo=3e-3, a0=1e-3)
    trajectory = simulate(constants, conditions, TIMES)
    expected = conditions.a0 * np.exp(-k_sink * TIMES)
    check("a pulse of aldehyde decays as exp(-k_sink t)",
          trajectory is not None
          and np.allclose(trajectory["A"], expected, rtol=1e-6, atol=1e-11),
          f"worst {np.max(np.abs(trajectory['A'] - expected)):.2e}"
          if trajectory else "no solution")
    check("...and benzoate absorbs exactly what the aldehyde loses",
          trajectory is not None
          and abs((conditions.a0 - trajectory["A"][-1]) - trajectory["BA"][-1]) < 1e-12)


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
    test_defaults_are_the_old_model()
    test_conservation_survives_the_extension()
    test_activation_makes_a_lag_with_r_below_one()
    test_substrate_saturation_order()
    test_peroxide_binding_saturates()
    test_sink_is_first_order()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
