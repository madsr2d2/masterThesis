"""
What the catalysed 4OMe-BnOH progress curves are doing before they start.

Every catalysed 4OMe run from 15 to 30 C begins slowly and speeds up over
thousands of seconds before it reaches the rate `temperature_series/ANALYSIS.md`
puts on an Arrhenius plot. That document names the induction, times it
(tau falls 6489 -> 3190 -> 945 -> 916 s from 15 to 30 C) and stops there.
`product_fate/` did the same job for the FALL at the other end of the curve.
This module is the attempt to name the RISE.

Four candidates, and each makes a different prediction about what the induction
waits for:

    seeding      the catalyst-free Cannizzaro loop needs product before it can
                 run (MECHANISM.md steps 1-3), so the induction ends when
                 enough A has been made.        -> needs no catalyst;
                                                   ends at a fixed PRODUCT
    scavenger    something in the reagents consumes the oxidant until it is
                 burned off.                    -> needs no catalyst;
                                                   ends at a fixed TURNOVER
    activation   the catalyst is not in its active form when the run starts
                 and converts into it.          -> needs the catalyst;
                                                   ends at a fixed TIME
    schedule     the induction is an artefact of when the operator started and
                 stopped recording.             -> tracks the run length and
                                                   nothing else

The first two end on the PRODUCT, the third on a CLOCK. That is the same
discrimination `slowdown.deceleration_drivers` makes for the fall, and it is
made the same way: within one curve the product only grows with time, so the
two can only be separated ACROSS curves whose rates differ. The substrate
ladder inside every 4OMe run moves the rate by half an order while holding the
schedule, the peroxide, the catalyst, the buffer and the temperature fixed,
which is exactly the lever this needs.

THE STATISTIC, AND WHY IT IS A LANDMARK HERE WHEN IT WAS NOT THERE.
`slowdown` withdrew a landmark -- when does the rate fall to three quarters of
its peak -- because a curve whose rate never falls that far has no landmark,
and dropping those curves biases the answer towards the clock. The rise has no
such problem, and for a reason that is structural rather than lucky: the
landmark here is a fraction of the curve's OWN maximum, and a curve always
reaches its own maximum. `induction_point` is therefore defined for every
curve in the archive, and a curve with no induction returns zero rather than
nothing. Nothing here is censored.

What it costs instead is a WINDOW. The rolling slope is read through a window
of a tenth of the run (`curve_metrics.LAG_WINDOW`), because a slow curve needs
a wide window to have a slope at all -- at 15 C the whole run is 0.04 AU -- and
a fixed window in seconds turns the 15 C curves into noise (a 300 s window puts
their induction at 529 s against 4289 s, and drags the block's activation
energy from 86 to 24 kJ/mol). A window that is a fraction of the run is safe
for comparisons WITHIN a run, where every cuvette shares the schedule, and
unsafe between runs, where it does not: across the 4OMe archive at 25 C the
induction time regresses on run length with an exponent of +0.75 +- 0.14 and on
pH, once run length is in the model, with -0.05 +- 0.09.

    SO: every concentration order here is measured within experiments.
    Between-run comparisons of `t_ind` are not evidence, and the one
    between-run quantity that is -- the temperature dependence -- is taken
    from `arrhenius`'s fitted `inverse_tau`, which is not windowed.

    python data/induction.py
"""
import sys
import os
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from curve_metrics import LAG_WINDOW, rolling_slope
from fit_dataset import source_floor
import arrhenius
import scope


# The fraction of the rise that marks "the reaction has started". Half, the
# same level curve_metrics.lag_time uses, so the two are the same landmark and
# differ only in the gate described below.
INDUCTION_LEVEL = 0.5

# The log floor for a curve with no induction at all. `t_ind` is a time and a
# curve that is fastest in its first window returns 0.0 honestly; the
# regressions take logs, so they need a floor rather than a dropped row --
# dropping is the censoring this statistic was built to avoid.
#
# 60 s is the coarsest sampling interval in the archive (they run 28-60 s), and
# the rolling slope's window centres are one interval apart, so 60 s is the
# shortest induction the readings can resolve: the smallest non-zero `t_ind`
# in the 4OMe archive is 49 s. Putting the floor there is what makes "no
# induction" and "an induction shorter than one reading" the same number,
# which is what they are. `report` prints the answer at floors from 1 to 300 s;
# the coefficient that carries this module's conclusion moves from +0.000 to
# -0.021 across that range, so nothing here rests on the choice.
INDUCTION_FLOOR = 60.0

# Below this the "induction" is a fraction of a percent of the peak rate and
# is not distinguishable from where the first window happened to land. It is a
# resolution threshold and not a gate: a curve below it gets `t_ind = 0`, which
# is a measurement, and keeps its row. An exact straight line needs it -- its
# rolling slopes are equal to within floating-point dust, the largest of them
# lands wherever the dust does, and without this a noiseless line reports an
# induction of 3890 s. `curve_metrics.lag_time` handles the same pathology with
# an acceleration gate, which is the thing this statistic may not have.
DEPTH_FLOOR = 0.01

# Where the two hypotheses put the slope of log(induction time) on log(rate).
# A clock fixes the TIME and lets the product land where it may; product
# control fixes the PRODUCT at which the induction ends, so the time to reach
# it is inversely proportional to the rate.
INDUCTION_CLOCK_SLOPE = 0.0
INDUCTION_PRODUCT_SLOPE = -1.0

# The floors `report` sweeps, to show what the answer owes to INDUCTION_FLOOR.
FLOOR_SWEEP = (1.0, 30.0, 60.0, 120.0, 300.0)


@dataclass(frozen=True)
class InductionPoint:
    """
    Where one curve's rate first reaches half of its own maximum.

    `t_ind` is that time in seconds from the first reading, `made` the
    absorbance built up by then, `depth` the fraction of the peak rate that is
    missing at the start, and `t_peak` where the maximum itself is. `depth` is
    a FRACTION on purpose: a catalyst that starts partly inactive gives the
    same fraction whatever the run's rate, and a build-up whose level depends
    on the substrate does not.

    `depth` IS A LOWER BOUND, not the amplitude. The first rolling slope is an
    average over a tenth of the run, so a curve that truly starts at zero rate
    reads back at `(tau/w)(1 - e^(-w/tau))` for a window `w`: a relaxation as
    long as the window reads 0.63 rather than 1.00, and only tau >> w recovers
    the whole of it. Every depth quoted anywhere is therefore an underestimate.
    The bias is towards zero for every curve and it is deepest where the
    induction is fastest, so the depths in this module understate how different
    the cold runs are from the warm ones, not the reverse.
    """
    t_ind: float
    made: float
    peak_rate: float
    start_rate: float
    depth: float
    t_peak: float
    points: int


def induction_point(curve, level=INDUCTION_LEVEL, window=LAG_WINDOW):
    """
    The landmark, read off the READINGS through `curve_metrics.rolling_slope`.

    No model and no extrapolation, the same discipline `slowdown.sink_fit`
    uses at the other end of the curve. The floor passed to the rolling slope
    is the curve's own -- `.rre` readings and `.txt` exports differ by a factor
    of 1096 and a floor left at the export's suppresses the very structure this
    is looking for (CLAUDE.md).

    This is `curve_metrics.lag_time` WITHOUT its acceleration gate, and the
    difference is the point: `lag_time` returns nan for a curve whose rise does
    not clear 3 sigma, which is the right answer when the question is "how many
    curves have a lag" and the wrong one when the question is "does the length
    of the lag track the rate", because the gate is passed more often by fast
    curves. `test_induction` checks that the two agree wherever both are
    defined.
    """
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    centres, slopes = rolling_slope(times, values, window,
                                    source_floor(curve.source))
    blank = InductionPoint(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0)
    if len(slopes) < 3:
        return blank
    top = int(np.argmax(slopes))
    peak, start = float(slopes[top]), float(slopes[0])
    if not peak > 0:
        return blank
    if peak <= start:
        # Fastest in its first window: no induction, and that is a measurement
        # rather than a failure. Every enzyme-free 4OMe curve lands here.
        return InductionPoint(0.0, 0.0, peak, start, 0.0, 0.0, len(slopes))
    if 1.0 - start / peak <= DEPTH_FLOOR:
        # A rise smaller than the resolution threshold. Same answer as no rise
        # at all, and for the same reason: nothing here is a measurement.
        return InductionPoint(0.0, 0.0, peak, start, 0.0, 0.0, len(slopes))
    threshold = start + level * (peak - start)
    index = int(np.flatnonzero(slopes >= threshold)[0])
    return InductionPoint(
        t_ind=float(centres[index] - centres[0]),
        made=float(np.interp(centres[index], times, values) - values[0]),
        peak_rate=peak, start_rate=start,
        depth=float(1.0 - start / peak),
        t_peak=float(centres[top] - centres[0]),
        points=int(len(slopes)))


def induction_table(experiments, frame=None):
    """
    `scope.frame`'s columns for these experiments, plus the five landmark ones.

    Built by joining onto the frame rather than by recomputing any of it, so
    every condition, rate and flag here is the same number the rest of the
    package uses.
    """
    import pandas as pd
    experiments = tuple(experiments)
    if frame is None:
        frame = scope.frame(experiments)
    frame = sign_table(frame)
    rows = []
    for curve in scope.curves(experiments):
        found = induction_point(curve)
        rows.append({"experiment": curve.experiment, "sample": curve.sample,
                     "t_ind": found.t_ind, "made": found.made,
                     "peak_rate": found.peak_rate,
                     "start_rate": found.start_rate, "depth": found.depth,
                     "t_peak": found.t_peak})
    landmarks = pd.DataFrame(rows)
    return frame.merge(landmarks, on=["experiment", "sample"], how="inner")


def substrate_lever(table):
    """
    How far the substrate ladder moves [S] INSIDE a run, and the rate with it.

    The whole of section 3 rests on this lever: a design that changes the rate
    while the schedule, the peroxide, the catalyst, the buffer, the pH and the
    temperature stay where they are. Quoted as the median over experiments,
    because two runs of the archive carry a single cuvette each.
    """
    live = table[table.live & (table.s0 > 0)]
    per_run = live.groupby("experiment").s0.agg(
        lambda column: float(column.max() / column.min()))
    laddered = per_run[per_run > 1]
    return {"experiments": int(len(per_run)),
            "laddered": int(len(laddered)),
            "median_lever": float(laddered.median()) if len(laddered) else np.nan,
            "largest_lever": float(per_run.max()),
            "curves": int(len(live))}


def induction_drivers(table, response="t_ind", rate="peak_rate",
                      floor=INDUCTION_FLOOR, fixed=True):
    """
    Regress how long the induction lasted on how FAST the curve went.

        log(induction time) = b . log(rate) + one offset per experiment

        b = 0    a clock: the induction takes the time it takes
        b = -1   product control: it ends at a fixed amount of product, so a
                 curve twice as fast gets there in half the time

    `fixed=True` puts one dummy per experiment, which absorbs temperature, pH,
    buffer, enzyme, cell, day AND the run length -- and the run length is the
    one that matters, because the window this landmark is read through is a
    fraction of it. Inside an experiment every cuvette shares a schedule, so
    the coefficient is carried only by the substrate ladder, which is what the
    ladder is for. Set `fixed=False` only to see how much worse the pooled
    answer is.

    The rate regressor is the landmark's own peak, not `vmax`, so the two sides
    of the regression come from the same windows; `report` prints the same fit
    against `vmax` and `v_peak` as a check that the answer is not an artefact
    of that choice.

    BEWARE THE DIRECTION OF THE BIAS THIS CANNOT REMOVE. Response and
    regressor are both measured off the same curve, so noise in the rate is
    errors-in-variables and attenuates `b` TOWARDS ZERO -- towards the clock,
    which is the answer this returns. A slope of zero here is therefore not by
    itself evidence. `order_ratio` is the route that does not have the problem:
    it replaces the measured rate with the composition the operator set, and
    it has to agree.
    """
    live = table[table.live & np.isfinite(table[response])
                 & (table[rate] > 0)].copy()
    if len(live) < 8:
        return {"points": int(len(live))}
    response_values = np.log(np.maximum(live[response].to_numpy(dtype=float),
                                        floor))
    columns = [np.log(live[rate].to_numpy(dtype=float))]
    if fixed:
        for experiment in sorted(live.experiment.unique()):
            columns.append((live.experiment == experiment).to_numpy(float))
    else:
        columns.append(np.ones(len(live)))
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, response_values, rcond=None)
    resid = response_values - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(live) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    spread = float(((response_values - response_values.mean()) ** 2).sum())
    return {"points": int(len(live)),
            "slope": float(beta[0]),
            "stderr": float(np.sqrt(max(covariance[0, 0], 0.0))),
            "r2": float(1 - float(resid @ resid) / spread) if spread else np.nan,
            "floored": int((live[response].to_numpy(dtype=float) < floor).sum()),
            "experiments": int(live.experiment.nunique())}


def order_ratio(table, response="t_ind", axis="s0", floor=INDUCTION_FLOOR,
                rate="v_peak"):
    """
    The same discrimination, done through the two substrate orders separately.

    `induction_drivers` regresses the induction time on the rate directly,
    which puts a quantity measured off the same windows on both sides.
    This route instead measures the order of the induction time in [S] and the
    order of the RATE in [S] -- two independent regressions against a
    composition the operator set -- and divides:

        d log(t_ind) / d log(rate) = order(t_ind) / order(rate)

    Product control predicts -1 and a clock predicts 0, exactly as above, and
    the answer must not depend on which route is taken. The error is the delta
    method on the ratio, with the two orders treated as independent: they are
    fitted to different responses on the same design, so their errors are
    correlated, and ignoring that makes this interval slightly too wide.

    Both orders come from `scope.orders`, which is the package's own
    within-experiment log-log fit; nothing about the design is re-derived here.
    """
    numerator = scope.orders(response, frame=table, floor=floor)
    denominator = scope.orders(rate, frame=table)
    top = numerator[f"order_{axis}"]
    bottom = denominator[f"order_{axis}"]
    ratio = top / bottom if bottom else np.nan
    error = (abs(ratio) * np.sqrt((numerator[f"stderr_{axis}"] / top) ** 2
                                  + (denominator[f"stderr_{axis}"] / bottom) ** 2)
             if top and bottom else np.nan)
    return {"axis": axis,
            "induction_order": top,
            "induction_stderr": numerator[f"stderr_{axis}"],
            "rate_order": bottom, "rate_stderr": denominator[f"stderr_{axis}"],
            "ratio": float(ratio), "ratio_stderr": float(error),
            "points": numerator["n"], "rate_points": denominator["n"]}


# How near a matched cell has to be in pH before two blocks may be compared.
# The 4OMe archive's nominal pH values are quoted to two decimals; 0.3 spans
# the 6.97/7.00 pair and the 7.50/7.53 pair and admits nothing else. Buffer
# identity and [H2O2] have to match exactly, not nearly: the buffer is a
# catalyst in this chemistry (temperature_series/ANALYSIS.md 3) and [H2O2] is
# the largest single lever on the rate there is.
MATCH_PH = 0.3
MATCH_TEMPERATURES = (25.0, 40.0)


def channel_contrast(frame, temperatures=MATCH_TEMPERATURES, tolerance=MATCH_PH):
    """
    The induction with the catalyst and without it, at matched composition.

    `differential` is the structural classification and not the filename: a
    run whose reference channel omits the ENZYME is a catalytic increment, one
    whose reference omits the H2O2 is the raw background (see `scope.frame`).
    Two of the three candidates -- a product-seeded loop and a scavenger in the
    reagents -- run in the cuvette with no catalyst in it, so a contrast at
    matched substrate, peroxide, buffer, pH and temperature is the whole test.

    The cell is centred on the ENZYME-FREE block, which is the scarce one: it
    exists at two temperatures and one buffer, and centring on the catalysed
    median instead lands at pH 8.5 where there is nothing to compare with.

    Returns one row per (temperature, channel) cell with the induction depth,
    the acceleration count and the longest run in it. The last column is there
    because the obvious objection to a null in the enzyme-free channel is that
    those runs were too short to show anything.
    """
    import pandas as pd
    four = frame[frame.substrate == "4OMe-BnOH"]
    rows = []
    for temperature in temperatures:
        cell = four[four.temperature == temperature]
        free = cell[~cell.differential]
        if not len(free):
            continue
        centre = float(free.pH.median())
        cell = cell[((cell.pH - centre).abs() <= tolerance)
                    & cell.buffer.isin(free.buffer.unique())
                    & cell.h2o2.isin(free.h2o2.unique())]
        for differential, block in cell.groupby("differential"):
            rows.append({
                "temperature": temperature,
                "channel": "catalysed" if differential else "enzyme-free",
                "curves": int(len(block)),
                "experiments": tuple(sorted(int(e) for e
                                            in block.experiment.unique())),
                "depth": float(block.depth.median()),
                "deep": int((block.depth > 0.25).sum()),
                "accelerates": int(block.accelerates.sum()),
                "longest_s": float(block.duration_s.max()),
                "pH": f"{block.pH.min():.2f}-{block.pH.max():.2f}",
            })
    return pd.DataFrame(rows)


def channel_summary(frame):
    """The same contrast over the whole 4OMe archive, matched on nothing."""
    four = frame[frame.substrate == "4OMe-BnOH"]
    out = {}
    for differential, block in four.groupby("differential"):
        name = "catalysed" if differential else "enzyme_free"
        out[name] = {
            "curves": int(len(block)),
            "median_depth": float(block.depth.median()),
            "deep": int((block.depth > 0.25).sum()),
            "accelerates": int(block.accelerates.sum()),
            "max_accel_z": float(block.accel_z.max()),
            "longest_s": float(block.duration_s.max()),
        }
    return out


# The one pair of runs in the block that holds the SCHEDULE fixed while the
# temperature moves: exps 19 and 14 were both recorded for 17934 s.
SCHEDULE_PAIR = (19, 14)


def schedule_control(frame, pair=SCHEDULE_PAIR):
    """
    Two runs of identical length whose induction times differ sevenfold.

    The induction's temperature dependence is the one between-run comparison
    this module makes, and the objection to it is that the operator gave the
    cold runs longer. These two did not differ: the same 17934 s at 15 C and at
    25 C. Whatever the run length does to `tau`, it does equally to both.
    """
    rows = frame[frame.experiment.isin(pair) & frame.tau_resolved]
    out = {}
    for experiment, block in rows.groupby("experiment"):
        out[int(experiment)] = {
            "temperature": float(block.temperature.iloc[0]),
            "duration_s": float(block.duration_s.iloc[0]),
            "curves": int(len(block)),
            "tau": float(block.tau.median()),
            "t_ind": float(block.t_ind.median()),
        }
    spans = {v["duration_s"] for v in out.values()}
    taus = [out[e]["tau"] for e in pair if e in out]
    out["same_span"] = len(spans) == 1
    out["tau_ratio"] = float(taus[0] / taus[1]) if len(taus) == 2 else np.nan
    return out


def activation_contrast():
    """
    The induction's activation parameters against the turnover's.

    Both come from `arrhenius.activation_parameters` and neither is refitted
    here. `inverse_tau` is the only parameter in that table whose Eyring
    entropy rests on no assumption about the rate law -- it is already a
    first-order relaxation rate in s^-1, so it needs neither an extinction
    coefficient nor an enzyme concentration -- which is what makes the entropy
    comparison below worth making at all.

    The gap that matters is not the enthalpy, which the two do not resolve
    apart, but the ENTROPY: a step that has to bring two solutes together in
    water pays 40-80 J/mol/K for it, and the induction pays nothing.
    """
    induction = arrhenius.activation_parameters("inverse_tau")
    turnover = arrhenius.activation_parameters("v_peak")
    entropy_gap = induction["entropy_J"] - turnover["entropy_J"]
    entropy_error = float(np.hypot(induction["entropy_stderr"],
                                   turnover["entropy_stderr"]))
    gibbs_gap = induction["gibbs_kJ"] - turnover["gibbs_kJ"]
    gibbs_error = float(np.hypot(induction["gibbs_stderr"],
                                 turnover["gibbs_stderr"]))
    return {
        "induction": induction, "turnover": turnover,
        "entropy_gap_J": float(entropy_gap),
        "entropy_gap_stderr": entropy_error,
        "gibbs_gap_kJ": float(gibbs_gap),
        "gibbs_gap_stderr": gibbs_error,
        # exp(-dG/RT): how many times faster the induction step is as a rate
        # constant at 298 K. A ratio of free energies, so the two prefactors
        # cancel and nothing about the rate law enters.
        "rate_ratio": float(np.exp(-gibbs_gap * 1000.0
                                   / (arrhenius.GAS_CONSTANT
                                      * arrhenius.REFERENCE_KELVIN))),
        "enthalpy_gap_kJ": float(induction["enthalpy_kJ"]
                                 - turnover["enthalpy_kJ"]),
        "enthalpy_gap_stderr": float(np.hypot(induction["enthalpy_stderr"],
                                              turnover["enthalpy_stderr"])),
    }


# The windows `landmark_window` compares, as a fraction of the run and in
# seconds. `curve_metrics.LAG_WINDOW` is the one everything else here uses.
WINDOW_SWEEP = (LAG_WINDOW, 0.05, 0.20)
ABSOLUTE_WINDOWS = (300.0, 900.0)
# The block whose schedule dependence is quoted: one temperature, so that the
# run length is the operator's choice and not the chemistry's.
SCHEDULE_TEMPERATURE = 25.0


def schedule_dependence(table, temperature=SCHEDULE_TEMPERATURE,
                        floor=INDUCTION_FLOOR):
    """
    How much of the induction time BETWEEN runs is the operator's schedule.

    The landmark's window is a tenth of the run, so a longer run reads its
    slope through a wider window and cannot resolve a short induction. This is
    the size of that problem, and it is why nothing in this module compares
    induction times across experiments.

    Fitted at one temperature, on the catalysed 4OMe runs, with pH beside run
    length because pH is the axis a between-run comparison would most want:
    it moves the rate by an order of magnitude across this block, and it is
    confounded with run length because the operator stopped the fast runs
    sooner.
    """
    block = table[(table.substrate == "4OMe-BnOH") & table.differential
                  & table.live & (table.temperature == temperature)
                  & np.isfinite(table.t_ind)]
    if len(block) < 8:
        return {"points": int(len(block))}
    response = np.log(np.maximum(block.t_ind.to_numpy(dtype=float), floor))
    design = np.column_stack([np.log(block.duration_s.to_numpy(dtype=float)),
                              block.pH.to_numpy(dtype=float),
                              np.ones(len(block))])
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    resid = response - design @ beta
    variance = float(resid @ resid) / max(1, len(block) - design.shape[1])
    covariance = variance * np.linalg.pinv(design.T @ design)
    spread = float(((response - response.mean()) ** 2).sum())
    return {"points": int(len(block)),
            "span": float(beta[0]),
            "span_stderr": float(np.sqrt(covariance[0, 0])),
            "pH": float(beta[1]),
            "pH_stderr": float(np.sqrt(covariance[1, 1])),
            "r2": float(1 - float(resid @ resid) / spread) if spread else np.nan,
            "experiments": int(block.experiment.nunique())}


def landmark_window(experiments=scope.TEMPERATURE_SERIES,
                    windows=WINDOW_SWEEP, absolute=ABSOLUTE_WINDOWS):
    """
    What the landmark and its activation energy do as the window changes.

    A window in SECONDS is the obvious way to make the statistic comparable
    between runs, and it does not survive contact with the cold end of this
    block: a 15 C run is 0.04 AU over five hours, so a 300 s window is reading
    a slope through noise, its maximum lands on an early excursion and the
    induction collapses. The activation energy goes with it. This is the
    measurement that decided `curve_metrics.LAG_WINDOW` -- a fraction of the
    run -- and with it the limitation stated at the top of this module.

    Returns one row per window: the geometric mean induction time at the
    coldest temperature, at the warmest, and the pooled activation energy of
    `1/t_ind` over all six temperatures with one offset per rung.
    """
    import pandas as pd
    frame = scope.frame(tuple(experiments))
    curves = list(scope.curves(tuple(experiments)))
    rows = []
    settings = ([(f"{w:.2f} of the run", w, None) for w in windows]
                + [(f"{a:.0f} s", None, a) for a in absolute])
    for label, fraction, seconds in settings:
        found = []
        for curve in curves:
            span = float(curve.times[-1] - curve.times[0])
            width = fraction if fraction is not None else min(seconds / span, 0.5)
            found.append({"experiment": curve.experiment,
                          "sample": curve.sample,
                          "t_ind": induction_point(curve, window=width).t_ind})
        table = frame.merge(pd.DataFrame(found), on=["experiment", "sample"])
        table = table[table.t_ind > 0].copy()
        table["rate_ind"] = 1.0 / table.t_ind
        fitted = arrhenius.pooled_arrhenius("rate_ind", frame=table,
                                            per_enzyme=False)
        cold = table[table.temperature == table.temperature.min()].t_ind
        warm = table[table.temperature == table.temperature.max()].t_ind
        rows.append({"window": label,
                     "curves": int(len(table)),
                     "cold_s": float(np.exp(np.log(cold).mean())),
                     "warm_s": float(np.exp(np.log(warm).mean())),
                     "activation_kJ": fitted["activation_kJ"],
                     "stderr_kJ": fitted["stderr_kJ"]})
    return pd.DataFrame(rows)


# The only cuvettes in the 4OMe archive that hold everything fixed and move
# the PEROXIDE: exps 127-131, 3.879 against 195.882 mM at five pH values, two
# cuvettes per level, 25 C, pyrophosphate, [S] fixed at 9.47 mM. Fifty-fold,
# and that is the whole lever this substrate has -- every other 4OMe run in the
# archive sits at 82.5 mM.
PEROXIDE_LEVER = (127, 128, 129, 130, 131)


def peroxide_lever(table, experiments=PEROXIDE_LEVER):
    """
    Does more peroxide shorten the induction, as forming an adduct with it must?

    `K + H2O2 <=> KP` relaxes at `1/tau = k_f[H2O2] + k_r`, so if what the
    induction times is the catalyst binding the oxidant, the induction time
    has an order in [H2O2] between 0 and -1 and can have no other sign. Any
    positive order falsifies that reading whatever else it means.

    The block is small and its cuvettes are not all live, so the answer comes
    with the signal-to-noise control beside it: the landmark is read off a
    rolling slope, a curve with no signal has a noisy one, and [H2O2] sets the
    signal. `signal_control` is the same test for the block this module's
    conclusion actually rests on.
    """
    block = table[table.experiment.isin(experiments)]
    induction_order = scope.orders("t_ind", frame=block, floor=INDUCTION_FLOOR)
    rate_order = scope.orders("v_peak", frame=block)
    return {"curves": int(len(block)),
            "live": int(block.live.sum()),
            "levels": sorted(float(v) for v in block.h2o2.unique()),
            "induction_order": induction_order["order_h2o2"],
            "induction_stderr": induction_order["stderr_h2o2"],
            "rate_order": rate_order["order_h2o2"],
            "rate_stderr": rate_order["stderr_h2o2"],
            "points": induction_order["n"],
            **signal_control(block)}


# ---------------------------------------------------------------------------
# What a peroxide adduct would require of the two orders at once.
#
# For `K + H2O2 <=> KP` with h = [H2O2] in 100-6000x excess over the catalyst,
# so that the forward leg is pseudo-first-order:
#
#     approach       1/tau = k_f h + k_r
#     destination    [KP]/E0 = K h / (1 + K h),      K = k_f / k_r
#
# and therefore, writing the log-log slopes in h,
#
#     d ln v   / d ln h = 1 / (1 + K h)        +1 unsaturated, 0 saturated
#     d ln tau / d ln h = -K h / (1 + K h)      0 unsaturated, -1 saturated
#
# THE TWO ARE LOCKED. Their difference is 1 identically, for every K and every
# h, so "is the rate first order in peroxide" and "does peroxide shorten the
# induction" are one question asked twice. That is what makes it testable
# without knowing where on the saturation curve the design sits, and it is
# tested here as ONE regression on log(v / t_ind), which also disposes of the
# correlation between two coefficients fitted to the same 15 curves.
#
# The other half of the scheme is worth stating separately: `1/tau` is
# monotonically increasing in h whatever the constants are, so a POSITIVE order
# on the induction time is not a statement about the regime -- it falsifies the
# scheme outright. `peroxide_lever` measures that order; this measures the
# constraint it violates.
PERHYDRATE_ORDER_GAP = 1.0

# Where a trap puts the same slope instead. If H2O2 parks the catalyst in a
# state that is NOT on the activation path, only the free fraction
# 1/(1 + K h) can activate, so 1/tau_obs = k_act/(1 + K h) and
#
#     d ln tau / d ln h = +K h / (1 + K h),   which lies in (0, +1)
#
# -- positive, bounded, and consistent with everything section 3 establishes,
# because k_act carries no concentration at all.
TRAP_ORDER_RANGE = (0.0, 1.0)

# The profile grid for the association constant, in 1/mM, and the F cutoff
# that turns a profile into a 95% interval. F(1, dof, 0.95) is 4.08 at dof 40
# and 3.94 at dof 200; 3.99 is the middle of the range this module fits in and
# the interval is not sensitive to the choice at three digits.
SATURATION_SPAN = (-4.0, 0.5)
SATURATION_GRID = 600
PROFILE_F = 3.99


def joint_peroxide_order(table, floor=INDUCTION_FLOOR):
    """
    Regress log(rate / induction time) on log[H2O2]: the scheme requires +1.

    One regression rather than two, because the two coefficients it replaces
    are fitted to the same curves on the same design and their errors are
    correlated; differencing them by hand would quote an error that is not
    theirs. A `log[S]` term joins the fit wherever the block moves the
    substrate as well, since both sides of the response carry a substrate
    order and the peroxide ladder is an L rather than a grid.

    Returns the coefficient, its error, and how many standard errors it sits
    from the +1 the adduct scheme demands.

    THE FLOOR MOVES THIS ONE, and which way depends on where the short
    inductions sit. `test_induction` plants an adduct fast enough to be clipped
    and shows the bias is downward there -- towards this module's own reading --
    so the sweep is printed by `report` rather than left to be assumed away. In
    this archive the short inductions are at LOW peroxide, so the clipping
    pushes the coefficient the other way, up towards +1, and the module's 60 s
    floor is the conservative end of the range rather than the flattering one.
    """
    live = table[table.live & np.isfinite(table.t_ind)
                 & (table.v_peak > 0) & (table.h2o2 > 0)].copy()
    if len(live) < 8 or live.h2o2.nunique() < 2:
        return {"points": int(len(live))}
    response = (np.log(live.v_peak.to_numpy(dtype=float))
                - np.log(np.maximum(live.t_ind.to_numpy(dtype=float), floor)))
    columns = [np.log(live.h2o2.to_numpy(dtype=float))]
    if live.s0.nunique() > 1:
        columns.append(np.log(live.s0.to_numpy(dtype=float)))
    for experiment in sorted(live.experiment.unique()):
        columns.append((live.experiment == experiment).to_numpy(float))
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    resid = response - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(live) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    slope = float(beta[0])
    stderr = float(np.sqrt(max(covariance[0, 0], 0.0)))
    return {"points": int(len(live)),
            "experiments": int(live.experiment.nunique()),
            "slope": slope, "stderr": stderr,
            "required": PERHYDRATE_ORDER_GAP,
            "sigma": float(abs(PERHYDRATE_ORDER_GAP - slope) / stderr)
            if stderr else np.nan}


def peroxide_ladder(table):
    """
    The cuvettes that move [H2O2] with everything else held, one row each.

    In exps 135-151 the seven cuvettes are an L: four step the substrate at the
    top peroxide, four step the peroxide at the top substrate. Only the second
    arm is a peroxide ladder, and taking the whole run instead puts the
    substrate ladder into a fit that has no term for it.
    """
    live = table[table.live & (table.vmax > 0) & (table.h2o2 > 0)]
    if not len(live):
        return live
    top = live.groupby("experiment").s0.transform("max")
    arm = live[np.isclose(live.s0, top)]
    counts = arm.groupby("experiment").h2o2.transform("nunique")
    return arm[counts >= 2]


def peroxide_saturation(table, parameter="vmax", grid=SATURATION_GRID,
                        span=SATURATION_SPAN, cutoff=PROFILE_F):
    """
    Is the rate first order in peroxide, and can one binding equilibrium do it?

    Two fits on the ladder, each with one free level per run so that pH, buffer,
    enzyme, cell and day are absorbed:

        free power      v ~ h^a                       `a` says what the order is
        the scheme      v ~ K h / (1 + K h)           exponent FIXED at 1, K
                                                      profiled on a log grid

    The exponent is fixed at 1 on purpose. Leaving it free lets the saturating
    form buy a fit with an exponent above 1 and a large K, which is a curve
    through the points and not the hypothesis: `v` is proportional to [KP], and
    [KP] is that expression with no exponent on it.

    Returns the free exponent with its error, the F statistic against a strict
    first order, and the profiled K with a 95% interval and the perhydrate
    fraction it implies at the archive's working 82.5 mM.
    """
    ladder = peroxide_ladder(table)
    if len(ladder) < 10:
        return {"points": int(len(ladder))}
    h = ladder.h2o2.to_numpy(dtype=float)
    y = np.log(ladder[parameter].to_numpy(dtype=float))
    dummies = np.column_stack([(ladder.experiment.to_numpy() == e).astype(float)
                               for e in sorted(ladder.experiment.unique())])
    degrees = max(1, len(y) - dummies.shape[1] - 1)

    def _levelled(offset):
        """Best SSE once each run is free to sit where it likes."""
        beta, *_ = np.linalg.lstsq(dummies, y - offset, rcond=None)
        resid = y - offset - dummies @ beta
        return float(resid @ resid)

    powers = np.linspace(0.0, 1.6, grid)
    power_sse = np.array([_levelled(a * np.log(h)) for a in powers])
    best = int(np.argmin(power_sse))
    # The profile interval on the exponent, read the same way as K's below.
    inside = powers[power_sse <= power_sse[best] * (1.0 + cutoff / degrees)]
    first_order = _levelled(np.log(h))

    constants = np.concatenate([[0.0], np.logspace(*span, grid)])
    scheme_sse = np.array([_levelled(np.log(k * h / (1.0 + k * h)))
                           if k > 0 else _levelled(np.log(h))
                           for k in constants])
    top = int(np.argmin(scheme_sse))
    allowed = constants[scheme_sse
                        <= scheme_sse[top] * (1.0 + cutoff / degrees)]
    working = 82.5
    return {
        "points": int(len(ladder)),
        "experiments": int(ladder.experiment.nunique()),
        "peroxide_low": float(h.min()), "peroxide_high": float(h.max()),
        "order": float(powers[best]),
        "order_low": float(inside.min()), "order_high": float(inside.max()),
        "first_order_f": float((first_order - power_sse[best])
                               / (power_sse[best] / degrees)),
        "first_order_sse": float(first_order),
        "power_sse": float(power_sse[best]),
        "scheme_sse": float(scheme_sse[top]),
        "constant": float(constants[top]),
        "constant_low": float(allowed.min()),
        "constant_high": float(allowed.max()),
        "bound_fraction": (float(constants[top] * working
                                 / (1.0 + constants[top] * working))),
        "bound_low": float(allowed.min() * working
                           / (1.0 + allowed.min() * working)),
        "bound_high": float(allowed.max() * working
                            / (1.0 + allowed.max() * working)),
        "working_mM": working,
    }


def trap_constant(order, stderr, peroxide):
    """
    Invert `d ln tau / d ln h = K h / (1 + K h)` for K, and turn it into a ΔG°.

    Only meaningful for an order inside `TRAP_ORDER_RANGE`; an order at or above
    1 has no finite K and one at or below 0 has no trap. The error is the delta
    method, dK/d(order) = 1 / (h (1 - order)^2), which blows up as the order
    approaches 1 -- correctly, because there the data stop constraining K.

    `peroxide` is the geometric-mean [H2O2] of the design the order was measured
    on: the order is a LOCAL log-log slope and K is only recoverable at the
    concentration it was measured at.

    ΔG° is quoted per mole with the association constant in M^-1, which is a
    thousand times the mM^-1 this package works in.
    """
    low, high = TRAP_ORDER_RANGE
    if not low < order < high or peroxide <= 0:
        return {"constant": np.nan, "stderr": np.nan,
                "free_energy_kJ": np.nan, "peroxide_mM": float(peroxide)}
    constant = order / (peroxide * (1.0 - order))
    error = stderr / (peroxide * (1.0 - order) ** 2)
    molar = constant * 1000.0
    return {"constant": float(constant), "stderr": float(error),
            "molar": float(molar),
            "free_energy_kJ": float(-arrhenius.GAS_CONSTANT
                                    * arrhenius.REFERENCE_KELVIN
                                    * np.log(molar) / 1000.0),
            "bound_fraction": float(constant * peroxide
                                    / (1.0 + constant * peroxide)),
            "peroxide_mM": float(peroxide)}


def peroxide_geometric_mean(table):
    """The [H2O2] a local order measured on this block belongs to."""
    live = table[table.live & (table.h2o2 > 0)]
    return float(np.exp(np.log(live.h2o2.to_numpy(dtype=float)).mean()))


def signal_control(table, floor=INDUCTION_FLOOR):
    """
    Regress the induction time on the curve's signal-to-noise, same design.

    The one artefact that could manufacture every result in this module: the
    landmark is the first crossing of half the LARGEST rolling slope, and on a
    curve with no signal the largest rolling slope is a noise excursion that
    can land anywhere, usually early. If `t_ind` tracked signal-to-noise the
    orders would be measuring the spectrophotometer.

    Signal-to-noise is `net/noise`, both of them `scope`'s own columns.
    """
    live = table[table.live & np.isfinite(table.t_ind)
                 & (table.net > 0) & (table.noise > 0)]
    if len(live) < 8:
        return {"signal_slope": np.nan, "signal_stderr": np.nan,
                "signal_points": int(len(live))}
    response = np.log(np.maximum(live.t_ind.to_numpy(dtype=float), floor))
    columns = [np.log((live.net / live.noise).to_numpy(dtype=float))]
    for experiment in sorted(live.experiment.unique()):
        columns.append((live.experiment == experiment).to_numpy(float))
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    resid = response - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(live) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    return {"signal_slope": float(beta[0]),
            "signal_stderr": float(np.sqrt(max(covariance[0, 0], 0.0))),
            "signal_points": int(len(live))}


# ---------------------------------------------------------------------------
# WHICH WAY THE INDUCTION POINTS.
#
# Everything above measures how LONG the induction is on curves that have one.
# The archive also contains curves that begin FAST and slow down, and until
# `progress_kind` reached `scope.frame` on 2026-09-02 nothing here could tell
# the two apart: `depth` is 1 - start/peak and cannot go below zero, so a curve
# that is fastest in its first window and one that starts a hair below its peak
# both read as "no induction".
#
# The sign comes from the form the curve earned, on the convention B > 0 means
# the rate starts BELOW its eventual value. The two-phase form has carried both
# signs since it was written -- `B1 > 0 > B2` is a lag then a fall, `B1 < 0` a
# burst -- so nothing needed refitting; what was missing was reading it.
LAG_FIRST_KINDS = ("lag", "lag then fall", "two lags")


def sign_table(frame):
    """
    The frame with `lag_first`: does this curve begin below its eventual rate?

    A BINARY, and deliberately. The obvious continuous statistic is the fast
    phase's amplitude, `B_fast`, normalised by the run's own signal -- and it
    is unusable. When the two exponentials are nearly degenerate the linear
    solve trades enormous opposite amplitudes between them (exp 135 sample 3
    comes back with B_fast = -241 and B_slow = +303 on a curve that moves 0.06
    AU), so `B_fast/net` has an interquartile range of 1.8 and a 10-90 range of
    35. The SIGN of that trade is stable where its size is not.
    """
    table = frame.copy()
    table["lag_first"] = table.progress_kind.isin(LAG_FIRST_KINDS)
    return table


def ladder_arms(table):
    """
    Split exps 135-151's L into the two arms that each move one thing.

    The seven cuvettes of an two-axis run are not a grid: four step the
    substrate at the run's top peroxide, and four step the peroxide at the run's
    top substrate, sharing the corner. Pooled, `log[S]` and `log[H2O2]`
    correlate about -0.5 by construction, so a banded table against either axis
    is reading the other one too. Splitting is not optional here -- the pooled
    bands say the sign tracks peroxide at 10% against 56%, and inside the
    peroxide arm alone that effect is gone.
    """
    live = table[table.live]
    top_substrate = live.groupby("experiment").s0.transform("max")
    top_peroxide = live.groupby("experiment").h2o2.transform("max")
    return {"substrate arm": live[np.isclose(live.h2o2, top_peroxide)],
            "peroxide arm": live[np.isclose(live.s0, top_substrate)]}


def sign_drivers(table, axis="s0", control=True):
    """
    Does the sign of the induction track `axis`, within runs?

    A LINEAR PROBABILITY MODEL: least squares of the 0/1 `lag_first` on
    `log(axis)` with one offset per experiment. Not a logit, and the reason is
    the design rather than taste -- 63 rows carrying 17 run offsets separate
    perfectly in several runs, where a logit's coefficient runs to infinity and
    its standard error with it. What is wanted here is a direction and a rough
    size, and for that the linear model is honest as long as it is read as one:
    the coefficient is a change in probability per e-fold of `axis`, and it is
    not constrained to keep the fitted values inside [0, 1].

    `control=True` adds each curve's own signal-to-noise. It has to be there:
    on a curve with little signal the two-phase fit has little to choose
    between the shapes, and the sign it returns leans burst. In-scope the share
    of lag-first curves rises from 0.29 to 0.50 across the signal-to-noise
    quartiles. What makes the substrate result below worth reporting is that it
    runs AGAINST that lean rather than with it.
    """
    live = table[table.live & np.isfinite(table[axis]) & (table[axis] > 0)
                 & (table.net > 0) & (table.noise > 0)].copy()
    if len(live) < 8 or live[axis].nunique() < 2:
        return {"points": int(len(live))}
    response = live.lag_first.to_numpy(dtype=float)
    columns = [np.log(live[axis].to_numpy(dtype=float))]
    names = [axis]
    if control:
        columns.append(np.log((live.net / live.noise).to_numpy(dtype=float)))
        names.append("signal")
    for experiment in sorted(live.experiment.unique()):
        columns.append((live.experiment == experiment).to_numpy(float))
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    resid = response - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(live) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    out = {"points": int(len(live)),
           "experiments": int(live.experiment.nunique()),
           "lag_first": float(live.lag_first.mean())}
    for index, name in enumerate(names):
        out[name] = float(beta[index])
        out[name + "_stderr"] = float(np.sqrt(max(covariance[index, index],
                                                  0.0)))
    return out


def composition_collinearity(table):
    """
    How far `[S]` and `[buf]` move together inside a run, per block.

    THE CAVEAT THIS FUNCTION EXISTS FOR. In every 4OMe run the substrate was
    added by volume and displaced buffer, so `[buf]` falls 80 -> 50 mM as `[S]`
    rises 1.85 -> 7.4: the two correlate at -0.96 in logs. That is already
    known to corrupt the substrate order of the RATE -- `temperature_series`
    3 corrects it with the buffer order measured on exps 32 and 34 -- and it
    corrupts the substrate order of the INDUCTION in exactly the same way.

    Exps 135-151 are the block where it does not: `[buf]` is constant across all
    seven cuvettes of all seventeen runs, so there `[S]` moves alone.

    Which is why the two blocks may disagree about the substrate without either
    being wrong, and why `induction_drivers` -- whose regressor is the measured
    RATE and not a concentration -- is the route in section 3 that this does
    not reach.
    """
    rows = []
    for experiment, block in table[table.live].groupby("experiment"):
        if block.s0.nunique() < 2:
            continue
        if block.buf.nunique() < 2:
            rows.append(0.0)
            continue
        rows.append(float(np.corrcoef(np.log(block.s0.to_numpy(dtype=float)),
                                      np.log(block.buf.to_numpy(dtype=float)))
                          [0, 1]))
    if not rows:
        return {"runs": 0}
    # And the SLOPE, not just the correlation: correcting an order needs
    # d log[buf] / d log[S], which is what a fit without a buffer term folds
    # into the substrate coefficient. `temperature_series` 3 uses -0.325 for
    # its own six runs and corrects the RATE's order with it.
    live = table[table.live & (table.s0 > 0) & (table.buf > 0)]
    laddered = live.groupby("experiment").s0.transform("nunique") > 1
    live = live[laddered]
    slope = np.nan
    if len(live) > 4 and live.buf.nunique() > 1:
        x = np.log(live.s0.to_numpy(dtype=float))
        y = np.log(live.buf.to_numpy(dtype=float))
        design = np.column_stack(
            [x] + [(live.experiment.to_numpy() == e).astype(float)
                   for e in sorted(live.experiment.unique())])
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        slope = float(beta[0])
    return {"runs": len(rows),
            "median": float(np.median(rows)),
            "slope": slope,
            "constant_buffer": int(sum(1 for value in rows if value == 0.0))}


def substrate_order_corrected(table, buffer_slope, buffer_stderr,
                              response="t_ind", floor=INDUCTION_FLOOR):
    """
    The substrate order with the buffer that rides on the ladder taken out.

    A fit without a buffer term returns `a' = a + d.g`, where `d` is the order
    in `[buf]` and `g = d log[buf] / d log[S]` is how steeply buffer rides on
    the substrate ladder. `temperature_series` 3 does exactly this to the
    RATE's substrate order, with d = +0.400 from exps 32 and 34; this does it
    to the INDUCTION's, with d from `buffer_order` on the same two runs.

    IT MATTERS, and not in the direction that flatters this module. The
    uncorrected order is -0.121 +- 0.148 against the -0.471 a product threshold
    requires, which excludes it; corrected it moves towards the threshold and
    stops excluding it. What is NOT reached by any of this is
    `induction_drivers`, whose regressor is the curve's own measured rate: it
    asks whether a faster cuvette's induction is shorter, and that question is
    well posed whatever is making the cuvette faster -- buffer included.
    """
    measured = scope.orders(response, frame=table, floor=floor)
    ladder = composition_collinearity(table)
    slope = ladder.get("slope", np.nan)
    correction = buffer_slope * slope
    corrected = measured["order_s0"] - correction
    rate = scope.orders("v_peak", frame=table)
    return {"measured": measured["order_s0"],
            "measured_stderr": measured["stderr_s0"],
            "buffer_slope": float(buffer_slope),
            "ladder_slope": float(slope),
            "correction": float(correction),
            "correction_stderr": float(abs(slope) * buffer_stderr),
            "corrected": float(corrected),
            "corrected_stderr": float(np.hypot(measured["stderr_s0"],
                                               abs(slope) * buffer_stderr)),
            # What a product threshold would require: the rate's own substrate
            # order with the sign reversed.
            "threshold": float(-rate["order_s0"]),
            "points": measured["n"]}


# The archive's one direct buffer lever on an induction: exps 32 and 34, 4OMe
# at 40 C and pH 7.00 with [S] at 8.251 mM and [H2O2] at 82.5 mM on all eight
# cuvettes, and [buf] stepped 3.125 -> 200 mM. Sixty-four fold, and two runs.
BUFFER_LEVER = (34, 32)

# The window this lever's landmark is read through, in SECONDS and not as a
# fraction of the run. Everywhere else in this module the window is a fraction,
# for the reason `landmark_window` measures: a cold slow curve needs a wide one.
# Here that rule inverts, because the two runs differ threefold in length
# (5280 s against 1767 s) and a fractional window would read them through
# windows that differ threefold too -- which is the whole reason the two runs
# looked like they disagreed. These are 40 C curves reaching 0.05-0.31 AU, so
# they can afford a fixed window where the 15 C curves cannot.
BUFFER_WINDOW = 450.0
BUFFER_WINDOW_SWEEP = (300.0, 450.0, 600.0, 900.0, 1200.0)


def buffer_landmark(curve, width=BUFFER_WINDOW):
    """`induction_point` with the window given in seconds rather than as a share."""
    span = float(curve.times[-1] - curve.times[0])
    return induction_point(curve, window=min(width / span, 0.5))


def buffer_lever(table, experiments=BUFFER_LEVER, width=BUFFER_WINDOW,
                 floor=INDUCTION_FLOOR):
    """
    Does buffer shorten the induction, as general base catalysis of E -> E* would?

    WHAT THIS FUNCTION GOT WRONG FIRST, because the mistake is the reason it is
    written this way. Its first version regressed `tau_fast` on `[buf]` inside
    each run and reported that the two runs disagreed in sign, +0.457 +- 0.097
    against -1.052 +- 0.469. They do not disagree. `tau_fast` is not the same
    quantity in the two runs: every curve of exp 34 earns the two-phase form
    (F = 71 to 819 against a threshold of 12) and every curve of exp 32 earns
    the one-phase form (F = 1.6 to 4.7), so one row was tau1 of a two-phase fit
    and the other tau of a one-phase fit. And the model form breaks exactly at
    the run boundary because the SCHEDULE does: exp 34 ran 5280 s and exp 32
    1767 s, so exp 34's runs are long enough to contain the slow fall and exp
    32's end before it. A comparison across that boundary was never between two
    measurements of one thing.

    What is comparable is a landmark read off the READINGS through a window
    common to both runs, which is what `buffer_landmark` does. On that footing
    the two runs agree, and they agree in the direction general base catalysis
    predicts.

    ONE LEVEL PER RUN, always. The rate -- `v_peak`, the one quantity defined
    identically on both forms -- FALLS 1.80x across the join, from 7.41e-5 at
    25 mM in exp 34 to 4.13e-5 at 50 mM in exp 32, while the buffer doubles.
    Whatever separates the two days is larger than the effect being measured,
    so a fit with a shared intercept is measuring the day. That is the same
    reason `temperature_series` 3 quotes the buffer order of the RATE as two
    range-specific numbers, +0.803 +- 0.173 below 25 mM and +0.400 +- 0.028
    above 50, rather than one pooled one.
    """
    import pandas as pd
    rows = []
    for curve in scope.curves(tuple(experiments)):
        found = buffer_landmark(curve, width)
        rows.append({"experiment": curve.experiment, "sample": curve.sample,
                     "buf": curve.buf, "t_ind": found.t_ind,
                     "depth": found.depth,
                     "span_s": float(curve.times[-1] - curve.times[0])})
    ladder = pd.DataFrame(rows)
    # `live`, `net` and `noise` come along so that `signal_control` runs on
    # THIS ladder unchanged -- the landmark here is the seconds-window one, and
    # a control read off the fractional-window `t_ind` would not be its control.
    joined = ladder.merge(table[["experiment", "sample", "phases",
                                 "two_phase_f", "v_peak", "progress_kind",
                                 "live", "net", "noise"]],
                          on=["experiment", "sample"])
    return joined


def buffer_order(ladder, response="t_ind", floor=INDUCTION_FLOOR):
    """
    The pooled log-log slope of the landmark on `[buf]`, one level per run.

    Returns the pooled slope and each run's own, so that "the two runs agree"
    is a statement the caller can check rather than one it has to believe.
    """
    values = ladder[response].to_numpy(dtype=float)
    y = np.log(np.maximum(values, floor if response == "t_ind" else DEPTH_FLOOR))
    x = np.log(ladder.buf.to_numpy(dtype=float))
    runs = sorted(ladder.experiment.unique())
    design = np.column_stack([x] + [(ladder.experiment.to_numpy() == run)
                                    .astype(float) for run in runs])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(y) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    out = {"points": int(len(y)), "runs": len(runs),
           "slope": float(beta[0]),
           "stderr": float(np.sqrt(max(covariance[0, 0], 0.0)))}
    for run in runs:
        block = ladder[ladder.experiment == run]
        if len(block) < 3:
            continue
        bx = np.log(block.buf.to_numpy(dtype=float))
        by = np.log(np.maximum(block[response].to_numpy(dtype=float),
                               floor if response == "t_ind" else DEPTH_FLOOR))
        single = np.column_stack([bx, np.ones(len(bx))])
        coefficients, *_ = np.linalg.lstsq(single, by, rcond=None)
        residual = by - single @ coefficients
        spread = float(residual @ residual) / max(1, len(bx) - 2)
        error = spread * np.linalg.pinv(single.T @ single)
        out[f"slope_{int(run)}"] = float(coefficients[0])
        out[f"stderr_{int(run)}"] = float(np.sqrt(max(error[0, 0], 0.0)))
    return out


def joint_buffer_order(ladder, floor=INDUCTION_FLOOR):
    """
    The same pre-equilibrium constraint as `joint_peroxide_order`, on `[buf]`.

    THE ALGEBRA DOES NOT CARE WHICH SPECIES IT IS. Section 4b's `+1` is not a
    fact about hydrogen peroxide; it is a fact about any scheme in which the
    catalyst is drawn into its active form by a species held in excess. Write
    that species X, at a concentration the run controls:

        E + X <-> E*        1/tau = k_f [X] + k_r,   [E*]/E0 = K[X]/(1 + K[X])

    so d ln v/d ln[X] = 1/(1 + K[X]) and d ln tau/d ln[X] = -K[X]/(1 + K[X]),
    and the DIFFERENCE is 1 for every K and every [X]. Put the buffer in that
    role -- as a general base, or through a peroxo adduct of the buffer itself,
    the algebra is identical -- and the constraint transfers unchanged.

    This is worth doing because the peroxide axis FAILS it (2.6 and 3.7 sigma)
    on blocks that also fail `signal_control`, and the buffer axis is the only
    other lever in the archive that moves a candidate activator. It is eight
    curves, so it is a weak test; it is a weak test of a parameter-free
    prediction, which is not the same as no test.

    Read on the same footing `buffer_lever` establishes: `t_ind` from a
    landmark whose window is in SECONDS, because exps 32 and 34 differ in
    length, and one free level per run, because the rate falls 1.80x across the
    join. One regression rather than two differenced by hand, so the covariance
    between the two orders is the fit's own.
    """
    live = ladder[np.isfinite(ladder.t_ind) & (ladder.v_peak > 0)
                  & (ladder.buf > 0)]
    if len(live) < 4 or live.buf.nunique() < 2:
        return {"points": int(len(live))}
    response = (np.log(live.v_peak.to_numpy(dtype=float))
                - np.log(np.maximum(live.t_ind.to_numpy(dtype=float), floor)))
    columns = [np.log(live.buf.to_numpy(dtype=float))]
    runs = sorted(live.experiment.unique())
    columns += [(live.experiment.to_numpy() == run).astype(float)
                for run in runs]
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    resid = response - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(live) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    slope = float(beta[0])
    stderr = float(np.sqrt(max(covariance[0, 0], 0.0)))
    return {"points": int(len(live)), "runs": len(runs),
            "slope": slope, "stderr": stderr,
            "required": PERHYDRATE_ORDER_GAP,
            "sigma": float(abs(PERHYDRATE_ORDER_GAP - slope) / stderr)
            if stderr else np.nan}


def buffer_join_step(ladder):
    """
    The between-run level step, which is why the two runs get separate levels.

    The top rung of the low-buffer run and the bottom rung of the high-buffer
    one differ twofold in `[buf]`, so the rate should RISE across them and it
    falls instead. That ratio is the size of whatever separates the two days.
    """
    low, high = ladder[ladder.experiment == BUFFER_LEVER[0]], \
        ladder[ladder.experiment == BUFFER_LEVER[1]]
    if not len(low) or not len(high):
        return {}
    top = low.loc[low.buf.idxmax()]
    bottom = high.loc[high.buf.idxmin()]
    return {"from_buf": float(top.buf), "to_buf": float(bottom.buf),
            "from_rate": float(top.v_peak), "to_rate": float(bottom.v_peak),
            "step": float(top.v_peak / bottom.v_peak)}


# The cuts this question needs, by the frame's own columns. `slowdown` has the
# same idea for the fall and a different set of blocks: this one separates the
# 4OMe archive by CHANNEL, because the catalyst is the variable under test, and
# keeps the BnOH scope whole, because there the question is whether the early
# curve is the same object at all.
def induction_blocks(frame):
    """The named cuts, from the frame's own columns."""
    four = frame.substrate == "4OMe-BnOH"
    return {
        "4OMe catalysed": frame[four & frame.differential],
        "4OMe catalysed, 25 C": frame[four & frame.differential
                                      & (frame.temperature == 25.0)],
        "4OMe enzyme-free": frame[four & ~frame.differential],
        "temperature series":
            frame[frame.experiment.isin(scope.TEMPERATURE_SERIES)],
        "BnOH two-axis (135-151)":
            frame[frame.experiment.isin(scope.TWO_AXIS_BLOCK)],
    }


WHOLE_ARCHIVE = tuple(range(1, 152))


def report(table=None):
    """Print the whole argument, in the order it has to be made."""
    if table is None:
        table = induction_table(WHOLE_ARCHIVE)
    named = induction_blocks(table)

    print("\n1. DOES THE INDUCTION NEED THE CATALYST")
    summary = channel_summary(table)
    for name in ("catalysed", "enzyme_free"):
        row = summary[name]
        print(f"   4OMe {name:11s} {row['curves']:3d} curves   "
              f"median depth {row['median_depth']:.3f}   "
              f"deep {row['deep']:3d}   accelerating {row['accelerates']:3d}   "
              f"max z {row['max_accel_z']:.2f}   "
              f"longest {row['longest_s']:.0f} s")
    print()
    print(channel_contrast(table).to_string(index=False))

    print("\n2. IS IT THE SCHEDULE")
    control = schedule_control(table)
    for experiment in SCHEDULE_PAIR:
        row = control[experiment]
        print(f"   exp {experiment}  {row['temperature']:.0f} C   "
              f"{row['duration_s']:.0f} s   tau {row['tau']:.0f} s   "
              f"t_ind {row['t_ind']:.0f} s")
    print(f"   same span: {control['same_span']}   "
          f"tau ratio {control['tau_ratio']:.1f}x")

    print("\n2a. WHAT THE STATISTIC OWES TO ITS WINDOW AND TO THE SCHEDULE")
    print(landmark_window().to_string(index=False))
    drift = schedule_dependence(table)
    print(f"   between runs at {SCHEDULE_TEMPERATURE:.0f} C, "
          f"{drift['points']} curves in {drift['experiments']} experiments:")
    print(f"   t_ind on run length {drift['span']:+.3f} +- "
          f"{drift['span_stderr']:.3f}   on pH {drift['pH']:+.3f} +- "
          f"{drift['pH_stderr']:.3f}")

    print("\n3. A CLOCK OR A PRODUCT THRESHOLD")
    lever = substrate_lever(named["4OMe catalysed"])
    print(f"   the lever: [S] moves {lever['median_lever']:.1f}x inside the "
          f"median run, on {lever['laddered']} of {lever['experiments']} "
          f"experiments")
    print("   log(induction time) on log(rate), one offset per experiment.")
    print(f"   clock predicts {INDUCTION_CLOCK_SLOPE:+.1f}, "
          f"product control {INDUCTION_PRODUCT_SLOPE:+.1f}")
    for name in ("4OMe catalysed", "temperature series",
                 "BnOH two-axis (135-151)"):
        block = named[name]
        for rate in ("peak_rate", "v_peak", "vmax"):
            fit = induction_drivers(block, rate=rate)
            if "slope" not in fit:
                continue
            print(f"   {name:26s} vs {rate:10s} "
                  f"{fit['slope']:+.3f} +- {fit['stderr']:.3f}   "
                  f"n={fit['points']:3d}  floored {fit['floored']:2d}")
    print("   and the same fit at other floors, to show it is not the floor:")
    for floor in FLOOR_SWEEP:
        fit = induction_drivers(named["4OMe catalysed"], floor=floor)
        print(f"   4OMe catalysed, floor {floor:5.0f} s        "
              f"{fit['slope']:+.3f} +- {fit['stderr']:.3f}")

    print("\n4. THE SAME QUESTION THROUGH THE SUBSTRATE ORDERS")
    print("   the route with no errors-in-variables: the regressor is the")
    print("   composition, not a rate measured off the same curve.")
    for name in ("4OMe catalysed", "temperature series"):
        ratio = order_ratio(named[name])
        print(f"   {name:22s} order(t_ind) "
              f"{ratio['induction_order']:+.3f} +- {ratio['induction_stderr']:.3f}"
              f"   order(v_peak) {ratio['rate_order']:+.3f} +- "
              f"{ratio['rate_stderr']:.3f}"
              f"   ratio {ratio['ratio']:+.2f} +- {ratio['ratio_stderr']:.2f}")
    print("   and the amplitude, which a product build-up would make "
          "composition dependent:")
    for name in ("4OMe catalysed", "temperature series"):
        depth = scope.orders("depth", frame=named[name], floor=DEPTH_FLOOR)
        print(f"   {name:22s} order(depth) in [S] "
              f"{depth['order_s0']:+.3f} +- {depth['stderr_s0']:.3f}  "
              f"n={depth['n']}")

    print("\n4a. THE ONE PEROXIDE LEVER THIS SUBSTRATE HAS")
    lever = peroxide_lever(table)
    print(f"   exps {PEROXIDE_LEVER} -- {lever['levels']} mM, "
          f"{lever['curves']} curves, {lever['live']} live")
    print(f"   order of t_ind in [H2O2]  "
          f"{lever['induction_order']:+.3f} +- {lever['induction_stderr']:.3f}"
          f"   (an adduct with H2O2 requires 0 to -1)")
    print(f"   order of v_peak in [H2O2] "
          f"{lever['rate_order']:+.3f} +- {lever['rate_stderr']:.3f}")
    print(f"   t_ind on signal-to-noise  "
          f"{lever['signal_slope']:+.3f} +- {lever['signal_stderr']:.3f}"
          f"   n={lever['signal_points']}")
    print("   the same control on the block the conclusion rests on:")
    for name in ("4OMe catalysed", "BnOH two-axis (135-151)"):
        got = signal_control(named[name])
        print(f"   {name:26s} {got['signal_slope']:+.3f} "
              f"+- {got['signal_stderr']:.3f}   n={got['signal_points']}")

    print("\n4b. WHAT AN ADDUCT WITH H2O2 WOULD REQUIRE OF BOTH ORDERS AT ONCE")
    print("   d ln v/d ln h - d ln tau/d ln h = 1 identically, for every K and")
    print("   every h, so the two questions are one. Measured as one fit:")
    peroxide_blocks = (("4OMe peroxide, exps 127-131",
                        table[table.experiment.isin(PEROXIDE_LEVER)]),
                       ("BnOH two-axis, exps 135-151",
                        named["BnOH two-axis (135-151)"]))
    for label, block in peroxide_blocks:
        joint = joint_peroxide_order(block)
        if "slope" not in joint:
            continue
        print(f"   {label:30s} order(v/t_ind) {joint['slope']:+.3f} +- "
              f"{joint['stderr']:.3f}   n={joint['points']:3d}   "
              f"{joint['sigma']:.1f} sigma from the required "
              f"{joint['required']:+.0f}")
    print("   and at the other floors, because a floor moves this one:")
    for label, block in peroxide_blocks:
        swept = [joint_peroxide_order(block, floor=floor)
                 for floor in FLOOR_SWEEP]
        print(f"   {label:30s} "
              + "  ".join(f"{floor:.0f}s {fit['slope']:+.3f}"
                          for floor, fit in zip(FLOOR_SWEEP, swept)))
    saturation = peroxide_saturation(named["BnOH two-axis (135-151)"])
    print(f"   and the rate's own order, on the {saturation['points']}-curve "
          f"ladder over {saturation['peroxide_low']:.2f}-"
          f"{saturation['peroxide_high']:.0f} mM:")
    print(f"   free power law      a = {saturation['order']:.3f} "
          f"({saturation['order_low']:.3f} to {saturation['order_high']:.3f})")
    print(f"   strict first order  rejected, F = "
          f"{saturation['first_order_f']:.1f}")
    print(f"   the scheme's own form fits "
          f"{'worse' if saturation['scheme_sse'] > saturation['power_sse'] else 'better'}"
          f" than a free power ({saturation['scheme_sse']:.2f} against "
          f"{saturation['power_sse']:.2f}), so the fractional order is not one")
    print("   binding equilibrium saturating.")
    print("   IF the positive order is real it is a TRAP, not the activation:")
    for label, block in peroxide_blocks:
        got = scope.orders("t_ind", frame=block, floor=INDUCTION_FLOOR)
        trap = trap_constant(got["order_h2o2"], got["stderr_h2o2"],
                             peroxide_geometric_mean(block))
        print(f"   {label:30s} K = {trap['constant']:.4f} +- "
              f"{trap['stderr']:.4f} /mM = {trap['molar']:.0f} /M   "
              f"dG = {trap['free_energy_kJ']:+.2f} kJ/mol")
    print(f"   {'the same K from the rates':30s} K = "
          f"{saturation['constant']:.4f} /mM "
          f"({saturation['constant_low']:.4f}-{saturation['constant_high']:.4f})"
          f"   {saturation['bound_fraction']:.0%} bound at "
          f"{saturation['working_mM']:.1f} mM")

    print("\n4c. WHICH WAY THE INDUCTION POINTS, AND WHAT MOVES IT")
    for name in ("4OMe catalysed", "temperature series",
                 "BnOH two-axis (135-151)", "4OMe enzyme-free"):
        block = named[name]
        block = block[block.live]
        counts = block.progress_kind.value_counts()
        print(f"   {name:26s} lag-first {int(block.lag_first.sum()):3d} of "
              f"{len(block):3d}   " + "  ".join(f"{k} {v}" for k, v
                                                in counts.items()))
    print("   [S] and [buf] move together inside a run, which is what makes")
    print("   a substrate order in the 4OMe block an order in the pair:")
    for name in ("4OMe catalysed", "temperature series",
                 "BnOH two-axis (135-151)"):
        got = composition_collinearity(named[name])
        print(f"   {name:26s} {got['runs']:2d} runs with a ladder, "
              f"median corr(log S, log buf) {got['median']:+.2f}, "
              f"{got['constant_buffer']} with [buf] constant")
    print("   P(lag first) per e-fold, within runs, signal-to-noise controlled:")
    arms = ladder_arms(named["BnOH two-axis (135-151)"])
    for label, block, axis in (("two-axis, substrate arm",
                                arms["substrate arm"], "s0"),
                               ("two-axis, peroxide arm",
                                arms["peroxide arm"], "h2o2"),
                               ("4OMe catalysed ([S] and [buf])",
                                named["4OMe catalysed"], "s0")):
        for control in (False, True):
            got = sign_drivers(block, axis=axis, control=control)
            tail = (f"   signal {got['signal']:+.3f} +- "
                    f"{got['signal_stderr']:.3f}" if control else "")
            print(f"   {label:32s} {'+ S/N' if control else '     '} "
                  f"{axis:5s} {got[axis]:+.3f} +- {got[axis + '_stderr']:.3f}"
                  f"   n={got['points']:3d}{tail}")
    print("   the one direct buffer lever, read on a common footing:")
    ladder = buffer_lever(table)
    print(ladder.to_string(index=False))
    step = buffer_join_step(ladder)
    print(f"   the rate falls {step['step']:.2f}x from {step['from_buf']:.0f} "
          f"to {step['to_buf']:.0f} mM, so the runs get separate levels")
    for response in ("t_ind", "depth"):
        got = buffer_order(ladder, response)
        print(f"   d log {response:6s}/d log[buf]  pooled "
              f"{got['slope']:+.3f} +- {got['stderr']:.3f}   "
              f"exp 34 {got['slope_34']:+.3f} +- {got['stderr_34']:.3f}   "
              f"exp 32 {got['slope_32']:+.3f} +- {got['stderr_32']:.3f}")
    print("   and across windows, none of which changes its sign:")
    print("   " + "  ".join(
        f"{width:.0f}s {buffer_order(buffer_lever(table, width=width))['slope']:+.3f}"
        for width in BUFFER_WINDOW_SWEEP))
    print("   4b's constraint again, on the buffer axis instead: the same +1,")
    print("   because the algebra is about a pre-equilibrium and not about H2O2.")
    print("   Each window with the signal control that decides whether to read it:")
    for width in BUFFER_WINDOW_SWEEP:
        rungs = buffer_lever(table, width=width)
        joint = joint_buffer_order(rungs)
        control = signal_control(rungs)
        verdict = ("passes" if abs(control["signal_slope"])
                   < 2 * control["signal_stderr"] else "FAILS")
        print(f"   {width:5.0f}s  order(v/t_ind) {joint['slope']:+.3f} +- "
              f"{joint['stderr']:.3f}  {joint['sigma']:.1f} sigma from +1   "
              f"S/N {control['signal_slope']:+.3f} +- "
              f"{control['signal_stderr']:.3f} {verdict}")
    print("   The windows that fail the control are the windows that overshoot,")
    print("   which is the direction the artefact predicts: more buffer, more")
    print("   signal, an earlier landmark, an inflated gap.")
    order = buffer_order(buffer_lever(table))
    fixed = substrate_order_corrected(named["4OMe catalysed"], order["slope"],
                                      order["stderr"])
    print(f"   what that does to section 3's route two: "
          f"{fixed['measured']:+.3f} +- {fixed['measured_stderr']:.3f} becomes "
          f"{fixed['corrected']:+.3f} +- {fixed['corrected_stderr']:.3f}, "
          f"against {fixed['threshold']:+.3f}")

    print("\n5. THE ACTIVATION PARAMETERS")
    gap = activation_contrast()
    for key in ("induction", "turnover"):
        row = gap[key]
        print(f"   {row['label'][:34]:34s} Ea {row['activation_kJ']:6.2f} "
              f"+- {row['activation_stderr']:5.2f}   "
              f"dH {row['enthalpy_kJ']:6.2f}   "
              f"dS {row['entropy_J']:+7.2f} +- {row['entropy_stderr']:5.2f}   "
              f"dG {row['gibbs_kJ']:6.2f}")
    print(f"   entropy gap {gap['entropy_gap_J']:+.1f} "
          f"+- {gap['entropy_gap_stderr']:.1f} J/mol/K")
    print(f"   enthalpy gap {gap['enthalpy_gap_kJ']:+.1f} "
          f"+- {gap['enthalpy_gap_stderr']:.1f} kJ/mol")
    print(f"   free energy gap {gap['gibbs_gap_kJ']:+.2f} "
          f"+- {gap['gibbs_gap_stderr']:.2f} kJ/mol "
          f"-> {gap['rate_ratio']:.0f}x faster at 298 K")

    print("\n6. THE OTHER SUBSTRATE")
    print("   exps 135-151 vary [H2O2] 30-fold inside every run, which the")
    print("   4OMe archive does not; the temperature series holds it at 82.5")
    print("   mM on all 24 curves and can carry no peroxide order at all.")
    scoped = named["BnOH two-axis (135-151)"]
    for parameter, floor in (("t_ind", INDUCTION_FLOOR),
                             ("depth", DEPTH_FLOOR), ("vmax", None)):
        got = scope.orders(parameter, frame=scoped, floor=floor)
        print(f"   order of {parameter:6s} in [S] {got['order_s0']:+.3f} "
              f"+- {got['stderr_s0']:.3f}   in [H2O2] "
              f"{got['order_h2o2']:+.3f} +- {got['stderr_h2o2']:.3f}   "
              f"n={got['n']}")
    return table


def main():
    report()


if __name__ == "__main__":
    main()
