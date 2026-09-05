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
import functools
from dataclasses import dataclass

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from curve_metrics import LAG_WINDOW, debubble, rolling_slope
from fit_dataset import source_floor
from summary_kinetics import fit_progress
from arrhenius import GAS_CONSTANT
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
        # THE SAME INTERVAL `t_ind` SPANS, and it was not until 2026-09-05:
        # `t_ind` is measured from the FIRST WINDOW CENTRE and `made` was
        # measured from the first READING, so it carried an extra half-window
        # of product -- 2000 s of it on a 40000 s run. That offset is
        # proportional to the rate, so it manufactured exactly the +1 slope
        # `product_at_landmark` tests for: a planted product THRESHOLD read
        # back as +1.03 instead of 0.00. Nothing published moved, because this
        # field was computed from the day the module was written and read by
        # nothing until the test that caught it.
        made=float(np.interp(centres[index], times, values)
                   - np.interp(centres[0], times, values)),
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


# Where the two hypotheses put the slope of log(product made by the landmark) on
# log(rate). This is `induction_drivers`' question asked in the units the
# PRODUCT hypothesis is stated in: a clock lets the product land wherever the
# rate puts it in a fixed time, so the product is proportional to the rate; a
# threshold fixes the product and lets the time land where it may.
PRODUCT_CLOCK_SLOPE = 1.0
PRODUCT_THRESHOLD_SLOPE = 0.0

# The geometries `product_recovery` plants on: (run length, tau at the middle
# of the ladder). Between them they span the archive's own range -- runs of
# 40000 down to 60000 s and time constants from a twentieth of the run to a
# third of it, which is the band the landmark can resolve at all.
PRODUCT_RECOVERY_PLANTINGS = ((40000.0, 1500.0), (60000.0, 6000.0),
                              (60000.0, 10000.0), (100000.0, 15000.0))

# The rate lever the plantings use, matched to the archive's own: the substrate
# ladder inside a 4OMe run moves the rate about four-fold (`substrate_lever`).
PRODUCT_RECOVERY_LEVER = 4.0
PRODUCT_RECOVERY_CURVES = 12


@functools.lru_cache(maxsize=4)
def product_recovery(plantings=PRODUCT_RECOVERY_PLANTINGS,
                     lever=PRODUCT_RECOVERY_LEVER,
                     curves=PRODUCT_RECOVERY_CURVES):
    """
    What `product_at_landmark` reads back when each hypothesis is TRUE.

    THE LANDMARK IS BIASED TOWARDS THE CLOCK AND THE BIAS IS NOT SMALL. The
    rolling slope's first centre sits half a window into the run, so part of the
    induction has already happened before the landmark's own clock starts --
    and how much depends on tau, which is exactly what a product threshold
    varies. Writing it out for a threshold planting with `v.tau = C` fixed:

        made = C(ln2 - 0.5 e^(-c0/tau))

    with `c0` the half-window. As tau falls the exponential dies and `made`
    rises by a factor of 3.6 across the band, which is a POSITIVE slope on the
    rate where the hypothesis says zero.

    So a planted threshold reads back near **+0.4**, not 0.0, and the sigma a
    result should be quoted against is the distance from that rather than from
    the nominal value. A planted clock reads back at +1.000 exactly, so the
    other end needs no correction.

    Returns the recovered slope for each hypothesis at each planting, and the
    worst (most clock-like) threshold reading, which is what
    `product_at_landmark` uses.
    """
    rows = []
    for run, middle in plantings:
        for name, tau_of in (
                ("clock", lambda rate, middle=middle: middle),
                ("threshold", lambda rate, middle=middle: middle * 2e-6
                 * np.sqrt(lever) / rate)):
            planted = []
            for index in range(curves):
                rate = 2e-6 * lever ** (index / (curves - 1.0))
                tau = tau_of(rate)
                times = np.linspace(0.0, run, 1000)
                values = rate * (times - tau * (1.0 - np.exp(-times / tau)))
                found = induction_point(_PlantedCurve(times, values))
                planted.append({"experiment": 1, "sample": index + 1,
                                "live": True, "t_ind": found.t_ind,
                                "made": found.made, "v_peak": found.peak_rate,
                                "epsilon": 1.0, "s0": 1.0})
            got = product_at_landmark(pd.DataFrame(planted),
                                      minimum_per_run=curves, calibrate=False)
            rows.append({"run_s": run, "tau_s": middle, "planted": name,
                         "recovered": got.get("slope", np.nan),
                         "spread_time": got.get("spread_time", np.nan),
                         "spread_product": got.get("spread_product", np.nan)})
    table = pd.DataFrame(rows)
    threshold = table[table.planted == "threshold"].recovered
    clock = table[table.planted == "clock"].recovered
    return {"table": table,
            "threshold_reads": float(threshold.max()),
            "threshold_range": (float(threshold.min()), float(threshold.max())),
            "clock_reads": float(clock.mean()),
            "clock_range": (float(clock.min()), float(clock.max()))}


class _PlantedCurve:
    """The three attributes `induction_point` reads, for a synthetic curve."""

    def __init__(self, times, absorbance, source="rre"):
        self.times = np.asarray(times, dtype=float)
        self.absorbance = np.asarray(absorbance, dtype=float)
        self.source = source


def product_at_landmark(table, rate="v_peak", minimum_per_run=3,
                        calibrate=True):
    """
    How much product had been made when the induction ended, and was it fixed?

    THE QUESTION IN ITS OWN UNITS. `induction_drivers` asks whether a faster
    cuvette's induction is SHORTER; this asks whether a faster cuvette's
    induction ends at the same PRODUCT, which is what "the induction waits for
    product" actually claims. `InductionPoint.made` -- the absorbance built up
    by the landmark -- has been computed on every curve since this module was
    written and read by nothing until 2026-09-05.

        d log(made) / d log(rate) = +1     a clock: fixed time, so the product
                                           is whatever the rate makes in it
                                         =  0     a threshold: fixed product,
                                           so the time is whatever it takes

    IT IS NOT INDEPENDENT OF `induction_drivers` and must not be quoted as
    though it were. Over the induction the product is roughly the rate times
    the time, so this slope is about `1 + (that one)` by construction: route
    one's -0.025 +/- 0.109 predicts +0.975 and this returns +0.979. What it
    adds is the units, and the second number below, which is not a regression
    at all.

    THE SPREAD COMPARISON, which owes nothing to any fit. Inside one run the
    schedule, pH, temperature, buffer and catalyst are fixed and only the
    substrate ladder moves, changing the rate about four-fold. So ask which of
    the two candidate constants is actually more nearly constant across those
    cuvettes: the pooled within-run standard deviation of log(t_ind) against
    that of log(made). A clock says the first is smaller and a threshold says
    the second is. On the catalysed 4OMe block they are 0.667 and 1.091.
    That comparison carries the same bias and `product_recovery` prices it: a
    planted threshold gives 0.44 against 0.19 rather than 0.44 against 0, so
    the two are not expected to separate cleanly -- and the archive's separate
    in the clock's direction.

    WHAT IT COSTS. A curve with no landmark has made no product by it, so this
    is the one statistic in the module that DROPS rows: 112 of the 147 live
    catalysed 4OMe curves carry a landmark. That is a selection towards curves
    that have an induction, which is the population the question is about, but
    it is a selection and the count is returned so it can be read.
    """
    live = table[table.live & (table.t_ind > 0) & (table.made > 0)
                 & (table[rate] > 0)]
    dropped = int(table.live.sum()) - len(live)
    if len(live) < 10:
        return {"curves": int(len(live)), "dropped": dropped}
    response = np.log(live.made.to_numpy(dtype=float))
    columns = [np.log(live[rate].to_numpy(dtype=float))]
    columns += [(live.experiment == run).to_numpy(float)
                for run in sorted(live.experiment.unique())]
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    resid = response - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(live) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    slope = float(beta[0])
    stderr = float(np.sqrt(max(covariance[0, 0], 0.0)))

    spreads = {"t_ind": [], "made": [], "conversion": [], rate: []}
    for _, group in live.groupby("experiment"):
        if len(group) < minimum_per_run:
            continue
        converted = (group.made / (group.epsilon * group.s0)).to_numpy(float)
        for name, values in (("t_ind", group.t_ind.to_numpy(float)),
                             ("made", group.made.to_numpy(float)),
                             ("conversion", converted),
                             (rate, group[rate].to_numpy(float))):
            if np.all(values > 0):
                spreads[name].append(float(np.std(np.log(values), ddof=1)))
    pooled = {name: float(np.sqrt(np.mean(np.square(values))))
              if values else np.nan for name, values in spreads.items()}
    # AGAINST WHAT A PLANTING ACTUALLY READS BACK, not against the nominal
    # value. A true product threshold does not read 0.0 through this landmark,
    # it reads about +0.4, because the rolling slope's first centre is half a
    # window into the run -- `product_recovery` derives it and plants it.
    # Quoting the nominal 0.0 would put the archive's answer 8.9 sigma from a
    # threshold when the honest distance is 5.3.
    reads = (product_recovery() if calibrate
             else {"threshold_reads": PRODUCT_THRESHOLD_SLOPE,
                   "clock_reads": PRODUCT_CLOCK_SLOPE})
    return {"curves": int(len(live)), "dropped": dropped,
            "experiments": int(live.experiment.nunique()),
            "runs_in_spread": len(spreads["t_ind"]),
            "slope": slope, "stderr": stderr,
            "clock_reads": reads["clock_reads"],
            "threshold_reads": reads["threshold_reads"],
            "clock_sigma": float(abs(slope - reads["clock_reads"]) / stderr)
            if stderr else np.nan,
            "threshold_sigma": float(abs(slope - reads["threshold_reads"])
                                     / stderr) if stderr else np.nan,
            "spread_time": pooled["t_ind"], "spread_product": pooled["made"],
            "spread_conversion": pooled["conversion"],
            "spread_rate": pooled[rate]}


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


def joint_order(table, axis, rate="v_peak", timescale="t_ind",
                covariates=(), floor=INDUCTION_FLOOR, minimum=8,
                live_only=True, gate=None, required=PERHYDRATE_ORDER_GAP):
    """
    `d ln(rate)/d ln(axis) - d ln(timescale)/d ln(axis)`, as ONE regression.

    THE ALGEBRA DOES NOT CARE WHICH SPECIES, WHICH RATE, OR WHICH CLOCK. A
    catalyst drawn into its active form by a species X held in excess gives

        E + X <-> E*    1/tau = k_f[X] + k_r,   [E*]/E0 = K[X]/(1 + K[X])

    so `d ln v/d ln[X] = 1/(1 + K[X])` and `d ln tau/d ln[X] = -K[X]/(1 +
    K[X])`, and their DIFFERENCE is 1 for every K and every [X]. That is a
    parameter-free prediction, and it holds for any timescale of the
    activation step and any axis the archive moves.

    ONE regression, not two differenced by hand. The two orders it replaces are
    fitted to the same curves on the same design, so their errors are
    correlated and a hand-differenced error is not theirs. `covariates` are the
    other axes that move, in logs; one free level per experiment absorbs pH,
    buffer, catalyst, cell and day.

    `floor` clips the timescale from below, which is right for a LANDMARK --
    `t_ind` is resolution-limited and a curve with no induction reads zero --
    and wrong for a FITTED time constant, which is either resolved or not.
    Pass `floor=None` with `gate="tau_slow_resolved"` for the latter: the gate
    is the honest filter, and clipping would put unresolved constants on the
    floor and call them fast.

    THE AXIS IS ALSO THE CONTROL. The +1 belongs to the species that activates
    the catalyst, so an axis whose species does NOT should miss it -- and on
    this archive's substrate axis it misses by 7 sigma and more. A test that
    only ever confirms is not a test; run the control axis beside the claim.

    `joint_peroxide_order` and `joint_buffer_order` are this function with the
    species filled in. Returns the coefficient, its error, and how many
    standard errors it sits from `required`.
    """
    live = table[table.live] if live_only else table
    if gate is not None:
        live = live[live[gate].astype(bool)]
    clock = live[timescale].to_numpy(dtype=float)
    keep = (np.isfinite(clock)
            & (live[rate].to_numpy(dtype=float) > 0)
            & (live[axis].to_numpy(dtype=float) > 0))
    if floor is None:
        keep = keep & (clock > 0)
    live = live[keep]
    if len(live) < minimum or live[axis].nunique() < 2:
        return {"points": int(len(live))}
    clock = live[timescale].to_numpy(dtype=float)
    if floor is not None:
        clock = np.maximum(clock, floor)
    response = np.log(live[rate].to_numpy(dtype=float)) - np.log(clock)
    columns = [np.log(live[axis].to_numpy(dtype=float))]
    for name in covariates:
        if live[name].nunique() > 1:
            columns.append(np.log(live[name].to_numpy(dtype=float)))
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
    return {"points": int(len(live)), "experiments": len(runs),
            "runs": len(runs), "slope": slope, "stderr": stderr,
            "required": required,
            "sigma": float(abs(required - slope) / stderr)
            if stderr else np.nan}


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
    return joint_order(table, axis="h2o2", rate="v_peak", timescale="t_ind",
                       covariates=("s0",), floor=floor, minimum=8)


# The clocks `joint_order` can be asked for, and the flag that says a fitted
# one was actually pinned. `t_ind` is the landmark and carries a floor; the
# other two come from the progress fit and carry a gate instead.
#
# THE DEFAULTS ARE THE REBUILT ONES, and on a peroxide axis that is not a
# detail. The O2 is made from peroxide, so leaving it in inflates the rate's
# peroxide order AND shortens the apparent clock, and both push
# `d ln v - d ln tau` towards the +1 this function tests. Asked of the readings
# the two-axis block's `tau_slow` row sits 0.3 sigma from +1; asked of the
# rebuilt curves it sits 1.4. `JOINT_CLOCKS_RAW` is the uncorrected set, kept
# so the difference can be shown rather than asserted.
# Each entry is (clock, resolution gate, is it windowed). The third field is
# carried rather than derived from the second: `lag_half_s` has NO gate, being
# a landmark on the fitted rate rather than a fitted constant, and NO window,
# and reading "windowed" off "gate is None" would have labelled it with the
# one property it was added to remove.
JOINT_CLOCKS_RAW = (("t_ind", None, True), ("tau", "tau_resolved", False),
                    ("tau_slow", "tau_slow_resolved", False))
JOINT_CLOCKS = (("t_ind", None, True),
                ("lag_half_s", None, False),
                ("tau_corrected", "tau_resolved_corrected", False),
                ("tau_slow_corrected", "tau_slow_resolved_corrected", False))


def joint_clocks(table, axis="h2o2", control="s0", rate="vmax_corrected",
                 clocks=JOINT_CLOCKS, floor=INDUCTION_FLOOR):
    """
    The +1 constraint on one axis, through every clock, beside its control.

    THE BLOCK THAT CANNOT USE A LANDMARK CAN STILL BE ASKED. `t_ind` is a
    rolling window a tenth of the run wide, so it is not comparable between
    runs of different length, and on the two-axis block `signal_control` fails
    -- the landmark there is partly measuring the spectrophotometer. `tau` and
    `tau_slow` come from the progress fit, carry no window, and are subject to
    the same identity. Where the two routes disagree the disagreement is the
    result, and it has to be visible in one table rather than argued.

    One row per clock per axis: the axis under test, the control axis, the
    curves each survives, and how far each sits from +1. READ THE CONTROL. The
    +1 belongs to the species that draws the catalyst into its active form, so
    an axis whose species does not -- the alcohol -- must MISS it, and a run of
    this table where both axes meet +1 is a regression that has stopped
    discriminating rather than a mechanism.

    `rate` is the fitted clocks' partner and `v_peak` the landmark's, which is
    the pairing each was defined with; the landmark's floor applies to it alone.

    THE LANDMARK ROW IS NOT GAS-CORRECTED and cannot be from this table --
    `t_ind` and `v_peak` are read off the curve as it stands. That row is the
    one the two-axis block rejects on other grounds anyway (`signal_control`),
    so it is carried as the comparator it is rather than repaired; do not read
    the gap between it and the fitted rows as though the two were like for like.
    """
    rows = []
    for clock, gate, windowed in clocks:
        if clock not in table.columns:
            continue
        if gate is not None and gate not in table.columns:
            continue
        for name, tested, other in (("axis", axis, control),
                                    ("control", control, axis)):
            got = joint_order(
                table, axis=tested,
                rate="v_peak" if clock == "t_ind" else rate,
                timescale=clock, covariates=(other,),
                floor=floor if gate is None else None, gate=gate)
            rows.append({"clock": clock, "role": name, "axis": tested,
                         "windowed": windowed,
                         "curves": got.get("points", 0),
                         "order": got.get("slope", np.nan),
                         "stderr": got.get("stderr", np.nan),
                         "sigma": got.get("sigma", np.nan)})
    return pd.DataFrame(rows).set_index(["clock", "role"])


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
    return joint_order(ladder, axis="buf", rate="v_peak", timescale="t_ind",
                       floor=floor, minimum=4, live_only=False)


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
# ---------------------------------------------------------------------------
# THE LAG WITHOUT A WINDOW, AND THE SEVEN AXES.
#
# Everything above is measured through `induction_point`, whose rolling window
# is a tenth of the run. That window is why this folder had, until 2026-09-05,
# measured the induction against exactly two of the seven variables the
# experiment moves: `[S]` and `[H2O2]`, both of which step INSIDE a run. The
# other five -- temperature, pH, `[buf]`, the buffer salt and the substrate --
# are one value per run in this archive, always, so every one of them is a
# BETWEEN-run comparison, and a window in run-fractions makes a between-run
# comparison a comparison of windows.
#
# `scope.frame` now carries `lag_depth` and `lag_half_s`, read off the fitted
# rate (`summary_kinetics.ProgressFit.lag_profile`) rather than a rolling one.
# They carry no window. What they still carry is the SCHEDULE -- the fit's tau
# grid runs from span/300 to 2.span, so a 1260 s run cannot report a 4000 s
# clock -- and the SIGNAL, because a curve with little signal gives the fit
# little to choose between the shapes. This section answers the first by
# refitting a block on a window every run in it shares, and the second by
# carrying the control beside every answer rather than in a footnote.
#
# THE THREE CONFOUNDS, and why the archive needs all four pH ladders. Across
# `scope.PH_LADDER_PHOSPHATE` pH correlates with log(run length) at -0.43 and
# with log(signal/noise) at +0.92; across `scope.PH_LADDER_BORIC` the same two
# are +0.71 and +0.67. The signs differ, so a coefficient that survives in both
# is not either confound. Read `lag_ladder`, never a raw regression on pH.

# Nodes of a rolling window are not involved anywhere below, so the only floor
# a lag clock needs is the one that makes "no lag" and "a lag shorter than one
# reading" the same number. It is INDUCTION_FLOOR, deliberately: two floors for
# one idea is how the lag statistic came to have two definitions in the first
# place.

# A window has to leave enough readings to fit four parameters and profile a
# time constant. Eight is `BURST_MINIMUM_POINTS` plus two.
LAG_WINDOW_MINIMUM_POINTS = 8

# The fractions of the shortest run in a block that `lag_window_sweep` walks.
# 1.0 is the longest window every run can supply; below it the block keeps all
# its runs and loses the slow ones's clocks, which is the trade the sweep
# exists to show.
LAG_WINDOW_SWEEP = (0.5, 0.75, 1.0)

# What the two bounded schemes allow, in the units each axis is read in.
# For a species X in pre-equilibrium OFF the activation path -- the trap of
# section 4a -- 1/tau = k_act/(1 + K[X]), so d ln tau/d ln[X] = K[X]/(1 + K[X])
# lies in (0, 1) and is the SATURATION FRACTION of that trap. For a species
# whose BOUND form activates, 1/tau = k_act.K[X]/(1 + K[X]) and the same
# quantity lies in (-1, 0). Neither bound involves a rate constant, so a
# measurement outside them falsifies the scheme rather than fitting it.
TRAP_BOUND = (0.0, 1.0)
ACTIVATING_BOUND = (-1.0, 0.0)

# One pH unit is one decade of every species whose formation consumes a proton,
# so a pH coefficient is converted to an order in that species by dividing by
# ln 10 before it is compared with the bounds above.
LN10 = float(np.log(10.0))


def common_window(experiments):
    """The longest window every run in `experiments` can supply, in seconds."""
    return min(float(np.asarray(curve.times, dtype=float)[-1])
               for curve in scope.curves(tuple(experiments)))


def lag_window_frame(experiments, window=None):
    """
    Refit every curve of `experiments` on a window they ALL share.

    One row per curve, with `lag_depth`, `lag_half_s`, `lag_peak`, the sign,
    and the design. Returns a COPY of a memoised frame, for the reason
    `scope.frame` does.

    WHY REFITTING RATHER THAN CONTROLLING. Run length is not a nuisance
    variable that a regression can absorb here, it is a bound on what the fit
    can report: the tau grid runs from span/300 to 2.span, so a run cannot
    measure a relaxation much longer than itself. In `scope.PH_LADDER_BORIC`
    the runs span 1260 to 17940 s and pH correlates with log(run length) at
    +0.71, so the raw pH coefficient on the clock (+1.376 +/- 0.318 per pH
    unit) is partly the schedule. Putting log(run length) in the regression
    hands the answer to a 9-point collinearity; truncating every run to 1260 s
    and refitting removes it by construction. What truncation costs is the
    clocks longer than the window, and `lag_window_sweep` is how that cost is
    shown rather than assumed.

    The same tool answers the same question for the temperature series, where
    it matters most: cold runs are the long ones (1/T against log(run length),
    r = +0.66), and the activation energy of the induction moves from
    83.7 +/- 8.9 kJ/mol on whole runs to 55.3 +/- 24.3 if log(run length) is
    put in the regression instead. At a common window it is 73-84 across every
    window that keeps all six temperatures. `lag_arrhenius`.

    The readings are gas-corrected first (`curve_metrics.debubble`), like every
    other clock in this package since 2026-09-04.
    """
    frame = _lag_window_frame(tuple(experiments),
                              None if window is None else float(window))
    return frame.copy()


@functools.lru_cache(maxsize=32)
def _lag_window_frame(experiments, window):
    """`lag_window_frame`'s work, memoised. Call `lag_window_frame`."""
    if window is None:
        window = common_window(experiments)
    rows = []
    for curve in scope.curves(experiments):
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        corrected, _ = debubble(times, values, curve.noise)
        keep = times <= window
        if int(keep.sum()) < LAG_WINDOW_MINIMUM_POINTS:
            continue
        times, corrected = times[keep], corrected[keep]
        span = float(times[-1] - times[0])
        fit = fit_progress(times, corrected)
        depth, half, peak, start = fit.lag_profile(span)
        net = float(corrected[-1] - corrected[0])
        rows.append({
            "experiment": curve.experiment, "sample": curve.sample,
            "substrate": curve.substrate, "buffer": curve.buffer,
            "temperature": curve.temperature, "kelvin": curve.temperature + 273.15,
            "pH": curve.pH, "buf": curve.buf,
            "s0": curve.conditions.s0, "h2o2": curve.conditions.h2o2,
            "e0": curve.conditions.e0, "hoo": curve.conditions.hoo,
            "window_s": span, "duration_s": float(np.asarray(curve.times)[-1]),
            "points": int(len(times)),
            "net": net, "noise": curve.noise,
            "live": bool(net > scope.LIVE_SIGNAL_NOISE_MULTIPLE * curve.noise),
            "lag_depth": depth, "lag_half_s": half,
            "lag_peak": peak, "lag_start": start,
            "phases": int(fit.phases), "progress_kind": fit.kind,
            "B_fast": fit.amplitudes[0],
            "lag_first": fit.kind in LAG_FIRST_KINDS,
        })
    return pd.DataFrame(rows)


# A regressor is identified only where the offsets cannot absorb it. This is
# the fraction of an axis's variance that survives projecting the offsets out;
# below it, the coefficient is whatever the pseudo-inverse's minimum-norm
# solution happens to split off, and its standard error is a fit to noise.
# It cost an afternoon on 2026-09-05: `lag_ph_ladders` first ran with one
# offset per EXPERIMENT on ladders where pH is one value per experiment, and
# reported +0.549 +/- 0.014 -- a forty-sigma pH effect that was the collinearity
# and nothing else. `scope._moves` guards the same mistake on the within-run
# orders; this guards it on the between-run ladders.
IDENTIFIED_SHARE = 1e-6


def _identified(values, offsets_matrix):
    """The share of `values`'s variance the offsets cannot absorb."""
    spread = float(((values - values.mean()) ** 2).sum())
    if spread <= 0:
        return 0.0
    beta, *_ = np.linalg.lstsq(offsets_matrix, values, rcond=None)
    resid = values - offsets_matrix @ beta
    return float(resid @ resid) / spread


def _lag_fit(table, response, terms, offsets=None, floor=INDUCTION_FLOOR):
    """
    Least squares of a lag response on `terms`, with one offset per group.

    `response` is "lag_half_s" (logged, floored), "lag_depth" (as it is, since
    it is already a fraction) or "lag_first" (0/1, a linear probability model
    for the reason `sign_drivers` gives). `offsets` is a column name -- usually
    "experiment", and "sample" on the two-axis pH ladders, where a cuvette is
    matched across runs and a run offset would absorb the pH axis itself.
    """
    if response == "lag_half_s":
        y = np.log(np.maximum(table.lag_half_s.to_numpy(dtype=float), floor))
    elif response == "lag_first":
        y = table.lag_first.to_numpy(dtype=float)
    else:
        y = table[response].to_numpy(dtype=float)
    columns = [np.asarray(term, dtype=float) for term in terms]
    if offsets is None:
        columns.append(np.ones(len(table)))
    else:
        for level in sorted(table[offsets].unique()):
            columns.append((table[offsets] == level).to_numpy(float))
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(y) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    return (float(beta[0]), float(np.sqrt(max(covariance[0, 0], 0.0))),
            int(len(y)))


def _partial_correlation(first, second, table, offsets):
    """Correlation of two columns with the offsets projected out of both."""
    if offsets is None:
        levels = np.ones((len(table), 1))
    else:
        levels = np.column_stack(
            [(table[offsets] == level).to_numpy(float)
             for level in sorted(table[offsets].unique())])
    out = []
    for column in (first, second):
        beta, *_ = np.linalg.lstsq(levels, column, rcond=None)
        out.append(column - levels @ beta)
    if out[0].std() == 0 or out[1].std() == 0:
        return np.nan
    return float(np.corrcoef(out[0], out[1])[0, 1])


LAG_RESPONSES = ("lag_half_s", "lag_depth", "lag_first")


LAG_ORDER_FLOORS = {"lag_half_s": INDUCTION_FLOOR, "lag_depth": DEPTH_FLOOR}


def lag_channel_table(frame=None):
    """
    The window-free lag by substrate and channel: four cells, and they differ.

    `channel_summary` above does this for the LANDMARK on the 4OMe block, where
    the contrast is the folder's first claim -- 0 of 49 enzyme-free curves have
    an induction. Run over the whole archive the same contrast holds on 4OMe
    and DOES NOT HOLD ON BnOH: 14 of 26 enzyme-free BnOH curves begin below
    their eventual rate, at a median depth of 0.138.

    That is not a counter-example to section 2, it is a second phenomenon. The
    BnOH lags sit in exps 3 and 6 -- pH 6.71 phosphate, 11025 and 17934 s, 8 of
    10 curves -- and in exp 65, whose four cuvettes share the synchronised
    break `scope.synchronised_break` measures. Steps 1-3 of `MECHANISM.md` are
    autocatalytic in the product and need no catalyst, and `product_fate`
    already finds them switched ON for benzaldehyde and OFF for the
    4-methoxy aldehyde, whose electron-rich ring makes it a worse hydride
    donor. An accelerating enzyme-free BnOH curve is what that predicts.

    So the archive separates the two: on 4OMe the lag needs the catalyst and is
    its activation; on BnOH a lag can also be the catalyst-free loop finding its
    product. Do not pool them.
    """
    if frame is None:
        frame = scope.frame(scope.archive())
    live = frame[frame.live]
    rows = []
    for (substrate, catalysed), group in live.groupby(
            ["substrate", "differential"]):
        rows.append({
            "substrate": substrate,
            "channel": "catalysed" if catalysed else "enzyme-free",
            "curves": int(len(group)),
            "with_a_lag": int((group.lag_depth > 0).sum()),
            "median_depth": float(group.lag_depth.median()),
            "median_clock_s": float(group.lag_half_s.median()),
            "experiments": int(group.experiment.nunique())})
    return pd.DataFrame(rows).set_index(["substrate", "channel"])


def lag_orders(table, terms=scope.ORDER_TERMS, live_only=True,
               floors=None):
    """
    Within-run orders of the window-free lag, with and without the signal.

    Delegates to `scope.orders`, which is the package's one within-experiment
    log-log fit and carries the `_moves` guard; this adds only the pairing that
    makes the answer readable -- the same fit run again with the curve's own
    log(net/noise) held as a covariate, and the partial correlation that says
    whether holding it was ever going to matter.

    FIT EVERY AXIS THE BLOCK MOVES. `terms` defaults to ("s0", "h2o2") and has
    to be ("s0", "h2o2", "buf") on the 4OMe archive, where substrate volume
    displaced buffer volume. One axis at a time is a different regression on an
    L: the two-axis block's cuvettes carry log[S] against log[H2O2] at about
    -0.5, and a substrate-only fit of the induction clock there reads
    -0.453 +/- 0.107 where the joint fit reads -0.225 +/- 0.115.

    THE PAIR IS THE RESULT. A fitted lag tracks signal-to-noise on every block
    that moves peroxide, because in those blocks the peroxide IS the signal:
    +0.92 +/- 0.26 on the two-axis block, +0.98 +/- 0.29 on exps 127-131. Where
    the held and unheld coefficients agree the axis is separated from the
    signal; where they do not, the design cannot tell them apart, and that is
    the finding rather than a caveat on one.

    The sign of the early curve is NOT here. It is a 0/1 and takes no
    logarithm, so it has its own estimator in `sign_drivers`, which section 6
    already reports.
    """
    live = table[table.live] if live_only else table
    live = live[(live.net > 0) & (live.noise > 0)].copy()
    live["signal"] = live.net / live.noise
    out = {"points": int(len(live)),
           "experiments": int(live.experiment.nunique()), "terms": tuple(terms)}
    for term in terms:
        out[f"signal_collinearity_{term}"] = _partial_correlation(
            np.log(live[term].to_numpy(dtype=float)),
            np.log(live.signal.to_numpy(dtype=float)), live, "experiment")
    for response, floor in (floors or LAG_ORDER_FLOORS).items():
        bare = scope.orders(response, frame=live, floor=floor, terms=terms)
        held = scope.orders(response, frame=live, floor=floor, terms=terms,
                            covariates=("signal",))
        # How many rows sat ON the floor. `lag_half_s` is 0 on a curve with no
        # lag and `lag_depth` is 0 with it, and those zeros are measurements --
        # but a response that is mostly floor is a response the log-log fit
        # cannot see, and the count is the reader's warning.
        values = live[response].to_numpy(dtype=float)
        row = {"n": bare["n"], "r2": bare["r2"],
               "floored": int((values < floor).sum()),
               "signal_order": held.get("held_signal", np.nan),
               "signal_stderr": held.get("held_stderr_signal", np.nan)}
        for term in terms:
            row[term] = {"order": bare[f"order_{term}"],
                         "stderr": bare[f"stderr_{term}"],
                         "controlled": held[f"order_{term}"],
                         "controlled_stderr": held[f"stderr_{term}"]}
        out[response] = row
    return out


def lag_order_floor_sweep(table, axis, terms=scope.ORDER_TERMS,
                          response="lag_half_s", floors=FLOOR_SWEEP):
    """
    One axis's order at every floor, because 44% of the block sits on it.

    A curve that begins at its fastest has `lag_half_s = 0` and `lag_depth = 0`
    honestly, and those zeros are 48 of the two-axis block's 110 live curves.
    They cannot be dropped -- dropping them is the censoring the whole
    statistic was built to avoid -- so they are floored, and the floor is a
    choice. This walks it, the way `report` walks INDUCTION_FLOOR for the
    landmark. A coefficient that changes sign across the sweep is the floor's.
    """
    rows = []
    for floor in floors:
        got = lag_orders(table, terms=terms, floors={response: floor})
        order = got[response][axis]
        rows.append({"floor": floor, "order": order["order"],
                     "stderr": order["stderr"],
                     "controlled": order["controlled"],
                     "controlled_stderr": order["controlled_stderr"],
                     "floored": got[response]["floored"]})
    return pd.DataFrame(rows).set_index("floor")


def lag_sign(table):
    """
    `lag_first` on a frame from either source, derived once.

    `lag_window_frame` builds the column itself; `scope.frame` carries the
    shape as `progress_kind_corrected`, from the same fit on the same rebuilt
    readings. Both mean "the fast phase is a lag", and LAG_FIRST_KINDS is the
    single definition of which shapes those are.
    """
    if "lag_first" in table.columns:
        return table
    table = table.copy()
    table["lag_first"] = table.progress_kind_corrected.isin(LAG_FIRST_KINDS)
    return table


def lag_signal_control(table, offsets="experiment"):
    """
    Does the WINDOW-FREE lag track the curve's own signal-to-noise?

    `signal_control` above asks this of the landmark; this asks it of the
    fitted lag, and the two do not have the same answer. The catalysed 4OMe
    block passes on both (+0.003 +/- 0.149 landmark, +0.171 +/- 0.138 clock).
    The two-axis block fails on both (+0.619 +/- 0.228, +0.832 +/- 0.252) and
    exps 127-131 fail harder on the fit (+0.989 +/- 0.285) -- so moving to a
    window-free statistic does NOT rescue the peroxide axis, and the folder's
    standing refusal to read a peroxide order off those blocks stands.

    It is not a gate. A block that fails it can still carry an axis that is
    not collinear with its own signal, which is exactly what the two-axis
    block's substrate arm is: the block's rate order in `[S]` is +0.01 +/- 0.04,
    so more substrate buys no more signal, and the substrate coefficient on the
    clock is unmoved by the control (-0.431 -> -0.444). Read `lag_orders`,
    which reports both.
    """
    table = lag_sign(table)
    live = table[table.live & np.isfinite(table.lag_half_s)
                 & (table.net > 0) & (table.noise > 0)]
    if len(live) < LAG_WINDOW_MINIMUM_POINTS:
        return {"points": int(len(live))}
    signal = np.log((live.net / live.noise).to_numpy(dtype=float))
    out = {"points": int(len(live))}
    for response in LAG_RESPONSES:
        slope, stderr, _ = _lag_fit(live, response, [signal], offsets)
        out[response] = {"slope": slope, "stderr": stderr}
    return out


def lag_ladder(experiments, axis="pH", window=None, offsets=None,
               logged=False):
    """
    A BETWEEN-run ladder, read at a window every run in it shares.

    `axis="pH"` is read per pH UNIT and not logged, because pH is already a
    logarithm; every concentration axis is logged. `offsets="sample"` on the
    two-axis ladders, where the seven compositions repeat across runs and a
    run offset would absorb the ladder.

    Returns the coefficient on each lag response, alone and with the signal
    control, plus the window used, the runs kept and the two collinearities
    that decide whether the answer means anything: `axis` against log(run
    length) and against log(signal/noise), over the runs.
    """
    experiments = tuple(experiments)
    if window is None:
        window = common_window(experiments)
    table = lag_window_frame(experiments, window)
    live = table[table.live & np.isfinite(table.lag_half_s)
                 & (table.net > 0) & (table.noise > 0)]
    if live.experiment.nunique() < 3:
        return {"axis": axis, "window_s": float(window),
                "runs": int(live.experiment.nunique()), "points": int(len(live))}
    values = (np.log(live[axis].to_numpy(dtype=float)) if logged
              else live[axis].to_numpy(dtype=float))
    signal = np.log((live.net / live.noise).to_numpy(dtype=float))
    # THE OFFSETS MAY NOT ABSORB THE LADDER. On a between-run ladder the axis
    # is one value per run, so `offsets="experiment"` is collinear with it by
    # construction and returns a minimum-norm split with a meaningless error.
    # Pass None (a single level) or "sample" -- the two-axis ladders repeat
    # seven compositions across runs, so a cuvette offset is matched and a run
    # offset is the axis.
    if offsets is None:
        levels = np.ones((len(live), 1))
    else:
        levels = np.column_stack(
            [(live[offsets] == level).to_numpy(float)
             for level in sorted(live[offsets].unique())])
    share = _identified(values, levels)
    if share < IDENTIFIED_SHARE:
        raise ValueError(
            f"{axis} is absorbed by one offset per {offsets}: "
            f"{share:.2e} of its variance survives. A between-run ladder needs "
            f"offsets=None, or offsets='sample' where cuvettes are matched.")
    # `level` and not `axis`: `axis` is a reserved keyword of DataFrame.agg.
    runs = live.groupby("experiment").agg(
        level=(axis, "first"), duration=("duration_s", "first"),
        signal=("net", "median"), noise=("noise", "median"))
    run_axis = (np.log(runs.level.to_numpy(float)) if logged
                else runs.level.to_numpy(float))
    out = {"axis": axis, "window_s": float(window),
           "runs": int(live.experiment.nunique()), "points": int(len(live)),
           "schedule_collinearity": float(np.corrcoef(
               run_axis, np.log(runs.duration.to_numpy(float)))[0, 1]),
           "signal_collinearity": float(np.corrcoef(
               run_axis, np.log((runs.signal / runs.noise).to_numpy(float)))[0, 1])}
    for response in LAG_RESPONSES:
        slope, stderr, _ = _lag_fit(live, response, [values], offsets)
        held, held_stderr, _ = _lag_fit(live, response, [values, signal], offsets)
        out[response] = {"slope": slope, "stderr": stderr,
                         "controlled": held, "controlled_stderr": held_stderr}
    return out


def lag_window_sweep(experiments, axis="pH", offsets=None,
                     fractions=LAG_WINDOW_SWEEP, logged=False):
    """`lag_ladder` at fractions of the block's common window. The systematic."""
    full = common_window(experiments)
    return [lag_ladder(experiments, axis=axis, window=full * fraction,
                       offsets=offsets, logged=logged)
            for fraction in fractions]


def lag_ph_ladders(window_fraction=1.0):
    """Every pH ladder in the archive, read the same way. `scope.PH_LADDERS`."""
    rows = []
    for name, experiments in scope.PH_LADDERS.items():
        # One offset per CUVETTE on the two-axis ladders, where the seven
        # compositions repeat across runs; a single level on the other two,
        # whose runs share nothing but the design. Never one per experiment:
        # pH is one value per run everywhere in this archive.
        offsets = "sample" if experiments in (
            scope.PH_LADDER_TWO_AXIS_LOW, scope.PH_LADDER_TWO_AXIS_HIGH) else None
        window = common_window(experiments) * window_fraction
        result = lag_ladder(experiments, axis="pH", window=window,
                            offsets=offsets)
        result["ladder"] = name
        result["offsets"] = offsets
        rows.append(result)
    return rows


def pooled_ladder(results, response="lag_half_s", controlled=True):
    """
    Combine several ladders' coefficients on one response, inverse-variance.

    Four pH ladders is not four times one pH ladder: they sit in three buffers
    and on two substrates, and their confounds have OPPOSITE signs -- pH
    against log(run length) runs -0.25, +0.71, -0.53, -0.79 across them and
    against log(signal/noise) +0.87, -0.65, +0.77, +0.79. So a pooled
    coefficient is worth more than its error suggests if the four agree, and
    `chi2` on `dof` is how that is checked rather than assumed.

    It is worth LESS than its error suggests in one respect, which
    `lag_window_sweep` prices and this cannot: the window is a choice, and on
    three of the four ladders the coefficient moves further across windows than
    its own error. Quote the pooled value with the sweep beside it.
    """
    key = "controlled" if controlled else "slope"
    error = "controlled_stderr" if controlled else "stderr"
    rows = [r[response] for r in results
            if response in r and np.isfinite(r[response][error])
            and r[response][error] > 0]
    if len(rows) < 2:
        return {"points": len(rows)}
    values = np.array([r[key] for r in rows], dtype=float)
    weights = 1.0 / np.array([r[error] for r in rows], dtype=float) ** 2
    mean = float((weights * values).sum() / weights.sum())
    stderr = float(1.0 / np.sqrt(weights.sum()))
    chi2 = float((weights * (values - mean) ** 2).sum())
    return {"ladders": len(rows), "pooled": mean, "stderr": stderr,
            "chi2": chi2, "dof": len(rows) - 1,
            "values": [float(v) for v in values],
            "errors": [float(e) for e in 1.0 / np.sqrt(weights)]}


def lag_arrhenius(experiments=None, window=None, response="lag_half_s"):
    """
    The induction's barrier, from the window-free clock at a common window.

    `arrhenius.activation_parameters("inverse_tau")` measures the same barrier
    from the one-phase fit's `tau` and gets 95.0 +/- 15.7 kJ/mol on 16 curves
    at FOUR temperatures -- 15 to 30 C, because above 32 C a decelerating curve
    has no lag for the one-phase form to find and tau lands on the top of its
    grid (`arrhenius.BURST_TRUSTWORTHY_BELOW_C`). The two-phase form does find
    it, so `lag_half_s` reaches all six temperatures.

    IT ALSO REACHES THE CONFOUND. Cold runs are long ones here, 1/T against
    log(run length) at +0.66, and the clock cannot exceed the run: on whole
    runs the barrier is 83.7 +/- 8.9 and putting log(run length) in the
    regression drops it to 55.3 +/- 24.3 on a six-point collinearity. At a
    window all six runs share it is 73.2 +/- 11.9, and `lag_window_sweep` on
    this block shows it stable from there to the longest run.

    Returns the Arrhenius slope as an activation energy in kJ/mol.
    """
    if experiments is None:
        experiments = scope.TEMPERATURE_SERIES
    table = lag_window_frame(experiments, window)
    live = table[table.live & (table.lag_half_s > 0)]
    if live.temperature.nunique() < 3:
        return {"temperatures": int(live.temperature.nunique()),
                "points": int(len(live))}
    inverse = 1.0 / live.kelvin.to_numpy(dtype=float)
    y = np.log(1.0 / live.lag_half_s.to_numpy(dtype=float))
    design = np.column_stack([inverse, np.ones(len(y))])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    covariance = (float(resid @ resid) / max(1, len(y) - 2)
                  * np.linalg.pinv(design.T @ design))
    return {"window_s": float(table.window_s.max()),
            "points": int(len(y)),
            "temperatures": int(live.temperature.nunique()),
            "activation_kj": float(-beta[0] * GAS_CONSTANT / 1000.0),
            "stderr_kj": float(np.sqrt(max(covariance[0, 0], 0.0))
                               * GAS_CONSTANT / 1000.0),
            "clock_by_temperature": {
                float(t): float(g.lag_half_s.median())
                for t, g in live.groupby("temperature")}}


# Multiples of the block's common window. 1.0 is the longest window every run
# can supply and is the one to quote; below it the cold curves' clocks are
# longer than the window and the barrier collapses, which is the censoring this
# whole section is about, shown rather than avoided. Above it the runs that are
# shorter fall back to their own length, so 4.0 -- past the 3.55x that separates
# the temperature series' longest run from its shortest -- IS whole runs.
ARRHENIUS_WINDOW_SWEEP = (0.3, 0.5, 0.75, 1.0, 2.0, 4.0)


def lag_arrhenius_sweep(experiments=None, fractions=ARRHENIUS_WINDOW_SWEEP):
    """The barrier at multiples of the common window; the last is whole runs."""
    if experiments is None:
        experiments = scope.TEMPERATURE_SERIES
    full = common_window(experiments)
    out = []
    for fraction in fractions:
        result = lag_arrhenius(experiments, window=full * fraction)
        result["fraction"] = fraction
        out.append(result)
    return out


def replicate_floor(experiments=None, window=None):
    """
    How far two runs that differ in NOTHING move the lag parameters.

    `scope.REPLICATE_RUNS` is four repeats of one composition -- the only
    four-fold repeat in the archive -- so it is the bar every between-run
    result in this section has to clear. Returns the spread of the per-run
    medians, as a ratio for the clock and a difference for the depth.
    """
    if experiments is None:
        experiments = scope.REPLICATE_RUNS
    table = lag_window_frame(experiments, window)
    live = table[table.live]
    runs = live.groupby("experiment").agg(
        clock=("lag_half_s", "median"), depth=("lag_depth", "median"),
        curves=("lag_depth", "size"))
    clocks = runs.clock.to_numpy(dtype=float)
    positive = clocks[clocks > 0]
    return {"runs": int(len(runs)), "curves": int(len(live)),
            "window_s": float(table.window_s.max()),
            "clock_ratio": float(positive.max() / positive.min())
                           if len(positive) > 1 else np.nan,
            "clock_range": (float(clocks.min()), float(clocks.max())),
            "depth_range": (float(runs.depth.min()), float(runs.depth.max())),
            "depth_spread": float(runs.depth.max() - runs.depth.min()),
            "table": runs}


def matched_pair(pair, window=None):
    """
    Two runs matched on everything but one variable, at a common window.

    The archive's substrate, buffer-salt and catalyst contrasts are all pairs
    (`scope.SUBSTRATE_PAIRS`, `scope.BUFFER_TYPE_PAIR`, `scope.ENZYME_PAIRS`),
    because none of those three was ever laddered. A pair cannot give an order;
    it gives a direction and a size, and `replicate_floor` says whether the size
    means anything -- which for the clock is a factor of 2.1 between runs that
    differ in nothing at all.
    """
    table = lag_window_frame(tuple(pair), window)
    live = table[table.live]
    runs = live.groupby("experiment").agg(
        substrate=("substrate", "first"), buffer=("buffer", "first"),
        pH=("pH", "first"), e0=("e0", "first"), temperature=("temperature", "first"),
        clock=("lag_half_s", "median"), depth=("lag_depth", "median"),
        peak=("lag_peak", "median"), lag_first=("lag_first", "sum"),
        curves=("lag_depth", "size"))
    return {"pair": tuple(pair), "window_s": float(table.window_s.max()),
            "table": runs.reindex([e for e in pair if e in runs.index])}


def saturation_fraction(slope, stderr, per_ph=True):
    """
    A measured lag coefficient as the saturation fraction it implies.

    Both bounded schemes of section 4a read the same way. For a trap -- a
    species X that holds the catalyst OFF the activation path --

        1/tau = k_act/(1 + K[X])   =>   d ln tau/d ln[X] = K[X]/(1 + K[X])

    and for a species whose BOUND form is what activates, the same expression
    with the sign reversed. So the coefficient IS the fraction of the catalyst
    held by X, between 0 and 1, with no rate constant anywhere in it. A pH
    coefficient is a coefficient in an axis that moves by a decade per unit, so
    it is divided by ln 10 first.

    Returns the fraction, its error, and whether it lies inside (0, 1) -- which
    is the whole test, since a value outside falsifies both schemes.
    """
    fraction = slope / LN10 if per_ph else slope
    error = stderr / LN10 if per_ph else stderr
    return {"fraction": float(fraction), "stderr": float(error),
            "inside_trap": bool(TRAP_BOUND[0] < fraction < TRAP_BOUND[1]),
            "inside_activating": bool(
                ACTIVATING_BOUND[0] < fraction < ACTIVATING_BOUND[1])}


LAG_AXES = ("temperature", "pH", "buf", "buffer", "s0", "substrate", "h2o2")


def lag_identifiability(frame=None):
    """
    For each of the seven variables: does the archive move it inside a run?

    THE TABLE THIS SECTION EXISTS FOR. A variable that moves only between runs
    can only be read through a statistic with no window in it, and then only
    against the schedule and the signal. Returns one row per axis with the
    number of runs that move it internally, the number of distinct values it
    takes across the archive, and the widest within-run span.
    """
    if frame is None:
        frame = scope.frame(scope.archive())
    live = frame[frame.live]
    rows = []
    for axis in LAG_AXES:
        groups = live.groupby("experiment")[axis]
        internal = int((groups.nunique() > 1).sum())
        # A ratio for a concentration and a difference for the two axes that
        # are already logarithms or already differences -- pH and temperature.
        # Reporting a pH "span" as max/min is how a constant-pH run comes back
        # as 1.0 and reads like a decade.
        if axis in ("buffer", "substrate"):
            span, unit = np.nan, "category"
        elif axis in ("pH", "temperature"):
            span = float(groups.apply(lambda s: s.max() - s.min()).max())
            unit = "units" if axis == "pH" else "C"
        else:
            span = float(groups.apply(
                lambda s: s.max() / s.min() if s.min() > 0 else np.nan).max())
            unit = "fold"
        rows.append({"axis": axis,
                     "runs_moving_it": internal,
                     "runs": int(live.experiment.nunique()),
                     "levels": int(live[axis].nunique()),
                     "widest_within_run": span, "unit": unit})
    return pd.DataFrame(rows).set_index("axis")


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

    print("\n3a. AND IN THE UNITS THE PRODUCT HYPOTHESIS IS STATED IN")
    recovery = product_recovery()
    print(f"   a planted clock reads back at {recovery['clock_reads']:+.3f}, "
          f"a planted threshold at {recovery['threshold_range'][0]:+.2f} to "
          f"{recovery['threshold_range'][1]:+.2f}")
    print(recovery["table"].to_string(index=False))
    for name in ("4OMe catalysed", "BnOH two-axis (135-151)",
                 "temperature series"):
        got = product_at_landmark(named[name])
        if "slope" not in got:
            continue
        print(f"   {name:26s} d log(made)/d log(v) "
              f"{got['slope']:+.3f} +- {got['stderr']:.3f}   "
              f"{got['clock_sigma']:.1f} sigma from a clock, "
              f"{got['threshold_sigma']:.1f} from a threshold   "
              f"{got['curves']} curves, {got['dropped']} dropped")
        print(f"   {'':26s} pooled sd: log t {got['spread_time']:.3f}, "
              f"log product {got['spread_product']:.3f}, "
              f"log conversion {got['spread_conversion']:.3f}")

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

    print("\n7. THE SEVEN AXES, AND WHERE THE ARCHIVE IDENTIFIES EACH")
    print(lag_identifiability().to_string())

    print("\n7a0. THE LAG BY SUBSTRATE AND CHANNEL, over the whole archive")
    print(lag_channel_table().to_string())

    print("\n7a. THE REPLICATE FLOOR -- what two identical runs do anyway")
    floor = replicate_floor()
    print(floor["table"].to_string())
    print(f"   clock {floor['clock_ratio']:.2f}x over {floor['runs']} runs, "
          f"depth spread {floor['depth_spread']:.3f}")

    print("\n7b. WITHIN RUNS: the axes the archive steps inside a cuvette set")
    archive = scope.frame(scope.archive())
    named = induction_blocks(archive)
    cases = (("4OMe catalysed", named["4OMe catalysed"], ("s0", "h2o2", "buf")),
             ("BnOH two-axis", named["BnOH two-axis (135-151)"], ("s0", "h2o2")),
             ("4OMe peroxide 127-131",
              archive[archive.experiment.isin(PEROXIDE_LEVER)], ("s0", "h2o2")),
             ("buffer titrations",
              archive[archive.experiment.isin(scope.BUFFER_TITRATIONS)],
              ("s0", "h2o2", "buf")))
    for label, block, terms in cases:
        got = lag_orders(block, terms=terms)
        print(f"   {label}  n={got['points']} runs={got['experiments']}")
        for response in LAG_ORDER_FLOORS:
            row = got[response]
            pieces = []
            for term in terms:
                order = row[term]
                if not np.isfinite(order["order"]):
                    continue
                pieces.append(
                    f"{term} {order['order']:+.3f}+-{order['stderr']:.3f}"
                    f" -> {order['controlled']:+.3f}+-"
                    f"{order['controlled_stderr']:.3f}"
                    f" (r={got['signal_collinearity_' + term]:+.2f})")
            print(f"      {response:11s} floored {row['floored']:3d}/{row['n']:3d}"
                  f"  " + "  ".join(pieces))
    print("\n7c. BETWEEN RUNS: the four pH ladders, at a window they share")
    for fraction in LAG_WINDOW_SWEEP:
        ladders = lag_ph_ladders(window_fraction=fraction)
        for row in ladders:
            clock = row.get("lag_half_s")
            if clock is None:
                continue
            print(f"   {row['ladder']:28s} f={fraction:.2f} "
                  f"W={row['window_s']:6.0f}s runs={row['runs']:2d} "
                  f"n={row['points']:3d}  "
                  f"r(pH,length) {row['schedule_collinearity']:+.2f} "
                  f"r(pH,signal) {row['signal_collinearity']:+.2f}  "
                  f"clock {clock['slope']:+.3f} +-{clock['stderr']:.3f} -> "
                  f"{clock['controlled']:+.3f} +-{clock['controlled_stderr']:.3f}")
        pooled = pooled_ladder(ladders, "lag_half_s")
        print(f"   {'POOLED':28s} f={fraction:.2f} "
              f"{pooled['pooled']:+.3f} +- {pooled['stderr']:.3f}  "
              f"chi2 {pooled['chi2']:.2f} on {pooled['dof']}")

    print("\n7d. TEMPERATURE, at a window every run in the series shares")
    for row in lag_arrhenius_sweep():
        print(f"   x{row['fraction']:.2f} W={row['window_s']:6.0f}s "
              f"n={row['points']:3d} T={row['temperatures']}  "
              f"Ea {row['activation_kj']:6.1f} +- {row['stderr_kj']:.1f} kJ/mol")

    print("\n7e. THE PAIRS: buffer salt, substrate, catalyst")
    for pair in (scope.BUFFER_TYPE_PAIR, *scope.SUBSTRATE_PAIRS,
                 *scope.ENZYME_PAIRS):
        got = matched_pair(pair)
        print(f"   {str(pair):10s} window {got['window_s']:6.0f} s")
        print(got["table"].to_string())
    return table


def main():
    report()


if __name__ == "__main__":
    main()
