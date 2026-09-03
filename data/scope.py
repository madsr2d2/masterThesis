"""
The one way to get at a block of experiments.

Fitting is scoped to exps 135-151 (see fit_dataset.TWO_AXIS_BLOCK and
FITTING.md). This module is the front door to them: it loads the curves, hangs
every derived per-curve quantity off them in one frame, and answers the
questions that keep getting asked about the block's design.

WHY THIS EXISTS. Every analysis in this repository's history has re-derived
its own initial rate, its own noise estimate, its own lag position -- and the
copies drifted, most damagingly in the lag statistic, whose two versions
disagreed on 96 of 402 curves. The measurements now live in curve_metrics and
the selection lives here. If an analysis needs a quantity this module does not
expose, ADD IT HERE rather than computing it in a script; that is the whole
point of the module.

    from scope import curves, frame, design

    cs = curves()                     # 119 Curve objects, exps 135-151
    df = frame()                      # one row per curve, derived columns filled
    design()                          # the block's design, as a table

    python data/scope.py              # print the summary
    python data/scope.py --design     # the per-experiment design table
"""
import itertools

import numpy as np
import pandas as pd

from curve_metrics import (ACCELERATION_SIGMA, BUBBLE_DROP_SIGMA,
                           INITIAL_WINDOW, LAG_THRESHOLD,
                           OUTLIER_SIGMA, acceleration, bubble_drops, bubble_load,
                           bubble_rate,
                           burst_amplitude, debubble, detachments,
                           initial_rate,
                           isolated_outliers, lag_time, local_outlier_z,
                           model_residual, monotone_bound, peak_position,
                           peak_rate,
                           quadratic_rate, segmented_fit, segment_selection,
                           SEGMENT_RATIO_STEEP,
                           whole_slope)
from fit_dataset import (TWO_AXIS_BLOCK, TWO_AXIS_GROUP, build_curves,
                         in_block, source_floor)
from solution_chemistry import dominant_buffer_pair
from summary_kinetics import fit_burst_bounded, fit_progress

# A run's own cuvettes have to move an axis by at least this much before that
# axis counts as measured inside the run rather than across experiments.
LADDER_MINIMUM = 2.0

# Net change below this multiple of a curve's own noise is not a measurement of
# anything. Exps 150 and 151 are mostly flat by this rule. That does NOT make
# them a background: every run in TWO_AXIS_BLOCK carries enzyme, and their
# cuvettes carry no concentration information either (concentration_agreement
# 0.61 and 0.005). The block has no enzyme-free control -- see
# DATA_VERIFICATION.md, 2026-08-31.
LIVE_SIGNAL_NOISE_MULTIPLE = 20.0


def curves_of(experiment):
    """Every curve of one experiment. A named alias so functions that
    take a `scope` argument can still reach the curve list."""
    return list(curves((int(experiment),)))


def curves(scope=TWO_AXIS_BLOCK):
    """The block's Curve objects, in (experiment, sample) order."""
    all_curves, _ = build_curves()
    return sorted(in_block(all_curves, scope),
                  key=lambda c: (c.experiment, c.sample))


def frame(scope=TWO_AXIS_BLOCK):
    """
    One row per curve of the block, with every derived quantity attached.

    Columns: experiment, source, sample, substrate, buffer, temperature, buf,
    pH, s0, h2o2, e0, hoo, duration_s,
    points, outliers, outliers_in_runs, first_point_z, first_point_flagged,
    noise, net, live, v0, v0_stderr, v0_rms, vmax, vmax_stderr, vmax_where,
    bubble_drops, bubble_load, vmax_corrected, vmax_monotone,
    gain, vmax_time_s, lag_time_s, conversion, peak, lags, accel_z, accel_where, accelerates,
    late_over_early.

    `v0` is the rate before the catalyst has built up and `vmax` the rate
    after; on this block they are different measurements with different
    orders, so pick deliberately rather than reaching for v0 by habit.

    `accelerates` is the autocatalysis verdict; prefer it to `lags` on this
    block, whose curves are small-amplitude enough that the point-wise
    gradient behind `lags` is part noise. See curve_metrics.acceleration.

    Reach for this before writing a loop over `curves()`. The columns are named
    for what they are, so a question like "does the substrate order hold at low
    peroxide" is a groupby on this frame rather than a new script.
    """
    rows = []
    for curve in curves(scope):
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        # curve.noise, not a fresh curve_noise call: build_curves floors it
        # by the curve's SOURCE, and a .rre curve floored at the .txt
        # export's quantisation reports 2.4x its real noise.
        noise = curve.noise
        # The same floor has to reach the RATES, not just the noise. Every one
        # of these divides by a standard error that line_fit floors, and until
        # 2026-09-01 that floor was hardcoded at the export's quantisation for
        # every curve -- suppressing the acceleration z on the .rre data this
        # scope is entirely made of. See fit_dataset.source_floor.
        floor = source_floor(curve.source)
        net = float(values[-1] - values[0])
        v0, v0_stderr, v0_rms = initial_rate(times, values, floor=floor)
        peak = peak_position(values, times)
        accel_z, accel_where = acceleration(times, values, floor=floor)
        vmax, vmax_stderr, vmax_where = peak_rate(times, values, floor=floor)
        # THE GAS. O2 from the catalysed decomposition of the peroxide grows on
        # the window and detaches, so the readings carry a sawtooth that no
        # kinetic form can hold. `debubble` splits the readings into a
        # non-decreasing chemistry and a gas made at a steady rate, and
        # `monotone_bound` brackets it from the assumption-free side;
        # `bubble_load` says which of the two a curve is entitled to.
        # curve_metrics.bubble_profile has the case against stitching, which
        # is the repair that suggests itself and the one that makes it worse.
        drops = bubble_drops(values, noise)
        events = detachments(values, noise)
        gas_rate = bubble_rate(times, values, events)
        corrected, _ = debubble(times, values, noise)
        vmax_corrected, _, _ = peak_rate(times, corrected, floor=floor)
        vmax_monotone, _, _ = peak_rate(
            times, monotone_bound(values), floor=floor)
        # Three more rate estimators, so that "does this conclusion depend on
        # how the rate was measured" is a groupby rather than an argument.
        # v0 uses the first 20% of the run, v0_whole every point with no bend
        # allowed, v0_quad every point with one bend allowed.
        v0_whole, v0_whole_stderr = whole_slope(times, values, floor=floor)
        v0_quad, v0_quad_stderr, curvature_t = quadratic_rate(
            times, values, floor=floor)
        # And the burst/lag form, the only estimator here that is a SHAPE
        # rather than a polynomial: A = c + v_ss t - B(1 - exp(-t/tau)), with
        # v0 = v_ss - B/tau. It is carried with three guards, because its v0
        # is the one number in this frame that can look confident while
        # meaning nothing:
        #
        #   v0_burst_bounded    the profile-likelihood interval over tau is
        #                       narrow enough to quote at all
        #   v0_burst_kind       burst / lag / clamped / unresolved. On a LAG
        #                       v0 is the INDUCTION rate, not the maximum, so
        #                       pooling the two into one column mixes two
        #                       different quantities.
        #   v0_burst_resid      residual RMS in units of the curve's own
        #                       noise. `bounded` asks whether the PARAMETER is
        #                       determined, not whether the MODEL fits, and
        #                       those come apart: exp 65's four cuvettes all
        #                       report bounded=True on fits that collapse to
        #                       B = 0 with tau at the floor of its grid --
        #                       a straight line wearing a four-parameter form
        #                       -- and sit at 7-8x noise.
        # The nested pair: one relaxation onto a steady rate, or two. The
        # one-phase form's rate is MONOTONE, so it cannot hold a curve whose
        # rate rises to a maximum and then falls -- 14 of the temperature
        # series' 24 do exactly that. `fit_progress` fits both and returns
        # whichever the curve earns on an F test. See summary_kinetics.
        # No `floor` argument: fit_progress's `floor` is the TAU GRID start as a
        # fraction of the run, not a noise floor. Passing the source floor here
        # set the grid to [span, 2*span] and made every curve look two-phase.
        progress = fit_progress(times, values)
        peak_fitted, peak_fitted_time = progress.peak_rate
        # How much the curve put on before it stopped rising, off the FITTED
        # curve rather than the readings -- curve_metrics.burst_amplitude says
        # why, and `burst_bounded` is the flag that keeps an unfinished rise
        # out of any comparison between runs of different length.
        burst_rise, burst_at, burst_bounded = burst_amplitude(
            times, progress.predict(times))
        burst = fit_burst_bounded(times, values, noise_floor=floor)
        burst_pred = (burst.c + burst.v_ss * times
                      - burst.B * (1.0 - np.exp(-times / burst.tau)))
        burst_resid = model_residual(values, burst_pred, 4, noise)
        # The quadratic scored the same way, so the two are comparable. Its
        # own docstring has warned since it was written that it misses by
        # 2-7x the noise where the deceleration is strong; this puts a number
        # on it per curve instead of leaving it as a caveat.
        quad_design = np.column_stack([np.ones(len(times)), times, times ** 2])
        quad_beta, *_ = np.linalg.lstsq(quad_design, values, rcond=None)
        v0_quad_resid = model_residual(values, quad_design @ quad_beta, 3, noise)
        # Suspect readings. `outliers` counts the ISOLATED ones -- a single
        # reading out of line with both neighbours, which nothing chemical can
        # produce at this sampling rate. `outliers_in_runs` counts flagged
        # readings with a flagged neighbour, which may be real structure and
        # are never treated as artefacts. Neither excludes anything.
        break_time, slope_before, slope_after, break_ratio = segmented_fit(
            times, values, floor=floor)
        # And the same search allowed TWO breakpoints. A rate that rises and
        # then falls has two, and a one-break search lands on whichever is
        # stronger and never reports the other -- which is how the early break
        # on the 40 C curves went unseen until it was spotted on a plot.
        segments = segment_selection(times, values)
        isolated, in_runs = isolated_outliers(times, values, noise)
        # The leading reading gets its own flag, taken from z[0] rather than
        # from membership of `isolated`: a bad leading point drags its
        # neighbour past the threshold and the pair then reads as a run. That
        # happened on 21 of the 86 flagged curves before the first-reading drop
        # and on 8 of the 56 after it.
        #
        # The flag is kept although its statistical case has weakened. The
        # instrument's own first reading was flagged on 21.4% of curves against
        # 14.7% for the last; since that reading is discarded the leading one
        # is flagged on 14.7% against 16.2%, so it is no longer the outlier
        # class it was. What survives is structural: t = 0 is where v0 is
        # extrapolated to, so a bad point there costs more than anywhere else,
        # and the run-masking above still hides it from `isolated`.
        outlier_z = local_outlier_z(times, values, noise)
        first_z = float(outlier_z[0]) if len(outlier_z) else np.nan
        buf_acid, buf_base, buf_pka = dominant_buffer_pair(
            curve.buffer, curve.pH, curve.buf)
        rows.append({
            "experiment": curve.experiment,
            "source": curve.source,
            "sample": curve.sample,
            # The block a curve belongs to, carried per row rather than
            # assumed: once `scope` is a free parameter a frame may span
            # several (substrate, temperature, buffer) cells, and a caller
            # that pools across one without meaning to has no way to notice.
            "substrate": curve.substrate,
            "buffer": curve.buffer,
            "temperature": curve.temperature,
            # Buffer CONCENTRATION, mM -- distinct from `buffer`, the salt.
            # In every enzyme-free titration this falls as `s0` rises, because
            # substrate volume displaced buffer volume; see BUFFER_CONFOUNDED.
            "buf": curve.buf,
            # The conjugate pair straddling this pH, in mM. General acid/base
            # catalysis and a peroxo-adduct route are both first order in a
            # SPECIES, not in the total -- see
            # solution_chemistry.dominant_buffer_pair, which also says why
            # this dataset cannot separate them.
            "buf_acid": buf_acid,
            "buf_base": buf_base,
            "buf_pka": buf_pka,
            "pH": curve.pH,
            "s0": curve.conditions.s0,
            "h2o2": curve.conditions.h2o2,
            "e0": curve.conditions.e0,
            "hoo": curve.conditions.hoo,
            "duration_s": float(times[-1] - times[0]),
            "points": len(times),
            "noise": noise,
            "net": net,
            "live": net > LIVE_SIGNAL_NOISE_MULTIPLE * noise,
            "v0": v0,
            "v0_stderr": v0_stderr,
            "v0_rms": v0_rms,
            "v0_whole": v0_whole,
            "v0_whole_stderr": v0_whole_stderr,
            "v0_quad": v0_quad,
            "v0_quad_stderr": v0_quad_stderr,
            "curvature_t": curvature_t,
            "v0_burst": burst.v0,
            "v0_burst_low": burst.v0_low,
            "v0_burst_high": burst.v0_high,
            "v0_burst_bounded": bool(burst.bounded),
            "v0_burst_kind": burst.kind,
            "v0_burst_resid": burst_resid,
            # How many relaxation phases the curve earned, and the evidence.
            "phases": int(progress.phases),
            "two_phase_f": float(progress.f_statistic),
            # WHICH WAY THE CURVE POINTS. The sign convention is B > 0 for a
            # LAG -- the rate starts below its eventual value and rises -- and
            # B < 0 for a burst. Both forms carry it and until 2026-09-02 this
            # frame carried neither: the time constants were here, the peak
            # rate was here, and nothing said whether the curve began slow or
            # began fast, so no analysis in the package could ask. The two-axis
            # block splits almost evenly between the two (46 lag-first against
            # 45 burst-first of 110 live), which is the thing that could not
            # be seen. See summary_kinetics.ProgressFit.kind.
            "progress_kind": progress.kind,
            "B_fast": progress.amplitudes[0],
            "B_slow": progress.amplitudes[1],
            # THE SIZE OF THE EARLY RISE, which `B_fast` is not: the two-phase
            # solve trades amplitude between its exponentials without moving
            # the curve, so the split is unstable where the prediction is not.
            # `turnovers` divides it by the catalyst, and is the whole question
            # of whether the rise is one turnover or many.
            "burst": burst_rise,
            "burst_time_s": burst_at,
            "burst_bounded": burst_bounded,
            # THE RATE TO USE ON AN ARRHENIUS PLOT. The largest rate the fitted
            # model reaches: v_ss for a one-phase lag, the interior maximum for
            # a two-phase curve. It is what `vmax` measures off the raw
            # readings with the truncation taken out -- the block statistic can
            # only find a maximum inside the window, and at 15 and 20 C the
            # window ends before the rate levels off.
            #
            # Do NOT use the two-phase form's `v_ss` instead: it is the
            # t -> infinity asymptote, unconstrained once a decay is in the
            # fit, and it comes out NEGATIVE on two of the 35 C curves.
            "v_peak": float(peak_fitted),
            "v_peak_time": float(peak_fitted_time),
            "tau_fast": float(progress.two.tau1 if progress.phases == 2
                              else burst.tau),
            "tau_slow": float(progress.two.tau2 if progress.phases == 2
                              else np.nan),
            "tau_slow_resolved": bool(progress.phases == 2
                                      and progress.two.resolved),
            "progress_resid": model_residual(
                values, progress.predict(times),
                6 if progress.phases == 2 else 4, noise),
            # The burst form's ASYMPTOTIC rate, after the lag is over. On a lag
            # curve this is the quantity `vmax` is trying to measure and
            # cannot when the run ends before the rate levels off -- which is
            # the whole 15-20 C problem in the temperature series. Carried
            # beside v0_burst (the INDUCTION rate) because on a lag curve the
            # two are the opposite ends of the same fit and confusing them is
            # the error curve_metrics.peak_rate exists to prevent.
            "v_ss": float(burst.v_ss),
            "v0_quad_resid": v0_quad_resid,
            "tau": burst.tau,
            "tau_resolved": bool(burst.tau_resolved),
            "outliers": len(isolated),
            "outliers_in_runs": len(in_runs),
            "first_point_z": first_z,
            "first_point_flagged": bool(np.isfinite(first_z)
                                        and abs(first_z) > OUTLIER_SIGMA),
            "vmax": vmax,
            "vmax_stderr": vmax_stderr,
            "vmax_where": vmax_where,
            # THE SAME RATE WITH THE GAS TAKEN OUT, and the bound that needs no
            # bubble model at all. Quote the pair: `vmax_corrected` is the
            # estimate and the gap to `vmax_monotone` is its systematic, the
            # way product_fate quotes the sink's activation energy. Neither is
            # a substitute for reading `bubble_load` first -- above 1 the
            # artefact carries more absorbance than the reaction and no repair
            # here recovers the rate.
            "bubble_drops": int(len(drops)),
            "bubble_events": int(len(events)),
            "bubble_load": bubble_load(values, drops),
            "gas_rate": gas_rate,
            "vmax_corrected": vmax_corrected,
            "vmax_monotone": vmax_monotone,
            # vmax/v0: how many times over the reaction sped up. accel_z
            # says whether the speed-up is real, this says how big it is --
            # z also carries the curve's noise and length, so z alone is not
            # comparable between cuvettes.
            "gain": vmax / v0 if np.isfinite(v0) and v0 > 0 else np.nan,
            # Fraction of the cuvette's substrate that has turned over by the
            # end of the run. This is the axis an autocatalysis driven by
            # PRODUCT lives on: a loop that feeds on product closes at a
            # given product/substrate ratio, not at a given absolute rate.
            # When the curve is steepest, in seconds. The FRACTION is what
            # `vmax_where` reports, and fractions are not comparable between
            # runs of 51 and 480 minutes; a mechanism predicts a time.
            "vmax_time_s": vmax_where * float(times[-1] - times[0]),
            "lag_time_s": lag_time(times, values, floor=floor),
            # The extinction coefficient the sheet declares for this run,
            # carried so a rate in AU/s can be turned into one in mM/s. That is
            # what an ABSOLUTE rate constant needs, and therefore what an
            # Eyring entropy needs: a slope survives any constant factor, an
            # intercept does not. `conversion` below already relies on it.
            "epsilon": float(curve.epsilon),
            "conversion": (net / (curve.epsilon * curve.conditions.s0)
                           if curve.epsilon > 0 and curve.conditions.s0 > 0
                           else np.nan),
            # THE EARLY RISE IN UNITS OF THE CATALYST. Below one, the run has
            # not turned the catalyst over even once and the rise is a
            # stoichiometric burst; well above it, the catalyst is being
            # regenerated. Both epsilon and [enz] come off the sheet, so this
            # is as good as the weighing -- `verify_enzyme_stock.py` is the
            # independent check on the second of them.
            "turnovers": (burst_rise / (curve.epsilon * curve.conditions.e0)
                          if curve.epsilon > 0 and curve.conditions.e0 > 0
                          else np.nan),
            "peak": peak,
            "lags": bool(peak > LAG_THRESHOLD) if np.isfinite(peak) else False,
            "accel_z": accel_z,
            "accel_where": accel_where,
            "accelerates": bool(accel_z > ACCELERATION_SIGMA)
            if np.isfinite(accel_z) else False,
            "late_over_early": _late_over_early(times, values),
            # WHETHER THIS CURVE CAN SHOW A BACKGROUND FEATURE AT ALL.
            # Every run is double-beam and what the reference channel omits
            # decides what the curve means (kinetics_io, DATA_VERIFICATION.md
            # 2026-08-31): an enzyme run's reference omits the ENZYME, so the
            # background appears in both beams and CANCELS, and the curve is a
            # catalytic increment. A background run's reference omits the
            # H2O2, so the curve is the raw reaction with nothing subtracted.
            # A background shape can therefore only be looked for on
            # `differential == False` rows, and comparing the two populations
            # for one is a category error -- which is how the boric probe's
            # first control set was chosen wrongly on 2026-09-01.
            "differential": bool(curve.conditions.e0 > 0),
            # The two-line split. Every other shape column here compares the
            # curve's start to its end, and a curve that breaks upward in the
            # MIDDLE and then plateaus defeats all of them -- exp 65 sat
            # mid-pack on `late_over_early` while steepening 1.8-15.9x across
            # a break its four cuvettes share to 56 s. `break_ratio` is the
            # one to read; see curve_metrics.segmented_fit and
            # `synchronised_break`.
            "break_time": break_time,
            "slope_before": slope_before,
            "slope_after": slope_after,
            "break_ratio": break_ratio,
            # Every breakpoint the curve earns, and what the slopes DO across
            # them. Read `break_pattern`, not `break_count`: three lines fit a
            # smooth bend better than two whether or not anything happened, so
            # the count is about approximation and the pattern is about the
            # curve. See curve_metrics.segment_selection.
            "break_count": int(segments["breaks"]),
            "break_times": tuple(round(float(v), 1) for v in segments["times"]),
            "break_pattern": segments["pattern"],
            "break_f": float(segments["f_statistic"]),
        })
    return pd.DataFrame(rows)


def _late_over_early(times, values, fraction=0.2):
    """
    Late-window slope divided by early-window slope.

    Above 1 the curve is accelerating, which is the autocatalysis signature.
    Kept separate from `peak_position` because it answers a different question:
    peak_position asks *where* the curve is steepest, this asks whether it is
    still getting steeper at the end.
    """
    count = max(4, int(len(times) * fraction))
    early = np.polyfit(times[:count], values[:count], 1)[0]
    late = np.polyfit(times[-count:], values[-count:], 1)[0]
    return float(late / early) if early > 0 else np.nan


def ladder(axis, scope=TWO_AXIS_BLOCK):
    """
    How far each run moves `axis` within its own cuvettes, as a factor.

    `axis` is a column of `frame()` -- "s0", "h2o2", "hoo", "e0".
    """
    data = frame(scope)
    return data.groupby("experiment")[axis].agg(
        lambda values: float(values.max() / max(values.min(), 1e-12)))


def within_experiment_share(axis, scope=TWO_AXIS_BLOCK):
    """
    The fraction of an axis's log-variance that lives inside experiments.

    This is the number the scope was chosen on. Near 1 means the order in that
    axis is measured within runs and cannot be absorbed by a per-experiment
    offset; near 0 means it rests entirely on comparing one run to another.
    """
    data = frame(scope)
    logged = np.log(np.maximum(data[axis].to_numpy(dtype=float), 1e-12))
    total = logged.var()
    if total <= 0:
        return 0.0
    groups = data.experiment.to_numpy()
    within = np.average(
        [logged[groups == e].var() for e in np.unique(groups)],
        weights=[int((groups == e).sum()) for e in np.unique(groups)])
    return float(within / total)


def design(scope=TWO_AXIS_BLOCK):
    """One row per experiment: what it varies, how long it ran, what it saw."""
    data = frame(scope)
    rows = []
    for experiment, group in data.groupby("experiment"):
        rows.append({
            "experiment": experiment,
            "pH": group.pH.iloc[0],
            "hoo_mM": group.hoo.iloc[0],
            "cuvettes": len(group),
            "s0_ladder": group.s0.max() / max(group.s0.min(), 1e-12),
            "h2o2_ladder": group.h2o2.max() / max(group.h2o2.min(), 1e-12),
            "duration_min": group.duration_s.max() / 60.0,
            "live": int(group.live.sum()),
            "lagging": int(group.lags.sum()),
            # Live curves only. A dead curve's "acceleration" is its
            # quantisation staircase happening to step late: exp 151 has
            # one live curve and would otherwise report two accelerating.
            "accelerating": int(group.loc[group.live, "accelerates"].sum()),
            "median_late_over_early": float(group.loc[group.live,
                                                      "late_over_early"].median()),
        })
    return pd.DataFrame(rows).set_index("experiment")


def blocks(scope=TWO_AXIS_BLOCK):
    """
    The (substrate, temperature, buffer) cells a scope spans, with counts.

    Rate constants may be pooled only within one cell (FITTING.md F7), so a
    scope that returns more than one row here is a scope no fit may be run on
    as a unit. TWO_AXIS_BLOCK returns exactly one; the enzyme-free BnOH set
    returns two, phosphate and boric, which is why it is a background
    characterisation and not a fit.
    """
    data = frame(scope)
    return (data.groupby(["substrate", "temperature", "buffer"])
            .agg(curves=("experiment", "size"),
                 experiments=("experiment", "nunique"))
            .sort_values("curves", ascending=False))


def summary(scope=TWO_AXIS_BLOCK):
    """The scope in one paragraph of numbers, for printing."""
    data = frame(scope)
    cells = blocks(scope)
    return {
        # Read off the curves, not assumed: with `scope` a free parameter this
        # was reporting TWO_AXIS_GROUP for every scope it was handed.
        "block": cells.index[0] if len(cells) == 1 else tuple(cells.index),
        "experiments": int(data.experiment.nunique()),
        "curves": len(data),
        "live_curves": int(data.live.sum()),
        "pH_range": (float(data.pH.min()), float(data.pH.max())),
        "hoo_decades": float(np.log10(data.hoo.max() / max(data.hoo.min(), 1e-12))),
        "within_experiment_s0": within_experiment_share("s0", scope),
        "within_experiment_h2o2": within_experiment_share("h2o2", scope),
        "lagging": int(data.lags.sum()),
        "accelerating": int(data.accelerates.sum()),
        "accelerating_live": int(data.loc[data.live, "accelerates"].sum()),
    }


# The other axis of a two-axis run is held constant along each ladder, so a
# ladder is the smallest set of cuvettes that isolates one concentration.
AXIS_PARTNER = {"s0": "h2o2", "h2o2": "s0"}


def ladder_groups(axis, scope=TWO_AXIS_BLOCK, live_only=True, minimum=3):
    """
    The runs of cuvettes that vary `axis` alone, as (label, sub-frame) pairs.

    Within one experiment, cuvettes sharing a value of the partner axis differ
    only in `axis` -- and in nothing else at all, since pH, [HOO-], enzyme,
    cell and day are properties of the run. A slope read along one of these is
    therefore an order, not a correlation. This is the same argument the scope
    rests on, applied to a single ladder instead of the whole block.

    `minimum` is the fewest cuvettes a ladder must have to be worth a slope.
    """
    partner = AXIS_PARTNER[axis]
    data = frame(scope)
    if live_only:
        data = data[data.live]
    groups = []
    for (experiment, held), group in data.groupby(["experiment", partner]):
        group = group.sort_values(axis)
        if len(group) < minimum or group[axis].nunique() < minimum:
            continue
        groups.append((f"exp {experiment}, {partner}={held:.4g}", group))
    return groups


def ladder_trend(parameter, axis, scope=TWO_AXIS_BLOCK):
    """
    A parameter's median value at the bottom and top rung of the ladders.

    For quantities that cannot be logged -- `accel_z` is a signed statistic, not
    a rate -- an order is meaningless, but the question "does this track the
    axis" still has an answer. Taking the bottom and top rung of each ladder
    keeps the comparison inside a run, so it carries the same protection from
    between-run confounding that `orders` gets from its offsets.

    Returns (low, high, n_ladders).
    """
    low, high = [], []
    for _, group in ladder_groups(axis, scope):
        ordered = group.sort_values(axis)
        low.append(float(ordered[parameter].iloc[0]))
        high.append(float(ordered[parameter].iloc[-1]))
    if not low:
        return np.nan, np.nan, 0
    # nanmedian, not median: `gain` is nan wherever v0 was not positive, and
    # one such rung would otherwise erase the whole trend.
    return float(np.nanmedian(low)), float(np.nanmedian(high)), len(low)


def local_orders(parameter, axis, scope=TWO_AXIS_BLOCK):
    """
    The log-log slope between each ADJACENT pair of rungs, ladder by ladder.

    `orders` fits one slope to the whole block and so cannot tell saturation
    from inhibition: both flatten the average. This returns the order as a
    function of concentration instead. Saturation approaches zero and stays
    there; inhibition goes negative at the top of the range. That distinction
    is the whole question when a curve stops responding to its substrate.

    Returns a frame with columns: experiment, axis value (geometric mean of the
    pair), order, and the pH of the run.
    """
    rows = []
    for _, group in ladder_groups(axis, scope):
        group = group.sort_values(axis)
        x = np.log(group[axis].to_numpy(dtype=float))
        y = group[parameter].to_numpy(dtype=float)
        for i in range(len(group) - 1):
            if not (y[i] > 0 and y[i + 1] > 0) or x[i + 1] <= x[i]:
                continue
            rows.append({
                "experiment": int(group.experiment.iloc[i]),
                axis: float(np.exp(0.5 * (x[i] + x[i + 1]))),
                "order": float((np.log(y[i + 1]) - np.log(y[i]))
                               / (x[i + 1] - x[i])),
                "pH": float(group.pH.iloc[i]),
            })
    return pd.DataFrame(rows)


def concentration_agreement(scope=TWO_AXIS_BLOCK, parameter="vmax"):
    """
    Per experiment: does the rate follow the concentrations across its cuvettes?

    A run's cuvettes share pH, [HOO-], enzyme, cell and day, so the only thing
    that may move `parameter` between them is [S] and [H2O2]. This correlates
    each run's observed log rate against the log rate the block's own fitted
    orders predict for its cuvettes. Near +1 the run obeys the rate law; near 0
    it does not, and whatever is moving its cuvettes is not the reaction.

    THIS IS THE SCREEN FOR DRIFT-DOMINATED RUNS. A cuvette measured at the level
    of the cell's own wander still produces a rate and a standard error, and
    nothing else in this package would call it out.

    Mildly circular: the orders are fitted over the same curves, so a block of
    drift-dominated runs would flatten the orders they are then judged against.
    Read it as a ranking between runs, not as an absolute score.

    Returns a frame indexed by experiment: hoo, median rate, agreement, n.
    """
    data = frame(scope)
    data = data[data.live & (data[parameter] > 0)]
    fitted = orders(parameter, scope)
    rows = []
    for experiment, group in data.groupby("experiment"):
        if len(group) < 4:
            continue
        predicted = (fitted["order_s0"] * np.log(group.s0)
                     + fitted["order_h2o2"] * np.log(group.h2o2))
        observed = np.log(group[parameter].to_numpy(dtype=float))
        if predicted.std() <= 0 or observed.std() <= 0:
            continue
        rows.append({"experiment": experiment,
                     "hoo": float(group.hoo.iloc[0]),
                     "median_rate": float(group[parameter].median()),
                     "agreement": float(np.corrcoef(observed, predicted)[0, 1]),
                     "cuvettes": len(group)})
    return pd.DataFrame(rows).set_index("experiment")


# The per-curve quantities it is meaningful to ask an order of. `v0` and `vmax`
# are rates; `net` is an extent, whose "order" is a shape statement about how
# far a cuvette gets, not a rate law; `gain` is the speed-up factor, whose
# order says how the autocatalysis itself depends on each concentration.
ORDER_PARAMETERS = ("v0", "vmax", "net", "gain")

# A run has to predict its own cuvettes this well before its orders are worth
# reading. concentration_agreement correlates each run's observed log vmax
# against the log rate its own cuvette concentrations imply; a run that scores
# low is telling you its cuvettes differ by something other than what was put
# in them -- drift, a bad blank, or a signal too small to carry the ladder.
#
# The archive separates cleanly here, which is why the threshold is not fine-
# tuned: the eleven runs above it score 0.724 to 0.974, and the five below
# score 0.609 down to 0.005. Exp 151 does not appear at all -- its cuvettes
# scatter 234-fold with two negative rates.
AGREEMENT_FLOOR = 0.70


def strong_runs(scope=TWO_AXIS_BLOCK, floor=AGREEMENT_FLOOR):
    """
    The experiments whose own cuvettes predict their own rates.

    Returns a sorted tuple of experiment numbers. Orders quoted in MECHANISM.md
    are measured over these; quoting them over all 17 runs moves the substrate
    order of vmax from +0.01 to +0.11 and drops the fit's R2 from 0.88 to 0.81,
    because the runs that fail this test contribute scatter and no signal.
    """
    table = concentration_agreement(scope)
    return tuple(sorted(int(e) for e in
                        table.index[table.agreement >= floor]))


def _moves(column, groups, within):
    """Is this regressor identified alongside the fit's intercepts?"""
    if not within:
        return float(np.ptp(column)) > 0
    return any(float(np.ptp(column[groups == g])) > 0 for g in np.unique(groups))


def orders(parameter="v0", scope=TWO_AXIS_BLOCK, within=True, live_only=True,
           frame=None, floor=None):
    """
    Apparent reaction orders in [S] and [H2O2], from a log-log regression.

    Fits log(parameter) = intercept + a*log[S] + b*log[H2O2], with a free
    offset per experiment when `within` is True.

    `frame` lets a caller measure the orders of a column this module does not
    build -- `induction` passes its own table in, the way `arrhenius` accepts
    one -- and then `scope` is ignored. `floor` is the response's log floor:
    the default None drops non-positive rows, which is right for a rate and
    wrong for a quantity whose zero is a measurement rather than a failure.

    THE OFFSETS ARE THE POINT. A run holds pH, [HOO-], enzyme batch, cell and
    day constant across its own cuvettes, so a per-experiment offset absorbs
    all of them and leaves `a` and `b` measured only from contrast between
    cuvettes of the same run. That is exactly what this scope was chosen to
    provide -- 100.0% of its log[S] variance and 94.1% of its log[H2O2]
    variance is within-experiment. Set `within=False` to see what the same
    numbers look like when between-run differences are allowed to carry the
    fit; the gap between the two is the confounding the scope exists to avoid.

    Returns a dict: order_s0, stderr_s0, order_h2o2, stderr_h2o2, n, r2.
    """
    data = globals()["frame"](scope) if frame is None else frame
    if live_only:
        data = data[data.live]
    data = data[np.isfinite(data[parameter])]
    data = data[data[parameter] > 0] if floor is None else data[
        data[parameter] >= 0]
    if len(data) < 4:
        return {"order_s0": np.nan, "stderr_s0": np.nan, "order_h2o2": np.nan,
                "stderr_h2o2": np.nan, "n": len(data), "r2": np.nan}

    y = data[parameter].to_numpy(dtype=float)
    y = np.log(y if floor is None else np.maximum(y, floor))
    columns = [np.log(data.s0.to_numpy(dtype=float)),
               np.log(data.h2o2.to_numpy(dtype=float))]
    # An axis the block does not move is not identified, and saying so is not
    # optional: with `within=True` a constant axis is collinear with the
    # experiment indicators, `lstsq` splits the coefficient between them
    # through the pseudo-inverse, and the number that comes back looks like a
    # measurement. Exps 127-131 hold [S] at 9.47 mM on every cuvette and
    # returned an order of +2.14 +- 0.18 in it before this guard existed.
    #
    # WITH OFFSETS, "DOES IT MOVE" MEANS "DOES IT MOVE INSIDE A RUN". Testing
    # the whole column instead was the same failure one level down, and it was
    # live until 2026-09-03: the two-axis block's L splits into an arm that
    # ladders [S] at fixed [H2O2] and an arm that ladders [H2O2] at fixed [S]
    # (induction.ladder_arms), and inside each arm the OTHER axis is constant
    # per run and varies only between runs -- 73.4 against 35.2 mM, 10.816
    # against 8.759. That clears a global ptp test, so the offsets absorbed it
    # and the pseudo-inverse handed back a substrate order of -6.067 +- 0.075
    # for the peroxide arm: tight, confident and meaningless.
    identified = [_moves(column, data.experiment.to_numpy(), within)
                  for column in columns]
    columns = [column for column, keep in zip(columns, identified) if keep]
    if within:
        # One indicator per experiment, and no separate intercept -- the
        # indicators already span it. Adding both would make X rank-deficient.
        labels = np.unique(data.experiment.to_numpy())
        columns += [(data.experiment.to_numpy() == e).astype(float)
                    for e in labels]
    else:
        columns.append(np.ones(len(data)))
    design_matrix = np.column_stack(columns)

    if not columns:
        return {"order_s0": np.nan, "stderr_s0": np.nan, "order_h2o2": np.nan,
                "stderr_h2o2": np.nan, "n": len(data), "r2": np.nan}
    coefficients, *_ = np.linalg.lstsq(design_matrix, y, rcond=None)
    residual = y - design_matrix @ coefficients
    degrees = max(1, len(y) - np.linalg.matrix_rank(design_matrix))
    variance = float(residual @ residual) / degrees
    covariance = variance * np.linalg.pinv(design_matrix.T @ design_matrix)
    total = float(((y - y.mean()) ** 2).sum())
    out = {"order_s0": np.nan, "stderr_s0": np.nan,
           "order_h2o2": np.nan, "stderr_h2o2": np.nan,
           "n": int(len(y)),
           "r2": float(1 - (residual @ residual) / total) if total > 0
           else np.nan}
    position = 0
    for name, keep in zip(("s0", "h2o2"), identified):
        if not keep:
            continue
        out[f"order_{name}"] = float(coefficients[position])
        out[f"stderr_{name}"] = float(np.sqrt(covariance[position, position]))
        position += 1
    return out


def order_table(scope=TWO_AXIS_BLOCK, parameters=ORDER_PARAMETERS):
    """Orders for each parameter, within-experiment and pooled, side by side."""
    rows = []
    for parameter in parameters:
        for within in (True, False):
            result = orders(parameter, scope, within=within)
            rows.append({"parameter": parameter,
                         "fit": "within-experiment" if within else "pooled",
                         **result})
    return pd.DataFrame(rows).set_index(["parameter", "fit"])




# ---------------------------------------------------------------------------
# The pH ladder inside the two-axis block.
#
# The block is usually described by what a run varies -- both concentration
# axes, seven cuvettes, an L rather than a grid (induction.ladder_arms). That
# describes ONE run. What the seventeen runs are to each other is a second
# design, and it is the stronger one. The block holds only TWO composition sets
# in sixteen runs: exps 136-142 and exps 143-151, which share the substrate
# ladder exactly (0.216, 0.865, 3.028, 10.816 mM) and differ only in their four
# peroxide levels. So inside either set a cuvette is matched one for one across
# the runs of a pH ladder running 5.47 to 9.73.
#
# Which makes the pH axis measurable the way the concentration axes are. An
# order in [HOO-] read across runs is normally hostage to everything else that
# differs between them; here the cuvettes are matched one for one, so a per-
# CUVETTE offset absorbs the composition exactly as a per-experiment offset
# absorbs the day. log[HOO-] = log[H2O2] + f(pH), and with the cuvette's own
# [H2O2] held in its offset what is left driving the regressor is f(pH) alone.
#
# The enzyme is what forces the split. It is not constant over the block --
# 0.069, 0.034, 0.014 and 0.021 mM -- and it steps between runs, exactly where
# pH does, so a per-cuvette offset cannot touch it. Exps 141 and 142 sit at
# 0.014 mM against exps 136-140's 0.034 on the same composition and at the top
# of that group's pH range, so pooling them would push a rate that scales with
# enzyme downward precisely at high pH and bias the slope down with it. A
# ladder is therefore runs sharing a composition AND a loading.
PH_LADDER_MINIMUM = 3
PH_AXIS = "hoo"


def _composition(group):
    """A run's cuvette set, as a hashable signature. Rounded, not exact."""
    return tuple(sorted((round(float(s), 6), round(float(h), 6))
                        for s, h in zip(group.s0, group.h2o2)))


def ph_ladders(scope=TWO_AXIS_BLOCK, minimum=PH_LADDER_MINIMUM):
    """
    Groups of runs that share a composition and an enzyme loading, over pH.

    Returns {label: DataFrame}, keyed by the loading, with `minimum` distinct
    pH values required before a group counts as a ladder. Derived from the
    frame -- no experiment number appears here -- so a run that stops matching
    stops being in a ladder.
    """
    data = frame(scope).copy()
    signatures = {experiment: _composition(group)
                  for experiment, group in data.groupby("experiment")}
    data["design"] = [signatures[e] for e in data.experiment]
    found = {}
    for (loading, _), group in data.groupby([data.e0.round(6), "design"],
                                            sort=False):
        if group.pH.nunique() >= minimum:
            found[f"{loading:g} mM chemzyme"] = group.drop(columns="design")
    return dict(sorted(found.items(),
                       key=lambda item: -item[1].experiment.nunique()))


def ph_order(parameter="vmax", scope=TWO_AXIS_BLOCK, axis=PH_AXIS,
             minimum=PH_LADDER_MINIMUM, ladders=None):
    """
    d log(`parameter`) / d log[HOO-] at fixed composition, one offset per cuvette.

    One row per ladder plus a pooled row carrying one offset per (ladder,
    cuvette). The offsets are what make this a within-design measurement: the
    cuvette is matched across every run of a ladder, so its substrate, its
    peroxide and its position in the cell holder all sit in its own intercept
    and the slope is read from pH alone.
    """
    groups = ph_ladders(scope, minimum) if ladders is None else ladders
    tagged = {}
    for label, group in groups.items():
        group = group.copy()
        # The offset key. Its cuvette, within its ladder: a cuvette of one
        # ladder is not the same cuvette as its namesake in the other, which
        # runs a different peroxide and a different enzyme.
        group["cuvette"] = [f"{label}/{s}" for s in group["sample"]]
        tagged[label] = group
    rows = []
    pooled = ([("pooled", pd.concat(tagged.values()))]
              if len(tagged) > 1 else [])
    for label, group in list(tagged.items()) + pooled:
        live = group[group.live & (group[parameter] > 0)
                     & np.isfinite(group[parameter]) & (group[axis] > 0)]
        keys = list(live.cuvette)
        if live.pH.nunique() < minimum or len(live) < minimum + 2:
            continue
        y = np.log(live[parameter].to_numpy(dtype=float))
        x = np.log(live[axis].to_numpy(dtype=float))
        cuvettes = sorted(set(keys))
        design_matrix = np.column_stack(
            [x] + [(np.array(keys) == key).astype(float) for key in cuvettes])
        beta, *_ = np.linalg.lstsq(design_matrix, y, rcond=None)
        residual = y - design_matrix @ beta
        rank = int(np.linalg.matrix_rank(design_matrix))
        variance = float(residual @ residual) / max(1, len(y) - rank)
        covariance = variance * np.linalg.pinv(design_matrix.T @ design_matrix)
        total = float(((y - y.mean()) ** 2).sum())
        rows.append({"ladder": label,
                     "runs": int(live.experiment.nunique()),
                     "curves": int(len(live)),
                     "cuvettes": len(cuvettes),
                     "pH_low": float(live.pH.min()),
                     "pH_high": float(live.pH.max()),
                     "order": float(beta[0]),
                     "stderr": float(np.sqrt(max(covariance[0, 0], 0.0))),
                     "r2": float(1 - (residual @ residual) / total)
                     if total > 0 else np.nan})
    return pd.DataFrame(rows).set_index("ladder")


ACCELERATION_SPLIT = 9.0


def acceleration_by_ph(scope=TWO_AXIS_BLOCK, split=ACCELERATION_SPLIT):
    """
    The autocatalytic acceleration, banded either side of a pH.

    The block's long runs are the ones that DECELERATE -- exps 138, 146 and
    148-151 run eight hours because they are slow -- so "the acceleration
    builds up over a long run" reads the design backwards. Banding by pH is the
    statement that survives: what the acceleration tracks is [HOO-], and a run
    at pH 9 shows it in an hour.

    Live curves only. A dead curve's `accelerates` is its quantisation
    staircase stepping late, which is why `design` bands the same way.
    """
    data = frame(scope)
    live = data[data.live]
    bands = np.where(live.pH >= split, f"pH >= {split:g}", f"pH < {split:g}")
    table = (live.assign(band=bands).groupby("band")
             .agg(curves=("live", "size"),
                  accelerating=("accelerates", "sum"),
                  median_late_over_early=("late_over_early", "median"))
             .sort_index(ascending=False))
    table["share"] = table.accelerating / table.curves
    return table


def run_dates(scope=TWO_AXIS_BLOCK):
    """
    When each run was collected, from the instrument's own export header.

    An INDEPENDENT source, in the sense verify_instrument.py means it: the
    workbook records no date at all, so this is the only account of when a run
    happened and it is written by the spectrophotometer rather than by hand.
    `kinetics_io.parse_experiment_data` already reads the field; nothing else
    in the pipeline had asked for it.

    Returns a frame indexed by experiment with `date` and `order`, the rank of
    that date, so a caller can ask what tracks the schedule without assuming
    experiment number is chronological. Over exps 135-151 it happens to be:
    3 to 14 September 2010, strictly increasing.
    """
    import datetime

    from kinetics_io import parse_experiment_data
    from verify_instrument import export_path

    rows = []
    for experiment in sorted(int(e) for e in scope):
        path = export_path(experiment)
        parsed = parse_experiment_data(path) if path else None
        text = (parsed or {}).get("date")
        stamp = None
        if text:
            for pattern in ("%m/%d/%Y", "%d/%m/%Y", "%Y-%m-%d"):
                try:
                    stamp = datetime.datetime.strptime(text.strip(), pattern)
                    break
                except ValueError:
                    continue
        rows.append({"experiment": experiment, "date": stamp})
    table = pd.DataFrame(rows).set_index("experiment")
    table["order"] = table.date.rank(method="dense")
    return table


def ph_schedule_control(parameter="vmax", scope=None):
    """
    Could the pH order be the stock ageing instead? The two ladders answer it.

    The pH ladders were not run in the same direction. Exps 136-140 climb
    6.95 -> 9.22 over 3-4 September; exps 143-151 descend 9.73 -> 5.47 over
    5-14 September. So pH correlates with the schedule POSITIVELY in one ladder
    and NEGATIVELY in the other, and anything that drifts monotonically with
    the schedule -- an enzyme stock losing activity over twelve days, a lamp,
    a stock of peroxide decomposing -- enters the two ladders' pH slopes with
    OPPOSITE signs.

    That makes the control free and it needs no assumption about what drifts or
    how fast. If both ladders return the same sign, no monotone schedule effect
    produced it; the sign the two agree on is the sign of pH.

    Returns the per-ladder slopes, each ladder's pH-against-schedule
    correlation, and whether the signs agree.
    """
    scope = strong_runs() if scope is None else scope
    schedule = run_dates(scope)
    table = ph_order(parameter, scope=scope)
    rows = []
    for label, group in ph_ladders(scope).items():
        if label not in table.index:
            continue
        runs = sorted(group.experiment.unique())
        pH = [float(group.loc[group.experiment == e, "pH"].iloc[0])
              for e in runs]
        order = [float(schedule.loc[e, "order"]) for e in runs]
        rows.append({"ladder": label,
                     "runs": len(runs),
                     "first": schedule.loc[runs[0], "date"],
                     "last": schedule.loc[runs[-1], "date"],
                     "pH_vs_schedule": float(np.corrcoef(pH, order)[0, 1]),
                     "order": float(table.loc[label, "order"]),
                     "stderr": float(table.loc[label, "stderr"])})
    out = pd.DataFrame(rows).set_index("ladder")
    signs = set(np.sign(out.order))
    directions = set(np.sign(out.pH_vs_schedule))
    return out, {"ladders": len(out),
                 "opposed_schedules": len(directions) > 1,
                 "orders_agree_in_sign": len(signs) == 1}


# ---------------------------------------------------------------------------
# The size of the early rise, and whether it is one turnover or many.
#
# Exp 150 cuvette 6 is the curve that prompted this: a rise to 0.0120 AU by
# 7900 s and then a steady fall at -2.6e-7 AU/s over the remaining six hours,
# on 0.057% of its substrate. Divided by the extinction coefficient and the
# catalyst it is HALF AN EQUIVALENT -- one turnover and stop, which is what a
# catalyst whose regeneration has failed looks like, and at pH 6.26 the step
# that regenerates it needs a [HOO-] four decades down.
#
# WHAT IS AND IS NOT COMPARABLE. `burst` is the fitted curve's rise above its
# own start, up to the fit's maximum. When that maximum sits at the end of the
# run the curve had not stopped rising and the value is a LOWER BOUND, which
# `burst_bounded` records -- and only 17 of the block's 110 live curves are
# bounded, so most of the block never finishes its rise.
#
# The bound bites BETWEEN runs and not within them. Every cuvette of one run
# shares its length, so a per-experiment offset absorbs the truncation exactly
# and an order in [S] or [H2O2] may be read off all the live curves;
# `burst_drivers` does that. Comparing one run's burst with another's may not,
# and the block spans 3000 to 28740 s, so `enzyme_pair` -- which has to
# compare runs, because [enz] never moves inside one -- reads both members at a
# window in SECONDS common to both instead.
# The window to read both runs at, as a fraction of the SHORTER run's length.
# One, not a half: the largest window both runs cover is the one choice here
# that is not arbitrary, it uses the most signal, and it needs no extrapolation
# of either fit beyond its own data. `enzyme_pair_sensitivity` sweeps it,
# because the number does move -- the apparent order runs 0.77 to 1.35 over
# 945-3780 s, a systematic the size of the statistical error, and every window
# in the sweep is consistent with first order and excludes no dependence.
ENZYME_PAIR_FLOOR = 1.0


def burst_table(scope=TWO_AXIS_BLOCK, live_only=True):
    """One row per curve: the early rise, when it peaked, and in catalyst units."""
    data = frame(scope)
    if live_only:
        data = data[data.live]
    columns = ["experiment", "sample", "pH", "e0", "s0", "h2o2", "hoo",
               "duration_s", "progress_kind", "burst", "burst_time_s",
               "burst_bounded", "turnovers"]
    return data[columns].sort_values(["pH", "experiment", "sample"])


def burst_drivers(scope=TWO_AXIS_BLOCK):
    """
    What the early rise scales with -- within runs, where the axes move.

    Returns the orders of `burst` in [S] and [H2O2] against one offset per
    experiment, and `enzyme_identified` to say what it cannot answer: [enz] is
    constant across every cuvette of a run and moves only between runs, so the
    offsets absorb it and `_moves` refuses it. `enzyme_pair` is the route to
    that one.
    """
    data = frame(scope)
    live = data[data.live & np.isfinite(data.burst) & (data.burst > 0)]
    result = orders("burst", frame=live)
    within = [float(np.ptp(block.e0)) > 0
              for _, block in live.groupby("experiment")]
    return {**result,
            "curves": int(len(live)),
            "bounded": int(live.burst_bounded.sum()),
            "enzyme_identified": any(within)}


def enzyme_pair(scope=TWO_AXIS_BLOCK, floor=ENZYME_PAIR_FLOOR):
    """
    Does the early rise scale with the CATALYST? The block's one lever on it.

    [enz] never moves inside a run, so this has to compare runs, and only one
    comparison in the block is worth making: two runs that share a composition,
    differ in loading, and sit as close in pH as the block allows. It is
    derived here rather than listed -- the pair is whichever qualifying pair
    has the smallest pH gap -- and comes out as exps 140 and 141, 0.034 against
    0.014 mM, pH 9.22 against 9.15.

    READ AT A WINDOW IN SECONDS COMMON TO BOTH, not at each run's own maximum,
    and at the largest such window rather than a fraction chosen to taste --
    `ENZYME_PAIR_FLOOR`, swept by `enzyme_pair_sensitivity`.
    The runs are 3780 s and 4680 s and `burst` stops at whichever time each
    fit peaks, so comparing those two numbers would be comparing two windows --
    the error `induction.buffer_lever` exists to avoid. Reading both fitted
    curves at the same absolute time gives all seven cuvette pairs instead of
    the one that happens to be bounded in both.

    AND CORRECTED FOR THE PH GAP THAT IS LEFT, using this block's own measured
    pH order (`ph_order`), because 0.07 pH units at pH 9.2 is about a 9%
    effect and the enzyme step being tested is only 2.4x.

    Returns (per-cuvette table, verdict). The verdict's `ratio` is the observed
    rise ratio after correction and `expected` the loading ratio: equal means
    the rise is proportional to catalyst.
    """
    from summary_kinetics import fit_progress

    data = frame(scope)
    signatures = {experiment: _composition(group)
                  for experiment, group in data.groupby("experiment")}
    candidates = []
    for high, low in itertools.permutations(sorted(signatures), 2):
        if signatures[high] != signatures[low]:
            continue
        rows = {e: data[data.experiment == e] for e in (high, low)}
        if rows[high].e0.iloc[0] <= rows[low].e0.iloc[0]:
            continue
        candidates.append((abs(rows[high].pH.iloc[0] - rows[low].pH.iloc[0]),
                           high, low))
    if not candidates:
        return pd.DataFrame(), {"pairs": 0}
    gap, high, low = min(candidates)

    lookup = {(c.experiment, c.sample): c for c in curves(scope)}
    window = floor * min(data[data.experiment == e].duration_s.max()
                         for e in (high, low))
    rises = {}
    for experiment in (high, low):
        for row in data[data.experiment == experiment].itertuples():
            curve = lookup.get((experiment, row.sample))
            if curve is None:
                continue
            times = np.asarray(curve.times, dtype=float)
            values = np.asarray(curve.absorbance, dtype=float)
            fitted = fit_progress(times, values)
            grid = np.array([times[0], window], dtype=float)
            predicted = fitted.predict(grid)
            rises[(experiment, row.sample)] = float(predicted[1]
                                                    - predicted[0])

    order = ph_order("vmax", scope=scope)
    slope = float(order.loc["pooled", "order"] if "pooled" in order.index
                  else order.iloc[0].order)
    high_rows = data[data.experiment == high]
    low_rows = data[data.experiment == low]
    # [HOO-] at the two pH values, at matched composition: its ratio is the
    # correction the pH order is applied to.
    hoo_ratio = float(high_rows.hoo.iloc[0] / low_rows.hoo.iloc[0])
    correction = hoo_ratio ** slope

    rows = []
    for sample in sorted(set(high_rows["sample"]) & set(low_rows["sample"])):
        top, bottom = rises.get((high, sample)), rises.get((low, sample))
        if not top or not bottom or bottom <= 0 or top <= 0:
            continue
        rows.append({"sample": int(sample),
                     "s0": float(high_rows.loc[high_rows["sample"] == sample,
                                               "s0"].iloc[0]),
                     "h2o2": float(high_rows.loc[high_rows["sample"] == sample,
                                                 "h2o2"].iloc[0]),
                     f"rise_{high}": top, f"rise_{low}": bottom,
                     "ratio": top / bottom,
                     "corrected": (top / bottom) / correction})
    table = pd.DataFrame(rows).set_index("sample")
    expected = float(high_rows.e0.iloc[0] / low_rows.e0.iloc[0])
    observed = float(np.exp(np.log(table.corrected.to_numpy()).mean())) \
        if len(table) else np.nan
    spread = float(np.exp(np.log(table.corrected.to_numpy()).std(ddof=1)
                          / np.sqrt(len(table)))) if len(table) > 1 else np.nan
    # In logs, because the quantity is a ratio: the spread is a FACTOR and a
    # symmetric error bar on a ratio is the wrong shape.
    log_error = float(np.log(spread)) if np.isfinite(spread) else np.nan
    return table, {"high": int(high), "low": int(low),
                   "pH_gap": float(gap), "window_s": float(window),
                   "cuvettes": int(len(table)),
                   "pH_correction": correction,
                   "expected": expected, "ratio": observed,
                   "stderr_factor": spread,
                   "order": float(np.log(observed) / np.log(expected))
                   if observed > 0 and expected > 1 else np.nan,
                   # Against the two hypotheses worth naming: the rise is
                   # proportional to catalyst, and the rise does not depend on
                   # catalyst at all.
                   "sigma_first_order": abs(np.log(observed)
                                            - np.log(expected)) / log_error
                   if log_error > 0 else np.nan,
                   "sigma_no_dependence": abs(np.log(observed)) / log_error
                   if log_error > 0 else np.nan}


ENZYME_PAIR_SWEEP = (0.25, 0.5, 0.75, 1.0)


def enzyme_pair_sensitivity(scope=TWO_AXIS_BLOCK, floors=ENZYME_PAIR_SWEEP):
    """
    The pair's ratio at four windows, because the window was a choice.

    `enzyme_pair` reads both runs at one absolute time and that time is not
    handed down by anything -- so the number is only worth quoting if it does
    not move much when the window does. This is the same discipline as
    `slowdown.sink_window_sensitivity` and `induction.BUFFER_WINDOW_SWEEP`.
    """
    rows = []
    for floor in floors:
        _, verdict = enzyme_pair(scope, floor=floor)
        rows.append({"floor": floor, "window_s": verdict.get("window_s"),
                     "ratio": verdict.get("ratio"),
                     "stderr_factor": verdict.get("stderr_factor"),
                     "order": verdict.get("order"),
                     "cuvettes": verdict.get("cuvettes")})
    return pd.DataFrame(rows).set_index("floor")


# Above this the detachments carry more absorbance than the reaction does, and
# no repair in curve_metrics recovers the rate: against planted sawtooths the
# subtraction is unbiased to a load of about 0.5 (1.01), leaves a tenth by 1
# (1.13) and half by 2 (1.52). It is a CEILING ON USE, not an exclusion --
# every curve stays in the frame, on the page and in the counts.
BUBBLE_LOAD_CEILING = 1.0

# The peroxide bands the drop rate is reported over. Edges, in mM.
BUBBLE_PEROXIDE_BANDS = (0.0, 5.0, 10.0, 25.0, 40.0, 80.0, 200.0)


def bubble_table(scope=TWO_AXIS_BLOCK, live_only=True):
    """
    Every curve's gas load, and what its rate looks like under each repair.

    Columns: experiment, sample, s0, h2o2, e0, bubble_drops, bubble_events,
    bubble_load, gas_rate, vmax, vmax_corrected, vmax_monotone, repairable.

    `bubble_drops` counts falling READINGS and `bubble_events` the bubbles
    that made them -- a fall may span two readings, and 50 detaching curves
    carry 220 falls in 206 events.

    `repairable` is `bubble_load <= BUBBLE_LOAD_CEILING`. Read it before
    quoting a rate off any curve in this block; it is the flag the folder
    documents cite and the one `check_numbers` counts.
    """
    data = frame(scope)
    if live_only:
        data = data[data.live]
    out = data[["experiment", "sample", "s0", "h2o2", "e0", "bubble_drops",
                "bubble_events", "bubble_load", "gas_rate", "vmax",
                "vmax_corrected", "vmax_monotone"]].copy()
    out["repairable"] = out.bubble_load <= BUBBLE_LOAD_CEILING
    return out.sort_values("bubble_load", ascending=False).reset_index(
        drop=True)


def bubble_ladder(scope=TWO_AXIS_BLOCK, bands=BUBBLE_PEROXIDE_BANDS):
    """
    Detachments against [H2O2]: the evidence that the chop is gas.

    One row per peroxide band, with the number of live curves, how many carry
    a detachment, and the mean absorbance lost to them. The count rises
    monotonically across the block's six bands, from none below 5 mM to every
    curve above 80 mM.

    PEROXIDE ALONE IS NOT ENOUGH -- see `bubble_turnover_control`, which is the
    other half of this argument and the reason the cause is the CATALYSED
    decomposition rather than the peroxide sitting in a cuvette.
    """
    data = frame(scope)
    data = data[data.live].copy()
    data["band"] = pd.cut(data.h2o2, list(bands))
    grouped = data.groupby("band", observed=True)
    return pd.DataFrame({
        "curves": grouped.size(),
        "with_drops": grouped.bubble_drops.apply(lambda s: int((s > 0).sum())),
        "mean_drops": grouped.bubble_drops.mean(),
        "mean_lost": grouped.apply(
            lambda g: float((g.bubble_load * g.net).sum() / len(g)),
            include_groups=False),
    })


def bubble_turnover_control(scope=TWO_AXIS_BLOCK):
    """
    Runs at the block's top peroxide that carry no detachment at all.

    The control on `bubble_ladder`. If gas came from peroxide standing in a
    cuvette, every run at 73.4 mM would chop; exps 136 and 137 sit there and
    carry none, and they are the block's two weakest runs on
    `concentration_agreement` (0.25 and 0.21) with median rates thirty times
    below exp 140's. So the gas needs TURNOVER as well as peroxide, which is
    what makes it the ketone-catalysed decomposition rather than the peroxide.

    One row per run, sorted by the peroxide it reaches: top_h2o2, drops,
    agreement, median_rate. Read the two runs that sit at the block's own top
    peroxide with a drop count of zero.
    """
    data = frame(scope)
    data = data[data.live]
    if data.empty:
        return pd.DataFrame()
    agreement = concentration_agreement(scope)
    rows = []
    for experiment, group in data.groupby("experiment"):
        rows.append({
            "experiment": int(experiment),
            "top_h2o2": float(group.h2o2.max()),
            "cuvettes": len(group),
            "drops": int(group.bubble_drops.sum()),
            "agreement": float(agreement.agreement.get(experiment, np.nan)),
            "median_rate": float(agreement.median_rate.get(experiment,
                                                           np.nan)),
        })
    return pd.DataFrame(rows).set_index("experiment").sort_values(
        "top_h2o2", ascending=False)


def bubble_synchrony(scope=TWO_AXIS_BLOCK):
    """
    Do the cuvettes of one run detach together? They do not.

    The instrument control. A lamp flicker, a shutter or the carousel would hit
    all seven cuvettes of a run at the same reading; a bubble on one window
    cannot. Counts detachment times shared by a pair of cuvettes in the same
    run against what independence predicts -- for a pair with `a` and `b`
    detachments over `n` intervals, `a * b / n`.

    Returns a dict: `observed`, `expected`, `pairs`, `runs`. Over the two-axis
    block the two agree, which leaves the light path itself as the only place
    the absorbance can be going.
    """
    observed = expected = 0.0
    pairs = runs = 0
    for experiment in sorted({c.experiment for c in curves(scope)}):
        group = [c for c in curves_of(experiment)]
        if len(group) < 2:
            continue
        runs += 1
        marks = []
        for curve in group:
            values = np.asarray(curve.absorbance, dtype=float)
            times = np.asarray(curve.times, dtype=float)
            drops = bubble_drops(values, curve.noise)
            marks.append((set(np.round(times[drops], 3)), len(values) - 1))
        for first in range(len(marks)):
            for second in range(first + 1, len(marks)):
                (one, count_one), (two, count_two) = marks[first], marks[second]
                intervals = min(count_one, count_two)
                if intervals <= 0:
                    continue
                pairs += 1
                observed += len(one & two)
                expected += len(one) * len(two) / intervals
    return {"observed": int(observed), "expected": float(expected),
            "pairs": int(pairs), "runs": int(runs)}


def bubble_mass_balance(scope=TWO_AXIS_BLOCK):
    """
    What each repair claims the substrate produced, as a fraction of what it
    could.

    THE TEST THAT REFUTES STITCHING. `stitched` is the net a curve would show
    if every detachment were added back to the readings after it. On exp 135
    sample 4 that is 1.26 -- more absorbance than 0.219 mM of substrate can
    make -- and four more curves land between 0.49 and 0.86 while their ramps
    STEEPEN (`ramp_gain` below 1 is a reaction spending its substrate; these
    read 1.07 to 23.5). Subtracting the ramp instead keeps every curve under 1.

    Columns: experiment, sample, s0, h2o2, ceiling, raw, stitched, corrected,
    ramp_gain. `ceiling` is `epsilon * s0` in AU.
    """
    rows = []
    for curve in curves(scope):
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        ceiling = float(curve.epsilon * curve.conditions.s0)
        if ceiling <= 0:
            continue
        steps = np.diff(values)
        drops = bubble_drops(values, curve.noise)
        corrected, _ = debubble(times, values, curve.noise)
        # The ramps' own slope, early against late, over rising steps only, so
        # that the detachments do not enter the comparison.
        quarter = max(1, len(steps) // 4)
        gains = []
        for span in (slice(0, quarter), slice(-quarter, None)):
            rising = steps[span] > 0
            gains.append(steps[span][rising].sum() / np.diff(times)[span][
                rising].sum() if rising.any() else np.nan)
        rows.append({
            "experiment": curve.experiment, "sample": curve.sample,
            "s0": curve.conditions.s0, "h2o2": curve.conditions.h2o2,
            "ceiling": ceiling,
            "raw": float(values[-1] - values[0]) / ceiling,
            "stitched": float(values[-1] - values[0]
                              - steps[drops].sum()) / ceiling,
            "corrected": float(corrected[-1] - corrected[0]) / ceiling,
            "ramp_gain": float(gains[1] / gains[0]) if gains[0] else np.nan,
        })
    return pd.DataFrame(rows)


# A step this far below zero, in units of the curve's own noise, is a fall the
# noise cannot explain. It is the detector's own threshold, reused: a
# reconstruction that still carries one has not removed the artefact it was
# built for.
# The severities the recovery table is reported at: the artefact's total rise
# as a multiple of the chemistry's. 2.0 is past anything in the block and is
# there because that is where the segment ramp broke.
RECOVERY_SEVERITIES = (0.25, 0.5, 1.0, 2.0)

# The seed the planting draws its detachment times from. Fixed so the table in
# the document and the table the test asserts are the same table.
RECOVERY_SEED = 7

# How many plantings each donor curve carries at each severity.
RECOVERY_REPEATS = 4


def bubble_recovery(severities=RECOVERY_SEVERITIES, emptying=True,
                    scope=TWO_AXIS_BLOCK, seed=RECOVERY_SEED,
                    repeats=RECOVERY_REPEATS):
    """
    Recovered `vmax` over true `vmax`, for each repair, at a known truth.

    The donors are the block's OWN clean curves, so the chemistry being
    recovered is real chemistry with real noise; only the artefact is planted,
    and the rate measured before planting is the truth. 1.00 is exact and
    above 1.00 is a rate that the gas inflated.

    `emptying` chooses the planting, and the choice is the argument. With
    `True` each bubble leaves completely at its detachment, which is what
    `bubble_ramp` assumed until 2026-09-03. With `False` each leaves only 40 to
    100% of what it holds and the rest carries over, which is what
    `bubble_record(141, 3)` shows -- a curve shedding its LARGEST drop after
    its SHORTEST growth window cannot be one bubble emptying each time. A
    repair has to survive both, and the segment ramp survived only the first:
    1.01, 1.01, 1.12, 1.60 emptying, against 1.05, 1.10, 1.28, 1.70 not. It
    is also what `curve_metrics.bubble_ceiling` is measured rather than
    assumed for: capping the stretches between detachments as well as the
    tail costs 1.11 and 1.34 here at 1x and 2x, while changing nothing
    under the emptying planting.

    Columns: severity, raw, stitched, rebuilt, n.
    """
    generator = np.random.default_rng(seed)
    donors = [c for c in curves(scope)
              if not len(bubble_drops(np.asarray(c.absorbance, dtype=float),
                                      c.noise))
              and float(np.ptp(np.asarray(c.absorbance, dtype=float))) > 0.02]
    rows = []
    for severity in severities:
        raw, stitched, rebuilt = [], [], []
        for curve in donors:
            times = np.asarray(curve.times, dtype=float)
            values = np.asarray(curve.absorbance, dtype=float)
            floor = source_floor(curve.source)
            truth = peak_rate(times, values, floor=floor)[0]
            if not np.isfinite(truth) or truth <= 0:
                continue
            for _ in range(repeats):
                edges = np.sort(generator.choice(
                    np.arange(5, len(times) - 5), size=5, replace=False))
                rate = severity * (values[-1] - values[0]) / times[-1]
                artefact = _planted_gas(times, edges, rate, emptying,
                                        generator)
                spoilt = values + artefact
                drops = bubble_drops(spoilt, curve.noise)
                joined = spoilt.copy()
                for index in drops:
                    joined[index + 1:] -= np.diff(spoilt)[index]
                fixed, _ = debubble(times, spoilt, curve.noise)
                for bag, series in ((raw, spoilt), (stitched, joined),
                                    (rebuilt, fixed)):
                    bag.append(peak_rate(times, series, floor=floor)[0] / truth)
        rows.append({"severity": severity,
                     "raw": float(np.median(raw)),
                     "stitched": float(np.median(stitched)),
                     "rebuilt": float(np.median(rebuilt)),
                     "n": len(raw)})
    return pd.DataFrame(rows).set_index("severity")


def _planted_gas(times, edges, rate, emptying, generator):
    """One planted artefact: gas made at `rate`, shed at each of `edges`."""
    if emptying:
        artefact = np.zeros(len(times))
        start = 0
        for edge in edges:
            artefact[start:edge + 1] = rate * (
                times[start:edge + 1] - times[start])
            start = edge + 1
        artefact[start:] = rate * (times[start:] - times[start])
        return artefact
    held = rate * (times - times[0])
    released = np.zeros(len(times))
    for edge in edges:
        released[edge + 1:] += (held[edge] - released[edge]) * (
            generator.uniform(0.4, 1.0))
    return held - released


REBUILD_STEP_SIGMA = BUBBLE_DROP_SIGMA


def rebuild_smoothness(scope=TWO_AXIS_BLOCK, sigma=REBUILD_STEP_SIGMA,
                       live_only=True):
    """
    Does a repaired curve look like a curve that never bubbled?

    THE TEST THE REPAIR HAS TO PASS, and the one the segment ramp failed. A
    reconstruction is only a reconstruction if what comes out carries no fall
    the noise cannot explain -- the curves that never bubbled set the standard,
    and the repaired ones have to meet it rather than merely improve on where
    they started.

    One row per curve: experiment, sample, bubble_events, bubble_load,
    raw_worst (the steepest single fall in the READINGS, in units of the
    curve's own noise), rebuilt_worst (the same for the reconstruction),
    `clean` for the curves with no detachment at all, and three columns for
    the second way a repair can go wrong: gas_held (the most the profile ever
    puts in the beam), biggest_bubble (the largest single detachment), and
    rebuilt_net.

    READ `gas_held` AGAINST `biggest_bubble`. A monotone reconstruction can
    still be a wrong one -- until 2026-09-03 this table would have shown 12 of
    49 curves holding more than twice their own largest bubble, one at 26.6x,
    and three finishing BELOW ZERO, because the production rate was
    extrapolated across hours in which nothing detached. Every one of those
    curves passed the smoothness test above, which is why this is measured
    separately: the fault was not roughness, it was level.

    Over the two-axis block the repaired curves' worst step goes from -260.4
    sigma to -9.6 sigma, and the only one still past `sigma` is exp 135
    cuvette 6 -- the first-interval detachment `debubble` returns untouched.
    The curves that never bubbled reach -5.8, so one curve separates the two
    populations and it is the one the model declines to touch.
    """
    keep = set()
    if live_only:
        data = frame(scope)
        keep = set(zip(data[data.live].experiment, data[data.live]["sample"]))
    rows = []
    for curve in curves(scope):
        if live_only and (curve.experiment, curve.sample) not in keep:
            continue
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        rebuilt, events = debubble(times, values, curve.noise)
        steps = np.diff(values)
        rows.append({
            "experiment": curve.experiment,
            "sample": curve.sample,
            "bubble_events": len(events),
            "bubble_load": bubble_load(values, bubble_drops(
                values, curve.noise)),
            "raw_worst": float(steps.min() / curve.noise),
            "rebuilt_worst": float(np.diff(rebuilt).min() / curve.noise),
            "clean": not len(events),
            "gas_held": float((values - rebuilt).max()),
            "biggest_bubble": max(
                (float(values[start] - values[stop])
                 for start, stop in events), default=np.nan),
            "rebuilt_net": float(rebuilt[-1] - rebuilt[0]),
        })
    return pd.DataFrame(rows)


def gas_rate_drivers(scope=TWO_AXIS_BLOCK):
    """
    What sets the rate `bubble_rate` fits -- and the check that it is O2.

    THE FIT NEVER SEES A CONCENTRATION. `bubble_rate` is read off the timing
    and size of the detachments alone, so its dependence on the composition is
    a prediction the gas argument makes and this measures. If the gas is the
    ketone-catalysed decomposition of the peroxide, the rate belongs to the
    PEROXIDE and the catalyst, not to the alcohol.

    IT IS FIRST ORDER IN PEROXIDE: +1.203 +/- 0.221 over the 49 curves that
    carry a rate, which is the catalysed decomposition of H2O2 measured
    without ever fitting a concentration to it. That is the strongest
    independent support the gas argument has -- `bubble_ladder` shows the
    drops get MORE COMMON with peroxide, and this shows the gas is made
    FASTER, in proportion, which is a different measurement of the same claim.

    The substrate carries -0.250 +/- 0.094 with one offset per experiment.
    Read it as a weak negative and not as a null: it is 2.6 sigma, it is what
    a substrate competing with the decomposition for the same catalyst would
    give, and it is measured on a derived quantity over the high-peroxide
    curves alone. What matters for the diagnosis is that it is not the
    substrate order of a rate -- the alcohol is not what is turning into gas.

    The peroxide is identified only BETWEEN runs here, so it is reported
    pooled and carries the usual between-run caveat; the substrate moves
    within runs and is read the way every other order in this module is.

    Returns the `orders` dictionary, plus `n_pooled` and the pooled peroxide
    coefficient for the axis the offsets cannot see.
    """
    data = frame(scope)
    data = data[data.live & (data.bubble_events > 0)
                & np.isfinite(data.gas_rate) & (data.gas_rate > 0)]
    within = orders("gas_rate", scope=scope, frame=data)
    pooled = orders("gas_rate", scope=scope, frame=data, within=False)
    return {"order_s0": within["order_s0"], "stderr_s0": within["stderr_s0"],
            "n": within["n"],
            "pooled_s0": pooled["order_s0"],
            "pooled_stderr_s0": pooled["stderr_s0"],
            "pooled_h2o2": pooled["order_h2o2"],
            "pooled_stderr_h2o2": pooled["stderr_h2o2"],
            "n_pooled": pooled["n"]}


BUBBLE_TREATMENTS = ("vmax", "vmax_corrected", "vmax_monotone")


def bubble_sensitivity(scope=TWO_AXIS_BLOCK, treatments=BUBBLE_TREATMENTS,
                       ceiling=BUBBLE_LOAD_CEILING):
    """
    Every concentration order this block reports, under each repair.

    THE ANSWER TO "IS vmax MISLEADING". Curve by curve, on the worst curves,
    badly -- the bound cuts `vmax` to a sixth on one of them. Block by block,
    no: over the strong runs every order here moves by less than its own
    standard error, and the substrate order moves AWAY from zero rather than
    toward it, which is the direction that matters because a substrate-blind,
    peroxide-driven artefact is exactly what would manufacture a flat one.

    Rows: each treatment, plus `vmax` over the repairable curves alone.
    """
    data = frame(scope)
    rows = []
    for treatment in treatments:
        rows.append({"treatment": treatment, "curves": "all live",
                     **orders(treatment, frame=data)})
    kept = data[data.bubble_load.fillna(0.0) <= ceiling]
    rows.append({"treatment": "vmax", "curves": f"load <= {ceiling:g}",
                 **orders("vmax", frame=kept)})
    return pd.DataFrame(rows).set_index(["treatment", "curves"])

# A step this many of a curve's own noise is "large" for the purpose of asking
# whether large steps rise or fall. Well above the detachment threshold, so the
# comparison is between the two TAILS and not between a tail and the noise.
BUBBLE_STEP_SIGMA = 20.0

# The curve section 5 walks through detachment by detachment. Exp 141 cuvette 3
# because it is the clearest case of the thing the repair cannot do: four
# detachments, the largest after the shortest growth window, and a level after
# the last one below every earlier level -- so more than one bubble was growing
# at a time. Named here so the document and its checker read the same curve.
BUBBLE_WORKED_EXAMPLE = (141, 3)


def bubble_step_asymmetry(scope=TWO_AXIS_BLOCK, sigma=BUBBLE_STEP_SIGMA):
    """
    Do the block's large steps rise or fall? They fall, about five to one.

    WHICH BEAM THE GAS IS IN, read off the shape alone. A bubble growing in the
    SAMPLE beam scatters light out of the aperture, so absorbance climbs while
    it grows and drops when it detaches: slow up, sudden down. A bubble in the
    REFERENCE beam does the same to the reference, and the difference inverts
    it: slow down, sudden up. The two populations are therefore separable by
    the sign of the large steps, and they are not balanced.

    Returns a dict: `steps`, `rises`, `falls` (counts beyond `sigma`),
    `largest_rise`, `largest_fall` (in noise), `ratio` = falls / rises.

    The reference cuvette omits only the ENZYME, so it holds the same peroxide
    as the sample and would bubble too if the gas came from peroxide standing
    in solution. It barely does, which is the same conclusion
    `bubble_turnover_control` reaches from the other side.
    """
    rises = falls = steps = 0
    largest_rise = largest_fall = 0.0
    for curve in curves(scope):
        values = np.asarray(curve.absorbance, dtype=float)
        if len(values) < 2 or not np.isfinite(curve.noise) or curve.noise <= 0:
            continue
        z = np.diff(values) / curve.noise
        steps += len(z)
        rises += int((z > sigma).sum())
        falls += int((z < -sigma).sum())
        largest_rise = max(largest_rise, float(z.max()))
        largest_fall = min(largest_fall, float(z.min()))
    return {"steps": steps, "rises": rises, "falls": falls,
            "largest_rise": largest_rise, "largest_fall": largest_fall,
            "ratio": (falls / rises) if rises else np.nan}


def bubble_record(experiment, sample, sigma=BUBBLE_DROP_SIGMA):
    """
    One curve's detachments, one row each: when, how much, and what preceded it.

    Columns: time_s, lost, sigma, grew_s (seconds since the previous
    detachment, or since t = 0), rose (how far the curve climbed over that
    window), ramp (that climb per second), after (the reading immediately
    following the detachment).

    `lost` and not `drop`, because `DataFrame.drop` is a method and
    `record.drop.max()` silently reaches the method rather than the column.

    READ `grew_s` AGAINST `lost`. One bubble growing steadily on one site would
    shed in proportion to the time it had to grow. Exp 141 cuvette 3 does not:
    600 s of growth sheds 0.0215, another 600 s sheds 0.0153, and 300 s sheds
    0.0303 -- the largest drop after the shortest window. `after` says the same
    thing from the other end, rising +0.0472, +0.0558, +0.0647 and then falling
    to +0.0430, which a monotone chemistry under a single bubble cannot do.
    At least two bubbles were growing on different sites at once, and the trace
    is their sum. THAT IS WHY THE GAS CARRIES OVER. `bubble_profile` reduces
    the held gas by a detachment instead of resetting it to zero, and
    `bubble_rate` answers the cumulative demand of the drops so far rather
    than each drop in turn -- because dating every bubble from the previous
    detachment, which is what `bubble_ramp` did until 2026-09-03, has to
    explain this curve's 0.0303 AU with 300 s of growth and can only do it by
    subtracting a ramp steeper than the curve rises.
    """
    curve = {c.sample: c for c in curves_of(experiment)}[int(sample)]
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    drops = bubble_drops(values, curve.noise, sigma=sigma)
    steps = np.diff(values)
    rows = []
    start = 0
    for index in drops:
        window = times[index] - times[start]
        climb = float(values[index] - values[start])
        rows.append({
            "time_s": float(times[index]),
            "lost": float(-steps[index]),
            "sigma": float(steps[index] / curve.noise),
            "grew_s": float(window),
            "rose": climb,
            "ramp": climb / window if window > 0 else np.nan,
            "after": float(values[index + 1]),
        })
        start = index + 1
    return pd.DataFrame(rows)

def hoo_consistency(parameter="vmax", scope=None):
    """
    Is [HOO-] the reactant? Move it two ways and ask whether the order agrees.

    THE ONE TEST ONLY THIS BLOCK CAN RUN. [HOO-] = [H2O2] * Ka / ([H+] + Ka),
    so there are two independent levers on it, and the block moves both:

      * WITHIN a run, at fixed pH, the peroxide arm ladders [H2O2] and [HOO-]
        follows it proportionally. d ln v / d ln[HOO-] is then exactly the
        peroxide order, measured against a per-experiment offset.
      * BETWEEN runs, at matched composition, `ph_order` ladders pH while every
        cuvette's [H2O2] stays put, measured against a per-cuvette offset.

    Nothing is shared between the two: different contrast, different offsets,
    different runs carrying the signal. If the hydroperoxide anion is what the
    chemistry consumes, the two orders are the same number. If the rate
    responded to pH through something else -- the buffer's speciation, the
    catalyst's own ionisation, the substrate -- the pH route would carry that
    too and the two would part.

    Returns the two orders, their difference and its size in sigma. `scope`
    defaults to `strong_runs()`, because a run whose own cuvettes do not
    predict its own rates cannot carry either side (AGREEMENT_FLOOR).
    """
    from induction import ladder_arms
    scope = strong_runs() if scope is None else scope
    arm = ladder_arms(frame(scope))["peroxide arm"]
    within = orders(parameter, frame=arm)
    across = ph_order(parameter, scope=scope)
    pooled = across.loc["pooled"] if "pooled" in across.index else across.iloc[0]
    gap = float(within["order_h2o2"] - pooled.order)
    spread = float(np.hypot(within["stderr_h2o2"], pooled.stderr))
    # WHERE ON THE CURVE EACH CONTRAST SITS. A rate that saturates in [HOO-]
    # has a local order that falls as [HOO-] rises, so two contrasts centred on
    # different levels may differ without anything but saturation happening.
    # Reporting the geometric mean of each is what lets a reader tell that
    # explanation from the other one.
    ladder = pd.concat(ph_ladders(scope).values())
    centres = {}
    for name, rows in (("within", arm), ("across", ladder)):
        live = rows[rows.live & (rows.hoo > 0)]
        centres[name] = float(np.exp(np.log(live.hoo.to_numpy(float)).mean()))
    return {"parameter": parameter,
            "within_hoo": centres["within"],
            "across_hoo": centres["across"],
            "within_order": float(within["order_h2o2"]),
            "within_stderr": float(within["stderr_h2o2"]),
            "within_curves": int(within["n"]),
            "across_order": float(pooled.order),
            "across_stderr": float(pooled.stderr),
            "across_curves": int(pooled.curves),
            "gap": gap,
            "sigma": abs(gap) / spread if spread > 0 else np.nan}


def arm_orders(scope=TWO_AXIS_BLOCK, parameters=ORDER_PARAMETERS):
    """
    Each order read from the arm that moves only its own axis.

    The joint fit of `orders` reads both coefficients from all seven cuvettes
    at once, which is right and is also an extrapolation: the L has no interior
    point, so no cuvette moves both axes and the fit's separation of the two
    rests on the model being additive in logs. Each arm holds the other axis
    fixed within a run, so its slope needs no such assumption. Agreement
    between the two is the check; it is not automatic.
    """
    from induction import ladder_arms
    arms = ladder_arms(frame(scope))
    axis_of = {"substrate arm": "s0", "peroxide arm": "h2o2"}
    rows = []
    for parameter in parameters:
        joint = orders(parameter, scope)
        for arm, table in arms.items():
            axis = axis_of[arm]
            alone = orders(parameter, frame=table)
            rows.append({
                "parameter": parameter, "arm": arm, "axis": axis,
                "curves": alone["n"],
                "arm_order": alone[f"order_{axis}"],
                "arm_stderr": alone[f"stderr_{axis}"],
                "joint_order": joint[f"order_{axis}"],
                "joint_stderr": joint[f"stderr_{axis}"]})
    table = pd.DataFrame(rows)
    table["sigma"] = ((table.arm_order - table.joint_order).abs()
                      / np.hypot(table.arm_stderr, table.joint_stderr))
    return table.set_index(["parameter", "arm"])


# ---------------------------------------------------------------------------
# The +/- chemzyme controls.
#
# Exps 65-71 are consecutive runs from June 2010 in which the same substrate
# ladder, [H2O2], buffer, pH and temperature were run twice, once without the
# chemzyme and once with it at 0.028 mM. They are the only paired controls in
# the archive on BnOH: every other enzyme-free BnOH run (exps 3 and 6) sits at
# a pH and [H2O2] no catalysed run shares.
#
# They are NOT in TWO_AXIS_BLOCK and cannot be pooled with it -- phosphate and
# boric buffer, and only one value of [H2O2] per run, so they carry no
# peroxide order. What they carry is the one comparison the two-axis block
# cannot make at all: the same chemistry with and without the catalyst.
#
# (enzyme-free experiments, catalysed partner, label)
PAIRED_CONTROLS = (
    ((65,), 66, "boric pH 8.51"),
    ((67,), 68, "phosphate pH 8.01"),
    ((69, 70), 71, "phosphate pH 8.01, low rungs"),
)

# Enzyme-free BnOH runs whose [buf] is CONSTANT along the substrate ladder,
# and so the only ones from which a substrate order may be read.
FREE_BNOH = (65, 67, 69, 70)

# The archive holds two more enzyme-free BnOH runs, and they are a trap. Exps 3
# and 6 are buffer titrations: [buf] falls 85 -> 25 mM as [sub] rises 1.28 ->
# 8.98 mM, r = -0.91 and -0.98 against log[sub]. Their rate falls with
# substrate, which looks like the turnover the catalysed block shows and is a
# buffer effect wearing substrate's clothes -- FITTING.md F1 has said so since
# 2026-08-29. Pooling them into FREE_BNOH turns 2 clean rung-pairs above 3 mM
# into 7 and swings the median order from -0.245 to -0.431. Do not.
FREE_BNOH_BUFFER_TITRATIONS = (3, 6)

# Rungs count as the same rung if their [BnOH] agrees to this relative
# tolerance. The ladders were made from one dilution series, so matching is
# exact in practice; the tolerance only absorbs rounding in the sheets.
RUNG_TOLERANCE = 0.02


def paired_controls(controls=PAIRED_CONTROLS):
    """
    One row per substrate rung that exists both with and without the chemzyme.

    Columns: pair, s0, v0_free, v0_enz, vmax_free, vmax_enz, ratio (the
    catalysed vmax over the enzyme-free one), accel_free, accel_enz, live --
    where `live` is True only if BOTH sides carry a signal. Read `ratio` only
    on live rows: a ratio against a dead curve is a ratio against noise, and
    three of the twelve catalysed cuvettes are dead where their enzyme-free
    partner is alive.

    For scale when reading `ratio`: exps 69 and 70 are the same experiment run
    twice, and their vmax disagrees by up to 1.55x rung for rung. Nothing
    inside that factor is a measurement of anything.
    """
    scope = tuple(sorted({e for free, cat, _ in controls
                          for e in (*free, cat)}))
    data = frame(scope)
    rows = []
    for free, catalysed, label in controls:
        left = data[data.experiment.isin(free)]
        right = data[data.experiment == catalysed]
        for s0 in sorted(right.s0.unique()):
            a = left[np.isclose(left.s0, s0, rtol=RUNG_TOLERANCE)]
            b = right[np.isclose(right.s0, s0, rtol=RUNG_TOLERANCE)]
            if not len(a) or not len(b):
                continue
            b = b.iloc[0]
            both_live = bool(a.live.all() and b.live)
            rows.append({
                "pair": label, "s0": float(s0),
                "v0_free": float(a.v0.median()), "v0_enz": float(b.v0),
                "vmax_free": float(a.vmax.median()), "vmax_enz": float(b.vmax),
                "ratio": float(b.vmax / a.vmax.median()) if both_live else np.nan,
                "accel_free": float(a.accel_z.max()),
                "accel_enz": float(b.accel_z),
                "live": both_live,
            })
    return pd.DataFrame(rows)


# The same experiment run twice: exps 69 and 70 share every declared
# condition. Their disagreement is the archive's only direct measure of
# run-to-run reproducibility on BnOH, and so the yardstick any ratio has to
# clear before it means anything.
REPLICATE_PAIR = (69, 70)


def catalytic_effect(controls=PAIRED_CONTROLS, replicate=REPLICATE_PAIR):
    """
    What adding 0.028 mM chemzyme does to BnOH oxidation, on the paired runs.

    Returns the median vmax ratio over the live matched rungs, its range, the
    count, and -- as the yardstick that decides whether the ratio means
    anything -- the largest vmax disagreement between exps 69 and 70, which
    are the same experiment run twice.
    """
    table = paired_controls(controls)
    live = table[table.live]
    first, second = replicate
    data = frame(replicate)
    repeats = []
    for s0 in sorted(data[data.experiment == first].s0.unique()):
        a = float(data[(data.experiment == first) & (data.s0 == s0)].vmax.iloc[0])
        b = float(data[(data.experiment == second) & (data.s0 == s0)].vmax.iloc[0])
        repeats.append(max(a, b) / min(a, b))
    return {
        "rungs": int(len(live)),
        "median_ratio": float(live.ratio.median()),
        "ratio_range": (float(live.ratio.min()), float(live.ratio.max())),
        "replicate_scatter": float(max(repeats)),
        "dead_with_enzyme": int((~table.live).sum()),
    }


# ---------------------------------------------------------------------------
# The literature's own kinetics for this reaction, for scale.
#
# From MECHANISM.md reference 4 (ChemCatChem 2025), the only Bols-group paper
# retrieved in full. Its catalyst (diketone 8) is not necessarily the
# "a-diesterketon" exp 66's sheet names, so this is an order-of-magnitude
# comparison and nothing finer.
LITERATURE = {
    "source": "ChemCatChem 2025, diketone 8 (MECHANISM.md ref 4)",
    "kcat_per_s": 44e-5,
    "km_mM": 1.25,
    "kcat_over_kuncat": 28000,
    "catalyst_mM": 0.4,
    "h2o2_mM": 72.0,
    "pH": 7.0,
}

# The sheets' extinction coefficient for BnOH at 285 nm, used only to turn our
# AU/s into mM/s. What the absorbance actually measures is MECHANISM.md's
# leading open question, so treat any rate in mM/s as uncertain by whatever
# that answer turns out to be.
BNOH_EPSILON = 1.23

# The enzyme-free BnOH runs at NEAR-NEUTRAL pH -- the only ones comparable to
# the literature's pH 7 without extrapolating across a decade of [HOO-].
#
# These are exps 3 and 6, the same buffer titrations FREE_BNOH excludes. The
# collinearity between [buf] and [sub] destroys any SUBSTRATE ORDER read from
# them; it does not destroy the RATE of an individual cuvette, which is what
# this comparison needs. Do not use them for an order. See FITTING.md F1.
FREE_BNOH_NEUTRAL = (3, 6)


# Every enzyme-free BnOH run the manifest keeps: the four with a constant-[buf]
# substrate ladder and the two buffer titrations. Exp 64 is NOT here -- it is
# excluded as an aborted run (7 minutes at dt = 28 s), which is why the boric
# pair 64/65 cannot supply a within-pair [H2O2] contrast.
FREE_BNOH_ALL = FREE_BNOH + FREE_BNOH_NEUTRAL

# The boric-buffer runs, isolated so that "what does borate carry" is a scope
# rather than an argument. In the enzyme-free BnOH data this is exp 65 alone.
#
# MECHANISM.md says to treat boric points as suspect and gives three separate
# reasons, all of which bite hardest exactly where exp 65 sits (pH 8.51):
# borate forms PEROXOBORATE with H2O2 (K = 2.0e-8, significant above pH ~7.7),
# a much faster oxidant than H2O2 itself; it generates DIOXABORIRANE, a
# competing electrophilic oxidant with a 2.8 kcal/mol barrier; and boric acid
# CATALYSES PEROXYACID HYDROLYSIS ~12-fold with a maximum at pH 8.4-9, which
# would be actively destroying the intermediate the mechanism runs through.
#
# The data agreed before that was consulted. Exp 65 is the only run in the set
# that NEITHER the quadratic NOR the burst/lag form fits -- 4.8-5.7x noise and
# 6.9-8.3x respectively -- and its samples 3 and 4 return a negative v0_quad,
# so they silently leave any log-log fit. See DATA_VERIFICATION.md 2026-09-01.
#
# EXCLUDING IT IS NOT FREE. Exp 65 is the only run at pH 8.51 and carries the
# top of the [HOO-] range on its own (0.089 mM against 0.041 and 0.0012), so
# without it the pH axis has two levels and a [HOO-] order taken across them
# is a two-point line that cannot be checked for curvature. That is why this
# is a named alternative scope and not a deletion.
BORIC_BUFFER = (65,)
FREE_BNOH_PHOSPHATE = tuple(e for e in FREE_BNOH_ALL if e not in BORIC_BUFFER)

# EXP 65 HAS NO USABLE RATE, and every default in this module that fits one now
# excludes it. Ruled 2026-09-01. The reason is not that boric is a different
# buffer -- that would be an argument for a sensitivity, which is what
# `boric_sensitivity` was -- it is that the run's curves cannot be reduced to a
# rate at all:
#
#   * all four cuvettes break upward mid-run, at 504-560 s, and steepen by
#     1.82-15.94x across the break (`synchronised_break`). So `v0` measures the
#     pre-break stretch, `vmax` the post-break one, and they are not estimates
#     of one quantity that disagree -- they are two different quantities.
#   * `v0_quad` returns a NEGATIVE rate on two of the four, which silently
#     drops them from any log-log fit, so the "all six runs" law rested on two
#     of exp 65's four cuvettes without saying so.
#   * neither rate form fits it: 6.9-8.3x noise for the burst, 4.8-5.7x for the
#     quadratic, against a median 1.08-1.18x over the rest of the block.
#   * all four burst fits collapse to B = 0 with tau at its grid floor while
#     still reporting bounded = True.
#
# There is no defensible choice of estimator, so there is no defensible number,
# and a sensitivity that reports the fit both ways implies one of them is right.
#
# THIS IS NOT A KNOWN_EXCLUSIONS ENTRY, deliberately. That would drop exp 65
# from the dataset, and its SHAPE is the most informative thing in the boric
# block: it is the only background run in the archive that breaks, which is the
# only direct evidence bearing on whether the buffer makes an oxidant. The
# shape is evidence; the rate is not. `FREE_BNOH_ALL` therefore still exists
# and still contains exp 65, for shape work and for `boric_sensitivity`, which
# is now a record of what the exclusion bought rather than a live alternative.
BORIC_RATE_UNUSABLE = BORIC_BUFFER


def background_orders(scope=FREE_BNOH_PHOSPHATE, terms=("s0", "h2o2", "hoo"),
                      within=False, parameter="vmax", drop_accelerating=False):
    """
    How the enzyme-free rate depends on substrate, peroxide and pH.

    Fitted across all six enzyme-free BnOH runs, which span pH 6.71 to 8.51.
    Returns a dict of orders with standard errors.

    The [HOO-] order is the one that matters here: it comes out near +0.84,
    i.e. the uncatalysed reaction is close to first order in the peroxide
    anion, so its rate climbs about tenfold per pH unit and a background
    measured at pH 8 says almost nothing about a background at pH 7. Ignoring
    that is how this module first reported an 876-fold discrepancy against the
    literature where the honest figure is nearer 34.

    [H2O2] and [HOO-] are collinear (hoo = h2o2 * f(pH)), so read their sum as
    the peroxide dependence and the [HOO-] term as the pH part.

    `terms` is the list of axes to fit. IT DOES NOT DEFAULT TO EVERYTHING, and
    the omission that matters is `buf`. Exps 3 and 6 are buffer titrations --
    [buf] falls 85 -> 25 mM as [sub] rises -- so a fit without a `buf` term
    does not drop the buffer effect, it RELABELS it as a substrate order. See
    BUFFER_CONFOUNDED and `buffer_dependence`. Every returned order carries a
    `vif_` alongside it for the same reason: an order whose variance is
    inflated tenfold is a number the design cannot support, and it looks
    exactly like one that can.
    """
    data = frame(scope)
    data = data[data.live]
    # A log-log order is undefined for a non-positive rate, and some
    # estimators produce them: `v0_quad` extrapolates to t = 0 and returns a
    # negative rate on the two exp 65 cuvettes whose curvature the quadratic
    # cannot hold. Dropping them here rather than propagating a nan means the
    # `n` this function reports is the count it actually fitted.
    data = data[(data[parameter] > 0) & np.isfinite(data[parameter])]
    if drop_accelerating:
        # An "initial rate" read off a curve whose rate is still RISING is the
        # induction rate, not the reaction at the stated concentrations -- the
        # same distinction curve_metrics.peak_rate draws between v0 and vmax.
        # Four enzyme-free BnOH curves accelerate and all four are in the
        # titrations, so this is a sensitivity worth reporting rather than a
        # default: it removes 2 of the 10 live curves the buffer order's
        # titration arm rests on.
        data = data[~data.accelerates]
    y = np.log(data[parameter].to_numpy(dtype=float))
    axes = [np.log(data[c].to_numpy(dtype=float)) for c in terms]
    if within:
        # One indicator per experiment instead of a single intercept, as in
        # `orders`. Anything constant across a run -- pH, [HOO-], [H2O2], the
        # cell, the day -- is absorbed, so the remaining orders are measured
        # only from contrast BETWEEN CUVETTES OF THE SAME RUN. A term that is
        # itself constant within every run is then unidentifiable and must not
        # be in `terms`; it would be collinear with the indicators.
        labels = np.unique(data.experiment.to_numpy())
        design_matrix = np.column_stack(
            axes + [(data.experiment.to_numpy() == e).astype(float)
                    for e in labels])
    else:
        design_matrix = np.column_stack([np.ones(len(data))] + axes)
    coefficients, *_ = np.linalg.lstsq(design_matrix, y, rcond=None)
    residual = y - design_matrix @ coefficients
    dof = max(1, len(data) - np.linalg.matrix_rank(design_matrix))
    variance = float(residual @ residual) / dof
    stderr = np.sqrt(np.diag(variance
                             * np.linalg.pinv(design_matrix.T @ design_matrix)))
    result = {"n": int(len(data)), "terms": tuple(terms), "within": bool(within),
              "dof": int(dof),
              "r2": float(1 - (residual @ residual)
                          / ((y - y.mean()) ** 2).sum())}
    # The orders are the LEADING coefficients in both layouts: `within` puts
    # the indicators after them, the pooled fit puts the intercept first.
    offset = 0 if within else 1
    inflation = variance_inflation(data, terms, within=within)
    for name, value, error in zip(terms, coefficients[offset:], stderr[offset:]):
        result[f"order_{name}"] = float(value)
        result[f"stderr_{name}"] = float(error)
        result[f"vif_{name}"] = inflation[name]
    return result


def variance_inflation(data, terms, within=False):
    """
    Each term's variance inflation factor against the others, on log axes.

    VIF = 1/(1 - R2) where R2 is from regressing one term on the rest. It is
    the factor by which collinearity widens that coefficient's standard error,
    so it says whether an order is measured or merely reported. Above about 10
    the coefficient is arithmetic, not evidence.
    """
    logged = {name: np.log(data[name].to_numpy(dtype=float)) for name in terms}
    result = {}
    for name in terms:
        others = [logged[o] for o in terms if o != name]
        if within:
            labels = np.unique(data.experiment.to_numpy())
            others = others + [(data.experiment.to_numpy() == e).astype(float)
                               for e in labels]
        elif not others:
            result[name] = 1.0
            continue
        matrix = np.column_stack([np.ones(len(data))] + others)
        target = logged[name]
        coefficients, *_ = np.linalg.lstsq(matrix, target, rcond=None)
        residual = target - matrix @ coefficients
        total = float(((target - target.mean()) ** 2).sum())
        r2 = 1 - float(residual @ residual) / total if total > 0 else 0.0
        result[name] = float(1.0 / (1.0 - r2)) if r2 < 1 else float("inf")
    return result


# The two enzyme-free BnOH designs, named by what they can and cannot measure.
#
# BUFFER_CONFOUNDED (exps 3, 6) move [buf] 85 -> 25 mM DOWN as [sub] moves
# 1.28 -> 8.98 mM UP, because substrate volume displaced buffer volume in the
# cuvette. BUFFER_FIXED (exps 67, 69, 70; exp 65 left it on 2026-09-01,
# see BORIC_RATE_UNUSABLE) holds [buf] at 85 mM along
# their substrate ladder. Neither design varies [buf] at constant [sub] --
# no enzyme-free run in the archive does, in any block -- so the buffer order
# is not directly measurable and has to be recovered from the disagreement
# between these two, which is what `buffer_dependence` does.
BUFFER_CONFOUNDED = FREE_BNOH_NEUTRAL
BUFFER_FIXED = tuple(e for e in FREE_BNOH if e not in BORIC_RATE_UNUSABLE)


def buffer_dependence(anchor=BUFFER_FIXED, titration=BUFFER_CONFOUNDED,
                      parameter="vmax", drop_accelerating=False):
    """
    The order in buffer concentration, from the substrate order it corrupts.

    THE DESIGN PROBLEM. No enzyme-free run varies [buf] at constant [sub], so
    no single run measures a buffer order. Two things are measurable instead:

      a   the substrate order where [buf] is CONSTANT (`anchor`), and
      a'  the substrate order where [buf] falls as [sub] rises (`titration`),

    both read within experiments, so both are free of pH, [H2O2], cell and day.
    If the rate goes as [sub]^a [buf]^d, and within the titration runs
    log[buf] = g log[sub] + constant, then fitting the titrations without a
    buffer term returns a' = a + d g, and so

        d = (a' - a) / g.

    WHY NOT JUST FIT BOTH TERMS. Because the titrations cannot carry it: on
    those eight live curves [sub] and [buf] run at VIF 11-14, and the joint fit
    returns d = +0.73 +/- 0.62, consistent with anything from zero to two.
    Pooling all six runs instead does return a tight d, but there [buf] is 85+
    in every pH 8.0-8.5 run and sweeps only in the pH 6.71 ones, so [buf] is
    partly a label for pH -- and the [HOO-] order drops from +0.84 to +0.74
    when the buffer term is added, which is that theft made visible. This
    route uses only within-run contrast on both sides and never asks one
    regression to separate the two.

    THE ASSUMPTION IT RESTS ON, stated because it is not testable here: that
    the substrate order is the same at the titrations' pH 6.71 as at the
    anchor's pH 8.01-8.51. The archive has no run that could check it.

    Returns a dict: order_s0_fixed, order_s0_titration, coupling (g),
    order_buf and its standard error, and the counts behind each.
    """
    clean = background_orders(anchor, terms=("s0",), within=True,
                              parameter=parameter,
                              drop_accelerating=drop_accelerating)
    dirty = background_orders(titration, terms=("s0",), within=True,
                              parameter=parameter,
                              drop_accelerating=drop_accelerating)

    # g: how [buf] tracks [sub] inside the titration runs, on log axes, with a
    # free offset per run -- the same within-experiment contrast the orders use.
    data = frame(titration)
    data = data[data.live]
    labels = np.unique(data.experiment.to_numpy())
    design_matrix = np.column_stack(
        [np.log(data.s0.to_numpy(dtype=float))]
        + [(data.experiment.to_numpy() == e).astype(float) for e in labels])
    coefficients, *_ = np.linalg.lstsq(
        design_matrix, np.log(data.buf.to_numpy(dtype=float)), rcond=None)
    coupling = float(coefficients[0])

    gap = dirty["order_s0"] - clean["order_s0"]
    # g is measured from concentrations that were pipetted, not from a noisy
    # observable, so its own error is negligible beside the two rate orders'.
    error = float(np.hypot(clean["stderr_s0"], dirty["stderr_s0"]) / abs(coupling))
    return {
        "order_s0_fixed": clean["order_s0"],
        "stderr_s0_fixed": clean["stderr_s0"],
        "n_fixed": clean["n"],
        "order_s0_titration": dirty["order_s0"],
        "stderr_s0_titration": dirty["stderr_s0"],
        "n_titration": dirty["n"],
        "coupling": coupling,
        "order_buf": float(gap / coupling),
        "stderr_buf": error,
    }


# The archive's ONLY temperature series: 4OMe-BnOH in 65 mM phosphate at
# pH 7.00 and 82.5 mM H2O2, run at 15, 20, 25, 30, 35 and 40 C with the same
# four-rung substrate ladder (1.850, 3.700, 5.549, 7.399 mM) in every one.
# Listed in temperature order rather than experiment order, because the
# experiment numbers do not run with temperature: 14 is 25 C and 16 is 40 C.
TEMPERATURE_SERIES = (19, 18, 14, 17, 15, 16)

# The 4OMe-BnOH / 40 C enzyme-free runs. Exp 31 is deliberately NOT here: it
# is the same design at 35 C, and temperature moves every rate constant through
# Arrhenius, so including it would pool two cells (FITTING.md F7).
FREE_4OME_40C = (23, 24, 25, 26, 27, 28, 29, 30, 38, 39)


def boric_sensitivity(estimators=("v0_quad", "v0_burst", "vmax", "v0",
                                  "v0_whole"),
                      terms=("s0", "h2o2", "hoo")):
    """
    Every order in the enzyme-free BnOH set, with and without the boric run.

    Returns a DataFrame indexed by (estimator, term) with columns `all`,
    `all_stderr`, `phosphate`, `phosphate_stderr`, plus a `buf` pair from
    `buffer_dependence` -- the buffer order is computed differently, from the
    contrast between the fixed-buffer and titration designs, so it cannot come
    out of the same regression.

    THE POINT IS THE SPREAD ACROSS ESTIMATORS, not any single row. Five ways of
    measuring a rate should give one order; where they do not, the disagreement
    is the measurement's, not the chemistry's. On the substrate order the five
    span 0.24 in log units with exp 65 in and 0.06 without it -- so the boric
    run was the disagreement. The buffer order, by contrast, barely notices it,
    which is the robustness the headline actually rests on.

    See BORIC_BUFFER for why excluding it is defensible and what it costs.
    """
    rows = []
    for estimator in estimators:
        both = {}
        for label, scope_ in (("all", FREE_BNOH_ALL),
                              ("phosphate", FREE_BNOH_PHOSPHATE)):
            both[label] = background_orders(scope_, terms=terms,
                                            parameter=estimator)
        for term in terms:
            rows.append({
                "estimator": estimator, "term": term,
                "all": both["all"].get(f"order_{term}", np.nan),
                "all_stderr": both["all"].get(f"stderr_{term}", np.nan),
                "phosphate": both["phosphate"].get(f"order_{term}", np.nan),
                "phosphate_stderr": both["phosphate"].get(f"stderr_{term}",
                                                          np.nan),
            })
        # The buffer order comes from the two-design contrast, so its two
        # versions are the anchor WITH and WITHOUT exp 65. Both anchors are
        # named explicitly: BUFFER_FIXED stopped containing exp 65 on
        # 2026-09-01, and reading the "all" column off the default silently
        # made both columns the same number.
        full = buffer_dependence(anchor=FREE_BNOH, parameter=estimator)
        cut = buffer_dependence(anchor=BUFFER_FIXED, parameter=estimator)
        rows.append({
            "estimator": estimator, "term": "buf",
            "all": full["order_buf"], "all_stderr": full["stderr_buf"],
            "phosphate": cut["order_buf"], "phosphate_stderr": cut["stderr_buf"],
        })
    table = pd.DataFrame(rows).set_index(["estimator", "term"])
    return table


def boric_spread(table=None):
    """
    How far the five estimators disagree on each order, with boric and without.

    Returns a DataFrame indexed by term with columns `all`, `phosphate` --
    the max-minus-min across estimators. This is the summary that decides
    whether excluding boric helped: a spread that shrinks means the estimators
    were disagreeing about exp 65 rather than about the chemistry.
    """
    if table is None:
        table = boric_sensitivity()
    spread = table.groupby("term").agg(
        all=("all", lambda v: float(v.max() - v.min())),
        phosphate=("phosphate", lambda v: float(v.max() - v.min())))
    return spread


def buffer_cross_check(scope=FREE_4OME_40C):
    """
    The buffer order again, on the 4OMe-BnOH / 40 C block, independently.

    That block is the only other place with buffer contrast: three substrate
    values recur at two buffer levels each across its experiments (0.38 mM at
    75 and 90, 2.06 at 80 and 85, 4.12 at 70 and 75). The contrast is BETWEEN
    experiments, so this fit takes no per-experiment offsets -- offsets would
    absorb the very thing being measured. pH (6.97-7.00) and [H2O2] (82.5 mM)
    are constant across the block, so there is no pH for [buf] to proxy.

    Different substrate and temperature, so this checks the SIGN and rough
    SIZE of `buffer_dependence`, not its value.
    """
    return background_orders(scope, terms=("s0", "buf"), within=False)

# The one matched pair in the enzyme-free archive that changes the buffer SALT
# and almost nothing else. Exps 65 and 67 ran the same substrate ladder
# (7.310 / 3.655 / 1.827 / 0.365 mM) at the same [H2O2] (122.426 mM), the same
# temperature, the same instrument, and the same .rre source, at 87.5 against
# 85.0 mM buffer. Only the salt -- boric against phosphate -- and the pH
# (8.51 against 8.01) differ.
PEROXO_PAIR = (65, 67)
# The SAME pair with enzyme, and the reason it is needed: exp 65's curves have
# a synchronised mid-run break that no other run in the block has, so its rate
# numbers do not describe one process (`synchronised_break`,
# DATA_VERIFICATION.md 2026-09-01). Exps 66 and 68 match on enzyme (0.028 mM),
# buffer (85.0 mM), peroxide (122.426 mM), substrate (2.741 mM BnOH) and
# temperature, differing only in salt and pH, and both run smooth.
#
# WITHDRAWN AS A PEROXO PROBE, same day, one hour later. Both runs are
# catalysed, so their reference channels omit the ENZYME and the background is
# subtracted out of both. A buffer-made oxidant acts on the background, which
# means it is present in BOTH beams of each run and CANCELS. This pair cannot
# detect the thing it was added to detect -- not confounded, blind. It is kept
# only as the catalysed boric-vs-phosphate comparison it actually is, which is
# a statement about the CATALYSED reaction in the two buffers.
CATALYSED_PEROXO_PAIR = (66, 68)


# The archive's only pair of runs that differ in [enz] and in nothing else that
# matters: BnOH, boric, pH 8.51, 25 C, 122.426 mM H2O2, the same substrate
# ladder, [enz] 0.028 against 0.014 -- exactly 2.000x, from 0.01 ml of a
# 5.596 mM stock against 0.05 ml of a 0.5596 mM one. Buffer differs 77.0 vs
# 75.0 mM, a 2.6% gap that a buffer order of ~0.4 turns into about 1%.
SELWYN_PAIR = (59, 60)


def selwyn_test(pair=SELWYN_PAIR, levels=(4, 8, 12, 16, 18)):
    """
    Is the curvature in these progress curves catalyst INACTIVATION?

    THE CLASSICAL TEST (Selwyn 1965). Whenever the rate is proportional to
    active catalyst and the departure from linearity comes from the SOLUTION --
    substrate depletion, product inhibition -- product is a function of
    [E]0 x t alone, so progress curves at different [E]0 superimpose on that
    axis. Catalyst inactivation breaks the superposition, because inactivation
    runs on real time rather than on [E]0 t: the low-[E] run needs twice the
    real time to reach the same [E]0 t and its catalyst has decayed twice as
    long, so it lands BELOW.

    So the ratio P(low [E]) / P(high [E]) at matched [E]0 t is the statistic:
    1.00 superimposes, below 1.00 is inactivation, above 1.00 is neither and
    means the rate is sub-first-order in catalyst.

    ONLY THE OVERLAP IS USABLE. The low-[E] run reaches half the [E]0 t of the
    high-[E] one in the same wall-clock time, so the comparison stops at the
    shorter axis -- 18.4 mM s here. Interpolating past it silently compares a
    curve with its own last reading, which is what the first draft of this did.

    Returns a DataFrame indexed by the rung's [S], one column per level.
    """
    low_experiment, high_experiment = (
        sorted(pair, key=lambda e: frame((e,)).e0.median()))
    curves = {}
    for experiment in pair:
        for curve in curves_of(experiment):
            curves[(experiment, round(float(curve.conditions.s0), 3))] = curve
    rows = []
    for s0 in sorted({key[1] for key in curves}):
        high = curves.get((high_experiment, s0))
        low = curves.get((low_experiment, s0))
        if high is None or low is None:
            continue
        high_axis = np.asarray(high.times) * high.conditions.e0
        low_axis = np.asarray(low.times) * low.conditions.e0
        ceiling = min(high_axis.max(), low_axis.max())
        row = {"s0": float(s0)}
        for level in levels:
            row[level] = (
                float(np.interp(level, low_axis, low.absorbance)
                      / np.interp(level, high_axis, high.absorbance))
                if level <= ceiling else np.nan)
        rows.append(row)
    return pd.DataFrame(rows).set_index("s0")


def catalyst_order(pair=SELWYN_PAIR, levels=(4, 8, 12, 16, 18)):
    """
    The order in catalyst implied by `selwyn_test`, and the inactivation verdict.

    At matched [E]0 t, a rate going as [E]^n gives P proportional to
    [E]^(n-1), so the ratio between a pair differing twofold is 0.5^(n-1) and

        n = 1 + ln(ratio) / ln(0.5).

    Returns a dict with the median ratio, n, the ratio's range, and
    `inactivation`, which is True only if the ratio falls below 1.
    """
    table = selwyn_test(pair, levels)
    values = table.to_numpy(dtype=float).ravel()
    values = values[np.isfinite(values)]
    median = float(np.median(values))
    return {"median_ratio": median,
            "lowest_ratio": float(values.min()),
            "highest_ratio": float(values.max()),
            "order_in_catalyst": float(1 + np.log(median) / np.log(0.5)),
            "inactivation": bool(values.min() < 1.0),
            "n": int(len(values)), "pair": tuple(pair)}


def synchronised_break(scope=FREE_BNOH_ALL):
    """
    Per run: do its cuvettes break at the same TIME, and do they steepen?

    A run's four cuvettes differ only in substrate -- same buffer, same
    peroxide, same cell, same day. So a breakpoint they SHARE cannot be driven
    by the substrate, and a break that is also a STEEPENING cannot be the
    reaction decelerating toward conversion. The two together are the
    signature this function reports.

    Returns a DataFrame indexed by experiment: `span` (max - min break time,
    in seconds), `median_break`, `max_ratio`, `steep` (how many cuvettes
    exceed SEGMENT_RATIO_STEEP), `n`, and the sorted break times and ratios.

    ONLY ASK IT OF BACKGROUND RUNS. A shape in the background is invisible in a
    catalysed run BY CONSTRUCTION: an enzyme run's reference channel omits the
    enzyme, so the background is in both beams and cancels, and the curve is a
    catalytic increment (`frame`'s `differential` column, kinetics_io,
    DATA_VERIFICATION.md 2026-08-31). This function raises if `scope` mixes the
    two, because the first control set chosen for the boric probe was four
    catalysed runs and the conclusion drawn from their smoothness -- that exp
    65's break was not borate chemistry -- did not follow.

    WHAT IT FOUND, on the 20 background experiments -- the whole un-subtracted
    population -- 18 of which have live curves. Seventeen are phosphate and
    none of them breaks: ratios 0.22-1.23, with one isolated cuvette of exp
    3's six at 2.44. Exactly one boric run has live curves -- exp 65 -- and
    all four of its cuvettes steepen, by 1.82,
    2.04, 5.59 and 15.94, across breaks spanning 56 s, two of its 28 s
    sampling intervals. The steepening is LARGEST at the lowest substrate
    (15.94 at 0.365 mM against 1.82 at 7.310 mM), so the clock is not the
    substrate.

    That is 1 boric run of 1 showing it against 17 of 17 not, which is
    consistent with borate chemistry and rests on a single run. The other
    boric background run, exp 64, was aborted at 448 s -- BEFORE exp 65's
    break -- and is dead besides, so it cannot test it. The missing experiment
    is a repeat of exp 65. See DATA_VERIFICATION.md 2026-09-01.
    """
    data = frame(scope)
    data = data[data.live & np.isfinite(data.break_time)]
    if data.differential.nunique() > 1:
        mixed = sorted(data[data.differential].experiment.unique())
        raise ValueError(
            "synchronised_break was given both catalysed and background runs "
            f"({mixed} are catalysed). A catalysed curve is an increment whose "
            "reference channel already subtracted the background, so it cannot "
            "show a background shape and its smoothness is not evidence about "
            "one. Pass background runs only.")
    rows = []
    for experiment, group in data.groupby("experiment"):
        order = group.sort_values("break_time")
        rows.append({
            "experiment": int(experiment),
            "n": int(len(group)),
            "span": float(group.break_time.max() - group.break_time.min()),
            "median_break": float(group.break_time.median()),
            "max_ratio": float(group.break_ratio.max()),
            "steep": int((group.break_ratio > SEGMENT_RATIO_STEEP).sum()),
            "breaks": [int(v) for v in order.break_time],
            "ratios": [round(float(v), 2) for v in order.break_ratio],
        })
    return pd.DataFrame(rows).set_index("experiment")


def peroxo_buffer_test(pair=PEROXO_PAIR, orders_scope=FREE_BNOH_PHOSPHATE,
                       estimators=("v0", "vmax", "v0_whole", "v0_quad")):
    """
    Does a buffer that DOES form a peroxo species run faster than the law?

    THE QUESTION. The enzyme-free rate is first order in buffer
    (`buffer_dependence`), and two mechanisms give that: general acid/base
    catalysis, or the buffer making an oxidant -- phosphate + H2O2 ->
    peroxomonophosphate. Within the phosphate runs those are indistinguishable
    (see background_reaction/ANALYSIS.md section 6b: log[buf], log[H2PO4-] and
    log[HPO4^2-] are the same variable, correlation 1.000000).

    THE WAY ROUND IT. Borate is the buffer where the peroxo route is not a
    hypothesis. MECHANISM.md item 39 has B(OH)3 + H2O2 -> peroxoborate with
    K = 2.0e-8, "significant above pH ~ 7.7", and the anionic peroxoborates
    are much faster oxidants than H2O2 itself. Exp 65 is boric buffer at
    pH 8.51 -- above that threshold, at 122 mM H2O2 -- so a substantial part
    of its boron is peroxoborate. If a buffer-derived peroxo oxidant is what
    carries a first-order buffer term, exp 65 must run far above a rate law
    fitted without one. It does not.

    READ THIS BEFORE QUOTING THAT. On the default pair the test is WEAK, and
    it is weak for a reason found after it was written: exp 65's four cuvettes
    share a mid-run breakpoint at 504-560 s across which every one of them
    STEEPENS, by 1.82 to 15.94x, most at the LOWEST substrate. No other run in
    the block does this -- every other one decelerates -- so a single rate
    number does not describe exp 65 at all: `vmax` reads the post-break
    stretch and `v0` the pre-break one, and they are not the same process. The
    excesses below are therefore a comparison between two different things,
    and the right reading of them is "borate is nowhere fast", not "borate
    matches the law".

    AND THERE IS NO SECOND PROBE. CATALYSED_PEROXO_PAIR was added as one and
    withdrawn the same day: both its runs are catalysed, so a buffer-made
    oxidant sits in both beams and cancels. Among the 21 enzyme-free
    background experiments -- the entire population in which a background
    feature is even visible -- exactly one is boric, and it is exp 65. So this test rests on
    one run and that run is the one with the break. See `synchronised_break`
    and DATA_VERIFICATION.md 2026-09-01. The missing experiment is a repeat of
    exp 65: enzyme-free, boric, run long.

    HOW THE PREDICTION IS MADE. The rate law is fitted on `orders_scope`,
    which excludes the boric run, so exp 65 is out of sample. Because the two
    runs share [S] and [H2O2] exactly, the predicted ratio depends only on the
    [buf] and [HOO-] orders; the substrate and peroxide orders, the two worst
    determined here, drop out. Matching is per cuvette on `s0`, not on run
    medians.

    WHAT IT DOES NOT SETTLE. pH is not matched, so the prediction leans on the
    [HOO-] order holding from 8.01 up to 8.51 -- an extrapolation, and exp 65
    is the only run there. It is one run of four cuvettes against one run of
    four, on a day and a cell that are not controlled. Exp 65 is also the run
    neither rate form fits (section 3a) and its noise runs 1.5-2.8x exp 67's.
    A null result here is evidence against the peroxo route, not proof, and it
    says nothing directly about phosphate: it says the mechanism does not show
    up where it certainly operates.

    Returns a DataFrame indexed by estimator: `predicted` (boric/phosphate
    from the phosphate-only law), `observed` (median over matched cuvettes),
    `excess` = observed / predicted, `n` matched cuvettes, and the per-cuvette
    ratios in `ratios`.
    """
    boric_experiment, phosphate_experiment = pair
    data = frame(tuple(sorted(set(pair))))
    data = data[data.live]
    left = data[data.experiment == boric_experiment]
    right = data[data.experiment == phosphate_experiment]

    rows = []
    for estimator in estimators:
        if orders_scope is None:
            # UNCORRECTED, for the catalysed pair. The enzyme-free rate law
            # does not describe a catalysed comparison and there is no
            # catalysed law in this buffer to put in its place, so quote the
            # raw ratio against the [HOO-] ratio and let the reader correct:
            # boric carries 2.19x the hydroperoxide, so anything at or below
            # 1.00 is boric running SLOWER than pH alone would give it.
            law = {"order_buf": 0.0, "order_hoo": 0.0}
        else:
            law = background_orders(orders_scope, parameter=estimator,
                                    terms=("s0", "h2o2", "hoo", "buf"))
        # [S] and [H2O2] are identical between the two runs, so only these two
        # terms survive the ratio. Asserted rather than assumed: a pair that
        # did not match on them would need the other two orders as well.
        predicted = ((float(left.buf.median()) / float(right.buf.median()))
                     ** law["order_buf"]
                     * (float(left.hoo.median()) / float(right.hoo.median()))
                     ** law["order_hoo"])
        ratios = {}
        for s0 in sorted(set(left.s0) & set(right.s0), reverse=True):
            # A log-log order is undefined for a non-positive rate, and
            # `v0_quad` returns one on two of exp 65's cuvettes: skip the
            # cuvette rather than propagate a nan into the median.
            a = left[left.s0 == s0][estimator].to_numpy(dtype=float)
            b = right[right.s0 == s0][estimator].to_numpy(dtype=float)
            if len(a) != 1 or len(b) != 1 or not (a[0] > 0 and b[0] > 0):
                continue
            ratios[float(s0)] = float(a[0] / b[0])
        observed = float(np.median(list(ratios.values()))) if ratios else np.nan
        rows.append({"estimator": estimator,
                     "order_buf": law["order_buf"],
                     "order_hoo": law["order_hoo"],
                     "predicted": float(predicted),
                     "observed": observed,
                     "excess": observed / float(predicted),
                     "n": len(ratios),
                     "ratios": ratios})
    return pd.DataFrame(rows).set_index("estimator")


def literature_comparison(scope=FREE_BNOH_NEUTRAL,
                          orders_scope=FREE_BNOH_PHOSPHATE,
                          orders_terms=("s0", "h2o2", "hoo")):
    """
    Our enzyme-free background against the literature's uncatalysed rate.

    Both are put at the literature's conditions -- pH 7.0, 72 mM H2O2 -- by
    scaling our rates with the orders `background_orders` measures. Using the
    pH 8.0-8.5 runs instead and correcting only [H2O2] inflates the answer
    about 25-fold, which is a statement about [HOO-] and not about our cuvettes.

    Returns a frame with one row per live cuvette plus a `summary` attribute
    holding the median excess and the enhancement the literature's kcat would
    produce at OUR catalyst loading -- which is the number that explains why no
    enhancement is visible anywhere in this archive.
    """
    orders = background_orders(orders_scope, terms=orders_terms)
    h2o2_order, hoo_order = orders["order_h2o2"], orders["order_hoo"]
    kuncat = LITERATURE["kcat_per_s"] / LITERATURE["kcat_over_kuncat"]

    data = frame(scope)
    rows = []
    for _, row in data[data.live].iterrows():
        # to the literature's [H2O2], then across the pH gap through [HOO-]
        factor = ((LITERATURE["h2o2_mM"] / row.h2o2) ** h2o2_order
                  * (10 ** (LITERATURE["pH"] - row.pH)) ** hoo_order)
        ours = float(row.vmax) / BNOH_EPSILON * factor
        theirs = kuncat * float(row.s0)
        rows.append({"experiment": int(row.experiment), "pH": float(row.pH),
                     "s0": float(row.s0), "ours_mM_s": ours,
                     "literature_mM_s": theirs, "excess": ours / theirs})
    table = pd.DataFrame(rows)
    return table


def background_model(scope=FREE_BNOH_NEUTRAL,
                     orders_scope=FREE_BNOH_PHOSPHATE,
                     orders_terms=("s0", "h2o2", "hoo")):
    """
    An amplitude and three orders that predict the enzyme-free rate, in mM/s.

    The orders come from all six enzyme-free runs (`background_orders`); the
    amplitude is anchored on the near-neutral ones, so the model is pinned
    where it is compared to the literature and extrapolated -- not fitted --
    across the pH gap to the catalysed runs.
    """
    orders = background_orders(orders_scope, terms=orders_terms)
    data = frame(scope)
    data = data[data.live]
    exponents = (orders["order_s0"], orders["order_h2o2"], orders["order_hoo"])
    predicted = (data.s0 ** exponents[0] * data.h2o2 ** exponents[1]
                 * data.hoo ** exponents[2])
    amplitude = float(np.median(data.vmax / BNOH_EPSILON / predicted))
    return amplitude, exponents


# The catalysed runs that have an enzyme-free counterpart to be judged against:
# the three paired-control partners plus exps 73 and 83.
CATALYSED_WITH_BACKGROUND = (66, 68, 71, 73, 83)


def predicted_enhancement(scope=CATALYSED_WITH_BACKGROUND,
                          background_scope=FREE_BNOH_NEUTRAL,
                          orders_scope=FREE_BNOH_PHOSPHATE):
    """
    What the literature's kcat would show at THIS archive's catalyst loading.

    For every live cuvette of the catalysed runs that have an enzyme-free
    counterpart -- exps 66, 68, 71 (the paired controls) and 73, 83 -- this
    predicts the background at that cuvette's own conditions, adds the
    catalytic contribution kcat*E0*S/(Km+S), and reports the ratio to the
    background alone. That ratio is what the experiment would have had to
    resolve.

    It comes out at a median 1.3x, range about 1.15-1.9x. Exps 69 and 70 are
    the SAME experiment run twice and their vmax disagrees by up to 1.55x. So
    the enhancement these runs were capable of detecting is smaller than their
    own reproducibility, and the observed 0.63x is not evidence about the
    catalyst. At the literature's own 0.4 mM the predicted ratio is above 40x,
    which nothing could miss; no BnOH run in this archive exceeds 0.069 mM.
    """
    amplitude, (a, b, c) = background_model(background_scope, orders_scope)
    rows = []
    for experiment in scope:
        data = frame((experiment,))
        for _, row in data[data.live].iterrows():
            background = amplitude * row.s0 ** a * row.h2o2 ** b * row.hoo ** c
            catalysed = (LITERATURE["kcat_per_s"] * row.e0 * row.s0
                         / (LITERATURE["km_mM"] + row.s0))
            rows.append({
                "experiment": experiment, "pH": float(row.pH),
                "s0": float(row.s0), "e0": float(row.e0),
                "observed_mM_s": float(row.vmax) / BNOH_EPSILON,
                "background_mM_s": float(background),
                "catalysed_mM_s": float(catalysed),
                "expected_ratio": float((background + catalysed) / background),
                "at_literature_loading": float(
                    (background + catalysed * LITERATURE["catalyst_mM"] / row.e0)
                    / background),
            })
    return pd.DataFrame(rows)


# The scopes worth a name. Anything else is spelled out on the command line as
# experiment numbers, so a one-off question does not need a constant.
NAMED_SCOPES = {
    "two-axis": TWO_AXIS_BLOCK,
    "free-bnoh": FREE_BNOH,
    "free-bnoh-all": FREE_BNOH_ALL,
    "free-bnoh-neutral": FREE_BNOH_NEUTRAL,
    "free-bnoh-phosphate": FREE_BNOH_PHOSPHATE,
    "boric": BORIC_BUFFER,
    "temperature-series": TEMPERATURE_SERIES,
    "paired": tuple(sorted({e for free, cat, _ in PAIRED_CONTROLS
                            for e in (*free, cat)})),
}


def parse_scope(text):
    """
    A scope from a name, or from experiment numbers: "3,6" or "135-151".

    Returns a frozenset. Raises ValueError on anything it cannot read, rather
    than quietly returning an empty scope -- an empty scope produces an empty
    frame, and an empty frame is a table of zeroes that looks like a result.
    """
    if text in NAMED_SCOPES:
        return frozenset(NAMED_SCOPES[text])
    experiments = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part.lstrip("-"):
            low, high = part.split("-", 1)
            experiments.update(range(int(low), int(high) + 1))
        else:
            experiments.add(int(part))
    if not experiments:
        raise ValueError(f"empty scope: {text!r}")
    return frozenset(experiments)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--scope", default="two-axis",
                        help="a named scope (%s) or experiment numbers "
                             "(\"3,6\", \"135-151\"); default two-axis"
                             % ", ".join(NAMED_SCOPES))
    parser.add_argument("--design", action="store_true",
                        help="print the per-experiment design table")
    parser.add_argument("--orders", action="store_true",
                        help="print the apparent reaction orders")
    parser.add_argument("--controls", action="store_true",
                        help="print the +/- chemzyme paired controls")
    parser.add_argument("--literature", action="store_true",
                        help="compare the paired controls against the literature")
    parser.add_argument("--buffer", action="store_true",
                        help="the enzyme-free rate's dependence on [buf]")
    arguments = parser.parse_args()
    chosen = parse_scope(arguments.scope)

    if arguments.design:
        table = design(chosen)
        with pd.option_context("display.width", 200, "display.max_columns", 20):
            print(table.to_string(float_format=lambda v: f"{v:.3g}"))
        return 0

    if arguments.orders:
        with pd.option_context("display.width", 200):
            print(order_table(chosen).to_string(
                float_format=lambda v: f"{v:.3f}"))
        return 0

    if arguments.controls:
        with pd.option_context("display.width", 200):
            print(paired_controls().to_string(
                index=False, float_format=lambda v: f"{v:.4g}"))
        effect = catalytic_effect()
        print(f"\nvmax(+chemzyme)/vmax(-chemzyme) over "
              f"{effect['rungs']} live matched rungs: "
              f"median {effect['median_ratio']:.2f}x, range "
              f"{effect['ratio_range'][0]:.2f}-{effect['ratio_range'][1]:.2f}x")
        print(f"same experiment run twice (exps 69 vs 70) disagrees by up to "
              f"{effect['replicate_scatter']:.2f}x -- the ratio is not "
              f"resolved from no effect")
        return 0

    if arguments.buffer:
        print("enzyme-free BnOH: the substrate order depends on what the "
              "BUFFER was doing\n")
        result = buffer_dependence()
        print(f"  [buf] held constant   exps {BUFFER_FIXED}, "
              f"n={result['n_fixed']:2d}:  order in [sub] "
              f"{result['order_s0_fixed']:+.3f} +/- {result['stderr_s0_fixed']:.3f}")
        print(f"  [buf] falling         exps {BUFFER_CONFOUNDED}, "
              f"n={result['n_titration']:2d}:  order in [sub] "
              f"{result['order_s0_titration']:+.3f} +/- "
              f"{result['stderr_s0_titration']:.3f}")
        print(f"\n  The same reaction reads a POSITIVE substrate order where "
              f"[buf] is held and a\n  NEGATIVE one where [buf] falls as "
              f"[sub] rises. The difference is the buffer.\n")
        print(f"  coupling  dlog[buf]/dlog[sub] within the titrations: "
              f"{result['coupling']:+.3f}")
        print(f"  => order in [buf] = {result['order_buf']:+.2f} +/- "
              f"{result['stderr_buf']:.2f}   (approximately FIRST order)\n")
        check = buffer_cross_check()
        print(f"  independent cross-check, 4OMe-BnOH / 40 C / phosphate, "
              f"n={check['n']} (between-run\n  contrast at fixed pH and "
              f"[H2O2], different substrate and temperature):")
        print(f"      order in [buf]  {check['order_buf']:+.2f} +/- "
              f"{check['stderr_buf']:.2f}   (VIF {check['vif_buf']:.1f})")
        print(f"      order in [sub]  {check['order_s0']:+.2f} +/- "
              f"{check['stderr_s0']:.2f}")
        with_boric = buffer_dependence(anchor=FREE_BNOH)
        print(f"\n  The anchor is exps {BUFFER_FIXED} -- PHOSPHATE ONLY. Exp 65 "
              f"left it on 2026-09-01:\n  its curves break mid-run and have no "
              f"usable rate (BORIC_RATE_UNUSABLE,\n  `synchronised_break`). "
              f"With it in, this read {with_boric['order_buf']:+.2f} +/- "
              f"{with_boric['stderr_buf']:.2f}.")

        print(f"\n  WHAT MAKES IT FIRST ORDER is NOT SETTLED. Catalysis by a "
              f"buffer species, or the\n  buffer making an oxidant? Both are "
              f"first order in a buffer species and the\n  phosphate design "
              f"cannot separate them. Borate is where the second is not a\n  "
              f"hypothesis, and exp {PEROXO_PAIR[0]} is the only boric "
              f"BACKGROUND run in the archive -- so the\n  test rests on the "
              f"one run whose curves are unusable. Predicting it from a law\n  "
              f"fitted on phosphate alone ({PEROXO_PAIR[0]} against "
              f"{PEROXO_PAIR[1]}, cuvette for cuvette):\n")
        peroxo = peroxo_buffer_test()
        for estimator, row in peroxo.iterrows():
            print(f"      {estimator:9s} predicted {row['predicted']:.2f}x   "
                  f"observed {row['observed']:.2f}x   "
                  f"excess {row['excess']:.2f}x   (n={int(row['n'])})")
        print(f"\n  Read as 'borate is nowhere fast', not as 'borate matches "
              f"the law': exp 65's\n  vmax reads the post-break stretch and "
              f"its v0 the pre-break one. The catalysed\n  runs cannot stand "
              f"in -- their reference channel already subtracted the\n  "
              f"background, so a buffer-made oxidant cancels. See "
              f"background_reaction/\n  ANALYSIS.md section 6b; the decisive "
              f"tests are 31P NMR and a repeat of exp 65.")

        print(f"\n  The two-axis block is UNAFFECTED: [buf] = 75.013 mM in "
              f"all 119 of its curves,\n  so no buffer variation can reach "
              f"its substrate order.")
        return 0

    if arguments.literature:
        orders = background_orders()
        print(f"enzyme-free BnOH background, {orders['n']} live curves "
              f"(R2 {orders['r2']:.2f}):")
        for name in ("s0", "h2o2", "hoo"):
            print(f"  order in {name:5s} {orders['order_' + name]:+.2f} "
                  f"+/- {orders['stderr_' + name]:.2f}")
        print(f"\n  the [HOO-] order is why pH matters: the background climbs "
              f"about tenfold per pH unit,\n  so a background measured at pH 8 "
              f"says little about one at pH 7.\n")

        table = literature_comparison()
        print(f"at the literature's pH {LITERATURE['pH']:.1f} / "
              f"{LITERATURE['h2o2_mM']:.0f} mM H2O2, our near-neutral "
              f"enzyme-free runs {FREE_BNOH_NEUTRAL} against its uncatalysed rate:")
        with pd.option_context("display.width", 200):
            print(table.to_string(index=False,
                                  float_format=lambda v: f"{v:.4g}"))
        print(f"\n  median excess: {table.excess.median():.0f}x  "
              f"(range {table.excess.min():.0f}-{table.excess.max():.0f}x)")

        enhancement = predicted_enhancement()
        print(f"\nwhat the literature's kcat would show at this archive's "
              f"loading:")
        print(f"  predicted enhancement over background: median "
              f"{enhancement.expected_ratio.median():.2f}x "
              f"(range {enhancement.expected_ratio.min():.2f}-"
              f"{enhancement.expected_ratio.max():.2f}x)")
        print(f"  observed:                              median "
              f"{(enhancement.observed_mM_s / enhancement.background_mM_s).median():.2f}x")
        print(f"  at the literature's {LITERATURE['catalyst_mM']} mM loading:   "
              f"    median "
              f"{enhancement.at_literature_loading.median():.0f}x")
        print(f"\n  exps 69 and 70 are the same experiment run twice and "
              f"disagree by up to "
              f"{catalytic_effect()['replicate_scatter']:.2f}x, so the "
              f"predicted\n  enhancement is smaller than the reproducibility. "
              f"These runs could not have seen it.")
        print(f"\n  source: {LITERATURE['source']}")
        return 0

    facts = summary(chosen)
    cells = blocks(chosen)
    span = (f"exps {min(chosen)}-{max(chosen)}"
            if sorted(chosen) == list(range(min(chosen), max(chosen) + 1))
            else "exps " + ",".join(str(e) for e in sorted(chosen)))
    print(f"scope        {span}")
    for (substrate, temperature, buffer_name), row in cells.iterrows():
        print(f"block        {substrate} / {temperature:.0f} C / {buffer_name}"
              f"  ({row.curves} curves, {row.experiments} experiments)")
    print(f"curves       {facts['curves']} over {facts['experiments']} experiments, "
          f"{facts['live_curves']} with a live signal")
    print(f"pH           {facts['pH_range'][0]:.2f} to {facts['pH_range'][1]:.2f}, "
          f"{facts['hoo_decades']:.1f} decades of [HOO-]")
    print(f"contrast     log[S] {100 * facts['within_experiment_s0']:.1f}% "
          f"within-experiment, log[H2O2] "
          f"{100 * facts['within_experiment_h2o2']:.1f}%")
    print(f"lag          {facts['lagging']} curves peak after "
          f"{LAG_THRESHOLD:.0%} of the run")
    print(f"acceleration {facts['accelerating_live']} of "
          f"{facts['live_curves']} live curves are steeper later than at "
          f"the start, by >{ACCELERATION_SIGMA:.0f} sigma")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
