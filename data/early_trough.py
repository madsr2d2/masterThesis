"""
Whether the catalyst transiently starves its own reference-subtracted signal.

A catalysed curve's absorbance is already reference-subtracted: the recorded
value is sample minus an enzyme-free reference cuvette holding the same
composition, so whatever reaction the two cuvettes share -- the uncatalysed
background chemistry `background_reaction/` measures directly -- cancels
between them. It cancels only while the two cuvettes see the same free
concentration of whatever the background reaction runs on. If the catalyst
engages one reactant fast enough to measurably deplete it in the SAMPLE
cuvette alone, the sample's own share of the shared background reaction runs
slower than the reference's for a while, and the reported curve dips BELOW
zero before the catalysed rate overtakes it and the curve turns positive.

Two reactants can be scarce enough, relative to the catalyst's own fixed
concentration, for this to be visible: the hydroperoxide anion, which is
nanomolar across most of the archive's pH range and the dominant driver here
(`dominance_correlation` finds it significant in BOTH substrates
independently -- and specifically HOO-, not total H2O2: the pH-weighted
anion ratio beats the un-weighted total-peroxide ratio in both substrates,
decisively so in 4OMe-BnOH, where the total-peroxide version is not even
correctly signed), and the substrate itself, on the rare curves where the
substrate rung is low enough that the catalyst is a non-negligible fraction
of it. Neither needs anything not already established elsewhere in this
project: a real background reaction (`background_reaction/`), reference
subtraction (the instrument design `BUBBLES.md` already relies on for the
gas argument), and a catalyst that engages the peroxide non-productively
(`MECHANISM.md` S4, `BUBBLES.md`'s gas). This is that same engagement, read
off the OTHER cuvette's signal instead of the sample's own gas.

The trough's own relaxation time (`tau_fast`, already fitted for every curve
by `summary_kinetics.fit_progress` -- a trough IS that functional form's
"lag" shape taken far enough that `B/tau` exceeds `v_ss`) gives a
pseudo-first-order rate constant for the step, `binding_rates`, with no new
fitting. `arrhenius_check` is the honest temperature check on it -- and the
honest answer is that the archive cannot resolve one: the scatter among
curves sharing a single temperature is nearly as wide as the whole
15-40 C range. `buffer_comparison` finds the apparent rate constant DOES
depend on which buffer is present, pyrophosphate running roughly 6x the
phosphate curves' geometric mean even after normalising by [enz] the same
way throughout -- evidence that what is measured is not a clean elementary
step, more likely a fast, cuvette-symmetric buffer-HOO- pre-equilibrium
setting the size of the pool the catalyst draws from, on top of the
catalyst's own, slower engagement with it.

    from early_trough import trough_table, binding_rates, arrhenius_check
"""
import functools

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import arrhenius
import curve_metrics
import scope

# A candidate worth testing further: the smoothed trough clears this many
# noise units below zero. Four, not the module's usual five or six, because
# the three-part screen below (sustained, no bubble overlap, survives
# debubble) is what does the real work of separating a genuine dip from
# noise -- see `trough_table`.
EARLY_TROUGH_CANDIDATE_Z = -4.0

# A candidate must still clear this after `curve_metrics.debubble` removes
# every detected O2 event, or the dip could be (part of) the gas artefact
# rather than the phenomenon this module is about. Set equal to the
# candidate threshold on purpose: nothing found in the archive needed a
# looser bar here, and two of the strongest candidates (141.4, 142.4) get
# MORE negative after correction, not less -- see DATA_VERIFICATION.md.
EARLY_TROUGH_SURVIVE_Z = -4.0

MINIMUM_READINGS = 15


def dominance(e0, s0, hoo):
    """
    Whichever reactant the catalyst's own fixed concentration comes closest
    to overwhelming: max([enz]/[S], [enz]/[HOO-]).

    NOT the minimum. [enz] is fixed and tiny against most of what this
    archive puts in a cuvette, so the reactant that matters is whichever one
    is scarce enough for [enz] to be a large FRACTION of it, and that is the
    larger of the two ratios, not the smaller. A curve with [enz]/[S] = 0.002
    and [enz]/[HOO-] = 4000 is dominated by the peroxide side, not the
    substrate side, and `max` is what selects that.
    """
    ratios = []
    if s0 > 0:
        ratios.append(e0 / s0)
    if hoo > 0:
        ratios.append(e0 / hoo)
    return max(ratios) if ratios else np.nan


def cluster(row):
    """
    Which starved reactant a genuine curve's dominance ratio points to.

    Purely descriptive of `row.e0_s0` against `row.e0_hoo` -- not a claim
    about which mechanism is correct, since the archive cannot separate a
    catalyst engaging the peroxide from one engaging the substrate on a
    single curve. Used to label the two clusters throughout, and computed
    for every scanned curve (not just the genuine ones) in `trough_table`.
    """
    if not (np.isfinite(row.e0_s0) or np.isfinite(row.e0_hoo)):
        return "unknown"
    return "substrate" if row.e0_s0 >= row.e0_hoo else "oxidant"


def _overlaps_event(events, start, stop):
    """Does any detected O2 event fall inside the reading range [start, stop)?"""
    return any(not (estop <= start or estart >= stop) for estart, estop in events)


def trough_table(block=None):
    """
    Every live catalysed curve's early trough, screened for three artefacts.

    Runs `curve_metrics.early_trough` over every live, catalysed curve in
    `block` (the whole archive by default -- this is a cross-archive
    statistic, not a single experimental design, so there is no narrower
    scope that would not simply throw curves away) and classifies each
    candidate (smoothed trough below `EARLY_TROUGH_CANDIDATE_Z`) against
    three independent ways it could be spurious rather than genuine:

      SUSTAINED   `early_trough`'s own leave-one-out test: drop the single
                  worst reading in the trough window and require the
                  remaining 8 to still average below its own threshold.
                  A first version instead required a fixed COUNT of the
                  window's readings to individually clear a per-reading
                  depth bar, and it wrongly rejected two real curves
                  (exps 4.1, 22.2) whose decline is genuine but shallower
                  per reading than the block's more dramatic examples --
                  9 of 9 readings negative in both, and removing the worst
                  one barely moves the mean (-7.3 to -6.1 sigma; -6.2 to
                  -5.1). See `curve_metrics.early_trough` and
                  DATA_VERIFICATION.md for the correction.

      NO BUBBLE   The trough window must not overlap a detected O2
                  detachment. A GROWING bubble raises absorbance (it
                  scatters light out of the detector), so it cannot produce
                  a dip by itself -- but a detachment is a sudden FALL, and
                  one landing inside the smoothing window could be read as
                  part of a decline that is actually the unrelated gas
                  artefact. None of the archive's candidates overlap one.

      SURVIVES    The trough recomputed on the `debubble`-corrected curve
      DEBUBBLE    must still clear `EARLY_TROUGH_SURVIVE_Z`. Confirms the
                  dip is not an O2 correction that happens to land in this
                  window even without a detected event inside it.

    `genuine` is all three at once: **17 of 17** candidates in the archive,
    with the corrected screen -- there are currently no known false
    positives left to exclude. Every column needed for the correlation and
    the write-up is here: `e0_s0`, `e0_hoo`, `e0_h2o2`, `dominance`,
    `cluster`.

    MEMOISED, and hands out a COPY -- the same convention `scope.frame` uses
    and for the same reason: an `lru_cache` handing out a shared DataFrame is
    one in-place edit from a silent wrong answer in whichever caller ran
    second.
    """
    return _trough_table(scope.archive() if block is None else block).copy()


@functools.cache
def _trough_table(scope_arg):
    frame = scope.frame(scope_arg)
    frame = frame[frame.live & (frame.e0 > 0)]
    lookup = {(c.experiment, c.sample): c for c in scope.curves(scope_arg)}

    rows = []
    for r in frame.itertuples():
        curve = lookup.get((r.experiment, r.sample))
        if curve is None or len(curve.times) < MINIMUM_READINGS:
            continue
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        z, t_trough, start, sustained = curve_metrics.early_trough(
            times, values, curve.noise)
        if not np.isfinite(z):
            continue

        stop = start + curve_metrics.EARLY_TROUGH_WINDOW
        events = curve_metrics.detachments(values, curve.noise)
        overlap = _overlaps_event(events, start, stop)

        corrected, _ = curve_metrics.debubble(times, values, curve.noise)
        z_corrected, _, _, _ = curve_metrics.early_trough(
            times, corrected, curve.noise)
        survives = (np.isfinite(z_corrected)
                   and z_corrected <= EARLY_TROUGH_SURVIVE_Z)

        candidate = z <= EARLY_TROUGH_CANDIDATE_Z
        genuine = candidate and sustained and not overlap and survives

        rows.append(dict(
            experiment=int(r.experiment), sample=int(r.sample),
            substrate=r.substrate, pH=float(r.pH), s0=float(r.s0),
            h2o2=float(r.h2o2), hoo=float(r.hoo), e0=float(r.e0),
            buffer=r.buffer, temperature=float(r.temperature),
            buf=float(r.buf),
            e0_s0=float(r.e0 / r.s0) if r.s0 > 0 else np.nan,
            e0_hoo=float(r.e0 / r.hoo) if r.hoo > 0 else np.nan,
            e0_h2o2=float(r.e0 / r.h2o2) if r.h2o2 > 0 else np.nan,
            dominance=dominance(r.e0, r.s0, r.hoo),
            z=z, t_trough=t_trough, candidate=candidate, sustained=sustained,
            bubble_overlap=overlap, z_corrected=z_corrected,
            survives=survives, genuine=genuine))
    table = pd.DataFrame(rows)
    table["cluster"] = table.apply(cluster, axis=1)
    return table


def dominance_correlation(table=None):
    """
    Spearman rank correlation of the trough depth against log-dominance,
    by substrate, over EVERY scanned curve (not just the genuine ones --
    this is the archive-wide test that the effect tracks the ratio rather
    than a hand-picked list of examples).

    `e0_h2o2` (the catalyst against TOTAL, un-weighted peroxide) is the
    control that distinguishes the anion from the neutral molecule: if the
    effect tracked peroxide in general rather than HOO- specifically, this
    column would correlate as well as `e0_hoo` does. It does not -- in
    4OMe-BnOH it is not even correctly signed (+0.160, not significant),
    against `e0_hoo`'s -0.621 at p = 5e-17.

    Returns {substrate: {"e0_s0": (rho, p, n), "e0_hoo": (rho, p, n),
    "e0_h2o2": (rho, p, n)}}.
    """
    table = trough_table() if table is None else table
    out = {}
    for substrate, group in table.groupby("substrate"):
        out[substrate] = {}
        for column in ("e0_s0", "e0_hoo", "e0_h2o2"):
            valid = group[column].notna() & (group[column] > 0)
            if valid.sum() < 10:
                out[substrate][column] = (np.nan, np.nan, int(valid.sum()))
                continue
            rho, p = spearmanr(np.log(group[column][valid]), group.z[valid])
            out[substrate][column] = (float(rho), float(p), int(valid.sum()))
    return out


def binding_rates(table=None):
    """
    A pseudo-first-order rate constant for each genuine curve's trough,
    from the SAME progress fit already computed for every curve in the
    archive -- no new fitting.

    A trough IS the "lag" shape of `summary_kinetics.fit_progress`'s
    one/two-phase functional form taken far enough that `B/tau` exceeds
    `v_ss`: `rate(0) = v_ss - B/tau < 0`. `tau_fast` (`scope.frame`'s own
    column -- `progress.two.tau1` for a two-phase curve, the one-phase
    burst's own `tau` otherwise) is therefore already the relaxation time
    of whatever produces the trough, and needs no re-fitting.

    Pseudo-first-order: `k_obs = 1/tau_fast = k_on * [excess reagent]`,
    where the excess reagent is whichever the curve's own `cluster` says
    the catalyst is closest to overwhelming -- `[enz]` for the oxidant
    cluster (ratio 95-4263 there, a clean excess) and `[S]` for the
    substrate cluster (ratio only ~10-15, a rougher approximation, since
    that is not a strong excess). `k_on` is returned in M-1 s-1.

    Returns `trough_table`'s genuine rows plus `k_obs` (1/s) and `k_on`
    (M-1 s-1).
    """
    table = trough_table() if table is None else table
    genuine = table[table.genuine].copy()
    frame = scope.frame(scope.archive()).set_index(["experiment", "sample"])
    k_obs, k_on = [], []
    for row in genuine.itertuples():
        fit_row = frame.loc[(row.experiment, row.sample)]
        obs = 1.0 / fit_row.tau_fast
        excess_mM = row.e0 if row.cluster == "oxidant" else row.s0
        k_obs.append(float(obs))
        k_on.append(float(obs / (excess_mM * 1e-3)))
    genuine["k_obs"] = k_obs
    genuine["k_on"] = k_on
    return genuine


def arrhenius_check(rates=None, cluster_name="oxidant"):
    """
    `arrhenius.arrhenius_fit` -- the same ln(rate)-against-1/T regression
    and the same honest standard error every temperature-series result in
    this project already uses -- applied to `k_on`, over every genuine
    curve in `cluster_name` and its own actual temperature. Not a two-point
    estimate: every available point, not just the two extremes.

    THE POINT IS THE SCATTER CHECK AS MUCH AS THE SLOPE. A two-point
    estimate (exp 19.1 at 15 C against exp 34.4 at 40 C, the obvious thing
    to compute by hand) gives `activation_kJ` = 85.7 -- but solving for the
    IMPLIED PRE-EXPONENTIAL FACTOR from those same two points gives
    2.1e15 M-1 s-1, about 2e5 times the diffusion limit, and the
    equivalent Eyring `dS_double_dagger` is +40 J/mol/K where a genuine
    bimolecular association (losing translational/rotational freedom) must
    be NEGATIVE. That two-point number is not a measurement.

    The proper fit, over all 14 oxidant-cluster curves and their own
    actual temperatures, gives `activation_kJ` = 78.9 +/- 48.6 -- not the
    same number tightened, but a demonstration that the slope is not
    resolved at all (`t_statistic` = 1.6, below the ~2 needed for
    significance). The reason is `same_temperature_spread`: among the 12
    curves that share the single most common temperature (298.15 K),
    `k_on` itself spans 1.43 orders of magnitude -- almost as wide as the
    1.74 spanned across every temperature from 15-40 C combined. Most of
    the apparent "trend" is ordinary curve-to-curve scatter a single point
    per temperature cannot distinguish from a real one.

    Adds `t_statistic`, `same_temperature_n` and `same_temperature_spread`
    (log10 units) to `arrhenius.arrhenius_fit`'s own returned dict
    (`activation_kJ`, `stderr_kJ`, `slope`, `slope_stderr`, `n`, `rms`, ...).
    """
    rates = binding_rates() if rates is None else rates
    block = rates[rates.cluster == cluster_name]
    kelvin = block.temperature.to_numpy() + 273.15
    fit = dict(arrhenius.arrhenius_fit(kelvin, block.k_on.to_numpy()))

    common_T = block.temperature.mode().iloc[0]
    same = block[block.temperature == common_T]
    fit["t_statistic"] = float(abs(fit["slope"] / fit["slope_stderr"]))
    fit["same_temperature_n"] = int(len(same))
    fit["same_temperature_spread"] = float(
        np.log10(same.k_on.max() / same.k_on.min()))
    return fit


def buffer_comparison(rates=None, cluster_name="oxidant"):
    """
    Geometric-mean `k_on` by buffer identity, over `cluster_name`.

    If the reaction driving the trough were a clean elementary step between
    the catalyst and free HOO-, buffer identity should not matter once
    `[enz]` is accounted for -- `k_on` should come out the same regardless
    of which buffer holds the pH. It does not: pyrophosphate's geometric
    mean runs roughly 6x phosphate's, even though both are normalised by
    `[enz]` identically. The more likely reading is that the buffer runs
    its own fast, cuvette-symmetric pre-equilibrium with HOO- (a buffer
    perhydrate, analogous to `buffer/ANALYSIS.md`'s own open question about
    what carries the oxidant) that sets the size of the pool the catalyst
    actually draws from -- a fast equilibrium present identically in both
    cuvettes cancels in the reference subtraction on its own, so it would
    not by itself explain the trough, but it would change how much
    reactive material is available for the catalyst's own, slower,
    genuinely asymmetric engagement to draw down.

    Returns a DataFrame: buffer, n, geometric_mean, min, max.
    """
    rates = binding_rates() if rates is None else rates
    block = rates[rates.cluster == cluster_name]
    rows = []
    for buf, group in block.groupby("buffer"):
        geo = float(np.exp(np.mean(np.log(group.k_on))))
        rows.append(dict(buffer=buf, n=len(group), geometric_mean=geo,
                          min=float(group.k_on.min()),
                          max=float(group.k_on.max())))
    return pd.DataFrame(rows).sort_values("buffer").reset_index(drop=True)
