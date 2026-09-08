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
independently), and the substrate itself, on the rare curves where the
substrate rung is low enough that the catalyst is a non-negligible fraction
of it. Neither needs anything not already established elsewhere in this
project: a real background reaction (`background_reaction/`), reference
subtraction (the instrument design `BUBBLES.md` already relies on for the
gas argument), and a catalyst that engages the peroxide non-productively
(`MECHANISM.md` S4, `BUBBLES.md`'s gas). This is that same engagement, read
off the OTHER cuvette's signal instead of the sample's own gas.

    from early_trough import trough_table, GENUINE, dominance
"""
import functools

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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

      SUSTAINED   `early_trough`'s own consecutive-reading test: at least
                  5 of the 9 readings inside the trough window are
                  individually below 2 sigma, not just the smoothed mean.
                  Rejects a single deep outlier reading masquerading as a
                  decline -- exps 4.1 and 22.2 both clear -6 sigma smoothed
                  and neither has more than 4 of 9.

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

    `genuine` is all three at once. Every column needed for the correlation
    and the write-up is here: `e0_s0`, `e0_hoo`, `dominance`.

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
            buffer=r.buffer,
            e0_s0=float(r.e0 / r.s0) if r.s0 > 0 else np.nan,
            e0_hoo=float(r.e0 / r.hoo) if r.hoo > 0 else np.nan,
            dominance=dominance(r.e0, r.s0, r.hoo),
            z=z, t_trough=t_trough, candidate=candidate, sustained=sustained,
            bubble_overlap=overlap, z_corrected=z_corrected,
            survives=survives, genuine=genuine))
    return pd.DataFrame(rows)


def dominance_correlation(table=None):
    """
    Spearman rank correlation of the trough depth against log-dominance,
    by substrate, over EVERY scanned curve (not just the genuine ones --
    this is the archive-wide test that the effect tracks the ratio rather
    than a hand-picked list of examples).

    Returns {substrate: {"e0_s0": (rho, p, n), "e0_hoo": (rho, p, n)}}.
    """
    table = trough_table() if table is None else table
    out = {}
    for substrate, group in table.groupby("substrate"):
        out[substrate] = {}
        for column in ("e0_s0", "e0_hoo"):
            valid = group[column].notna() & (group[column] > 0)
            if valid.sum() < 10:
                out[substrate][column] = (np.nan, np.nan, int(valid.sum()))
                continue
            rho, p = spearmanr(np.log(group[column][valid]), group.z[valid])
            out[substrate][column] = (float(rho), float(p), int(valid.sum()))
    return out


def cluster(row):
    """
    Which starved reactant a genuine curve's dominance ratio points to.

    Purely descriptive of `row.e0_s0` against `row.e0_hoo` -- not a claim
    about which mechanism is correct, since the archive cannot separate a
    catalyst engaging the peroxide from one engaging the substrate on a
    single curve. Used only to label the two clusters in the write-up.
    """
    if not (np.isfinite(row.e0_s0) or np.isfinite(row.e0_hoo)):
        return "unknown"
    return "substrate" if row.e0_s0 >= row.e0_hoo else "oxidant"
