"""
Which curves may be used for what, and which look broken.

Two separate questions, kept apart because conflating them is how real
chemistry gets deleted.

ELIGIBILITY is about measurement power. A curve whose fitting window climbs
less than a few instrument quanta cannot constrain a rate -- not because the
cuvette failed, but because absorbance is reported to three decimals. That is
true of experiment 25's dead sample 2 (0.8 quanta) and equally true of the
legitimate bottom rung of a titration (0.2 quanta). Both are useless for a
substrate order; only one is broken. So eligibility attaches to a USE, never to
the curve:

    rate    orders, Km, anything fitted to v0    -- needs a measurable slope
    shape   lag fraction, burst amplitude, tau   -- needs amplitude, not slope

A curve can be ineligible for the first and perfectly good for the second.

The rate cut is applied PER BLOCK and only where it demonstrably moves no
answer. `validate_power_cut` refits each group with and without the low-power
curves; if a fitted order shifts by more than a standard error the cut is
selecting on the outcome rather than on power, and `unsafe_blocks` says so.
BnOH/25/Pyrophosphate fails exactly that way.

DEFECTS are about the cuvette failing, and are screened in three layers with
very different reach:

    1  condition-free   truncated runs, physical impossibilities        100%
    2  ladder           a rung far below or above BOTH its neighbours    68%
    3  peer             disagreement with cuvettes at matched conditions  6%

Layer 2 NOMINATES and never convicts. In an alternating ladder it flags the
survivors too: experiment 25's sample 3 is a sound cuvette that reads as a 10x
spike purely because both its neighbours are dead. Only layer 3 assigns blame,
and it reaches 6% of the data -- there is no replication behind the rest, which
is a fact about the experimental design and not something a rule can fix.

Nothing here writes an exclusion. Convictions go into
`build_manifest.KNOWN_SAMPLE_EXCLUSIONS` by hand, with their evidence.

CURVE SHAPE IS NEVER A DEFECT. Initial dips, lags and bursts are reported and
never screened out: at 49-60 s sampling every real kinetic feature is slow and
smooth, so there is no timescale separating chemistry from artefact for shape
to exploit, and the striking shapes in this dataset are live hypotheses. A dip
may be substrate sequestered into the cyclodextrin cavity before turnover
begins; see DATA_VERIFICATION.md.

    python data/curve_screen.py                 # the screen, all groups
    python data/curve_screen.py --candidates    # only what two layers agree on
    python data/curve_screen.py --validate      # does the power cut move any order?
"""
import argparse
from dataclasses import dataclass

import numpy as np
import pandas as pd

from fit_dataset import ABSORBANCE_QUANTUM, build_curves
from summary_kinetics import (INITIAL_WINDOW, REPLICATE_RSD, buffer_concentrations,
                              profile_km, regress, summarise, to_frame)

# A window that climbs fewer than this many 0.001 AU steps cannot support a
# slope. Three is where the criterion stops being sensitive: on the best-
# designed block it removes one curve of 57 and moves no fitted order at all.
MINIMUM_WINDOW_QUANTA = 3.0

# Below this many points there is no curve to speak of. Deliberately far below
# the 20 that `build_dossier` uses as a defect: experiment 26, the only true
# replicate set in the enzyme-free data, has ten-point curves and is sound.
MINIMUM_POINTS = 8

# A curve needs to have moved this many times its own noise before its SHAPE is
# worth discussing, whatever its rate does.
MINIMUM_SHAPE_SIGMA = 10.0

# Physical impossibilities. Neither fires on the present dataset; they are here
# as cheap guards against a future data-entry accident, not as live filters.
MAXIMUM_ABSORBANCE = 2.0      # beyond the detector's linear range
MAXIMUM_CONVERSION = 1.05     # more product than there was substrate

# A ladder rung this far from BOTH its neighbours is inconsistent with any
# smooth dependence on the conditions. Flat between 2.5x and 4x, so this is a
# stable choice rather than a tuned one.
LADDER_FACTOR = 3.0

# How close two cuvettes must be in each axis to count as peers. [HOO-] and
# [buf] get loose tolerances because their orders are absorbed by the
# per-experiment offsets and so cannot be corrected for.
PEER_TOLERANCE = {"s0": 0.02, "h2o2": 0.02, "hoo": 0.25, "buf": 0.25, "e0": 0.10}
MINIMUM_PEERS = 3
PEER_Z = 4.0

# How far a fitted order may move when the power cut is applied before the cut
# is judged unsafe for that block, in pooled standard errors. A power criterion
# that changes an answer is not measuring power -- it is selecting on the
# outcome. BnOH/25/Pyrophosphate fails this: the cut removes 12 of the 16
# curves below the 10th percentile of [sub] there, truncating the bottom of the
# substrate ladder, and the order falls from +0.41 to +0.11. The same cut on
# 4OMe-BnOH/40/Phosphate removes one curve of 61, at HIGH [sub], and moves
# nothing. So the cut is applied per block, never globally.
POWER_SHIFT_LIMIT = 1.0


@dataclass(frozen=True)
class Eligibility:
    """What one curve may be used for, and why not otherwise."""
    rate: bool
    shape: bool
    reasons: tuple = ()

    @property
    def summary(self):
        parts = [use for use, ok in (("rate", self.rate), ("shape", self.shape)) if ok]
        return "+".join(parts) if parts else "neither"


def eligibility(row):
    """
    Per-use eligibility for one curve summary (a CurveSummary or a mapping).

    Never consults other curves, so it cannot be confounded by a bad neighbour
    or a missing peer set, and applies to 100% of the data.
    """
    get = row.get if hasattr(row, "get") else (lambda k, d=None: getattr(row, k, d))
    rate_ok, shape_ok, reasons = True, True, []

    points = get("points", 0) or 0
    if points < MINIMUM_POINTS:
        rate_ok = shape_ok = False
        reasons.append(("both", f"only {points} points"))

    amplitude, noise = get("amplitude", np.nan), get("noise", np.nan)
    if np.isfinite(amplitude) and np.isfinite(noise) and noise > 0:
        if amplitude > MAXIMUM_ABSORBANCE:
            rate_ok = shape_ok = False
            reasons.append(("both", f"absorbance {amplitude:.2f} AU is past the "
                                    f"detector's linear range"))
        if amplitude < MINIMUM_SHAPE_SIGMA * noise:
            shape_ok = False
            reasons.append(("shape", f"moved only {amplitude / noise:.0f}x its noise, "
                                     f"so there is no shape to read"))
    conversion = get("conversion", np.nan)
    if np.isfinite(conversion) and conversion > MAXIMUM_CONVERSION:
        rate_ok = shape_ok = False
        reasons.append(("both", f"conversion {conversion:.0%} exceeds the substrate "
                                f"present"))

    v0 = get("v0", np.nan)
    if not np.isfinite(v0):
        rate_ok = False
        reasons.append(("rate", "no initial slope could be fitted"))
    quanta = get("window_quanta", np.nan)
    if np.isfinite(quanta) and quanta < MINIMUM_WINDOW_QUANTA:
        rate_ok = False
        reasons.append(("rate", f"the fitting window rises {quanta:.1f} quanta, "
                                f"below the {MINIMUM_WINDOW_QUANTA:.0f} needed to "
                                f"measure a slope at 0.001 AU resolution"))
    return Eligibility(rate=rate_ok, shape=shape_ok, reasons=tuple(reasons))


def add_eligibility(frame):
    """Attach `eligible_rate`, `eligible_shape` and `eligibility_reasons`."""
    flags = [eligibility(row) for row in frame.to_dict("records")]
    frame = frame.copy()
    frame["eligible_rate"] = [f.rate for f in flags]
    frame["eligible_shape"] = [f.shape for f in flags]
    frame["eligibility"] = [f.summary for f in flags]
    frame["eligibility_reasons"] = [f.reasons for f in flags]
    return frame


def ladder_anomalies(frame, factor=LADDER_FACTOR):
    """
    Rungs that sit far below or far above BOTH their neighbours.

    Symmetric on purpose: a rung above both neighbours is as impossible in a
    monotone titration as one below, and has its own causes -- over-pipetted
    substrate, a bubble, contamination.

    NOMINATES ONLY. Where two adjacent rungs fail, the sound cuvette between
    them reads as a spike; experiment 25 sample 3 does exactly that. A hit here
    means "this ladder is internally inconsistent", never "this cuvette is
    broken".
    """
    hits = []
    for experiment, block in frame.groupby("experiment"):
        if len(block) < 3:
            continue
        axis = next((k for k in ("s0", "hoo", "buf")
                     if block[k].nunique() == len(block)), None)
        if axis is None:
            continue
        block = block.sort_values(axis).reset_index(drop=True)
        for i in range(len(block)):
            here = block.v0[i]
            neighbours = [block.v0[j] for j in (i - 1, i + 1)
                          if 0 <= j < len(block)
                          and np.isfinite(block.v0[j]) and block.v0[j] > 0]
            if len(neighbours) < 2 or not np.isfinite(here) or here <= 0:
                continue
            below = min(n / here for n in neighbours)
            above = min(here / n for n in neighbours)
            if below > factor:
                kind, size = "dip", below
            elif above > factor:
                kind, size = "spike", above
            else:
                continue
            hits.append(dict(experiment=int(experiment),
                             sample=int(block["sample"][i]),
                             kind=kind, factor=float(size), axis=axis,
                             v0=float(here),
                             neighbours=float(np.median(neighbours))))
    return hits


def peer_sets(frame, tolerance=None):
    """{row label: [labels of cuvettes at matched conditions]} within a group."""
    tolerance = tolerance or PEER_TOLERANCE

    def close(a, b, tol):
        if a == 0 and b == 0:
            return True
        if a == 0 or b == 0 or not (np.isfinite(a) and np.isfinite(b)):
            return False
        return abs(np.log(a / b)) <= np.log(1 + tol)

    peers = {}
    for _, block in frame.groupby("group"):
        records = block.to_dict("index")
        for label, row in records.items():
            peers[label] = [other for other, q in records.items()
                            if other != label
                            and all(close(row[k], q[k], t) for k, t in tolerance.items())]
    return peers


def peer_scores(frame, tolerance=None, floor=REPLICATE_RSD, minimum=MINIMUM_PEERS):
    """
    z of each curve against cuvettes at matched conditions elsewhere.

    Uses no fitted orders, so nothing an absorbed axis can corrupt, and no
    per-experiment offset, so an experiment cannot hide its own failures behind
    one. The price is reach: most curves have no peer set and come back NaN,
    which is reported as unjudgeable rather than as passing.
    """
    peers = peer_sets(frame, tolerance)
    floor_log = np.log10(1 + floor)
    z, counts, centres = {}, {}, {}
    for label, group in peers.items():
        rates = frame.loc[group, "v0"] if group else pd.Series(dtype=float)
        rates = rates[np.isfinite(rates) & (rates > 0)]
        here = frame.at[label, "v0"]
        counts[label] = len(rates)
        if len(rates) < minimum or not np.isfinite(here) or here <= 0:
            z[label], centres[label] = np.nan, np.nan
            continue
        logs = np.log10(rates.to_numpy())
        centre = float(np.median(logs))
        scale = max(1.4826 * float(np.median(np.abs(logs - centre))), floor_log)
        z[label] = (np.log10(here) - centre) / scale
        centres[label] = 10 ** centre
    return (pd.Series(z, name="peer_z"), pd.Series(counts, name="n_peers"),
            pd.Series(centres, name="peer_v0"))


def screen(frame, factor=LADDER_FACTOR, peer_threshold=PEER_Z):
    """
    Run all three layers plus eligibility. Returns the frame with the verdicts.

    `layers` counts only the layers that can convict -- 1 and 3. Layer 2 is
    recorded in `ladder` and deliberately does not count, because it cannot
    tell a failed rung from a sound one between two failures.
    """
    frame = add_eligibility(frame)
    z, counts, centres = peer_scores(frame)
    frame["peer_z"], frame["n_peers"], frame["peer_v0"] = z, counts, centres

    frame["ladder"] = ""
    frame["ladder_factor"] = np.nan
    for hit in ladder_anomalies(frame, factor):
        mask = ((frame.experiment == hit["experiment"])
                & (frame["sample"] == hit["sample"]))
        frame.loc[mask, "ladder"] = hit["kind"]
        frame.loc[mask, "ladder_factor"] = hit["factor"]

    layer1 = ~frame.eligible_rate & ~frame.eligible_shape      # condition-free fault
    layer3 = frame.peer_z.notna() & (frame.peer_z < -peer_threshold)
    frame["layer1"], frame["layer3"] = layer1, layer3
    frame["layers"] = layer1.astype(int) + layer3.astype(int)
    return frame


def validate_power_cut(frame, quanta=MINIMUM_WINDOW_QUANTA):
    """
    Does dropping low-power curves move any fitted order, or only the scatter?

    The check that has to pass before the cut is applied to a block. A power
    criterion that changes an answer is not measuring power.
    """
    report = []
    for name, block in frame.groupby("group"):
        kept = block[block.window_quanta >= quanta]
        row = dict(group=name, n=len(block), dropped=len(block) - len(kept))
        for label, data in (("before", block), ("after", kept)):
            try:
                fitted = regress(data, "v0", None, per_experiment=True)
                order, stderr = fitted.coefficient("log[S]")
                row[f"order_{label}"] = order
                row[f"stderr_{label}"] = stderr
                row[f"scatter_{label}"] = fitted.residual_scatter
                row[f"km_{label}"] = profile_km(data, "v0", True).km
            except (ValueError, np.linalg.LinAlgError):
                row[f"order_{label}"] = np.nan
                row[f"stderr_{label}"] = np.nan
                row[f"scatter_{label}"] = np.nan
                row[f"km_{label}"] = np.nan
        shift = row.get("order_after", np.nan) - row.get("order_before", np.nan)
        pooled = np.hypot(row.get("stderr_before", np.nan), row.get("stderr_after", np.nan))
        row["shift_sigma"] = shift / pooled if pooled and np.isfinite(pooled) else np.nan
        report.append(row)
    return pd.DataFrame(report)


def unsafe_blocks(frame, quanta=MINIMUM_WINDOW_QUANTA, limit=POWER_SHIFT_LIMIT):
    """
    Groups where the power cut moves a fitted order rather than only the noise.

    Filtering on `eligible_rate` inside one of these biases the answer, because
    the low-power curves there are the low-substrate rungs and dropping them
    truncates the ladder. Callers must consult this before applying the cut.
    """
    report = validate_power_cut(frame, quanta)
    moved = report[report.shift_sigma.abs() > limit]
    return set(moved.group), report


def report_eligibility(frame):
    print("\nELIGIBILITY -- what each curve may be used for")
    print(f"  {'group':38s} {'n':>4s} {'rate':>6s} {'shape':>6s} {'neither':>8s}")
    for name, block in frame.groupby("group"):
        print(f"  {name:38s} {len(block):4d} "
              f"{block.eligible_rate.sum():6d} {block.eligible_shape.sum():6d} "
              f"{(~block.eligible_rate & ~block.eligible_shape).sum():8d}")
    print(f"  {'TOTAL':38s} {len(frame):4d} "
          f"{frame.eligible_rate.sum():6d} {frame.eligible_shape.sum():6d} "
          f"{(~frame.eligible_rate & ~frame.eligible_shape).sum():8d}")
    only_shape = frame[~frame.eligible_rate & frame.eligible_shape]
    print(f"\n  {len(only_shape)} curves carry no measurable rate but a readable shape.")
    print("  Those are not broken -- they are the slow rungs, and they are exactly")
    print("  the curves a lag or burst study still needs.")

    unsafe, report = unsafe_blocks(frame)
    checked = report.shift_sigma.notna().sum()
    print(f"\n  POWER CUT SAFETY: {checked}/{len(report)} blocks could be checked "
          f"(the rest have too few curves to regress).")
    if unsafe:
        print("  DO NOT filter on eligible_rate in these blocks -- the cut moves the")
        print("  answer, not just the noise:")
        for name in sorted(unsafe):
            row = report[report.group == name].iloc[0]
            print(f"    {name}: order {row.order_before:+.3f} -> {row.order_after:+.3f} "
                  f"({row.shift_sigma:+.1f} sigma), {row.dropped}/{row.n} dropped")
    else:
        print("  Every checkable block is safe.")


def report_screen(frame, peer_threshold=PEER_Z):
    print("\nDEFECT SCREEN")
    print(f"  layer 1 (condition-free fault)      {frame.layer1.sum():4d} curves")
    print(f"  layer 2 (ladder dip or spike)       "
          f"{(frame.ladder != '').sum():4d} curves  -- nominates only")
    judged = frame.peer_z.notna()
    print(f"  layer 3 (peer disagreement)         {frame.layer3.sum():4d} curves"
          f"   [{judged.sum()}/{len(frame)} judgeable, {100 * judged.mean():.0f}%]")
    both = frame[frame.layers >= 2]
    print(f"\n  CONVICTED by two convicting layers: {len(both)}")
    for row in both.itertuples():
        print(f"    exp {row.experiment} sample {row.sample}  "
              f"peer z={row.peer_z:+.1f} on {int(row.n_peers)} peers")
    print("\n  NOMINATED by the ladder (needs a second layer, or your eye):")
    for row in frame[frame.ladder != ""].itertuples():
        agree = "  <-- layer 3 agrees" if row.layer3 else ""
        z = f"{row.peer_z:+.1f}" if np.isfinite(row.peer_z) else "no peers"
        print(f"    exp {row.experiment} sample {row.sample}  {row.ladder:5s} "
              f"{row.ladder_factor:5.1f}x vs both neighbours   peer z={z}{agree}")
    lonely = frame[frame.layer3 & (frame.layers < 2)]
    if len(lonely):
        print("\n  Peer disagreement alone (look, do not exclude):")
        for row in lonely.itertuples():
            print(f"    exp {row.experiment} sample {row.sample}  z={row.peer_z:+.1f} "
                  f"on {int(row.n_peers)} peers, v0={row.v0:.2e} vs {row.peer_v0:.2e}")


def report_validation(report):
    print("\nPOWER CUT VALIDATION -- does dropping low-power curves move an order?")
    print(f"  {'group':38s} {'drop':>9s} {'order before':>13s} {'order after':>12s} "
          f"{'shift':>7s} {'scatter':>15s}")
    for row in report.itertuples():
        shift = f"{row.shift_sigma:+.2f}s" if np.isfinite(row.shift_sigma) else "   --"
        before = (f"{row.order_before:+.3f}" if np.isfinite(row.order_before) else "  --")
        after = (f"{row.order_after:+.3f}" if np.isfinite(row.order_after) else "  --")
        scatter = (f"{row.scatter_before:.2f}x -> {row.scatter_after:.2f}x"
                   if np.isfinite(row.scatter_before) and np.isfinite(row.scatter_after)
                   else "      --")
        print(f"  {row.group:38s} {row.dropped:3d}/{row.n:<5d} {before:>13s} {after:>12s} "
              f"{shift:>7s} {scatter:>15s}")
    moved = report[report.shift_sigma.abs() > 1.0]
    if len(moved):
        print(f"\n  {len(moved)} block(s) move by more than 1 sigma. The cut is NOT safe")
        print("  there and must not be applied to them without a reason:")
        for row in moved.itertuples():
            print(f"    {row.group}: {row.order_before:+.3f} -> {row.order_after:+.3f} "
                  f"({row.shift_sigma:+.1f} sigma)")
    else:
        print("\n  No block moves by more than 1 sigma: the cut removes noise, not answer.")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--enzyme-free", action="store_true")
    parser.add_argument("--catalysed", action="store_true")
    parser.add_argument("--candidates", action="store_true",
                        help="only what two convicting layers agree on")
    parser.add_argument("--validate", action="store_true",
                        help="check the power cut moves no fitted order")
    parser.add_argument("--quanta", type=float, default=MINIMUM_WINDOW_QUANTA)
    parser.add_argument("--save", default=None)
    arguments = parser.parse_args()

    curves, _ = build_curves()
    if arguments.enzyme_free:
        curves = [c for c in curves if c.conditions.e0 == 0]
    elif arguments.catalysed:
        curves = [c for c in curves if c.conditions.e0 > 0]
    frame = to_frame(summarise(curves, buffer_concentrations()))
    frame = screen(frame)

    print(f"{len(frame)} curves, {frame.experiment.nunique()} experiments, "
          f"{frame.group.nunique()} group(s); window {INITIAL_WINDOW:.0%}, "
          f"power floor {arguments.quanta:.0f} quanta of {ABSORBANCE_QUANTUM} AU")
    if arguments.candidates:
        report_screen(frame)
    elif arguments.validate:
        report_validation(validate_power_cut(frame, arguments.quanta))
    else:
        report_eligibility(frame)
        report_screen(frame)
    if arguments.save:
        frame.drop(columns=["eligibility_reasons"]).to_csv(arguments.save, index=False)
        print(f"\nscreen -> {arguments.save}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
