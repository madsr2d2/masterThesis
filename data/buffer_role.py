"""
What the buffer does in this chemistry, and to what.

The buffer turns up in three places in this archive and has never been looked
at as one question. It is a REAGENT -- the uncatalysed reaction is first order
in it (`background_reaction/`), the catalysed turnover has an order in it that
saturates (`temperature_series/` 3), and the induction shortens with it
(`induction/` 6). It is also a CONFOUND -- substrate was added by volume and
displaced buffer, so `[buf]` falls 80 -> 50 mM as `[S]` rises in every 4OMe run
and a substrate order measured there is an order in the pair. And it is the
leading candidate for the step `induction/` is trying to name: general
acid/base catalysis of E -> E*.

WHAT A BUFFER TITRATION CAN AND CANNOT ASK. At a fixed pH the acid, the base
and the total are all proportional to each other -- the ratio is set by the pH
and the pKa and nothing else -- so an order in `[buf]` is an order in TOTAL
buffer and says nothing about which species does the work. Every number in this
module is such an order. The only escape is to titrate the buffer at more than
one pH and watch the CATALYTIC COEFFICIENT move: general base predicts it rises
with the base fraction, general acid that it falls with the acid fraction, and
the two predictions differ by a factor of three over the pH range available.

The archive has that design and it is not strong enough to use. See
`catalytic_coefficient`.

    python data/buffer_role.py

Names here are unique to this module; test_curve_metrics.test_no_duplicate_
definitions enforces that.
"""
import sys
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import arrhenius
import induction
import scope
import solution_chemistry


# Every run in the archive that steps [buf] with [S] held fixed inside it.
# Five, not the two that `temperature_series` 3 uses: 4OMe-BnOH at 40 C in
# phosphate, [H2O2] at 82.5 mM, [buf] over 3.125-200 mM, at THREE pH values --
# which is what makes the acid/base question askable at all.
TITRATIONS = (32, 34, 35, 36, 37)

# The pH split those five fall into. 7.20 is phosphate's second pKa, so the two
# groups sit either side of it and the base fraction moves from 0.387 to 0.681.
TITRATION_PH_SPLIT = 7.2

# The runs that share a buffer RANGE, so their catalytic coefficients are
# comparable without extrapolating one of them: 50-200 mM. Exp 34 covers
# 3.125-25 and is the only measurement of the low end, at one pH.
WIDE_RANGE = (32, 35, 36, 37)

# The substrate order the coefficient is normalised by before two runs at
# different [S] are compared. It is `induction.order_ratio`'s denominator on
# the catalysed 4OMe block, measured on the rate and not on anything here.
SUBSTRATE_ORDER = 0.471


def titration_table(experiments=TITRATIONS, parameter="v_peak"):
    """
    One row per buffer titration: its order in [buf], and its linear split.

    The order comes from `arrhenius.catalysed_buffer_order`, which is the
    package's own within-run log-log fit and is not re-derived here. The linear
    split is the classic buffer plot instead -- `v = a + b.[buf]`, so `a` is
    whatever the rate would be with no buffer and `b` is the catalytic
    coefficient. The two say different things and both are wanted: the order
    says how much of the rate the buffer carries, the coefficient says how
    fast the buffer-catalysed route runs.
    """
    import pandas as pd
    rows = []
    for experiment in experiments:
        block = scope.frame((experiment,))
        block = block[block.live & (block[parameter] > 0)].sort_values("buf")
        if len(block) < 3:
            continue
        order = arrhenius.catalysed_buffer_order((experiment,),
                                                 parameter=parameter)
        buffer_mM = block.buf.to_numpy(dtype=float)
        rate = block[parameter].to_numpy(dtype=float)
        design = np.column_stack([buffer_mM, np.ones(len(buffer_mM))])
        beta, *_ = np.linalg.lstsq(design, rate, rcond=None)
        resid = rate - design @ beta
        variance = float(resid @ resid) / max(1, len(rate) - 2)
        covariance = variance * np.linalg.pinv(design.T @ design)
        spread = float(((rate - rate.mean()) ** 2).sum())
        acid, base, pka = solution_chemistry.dominant_buffer_pair(
            str(block.buffer.iloc[0]), float(block.pH.iloc[0]), 100.0)
        rows.append({
            "experiment": int(experiment),
            "pH": float(block.pH.iloc[0]),
            "s0": float(block.s0.iloc[0]),
            "buffer_low": float(buffer_mM.min()),
            "buffer_high": float(buffer_mM.max()),
            "curves": int(len(block)),
            "order": order["order_buf"], "order_stderr": order["stderr_buf"],
            "order_r2": order["r2"],
            "coefficient": float(beta[0]),
            "coefficient_stderr": float(np.sqrt(max(covariance[0, 0], 0.0))),
            "intercept": float(beta[1]),
            "linear_r2": float(1 - float(resid @ resid) / spread)
            if spread else np.nan,
            "base_fraction": float(base / (acid + base)),
            "pka": float(pka),
        })
    return pd.DataFrame(rows)


def species_prediction(low_pH, high_pH, buffer_name="Phosphate", buffer_mM=100.0):
    """
    What each species would do to the catalytic coefficient between two pHs.

    If the buffer catalyses as its BASE, the coefficient per total buffer goes
    as the base fraction and so rises by `base(high)/base(low)`. If as its
    ACID, it falls by `acid(high)/acid(low)`. If the buffer is a spectator and
    the apparent order is something else -- ionic strength, say -- the ratio is
    1 and the design cannot tell that from either.
    """
    low_acid, low_base, pka = solution_chemistry.dominant_buffer_pair(
        buffer_name, low_pH, buffer_mM)
    high_acid, high_base, _ = solution_chemistry.dominant_buffer_pair(
        buffer_name, high_pH, buffer_mM)
    return {"pka": float(pka),
            "low_base_fraction": float(low_base / (low_acid + low_base)),
            "high_base_fraction": float(high_base / (high_acid + high_base)),
            "general_base": float(high_base / low_base),
            "general_acid": float(high_acid / low_acid),
            "spectator": 1.0}


def catalytic_coefficient(experiments=WIDE_RANGE, drop=(),
                          split=TITRATION_PH_SPLIT, parameter="v_peak",
                          substrate_order=SUBSTRATE_ORDER, frame=None):
    """
    Does the buffer's catalytic coefficient track the base, or the acid?

    Fits `v/[S]^n = a_run + b(pH group).[buf]` over the runs that share the
    50-200 mM range: one free level per run, and one slope for each side of the
    buffer's pKa. The ratio of the two slopes is the whole question.

    THE NORMALISATION IS NOT COSMETIC. The two pH groups do not sit at the same
    substrate -- 8.251 mM at pH 7.00 against 12.228 and 57.900 at pH 7.53 --
    and a free intercept per run absorbs that from the LEVEL but not from the
    SLOPE, because the buffer-catalysed route presumably carries the same
    substrate order the rest of the rate does. Dividing by `[S]^n` first is the
    assumption that it does, and it is an assumption: `report` prints the
    unnormalised fit beside it so the reader can see what it is worth.

    THE ANSWER IS THAT THE ARCHIVE CANNOT SAY. The ratio comes out near 1 with
    an error that spans both hypotheses, and the reason is visible in the
    inputs rather than only in the output: at pH 7.5 the buffer-independent
    route is eight times larger than at pH 7.00, so the buffer term is a much
    smaller share of a much bigger number and its coefficient is 44-127%
    uncertain there.
    """
    frame = scope.frame(tuple(experiments)) if frame is None else frame
    frame = frame[frame.live & (frame[parameter] > 0)]
    if drop:
        frame = frame[~frame.experiment.isin(drop)]
    runs = sorted(frame.experiment.unique())
    if len(runs) < 2:
        return {"curves": int(len(frame))}
    high = (frame.pH.to_numpy(dtype=float) > split)
    scale = frame.s0.to_numpy(dtype=float) ** substrate_order
    response = frame[parameter].to_numpy(dtype=float) / scale
    buffer_mM = frame.buf.to_numpy(dtype=float)
    columns = [buffer_mM * ~high, buffer_mM * high]
    columns += [(frame.experiment.to_numpy() == run) / scale for run in runs]
    design = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    resid = response - design @ beta
    rank = int(np.linalg.matrix_rank(design))
    variance = float(resid @ resid) / max(1, len(response) - rank)
    covariance = variance * np.linalg.pinv(design.T @ design)
    low, top = float(beta[0]), float(beta[1])
    low_error = float(np.sqrt(max(covariance[0, 0], 0.0)))
    top_error = float(np.sqrt(max(covariance[1, 1], 0.0)))
    ratio = top / low if low else np.nan
    # Delta method on a ratio, keeping the covariance: the two slopes share the
    # run levels, so treating them as independent would overstate the error.
    error = (abs(ratio) * np.sqrt((top_error / top) ** 2
                                  + (low_error / low) ** 2
                                  - 2 * covariance[0, 1] / (top * low))
             if top and low else np.nan)
    # The run LEVELS come back too. The design is `v/[S]^n = a_run/[S]^n +
    # b.[buf]`, so `a_run` is in raw rate units and is the run's rate
    # extrapolated to no buffer -- which is what `free_route_order` reads the
    # pH dependence of the buffer-INDEPENDENT route from.
    levels = {int(run): float(beta[2 + index])
              for index, run in enumerate(runs)}
    return {"curves": int(len(response)), "runs": len(runs),
            "low_pH": float(frame.pH[~high].max()) if (~high).any() else np.nan,
            "high_pH": float(frame.pH[high].min()) if high.any() else np.nan,
            "low": low, "low_stderr": low_error,
            "high": top, "high_stderr": top_error,
            "ratio": float(ratio), "ratio_stderr": float(error),
            "levels": levels,
            "normalised": substrate_order != 0.0}


# The pair of titrations that sit at the SAME substrate, so the level of one
# can be divided by the level of the other with nothing to normalise away:
# exp 32 at pH 7.00 and exp 35 at pH 7.50, both at 8.251 mM 4OMe-BnOH, both
# over 50-200 mM buffer. Exp 35's own SLOPE is unusable (R^2 0.176) and is
# dropped from `catalytic_coefficient`; its LEVEL is essentially the run's mean
# rate and survives a non-monotone slope, which is why the two uses differ.
MATCHED_LEVEL_PAIR = (32, 35)


def free_route_order(pair=MATCHED_LEVEL_PAIR, frame=None,
                     substrate_order=SUBSTRATE_ORDER):
    """
    How steeply does the buffer-INDEPENDENT route rise with pH?

    Section 2 of the folder's analysis says the buffer's ORDER collapses above
    the pKa not because the buffer term shrinks but because the rest of the
    rate grows, and attributes that growth to `[HOO-]`. That attribution is a
    quantitative claim and this checks it: `[HOO-]` is first order in `[OH-]`
    anywhere far below H2O2's own pKa, so the buffer-free level should rise by
    exactly `10^dpH` and no more.

    It rises by more. What that leaves is a route with a SECOND hydroxide in
    it, and the mechanistic reading is in the analysis: making a dioxirane out
    of plain H2O2 has to expel hydroxide from the Criegee adduct, and a step
    that poor is why the same chemistry elsewhere uses peroxymonosulfate or a
    peracid instead -- a peroxide that already carries a leaving group.

    Read on `MATCHED_LEVEL_PAIR` first, where the two runs sit at the same
    substrate and the ratio needs no normalisation at all, and then on every
    run through `catalytic_coefficient`'s levels, which do.
    """
    frame = scope.frame(TITRATIONS) if frame is None else frame
    rows = {}
    for experiment in pair:
        block = frame[(frame.experiment == experiment) & frame.live]
        if not len(block):
            return {}
        rows[experiment] = {
            "pH": float(block.pH.iloc[0]), "s0": float(block.s0.iloc[0]),
            "hoo": float(block.hoo.median()),
        }
    low, high = pair
    fit = catalytic_coefficient(experiments=TITRATIONS, frame=frame,
                                substrate_order=substrate_order)
    levels = fit["levels"]
    delta_pH = rows[high]["pH"] - rows[low]["pH"]
    level_ratio = levels[high] / levels[low] if levels.get(low) else np.nan
    hoo_ratio = rows[high]["hoo"] / rows[low]["hoo"]
    return {"low": low, "high": high,
            "low_pH": rows[low]["pH"], "high_pH": rows[high]["pH"],
            "delta_pH": float(delta_pH),
            "matched_s0": bool(np.isclose(rows[low]["s0"], rows[high]["s0"])),
            "s0": rows[low]["s0"],
            "level_low": levels.get(low, np.nan),
            "level_high": levels.get(high, np.nan),
            "level_ratio": float(level_ratio),
            "hoo_ratio": float(hoo_ratio),
            "hoo_order": float(np.log10(hoo_ratio) / delta_pH),
            "apparent_order": float(np.log10(level_ratio) / delta_pH)
            if level_ratio > 0 else np.nan,
            "levels": levels}


def peroxide_crossing(experiments=None, frame=None):
    """
    Does any run in this archive move `[buf]` and `[H2O2]` at once? None does.

    This is the design fact that decides how far the buffer question can be
    taken, and it is worth a function because the answer is not obvious from
    any table already here. A buffer acting as a general base is a term in
    `[buf]` alone. A buffer acting through a PEROXO ADDUCT OF ITSELF -- the
    catalyst oxidised by a buffer perhydrate rather than by H2O2 -- puts the
    two concentrations in the same term, because the adduct's concentration is
    set by the product `[buf][H2O2]`. The two schemes differ by an INTERACTION
    and by nothing else at a single pH, so separating them needs a run, or a
    pair of runs, that crosses the two ladders.

    Returns every run that steps `[buf]`, with the number of distinct `[H2O2]`
    it holds while doing so, and the same for the runs that step `[H2O2]`.
    """
    frame = scope.frame(tuple(experiments) if experiments
                        else tuple(range(1, 152))) if frame is None else frame
    live = frame[frame.live]
    grouped = live.groupby("experiment")
    rows = []
    for experiment, block in grouped:
        rows.append({"experiment": int(experiment),
                     "buffer_levels": int(block.buf.nunique()),
                     "peroxide_levels": int(block.h2o2.nunique())})
    import pandas as pd
    table = pd.DataFrame(rows)
    steps_buffer = table[table.buffer_levels > 1]
    steps_peroxide = table[table.peroxide_levels > 1]
    both = table[(table.buffer_levels > 1) & (table.peroxide_levels > 1)]
    # Across runs rather than within one: do the blocks that move the buffer
    # sit at more than one peroxide between them?
    buffer_runs = live[live.experiment.isin(steps_buffer.experiment)]
    titration_peroxides = sorted(
        float(value) for value in
        live[live.experiment.isin(TITRATIONS)].h2o2.unique())
    return {"runs": int(len(table)),
            "steps_buffer": int(len(steps_buffer)),
            "steps_peroxide": int(len(steps_peroxide)),
            "steps_both": int(len(both)),
            "buffer_run_peroxides": int(buffer_runs.h2o2.nunique()),
            "titration_peroxides": titration_peroxides}


def separable(coefficient, prediction, sigma=2.0):
    """
    Which of the three hypotheses the measured ratio excludes. Usually none.

    A hypothesis is excluded when its predicted ratio is more than `sigma`
    standard errors from the measured one. Reporting the survivors rather than
    a p-value is the point: the useful output of this design is the LIST of
    things still standing.
    """
    out = {}
    for name in ("general_base", "general_acid", "spectator"):
        gap = abs(coefficient["ratio"] - prediction[name])
        out[name] = {"predicted": float(prediction[name]),
                     "sigma": float(gap / coefficient["ratio_stderr"])
                     if coefficient["ratio_stderr"] else np.nan,
                     "excluded": bool(gap > sigma * coefficient["ratio_stderr"])}
    out["survivors"] = tuple(name for name in
                             ("general_base", "general_acid", "spectator")
                             if not out[name]["excluded"])
    return out


def identity_overlap(frame=None):
    """
    Can the buffer's IDENTITY be separated from the pH it was used at?

    Returns one row per (substrate, channel, buffer) with its pH range, so that
    "phosphate is slower than pyrophosphate" can be checked for whether the two
    were ever run at the same pH. Mostly they were not: the three buffers were
    each chosen for the pH they hold, which is what a buffer is for and what
    makes an identity comparison a pH comparison wearing a label.
    """
    import pandas as pd
    if frame is None:
        frame = scope.frame(tuple(range(1, 152)))
    live = frame[frame.live]
    rows = []
    for (substrate, differential, buffer_name), block in live.groupby(
            ["substrate", "differential", "buffer"]):
        rows.append({"substrate": substrate,
                     "channel": "catalysed" if differential else "enzyme-free",
                     "buffer": buffer_name,
                     "curves": int(len(block)),
                     "experiments": int(block.experiment.nunique()),
                     "pH_low": float(block.pH.min()),
                     "pH_high": float(block.pH.max()),
                     "buf_low": float(block.buf.min()),
                     "buf_high": float(block.buf.max()),
                     # Carried because an overlap in pH is only usable if the
                     # cells also share a peroxide: the 4OMe phosphate and
                     # pyrophosphate blocks overlap by two pH units and sit at
                     # 82.5 mM and 3.879/195.882 mM respectively.
                     "h2o2_low": float(block.h2o2.min()),
                     "h2o2_high": float(block.h2o2.max()),
                     # The VALUES, not the range. Ranges overlap where values
                     # never meet: the pyrophosphate block sits at exactly
                     # 3.879 and 195.882 mM and the phosphate block at exactly
                     # 82.5, so a range test calls those a match and they are
                     # not one.
                     "h2o2_levels": tuple(sorted(
                         round(float(v), 3) for v in block.h2o2.unique()))})
    return pd.DataFrame(rows)


def overlap_width(table):
    """
    The widest pH range two different buffers share inside one substrate/channel.

    Zero means the identity comparison cannot be made anywhere: whatever
    differs between two buffers here also differs in pH.
    """
    best = {"width": 0.0, "cell": None, "buffers": ()}
    for (substrate, channel), block in table.groupby(["substrate", "channel"]):
        if len(block) < 2:
            continue
        for first in range(len(block)):
            for second in range(first + 1, len(block)):
                one, two = block.iloc[first], block.iloc[second]
                low = max(one.pH_low, two.pH_low)
                high = min(one.pH_high, two.pH_high)
                if high - low > best["width"]:
                    shared = bool(set(one.h2o2_levels) & set(two.h2o2_levels))
                    best = {"width": float(high - low),
                            "cell": f"{substrate}, {channel}",
                            "buffers": (one.buffer, two.buffer),
                            "range": (float(low), float(high)),
                            "shares_peroxide": bool(shared),
                            "peroxide": (one.h2o2_levels, two.h2o2_levels)}
    return best


def report():
    """Print the whole argument, in the order it has to be made."""
    print("\n1. WHAT A TITRATION AT ONE pH CAN ASK")
    prediction = species_prediction(7.00, 7.53)
    print(f"   phosphate's pKa here is {prediction['pka']:.2f}, so the base "
          f"fraction is {prediction['low_base_fraction']:.3f} at pH 7.00 and "
          f"{prediction['high_base_fraction']:.3f} at 7.53.")
    print("   At ONE pH the acid, the base and the total are proportional, so")
    print("   an order in [buf] is an order in TOTAL buffer and nothing more.")
    print(f"   Across those two pHs the coefficient would rise "
          f"{prediction['general_base']:.2f}x for general base, "
          f"{prediction['general_acid']:.2f}x for general acid.")

    print("\n2. THE FIVE TITRATIONS, WHICH IS MORE THAN WAS BEING USED")
    table = titration_table()
    print(table.to_string(index=False))

    print("\n3. DOES THE COEFFICIENT TRACK THE BASE OR THE ACID")
    for label, drop, order in (("all four 50-200 mM runs", (), SUBSTRATE_ORDER),
                               ("without exp 35", (35,), SUBSTRATE_ORDER),
                               ("unnormalised", (), 0.0)):
        got = catalytic_coefficient(drop=drop, substrate_order=order)
        print(f"   {label:26s} b(low) {got['low']:.3g} +- {got['low_stderr']:.2g}"
              f"   b(high) {got['high']:.3g} +- {got['high_stderr']:.2g}"
              f"   ratio {got['ratio']:+.2f} +- {got['ratio_stderr']:.2f}")
    verdict = separable(catalytic_coefficient(drop=(35,)), prediction)
    for name in ("general_base", "general_acid", "spectator"):
        row = verdict[name]
        print(f"   {name:14s} predicts {row['predicted']:.2f}   "
              f"{row['sigma']:.1f} sigma away   "
              f"{'EXCLUDED' if row['excluded'] else 'survives'}")
    print(f"   survivors: {verdict['survivors']}")

    print("\n4. THE BUFFER AS A CONFOUND")
    frame = scope.frame(tuple(range(1, 152)))
    for name, block in induction.induction_blocks(frame).items():
        got = induction.composition_collinearity(block)
        if not got.get("runs"):
            continue
        print(f"   {name:26s} {got['runs']:2d} runs with a substrate ladder, "
              f"median corr(log S, log buf) {got['median']:+.2f}, "
              f"slope {got.get('slope', float('nan')):+.3f}".replace("+nan", "  --")
              + ", "
              f"{got['constant_buffer']} with [buf] constant")

    print("\n4a. IS THE BUFFER A BASE, OR IS IT CARRYING THE PEROXIDE")
    crossing = peroxide_crossing(frame=frame)
    print(f"   of {crossing['runs']} runs, {crossing['steps_buffer']} step "
          f"[buf] and {crossing['steps_peroxide']} step [H2O2];")
    print(f"   {crossing['steps_both']} step both, and the five titrations sit "
          f"at {crossing['titration_peroxides']} mM.")
    print("   A buffer perhydrate puts [buf] and [H2O2] in ONE term, so the two")
    print("   schemes differ by an interaction the archive never varies.")
    free = free_route_order(frame=scope.frame(TITRATIONS))
    print(f"   what the buffer-free level does instead, exps "
          f"{free['low']} -> {free['high']} at matched [S] = {free['s0']} mM:")
    print(f"   level {free['level_low']:.3g} -> {free['level_high']:.3g} over "
          f"dpH {free['delta_pH']:.2f}   = {free['level_ratio']:.2f}x")
    print(f"   [HOO-]                                        "
          f"= {free['hoo_ratio']:.2f}x   (order {free['hoo_order']:+.2f} "
          f"in [OH-], as it must be)")
    print(f"   so the level's apparent order in pH is "
          f"{free['apparent_order']:+.2f}, not {free['hoo_order']:+.2f}: the "
          f"buffer-free route")
    print(f"   carries {free['level_ratio'] / free['hoo_ratio']:.2f}x more "
          f"than one hydroperoxide accounts for.")

    print("\n5. BUFFER IDENTITY AGAINST pH")
    identity = identity_overlap(frame)
    print(identity.to_string(index=False))
    widest = overlap_width(identity)
    print(f"   widest pH range two buffers share in one cell: "
          f"{widest['width']:.2f} units  {widest['buffers']}  "
          f"{widest.get('range')}")
    print(f"   and do those two cells share a peroxide? "
          f"{widest.get('shares_peroxide')}   {widest.get('peroxide')}")
    return table


def main():
    report()


if __name__ == "__main__":
    main()
