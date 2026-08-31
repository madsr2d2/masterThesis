"""
Assembles the fitting dataset: metadata from experiment_data.csv, progress
curves from the instrument exports, and the [HOO-] each cuvette was run at.

Kept separate from `fit_kinetics.py` so that "which rows are fittable, and what
exactly is each curve" is answerable without running an optimiser, and so the
selection rule has one implementation rather than one per caller.

Row selection is `build_manifest`'s declared data, not a list retyped here:
KNOWN_EXCLUSIONS (whole experiments), KNOWN_SAMPLE_EXCLUSIONS (single
cuvettes) and EXCLUDED_BUFFERS (carbonate). That is the same selection the
notebook's clean_experiment_dataframe makes, and `test_fit_dataset.py` pins the
resulting counts so the two cannot drift apart silently.

    python data/fit_dataset.py            # summarise the fittable curves
    python data/fit_dataset.py --group    # break them down by fit group
"""
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

from build_manifest import (EXCLUDED_BUFFERS, KNOWN_EXCLUSIONS,
                            KNOWN_SAMPLE_EXCLUSIONS)
from curve_metrics import (ABSORBANCE_QUANTUM, QUANTISATION_SIGMA,
                           curve_noise)
from kinetic_model import Conditions
from kinetics_io import parse_experiment_data
from read_rre import ARCHIVE_DIR as RRE_DIRECTORY, RRE_SIGMA
from read_rre import read_all as read_all_rre
from solution_chemistry import add_solution_columns

DATASET_PATH = "data/experiment_data.csv"
CURVE_DIRECTORY = "data/data"

# Re-exported from curve_metrics so existing importers keep working. Curve
# shape is measured in one module; see curve_metrics for why.
ABSORBANCE_QUANTUM = ABSORBANCE_QUANTUM
QUANTISATION_SIGMA = QUANTISATION_SIGMA
curve_noise = curve_noise

# How many leading points the baseline is taken from. The model's signal is
# zero at t = 0 by construction, so the data has to be put on the same footing;
# a median over a few points rather than the first point alone keeps one noisy
# reading from shifting a whole curve.
BASELINE_POINTS = 5


# The block the fitting effort is scoped to, decided 2026-08-31.
#
# Exps 135-151 are the only runs in the archive that carry BOTH a substrate
# ladder and a peroxide ladder *inside a single run*: 100.0% of the scope's
# log[S] variance and 94.1% of its log[H2O2] variance is within-experiment, so
# neither order can be absorbed by a per-experiment offset. Every other design
# in the archive puts at most one axis inside the run, which is how the
# 4OMe-BnOH/40 C block came to rest its substrate order on 12.9%
# within-experiment contrast. They also span 19 pH values from 5.47 to 9.73 --
# 5.1 decades of [HOO-] -- in one (substrate, temperature, buffer) cell, and
# carry no exclusions, no accepted deviations and no open questions.
#
# This is NOT the same set as the hand-sorted data/Mads/'good data BnOH'
# folder, which also holds exp 50 (already excluded), exp 51 (a 4-cuvette
# borate run) and exp 134 (a sheet with no instrument export). Exp 50 survived
# earlier passes precisely because that folder looked authoritative, so the
# scope is defined by the design and re-derived by test_fit_kinetics, never by
# where a file was filed.
#
# Exps 75 and 76 share the block's (BnOH, 25 C, Pyrophosphate) key but are
# excluded from the scope: they carry the unresolved hexametaphosphate
# speciation question (DATA_VERIFICATION.md 2026-08-31).
#
# One condition must be met before a fit on this scope is quotable, recorded in
# FITTING.md: the high-peroxide cuvettes of exps 135 and 138 need a ruling,
# since they backtrack up to 0.35 AU and that lands straight in the residuals.
#
# The two-salt buffer (Na4P2O7 with Na2HPO4, 143-151 NaH2PO4, equimolar) is
# still treated as pure pyrophosphate and is worth correcting -- it moves 6 of
# the 17 runs from above Davies' 500 mM ceiling to below it, and makes the
# unrecorded titrant recoverable by electroneutrality -- but it is not a
# blocker: across this scope Davies' pKa_eff varies only 11.478 to 11.494, so
# the entire ionic-strength apparatus is worth under 3.2% in [HOO-] against a
# 129000x span driven by pH and [H2O2] alone.
PRIMARY_SCOPE = frozenset(range(135, 152))
PRIMARY_SCOPE_BLOCK = ("BnOH", 25.0, "Pyrophosphate")


def in_scope(curves, scope=PRIMARY_SCOPE):
    """The curves of `curves` whose experiment is in `scope`."""
    return [curve for curve in curves if curve.experiment in scope]


@dataclass
class Curve:
    """One cuvette: what was in it, and what the instrument recorded."""
    experiment: int
    sample: int
    substrate: str
    buffer: str
    pH: float
    temperature: float
    epsilon: float          # mM^-1 cm^-1, the substrate's extinction coefficient
    times: np.ndarray       # s, starting at 0
    absorbance: np.ndarray  # baseline-subtracted, as recorded
    baseline: float         # what was subtracted
    noise: float            # absorbance units, 1 sigma
    conditions: Conditions
    source: str = "txt"     # "rre" where the instrument file was read

    @property
    def group(self):
        """
        Rate constants are only shared where they can be. Temperature changes
        every constant through Arrhenius and the two substrates are different
        molecules, so a fit may pool curves only within one (substrate,
        temperature, buffer) cell.
        """
        return (self.substrate, self.temperature, self.buffer)

    def __len__(self):
        return len(self.times)


def select_fittable(data):
    """
    The rows a fit may use, by build_manifest's declared exclusions.

    Returns (selected, report) where report counts what each rule removed, so a
    caller can print why 454 rows became 402 instead of asserting it.
    """
    report = {"rows_in": len(data), "experiments_in": data.experiment.nunique()}

    selected = data[~data.experiment.isin(KNOWN_EXCLUSIONS)]
    report["excluded_experiments"] = report["rows_in"] - len(selected)

    pairs = list(KNOWN_SAMPLE_EXCLUSIONS)
    mask = np.zeros(len(selected), dtype=bool)
    for experiment, sample in pairs:
        mask |= ((selected.experiment == experiment) & (selected["sample"] == sample)).to_numpy()
    before = len(selected)
    selected = selected[~mask]
    report["excluded_samples"] = before - len(selected)

    before = len(selected)
    selected = selected[~selected.buffer.isin(EXCLUDED_BUFFERS)]
    report["excluded_buffers"] = before - len(selected)

    report["rows_out"] = len(selected)
    report["experiments_out"] = selected.experiment.nunique()
    return selected, report


def read_all_curves(directory=CURVE_DIRECTORY, rre_directory=RRE_DIRECTORY):
    """
    Every curve once, as {experiment: [(time, values, source)]} in sample order.

    Values come from the instrument's own .rre where one exists and from the
    .txt export otherwise. The .rre is the same measurement at ~1000x the
    resolution -- the .txt is rounded to 0.001 AU, which zeroes the measured
    noise on 67 of the 119 in-scope curves -- so it is preferred wherever it
    is available and agrees. See rre_io.

    A .rre is used only when its sample has the same number of points as the
    export and tracks it to within the export's own rounding step. Anything
    else falls back to the .txt and is counted, never silently substituted:
    the two files are different formats written by different code paths and
    a misalignment would be invisible in the result.

    kinetics_io.load_experiment rescans the whole directory per experiment,
    which is fine for plotting one run and quadratic for fitting a hundred.
    Sample order is `parse_experiment_data`'s dict order, which is the order of
    the export's own header row -- the same mapping load_experiment uses, so
    sample numbers agree with the compiled dataset.
    """
    high_precision = read_all_rre(rre_directory)
    curves = {}
    for filename in sorted(os.listdir(directory)):
        if not filename.endswith(".txt"):
            continue
        parsed = parse_experiment_data(os.path.join(directory, filename))
        if parsed is None or parsed.get("num") is None:
            continue
        number = parsed["num"]
        better = high_precision.get(number, {})
        series = []
        for index, sample in enumerate(parsed["samples"].values(), start=1):
            times = np.asarray(sample["time"], dtype=float)
            values = np.asarray(sample["values"], dtype=float)
            series.append((times, *_prefer_rre(values, better.get(index))))
        curves[number] = series
    return curves


def _prefer_rre(exported, instrument):
    """(values, source) -- the .rre when it agrees with the export."""
    if instrument is None or len(instrument) != len(exported):
        return exported, "txt"
    # Both are compared as changes from their own first point: the export is
    # rounded absorbance and the .rre is -log10(%T/100) against a baseline the
    # instrument stored, so their offsets need not agree, only their shapes.
    drift = np.abs((instrument - instrument[0]) - (exported - exported[0])).max()
    if drift > ABSORBANCE_QUANTUM:
        return exported, "txt"
    return instrument, "rre"


def build_curves(dataset_path=DATASET_PATH, directory=CURVE_DIRECTORY,
                 minimum_points=10):
    """
    Every fittable cuvette, as Curve objects.

    Returns (curves, report). Rows whose export is missing, whose curve is too
    short to constrain anything, or whose extinction coefficient is unknown are
    dropped and counted in the report rather than silently skipped.
    """
    data = pd.read_csv(dataset_path)
    selected, report = select_fittable(data)
    selected = add_solution_columns(selected)
    exports = read_all_curves(directory)

    curves, dropped = [], {"no_export": 0, "no_curve": 0, "too_short": 0, "no_epsilon": 0}
    # Indexed by column name rather than by itertuples position: 'sample'
    # collides with a namedtuple method and would have to be reached
    # positionally, which breaks silently the moment a column is added.
    for row in selected.to_dict("records"):
        experiment = int(row["experiment"])
        sample = int(row["sample"])
        if experiment not in exports:
            dropped["no_export"] += 1
            continue
        samples = exports[experiment]
        if not 1 <= sample <= len(samples):
            dropped["no_curve"] += 1
            continue
        times, values, source = samples[sample - 1]
        if len(times) < minimum_points:
            dropped["too_short"] += 1
            continue
        epsilon = float(row["e"] or 0.0)
        if epsilon <= 0:
            dropped["no_epsilon"] += 1
            continue

        baseline = float(np.median(values[:max(1, min(BASELINE_POINTS, len(values) // 10))]))
        curves.append(Curve(
            experiment=experiment,
            sample=sample,
            substrate=str(row["substrate"]),
            buffer=str(row["buffer"]),
            pH=float(row["pH"]),
            temperature=float(row["T"]),
            epsilon=epsilon,
            times=times - times[0],
            absorbance=values - baseline,
            baseline=baseline,
            # The floor belongs to the source: a .rre curve floored at the
            # export's 0.001 AU quantisation would report 2.4x its real noise.
            noise=curve_noise(values,
                              RRE_SIGMA if source == "rre"
                              else QUANTISATION_SIGMA),
            source=source,
            conditions=Conditions(
                s0=float(row["[sub]"]),
                h2o2=float(row["[h2o2]"]),
                e0=float(row["[enz]"]),
                hoo=float(row["[HOO-]"]),
            ),
        ))
    report["dropped"] = dropped
    report["curves"] = len(curves)
    return curves, report


def group_curves(curves, enzyme_free=None):
    """
    Buckets curves by (substrate, temperature, buffer).

    enzyme_free=True keeps only E0 = 0 curves, False only E0 > 0, None keeps
    both. The two stages of the sequential fit are exactly these two subsets.
    """
    grouped = {}
    for curve in curves:
        if enzyme_free is True and curve.conditions.e0 != 0:
            continue
        if enzyme_free is False and curve.conditions.e0 == 0:
            continue
        grouped.setdefault(curve.group, []).append(curve)
    return grouped


def _summarise(curves, title):
    print(f"\n{title}")
    print(f"  {'substrate':11s} {'T':>4s} {'buffer':11s} {'curves':>6s} {'exps':>5s} "
          f"{'pH':>13s} {'[S] mM':>15s} {'points':>7s}")
    for key in sorted(group_curves(curves), key=str):
        block = group_curves(curves)[key]
        substrate, temperature, buffer_name = key
        pH = [c.pH for c in block]
        substrate_conc = [c.conditions.s0 for c in block]
        print(f"  {substrate:11s} {temperature:4.0f} {buffer_name:11s} {len(block):6d} "
              f"{len({c.experiment for c in block}):5d} "
              f"{min(pH):5.2f}-{max(pH):5.2f}  {min(substrate_conc):6.2f}-{max(substrate_conc):6.2f}  "
              f"{sum(len(c) for c in block):7d}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Summarise the fitting dataset.")
    parser.add_argument("--dataset", default=DATASET_PATH)
    parser.add_argument("--directory", default=CURVE_DIRECTORY)
    arguments = parser.parse_args()

    built, built_report = build_curves(arguments.dataset, arguments.directory)
    print(f"{built_report['rows_in']} rows / {built_report['experiments_in']} experiments compiled")
    print(f"  -{built_report['excluded_experiments']:3d} rows  excluded experiments "
          f"({len(KNOWN_EXCLUSIONS)} of them)")
    print(f"  -{built_report['excluded_samples']:3d} rows  excluded samples "
          f"({len(KNOWN_SAMPLE_EXCLUSIONS)} of them)")
    print(f"  -{built_report['excluded_buffers']:3d} rows  excluded buffers "
          f"({', '.join(sorted(EXCLUDED_BUFFERS))})")
    print(f"  ={built_report['rows_out']:4d} rows / {built_report['experiments_out']} "
          f"experiments fittable")
    if any(built_report["dropped"].values()):
        print(f"  dropped when loading curves: {built_report['dropped']}")
    print(f"  {built_report['curves']} curves carrying "
          f"{sum(len(c) for c in built):,} points")

    _summarise([c for c in built if c.conditions.e0 == 0],
               "ENZYME-FREE (stage 1: k_can, k3, k0, r)")
    _summarise([c for c in built if c.conditions.e0 > 0],
               "CATALYSED (stage 2: k5, k6)")
