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
from kinetic_model import Conditions
from kinetics_io import parse_experiment_data
from solution_chemistry import add_solution_columns

DATASET_PATH = "data/experiment_data.csv"
CURVE_DIRECTORY = "data/data"

# Absorbance is recorded to three decimals, so a reading's quantisation alone
# contributes this much standard deviation. Used as a floor on the per-curve
# noise estimate, exactly as in build_dossier.curve_flags -- a curve that sits
# on three or four distinct levels otherwise reports zero noise, and a zero
# denominator would give it infinite weight in the fit.
ABSORBANCE_QUANTUM = 0.001
QUANTISATION_SIGMA = ABSORBANCE_QUANTUM / np.sqrt(12)

# How many leading points the baseline is taken from. The model's signal is
# zero at t = 0 by construction, so the data has to be put on the same footing;
# a median over a few points rather than the first point alone keeps one noisy
# reading from shifting a whole curve.
BASELINE_POINTS = 5


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


def curve_noise(values):
    """
    Point-to-point noise from the median absolute second difference.

    The second difference annihilates any linear trend, so this measures a
    progress curve's scatter without being inflated by the progress itself.
    1.4826 converts a median absolute deviation to a standard deviation and the
    sqrt(6) undoes the variance the second difference introduces. Identical to
    build_dossier.curve_noise; duplicated rather than imported because that
    module pulls in matplotlib and builds a 6 MB page on import.
    """
    values = np.asarray(values, dtype=float)
    if len(values) < 5:
        return QUANTISATION_SIGMA
    curvature = values[2:] - 2 * values[1:-1] + values[:-2]
    estimate = 1.4826 * np.median(np.abs(curvature)) / np.sqrt(6)
    return max(float(estimate), QUANTISATION_SIGMA)


def select_fittable(data):
    """
    The rows a fit may use, by build_manifest's declared exclusions.

    Returns (selected, report) where report counts what each rule removed, so a
    caller can print why 454 rows became 404 instead of asserting it.
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


def read_all_curves(directory=CURVE_DIRECTORY):
    """
    Parses every instrument export once, returning {experiment: [(time, values)]}
    indexed by sample number.

    kinetics_io.load_experiment rescans the whole directory per experiment,
    which is fine for plotting one run and quadratic for fitting a hundred.
    Sample order is `parse_experiment_data`'s dict order, which is the order of
    the export's own header row -- the same mapping load_experiment uses, so
    sample numbers agree with the compiled dataset.
    """
    curves = {}
    for filename in sorted(os.listdir(directory)):
        if not filename.endswith(".txt"):
            continue
        parsed = parse_experiment_data(os.path.join(directory, filename))
        if parsed is None or parsed.get("num") is None:
            continue
        curves[parsed["num"]] = [
            (np.asarray(sample["time"], dtype=float),
             np.asarray(sample["values"], dtype=float))
            for sample in parsed["samples"].values()
        ]
    return curves


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
        times, values = samples[sample - 1]
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
            noise=curve_noise(values),
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
