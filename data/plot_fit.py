"""
Plots a fit against the curves it was fitted to.

A residual number says a fit is bad; a plot says how. This draws three figures
from a results file written by `fit_kinetics.py --save`:

    <block>_stage1.png     enzyme-free curves, data vs model, one panel per experiment
    <block>_stage2.png     the same for the catalysed curves
    <block>_diagnostics.png  parity, residuals, shape, and the hidden species

Nothing is refitted -- the constants come from the results file -- so this is
seconds rather than the half-hour a fit takes.

    python data/fit_kinetics.py --substrate BnOH --save fit.json
    python data/plot_fit.py fit.json

Model curves are put through `fit_kinetics.baseline_like_data`, the same
transformation the measured curves get, because otherwise the picture would show
a disagreement the fit was never asked to remove.
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

from fit_dataset import DATASET_PATH, build_curves
from fit_kinetics import baseline_like_data
from kinetic_model import Conditions, RateConstants, observable, simulate

_CMAP = plt.get_cmap("tab10")

# Drawn on the diagnostics figure as the band a fit would have to reach to be
# believable: the model within a few times each curve's own reading noise.
CREDIBLE_SIGMA = 3.0


def _colour(index):
    return _CMAP(index % 10)


def _constants(record):
    return RateConstants(**{name: record["constants"][name]
                            for name in ("k_can", "k3", "k0", "k5", "k6", "r")})


def _model(curve, constants):
    """The modelled absorbance for one curve, on the data's own footing."""
    signal = observable(constants, curve.conditions, curve.times)
    if signal is None:
        return None
    return baseline_like_data(curve.epsilon * signal)


def _label(curve):
    parts = [f"[S]={curve.conditions.s0:.2f}"]
    if curve.conditions.e0:
        parts.append(f"[E]={curve.conditions.e0:.3f}")
    return " ".join(parts)


def plot_stage(curves, constants, title, path):
    """One panel per experiment; every cuvette in it, data as points, model as a line."""
    experiments = sorted({curve.experiment for curve in curves})
    columns = min(3, len(experiments))
    rows = int(np.ceil(len(experiments) / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(5.2 * columns, 3.9 * rows),
                                squeeze=False)

    for position, experiment in enumerate(experiments):
        axis = axes[position // columns][position % columns]
        block = [c for c in curves if c.experiment == experiment]
        for index, curve in enumerate(sorted(block, key=lambda c: c.sample)):
            colour = _colour(index)
            axis.plot(curve.times / 60.0, curve.absorbance, ".", color=colour,
                      markersize=3, alpha=0.55, label=_label(curve))
            model = _model(curve, constants)
            if model is not None:
                axis.plot(curve.times / 60.0, model, "-", color=colour, linewidth=1.6)
        first = block[0]
        axis.set_title(f"exp {experiment}   pH {first.pH:g}, "
                       f"[H2O2] {first.conditions.h2o2:.0f} mM", fontsize=9)
        axis.set_xlabel("time (min)", fontsize=8)
        axis.set_ylabel("absorbance (baseline-subtracted)", fontsize=8)
        axis.tick_params(labelsize=7)
        axis.legend(fontsize=6, loc="upper left")
        axis.axhline(0.0, color="0.8", linewidth=0.6, zorder=0)

    for position in range(len(experiments), rows * columns):
        axes[position // columns][position % columns].axis("off")

    figure.suptitle(title + "      points = measured, lines = model", fontsize=11)
    figure.tight_layout(rect=[0, 0, 1, 0.96])
    figure.savefig(path, dpi=140)
    plt.close(figure)
    return path


def _peak_position(values, times, smooth_fraction=0.05):
    """
    Where the steepest point sits, as a fraction of the run, smoothed over ~5%
    of the run first -- the same method MECHANISM.md uses, so the numbers on this
    figure are comparable to the ones in that document.
    """
    values = np.asarray(values, dtype=float)
    window = max(3, int(len(values) * smooth_fraction) | 1)
    if len(values) < window + 4:
        return np.nan
    smoothed = np.convolve(values, np.ones(window) / window, mode="valid")
    centres = times[window // 2:window // 2 + len(smoothed)]
    slope = np.gradient(smoothed, centres)
    if slope[0] <= 0 or slope.max() <= 1.05 * slope[0]:
        return 0.0
    return float((centres[np.argmax(slope)] - centres[0]) / (centres[-1] - centres[0]))


def initial_rate(times, values, fraction=0.2):
    """Slope over the first `fraction` of a run, in absorbance per second."""
    count = max(5, int(len(times) * fraction))
    return float(np.polyfit(times[:count], values[:count], 1)[0])


def substrate_orders(curves, constants, buffer_by_sample):
    """
    Effective reaction order in [S] per experiment, for data and for model.

    Experiments whose [buf] changes along the cuvette ladder are reported
    separately and never pooled: in those series buffer and substrate were
    varied together, so a slope against [S] is partly a buffer effect. That
    collinearity is a standing finding in DATA_VERIFICATION.md, and exps 3 and 6
    here are exactly it -- they give a NEGATIVE apparent substrate order.
    """
    clean, confounded = [], []
    for experiment in sorted({c.experiment for c in curves}):
        block = [c for c in curves if c.experiment == experiment]
        if len(block) < 3:
            continue
        substrate = np.array([c.conditions.s0 for c in block])
        buffers = np.array([buffer_by_sample.get((c.experiment, c.sample), np.nan)
                            for c in block])
        measured = np.array([initial_rate(c.times, c.absorbance) for c in block])
        modelled = []
        for curve in block:
            model = _model(curve, constants)
            modelled.append(np.nan if model is None
                            else initial_rate(curve.times, model))
        modelled = np.array(modelled)
        varies = (np.nanmax(buffers) - np.nanmin(buffers)) / np.nanmax(buffers) > 0.05
        (confounded if varies else clean).append(
            (experiment, substrate, measured, modelled))
    return clean, confounded


def _log_slope(x, y):
    keep = (np.asarray(x) > 0) & (np.asarray(y) > 0) & np.isfinite(y)
    if keep.sum() < 3:
        return np.nan
    return float(np.polyfit(np.log(np.asarray(x)[keep]),
                            np.log(np.asarray(y)[keep]), 1)[0])


def plot_diagnostics(stages, title, path):
    """
    Six panels that say how the fit fails rather than by how much:

      a  net signal, measured against modelled, with the 1:1 line
      b  residuals in units of each curve's own noise, against time
      c  where the steepest point sits, measured against modelled
      d  the species the model is proposing, for one representative curve
      e  initial rate against [S] -- the model is first order, the data is not
      f  the same as an order per experiment, with the buffer titrations flagged
    """
    figure, axes = plt.subplots(2, 3, figsize=(18.0, 9.5))
    (parity, residual, shape), (species, orders, per_experiment) = axes

    for stage_index, (name, curves, constants) in enumerate(stages):
        colour = _colour(stage_index)
        nets_measured, nets_modelled, positions = [], [], []
        for curve in curves:
            model = _model(curve, constants)
            if model is None:
                continue
            nets_measured.append(curve.absorbance[-1])
            nets_modelled.append(model[-1])
            positions.append((_peak_position(curve.absorbance, curve.times),
                              _peak_position(model, curve.times)))
            residual.plot(curve.times / 60.0,
                          (model - curve.absorbance) / curve.noise,
                          "-", color=colour, linewidth=0.7, alpha=0.5)
        parity.plot(nets_measured, nets_modelled, "o", color=colour,
                    markersize=5, alpha=0.75, label=name)
        measured, modelled = zip(*positions) if positions else ((), ())
        shape.plot(np.array(measured) * 100, np.array(modelled) * 100, "o",
                   color=colour, markersize=5, alpha=0.75, label=name)

    limit = max(parity.get_xlim()[1], parity.get_ylim()[1])
    parity.plot([0, limit], [0, limit], "k--", linewidth=1, label="1:1")
    parity.set_xlabel("measured net signal (AU)")
    parity.set_ylabel("modelled net signal (AU)")
    parity.set_title("a. the model reproduces the wrong amount of signal", fontsize=10)
    parity.legend(fontsize=8)

    residual.axhspan(-CREDIBLE_SIGMA, CREDIBLE_SIGMA, color="0.85", zorder=0,
                     label=f"+/-{CREDIBLE_SIGMA:g} sigma")
    residual.axhline(0.0, color="0.4", linewidth=0.8)
    residual.set_xlabel("time (min)")
    residual.set_ylabel("(model - data) / that curve's noise")
    residual.set_title("b. residuals, in units of each curve's own noise", fontsize=10)
    residual.legend(fontsize=8)

    shape.plot([0, 100], [0, 100], "k--", linewidth=1)
    shape.set_xlabel("measured position of peak slope (% into run)")
    shape.set_ylabel("modelled position (% into run)")
    shape.set_title("c. the model lags where the data does not", fontsize=10)
    shape.legend(fontsize=8)

    # The species panel: what the model is actually claiming happens, for the
    # median-substrate enzyme-free curve. This is where a reader can see that
    # the aldehyde pool is a trace and benzoate carries the signal.
    name, curves, constants = stages[0]
    chosen = sorted(curves, key=lambda c: c.conditions.s0)[len(curves) // 2]
    trajectory = simulate(constants, chosen.conditions, chosen.times)
    if trajectory is not None:
        minutes = chosen.times / 60.0
        for key, style in (("A", "-"), ("PBA", "--"), ("BA", "-.")):
            species.plot(minutes, trajectory[key], style, linewidth=1.6, label=f"[{key}]")
        species.plot(minutes, chosen.absorbance / chosen.epsilon, ".", color="0.35",
                     markersize=3, label="measured / eps")
        species.set_yscale("log")
        species.set_xlabel("time (min)")
        species.set_ylabel("concentration (mM), log scale")
        species.set_title(f"d. what the model proposes is happening "
                          f"(exp {chosen.experiment} sample {chosen.sample}, "
                          f"[S]={chosen.conditions.s0:.2f} mM)", fontsize=10)
        species.legend(fontsize=8)

    # Panels e and f: the reaction order in substrate. This is the sharpest
    # single statement of the mismatch, so it gets two panels -- one showing the
    # raw rate-vs-[S] scaling, one reducing each experiment to a number.
    buffer_by_sample = {(int(row["experiment"]), int(row["sample"])): float(row["[buf]"])
                        for row in pd.read_csv(DATASET_PATH).to_dict("records")}
    data_slopes, model_slopes, flagged = [], [], []
    for stage_index, (name, curves, constants) in enumerate(stages):
        colour = _colour(stage_index)
        clean, confounded = substrate_orders(curves, constants, buffer_by_sample)
        for experiment, substrate, measured, modelled in clean:
            orders.loglog(substrate, np.maximum(measured, 1e-12), "o",
                          color=colour, markersize=5, alpha=0.8)
            orders.loglog(substrate, np.maximum(modelled, 1e-12), "x",
                          color=colour, markersize=6, alpha=0.8)
            data_slopes.append((experiment, _log_slope(substrate, measured), colour))
            model_slopes.append(_log_slope(substrate, modelled))
        for experiment, substrate, measured, modelled in confounded:
            flagged.append((experiment, _log_slope(substrate, measured), colour))

    orders.set_xlabel("[S] (mM)")
    orders.set_ylabel("initial rate (AU/s)")
    orders.set_title("e. circles = measured, crosses = model\n"
                     "(only experiments with [buf] held constant)", fontsize=10)

    positions = np.arange(len(data_slopes))
    for position, (experiment, slope, colour) in zip(positions, data_slopes):
        per_experiment.bar(position, slope, color=colour, alpha=0.85)
        per_experiment.text(position, slope + 0.03, str(experiment),
                            ha="center", fontsize=7)
    for offset, (experiment, slope, colour) in enumerate(flagged):
        position = len(data_slopes) + offset
        per_experiment.bar(position, slope, color=colour, alpha=0.35, hatch="//")
        per_experiment.text(position, slope + 0.03, str(experiment),
                            ha="center", fontsize=7)
    clean_mean = np.nanmean([slope for _, slope, _ in data_slopes])
    per_experiment.axhline(1.0, color="k", linestyle="--", linewidth=1.2,
                           label=f"model: {np.nanmean(model_slopes):.2f}")
    per_experiment.axhline(clean_mean, color="crimson", linewidth=1.2,
                           label=f"data mean: {clean_mean:.2f}")
    per_experiment.axhline(0.0, color="0.6", linewidth=0.8)
    per_experiment.set_ylabel("effective order in [S]")
    per_experiment.set_xticks([])
    per_experiment.set_title("f. order per experiment (hatched = [buf] varies too,\n"
                             "so substrate and buffer are confounded)", fontsize=10)
    per_experiment.legend(fontsize=8)

    figure.suptitle(title, fontsize=12)
    figure.tight_layout(rect=[0, 0, 1, 0.95])
    figure.savefig(path, dpi=140)
    plt.close(figure)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("results", help="a JSON file from fit_kinetics.py --save")
    parser.add_argument("--outdir", default="figures")
    arguments = parser.parse_args()

    with open(arguments.results) as handle:
        results = json.load(handle)
    os.makedirs(arguments.outdir, exist_ok=True)

    substrate, temperature, buffer_name = [part.strip()
                                           for part in results["block"].split(",")]
    key = (substrate, float(temperature.rstrip(" C")), buffer_name)
    slug = f"{substrate}_{temperature.replace(' ', '')}_{buffer_name}".replace("/", "-")

    curves, _ = build_curves()
    block = [curve for curve in curves if curve.group == key]
    if not block:
        print(f"no curves found for block {key}")
        return 1

    stages, written = [], []
    for stage, label, enzyme_free in (("stage_1", "STAGE 1 enzyme-free", True),
                                      ("stage_2", "STAGE 2 catalysed", False)):
        if stage not in results:
            continue
        constants = _constants(results[stage])
        subset = [c for c in block if (c.conditions.e0 == 0) == enzyme_free]
        stages.append((label, subset, constants))
        written.append(plot_stage(
            subset, constants,
            f"{results['block']}   {label}   "
            f"({results[stage]['rms_sigma']:.0f}x the curves' own noise)",
            os.path.join(arguments.outdir, f"{slug}_{stage}.png")))

    if stages:
        written.append(plot_diagnostics(
            stages, f"{results['block']} -- how the fit fails",
            os.path.join(arguments.outdir, f"{slug}_diagnostics.png")))

    for path in written:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
