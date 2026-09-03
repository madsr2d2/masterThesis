"""
Per-curve parameters against [S] and against [H2O2], one panel each.

Every point is one cuvette. Points are joined into LADDERS -- cuvettes of the
same experiment that share a value of the other axis -- because a ladder is the
only comparison in this dataset that varies one concentration and nothing else.
Two cuvettes from different runs differ in pH, [HOO-], enzyme batch, cell and
day as well, and a slope read across them is a correlation rather than an
order. The scope was chosen so that these ladders exist: 100.0% of its log[S]
variance and 94.1% of its log[H2O2] variance is within-experiment.

The slope quoted on each panel is the common order from
`scope.orders(..., within=True)`, which fits every ladder at once with a free
offset per experiment. The pooled figure is printed beside it; where the two
disagree, the difference is the between-run confounding.

v0 and vmax are both plotted deliberately. v0 is the rate before any catalyst
has built up -- on an accelerating curve that is the induction period, not the
reaction -- and vmax is the rate after. They do not carry the same order.

    python data/plot_orders.py                      # -> figures/
    python data/plot_orders.py --outdir /tmp/x
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from curve_metrics import ACCELERATION_SIGMA
from scope import (AXIS_PARTNER, TWO_AXIS_BLOCK, frame, ladder_groups,
                   ladder_trend, orders)

# label, column, log y?  -- the parameters worth an order.
PANELS = [
    ("initial rate $v_0$ (AU/s)", "v0", True),
    ("developed rate $v_{max}$ (AU/s)", "vmax", True),
    ("total change (AU)", "net", True),
    ("speed-up $v_{max}/v_0$", "gain", True),
    ("acceleration ($\\sigma$)", "accel_z", False),
]
AXES = [("s0", "[BnOH] (mM)"), ("h2o2", "[H$_2$O$_2$] (mM)")]


def _panel(axis_key, axis_label, parameter, parameter_label, log_y, ax, data):
    """One panel: every ladder in `axis_key`, coloured by its run's pH."""
    norm = plt.Normalize(data.pH.min(), data.pH.max())
    cmap = plt.get_cmap("viridis")
    for label, group in ladder_groups(axis_key):
        y = group[parameter].to_numpy(dtype=float)
        if log_y:
            # A ladder with a non-positive rung cannot be drawn on a log axis.
            # Dropping the rung keeps the rest of the ladder rather than the
            # whole ladder being silently absent from the figure.
            keep = y > 0
            if keep.sum() < 2:
                continue
            group, y = group[keep], y[keep]
        ax.plot(group[axis_key], y, "-o", ms=4, lw=1.1, alpha=0.85,
                color=cmap(norm(group.pH.iloc[0])))
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    else:
        ax.axhline(0, color="#888", lw=0.8)
        ax.axhline(ACCELERATION_SIGMA, color="#c0522a", lw=0.8, ls="--")
    ax.set_xlabel(axis_label)
    ax.set_ylabel(parameter_label)
    ax.grid(alpha=0.25, which="both")

    if log_y:
        within = orders(parameter, within=True)
        pooled = orders(parameter, within=False)
        key = "order_s0" if axis_key == "s0" else "order_h2o2"
        err = "stderr_s0" if axis_key == "s0" else "stderr_h2o2"
        ax.set_title(f"order {within[key]:+.2f} $\\pm$ {within[err]:.2f}"
                     f"   (pooled {pooled[key]:+.2f})", fontsize=10)
    else:
        low, high, count = ladder_trend(parameter, axis_key)
        ax.set_title(f"median {low:+.1f}$\\sigma$ at the bottom rung, "
                     f"{high:+.1f}$\\sigma$ at the top ({count} ladders)",
                     fontsize=10)
    return norm, cmap


def build(outdir="figures", scope=TWO_AXIS_BLOCK):
    """Write the figure. Returns its path."""
    os.makedirs(outdir, exist_ok=True)
    data = frame(scope)
    data = data[data.live]

    fig, axes = plt.subplots(len(PANELS), 2, figsize=(12, 4.4 * len(PANELS)))
    # The panels carry their own x labels, so they need room between rows;
    # the default leaves a label sitting on the title beneath it.
    fig.subplots_adjust(hspace=0.42, top=0.94)
    for row, (parameter_label, parameter, log_y) in zip(axes, PANELS):
        for ax, (axis_key, axis_label) in zip(row, AXES):
            norm, cmap = _panel(axis_key, axis_label, parameter,
                                parameter_label, log_y, ax, data)

    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    bar = fig.colorbar(mappable, ax=axes, fraction=0.025, pad=0.02)
    bar.set_label("pH of the run")
    fig.suptitle("In-scope per-cuvette parameters against each concentration.\n"
                 "A line joins cuvettes of one run that differ in that axis "
                 "alone; the quoted order is the common within-experiment slope.",
                 fontsize=12)
    path = os.path.join(outdir, "scope_orders.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="figures")
    arguments = parser.parse_args()
    print(f"wrote {build(arguments.outdir)}")
