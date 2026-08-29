"""
Consistent plotting for the kinetics dataset, built on kinetics_io.load_experiment.

Usage:
    python data/plot_kinetics.py 128
    python data/plot_kinetics.py 128 --mark 2:excluded --mark 3:excluded --mark 5:ref-channel
    python data/plot_kinetics.py 128 --save fig.png

Or from a notebook/script:
    from plot_kinetics import plot_experiment
    plot_experiment(128, mark_samples={2: "excluded", 3: "excluded", 5: "ref channel"})
"""
import argparse
import matplotlib.pyplot as plt

from kinetics_io import load_experiment

# Fixed color-per-sample-index convention so the same sample number always
# gets the same color across different plots/experiments.
_CMAP = plt.get_cmap("tab10")


def _sample_color(sample_index_zero_based):
    return _CMAP(sample_index_zero_based % 10)


def _sample_label(s):
    parts = [f"#{s['sample']}"]
    if s["[enz]"]:
        parts.append(f"[E]={s['[enz]']:.3g}")
    if s["[sub]"]:
        parts.append(f"[S]={s['[sub]']:.3g}")
    if s["[h2o2]"]:
        parts.append(f"[H2O2]={s['[h2o2]']:.3g}")
    if s["[buf]"]:
        parts.append(f"[buf]={s['[buf]']:.3g}")
    return " ".join(parts)


def plot_experiment(experiment_number, directory="data/data", mark_samples=None,
                     save_path=None, show=True):
    """
    Plots every sample in one experiment: raw absorbance vs time (left) and
    the Beer-Lambert-converted [P] vs time (right, when the substrate's
    extinction coefficient is known).

    mark_samples: optional {sample_number: annotation_text}. Matching samples
    are drawn dashed and get the annotation appended to their legend label --
    use this to flag samples known/suspected to be excluded from the fitting
    dataset (e.g. a reference channel, or a backwards/removed curve), without
    this module needing to know anything about clean_experiment_dataframe.
    """
    mark_samples = mark_samples or {}
    exp = load_experiment(experiment_number, directory=directory)
    if exp is None:
        raise ValueError(f"Could not load experiment {experiment_number}")

    fig, (ax_raw, ax_p) = plt.subplots(1, 2, figsize=(13, 5))

    for i, s in enumerate(exp["samples"]):
        color = _sample_color(i)
        marked = s["sample"] in mark_samples
        linestyle = "--" if marked else "-"
        label = _sample_label(s)
        if marked:
            label += f"  ({mark_samples[s['sample']]})"

        t = s["time"]
        v = s["values"]
        ax_raw.plot(t, v, linestyle, color=color, label=label, alpha=0.5 if marked else 1.0)

        if s["e"]:
            p = [x / s["e"] for x in v]
            ax_p.plot(t, p, linestyle, color=color, alpha=0.5 if marked else 1.0)

    ax_raw.set_xlabel("time (s)")
    ax_raw.set_ylabel("raw absorbance")
    ax_raw.set_title("Raw signal")
    ax_raw.legend(fontsize=7, loc="best")

    ax_p.set_xlabel("time (s)")
    ax_p.set_ylabel("[P] (mM)")
    ax_p.set_title("Beer-Lambert converted")

    title = (f"Experiment {exp['experiment']}  |  {exp['substrate']}, {exp['buffer']} buffer, "
             f"pH {exp['pH']}, T {exp['T']}\n{exp['txt_file']}  /  {exp['xls_file']}")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved to {save_path}")
    if show:
        plt.show()
    return fig


def _parse_mark_arg(values):
    marks = {}
    for v in values or []:
        sample_str, _, note = v.partition(":")
        marks[int(sample_str)] = note or "marked"
    return marks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot one experiment's kinetics curves.")
    parser.add_argument("experiment", type=int)
    parser.add_argument("--directory", default="data/data")
    parser.add_argument("--mark", action="append", help="sample:note, e.g. 5:ref-channel. Repeatable.")
    parser.add_argument("--save", default=None, help="Path to save the figure instead of / in addition to showing it.")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    plot_experiment(
        args.experiment,
        directory=args.directory,
        mark_samples=_parse_mark_arg(args.mark),
        save_path=args.save,
        show=not args.no_show,
    )
