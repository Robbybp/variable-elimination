#  ___________________________________________________________________________
#
#  Variable Elimination: Research code for variable elimination in NLPs
#
#  Copyright (c) 2023. Triad National Security, LLC. All rights reserved.
#
#  This software is distributed under the 3-clause BSD license.
#  ___________________________________________________________________________

import itertools
import os

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

import var_elim.scripts.config as config

CMAP = ListedColormap([(0.95, 0.55, 0.55), (0.05, 0.05, 0.35)])
MODEL_NAMES = ["mb-steady", "distill", "pipeline"]
PARAMETER_LOOKUP = {
    "mb-steady": (
        "fs.moving_bed.solid_phase.properties[0,1].temperature",
        "fs.moving_bed.solid_phase.properties[0,1].flow_mass",
    ),
    "distill": ("vol", "x_Feed"),
    "pipeline": (
        "fs.nodes[0].state[*].temperature",
        "fs.nodes[0].state[*].pressure",
    ),
}
PARAMETER_LABEL_LOOKUP = {
    "mb-steady": ("Temperature (K)", "Flow rate (kg/s)"),
    "distill": ("Relative volatility", "Feed mole fraction"),
    "pipeline": ("Gas temperature (K)", "Supply pressure (bar)"),
}
ELIM_NAMES = ("no-elim", "d1", "ecd2", "linear-d2", "d2", "greedy", "matching")
METHOD_SUBSETS = [
    ("No elimination", ("no-elim",)),
    ("Virtual best: No elimination + Linear degree-2", ("no-elim", "linear-d2")),
    ("Virtual best: No elimination + Linear degree-2 + Greedy", ("no-elim", "linear-d2", "greedy")),
    ("Virtual best: All algorithms", ELIM_NAMES),
]
TITLE_LOOKUP = {
    "mb-steady": "Moving bed reactor",
    "distill": "Distillation",
    "pipeline": "Pipeline",
}


def _success_series(df):
    success = df
    if success.dtype == bool:
        return success
    return success.astype(str).str.lower().isin(("true", "1"))


def load_virtual_best(results_dir, model_name, method_names, suffix):
    parameter_names = PARAMETER_LOOKUP[model_name]
    suffix_str = "" if suffix is None else f"-{suffix}"
    successes = []
    for method_name in method_names:
        filename = f"{model_name}-{method_name}-sweep{suffix_str}.csv"
        fpath = os.path.join(results_dir, filename)
        df = pd.read_csv(fpath)
        required_columns = set(parameter_names) | {"success"}
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            raise ValueError(f"{fpath} is missing columns: {sorted(missing_columns)}")
        successes.append(
            df.set_index(list(parameter_names))["success"].pipe(_success_series).rename(method_name)
        )

    success_df = pd.concat(successes, axis=1, join="outer")
    if success_df.isna().any().any():
        raise ValueError(
            f"Sweep results for {model_name} do not contain the same parameter instances"
        )
    success_df["success"] = success_df.any(axis=1)
    return success_df.reset_index()[list(parameter_names) + ["success"]]


def plot_convergence(ax, df, parameter_names, parameter_labels, title):
    parameters = [list(sorted(set(df[name]))) for name in parameter_names]
    param_index_maps = [{p: i for i, p in enumerate(params)} for params in parameters]
    success_lookup = {
        tuple(row[name] for name in parameter_names): int(row["success"])
        for _, row in df.iterrows()
    }
    n_success = sum(success_lookup.values())
    print(f"{title}: converged {n_success} / {len(df)} instances")

    convergence_array = np.zeros(tuple(len(params) for params in parameters))
    for params in itertools.product(*parameters):
        indices = tuple(idx_map[p] for idx_map, p in zip(param_index_maps, params))
        convergence_array[indices] = success_lookup[params]

    ax.imshow(
        convergence_array,
        aspect="equal",
        origin="lower",
        cmap=CMAP,
        vmin=0,
        vmax=1,
    )

    x_ticks = list(range(len(parameters[1])))
    x_tick_labels = [str(round(parameters[1][i])) if i % 2 else "" for i in x_ticks]
    non_blank = [label for label in x_tick_labels if label]
    if len(set(non_blank)) != len(non_blank):
        x_tick_labels = ["%0.1f" % parameters[1][i] if i % 2 else "" for i in x_ticks]
    ax.set_xticks(x_ticks, labels=x_tick_labels)

    y_ticks = list(range(len(parameters[0])))
    y_tick_labels = [str(round(parameters[0][i])) if i % 2 else "" for i in y_ticks]
    non_blank = [label for label in y_tick_labels if label]
    if len(set(non_blank)) != len(non_blank):
        y_tick_labels = ["%0.1f" % parameters[0][i] if i % 2 else "" for i in y_ticks]
    ax.set_yticks(y_ticks, labels=y_tick_labels)

    ax.tick_params(length=0)
    ax.grid(which="minor", linestyle="-", linewidth=1.5)
    ax.set_xticks(np.arange(-0.5, len(parameters[1]), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(parameters[0]), 1), minor=True)
    ax.set_xlabel(parameter_labels[1])
    ax.set_ylabel(parameter_labels[0])
    ax.set_title(title)


def main(args):
    results_dir = args.results_dir
    if os.path.basename(os.path.normpath(results_dir)) != "sweep":
        results_dir = os.path.join(results_dir, "sweep")

    plt.rcParams["font.size"] = 22
    plt.rcParams["font.family"] = "serif"
    fig, axes = plt.subplots(len(METHOD_SUBSETS), len(MODEL_NAMES), figsize=(18, 24))

    subset_titles = []
    for row, (subset_name, method_names) in enumerate(METHOD_SUBSETS):
        n_success = 0
        n_total = 0
        for col, model_name in enumerate(MODEL_NAMES):
            virtual_best = load_virtual_best(
                results_dir, model_name, method_names, args.suffix
            )
            n_success += virtual_best["success"].sum()
            n_total += len(virtual_best)
            plot_convergence(
                axes[row, col],
                virtual_best,
                PARAMETER_LOOKUP[model_name],
                PARAMETER_LABEL_LOOKUP[model_name],
                TITLE_LOOKUP[model_name],
            )
        if subset_name == "No elimination":
            space = "        "
        else:
            space = ""
        subset_titles.append(f"{subset_name} ({100 * n_success / n_total:.0f}%){space}")
        #axes[row, 1].annotate(
        #    subset_name,
        #    xy=(0.5, 1.15),
        #    xycoords="axes fraction",
        #    ha="center",
        #    va="bottom",
        #    fontsize=32,
        #    fontweight="bold",
        #)

    if not args.no_legend:
        fig.legend(
            handles=[
                Patch(color=CMAP(0), label="Unsuccessful"),
                Patch(color=CMAP(1), label="Successful"),
            ],
            loc="upper right",
            ncol=2,
        )

    y_positions = [0.24, 0.49, 0.75, 1.0]
    y_positions.reverse()
    for row, subset_title in enumerate(subset_titles):
        row_top = axes[row, 0].get_position().y1
        fig.text(
            0.53,
            #row_top + 0.01,
            y_positions[row],
            subset_title,
            ha="center",
            va="top",
            fontsize=30,
            fontweight="bold",
        )

    fig.tight_layout(h_pad=5.0, w_pad=-5.0, rect=(0.0, 0.0, 1.0, 0.985))
    if not args.no_save:
        suffix_str = "" if args.suffix is None else f"-{args.suffix}"
        fpath = os.path.join(
            args.image_dir, f"virtual-best-sweep-convergence{suffix_str}.pdf"
        )
        print(f"Saving figure to {fpath}")
        fig.savefig(fpath, transparent=not args.opaque)

    if args.show:
        plt.show()


if __name__ == "__main__":
    argparser = config.get_plot_argparser()
    args = argparser.parse_args()
    main(args)
