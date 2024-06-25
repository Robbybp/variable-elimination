#  ___________________________________________________________________________
#
#  Variable Elimination: Research code for variable elimination in NLPs
#
#  Copyright (c) 2023. Triad National Security, LLC. All rights reserved.
#
#  This program was produced under U.S. Government contract 89233218CNA000001
#  for Los Alamos National Laboratory (LANL), which is operated by Triad
#  National Security, LLC for the U.S. Department of Energy/National Nuclear
#  Security Administration. All rights in the program are reserved by Triad
#  National Security, LLC, and the U.S. Department of Energy/National Nuclear
#  Security Administration. The Government is granted for itself and others
#  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
#  in this material to reproduce, prepare derivative works, distribute copies
#  to the public, perform publicly and display publicly, and to permit others
#  to do so.
#
#  This software is distributed under the 3-clause BSD license.
#  ___________________________________________________________________________

import os
import itertools
import var_elim.scripts.config as config
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np

# These parameters depend on the problem we are trying to solve.
# TODO: How to determine which problem we are solving? Do we have to pass some
# metadata to this script as well?
PARAMETER_NAMES = [
    (
        "fs.moving_bed.solid_phase.properties[0,1].temperature",
        "fs.moving_bed.solid_phase.properties[0,1].flow_mass",
    ),
    ("vol", "x_Feed"),
    ("fs.nodes[0].state[*].temperature", "fs.nodes[0].state[*].pressure")
]

PARAMETER_LABELS = [
    ("Temperature (K)", "Flow rate (kg/s)"),
    ("Relative volatility", "Feed mole fraction"),
    ("Gas temperature", "Supply pressure"),
]

MODEL_NAMES = [
    "mb-steady",
    "distill",
    "pipeline",
]

PARAMETER_LOOKUP = dict(zip(MODEL_NAMES, PARAMETER_NAMES))
PARAMETER_LABEL_LOOKUP = dict(zip(MODEL_NAMES, PARAMETER_LABELS))

CMAP = ListedColormap([(0.95, 0.55, 0.55), (0.05, 0.05, 0.35)])


def plot_convergence(
    df,
    parameter_names,
    parameter_labels=None,
    legend=True,
):
    plt.rcParams["font.size"] = 16
    plt.rcParams["font.family"] = "serif"

    # If we are going to plot parameter sweep results, we better have two parameters...
    parameters = [list(sorted(set(df[name]))) for name in parameter_names]
    param_index_maps = [{p: i for i, p in enumerate(params)} for params in parameters]
    success_lookup = {}
    for i in range(len(df)):
        key = tuple(df[name][i] for name in parameter_names)
        success_lookup[key] = int(df["success"][i])
    n_success = list(success_lookup.values()).count(1)
    n_total = len(df)
    print(f"Converged {n_success} / {n_total} instances")
    convergence_array = np.zeros(tuple(len(params) for params in parameters))
    for params in itertools.product(*parameters):
        indices = tuple(idx_map[p] for idx_map, p in zip(param_index_maps, params))
        convergence_array[indices] = success_lookup[params]

    fig, ax = plt.subplots()

    ax.imshow(
        convergence_array,
        aspect="equal",
        origin="lower",
        cmap=CMAP,
        vmin=0,
        vmax=1,
    )
    # Since we are plotting on a 2D grid, it is unclear if the more general code
    # above to parse an arbitrary number of parameters has any value.
    assert len(parameters) == 2

    # Label every grid cell; turn off (major) tick marks
    x_ticks = [i for i in range(len(parameters[1]))]
    # TODO: If any two labels are the same (rounded to nearest int), then
    # add a decimal place.
    x_tick_labels = [str(round(parameters[1][i])) if i%2 else "" for i in x_ticks]
    ax.set_xticks(x_ticks, labels=x_tick_labels)

    y_ticks = [i for i in range(len(parameters[0]))]
    y_tick_labels = [str(round(parameters[0][i])) if i%2 else "" for i in y_ticks]
    ax.set_yticks(y_ticks, labels=y_tick_labels)

    # Turn off (major) tick marks for both axes
    ax.tick_params(length=0)

    # Set gridlines
    ax.grid(which="minor", linestyle="-", linewidth=1.5)
    ax.tick_params(length=0)
    ax.set_xticks(np.arange(-0.5, len(parameters[1]), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(parameters[0]), 1), minor=True)

    if parameter_labels is None:
        parameter_labels = parameter_names
    xlabel = parameter_labels[1]
    ylabel = parameter_labels[0]
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if legend:
        # Add labels
        colors = [CMAP(0), CMAP(1)]
        labels = ["Unsuccessful", "Successful"]
        patches = [
            Patch(color=color, label=label) for color, label in zip(colors, labels)
        ]
        ax.legend(
            handles=patches,
            loc=(1.05, 0.8),
        )
        w, h = fig.get_size_inches()
        fig.set_size_inches(1.2*w, h)

    return fig, ax


def main(args):
    model_name = None
    if args.model is None:
        for name in MODEL_NAMES:
            if name in os.path.basename(args.fpath):
                model_name = name
                break
        if model_name is None:
            raise RuntimeError(
                "model was not provided and cannot be inferred from file name"
            )
    elif args.model in MODEL_NAMES:
        model_name = args.model
    else:
        raise ValueError("Unrecognized model name provided")

    df = pd.read_csv(args.fpath)
    parameter_names = PARAMETER_LOOKUP[model_name]
    parameter_labels = PARAMETER_LABEL_LOOKUP[model_name]
    fig, ax = plot_convergence(
        df, parameter_names, parameter_labels=parameter_labels
    )

    if args.title is not None:
        ax.set_title(args.title)

    fig.tight_layout()

    if not args.no_save:
        input_basename = os.path.basename(args.fpath)
        if args.suffix is not None:
            # Since we just replace the extension of the input file, we don't allow
            # adding a suffix here. The suffix should just be added to the input
            # file name.
            raise ValueError("--suffix is not supported when plotting convergence")
        if input_basename.endswith(".csv"):
            plot_fname = input_basename.replace(".csv", "-convergence.pdf")
        elif input_basename.endswith(".CSV"):
            plot_fname = input_basename.replace(".CSV", "-convergence.pdf")
        plot_fpath = os.path.join(args.image_dir, plot_fname)

        fig.savefig(plot_fpath, transparent=not args.opaque)

    if args.show:
        plt.show()


if __name__ == "__main__":
    argparser = config.get_plot_argparser()
    argparser.add_argument("fpath", help="CSV file with parameter sweep results to plot")
    args = argparser.parse_args()
    if not args.fpath.endswith(".csv") and not args.fpath.endswith(".CSV"):
        raise ValueError("fpath must end with .csv or .CSV")
    main(args)
