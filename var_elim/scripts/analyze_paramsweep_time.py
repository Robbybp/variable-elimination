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

import itertools
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import var_elim.scripts.config as config
import os 
import pandas as pd

# This script takes in a data frame containing 
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
    ("Gas temperature (K)", "Supply pressure (bar)"),
]

MODEL_NAMES = [
    "mb-steady",
    "distill",
    "pipeline",
]

PARAMETER_LOOKUP = dict(zip(MODEL_NAMES, PARAMETER_NAMES))
PARAMETER_LABEL_LOOKUP = dict(zip(MODEL_NAMES, PARAMETER_LABELS))

def analyze_sweep(df, 
                  parameter_names,
                  parameter_labels=None,
                  title=None):
    
    plt.rcParams["font.size"] = 20
    plt.rcParams["font.family"] = "serif"
    
    # If we are going to plot parameter sweep results, we better have two parameters...
    parameters = [list(sorted(set(df[name]))) for name in parameter_names]
    param_index_maps = [{p: i for i, p in enumerate(params)} for params in parameters]
    success_lookup = {}
    solve_time = {}
    for i in range(len(df)):
        key = tuple(df[name][i] for name in parameter_names)
        success_lookup[key] = int(df["success"][i])
        solve_time[key] = round(df["solve-time"][i])
    n_success = list(success_lookup.values()).count(1)
    n_total = len(df)
    print(f"Converged {n_success} / {n_total} instances")
    convergence_array = np.zeros(tuple(len(params) for params in parameters))
    
    solve_time_array = np.zeros(tuple(len(params) for params in parameters))
    for params in itertools.product(*parameters):
        indices = tuple(idx_map[p] for idx_map, p in zip(param_index_maps, params))
        convergence_array[indices] = success_lookup[params]
        solve_time_array[indices] = solve_time[params]
    
    fig, ax = plt.subplots()
    ax.imshow(
        solve_time_array,
        aspect="equal",
        origin="lower"
    )
    
    for i in range(convergence_array.shape[0]):
        for j in range(convergence_array.shape[1]):
            if convergence_array[i, j] == 1:
                t = solve_time_array[i, j]
                txt = f"{t:.0f}"
    
            else:
                txt = "×"
    
            ax.text(
                j, i,
                txt,
                ha="center",
                va="center",
                fontsize=14,
                color="white",
            )
    
    set_figure_axes(parameters, 
                        parameter_labels, 
                        parameter_names, 
                        title, 
                        ax)
    

    return fig, ax

def set_figure_axes(parameters, 
                    parameter_labels, 
                    parameter_names, 
                    title, 
                    ax):
    
    # Label every grid cell; turn off (major) tick marks
    x_ticks = [i for i in range(len(parameters[1]))]
    # TODO: If any two labels are the same (rounded to nearest int), then
    # add a decimal place.
    x_tick_labels = [str(round(parameters[1][i])) if i%2 else "" for i in x_ticks]
    non_blank = [l for l in x_tick_labels if l != ""]
    if len(set(non_blank)) != len(non_blank):
        # There are duplicates!
        x_tick_labels = ["%0.1f" % parameters[1][i] if i%2 else "" for i in x_ticks]
    ax.set_xticks(x_ticks, labels=x_tick_labels)

    y_ticks = [i for i in range(len(parameters[0]))]
    y_tick_labels = [str(round(parameters[0][i])) if i%2 else "" for i in y_ticks]
    non_blank = [l for l in y_tick_labels if l != ""]
    if len(set(non_blank)) != len(non_blank):
        # There are duplicates!
        y_tick_labels = ["%0.1f" % parameters[1][i] if i%2 else "" for i in y_ticks]
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

    if title is not None:
        ax.set_title(title)


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

    from var_elim.scripts.config import ELIM_NAMES
    ELIM_NAME_TO_TITLE = {
        "no-elim": "Original model",
        "d1": "Linear degree-one",
        "trivial": "Equal-coefficient degree-two",
        "linear-d2": "Linear degree-two",
        "d2": "Degree two",
        "ampl": "Greedy",
        "matching": "Linear matching",
    }
    title = None
    if args.title is None:
        print("--title not provided. Attempting to infer from filename")
        for name in ELIM_NAMES:
            if name in os.path.basename(args.fpath):
                print(f"Method recognized as {name}. Setting title")
                if title is not None:
                    if "linear" in os.path.basename(args.fpath):
                        # This is a hack so we don't recognize linear-d2 as d2
                        print("'linear' in method name. not resetting title")
                        continue
                    print("WARNING: multiple matches. Resetting title to None")
                    title = None
                    break
                else:
                    title = ELIM_NAME_TO_TITLE[name]
        if title is None:
            print("Could not infer title from filename")

    fig, ax = analyze_sweep(
        df,
        parameter_names,
        parameter_labels=parameter_labels,
        title=title,
    )

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

