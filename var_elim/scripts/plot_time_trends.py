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
from var_elim.scripts.config import get_plot_argparser
import matplotlib.pyplot as plt
import collections
import pandas as pd


SolveTimeTuple = collections.namedtuple(
    "SolveTimeTuple", ["function", "jacobian", "hessian", "other"]
)


COLORS = ["blue", "orange", "green", "brown"]
MARKERS = ["o", "s", "^", "*", "d"]


XVAL_GETTER = {
    "nodecount":           lambda df, times, model: [row["nnode-nl-linear"] + row["nnode-nl-nonlinear"] for _, row in df.iterrows() if row["model"] == model],
    "nonlinear-nodecount": lambda df, times, model: [row["nnode-nl-nonlinear"] for _, row in df.iterrows() if row["model"] == model],
    "function":            lambda df, times, model: [times[row["model"], row["method"]].function for _, row in df.iterrows() if row["model"] == model],
    "jacobian":            lambda df, times, model: [times[row["model"], row["method"]].jacobian for _, row in df.iterrows() if row["model"] == model],
    "nvarcon":             lambda df, times, model: [row["nvar"]+row["nvar"] for _, row in df.iterrows() if row["model"] == model],
    "nvarnode":            lambda df, times, model: [row["nvar"]*row["nnode-nl-nonlinear"] for _, row in df.iterrows() if row["model"] == model],
    "nvarfunction":        lambda df, times, model: [row["nvar"]*times[row["model"], row["method"]].function for _, row in df.iterrows() if row["model"] == model],
    "nnz-total":           lambda df, times, model: [2*row["nnz"] + row["nnz-hessian"] for _, row in df.iterrows() if row["model"] == model],
}


XTYPE_GETTER = {
    "function": "nodecount",
    "jacobian": "nodecount",
    #"hessian": "nvarnode",
    "hessian": "nonlinear-nodecount",
    #"other": "nvarcon",
    "other": "nnz-total",
}


def plot_evaluation_time(
    structure_df,
    solvetime_df,
    which="function",
    models=None,
    normalize=False,
    fig_ax=None,
    markers=None,
    colors=None,
    labels=None,
    xtype=None,  # See XVAL_GETTER
):
    if fig_ax is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = fig_ax
    if markers is None:
        markers = MARKERS
    if colors is None:
        colors = COLORS
    if xtype is None:
        xtype = XTYPE_GETTER[which]

    solvetime_lookup = {
        (row["model"], row["method"]): SolveTimeTuple(
            row["function-per100"],
            row["jacobian-per100"],
            row["hessian-per100"],
            row["other-per100"],
        )
        for i, row in solvetime_df.iterrows()
    }

    if models is None:
        model_set = set(structure_df["model"])
        models = list(sorted(model_set))
    else:
        model_set = set(models)
    if labels is None:
        labels = models

    #nodecounts = {}
    xvalues = {}
    xreference = {}
    yreference = {}
    times = {}
    for model in models:
        # TODO: Unclear that we should count linear and nonlinear equally, but
        # also unclear how to combine them. Maybe just use the nonlinear count?
        #nodecounts[model] = [
        #    row["nnode-nl-linear"] + row["nnode-nl-nonlinear"]
        #    for i, row in structure_df.iterrows()
        #    if row["model"] == model
        #]
        xvalues[model] = XVAL_GETTER[xtype](structure_df, solvetime_lookup, model)
        # For each (model, method) present in the structure table, assume it exists
        # in the solvetime table and look up.
        solvetimes = [
            solvetime_lookup[row["model"], row["method"]]
            for i, row in structure_df.iterrows()
            if row["model"] == model
        ]
        times[model] = [getattr(time, which) for time in solvetimes]

        # NOTE: Assuming that the first item corresponds to no-elim
        xreference[model] = xvalues[model][0]
        yreference[model] = getattr(solvetime_lookup[model, "no-elim"], which)

    for i, model in enumerate(models):
        if normalize:
            xref = xreference[model]
            yref = yreference[model]
            norm_xvalues = [x/xref for x in xvalues[model]]
            norm_times = [y/yref for y in times[model]]
        else:
            norm_xvalues = xvalues[model]
            norm_times = times[model]
        ax.scatter(
            norm_xvalues,
            norm_times,
            marker=markers[i],
            color=colors[i],
            label=labels[i],
        )
    ax.legend()
    return fig, ax


def main(args):
    structure_df = pd.read_csv(args.structure_fpath)
    solvetime_df = pd.read_csv(args.solvetime_fpath)
    if args.model is None:
        models = None
    else:
        models = [args.model]

    fig, ax = plot_evaluation_time(
        structure_df,
        solvetime_df,
        which=args.which,
        models=models,
        normalize=args.normalize
    )

    if not args.no_save:
        input_basename = os.path.basename(args.structure_fpath)
        if args.suffix is not None:
            # Since we just replace the extension of the input file, we don't allow
            # adding a suffix here. The suffix should just be added to the input
            # file name.
            raise ValueError("--suffix is not supported when plotting evaluation times")

        if input_basename.endswith(".csv"):
            plot_fname = input_basename.replace(".csv", "-functiontime.pdf")
        elif input_basename.endswith(".CSV"):
            plot_fname = input_basename.replace(".CSV", "-functiontime.pdf")
        plot_fpath = os.path.join(args.image_dir, plot_fname)

        fig.savefig(plot_fpath, transparent=not args.opaque)

    if args.show:
        plt.show()


if __name__ == "__main__":
    argparser = get_plot_argparser()
    argparser.add_argument("structure_fpath", help="CSV file of structural results")
    argparser.add_argument("solvetime_fpath", help="CSV file of solve-time results")
    argparser.add_argument("--which", default="function", help="function, jacobian, hessian, or other. Default=function")
    argparser.add_argument("--normalize", action="store_true")
    args = argparser.parse_args()
    main(args)
