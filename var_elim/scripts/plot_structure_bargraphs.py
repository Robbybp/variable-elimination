import os
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np

from var_elim.scripts.config import get_plot_argparser


plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 20
plt.rcParams["font.family"] = "serif"


method_ord = {"no-elim": 0, "d1": 1, "trivial": 2, "linear-d2": 3, "d2": 4, "ampl": 5, "matching": 6}
method_to_label = {"no-elim": "--", "d1": "LD1", "trivial": "JP", "linear-d2": "LD2", "d2": "D2", "ampl": "GR", "matching": "LM"}
model_to_label = {"distill": "Distillation", "mb-steady": "Moving bed", "opf": "OPF", "pipeline": "Pipeline"}


def _plot_fraction_eliminated(df):
    # Models alphabetically are the same order as in the paper table, luckily.
    models = list(sorted(set(df["model"])))
    methods = list(sorted(set(df["method"]), key=lambda m: method_ord[m]))
    nvar_by_model = {}
    nelim_by_model_method = {}

    for i, row in df.iterrows():
        if row["method"] == "no-elim":
            nvar_by_model[row["model"]] = row["nvar"]
        nelim_by_model_method[row["model"], row["method"]] = row["n-elim"]

    nmodel = len(models)
    nmethod = len(methods)
    inter_model_spacing = 0.5
    width = 1.0
    tickpos_by_model = {
        model: (nmethod + inter_model_spacing) * i + nmethod / 2
        for i, model in enumerate(models)
    }
    offset_by_method = {method: method_ord[method] - 3.5 for method in methods}

    to_omit = ("no-elim",)
    model_methods = [(mod, met) for mod in models for met in methods if met not in to_omit]

    x_array = np.array([
        tickpos_by_model[model] + offset_by_method[method]
        for model, method in model_methods
    ])
    frac_elim_array = np.array([
        nelim_by_model_method[model, method] / nvar_by_model[model]
        for model, method in model_methods
    ])
    percent_elim_array = 100 * frac_elim_array

    fig, ax = plt.subplots()
    # Using different colors for each method doesn't look good
    #colors = ["blue", "orange", "green", "purple", "brown", "red"]
    ax.bar(
        x_array,
        percent_elim_array,
        width=width,
        align="center",
        edgecolor="black",
    )

    w, h = fig.get_size_inches()
    fig.set_size_inches(1.2*w, h)
    ax.set_ylabel(
        "Percent of\\\\\nvariables\\\\\neliminated",
        rotation=0,
        labelpad=50,
    )
    if args.title is not None:
        ax.set_title(args.title)

    # A little hacky. Here we're adding extra tick labels for the first instance
    # of each method. Assumming x array is ordered by method, and that we omit exactly
    # one method.
    method_tickpos = list(x_array[0:nmethod-1])
    method_ticklabels = [method_to_label[method] for method in methods if method not in to_omit]
    model_tickpos = [tickpos_by_model[model] for model in models]
    model_ticklabels = [
        "\n\n" + model_to_label[model] if i == 0 else "\n" + model_to_label[model]
        for i, model in enumerate(models)
    ]
    ax.set_xticks(model_tickpos, model_ticklabels, rotation=0, minor=False)
    ax.set_xticks(
        method_tickpos,
        method_ticklabels,
        rotation=90,
        minor=True,
        fontsize=14,
        ha="center",
    )
    ax.xaxis.set_tick_params(which="minor", length=0)

    fig.tight_layout()
    
    return fig, ax

def _plot_nnz_per_constraint(df):
    # Models alphabetically are the same order as in the paper table, luckily.
    models = list(sorted(set(df["model"])))
    methods = list(sorted(set(df["method"]), key=lambda m: method_ord[m]))
    nonlin_nnz_per_con = {}
    lin_nnz_per_con = {}

    for i, row in df.iterrows():
        nonlin_nnz_per_con[row["model"], row["method"]] = (row["nnz"] -  row["nnz-linear"])/row["ncon"]
        lin_nnz_per_con[row["model"], row["method"]] = row["nnz-linear"]/row["ncon"]

    nmodel = len(models)
    nmethod = len(methods)
    inter_model_spacing = 1.5
    width = 1.0
    tickpos_by_model = {
        model: (nmethod + inter_model_spacing) * i + nmethod / 2
        for i, model in enumerate(models)
    }
    offset_by_method = {method: method_ord[method] - 3.5 for method in methods}

    to_omit = ()
    model_methods = [(mod, met) for mod in models for met in methods if met not in to_omit]

    x_array = np.array([
        tickpos_by_model[model] + offset_by_method[method]
        for model, method in model_methods
    ])
    
    fig, ax = plt.subplots()
    # Using different colors for each method doesn't look good
    #colors = ["blue", "orange", "green", "purple", "brown", "red"]
    
    nonlin_nnz_per_con_array = np.array([
        nonlin_nnz_per_con[model, method]
        for model, method in model_methods
    ])
    
    lin_nnz_per_con_array = np.array([
        lin_nnz_per_con[model, method]
        for model, method in model_methods
    ])
    ax.bar(
        x_array,
        nonlin_nnz_per_con_array,
        width=width,
        align="center",
        edgecolor="black",
        label = 'Nonlinear NNZ/con'
    )
    
    ax.bar(
        x_array,
        lin_nnz_per_con_array,
        bottom = nonlin_nnz_per_con_array,
        width=width,
        align="center",
        edgecolor="black",
        label = 'Linear NNZ/con'
    )

    w, h = fig.get_size_inches()
    fig.set_size_inches(1.3*w, h)
    ax.set_ylabel(
        "Number of\\\\\nnon-zeros per\\\\\nconstraint",
        rotation=0,
        labelpad=60,
    )
    if args.title is not None:
        ax.set_title(args.title)

    # A little hacky. Here we're adding extra tick labels for the first instance
    # of each method. Assumming x array is ordered by method, and that we omit exactly
    # one method.
    method_tickpos = list(x_array[0:nmethod-len(to_omit)])
    method_ticklabels = [method_to_label[method] for method in methods if method not in to_omit]
    model_tickpos = [tickpos_by_model[model] for model in models]
    model_ticklabels = [
        "\n\n" + model_to_label[model] if i == 0 else "\n" + model_to_label[model]
        for i, model in enumerate(models)
    ]
    
    ax.set_xticks(model_tickpos, model_ticklabels, rotation=0, minor=False)
    ax.set_xticks(
        method_tickpos,
        method_ticklabels,
        rotation=90,
        minor=True,
        fontsize=14,
        ha="center",
    )
    ax.xaxis.set_tick_params(which="minor", length=0)
    plt.legend(loc = "upper left")
    fig.tight_layout()
    
    return fig, ax

if __name__ == "__main__":
    argparser = get_plot_argparser()
    argparser.add_argument(
        "structure_results",
        help="CSV file containing structural results",
    )
    args = argparser.parse_args()

    if not (args.structure_results.endswith(".csv") or args.structure_results.endswith(".CSV")):
        raise ValueError("structure results file must end with '.csv' or '.CSV'")
    if args.model is not None:
        raise ValueError("--model argument cannot be used")
    if args.method is not None:
        raise ValueError("--method argument cannot be used")

    df = pd.read_csv(args.structure_results)
    fig1, ax1 = _plot_fraction_eliminated(df)
    fig2, ax2 = _plot_nnz_per_constraint(df)
    if not args.no_save:
        suff_str = "" if args.suffix is None else f"-{args.suffix}"
        fname = f"fraction-elim{suff_str}.pdf"
        fpath = os.path.join(args.image_dir, fname)
        fig1.savefig(fpath, transparent=not args.opaque)
        
        suff_str = "" if args.suffix is None else f"-{args.suffix}"
        fname = f"nnz-per-con{suff_str}.pdf"
        fpath = os.path.join(args.image_dir, fname)
        fig2.savefig(fpath, transparent=not args.opaque)

    if args.show:
        plt.show()
