import os
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np

from var_elim.scripts.config import get_plot_argparser


plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 20
plt.rcParams["font.family"] = "serif"


method_ord = {"no-elim": 0, "d1": 1, "ecd2": 2, "linear-d2": 3, "d2": 4, "greedy": 5, "matching": 6}
method_to_label = {"no-elim": "--", "d1": "LD1", "ecd2": "ECD2", "linear-d2": "LD2", "d2": "D2", "greedy": "GR", "matching": "LM"}
model_to_label = {"distill": "Distillation", "mb-steady": "Moving bed", "opf": "OPF", "pipeline": "Pipeline"}

def _plot_solve_time_fractions(df):
    # Models alphabetically are the same order as in the paper table, luckily.
    models = list(sorted(set(df["model"])))
    methods = list(sorted(set(df["method"]), key=lambda m: method_ord[m]))
    total_time_by_model_method = {}
    smallest_total_time = {}
    smallest_time_per_iteration = {}
    func_eval_time_by_model_method = {}
    jac_eval_time_by_model_method = {}
    hess_eval_time_by_model_method = {}
    other_eval_time_by_model_method = {}
    iter_by_method = {}

    for i, row in df.iterrows():
        total_time_by_model_method[row["model"], row["method"]] = row["solve-time"]
        func_eval_time_by_model_method[row["model"], row["method"]] = row["function-time"]
        jac_eval_time_by_model_method[row["model"], row["method"]] = row["jacobian-time"]
        hess_eval_time_by_model_method[row["model"], row["method"]] = row["hessian-time"]
        other_eval_time_by_model_method[row["model"], row["method"]] = row["solve-time"] - (row["function-time"] 
                                                                                            + row["jacobian-time"]
                                                                                            + row["hessian-time"])
        iter_by_method[row["model"], row["method"]] = row['n-iter']

    for model in models:
        smallest_total_time[model] = min(value for (key1, key2), value in total_time_by_model_method.items() if key1 == model)
        smallest_time_per_iteration[model] = min(value/iter_by_method[key1, key2] for (key1, key2), value in total_time_by_model_method.items() if key1 == model)
    
    nmodel = len(models)
    nmethod = len(methods)
    inter_model_spacing = 1
    width = 1
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
    
    normalized_fun_eval_time_array = np.array([
        func_eval_time_by_model_method[model, method] / iter_by_method[model,method]/smallest_time_per_iteration[model]
        for model, method in model_methods
    ])
    
    normalized_jac_eval_time_array = np.array([
        jac_eval_time_by_model_method[model, method] / iter_by_method[model,method]/smallest_time_per_iteration[model]
        for model, method in model_methods
    ])
    
    normalized_hess_eval_time_array = np.array([
        hess_eval_time_by_model_method[model, method] / iter_by_method[model,method]/smallest_time_per_iteration[model]
        for model, method in model_methods
    ])
    
    normalized_other_eval_time_array = np.array([
        other_eval_time_by_model_method[model, method] / iter_by_method[model,method]/smallest_time_per_iteration[model]
        for model, method in model_methods
    ])
    
    

    fig, ax = plt.subplots()
    # Using different colors for each method doesn't look good
    #colors = ["blue", "orange", "green", "purple", "brown", "red"]
    ax.bar(
        x_array,
        normalized_other_eval_time_array,
        width=width,
        align="center",
        edgecolor="black",
        label = 'Other'
    )
    
    ax.bar(
        x_array,
        normalized_hess_eval_time_array,
        width=width,
        align="center",
        edgecolor="black",
        bottom = normalized_other_eval_time_array,
        label = 'Hessian eval'
    )
    
    ax.bar(
        x_array,
        normalized_jac_eval_time_array,
        width=width,
        align="center",
        edgecolor="black",
        bottom = normalized_other_eval_time_array + normalized_hess_eval_time_array,
        label = 'Jacobian eval'
    )
    
    ax.bar(
        x_array,
        normalized_fun_eval_time_array,
        width=width,
        align="center",
        edgecolor="black",
        bottom = normalized_other_eval_time_array + normalized_hess_eval_time_array + normalized_jac_eval_time_array,
        label = 'Function eval'
    )
    

    #"solvetime nper iter/fastest solvetime per iter"
    w, h = fig.get_size_inches()
    fig.set_size_inches(1.5*w, h)
    ax.set_ylabel(
        "Normalized\\\\\nsolvetime per\\\\\niteration",
        rotation=0,
        labelpad=50,
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
    plt.legend()
    fig.tight_layout()
    
    return fig, ax

if __name__ == "__main__":
    argparser = get_plot_argparser()
    argparser.add_argument(
        "timing_results",
        help="CSV file containing structural results",
    )
    args = argparser.parse_args()

    if not (args.timing_results.endswith(".csv") or args.structure_results.endswith(".CSV")):
        raise ValueError("timing results file must end with '.csv' or '.CSV'")
    if args.model is not None:
        raise ValueError("--model argument cannot be used")
    if args.method is not None:
        raise ValueError("--method argument cannot be used")

    df = pd.read_csv(args.timing_results)
    fig, ax = _plot_solve_time_fractions(df)
 
    if not args.no_save:
        suff_str = "" if args.suffix is None else f"-{args.suffix}"
        fname = f"fraction-solvetime{suff_str}.pdf"
        fpath = os.path.join(args.image_dir, fname)
        print(f"Writing figure to {fpath}")
        fig.savefig(fpath, transparent=not args.opaque)
        

    if args.show:
        plt.show()
