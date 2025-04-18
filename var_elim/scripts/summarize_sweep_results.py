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

import pyomo.environ as pyo
from pyomo.common.timing import HierarchicalTimer
import pselib
from idaes.core.util.parameter_sweep import (
    ParameterSweepSpecification,
    SequentialSweepRunner,
    UniformSampling,
)
import pandas as pd

import json
from collections import namedtuple

from var_elim.algorithms.validate import validate_solution
from var_elim.elimination_callbacks import (
    no_elim_callback,
    trivial_elim_callback,
    linear_d2_elim_callback,
    d2_elim_callback,
    matching_elim_callback,
)


import os
import var_elim.scripts.config as config

SolveResults = namedtuple(
    "SolveResults",
    ["objective", "feasible", "timer"],
)


def summarize_single_model(args, elim_names):
    sweep_data = {
        "model": [],
        "method": [],
        "n-success": [],
        "n-total": [],
        "percent-success": [],
        "ave-elim-time": [],
        "ave-solve-time": [],
    }
    # NOTE: Code is repeated below. TODO: Consolidate if this changes.
    suff_str = "-" + args.suffix if args.suffix is not None else ""
    for elimname in elim_names:
        fname = args.model + "-" + elimname + "-sweep" + suff_str + ".csv"
        fpath = os.path.join(args.results_dir, fname)
        df = pd.read_csv(fpath)

        n_success = list(df["success"]).count(True)
        n_total = len(df["success"])
        ave_solvetime = sum(df["solve-time"])/len(df["solve-time"])
        ave_elimtime = sum(df["elim-time"])/len(df["elim-time"])

        sweep_data["model"].append(args.model)
        sweep_data["method"].append(elimname)
        sweep_data["n-success"].append(n_success)
        sweep_data["n-total"].append(n_total)
        sweep_data["percent-success"].append(n_success / n_total * 100)
        sweep_data["ave-elim-time"].append(ave_elimtime)
        sweep_data["ave-solve-time"].append(ave_solvetime)

    output_df = pd.DataFrame(sweep_data)
    return output_df


SWEEP_MODELS = ["distill", "mb-steady", "pipeline"]


def summarize_all_models(args, elim_names):
    data = dict(method=[])
    for model in SWEEP_MODELS:
        data[model] = []
    # Insert this last so table is in correct order
    data["total"] = []

    orig_model_arg = args.model
    model_dfs = dict()
    for model in SWEEP_MODELS:
        # HACK: temporarily override args.model
        args.model = model
        model_dfs[model] = summarize_single_model(args, elim_names)
        args.model = orig_model_arg

    model_method_rows = {}
    for model in SWEEP_MODELS:
        df = model_dfs[model]
        for i, row in df.iterrows():
            model_method_rows[model, row["method"]] = row

    for method in elim_names:
        data["method"].append(method)
        total_instances = 0
        total_success = 0
        for model in SWEEP_MODELS:
            # Should we store the percent successful, or the raw number?
            row = model_method_rows[model, method]
            data[model].append(row["percent-success"])
            total_instances += row["n-total"]
            total_success += row["n-success"]
        total_percent = 100.0 * total_success / total_instances
        data["total"].append(total_percent)

    output_df = pd.DataFrame(data)
    return output_df


def main(args):
    if args.method is None:
        elimination_callbacks = config.ELIM_CALLBACKS
    else:
        elimination_callbacks = [(args.method, config.ELIM_LOOKUP[args.method])]
    elim_names = [name for name, _ in elimination_callbacks]

    if args.model is None:
        output_df = summarize_all_models(args, elim_names)
    else:
        output_df = summarize_single_model(args, elim_names)

    if args.method is None:
        method_str = ""
    else:
        method_str = f"-{args.method}"

    method_str = "" if args.method is None else f"{args.method}-"
    model_str = "" if args.model is None else f"{args.model}-"
    suff_str = "" if args.suffix is None else f"-{args.suffix}"

    output_fname = model_str + method_str + "sweep-summary" + suff_str + ".csv"

    if args.method is None:
        # If this is the summary for a particular model, we put it in the top-level
        # results dir.
        if os.path.basename(args.results_dir) == "sweep":
            resdir = os.path.dirname(args.results_dir)
        else:
            resdir = args.results_dir
        output_fpath = os.path.join(resdir, output_fname)
    else:
        output_fpath = os.path.join(args.results_dir, output_fname)
    if not args.no_save:
        print(f"Writing sweep summary to {output_fpath}")
        output_df.to_csv(output_fpath)

    print(output_df)


if __name__ == "__main__":
    argparser = config.get_argparser()

    # HACK: We change the default of the argparser so we can handle it specially
    # if --method or --model are used.
    # It's unclear whether this hack will be worth the convenience, but let's try it.
    argparser.set_defaults(results_dir=None)

    args = argparser.parse_args()

    if args.results_dir is None:
        if args.method is None and args.model is None:
            # If neither method nor model is used (we are collecting all results)
            # we put results in the top-level results directory.
            args.results_dir = config.get_results_dir()
        else:
            # If either method or model is used, we put the results in the
            # results/structure subdirectory. This is because we don't want the
            # top-level results getting polluted with a bunch of files.
            resdir = os.path.join(config.get_results_dir(), "sweep")
            config.validate_dir(resdir)
            args.results_dir = resdir

    main(args)
