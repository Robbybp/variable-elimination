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

import pandas as pd
import var_elim.scripts.config as config
import os


FILEDIR = os.path.dirname(__file__)

script_lookup = dict(
    structure="analyze_structure.py",
    solvetime="analyze_solvetime.py",
    sweep="run_param_sweep.py",
)
script_lookup["plot-sweep"] = "plot_sweep_results.py"


def main(args):
    cl_lists = []

    if args.by == "both" or args.by == "model":
        model_names = list(config.MODEL_NAMES)
    else:
        model_names = []

    if args.by == "both" or args.by == "method":
        elim_names = list(config.ELIM_NAMES)
    else:
        elim_names = []

    scriptname = script_lookup[args.result_type]
    # Script name is WRT file directory
    if scriptname not in os.listdir(os.getcwd()):
        # Prepend FILEDIR relative to CWD
        relative_dirname = os.path.relpath(FILEDIR, os.getcwd())
        scriptname = os.path.join(relative_dirname, scriptname)

    # Does it make sense to override model/method names here?
    if args.model is not None:
        model_names = [args.model]
    if args.method is not None:
        elim_names = [args.method]

    if args.result_type != "plot-sweep":
        if args.by == "both":
            cl_lists = [
                ["python", scriptname, f"--model={mname}", f"--method={ename}"]
                for mname in model_names
                for ename in elim_names
            ]
        elif args.by == "model":
            cl = ["python", scriptname]
            if args.method is not None:
                cl.append(f"--method={args.method}")
            cl_lists = [cl + [f"--model={mname}"] for mname in model_names]
        elif args.by == "method":
            cl = ["python", scriptname]
            if args.model is not None:
                cl.append(f"--model={args.model}")
            cl_lists = [cl + [f"--method={ename}"] for ename in elim_names]

        if args.suffix is not None:
            for cl in cl_lists:
                cl.append(f"--suffix={args.suffix}")
        if args.results_dir != config.get_results_dir():
            for cl in cl_lists:
                cl.append(f"--results-dir={args.results_dir}")
        if args.tee:
            for cl in cl_lists:
                cl.append("--tee")

    else:
        results_dir = os.path.join(os.path.dirname(__file__), "results", "sweep")
        suff_str = "" if args.suffix is None else f"-{args.suffix}"
        results_fnames = [
            f"{mname}-{ename}-sweep{suff_str}.csv"
            for mname in model_names
            for ename in elim_names
        ]
        results_fpaths = [os.path.join(results_dir, fname) for fname in results_fnames]
        # TODO: Allow more arguments to propagate to the command lines we write.
        # Maybe this should be handled in its own script?
        cl_lists = [["python", scriptname, fpath] for fpath in results_fpaths]

    command_lines = [" ".join(cl) for cl in cl_lists]

    fname = args.result_type + "-commands.txt"
    fpath = os.path.join(config.validate_dir(args.commands_dir), fname)

    print(f"Writing the following commands to {fpath}")
    print()
    for cl in command_lines:
        print(cl)
    print()

    command_lines = [cl + "\n" for cl in command_lines]
    if not args.no_save:
        with open(fpath, "w") as f:
            f.write("".join(command_lines))


if __name__ == "__main__":
    argparser = config.get_argparser()
    argparser.add_argument("result_type", help="'structure', 'solvetime', or 'sweep'")
    argparser.add_argument(
        "--commands-dir",
        help="Directory to store files full of command lines",
        default=config.get_commands_dir(),
    )
    argparser.add_argument(
        "--by",
        help="Parallelize by model, method, or both. Default=both",
        default="both",
    )
    argparser.add_argument("--tee", action="store_true", help="Solver log to stdout")
    args = argparser.parse_args()
    main(args)
