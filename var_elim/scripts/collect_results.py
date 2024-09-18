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
import pandas as pd
import var_elim.scripts.config as config


def main(args):
    if args.result_type != "structure" and args.result_type != "solvetime":
        raise ValueError("result_type must be 'structure' or 'solvetime'")
    else:
        basename = args.result_type

    if args.by == "both":
        model_names = config.MODEL_NAMES
        elim_names = config.ELIM_NAMES
        if args.model is not None:
            model_names = [args.model]
        if args.method is not None:
            elim_names = [args.method]
        suff_iter = itertools.product(model_names, elim_names)
        suffixes = ["-".join(elim_model_name) for elim_model_name in suff_iter]
    elif args.by == "model":
        if args.model is not None or args.method is not None:
            raise RuntimeError("--model and --method require --by=both")
        suffixes = list(config.MODEL_NAMES)
    elif args.by == "method":
        if args.model is not None or args.method is not None:
            raise RuntimeError("--model and --method require --by=both")
        suffixes = list(config.ELIM_NAMES)

    if args.suffix is not None:
        suffixes = [s + f"-{args.suffix}" for s in suffixes]

    fnames = [basename + "-" + s + ".csv" for s in suffixes]
    filedir = os.path.dirname(__file__)
    fpaths = [
        os.path.join(filedir, "results", args.result_type, fname)
        for fname in fnames
    ]
    print()
    print("Collecting results from files:")
    print()
    for fpath in fpaths:
        print(fpath)
    print()

    suff_str = "" if args.suffix is None else f"-{args.suffix}"
    output_fname = args.result_type + suff_str + ".csv"
    output_fpath = os.path.join(config.get_results_dir(), output_fname)

    df = pd.read_csv(files[0])
    for fname in files[1:]:
        file_df = pd.read_csv(fname)
        # - make sure dataframes have same columns
        # - concatenate rows of dataframes
    #df.to_csv(output_fname)


if __name__ == "__main__":
    argparser = config.get_argparser()
    argparser.add_argument("result_type", help="'structure' or 'solvetime'")
    argparser.add_argument(
        "--by",
        help="Filenames specified by model, method, or both",
        default="both",
    )

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
            resdir = os.path.join(config.get_results_dir(), "solvetime")
            config.validate_dir(resdir)
            args.results_dir = resdir

    main(args)
