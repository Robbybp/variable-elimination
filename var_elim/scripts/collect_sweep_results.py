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
    if args.model is None or args.method is None:
        raise RuntimeError("--model and --method must be provided")

    # NOTE: We assume that we are sweeping over two parameters
    nparameters = 2
    ninstances = args.nsgreedyes ** nparameters

    fnames = [
        f"{args.model}-{args.method}-sweep-{i}of{ninstances}"
        for i in range(1, ninstances + 1)
    ]
    if args.suffix is not None:
        fnames = [fname + f"-{args.suffix}" for fname in fnames]
    fnames = [fname + ".csv" for fname in fnames]

    filedir = os.path.dirname(__file__)
    fpaths = [os.path.join(args.results_dir, fname) for fname in fnames]
    print()
    print("Collecting results from files:")
    print()
    for fpath in fpaths:
        print(fpath)
    print()

    valid_fpaths = []
    for fpath in fpaths:
        if os.path.isfile(fpath):
            valid_fpaths.append(fpath)
        else:
            print(f"WARNING: {fpath} does not exist or is not a file. Skipping...")
    fpaths = valid_fpaths

    suff_str = ""
    if args.suffix is not None:
        suff_str = f"-{args.suffix}"
    if args.output_suffix is not None:
        # output-suffix overrides suffix
        suff_str = f"-{args.output_suffix}"
    output_fname = f"{args.model}-{args.method}-sweep" + suff_str + ".csv"
    output_fpath = os.path.join(args.results_dir, output_fname)

    dfs = [pd.read_csv(fpath) for fpath in fpaths]
    output_df = pd.concat(dfs, ignore_index=True, join="inner")
    if not all(all(output_df.columns == df.columns) for df in dfs):
        raise RuntimeError(
            "Columns changed as a result of inner join."
            " At least one dataframe is missing a column"
        )
    if not args.no_save:
        print(f"Writing combined dataframe to {output_fpath}")
        output_df.to_csv(output_fpath)


if __name__ == "__main__":
    argparser = config.get_sweep_argparser()
    argparser.add_argument(
        "--output-suffix",
        help=(
            "Suffix to apply to the collected result file. This is useful when we don't"
            " want to overwrite an existing result file, e.g. that was run in serial"
        ),
        default=None,
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
            resdir = os.path.join(config.get_results_dir(), "sweep")
            config.validate_dir(resdir)
            args.results_dir = resdir

    main(args)
