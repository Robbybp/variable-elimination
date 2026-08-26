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

"""Write one parameter sweep command per model-method pair for a single solver

Each command runs a full sweep, so that it may be submitted to its own node. This
is preferable to splitting a sweep across the cores of a node with GNU parallel
when the solver is itself multithreaded, e.g. Gurobi.

Commands are written to {solver}-sweep-commands.txt, and results are written to a
per-solver subdirectory, so that this script may be called once per solver without
overwriting anything.

"""

import argparse
import os

import var_elim.scripts.config as config

FILEDIR = os.path.dirname(__file__)

SWEEP_SCRIPT = "run_param_sweep.py"


def get_script_path(scriptname):
    """Path of a script, relative to the current working directory"""
    if scriptname in os.listdir(os.getcwd()):
        return scriptname
    return os.path.join(os.path.relpath(FILEDIR, os.getcwd()), scriptname)


def get_argparser():
    argparser = argparse.ArgumentParser(description=__doc__)
    solver_list_str = ", ".join(config.SOLVER_NAMES)
    argparser.add_argument(
        "--solver",
        default="ipopt",
        help=f"Solver preset to use. Options are: {solver_list_str}",
    )
    argparser.add_argument(
        "--solver-options",
        default=None,
        help=(
            "Comma-separated solver options, e.g. Presolve=0,ScaleFlag=2. These"
            " override the options of the --solver preset."
        ),
    )
    model_list_str = ", ".join(config.TESTPROBLEM_LOOKUP)
    argparser.add_argument(
        "--models",
        default=None,
        help=(
            f"Comma-separated models to sweep. Options are: {model_list_str}."
            " Default is all of them."
        ),
    )
    elim_list_str = ", ".join(config.ELIM_NAMES)
    argparser.add_argument(
        "--methods",
        default=None,
        help=(
            f"Comma-separated elimination methods to sweep. Options are:"
            f" {elim_list_str}. Default is all of them."
        ),
    )
    argparser.add_argument(
        "--nsamples",
        type=int,
        default=11,
        help="Number of samples per parameter in each parameter sweep",
    )
    argparser.add_argument(
        "--results-dir",
        default=None,
        help=(
            "Base directory for sweep results. The commands write to the"
            " SOLVER subdirectory of this directory."
        ),
    )
    argparser.add_argument(
        "--commands-dir",
        default=None,
        help="Directory to store the file of command lines",
    )
    argparser.add_argument(
        "--feastol",
        type=float,
        default=None,
        help="Tolerance for checking feasibility",
    )
    argparser.add_argument(
        "--suffix",
        default=None,
        help="Suffix to append to result file names",
    )
    argparser.add_argument(
        "--no-save",
        action="store_true",
        help="Print the commands without writing them to a file",
    )
    return argparser


def main(args):
    if args.solver not in config.SOLVER_LOOKUP:
        raise ValueError(
            f"Unrecognized solver '{args.solver}'. Options are:"
            f" {', '.join(config.SOLVER_NAMES)}"
        )

    if args.models is None:
        mnames = list(config.TESTPROBLEM_LOOKUP)
    else:
        mnames = args.models.split(",")
        for mname in mnames:
            if mname not in config.TESTPROBLEM_LOOKUP:
                raise ValueError(f"Model '{mname}' does not have a test problem")

    if args.methods is None:
        enames = list(config.ELIM_NAMES)
    else:
        enames = args.methods.split(",")
        for ename in enames:
            if ename not in config.ELIM_LOOKUP:
                raise ValueError(
                    f"Unrecognized method '{ename}'. Options are:"
                    f" {', '.join(config.ELIM_NAMES)}"
                )

    # Results for different solvers are distinguished by their directory, so that
    # we don't have to encode the solver in every result file name.
    results_dir = os.path.join(args.results_dir, args.solver)

    scriptname = get_script_path(SWEEP_SCRIPT)
    command_lines = []
    for mname in mnames:
        for ename in enames:
            cmd = [
                "python",
                scriptname,
                f"--model={mname}",
                f"--method={ename}",
                f"--nsamples={args.nsamples}",
                f"--solver={args.solver}",
                f"--results-dir={results_dir}",
            ]
            if args.solver_options is not None:
                cmd.append(f"--solver-options={args.solver_options}")
            if args.feastol is not None:
                cmd.append(f"--feastol={args.feastol}")
            if args.suffix is not None:
                cmd.append(f"--suffix={args.suffix}")
            command_lines.append(" ".join(cmd))

    fname = f"{args.solver}-sweep-commands.txt"
    fpath = os.path.join(config.validate_dir(args.commands_dir), fname)

    print(f"Writing the following commands to {fpath}")
    print()
    for cl in command_lines:
        print(cl)
    print()
    print(f"Results will be written to {results_dir}")

    if not args.no_save:
        with open(fpath, "w") as f:
            f.write("".join(cl + "\n" for cl in command_lines))


if __name__ == "__main__":
    args = get_argparser().parse_args()
    if args.results_dir is None:
        args.results_dir = os.path.join(config.get_results_dir(), "sweep")
    if args.commands_dir is None:
        args.commands_dir = config.get_commands_dir()
    main(args)
