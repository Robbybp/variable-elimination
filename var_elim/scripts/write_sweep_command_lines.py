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
    if args.model is None:
        mnames = config.MODEL_NAMES
    else:
        mnames = [args.model]
    if args.method is None:
        enames = config.ELIM_NAMES
    else:
        enames = [args.method]

    # NOTE: Assumming two parameters!
    NPARAM = 2
    nsamples_total = args.nsamples ** NPARAM

    parallel_sweep_commands = []
    suff_str = "" if args.suffix is None else f"-{args.suffix}"
    for mname in mnames:
        for ename in enames:
            # TODO: Write sweep commands to sweep-commands-mname-ename-suffix.txt
            sample_commands = [
                [
                    "python",
                    "run_param_sweep.py",
                    f"--model={mname}",
                    f"--method={ename}",
                    f"--sample={i}",
                ]
                # Note that sample arg (and output filename) is base-1
                for i in range(1, nsamples_total + 1)
            ]
            if args.suffix is not None:
                for cmd in sample_commands:
                    cmd.append(f"--suffix={args.suffix}")

            sample_commands_str = [" ".join(cmd) + "\n" for cmd in sample_commands]

            sweep_command_fname = f"sweep-commands-{mname}-{ename}{suff_str}.txt"
            sweep_command_fpath = os.path.join(args.commands_dir, sweep_command_fname)

            if not args.no_save:
                print(f"Writing sample commands for {mname}-{ename} to {sweep_command_fpath}")
                with open(sweep_command_fpath, "w") as f:
                    for cmd in sample_commands_str:
                        f.write(cmd)

            parallel_sweep_commands.append(["parallel", "-a", sweep_command_fpath])

    if args.model is None or args.method is None:
        # We only write this file full of commands if we are looping over something,
        # so we don't pollute our directory with useless models.
        parallel_commands_fname = "parallel-sweep-commands"
        if args.model is not None:
            parallel_commands_fname += f"-{args.model}"
        if args.method is not None:
            parallel_commands_fname += f"-{args.method}"
        if args.suffix is not None:
            parallel_commands_fname += f"-{args.suffix}"
        parallel_commands_fname += ".txt"


if __name__ == "__main__":
    argparser = config.get_sweep_argparser()
    argparser.add_argument("--commands-dir", default=None)
    args = argparser.parse_args()
    if args.commands_dir is None:
        args.commands_dir = config.get_commands_dir()
    main(args)
