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


def validate_scriptname(scriptname):
    if scriptname not in os.listdir(os.getcwd()):
        # Prepend FILEDIR relative to CWD
        relative_dirname = os.path.relpath(FILEDIR, os.getcwd())
        scriptname = os.path.join(relative_dirname, scriptname)
    return scriptname


def main(args):
    if args.model is None:
        mnames = list(config.TESTPROBLEM_LOOKUP.keys())
    else:
        if args.model not in config.TESTPROBLEM_LOOKUP:
            raise ValueError("--model must have a test problem")
        mnames = [args.model]
    if args.method is None:
        enames = config.ELIM_NAMES
    else:
        enames = [args.method]

    # NOTE: Assumming two parameters!
    NPARAM = 2
    nsamples_total = args.nsamples ** NPARAM

    parallel_sweep_commands = []
    collect_commands = []
    plot_commands = []
    suff_str = "" if args.suffix is None else f"-{args.suffix}"
    for mname in mnames:
        for ename in enames:
            # Here we validate that the script we are running is in the working directory
            scriptname = validate_scriptname("run_param_sweep.py")
            sample_commands = [
                [
                    "python",
                    scriptname,
                    f"--model={mname}",
                    f"--method={ename}",
                    f"--sample={i}",
                ]
                # Note that sample arg (and output filename) is base-1
                for i in range(1, nsamples_total + 1)
            ]
            for cmd in sample_commands:
                if args.suffix is not None:
                    cmd.append(f"--suffix={args.suffix}")
                if args.results_dir is not None:
                    cmd.append(f"--results-dir={args.results_dir}")

            sample_commands_str = [" ".join(cmd) + "\n" for cmd in sample_commands]

            sweep_command_fname = f"sweep-commands-{mname}-{ename}{suff_str}.txt"
            sweep_command_fpath = os.path.join(args.commands_dir, sweep_command_fname)

            if not args.no_save:
                print(f"Writing sample commands for {mname}-{ename} to {sweep_command_fpath}")
                with open(sweep_command_fpath, "w") as f:
                    for cmd in sample_commands_str:
                        f.write(cmd)

            parallel_sweep_commands.append(["parallel", "-a", sweep_command_fpath])
            scriptname = validate_scriptname("collect_sweep_results.py")
            collect_cmd = [
                "python",
                scriptname,
                f"--model={mname}",
                f"--method={ename}",
            ]
            if args.suffix is not None:
                collect_cmd.append(f"--suffix={args.suffix}")
            if args.output_suffix is not None:
                # Here, we potentially send both --suffix and --output-suffix.
                # suffix will be used to identify the right input files, while
                # output-suffix will be used for the output file.
                collect_cmd.append(f"--output_suffix={args.output_suffix}")
            if args.results_dir is not None:
                collect_cmd.append(f"--results-dir={args.results_dir}")
            collect_commands.append(collect_cmd)

            # TODO: Maybe send other arguments to the plotting command?
            # Here, --output-suffix overrides suffix, as it does in the output file
            # from collect_sweep_results
            output_suffix_str = suff_str if args.output_suffix is None else f"-{args.output_suffix}"
            result_fname = f"{mname}-{ename}-sweep{output_suffix_str}.csv"
            # Why did I used to hardcode this sweep directory? the collect script, above, used
            # a different results_dir then hard-coded the sweep subdirectory?
            sweepresultsdir = args.results_dir if args.results_dir.endswith("sweep") else os.path.join(args.results_dir, "sweep")
            result_fpath = os.path.join(sweepresultsdir, result_fname)
            scriptname = validate_scriptname("plot_sweep_results.py")
            plot_cmd = [
                "python",
                scriptname,
                result_fpath,
                f"--image-dir={args.image_dir}",
            ]
            plot_commands.append(plot_cmd)

    parallel_sweep_commands_str = [" ".join(cmd) + "\n" for cmd in parallel_sweep_commands]
    collect_commands_str = [" ".join(cmd) + "\n" for cmd in collect_commands]
    plot_commands_str = [" ".join(cmd) + "\n" for cmd in plot_commands]
    collect_commands_fname = "collect-sweep-commands"
    plot_commands_fname = "plot-sweep-commands"

    def update_fname_from_args(args, fname):
        if args.model is not None:
            fname += f"-{args.model}"
        if args.method is not None:
            fname += f"-{args.method}"
        if args.suffix is not None:
            fname += f"-{args.suffix}"
        return fname + ".txt"

    collect_commands_fname = update_fname_from_args(args, collect_commands_fname)
    collect_commands_fpath = os.path.join(args.commands_dir, collect_commands_fname)
    plot_commands_fname = update_fname_from_args(args, plot_commands_fname)
    plot_commands_fpath = os.path.join(args.commands_dir, plot_commands_fname)

    if not args.no_save:
        print(f"Writing collect-sweep-results commands to {collect_commands_fpath}")
        with open(collect_commands_fpath, "w") as f:
            for cmd in collect_commands_str:
                f.write(cmd)
        print(f"Writing plot-sweep commands to {plot_commands_fpath}")
        with open(plot_commands_fpath, "w") as f:
            for cmd in plot_commands_str:
                f.write(cmd)

    if args.model is None or args.method is None:
        # We only write this file full of commands if we are looping over something,
        # so we don't pollute our directory with useless files.
        parallel_commands_fname = "parallel-sweep-commands"
        parallel_commands_fname = update_fname_from_args(args, parallel_commands_fname)
        parallel_commands_fpath = os.path.join(args.commands_dir, parallel_commands_fname)

        if not args.no_save:
            print(f"Writing `parallel` commands for parameter sweeps to {parallel_commands_fpath}")
            with open(parallel_commands_fpath, "w") as f:
                for cmd in parallel_sweep_commands_str:
                    f.write(cmd)


if __name__ == "__main__":
    argparser = config.get_sweep_argparser()
    argparser.add_argument("--commands-dir", default=None)
    argparser.add_argument("--image-dir", default=None)
    argparser.add_argument(
        "--output-suffix",
        help=(
            "Suffix to apply to the collected result file. This is useful when we don't"
            " want to overwrite an existing result file, e.g. that was run in serial"
        ),
        default=None,
    )
    args = argparser.parse_args()
    if args.commands_dir is None:
        args.commands_dir = config.get_commands_dir()
    if args.image_dir is None:
        args.image_dir = config.get_image_dir()
    main(args)
