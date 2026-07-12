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

import argparse
import datetime
import os
import shlex
import subprocess


def default_run_dir():
    date = datetime.date.today().strftime("%Y%m%d")
    return os.path.join("runs", date)


def script_path(script_name):
    return os.path.join("var_elim", "scripts", script_name)


def require_repo_root():
    expected = [
        "pyproject.toml",
        os.path.join("var_elim", "scripts"),
    ]
    missing = [path for path in expected if not os.path.exists(path)]
    if missing:
        raise RuntimeError(
            "reproduce.py must be run from the repository root. "
            f"Missing expected path: {missing[0]}"
        )


def run_command(cmd, dry_run=False):
    print("+ " + shlex.join(cmd))
    if not dry_run:
        subprocess.run(cmd, check=True)


def make_dir(dirname, dry_run=False):
    print("+ mkdir -p " + shlex.quote(dirname))
    if not dry_run:
        os.makedirs(dirname, exist_ok=True)


def add_common_args(parser):
    parser.add_argument(
        "--run-dir",
        default=None,
        help="Run directory. Defaults to runs/YYYYMMDD.",
    )
    parser.add_argument(
        "--results-dir",
        default=None,
        help="Directory for CSV and table outputs. Defaults to RUN_DIR/results.",
    )
    parser.add_argument(
        "--image-dir",
        default=None,
        help="Directory for plot outputs. Defaults to RUN_DIR/images.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running them.",
    )


def add_structure_args(parser):
    add_common_args(parser)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run only the distillation structural workflow for a quicker check.",
    )
    parser.add_argument(
        "--skip-tables",
        action="store_true",
        help="Skip LaTeX table generation.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip plot generation.",
    )


def resolve_output_dirs(args):
    run_dir = args.run_dir if args.run_dir is not None else default_run_dir()
    results_dir = (
        args.results_dir
        if args.results_dir is not None
        else os.path.join(run_dir, "results")
    )
    image_dir = (
        args.image_dir
        if args.image_dir is not None
        else os.path.join(run_dir, "images")
    )
    return run_dir, results_dir, image_dir


def run_structure(args):
    _, results_dir, image_dir = resolve_output_dirs(args)

    make_dir(results_dir, dry_run=args.dry_run)
    make_dir(image_dir, dry_run=args.dry_run)

    analyze_cmd = [
        "python",
        script_path("analyze_structure.py"),
        "--results-dir",
        results_dir,
    ]
    if args.smoke:
        analyze_cmd.extend(["--model", "distill"])

    run_command(analyze_cmd, dry_run=args.dry_run)

    csv_name = "structure-distill.csv" if args.smoke else "structure.csv"
    structure_csv = os.path.join(results_dir, csv_name)

    if not args.skip_tables:
        table_cmd = [
            "python",
            script_path("write_latex_table.py"),
            structure_csv,
            "--results-dir",
            results_dir,
        ]
        run_command(table_cmd, dry_run=args.dry_run)

        matching_table_cmd = table_cmd + ["--which", "matching-bounds"]
        run_command(matching_table_cmd, dry_run=args.dry_run)

    if not args.skip_plots:
        plot_cmd = [
            "python",
            script_path("plot_structure_bargraphs.py"),
            structure_csv,
            "--image-dir",
            image_dir,
        ]
        run_command(plot_cmd, dry_run=args.dry_run)

        if not args.smoke:
            sparsity_cmd = [
                "python",
                script_path("plot_sparsity.py"),
                "--model",
                "mb-steady",
                "--image-dir",
                image_dir,
            ]
            run_command(sparsity_cmd, dry_run=args.dry_run)


def main():
    require_repo_root()

    parser = argparse.ArgumentParser(
        description="Run end-to-end result reproduction workflows."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    structure_parser = subparsers.add_parser(
        "structure",
        help="Run structural analysis, tables, and plots.",
    )
    add_structure_args(structure_parser)
    structure_parser.set_defaults(func=run_structure)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
