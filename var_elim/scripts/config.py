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
import argparse
from var_elim.elimination_callbacks import (
    no_elim_callback,
    linear_degree1_elim_callback,
    equalcoef_degree2_elim_callback,
    linear_degree2_elim_callback,
    degree2_elim_callback,
    matching_elim_callback,
    greedy_elim_callback,
)
from var_elim.models.testproblems import (
    DistillationTestProblem,
    PipelineTestProblem
)
from var_elim.models.opf.opf_model import make_model as create_opf
from var_elim.models.gas_pipelines.gas_network_model import make_dynamic_model as create_pipeline
import pselib
from var_elim.cyipopt import TimedPyomoCyIpoptSolver, Callback
import pyomo.environ as pyo

filedir = os.path.dirname(__file__)

ELIM_CALLBACKS = [
    # These names are used in the "method" column in dataframes, or in the
    # filename of parameter sweep results.
    ("no-elim", no_elim_callback),
    ("d1", linear_degree1_elim_callback),
    ("ecd2", equalcoef_degree2_elim_callback),
    ("linear-d2", linear_degree2_elim_callback),
    ("d2", degree2_elim_callback),
    ("greedy", greedy_elim_callback),
    ("matching", matching_elim_callback),
]
ELIM_LOOKUP = dict(ELIM_CALLBACKS)
ELIM_NAMES = [name for name, _ in ELIM_CALLBACKS]

# Each preset maps a name to the Pyomo solver interface and the options it is
# constructed with. These names are used in the "solver" column in dataframes.
# Note that OptimalityTarget=1 is necessary for the Gurobi presets: without it,
# Gurobi attempts to solve these nonconvex problems to global optimality.
SOLVER_PRESETS = [
    ("ipopt",            ("cyipopt",             {"print_user_options": "yes", "max_iter": 3000})),
    ("ipopt-tol6",       ("ipopt",               {"print_user_options": "yes", "max_iter": 3000, "tol": 1e-6})),
    ("gurobi",           ("gurobi_direct_minlp", {"OptimalityTarget": 1, "NLBarIterLimit": 3000})),
    ("gurobi-Presolve0", ("gurobi_direct_minlp", {"OptimalityTarget": 1, "NLBarIterLimit": 3000, "Presolve": 0})),
    ("gurobi-Presolve1", ("gurobi_direct_minlp", {"OptimalityTarget": 1, "NLBarIterLimit": 3000, "Presolve": 1})),
    ("gurobi-Presolve2", ("gurobi_direct_minlp", {"OptimalityTarget": 1, "NLBarIterLimit": 3000, "Presolve": 2})),
]
SOLVER_LOOKUP = dict(SOLVER_PRESETS)
SOLVER_NAMES = [name for name, _ in SOLVER_PRESETS]


def _patch_gurobi_status_map():
    """Add the LOCALLY_* Gurobi statuses to Pyomo's termination condition map

    Pyomo's map predates these statuses, so they fall through to `unknown`. This
    makes every locally optimal solve look like a failure to
    pyo.check_optimal_termination. Note that the new-style TerminationCondition
    enum does not distinguish local from global optimality.

    This can be removed once Pyomo's map includes these statuses.

    """
    from gurobipy import GRB
    from pyomo.contrib.solver.common.results import TerminationCondition
    from pyomo.contrib.solver.solvers.gurobi.gurobi_direct import GurobiDirectBase

    # The map is built lazily and cached on the class, so we have to construct it
    # before we can update it.
    tc_map = GurobiDirectBase._get_tc_map(None)
    tc_map[GRB.LOCALLY_OPTIMAL] = TerminationCondition.convergenceCriteriaSatisfied
    tc_map[GRB.LOCALLY_INFEASIBLE] = TerminationCondition.locallyInfeasible


_patch_gurobi_status_map()

MODEL_NAMES = [
    # These names are used in the "model" column in dataframes, or in the
    # filename of parameter sweep results.
    "distill",
    "mb-steady",
    "opf",
    "pipeline",
]

def mb_steady_constructor():
    model = pselib.get_problem("MBCLC-METHANE-STEADY").create_instance()
    pyo.TransformationFactory("core.scale_model").apply_to(model)
    return model

# NOTE: Via the construction of CONSTRUCTOR_LOOKUP below, we rely on
# MODEL_CONSTRUCTORS and MODEL_NAMES having the same order.
MODEL_CONSTRUCTORS = [
    DistillationTestProblem().create_instance,
    # Note that the mb-steady constructor already scales, it *should not*
    # be used for parameter sweeps. This should be handled by the TestProblem
    # constructor below.
    mb_steady_constructor,
    create_opf,
    PipelineTestProblem().create_instance,
]
CONSTRUCTOR_LOOKUP = dict(zip(MODEL_NAMES, MODEL_CONSTRUCTORS))

TESTPROBLEM_LOOKUP = {
    "distill": DistillationTestProblem(),
    "mb-steady": pselib.get_problem("MBCLC-METHANE-STEADY"),
    "pipeline": PipelineTestProblem(),
}


def get_optimization_solver():
    options = {
        "print_user_options": "yes",
        "print_timing_statistics": "yes",
        "max_iter": 3000,
    }
    # Simple callback to get the number of iterations
    #callback = Callback()
    solver = TimedPyomoCyIpoptSolver(
        options=options,
        # Note that we don't set the intermediate callback here, as we don't want
        # to accidentally use it for multiple solves.
        # The alternative is just to construct a new solver each time it needs
        # to be used.
        #intermediate_callback=callback,
    )
    return solver


def parse_solver_options(options_str):
    """Parse a comma-separated string of option=value pairs into a dict

    Values are left as strings. Both the Gurobi and Ipopt interfaces accept
    string-valued options, so we don't need to know each option's type.

    """
    if options_str is None:
        return {}
    options = {}
    for pair in options_str.split(","):
        key, sep, value = pair.partition("=")
        if not sep:
            raise ValueError(
                f"Could not parse solver option '{pair}'. Options must be"
                " comma-separated key=value pairs, e.g. Presolve=0,ScaleFlag=2"
            )
        options[key.strip()] = value.strip()
    return options


def get_solver(solvername, options_str=None):
    """Construct a solver from a preset, overriding options from a string

    Arguments
    ---------
    solvername: Key of SOLVER_LOOKUP
    options_str: Optional comma-separated string of option=value pairs. These
        override the options of the preset.

    """
    if solvername not in SOLVER_LOOKUP:
        raise ValueError(
            f"Unrecognized solver '{solvername}'. Options are:"
            f" {', '.join(SOLVER_NAMES)}"
        )
    interface, options = SOLVER_LOOKUP[solvername]
    # Don't modify the preset's dict
    options = dict(options)
    options.update(parse_solver_options(options_str))
    return pyo.SolverFactory(interface, options=options)


def get_basename(basename, *suffixes):
    if not basename.endswith(".csv") and not basename.endswith(".CSV"):
        raise ValueError("basename must end with .csv or .CSV")
    name_noext = basename[:-4]
    ext = basename[-4:]
    suffix = "".join(["-" + s for s in suffixes if s is not None])
    name = name_noext + suffix + ext
    return name


def validate_dir(name):
    if os.path.isfile(name):
        raise OSError(f"Default directory {name} is already a file")
    elif not os.path.isdir(name):
        os.makedirs(name)
    return name


def get_results_dir():
    resdir = os.path.join(filedir, "results")
    validate_dir(resdir)
    return resdir

def get_commands_dir():
    comdir = os.path.join(filedir, "commands")
    validate_dir(comdir)
    return comdir

def get_image_dir():
    plotdir = os.path.join(filedir, "images")
    validate_dir(plotdir)
    return plotdir

def get_argparser():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--results-dir", default=get_results_dir())
    model_list_str = ", ".join(MODEL_NAMES)
    argparser.add_argument(
        "--model",
        default=None,
        help=(
            "Model to analyze or use for parameter sweep (optional). Options are: "
            + model_list_str
        ),
    )
    elim_list_str = ", ".join(ELIM_NAMES)
    argparser.add_argument(
        "--method",
        default=None,
        help=(
            f"Method to apply. Options are: {elim_list_str}. Default (None) loops"
            " over all methods."
        ),
    )
    solver_list_str = ", ".join(SOLVER_NAMES)
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
    argparser.add_argument("--no-save", action="store_true", help="Don't save results")
    argparser.add_argument(
        "--suffix",
        default=None,
        help=(
            "Suffix to append to result file names. This is useful for storing"
            " preliminary results created with different dependency versions."
        ),
    )
    argparser.add_argument("--feastol", type=float, default=1e-4, help="Tolerance for checking feasibility")
    return argparser

def get_sweep_argparser():
    argparser = get_argparser()
    argparser.add_argument(
        "--nsamples",
        type=int,
        default=11,
        help="Number of samples per parameter in parameter sweep"
    )
    argparser.add_argument(
        "--sample",
        type=int,
        default=None,
        help="Index (base-1) of sample to run (default, None, runs all samples)",
    )
    return argparser

def get_plot_argparser():
    argparser = get_argparser()
    argparser.add_argument("--image-dir", default=get_image_dir())
    argparser.add_argument("--show", action="store_true")
    argparser.add_argument("--opaque", action="store_true")
    argparser.add_argument("--title", default=None)
    argparser.add_argument("--no-legend", action="store_true")
    return argparser
