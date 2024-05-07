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
    d1_elim_callback,
    trivial_elim_callback,
    linear_d2_elim_callback,
    d2_elim_callback,
    matching_elim_callback,
)
from var_elim.models.testproblems import (
    DistillationTestProblem,
)
from var_elim.models.opf.opf_model import make_model as create_opf
from var_elim.models.gas_pipelines.gas_network_model import make_dynamic_model as create_pipeline
import pselib
from var_elim.cyipopt import TimedPyomoCyIpoptSolver, Callback

filedir = os.path.dirname(__file__)

ELIM_CALLBACKS = [
    # These names are used in the "method" column in dataframes, or in the
    # filename of parameter sweep results.
    ("no-elim", no_elim_callback),
    ("d1", d1_elim_callback),
    ("trivial", trivial_elim_callback),
    ("linear-d2", linear_d2_elim_callback),
    ("d2", d2_elim_callback),
    ("matching", matching_elim_callback),
]

MODEL_NAMES = [
    # These names are used in the "model" column in dataframes, or in the
    # filename of parameter sweep results.
    "distill",
    "mb-steady",
    "opf",
    "pipeline",
]

MODEL_CONSTRUCTORS = [
    DistillationTestProblem().create_instance,
    pselib.get_problem("MBCLC-METHANE-STEADY").create_instance,
    create_opf,
    create_pipeline,
]

CONSTRUCTOR_LOOKUP = dict(zip(MODEL_NAMES, MODEL_CONSTRUCTORS))

TESTPROBLEM_LOOKUP = {
    "distill": DistillationTestProblem(),
    "mb-steady": pselib.get_problem("MBCLC-METHANE-STEADY"),
}


def get_optimization_solver():
    options = {
        "print_user_options": "yes",
        "max_iter": 1000,
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


def get_results_dir():
    resdir = os.path.join(filedir, "results")
    if os.path.isfile(resdir):
        raise OSError(f"Default results directory {resdir} is already a file")
    elif not os.path.isdir(resdir):
        os.mkdir(resdir)
    return resdir

def get_image_dir():
    plotdir = os.path.join(filedir, "images")
    if os.path.isfile(plotdir):
        raise OSError(f"Default image directory {plotdir} is already a file")
    elif not os.path.isdir(plotdir):
        os.mkdir(plotdir)
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
    argparser.add_argument("--no-save", action="store_true", help="Don't save results")
    # TODO: feastol argument, include max infeasibility in parameter sweep results.
    #argparser.add_argument("--feastol", type=float, default=1e-5, help="Tolerance for checking feasibility")
    return argparser

def get_sweep_argparser():
    argparser = get_argparser()
    argparser.add_argument(
        "--nsamples",
        type=int,
        default=11,
        help="Number of samples per parameter in parameter sweep"
    )
    return argparser

def get_plot_argparser():
    argparser = get_argparser()
    argparser.add_argument("--image-dir", default=get_image_dir())
    argparser.add_argument("--show", action="store_true")
    argparser.add_argument("--no-save", action="store_true")
    argparser.add_argument("--opaque", action="store_true")
    argparser.add_argument("--title", default=None)
    return argparser
