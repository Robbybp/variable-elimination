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

import time
import pyomo.environ as pyo
from pyomo.common.collections import ComponentMap
from pyomo.common.timing import TicTocTimer, HierarchicalTimer
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.config import IncidenceMethod
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix

from collections import namedtuple
import itertools

import scipy.sparse as sps
import matplotlib.pyplot as plt

from var_elim.models.distillation.distill import create_instance as create_distill
from var_elim.models.opf.opf_model import make_model as create_opf
from var_elim.models.gas_pipelines.gas_network_model import make_dynamic_model as create_pipeline
from var_elim.heuristics.matching import (
    generate_elimination_via_matching,
    define_elimination_order,
)
from var_elim.heuristics.trivial_elimination import (
    get_degree_one_elimination,
    get_degree_two_elimination,
    get_linear_degree_two_elimination,
    get_trivial_constraint_elimination,
    filter_constraints,
)
from var_elim.algorithms.replace import eliminate_variables 
from var_elim.algorithms.expr import (
    count_nodes, count_model_nodes, count_amplrepn_nodes
)
from var_elim.algorithms.validate import validate_solution

import var_elim.scripts.config as config
from var_elim.elimination_callbacks import matching_elim_callback
# We need to import the callback here so we can re-construct it for every model
from var_elim.cyipopt import Callback

import os
import pandas as pd


USE_NAMED_EXPRESSIONS = True


def solve_original(m, tee=True):
    solver = pyo.SolverFactory("ipopt")
    #solver.options["print_timing_statistics"] = "yes"
    solver.solve(m, tee=tee)
    return m


def solve_reduced(m, tee=True, callback=matching_elim_callback):
    callback(m)
    solver = pyo.SolverFactory("ipopt")
    #solver.options["print_timing_statistics"] = "yes"
    solver.solve(m, tee=tee)
    return m


def main(args):
    if args.model is not None:
        models = [(args.model, config.CONSTRUCTOR_LOOKUP[args.model])]
    else:
        # Note that this is a hard-coded subset of models that are the default
        # to use for the timing analysis
        # TODO: Add 100-discr version of mb-steady
        model_names = ["distill", "mb-steady", "opf", "pipeline"]
        models = [(name, config.CONSTRUCTOR_LOOKUP[name]) for name in model_names]

    if args.method is None:
        elim_callbacks = config.ELIM_CALLBACKS
    else:
        elim_callbacks = [(args.method, config.ELIM_LOOKUP[args.method])]

    solvers = []
    solver = config.get_optimization_solver()
    solvers.append(solver)

    model_cb_elim_solver_prod = list(itertools.product(models, elim_callbacks, solvers))

    model_elim_solver_prod = []
    model_build_time_lookup = {}
    for (mname, model_cb), (ename, elim_cb), solver in model_cb_elim_solver_prod:
        _t = time.time()
        model = model_cb()
        model_build_time_lookup[mname, ename] = time.time() - _t
        model_elim_solver_prod.append(((mname, model_cb, model), (ename, elim_cb), solver))

    data = {
        "model": [],
        "method": [],
        "success": [],
        "feasible": [],
        "elim-time": [],
        "solve-time": [],
        "init-time": [],
        "build-time": [],
        # I don't see any reason to distinguish between e.g. constraint and objective
        # time here.
        "function-time": [],
        "jacobian-time": [],
        "hessian-time": [],
        "n-iter": [],
        "ave-ls-trials": [],
        "function-per100": [],
        "jacobian-per100": [],
        "hessian-per100": [],
        "other-per100": [],
    }

    timer = TicTocTimer()
    htimer = HierarchicalTimer()
    htimer.start("root")
    for (mname, model_cb, model), (elim_name, elim_callback), solver in model_elim_solver_prod:
        print(mname, elim_name)
        timer.tic()
        elim_res = elim_callback(model)
        elim_time = timer.toc("Apply elimination")

        label = "-".join((mname, elim_name))
        htimer.start(label)
        # We need to re-set the callback each time we solve a model
        solver.config.intermediate_callback = Callback()
        res = solver.solve(model, tee=False, timer=htimer)
        htimer.stop(label)

        timer.toc("Solve model")

        try:
            # For this particular experiment, it actually doesn't matter whether
            # we solve the problem. We are just interested in timing the iterations.
            # Obviously, all else equal, we would like to solve the problem (lest
            # these iterations be called unrepresentative), but we can still
            # collect what I believe to be meaningful data if we don't.
            #
            # To collect meaninful data on the per-iteration (percall) level,
            # all we need to do is take sufficiently many iterations that noise
            # is negligible.
            pyo.assert_optimal_termination(res)
            success = pyo.check_optimal_termination(res)
        except RuntimeError as err:
            print(err)
            print("ERROR: BAD SOLVER STATUS. CONTINUING ANYWAY.")
        valid, violations = validate_solution(
            # TODO: Use args.feastol here?
            model, elim_res.var_expressions, elim_res.constraints, tolerance=1e-6
        )
        if not valid:
            print("WARNING: Result is not valid!")
        timer.toc("Validate result")

        # Extract timing information from the hierarchical timer

        # We combine per-call function evaluation for the objective and
        # constraints.
        con_timer = htimer.timers["root"].timers[label].timers["solve"].timers["constraints"]
        obj_timer = htimer.timers["root"].timers[label].timers["solve"].timers["objective"]
        obj_percall = obj_timer.total_time / obj_timer.n_calls
        con_percall = con_timer.total_time / con_timer.n_calls
        function_eval_percall = obj_percall + con_percall
        function_eval_per100 = 100 * function_eval_percall

        # We combine per-call function evaluation for the objective gradient and
        # constraint Jacobian
        conjac_timer = htimer.timers["root"].timers[label].timers["solve"].timers["constraint-jac"]
        objgrad_timer = htimer.timers["root"].timers[label].timers["solve"].timers["objective-grad"]
        objgrad_percall = objgrad_timer.total_time / objgrad_timer.n_calls
        conjac_percall = conjac_timer.total_time / conjac_timer.n_calls
        jacobian_eval_percall = objgrad_percall + conjac_percall
        jacobian_eval_per100 = 100 * jacobian_eval_percall

        laghess_timer = htimer.timers["root"].timers[label].timers["solve"].timers["lagrangian-hess"]
        laghess_eval_percall = laghess_timer.total_time / laghess_timer.n_calls
        laghess_eval_per100 = 100 * laghess_eval_percall

        # Unclear exactly how to count iterations. Sometimes the same iteration
        # number comes up multiple times in the callback.
        #
        # This accesses the second element of the last tuple of iterate data.
        # We add 1 to account for iteration zero.
        # (Unclear if we should actually include iteration zero...)
        n_iter = solver.config.intermediate_callback.iterate_data[-1][1] + 1
        solve_time = htimer.timers["root"].timers[label].timers["solve"].total_time
        ls_trials = [data[-1] for data in solver.config.intermediate_callback.iterate_data]
        ave_ls_trials = sum(ls_trials)/len(ls_trials)
        print(f"Time to build the model was {model_build_time_lookup[mname, elim_name]} s")
        print(f"Ipopt took {n_iter} iterations (including iteration 0)")
        print(f"Ipopt took {solve_time} s to solve the problem")
        print(f"Ipopt took, on average, {ave_ls_trials} line search trials per iteration")

        # These three categories sum to the total evaluation time for an "ideal"
        # iteration. In reality, iterations can involve more constraint and objective
        # evaluations due to failed line search trials.
        print(f"Time per 100 function evaluations: {function_eval_per100}")
        print(f"Time per 100 Jacobian evaluations: {jacobian_eval_per100}")
        print(f"Time per 100 Hessian evaluations:  {laghess_eval_per100}")

        solve_timer = htimer.timers["root"].timers[label].timers["solve"]
        # This is really "time spent in Ipopt"
        total_solve_time = solve_timer.total_time
        # This is time in Ipopt, not measured by our callbacks. We assume this
        # is dominated by KKT matrix factorization.
        other_solve_time = total_solve_time
        for subtimer in solve_timer.timers.values():
            other_solve_time -= subtimer.total_time
        other_solve_time_periter = other_solve_time / n_iter
        other_solve_time_per100 = 100 * other_solve_time_periter
        print(f"Other time per 100 iterations:     {other_solve_time_per100}")
        # This is time spent in the solver.solve method, *not* spent in Ipopt.
        # This includes .nl file write and initialization of the PyomoNLP.
        init_time = htimer.timers["root"].timers[label].total_time - total_solve_time
        print(f"Time spent initializing solver:    {init_time}")
        print()

        data["model"].append(mname)
        data["method"].append(elim_name)
        data["elim-time"].append(elim_time)
        data["success"].append(success)
        data["feasible"].append(valid)
        data["solve-time"].append(solve_time)
        data["init-time"].append(init_time)
        # This is time to build the model, and has nothing to do with the elimination
        data["build-time"].append(model_build_time_lookup[mname, elim_name])
        function_time = con_timer.total_time + obj_timer.total_time
        data["function-time"].append(function_time)
        jacobian_time = conjac_timer.total_time + objgrad_timer.total_time
        data["jacobian-time"].append(jacobian_time)
        hessian_time = laghess_timer.total_time
        data["hessian-time"].append(hessian_time)
        data["n-iter"].append(n_iter)
        data["ave-ls-trials"].append(ave_ls_trials)
        data["function-per100"].append(function_eval_per100)
        data["jacobian-per100"].append(jacobian_eval_per100)
        data["hessian-per100"].append(laghess_eval_per100)
        data["other-per100"].append(other_solve_time_per100)

    htimer.stop("root")
    print(htimer)

    df = pd.DataFrame(data)
    suffixes = [args.model, args.method, args.suffix]
    fname = config.get_basename("solvetime.csv", *suffixes)
    fpath = os.path.join(args.results_dir, fname)
    print(df)
    if not args.no_save:
        fpath = os.path.join(args.results_dir, fname)
        print(f"Writing results to {fpath}")
        df.to_csv(fpath)


if __name__ == "__main__":
    argparser = config.get_argparser()
    argparser.add_argument(
        "--fname",
        default=None,
        help="Basename for file to write results to",
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
