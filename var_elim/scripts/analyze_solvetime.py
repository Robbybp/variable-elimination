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


def main():
    horizon = 300
    nfe = 300
    models = [
        ("Distill", lambda : create_distill(horizon=horizon, nfe=nfe)),
        #("OPF-4917", create_opf),
        #("Pipeline", create_pipeline),
    ]
    #m1 = its.get_problem("MBCLC-METHANE-STEADY").create_instance()
    #m2 = its.get_problem("MBCLC-METHANE-STEADY").create_instance()

    elim_callbacks = config.ELIM_CALLBACKS

    solvers = []
    # Note that these options may not be applicable for all models.
    # (i.e. 100 iterations may not be sufficient to solve)
    options = {
        "print_user_options": "yes",
        "max_iter": 100,
    }

    from var_elim.cyipopt import TimedPyomoCyIpoptSolver, Callback

    # Note that we allow any solver to be used here. However, below we assume
    # that the HierarchicalTimer will get populated with timing categories that
    # are specific to TimedPyomoCyIpoptSolver
    solvers.append(TimedPyomoCyIpoptSolver(options=options))

    model_cb_elim_solver_prod = list(itertools.product(models, elim_callbacks, solvers))

    model_elim_solver_prod = []
    for (mname, model_cb), (ename, elim_cb), solver in model_cb_elim_solver_prod:
        model = model_cb()
        model_elim_solver_prod.append(((mname, model_cb, model), (ename, elim_cb), solver))

    timer = TicTocTimer()
    htimer = HierarchicalTimer()
    htimer.start("root")
    for (mname, model_cb, model), (elim_name, elim_callback), solver in model_elim_solver_prod:
        print(mname, elim_name)
        timer.tic()
        elim_res = elim_callback(model)
        timer.toc("Apply elimination")
        
        label = "-".join((mname, elim_name))
        htimer.start(label)
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
        except RuntimeError as err:
            print(err)
            print("ERROR: BAD SOLVER STATUS. CONTINUING ANYWAY.")
        valid, violations = validate_solution(
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
        n_iter = solver.config.intermediate_callback.iterate_data[-1][1] + 1
        solve_time = htimer.timers["root"].timers[label].timers["solve"].total_time
        ls_trials = [data[-1] for data in solver.config.intermediate_callback.iterate_data]
        ave_ls_trials = sum(ls_trials)/len(ls_trials)
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
        total_solve_time = solve_timer.total_time
        other_solve_time = total_solve_time
        for subtimer in solve_timer.timers.values():
            other_solve_time -= subtimer.total_time
        other_solve_time_periter = other_solve_time / n_iter
        print(f"Other time per 100 iterations:     {100*other_solve_time_periter}")
        print()

    htimer.stop("root")
    print(htimer)


if __name__ == "__main__":
    main()
