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
from pyomo.common.timing import HierarchicalTimer
import pselib
from idaes.core.util.parameter_sweep import (
    ParameterSweepSpecification,
    SequentialSweepRunner,
    UniformSampling,
)
import pandas as pd

import json
from collections import namedtuple

from var_elim.algorithms.validate import validate_solution
from var_elim.elimination_callbacks import (
    no_elim_callback,
    trivial_elim_callback,
    linear_d2_elim_callback,
    d2_elim_callback,
    matching_elim_callback,
)
from var_elim.models.testproblems import (
    DistillationTestProblem,
)


import os
import var_elim.scripts.config as config

SolveResults = namedtuple(
    "SolveResults",
    ["objective", "feasible", "timer"],
)

#N_SAMPLES = 11

testset = pselib.TestSet()


def main(args):

    elimination_callbacks = [
        ("no-elim", no_elim_callback),
        ("trivial-elim", trivial_elim_callback),
        ("linear-d2", linear_d2_elim_callback),
        ("d2", d2_elim_callback),
        ("matching", matching_elim_callback),
    ]

    # Instead of defining model construction callbacks, we will define TestProblems,
    # which define parameters and ranges as well as construction callbacks
    problems = [
        #("mb-steady", pselib.get_problem("MBCLC-METHANE-STEADY")),
        ("distill", DistillationTestProblem()),
    ]
    scale_problem = dict([("mb-steady", True), ("distill", False)])

    # We will store the results in a dict mapping callback name and problem
    # name to the results for a specific parameter sweep. We do this instead
    # of storing the results in a flattened list so sweep result files for
    # individual problem-method pairs are easier to generate.
    sweep_results_lookup = {}

    for problem_name, problem in problems:
        for elim_name, elim_cb in elimination_callbacks:
            print(f"Performing parameter sweep for problem: {problem.uid}")

            sweep = ParameterSweepSpecification()
            sweep.set_sampling_method(UniformSampling)
            sample_sizes = []
            for param in problem.parameters:
                sample_sizes.append(args.nsamples)
                lb, ub = problem.parameter_ranges[param]
                sweep.add_sampled_input(
                    str(param),
                    pyo.value(lb),
                    pyo.value(ub),
                )
            sweep.set_sample_size(sample_sizes)
            sweep.generate_samples()

            print("Samples:")
            print(sweep.samples)

            # Why is this necessary? It seems like a reasonable default could
            # just return results
            def build_outputs(model, results):
                # This is what we see in the runner.results["results"] column
                return results

            options = {
                "print_user_options": "yes",
                "max_iter": 1000,
            }
            # Note that if we want to collect more detailed information, we could
            # use TimedCyIpoptSolver
            solver = pyo.SolverFactory("cyipopt", options=options)

            def run_model(model, solver):
                timer = HierarchicalTimer()

                # Optional: re-initialize
                # model.fs.moving_bed.initialize()
                if scale_problem[problem_name]:
                    timer.start("scale")
                    # We scale the model here so we don't have to convert parameters into
                    # the scaled space (parameters are set before this function is called).
                    # This is really inconvenient... TODO: Re-consider scaling
                    pyo.TransformationFactory("core.scale_model").apply_to(model)
                    timer.stop("scale")

                # Run the elimination method
                timer.start("elimination")
                elim_results = elim_cb(model)
                timer.stop("elimination")

                timer.start("solve")
                results = solver.solve(model, tee=True)
                solved = pyo.check_optimal_termination(results)
                timer.stop("solve")

                timer.start("validate")
                # TODO: Should we exit if the solver didn't converge optimal?
                valid, violations = validate_solution(
                    model,
                    elim_results.var_expressions,
                    elim_results.constraints,
                    tolerance=1e-5,
                )
                timer.stop("validate")

                obj = next(iter(model.component_data_objects(pyo.Objective, active=True)))

                # Why is "solved" returned separately from the rest of the results
                results = SolveResults(pyo.value(obj), valid, timer)

                # Should runner check that solved is bool?
                return solved, results

            runner = SequentialSweepRunner(
                build_model=problem.create_instance,
                run_model=run_model,
                input_specification=sweep,
                # If we are implementing a run_model method, there is no reason to
                # provide a solver (other than that the solver is then global wrt
                # the run_model function)
                solver=solver,
                build_outputs=build_outputs,
            )

            runner.execute_parameter_sweep()

            sweep_results_lookup[problem_name, elim_name] = (
                sweep.samples, runner.results
            )

    n_converged_lookup = {}
    for problem_name, problem in problems:
        for elim_name, _ in elimination_callbacks:
            samples, results = sweep_results_lookup[problem_name, elim_name]
            results_df = pd.DataFrame(results).transpose()

            n_converged_lookup[problem_name, elim_name] = list(results_df["success"]).count(True)

            sweep_data_df = {}
            for param in problem.parameters:
                sweep_data_df[str(param)] = samples[str(param)]
            sweep_data_df["success"] = results_df["success"]
            sweep_data_df["error"] = results_df["error"]
            sweep_data_df["feasible"] = [res.feasible if res is not None else False for res in results_df["results"]]
            sweep_data_df["objective"] = [res.objective if res is not None else None for res in results_df["results"]]
            sweep_data_df["solve-time"] = [res.timer.timers["solve"].total_time if res is not None else None for res in results_df["results"]]
            sweep_data_df["elim-time"] = [res.timer.timers["elimination"].total_time if res is not None else None for res in results_df["results"]]
            # TODO: Add values for degree-of-freedom variables?
            # What do I want these for? A later evaluation step?

            sweep_data_df = pd.DataFrame(sweep_data_df)

            print(f"Sweep data for {problem_name}-{elim_name}:")
            print(sweep_data_df)

            fname = f"{problem_name}-{elim_name}-sweep.csv"
            fpath = os.path.join(args.results_dir, fname)
            sweep_data_df.to_csv(fpath)

    for problem_name, problem in problems:
        n_instances = N_SAMPLES ** len(problem.parameters)
        for elim_name, _ in elimination_callbacks:
            n_converged = n_converged_lookup[problem_name, elim_name]
            print(f"{problem_name}-{elim_name} converged {n_converged} / {n_instances} instances")


if __name__ == "__main__":
    argparser = config.get_sweep_argparser()
    args = argparser.parse_args()
    main(args)
