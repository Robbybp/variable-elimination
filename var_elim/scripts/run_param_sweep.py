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
    ["objective", "feasible", "max_infeasibility", "timer"],
)

testset = pselib.TestSet()


def main(args):

    if args.method is None:
        elimination_callbacks = config.ELIM_CALLBACKS
    else:
        elimination_callbacks = [(args.method, config.ELIM_LOOKUP[args.method])]

    if args.model is None:
        raise RuntimeError("--model argument must be provided for parameter sweep")
    # For now, this script is only set up to run with a single test problem at a time
    problems = [(args.model, config.TESTPROBLEM_LOOKUP[args.model])]

    print()
    print("Running parameter sweep with following arguments")
    print(f"  model={args.method}")
    print(f"  method={args.method}")
    print(f"  nsamples={args.nsamples}")
    print(f"  sample={args.sample}")

    # Annoyingly, the MB problem requires scaling. Because the input parameters
    # will be scaled, we need to do this scaling *after* setting these input
    # parameters. We don't want to unnecessarily try to scale models if it is
    # not necessary (or actively harmful), so we set this flag if the model needs
    # to be scaled.
    # ... This should really be handled better in pselib...
    scale_problem = dict([("mb-steady", True), ("distill", False), ("pipeline", False)])

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

            print("Samples: (indices are base-0)")
            if args.sample is None:
                print(sweep.samples)
            else:
                print(pd.DataFrame(sweep.samples.loc[args.sample-1]).transpose())

            # Why is this necessary? It seems like a reasonable default could
            # just return results
            def build_outputs(model, results):
                # This is what we see in the runner.results["results"] column
                return results

            options = {
                "print_user_options": "yes",
                "max_iter": 3000,
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
                    tolerance=args.feastol,
                )
                viol_cons_reduced, viol_bounds, viol_elim_cons = violations
                max_con_viol = max(viol_cons_reduced, default=(0.0,), key=lambda item: abs(item[-1]))
                max_bound_viol = max(viol_bounds, default=(0.0,), key=lambda item: abs(item[-1]))
                max_elim_viol = max(viol_elim_cons, default=(0.0,), key=lambda item: abs(item[-1]))
                max_infeas = max(abs(max_con_viol[-1]), abs(max_bound_viol[-1]), abs(max_elim_viol[-1]))
                timer.stop("validate")

                obj = next(iter(model.component_data_objects(pyo.Objective, active=True)))

                # Why is "solved" returned separately from the rest of the results
                results = SolveResults(pyo.value(obj), valid, max_infeas, timer)
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

            if args.sample is None:
                # If a sample was not provided, run the entire sweep
                runner.execute_parameter_sweep()
                samples = sweep.samples
                results = runner.results
            else:
                # If a sample was provided, only run the sample.
                s = args.sample - 1
                sresults, success, error = runner.execute_single_sample(s)
                # Mock up our own samples/results data structures
                samples = pd.DataFrame(sweep.samples.loc[s]).transpose()
                results = {s: dict(results=sresults, success=success, error=error)}

            sweep_results_lookup[problem_name, elim_name] = (samples, results)

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
            sweep_data_df["max-infeasibility"] = [res.max_infeasibility if res is not None else float("nan") for res in results_df["results"]]
            sweep_data_df["objective"] = [res.objective if res is not None else None for res in results_df["results"]]
            sweep_data_df["solve-time"] = [res.timer.timers["solve"].total_time if res is not None else None for res in results_df["results"]]
            sweep_data_df["elim-time"] = [res.timer.timers["elimination"].total_time if res is not None else None for res in results_df["results"]]
            # TODO: Add values for degree-of-freedom variables?
            # What do I want these for? A later evaluation step?

            sweep_data_df = pd.DataFrame(sweep_data_df)

            print(f"Sweep data for {problem_name}-{elim_name} (indices are base-0):")
            print(sweep_data_df)

            suffix = ""
            if args.sample is not None:
                n_instances_total = args.nsamples ** len(problem.parameters)
                # Note that, to collect these files, I will need to know how many
                # samples I expect.
                suffix += f"-{args.sample}of{n_instances_total}"
            if args.suffix is not None:
                suffix += f"-{args.suffix}"
            # TODO: Update this naming convention for consistency with structure
            # and solvetime. This will require an update to the plotting script.
            fname = f"{problem_name}-{elim_name}-sweep{suffix}.csv"
            fpath = os.path.join(args.results_dir, fname)
            if not args.no_save:
                print(f"Writing parameter sweep results to {fpath}")
                sweep_data_df.to_csv(fpath)

    # TODO: Write a sweep summary dataframe.
    # Should we only do this if args.sample is None?
    for problem_name, problem in problems:
        for elim_name, _ in elimination_callbacks:
            samples, results = sweep_results_lookup[problem_name, elim_name]
            n_instances = len(samples)
            n_converged = n_converged_lookup[problem_name, elim_name]
            print(f"{problem_name}-{elim_name} converged {n_converged} / {n_instances} instances")


if __name__ == "__main__":
    argparser = config.get_sweep_argparser()

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
