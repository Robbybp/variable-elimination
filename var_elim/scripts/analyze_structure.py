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
from pyomo.core.expr import EqualityExpression
from pyomo.common.collections import ComponentMap
from pyomo.common.timing import TicTocTimer, HierarchicalTimer
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.config import IncidenceMethod
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix

import itertools
from collections import namedtuple

import scipy.sparse as sps
import matplotlib.pyplot as plt

from var_elim.models.distillation.distill import create_instance as create_distill
#from var_elim.models.opf.opf_model import make_model as create_opf
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
    count_nodes,
    count_model_nodes,
    count_amplrepn_nodes,
    count_model_amplrepn_nodes,
)
from var_elim.elimination_callbacks import (
    get_equality_constraints,
    matching_elim_callback,
    trivial_elim_callback,
    d1_elim_callback,
    d2_elim_callback,
    linear_d2_elim_callback,
)

import os
import var_elim.scripts.config as config
import pandas as pd
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP

import pselib


IncStructure = namedtuple(
    "IncStructure",
    [
        "nvar",
        "ncon",
        "nnz",
        "nnz_linear",
        "nnz_hessian",
        "nnode",
        "n_nonlinear_node",
        "n_linear_node",
    ],
)
ElimResults = namedtuple("ElimResults", ["upper_bound"])
StructuralResults = namedtuple(
    "StructuralResults",
    ["orig", "reduced", "elim"],
)


# TODO: Command line argument
USE_NAMED_EXPRESSIONS = True


def get_nnz_in_hessian(model, nlp=None):
    if nlp is None:
        nlp = PyomoNLP(model)
    hessian = nlp.evaluate_hessian_lag()
    return hessian.nnz


def get_structural_results(model, elim_callback, htimer=None):
    htimer = HierarchicalTimer() if htimer is None else htimer
    timer = TicTocTimer()

    timer.tic()
    orig_igraph = IncidenceGraphInterface(
        model, include_inequality=True, method=IncidenceMethod.ampl_repn
    )
    timer.toc("Original incidence graph")
    orig_linear_igraph = IncidenceGraphInterface(
        model,
        include_inequality=True,
        linear_only=True,
        method=IncidenceMethod.ampl_repn,
    )
    timer.toc("Linear subgraph")
    eq_cons = get_equality_constraints(model)
    orig_linear_eq_igraph = orig_linear_igraph.subgraph(
        orig_linear_igraph.variables, eq_cons
    )
    timer.toc("Linear eq-only subgraph")

    orig_nnode = count_model_nodes(model)
    timer.toc("Count Pyomo nodes")

    orig_nodecount = count_model_amplrepn_nodes(model)
    #orig_nlnode = count_model_nodes(model, amplrepn=True)
    timer.toc("Count nl nodes")
    #orig_linear_nlnode = count_model_nodes(model, amplrepn=True, linear_only=True)
    #timer.toc("Count linear nl nodes")

    orig_nvar = len(orig_igraph.variables)
    orig_ncon = len(orig_igraph.constraints)
    orig_nnz = orig_igraph.n_edges
    orig_nnz_linear = orig_linear_igraph.n_edges

    orig_nnz_hessian = get_nnz_in_hessian(model)

    elim_res = elim_callback(
        model, igraph=orig_igraph, linear_igraph=orig_linear_eq_igraph, timer=htimer
    )
    timer.toc("Perform elimination")

    reduced_nnode = count_model_nodes(model)
    timer.toc("Count reduced Pyomo nodes")

    reduced_nodecount = count_model_amplrepn_nodes(model)
    #reduced_nlnode = count_model_nodes(model, amplrepn=True)
    timer.toc("Count reduced nl nodes")
    #reduced_linear_nlnode = count_model_nodes(model, amplrepn=True, linear_only=True)
    #timer.toc("Count reduced linear nl nodes")

    reduced_igraph = IncidenceGraphInterface(
        model, include_inequality=True, method=IncidenceMethod.ampl_repn
    )
    timer.toc("Reduced incidence graph")
    reduced_linear_igraph = IncidenceGraphInterface(
        model,
        include_inequality=True,
        linear_only=True,
        method=IncidenceMethod.ampl_repn,
    )
    timer.toc("Reduced linear subgraph")

    reduced_nvar = len(reduced_igraph.variables)
    reduced_ncon = len(reduced_igraph.constraints)
    reduced_nnz = reduced_igraph.n_edges
    reduced_nnz_linear = reduced_linear_igraph.n_edges

    reduced_nnz_hessian = get_nnz_in_hessian(model)

    orig_struc = IncStructure(
        orig_nvar,
        orig_ncon,
        orig_nnz,
        orig_nnz_linear,
        orig_nnz_hessian,
        orig_nnode,
        orig_nodecount.nonlinear,
        orig_nodecount.linear,
        #orig_nlnode,
        #orig_linear_nlnode,
    )
    reduced_struc = IncStructure(
        reduced_nvar,
        reduced_ncon,
        reduced_nnz,
        reduced_nnz_linear,
        reduced_nnz_hessian,
        reduced_nnode,
        reduced_nodecount.nonlinear,
        reduced_nodecount.linear,
        #reduced_nlnode,
        #reduced_linear_nlnode,
    )

    results = StructuralResults(orig_struc, reduced_struc, elim_res)
    return results


def solve_original(m, tee=True):
    solver = pyo.SolverFactory("ipopt")
    #solver.options["print_timing_statistics"] = "yes"
    solver.solve(m, tee=tee)
    return m


def solve_reduced(m, tee=True):
    timer = TicTocTimer()
    timer.tic()
    var_elim, con_elim = generate_elimination_via_matching(m)
    var_exprs, var_lbs, var_ubs = eliminate_variables(m, var_elim, con_elim, use_named_expressions=USE_NAMED_EXPRESSIONS)
    timer.toc("eliminate variables")
    solver = pyo.SolverFactory("ipopt")
    #solver.options["print_timing_statistics"] = "yes"
    solver.solve(m, tee=tee)
    timer.toc("solve")
    from var_elim.algorithms.validate import validate_solution
    validate_solution(m, var_exprs, con_elim, tolerance=1e-6)
    timer.toc("validate")
    return m


def main(args):
    if args.model is not None:
        models = [(args.model, config.CONSTRUCTOR_LOOKUP[args.model])]
    else:
        models = list(zip(config.MODEL_NAMES, config.MODEL_CONSTRUCTORS))

    if args.method is not None:
        elim_callbacks = [(args.method, config.ELIM_LOOKUP[args.method])]
    else:
        elim_callbacks = config.ELIM_CALLBACKS

    model_cb_elim_prod = list(itertools.product(models, elim_callbacks))
    model_elim_prod = []
    for (mname, model_cb), (ename, elim_cb) in model_cb_elim_prod:
        model = model_cb()
        model_elim_prod.append(((mname, model), (ename, elim_cb)))

    data = {
        "model": [],
        "method": [],
        "nvar": [],
        "ncon": [],
        "n-elim": [],
        "n-elim-ub": [],
        "n-elim-lb": [],
        "max-block-size": [],
        "nnz": [],
        "nnz-linear": [],
        "nnz-hessian": [],
        "nnode-pyomo": [],
        "nnode-nl-linear": [],
        "nnode-nl-nonlinear": [],
    }

    timer = HierarchicalTimer()
    timer.start("root")

    for i in range(len(model_cb_elim_prod)):
        mname, model = model_elim_prod[i][0]
        elim_name, elim_callback = model_elim_prod[i][1]
        nchar = len(mname) + len(elim_name) + 5
        timer.start(f"{mname}-{elim_name}")
        print()
        print(f"{mname} -- {elim_name}")
        print("-"*nchar)
        results = get_structural_results(model, elim_callback, htimer=timer)

        orig_nnz_per_con = results.orig.nnz / results.orig.ncon
        reduced_nnz_per_con = results.reduced.nnz / results.reduced.ncon
        n_elim = results.orig.nvar - results.reduced.nvar
        ncon_diff = results.orig.ncon - results.reduced.ncon
        if results.elim.upper_bound is not None:
            elimination_gap = (results.elim.upper_bound - n_elim) / results.elim.upper_bound
        else:
            elimination_gap = None

        #pyo.SolverFactory("ipopt").solve(model, tee=True)
        # Not printing results for now, as the data structure is quite verbose.
        #print(results)
        print()
        print(f"N. var eliminated: {n_elim}")
        # Note the ncon here includes inequalities, which we could be adding
        print(f"N. con eliminated: {ncon_diff}")

        print(f"Original total NNZ: {results.orig.nnz}")
        print(f"Reduced total NNZ: {results.reduced.nnz}")
        print(f"Original linear NNZ: {results.orig.nnz_linear}")
        print(f"Reduced linear NNZ: {results.reduced.nnz_linear}")
        print(f"Original Hessian NNZ: {results.orig.nnz_hessian}")
        print(f"Reduced Hessian NNZ: {results.reduced.nnz_hessian}")

        print(f"Original NNZ/con: {orig_nnz_per_con}")
        print(f"Reduced NNZ/con: {reduced_nnz_per_con}")

        print(f"Size of linear matching: {results.elim.upper_bound}")
        print(f"Gap to maximum elimination: {elimination_gap}")
        print(f"Original n. nodes: {results.orig.nnode}")
        print(f"Reduced n. nodes: {results.reduced.nnode}")
        print(f"Original n. nl nodes (total): {results.orig.n_nonlinear_node}")
        print(f"Reduced n. nl nodes (total): {results.reduced.n_nonlinear_node}")
        # Note that these are number of linear nodes in the nl file. We can't
        # (easily) extract the linear nodes directly from the Pyomo expression tree.
        print(f"Original n. linear nodes: {results.orig.n_linear_node}")
        print(f"Reduced n. linear nodes: {results.reduced.n_linear_node}")
        orig_nonlin_nodes = results.orig.n_nonlinear_node
        reduced_nonlin_nodes = results.reduced.n_nonlinear_node
        print(f"Original n. nonlinear nodes: {orig_nonlin_nodes}")
        print(f"Reduced n. nonlinear nodes: {reduced_nonlin_nodes}")

        data["model"].append(mname)
        data["method"].append(elim_name)
        data["nvar"].append(results.reduced.nvar)
        data["ncon"].append(results.reduced.ncon)
        data["n-elim"].append(n_elim)
        data["n-elim-ub"].append(results.elim.upper_bound)
        data["n-elim-lb"].append(results.elim.lower_bound)
        data["max-block-size"].append(results.elim.max_block_size)
        data["nnz"].append(results.reduced.nnz)
        data["nnz-linear"].append(results.reduced.nnz_linear)
        data["nnz-hessian"].append(results.reduced.nnz_hessian)
        data["nnode-pyomo"].append(results.reduced.nnode)
        data["nnode-nl-linear"].append(results.reduced.n_linear_node)
        data["nnode-nl-nonlinear"].append(reduced_nonlin_nodes)

        timer.stop(f"{mname}-{elim_name}")

    df = pd.DataFrame(data)
    suffixes = [args.model, args.method, args.suffix]
    fname = config.get_basename("structure.csv", *suffixes)
    fpath = os.path.join(args.results_dir, fname)
    if not args.no_save:
        print(f"Writing results to {fpath}")
        df.to_csv(fpath)
    print(df)

    timer.stop("root")
    print(timer)


if __name__ == "__main__":
    argparser = config.get_argparser()
    argparser.add_argument(
        "--fname", default=None, help="Basename for output file (optional)"
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
            resdir = os.path.join(config.get_results_dir(), "structure")
            config.validate_dir(resdir)
            args.results_dir = resdir

    main(args)
