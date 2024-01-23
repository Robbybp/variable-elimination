import pyomo.environ as pyo
from pyomo.common.collections import ComponentMap
from pyomo.common.timing import TicTocTimer
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.config import IncidenceMethod
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix

import itertools
from collections import namedtuple

import scipy.sparse as sps
import matplotlib.pyplot as plt

from var_elim.models.distillation.distill import create_instance as create_distill
#from var_elim.models.opf.opf_model import make_model as create_opf
#from var_elim.models.gas_pipelines.gas_network_model import make_dynamic_model as create_pipeline
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


IncStructure = namedtuple(
    "IncStructure",
    ["nvar", "ncon", "nnz", "nnz_linear", "nnode", "n_nl_node", "n_linear_node"],
)
ElimResults = namedtuple("ElimResults", ["upper_bound"])
StructuralResults = namedtuple(
    "StructuralResults",
    ["orig", "reduced", "elim"],
)


def get_structural_results(model, elim_callback):
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
    # This is too matching-specific. This should be handled by the elimination callback
    #matching = orig_linear_igraph.maximum_matching()
    #n_linear_match = len(matching)

    orig_nnode = count_model_nodes(model)
    timer.toc("Count Pyomo nodes")
    orig_nlnode = count_model_nodes(model, amplrepn=True)
    timer.toc("Count nl nodes")
    orig_linear_nlnode = count_model_nodes(model, amplrepn=True, linear_only=True)
    timer.toc("Count linear nl nodes")

    elim_res = elim_callback(model)
    timer.toc("Perform elimination")

    reduced_nnode = count_model_nodes(model)
    timer.toc("Count reduced Pyomo nodes")
    reduced_nlnode = count_model_nodes(model, amplrepn=True)
    timer.toc("Count reduced nl nodes")
    reduced_linear_nlnode = count_model_nodes(model, amplrepn=True, linear_only=True)
    timer.toc("Count reduced linear nl nodes")

    reduced_igraph = IncidenceGraphInterface(
        model, include_inequality=True, method=IncidenceMethod.ampl_repn
    )
    reduced_linear_igraph = IncidenceGraphInterface(
        model,
        include_inequality=True,
        linear_only=True,
        method=IncidenceMethod.ampl_repn,
    )
    timer.toc("Reduced incidence graph")

    orig_nvar = len(orig_igraph.variables)
    orig_ncon = len(orig_igraph.constraints)
    orig_nnz = orig_igraph.n_edges
    orig_nnz_linear = orig_linear_igraph.n_edges

    reduced_nvar = len(reduced_igraph.variables)
    reduced_ncon = len(reduced_igraph.constraints)
    reduced_nnz = reduced_igraph.n_edges
    reduced_nnz_linear = reduced_linear_igraph.n_edges

    orig_struc = IncStructure(
        orig_nvar,
        orig_ncon,
        orig_nnz,
        orig_nnz_linear,
        orig_nnode,
        orig_nlnode,
        orig_linear_nlnode,
    )
    reduced_struc = IncStructure(
        reduced_nvar,
        reduced_ncon,
        reduced_nnz,
        reduced_nnz_linear,
        reduced_nnode,
        reduced_nlnode,
        reduced_linear_nlnode,
    )

    results = StructuralResults(orig_struc, reduced_struc, elim_res)
    return results


def matching_elim_callback(model):
    linear_igraph = IncidenceGraphInterface(
        model,
        include_inequality=False,
        linear_only=True,
        method=IncidenceMethod.ampl_repn,
    )
    matching = linear_igraph.maximum_matching()
    ub = len(matching)
    var_elim, con_elim = generate_elimination_via_matching(model)
    var_elim, con_elim = define_elimination_order(var_elim, con_elim)
    eliminate_variables(model, var_elim, con_elim)
    results = ElimResults(ub)
    return results


def d1_elim_callback(model):
    while True:
        var_elim, con_elim = get_degree_one_elimination(model)
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            eliminate_variables(model, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            continue
        break
    return ElimResults(None)


def d2_elim_callback(model):
    while True:
        var_elim, con_elim = get_degree_one_elimination(model)
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            eliminate_variables(model, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            continue

        var_elim, con_elim = get_degree_two_elimination(model)
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            eliminate_variables(model, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            continue

        # No d1 cons and no d2 cons
        break

    return ElimResults(None)


def trivial_elim_callback(model):
    while True:
        var_elim, con_elim = get_degree_one_elimination(model)
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            eliminate_variables(model, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            continue

        var_elim, con_elim = get_trivial_constraint_elimination(model, allow_affine=True)
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            eliminate_variables(model, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            continue

        # No d1 cons and no d2 cons
        break

    return ElimResults(None)


def linear_d2_elim_callback(model):
    while True:
        var_elim, con_elim = get_degree_one_elimination(model)
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            eliminate_variables(model, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            continue

        var_elim, con_elim = get_linear_degree_two_elimination(
            model, allow_affine=True
        )
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            eliminate_variables(model, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            continue

        # no d1 cons or d2 cons
        break

    return ElimResults(None)


def solve_original(m, tee=True):
    solver = pyo.SolverFactory("ipopt")
    #solver.options["print_timing_statistics"] = "yes"
    solver.solve(m, tee=tee)
    return m


def solve_reduced(m, tee=True):
    var_elim, con_elim = generate_elimination_via_matching(m)
    var_elim, con_elim = define_elimination_order(var_elim, con_elim)
    eliminate_variables(m, var_elim, con_elim)
    solver = pyo.SolverFactory("ipopt")
    #solver.options["print_timing_statistics"] = "yes"
    solver.solve(m, tee=tee)
    return m


def main():
    horizon = 300
    nfe = 300
    models = [
        ("Distill", lambda : create_distill(horizon=horizon, nfe=nfe)),
        #("OPF-4917", create_opf()),
        #("Pipeline", create_pipeline()),
    ]

    elim_callbacks = [
        ("Degree=1", d1_elim_callback),
        ("Trivial", trivial_elim_callback),
        ("Linear, degree=2", linear_d2_elim_callback),
        ("Degree=2", d2_elim_callback),
        ("Matching", matching_elim_callback),
    ]
    model_cb_elim_prod = list(itertools.product(models, elim_callbacks))
    model_elim_prod = []
    for (mname, model_cb), (ename, elim_cb) in model_cb_elim_prod:
        model = model_cb()
        model_elim_prod.append(((mname, model), (ename, elim_cb)))

    #m1 = create_distill(horizon=horizon, nfe=nfe)
    #m2 = create_distill(horizon=horizon, nfe=nfe)
    #solve_original(m1, tee=True)
    #solve_reduced(m2, tee=True)

    for i in range(len(model_cb_elim_prod)):
        mname, model = model_elim_prod[i][0]
        elim_name, elim_callback = model_elim_prod[i][1]
        nchar = len(mname) + len(elim_name) + 5
        print()
        print(f"{mname} -- {elim_name}")
        print("-"*nchar)
        #print_statistics(model)
        results = get_structural_results(model, elim_callback)

        orig_nnz_per_con = results.orig.nnz / results.orig.ncon
        reduced_nnz_per_con = results.reduced.nnz / results.reduced.ncon
        n_elim = results.orig.nvar - results.reduced.nvar
        ncon_diff = results.orig.ncon - results.reduced.ncon
        if results.elim.upper_bound is not None:
            elimination_gap = (results.elim.upper_bound - n_elim) / results.elim.upper_bound
        else:
            elimination_gap = None

        #pyo.SolverFactory("ipopt").solve(model, tee=True)
        print(results)
        print()
        print(f"N. var eliminated: {n_elim}")
        # Note the ncon here includes inequalities, which we could be adding
        print(f"N. con eliminated: {ncon_diff}")

        print(f"Original total NNZ: {results.orig.nnz}")
        print(f"Reduced total NNZ: {results.reduced.nnz}")
        print(f"Original linear NNZ: {results.orig.nnz_linear}")
        print(f"Reduced linear NNZ: {results.reduced.nnz_linear}")

        print(f"Original NNZ/con: {orig_nnz_per_con}")
        print(f"Reduced NNZ/con: {reduced_nnz_per_con}")

        print(f"Size of linear matching: {results.elim.upper_bound}")
        print(f"Gap to maximum elimination: {elimination_gap}")
        print(f"Original n. nodes: {results.orig.nnode}")
        print(f"Reduced n. nodes: {results.reduced.nnode}")
        print(f"Original n. nl nodes (total): {results.orig.n_nl_node}")
        print(f"Reduced n. nl nodes (total): {results.reduced.n_nl_node}")
        # Note that these are number of linear nodes in the nl file. We can't
        # (easily) extract the linear nodes directly from the Pyomo expression tree.
        print(f"Original n. linear nodes: {results.orig.n_linear_node}")
        print(f"Reduced n. linear nodes: {results.reduced.n_linear_node}")
        orig_nonlin_nodes = results.orig.n_nl_node - results.orig.n_linear_node
        reduced_nonlin_nodes = results.reduced.n_nl_node - results.reduced.n_linear_node
        print(f"Original n. nonlinear nodes: {orig_nonlin_nodes}")
        print(f"Reduced n. nonlinear nodes: {reduced_nonlin_nodes}")


if __name__ == "__main__":
    main()
