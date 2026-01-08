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
from pyomo.common.timing import TicTocTimer, HierarchicalTimer
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.config import IncidenceMethod

from collections import namedtuple

from var_elim.heuristics.matching import (
    generate_elimination_via_matching,
    define_elimination_order,
)
from var_elim.heuristics.trivial_elimination import (
    get_degree_one_elimination,
    get_degree_two_elimination,
    get_linear_degree_two_elimination,
    get_equalcoef_degree_two_elimination,
    filter_constraints,
)
from var_elim.heuristics.ampl_heuristic import identify_vars_for_elim_ampl
from var_elim.algorithms.replace import eliminate_variables, eliminate_nodes_from_graph


"""elimination_callbacks.py

Functions for performing different types of variable eliminations, with
consistent call signatures. These are implemented in one place for consistency
when running across multiple models and experiments.

"""


# TODO: Potentially move these to their own "data.py" file
IncStructure = namedtuple(
    "IncStructure",
    ["nvar", "ncon", "nnz", "nnz_linear", "nnode", "n_nl_node", "n_linear_node"],
)
ElimResults = namedtuple(
    "ElimResults",
    ["upper_bound", "lower_bound", "constraints", "var_expressions", "max_block_size"],
)
StructuralResults = namedtuple(
    "StructuralResults",
    ["orig", "reduced", "elim"],
)


# TODO: Command line argument
USE_NAMED_EXPRESSIONS = True


def get_equality_constraints(model):
    eq_cons = []
    for con in model.component_data_objects(pyo.Constraint, active=True):
        if isinstance(con.expr, EqualityExpression):
            eq_cons.append(con)
    return eq_cons


def matching_elim_callback(model, **kwds):
    timer = kwds.pop("timer", HierarchicalTimer())
    # In case we don't want to modify the model in-place
    eliminate = kwds.pop("eliminate", True)

    igraph = kwds.pop("igraph", None)
    linear_igraph = kwds.pop("linear_igraph", None)

    # Flag we can set to avoid the unnecessary maximum matching and block triangularize
    # that we use to get structural results.
    skip_extra_info = kwds.pop("skip_extra_info", False)

    # Construct all the incidence graphs we will need for this analysis
    #
    # The full incidence graph is used to determine variable/constraint incidence.
    # It does not necessarily need to use ampl_repn, although "unnecessary
    # replacements" may lead to spurious variables in the NL file (that only
    # participate with coefficients of zero). (Aside: I'm not sure this is
    # true. If AMPLRepn is always built root-to-leaf, then I think we should
    # be fine. Where exactly was the bug that led to me using ampl_repn?)
    if igraph is None:
        timer.start("igraph")
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
        timer.stop("igraph")
    eq_cons = get_equality_constraints(model)

    # We need an incidence graph on only equality constraints to enforce lower
    # triangularity of the eliminated variables and constraints.
    timer.start("subgraph")
    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)
    timer.stop("subgraph")

    # We need a linear, equality-only incidence graph to actually select
    # variable/constraint elimination pairs. This *could potentially* be
    # constructed as an edge-subgraph, if we stored a linear/nonlinear
    # indicator for each edge.
    timer.start("linear_igraph")
    if linear_igraph is None:
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
    timer.stop("linear_igraph")

    if skip_extra_info:
        ub = None
        lb = None
        max_block_size = None
    else:
        timer.start("maximum_matching")
        matching = linear_igraph.maximum_matching()
        timer.stop("maximum_matching")
        ub = len(matching)
        timer.start("subgraph")
        matching_subgraph = eq_igraph.subgraph(list(matching.values()), list(matching.keys()))
        timer.stop("subgraph")
        timer.start("block_triangularize")
        vblocks, cblocks = matching_subgraph.block_triangularize()
        timer.stop("block_triangularize")
        lb = len(vblocks)
        max_block_size = max(len(b) for b in vblocks)

    timer.start("generate_elimination")
    var_elim, con_elim = generate_elimination_via_matching(
        model,
        linear_igraph=linear_igraph,
        igraph=eq_igraph,
        timer=timer,
    )
    timer.stop("generate_elimination")

    #timer.start("define_order")
    # Not necessary for this algorithm
    #elim_subgraph = eq_igraph.subgraph(con_elim, var_elim)
    #var_elim, con_elim = define_elimination_order(
    #    var_elim, con_elim, igraph=igraph,
    #)
    #timer.stop("define_order")

    if eliminate:
        var_exprs, var_lb_map, var_ub_map = eliminate_variables(
            model,
            var_elim,
            con_elim,
            igraph=igraph,
            use_named_expressions=USE_NAMED_EXPRESSIONS,
            timer=timer,
        )
    else:
        # TODO: Find a way to communicate this information other than
        # overloading the var_expressions field
        var_exprs = var_elim

    results = ElimResults(ub, lb, con_elim, var_exprs, max_block_size)
    return results

def linear_degree1_elim_callback(model, **kwds):
    igraph = kwds.pop("igraph", None)
    linear_igraph = kwds.pop("linear_igraph", None)
    timer = kwds.pop("timer", HierarchicalTimer())

    if igraph is None:
        timer.start("igraph")
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
        timer.stop("igraph")
    if linear_igraph is None:
        timer.start("linear-igraph")
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
        timer.stop("linear-igraph")
    eq_cons = get_equality_constraints(model)

    timer.start("subgraph")
    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)
    timer.stop("subgraph")

    total_con_elim = []
    total_var_exprs = []
    while True:
        timer.start("get-d1")
        var_elim, con_elim = get_degree_one_elimination(
            model,
            linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        timer.stop("get-d1")
        if var_elim:
            #timer.start("order")
            #var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            #timer.stop("order")
            timer.start("eliminate")
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            timer.stop("eliminate")
            # We need to perform the elimination on the linear and equality-only
            # incidence graphs so we can choose new degree-1 constraints that
            # were potentially revealed.
            #
            # We don't have the potential bug where we add a nonlinear edge
            # to the linear incidence graph because eliminating a degree-1
            # constraint never adds edges. This allows us to reuse the linear
            # incidence graph and be slightly more efficient here.
            timer.start("eliminate-nodes")
            eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            timer.stop("eliminate-nodes")
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue
        break
    return ElimResults(None, None, total_con_elim, total_var_exprs, None)

def degree2_elim_callback(model, **kwds):
    igraph = kwds.pop("igraph", None)
    linear_igraph = kwds.pop("linear_igraph", None)
    timer = kwds.pop("timer", HierarchicalTimer())

    if igraph is None:
        timer.start("igraph")
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
        timer.stop("igraph")
    if linear_igraph is None:
        timer.start("linear-igraph")
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
        timer.stop("linear-igraph")
    eq_cons = get_equality_constraints(model)
    timer.start("subgraph")
    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)
    timer.stop("subgraph")

    total_con_elim = []
    total_var_exprs = []
    while True:
        timer.start("get-d1")
        var_elim, con_elim = get_degree_one_elimination(
            model,
            # I think we can't reuse this as we don't have a good way of modifying
            # the linear-only incidence graph.
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        timer.stop("get-d1")
        if var_elim:
            #timer.start("order")
            #var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            #timer.stop("order")
            timer.start("eliminate")
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            timer.stop("eliminate")
            timer.start("eliminate-nodes")
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            timer.stop("eliminate-nodes")
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        timer.start("get-d2")
        var_elim, con_elim = get_degree_two_elimination(
            model,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        timer.stop("get-d2")
        if var_elim:
            #timer.start("order")
            #var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            #timer.stop("order")
            timer.start("eliminate")
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            timer.stop("eliminate")
            timer.start("eliminate-nodes")
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            timer.stop("eliminate-nodes")
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        # No d1 cons and no d2 cons
        break

    return ElimResults(None, None, total_con_elim, total_var_exprs, None)

def equalcoef_degree1_elim_callback(model, **kwds):
    igraph = kwds.pop("igraph", None)
    linear_igraph = kwds.pop("linear_igraph", None)
    timer = kwds.pop("timer", HierarchicalTimer())

    if igraph is None:
        timer.start("igraph")
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
        timer.stop("igraph")
    if linear_igraph is None:
        timer.start("linear-igraph")
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
        timer.stop("linear-igraph")
    eq_cons = get_equality_constraints(model)
    timer.start("subgraph")
    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)
    timer.stop("subgraph")

    total_con_elim = []
    total_var_exprs = []
    while True:
        timer.start("get-d1")
        var_elim, con_elim = get_degree_one_elimination(
            model,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        timer.stop("get-d1")
        if var_elim:
            #timer.start("order")
            #var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            #timer.stop("order")
            timer.start("eliminate")
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            timer.stop("eliminate")
            timer.start("eliminate-nodes")
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            timer.stop("eliminate-nodes")
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        timer.start("get-d2")
        var_elim, con_elim = get_equalcoef_degree_two_elimination(
            model,
            allow_affine=True,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        timer.stop("get-d2")
        if var_elim:
            #timer.start("order")
            #var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            #timer.stop("order")
            timer.start("eliminate")
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            timer.stop("eliminate")
            timer.start("eliminate-nodes")
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            timer.stop("eliminate-nodes")
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        # No d1 cons and no d2 cons
        break

    return ElimResults(None, None, total_con_elim, total_var_exprs, None)

def linear_degree2_elim_callback(model, **kwds):
    igraph = kwds.pop("igraph", None)
    linear_igraph = kwds.pop("linear_igraph", None)
    timer = kwds.pop("timer", HierarchicalTimer())
    if igraph is None:
        timer.start("igraph")
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
        timer.stop("igraph")
    if linear_igraph is None:
        timer.start("linear-igraph")
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
        timer.stop("linear-igraph")
    eq_cons = get_equality_constraints(model)
    timer.start("subgraph")
    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)
    timer.stop("subgraph")

    total_con_elim = []
    total_var_exprs = []
    while True:
        timer.start("get-d1")
        var_elim, con_elim = get_degree_one_elimination(
            model,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        timer.stop("get-d1")
        if var_elim:
            #timer.start("order")
            #var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            #timer.stop("order")
            timer.start("eliminate")
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            timer.stop("eliminate")
            timer.start("eliminate-nodes")
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            timer.stop("eliminate-nodes")
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        timer.start("get-d2")
        var_elim, con_elim = get_linear_degree_two_elimination(
            model,
            allow_affine=True,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        timer.stop("get-d2")
        if var_elim:
            #timer.start("order")
            #var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            #timer.stop("order")
            timer.start("eliminate")
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            timer.stop("eliminate")
            timer.start("eliminate-nodes")
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            timer.stop("eliminate-nodes")
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        # no d1 cons or d2 cons
        break

    return ElimResults(None, None, total_con_elim, total_var_exprs, None)


def no_elim_callback(model, **kwds):
    return ElimResults(None, None, [], [], None)

def greedy_elim_callback(model, **kwds):
    igraph = kwds.pop('igraph', None)
    timer = kwds.pop("timer", HierarchicalTimer())
    if igraph is None:
        timer.start("igraph")
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
        timer.stop("igraph")

    randomize = kwds.pop('randomize', False)
    eliminate_bounded_vars = kwds.pop('eliminate_bounded_vars', True)
    eliminate_linear_cons_only = kwds.pop('eliminate_linear_cons_only', False)
    first_variable_only = kwds.pop("first_variable_only", False)

    timer.start("iden-ampl")
    var_elim, con_elim = identify_vars_for_elim_ampl(
        model,
        randomize=randomize,
        eliminate_bounded_vars=eliminate_bounded_vars,
        eliminate_linear_cons_only=eliminate_linear_cons_only,
        first_variable_only=first_variable_only,
    )
    timer.stop("iden-ampl")
    if var_elim:
        #timer.start("order")
        #var_elim, con_elim = define_elimination_order(var_elim, con_elim)
        #timer.stop("order")
        timer.start("eliminate")
        var_exprs, _, _ = eliminate_variables(
            model,
            var_elim,
            con_elim,
            igraph = igraph,
            use_named_expressions=USE_NAMED_EXPRESSIONS,
        )
        timer.stop("eliminate")
    else:
        con_elim = []
        var_exprs = []
    return ElimResults(None, None, con_elim, var_exprs, None)
