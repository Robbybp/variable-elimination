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
    get_trivial_constraint_elimination,
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
    ["upper_bound", "constraints", "var_expressions"],
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

    timer.start("maximum_matching")
    matching = linear_igraph.maximum_matching()
    timer.stop("maximum_matching")
    ub = len(matching)

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

    results = ElimResults(ub, con_elim, var_exprs)
    return results

def d1_elim_callback(model, **kwds):
    igraph = kwds.pop("igraph", None)
    linear_igraph = kwds.pop("linear_igraph", None)

    if igraph is None:
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
    if linear_igraph is None:
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
    eq_cons = get_equality_constraints(model)

    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)

    total_con_elim = []
    total_var_exprs = []
    while True:
        var_elim, con_elim = get_degree_one_elimination(model, linear_igraph = linear_igraph, eq_igraph = eq_igraph)
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            # We need to perform the elimination on the linear and equality-only
            # incidence graphs so we can choose new degree-1 constraints that
            # were potentially revealed.
            #
            # We don't have the potential bug where we add a nonlinear edge
            # to the linear incidence graph because eliminating a degree-1
            # constraint never adds edges. This allows us to reuse the linear
            # incidence graph and be slightly more efficient here.
            eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue
        break
    return ElimResults(None, total_con_elim, total_var_exprs)

def d2_elim_callback(model, **kwds):
    igraph = kwds.pop("igraph", None)

    linear_igraph = kwds.pop("linear_igraph", None)

    if igraph is None:
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
    if linear_igraph is None:
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
    eq_cons = get_equality_constraints(model)
    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)

    total_con_elim = []
    total_var_exprs = []
    while True:
        var_elim, con_elim = get_degree_one_elimination(
            model,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        var_elim, con_elim = get_degree_two_elimination(
            model,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        # No d1 cons and no d2 cons
        break

    return ElimResults(None, total_con_elim, total_var_exprs)

def trivial_elim_callback(model, **kwds):
    igraph = kwds.pop("igraph", None)
    linear_igraph = kwds.pop("linear_igraph", None)

    if igraph is None:
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
    if linear_igraph is None:
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
    eq_cons = get_equality_constraints(model)
    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)

    total_con_elim = []
    total_var_exprs = []
    while True:
        var_elim, con_elim = get_degree_one_elimination(
            model,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        var_elim, con_elim = get_trivial_constraint_elimination(
            model,
            allow_affine=True,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        # No d1 cons and no d2 cons
        break

    return ElimResults(None, total_con_elim, total_var_exprs)

def linear_d2_elim_callback(model, **kwds):
    igraph = kwds.pop("igraph", None)
    linear_igraph = kwds.pop("linear_igraph", None)
    if igraph is None:
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )
    if linear_igraph is None:
        linear_igraph = IncidenceGraphInterface(
            model, linear_only=True, include_inequality=False
        )
    eq_cons = get_equality_constraints(model)
    eq_igraph = igraph.subgraph(igraph.variables, eq_cons)

    total_con_elim = []
    total_var_exprs = []
    while True:
        var_elim, con_elim = get_degree_one_elimination(
            model,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 1")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        var_elim, con_elim = get_linear_degree_two_elimination(
            model,
            allow_affine=True,
            #linear_igraph=linear_igraph,
            eq_igraph=eq_igraph,
        )
        if var_elim:
            var_elim, con_elim = define_elimination_order(var_elim, con_elim)
            var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph=igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
            )
            #eliminate_nodes_from_graph(linear_igraph, var_elim, con_elim)
            eliminate_nodes_from_graph(eq_igraph, var_elim, con_elim)
            print(f"Eliminated {len(var_elim)} constraints of degree 2")
            total_con_elim.extend(con_elim)
            total_var_exprs.extend(var_exprs)
            continue

        # no d1 cons or d2 cons
        break

    return ElimResults(None, total_con_elim, total_var_exprs)


def no_elim_callback(model, **kwds):
    return ElimResults(None, [], [])

def ampl_elim_callback(model, **kwds):
    igraph = kwds.pop('igraph', None)
    if igraph is None:
        igraph = IncidenceGraphInterface(
            model,
            linear_only=False,
            include_inequality=True,
            method=IncidenceMethod.ampl_repn,
        )

    randomize = kwds.pop('randomize', False)
    eliminate_bounded_vars =kwds.pop('eliminate_bounded_vars', False)
    eliminate_linear_cons_only =kwds.pop('eliminate_lineaR_cons_only', False)
    
    var_elim, con_elim = identify_vars_for_elim_ampl(model, randomize, eliminate_bounded_vars, eliminate_linear_cons_only)
    if var_elim:
        var_elim, con_elim = define_elimination_order(var_elim, con_elim)
        var_exprs, _, _ = eliminate_variables(
                model,
                var_elim,
                con_elim,
                igraph = igraph,
                use_named_expressions=USE_NAMED_EXPRESSIONS,
                )
    else:
        con_elim = []
        var_exprs = []
    return ElimResults(None, con_elim, var_exprs)
