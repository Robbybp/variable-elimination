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

import enum
from pyomo.util.subsystems import create_subsystem_block, TemporarySubsystemManager
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.config import IncidenceMethod
from var_elim.algorithms.replace import define_elimination_order


from pyomo.common.timing import HierarchicalTimer


class TearMethod(enum.Enum):
    greedy = 0
    # TODO: more sophisticated heuristic, e.g. Elmqvist & Otter


def break_algebraic_loop_greedy(igraph, matching):
    # We assume that igraph does not decompose
    reduced_system = []
    # Need a nonzero diagonal of linear incidence
    used_constraints = set()
    var_set = set(id(var) for var in igraph.variables)
    for con in igraph.constraints:
        # We assume that this (var, con) edge is linear
        var = matching[con]

        # My math says that con's matched variable, in any perfect matching, must
        # be part of this diagonal block. Assert this in case I am wrong.
        assert id(var) in var_set

        if not any(
            # If the variable participates in some constraint that we've seen
            # before, don't use it. This guarantees that incidence matrix is
            # lower triangularizable.
            id(adj_con) in used_constraints
            for adj_con in igraph.get_adjacent_to(var)
        ):
            used_constraints.add(id(con))
            reduced_system.append((var, con))
    # Since each variable was specified by the matching, it is linear in
    # its associated constraint.
    variables = [var for var, con in reduced_system]
    constraints = [con for var, con in reduced_system]
    return variables, constraints


_dispatcher = {
    TearMethod.greedy: break_algebraic_loop_greedy,
}


def break_algebraic_loop(igraph, matching, method=TearMethod.greedy):
    # TODO: break_algebraic_loops function that allows decomposable systems
    return _dispatcher[method](igraph, matching)


def generate_elimination_via_matching(
    m,
    linear_igraph=None,
    igraph=None,
    timer=None,
):
    """
    Parameters
    ----------
    m : model
    linear_igraph : Incidence graph with edges corresponding to linear
        variable-constraint incidence. Should *not* include inequalities.
    igraph : Full incidence graph. Should *not* include inequalities.

    """
    if timer is None:
        timer = HierarchicalTimer()
    timer.start("linear_igraph")
    if linear_igraph is None:
        linear_igraph = IncidenceGraphInterface(
            m,
            active=True,
            include_fixed=False,
            include_inequality=False,
            linear_only=True,
        )
    timer.stop("linear_igraph")
    timer.start("maximum_matching")
    matching = linear_igraph.maximum_matching()
    timer.stop("maximum_matching")

    con_list = list(matching.keys())
    var_list = list(matching.values())

    timer.start("igraph")
    if igraph is None:
        # NOTE: Use ampl_repn here as we don't want spurious nonzeros to contribute
        # to algebraic loops.
        igraph = IncidenceGraphInterface(
            m,
            active=True,
            include_fixed=False,
            include_inequality=False,
            linear_only=False,
        )
    timer.stop("igraph")
    timer.start("block_triang")
    # Note that this already uses an efficient subgraph, so igraph does not
    # need to be restricted to any subset of variables/constraints.
    matched_subgraph = igraph.subgraph(var_list, con_list)
    var_blocks, con_blocks = matched_subgraph.block_triangularize()
    timer.stop("block_triang")

    var_order = []
    con_order = []
    for vb, cb in zip(var_blocks, con_blocks):
        if len(vb) == 1:
            # If v, c is in a 1x1 block in this matching, we know that it is
            # linear:
            # - (v, c) is in any perfect matching on this graph
            # Claim: zip(var_list, con_list) is a perfect matching on this graph.
            # - var_list, con_list is a perfect matching on the linear subgraph
            #   on var_list, con_list (almost by definition - it contains every
            #   node in var_list, con_list)
            # - Therefore var_list, con_list is a perfect matching on the full
            #   (linear & nonlinear) subgraph on var_list, con_list. (All we
            #   did was add edges, we have the same nodes)
            # Therefore, (v, c) is in var_list, con_list, meaning that (v, c)
            # is linear.
            var_order.extend(vb)
            con_order.extend(cb)
        else:
            # We need to ensure that the reduced subsystem is linear.
            # Is this guaranteed by the fact that the matching is linear?
            # Yes. We only return variables and constraints that are paired
            # according to this matching.
            #
            # We assume that, for any con in this block, matching[con] is a
            # variable in this block. I believe this is true (TODO: prove
            # -- this should follow from uniqueness of strongly connected
            # components).
            timer.start("break_loop")
            timer.start("subgraph")
            subgraph = matched_subgraph.subgraph(vb, cb)
            timer.stop("subgraph")
            reduced_vb, reduced_cb = break_algebraic_loop(subgraph, matching)
            timer.stop("break_loop")
            var_order.extend(reduced_vb)
            con_order.extend(reduced_cb)

    return var_order, con_order


if __name__ == "__main__":
    from var_elim.models.distillation.distill import create_instance

    m = create_instance()
    var_order, con_order = generate_elimination_via_matching(m)
    var_order, con_order = define_elimination_order(var_order, con_order)
