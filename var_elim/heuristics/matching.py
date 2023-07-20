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
from pyomo.util.subsystems import (
    create_subsystem_block, TemporarySubsystemManager
)
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from var_elim.algorithms.replace import define_elimination_order


class TearMethod(enum.Enum):
    greedy = 0
    # TODO: more sophisticated heuristic, e.g. Elmqvist & Otter


def break_algebraic_loop_greedy(igraph, matching):
    # We assume that igraph does not decompose
    reduced_system = []
    # Need a nonzero diagonal of linear incidence
    used_constraints = set()
    for con in igraph.constraints:
        # We assume that this (var, con) edge is linear
        var = matching[con]
        if not any(
            # If the variable participates in some constraint that we've seen
            # before, don't use it. This guarantees that incidence matrix is
            # lower triangularizable.
            id(adj_con) in used_constraints
            for adj_con in igraph.get_adjacent_to(var)
        ):
            used_constraints.add(id(con))
            reduced_system.append((var, con))
    variables = [var for var, con in reduced_system]
    constraints = [con for var, con in reduced_system]
    return variables, constraints


_dispatcher = {
    TearMethod.greedy: break_algebraic_loop_greedy,
}


def break_algebraic_loop(
    variables, constraints, matching, method=TearMethod.greedy
):
    # TODO: Optional IncidenceGraphInterface argument
    # TODO: break_algebraic_loops function that allows decomposable systems
    subsystem = create_subsystem_block(constraints, variables)
    to_fix = list(subsystem.input_vars[:])
    with TemporarySubsystemManager(to_fix=to_fix):
        igraph = IncidenceGraphInterface(subsystem)
    var_blocks, con_blocks = igraph.block_triangularize()
    if len(var_blocks) != 1:
        # The incidence matrix does not satisfy the strong Hall property.
        raise RuntimeError(
            f"break_algebraic_loop only accepts systems that do not decompose."
            f"Got {len(var_blocks)} strongly connected components."
        )
    return _dispatcher[method](igraph, matching)
    


def generate_elimination_via_matching(m):
    # TODO: Optional IncidenceGraphInterface argument
    linear_igraph = IncidenceGraphInterface(
        m, active=True, include_fixed=False, linear_only=True
    )
    matching = linear_igraph.maximum_matching()

    con_list = list(matching.keys())
    var_list = list(matching.values())

    igraph = IncidenceGraphInterface(m)
    var_blocks, con_blocks = igraph.block_triangularize(var_list, con_list)

    var_order = []
    con_order = []
    for vb, cb in zip(var_blocks, con_blocks):
        if len(vb) == 1:
            var_order.extend(vb)
            con_order.extend(cb)
        else:
            reduced_vb, reduced_cb = break_algebraic_loop(vb, cb, matching)
            var_order.extend(reduced_vb)
            con_order.extend(reduced_cb)

    return var_order, con_order


if __name__ == "__main__":
    from var_elim.models.distillation.distill import create_instance
    m = create_instance()
    var_order, con_order = generate_elimination_via_matching(m)
    var_order, con_order = define_elimination_order(var_order, con_order)
