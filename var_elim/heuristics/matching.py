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

from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from var_elim.algorithms.replace import define_elimination_order

def generate_elimination_via_matching(m):
    # TODO: Optional IncidenceGraphInterface argument
    linear_igraph = IncidenceGraphInterface(
        m, active=True, include_fixed=False, linear_only=True
    )
    matching = linear_igraph.maximum_matching()
    import pdb; pdb.set_trace()
    
    con_list = list(matching.keys())
    var_list = list(matching.values())

    var_order, con_order = define_elimination_order(var_list, con_list)

    # TODO: If we are not lower triangular in the full incidence graph,
    # generate a "reduced elimination order" that is lower triangular.
    # Our assumption is that generated orderings will be "almost lower
    # triangular".

    return var_order, con_order


if __name__ == "__main__":
    from var_elim.models.distillation.distill import create_instance
    m = create_instance()
    var_order, con_order = generate_elimination_via_matching(m)
