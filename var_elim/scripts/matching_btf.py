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

import var_elim.scripts.config as config
from var_elim.heuristics.matching import break_algebraic_loop_greedy
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface


def main(args):
    model = config.CONSTRUCTOR_LOOKUP[args.model]()
    igraph = IncidenceGraphInterface(model, include_inequality=False)
    matching = igraph.maximum_matching()
    variables = list(matching.values())
    constraints = list(matching.keys())
    subg = igraph.subgraph(variables, constraints)
    vblocks, cblocks = subg.block_triangularize()

    lt_matching = []
    for vb, cb in zip(vblocks, cblocks):
        if len(vb) == 1:
            lt_matching.append((vb[0], cb[0]))
        else:
            sg = subg.subgraph(vb, cb)
            variables, constraints = break_algebraic_loop_greedy(sg, matching)
            lt_matching.extend(zip(variables, constraints))

    nvar = len(igraph.variables)
    ncon = len(igraph.constraints)
    nmatch = len(matching)
    nblocks = len(vblocks)
    max_block_size = max(len(b) for b in vblocks)
    lt_matching_size = len(lt_matching)

    print(f"N. variables                      : {nvar}")
    print(f"N. eq. constraints                : {ncon}")
    print(f"|M|                               : {nmatch}")
    print(f"N. diagonal blocks                : {nblocks}")
    print(f"Max block size                    : {max_block_size}")
    print(f"Size of lower triangular matching : {lt_matching_size}")


if __name__ == "__main__":
    argparser = config.get_argparser()
    args = argparser.parse_args()
    if args.model is None:
        raise ValueError("--model must be provided.")
    main(args)
