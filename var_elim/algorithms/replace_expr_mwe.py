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
from pyomo.core.expr.visitor import replace_expressions


def main():
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3])
    expr = m.x[1] + 2*m.x[2] + 3*m.x[3]

    # replace_expressions accepts a map from expression (including just a
    # variable) ids to the "target" expression that the original should
    # be replaced with.
    substitution_map = {
        id(m.x[2]): m.x[1]**2 + 1/m.x[3]
    }

    new_expr = replace_expressions(expr, substitution_map)

    print("Before replacement") 
    print("------------------") 
    print(expr)

    print("After replacement") 
    print("-----------------") 
    print(new_expr)


if __name__ == "__main__":
    main()
