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
from pyomo.environ import Constraint
from pyomo.core.expr.visitor import identify_variables
from pyomo.repn import generate_standard_repn
from pyomo.core.base.set import Integers, Binary
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.current import EqualityExpression
import random
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface


def identify_vars_for_elim_ampl(m, 
                                randomize = False, 
                                eliminate_bounded_vars = False, 
                                eliminate_linear_cons_only = False,
                                eliminate_linear_vars_only = False):
    """Identify defined variables and defining constraints via the heuristic
    from the AMPL preprocessor

    Parameters
    ----------
    m : Pyomo Model

    Returns
    -------
    var_list : List of variables to be eliminated
    con_list : List of constraints used to eliminate the variables

    """
    #Get active equality constraints from the model
    cons = []
    for c in m.component_data_objects(Constraint, active=True):
        if isinstance(c.expr, EqualityExpression):
            cons.append(c)
    
    #Shuffle constraints for randomizing constraint order
    if randomize == True:
        random.shuffle(cons)
   
    #Generate full incidence graph and linear incidence graph only if 
    #we need to eliminate only linear vars
    if eliminate_linear_vars_only:
        igraph = IncidenceGraphInterface(m, active=True, include_inequality= True)
        linear_igraph = IncidenceGraphInterface(m, active = True, linear_only=True, include_inequality = False)
    
    # Identify variables of the type ===> coef*v (+/-) expr = 0
    var_list = []
    con_list = []
    defining_var_ids = set()
    for c in cons:
        # gets all vars in order from left to right in the constraint expression
        expr_vars = list(identify_variables(c.expr))

        # gets linear vars in the constraint expression
        repn = generate_standard_repn(c.body, compute_values=False, quadratic=False)
        linear_vars = ComponentSet(repn.linear_vars)

        nonlinear_vars = ComponentSet(repn.nonlinear_vars)

        if expr_vars[0].domain is Integers or expr_vars[0].domain is Binary:
            pass
        elif not eliminate_bounded_vars and (expr_vars[0].lb is not None or expr_vars[0].ub is not None):
            pass
        elif expr_vars[0] in nonlinear_vars:
            pass
        elif eliminate_linear_cons_only and len(nonlinear_vars) != 0:
            pass
        elif eliminate_linear_vars_only and (len(igraph.get_adjacent_to(expr_vars[0])) != len(linear_igraph.get_adjacent_to(expr_vars[0]))):
            pass
        elif id(expr_vars[0]) not in defining_var_ids:
            if expr_vars[0] in linear_vars:
                var_list.append(expr_vars[0])
                con_list.append(c)
                for var in expr_vars:
                    defining_var_ids.add(id(var))

                # This will add all vars from the expression to the defining vars list
                # This works because we anyways don't want to eliminate the same var twice
                # using 2 different constraints

    return var_list, con_list
