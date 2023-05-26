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

from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.repn import generate_standard_repn
from pyomo.core.expr.relational_expr import EqualityExpression
from pyomo.core.expr.visitor import replace_expressions
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface


def define_variable_from_constraint(variable, constraint):
    """Get the expression that defines the variable according to the
    constraint.

    This only works if the variable participates linearly in the constraint.
    This should be used to generate substitution maps for
    ``replace_expressions``.

    Returns
    -------
    NumericExpression
        Defines variable using the expression in the constraint

    """
    if not isinstance(constraint.expr, EqualityExpression):
        raise RuntimeError(
            f"{constraint.name} does not contain an EqualityExpression."
            f" Note that ranged inequalities with lower==upper are not"
            f" supported at this time."
        )
    # Generate standard repn 
    # - compute_values=False so we identify linear coefficients in terms of
    #   whatever parameters (or fixed variables) are present. This is important
    #   for reconstructing the full expression
    # - quadratic=False to avoid unnecessary complexity
    repn = generate_standard_repn(
        constraint.body, compute_values=False, quadratic=False
    )
    var_coef_list = zip(repn.linear_vars, repn.linear_coefs)
    linear_coef_map = ComponentMap(var_coef_list)
    nonlinear_vars = ComponentSet(repn.nonlinear_vars)
    if variable not in linear_coef_map.keys():
        raise RuntimeError(
            f"{variable.name} does not participate linearly in {constraint.name}!"
        )
    elif variable in nonlinear_vars:
        raise RuntimeError(
            f"{variable.name} participates nonlinearly in {constraint.name}!"
        )
    coef = linear_coef_map[variable]
    # coef will never be a constant of zero, but it could be an expression of
    # parameters that evaluates to zero. TODO: raise an error if this happens.

    other_linear_vars = list(filter(lambda v: v is not variable, repn.linear_vars))
    linear_subexpr = sum(linear_coef_map[v] * v for v in other_linear_vars)

    # constraint: const + coef * var + sum(coef * other_vars) + nonlinear == upper
    # =>
    # var := (1/coef) * (upper - const - sum(coef * other_vars) - nonlinear)
    
    nonlinear_expr = repn.nonlinear_expr if repn.nonlinear_expr is not None else 0.0
   
    var_expr = (1 / coef) * (
        constraint.upper  # Either constraint.upper or constraint.lower would work
        - repn.constant
        - linear_subexpr
        - nonlinear_expr
    )
   
    return var_expr

def define_elimination_order(var_list, con_list, igraph = None):
    """
    Returns elimination order using block triangularize from incidence graph interface
    
    """
    if igraph is None:
        igraph = IncidenceGraphInterface()
    
    var_dmp, _ = igraph.dulmage_mendelsohn(var_list, con_list)
    assert not var_dmp.unmatched
    
    var_blocks, con_blocks = igraph.block_triangularize(var_list, con_list)
    for vb, cb in zip(var_blocks, con_blocks):
        assert len(vb) == 1
        assert len(cb) == 1 
    var_order = sum(var_blocks, [])
    con_order = sum(con_blocks, [])
    return var_order, con_order


def eliminate_variables(m, var_order, con_order, igraph = None):
    """
    Does the actual elimination by defining variable from constraint, deactivating
    the constraint used for variable definition, and replacing the variable in
    every adjacent constraint

    Returns
    -------
    Reduced Model

    """
    if igraph is None:
        igraph = IncidenceGraphInterface(m, include_inequality = False)
        
    for var, con in zip(var_order, con_order):
        #Get expression for the variable from constraint
        var_expr = define_variable_from_constraint(var, con)
        con.deactivate()
        
        #Build substitution map
        substitution_map = {id(var): var_expr}
        
        #Get constraints in which the variable appears
        #This will have the deactivated constraints too
        adj_cons = igraph.get_adjacent_to(var)
        for ad_con in adj_cons:
            if ad_con is not con: 
                new_expr = replace_expressions(ad_con.expr, substitution_map)
                ad_con.set_value(new_expr)
    return m



if __name__ == "__main__":
    import pyomo.environ as pyo
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3])
    m.p = pyo.Param([1, 2], initialize={1: 10, 2: 20}, mutable=True)
    expr = (
        3*m.p[1]*m.x[1] + 4*m.p[1]*m.x[2] - 5*m.p[2]*m.x[3]
        + m.p[2]*m.x[1]**2 + m.x[3]*pyo.exp(-m.p[2]/m.x[1])
    )
    m.con = pyo.Constraint(expr=expr == (m.x[2] + m.p[1] - m.p[2])/2)

    x2_expr = define_variable_from_constraint(m.x[2], m.con)
    print("Constraint expression")
    print("---------------------")
    print(m.con.expr)
    print()
    print("Expression defining x[2]")
    print("------------------------")
    print(x2_expr)
