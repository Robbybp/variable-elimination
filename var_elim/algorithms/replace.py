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

from pyomo.core.base.objective import Objective
from pyomo.core.base.constraint import Constraint
from pyomo.core.base.expression import Expression
from pyomo.core.base.set import Set, Integers, Binary
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.repn import generate_standard_repn
from pyomo.core.expr.relational_expr import EqualityExpression
from pyomo.core.expr.visitor import (
    replace_expressions,
    identify_variables,
    ExpressionReplacementVisitor,
)
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.config import IncidenceMethod
from pyomo.common.modeling import unique_component_name
from pyomo.common.timing import HierarchicalTimer


def define_variable_from_constraint(variable, constraint, timer=None):
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
    timer = HierarchicalTimer() if timer is None else timer
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
    timer.start("standard_repn")
    repn = generate_standard_repn(
        constraint.body, compute_values=False, quadratic=False
    )
    timer.stop("standard_repn")
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


def define_elimination_order(var_list, con_list, igraph=None):
    """
    Returns elimination order using block triangularize from incidence graph interface

    """
    if igraph is None:
        # Use ampl_repn here as we don't want spurious nonzeros to mistakenly
        # induce an algebraic loop
        igraph = IncidenceGraphInterface(method=IncidenceMethod.ampl_repn)

    var_blocks, con_blocks = igraph.block_triangularize(var_list, con_list)
    for vb, cb in zip(var_blocks, con_blocks):
        assert len(vb) == 1
        assert len(cb) == 1
    var_order = sum(var_blocks, [])
    con_order = sum(con_blocks, [])
    return var_order, con_order


def add_bounds_to_expr(var, var_expr):
    #Need to change doc string of this function
    
    #We cannot write expressions which lead to trivial booleans
    # Eg. 0 >= -2 or 0 <= 2
    #Ths can happen when var_expr is a constant and the var lb and ub contain that constant
    """
    This function takes in a variable, the expression for variable replacement
    and an indexed constraint list - bound_cons. It updates the list with inequality
    constraints on the expression if the variables replaced were bounded

    Each constraint added to the bound_cons list is indexed by var_name_ub or
    var_name_lb depending upon whihc bound it adds to the expression
    """
    if var.ub is None and var.lb is None:
        lb_expr = None
        ub_expr = None
    elif var.lb is not None and var.ub is None:
        lb_expr = var_expr >= var.lb
        ub_expr = None
    elif var.ub is not None and var.lb is None:
        lb_expr = None
        ub_expr = var_expr <= var.ub
    else:
        lb_expr = var_expr >= var.lb
        ub_expr = var_expr <= var.ub
    
    return lb_expr, ub_expr


def _get_elimination_map(
    m,
    igraph,
    variables,
    constraints,
    use_named_expressions=False,
):
    subgraph = igraph.subgraph(variables, constraints)

    elim_var_set = Set(initialize=[])
    m.add_component(
        unique_component_name(m, "replaced_variable_set"), elim_var_set
    )
    elim_var_expr = Expression(elim_var_set)
    m.add_component(
        unique_component_name(m, "eliminated_variable_expressions"), elim_var_expr
    )

    # Assume variables/constraints are in lower triangular order
    substitution_map = {}
    to_replace = set()
    for var, con in zip(variables, constraints):
        # If necessary, replace expressions in constraint, before
        # using it to generate an expression defining the var
        if id(con) in to_replace:
            to_replace.remove(id(con))
            new_expr = replace_expressions(
                con.expr,
                substitution_map,
                descend_into_named_expressions=True,
                remove_named_expressions=False,
            )
            con.set_value(new_expr)

        # Define variable from the constraint and update the substitution map
        var_expr = define_variable_from_constraint(var, con)
        elim_var_set.add(var.name)
        elim_var_expr[var.name] = var_expr
        if use_named_expressions:
            var_expr = elim_var_expr[var.name]
        substitution_map[id(var)] = var_expr

        # Add other constraints in which the variable appears to the set of
        # constraints we must process
        for adj_con in subgraph.get_adjacent_to(var):
            if adj_con is not con:
                to_replace.add(id(adj_con))
    return substitution_map, elim_var_set, elim_var_expr


def eliminate_nodes_from_graph(igraph, variables, constraints, timer=None):
    """Update an incidence graph to account for eliminated variables and
    constraints

    Nodes corresponding to eliminated variables and constraints are removed from
    the graph. For every eliminated variable, constraint pair, an edge is drawn
    between all pairs of adjacent (non-eliminated) constraints and variables.

    Parameters
    ----------

    igraph: IncidenceGraphInterface
        Graph to modify in-place.

    variables: list of VarData
        Variables to eliminated

    constraints: list of ConData
        Constraints to eliminate

    """
    if timer is None:
        timer = HierarchicalTimer()
    if len(variables) != len(constraints):
        raise RuntimeError("Dimension mismatch")
    timer.start("add_edge")
    for var, con in zip(variables, constraints):
        adj_vars = igraph.get_adjacent_to(con)
        adj_cons = igraph.get_adjacent_to(var)
        for adjvar in adj_vars:
            for adjcon in adj_cons:
                igraph.add_edge(adjvar, adjcon)
    timer.stop("add_edge")
    timer.start("remove_nodes")
    igraph.remove_nodes(variables, constraints)
    timer.stop("remove_nodes")


def eliminate_variables(
    m,
    var_order,
    con_order,
    igraph=None,
    linear_igraph = None,
    eq_igraph = None,
    use_named_expressions=False,
    timer=None,
):
    """
    Does the actual elimination by defining variable from constraint, deactivating
    the constraint used for variable definition, and replacing the variable in
    every adjacent constraint

    Parameters
    ----------

    m: ConcreteModel
        Model whose variables and constraints will be eliminated

    var_order: list of VarData
        Variables to eliminate in the order provided

    con_order: list of ConData
        Equality constraints used to eliminate the variables, in order.
        These constraints will also be eliminated.

    igraph: IncidenceGraphInterface
        Incidence graph of the model. This graph should include inequality
        constraints. It is used to determine which constraints contain the
        variables that have been eliminated, so these variables can be replaced
        by the corresponding expressions.

    Returns
    -------
    Reduced Model

    """
    if timer is None:
        timer = HierarchicalTimer()
    timer.start("eliminate_variables")
    for var in var_order:
        if var.domain is Integers or var.domain is Binary:
            raise RuntimeError(
                f"Cannot eliminate discrete variable {var.name}"
            )

    # TODO: This would not be necessary if IncidenceGraphInterface supported
    # objectives as nodes.
    #
    # This would get slightly simpler if we only support single-objective
    # problems...
    objectives = list(m.component_data_objects(Objective, active=True))
    var_obj_map = ComponentMap()
    for obj in objectives:
        # We do not need include_fixed here. These variables do not determine
        # the variables that will be replaced. Variables in the objective
        # with a coefficient of zero will not be a problem either (the replaced
        # expression will just have a coefficient of zero).
        for var in identify_variables(obj.expr):
            if var in var_obj_map:
                var_obj_map[var].append(obj)
            else:
                var_obj_map[var] = [obj]

    # Set that will store names of bounds on replacement expressions
    bound_con_set = Set(initialize=[])
    m.add_component(
        unique_component_name(m, "replaced_variable_bounds_set"), bound_con_set
    )
    # Constraint that will store bounds on replacement expressions
    bound_con = Constraint(bound_con_set)
    m.add_component(
        unique_component_name(m, "replaced_variable_bounds"), bound_con
    )

    var_lb_map = ComponentMap()
    var_ub_map = ComponentMap()
    var_exprs = []

    # Including inequalities in this incidence graph replaces variables in the
    # adjacent inequality constraints too. If the user supplies an igraph,
    # it needs to have the inequality constraints included
    timer.start("igraph")
    igraph_provided = igraph is not None
    if not igraph_provided:
        igraph = IncidenceGraphInterface(
            m,
            include_inequality=True,
            # Use ampl_repn as we don't want to do extra work due to spurious
            # nonzeros (and introduce more spurious nonzeros).
            method=IncidenceMethod.ampl_repn,
        )
    timer.stop("igraph")

    substitution_map, elim_var_set, elim_var_expr = _get_elimination_map(
        m, igraph, var_order, con_order, use_named_expressions=use_named_expressions
    )

    # Deactivate constraints that define variables
    elim_con_set = set(id(con) for con in con_order)
    for con in con_order:
        con.deactivate()

    # Update data structures for eliminated variables
    for var in var_order:
        var_expr = substitution_map[id(var)]
        lb_expr, ub_expr = add_bounds_to_expr(var, var_expr)
        lb_name = var.name + "_lb"
        ub_name = var.name + "_ub"
        if lb_expr is not None and type(lb_expr) is not bool:
            if lb_expr is False:
                raise RuntimeError("Lower bound resolved to trivial infeasible constraint")
            bound_con_set.add(lb_name)
            bound_con[lb_name] = lb_expr
            var_lb_map[var] = bound_con[lb_name]
        if ub_expr is not None and type(ub_expr) is not bool:
            if ub_expr is False:
                raise RuntimeError("Upper bound resolved to trivial infeasible constraint")
            bound_con_set.add(ub_name)
            bound_con[ub_name] = ub_expr
            var_ub_map[var] = bound_con[ub_name]
        var_exprs.append((var, var_expr))

    # Replace eliminated variables in the rest of the constraints
    replacement_visitor = ExpressionReplacementVisitor(
        substitute=substitution_map,
        descend_into_named_expressions=True,
        remove_named_expressions=False,
    )
    # TODO: There is a potential "sparse" implementation of this loop, where we
    # iterate over variables-to-eliminate, and perform the substitution for each
    # (new) constraint that they're adjacent to. This way we don't check every
    # constraint if we're only eliminating a small number of variables.
    for con in igraph.constraints:
        if (
            id(con) not in elim_con_set
            and any(id(var) in substitution_map for var in igraph.get_adjacent_to(con))
        ):
            new_expr = replacement_visitor.walk_expression(con.expr)
            con.set_value(new_expr)

    for obj in m.component_data_objects(Objective, active=True):
        new_expr = replacement_visitor.walk_expression(obj.expr)
        obj.set_value(new_expr)

    if igraph_provided:
        eliminate_nodes_from_graph(igraph, var_order, con_order)

    timer.stop("eliminate_variables")
    return var_exprs, var_lb_map, var_ub_map


if __name__ == "__main__":
    import pyomo.environ as pyo

    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3])
    m.p = pyo.Param([1, 2], initialize={1: 10, 2: 20}, mutable=True)
    expr = (
        3 * m.p[1] * m.x[1]
        + 4 * m.p[1] * m.x[2]
        - 5 * m.p[2] * m.x[3]
        + m.p[2] * m.x[1] ** 2
        + m.x[3] * pyo.exp(-m.p[2] / m.x[1])
    )
    m.con = pyo.Constraint(expr=expr == (m.x[2] + m.p[1] - m.p[2]) / 2)

    x2_expr = define_variable_from_constraint(m.x[2], m.con)
    print("Constraint expression")
    print("---------------------")
    print(m.con.expr)
    print()
    print("Expression defining x[2]")
    print("------------------------")
    print(x2_expr)
