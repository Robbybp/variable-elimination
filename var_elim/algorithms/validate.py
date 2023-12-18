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

from pyomo.core.expr import value as pyo_value
from pyomo.core.base.var import Var
from idaes.core.util.model_statistics import large_residuals_set


def validate_solution(
    model,
    eliminated_var_exprs,
    eliminated_constraints,
    tolerance=0.0,
):
    violated_cons_reduced = large_residuals_set(model, tol=tolerance)

    # Set variables to value defined by elimination expression
    # We assume these expressions are in terms of reduced-space variables.
    for var, expr in eliminated_var_exprs:
        var.set_value(pyo_value(expr))

    vars_violating_bounds = []
    for var in model.component_data_objects(Var):
        if var.value is None:
            continue
        if var.ub is not None:
            ub_diff = pyo_value(var.value - var.ub)
            if ub_diff > tolerance:
                vars_violating_bounds.append((var, var.ub, ub_diff))
        if var.lb is not None:
            lb_diff = pyo.value(var.value - var.lb)
            if lb_diff < - tolerance:
                vars_violating_bounds.append((var, var.lb, lb_diff))

    violated_eliminated_cons = []
    for con in eliminated_constraints:
        # Should be all equality constraints
        resid = pyo_value(con.body - con.upper)
        if resid > tolerance:
            violated_eliminated_cons.append((con, resid))

    if violated_cons_reduced:
        print("WARNING: Constraints in the reduced-space model are violated")

    if vars_violating_bounds:
        print("WARNING: There are variables violating their bounds")

    if violated_eliminated_cons:
        print("WARNING: Eliminated constraints are violated")

    violations = (violated_cons_reduced, vars_violating_bounds, violated_eliminated_cons)
    valid = not any(violations)

    return valid, violations
