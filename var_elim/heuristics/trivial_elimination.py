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

from pyomo.core.base.constraint import Constraint
from pyomo.repn import generate_standard_repn
from pyomo.core.expr import EqualityExpression
from pyomo.util.subsystems import create_subsystem_block
from var_elim.heuristics.matching import generate_elimination_via_matching


def expr_filter(
    expr,
    degree=None,
    linear=None,
    affine=None,
    equal_coefficients=None,
):
    # This function does not actually need the expressions generated
    # by standard repn. TODO: explore using ampl_repn for consistency.
    repn = generate_standard_repn(
        expr,
        # We compute values here to try to filter out zero terms in nonlinear
        # subexpression.
        compute_values=True,
        quadratic=False,
    )
    con_is_affine = (len(repn.nonlinear_vars) == 0)
    if linear is not None:
        # We filter according to whether the constraint's linearity
        # matches the provided flag.
        # Note that we follow the typical math programming convention where
        # "linear" means "linear or affine".
        if linear != con_is_affine:
            # If we fail the linearity check, return immediately.
            # Otherwise, we continue with other checks.
            return False
    if affine is not None:
        # This specifically checks whether the constraint is affine,
        # i.e. "affine and not linear"
        con_is_linear = (con_is_affine and (repn.constant == 0))
        if affine != (con_is_affine and not con_is_linear):
            # We fail the affine check
            return False
    if equal_coefficients is not None and len(repn.linear_coefs) >= 1:
        # Note that this only makes sense if we have linear variables.
        # It also arguably doesn't make sense if we have any nonlinear variables.
        # Note that we only care whether the *magnitude* of the coefficients is
        # equal. This flag could probably use a better name...
        coef = repn.linear_coefs[0]
        if (
            equal_coefficients
            != all(abs(c) == abs(coef) for c in repn.linear_coefs)
        ):
            return False
    if degree is not None:
        # TODO: When/if this goes into Pyomo, we should make sure that this
        # method of determining degree is consistent with get_incidence_variables.
        linear_var_ids = set(id(var) for var in repn.linear_vars)
        nonlinear_var_ids = set(id(var) for var in repn.nonlinear_vars)
        var_ids = (linear_var_ids | nonlinear_var_ids)
        con_degree = len(var_ids)
        if con_degree != degree:
            return False
    return True


def filter_constraints(
    model,
    active=True,
    include_inequality=False,
    degree=None,
    linear=None,
    affine=None,
    equal_coefficients=None,
):
    return list(filter(
        lambda constraint: (
            include_inequality or isinstance(constraint.expr, EqualityExpression)
            and expr_filter(
                constraint.body,
                degree=degree,
                linear=linear,
                affine=affine,
                equal_coefficients=equal_coefficients,
            )
        ),
        model.component_data_objects(Constraint, active=active),
    ))


def get_trivial_constraint_elimination(model):
    trivial_cons = filter_constraints(
        model,
        degree=2,
        linear=True,
        affine=False,
        equal_coefficients=True,
    )
    temp_block = create_subsystem_block(trivial_cons)
    return generate_elimination_via_matching(temp_block)


if __name__ == "__main__":
    from var_elim.models.distillation.distill import create_instance
    from var_elim.algorithms.replace import define_elimination_order

    m = create_instance()
    var_order, con_order = get_trivial_constraint_elimination(m)
    var_order, con_order = define_elimination_order(var_order, con_order)
