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

from pyomo.core.expr.visitor import StreamBasedExpressionVisitor
from pyomo.core.base.constraint import Constraint
from pyomo.core.base.expression import Expression
from pyomo.repn.plugins.nl_writer import AMPLRepnVisitor, AMPLRepn, text_nl_template
from pyomo.repn.util import FileDeterminism, FileDeterminism_to_SortComponents


class NodeCounter(StreamBasedExpressionVisitor):

    def initializeWalker(self, expr):
        self._count = 0
        return True, expr

    def enterNode(self, node):
        self._count += 1
        return None

    def finalizeResult(self, result):
        return self._count


def count_nodes(expr):
    visitor = NodeCounter()
    return visitor.walk_expression(expr)


def count_model_nodes(
    model,
    amplrepn=False,
    **kwds,
):
    if kwds and not amplrepn:
        raise RuntimeError("kwds not supported with amplrepn=False")
    count = 0
    for con in model.component_data_objects(Constraint, active=True):
        if amplrepn:
            expr_cache = {}
            count += count_amplrepn_nodes(con.body, expression_cache=expr_cache, **kwds)
        else:
            count += count_nodes(con.body)
    
    if amplrepn:
        expr_ids = list(expr_cache)
        while expr_ids:
            e_id = expr_ids.pop()
            # Why is there sometimes a third object in expr_cache?
            e_obj, repn, _ = expr_cache[e_id]

            new_expr_cache = {}
            # NOTE: The named expression subtree will replace at least one node
            # in each nonlinear constraint where it appears. We don't attempt
            # to account for this.
            count += count_amplrepn_nodes(
                e_obj.expr, expression_cache=new_expr_cache, **kwds
            )
            for new_e_id in new_expr_cache:
                if new_e_id not in expr_cache:
                    # Update "global" dict with newly discovered
                    # subexpressions
                    expr_cache[new_e_id] = new_expr_cache[new_e_id]
                    # Push to the top of our stack
                    expr_ids.append(new_e_id)

    return count


def count_amplrepn_nodes(
    expr,
    export_defined_variables=True,
    expression_cache=None,
    linear_only=False,
):
    """
    Use export_defined_variables=False to descend into named expressions while
    computing nodes.

    Side effect: If expression_cache (dict) is provided, it will be updated
    with the named expressions found in the provided expr.

    """
    local_subexpression_cache = {}
    subexpression_order = []
    external_functions = {}
    var_map = {}
    used_named_expressions = set()
    symbolic_solver_labels = False
    sorter = FileDeterminism_to_SortComponents(FileDeterminism.ORDERED)
    visitor = AMPLRepnVisitor(
        text_nl_template,
        local_subexpression_cache,
        subexpression_order,
        external_functions,
        var_map,
        used_named_expressions,
        symbolic_solver_labels,
        export_defined_variables,
        sorter,
    )
    AMPLRepn.ActiveVisitor = visitor
    try:
        repn = visitor.walk_expression((expr, None, 0, 1.0))
    finally:
        AMPLRepn.ActiveVisitor = None

    count = 0

    if repn.const != 0.0:
        # Add the nodes for the +(const) operation
        count += 2

    # We model each linear term as +(*(coef, var)), which is four nodes
    count += max(
        0,
        # Subtract one as, for n terms, we only need n-1 (+) operations
        4 * len([vid for vid, coef in repn.linear.items() if coef != 0.0]) - 1,
    )

    if not linear_only:
        # If we only want to count linear nodes, we ignore the nonlinear
        # subexpression. We also ignore named subexpressions, as these are only
        # ever referenced in the nonlinear portion of the constraint expression.

        if repn.nonlinear is not None:
            # Add one node for the (+) that connects the linear and nonlinear
            # subexpressions
            count += 1
            # Each line is a new node in the nonlinear expression
            count += repn.nonlinear[0].count("\n")

        if (
            expression_cache is not None
            and export_defined_variables
            and repn.named_exprs # Is not an empty set
        ):
            expression_cache.update(local_subexpression_cache)

    return count
