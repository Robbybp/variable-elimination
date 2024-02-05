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

from pyomo.core.expr.numvalue import NumericValue
from pyomo.core.expr.visitor import StreamBasedExpressionVisitor
from pyomo.core.base.constraint import Constraint
from pyomo.core.base.objective import Objective
from pyomo.core.base.expression import Expression
from pyomo.repn.plugins.nl_writer import (
    AMPLRepnVisitor, AMPLRepn, text_nl_template, NLFragment
)
from pyomo.repn.util import FileDeterminism, FileDeterminism_to_SortComponents


class NodeCounter(StreamBasedExpressionVisitor):

    def __init__(self, descend_into_named_expressions=True):
        super().__init__()
        self._descend_into_named_expressions = descend_into_named_expressions
        # Note that these are only used if we are not descending
        # These are intended to store the named expression objects so we can
        # count their nodes later.
        self.named_expressions = []
        self.named_expr_map = {}

    def initializeWalker(self, expr):
        self.count = 0
        return True, expr

    def beforeChild(self, parent, child, idx):
        if (
            not self._descend_into_named_expressions
            and isinstance(child, NumericValue)
            and child.is_named_expression_type()
        ):
            # Because we will not enter the child node, we need to update
            # the count here.
            self.count += 1
            if id(child) not in self.named_expr_map:
                self.named_expr_map[id(child)] = child
                self.named_expressions.append(child)
            return False, None
        else:
            return True, None

    def enterNode(self, node):
        self.count += 1
        return None

    def finalizeResult(self, result):
        return self.count


def count_nodes(expr, **kwds):
    visitor = NodeCounter(**kwds)
    return visitor.walk_expression(expr)


def count_model_nodes(
    model,
    amplrepn=False,
    **kwds,
):
    # TODO: Separate functions for amplrepn vs. Pyomo (that can be called by this
    # function for convenience).
    #
    # When counting nodes for the entire model, I think we always want to consider
    # named expressions independently.
    #descend_into_named_expressions = kwds.pop("descend_into_named_expressions", True)
    if kwds and not amplrepn:
        raise RuntimeError("kwds not supported with amplrepn=False")
    visitor = NodeCounter(descend_into_named_expressions=False)

    count = 0
    component_exprs = (
        [con.body for con in model.component_data_objects(Constraint, active=True)]
        + [obj.expr for obj in model.component_data_objects(Objective, active=True)]
    )
    for expr in component_exprs:
        count += visitor.walk_expression(expr)

    # This is the stack of expressions we still need to process.
    expr_stack = list(visitor.named_expressions)
    # This is the set of all expressions that have ever been added to the
    # stack.
    encountered_expr_ids = set(visitor.named_expr_map)
    while expr_stack:
        named_expr = expr_stack.pop()

        # Walk expression to count nodes in this named expression.
        # Additionally, this adds any new named expressions
        #
        # Clear named expressions so we know which ones were discovered
        # this iteration. This avoids quadratic scaling as we build up
        # lots of named expressions.
        visitor.named_expr_map.clear()
        visitor.named_expressions = []
        count += visitor.walk_expression(named_expr.expr)

        # Add new expressions to the "global set"
        for e_id, new_expr in visitor.named_expr_map.items():
            if e_id not in encountered_expr_ids:
                # encountered_expr_ids stays in-sync with the visitor's set
                # of expressions, but we need to maintain these two sets
                # so we know which expressions to add to the stack.
                # We could just re-construct (or clear) the visitor, but then
                # we need to maintain a separate count. Can consider this for
                # performance later.
                encountered_expr_ids.add(e_id)
                # need to link the id to the actual expression
                expr_stack.append(new_expr)

    return count
