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

from pyomo.common.timing import TicTocTimer
from pyomo.core.expr.numvalue import NumericValue
from pyomo.core.expr.visitor import StreamBasedExpressionVisitor
from pyomo.core.base.constraint import Constraint
from pyomo.core.base.objective import Objective
from pyomo.core.base.expression import ExpressionData
from pyomo.repn.ampl import AMPLRepnVisitor, NLFragment
from pyomo.repn.plugins.nl_writer import (
    AMPLRepn,
)
from pyomo.repn.util import FileDeterminism, FileDeterminism_to_SortComponents

import pyomo
pyomo_version = pyomo.version.version_info[:3]
pyomo_ge_673 = pyomo_version >= (6, 7, 3)
pyomo_ge_680 = pyomo_version >= (6, 8, 0)

from collections import namedtuple
# I don't really care about the constant node. We either have zero or one
# constant nodes... not a huge difference.
NodeCount = namedtuple("NodeCount", ["linear", "nonlinear"])


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
        if (
            expr.is_named_expression_type()
            and not self._descend_into_named_expressions
        ):
            self.count += 1
            # self.count is returned as the result of this walk. We increment
            # here as enterNode is not called.
            return False, self.count
        else:
            return True, None

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
    if amplrepn:
        subexpression_cache = {}
        subexpression_order = []
        external_functions = {}
        var_map = {}
        used_named_expressions = set()
        symbolic_solver_labels = False
        export_defined_variables = True
        sorter = FileDeterminism_to_SortComponents(FileDeterminism.ORDERED)
        if pyomo_ge_680:
            visitor_args = (
                subexpression_cache,
                external_functions,
                var_map,
                used_named_expressions,
                symbolic_solver_labels,
                export_defined_variables,
                sorter,
            )
        elif pyomo_ge_673:
            visitor_args = (
                text_nl_template,
                subexpression_cache,
                subexpression_order,
                external_functions,
                var_map,
                used_named_expressions,
                symbolic_solver_labels,
                export_defined_variables,
                sorter,
            )
        else:
            visitor_args = (
                text_nl_template,
                subexpression_cache,
                external_functions,
                var_map,
                used_named_expressions,
                symbolic_solver_labels,
                export_defined_variables,
                sorter,
            )
        visitor = AMPLRepnVisitor(*visitor_args)
    else:
        visitor = NodeCounter(descend_into_named_expressions=False)

    count = 0
    component_exprs = (
        [con.body for con in model.component_data_objects(Constraint, active=True)]
        + [obj.expr for obj in model.component_data_objects(Objective, active=True)]
    )
    for expr in component_exprs:
        if amplrepn:
            expr_cache = subexpression_cache
            count += count_amplrepn_nodes(
                expr,
                visitor=visitor,
                **kwds,
            )
        else:
            count += visitor.walk_expression(expr)

    if amplrepn:
        expr_ids = list(expr_cache)
        while expr_ids:
            e_id = expr_ids.pop()

            e_obj, repn, _ = expr_cache[e_id]
            if isinstance(e_obj, NLFragment):
                # NLFragment objects store the nonlinear portion of a named
                # expression, as linear and nonlinear portions are written
                # separately in the nl file. We skip this, as we will encounter
                # and process the full named expression later.
                #
                # Add 1 to account for the indirection that occurs between
                # the first and second portions of the subexpression.
                # This is for done for consistency with how we count subexpressions
                # elsewhere.
                count += 1
                continue

            if isinstance(e_obj, NLFragment):
                count += 1
                continue

            # NOTE: The named expression subtree will replace at least one node
            # in each nonlinear constraint where it appears. We don't attempt
            # to account for this.

            # Re-set subexpression cache so we know which expressions are the
            # new ones.
            visitor.subexpression_cache = {}
            new_expr_cache = visitor.subexpression_cache
            count += count_amplrepn_nodes(
                e_obj.expr,
                visitor=visitor,
                **kwds,
            )
            for new_e_id in new_expr_cache:
                if new_e_id not in expr_cache:
                    # Update "global" dict with newly discovered
                    # subexpressions
                    expr_cache[new_e_id] = new_expr_cache[new_e_id]
                    # Push to the top of our stack
                    expr_ids.append(new_e_id)
    else:
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
            # We "already counted" the named expression node itself, so we just
            # count the expression that defines it here. Since we have
            # visitor._descend_into_named_expressions==False, counting the
            # expression itself would give us a trivial result.
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
        # Are we double-counting expressions here?

    return count


def count_amplrepn_nodes(
    expr,
    export_defined_variables=True,
    expression_cache=None,
    visitor=None
):
    """
    Use export_defined_variables=False to descend into named expressions while
    computing nodes.

    Side effect: If expression_cache (dict) is provided, it will be updated
    with the named expressions found in the provided expr.

    """
    if visitor is None:
        local_subexpression_cache = {}
        subexpression_order = []
        external_functions = {}
        var_map = {}
        used_named_expressions = set()
        symbolic_solver_labels = False
        sorter = FileDeterminism_to_SortComponents(FileDeterminism.ORDERED)
        if pyomo_ge_680:
            visitor_args = (
                local_subexpression_cache,
                external_functions,
                var_map,
                used_named_expressions,
                symbolic_solver_labels,
                export_defined_variables,
                sorter,
            )
        elif pyomo_ge_673:
            visitor_args = (
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
        else:
            visitor_args = (
                text_nl_template,
                local_subexpression_cache,
                external_functions,
                var_map,
                used_named_expressions,
                symbolic_solver_labels,
                export_defined_variables,
                sorter,
            )
        visitor = AMPLRepnVisitor(*visitor_args)
    AMPLRepn.ActiveVisitor = visitor
    try:
        repn = visitor.walk_expression((expr, None, 0, 1.0))
        local_subexpression_cache = visitor.subexpression_cache
    finally:
        AMPLRepn.ActiveVisitor = None

    n_linear_nodes = len(repn.linear)

    if repn.nonlinear is not None:
        # We count nonlinear nodes by counting the number of newlines in the
        # nl string.
        #
        # NOTE: As of Pyomo 6.7.3, variables nodes in the nl template string
        # (i.e. "%s") no longer have newlines appended. We now need to count these
        # substrings as well.
        n_nonlinear_nodes = (
            repn.nonlinear[0].count("\n")
            + repn.nonlinear[0].count("%s")
        )
    else:
        n_nonlinear_nodes = 0

    if (
        expression_cache is not None
        and export_defined_variables
        and repn.named_exprs # Is not None or an empty set
    ):
        # NOTE: To deal with the splitting of subexpressions between NLFragments
        # and linear parts, we need to only update the cache with subexpressions
        # that are actually used.
        # The possibilities are:
        # - the subexpression participates nonlinearly => add NLFragment and expression
        # - The subexpression participates linearly => add only NLFragment
        used_local_subexpressions = {}
        for e_id in repn.named_exprs:
            # Is e_id guaranteed to exist in this cache? I believe so.
            data = local_subexpression_cache[e_id]
            if isinstance(data[0], ExpressionData):
                used_local_subexpressions[e_id] = data
            elif isinstance(data[0], NLFragment):
                # A new NLFragment is constructed every time we encounter
                # the same subexpression, so using the id as a key will result
                # in this expression be counted multiple times. We need a key
                # that is unique to this expression, so we use the NLFragment's
                # name, which is e.g. "nl(subexpression[1])"
                # The key doesn't really matter here, as we'll iterate over
                # the values later to count the subexpression nodes.
                used_local_subexpressions[data[0].name] = data
            else:
                raise TypeError(f"Unrecognized expression type {data[0].__class__}")
        expression_cache.update(used_local_subexpressions)

    nodecount = NodeCount(n_linear_nodes, n_nonlinear_nodes)
    return nodecount


def count_model_amplrepn_nodes(model, **kwds):
    # Initialize walker
    subexpression_cache = {}
    subexpression_order = []
    external_functions = {}
    var_map = {}
    used_named_expressions = set()
    symbolic_solver_labels = False
    export_defined_variables = True
    sorter = FileDeterminism_to_SortComponents(FileDeterminism.ORDERED)
    if pyomo_ge_680:
        visitor_args = (
            subexpression_cache,
            external_functions,
            var_map,
            used_named_expressions,
            symbolic_solver_labels,
            export_defined_variables,
            sorter,
        )
    elif pyomo_ge_673:
        visitor_args = (
            text_nl_template,
            subexpression_cache,
            subexpression_order,
            external_functions,
            var_map,
            used_named_expressions,
            symbolic_solver_labels,
            export_defined_variables,
            sorter,
        )
    else:
        visitor_args = (
            text_nl_template,
            subexpression_cache,
            external_functions,
            var_map,
            used_named_expressions,
            symbolic_solver_labels,
            export_defined_variables,
            sorter,
        )
    visitor = AMPLRepnVisitor(*visitor_args)

    component_exprs = (
        [con.body for con in model.component_data_objects(Constraint, active=True)]
        + [obj.expr for obj in model.component_data_objects(Objective, active=True)]
    )
    n_linear_nodes = 0
    n_nonlinear_nodes = 0

    # We'll build up a new cache of subexpressions, out of only the subexpressions
    # that are actually used
    used_expr_cache = {}
    for expr in component_exprs:
        nodecount = count_amplrepn_nodes(
            expr,
            # NOTE: This is too slow if I don't pass in the visitor here
            visitor=visitor,
            export_defined_variables=True,
            expression_cache=used_expr_cache,
            **kwds,
        )
        # Note that these linear nodes include any linear nodes from subexpressions
        # that participate linearly in this component.
        n_linear_nodes += nodecount.linear
        n_nonlinear_nodes += nodecount.nonlinear

    expr_cache = used_expr_cache
    expr_ids = list(expr_cache)

    for data in expr_cache.values():
        e_obj, repn, _ = data
        if isinstance(e_obj, NLFragment):
            # NLFragment should never have a linear part
            assert repn.linear is None
            # NLFragment should always have a nonlinear part
            assert repn.nonlinear is not None
            n_nonlinear_nodes += (
                repn.nonlinear[0].count("\n")
                + repn.nonlinear[0].count("%s")
            )
        elif isinstance(e_obj, ExpressionData):
            # We are using the full expression. This means it appears nonlinearly
            # somewhere.
            # We count linear terms and nonlinear nodes equally as nonlinear nodes
            # in the parent expression
            if repn.linear is not None:
                n_nonlinear_nodes += len(repn.linear)
            if repn.nonlinear is not None:
                n_nonlinear_nodes += (
                    repn.nonlinear[0].count("\n")
                    + repn.nonlinear[0].count("%s")
                )
        else:
            raise TypeError(f"Unrecognized expression type {e_obj.__class__}")

    model_nodecount = NodeCount(n_linear_nodes, n_nonlinear_nodes)
    return model_nodecount
