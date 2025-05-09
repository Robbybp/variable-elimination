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

import pytest
import pyomo.environ as pyo
from var_elim.algorithms.expr import (
    NodeCounter,
    count_nodes,
    count_model_nodes,
    count_amplrepn_nodes,
    count_model_amplrepn_nodes,
)
from var_elim.models.distillation.distill import create_instance as create_distill_instance
from var_elim.models.opf.opf_model import make_model as make_opf_model
from var_elim.algorithms.replace import (
    eliminate_variables, define_variable_from_constraint
)


class TestNodeCounter:

    def test_count_nodes_simple(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        expr = m.x[1] + m.x[2] + 2*m.x[3]

        n_nodes = count_nodes(expr)
        # Would have expected 6 nodes, but looks like x[1] and x[2] terms are
        # represented as MonomialTermExpressions...
        # This has been updated somewhere between Pyomo 6.7.2 and 6.7.3
        assert n_nodes == 6

    def test_count_nodes_distill_vol_expr(self):
        m = create_distill_instance(horizon=20, nfe=4)
        n_nodes = count_nodes(m.mole_frac_balance[1, 1].body)
        assert n_nodes == 14

    def test_count_nodes_opf_pf_expr(self):
        m = make_opf_model("ieee118")
        n_nodes = count_nodes(m.eq_pf_branch["1"].body)
        assert n_nodes == 21

    def test_count_nodes_no_descend(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        m.eq = pyo.Constraint(pyo.Integers)
        m.subexpr = pyo.Expression(pyo.Integers)

        # 6 nodes
        m.subexpr[1] = m.x[2] * pyo.exp(3*m.x[1])
        # 13 nodes in body: sum(linear(monomials), negation(expr))
        m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[1] == 1.5
        # 9 nodes
        m.eq[2] = 4*m.x[2] + m.x[3]**3 * m.subexpr[1] == 0.0

        n_nodes = count_nodes(m.eq[1].body, descend_into_named_expressions=False)
        assert n_nodes == 9

        # Note that the named subexpression still occupys a node even when
        # we descend into it.
        n_nodes = count_nodes(m.eq[1].body, descend_into_named_expressions=True)
        assert n_nodes == 15

        n_nodes = count_nodes(m.eq[2].body, descend_into_named_expressions=False)
        assert n_nodes == 9

        n_nodes = count_nodes(m.eq[2].body, descend_into_named_expressions=True)
        assert n_nodes == 15

    def test_count_nodes_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        m.eq = pyo.Constraint(pyo.Integers)
        m.subexpr = pyo.Expression(pyo.Integers)

        # 6 nodes
        m.subexpr[1] = m.x[2] * pyo.exp(3*m.x[1])
        # 13 nodes in body: sum(linear(monomials), negation(expr))
        m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[1] == 1.5
        # 9 nodes
        m.eq[2] = 4*m.x[2] + m.x[3]**3 * m.subexpr[1] == 0.0

        n_nodes = count_model_nodes(m, amplrepn=False)
        assert n_nodes == 24

    def test_count_nodes_nested_expr(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        m.eq = pyo.Constraint(pyo.Integers)
        m.subexpr = pyo.Expression(pyo.Integers)
        # 6 nodes
        m.subexpr[1] = m.x[2] * pyo.exp(3*m.x[1])
        # 9 nodes: product(power(linear(1, monomial), 2), expr)
        m.subexpr[2] = (1 - m.x[3])**2 * m.subexpr[1]
        # 13 nodes
        m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[2] == 1.5
        # 10 nodes: sum(monomial, product(pow, expr), expr)
        m.eq[2] = 4*m.x[2] + m.x[3]**3 * m.subexpr[2] + m.subexpr[1] == 0.0

        n_nodes = count_model_nodes(m, amplrepn=False)
        assert n_nodes == 34

    def test_count_nodes_of_named_expr(self):
        # If the root node of an expression *is* a named expression, we appear
        # to not be handling it properly
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        m.subexpr = pyo.Expression(expr=m.x[1] * m.x[2] + 2 * m.x[3])
        m.eq = pyo.Constraint(expr=m.subexpr == 3)

        visitor = NodeCounter(descend_into_named_expressions=False)
        # The root node of eq.body is a named expression. This should have
        # a node-count of 1 if we do not descend.
        assert m.eq.body.is_named_expression_type()
        count = visitor.walk_expression(m.eq.body)
        assert count == 1

    def test_model_nodecounter_replace_bounds(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3], bounds=(-1, 10))
        m.eq = pyo.Constraint(pyo.PositiveIntegers)
        m.eq[1] = m.x[1] + m.x[2] * m.x[3] == 2
        m.eq[2] = m.x[1] * m.x[2] * m.x[3] == 7
        m.eq[3] = m.x[1] + m.x[3] == 2
        var_elim = [m.x[1]]
        con_elim = [m.eq[1]]
        eliminate_variables(m, var_elim, con_elim, use_named_expressions=True)
        nnode = count_model_nodes(m, amplrepn=False)
        # Expresion for x: 6 nodes
        # eq[2]: 5 nodes
        # eq[3]: 3 nodes
        # x[1]_lb: 1
        # x[1]_ub: 1
        assert nnode == 16

    def test_defining_expression_with_subexpression(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3], bounds=(0, 10))
        m.eq = pyo.Constraint(pyo.PositiveIntegers)
        m.subexpr = pyo.Expression(pyo.PositiveIntegers)
        m.subexpr[1] = sum(
            (i + j/2) * m.x[1]**i * m.x[2]**j
            for i in range(1, 6) for j in range(1, 6)
        )
        m.eq[1] = 2*m.x[3] + m.subexpr[1] == 5
        expr = define_variable_from_constraint(m.x[3], m.eq[1])
        n_nodes_subexpr = count_nodes(m.subexpr[1].expr, descend_into_named_expressions=True)
        n_nodes_total = count_nodes(expr, descend_into_named_expressions=True)
        n_nodes_nodescend = count_nodes(expr, descend_into_named_expressions=False)

        # NOTE: Standard repn breaks apart named expression in order to separate linear
        # and nonlinear parts. This this way, we cannot exploit nested defined variables.
        # This is why we see large jumps in Pyomo node counts with matching-based
        # elimination sometimes.
        # Standard repn does not attempt to exploit a "subexpression fragment" that
        # represents the nonlinear portion, as AMPLRepn does.
        # We do not necessarily rely on this behavior (in fact, we would like nonlinear
        # subexpression fragments to be exploited...), but this test is here to remind
        # us that this behavior exists, as it significantly impacts our structural
        # results.
        assert n_nodes_total == n_nodes_nodescend
        assert n_nodes_total == 211


class TestAmplNodeCounter:

    def test_count_nodes_simple(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        expr = m.x[1] + m.x[2] + 2*m.x[3] - m.x[2] * pyo.exp(3*m.x[1])

        nodecount = count_amplrepn_nodes(expr)
        assert nodecount.linear == 3
        assert nodecount.nonlinear == 7

    def test_count_nodes_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        m.eq = pyo.Constraint(pyo.Integers)
        m.subexpr = pyo.Expression(pyo.Integers)

        # 6 nodes
        m.subexpr[1] = m.x[2] * pyo.exp(3*m.x[1])

        # 14 nodes in body (somehow)
        m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[1] == 1.5

        # 9 nodes
        m.eq[2] = 4*m.x[2] + m.x[3]**3 * m.subexpr[1] == 0.0

        n_nodes = count_model_amplrepn_nodes(m)
        assert n_nodes.linear == 4
        assert n_nodes.nonlinear == 6 + 2 + 5

        m.obj = pyo.Objective(expr=m.x[3]**2 + m.subexpr[1]**2)
        n_nodes = count_model_amplrepn_nodes(m)
        assert n_nodes.nonlinear == 6 + 2 + 5 + 7

    def test_count_nodes_nested_expr(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        m.eq = pyo.Constraint(pyo.Integers)
        m.subexpr = pyo.Expression(pyo.Integers)

        # 6 nodes
        m.subexpr[1] = m.x[2] * pyo.exp(3*m.x[1])

        # I count 7 nodes, all nonlinear
        # The (- x[3]) seems to get expanded into -1 * x[3]
        m.subexpr[2] = (1 - m.x[3])**2 * m.subexpr[1]

        # 3 linear terms + 2 nonlinear nodes
        m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[2] == 1.5

        # 1 linear term + 5 nonlinear nodes
        m.eq[2] = 4*m.x[2] + m.x[3]**3 * m.subexpr[2] == 0.0

        n_nodes = count_model_amplrepn_nodes(m)
        assert n_nodes.linear == 4
        assert n_nodes.nonlinear == 22

    def test_count_nodes_nlfragment(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        m.eq = pyo.Constraint(pyo.Integers)
        m.subexpr = pyo.Expression(pyo.Integers)
        # 2 linear terms + 8 nonlinear nodes
        m.subexpr[1] = m.x[3] + 2*m.x[1] - m.x[2] * pyo.exp(3*m.x[1])
        # To maximize the size of the linear subexpression, the subexpr[1] term
        # is split into linear(subexpr[1]) and nonlinear(subexpr[1]) (the latter
        # is the NLFragment). The linear terms are then inserted directly into
        # this expression (as linear subexpressions can't contain defined variables).
        # This increases node count, but decreases the work that must be done
        # to take derivatives.
        #
        # 1 linear term + 2 linear terms from subexpr[1] + 11 nonlinear nodes
        # sub[1] + (1 + -1*x[3])**2 * sub[1]
        m.subexpr[2] = m.subexpr[1] + 3*m.x[2] + (1 - m.x[3])**2 * m.subexpr[1]

        # Similarly, linear portions of common subexpressions are "lifted" into
        # the constraints.

        # 3 linear terms, 2 nonlinear nodes (Maybe 3 nl nodes?)
        m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[2] == 1.5

        # 7 nonlinear nodes, ... 2 linear terms from subexpr[1]?
        m.eq[2] = m.x[3]**3 * m.subexpr[2] + m.subexpr[1] == 0.0

        n_nodes = count_model_amplrepn_nodes(m)
        assert n_nodes.linear == 2 + 3
        # m.eq[2] gives 7 nonlinear node
        # eq[1] gives 3 nonlinear nodes
        # sub[1] fragment gives 8 nonlinear nodes
        # sub[2] fragment gives 11 nonlinear nodes
        # sub[1] expr (in sub[2]) gives 3 nodes
        # sub[2] expr (in eq[2]) gives 2 nodes
        # 34 = 7 + 3 + 8 + 11 * 3 + 2
        assert n_nodes.nonlinear == 34

    def test_fragment_used_but_not_expression(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        m.eq = pyo.Constraint(pyo.Integers)
        m.subexpr = pyo.Expression(pyo.Integers)
        m.subexpr[1] = m.x[1] + 2*m.x[2]**2 + m.x[1]*pyo.exp(m.x[3])
        m.eq[1] = 2*m.x[2] - m.subexpr[1] == 1
        m.eq[2] = m.x[3] + m.subexpr[1] == 2
        nnode = count_model_amplrepn_nodes(m)

        # The linear part of subexpr[1] gets inserted directly into eq[1] and
        # eq[2].
        assert nnode.linear == 4
        # The subexpresion has 10 nodes. The +subexpr[1] terms is one node,
        # while the -subexpr[1] term is two nodes.
        # This tests that we don't count the nodes in the full subexpression,
        # as this is not used anywhere.
        assert nnode.nonlinear == 13 # 3 + 10


if __name__ == "__main__":
    #pytest.main([__file__])
    TestNodeCounter().test_defining_expression_with_subexpression()
