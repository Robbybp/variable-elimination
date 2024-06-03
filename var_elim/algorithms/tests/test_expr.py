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
)
from var_elim.models.distillation.distill import create_instance as create_distill_instance
from var_elim.models.opf.opf_model import make_model as make_opf_model


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

class TestAmplNodeCounter:

    def test_count_nodes_simple(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3])
        expr = m.x[1] + m.x[2] + 2*m.x[3] - m.x[2] * pyo.exp(3*m.x[1])

        nodecount = count_amplrepn_nodes(expr)
        assert nodecount.linear == 3
        assert nodecount.nonlinear == 5
        #assert n_nodes == 17

    #def test_count_nodes_model(self):
    #    m = pyo.ConcreteModel()
    #    m.x = pyo.Var([1, 2, 3])
    #    m.eq = pyo.Constraint(pyo.Integers)
    #    m.subexpr = pyo.Expression(pyo.Integers)

    #    # 6 nodes (somehow)
    #    m.subexpr[1] = m.x[2] * pyo.exp(3*m.x[1])

    #    # 14 nodes in body (somehow)
    #    m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[1] == 1.5

    #    # 9 nodes
    #    m.eq[2] = 4*m.x[2] + m.x[3]**3 * m.subexpr[1] == 0.0

    #    n_nodes = count_model_nodes(m, amplrepn=True)
    #    assert n_nodes == 24

    #    n_linear_nodes = count_model_nodes(m, amplrepn=True, linear_only=True)
    #    assert n_linear_nodes == 14

    #    m.obj = pyo.Objective(expr=m.x[3]**2 + m.subexpr[1]**2)
    #    n_nodes = count_model_nodes(m, amplrepn=True)
    #    # I count 7 nodes in the objective, but somehow the nl repn contains 8 nodes
    #    assert n_nodes == 29

    #def test_count_nodes_nested_expr(self):
    #    m = pyo.ConcreteModel()
    #    m.x = pyo.Var([1, 2, 3])
    #    m.eq = pyo.Constraint(pyo.Integers)
    #    m.subexpr = pyo.Expression(pyo.Integers)

    #    # 6 nodes
    #    m.subexpr[1] = m.x[2] * pyo.exp(3*m.x[1])

    #    # I count 7 nodes, all nonlinear
    #    m.subexpr[2] = (1 - m.x[3])**2 * m.subexpr[1]

    #    # 14 nodes in body (somehow)
    #    m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[2] == 1.5

    #    # 9 nodes
    #    m.eq[2] = 4*m.x[2] + m.x[3]**3 * m.subexpr[2] == 0.0

    #    n_nodes = count_model_nodes(m, amplrepn=True)
    #    assert n_nodes == 31

    #    n_linear_nodes = count_model_nodes(m, amplrepn=True, linear_only=True)
    #    assert n_linear_nodes == 14

    #def test_count_nodes_nlfragment(self):
    #    m = pyo.ConcreteModel()
    #    m.x = pyo.Var([1, 2, 3])
    #    m.eq = pyo.Constraint(pyo.Integers)
    #    m.subexpr = pyo.Expression(pyo.Integers)
    #    # 15 nodes
    #    m.subexpr[1] = m.x[3] + 2*m.x[1] - m.x[2] * pyo.exp(3*m.x[1])
    #    # To maximize the size of the linear subexpression, the subexpr[1] term
    #    # is split into linear(subexpr[1]) and nonlinear(subexpr[1]) (the latter
    #    # is the NLFragment). The linear terms are then inserted directly into
    #    # this expression (as linear subexpressions can't contain defined variables).
    #    # This increases node count, but decreases the work that must be done
    #    # to take derivatives.
    #    #
    #    # Linear: 11 nodes, nonlinear: 11 nodes (and the plus to combine them)
    #    # Then we add 1 to account for the indirection required to access the
    #    # NLFragment.
    #    # Total: 24 nodes
    #    m.subexpr[2] = m.subexpr[1] + 3*m.x[2] + (1 - m.x[3])**2 * m.subexpr[1]
    #    # Similarly, linear portions of common subexpressions are "lifted" into
    #    # the constraints.
    #    # Linear: 11 nodes. Nonlinear: 2 nodes.
    #    # Total: 15 nodes
    #    m.eq[1] = m.x[1] + m.x[2] + 2*m.x[3] - m.subexpr[2] == 1.5
    #    # Linear: 7 nodes. Nonlinear: 7 nodes.
    #    # Total: 16 nodes
    #    m.eq[2] = m.x[3]**3 * m.subexpr[2] + m.subexpr[1] == 0.0

    #    n_nodes = count_model_nodes(m, amplrepn=True)
    #    assert n_nodes == 61


if __name__ == "__main__":
    pytest.main([__file__])
