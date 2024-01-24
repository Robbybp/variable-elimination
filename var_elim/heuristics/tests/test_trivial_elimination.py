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
import math
import pyomo.environ as pyo
from pyomo.common.collections import ComponentSet
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from var_elim.heuristics.trivial_elimination import (
    expr_filter,
    filter_constraints,
    get_trivial_constraint_elimination,
)
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables,
)


ipopt_avail = pyo.SolverFactory("ipopt").available()


def _make_model():
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3], initialize=1.0)
    m.y = pyo.Var(initialize = 0)
    m.y.fix(0)
    m.eq1 = pyo.Constraint(expr=m.x[1] == m.x[2])
    m.eq2 = pyo.Constraint(expr=m.x[1] == -m.x[2] + 1)
    m.eq3 = pyo.Constraint(expr=m.x[1] == 2*m.x[3])
    m.eq4 = pyo.Constraint(expr=m.x[1] + m.x[2] + m.x[3] == 0)
    m.eq5 = pyo.Constraint(expr=m.x[1] + 0.5*m.x[2] + m.x[3] == 0)
    m.eq6 = pyo.Constraint(expr=m.x[1] - m.x[1] + 2*m.x[3] == 0)
    m.eq7 = pyo.Constraint(expr=m.x[2] + m.x[2]**3 + m.x[1] == 0)
    m.eq8 = pyo.Constraint(expr=0*m.x[2]**3 + m.x[1] == 0)
    m.eq9 = pyo.Constraint(expr= m.y*m.x[1] == m.x[2]**2)
    return m


class TestExprFilter:
    def test_degree(self):
        m = _make_model()
        assert expr_filter(m.eq1.body, degree=2)
        assert not expr_filter(m.eq1.body, degree=1)
        assert expr_filter(m.eq4.body, degree=3)
        assert not expr_filter(m.eq4.body, degree=4)
        assert expr_filter(m.eq6.body, degree=1)

    def test_linear(self):
        # NOTE: linear in this context means "linear or affine"
        m = _make_model()
        assert expr_filter(m.eq1.body, linear=True)
        assert expr_filter(m.eq2.body, linear=True)
        assert not expr_filter(m.eq7.body, linear=True)
        assert expr_filter(m.eq8.body, linear=True)

    def test_affine(self):
        # Here affine means "affine and not linear"
        m = _make_model()
        assert not expr_filter(m.eq1.body, affine=True)
        assert expr_filter(m.eq2.body, affine=True)
        assert not expr_filter(m.eq7.body, affine=True)

    def test_linear_and_not_affine(self):
        m = _make_model()
        assert expr_filter(m.eq1.body, linear=True, affine=False)
        assert not expr_filter(m.eq2.body, linear=True, affine=False)

    def test_equal_coefficients(self):
        m = _make_model()
        assert expr_filter(m.eq2.body, equal_coefficients=True)
        assert expr_filter(m.eq4.body, equal_coefficients=True)
        assert not expr_filter(m.eq5.body, equal_coefficients=True)
        assert expr_filter(m.eq5.body, equal_coefficients=False)


class TestFilterConstraints:
    def test_get_trivial_equality(self):
        m = _make_model()
        trivial_eqs = filter_constraints(
            m,
            degree=2,
            linear=True,
            affine=False,
            equal_coefficients=True,
        )
        assert len(trivial_eqs) == 1
        assert trivial_eqs[0] is m.eq1

        trivial_affine_eqs = filter_constraints(
            m,
            degree=2,
            linear=True,
            affine=True,
            equal_coefficients=True,
        )
        assert len(trivial_affine_eqs) == 1
        assert trivial_affine_eqs[0] is m.eq2

        trivial_linear_or_affine = filter_constraints(
            m,
            degree=2,
            linear=True,
            equal_coefficients=True,
        )
        assert len(trivial_linear_or_affine) == 2
        assert (
            ComponentSet(trivial_linear_or_affine)
            == ComponentSet([m.eq1, m.eq2])
        )


class TestTrivialElimination:
    def _make_simple_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == m.y[1])
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3 * m.y[2] ** 3)
        m.eq3 = pyo.Constraint(expr=m.x[1] * m.x[2] == 1.0)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1] ** 2 + m.y[2] ** 2)

        return m

    def _make_model_for_nested_replacement(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2, 3, 4], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1, bounds=(0, 10))

        m.eq1 = pyo.Constraint(expr=m.x[4] == m.x[2] + m.x[1] - m.y[1]**2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == m.y[1]*m.y[2])
        m.eq3 = pyo.Constraint(expr=m.x[3] == m.x[1])
        m.eq4 = pyo.Constraint(expr=m.x[1] == m.y[2])
        # Note that eq5 cannot be used for elimination.
        m.eq5 = pyo.Constraint(expr=m.x[1] + m.x[2] + m.x[3] * m.x[4] == 2)

        # Need this initialization for convergence of the reduced-space problem
        calculate_variable_from_constraint(m.x[2], m.eq2)
        calculate_variable_from_constraint(m.x[1], m.eq4)
        calculate_variable_from_constraint(m.x[4], m.eq1)
        calculate_variable_from_constraint(m.x[3], m.eq3)

        m.obj = pyo.Objective(expr=m.x[1]**2 + 2*m.x[2]**2 + 3*m.x[3]**2 + 4*m.x[4]**2)

        return m

    def test_replacement_vars_simple(self):
        m = self._make_simple_model()
        vars_to_elim, cons_to_elim = get_trivial_constraint_elimination(m)
        # Make sure correct variables are eliminated
        assert len(vars_to_elim) == 1
        # Note that y[1] could also be eliminated. This test relies on the
        # order in which variables are encountered.
        assert m.x[1] is vars_to_elim[0]

    @pytest.mark.skipif(not ipopt_avail, reason="ipopt is not available")
    def test_same_solution_simple(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model()

        vars_to_elim, cons_to_elim = get_trivial_constraint_elimination(m2)

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)

    def test_nested_replacement(self):
        m = self._make_model_for_nested_replacement()
        vars_to_elim, cons_to_elim = get_trivial_constraint_elimination(m)

        # This is just the elimination we happen to arrive at given (a) the
        # maximum matching and (b) the tearing method. This could change if
        # either of these change.
        #
        # Since "trivial equality constraints" are symmetric, y[2] could
        # be eliminated instead of either of these. This is just what we
        # land on given how the expressions are written.
        assert len(vars_to_elim) == 2
        assert ComponentSet(vars_to_elim) == ComponentSet([m.x[1], m.x[3]])
        assert ComponentSet(cons_to_elim) == ComponentSet([m.eq3, m.eq4])

    @pytest.mark.skipif(not ipopt_avail, reason="ipopt is not available")
    def test_same_solution_nested_replacement(self):
        m1 = self._make_model_for_nested_replacement()
        solver = pyo.SolverFactory("ipopt")
        #res = solver.solve(m1, tee=True)
        #pyo.assert_optimal_termination(res)

        m2 = self._make_model_for_nested_replacement()
        vars_to_elim, cons_to_elim = get_trivial_constraint_elimination(m2)
        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)
        solver = pyo.SolverFactory("ipopt")

        # With default initialization, these solves converge to different local
        # solutions.
        res = solver.solve(m2, tee=True)
        pyo.assert_optimal_termination(res)

        # Initialize the full-space system so it converges to the same solution
        # that the reduced-space system arives at.
        m1.x[1] = m2.x[1]
        m1.x[2] = m2.x[2]
        m1.y[1] = m2.y[1]
        calculate_variable_from_constraint(m1.x[4], m1.eq1)
        calculate_variable_from_constraint(m1.x[3], m1.eq3)
        calculate_variable_from_constraint(m1.y[2], m1.eq4)

        solver.solve(m1, tee=True)
        pyo.assert_optimal_termination(res)

        # Now assert the solutions are the same.
        varnames = ["x[2]", "x[4]", "y[1]", "y[2]"]
        for name in varnames:
            v1 = m1.find_component(name)
            v2 = m2.find_component(name)
            assert math.isclose(v1.value, v2.value)


if __name__ == "__main__":
    pytest.main()
