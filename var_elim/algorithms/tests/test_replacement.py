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
from pyomo.core.expr.visitor import identify_variables
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables,
)
from var_elim.algorithms.validate import validate_solution


ipopt_avail = pyo.SolverFactory("ipopt").available()


class TestReplacementSimpleModel:
    def _make_simple_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2 * m.y[1] ** 2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3 * m.y[2] ** 3)
        m.eq3 = pyo.Constraint(expr=m.x[1] * m.x[2] == 1.0)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1] ** 2 + m.y[2] ** 2)

        return m

    def test_simple_replacement(self):
        m = self._make_simple_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m, var_order, con_order)

        new_igraph = IncidenceGraphInterface(m)
        # Make sure new model has the correct size
        assert len(new_igraph.constraints) == 1
        assert len(new_igraph.variables) == 2

        assert new_igraph.constraints[0] is m.eq3

        # Make sure proper replacement happened here
        assert ComponentSet(identify_variables(m.eq3.expr)) == ComponentSet(m.y[:])
        # Make sure no replacement happened here
        assert ComponentSet(identify_variables(m.obj)) == ComponentSet(m.y[:])

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model()

        vars_to_elim = [m2.x[1], m2.x[2]]
        cons_to_elim = [m2.eq1, m2.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)


class TestReplacementInObjective:
    def _make_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2 * m.y[1] ** 2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3 * m.y[2] ** 3)
        m.eq3 = pyo.Constraint(expr=m.x[1] * m.x[2] == 1.0)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.x[1] ** 2 + m.x[2] ** 2)

        return m

    def test_simple_replacement(self):
        m = self._make_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m, var_order, con_order)

        new_igraph = IncidenceGraphInterface(m)
        # Make sure new model has the correct size
        assert len(new_igraph.constraints) == 1
        assert len(new_igraph.variables) == 2

        assert new_igraph.constraints[0] is m.eq3

        # Make sure proper replacement happened here
        assert ComponentSet(identify_variables(m.eq3.expr)) == ComponentSet(m.y[:])

        # Make sure proper replacement happened in objective
        assert ComponentSet(identify_variables(m.obj)) == ComponentSet(m.y[:])

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution(self):
        m1 = self._make_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_model()

        vars_to_elim = [m2.x[1], m2.x[2]]
        cons_to_elim = [m2.eq1, m2.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        assert math.isclose(pyo.value(m1.obj), pyo.value(m2.obj), rel_tol=1e-6)


class TestReplacementWithBounds:
    def _make_simple_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1, bounds=(-5, 5))
        m.y = pyo.Var([1, 2], initialize=1)

        # Set bounds on eliminated variables, one of which will be active.
        m.x[2].setlb(2)
        m.x[2].setub(5)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2 * m.y[1] ** 2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3 * m.y[2] ** 3)
        m.eq3 = pyo.Constraint(expr=m.x[1] * m.x[2] == 1.0)

        m.y[1].setlb(0.1)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1] ** 2 + m.y[2] ** 2)

        return m

    def test_constraints_are_added(self):
        m = self._make_simple_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m, var_order, con_order)

        # Make sure new model has correct number of constraints
        new_igraph = IncidenceGraphInterface(m, include_inequality=True)

        assert len(new_igraph.constraints) == 5
        assert len(new_igraph.variables) == 2

        # I want to check that the correct expressions have the correct bounds

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        # Assert that x[2] is near its bound.
        assert math.isclose(m1.x[2].value, m1.x[2].lb)

        m2 = self._make_simple_model()

        vars_to_elim = [m2.x[1], m2.x[2]]
        cons_to_elim = [m2.eq1, m2.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        _, var_exprs, var_lb_map, var_ub_map = eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)

        x2_lb_con = var_lb_map[m2.x[2]]
        assert math.isclose(
            pyo.value(x2_lb_con.body),
            pyo.value(x2_lb_con.lower),
        )


class TestReplacementInInequalities:
    def _make_simple_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1, bounds=(-5, 5))
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2 * m.y[1] ** 2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3 * m.y[2] ** 3)
        m.eq3 = pyo.Constraint(expr=m.x[1] * m.x[2] == 1.0)
        m.eq4 = pyo.Constraint(expr=m.x[1] + m.x[2] >= 3)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1] ** 2 + m.y[2] ** 2)

        return m

    def test_simple_replacement(self):
        m = self._make_simple_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m, var_order, con_order)

        # Make sure new model has correct number of constraints
        new_igraph = IncidenceGraphInterface(m, include_inequality=True)

        assert len(new_igraph.constraints) == 6

        # Make sure proper replacement happened here
        assert ComponentSet(identify_variables(m.eq4.expr)) == ComponentSet(m.y[:])

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model()

        vars_to_elim = [m2.x[1], m2.x[2]]
        cons_to_elim = [m2.eq1, m2.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(pyo.value(m1.eq4.body), pyo.value(m1.eq4.lb), rel_tol=1e-6)
        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)


class TestReplaceWithNamedExpressions:
    def _make_simple_model(self, named_expressions=False):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1, bounds=(-5, 5))
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2 * m.y[1] ** 2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3 * m.y[2] ** 3)

        if named_expressions:
            m.subexpr = pyo.Expression(pyo.Integers)
            m.subexpr[1] = m.x[1] * m.x[2]
            m.subexpr[2] = m.x[1] + m.x[2]
            m.eq3 = pyo.Constraint(expr=m.subexpr[1] == 1.0)
            m.eq4 = pyo.Constraint(expr=m.subexpr[2] >= 3.0)
        else:
            m.eq3 = pyo.Constraint(expr=m.x[1] * m.x[2] == 1.0)
            m.eq4 = pyo.Constraint(expr=m.x[1] + m.x[2] >= 3)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1] ** 2 + m.y[2] ** 2)

        return m

    def test_defining_expressions(self):
        m = self._make_simple_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        m, var_exprs, var_lb_map, var_ub_map = eliminate_variables(
            m, var_order, con_order, use_named_expressions=True
        )

        assert isinstance(m.eliminated_variable_expressions, pyo.Expression)

        # Make sure the common subexpressions are being used in the constraint body
        assert m.eq3.body.args[0].ctype is pyo.Expression
        assert m.eq3.body.args[1].ctype is pyo.Expression
        assert m.eq4.body.args[0].ctype is pyo.Expression
        assert m.eq4.body.args[1].ctype is pyo.Expression

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_solve_and_validate(self):
        m = self._make_simple_model()
        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]
        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        m, var_exprs, var_lb_map, var_ub_map = eliminate_variables(
            m, var_order, con_order, use_named_expressions=True
        )

        solver = pyo.SolverFactory("ipopt")
        print(solver.executable())
        res = solver.solve(m, tee=True)
        pyo.assert_optimal_termination(res)

        valid, violations = validate_solution(
            m, var_exprs, cons_to_elim, tolerance=1e-6
        )
        violated_cons, violated_bounds, violated_elim_cons = violations
        for con in violated_cons:
            print(con.name)
        assert valid

    def test_replace_in_named_expressions(self):
        m = self._make_simple_model(named_expressions=True)

        # Sanity check that these are named expressions before
        # any elimination
        assert m.eq3.body.ctype is pyo.Expression
        assert m.eq4.body.ctype is pyo.Expression

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]
        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m, var_order, con_order)

        # Make sure named expressions are still used by constraints
        assert m.eq3.body.ctype is pyo.Expression
        assert m.eq4.body.ctype is pyo.Expression

        # Make sure replacement has happened in the named expressions
        assert (
            ComponentSet(identify_variables(m.subexpr[1].expr))
            == ComponentSet([m.y[1], m.y[2]])
        )
        assert (
            ComponentSet(identify_variables(m.subexpr[2].expr))
            == ComponentSet([m.y[1], m.y[2]])
        )


class TestExceptions:
    def test_discrete_variable(self):
        m = pyo.ConcreteModel()
        m.y = pyo.Var([1, 2], domain=pyo.Binary)
        m.eq = pyo.Constraint(expr=m.y[1] == 1 - m.y[2])
        msg = "Cannot eliminate discrete variable"
        with pytest.raises(RuntimeError, match=msg):
            eliminate_variables(m, [m.y[1]], [m.eq])


if __name__ == "__main__":
    pytest.main([__file__])
