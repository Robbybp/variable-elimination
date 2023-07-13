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
from var_elim.mip_formulations.mip_elim import identify_vars_for_elim_mip

ipopt_avail = pyo.SolverFactory("ipopt").available()
mip_solvers = [
    "gurobi",
    "cbc",
    "glpk",
]
solver_avail = False
for solver_name in mip_solvers:
    if pyo.SolverFactory(solver_name).available():
        solver_avail = pyo.SolverFactory(solver_name).available()
        solver_to_use = solver_name
        break


@pytest.mark.skipif(not solver_avail, reason="MIP solver is not available")
class TestMipFormulation1:
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

    def _make_complex_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)
        m.z = pyo.Var(initialize=1)

        m.eq1 = pyo.Constraint(expr=m.z == m.x[1] + m.x[2])
        m.eq2 = pyo.Constraint(expr=m.x[1] == -2 * m.y[1] ** 2 + m.x[2] ** 3)
        m.eq3 = pyo.Constraint(expr=m.x[2] == 1/2 * m.y[2] ** 3 + m.x[1])
        m.eq4 = pyo.Constraint(expr=m.x[1] * m.x[2] == 1.0)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)
        m.z.setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1] ** 2 + m.y[2] ** 2)

        return m

    def _make_simple_model2(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2 * m.y[1] ** 2)
        m.eq2 = pyo.Constraint(
            expr=m.x[2] == 2 * m.y[1] ** 2 + 3 * m.y[2] ** 3 - 4 * m.x[1]
        )
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
        m.eq3 = pyo.Constraint(expr=m.x[3] == m.x[1]*m.x[2] + m.x[4])
        m.eq4 = pyo.Constraint(expr=m.x[1] == m.x[2] - m.y[2])
        # Note that eq5 cannot be used for elimination.
        m.eq5 = pyo.Constraint(expr=m.x[1] + m.x[2] + m.x[3] * m.x[4] == 2)

        m.obj = pyo.Objective(expr=m.x[1]**2 + m.x[2]**2 + m.x[3]**2 + m.x[4]**2)

        return m

    def test_replacement_vars_simple(self):
        m = self._make_simple_model()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_mip(m, solver_name = solver_to_use)

        # Make sure correct variables are eliminated
        assert len(vars_to_elim) == 2
        assert m.x[1] is vars_to_elim[0]
        assert m.x[2] is vars_to_elim[1]

    def test_replacement_vars_complex(self):
        m = self._make_complex_model()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_mip(m, solver_name = solver_to_use)

        # Make sure correct variables are eliminated
        assert len(vars_to_elim) == 2
        assert m.z in ComponentSet(vars_to_elim)
        assert m.x[1] in ComponentSet(vars_to_elim) or m.x[2] in ComponentSet(vars_to_elim)
        

    def test_replacement_vars_simple2(self):
        m = self._make_simple_model2()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_mip(m, solver_name = solver_to_use)

        assert len(vars_to_elim) == 2
        assert m.x[1] is vars_to_elim[0]
        assert m.x[2] is vars_to_elim[1]

    @pytest.mark.skipif(not ipopt_avail, reason="ipopt is not available")
    def test_same_solution_simple(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_mip(m2, solver_name = solver_to_use)

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)

    @pytest.mark.skipif(not ipopt_avail, reason="ipopt is not available")
    def test_same_solution_simple2(self):
        m1 = self._make_simple_model2()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model2()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_mip(m2, solver_name = solver_to_use)

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        
    @pytest.mark.skipif(not ipopt_avail, reason="ipopt is not available")
    def test_same_solution_complex(self):
        m1 = self._make_complex_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=True)
        m1.display()
        pyo.assert_optimal_termination(res)

        m2 = self._make_complex_model()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_mip(m2, solver_name = solver_to_use)

        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        #Can't assert for x[1] or x[2] as either one of them can be eliminated using the mip but ot both 
        
    

    def test_nested_replacement(self):
        m = self._make_model_for_nested_replacement()
        vars_to_elim, cons_to_elim = identify_vars_for_elim_mip(m, solver_name = solver_to_use)

        assert ComponentSet(vars_to_elim) == ComponentSet(m.x[:])
        assert ComponentSet(cons_to_elim) == ComponentSet([m.eq1, m.eq2, m.eq3, m.eq4])

    @pytest.mark.skipif(not ipopt_avail, reason="ipopt is not available")
    def test_same_solution_nested_replacement(self):
        m1 = self._make_model_for_nested_replacement()
        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=True)
        pyo.assert_optimal_termination(res)

        m2 = self._make_model_for_nested_replacement()
        vars_to_elim, cons_to_elim = identify_vars_for_elim_mip(m2, solver_name = solver_to_use)
        var_order, con_order = define_elimination_order(vars_to_elim, cons_to_elim)
        eliminate_variables(m2, var_order, con_order)
        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m2, tee=True)
        pyo.assert_optimal_termination(res)

        for v1, v2 in zip(m1.y[:], m2.y[:]):
            assert math.isclose(v1.value, v2.value)


if __name__ == "__main__":
    pytest.main()
