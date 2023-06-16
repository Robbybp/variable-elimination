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
from var_elim.heuristics.min_degree_heuristic import identify_vars_for_elim_min_degree
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables,
)

ipopt_avail = pyo.SolverFactory("ipopt").available()

class TestMinDegreeHeuristic:
    def _make_simple_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2*m.y[1]**2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3*m.y[2]**3)
        m.eq3 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)

        return m
    
    def _make_simple_model2(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == m.x[1]**2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3*m.y[2]**3 + m.y[1]**2)
        m.eq3 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)

        m.y[1].setlb(0.2)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)

        return m
        
    
    def _make_simple_model3(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2*m.y[1]**2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 2*m.y[1]**2 + 3*m.y[2]**3 - 4*m.x[1])
        m.eq3 = pyo.Constraint(expr=m.x[2] == m.y[2]**2)
        m.eq4 = pyo.Constraint(expr = m.x[2]*m.y[1] == 2)
        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)
        
        return m
    
    def test_replacement_vars_simple3(self):
        m = self._make_simple_model3()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_min_degree(m)
        
        assert len(vars_to_elim) == 2
        assert len(cons_to_elim) == 2
        assert m.x[1] in ComponentSet(vars_to_elim)
        assert m.x[2] in ComponentSet(vars_to_elim)
        assert m.eq1 in ComponentSet(cons_to_elim)
        assert m.eq3 in ComponentSet(cons_to_elim)
    
    
    def test_replacement_vars_simple(self):
        m = self._make_simple_model()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_min_degree(m)
        
        #Make sure correct variables and constraints are eliminated
        assert len(vars_to_elim) == 2
        assert len(cons_to_elim) == 2
        assert m.x[1] in ComponentSet(vars_to_elim)
        assert m.x[2] in ComponentSet(vars_to_elim)
        assert m.eq1 in ComponentSet(cons_to_elim)
        assert m.eq2 in ComponentSet(cons_to_elim)
    
    
    def test_replacement_vars_simple2(self):
        m = self._make_simple_model2()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_min_degree(m)
        
        assert len(vars_to_elim) == 1
        assert len(cons_to_elim) == 1
        assert m.x[2] in ComponentSet(vars_to_elim)        
        
        
    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution_simple(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_min_degree(m2)

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        
        
    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution_simple2(self):
        m1 = self._make_simple_model2()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model2()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_min_degree(m2)

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        assert math.isclose(m1.x[1].value, m2.x[1].value)
        
if __name__ == "__main__":
    pytest.main()
