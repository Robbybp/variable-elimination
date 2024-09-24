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
from var_elim.heuristics.ampl_heuristic import identify_vars_for_elim_ampl

ipopt_avail = pyo.SolverFactory("ipopt").available()


class TestAmplHeuristic:
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
    
    def _make_complex_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)
        m.z = pyo.Var(initialize= 1)

        m.eq1 = pyo.Constraint(expr=m.z == m.x[1] + m.x[2])
        m.eq2 = pyo.Constraint(expr=m.x[1] == 2*m.y[1]**2 + m.x[2]**3 - 4*m.x[1])
        m.eq3 = pyo.Constraint(expr=m.x[2] == 3*m.y[2]**3)
        m.eq4 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)
        

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)
        m.z.setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)

        return m
    
    def _make_simple_model2(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2*m.y[1]**2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 2*m.y[1]**2 + 3*m.y[2]**3 - 4*m.x[1])
        m.eq3 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)
        
        return m
    
    def _make_model(self):
        # This model is to test that the variable will be eliminated even if it 
        # is no the first variable in the expression
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1)
        m.y = pyo.Var([1, 2], initialize=1)
        
        m.eq1 = pyo.Constraint(expr=2*m.y[1]**2 == m.x[1])
        m.eq2 = pyo.Constraint(expr=2*m.y[1]**2 + 3*m.y[2]**3 - 4*m.x[1]- m.x[2] == 0)
        m.eq3 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)
        
        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)
        
        m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)
        
        return m

    def test_replacement_vars_simple(self):
        m = self._make_simple_model()

        vars_to_elim, cons_to_elim = identify_vars_for_elim_ampl(m)
        
        #Make sure correct variables are eliminated
        assert len(vars_to_elim) == 2
        assert m.x[1] is vars_to_elim[0]
        assert m.x[2] is vars_to_elim[1]
        
    def test_replacement_vars_complex(self):
        m = self._make_complex_model()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_ampl(m)
        
        #Make sure correct variables are eliminated
        assert len(vars_to_elim) == 1
        assert m.x[1] is vars_to_elim[0]
        
    def test_replacement_vars_simple2(self):
        m = self._make_simple_model2()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_ampl(m)
        
        assert len(vars_to_elim) == 2
        assert m.x[1] is vars_to_elim[0]
        assert m.x[2] is vars_to_elim[1]
        
    def test_replacement_vars_unordered(self):
        m = self._make_model()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_ampl(m)
        
        assert len(vars_to_elim) == 2
        assert m.x[1] is vars_to_elim[0]
        assert m.x[2] is vars_to_elim[1]
        

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution_simple(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_ampl(m2)

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        
    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution_complex(self):
        m1 = self._make_complex_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_complex_model()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_ampl(m2)

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        assert math.isclose(m1.x[2].value, m2.x[2].value)
        assert math.isclose(m1.z.value, m2.z.value)

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution_simple2(self):
        m1 = self._make_simple_model2()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model2()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_ampl(m2)

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        
    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution_unordered(self):
        m1 = self._make_model()
        
        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)
        
        m2 = self._make_model()
        
        vars_to_elim, cons_to_elim = identify_vars_for_elim_ampl(m2)

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        
 
if __name__ == "__main__":
    pytest.main()