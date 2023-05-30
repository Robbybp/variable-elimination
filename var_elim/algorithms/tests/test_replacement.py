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


ipopt_avail = pyo.SolverFactory("ipopt").available()


class TestReplacementSimpleModel:
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

    def test_simple_replacement(self):
        m = self._make_simple_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m, var_order, con_order)

        new_igraph = IncidenceGraphInterface(m)
        # Make sure new model has the correct size
        assert len(new_igraph.constraints) == 1
        assert len(new_igraph.variables) == 2

        assert new_igraph.constraints[0] is m.eq3

        # Make sure proper replacement happened here
        assert (
            ComponentSet(identify_variables(m.eq3.expr))
            == ComponentSet(m.y[:])
        )
        # Make sure no replacement happened here
        assert (
            ComponentSet(identify_variables(m.obj))
            == ComponentSet(m.y[:])
        )

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model()

        vars_to_elim = [m2.x[1], m2.x[2]]
        cons_to_elim = [m2.eq1, m2.eq2]

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
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

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2*m.y[1]**2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3*m.y[2]**3)
        m.eq3 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.x[1]**2 + m.x[2]**2)

        return m

    def test_simple_replacement(self):
        m = self._make_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m, var_order, con_order)

        new_igraph = IncidenceGraphInterface(m)
        # Make sure new model has the correct size
        assert len(new_igraph.constraints) == 1
        assert len(new_igraph.variables) == 2

        assert new_igraph.constraints[0] is m.eq3

        # Make sure proper replacement happened here
        assert (
            ComponentSet(identify_variables(m.eq3.expr))
            == ComponentSet(m.y[:])
        )
      
        # Make sure proper replacement happened in objective
        assert (
            ComponentSet(identify_variables(m.obj))
            == ComponentSet(m.y[:])
        )

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution(self):
        m1 = self._make_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_model()

        vars_to_elim = [m2.x[1], m2.x[2]]
        cons_to_elim = [m2.eq1, m2.eq2]

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        assert math.isclose(pyo.value(m1.obj), pyo.value(m2.obj))

class TestReplacementWithBounds:
    def _make_simple_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1, bounds = (-5, 5))
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2*m.y[1]**2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3*m.y[2]**3)
        m.eq3 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)

        return m
    
    def test_constraints_are_added(self):
        m = self._make_simple_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m, var_order, con_order)
        
        #Make sure new model has correct number of constraints
        new_igraph = IncidenceGraphInterface(m, include_inequality=True)
        
        assert len(new_igraph.constraints) == 5
        
        #I want to check that the correct expressions have the correct bounds 
        
    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution(self):
        m1 = self._make_simple_model()

        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)

        m2 = self._make_simple_model()

        vars_to_elim = [m2.x[1], m2.x[2]]
        cons_to_elim = [m2.eq1, m2.eq2]

        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m2, var_order, con_order)

        solver.solve(m2, tee=False)
        pyo.assert_optimal_termination(res)

        assert math.isclose(m1.y[1].value, m2.y[1].value)
        assert math.isclose(m1.y[2].value, m2.y[2].value)
        

class TestReplacementInInequalities:
    def _make_simple_model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var([1, 2], initialize=1, bounds = (-5, 5))
        m.y = pyo.Var([1, 2], initialize=1)

        m.eq1 = pyo.Constraint(expr=m.x[1] == 2*m.y[1]**2)
        m.eq2 = pyo.Constraint(expr=m.x[2] == 3*m.y[2]**3)
        m.eq3 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)
        m.eq4 = pyo.Constraint(expr = m.x[1]**2 + m.x[2]**2 >= 1)

        m.y[1].setlb(1.0)
        m.y[2].setlb(0.5)

        m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)

        return m
    
    def test_simple_replacement(self):
        m = self._make_simple_model()

        vars_to_elim = [m.x[1], m.x[2]]
        cons_to_elim = [m.eq1, m.eq2]
        
        var_order, con_order = define_elimination_order(
            vars_to_elim, cons_to_elim
        )
        eliminate_variables(m, var_order, con_order)
        
        #Make sure new model has correct number of constraints
        new_igraph = IncidenceGraphInterface(m, include_inequality=True)
        
        assert len(new_igraph.constraints) == 6
        
        # Make sure proper replacement happened here
        assert (
            ComponentSet(identify_variables(m.eq4.expr))
            == ComponentSet(m.y[:])
        )
      

    @pytest.mark.skipif(not ipopt_avail, reason="Ipopt is not available")
    def test_same_solution(self):
        m1 = self._make_simple_model()
   
        solver = pyo.SolverFactory("ipopt")
        res = solver.solve(m1, tee=False)
        pyo.assert_optimal_termination(res)
   
        m2 = self._make_simple_model()
   
        vars_to_elim = [m2.x[1], m2.x[2]]
        cons_to_elim = [m2.eq1, m2.eq2]
   
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
