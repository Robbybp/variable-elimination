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

import pyomo.environ as pyo
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix
from var_elim.models.distillation.distill import create_instance
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables
)

from var_elim.heuristics.ampl_heuristic import identify_vars_for_elim_ampl
import time
    
def main():
    m = create_instance(horizon=1500, nfe=400)
    
    #Using AMPL heuristic for elimination
    t0 = time.time()
    var_list, con_list = identify_vars_for_elim_ampl(m, randomize=True)
    t1 = time.time() - t0
    elim_vars = len(var_list)
    elim_cons = len(con_list)
    print("Time for ampl heuristic = ", t1)
    
    #Make incidece graph
    igraph = IncidenceGraphInterface(m, include_inequality = False)
    
    #Get ordered variable and corresponding constraints list
    t0 = time.time()
    var_order, con_order = define_elimination_order(var_list, con_list, igraph = igraph)
    t1 = time.time() - t0

    print("Time for getting elimination order = ", t1)    
    #Variable elimination
    t0 = time.time()
    m_reduced, _, _ = eliminate_variables(m, var_order, con_order, igraph = igraph)
    #new_var_bounds = fbbt(m_reduced)
    t1 = time.time() - t0
    
    print("Time to eliminate the variables = ", t1)
    ipopt = pyo.SolverFactory('ipopt')
    ipopt.options["print_timing_statistics"] = "yes"
    ipopt.solve(m_reduced, tee= True)
    return m_reduced, elim_vars, elim_cons
    
    
if __name__ == "__main__":
    variables_eliminated = []
    constraints_eliminated = []
    for i in range(0,10):
        m_reduced, elim_vars, elim_cons = main()
        variables_eliminated.append(elim_vars)
        constraints_eliminated.append(elim_cons)
 

