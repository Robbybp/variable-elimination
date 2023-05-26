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

"""
Created on Wed May 24 09:28:28 2023

@author:  sakshi
"""

import pyomo.environ as pyo
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix
from var_elim.distill import create_instance
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables
)
                                         
import matplotlib.pyplot as plt

def main():
    m = create_instance()
    
    #Cuids of variables and constraints that can be used for elimination
    var_list= []
    con_list = []
    for t in m.t:
        var_list.append(m.rr[t])
        var_list.append(m.L[t])
        var_list.append(m.FL[t])
        con_list.append(m.reflux_ratio[t])
        con_list.append(m.vapor_column[t])
        con_list.append(m.flowrate_stripping[t])
        
        if t!= 1:
            for n in m.S_TRAYS:
                var_list.append(m.dx[n, t])
                con_list.append(m.diffeq[n,t])
    
    #Creating the incidence graph
    #NOTE: We cannon deactivate the constraints before making the igraph
    igraph = IncidenceGraphInterface(m, include_inequality = False)
    
    #Get ordered variable and corresponding constraints list
    var_order, con_order = define_elimination_order(var_list, con_list, igraph = igraph)
    
    #Sanity check: Plots the ordered matrix 
    imat = get_structural_incidence_matrix(var_order, con_order)
    plt.spy(imat)
    
    #Variable elimination
    m_reduced = eliminate_variables(m, var_order, con_order, igraph = igraph)
    ipopt = pyo.SolverFactory('ipopt')
    ipopt.solve(m_reduced, tee= True)
    
if __name__ == "__main__":
    main()
