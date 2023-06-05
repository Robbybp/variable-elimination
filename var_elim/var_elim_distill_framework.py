#  ___________________________________________________________________________
#
#  Variable Elimination: Research code for variable elimination in NLPs
#
#  Copyright (c) 2023. Triad National Security, LLC. All rights reserved.
#
#  This program was produced under U.S. Government contract 89233218CNA000001
#  for Los Alamos National Laboratory (LANL), which i"s operated by Triad
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
from pyomo.common.collections import ComponentSet, ComponentMap
from var_elim.distill import create_instance
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables
)
#from var_elim.algorithms.ampl_heuristic import identify_vars_for_elim_ampl
                                        
import matplotlib.pyplot as plt

def main():
    m = create_instance()

    igraph = IncidenceGraphInterface(m)
    print("Before doing anything")
    print(f"N. variables: {len(igraph.variables)}")
    print(f"N. constraints: {len(igraph.constraints)}")

    print()
    print(f"Time set: {list(m.t)}")

    print()
    var = m.u1[1]
    print(f"Constraints adjacent to {var.name}:")
    for con in igraph.get_adjacent_to(var):
        print(f"  {con.name}: {con.expr.to_string()}")
    var = m.rr[1]
    print()
    print(f"Constraints adjacent to {var.name}:")
    for con in igraph.get_adjacent_to(var):
        print(f"  {con.name}: {con.expr.to_string()}")
    var = m.L[1]
    print()
    print(f"Constraints adjacent to {var.name}:")
    for con in igraph.get_adjacent_to(var):
        print(f"  {con.name}: {con.expr.to_string()}")
    var = m.V[1]
    print()
    print(f"Constraints adjacent to {var.name}:")
    for con in igraph.get_adjacent_to(var):
        print(f"  {con.name}: {con.expr.to_string()}")
    var = m.FL[1]
    print()
    print(f"Constraints adjacent to {var.name}:")
    for con in igraph.get_adjacent_to(var):
        print(f"  {con.name}: {con.expr.to_string()}")

    #Cuids of variables and constraints that can be used for elimination
    var_list= []
    con_list = []
    for t in m.t:
        #These set of variables when eliminated cause the degrees of freedom to decrease by 1
        var_list.append(m.rr[t])
        var_list.append(m.L[t])
        var_list.append(m.V[t])
        var_list.append(m.FL[t])
        con_list.append(m.reflux_ratio[t])
        con_list.append(m.flowrate_rectification[t])
        con_list.append(m.vapor_column[t])
        con_list.append(m.flowrate_stripping[t])

        # for n in m.S_TRAYS:
        #     var_list.append(m.y[n,t])
        #     con_list.append(m.mole_frac_balance[n,t])
        #     if t!= 1:
        #         var_list.append(m.dx[n,t])
        #         con_list.append(m.diffeq[n,t])

    print()
    print(f"N. variables to elim: {len(var_list)}")
    print(f"N. constraints to elim: {len(var_list)}")
    print(f"Predicted n. variables: {len(igraph.variables) - len(var_list)}")
    print(f"Predicted n. constraints: {len(igraph.constraints) - len(con_list)}")

    elim_var_set = ComponentSet(var_list)
    pred_var = [var for var in igraph.variables if var not in elim_var_set]

    print()
    var = m.rr[1]
    in_set = (var in elim_var_set)
    print(f"{var.name} eliminated: {in_set}")
    var = m.u1[1]
    in_set = (var in elim_var_set)
    print(f"{var.name} eliminated: {in_set}")
    var = m.L[1]
    in_set = (var in elim_var_set)
    print(f"{var.name} eliminated: {in_set}")

    #Creating the incidence graph
    #NOTE: We cannon deactivate the constraints before making the igraph
    igraph = IncidenceGraphInterface(m, include_inequality = False)

    #Get ordered variable and corresponding constraints list
    var_order, con_order = define_elimination_order(var_list, con_list, igraph = igraph)

    print()
    print(f"len(var_order) = {len(var_order)}")
    var_ord_map = ComponentMap([(var, i) for i, var in enumerate(var_order)])
    con_ord_map = ComponentMap([(con, i) for i, var in enumerate(con_order)])

    #print()
    #var = m.u1[1]
    #print(f"{var.name} eliminated at position {var_ord_map[var]} by constraint {con_order[var_ord_map[var]].name}")
    var = m.rr[1]
    print(f"{var.name} eliminated at position {var_ord_map[var]} by constraint {con_order[var_ord_map[var]].name}")
    var = m.L[1]
    print(f"{var.name} eliminated at position {var_ord_map[var]} by constraint {con_order[var_ord_map[var]].name}")

    #Sanity check: Plots the ordered matrix 
    imat = get_structural_incidence_matrix(var_order, con_order)
    mat = imat.tocsr()
    plt.spy(mat[0:20, 0:20])

    #Variable elimination
    m_reduced, var_expr_map, _, _ = eliminate_variables(m, var_order, con_order, igraph = igraph)

    var = m.L[1]
    print(f"{var.name} eliminated using expression {var_expr_map[var].to_string()}")

    igraph = IncidenceGraphInterface(m_reduced)
    print()
    print("After elimination")
    print(f"N. variables: {len(igraph.variables)}")
    print(f"N. constraints: {len(igraph.constraints)}")

    reduced_var_set = ComponentSet(igraph.variables)
    print()
    print("Not in the set of reduced variables:")
    for var in pred_var:
        if var not in reduced_var_set:
            print(f"  {var.name}")

    m_reduced.u1[:].fix()

    print()
    print("After elimination and fixing DOF")
    igraph = IncidenceGraphInterface(m_reduced)
    print(f"N. variables: {len(igraph.variables)}")
    print(f"N. constraints: {len(igraph.constraints)}")

    #ipopt = pyo.SolverFactory('ipopt')
    #ipopt.solve(m_reduced, tee= True)
    return m_reduced
    
    
if __name__ == "__main__":
    m_reduced = main()
