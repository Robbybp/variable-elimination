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
#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

"""
Created on Wed May 24 09:28:28 2023

@author:  sakshi
"""

import pyomo.environ as pyo
from pyomo.core.expr.visitor import replace_expressions
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix

from distill import create_instance
from algorithms.replace import define_variable_from_constraint

import matplotlib.pyplot as plt

import copy


def choose_var_con_pair(m, cuid_var, cuid_con):
    """
    Pick some logical variables and constraints which can be eliminated
    """
    var_list = []
    con_list = []
    
    for cuid in cuid_var:
        var_list.append(m.find_component(cuid))
    
    for cuid in cuid_con:
        con_list.append(m.find_component(cuid))
        
    assert len(var_list) == len(con_list)
    
    #Possible issues we could run into in terms of error throwing:::
    #The replace function asserts if the variable participates linearly in the constraint or not 
    #But we do block triangularization first to find the elimination order 
    #and that won't work if there is no perfect matching
    return var_list, con_list


def define_elimination_order(igraph, var_list, con_list):
    """
    Finds elimination order using block triangularize from incidence graph interface
    """
    
    var_blocks, con_blocks = igraph.block_triangularize(var_list, con_list)
    
    for vb, cb in zip(var_blocks, con_blocks):
        assert len(vb) == 1
        assert len(cb) == 1 
    
    var_order = sum(var_blocks, [])
    con_order = sum(con_blocks, [])
    return var_order, con_order


def var_elimination_routine(m, igraph, var_order, con_order):
    """
    Does the actual elimination by defining expression from constraint, defines 
    susbtitution map and replaces the variable in every adjacent constraint
    """
    
    for var, con in zip(var_order, con_order):
        #Get expression for the variable from constraint
        var_expr = define_variable_from_constraint(var, con)
        
        #Build substitution map
        substitution_map = {id(var): var_expr}
        
        #Get constraints in which the variable appears
        #This will have the deactivated constraints too
        
        adj_cons = igraph.get_adjacent_to(var)
        for ad_con in adj_cons:
            if ad_con.name != con.name: 
                new_expr = replace_expressions(ad_con.expr, substitution_map)
                ad_con.set_value(new_expr)
    return m


def main():
    m = create_instance()
    
    #Cuids of variables and constraints that can be used for elimination
    cuid_var = []
    cuid_con = []
    for t in m.t:
        cuid_var.append('rr[{}]'.format(t))
        cuid_var.append('L[{}]'.format(t))
        cuid_var.append('FL[{}]'.format(t))
        cuid_con.append('reflux_ratio[{}]'.format(t))
        cuid_con.append('vapor_column[{}]'.format(t))
        cuid_con.append('flowrate_stripping[{}]'.format(t))
        
        if t!= 1:
            for n in m.S_TRAYS:
                cuid_var.append('dx[{}, {}]'.format(n,t))
                cuid_con.append('diffeq[{},{}]'.format(n, t))
            
    #Actual variables and constraints from the model
    var_list, con_list = choose_var_con_pair(m, cuid_var, cuid_con)
    
    #Creating the incidence graph
    #NOTE: We cannon deactivate the constraints before making the igraph
    igraph = IncidenceGraphInterface(m, include_inequality = False)
    
    #Deactivating constraints which are used for variable elimination
    for con in con_list:
        con.deactivate()
    
    #Get ordered variable and corresponding constraints list
    var_order, con_order = define_elimination_order(igraph, var_list, con_list)
    
    #Sanity check: Plots the ordered matrix 
    imat = get_structural_incidence_matrix(var_order, con_order)
    plt.spy(imat)
    
    #Variable elimination
    m_reduced = var_elimination_routine(m, igraph, var_order, con_order)
    ipopt = pyo.SolverFactory('ipopt')
    ipopt.solve(m_reduced, tee= True)
    
if __name__ == "__main__":
    main()
