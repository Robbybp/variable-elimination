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
Created on Tue May 23 10:39:34 2023

@author: sakshi
"""

import pyomo.environ as pyo
from distill import create_instance
from pyomo.core.expr.visitor import replace_expressions
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix
import matplotlib.pyplot as plt

def replacement_of_dx_vars(m, n, t):
    #A function to give the right substitution map from the diffeq constraints
    if n == 1:
        substitution_map = {id(m.dx[n, t]):1 / m.acond * m.V[t] * (
            m.y[n + 1, t] - m.x[n, t])
        }
        return substitution_map
    elif n in m.S_RECTIFICATION:
        substitution_map = {id(m.dx[n, t]):1 / m.atray * (
            m.L[t] * (m.x[n - 1, t] - m.x[n, t])
            - m.V[t] * (m.y[n, t] - m.y[n + 1, t])
        )
        }
        return substitution_map
    elif n == 17:
        substitution_map = {id(m.dx[n, t]):1 / m.atray * (
            m.Feed * m.x_Feed
            + m.L[t] * m.x[n - 1, t]
            - m.FL[t] * m.x[n, t]
            - m.V[t] * (m.y[n, t] - m.y[n + 1, t])
        )}
        return substitution_map
    elif n in m.S_STRIPPING:
        substitution_map = {id(m.dx[n, t]):1 / m.atray * (
            m.FL[t] * (m.x[n - 1, t] - m.x[n, t])
            - m.V[t] * (m.y[n, t] - m.y[n + 1, t])
        )}
        return substitution_map
    else:
        substitution_map = {id(m.dx[n, t]): 1 / m.areb * (
            m.FL[t] * m.x[n - 1, t]
            - (m.Feed - m.D) * m.x[n, t]
            - m.V[t] * m.y[n, t]
        )}
        return substitution_map


def var_elem():
    #Create distillation column instance and solve without var elem
    m = create_instance()
    ipopt = pyo.SolverFactory('ipopt')
    
    
    #Generate igraph
    igraph = IncidenceGraphInterface(m, include_inequality=False)
    cons = [con.name for con in igraph.get_adjacent_to(m.rr[1])]
    print(cons)
   
    for t in m.t:
        #Eliminate rr variable --> replace with u1
        substitution_map = {id(m.rr[t]): m.u1[t]}
        new_expr = replace_expressions(m.flowrate_rectification[t].expr, substitution_map)
        m.flowrate_rectification[t] = new_expr
        
        #Eliminate dx variable  with its RHS from the discretized equations dx_disc_eq
        #These equations look like dx[n, t] = x[n, t] - x[n, t-1]
        for n in m.S_TRAYS:
            if t != 1:
                substitution_map = replacement_of_dx_vars(m, n, t)
                new_expr = replace_expressions(m.dx_disc_eq[n,t].expr, substitution_map)
                m.dx_disc_eq[n,t] = new_expr
                
        #Replace L with V - D from  
        substitution_map = {id(m.L[t]): m.V[t] - m.D}
        new_expr = replace_expressions(m.flowrate_rectification[t].expr, substitution_map)
        m.flowrate_rectification[t] = new_expr
        
        new_expr = replace_expressions(m.flowrate_stripping[t].expr, substitution_map)
        m.flowrate_stripping[t] = new_expr
        
        for n in m.S_TRAYS:
            if t != 1:
                new_expr = replace_expressions(m.dx_disc_eq[n, t].expr, substitution_map)
                m.dx_disc_eq[n, t] = new_expr
            
        
   
    #Deactivate the reflux ratio constraint
    m.reflux_ratio.deactivate()
   
    #Deactivate the diffeq equation
    m.diffeq.deactivate()
       
    #Deactivate vapor column rule
    m.vapor_column.deactivate()
    
    #Solve model again with the eliminated variables
    ipopt.solve(m, tee= True)
   
    #A list of vars and cons for elimination
    var_list = [m.rr[2], m.dx[2,2], m.L[2]]
    con_list = [m.diffeq[2,2], m.vapor_column[2], m.reflux_ratio[2]]
    
    var_blocks, con_blocks = igraph.block_triangularize(var_list, con_list)
    for vb, cb in zip(var_blocks, con_blocks):
        assert len(vb) == 1
        assert len(cb) == 1 
    
    var_order = sum(var_blocks, [])
    con_order = sum(con_blocks, [])
    imat = get_structural_incidence_matrix(var_order, con_order)
    plt.spy(imat)
    
    print([var.name for var in var_order])
    print([con.name for con in con_order])
    return m, var_list, con_list
m, var_list, con_list = var_elem()
   
   
   

