#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:57:13 2023

@author:  sakshi
"""

from pyomo.environ import Constraint
from var_elim.distill import create_instance
from pyomo.core.expr.visitor import identify_variables
from pyomo.repn import generate_standard_repn

#Component data objects to get a list of all constraints

#In each constraint see if the form is coef*var - expr = 0
#Identify variables in order using general_repn and see if a linear var comes at [0] 

def identify_vars_for_elim_ampl(m):
    """
    Implements the ampl heuristic from the ampl book and returns a list 
    of variables to eliminate and the corresponding constraint using which the
    variable is eliminated

    """
    #Get constraint data in the order in which constraints are written
    cons = list(m.component_data_objects(Constraint))
    
    #Identify variables of the type ===> coef*v (+/-) expr = 0
    var_list = []
    con_list = []
    defining_var_ids = []
    for c in cons:
        #gets all vars in order from left to right in the constraint expression
        expr_vars = list(identify_variables(c.expr))
        
        #gets linear vars in the constraint expression
        repn = generate_standard_repn(c.body, compute_values=False, quadratic=False)
        linear_vars = repn.linear_vars
        
        #If the first variable in the constraint expression hasn't appeared in 
        #the rhs of a defining constraint and is linear, add it to var_list

        if id(expr_vars[0]) not in defining_var_ids:
            if expr_vars[0] is linear_vars[0]:
                var_list.append(expr_vars[0])
                con_list.append(c)
                for var in expr_vars:
                    defining_var_ids.append(id(var))
                    
                    #This will add all vars from the expression to the defining vars list
                    #This works because we anyways don't want to eliminate the same var twice 
                    #using 2 different constraints
        
    return var_list, con_list
        