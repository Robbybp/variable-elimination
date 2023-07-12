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

from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.repn import generate_standard_repn

def identify_vars_for_elim_min_degree(m, 
                                      major_elim = 'Variables', 
                                      eliminate_bounded_vars = False,
                                      eliminate_linear_cons_only = False,
                                      eliminate_linear_vars_only = False):
    if major_elim == 'Variables':
        var_list, con_list = var_major_elimination(m, 
                                                   eliminate_bounded_vars = eliminate_bounded_vars,
                                                   eliminate_linear_cons_only=eliminate_linear_cons_only, 
                                                   eliminate_linear_vars_only = eliminate_linear_vars_only)
    elif major_elim == 'Constraints':
        var_list, con_list = con_major_elimination(m, 
                                                   eliminate_bounded_vars = eliminate_bounded_vars,
                                                   eliminate_linear_cons_only=eliminate_linear_cons_only,
                                                   eliminate_linear_vars_only=eliminate_linear_vars_only)
    else:
        raise ValueError("major_elim must be 'Variables' or 'Constraints'")
        
    return var_list, con_list

def var_major_elimination(m,
                          eliminate_bounded_vars = False,
                          eliminate_linear_cons_only = False,
                          eliminate_linear_vars_only = False):
    """
    Identify variables for elimination and constraints to eliminate the variables
    using a minimum degree heuristic with first sorting on variable degree

    Parameters
    ----------
    m : Pyomo Model

    Returns
    -------
    var_list : List of variables to be eliminated
    con_list : List of constraints used to eliminate the variables

    """
    var_list = []
    con_list = []
    defining_var_ids = set()
    defining_con_ids = set()

    # A map which holds adjacent constraints to variables
    adj_cons_map = ComponentMap()

    # A map which holds adjacent variables to constraints
    adj_vars_map = ComponentMap()

    # Degree map for variables
    degree_map_var = ComponentMap()

    # Degree map for constraints
    degree_map_con = ComponentMap()

    # Generate incidence graph with active constraints
    igraph = IncidenceGraphInterface(m, active=True, include_inequality= True)

    # Generate linear incidence graph to identify variables appearing linearly
    # in the constraints
    # This is simply to narrow down the variables and constraints that are
    # candidates for replacement so we potentially have to do fewer checks
    # of standard_repn below.
    linear_igraph = IncidenceGraphInterface(m, active = True, linear_only=True, include_inequality = False)
    linear_vars = ComponentSet(linear_igraph.variables)
    linear_cons = ComponentSet(linear_igraph.constraints)

    # Get the degree of linear variables from the full graph
    #We should just look at vars and cons in the linear graph 
    #but get adjacency from the full graph
    for var in linear_vars:
        adj_cons = igraph.get_adjacent_to(var)
        adj_cons_map[var] = adj_cons
        degree_map_var[var] = len(adj_cons)

    # Get the degree of linear constraints from the full graph
    for con in linear_cons:
        adj_vars = igraph.get_adjacent_to(con)
        adj_vars_map[con] = adj_vars
        degree_map_con[con] = len(adj_vars)

    # Sort the degree map of the variables
    sorted_vars = ComponentMap(sorted(degree_map_var.items(), key=lambda item: item[1]))

    for var in sorted_vars:
        if not eliminate_bounded_vars and (var.lb is not None or var.ub is not None):
            pass
        elif eliminate_linear_vars_only and (len(adj_cons_map[var]) != len(linear_igraph.get_adjacent_to(var))):
            pass
        elif id(var) not in defining_var_ids and var in linear_vars:
            # This maps constraints that are valid for elimination to their degree.
            # It is used to both define the valid constraints store their degrees.
            degree_adj_cons = ComponentMap()
            for con in adj_cons_map[var]:
                if id(con) not in defining_con_ids and con in linear_cons:
                    # We still need a check that the var doesn't appear nonlinearly in the constraint
                    # Generating the repn here will call it much fewer times than anywhere else
                    repn = generate_standard_repn(
                        con.body, compute_values=False, quadratic=False
                    )
                    
                    if eliminate_linear_cons_only:
                        if len(repn.nonlinear_vars) == 0:
                            degree_adj_cons[con] = degree_map_con[con]
                    else:
                        if var not in ComponentSet(repn.nonlinear_vars):
                            degree_adj_cons[con] = degree_map_con[con]

                    
            # If variable appears linearly atleast in one constraint then find
            # the defining constraint
            if len(degree_adj_cons) != 0:
                defining_con = min(degree_adj_cons.items(), key=lambda item: item[1])[0]
                # Identify variables in the defining constraint to add to the list
                # of variables that cannot be eliminated
                expr_vars = list(identify_variables(defining_con.expr))
                for v in expr_vars:
                    defining_var_ids.add(id(v))

                defining_con_ids.add(id(defining_con))
                var_list.append(var)
                con_list.append(defining_con)
               
    return var_list, con_list

def con_major_elimination(m,
                          eliminate_bounded_vars = False,
                          eliminate_linear_cons_only = False,
                          eliminate_linear_vars_only = False):
    """
    Identify variables for elimination and constraints to eliminate the variables
    using a minimum degree heuristic with first sorting on constraint degree

    Parameters
    ----------
    m : Pyomo Model

    Returns
    -------
    var_list : List of variables to be eliminated
    con_list : List of constraints used to eliminate the variables

    """
    var_list = []
    con_list = []
    defining_var_ids = set()
    defining_con_ids = set()

    # A map which holds adjacent constraints to variables
    adj_cons_map = ComponentMap()

    # A map which holds adjacent variables to constraints
    adj_vars_map = ComponentMap()

    # Degree map for variables
    degree_map_var = ComponentMap()

    # Degree map for constraints
    degree_map_con = ComponentMap()

    # Generate incidence graph with active constraints
    igraph = IncidenceGraphInterface(m, active=True, include_inequality= True)
    
    # Generate linear incidence graph to identify variables appearing linearly
    # in the constraints
    linear_igraph = IncidenceGraphInterface(m, active = True, linear_only=True, include_inequality = False)
    linear_vars = ComponentSet(linear_igraph.variables)
    linear_cons = ComponentSet(linear_igraph.constraints)

    # Get the degree of linear variables from the full graph
    #We should just look at vars and cons in the linear graph 
    #but get adjacency from the full graph
    for var in linear_vars:
        adj_cons = igraph.get_adjacent_to(var)
        adj_cons_map[var] = adj_cons
        degree_map_var[var] = len(adj_cons)

    # Get the degree of linear constraints from the full graph
    for con in linear_cons:
        adj_vars = igraph.get_adjacent_to(con)
        adj_vars_map[con] = adj_vars
        degree_map_con[con] = len(adj_vars)
        
    # Sort the degree map of the constraints
    sorted_cons = ComponentMap(sorted(degree_map_con.items(), key=lambda item: item[1]))
    
    for con in sorted_cons:
        if id(con) not in defining_con_ids:
            #Generating repn here is the best since it'll generate it only once for each constraint
            repn = generate_standard_repn(
                con.body, compute_values=False, quadratic=False
            )
            if eliminate_linear_cons_only and len(repn.nonlinear_vars) != 0:
                pass
            else:
                degree_adj_vars = ComponentMap()
                for var in adj_vars_map[con]:
                    if not eliminate_bounded_vars and (var.lb is not None or var.ub is not None):
                        pass
                    elif eliminate_linear_vars_only and (len(adj_cons_map[var]) != len(linear_igraph.get_adjacent_to(var))):
                        pass
                    elif id(var) not in defining_var_ids and var in linear_vars:
                        if var not in ComponentSet(repn.nonlinear_vars):
                            degree_adj_vars[var] = degree_map_var[var]
            
                if len(degree_adj_vars) != 0:
                    defining_var = min(degree_adj_vars.items(), key=lambda item: item[1])[0]
                    
                    # Identify variables in the defining constraint to add to the list
                    # of variables that cannot be eliminated
                    expr_vars = list(identify_variables(con.expr))
                    for v in expr_vars:
                        defining_var_ids.add(id(v))
    
                    defining_con_ids.add(id(con))
                    var_list.append(defining_var)
                    con_list.append(con)
        
    return var_list, con_list
