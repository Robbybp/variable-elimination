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
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.repn import generate_standard_repn
from var_elim.models.distillation.distill import create_instance

def get_components_from_model(m):
    """Return a map of variables to their indices, constraints to their indices, 
    set of edges and set of linear edges

    Parameters
    ----------
    m : Pyomo Model
    
    Returns
    -------
    var_idx_map: Map from variables to their indices
    con_idx_map: Map from constraints to their indices
    edge_set: List of tuples mapping constraints to variables
    linear_edge_Set: List of tuples mapping constraints to linear variables

    """
    #Get the incidence graph from the model to find the number of vars and cons
    igraph = IncidenceGraphInterface(m, active = True, include_inequality= False)
    var_idx_map = ComponentMap((v, i) for i, v in enumerate(igraph.variables))
    con_idx_map = ComponentMap((c, i) for i, c in enumerate(igraph.constraints))
    
    #Get the set of edges necessary
    #we can use linear igraph to get linear edges after Robbys PR is merged in 
    #pyomo main
    edge_set = []
    linear_edge_set = []
    for var in var_idx_map:
        adj_cons = igraph.get_adjacent_to(var)
        for con in adj_cons:
            edge_set.append((con_idx_map[con], var_idx_map[var]))
            repn = generate_standard_repn(
                con.body, compute_values=False, quadratic=False
            )
            if var in ComponentSet(repn.nonlinear_vars):
                pass
            else:
                linear_edge_set.append((con_idx_map[con], var_idx_map[var]))
    return var_idx_map, con_idx_map, edge_set, linear_edge_set

def get_var_con_pairs(m, var_idx_map, con_idx_map):
    """Identify defined variables and defining constraints from indices 
    from the mip solved model m

    Parameters
    ----------
    m : Pyomo Model
    var_idx_map: Map from variables to their indices
    con_idx_map: Map from constraints to their indices

    Returns
    -------
    var_list : List of variables to be eliminated
    con_list : List of constraints used to eliminate the variables

    """
    var_list = []
    con_list = []
    
    var_indices = []
    con_indices = []
    for i in range(0, len(con_idx_map)):
        for j in range(0, len(var_idx_map)):
            if pyo.value(m.y[i,j]) == 1:
                var_indices.append(j)
                con_indices.append(i)
                
    for var in var_idx_map:
        if var_idx_map[var] in var_indices:
            var_list.append(var)
            
    for con in con_idx_map:
        if con_idx_map[con] in con_indices:
            con_list.append(con)
    return var_list, con_list
    

def identify_vars_for_elim_mip(model):
    """Identify defined variables and defining constraints using a mip 
    formulation 

    Parameters
    ----------
    m : Pyomo Model

    Returns
    -------
    var_list : List of variables to be eliminated
    con_list : List of constraints used to eliminate the variables

    """
    
    var_idx_map, con_idx_map, edge_set, linear_edge_set = get_components_from_model(model)
    
    m = pyo.ConcreteModel()
    
    #Sets
    m.I = pyo.Set(initialize = range(0, len(con_idx_map) ))
    m.J = pyo.Set(initialize = range(0, len(var_idx_map)))
    
    m.E = pyo.Set(initialize = edge_set)
    m.LE = pyo.Set(initialize = linear_edge_set)
    
    if len(m.J) >= len(m.I):
        m.K = pyo.Set(initialize = range(0, len(m.J)))
    else:
        m.K = pyo.Set(initialize = range(0, len(m.I)))

    #Parameters
    #P_k needs to go from 1 to K 
    def _param_definition(m, k):
        return k + 1
    m.P = pyo.Param(m.K, rule = _param_definition)
    m.big_M = pyo.Param(initialize = len(m.K) + 1)
    
    #Variables
    m.Z = pyo.Var(domain = pyo.Integers)
    m.y = pyo.Var(m.I, m.J, domain = pyo.Binary)
    m.q = pyo.Var(m.I, m.I, domain = pyo.Binary)
    m.p = pyo.Var(m.J, m.J, domain = pyo.Binary)
    m.z = pyo.Var(m.K, domain = pyo.Binary)

    #Constraints
    #Define the objective function
    def _obj_linking_con(m):
        return m.Z == sum(sum(m.y[i, j] for i in m.I) for j in m.J)
    m.obj_linking_con = pyo.Constraint(rule = _obj_linking_con)
    
    #Each variable can be matched only to one constraint
    def _match_var_con1(m, j):
        return sum(m.y[i, j] for i in m.I) <= 1
    m.match_var_con1 = pyo.Constraint(m.J, rule = _match_var_con1)
    
    #Each constraint can be matched only to one variable
    def _match_var_con2(m, i):
        return sum(m.y[i, j] for j in m.J) <= 1
    m.match_var_con2 = pyo.Constraint(m.I, rule = _match_var_con2)
    
    #y is 0 if edge is not in th set of linear edges
    def _non_linear_edge_con(m, i, j):
        if (i,j) not in m.LE:
            return m.y[i, j] == 0
        else:
            return pyo.Constraint.Skip
    m.non_linear_edge_con = pyo.Constraint(m.I, m.J, rule = _non_linear_edge_con)
   
    #Variable reordering constraints
    def _var_reordering_con1(m, k):
        return m.Z+1 <= m.big_M*m.z[k] + m.P[k]
    m.var_reordering_con1 = pyo.Constraint(m.K, rule = _var_reordering_con1)
    
    def _var_reordering_con2(m, k):
        return m.P[k] <= m.big_M*(1-m.z[k]) + m.Z
    m.var_reordering_con2 = pyo.Constraint(m.K, rule =  _var_reordering_con2)
    
    #Assignment constraints
    def _assignment_con1(m, ip):
        return sum(m.q[i , ip] for i in m.I) == m.z[ip]
    m.assignment_con1 = pyo.Constraint(m.I, rule = _assignment_con1)
    
    def _assignment_con2(m, i):
        return sum(m.q[i, ip] for ip in m.I) == sum(m.y[i, j] for j in m.J)
    m.assignment_con2 = pyo.Constraint(m.I, rule = _assignment_con2)
    
    def _assignment_con3(m, jp):
        return sum(m.p[j, jp] for j in m.J) == m.z[jp]
    m.assignment_con3 = pyo.Constraint(m.J, rule = _assignment_con3)
    
    def _assignment_con4(m, j):
        return sum(m.p[j, jp] for jp in m.J) == sum(m.y[i, j] for i in m.I)
    m.assignment_con4 = pyo.Constraint(m.J, rule = _assignment_con4)
    
    #Lower traingular constraint
    def _lower_triangular_con(m, i, j):
        return sum(m.P[jp]*m.p[j, jp] for jp in m.J) <= sum(m.P[ip]*m.q[i, ip] for ip in m.I) + m.big_M*(2- sum(m.y[k,j]for k in m.I) - sum(m.y[i, l] for l in m.J))
    m.lower_triangular_con = pyo.Constraint(m.E, rule = _lower_triangular_con)
    
    #Objective function
    m.obj = pyo.Objective(expr = -m.Z, sense = pyo.minimize)
    solver = pyo.SolverFactory('glpk')
    solver.solve(m, tee = True)
    
    var_list, con_list = get_var_con_pairs(m, var_idx_map, con_idx_map)
    return var_list, con_list