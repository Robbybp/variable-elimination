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
from pyomo.core.expr.visitor import replace_expressions, identify_variables
from var_elim.models.distillation.distill import create_instance
import itertools
import warnings

def get_components_from_model(m):
    #All sets are stored as indices of vars and constraints
    #Set C - all constraints with atleast 1 linear variable
    linear_cons = set()
    
    #Set X - all vars which appear linearly atleast once in the constraints
    #We don't want vars repeated when we append. Hence I chose sets here. 
    #Can also do list and append only if the var hasn't appeared before
    var_atleast_once_linear = set()
    
    #Set C_x - all constraints where x appears linearly
    cons_where_var_linear = {}
    
    #Set X_c - linear vars in constraint c
    linear_vars_in_cons = {}
    
    #Set X_C - all vars in constraint c which appear in X
    #aol stands for atleast once linear somewhere in the model
    all_vars_in_cons_aol = {}
    
   
    #Get the incidence graph from the model to find the number of vars and cons
    igraph = IncidenceGraphInterface(m, active = True, include_inequality= False)
    var_idx_map = ComponentMap((v, i) for i, v in enumerate(igraph.variables))
    con_idx_map = ComponentMap((c, i) for i, c in enumerate(igraph.constraints))
    
    for con in igraph.constraints:
        repn = generate_standard_repn(
            con.body, compute_values=False, quadratic=False
        )
        if len(repn.linear_vars) != 0:
            linear_cons.add(con_idx_map[con])
            linear_vars_in_cons[con_idx_map[con]] = [var_idx_map[var] for var in repn.linear_vars]
            for var in ComponentSet(repn.linear_vars):
                try:
                    var_ls = cons_where_var_linear[var_idx_map[var]]
                except:
                    var_ls = []
                var_atleast_once_linear.add(var_idx_map[var])
                
                var_ls.append(con_idx_map[con])
                cons_where_var_linear[var_idx_map[var]] = var_ls
    
    for con in igraph.constraints:
        all_vars_in_cons_aol[con_idx_map[con]] =[var_idx_map[var] for var in identify_variables(con.expr) if var_idx_map[var] in var_atleast_once_linear]
        
    #Too many returns maybe I should write this a a class(deal with it later)
    return linear_cons, var_atleast_once_linear, cons_where_var_linear, linear_vars_in_cons, all_vars_in_cons_aol, var_idx_map, con_idx_map

def get_var_con_pairs(m, var_idx_map, con_idx_map):
    #This function should probably go into a util folder since its general
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

def identify_vars_for_elim_mip2(model, solver_name = 'gurobi', tee = True):       
    """Identify defined variables and defining constraints using a mip 
    formulation (Russells)

    Parameters
    ----------
    m : Pyomo Model

    Returns
    -------
    var_list : List of variables to be eliminated
    con_list : List of constraints used to eliminate the variables

    """
    #Picking the solver
    solver_name_original = solver_name
    if not pyo.SolverFactory(solver_name).available():
        solver_name = None
        mip_solvers = ["gurobi", "cbc", "glpk"]
        for name in mip_solvers:
            if pyo.SolverFactory(name).available():
                solver_name = name
                break
    if solver_name is None: 
        raise RuntimeError(
            "MIP solver is not available"
        )
    if solver_name != solver_name_original:
        warnings.warn(
            f"{solver_name_original} not found using {solver_name} instead for the MIP solve"
            ) 
    
    var_list = []
    con_list = []
    linear_cons, var_atleast_once_linear, cons_where_var_linear, linear_vars_in_cons, all_vars_in_cons_aol, var_idx_map, con_idx_map = get_components_from_model(model)
    
    #All directed edges between variable nodes
    directed_edges_vars = list(itertools.product(var_idx_map.values(), var_idx_map.values()))
    
    m = pyo.ConcreteModel()
    
    #Sets
    m.C = pyo.Set(initialize = list(linear_cons))
    m.X = pyo.Set(initialize = list(var_atleast_once_linear))
    m.Cx = pyo.Set(cons_where_var_linear.keys(), initialize = cons_where_var_linear)
    m.Xc = pyo.Set(linear_vars_in_cons.keys(), initialize = linear_vars_in_cons)
    m.XC = pyo.Set(all_vars_in_cons_aol.keys(), initialize = all_vars_in_cons_aol)
    m.directed_edges = pyo.Set(initialize = directed_edges_vars)
    m.dummy_edges = pyo.Set(initialize = list(var_idx_map.values()))
    
    #Parameter
    m.X_size = pyo.Param(initialize=len(var_atleast_once_linear))
    
    #Variables
    m.y = pyo.Var(m.X, m.C, initialize = 1, domain = pyo.Binary)
    m.z = pyo.Var(m.directed_edges, domain = pyo.Binary)
    m.z_phi = pyo.Var(m.dummy_edges, domain = pyo.Binary)
    m.f = pyo.Var(m.directed_edges, domain = pyo.Binary)
    m.f_phi = pyo.Var(m.dummy_edges, domain = pyo.Binary)
    
    #Constraints
    #Atmost one constraint can be used to eliminate a variable
    def _match_var_to_con1(m, x):
        return sum(m.y[x, c] for c in m.Cx[x]) <= 1
    m.match_var_to_con1 = pyo.Constraint(m.X, rule = _match_var_to_con1)
    

    #A constraint eliminates atmost one variable
    def _match_var_to_con2(m, c):
        return sum(m.y[x, c] for x in m.Xc[c]) <= 1
    m.match_var_to_con2 = pyo.Constraint(m.C, rule = _match_var_to_con2)
    
    #We can write this constraint either on all linear vars or on all linear constraints
    #Here it is written on all linear vars first
    #Logic: for each var that appears linearly, look at each constraint in which it
    #appear linearly and look at the other variables that appear in that constraint.
    
    m.map_graph_variables_con = pyo.ConstraintList()
    
    for i in m.X:
        for c in m.Cx[i]:
            for j in m.XC[c]:
                m.map_graph_variables_con.add(m.z[i,j] >= m.y[i,c])
                
    #All constraints corresponding to the flow model
    def _flow_con1(m):
        return sum(m.z[i,j] for (i,j) in m.directed_edges) + sum(m.z_phi[i] for i in m.dummy_edges) == m.X_size
    m.flow_con1 = pyo.Constraint(rule = _flow_con1)
     
    def _flow_con2(m, i, j):
        return m.f[i,j] <= m.X_size*m.z[i,j]
    m.flow_con2 = pyo.Constraint(m.directed_edges, rule = _flow_con2)
    
    def _flow_con3(m, i):
        return m.f_phi[i] <= m.X_size*m.z_phi[i]
    m.flow_con3 = pyo.Constraint(m.dummy_edges, rule = _flow_con3)    
    
    def _flow_con4(m):
        return sum(m.f_phi[i] for i in m.dummy_edges) == m.X_size
    m.flow_con4 = pyo.Constraint(rule = _flow_con4)
    
    def _flow_con5(m, i):
        return m.f_phi[i] + sum(m.f[j, i] for j in m.dummy_edges) - sum(m.f[i, j] for j in m.dummy_edges) == 1
    m.flow_con5 = pyo.Constraint(m.dummy_edges, rule = _flow_con5)
    
    m.obj = pyo.Objective(expr = sum(sum(m.y[x, c] for c in m.Cx[x]) for x in m.X))
    
    solver = pyo.SolverFactory(solver_name)
    solver.solve(m, tee = tee, sense = 'Maximize')
    m.y.pprint()
    m.f_phi.pprint()
    m.f.pprint()
    m.z.pprint()
    import pdb;pdb.set_trace()
    #var_list, con_list = get_var_con_pairs(m, var_idx_map, con_idx_map)
    return var_list, con_list

    
m = pyo.ConcreteModel()
m.x = pyo.Var([1, 2], initialize=1)
m.y = pyo.Var([1, 2], initialize=1)

m.eq1 = pyo.Constraint(expr=m.x[1] + m.x[2] == 2*m.y[1]**2)
m.eq2 = pyo.Constraint(expr=m.x[2] == 3*m.y[2]**3)
m.eq3 = pyo.Constraint(expr=m.x[1]*m.x[2] == 1.0)

m.y[1].setlb(1.0)
m.y[2].setlb(0.5)

m.obj = pyo.Objective(expr=m.y[1]**2 + m.y[2]**2)
identify_vars_for_elim_mip2(m)
                
    