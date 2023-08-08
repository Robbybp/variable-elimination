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
from pyomo.core.base.componentuid import ComponentUID
from var_elim.models.distillation.distill import create_instance
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables
)
from pyomo.environ import Constraint
from pyomo.core.expr.current import EqualityExpression
from var_elim.heuristics.min_degree_heuristic import identify_vars_for_elim_min_degree
from var_elim.heuristics.ampl_heuristic import identify_vars_for_elim_ampl
import time
from var_elim.util.ipopt_helper import ipopt_solver_timing_records
import json
import sys
import random

path = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/distillation/'
path_ampl = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/distillation/randomized_ampl/'
path_vm = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/distillation/randomized_var_major/'
path_cm = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/distillation/randomized_con_major/'
for i in range(0,5):
    timing_data = {}
    #nfe = 2000, ncp =5
    m_ampl = create_instance(horizon = 1500, nfe = 400)
    
    #Get constraints and shuffle them
    #Get active equality constraints from the model
    cons = []
    
    for c in m_ampl.component_data_objects(Constraint, active=True):
        if isinstance(c.expr, EqualityExpression):
            cons.append(c)
    
    ids = list(range(0, len(cons)))
    random.shuffle(ids)
    
    with open(path + "constraint_ordering_" + str(i) + ".txt", 'w') as f:
        f.write(json.dumps(ids))
    
    #Ampl elimination
    t0 = time.time()
    var_list, con_list, cons_ampl = identify_vars_for_elim_ampl(m_ampl, eliminate_bounded_vars=False, eliminate_linear_cons_only=False, constraint_ordering= ids)
    t1 = time.time() - t0
    
    igraph = IncidenceGraphInterface(m_ampl, include_inequality = False)
    var_order, con_order = define_elimination_order(var_list, con_list, igraph = igraph)
    
    igraph = IncidenceGraphInterface(m_ampl, include_inequality = True)
    m_reduced, _, _ = eliminate_variables(m_ampl, var_order, con_order, igraph = igraph)
    
    nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m_reduced, tee=True)
    
    timing_data[str(i) + ', num_vars'] = nv
    timing_data[str(i) + ', num_cons'] = nc
    timing_data[str(i) + ', num_ineq'] = ni
    timing_data[str(i) + ', ipopt_linsolve'] = linsolvetime
    timing_data[str(i) + ', nlp_eval_time'] = nlptime
    timing_data[str(i) + ', heuristic_time'] = t1
    timing_data[str(i) + ', NZ_eq_jac'] = nze
    timing_data[str(i) + ', NZ_ineq_jac'] = nzie
    timing_data[str(i) + ', NZ_lag_hess'] = nzlh
    timing_data[str(i) + ', iter'] = niter
    
    with open(path_ampl + "timing_data_" + str(i) + ".txt", 'w') as f:
        f.write(json.dumps(timing_data))
    
    #Var major elimination
    #Make model again to make sure it is not initialized to the previous solution
    timing_data = {}
    m_var = create_instance(horizon = 1500, nfe = 400)
    
    t0 = time.time()
    var_list, con_list, cons_var = identify_vars_for_elim_min_degree(m_var, major_elim = "Variables", eliminate_bounded_vars=False, eliminate_linear_cons_only=False, constraint_ordering= ids)
    t1 = time.time() - t0
    
    igraph = IncidenceGraphInterface(m_var, include_inequality = False)
    var_order, con_order = define_elimination_order(var_list, con_list, igraph = igraph)
    
    igraph = IncidenceGraphInterface(m_var, include_inequality = True)
    m_reduced, _, _ = eliminate_variables(m_var, var_order, con_order, igraph = igraph)
    
    nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m_reduced, tee=True)
    
    timing_data[str(i) + ', num_vars'] = nv
    timing_data[str(i) + ', num_cons'] = nc
    timing_data[str(i) + ', num_ineq'] = ni
    timing_data[str(i) + ', ipopt_linsolve'] = linsolvetime
    timing_data[str(i) + ', nlp_eval_time'] = nlptime
    timing_data[str(i) + ', heuristic_time'] = t1
    timing_data[str(i) + ', NZ_eq_jac'] = nze
    timing_data[str(i) + ', NZ_ineq_jac'] = nzie
    timing_data[str(i) + ', NZ_lag_hess'] = nzlh
    timing_data[str(i) + ', iter'] = niter
    
    with open(path_vm + "timing_data_" + str(i) + ".txt", 'w') as f:
        f.write(json.dumps(timing_data))
    #Con major elimination
    #Make model again to make sure it is not initialized to the previous solution
    timing_data = {}
    m_con = create_instance(horizon = 1500, nfe = 400)
    
    t0 = time.time()
    var_list, con_list, cons_con = identify_vars_for_elim_min_degree(m_con, major_elim = "Constraints", eliminate_bounded_vars=False, eliminate_linear_cons_only=False, constraint_ordering= ids)
    t1 = time.time() - t0
    
    igraph = IncidenceGraphInterface(m_con, include_inequality = False)
    var_order, con_order = define_elimination_order(var_list, con_list, igraph = igraph)
    
    igraph = IncidenceGraphInterface(m_con, include_inequality = True)
    m_reduced, _, _ = eliminate_variables(m_con, var_order, con_order, igraph = igraph)
    
    nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m_reduced, tee=True)
    timing_data[str(i) + ', num_vars'] = nv
    timing_data[str(i) + ', num_cons'] = nc
    timing_data[str(i) + ', num_ineq'] = ni
    timing_data[str(i) + ', ipopt_linsolve'] = linsolvetime
    timing_data[str(i) + ', nlp_eval_time'] = nlptime
    timing_data[str(i) + ', heuristic_time'] = t1
    timing_data[str(i) + ', NZ_eq_jac'] = nze
    timing_data[str(i) + ', NZ_ineq_jac'] = nzie
    timing_data[str(i) + ', NZ_lag_hess'] = nzlh
    timing_data[str(i) + ', iter'] = niter
    
    with open(path_cm + "timing_data_" + str(i) + ".txt", 'w') as f:
        f.write(json.dumps(timing_data))
        
    for i in range(0, len(cons)):
        assert cons_ampl[i].name == cons_var[i].name
        assert cons_ampl[i].name == cons_con[i].name