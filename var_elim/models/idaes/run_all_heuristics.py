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
from var_elim.models.idaes.pid_steam_tank import create_model
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables
)

from var_elim.heuristics.min_degree_heuristic import identify_vars_for_elim_min_degree
from var_elim.heuristics.ampl_heuristic import identify_vars_for_elim_ampl
import time
from var_elim.util.ipopt_helper import ipopt_solver_timing_records
import json
import sys
import dill

#Heuristics are ordered as 
# 0) full model - Non heuristic

# 1) ampl - eliminate_bounded_vars = False, eliminate_linear_only = False
# 2) ampl - eliminate_bounded_vars = True, eliminate_linear_only = False
# 3) ampl - eliminate_bounded_vars = False, eliminate_linear_only = True
# 4) ampl - eliminate_bounded_vars = True, eliminate_linear_only = True

# 5) var_major - eliminate_bounded_vars = False, eliminate_linear_only = False
# 6) var_major - eliminate_bounded_vars = True, eliminate_linear_only = False
# 7) var_major - eliminate_bounded_vars = False, eliminate_linear_only = True
# 8) var_major - eliminate_bounded_vars = True, eliminate_linear_only = True

# 9) con_major - eliminate_bounded_vars = False, eliminate_linear_only = False
# 10) con_major - eliminate_bounded_vars = True, eliminate_linear_only = False
# 11) con_major - eliminate_bounded_vars = False, eliminate_linear_only = True
# 12) con_major - eliminate_bounded_vars = True, eliminate_linear_only = True
timing_data = {}
sys.setrecursionlimit(10000)
path = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/models/idaes/reduced_models'
for i in range(0,13):
    print(i)
    m = create_model(time_set =[0, 100])
    if i == 0:
        nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m, tee=True)
        t1 = None
    elif i == 1:
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_ampl(m, eliminate_bounded_vars=False, eliminate_linear_cons_only=False)
        t1 = time.time() - t0
    elif i == 2: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_ampl(m, eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
        t1 = time.time() - t0    
    elif i == 3: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_ampl(m, eliminate_bounded_vars=False, eliminate_linear_cons_only=True)
    elif i == 4: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_ampl(m, eliminate_bounded_vars=True, eliminate_linear_cons_only=True)
        t1 = time.time() - t0
        t1 = time.time() - t0
    elif i == 5: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = "Variables", eliminate_bounded_vars=False, eliminate_linear_cons_only=False)
        t1 = time.time() - t0
    elif i == 6: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = "Variables", eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
        t1 = time.time() - t0
    elif i == 7: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = "Variables", eliminate_bounded_vars=False, eliminate_linear_cons_only=True)
        t1 = time.time() - t0
    elif i == 8: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = "Variables", eliminate_bounded_vars=True, eliminate_linear_cons_only=True)
        t1 = time.time() - t0
    elif i == 9: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = "Constraints", eliminate_bounded_vars=False, eliminate_linear_cons_only=False)
        t1 = time.time() - t0
    elif i == 10: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = "Constraints", eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
        t1 = time.time() - t0
    elif i == 11: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = "Constraints", eliminate_bounded_vars=False, eliminate_linear_cons_only=True)
        t1 = time.time() - t0
    elif i == 12: 
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = "Constraints", eliminate_bounded_vars=True, eliminate_linear_cons_only=True)
        t1 = time.time() - t0
    if i in range(1, 13):
        igraph = IncidenceGraphInterface(m, include_inequality = False)
        var_order, con_order = define_elimination_order(var_list, con_list, igraph = igraph)
        igraph = IncidenceGraphInterface(m, include_inequality = True)
        m_reduced, _, _ = eliminate_variables(m, var_order, con_order, igraph = igraph)
        file_name = '/model_' + str(i)
        filepath = path + file_name
        with open(filepath , mode='wb') as file:
            dill.dump(m_reduced, file)
        
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
    with open("run_all_heuristics_data_new.txt", 'w') as f:
        f.write(json.dumps(timing_data))
#{"0, num_vars": 151367, "0, num_cons": 149868, "0, num_ineq": 0, "0, ipopt_linsolve": 16.467, "0, nlp_eval_time": 3.573, "0, heuristic_time": null, "1, num_vars": 49466, "1, num_cons": 47968, "1, num_ineq": 0, "1, ipopt_linsolve": 10.801, "1, nlp_eval_time": 4.322, "1, heuristic_time": 6.027541637420654, "2, num_vars": 49466, "2, num_cons": 47968, "2, num_ineq": 0, "2, ipopt_linsolve": 10.446, "2, nlp_eval_time": 4.041, "2, heuristic_time": 5.777081489562988, "3, num_vars": 97402, "3, num_cons": 95904, "3, num_ineq": 0, "3, ipopt_linsolve": 7.55, "3, nlp_eval_time": 3.543, "3, heuristic_time": 5.777081489562988, "4, num_vars": 97402, "4, num_cons": 95904, "4, num_ineq": 0, "4, ipopt_linsolve": 7.447, "4, nlp_eval_time": 3.54, "4, heuristic_time": 5.596911191940308, "5, num_vars": 50966, "5, num_cons": 49467, "5, num_ineq": 0, "5, ipopt_linsolve": 10.766, "5, nlp_eval_time": 4.009, "5, heuristic_time": 79.6095917224884, "6, num_vars": 50966, "6, num_cons": 49467, "6, num_ineq": 2998, "6, ipopt_linsolve": 10.769, "6, nlp_eval_time": 3.86, "6, heuristic_time": 79.85096859931946, "7, num_vars": 98934, "7, num_cons": 97435, "7, num_ineq": 0, "7, ipopt_linsolve": 8.016, "7, nlp_eval_time": 3.569, "7, heuristic_time": 79.06965494155884, "8, num_vars": 98934, "8, num_cons": 97435, "8, num_ineq": 2998, "8, ipopt_linsolve": 8.109, "8, nlp_eval_time": 3.612, "8, heuristic_time": 80.69568419456482, "9, num_vars": 49434, "9, num_cons": 47936, "9, num_ineq": 0, "9, ipopt_linsolve": 10.334, "9, nlp_eval_time": 4.078, "9, heuristic_time": 72.42805910110474, "10, num_vars": 49435, "10, num_cons": 47936, "10, num_ineq": 2998, "10, ipopt_linsolve": 10.458, "10, nlp_eval_time": 4.033, "10, heuristic_time": 71.6800754070282}
    
    

