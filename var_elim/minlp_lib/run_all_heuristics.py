import os
import sys
import pyomo.environ as pyo
from var_elim.heuristics.min_degree_heuristic import identify_vars_for_elim_min_degree
import time
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables
)
import pyomo.environ as pyo
from var_elim.heuristics.ampl_heuristic import identify_vars_for_elim_ampl
import json
from var_elim.util.ipopt_helper import ipopt_solver_timing_records
from pyomo.util.subsystems import TemporarySubsystemManager
from pyomo.contrib import appsi


path = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/minlp_lib/minlplib_py/minlplib/py"
sys.path.append(path)
sys.setrecursionlimit(1000000)
dir_list = os.listdir(path)
nlp_files = []
module_names =[]
done_names = ["waterno1_12", "squfl040-080persp", "supplychainp1_022020", "saa_2", "sepasequ_complex", "sfacloc1_4_80", "unitcommit_50_20_2_mod_8", "parabol_p", "qspp_0_11_0_1_10_1", "squfl030-100persp", "nuclear49b", "knp5-44", "ringpack_20_1", "truck"]
for file in dir_list:
    if file [-3:] == '.py':
        module_names.append(file[:-3])
path2 = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/minlp_lib/run_all_heuristics/"
path3 = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/minlp_lib/"
with open(path3 + "acceptable_module_names.txt", 'r') as fp:
    data = fp.read()
acceptable_module_names = json.loads(data)
ipopt = pyo.SolverFactory('ipopt')
for name in module_names:
    #if name in acceptable_module_names and name not in done_names:
    if name in ['junkturn']:
        timing_data = {}
        module = __import__(name)
        m2 = module.model
        for i in range(5,13):
            print(i)
            m = m2.clone()
            xfrm = pyo.TransformationFactory('core.relax_integer_vars')
            xfrm.apply_to(m)
            igraph = IncidenceGraphInterface(m)
           
            if i == 0:
                nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m, tee=True)
                t1 = None
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
                # with open(path2 + name + "_data.txt", 'w') as f:
                #     f.write(json.dumps(timing_data))
                import pdb;pdb.set_trace()
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
                t1 = time.time() - t0
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
                try:
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
                    with open(path2 + name + "_data.txt", 'w') as f:
                        f.write(json.dumps(timing_data))
                except:
                    timing_data[str(i) + ', num_vars'] = None
                    timing_data[str(i) + ', num_cons'] = None
                    timing_data[str(i) + ', num_ineq'] = None
                    timing_data[str(i) + ', ipopt_linsolve'] = None
                    timing_data[str(i) + ', nlp_eval_time'] = None
                    timing_data[str(i) + ', heuristic_time'] = None
                    timing_data[str(i) + ', NZ_eq_jac'] = None
                    timing_data[str(i) + ', NZ_ineq_jac'] = None
                    timing_data[str(i) + ', NZ_lag_hess'] = None
                    timing_data[str(i) + ', iter'] = None
                    # with open(path2 + name + "_data.txt", 'w') as f:
                    #     f.write(json.dumps(timing_data))
                    
            
           