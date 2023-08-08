import json
from var_elim.heuristics.ampl_heuristic import identify_vars_for_elim_ampl
from var_elim.heuristics.min_degree_heuristic import identify_vars_for_elim_min_degree
from var_elim.mip_formulations.mip_elim import identify_vars_for_elim_mip
from var_elim.algorithms.replace import (
    define_elimination_order,
    eliminate_variables
)
from var_elim.util.ipopt_helper import ipopt_solver_timing_records
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
import sys
import pyomo.environ as pyo
import time

path = r"/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/models/cute/"
with open(path + "acceptable_cute_problems.txt",'r') as f:
    acceptable_module_names = json.load(f)
print(acceptable_module_names)
import pdb; pdb.set_trace()
path = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/models/cute/pyomo_model_libraries_main/cute"
sys.path.append(path)
path2 = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/mip_compare_data/"

ipopt = pyo.SolverFactory('ipopt')

methods = ['no_elim','ampl', 'var_major', 'con_major', 'mip']
for name in acceptable_module_names:
    timing_data = {}
    module = __import__(name)
    print(name)
    
    for method in methods:
        if method == 'no_elim':
            try:
                m = module.model
                nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m, tee=True)
                t1 = None
            except:
                import pdb;pdb.set_trace()
            
        elif method == 'ampl':
            m = module.model
            t0 = time.time()
            var_list, con_list = identify_vars_for_elim_ampl(m , eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
            t1 = time.time() - t0
        elif method == 'var_major':
            m = module.model
            t0 = time.time()
            var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = 'Variables', eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
            t1 = time.time() - t0
        elif method == 'con_major':
            m = module.model
            t0 = time.time()
            var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim ='Constraints', eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
            t1 = time.time() - t0
        elif method == 'mip':
            m = module.model
            t0 = time.time()
            var_list, con_list = identify_vars_for_elim_mip(m)
            t1 = time.time() - t0
        
        if method != 'no_elim':
            igraph = IncidenceGraphInterface(m, include_inequality = False)
            var_order, con_order = define_elimination_order(var_list, con_list, igraph = igraph)
            igraph = IncidenceGraphInterface(m, include_inequality = True)
            m_reduced, _, _ = eliminate_variables(m, var_order, con_order, igraph = igraph)
            
            
            nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m_reduced, tee=True)
        timing_data[str(method) + ', num_vars'] = nv
        timing_data[str(method) + ', num_cons'] = nc
        timing_data[str(method) + ', num_ineq'] = ni
        timing_data[str(method) + ', ipopt_linsolve'] = linsolvetime
        timing_data[str(method) + ', nlp_eval_time'] = nlptime
        timing_data[str(method) + ', heuristic_time'] = t1
        timing_data[str(method) + ', NZ_eq_jac'] = nze
        timing_data[str(method) + ', NZ_ineq_jac'] = nzie
        timing_data[str(method) + ', NZ_lag_hess'] = nzlh
        timing_data[str(method) + ', iter'] = niter
    with open(path2+"mip_compare_" +name+".txt", 'w') as f:
        f.write(json.dumps(timing_data))
        
    
    
    
