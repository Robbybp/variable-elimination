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
from var_elim.models.pyomo_dae.optimal_control_example.optimal_control_pyomo_dae.optimal_control import make_model


path2 = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/mip_compare_cheme/"

methods = ['no_elim','ampl', 'var_major', 'con_major', 'mip']

timing_data = {}
nfe = 7
ncp = 3

for method in methods:
    if method == 'no_elim':
        try:
            m = make_model(nfe = nfe, ncp = ncp)
            nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m, tee=True)
            t1 = None
        except:
            import pdb;pdb.set_trace()
        
    elif method == 'ampl':
        m = make_model(nfe = nfe, ncp = ncp)
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_ampl(m , eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
        t1 = time.time() - t0
    elif method == 'var_major':
        m = make_model(nfe = nfe, ncp = ncp)
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim = 'Variables', eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
        t1 = time.time() - t0
    elif method == 'con_major':
        m = make_model(nfe = nfe, ncp = ncp)
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_min_degree(m, major_elim ='Constraints', eliminate_bounded_vars=True, eliminate_linear_cons_only=False)
        t1 = time.time() - t0
    elif method == 'mip':
        m = make_model(nfe = nfe, ncp = ncp)
        t0 = time.time()
        var_list, con_list = identify_vars_for_elim_mip(m, warm_start = True, warm_start_var_order= var_order, warm_start_con_order= con_order, tee= True)
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

print(timing_data)

with open(path2+"mip_compare_opt_con.txt", 'w') as f:
    f.write(json.dumps(timing_data))
