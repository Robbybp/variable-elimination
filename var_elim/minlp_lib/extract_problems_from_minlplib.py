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
path = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/minlp_lib/minlplib_py/minlplib/py"
sys.path.append(path)
sys.setrecursionlimit(10000)
dir_list = os.listdir(path)
nlp_files = []
module_names =[]
for file in dir_list:
    if file [-3:] == '.py':
        module_names.append(file[:-3])
acceptable_module_names = []
ipopt = pyo.SolverFactory('ipopt')
for name in module_names:
    print(name)
 
    try:
        module = __import__(name)
        m = module.model
        xfrm = pyo.TransformationFactory('core.relax_integer_vars')
        xfrm.apply_to(m)
        igraph = IncidenceGraphInterface(m, linear_only = True)
        if len(igraph.variables) >= 1 and len(igraph.constraints) >= 1:
            nv, nc, ni, linsolvetime, nlptime, nze, nzie, nzlh, niter = ipopt_solver_timing_records(m, tee=True)
            if linsolvetime+nlptime >= 2:
                acceptable_module_names.append(name)
                with open("acceptable_module_names.txt", 'w') as fp:
                    fp.write(json.dumps(acceptable_module_names))
                
    except:
        pass

print(acceptable_module_names)    