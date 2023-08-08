from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix
from var_elim.models.gas_pipelines.gas_network_model import make_dynamic_model
from var_elim.models.distillation.distill import create_instance
from var_elim.models.pyomo_dae.pde_mathworks.pde_ex_mathworks import make_model
import matplotlib.pyplot as plt
from pyomo.core.expr.current import EqualityExpression
from var_elim.heuristics.ampl_heuristic import identify_vars_for_elim_ampl
import json
import pyomo.environ as pyo
path =  r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/pde_heat/'
m = make_model(nfe = 100, ncp = 3)
igraph = IncidenceGraphInterface(m)
matrix = get_structural_incidence_matrix(igraph.variables, igraph.constraints)
plt.figure()
plt.spy(matrix, markersize=0.1)


with open(path + 'constraint_ordering_0.txt') as f:
    constraint_ordering = json.load(f)
  
cons_orig = []
for c in m.component_data_objects(pyo.Constraint, active=True):
            if isinstance(c.expr, EqualityExpression):
                cons_orig.append(c)
        
cons = [c for _,c in sorted(zip(constraint_ordering, cons_orig))]

matrix2 = get_structural_incidence_matrix(igraph.variables, cons)
plt.figure()
plt.spy(matrix2, markersize= 0.2)
plt.show()