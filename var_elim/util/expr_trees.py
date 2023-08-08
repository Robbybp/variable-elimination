import dill
import pyomo.environ as pyo
import os
from var_elim.algorithms.expr import count_nodes
import json

path = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/models/opf/reduced_models/'
expr_nodes = {}
for file in os.listdir(path):
    print(file)
    with open(path + file, mode='rb') as f:
        m_reduced = dill.load(f)
        key = int(file.split('_')[-1])
        total_nodes = 0
        for c in m_reduced.component_data_objects(pyo.Constraint, active=True):
            total_nodes+=count_nodes(c.body)
        expr_nodes[key] = total_nodes

path_2 =  r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/models/opf/'
with open(path_2 + "expr_nodes.txt", 'w') as f:
    f.write(json.dumps(expr_nodes))
    