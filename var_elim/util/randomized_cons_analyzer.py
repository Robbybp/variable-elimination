import json
import os
import numpy as np
import pandas as pd
path_ampl = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/path_constraints/randomized_ampl/'
path_vm = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/path_constraints/randomized_var_major/'
path_cm = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/path_constraints/randomized_con_major/'

lin_solve_ampl = []
nlp_eval_ampl = []
total_ampl = []
lin_solve_vm = []
nlp_eval_vm = []
total_vm = []
lin_solve_cm = []
nlp_eval_cm = []
total_cm = []


for file in os.listdir(path_ampl):
    with open(path_ampl + file) as f:
        timing_data = json.load(f)
        lin_solve_ampl.append(timing_data[file[-5] + ', ipopt_linsolve'])
        nlp_eval_ampl.append(timing_data[file[-5] +', nlp_eval_time'])
        total_ampl.append(timing_data[file[-5] + ', ipopt_linsolve'] + timing_data[file[-5] +', nlp_eval_time'])

print(total_ampl)
print(np.mean(total_ampl), np.var(total_ampl))



for file in os.listdir(path_ampl):
    with open(path_vm + file) as f:
        timing_data = json.load(f)
        lin_solve_vm.append(timing_data[file[-5] + ', ipopt_linsolve'])
        nlp_eval_vm.append(timing_data[file[-5] +', nlp_eval_time'])
        total_vm.append(timing_data[file[-5] + ', ipopt_linsolve'] + timing_data[file[-5] +', nlp_eval_time'])

print(total_vm)
print(np.mean(total_vm), np.var(total_vm))

for file in os.listdir(path_ampl):
    with open(path_cm + file) as f:
        timing_data = json.load(f)
        lin_solve_cm.append(timing_data[file[-5] + ', ipopt_linsolve'])
        nlp_eval_cm.append(timing_data[file[-5] +', nlp_eval_time'])
        total_cm.append(timing_data[file[-5] + ', ipopt_linsolve'] + timing_data[file[-5] +', nlp_eval_time'])
        
print(total_cm)
print(np.mean(total_cm), np.var(total_cm))

# path = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/"
# save_file = "data_randomized.xlsx"

# names = [str(file) for file in os.listdir(path) if file [-6:] != '.xlsx#']
# print(names)
# import pdb;pdb.set_trace()
# i = 0
# for name in names:
#     print(name)
#     dictn = {}
#     path_ampl = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/" + name + "/randomized_ampl/"
#     path_vm = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/" + name + "/randomized_var_major/"
#     path_cm = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/" + name + "/randomized_con_major/"
#     for j in range(0,15):
#         if j<=4:
#             path = path_ampl
#             dict_key = "ampl_"
#         elif j<=9:
#             path = path_vm
#             dict_key = "vm_"
#         else:
#             path = path_cm
#             dict_key = "cm_"
#         for file in os.listdir(path):
#             with open(path + file) as f:
#                 data = f.read()
#             js = json.loads(data)
            
#             entries = []
#             for key in js.keys():
#                 if key.startswith("{}, ".format(str(file[-5]))):
#                     entries.append(js[key])
#             dictn[dict_key + str(file[-5])] = entries
        
#             df = pd.DataFrame.from_dict(dictn, orient ='index')
#     #df.rename(mapper= {0:"num_vars", 1:"num_columns", 2: "num_ineq", 3:"ipopt_linsolve", 4:"nlp_eval_time", 5:"heuristic_time"}, axis= 'columns')
#     df.columns = ["num_vars", "num_columns", "num_ineq", "lin_solve_time(s)", "nlp_eval_time(s)", "heuristic_time(s)", "NZ_eq_jac", "NZ_ineq_jac", "NZ_lag_hess", "Iterations"]
#     if i == 0:
#         with pd.ExcelWriter(save_file) as writer:  
#             df.to_excel(writer, sheet_name=name)
#     else:
#         with pd.ExcelWriter(save_file,engine='openpyxl', mode='a') as writer:  
#             df.to_excel(writer, sheet_name=name)
#     i = i+1
