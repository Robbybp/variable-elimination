import json
import os 
import sys
import matplotlib.pyplot as plt
import pandas as pd
path = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/new_data/"
#sys.path.append(path)
# reading the data from the file
i = 0
save_file = "data_pandas_new.xlsx"
for file in os.listdir(path):
    if file[-8:] =="data.txt":
        with open(path + file) as f:
            data = f.read()
        js = json.loads(data)
        dictn = {}
        for h in range(0,13):
            entries = []
            for key in js.keys():
                if key.startswith("{}, ".format(h)):
                    entries.append(js[key])
            dictn[h] = entries
        df = pd.DataFrame.from_dict(dictn, orient ='index')
        #df.rename(mapper= {0:"num_vars", 1:"num_columns", 2: "num_ineq", 3:"ipopt_linsolve", 4:"nlp_eval_time", 5:"heuristic_time"}, axis= 'columns')
        df.columns = ["num_vars", "num_columns", "num_ineq", "lin_solve_time(s)", "nlp_eval_time(s)", "heuristic_time(s)", "NZ_eq_jac", "NZ_ineq_jac", "NZ_lag_hess", "Iterations"]
        if i == 0:
            with pd.ExcelWriter(save_file) as writer:  
                df.to_excel(writer, sheet_name=str(file[:-4]))
        else:
            with pd.ExcelWriter(save_file,engine='openpyxl', mode='a') as writer:  
                df.to_excel(writer, sheet_name=str(file[:-4]))
        i = i+1

    