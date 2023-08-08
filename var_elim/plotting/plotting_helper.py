import pandas as pd

path = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/new_data/data_pandas_new.xlsx'

xls = pd.ExcelFile(path)
sheets = xls.sheet_names
full_model = []
ampl = []
var_major = []
con_major = []
best_time = []
avg_time = []
ub_l = [3, 7, 11]
for name in sheets:
    df = pd.read_excel(xls, name)
    time = df['Total time ']
    full_model.append(time[0])
    ampl.append(time[3])
    var_major.append(time[7])
    con_major.append(time[11])    
    best_time.append(min(time))
    avg_time.append(min(time[i] for i in ub_l))
import matplotlib.pyplot as plt
import numpy as np


timing_data ={
    'Full model': full_model[:-1],
    'AMPL': ampl[:-1],
    'Min degree(Var major)': var_major[:-1],
    'Min degree(Con major)': con_major[:-1]}

x = np.arange(len(sheets[:-1]))  # the label locations
width = 0.15  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')
colors = {'Full model': 'red', 'AMPL': 'blue', 'Min degree(Var major)': 'green', 'Min degree(Con major)':'gold'}

for attribute, measurement in timing_data.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label = attribute, color = colors[attribute])
    #ax.bar_label(rects, padding=4)
    multiplier += 1
    
print(sheets)
sheet_names = ['path_con', "spps", "distill", 'pde_heat', 'gas_net', 'opt_con', 'opf']


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('IPOPT total solve time (secs)')
ax.set_title('Eliminating only linear constraints and unbounded variables')
ax.set_xticks(x + width, sheet_names)
ax.legend(loc='upper right', ncols=2)
ax.set_ylim(0, 55)
plt.savefig('elim_ub_linear.png', dpi = 300)
plt.show()

timing_data_spec = {'min_ub_linear_time': avg_time[:-1],
                    'best_time': best_time[:-1]}

x = np.arange(len(sheets[:-1]))  # the label locations
width = 0.15  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')
colors = {'min_ub_linear_time': 'red', 'best_time': 'blue'}

for attribute, measurement in timing_data_spec.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label = attribute, color = colors[attribute])
    #ax.bar_label(rects, padding=4)
    multiplier += 1
    
print(sheets)
sheet_names = ['path_con', "spps", "distill", 'pde_heat', 'gas_net', 'opt_con', 'opf']


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('IPOPT total solve time (secs)')
ax.set_title('Eliminating only linear constraints and unbounded variables')
ax.set_xticks(x + width, sheet_names)
ax.legend(loc='upper left', ncols=4)
ax.set_ylim(0, 30)

plt.show()