import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
path = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/randomized_data/data_randomized.xlsx'

xls = pd.ExcelFile(path)
sheets = xls.sheet_names

ampl = []
vm = []
cm = []
for name in sheets:
    df = pd.read_excel(xls, name)
    dev = df['Percent_deviation_from_average']
    ampl.append(max(abs(dev[0:5])))
    vm.append(max(abs(dev[5:10])))
    cm.append(max(abs(dev[10:15])))
    
df = pd.DataFrame({'ampl':ampl,
                   'vm':vm,
                   'cm':cm})


fig, ax = plt.subplots(layout='constrained')
colors = {'ampl': 'crimson', 'vm': 'darkgreen', 'cm':'gold'}
width = 0.15
multiplier = 0
x = np.arange(len(sheets))  # the label locations
for attribute, measurement in df.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label = attribute, color = colors[attribute])
    #ax.bar_label(rects, padding=4)
    multiplier += 1
    
sheet_names = ['path_con', 'distill', 'gas_net', 'pde_heat', 'opt_con', 'spps']
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Maximum percentage deviation from the mean')
ax.set_title('Randomizing constraint orders for heuristics')
ax.set_xticks(x + width, sheet_names)
ax.legend(loc='upper right', ncols=3)
ax.set_ylim(0, 140)
plt.savefig("randomized_cons.png", )
plt.show()
    
