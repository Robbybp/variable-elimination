import os
import json
import numpy as np
import matplotlib.pyplot as plt
path = "/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/mip_compare_data/"

name = []
total_vars = []
ampl_vars = []
vm_vars = []
cm_vars = []
mip_vars =[]

for file in os.listdir(path):
    with open(path + file, 'r') as f:
        data = f.read()
    timing_data = json.loads(data)
    name.append(str(file))
    total_vars.append(timing_data["no_elim, num_vars"] )
    ampl_vars.append(timing_data["ampl, num_vars"])
    vm_vars.append(timing_data["var_major, num_vars"])
    cm_vars.append(timing_data["con_major, num_vars"])
    mip_vars.append(timing_data["mip, num_vars"])
        
print(total_vars)
print(ampl_vars)
print(vm_vars)
print(cm_vars)
print(mip_vars)

timing_data ={
    'Full model': total_vars,
    'AMPL': ampl_vars,
    'var_major': vm_vars,
    'con_major': cm_vars,
    'mip_vars': mip_vars}

x = np.arange(len(name))  # the label locations
width = 0.13 # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')
colors = {'Full model': 'crimson', 'AMPL': 'darkblue', 'var_major': 'green', 'con_major':'coral', 'mip_vars':'cornflowerblue'}

for attribute, measurement in timing_data.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label = attribute, color = colors[attribute], align = 'center')
    #ax.bar_label(rects, padding=4)
    multiplier += 1
    


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Variables in the problems')
ax.set_title('Eliminating using heuristics vs MIP')
ax.set_xticks(x + width, range(1, len(name)+1))
ax.legend(loc='upper left', ncols=3)
ax.set_ylim(0, 30)
plt.savefig("mip_comparison.png", dpi = 300)
plt.show()

