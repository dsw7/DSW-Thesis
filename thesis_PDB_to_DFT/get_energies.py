"""
This script gets all energies from .out files.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

out_dirs = os.listdir('orca_workdir')
out_dirs = [i for i in out_dirs if i[-4:] == '.out']

out_dict = {}
for out_deck in out_dirs:

    with open(os.path.join('orca_workdir', out_deck)) as f:
        lines = f.read().split('\n')
    
    try:
        # get index of where orbital energies start
        idx =  lines.index('ORBITAL ENERGIES')
    except ValueError:
        continue
    
    energies = lines[idx + 4: idx + 104]  # get the first 100 lines
    
    energies = [item.split(' ') for item in energies]
    energies = [[i for i in item if i != ''] for item in energies]
    
    df = pd.DataFrame(energies)
    df.columns = ['N0', '0CC', 'E(Eh)', 'E(eV)']  # no idea what N0 and 0CC mean
    
    df = df[['0CC', 'E(eV)']]
    
    for key, gb_obj in df.groupby('0CC'):
        if int(float(key)) == 2:
            HOMO = gb_obj.tail(1).astype(float).iloc[0][1]
        else:
            LUMO = gb_obj.head(1).astype(float).iloc[0][1]
            
    key_deck = float(out_deck.split('_')[2][:-4])
    out_dict[key_deck] = [HOMO, LUMO]

all_keys = sorted(list(out_dict.keys()))

y_1 = [out_dict.get(i)[0] for i in all_keys]
y_2 = [out_dict.get(i)[1] for i in all_keys]

FLIP_X = False
SAVE_FIG = True

plt.figure(figsize=(10, 8))
plt.xlabel('Translation (Angstroms)', size=17)
plt.ylabel('Energy (eV)', size=17)
if FLIP_X:
    plt.plot(all_keys[::-1], y_1, label='HOMO')
    plt.plot(all_keys[::-1], y_2, label='LUMO')
else:
    plt.plot(all_keys, y_1, 'ro', label='HOMO')
    plt.plot(all_keys, y_2, 'bx', label='LUMO')
plt.legend(fontsize=17)
if SAVE_FIG:
    plt.savefig('dft_plot.png', dpi=600)
plt.show()
