"""
Written by David S. Weber
dsw7@sfu.ca
"""

from numpy import arange
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.dpi'] = 100

# ------------------------------------------------------------------------------
# data obtained from analyze.py

data = {
    '0': {'PHE-PHE': 1952, 'TYR-TYR': 753, 'TRP-TRP': 172, 'TYR-TRP': 774, 'PHE-TYR': 2373, 'PHE-TRP': 1078},
    '1': {'PHE-PHE': 423, 'TYR-TYR': 170, 'TRP-TRP': 56, 'TYR-TRP': 175, 'PHE-TYR': 480, 'PHE-TRP': 267},
    '2': {'PHE-PHE': 560, 'TYR-TYR': 258, 'TRP-TRP': 66, 'TYR-TRP': 177, 'PHE-TYR': 693, 'PHE-TRP': 281},
    '3': {'PHE-PHE': 645, 'TYR-TYR': 332, 'TRP-TRP': 109, 'TYR-TRP': 368, 'PHE-TYR': 833, 'PHE-TRP': 430},
    '4': {'PHE-PHE': 148, 'TYR-TYR': 102, 'TRP-TRP': 8, 'TYR-TRP': 60, 'PHE-TYR': 196, 'PHE-TRP': 78},
    '5': {'PHE-PHE': 149, 'TYR-TYR': 42, 'TRP-TRP': 6, 'TYR-TRP': 26, 'PHE-TYR': 114, 'PHE-TRP': 45},
    '6': {'PHE-PHE': 103, 'TYR-TYR': 50, 'TRP-TRP': 12, 'TYR-TRP': 31, 'PHE-TYR': 157, 'PHE-TRP': 56}
}

types = list(data.get('0').keys())

# ------------------------------------------------------------------------------

# reshape data
ECs = []
counts = []
for ec in data:
    ECs.append(ec)
    row_EC = data.get(ec)
    row = []
    for key in row_EC:
        row.append(row_EC.get(key))
    counts.append(row)

fig, ax = plt.subplots(figsize=(4, 4))
ax.imshow(counts, cmap='coolwarm')
ax.set_xticks(arange(len(types)))
ax.set_yticks(arange(len(ECs)))
ax.set_xticklabels(types)
ax.set_yticklabels(ECs)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", va="center", rotation_mode="anchor")

# add white grid
for line in range(0, len(types)):
    plt.axvline(line - 0.5, lw=3, c='white')

for line in range(0, len(data)):
    plt.axhline(line - 0.5, lw=3, c='white')

# remove border
plt.gca().set_frame_on(False)

# annotate pixels
for i in range(len(ECs)):
    for j in range(len(types)):
        text = ax.text(j, i, counts[i][j], ha="center", va="center", color="k")

# show results
fig.tight_layout()
plt.show()

