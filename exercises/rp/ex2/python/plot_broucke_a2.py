import os
import numpy as np
import matplotlib.pyplot as plt

ex_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

dir_data = ex_dir + '/fortran/not_parallel/results/broucke_a2/'
files = os.listdir(dir_data)
files = [dir_data+file for file in files]

M = len(files)
pos = []

for m in range(M):
    snap_data = np.loadtxt(files[m])
    pos.append(snap_data)
pos.append(pos[0])

pos = np.array(pos)

fig, ax = plt.subplots()

colors = ['deepskyblue', 'indianred', 'gold']
for i in range(pos.shape[1]):
    ax.plot(pos[:, i, 0], pos[:, i, 1], c=colors[i])

ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax.set_aspect('equal')

bc = '#1c1c1c'
fig.patch.set_facecolor(bc)
ax.set_facecolor(bc)       

for spine in ax.spines.values():
    spine.set_visible(False)

fig.savefig(ex_dir+'/figures/broucke_a2.png', dpi=300, bbox_inches='tight')