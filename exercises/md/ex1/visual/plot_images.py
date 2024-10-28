import os
import sys
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib.pyplot as plt

import utils

plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['savefig.dpi'] = 200
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Charter'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{bm}'
plt.rcParams['font.size'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['legend.edgecolor'] = 'black'
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.framealpha'] = 1
plt.rcParams['legend.fancybox'] = False
plt.rcParams['text.antialiased'] = True
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['lines.antialiased'] = True
plt.rcParams['text.antialiased'] = True

if len(sys.argv) not in [3, 4]:
    print(
        "Use: python plot_images.py <output_file> <input_file> (optional):<num_cores>"
    )
    sys.exit(1)

output_file = str(sys.argv[1])
input_file = str(sys.argv[2])
if len(sys.argv) == 4:
    num_cores = int(sys.argv[3])
else:
    num_cores = os.cpu_count()

pos, ids, mass = utils.read_data(output_file)
sim_info = utils.read_input_file(input_file)
dt = sim_info['time_step']
dt_out = sim_info['output_time_step']

x_min, x_max = np.min(pos[:, 0]) - 0.2 * np.abs(np.min(pos[:, 0])), np.max(
    pos[:, 0]
) + 0.2 * np.abs(np.max(pos[:, 0]))
y_min, y_max = np.min(pos[:, 1]) - 0.2 * np.abs(np.min(pos[:, 1])), np.max(
    pos[:, 1]
) + 0.2 * np.abs(np.max(pos[:, 1]))
z_min, z_max = np.min(pos[:, 2]) - 0.2 * np.abs(np.min(pos[:, 2])), np.max(
    pos[:, 2]
) + 0.2 * np.abs(np.max(pos[:, 2]))

xyz = [(0, 1), (0, 2), (1, 2)]
colors = utils.generate_colors(len(ids), cmap='gnuplot_r')


def save_figure(t):
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f'$t={t* dt_out:.3f}$', fontsize=20)
    for i, particle in enumerate(ids):
        dist_in_fin = np.linalg.norm(pos[particle, :][0] - pos[particle, :][-1])
        for j in range(3):
            if dist_in_fin < 1e-3:
                ax[j].plot(
                    pos[particle, xyz[j][0]][t],
                    pos[particle, xyz[j][1]][t],
                    'o',
                    color=colors[i],
                    ms=20,
                )
            else:
                ax[j].plot(
                    pos[particle, xyz[j][0]][:t],
                    pos[particle, xyz[j][1]][:t],
                    '-',
                    color=colors[i],
                    alpha=0.6,
                )
                ax[j].plot(
                    pos[particle, xyz[j][0]][t],
                    pos[particle, xyz[j][1]][t],
                    'o',
                    color=colors[i],
                    ms=10,
                )
    fig, ax = utils.labels_plots(fig, ax)

    ax[0].set_xlim(x_min, x_max)
    ax[0].set_ylim(y_min, y_max)
    ax[1].set_xlim(x_min, x_max)
    ax[1].set_ylim(z_min, z_max)
    ax[2].set_xlim(y_min, y_max)
    ax[2].set_ylim(z_min, z_max)
    fig.subplots_adjust(top=0.9, bottom=0.15, left=0.07, right=0.99, wspace=0.25)
    output_dir = './output/images_' + sim_info['simulation_name']
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, f'im_{t:01d}.png'))
    fig.close()


if __name__ == '__main__':
    time_indexs = range(int(sim_info['final_time'] // dt_out))

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [executor.submit(save_figure, t) for t in time_indexs]
        for _ in tqdm(as_completed(futures), total=len(time_indexs)):
            pass
