import os
import sys
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib.pyplot as plt

import utils

# This script generates a set of images from the output file
# These images can be used to generate a video with the script image_to_video.py
# It uses parallel processing to generate images with different threads
# I mean, each image can be generated in a different thread, and not in sequence
# The main thing is to numerate the images in the correct order

if len(sys.argv) not in [3, 4, 6]:
    print(
        "Use: python plot_images.py <output_file> <input_file> (optional):<num_cores> (optional):<(lim_min, lim_max)>"
    )
    sys.exit(1)

# Reading the output file and the input file
output_file = str(sys.argv[1])
input_file = str(sys.argv[2])
if len(sys.argv) == 4:
    # Number of cores (threads) to use
    if str(sys.argv[3]) == 'None':
        num_cores = os.cpu_count()
    elif sys.argv[3].isnumeric() and int(sys.argv[3]) <= os.cpu_count():
        num_cores = int(sys.argv[3])
    else:
        print("Core number is invalid. Try it again")
        sys.exit(1)
elif len(sys.argv) == 6:
    # Limits of the plots
    lims = (sys.argv[4], sys.argv[5])
    num_cores = int(sys.argv[3])
else:
    num_cores = os.cpu_count()


# Positions, identification of the particles and mass
pos, ids, mass, time = utils.read_data(output_file)
# Simulation information
sim_info = utils.read_input_file(input_file)
dt = sim_info['time_step']
dt_out = sim_info['output_time_step']

# To set the limits of the plots, all the images must have the same limits
if 'lims' in locals():
    x_min, x_max = float(lims[0]), float(lims[1])
    y_min, y_max = float(lims[0]), float(lims[1])
    z_min, z_max = float(lims[0]), float(lims[1])
else:
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
if len(ids) < 15:
    colors = utils.generate_colors(len(ids), cmap='gnuplot_r')


# Function to save one image at an specific time of the simulation
def save_figure(t):
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f'$t={t* dt_out:.3f}$', fontsize=20)
    # Loop over all the particles types
    for i, particle in enumerate(ids):
        # Distance between the initial and final position of the particle
        dist_in_fin = np.linalg.norm(pos[particle, :][0] - pos[particle, :][-1])
        for j in range(3):
            # If the particle is at the same position at the beginning and at the
            # end of the simulation, don't plot the trajectory
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
    plt.close(fig)


def save_figure_nolines(t):
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f'$t={t* dt_out:.3f}$', fontsize=20)

    # Loop over all the particles types
    for i, particle in enumerate(ids):
        for j in range(3):
            ax[j].plot(
                pos[particle, xyz[j][0]][t],
                pos[particle, xyz[j][1]][t],
                marker='.',
                color='red',
                alpha=0.6,
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

    plt.close(fig)


if __name__ == '__main__':
    time_indexs = range(int(sim_info['final_time'] // dt_out))
    if len(ids) > 15:
        save_figure = save_figure_nolines
    # Parallelization: each image is generated in a different thread
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [executor.submit(save_figure, t) for t in time_indexs]
        for _ in tqdm(as_completed(futures), total=len(time_indexs)):
            pass
