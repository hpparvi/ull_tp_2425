import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_hex


def generate_colors(N, cmap='jet'):
    values = np.linspace(0, 1, N)
    cmap = plt.get_cmap(cmap)
    colors = cmap(values)
    colors_hex = [to_hex(color) for color in colors]

    return colors_hex


def read_data(file):
    if '.txt' not in file:
        file += '.txt'
    if 'output' not in file:
        file = 'output_' + file
    data = np.loadtxt(file)
    unique_ids = np.unique(data[:, 0])
    ids = [data[:, 0] == u_id for u_id in unique_ids]
    pos = data[:, 1:4]
    mass = data[:, 4]
    return pos, ids, mass


def read_input_file(filename):
    data = np.genfromtxt(
        filename, comments='!', dtype=None, encoding='utf-8', delimiter='\n'
    )

    data = [line for line in data if line.strip() != '']

    sim_info = {}
    sim_info['simulation_name'] = str(data[0].strip())
    sim_info['time_step'] = float(data[1].strip())
    sim_info['output_time_step'] = float(data[2].strip())
    sim_info['final_time'] = float(data[3].strip())
    sim_info['number_of_particles'] = int(data[4].strip())

    return sim_info


def labels_plots(fig, ax):
    # XY plot
    ax[0].set_xlabel('$X$')
    ax[0].set_ylabel('$Y$')
    ax[0].axis('equal')

    # XZ plot
    ax[1].set_xlabel('$X$')
    ax[1].set_ylabel('$Z$')
    ax[1].axis('equal')

    # YZ plot
    ax[2].set_xlabel('$Y$')
    ax[2].set_ylabel('$Z$')
    ax[2].axis('equal')

    fig.tight_layout()
    return fig, ax
