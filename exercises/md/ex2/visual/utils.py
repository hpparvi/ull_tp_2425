import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_hex


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


# To generate a list of colors, based in a colormap
def generate_colors(N, cmap='jet'):
    values = np.linspace(0, 1, N)
    cmap = plt.get_cmap(cmap)
    colors = cmap(values)
    colors_hex = [to_hex(color) for color in colors]

    return colors_hex


# To read the data from the output file
def read_data(file):
    if '.dat' not in file:
        file += '.dat'
    if 'output' not in file:
        file = '../output/' + file
    data = np.loadtxt(file)
    with open(file, 'r') as f:
        exec_time = f.readlines()[-1].split()[-1]
    unique_ids = np.unique(data[:, 0])
    ids = [data[:, 0] == u_id for u_id in unique_ids]
    pos = data[:, 1:4]
    mass = data[:, 4]
    t = data[:, 5]
    return pos, ids, mass, t, exec_time


# To read the input file
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
    sim_info['theta'] = float(data[5].strip())

    return sim_info


# To label the plots
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
