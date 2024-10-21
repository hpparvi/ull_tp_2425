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
    data = np.loadtxt('../output/' + file + '.txt')
    unique_ids = np.unique(data[:, 0])
    ids = [data[:, 0] == u_id for u_id in unique_ids]
    pos = data[:, 1:4]
    mass = data[:, 4]
    return pos, ids, mass


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
