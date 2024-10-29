import os
import numpy as np

ex1_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ics_dir = ex1_dir + "/data/ics/"
ics_filename = "ics_broucke_a2.dat"

N = 3  # Number of particles
T = 7.702408  # Integration time
DT = 0.001  # Time-step
EPSILON = 0.001  # Smoothing length

N_SNAPSHOTS = 100

masses = [1, 1, 1]
positions = [[0.336130095, 0, 0], [0.7699893804, 0, 0], [-1.1061194753, 0, 0]]
velocities = [[0, 1.532431537, 0], [0, -0.6287350978, 0], [0, -0.9036964391, 0]]


def generate_ics():

    if not os.path.exists(ics_dir):
        os.makedirs(ics_dir)

    ics_path = ics_dir + ics_filename
    with open(ics_path, "w") as f:
        f.write(f"{N}\n")
        f.write(f"{T}\n")
        f.write(f"{DT}\n")
        f.write(f"{EPSILON}\n")
        f.write(f"{N_SNAPSHOTS}\n\n")
        for i in range(N):
            f.write(f"{masses[i]} ")
            f.write(f"{positions[i][0]} {positions[i][1]} {positions[i][2]} ")
            f.write(f"{velocities[i][0]} {velocities[i][1]} {velocities[i][2]} ")
            f.write('\n')

if __name__ == "__main__":
    generate_ics()
