import os
import numpy as np

ex_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ics_dir = ex_dir + "/ics/"
ics_filename = "uniform_sphere.dat"


N = 1e3  #Number of particles (approximate)
MASS_TOTAL = 100
R = 1  #Radius sphere

Q = 6/np.pi  # ratio between volume of a cube and sphere
N_CUBE = int((Q*N) * 1)

pos = np.random.uniform(size=(N_CUBE, 3)) * 2*R - R

r = np.linalg.norm(pos, axis=1)
r_bool = (np.abs(r)<R)

pos = pos[r_bool, :]
N = pos.shape[0]

vel = np.zeros((N, 3))

MASS = MASS_TOTAL/N

def generate_ics():

    if not os.path.exists(ics_dir):
        os.makedirs(ics_dir)

    print(f'ICS generated for {N} particles')

    ics_path = ics_dir + ics_filename
    with open(ics_path, "w") as f:
        f.write(f"{N}\n")
        for i in range(N):
            f.write(f"{MASS} ")
            f.write(f"{pos[i][0]} {pos[i][1]} {pos[i][2]} ")
            f.write(f"{vel[i][0]} {vel[i][1]} {vel[i][2]} ")
            f.write("\n")


def print_freefall_time(G):

    V_SPH = 4/3 * np.pi * R**3
    density = MASS_TOTAL / (V_SPH)

    t_ff = np.sqrt((3*np.pi) / (32 * G * density))

    print(f'Free-fall time: {t_ff:0.3e}')



if __name__ == "__main__":
    generate_ics()
    print_freefall_time(G=1.)






###
# Check sphere is uniform
###
# import matplotlib.pyplot as plt
# plt.scatter(pos[:,0], pos[:,1], s=3, lw=0, alpha=0.5)
# plt.show()
