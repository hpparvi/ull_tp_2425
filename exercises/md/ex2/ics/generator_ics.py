import numpy as np

G = 4.302e-9  # kpc km^2/s^2
G = 1


class ic:
    def __init__(
        self,
        N,
        file_name,
        sim_name="test",
        dt=0.01,
        dt_out=0.1,
        t_end=15,
        epsilon=None,
        theta=1,
    ):
        self.N = int(N)
        self.file_name = file_name
        self.sim_name = sim_name
        self.dt = dt
        self.dt_out = dt_out
        self.t_end = t_end
        self.theta = theta
        self.epsilon = epsilon

        self.X = np.empty((0, 7))

        self.components = True

    def uni_sphere(self, N, R, M_total, V_max):
        r = np.random.uniform(0, R, size=N)
        theta = np.random.uniform(0, np.pi, size=N)
        phi = np.random.uniform(0, 2 * np.pi, size=N)

        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        m = np.ones_like(r) * M_total / len(r)
        M = np.empty_like(r)
        for i in range(len(r)):
            M[i] = np.sum(m[r <= r[i]])

        v_x = -V_max * np.sin(theta) * np.cos(phi) / 1e8
        v_y = -V_max * np.sin(theta) * np.sin(phi) / 1e8
        v_z = -V_max * np.cos(theta) / 1e8

        X = np.column_stack([m, x, y, z, v_x, v_y, v_z])

        self.X = np.vstack([self.X, X])
        self.L = 2 * R  # Characteristic length scale

    def exp_disk(self, N, R_d, z_scale, M_disk):
        r = np.random.exponential(scale=R_d, size=N)
        angles = np.random.uniform(0, 2 * np.pi, size=N)
        x = r * np.cos(angles)
        y = r * np.sin(angles)
        z = np.random.normal(scale=z_scale, size=N)

        m = np.ones_like(r) * M_disk / len(r)
        M = np.empty_like(r)
        for i in range(len(r)):
            M[i] = np.sum(m[r <= r[i]])

        v_circ = np.sqrt(G * M_disk / r)

        v_x = -v_circ * np.sin(angles) / 1e8
        v_y = v_circ * np.cos(angles) / 1e8
        v_z = np.random.normal(scale=10, size=N) / 1e8

        X = np.column_stack([m, x, y, z, v_x, v_y, v_z])

        self.X = np.vstack([self.X, X])
        self.L = 2 * R_d

    def save_file(self):
        if self.epsilon is None:
            factor = 0.1
            self.epsilon = factor * self.L / self.N ** (1 / 3)
        self.N = len(self.X)
        with open(self.file_name, 'w') as f:
            f.write("! Simulation name:\n")
            f.write(f"{self.sim_name}\n\n")

            f.write("! Time step:\n")
            f.write(f"{self.dt}\n\n")

            f.write("! Output time step:\n")
            f.write(f"{self.dt_out}\n\n")

            f.write("! Final time:\n")
            f.write(f"{self.t_end}\n\n")

            f.write("! Number of particles:\n")
            f.write(f"{self.N}\n\n")

            f.write("! Theta:\n")
            f.write(f"{self.theta}\n\n")

            f.write("! Epsilon (softening length):\n")
            f.write(f"{self.epsilon}\n\n")

            f.write("! Positions (x,y,z) (each line, a particle):\n")
            positions = self.X[:, 1:4]
            np.savetxt(f, positions, fmt="%.10f", delimiter=", ")
            f.write("\n")

            f.write("! Velocity (vx, vy, vz) (each line, a particle):\n")
            velocities = self.X[:, 4:7]
            np.savetxt(f, velocities, fmt="%.10f", delimiter=", ")
            f.write("\n")

            f.write("! Mass (each line, a particle):\n")
            masses = self.X[:, 0]
            np.savetxt(f, masses, fmt="%.10f")


if __name__ == "__main__":
    ic1 = ic(3000, "ics/ic_test.txt", theta=0.0001, dt=1e-3, dt_out=1e-1, t_end=10)
    ic1.uni_sphere(3000, 12, 1e12, 0)
    # ic1.exp_disk(100, 10.0, 0.3, 5e10)
    # ic1.uni_sphere(100, 1e9, 20.0, 1e11, N=100)
    # ic1.NFW(5932371.0, 20.0, 5e10, 0.3, 10, N=100)
    ic1.save_file()
