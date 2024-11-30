import numpy as np

G = 4.302e-6  # kpc km^2/s^2
G = 1


class ic:
    def __init__(
        self,
        file_name,
        sim_name="test",
        dt=0.01,
        dt_out=0.1,
        t_end=15,
        epsilon=None,
        theta=1,
    ):

        self.file_name = file_name
        self.sim_name = sim_name
        self.dt = dt
        self.dt_out = dt_out
        self.t_end = t_end
        self.theta = theta
        self.epsilon = epsilon

        self.X = np.empty((0, 7))

        self.components = True

    def hernquist_sphere(self, N, M_total, a, R_max):
        """
        Genera partículas distribuidas siguiendo un perfil de Hernquist.

        Parámetros:
            N (int): Número de partículas.
            M_total (float): Masa total del sistema (M_sun).
            a (float): Parámetro de escala del perfil de Hernquist (kpc).
            R_max (float): Radio máximo de las posiciones (kpc).

        Retorna:
            None. Agrega las partículas al atributo `self.X`.
        """

        def inverse_cdf(fraction, a):
            """
            Inversa de la CDF del perfil de Hernquist.
            """
            return a * fraction / (1 - fraction)

        # Generar fracciones de masa acumulada uniformemente distribuidas
        mass_fractions = np.random.uniform(0, 1, N)

        # Calcular radios usando la inversa de la CDF
        radii = inverse_cdf(mass_fractions, a)

        # Filtrar partículas fuera de R_max
        radii = radii[radii <= R_max]

        # Generar más partículas si el filtro elimina demasiadas
        while len(radii) < N:
            extra_fractions = np.random.uniform(0, 1, N - len(radii))
            extra_radii = inverse_cdf(extra_fractions, a)
            radii = np.concatenate((radii, extra_radii[extra_radii <= R_max]))
        radii = radii[:N]

        # Generar direcciones aleatorias en 3D
        theta = np.random.uniform(0, np.pi, N)
        phi = np.random.uniform(0, 2 * np.pi, N)

        x = radii * np.sin(theta) * np.cos(phi)
        y = radii * np.sin(theta) * np.sin(phi)
        z = radii * np.cos(theta)

        # Asignar masas iguales a todas las partículas
        m = np.ones_like(radii) * M_total / N

        # Asignar velocidades iniciales (inicialmente cero)
        v_x = np.zeros(N)
        v_y = np.zeros(N)
        v_z = np.zeros(N)

        # Combinar en un array
        X = np.column_stack([m, x, y, z, v_x, v_y, v_z])

        # Agregar al atributo de posiciones
        self.X = np.vstack([self.X, X])
        self.L = 2 * R_max  # Escala característica del sistema

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

        v_x = -V_max * np.sin(theta) * np.cos(phi)
        v_y = -V_max * np.sin(theta) * np.sin(phi)
        v_z = -V_max * np.cos(theta)

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
            r_others = np.sqrt(
                self.X[:, 1] ** 2 + self.X[:, 2] ** 2 + self.X[:, 3] ** 2
            )
            M[i] = M[i] + np.sum(self.X[:, 0][r_others <= r[i]])

        v_circ = np.sqrt(G * M / r)

        v_x = -v_circ * np.sin(angles)
        v_y = v_circ * np.cos(angles)
        v_z = np.random.normal(scale=10, size=N)

        X = np.column_stack([m, x, y, z, v_x, v_y, v_z])

        self.X = np.vstack([self.X, X])
        self.L = 2 * R_d

    def save_file(self):
        self.N = len(self.X)
        if self.epsilon is None:
            factor = 0.1
            self.epsilon = factor * self.L / self.N ** (1 / 3)
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
    ic1 = ic("ics/ic_test.txt", theta=1, dt=1e-2, dt_out=1, t_end=100)
    ic1.uni_sphere(1000, 10, 1000, 0)
    # ic1.hernquist_sphere(N=1000, M_total=1e12, a=10, R_max=100)
    # ic1.exp_disk(2000, 20.0, 0.3, 5e10)
    # ic1.uni_sphere(100, 1e9, 20.0, 1e11, N=100)
    # ic1.NFW(5932371.0, 20.0, 5e10, 0.3, 10, N=100)
    ic1.save_file()
