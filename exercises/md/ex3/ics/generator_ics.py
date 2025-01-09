from PIL import Image, ImageDraw, ImageFont
import numpy as np


G = 4.302e-6  # kpc km^2/s^2
G = 1


class ic:
    def __init__(
        self,
        sim_name,
        dt=0.01,
        dt_out=0.1,
        t_end=15,
        epsilon=None,
        theta=1,
    ):

        self.file_name = 'ics/ic_' + sim_name + '.txt'
        self.sim_name = sim_name
        self.dt = dt
        self.dt_out = dt_out
        self.t_end = t_end
        self.theta = theta
        self.epsilon = epsilon

        self.X = np.empty((0, 7))

        self.components = True

    def uni_sphere(self, N, R, M_total, V_max, xyz_0):
        r = np.random.uniform(0, R, size=N)
        theta = np.random.uniform(0, np.pi, size=N)
        phi = np.random.uniform(0, 2 * np.pi, size=N)

        x = r * np.sin(theta) * np.cos(phi) + xyz_0[0]
        y = r * np.sin(theta) * np.sin(phi) + xyz_0[1]
        z = r * np.cos(theta) + xyz_0[2]

        m = np.ones_like(r) * M_total / len(r)

        v_x = -V_max * np.sin(theta) * np.cos(phi)
        v_y = -V_max * np.sin(theta) * np.sin(phi)
        v_z = -V_max * np.cos(theta)

        X = np.column_stack([m, x, y, z, v_x, v_y, v_z])

        self.X = np.vstack([self.X, X])
        self.L = 2 * R

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
    
    def word(self, text, N, font=None, font_size=500, box_size = (50, 50)):
        """
        Generates initial positions, velocities, and masses for particles in a text shape.
            text (str): The text to form (e.g., "MPI").
            N (int): Total number of particles.
            font (str): Path to the font file (default is None). For example, arial.ttf
            font_size (int): Font size for the text.
            box_size (tuple): Dimensions of the simulation box (width, height).
        """

        # Set image size and create a new image with a black background
        img_size = (1000, 1000) 
        img = Image.new("L", img_size, 0)
        draw = ImageDraw.Draw(img)
        
        # Load the font, use default if none is provided
        if font is None:
            font = ImageFont.load_default()
        else:
            font = ImageFont.truetype(font, font_size)
        
        # Calculate the bounding box of the text
        bbox = draw.textbbox((0, 0), text, font=font)
        text_size = (bbox[2] - bbox[0], bbox[3] - bbox[1])
        
        # Draw the text in the center of the image
        draw.text(((img_size[0] - text_size[0]) // 2, (img_size[1] - text_size[1]) // 2), 
            text, fill=255, font=font)
        
        # Get coordinates of white pixels
        img_array = np.rot90(np.array(img).T)
        y, x = np.where(img_array > 128)  # Coordinates of white pixels
        coords = np.stack((x, y), axis=-1)
        if coords.shape[0] < N:
            print(f'Not enough particles, using {coords.shape[0]}')
            N = coords.shape[0]
            
        # Scale and center the coordinates
        coords = coords - coords.mean(axis=0)  # Center
        coords = coords / np.max(coords) * box_size[0] / 2  # Scale
        
        # Select random particles
        indices = np.random.choice(len(coords), size=N, replace=False)
        positions = coords[indices]
        
        # Generate random z-coordinates within a small range
        z_coords = np.random.uniform(-box_size[0] / 10, box_size[0] / 10, size=N)
        positions = np.hstack((positions, z_coords[:, None]))
        
        # Initialize velocities to zero
        velocities = np.zeros_like(positions)
        
        masses = np.ones(N) * 0.5
        
        X = np.column_stack((masses, positions, velocities))
        
        self.X = np.vstack([self.X, X])
        self.L = box_size[0]

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

    file_name = input("Enter the file name: ")
    theta = float(input("Enter theta: "))
    dt = float(input("Enter dt: "))
    dt_out = float(input("Enter dt_out: "))
    t_end = float(input("Enter t_end: "))

    ics = ic(file_name, theta=theta, dt=dt, dt_out=dt_out, t_end=t_end)
    
    while True:
        choice = input("Do you want to generate text or uniform spheres? (Enter 'text' or 'spheres'): ").strip().lower()

        if choice == 'text':
            text = input("Enter the text to form: ")
            N = int(input("Enter the number of particles: "))
            font = input("Enter the font path (or leave blank for default): ")
            font_size = int(input("Enter the font size: "))
            box_size = (float(input("Enter the box width: ")), float(input("Enter the box height: ")))
            ics.word(text, N, font=font if font else None, font_size=font_size, box_size=box_size)
            break
        elif choice == 'spheres':
            num_unispheres = int(input("Enter the number of uniform spheres: "))
            for i in range(num_unispheres):
                print(f"Enter parameters for unisphere {i+1}:")
                N = int(input("  Enter N (particle numbers): "))
                R = float(input("  Enter R (maximum radius of the unisphere): "))
                M_total = float(input("  Enter M_total (total mass of the sphere): "))
                V_max = float(
                    input(
                        "  Enter V_max (maximum velocity in magnitude of a particle in the sphere): "
                    )
                )
                xyz_0 = (
                    float(input("  Enter x0: ")),
                    float(input("  Enter y0: ")),
                    float(input("  Enter z0: ")),
                )
                ics.uni_sphere(N, R, M_total, V_max, xyz_0)
            break
        else:
            print("Invalid choice. Please enter 'text' or 'spheres'.")

    ics.save_file()
