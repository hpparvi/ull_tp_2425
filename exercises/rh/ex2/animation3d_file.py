import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Import data of time and particle positions
data = np.loadtxt('output.dat', unpack = True)

# Parameters
num_particles = int((np.shape(data)[0]-1)/3)      # Number of particles
num_steps = np.shape(data)[1]         # Number of time steps
box_size = np.max(data[1:,:])            # Size of the simulation box

positions = np.zeros((num_steps, num_particles, 3))
time = data[1,:]
for i in range(num_particles):
	positions[:,i,0] = data[3*i+1,:]
	positions[:,i,1] = data[3*i+2,:]
	positions[:,i,2] = data[3*i+3,:]

# Create the figure and 3D axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-1, 2)
ax.set_ylim(-1, 2)
ax.set_zlim(-1, 2)

# Initialize particle markers and trajectory lines
particles, = ax.plot([], [], [], 'co', markersize=5)  # Particle markers
traces = [ax.plot([], [], [], '-', alpha=0.5)[0] for _ in range(num_particles)]  # Trajectories

def init():
    particles.set_data([], [])
    particles.set_3d_properties([])
    for trace in traces:
        trace.set_data([], [])
        trace.set_3d_properties([])
    return [particles, *traces]

def update(frame):
    # Update particle positions
    particles.set_data(positions[frame, :, 0], positions[frame, :, 1])
    particles.set_3d_properties(positions[frame, :, 2])

    # Update trajectories
    for i, trace in enumerate(traces):
        trace.set_data(positions[:frame + 1, i, 0], positions[:frame + 1, i, 1])
        trace.set_3d_properties(positions[:frame + 1, i, 2])

    return [particles, *traces]

# Create the animation
anim = FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=False, interval=50)

# Show animation
plt.show()
