import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

# Crear una figura y un eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Lee el archivo de texto
DATAPATH = '/home/methos/Documentos/Programacion2004/ull_tp_2425/exercises/gd/ex2/output.dat'
df = pd.read_csv(DATAPATH, delim_whitespace=True, header=None)

# Número de partículas
N = int((df.shape[1]-1)/3)  # Número de partículas
frames = len(df)  # Número de cuadros basado en las filas del archivo

# Configuración inicial
particles, = ax.plot([], [], [], 'bo', markersize=2, label="Partículas")  # Puntos para las partículas
trails = []  # Para las trayectorias
for _ in range(N):
    trail, = ax.plot([], [], [], alpha=1)  # Una línea para cada órbita
    trails.append(trail)

# Límites del gráfico
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_zlim(-2, 2)

# Etiquetas de los ejes
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')

# Función para inicializar la animación
def init():
    particles.set_data([], [])
    particles.set_3d_properties([])
    for trail in trails:
        trail.set_data([], [])
        trail.set_3d_properties([])
    return [particles] + trails

N_TRAIL = 25  # Número máximo de puntos en la trayectoria

def update(frame):
    positions = np.array([df.iloc[frame, 1 + 3 * i:4 + 3 * i].values for i in range(N)])  # Extraer posiciones
    particles.set_data(positions[:, 0], positions[:, 1])
    particles.set_3d_properties(positions[:, 2])
    
    # Actualizar trayectorias con rastro limitado
    for i, trail in enumerate(trails):
        start = max(0, frame - N_TRAIL)  # Empieza desde el cuadro más reciente limitado a N_TRAIL
        trail_data = df.iloc[start:frame, 1 + 3 * i:4 + 3 * i].values.T
        trail.set_data(trail_data[0], trail_data[1])
        trail.set_3d_properties(trail_data[2])
    return [particles] + trails

# Crear la animación
ani = FuncAnimation(fig, update, frames=frames, init_func=init, blit=True, interval=50)
ani.save('animacion_inf.mp4', writer='ffmpeg', fps=20)
# Mostrar la animación
plt.show()