# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 19:23:24 2024

@author: XVkou
"""

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Crear una figura y un eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Lee el archivo de texto
df = pd.read_csv('data_output.dat', delim_whitespace=True, header=None)


# Selecciona solo los elementos en posiciones impares
posiciones_p1 = df[1::3]  # Comienza desde el índice 1 y toma cada segundo elemento
posiciones_p2 = df[0::3]  # Comienza desde el índice 1 y toma cada segundo elemento
posiciones_p3 = df[2::3]  # Comienza desde el índice 1 y toma cada tercer elemento

# Plot de las tres trayectorias
ax.plot(posiciones_p1[0], posiciones_p1[1], posiciones_p1[2], c='r', label='Partícula 1')

ax.plot(posiciones_p2[0], posiciones_p2[1], posiciones_p2[2], c='b', label='Partícula 2')

ax.plot(posiciones_p3[0], posiciones_p3[1], posiciones_p3[2], c='y', label='Partícula 2')

# Etiquetas de los ejes
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')


# Añadir leyenda
ax.legend()

# Mostrar el gráfico
plt.show()

"""
#If you want to animate the data I used this code from another lecture, but it only works in google colab 
    
    
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
import os
from google.colab import drive
drive.mount('/content/drive')


path = '/content/drive/MyDrive/Fortran'
os.chdir(path)

# Cargar los datos de posiciones desde el archivo
df = pd.read_csv('data_output.dat', delim_whitespace=True, header=None)

# Dividir en tres grupos de partículas
posiciones_p1 = df[1::3].reset_index(drop=True)  # Posiciones de la partícula 1
posiciones_p2 = df[0::3].reset_index(drop=True)  # Posiciones de la partícula 2
posiciones_p3 = df[2::3].reset_index(drop=True)  # Posiciones de la partícula 3


# Número de frames (suponiendo que las partículas tienen la misma longitud de datos)
#n_frames = len(posiciones_p1)
n_frames = 100

# Crear figura y eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Configurar límites
ax.set_xlim(np.min(df[0]), np.max(df[0]))
ax.set_ylim(np.min(df[1]), np.max(df[1]))
ax.set_zlim(np.min(df[2]), np.max(df[2]))

# Trayectorias y estrellas
orbit_lines = [ax.plot([], [], [], color=c, linestyle='dashed', linewidth=1)[0] for c in ['r', 'b', 'y']]
stars = [ax.plot([], [], [], 'o', color=c)[0] for c in ['r', 'b', 'y']]

# Etiquetas de los ejes
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')

# Inicializar la animación
def init():
    for line, star in zip(orbit_lines, stars):
        line.set_data([], [])
        line.set_3d_properties([])
        star.set_data([], [])
        star.set_3d_properties([])
    return orbit_lines + stars

# Función de animación
def animate(i):
    # Para cada partícula, actualiza la línea de trayectoria y la estrella
    for j, (line, star, pos) in enumerate(zip(orbit_lines, stars, [posiciones_p1, posiciones_p2, posiciones_p3])):
        line.set_data(pos.iloc[:i, 0], pos.iloc[:i, 1])  # Trayectoria en X e Y
        line.set_3d_properties(pos.iloc[:i, 2])  # Trayectoria en Z
        star.set_data(pos.iloc[i, 0], pos.iloc[i, 1])  # Posición actual en X e Y
        star.set_3d_properties(pos.iloc[i, 2])  # Posición actual en Z
    return orbit_lines + stars

# Crear la animación
ani = FuncAnimation(fig, animate, frames=n_frames, init_func=init, interval=50, blit=True)

# Mostrar la animación
plt.close()
HTML(ani.to_html5_video())
"""