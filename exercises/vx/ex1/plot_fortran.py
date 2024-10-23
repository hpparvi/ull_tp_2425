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