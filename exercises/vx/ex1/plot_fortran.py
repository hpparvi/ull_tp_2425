# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 19:27:10 2024

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

# Imprime el contenido del DataFrame
#print(df)

# Selecciona solo los elementos en posiciones impares
posiciones_p1 = df[1::2]  # Comienza desde el índice 1 y toma cada segundo elemento
posiciones_p2 = df[0::2]  # Comienza desde el índice 1 y toma cada segundo elemento

ax.plot(posiciones_p1[0], posiciones_p1[1], posiciones_p1[2], c='r', label='Partícula 1')

# Diagrama de dispersión para la segunda partícula
ax.plot(posiciones_p2[0], posiciones_p2[1], posiciones_p2[2], c='b', label='Partícula 2')

# Limitar los ejes
#ax.set_xlim(-3, -2)
#ax.set_ylim(-1, 1)
#ax.set_zlim(-1, 1)

# Etiquetas de los ejes
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')

ax.set_xlim(5, 6.5)
ax.set_ylim(5, 6.5)
ax.set_zlim(5, 6.5)

# Añadir leyenda
ax.legend()

# Mostrar el gráfico
plt.show()