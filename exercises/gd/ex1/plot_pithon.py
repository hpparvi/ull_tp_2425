import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Crear una figura y un eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Lee el archivo de texto

DATAPATH = '/home/methos/Documentos/Programacion2004/ull_tp_2425/exercises/gd/ex1/orbit.dat'
df = pd.read_csv(DATAPATH, delim_whitespace=True, header=None)


# Imprime el contenido del DataFrame
#print(df)

# Selecciona solo los elementos en posiciones impares
p1 = np.array(df[[1,2,3]])  # Comienza desde el índice 1 y toma cada segundo elemento
p2 = np.array(df[[4,5,6]])  # Comienza desde el índice 1 y toma cada segundo elemento
p3 = np.array(df[[7,8,9]])  # Comienza desde el índice 1 y toma cada segundo elemento


ax.plot(p1.T[0], p1.T[1], p1.T[2], alpha = 0.2, c='r', label='Partícula 1')
# Diagrama de dispersión para la segunda partícula
ax.plot(p2.T[0], p2.T[1], p2.T[2], alpha = 0.2, c='b', label='Partícula 2')
ax.plot(p3.T[0], p3.T[1], p3.T[2], alpha = 0.2, c='black', label='Partícula 3')

# Limitar los ejes
#ax.set_xlim(-3, -2)
#ax.set_ylim(-1, 1)
#ax.set_zlim(-1, 1)

# Etiquetas de los ejes
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')

ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_zlim(-2, 2)

# Añadir leyenda
ax.legend()

# Mostrar el gráfico
plt.show()