import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Crear una figura y un eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Lee el archivo de texto

DATAPATH = '/home/methos/Documentos/Programacion2004/ull_tp_2425/exercises/gd/ex2/orbit.dat'
df = pd.read_csv(DATAPATH, delim_whitespace=True, header=None)


# Imprime el contenido del DataFrame
#print(df)

# Selecciona solo los elementos en posiciones impares
p1 = np.array(df[[1,2,3]])  # Comienza desde el índice 1 y toma cada segundo elemento
p2 = np.array(df[[4,5,6]])  # Comienza desde el índice 1 y toma cada segundo elemento
p3 = np.array(df[[7,8,9]])  # Comienza desde el índice 1 y toma cada segundo elemento
p4 = np.array(df[[10,11,12]])  # Comienza desde el índice 1 y toma cada segundo elemento
p5 = np.array(df[[13,14,15]])  # Comienza desde el índice 1 y toma cada segundo elemento
p6 = np.array(df[[16,17,18]])  # Comienza desde el índice 1 y toma cada segundo elemento
p7 = np.array(df[[19,20,21]])  # Comienza desde el índice 1 y toma cada segundo elemento
p8 = np.array(df[[22,23,24]])  # Comienza desde el índice 1 y toma cada segundo elemento
p9 = np.array(df[[25,26,27]])  # Comienza desde el índice 1 y toma cada segundo elemento
p10 = np.array(df[[28,29,30]])  # Comienza desde el índice 1 y toma cada segundo elemento


ax.plot(p1.T[0], p1.T[1], p1.T[2], alpha = 0.2, c='r', label='Partícula 1')
# Diagrama de dispersión para la segunda partícula
ax.plot(p2.T[0], p2.T[1], p2.T[2], alpha = 0.2, c='b', label='Partícula 2')
ax.plot(p3.T[0], p3.T[1], p3.T[2], alpha = 0.2, c='black', label='Partícula 3')
ax.plot(p4.T[0], p4.T[1], p4.T[2], alpha = 0.2, c='black', label='Partícula 4')
ax.plot(p5.T[0], p5.T[1], p5.T[2], alpha = 0.2, c='black', label='Partícula 5')
ax.plot(p6.T[0], p6.T[1], p6.T[2], alpha = 0.2, c='black', label='Partícula 6')
ax.plot(p7.T[0], p7.T[1], p7.T[2], alpha = 0.2, c='black', label='Partícula 7')
ax.plot(p8.T[0], p8.T[1], p8.T[2], alpha = 0.2, c='black', label='Partícula 8')
ax.plot(p9.T[0], p9.T[1], p9.T[2], alpha = 0.2, c='black', label='Partícula 9')
ax.plot(p10.T[0], p10.T[1], p10.T[2], alpha = 0.2, c='black', label='Partícula 10')

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


# Mostrar el gráfico
plt.show()