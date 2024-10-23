import numpy as np
import matplotlib.pyplot as plt

x1,y1,z1,x2,y2,z2,x3,y3,z3 = np.loadtxt('results.dat', unpack = True)

ax = plt.figure().add_subplot(projection='3d')
ax.plot(x1,y1,z1,'r')
ax.plot(x2,y2,z2,'b')
ax.plot(x3,y3,z3,'g')
plt.show()
