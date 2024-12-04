#python program to make a 3dimentional plot with the results of 3 particles
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('output.dat', unpack = True)



ax = plt.figure().add_subplot(projection='3d')


n = int((np.shape(data)[0]-1)/3)
print(n)
t = data[0,:]
for i in range(n):
	xx = data[3*i+1,:]
	yy = data[3*i+2,:]
	zz = data[3*i+3,:]
	ax.plot(xx,yy,zz,alpha = 0.5)

plt.show()

