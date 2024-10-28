import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir('/Users/oscarna/Documents/Compu/ull_tp_2425/exercises/so/ex1')

dt = 0.01
print_t = 1
final_t = 2
parts = ['1.0 .9700436 -.24308753 0.0 .466203685 0.43236573 0.0\n',
         '1.0 -.9700436 .24308753 0.0 .466203685 0.43236573 0.0\n',
         '1.0 0.0 0.0 0.0 -0.93249737 -0.86473146 0.0\n']

# parts = ['1.0 1.0 0.0 0.0 0.0 1.0 0.0\n',
#          '1.0 -1.0 0.0 0.0 0.0 -1.0 0.0\n']

with open('stars.txt', 'w') as file:
    file.write('%s %s %s\n'%(dt, print_t, final_t))
    file.writelines(parts)


os.system('make clean')
os.system('make')
os.system('./ex1 stars.txt')

data  = np.loadtxt('output.txt')
plt.close(1)

n_parts = 3

fig, ax = plt.subplots(num = 1)
for i in range(n_parts):
    ax.plot(data[2*i::2*n_parts,2],  data[2*i::2*n_parts,3], color = plt.colormaps['tab10'](i), marker = '+', ls = '', label = 'm%i, part'%(i+1), alpha = 0.5)
    ax.plot(data[2*i+1::2*n_parts,2],  data[2*i+1::2*n_parts,3], color = plt.colormaps['tab10'](i), marker = 'o', ls = '', label = 'm%i, check'%(i+1), alpha = 0.5)

ax.minorticks_on()
ax.grid()

ax.set_title('Representation of particle type method and simple method\nTotal time: %.2f\nTotal number of iterations: %i'%(final_t, final_t/dt))
ax.set_xlabel('x')
ax.set_ylabel('y')

ax.legend()
plt.tight_layout()
#%% Comparison of methods

fig, ax = plt.subplots(num = 2)
for i in range(n_parts):
    ax.plot(data[2*i::2*n_parts,2] - data[2*i+1::2*n_parts,2],  data[2*i::2*n_parts,3] - data[2*i+1::2*n_parts,3], color = plt.colormaps['tab10'](i), marker = '+', ls = '', label = 'm%i, diff'%(i+1), alpha = 0.5)
ax.legend()
ax.set_title('$r_{particle} - r_{simple}$\nTotal time: %.2f\nTotal number of iterations: %i'%(final_t, final_t/dt))
ax.set_xlabel('x')
ax.set_ylabel('y')

ax.minorticks_on()
ax.grid(alpha = 0.1)
plt.tight_layout()