import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


dt = 0.01
t_end = 100
n_it = 10000

bodies = np.array([1,
                   5,
                   10,
                   50,
                   100,
                   200,
                   400,
                   1000,
                   1500,
                   2000,
                   1,
                   5,
                   10,
                   50,
                   100,
                   200,
                   400,
                   1000,
                   1500,
                   2000,
                   6,
                   10,
                   50,
                   100,
                   200,
                   400,
                   1000,
                   1500,
                   2000])
                   
                   
                   
                   
parallel = np.array([True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     False,
                     True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     True,
                     ])

type_sim = np.array([0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     2,
                     2,
                     2,
                     2,
                     2,
                     2,
                     2,
                     2,
                     2,
                     ])
                     
                     
                     
data = np.array([[0.90,  0.85,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.83,  0.00,  0.00,  0.04, np.nan],
                 [0.94,  0.84,  0.01,  0.01,  0.03,  0.01,  0.00,  0.00,  0.78,  0.01,  0.00,  0.09, np.nan],
                 [1.18,  0.99,  0.01,  0.00,  0.14,  0.02,  0.02,  0.00,  0.79,  0.04,  0.01,  0.14, np.nan],
                 [2.53,  1.81,  0.02,  0.00,  0.49,  0.10,  0.09,  0.00,  1.10,  0.18,  0.02,  0.53, np.nan],
                 [3.98,  2.66,  0.04,  0.00,  0.84,  0.17,  0.14,  0.00,  1.46,  0.29,  0.05,  0.97, np.nan],
                 [8.47,  5.75,  0.05,  0.00,  1.55,  0.29,  0.24,  0.01,  3.62,  0.52,  0.14,  2.05, np.nan],
                 [16.78, 11.78,  0.07,  0.00,  3.14,  0.56,  0.47,  0.02,  7.51,  0.97,  0.25,  3.76, np.nan],
                 [38.92, 27.44,  0.15,  0.00,  8.22,  1.48,  1.20,  0.03, 16.35,  2.56,  0.56,  8.34, np.nan],
                 [61.45, 44.29,  0.24,  0.00, 12.31,  2.15,  1.83,  0.05, 27.70,  3.92,  0.80, 12.40, np.nan],
                 [89.64, 66.21,  0.34,  0.00, 17.96,  2.83,  2.73,  0.08, 42.25,  5.37,  1.07, 16.94, np.nan],
                 [0.05,  0.01,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.03, np.nan],
                 [0.14,  0.06,  0.01,  0.00,  0.02,  0.00,  0.01,  0.00,  0.01,  0.01,  0.01,  0.06, np.nan],
                 [0.36,  0.21,  0.00,  0.00,  0.11,  0.02,  0.01,  0.00,  0.06,  0.03,  0.01,  0.11, np.nan],
                 [2.79,  2.06,  0.02,  0.00,  0.53,  0.12,  0.07,  0.00,  1.31,  0.18,  0.02,  0.53, np.nan],
                 [6.04,  4.68,  0.03,  0.00,  0.94,  0.19,  0.12,  0.00,  3.39,  0.29,  0.06,  1.01, np.nan],
                 [12.96, 10.55,  0.04,  0.00,  1.65,  0.30,  0.23,  0.01,  8.32,  0.55,  0.09,  1.75, np.nan],
                 [31.24, 26.70,  0.07,  0.00,  3.30,  0.61,  0.50,  0.01, 22.21,  1.01,  0.17,  3.33, np.nan],
                 [93.86, 82.70,  0.14,  0.00,  7.91,  1.47,  1.27,  0.03, 71.87,  2.71,  0.41,  8.02, np.nan],
                 [165.98,150.32,  0.21,  0.00, 11.60,  1.97,  1.73,  0.06,134.74,  3.69,  0.57, 11.38, np.nan],
                 [245.87,223.88,  0.33,  0.00, 16.96,  2.88,  2.50,  0.08,201.12,  5.46,  0.78, 15.72, np.nan],
                 [  0.15,  0.10,  0.01,  0.00,  0.04,  0.01,  0.01,  0.00,  0.00,  0.01,  0.00,  0.00,  0.03],
                 [  0.36,  0.27,  0.01,  0.00,  0.14,  0.03,  0.02,  0.00,  0.01,  0.04,  0.00,  0.00,  0.03],
                 [  1.24,  1.05,  0.02,  0.00,  0.45,  0.09,  0.06,  0.00,  0.20,  0.10,  0.01,  0.01,  0.06],
                 [  2.39,  2.06,  0.03,  0.00,  0.90,  0.16,  0.12,  0.00,  0.51,  0.23,  0.01,  0.01,  0.07],
                 [  5.58,  4.75,  0.04,  0.00,  1.92,  0.33,  0.28,  0.01,  1.72,  0.63,  0.02,  0.04,  0.14],
                 [ 14.06, 12.29,  0.09,  0.00,  4.24,  0.73,  0.65,  0.02,  4.90,  1.42,  0.04,  0.05,  0.24],
                 [ 39.91, 35.39,  0.20,  0.00, 11.30,  2.12,  1.74,  0.05, 16.55,  3.87,  0.10,  0.12,  0.40],
                 [ 64.95, 58.00,  0.29,  0.00, 16.87,  3.45,  2.80,  0.07, 30.09,  6.06,  0.15,  0.16,  0.56],
                 [ 93.01, 83.26,  0.40,  0.00, 23.42,  5.06,  4.20,  0.09, 44.45,  8.66,  0.20,  0.21,  0.65]])

times = ['Total wall time', 'Making trees', 'Calculating ranges', 'Nullifying pointers', 'Finding and placing cells', 'Deleting empty leaves',
         'Calculating masses', 'Initializing accelerations', 'Calculating forces', 'Deleting trees', 'Updating particles', 'Writing results', 'Gathering data']
colors = plt.colormaps['tab20'](range(len(times)))

plt.close(1)
fig, ax = plt.subplots(num=1, figsize=(10,7))

N = np.linspace(bodies[0], bodies[9], 100)
ax.plot(N, N*np.log(N)*1.5e-2, color = 'orange', lw = 0.7)
ax.plot(N, N**2*6e-5, color = 'green', lw = 0.7)
ax.plot(N, N*np.log(N)*5e-3, color = 'orange', lw = 0.7)

legend_elements = []
parallel = np.where(type_sim == 0)[0]
non_parallel = np.where(type_sim == 1)[0]
mpi = np.where(type_sim == 2)[0]
for i_t, tim in enumerate(times):
    legend_elements.append(Line2D([0], [0], color=colors[i_t], lw=4, label=tim))
    ax.plot(bodies[parallel], data[:,i_t][parallel], color = colors[i_t], marker = 'o', ls = 'dashed')
    ax.plot(bodies[non_parallel], data[:,i_t][non_parallel], color = colors[i_t], marker = '*', ls = 'dotted')
    ax.plot(bodies[mpi], data[:,i_t][mpi], color = colors[i_t], marker = 'X', ls = 'dashdot')

    
legend_elements_2 = [Line2D([0], [0], marker='o', color='grey', label='OpenMP (12 threads)',
                            markerfacecolor='grey',ls='dashed', markersize=10),
                     Line2D([0], [0], marker='*', color='grey', label='Non-Parallel', 
                            markerfacecolor='grey', ls='dotted',  markersize=10),
                     Line2D([0], [0], marker='X', color='grey', label='MPI (6 processes)', 
                            markerfacecolor='grey', ls='dashdot',  markersize=10)]

legend_elements_3 = [Line2D([0], [0], color='orange', lw=1, label='$K\\cdot N\\cdot \\log N$'),
                     Line2D([0], [0], color='green', lw=1, label='$K\\cdot N^2$')]
legend1 = plt.legend(handles=legend_elements)
ax.add_artist(legend1)
legend2 = plt.legend(handles = legend_elements_2, loc = (0.35,0.855))
ax.add_artist(legend2)
plt.legend(handles = legend_elements_3, loc = (0.35, 0.755))

plt.title('Computational time variation with respect to number of bodies used\nSimulation time = %i;  $\\Delta t$ = %.0e;  Iterations = %i'%(t_end,dt,n_it))
plt.xlabel('Number of bodies')
plt.ylabel('Computational time [s]')
ax.minorticks_on()
plt.grid(which = 'both')

plt.show()
