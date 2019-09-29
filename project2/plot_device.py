import numpy as np
import matplotlib.pyplot as plt

# Plot timing
test = np.loadtxt('timedata')
n_values = np.loadtxt('n_values_timing')
plt.plot(n_values, test)
plt.show()


import pandas as pd
import seaborn as sns

# Plot error surface
rho_max = np.loadtxt('max_rho_values', delimiter = ',')
eigen_val_diff = np.loadtxt('mean_eigenvalue_error', delimiter = ',')
n_values = np.loadtxt('n_values_error_surface', delimiter = ',')

plt.pcolormesh(rho_max, n_values, np.log(eigen_val_diff), cmap = 'seismic')
plt.colorbar()
plt.show()


# 2 electron wave function
wave_functions = np.loadtxt("eigenvectors_2d", delimiter = ',')
rho = np.loadtxt('eigenvecs_rho')
for i in range(3):
    plt.plot(rho, wave_functions[:, i]**2)

plt.show()
