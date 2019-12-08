
import numpy as np
import matplotlib.pyplot as plt

psi1_params = np.genfromtxt('results/psi1_gridsearch', delimiter = ',')
plt.plot(psi1_params[1], psi1_params[0])
plt.show()

full_grid = np.genfromtxt('results/full_gridsearch_psi2', delimiter = ',')
axes = np.genfromtxt('results/full_gs_psi2_varparams', delimiter = ',')
min_alpha = axes[0,0]; min_beta = axes[1,0]
max_alpha = axes[0, -1]; max_beta = axes[1, -1]
plt.imshow(full_grid[30:], extent = [min_alpha, max_alpha, min_beta, max_beta] ,cmap = 'seismic')
plt.colorbar()
plt.xlabel('$\\alpha$', fontsize = 12)
plt.ylabel('$\\beta$', fontsize = 12)

plt.show()

print(np.min(full_grid))
