
import numpy as np
import matplotlib.pyplot as plt
# Set font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

# Grid search psi1
psi1_params = np.genfromtxt('results/psi1_gridsearch', delimiter = ',')
plt.plot(psi1_params[1], psi1_params[0], 'k--')
mina = np.argmin(psi1_params[0])
plt.plot(psi1_params[1,mina], psi1_params[0,mina], 'ro', markersize = 6)
print('Min a, Min E:', psi1_params[1,mina], psi1_params[0,mina])
plt.xlabel('$\\alpha$', fontsize = 15)
plt.ylabel('$\langle E_{T1}\\rangle$', fontsize = 15)
plt.show()


full_grid = np.genfromtxt('results/full_gridsearch_psi2_big', delimiter = ',')
axes = np.genfromtxt('results/full_gs_psi2_varparams_big', delimiter = ',')
min_alpha = axes[0,0]; min_beta = axes[1,0]
max_alpha = axes[0, -1]; max_beta = axes[1, -1]
plt.imshow(full_grid, extent = [min_alpha, max_alpha,max_beta, min_beta] ,cmap = 'seismic', origin = 'upper')
plt.colorbar()
plt.xlabel('$\\alpha$', fontsize = 12)
plt.ylabel('$\\beta$', fontsize = 12)

plt.show()

print(np.min(full_grid))

plt.figure()

n = 0
full_grid = np.genfromtxt('results/full_gridsearch_psi2', delimiter = ',')
axes = np.genfromtxt('results/full_gs_psi2_varparams', delimiter = ',')
min_alpha = axes[0,n]; min_beta = axes[1,n]
max_alpha = axes[0, -1]; max_beta = axes[1, -1]
plt.imshow(full_grid[n:], extent = [min_alpha, max_alpha, min_beta, max_beta], origin = 'upper')
plt.colorbar()
plt.xlabel('$\\alpha$', fontsize = 12)
plt.ylabel('$\\beta$', fontsize = 12)
print(np.min(full_grid))
plt.show()

# Decent: beta = 1.005, al = 0.23

avg_psi2 = np.genfromtxt('results/psi2_run', delimiter = ',')

var_en = avg_psi2[:, 1] - avg_psi2[:,0]**2

plt.plot(var_en)
plt.show()

plt.plot(avg_psi2[:,2])
plt.plot(avg_psi2[:,3])
plt.show()
