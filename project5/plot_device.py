
import numpy as np
import matplotlib.pyplot as plt
# Set font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

# Grid search psi1

fig, (ax1, ax2) = plt.subplots(1,2)

for w, style in zip(['0.01', '0.5', '1.00'],['--', ':', '-.']):

    psi1_params = np.genfromtxt('results/psi1/psi1_gridsearch_w=' + w, delimiter = ',')
    ax1.semilogy(psi1_params[1], psi1_params[0], style, label = '$\omega =$' + w, linewidth = 1)
    mina = np.argmin(psi1_params[0])
    ax1.plot(psi1_params[1,mina], psi1_params[0,mina], 'ro', markersize = 6)
    print('Min a, Min E:', psi1_params[1,mina], psi1_params[0,mina])
    minb = np.argmin(psi1_params[2]-psi1_params[0]**2)
    ax2.semilogy(psi1_params[1], psi1_params[2]-psi1_params[0]**2, style,  label = '$\omega =$' + w, linewidth = 1)
    ax2.semilogy(psi1_params[1, minb], psi1_params[2,minb]-psi1_params[0,minb]**2, 'ko')
ax1.legend(frameon = False)
ax1.set_xlabel('$\\alpha$', fontsize = 15)
ax1.set_ylabel('$\langle E_{T1}\\rangle$', fontsize = 15)
ax2.set_xlabel('$\\alpha$', fontsize = 15)
ax2.set_ylabel('$\langle \sigma_{E,T1}\\rangle$', fontsize = 17)
ax2.legend(frameon = False)
plt.tight_layout()
plt.show()
'''
for w in ['0.01', '0.5', '1.00']:
    psi1_avg = np.genfromtxt('results/psi1/psi1_opt_w=' + w, delimiter = ',')
    plt.semilogy(psi1_avg[:, 2], label = '$\omega =$' + w )

plt.legend(frameon = False)
plt.xlabel('$n$', fontsize = 15)
plt.ylabel('$\langle \sigma_{E,T1}\\rangle$', fontsize = 15)
plt.show()
'''
'''
full_grid = np.genfromtxt('results/psi2/full_gridsearch_psi2_big', delimiter = ',')
axes = np.genfromtxt('results/psi2/full_gs_psi2_varparams_big', delimiter = ',')
min_alpha = axes[0,0]; min_beta = axes[1,0]
max_alpha = axes[0, -1]; max_beta = axes[1, -1]
plt.imshow(full_grid.T, extent = [min_alpha, max_alpha, max_beta, min_beta] ,cmap = 'seismic', origin = 'upper', vmax = 3.8, vmin = 3.7)
plt.colorbar()
plt.xlabel('$\\alpha$', fontsize = 12)
plt.ylabel('$\\beta$', fontsize = 12)
plt.show()

ind = np.unravel_index(np.argmin(full_grid, axis=None), full_grid.shape)
print('Min energy:',  full_grid[ind])
print('Min at alpha, beta: %.3f, %.3f' %(axes[0,ind[0]], axes[1, ind[1]]))

plt.figure()
'''
n = 0
full_grid = np.genfromtxt('results/psi2/full_gridsearch_psi2', delimiter = ',')
axes = np.genfromtxt('results/psi2/full_gs_psi2_varparams', delimiter = ',')
min_alpha = axes[0,n]; min_beta = axes[1,n]
max_alpha = axes[0, -1]; max_beta = axes[1, -1]
plt.imshow(full_grid.T, extent = [min_alpha, max_alpha, max_beta, min_beta], cmap = 'jet', origin = 'upper', vmax = 3.8, vmin = 3.72)
plt.colorbar()
plt.xlabel('$\\alpha$', fontsize = 12)
plt.ylabel('$\\beta$', fontsize = 12)
plt.show()

ind = np.unravel_index(np.argmin(full_grid, axis=None), full_grid.shape)
print('Min energy:',  full_grid[ind])
print('Min at alpha, beta: %.3f, %.3f' %(axes[0,ind[0]], axes[1, ind[1]]))

avg_psi2 = np.genfromtxt('results/psi2/psi2_opt_w=1.00', delimiter = ',')
print(avg_psi2[-1,0])
var_en = avg_psi2[:, 1] - avg_psi2[:,0]**2

plt.plot(var_en)
plt.show()

plt.plot(avg_psi2[:,0])
plt.title(':)')
plt.show()
plt.plot(avg_psi2[:,2])
plt.plot(avg_psi2[:,3])
plt.show()


# Virial theorem
###########################################################################
ek_ep_omega = np.genfromtxt('results/psi2/virial', delimiter = ',')

plt.plot(ek_ep_omega[3], ek_ep_omega[0]/ek_ep_omega[1],'o')
plt.plot(ek_ep_omega[3], ek_ep_omega[0]/ek_ep_omega[2],'*')
plt.ylabel('$\langle T\\rangle/ \langle V \\rangle$', fontsize = 12)
plt.xlabel('$\omega$', fontsize = 12)
plt.show()
