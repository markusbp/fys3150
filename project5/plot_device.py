import numpy as np
import matplotlib.pyplot as plt
# Set font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

# Grid search psi1
################################################################################
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
ax1.set_ylabel('$\langle E_{T1}\\rangle$ [a.u.]', fontsize = 15)
ax2.set_xlabel('$\\alpha$', fontsize = 15)
ax2.set_ylabel('$\sigma_{E,T1}^2 $ [a.u.$^2$]', fontsize = 17)
ax2.legend(frameon = False)
plt.show()

print('Non-Interacting:')
for w, style in zip(['0.01', '0.5', '1.00'],['--', ':', '-.']):
    psi1_params = np.genfromtxt('results/psi1/psi1_gridsearch_w=' + w, delimiter = ',')
    plt.semilogy(psi1_params[1], psi1_params[3], style, label = '$\omega =$' + w, linewidth = 1)
    mina = np.argmin(psi1_params[3])
    plt.plot(psi1_params[1,mina], psi1_params[3,mina], 'ro', markersize = 6)
    print('Min a, Min E:', psi1_params[1,mina], psi1_params[3,mina])
plt.legend(frameon = False)
plt.xlabel('$\\alpha$', fontsize = 15)
plt.ylabel('$\langle E_{T1}\\rangle$ [a.u.]', fontsize = 15)
plt.legend(frameon = False)
plt.show()


# Psi2 grid search
################################################################################

full_grid = np.genfromtxt('results/psi2/full_gridsearch_psi2', delimiter = ',')
axes = np.genfromtxt('results/psi2/full_gs_psi2_varparams', delimiter = ',')

ind = np.unravel_index(np.argmin(full_grid, axis=None), full_grid.shape)
print('Min energy:',  full_grid[ind])
print('Min at alpha, beta: %.3f, %.3f' %(axes[0,ind[0]], axes[1, ind[1]]))

a, b = np.meshgrid(axes[0], axes[1])
plt.pcolor(b, a, full_grid, cmap = 'gist_earth', vmax = 5)
plt.xlabel('$\\alpha$', fontsize = 14)
plt.ylabel('$\\beta$', fontsize = 14)
plt.colorbar()
plt.show()

# Stability "analysis"
################################################################################

fig, axs = plt.subplots(1,3)
for w, i, n in zip(['0.01', '0.5', '1.00'], range(3), [2500, 5000, 10000]):
    psi1_avg = np.genfromtxt('results/psi1/psi1_opt_w=' + w, delimiter = ',') # alph = 0.88
    psi2_avg = np.genfromtxt('results/psi2/psi2_opt_w=' + w, delimiter = ',') # beta = 0.270, alph = 1.00
    axs[i].plot(psi1_avg[:n,0], '--',label = '$\psi_{T1}$')
    axs[i].plot(psi2_avg[:n,0], ':', label = '$\psi_{T2}$')
    axs[i].set_title('$\omega = $' + w)
    axs[i].set_xlabel('n', fontsize = 15)
    axs[i].legend(frameon = False)
    print('Particle separation for Psi1:', psi1_avg[-1, 4], 'omega', w)
    print('Particle separation for Psi2:', psi2_avg[-1, 4], 'omega', w)
    print('Energy for Psi1:', psi1_avg[-1, 0], 'omega', w)
    print('Energy for Psi2:', psi2_avg[-1, 0], 'omega', w)
    print('Energy var. for Psi1:', psi1_avg[-1, 1]- psi1_avg[-1, 0]**2, 'omega', w)
    print('Energy var. for Psi2:', psi2_avg[-1, 1]- psi2_avg[-1, 0]**2, 'omega', w)
    print('\n')
axs[0].set_ylabel('$\langle E \\rangle$ [a.u.]', fontsize = 15)

plt.legend(frameon = False)
plt.show()


# Virial theorem
###########################################################################
ek_ep_omega = np.genfromtxt('results/psi2/virial', delimiter = ',')
plt.plot(ek_ep_omega[3], ek_ep_omega[0]/ek_ep_omega[1],'o', markersize = 2)
plt.plot(ek_ep_omega[3], ek_ep_omega[4]/ek_ep_omega[6],'+', markersize = 4)
plt.ylabel('$\langle T\\rangle/ \langle V \\rangle$', fontsize = 14)
plt.xlabel('$\omega$', fontsize = 14)
plt.legend(['Interacting', 'Non-Interacting'], frameon = False)
plt.show()
