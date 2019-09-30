import numpy as np
import matplotlib.pyplot as plt

# Plot timing
time = np.loadtxt('./results/timedata')
n_values = np.loadtxt('./results/n_values_timing')
poly = np.polyfit(n_values, time, deg = 3)

# Fit polynomial of deg. 3 to execution time
res = np.polyval(poly, n_values)
plt.plot(n_values, time,'k', n_values, res, '--r')
name = '%.3f' %poly[0] + ' ' + '%.3f' %poly[1] + 'n '
for val in range(2, len(poly)):
    name += '+ %.3f' %poly[val] + 'n^' + str(val) + ' '
plt.legend(['Measured', name ], frameon = False)
plt.xlabel('N', fontsize = 12)
plt.ylabel('Execution Time (s)', fontsize = 12)
plt.show()

# Plot error surface
rho_max = np.loadtxt('./results/max_rho_values', delimiter = ',') # load csv
eigen_val_diff = np.loadtxt('./results/mean_eigenvalue_error', delimiter = ',')
n_values = np.loadtxt('./results/n_values_error_surface', delimiter = ',')
log_error = np.log(eigen_val_diff)
maxval = np.max(np.abs(log_error))
plt.pcolormesh(rho_max, n_values, log_error, cmap = 'seismic', vmin = -maxval, vmax = maxval)
plt.xlabel('$\\rho_{max}$', fontsize = 12)
plt.ylabel('n', fontsize = 12)
plt.colorbar()
plt.show()

# 2 electron wave function
rho = np.loadtxt('./results/eigenvecs_rho_2')
omegas = [0.01, 0.25, 0.5, 1, 5]
for omega in omegas:
    eigenvalues = np.loadtxt('./results/eigenvals_2_omega_' + str(omega))
    wave_functions = np.loadtxt("./results/eigenvecs_2_omega_" + str(omega), delimiter = ',')
    for i in range(1):
        plt.plot(rho, wave_functions[:, i]**2, label = 'n =%.2f'%i +' $\omega_r$ = %.2f'%omega)
        print('Lowest lying eigenvalue: \n', eigenvalues[0], 'for omega = ', omega)
        print('Norm of eigenstate: ', np.trapz(wave_functions[:, i]**2, rho))
plt.legend(frameon = False)
plt.xlabel('$\\rho$', fontsize = 12)
plt.ylabel('$|u(\\rho)|^2$', fontsize = 12)
plt.show()

# n = 0, compare with analytical solution
wave_functions = np.loadtxt('./results/eigenvecs_2_omega_0.25', delimiter = ',')
u2 = lambda r: r*np.exp(-r**2/8)*(1 + r/2) # Result from M. Taut
u2 = u2(rho)/np.sqrt(np.trapz(u2(rho)**2, rho))
plt.plot(rho, wave_functions[:, 0]**2, linewidth = 5)
plt.plot(rho, u2**2)
plt.legend(['$\omega_r$ = 0.25, n = 0', 'Analytical Solution'], frameon = False)
plt.xlabel('$\\rho$', fontsize = 12)
plt.ylabel('$|u(\\rho)|^2$', fontsize = 12)
plt.show()
