import numpy as np
import matplotlib.pyplot as plt

# Set font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

# plot initial wave_function
r = np.linspace(0, 10, 1e3)
alpha = 2
wave_function = np.exp(-alpha*r)
plt.plot(r, wave_function, 'k')
plt.xlabel('|$\mathbf{r}$|')
plt.ylabel('$\psi\ (\mathbf{r}$)')
plt.show()

exact_val = 5*np.pi**2/16**2

# Profile Gauss-Legendre Quadrature
legendre_integral = np.loadtxt('results/legendre_profile_results', delimiter = ',')
legendre_params = np.loadtxt('results/legendre_profile_params', delimiter = ',')
legendre_timing = np.loadtxt('results/legendre_timing', delimiter = ',')
lower = legendre_params[0]
upper = legendre_params[1]

n_legendre = legendre_params[2:]

# Profile Gauss-Laguerre Quadrature
laguerre_integral = np.loadtxt('results/laguerre_profile_results', delimiter = ',')
laguerre_params = np.loadtxt('results/laguerre_profile_params', delimiter = ',')
laguerre_timing = np.loadtxt('results/laguerre_timing', delimiter = ',')

n_laguerre = laguerre_params[2:]
plt.plot(n_legendre, legendre_integral, '-ko', linewidth = 0.5, markersize = 2)
plt.plot(n_laguerre, laguerre_integral, '-+', linewidth = 0.5)
plt.plot(n_legendre, np.ones(len(n_legendre))*exact_val, '--')
plt.legend(['Legendre', 'Laguerre/Legendre', 'Exact Value'], frameon = False)
plt.xlabel('n')
plt.ylabel('Integral')
plt.show()

plt.plot(10**np.log10(n_legendre), 10**np.log10(legendre_timing), '-ko',markersize = 2)
plt.plot(10**np.log10(n_laguerre), 10**np.log10(laguerre_timing), '-+', linewidth = 0.5)
#print(np.polyfit(np.log10(n_legendre), np.log10(legendre_timing), 1))
#print(np.polyfit(np.log10(n_laguerre), np.log10(laguerre_timing), 1))

plt.legend(['Legendre', 'Legendre/Laguerre'], frameon = False)
plt.xlabel('$\log_{10}(n)$')
plt.ylabel('$\log_{10}$(Execution Time) [s]')
plt.show()

# Error for Legendre/Laguerre
legendrerror = np.log10(np.abs(exact_val - legendre_integral))
laguerror = np.log10(np.abs(exact_val - laguerre_integral))

plt.plot(n_legendre, legendrerror, '-ko', linewidth = 0.5, markersize = 2)
plt.plot(n_laguerre, laguerror, '-+', linewidth = 0.5)
plt.legend(['Legendre', 'Legendre/Laguerre'], frameon = False)
plt.xlabel('n')
plt.ylabel('Absolute Error')
plt.show()

# Profile Brute force monte carlo
brute_timing = np.loadtxt('results/monte_carlo_brute_force_timing', delimiter = ',')
brute_integral = np.loadtxt('results/monte_carlo_brute_force_integral', delimiter = ',')
brute_n_values_monte_carlo = np.loadtxt('results/monte_carlo_brute_force_n_values', delimiter = ',')
brute_std = np.loadtxt("results/monte_carlo_brute_force_std", delimiter = ',')

plt.plot(brute_n_values_monte_carlo, brute_timing, '-', label = 'Serial', markersize = 4)
plt.xlabel('n')
plt.ylabel('Execution Time [s]')
plt.show()

# Profile importance monte carlo without parallel processes
timing = np.loadtxt('results/monte_carlo_importance_timing', delimiter = ',')
integral = np.loadtxt('results/monte_carlo_importance_integral', delimiter = ',')
n_values_monte_carlo = np.loadtxt('results/monte_carlo_n_values', delimiter = ',')
importance_std = np.loadtxt("results/monte_carlo_importance_std", delimiter = ',')


# Profiling Parallel Monte Carlo for 4!! parallel processes
timing_para = np.loadtxt('results/parallel_monte_carlo_timing', delimiter = ',')
integral_para = np.loadtxt('results/parallel_monte_carlo_integral', delimiter = ',')
n_values_para = np.loadtxt('results/parallel_monte_carlo_n_values', delimiter = ',')

plt.plot(n_values_monte_carlo, timing, '.', label = 'Serial', markersize = 4)
plt.plot(n_values_para[:, 0], np.mean(timing_para, axis = 1), '.',  label = 'Parallel', markersize = 4)
plt.legend(frameon = False)
plt.xlabel('n')
plt.ylabel('Execution time [s]')
plt.show()
