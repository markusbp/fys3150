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
legendre_integral = np.loadtxt('legendre_profile_results', delimiter = ',')
legendre_params = np.loadtxt('legendre_profile_params', delimiter = ',')
lower = legendre_params[0]
upper = legendre_params[1]

n_legendre = legendre_params[2:]

# Profile Gauss-Laguerre Quadrature
laguerre_integral = np.loadtxt('laguerre_profile_results', delimiter = ',')
laguerre_params = np.loadtxt('laguerre_profile_params', delimiter = ',')
n_laguerre = laguerre_params[2:]
plt.plot(n_legendre, legendre_integral, '-ko', linewidth = 0.5)
plt.plot(n_laguerre, laguerre_integral, '-+', linewidth = 0.5)
plt.plot(n_legendre, np.ones(len(n_legendre))*exact_val, '--')
plt.legend(['Legendre', 'Laguerre/Legendre', 'Analytical Value'])
plt.show()


legendrerror = np.log(np.abs(exact_val - legendre_integral))
laguerror = np.log(np.abs(exact_val - laguerre_integral))

plt.plot(n_legendre, legendrerror, '-o')
plt.plot(n_laguerre, laguerror, '-+')
plt.show()
