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

# Profiling Parallel Monte Carlo for 8 parallel processes
n_para = [40, 400, 4e3, 4e4, 4e5, 4e6, 4e7]
integral_para = np.array([
[0.207493, 0.240491, 0.172651, 0.175078, 0.207643],
[0.155835, 0.196483, 0.191916, 0.167415, 0.175165],
[0.189991, 0.187443, 0.184698, 0.183476, 0.188573],
[0.191491, 0.194402, 0.191319, 0.191336, 0.190610],
[0.192437, 0.192317, 0.192352, 0.192679, 0.192899],
[0.19274, 0.19270, 0.192798, 0.19276, 0.192722],
[0.19277, 0.192748, 0.192765, 0.192741, 0.192764]])
time_para = np.array([[0.149311, 0.148286, 0.164281,  0.155366, 0.160207],
[0.157271, 0.14939, 0.147956, 0.152362, 0.148212],
[0.158482, 0.156621, 0.149835, 0.147963, 0.156868],
[0.180678, 0.172137, 0.178925, 0.177415, 0.169764],
[0.34657, 0.347212, 0.34913, 0.424784, 0.401711],
[3.5256, 3.97431, 3.87037, 3.89731, 4.02818],
[42.8214, 40.4804, 40.6477, 42.5931, 43.4205]])

for i in range(5):
        plt.semilogx(n_para, integral_para[:,i], '.')
plt.show()

plt.plot(np.log10(n_para), np.log10(np.mean(time_para, axis = 1)), '-ko', linewidth = 0.5)
plt.xlabel('$\log_{10}(n)$')
plt.ylabel('Mean Execution Time [s]')
plt.show()
