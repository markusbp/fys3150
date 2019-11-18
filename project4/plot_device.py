import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Machinery for generating plots and plots and plots

# Set font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

'''
# Task 4b)
############################################################################
temp = 1.0
partition_func = 4*np.cosh(8/temp) + 12
analytical_abs_mag = 8/partition_func*(np.exp(8/temp) +2)/4
analytical_mean_energy = -32/partition_func*np.sinh(8/temp)/4
analytical_sus = 32/temp/partition_func*(np.exp(8/temp)+ 1)/4
analytical_heatcap = 256/(temp**2*partition_func**2)*(partition_func*np.cosh(8/temp) - 4*np.sinh(8/temp)**2)/4

results = './results/2x2/'

mean_mag = np.loadtxt(results + '_mean_abs_mag_t=1.00')
mean_en = np.loadtxt(results + '_mean_energy_t=1.00')
susc = np.loadtxt(results + '_susceptibility_t=1.00')
heatcap = np.loadtxt(results + '_heatcap_t=1.00')

mc_cycles = len(mean_mag)
n_vals = np.arange(mc_cycles)
fig, ax = plt.subplots(4, 1)
mag_error = np.log10(np.abs((mean_mag - analytical_abs_mag)/analytical_abs_mag))
ax[0].plot(n_vals, mag_error, 'slategrey', linewidth = 1)
ax[0].legend(['$\langle E \\rangle /J $'], frameon = False)

en_error = np.log10(np.abs((mean_en - analytical_mean_energy)/analytical_mean_energy))
ax[1].plot(n_vals, en_error,'firebrick', linewidth = 1)
ax[1].legend(['$\langle |M| \\rangle$'], frameon = False)

sus_error = np.log10(np.abs((susc - analytical_sus)/analytical_sus))
ax[2].plot(n_vals, sus_error, 'steelblue', linewidth = 1)
ax[2].legend(['$\hat{\chi} $'], frameon = False)

heatcap_error = np.log10(np.abs((heatcap - analytical_heatcap)/analytical_heatcap))
ax[3].plot(n_vals, heatcap_error,'k', linewidth = 1)
ax[3].legend(['$\hat{C}_V$'], frameon = False)
plt.tight_layout()
plt.xlabel('MC Cycles', fontsize = 12)
fig.text(0.025,0.5, "$\log_{10}$ Absolute Relative Error", ha="center", va="center", rotation=90)
plt.show()
'''
##########################################################################
'''

# Task 4c, energy/magnetization PER SPIN!
# t = 1.0
results = './results/20x20/ordered'
mean_mag_o_1 = np.loadtxt(results + '_mean_abs_mag_t=1.00')[:2000]
mean_en_o_1 = np.loadtxt(results + '_mean_energy_t=1.00')[:2500]
results = './results/20x20/random'
mean_mag_r_1 = np.loadtxt(results + '_mean_abs_mag_t=1.00')[:10000]
mean_en_r_1 = np.loadtxt(results + '_mean_energy_t=1.00')[:5000]

# t = 2.4
results = './results/20x20/ordered'
mean_mag_o_24 = np.loadtxt(results + '_mean_abs_mag_t=2.40')[:10000]
mean_en_o_24 = np.loadtxt(results + '_mean_energy_t=2.40')[:3000]
results = './results/20x20/random'
mean_mag_r_24 = np.loadtxt(results + '_mean_abs_mag_t=2.40')[:10000]
mean_en_r_24 = np.loadtxt(results + '_mean_energy_t=2.40')[:3000]

fig, ax = plt.subplots(2, 2)

mc_cycles = len(mean_mag_o_1)
n_vals = np.arange(mc_cycles)
ax[0,0].plot(n_vals, mean_mag_o_1,'k', label = '$\langle |M| \\rangle$')
ax[0,0].set_title('Ordered, $\hat{T} = 1.0$', fontsize = 14)
ax[0,0].legend(frameon = False, loc = 'upper right')

mc_cycles = len(mean_mag_r_1)
n_vals = np.arange(mc_cycles)
ax[0,1].plot(n_vals, mean_mag_r_1,'k', label = '$\langle |M| \\rangle$')
ax[0,1].set_title('Random, $\hat{T} = 1.0$', fontsize = 14)
ax[0,1].legend(frameon = False)

mc_cycles = len(mean_en_o_1)
n_vals = np.arange(mc_cycles)
ax[1,0].plot(n_vals, mean_en_o_1,'k', label = '$\langle E \\rangle/J$')
ax[1,0].legend(frameon = False, loc = 'upper right')

mc_cycles = len(mean_en_r_1)
n_vals = np.arange(mc_cycles)
ax[1,1].plot(n_vals, mean_en_r_1,'k', label = '$\langle E \\rangle/J$')
ax[1,1].legend(frameon = False)

fig.text(0.5, 0.005, 'MC Cycles', ha='center')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(2, 2)

mc_cycles = len(mean_en_o_24)
n_vals = np.arange(mc_cycles)
ax[0,0].plot(n_vals, mean_en_o_24,'k', label = '$\langle E \\rangle/J$')
ax[0,0].set_title('Ordered, $\hat{T} = 2.4$', fontsize = 14)
ax[0,0].legend(frameon = False)

mc_cycles = len(mean_en_r_24)
n_vals = np.arange(mc_cycles)
ax[0,1].plot(n_vals, mean_en_r_24,'k', label = '$\langle E \\rangle/J$')
ax[0,1].set_title('Random, $\hat{T} = 2.4$', fontsize = 14)
ax[0,1].legend(frameon = False)

mc_cycles = len(mean_mag_o_24)
n_vals = np.arange(mc_cycles)
ax[1,0].plot(n_vals, mean_mag_o_24,'k', label = '$\langle |M| \\rangle$')
ax[1,0].legend(frameon = False)

mc_cycles = len(mean_mag_r_24)
n_vals = np.arange(mc_cycles)
ax[1,1].plot(n_vals, mean_mag_r_24, 'k',label = '$\langle |M| \\rangle$')
ax[1,1].legend(frameon = False)

fig.text(0.5, 0.005, 'MC Cycles', ha='center')
plt.tight_layout()
#fig.text(0.04, 0.5, 'common Y', va='center', rotation='vertical')
plt.show()
'''
'''
accept_r_10 = np.genfromtxt('./results/20x20/random_accepted_flips_t=1.00', delimiter = ',')
accept_o_10 = np.genfromtxt('./results/20x20/ordered_accepted_flips_t=1.00', delimiter = ',')
accept_r_24 = np.genfromtxt('./results/20x20/random_accepted_flips_t=2.40', delimiter = ',')
accept_o_24 = np.genfromtxt('./results/20x20/ordered_accepted_flips_t=2.40', delimiter = ',')

n_values = np.arange(1,len(accept_r_10) +1)
plt.loglog(n_values, accept_r_10, '--', color = 'steelblue', linewidth = 1)
plt.loglog(n_values, accept_o_10, color = 'steelblue', linewidth = 1)
plt.loglog(n_values, accept_r_24, '--', color = 'tomato', linewidth = 1)
plt.loglog(n_values, accept_o_24, color = 'tomato', linewidth = 1)
plt.legend(['Random,  $\hat{T} = 1.0$','Ordered,  $\hat{T} = 1.0$','Random,  $\hat{T} = 2.4$', 'Ordered,  $\hat{T} = 2.4$'],frameon=  False, fontsize = 12)
plt.show()
'''
'''
# Task d
################################################################################

ll = 20*20
energies = np.loadtxt('./results/20x20/ordered_all_energies_t=1.00', delimiter = ',')[1000:]/ll
bins = np.arange(np.min(energies), np.max(energies), 2/ll) # Each energy gets bin
counts, bins = np.histogram(energies,bins)
plt.hist(bins[:-1], bins, weights=counts/len(energies), color='tab:blue') # Convert to probability
energies = np.loadtxt('./results/20x20/random_all_energies_t=2.40', delimiter = ',')[2000:]/ll
bins = np.arange(np.min(energies), np.max(energies), 2.5/ll) # Each energy gets bin
counts, bins = np.histogram(energies,bins)
plt.hist(bins[:-1], bins, weights=counts/len(energies), color='k') # Convert to probability
plt.legend(['$\hat{T} = 1.0$', '$\hat{T} = 2.4$'], frameon = False)
plt.xlabel('E/J', fontsize = 14)
plt.ylabel('P(E/J)', fontsize = 14)
plt.yscale('log')
plt.show()


'''
'''
# Task e/f
##############################################################################
temp_range = np.arange(2, 2.6, 0.025)

ylab = ['$\langle E \\rangle /(lJ) $' ,'$ \langle |M| \\rangle /l $' ,'$(J \chi)/l$' ,'$C_V/(lk_B)$' ]

ylab = ['$\langle E \\rangle /J $' ,'$ \langle |M| \\rangle $', '$\hat{\chi}$', '$\hat{C}_V$']
for i, yl in zip(range(4), ylab):
    f1 = np.genfromtxt(open("./results/parallel/40x40", "rb"), delimiter=",")
    plt.plot(temp_range,f1[i], '-+', color = 'slategrey', linewidth = 1)
    f2 = np.genfromtxt(open("./results/parallel/60x60", "rb"), delimiter=",")
    plt.plot(temp_range,f2[i], '-x', color = 'firebrick', linewidth = 1)
    f3 = np.genfromtxt(open("./results/parallel/80x80", "rb"), delimiter=",")
    plt.plot(temp_range,f3[i], '-*', color = 'steelblue', linewidth = 1)
    f4 = np.genfromtxt(open("./results/parallel/100x100", "rb"), delimiter=",")
    plt.plot(temp_range,f4[i], '-o', color = 'k', markersize = 4, linewidth = 1)
    plt.legend(["40x40", '60x60', '80x80', '100x100'], frameon = False)
    plt.xlabel('$\hat{T}$', fontsize = 16)
    plt.ylabel(yl, fontsize = 16)
    plt.tight_layout()
    plt.show()

tl = [temp_range[np.argmax(t[3])] for t in [f2, f3, f4]]
poly = np.polyfit([60, 80,100], tl, 1)
print('Estimated critical temperature: ', poly[1])
'''

cm = 'tab20b'
test = np.genfromtxt('./results/spin/spins_t=1.00', delimiter = ',') # n = 10**5
plt.imshow(test, cmap = cm)
plt.show()
test = np.genfromtxt('./results/spin/spins_t=2.00', delimiter = ',')
plt.imshow(test, cmap = cm)
plt.show()
test = np.genfromtxt('./results/spin/spins_t=2.269', delimiter = ',')
plt.imshow(test, cmap = cm)
plt.show()

test = np.genfromtxt('./results/spin/spins_t=2.30', delimiter = ',')
plt.imshow(test, cmap = cm)
plt.show()
test = np.genfromtxt('./results/spin/spins_t=2.50', delimiter = ',')
plt.imshow(test, cmap = cm)
plt.show()
