import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Project 4a
sa = np.array([-1,1])
sb = np.copy(sa)
sd = np.copy(sa)
sc = np.copy(sa)

degeneracy = {}

absmag = 0
for a in sa:
    for b in sb:
        for c in sc:
            for d in sd:
                magnetization = a + b + c + d
                absmag =+ abs(magnetization)
                energy = 2*(a*d+a*b+b*c+c*d)
                state = 'M = '+ str(magnetization) + '  E = ' + str(energy)
                if state in degeneracy:
                    degeneracy[state] += 1
                else:
                    degeneracy.update({state:1})
'''
test = np.loadtxt('./results/file', skiprows = 1)
infile = open('./results/file','r')
dims = infile.readline().strip('\n').split(' ')
test = test.reshape((-1, int(dims[0]), int(dims[1])))

for i in np.arange(0, 1000):
    plt.imshow(test[i])
    plt.pause(0.00001)
    plt.cla()

plt.show()
'''

energies = np.loadtxt('./results/20x20/all_energies', delimiter = ',')
print(np.max(np.abs(energies)))
plt.hist(energies, bins = 800)
print(energies)
plt.show()


'''
temp = 3
partition_func = 4*np.cosh(8/temp) + 12
analytical_abs_mag = 8/partition_func*(np.exp(8/temp) +2)
analytical_mean_energy = -32/partition_func*np.sinh(8/temp)

mean_mag = np.loadtxt('./results/ordered_mean_mag', delimiter = ',')
plt.plot(mean_mag)
plt.show()

mean_abs_mag = np.loadtxt('./results/ordered_mean_abs_mag', delimiter = ',')
plt.plot(mean_abs_mag*4)
plt.plot(np.arange(0, len(mean_abs_mag)), np.ones(len(mean_abs_mag))*analytical_abs_mag, '--')
plt.show()

mean_energy = np.loadtxt('./results/ordered_mean_energy', delimiter = ',')
plt.plot(mean_energy*4)
plt.plot(np.arange(0, len(mean_energy)), np.ones(len(mean_energy))*analytical_mean_energy, '--')
plt.xlabel('Monte Carlo Cycles')
plt.ylabel('Energy') # !!!!!!!!!!!!!!!!! UNITS
plt.show()

accepted_flips = np.loadtxt('./results/ordered_accepted_flips', delimiter = ',')
plt.show()

susceptibility = np.loadtxt('./results/ordered_susceptibility', delimiter = ',')
plt.plot(susceptibility, 'c')
plt.show()

heatcap = np.loadtxt('./results/ordered_heatcap', delimiter = ',')
plt.plot(heatcap)
plt.show()

mean_mag = np.loadtxt('./results/random_mean_mag', delimiter = ',')
plt.plot(mean_mag)
plt.show()

mean_abs_mag = np.loadtxt('./results/random_mean_abs_mag', delimiter = ',')
plt.plot(mean_abs_mag*4)
plt.plot(np.arange(0, len(mean_abs_mag)), np.ones(len(mean_abs_mag))*analytical_abs_mag, '--')
plt.show()

mean_energy = np.loadtxt('./results/random_mean_energy', delimiter = ',')
plt.plot(mean_energy*4)
plt.plot(np.arange(0, len(mean_energy)), np.ones(len(mean_energy))*analytical_mean_energy, '--')
plt.xlabel('Monte Carlo Cycles')
plt.ylabel('Energy') # !!!!!!!!!!!!!!!!! UNITS
plt.show()
'''
