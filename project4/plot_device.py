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
