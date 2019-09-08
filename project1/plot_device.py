import numpy as np
import matplotlib.pyplot as plt
import os

for root, dirs, files in os.walk('.'):
    filenames = files


for filename in filenames:
    if '0x1' in filename:
        approx = np.loadtxt(filename)
        h = 1/(len(approx)+1)
        x = np.array([h*(i+1) for i in range(len(approx))])
        n_value = filename.split('x')[0]
        plt.plot(x, approx, '--', label = 'n = %s' %n_value)

u = 1.0 - (1.0 - np.exp(-10.0))*x  - np.exp(-10.0*x);
plt.plot(x, u, label = 'Exact Solution', linewidth = 1)
plt.xlabel('x', fontsize = 12)
plt.ylabel('$v$', fontsize = 12)
plt.legend(frameon = False)
plt.show()
