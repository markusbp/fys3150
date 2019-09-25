import numpy as np
import matplotlib.pyplot as plt


test = np.fromfile('timedata')
nn = np.linspace(5, 5+10*20, 20)
plt.plot(nn, test, nn, 1/2*nn**2*20/nn[-1]**2)
plt.show()
