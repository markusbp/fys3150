
import numpy as np
import matplotlib.pyplot as plt

w = np.genfromtxt('w_vals', delimiter = ',')

plt.plot(w)
plt.show()
plt.hist(w, bins = 1000)
plt.show()
