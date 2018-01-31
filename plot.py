import numpy as np
from matplotlib import pyplot as plt
import sys

E = np.genfromtxt(sys.argv[1], delimiter=",")
plt.plot(np.sqrt(np.power(E[:,0],2) + np.power(E[:,1],2)))
plt.show()
