import sys
from itertools import groupby
from matplotlib import pyplot as plt
import numpy as np

data_list = [] # Initialize array of datasets. Each element is a time step.
with open(sys.argv[1]) as f:
    # Split datafile into sections separated by empty lines
    for k, g in groupby(f, lambda x: x == "\n"):
        if not k:
            # Split line of data after comma, striping whitespace.
            data_list.append(np.array([[float(x) for x in d.split(',')]
                                       for d in g if len(d.strip())]))
# Initialize figure
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# Position in Color Space
cpos = np.linspace(0, 0.8, len(data_list))
# Plot each timestep
for c, data_set in zip(cpos, data_list):
    re, im = data_set.T
    ax1.plot(np.sqrt(np.power(re,2) + np.power(im,2)), color=plt.cm.viridis(c))


for c, data_set in zip(cpos, data_list):
    re, im = data_set.T
    Ef = np.fft.fftshift(np.fft.fft(np.fft.fftshift(re+im*1j)))
    re = np.real(Ef)
    im = np.imag(Ef)
    ax2.plot(np.sqrt(np.power(re,2) + np.power(im, 2)), color=plt.cm.viridis(c))

plt.show()
