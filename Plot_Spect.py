import sys
from itertools import groupby
import yaml
from matplotlib import pyplot as plt
import numpy as np

def norm_spectrum(data_set):
    re, im = data_set.T
    Ef = np.fft.fftshift(np.fft.fft(np.fft.fftshift(re+im*1j)))
    re = np.real(Ef)
    im = np.imag(Ef)
    norm = np.sqrt(np.power(re,2) + np.power(im,2))
    return norm/max(norm)

folder = sys.argv[1]
data_list = [] # Initialize array of datasets. Each element is a time step.
with open(folder + "/E") as f:
    # Split datafile into sections separated by empty lines
    for k, g in groupby(f, lambda x: x == "\n"):
        if not k:
            # Split line of data after comma, striping whitespace.
            data_list.append(np.array([[float(x) for x in d.split(',')]
                                       for d in g if len(d.strip())]))

with open(folder + "/params") as f:
    p = yaml.load(f)
    for key, val in p.items():
        p[key] = float(val)

# x axese
dt     = p["tmax"]/(p["Nt"]-1)        # Time step
points = np.arange(-p["Nt"]/2,p["Nt"]/2,1)
t      = points*dt  # Time grid iterator
f      = points/p["tmax"]    # Frequency grid 0 centered
f      = f + 299792458/p["lambda"]
# frequency axis
# Initialize figure
spectra = np.array([norm_spectrum(data_set) for data_set in data_list])
plt.pcolormesh(spectra.T)
plt.axis([0,55,8000,9000])
plt.show()

