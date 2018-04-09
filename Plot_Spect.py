import sys
from itertools import groupby
import yaml
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

cf = 6.626e-34 * 6.242e18

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

# y axis
dt     = p["tmax"]/(p["Nt"]-1)        # Time step
points = np.arange(-p["Nt"]/2,p["Nt"]/2+1,1)
t      = points*dt  # Time grid iterator
f      = points/p["tmax"]    # Frequency grid 0 centered
f      = f + 299792458/p["lambda"]
g      = f[np.where(f > 0)]
ev     = g * cf
wl = 1239.84193/ev

X, Y = np.meshgrid(np.linspace(0,2.5, len(data_list)), wl[wl < 1200])
# frequency axis
# Initialize figure
data_list = [norm_spectrum(data_set) for data_set in data_list]
data_list = [spectrum[np.where(f>0)] for spectrum in data_list]
data_list = [spectrum[np.where(wl<1200)] for spectrum in data_list]
data_list = [spectrum * 1239.84193/np.power(ev[np.where(wl<1200)],2) for spectrum in data_list]
spectra = np.array([[spectrum] for spectrum in data_list])[:,0,:]
fig, ax = plt.subplots(1)
plt.pcolormesh(X, Y, spectra.T, figure=fig, cmap="plasma")
ax.set_ylim(564,1130)
#ax.set_ylim(500,1100)
ax.set_xlabel("Propagation Length (m)")
ax.set_ylabel("Wavelength (nm)")
plt.show()

