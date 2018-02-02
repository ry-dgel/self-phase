import sys
import yaml
from itertools import groupby
from matplotlib import pyplot as plt
import numpy as np

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

# x axese
dt     = p["tmax"]/(p["Nt"]-1)        # Time step
points = np.arange(-p["Nt"]/2,p["Nt"]/2,1)
t      = points*dt  # Time grid iterator
f      = points/p["tmax"]    # Frequency grid 0 centered
f      = f + 299792458/p["λ"]
# frequency axis
# Initialize figure
fig = plt.figure()
fig.suptitle(r"Simulation of: $\lambda = %.0f$ nm, $E = %.0f$ $\mu$J, $\Delta t_{fwhm} = %.0f$ fs, $P=%.0f$ Bar" % (p["λ"]*1E9, p["Energy"]*1E6, p["Tfwhm"]*1E15, p["pressure"]))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# Position in Color Space
cpos = np.linspace(0.2, 0.8, len(data_list))

if len(sys.argv) > 2 and sys.argv[2] == "last":
    data_list = [data_list[-1]]
# Plot each timestep
for c, data_set in zip(cpos, data_list):
    re, im = data_set.T
    ax1.plot(t*1e15,1e-4*np.sqrt(np.power(re,2) + np.power(im,2)),
             color=plt.cm.viridis(c))
ax1.set_title("Pulse Shape")
ax1.set_xlabel("t (fs)")
ax1.set_ylabel(r"Intensity (W$\cdot$cm$^{-2}$)")


for c, data_set in zip(cpos, data_list):
    re, im = data_set.T
    Ef = np.fft.fftshift(np.fft.fft(np.fft.fftshift(re+im*1j)))
    re = np.real(Ef)
    im = np.imag(Ef)
    ax2.plot(f*1e-12,np.sqrt(np.power(re,2) + np.power(im, 2)),
             color=plt.cm.viridis(c))
# THIS IS DIRTY, FIX THIS WHEN YOU'RE LESS LAZY
ax2.set_xlim(np.array([-100,100]) + 299792458/p["λ"]*1e-12)
ax2.set_title("Pulse Spectrum")
ax2.set_xlabel("f (THz)")
ax2.set_ylabel("Energy Density (a.u)")

plt.show()
