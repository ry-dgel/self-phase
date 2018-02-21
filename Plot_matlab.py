import sys
import yaml
from itertools import groupby
from matplotlib import pyplot as plt
import numpy as np

folder = sys.argv[1]
data_list = [] # Initialize array of datasets. Each element is a time step.
with open(folder + "/E") as f:
    # Split datafile into sections separated by empty lines
    for k, g in groupby(f, lambda x: x=="\n"):
        if not k:
            # Split line of data after comma, striping whitespace.
            data_list.append(np.array([[float(x) for x in d.split(',')]
                                        for d in g if len(d.strip()) ]))

with open(folder + "/params") as f:
    p = yaml.load(f)
    for key, val in p.items():
        p[key] = float(val)

matlablines = open("%s/inte0" % sys.argv[2]).read().splitlines()
for line in matlablines:
    line = float(line.strip())
I_init = 1e-4*np.asarray(matlablines,dtype='double')
matlablines = open("%s/inte" % sys.argv[2]).read().splitlines()
for line in matlablines:
    line = float(line.strip())
I_finl = 1e-4*np.asarray(matlablines,dtype='double')
matlablines = open("%s/spec0" % sys.argv[2]).read().splitlines()
for line in matlablines:
    line = float(line.strip())
If_init = 1e-4*np.asarray(matlablines,dtype='double')
matlablines = open("%s/spec" % sys.argv[2]).read().splitlines()
for line in matlablines:
    line = float(line.strip())
If_finl = 1e-4*np.asarray(matlablines,dtype='double')

# x axese
dt     = p["tmax"]/(p["Nt"]-1)        # Time step
points = np.arange(-p["Nt"]/2,p["Nt"]/2,1)
t      = points*dt  # Time grid iterator
f      = points/p["tmax"]    # Frequency grid 0 centered
f      = f + 299792458/p["lambda"]
# frequency axis
# Initialize figure
fig = plt.figure()
fig.suptitle(r"Simulation of: $\lambda = %.0f$ nm, $E = %.0f$ $\mu$J, $\Delta t_{fwhm} = %.0f$ fs, $P=%.0f$ Bar" % (p["lambda"]*1E9, p["Energy"]*1E6, p["Tfwhm"]*1E15, p["Pout"]))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)


# Only make sys.argv[2]# of plots in between the intitial and final.
# If set to only 1, prints the final plot
#data_list = [data_list[i] for i in map(int, np.linspace(len(data_list)-1,0,2))]
data_list = [data_list[i] for i in map(int, np.linspace(len(data_list)-1,0,2))]

# Position in Color Space
cpos = np.linspace(0.8, 0.2, len(data_list))

# Plot each timestep
re, im = data_list[0].T
ax3.plot(t*1e15,1E-4*(np.power(re,2) + np.power(im,2))[::-1],
             #color=plt.cm.viridis(c))
        )
re, im = data_list[1].T
ax1.plot(t*1e15,1E-4*(np.power(re,2) + np.power(im,2)),
             #color=plt.cm.viridis(c))
        )
ax1.plot(t*1e15, I_init)
ax3.plot(t*1e15, I_finl)
ax1.set_title("Pulse Shape")
ax1.set_xlabel("t (fs)")
ax1.set_ylabel(r"Intensity (W$\cdot$cm$^{-2}$)")

ax3.set_title("Pulse Shape")
ax3.set_xlabel("t (fs)")
ax3.set_ylabel(r"Intensity (W$\cdot$cm$^{-2}$)")

re, im = data_list[0].T
Ef = np.fft.fftshift(np.fft.fft(np.fft.fftshift(re+im*1j)))
re = np.real(Ef)
im = np.imag(Ef)
If = (np.power(re,2) + np.power(im,2))
ax4.plot(f*1e-12,(If/max(If))[::-1],
         #color=plt.cm.viridis(c))
         )
re, im = data_list[1].T
Ef = np.fft.fftshift(np.fft.fft(np.fft.fftshift(re+im*1j)))
re = np.real(Ef)
im = np.imag(Ef)
If = (np.power(re,2) + np.power(im,2))
ax2.plot(f*1e-12,If/max(If),
         #color=plt.cm.viridis(c))
         )
# THIS IS DIRTY, FIX THIS WHEN YOU'RE LESS LAZY
ax4.plot(f*1e-12, If_finl/max(If_finl))
ax2.plot(f*1e-12, If_init/max(If_init))
ax2.set_xlim(np.array([-60,60]) + 299792458/p["lambda"]*1e-12)
ax2.set_title("Pulse Spectrum")
ax2.set_xlabel("f (THz)")
ax2.set_ylabel("Energy Density (a.u)")

ax4.set_xlim(np.array([-60, 60]) + 299792458/p["lambda"]*1e-12)
ax4.set_title("Pulse Spectrum")
ax4.set_xlabel("f (THz)")
ax4.set_ylabel("Energy Density (a.u)")
plt.show()
