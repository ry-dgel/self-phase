import matplotlib.pyplot as plt
import numpy as np
from plot import *
import matplotlib

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

fibs = fs.load_data("/home/ry/Programming/self-phase/data_comp")
f270 = fibs.const_param("fiberD", 270E-6)
exp_data = np.genfromtxt("/home/ry/Documents/Thesis/expData.csv", delimiter=",").T

cf = 6.626e-34 * 6.242e18
p = f270.runs[0].params
dt     = p["tmax"]/(p["Nt"]-1)        # Time step
points = np.arange(-p["Nt"]/2,p["Nt"]/2+1,1)
t      = points*dt  # Time grid iterator
f      = points/p["tmax"]    # Frequency grid 0 centered
f      = f*1e15 + 299792458e9/p["lambda"]
g      = f[np.where(f > 0)]
ev     = g * cf
wl = 1239.84193/ev


fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharey=True, sharex=True)
fig.subplots_adjust(hspace=0)
ax1.set_xlim(550,950)
ax1.set_ylim(-0.1,1.1)

# First Plot
spectrum = f270.const_param("Energy", 400).runs[0].spectra()[-1][::-1]
spectrum = spectrum[np.where(f>0)]
spectrum = spectrum * 1239.84193/np.power(ev,2)
spectrum = spectrum/max(spectrum)

ax1.plot(exp_data[0,:], exp_data[2,:],'-', color="gray")
ax1.scatter(exp_data[0,:], exp_data[2,:], marker="x", s=16)
ax1.plot(wl,spectrum,color="red")
ax1.set_xticks([600,650,700,750,800,850,900])
ax1.tick_params(direction='in', width = 2, length = 8)

# Second Plot
spectrum = f270.const_param("Energy", 800).runs[0].spectra()[-1][::-1]
spectrum = spectrum[np.where(f>0)]
spectrum = spectrum * 1239.84193/np.power(ev,2)
spectrum = spectrum/max(spectrum)

ax2.plot(exp_data[3,:], exp_data[5,:],'-', dashes=[8,2], color="gray")
ax2.scatter(exp_data[3,:], exp_data[5,:], marker="x", s=16)
ax2.plot(wl,spectrum,color="red")
ax2.set_xticks([600,650,700,750,800,850,900])
ax2.tick_params(direction='in', width = 2, length = 8)

# Third Plot
spectrum = f270.const_param("Energy", 1600).runs[0].spectra()[-1][::-1]
spectrum = spectrum[np.where(f>0)]
spectrum = spectrum #* 1239.84193/np.power(ev[::-1],2)
spectrum = spectrum/max(spectrum)

ax3.plot(exp_data[6,:], exp_data[8,:],'-', dashes=[8,2], color="gray")
ax3.scatter(exp_data[6,:], exp_data[8,:], marker="o", s=16)
ax3.plot(wl,spectrum,color="red")
ax3.set_xticks([600,650,700,750,800,850,900])
ax3.tick_params(direction='in', width = 2, length = 8)

ax3.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("Intensity (A.U)")
print(f)
plt.show()
