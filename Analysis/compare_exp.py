import matplotlib.pyplot as plt
import numpy as np
from plot import *
import matplotlib

font = {'family' : 'normal',
        'size'   : 30}

matplotlib.rc('font', **font)

fibs = fs.load_data("/home/ry/Programming/self-phase/Analysis/compare_Exp/simData")
exp_data = np.genfromtxt("/home/ry/Documents/Thesis/expData.csv", delimiter=",").T



fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharey=True, sharex=True)
fig.subplots_adjust(hspace=0)
ax1.set_xlim(550,950)
ax1.set_ylim(-0.1,1.1)

# First Plot
run = fibs.runs[1]
wl = run.make_wavelength_scale()
spectrum = run.apply_jacob(True)[-1]
ax1.plot(exp_data[0,:], exp_data[2,:],'-', color="gray")
ax1.scatter(exp_data[0,:], exp_data[2,:], marker="x", s=16, label="Exp.")
ax1.plot(wl,spectrum,color="red", label="Sim.")
ax1.set_xticks([600,650,700,750,800,850,900])
ax1.text(0.05, 0.8,r"400$\mu$J")
ax1.tick_params(direction='in', width = 2, length = 8)
plt.legend(bbox_to_anchor=(0.62, 0.95), loc=2, borderaxespad=0.)

# Second Plot
run = fibs.runs[2]
wl = run.make_wavelength_scale()
spectrum = run.apply_jacob(True)[-1]
ax2.plot(exp_data[3,:], exp_data[5,:],'-', dashes=[8,2], color="gray")
ax2.scatter(exp_data[3,:], exp_data[5,:], marker="x", s=16)
ax2.plot(wl,spectrum,color="red")
ax2.set_xticks([600,650,700,750,800,850,900])
ax2.text(0.05, 0.8,r"800$\mu$J")
ax2.tick_params(direction='in', width = 2, length = 8)

# Third Plot
run = fibs.runs[0]
wl = run.make_wavelength_scale()
spectrum = run.apply_jacob(True)[-1]
ax3.plot(exp_data[6,:], exp_data[8,:],'-', dashes=[8,2], color="gray")
ax3.scatter(exp_data[6,:], exp_data[8,:], marker="o", s=16)
ax3.plot(wl,spectrum,color="red")
ax3.set_xticks([600,650,700,750,800,850,900])
ax3.text(0.05, 0.8,r"1600$\mu$J")
ax3.tick_params(direction='in', width = 2, length = 8)

ax3.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("Intensity (A.U)")
plt.show()
