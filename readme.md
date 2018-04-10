# Self-Phase Modulation
## Code for simulation to be translated from matlab -> julia
Done by Rigel and Joseph as part of undergraduate thesis.

```
julia -p <number of sims to run at once/cores to use> multi-thread.jl <Parameter File> <Number of times to save data>
```
## Parameters
Parameters are input as either single numbers, or a triplet defining a range.
In both cases, the parameter must be surrounded by square brackets '[]'.
Ranges are given as [<min value>, <max value>, <step size>].

All parameters in the following list must be defined, and no others are available:
* Energy  -  The total energy of the laser pulse, in Joules.
* Tfwhm   -  The full width half maximum pulse duration, in seconds.
* lambda  -  The central wavelength of the pulse, in meters.
* Pin     -  Entrance pressure in the fiber.
* Pout    -  Exit pressure in the fiber; set to equal to Pin for static pressure.
* Chirp   -  Chirp value of pulse.
* fiberD  -  Diameter of optical fiber.
* dz      -  The spatial step of the simulation, in seconds.
* zmax    -  The simulation length (i.e. the length of the fiber), in meters.
* Nt      -  Number of points in data grid.
* tmax    -  The max time before/after the central peak of the pulse to simulate; defines max values for spatial and Fourier space.
View params.yaml for an example.

## Analysis
We've provided some Python code for analyzing the data sets produced from multiple simulations.
These are locaded in the `Analysis` folder.

As is, the easiest way to import and use it all is to include the following in a script:
```
from plot import *
```

with a single directory containing one subfolder for each simulation run.
The following loads the entire data set:
```
fibs = fs.load(<path/to/data>)
```
This will load in the data set, and for each one, load in or compute some useful 
metrics:
* pwidth - Full width half maximum of pulse, in femtoseconds.
* bwidth - Tenth width half maximum of pulse spectrum, in Hertz.
* l_edge - Leftmost edge used to calculate bwidth, in Hertz. 
* r_edge - Rightmost edge used to calculate bwidth, in Hertz.
* rho_max - Maximum plasma density during simulation, in ions/cubic centimeter.
* power - Ratio of final and initial power.

## Todo
* ~~Translate Main Loop~~
* ~~Translate "PPT.m"~~
* ~~Turn chunks of main loop into functions~~
* Optimize (This will never be crossed out)
* ~~Test~~
* ~~Use plan fft/ifft to speed up all the transforms~~
* ~~Multi-Thread drifting~~
