import numpy as np
from numpy import fft as f 

def populate_metrics(fiber_run, fname):
	bandwidth = bandwidth(fiber_run.fields[-1])
	pulse_width = pulse_width(fiber_run.fields[-1])
	power_ratio = power_ratio(fiber_run.fields[-1], fiber_run.fields[0])
	rho_max = plas_denzel(fname)
	metrics = {"pwidth":pulse_width, "bwidth":bandwidth[0], "l_edge":bandwidth[1], "r_edge":bandwidth[2],
		"rho_max":rho_max, "power":power_ratio}

def plas_denzel(fname):
	plas_denzel = open(fname + "/PlasmaDensity").read().split("\n")
	for x in plas_denzel:
		x = float(x)
	rho_max = max(plas_denzel)
	return rho_max	

def bandwidth(E):
	Ef = f.fftshift(f.fft(f.fft(E)))
	re = np.real(Ef)
	im = np.im(Ef)
	If = re**2 + im**2
	tm = 0.1*max(If)
	right_edge = If[np.argmax(np.flip(If,0) > tm)]
	left_edge = If[np.argmax(If > tm)]
	bandwidth = right_edge - left_edge
	return [bandwidth, left_edge, right_edge]

def pulse_width(E):
	re = np.real(E)
	im = np.im(E)
	I = re**2 + im**2
	tm = 0.1*max(I)
	pulse_width = I[np.argmax(np.flip(I,0) > tm)] - I[np.argmax(I > tm)]
	return pulse_width

def power_ratio(E, E0):
	re = np.real(E)
	im = np.im(E)
	I = re**2 + im**2
	re = np.real(E0)
	im = np.im(E0)
	I0 = re**2 + im**2
	P0 = np.trapz(I0)
	Pf = np.trapz(I)
	return Pf/P0