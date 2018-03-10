import numpy as np
from numpy import fft as f 

# Given a fiber_run, computes each metric putting it in fiber_run.metrics
def populate_metrics(fiber_run, fname):
    params = fiber_run.params
    bw = bandwidth(params, fiber_run.fields[-1])
    pw = pulse_width(params, fiber_run.fields[-1])
    pr = power_ratio(fiber_run.fields[-1], fiber_run.fields[0])
    rho_max = plas_denzel(fname)
    metrics = {"pwidth":pw, "bwidth":bw[0], "l_edge":bw[1], "r_edge":bw[2],
                "rho_max":rho_max, "power":pr}
    return metrics

# Gets the maximum plasma density from a data set
def plas_denzel(fname):
    plas_denzel = np.genfromtxt(fname + "/PlasmaDensity")
    rho_max = max(plas_denzel)
    return rho_max	

# Computes the FWTM bandwidth of the final field
def bandwidth(params, E):
    Ef = f.fftshift(f.fft(f.fftshift(E)))
    If = np.power(np.abs(Ef),2)

    freq = np.arange(-params["Nt"]/2, params["Nt"]/2, 1)/params["tmax"]
    freq = freq + 299792458/params["lambda"]

    tm = 0.1*max(If)
    right_edge = freq[len(If)-np.argmax(np.flip(If,0) > tm)]
    left_edge = freq[np.argmax(If > tm)]
    bw = right_edge - left_edge
    return [bw, left_edge, right_edge]

# Computes the FWHM pulse width of the final field
def pulse_width(params, E):
    I = np.power(np.abs(E),2)
    
    dt   = params["tmax"]/params["Nt"]
    time = np.arange(-params["Nt"]/2, params["Nt"]/2, 1) * dt

    tm = 0.5*max(I)
    pw = time[len(I) - np.argmax(np.flip(I,0) > tm)] - time[np.argmax(I > tm)]
    return pw

# Computes the power ratio of the initial and final spectra
def power_ratio(E, E0):
    I = np.power(np.abs(E),2)
    I0 = np.power(np.abs(E0),2)
    P0 = np.trapz(I0)
    Pf = np.trapz(I)
    return Pf/P0
