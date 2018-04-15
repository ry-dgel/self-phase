import numpy as np
from numpy import fft as f 

cf = 6.626e-34 * 6.242e18

# Given a fiber_run, computes each metric putting it in fiber_run.metrics
def populate_metrics(fiber_run, fname):
    params = fiber_run.params
    bw = bandwidth(fiber_run)
    print(bw)
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

# Computes the FWTM bandwidth of the final field in nm
def bandwidth(fiber_run):
    wl = fiber_run.make_wavelength_scale() 
    spectrum = fiber_run.apply_jacob(True)[-1]

    # 0.1 since looking for FWTM of a normalized spectrum
    # made to be 0.0995 for a little leeway in rounding.
    lmin = min(np.where(spectrum >= 0.0995)[0])
    rmin = min(np.where(np.flip(spectrum,0) >= 0.0995)[0])
    
    # Probably over shot 0.1 by a bit, so take average with previous wl
    left_edge = np.mean(np.flip(wl,0)[rmin-1:rmin])
    right_edge = np.mean(wl[lmin-1:lmin])
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
