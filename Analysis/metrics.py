import numpy as np
from numpy import fft as f 

def populate_metrics(fiber_run, fname):
    bw = bandwidth(fiber_run.fields[-1])
    pw = pulse_width(fiber_run.fields[-1])
    pr = power_ratio(fiber_run.fields[-1], fiber_run.fields[0])
    rho_max = plas_denzel(fname)
    metrics = {"pwidth":pw, "bwidth":bw[0], "l_edge":bw[1], "r_edge":bw[2],
                "rho_max":rho_max, "power":pr}
    return metrics

def plas_denzel(fname):
    plas_denzel = np.genfromtxt(fname + "/PlasmaDensity")
    rho_max = max(plas_denzel)
    return rho_max	

def bandwidth(E):
    Ef = f.fftshift(f.fft(f.fftshift(E)))
    If = np.power(np.abs(Ef),2)
    tm = 0.1*max(If)
    right_edge = If[len(If)-np.argmax(np.flip(If,0) > tm)]
    left_edge = If[np.argmax(If > tm)]
    bw = right_edge - left_edge
    return [bw, left_edge, right_edge]

def pulse_width(E):
    I = np.power(np.abs(E),2)
    tm = 0.1*max(I)
    pw = I[len(I) - np.argmax(np.flip(I,0) > tm)] - I[np.argmax(I > tm)]
    return pw

def power_ratio(E, E0):
    I = np.power(np.abs(E),2)
    I0 = np.power(np.abs(E0),2)
    P0 = np.trapz(I0)
    Pf = np.trapz(I)
    return Pf/P0
