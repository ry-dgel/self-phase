import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fiber_set as fs

hbar_evpj = 6.626e-34 * 6.242e18  # Convert Hz to eV
hc = 1239.84193 # Convert eV to nm

def frame_spectra(fiber_set, plist=[], mlist=[]):
    wl = fiber_set.runs[0].make_wavelength_scale()
    wl_length = len(wl)
    wl = np.tile(wl, len(fiber_set.runs))
    I = np.concatenate([run.apply_jacob(normed = True)[-1] for run in fiber_set.runs])
    dic = {"wl": wl, "I" : I}

    #Make extra columns given by plist and mlist.
    dic.update({param : np.concatenate([np.repeat(run.params[param], wl_length) 
                                   for run in fiber_set.runs]) 
                for param in plist})
    dic.update({metric : np.concatenate([np.repeat(run.metrics[metric], wl_length)
                                    for run in fiber_set.runs])
                for metric in mlist})

    return pd.DataFrame(dic)

def frame_full_spectra(fiber_run):
    wl = fiber_run.make_wavelength_scale()
    Is = fiber_run.apply_jacob(normed = True)
    z = np.linspace(0, np.run.params["zmax"], len(fiber_run.fields))
    dic = {"wl" : np.tile(wl, len(fiber_run)), "I" : np.concatenate(Is)}
    return pd.DataFrame(dic)


def frame_pulse(fiber_set, plist, mlist):
    t = np.tile(fiber_set.runs[0].make_time_scale(), len(fiber.set.runs))
    I = np.concatenate([np.power(np.abs(run.fields[-1]),2) for run in fiber_set.runs])
    dic = {"t": t, "I" : I}

    #Make columns given by plist and mlist.
    dic.update({param : [run.params[param] for run in fiber_set.runs]
                for param in plist})
    dic.update({metric : [run.metrics[metric] for run in fiber_set.runs]
                for metric in mlist})

    return pd.DataFrame(dic)

def frame(fiber_set, plist, mlist):
    dic = {}
    #Make columns given by plist and mlist.
    dic.update({param : [run.params[param] for run in fiber_set.runs]
                for param in plist})
    dic.update({metric : [run.metrics[metric] for run in fiber_set.runs]
                for metric in mlist})
    return pd.DataFrame(dic)

def plot_spect_grid(fiber_set, p1, p2=None):
    if p2:
        df = frame_spectra(fiber_set, plist=[p1,p2], mlist=["l_edge", "r_edge"])
        g = sns.FacetGrid(df, row=p1, col=p2, aspect=1)
    else:
        df = frame_spectra(fiber_set, plist=[p1], mlist=["l_edge", "r_edge"])
        g = sns.FacetGrid(df, row=p1, aspect=1)
    xmin = 0.95 * df["l_edge"].min(axis=0)
    xmax = 1.05 * df["r_edge"].max(axis=0)
    ymin = 0
    ymax = 1.1
    g.map(plt.plot, "wl", "I")
    g.set(xlim=(xmin, xmax))
    g.set(ylim=(ymin, ymax))
    return g

def plot_full_spectrum(fiber_run):
    wl = fiber_run.make_wavelength_scale()
    Is = np.array(fiber_run.apply_jacob(normed = True))
    z = np.linspace(0, fiber_run.params["zmax"], len(fiber_run.fields))
    X, Y = np.meshgrid(z, wl)
    fig, ax = plt.subplots(1)
    ax.pcolormesh(X, Y, Is.T, figure=fig, cmap="inferno")
    ax.set_ylim(fiber_run.metrics["l_edge"] * 0.95,fiber_run.metrics["r_edge"] * 1.05)
    ax.set_xlabel("Propagation Length (m)")
    ax.set_ylabel("Wavelength (nm)")
    return ax

def metric_heatmap(fiber_set, p1, p2, metric):
    df = frame(fiber_set, plist=[p1,p2], mlist=[metric]).pivot(p1, p2, metric)
    return sns.heatmap(df)
