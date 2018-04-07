import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fiber_set as fs

cf = 6.626e-34 * 6.242e18  #converts frequency in Hz to energy in eV using h

def frame_pulse(fiber_set, plist=[], mlist=[]):
    # Make time axis
    Nt = fiber_set.runs[0].params["Nt"]
    tmax = fiber_set.runs[0].params["tmax"] 
    dt = tmax/Nt
    points = np.arange(-Nt/2, Nt/2, 1)
    t = points*dt

    dic = {} 
    # Repeat time axis for every run to make time column
    t = np.repeat(t, len(fiber_set.runs))
    # Concat fields together to make intensity column
    I = np.concatenate(list(map(lambda field: np.power(np.abs(field),2),
                           [run.fields[-1] for run in fiber_set.runs])))
    dic.update({"t": t, "I" : I})

    #Make extra columns given by plist and mlist.
    dic.update({param : np.concatenate([np.repeat(run.params[param], len(run.fields[-1])) 
                                   for run in fiber_set.runs]) 
                for param in plist})
    dic.update({metric : np.concatenate([np.repeat(run.metrics[metric], len(run.fields[-1]))
                                    for run in fiber_set.runs])
                for metric in mlist})
    return pd.DataFrame(dic)


def frame_spectra(fiber_set, plist=[], mlist=[]):
    # Make freq axis
    Nt = fiber_set.runs[0].params["Nt"]
    tmax = fiber_set.runs[0].params["tmax"]
    points = np.arange(-Nt/2, Nt/2+1, 1)
    f = np.concatenate([points/tmax + 299792458/run.params["lambda"] 
                        for run in fiber_set.runs])

    dic = {} 
    # Concat fields together to make intensity column
    I = np.concatenate(list(map(lambda spectrum: np.power(np.abs(spectrum),2),
                           [cf*run.spectra()[-1] for run in fiber_set.runs])))
    dic.update({"f": f, "I" : I})

    #Make extra columns given by plist and mlist.
    dic.update({param : np.concatenate([np.repeat(run.params[param], len(run.fields[-1])) 
                                   for run in fiber_set.runs]) 
                for param in plist})
    dic.update({metric : np.concatenate([np.repeat(run.metrics[metric], len(run.fields[-1]))
                                    for run in fiber_set.runs])
                for metric in mlist})
    return pd.DataFrame(dic)

def frame_full_spectra(fiber_set, plist=[], mlist=[]):
    # Make freq axis
    Nt = fiber_set.runs[0].params["Nt"]
    tmax = fiber_set.runs[0].params["tmax"]
    points = np.arange(-Nt/2, Nt/2, 1)

    f = np.concatenate([np.repeat(points/tmax + 299792458/run.params["lambda"],
                                  len(run.fields)) 
                        for run in fiber_set.runs])

    I = np.concatenate([np.concatenate([np.power(np.abs(spectrum), 2) 
                                        for spectrum in run.spectra()]) 
                        for run in fiber_set.runs])

    dic = {"f" : f, "I" : I}

    #Make columns given by plist and mlist.
    dic.update({param : np.concatenate([np.repeat(run.params[param], 
                                                  len(run.fields[-1] * len(run.fields))) 
                                   for run in fiber_set.runs]) 
                for param in plist})

    dic.update({metric : np.concatenate([np.repeat(run.metrics[metric], 
                                                   len(run.fields[-1] * len(run.fields)))
                                    for run in fiber_set.runs])
                for metric in mlist})

    return dic

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
    ymax = 1.1 * df["I"].max(axis=0)
    g.map(plt.plot, "f", "I")
    g.set(xlim=(xmin, xmax))
    g.set(ylim=(ymin, ymax))
    return g

def metric_heatmap(fiber_set, p1, p2, metric):
    df = frame(fiber_set, plist=[p1,p2], mlist=[metric]).pivot(p1, p2, metric)
    return sns.heatmap(df)

