import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fiber_set as fs

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
    points = np.arange(-Nt/2, Nt/2, 1)
    f = np.concatenate([points/tmax + 299792458/run.params["lambda"] 
                        for run in fiber_set.runs])

    dic = {} 
    # Concat fields together to make intensity column
    I = np.concatenate(list(map(lambda spectrum: np.power(np.abs(spectrum),2),
                           [run.spectra()[-1] for run in fiber_set.runs])))
    dic.update({"f": f, "I" : I})

    #Make extra columns given by plist and mlist.
    dic.update({param : np.concatenate([np.repeat(run.params[param], len(run.fields[-1])) 
                                   for run in fiber_set.runs]) 
                for param in plist})
    dic.update({metric : np.concatenate([np.repeat(run.metrics[metric], len(run.fields[-1]))
                                    for run in fiber_set.runs])
                for metric in mlist})
    return pd.DataFrame(dic)

def frame(plist, mlist):
    dic = {}

    #Make columns given by plist and mlist.
    dic.update({param : np.concatenate([np.repeat(run.params[param], len(run.fields[-1])) 
                                   for run in fiber_set.runs]) 
                for param in plist})
    dic.update({metric : np.concatenate([np.repeat(run.metrics[metric], len(run.fields[-1]))
                                    for run in fiber_set.runs])
                for metric in mlist})

    return pd.DataFrame(dic)

def plot_spect_grid(fiber_set, p1, p2=None):
    df = frame_spectra(fiber_set, plist=[p1,p2])
    if p2:
        g = sns.FacetGrid(df, row=p1, aspect=2)
    else:
        g = sns.FacetGrid(df, row=p1, col=p2, aspect=2)
    g.map(plt.plot, "f", "I")
    return g
