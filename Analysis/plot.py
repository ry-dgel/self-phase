import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import fiber_set as fs

font = {'size'   : 28}

matplotlib.rc('font', **font)

hbar_evpj = 6.626e-34 * 6.242e18  # Convert Hz to eV
hc = 1239.84193 # Convert eV to nm

def frame_spectra(fiber_set, plist=[], mlist=[]):
    wls = []
    Is = []
    params = {param : [] for param in plist}
    metrics = {metric : [] for metric in mlist}
    for run in fiber_set.runs:
        wl = run.make_wavelength_scale()
        peak_ind = np.argmin(np.abs(wl - run.params["lambda"]))
        wls.append(wl)
        I = run.apply_jacob(normed = True)[-1]
#        J = np.zeros(np.size(I))
#        for i,_ in enumerate(I):
#            J[i] = I[2*peak_ind - i]
        Is.append(I)
        for param in plist:
            params[param] = np.concatenate([params[param],np.repeat(run.params[param], len(wl))])
        for metric in mlist:
            metrics[metric] = np.concatenate([metrics[metric],np.repeat(run.metrics[metric], len(wl))])

    dic = {"wl" : np.concatenate(wls),
            "I" : np.concatenate(Is),
          }

    if plist:
        dic.update(params)
    if mlist:
        dic.update(metrics)
    return pd.DataFrame(dic)

def frame_full_spectra(fiber_run):
    wl = fiber_run.make_wavelength_scale()
    Is = fiber_run.apply_jacob(normed = True)
    dic = {"wl" : np.tile(wl, len(fiber_run)), "I" : np.concatenate(Is)}
    return pd.DataFrame(dic)


def frame_pulse(fiber_set, plist, mlist):
    t = np.tile(fiber_set.runs[0].make_time_scale(), len(fiber_set.runs))
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

def plot_before_final(fiber_run):
    t = fiber_run.make_time_scale()
    wl = fiber_run.make_wavelength_scale()

    fig, (ax1,ax2) = plt.subplots(1, 2, sharey=True)

    # Time Plot
    Ii = np.power(np.abs(fiber_run.fields[0]),2)
    Ii = Ii/max(Ii)
    If = np.power(np.abs(fiber_run.fields[-1]),2)
    If = If/max(If)
    ax1.plot(t, Ii)
    ax1.plot(t, If)
    tmin = -fiber_run.metrics["l_edge"] / 2 * 1.05
    tmax = fiber_run.metrics["l_edge"] / 2 * 1.05
    ax1.set_xlim(tmin,tmax)
    ax1.set_xlabel("Time (fs)")
    ax1.set_ylabel("Intensity (A.U.)")
    
    # Spectrum Plot
    spectra = fiber_run.apply_jacob(True)
    Fi = spectra[0]
    Ff = spectra[-1]
    ax2.plot(wl, Fi, label="Initial")
    ax2.plot(wl, Ff, label="Final")
    wlmin = 0.95 * fiber_run.metrics["l_edge"]
    wlmax = 1.05 * fiber_run.metrics["r_edge"]
    ax2.set_xlim(wlmin, wlmax)
    ax2.set_ylim(0,1.05)
    ax2.set_xlabel("Wavelength (nm)")
    plt.legend(bbox_to_anchor=(0.62, 0.95), loc=2, borderaxespad=0.)
    return fig

