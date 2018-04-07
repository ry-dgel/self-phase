import yaml
from itertools import groupby
import numpy as np
import numpy.fft as fft
import metrics

class fiber_run:
    def __init__(self, fname):
        self.fields  = []
        self.params  = {}
        self.metrics = {}

        # Read in fields from fname/E
        data_list = []
        with open(fname + "/E") as f:
            for k, g in groupby(f, lambda x: x=="\n"):
                if not k:
                    data_list.append(np.array([[float(x) for x in d.split(',')]
                                                for d in g if len(d.strip())]))

        # Turn raw field data in numpy arrays of complex numbers
        for data in data_list:
            re, im = data.T
            self.fields.append(re + 1j*im)

        # Load parameters
        with open(fname + "/params") as f:
            self.params = yaml.load(f)
            for key, val in self.params.items():
                self.params[key] = float(val)

        self.metrics = metrics.populate_metrics(self, fname)
        self.fix_units()
    
    # Returns a list of the spectra corresponding to each field
    def spectra(self):
        return [norm_spectra(field) for field in self.fields]

    def fix_units(self):
        self.params['lambda'] = self.params['lambda']*1e9
        self.params['Energy'] = self.params['Energy']*1e6
        self.params['Tfwhm'] = self.params['Tfwhm']*1e15
        self.params['tmax'] = self.params['tmax']*1e15

def norm_spectra(field):
    transform = np.power(np.abs(fft.fftshift(fft.fft(fft.fftshift(field)))),2)
    return transform/np.max(transform)

