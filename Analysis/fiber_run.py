import yaml
from itertools import groupby
import numpy as np
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

        # Load paramters
        with open(fname + "/params") as f:
            self.params = yaml.load(f)
            for key, val in self.params.items():
                self.params[key] = float(val)

        self.metrics = metrics.populate_metrics(self, fname)
