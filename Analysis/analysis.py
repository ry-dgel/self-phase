import os
import yaml
from itertools import groupby
import numpy as np

# Load in the data from a folder full of simulations
def load_data(folder):
    data_list = []
    for root, folders, files in os.walk(folder):
        for folder in folders:
            print("Loading data from %s" % folder, end='\r')
            data_folder = os.path.join(root, folder)
            data_list.append(read_folder(data_folder))
    return data_list
# Load in the data from an individual simulation run
def read_folder(folder):
    # Read in the parameters
    params = {}
    with open(os.path.join(folder, "params")) as f:
        params = yaml.load(f)
        for key, val in params.items():
            params[key] = float(val)
    # Read in the final field
    with open(os.path.join(folder, "E")) as f:
        data_list = []
        for k, g in groupby(f, lambda x: x=="\n"):
            if not k:
                data_list.append(np.array([[float(x) for x in d.split(',')]
                                            for d in g if len(d.strip())]))
        field = data_list[-1]
    # Read in the max plasma density
    with open(os.path.join(folder, "PlasmaDensity")) as f:
        maxDensity = max(np.genfromtxt(f))
    
    return params, field, maxDensity
