import seaborn as sns
import pandas as pd
import fiber_set as fs
import numpy as np
import matplotlib.pyplot as plt

def two_param_frame(fiber_set, p1, p2):
    col_labels = set(fiber_set.list_param(p1))
    row_labels = set(fiber_set.const_param(p1, next(iter(col_labels))).list_param(p2))
    
    cols = {str(col) : 
            {str(row) : 
             fiber_set.const_param(p1, col).const_param(p2, row).runs[0]
             for row in row_labels}
            for col in col_labels}

    return pd.DataFrame(cols)

def plot_grid(fiber_set, p1, p2):
    grid = sns.FacetGrid(two_param_frame(fiber_set, p1, p2))
    grid.map(lambda run: plt.plot(np.abs(run.fields[-1])**2), "Time", "Intensity")
    grid.fig.tight_layout(w_pad=1)
