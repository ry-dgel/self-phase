import os
import fiber_run

class fiber_set:
    def __init__(self, data_list):
        self.runs = data_list
    
    # Determines which fiber_run maximizes a metric, returns the metric and
    # the fiber_run object
    def max_metric(self, metric):
        maximum = max(map(lambda run : run.metrics[metric], self.runs))
        return maximum, next(run for run in self.runs if 
                            run.metrics[metric] == maximum)

    # Determines which fiber_run minimizes a metric, returns the metric and
    # the fiber_run object
    def min_metric(self, metric):
        minimum = min(map(lambda run : run.metrics[metric], self.runs))
        return minimum, next(run for run in self.runs if
                             run.metrics[metric] == minimum)
    
    # Returns a fiber_set object where each fiber_run has $value as $param
    def const_param(self, param, value):
        fibs = fiber_set([run for run in self.runs if 
                          run.params[param] == value])
        return fibs
    
    # Returns a list of the parameter dictionaries of each fiber_run
    def list_params(self):
        return [run.params for run in self.runs]
    
    # Returns a list of each value of a specific parameter for each fiber_run
    def list_param(self, param):
        return [run.params[param] for run in self.runs]
    
    # Returns a list of the metric dictionaries of each fiber_run
    def list_metrics(self):
        return [run.metrics for run in self.runs]

    # Returns a list of each value of a specific metric for each fiber_run
    def list_metric(self, metric):
        return [run.metrics[metric] for run in self.runs]

# Returns a fiber_set where each run is taken from a data folder.
def load_data(folder):
    data_list = []
    for root, folders, files, in os.walk(folder):
        for folder in folders:
            print("Loading data from %s" % folder, end='\r')
            data_folder = os.path.join(root, folder)
            data_list.append(fiber_run.fiber_run(data_folder))

    return fiber_set(data_list)
