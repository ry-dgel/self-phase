import os
import fiber_run

class fiber_set:
    def __init__(self, data_list):
        self.runs = data_list

    def max_metric(self, metric):
        maximum = max(map(lambda run : run.metrics[metric], self.runs))
        return maximum, next(run for run in self.runs if 
                            run.metrics[metric] == maximum)

    def min_metric(self, metric):
        minimum = min(map(lambda run : run.metrics[metric], self.runs))
        return minimum, next(run for run in self.runs if
                             run.metrics[metric] == minimum)

    def const_param(self, param, value):
        fibs = fiber_set([run for run in self.runs if 
                          run.params[param] == value])
        return fibs

    def list_params(self):
        return [run.params for run in self.runs]

    def list_param(self, param):
        return [run.params[param] for run in self.runs]

    def list_metrics(self):
        return [run.metrics for run in self.runs]

    def list_metric(self, metric):
        return [run.metrics[metric] for run in self.runs]

def load_data(folder):
    data_list = []
    for root, folders, files, in os.walk(folder):
        for folder in folders:
            print("Loading data from %s" % folder, end='\r')
            data_folder = os.path.join(root, folder)
            data_list.append(fiber_run.fiber_run(data_folder))

    return fiber_set(data_list)
