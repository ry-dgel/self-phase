import os
import fiber_run

class fiber_set:
    def __init__(self, data_folder):
        self.runs = []
        for root, folders, files, in os.walk(data_folder):
            for folder in folders:
                print("Loading data from %s" % folder, end='\r')
                data_folder = os.path.join(root, folder)
                self.runs.append(fiber_run.fiber_run(data_folder))

    def maxmetric(self, metric):
        maximum = max(map(lambda run : run.metrics[metric], self.runs))
        return maximum, next(run for run in self.runs if 
                            run.metrics[metric] == maximum)

    def minmetric(self, metric):
        minimum = min(map(lambda run : run.metrics[metric], self.runs))
        return minimum, next(run for run in self.runs if
                             run.metrics[metric] == minimum)

    def constparam(self, param, value):
        fibs = fiber_set()
        fibs.runs = [run for run in self.runs if 
                     self.runs.params[param] == value]
        return fibs
