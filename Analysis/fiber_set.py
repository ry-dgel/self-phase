import os
import fiber_run

class fiber_set:
    def __init__(self, data_folder):
        self.runs = []
        for root, folders, files, in os.walk(data_folder):
            for folder in folders:
                data_folder = os.join.path(root, folder)
                self.runs.append(fiber_run(data_folder))

    def maxmetric(self, metric):
        maximum = max(map(lambda run : run.metrics[metric], self.runs))
        return next(run for run in self.runs if 
                    self.runs.metrics[metric] == maximum)

    def minmetric(self, metric):
        minimum = min(map(lambda run : run.metrics[metric], self.runs))
        return next(run for run in self.runs if
                    self.runs.metrics[metric] == minimum)

    def constparam(self, param, value):
        return [run for run in self.runs if self.runs.params[param] == value]
