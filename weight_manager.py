import numpy as np

class WeightManager:
    def __init__(self):
        self.run2_luminosity = 138
    
    def get_weights(self, events, *weights, **kwargs):
        total_weight = np.prod([])
        for weight in weights:
            new_weight = getattr(self, weight)
            total_weight = total_weight * new_weight(events, **kwargs)
        return total_weight
    
    def xsec(self, events):
        return events.metadata["xsec"]
    
    def luminosity(self, events):
        return self.run2_luminosity
    
    def sum_genweight(self, events):
        print(1./events["n_primary"])
        return 1./events["n_primary"]