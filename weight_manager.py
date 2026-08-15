import numpy as np

from weight_variations import WEIGHT_INDEX, NOMINAL_VARIATION

class WeightManager:
    def __init__(self):
        self.run2_luminosity = 138
    
    def get_weights(self, events, *weights, **kwargs):
        total_weight = np.prod([])
        for weight in weights:
            new_weight = getattr(self, weight)
            total_weight = total_weight * new_weight(events, **kwargs)
        return total_weight
    
    def xsec(self, events, **kwargs):
        return events.metadata["xsec"]*1000
    
    def luminosity(self, events, **kwargs):
        return self.run2_luminosity    
    def sum_genweight(self, events, **kwargs):
        return 1./events.metadata["nevents"]

    def variation(self, events, variation=NOMINAL_VARIATION, **kwargs):
        """Per-event scale/PDF systematic weight, expressed as a ratio to
        the nominal (MUF=1, MUR=1, PDF=247000) weight.

        Multiplying this ratio into the usual xsec*luminosity/sum_genweight
        normalization reweights both the *shape* and the *rate* of the
        distribution for the requested variation: since the LHE sample is
        unweighted with event_norm=average, average(w_variation)/average(w_nominal)
        over the full sample equals xsec_variation/xsec_nominal, so no
        separate per-variation cross section needs to be looked up.

        `variation` must be one of the keys in weight_variations.WEIGHT_INDEX.
        Defaults to the nominal weight itself (ratio == 1) if not given.
        """
        w = events.Weight.Weight
        idx = WEIGHT_INDEX[variation]
        nominal_idx = WEIGHT_INDEX[NOMINAL_VARIATION]
        return w[:, idx] / w[:, nominal_idx]
