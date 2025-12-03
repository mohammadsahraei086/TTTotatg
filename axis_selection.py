import awkward as ak

def get_jet_multiplicity(events):
    n_jets = ak.num(events.GoodJets)
    n_bjets = ak.num(events.GoodBJets)
    
    categories = ak.Array(["other"] * len(events))
    
    categories = ak.where((n_jets == 1) & (n_bjets == 1), "1,1b", categories)
    categories = ak.where((n_jets == 2) & (n_bjets == 1), "2,1b", categories)
    categories = ak.where((n_jets >= 3) & (n_bjets == 1), ">=3,1b", categories)
    categories = ak.where((n_jets == 2) & (n_bjets == 2), "2,2b", categories)
    categories = ak.where((n_jets >= 3) & (n_bjets == 2), ">=3,2b", categories)
    categories = ak.where((n_jets >= 3) & (n_bjets >= 3), ">=3,>=3b", categories)  # etc.
    
    return ak.to_numpy(categories)