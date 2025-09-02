import awkward as ak


class ObjectSelector:
    def __init__(self, events):
        self.events = events
        self.particles = self.events.Particle
        self.particle_stable = self.particles[self.particles.Status == 1]
        self.particle_electrons = self.particles[(abs(self.particles.PID) == 11) & (self.particles.Status == 1)]
        self.particle_muons = self.particles[(abs(self.particles.PID) == 13) & (self.particles.Status == 1)]
        self.particle_photons = self.particles[(self.particles.PID == 22) & (self.particles.Status == 1)]       
        
    def selected_electrons(self):
        
        broadcast_photons = ak.ones_like(self.particle_electrons.pt)[:,:, None] * self.particle_photons[:, None,:]
        # ele_pho = ak.cartesian({"electron": self.particle_electrons, "photon": self.particle_photons}, nested=True)
        # broadcast_photons = ele_pho.photon
        dr_mask = self.particle_electrons.metric_table(self.particle_photons) < 0.1
        selected_photons = broadcast_photons[dr_mask]
        reco_electrons = self.particle_electrons + ak.with_name(ak.sum(selected_photons, axis=2), name="PtEtaPhiMLorentzVector")
        reco_electrons = ak.with_field(reco_electrons, self.particle_electrons.Charge, "charge")
        reco_electrons = ak.with_field(reco_electrons, "e", "flavor")
        selected_electrons = reco_electrons[(reco_electrons.pt > 15) & (abs(reco_electrons.eta) < 2.4)]
        
        return selected_electrons
    
    def selected_muons(self):
        
        broadcast_photons = ak.ones_like(self.particle_muons.pt)[:,:, None] * self.particle_photons[:, None,:]
        # mu_pho = ak.cartesian({"muon": self.particle_muons, "photon": self.particle_photons}, nested=True)
        # broadcast_photons = mu_pho.photon
        dr_mask = self.particle_muons.metric_table(self.particle_photons) < 0.1
        selected_photons = broadcast_photons[dr_mask]
        reco_muons = self.particle_muons + ak.with_name(ak.sum(selected_photons, axis=2), name="PtEtaPhiMLorentzVector")
        reco_muons = ak.with_field(reco_muons, self.particle_muons.Charge, "charge")
        reco_muons = ak.with_field(reco_muons, "mu", "flavor")
        selected_muons = reco_muons[(reco_muons.pt > 15) & (abs(reco_muons.eta) < 2.4)]
        
        return selected_muons
    
    def selected_photons(self, leptons):
        
        selected_photons = self.particle_photons[(self.particle_photons.pt > 20) & (abs(self.particle_photons.eta) < 1.44)]
        stable_particles = self.particle_stable[self.particle_stable.pt > 5]
        dr_mask = ak.all(
            (selected_photons.metric_table(stable_particles) > 0.1) |
            (selected_photons.metric_table(stable_particles) == 0),
            axis=2
        )
        selected_photons = selected_photons[dr_mask]
        selected_photons = selected_photons[ak.prod(selected_photons.metric_table(leptons) > 0.4, axis=2) == 1]
        
        return selected_photons
    
    def selected_jets(self, leptons, photons):
        
        selected_jets = self.events.GenJet[(self.events.GenJet.pt > 30) & (abs(self.events.GenJet.eta) < 2.4)]
        selected_jets = selected_jets[ak.all(selected_jets.metric_table(leptons) > 0.4, axis=2)]
        selected_jets = selected_jets[ak.all(selected_jets.metric_table(photons) > 0.1, axis=2)]
        
        return selected_jets
    
    def selected_b_jets(self, jets):
        
        gen_b = self.particles[(self.particles.Status == 23) & (abs(self.particles.PID) == 5)]
        selected_b_jets = jets[ak.any(jets.metric_table(gen_b) < 0.4, axis=2)]
        b_argmin = ak.argmin(selected_b_jets.metric_table(gen_b), axis=2)
        selected_b_jets = selected_b_jets[((selected_b_jets.pt-gen_b[b_argmin].pt)/selected_b_jets.pt) < 1]
        
        return selected_b_jets
    
    def select_good_objects(self):
        
        self.events["GoodElectrons"] = self.selected_electrons()
        self.events["GoodMuons"] = self.selected_muons()
        self.events["GoodLeptons"] = ak.concatenate((self.events["GoodElectrons"], self.events["GoodMuons"]), axis=1)
            #name='PtEtaPhiMLorentzVector',
        # )
        arg = ak.argsort(self.events.GoodLeptons.pt, ascending=False)
        self.events["GoodLeptons"] = self.events["GoodLeptons"][arg]
        self.events["GoodPhotons"] = self.selected_photons(self.events.GoodLeptons)
        self.events["GoodJets"] = self.selected_jets(self.events.GoodLeptons, self.events.GoodPhotons)
        self.events["GoodBJets"] = self.selected_b_jets(self.events.GoodJets)
        
    def count_good_objects(self):
        
        self.events["nGoodElectrons"] = ak.num(self.events.GoodElectrons)
        self.events["nGoodMuons"] = ak.num(self.events.GoodMuons)
        self.events["nGoodLeptons"] = ak.num(self.events.GoodLeptons)
        self.events["nGoodPhotons"] = ak.num(self.events.GoodPhotons)
        self.events["nGoodJets"] = ak.num(self.events.GoodJets)
        self.events["nGoodBJets"] = ak.num(self.events.GoodBJets)