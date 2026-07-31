from coffea.analysis_tools import PackedSelection
import awkward as ak

class EventSelector:
    def __init__(self, events):
        self.events = events

    def add_trigger_selection(self):
        pass

    def primary_skim(self):
        self.add_trigger_selection()

        return self._selection("trigger")
    
    def select_two_lep_events(self):
        selection = PackedSelection()
        selection.add("twoLep", self.events.nGoodLeptons == 2)
        twolep_events = self.events[selection.require(twoLep=True)]
        
        return twolep_events

    def define_variables_before_selection(self, events):
        
        t = events.Particle[(abs(events.Particle.PID)==9000017) & (events.Particle.Status==22)]
        events["TTbarMass"] = (t[:, 0]+t[:, 1]).mass

    def select_good_events(self, channel="total"):
        selection = PackedSelection()
        cutflow = {}
        weight = (self.events.metadata["xsec"]*1000 * 138)/self.events.metadata["nevents"]
        cutflow["primary"] = len(self.events)*weight

        if not "SM" in self.events.metadata["dataset"]:
            self.define_variables_before_selection(self.events)
            lambda_eff = 7071 # 1/sqrt(c_tg^2+c_tgamma^2)
            cutflow["eft_val_total"] = len(self.events[self.events.TTbarMass < lambda_eff])*weight

        selected_events = self.select_two_lep_events()
        cutflow["n_lep=2"] = len(selected_events)*weight

        selection.add("leadingLepPT", selected_events.GoodLeptons[:, 0].pt > 25)
        selection.add("OCLep", (selected_events.GoodLeptons[:, 0].charge + selected_events.GoodLeptons[:, 1].charge) == 0)
        selection.add("lepInvariantMass", (selected_events.GoodLeptons[:, 0] + selected_events.GoodLeptons[:, 1]).mass > 20)
        selection.add("onePhoton", selected_events.nGoodPhotons == 1)
        selection.add("atLeastOneBJet", selected_events.nGoodBJets >= 1)
        if not "SM" in self.events.metadata["dataset"]:
            selection.add("eft_val", selected_events.TTbarMass < lambda_eff)

        mask = ak.Array([True] * len(selected_events))
        for name in selection.names:
            new_mask = selection.all(name)
            cutflow[name] = ak.sum(mask & new_mask)*weight
            mask = mask & new_mask

        # Add selection for different channels
        selection.add("emu", selected_events.GoodLeptons.flavor[:, 0] != selected_events.GoodLeptons.flavor[:, 1])
        selection.add("ee", (selected_events.GoodLeptons.flavor[:, 0] == "e") & (selected_events.GoodLeptons.flavor[:, 1]=="e"))
        selection.add("mumu", (selected_events.GoodLeptons.flavor[:, 0] == "mu") & (selected_events.GoodLeptons.flavor[:, 1]=="mu"))
        
        if channel=="total":
            if not "SM" in self.events.metadata["dataset"]:
                selected_events = selected_events[selection.all("leadingLepPT", "OCLep", "lepInvariantMass", "onePhoton", "atLeastOneBJet", "eft_val")]
            else:
                selected_events = selected_events[selection.all("leadingLepPT", "OCLep", "lepInvariantMass", "onePhoton", "atLeastOneBJet")]
        else:
            if not "SM" in self.events.metadata["dataset"]:
                selected_events = selected_events[selection.all("leadingLepPT", "OCLep", "lepInvariantMass", "onePhoton", "atLeastOneBJet", f"{channel}", "eft_val")]
            else:
                selected_events = selected_events[selection.all("leadingLepPT", "OCLep", "lepInvariantMass", "onePhoton", "atLeastOneBJet", f"{channel}")]
            cutflow[channel] = ak.sum(mask & selection.all(channel))*weight
                            
        return selected_events, cutflow
