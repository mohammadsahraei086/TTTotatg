from coffea.analysis_tools import PackedSelection

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
        
    def select_good_events(self, channel="total"):
        selection = PackedSelection()
        selected_events = self.select_two_lep_events()
        
        selection.add("leadingLepPT", selected_events.GoodLeptons[:, 0].pt > 25)
        selection.add("OCLep", (selected_events.GoodLeptons[:, 0].charge + selected_events.GoodLeptons[:, 1].charge) == 0)
        selection.add("lepInvariantMass", (selected_events.GoodLeptons[:, 0] + selected_events.GoodLeptons[:, 1]).mass > 20)
        selection.add("onePhoton", selected_events.nGoodPhotons == 1)
        selection.add("atLeastOneBJet", selected_events.nGoodBJets >= 1)

        # Add selection for different channels
        selection.add("emu", selected_events.GoodLeptons.flavor[:, 0] != selected_events.GoodLeptons.flavor[:, 1])
        selection.add("ee", (selected_events.GoodLeptons.flavor[:, 0] == "e") & (selected_events.GoodLeptons.flavor[:, 1]=="e"))
        selection.add("mumu", (selected_events.GoodLeptons.flavor[:, 0] == "mu") & (selected_events.GoodLeptons.flavor[:, 1]=="mu"))
        
        if channel=="total":
            selected_events = selected_events[selection.all("leadingLepPT", "OCLep", "lepInvariantMass", "onePhoton", "atLeastOneBJet")]
        else:
            selected_events = selected_events[selection.all("leadingLepPT", "OCLep", "lepInvariantMass", "onePhoton", "atLeastOneBJet", f"{channel}")]
                            
        return selected_events