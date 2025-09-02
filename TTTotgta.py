import time
import numpy as np
import awkward as ak
import os

from coffea import processor
from coffea.nanoevents import DelphesSchema

import coffea

# from distributed import Client

class TTPairTotatg(processor.ProcessorABC):
    
    def __init__(self):
        self.categories = ["total", "emu", "ee", "mumu"]
        self.variables = ["pt", "eta"]
        self.luminosity = 138e3
    
    def define_output_layout(self, dataset):
        self.output = {}
        self.output["nEvents"] = {}
        self.output["nEvents"]["primary"] = {}
        self.output["nEvents"]["selected"] = {}
        self.output["hists"] = {}
        for cat in self.categories:
            self.output["nEvents"]["selected"][cat] = {}
            self.output["hists"][cat] = {}
            for var in self.variables:
                self.output["hists"][cat][var] = {}

    def process(self, events):
        dataset = events.metadata["dataset"]
        self.define_output_layout(dataset)
        self.events = events
        self.output["nEvents"]["primary"][dataset] = len(self.events)
        lum_weight = (events.metadata["xsec"]*self.luminosity)/len(self.events)

        object_selector = ObjectSelector(self.events)
        object_selector.select_good_objects()
        object_selector.count_good_objects()
        
        events_selector = EventSelector(self.events)
        hist_manager = HistManager()
        hist_manager.define_histograms()

        for cat in self.categories:
            selected_events = events_selector.select_good_events(cat)
            self.output["nEvents"]["selected"][cat][dataset] = len(selected_events)
            for var in self.variables:
                self.output["hists"][cat][var][dataset] = hist_manager.fill_histogram(selected_events, var, weight=lum_weight)
        
        return self.output

    def postprocess(self, accumulator):
        
        hist_manager = HistManager()
        hist_manager.define_histograms()
        for cat in self.categories:
            if cat == "total":
                continue
            hist_manager.plot_histograms(accumulator["hists"], channel=cat, signal=["Signal_500", "Signal_1000"])
            hist_manager.plot_histograms(accumulator["hists"], channel=cat, signal=["Signal_500", "Signal_1000"], normalize=True)


#####################################################################################################################
def main():
    # client = Client()

    fileset = {
        "Signal_500":{
            "files":{
                "/home/mohammad/Softwares/MG5_aMC_v3.6.3/MG5_aMC_v3_6_3/TTpairTotgta/TTpairTotgta/Events/run_02/tag_1_delphes_events.root":"Delphes"
                },
            "metadata":{
                "xsec": 3.959e-05
            }
        },
        "Signal_1000":{
            "files":{
                "/home/mohammad/Softwares/MG5_aMC_v3.6.3/MG5_aMC_v3_6_3/TTpairTotgta/TTpairTotgta/Events/run_03/tag_1_delphes_events.root":"Delphes"
            },
            "metadata":{
                "xsec": 6.409e-05
            }
        }
    }

#     dataset_runnable, dataset_updated = preprocess(
#         fileset,
#         align_clusters=False,
#         step_size=100_00,
#         files_per_batch=1,
#         save_form=False,
#     )

    tstart = time.time()
    
#     futures_run = processor.Runner(
#         executor = processor.FuturesExecutor(compression=None, workers=2),
#         schema=DelphesSchema,
#         maxchunks=10,
#     )

#     out = futures_run(
#         fileset,
#         "Delphes",
#         processor_instance=TTPairTotatg()
#     )
    
    iterative_run = processor.Runner(
        executor = processor.IterativeExecutor(compression=None),
        schema=DelphesSchema,
        maxchunks=10,
    )
    
    out = iterative_run(
        fileset,
        treename="Delphes",
        processor_instance=TTPairTotatg(),
    )
    print(out)
    
    # for i in out['Lepton']:
    #     print(i)
    # for i in out['arg']:
    #     print(i)
    
    elapsed = time.time() - tstart
    print(elapsed)

if __name__ == "__main__":
    main()