import time

from coffea import processor
from coffea.nanoevents import DelphesSchema

from hist_manager import HistManager
from object_selector import ObjectSelector
from event_selector import EventSelector

class TTPairTotatg(processor.ProcessorABC):
    
    def __init__(self):
        self.hist_manager = HistManager()
        self.hist_manager.define_axes()
        self.hist_manager.define_histograms()
        self.histograms = self.hist_manager.get_histograms()
        self.categories = ["total", "emu", "ee", "mumu"]
    
    def define_output_layout(self):
        output = {}
        output["nEvents"] = {}
        output["nEvents"]["primary"] = {}
        output["nEvents"]["selected"] = {}
        output["hists"] = {}
        for cat in self.categories:
            output["nEvents"]["selected"][cat] = {}
            output["hists"][cat] = {}
            for hist in self.histograms:
                output["hists"][cat][hist] = {}
        return output

    def process(self, events):
        dataset = events.metadata["dataset"]
        self.output = self.define_output_layout()
        self.events = events
        self.output["nEvents"]["primary"][dataset] = len(self.events)
        self.events["n_primary"] = len(self.events)

        object_selector = ObjectSelector(self.events)
        object_selector.select_good_objects()
        object_selector.count_good_objects()
        
        event_selector = EventSelector(self.events)

        for cat in self.categories:
            selected_events = event_selector.select_good_events(cat)
            self.output["nEvents"]["selected"][cat][dataset] = len(selected_events)
            for name, hist in self.histograms.items():
                self.output["hists"][cat][name][dataset] = hist.fill(selected_events)
        
        return self.output

    def postprocess(self, accumulator):
        pass
        # hist_manager = HistManager()
        # hist_manager.define_histograms()
        # for cat in self.categories:
        #     if cat == "total":
        #         continue
        #     hist_manager.plot_histograms(accumulator["hists"], channel=cat, signal=["Signal_500", "Signal_1000"])
        #     hist_manager.plot_histograms(accumulator["hists"], channel=cat, signal=["Signal_500", "Signal_1000"], normalize=True)


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
    
    elapsed = time.time() - tstart
    print(elapsed)

if __name__ == "__main__":
    main()