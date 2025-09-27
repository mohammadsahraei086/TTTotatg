import time

from coffea import processor
from coffea.util import save
from coffea.nanoevents import DelphesSchema

from hist_manager import HistManager
from object_selector import ObjectSelector
from event_selector import EventSelector
from fileset import *

fileset = fileset_pc_limit

class TTPairTotatg(processor.ProcessorABC):
    
    def __init__(self):
        self.hist_manager = HistManager()
        self.hist_manager.define_axes()
        self.hist_manager.define_histograms()
        self.histograms = self.hist_manager.get_histograms()
        self.categories = ["total", "emu", "ee", "mumu"]
    
    def define_output_layout(self):
        masses = set([])
        for item, dic in fileset.items():
            masses.add(dic["metadata"]["mass"])
        output = {}
        output["metadata"] = {}
        output["nEvents"] = {}
        output["nEvents"]["primary"] = {}
        output["nEvents"]["selected"] = {}
        output["hists"] = {}
        for mass in masses:
            output["metadata"][mass] = {}
        for cat in self.categories:
            output["nEvents"]["selected"][cat] = {}
            output["hists"][cat] = {}
            for hist in self.histograms:
                output["hists"][cat][hist] = {}
                for mass in masses:
                    output["nEvents"]["selected"][cat][mass] = {}
                    output["hists"][cat][hist][mass] = {}
                    
                    
        return output

    def process(self, events):
        dataset = events.metadata["dataset"]
        mass = events.metadata["mass"]
        self.output = self.define_output_layout()
        self.events = events
        self.output["metadata"][mass][dataset] = events.metadata
        self.output["nEvents"]["primary"][dataset] = len(self.events)
        self.events["n_primary"] = len(self.events)

        object_selector = ObjectSelector(self.events)
        object_selector.select_good_objects()
        object_selector.count_good_objects()
        
        event_selector = EventSelector(self.events)

        for cat in self.categories:
            selected_events = event_selector.select_good_events(cat)
            self.output["nEvents"]["selected"][cat][mass][dataset] = len(selected_events)
            for name, hist in self.histograms.items():
                hist.fill(selected_events)
                self.output["hists"][cat][name][mass][dataset] = hist.get_histogram()
        
        return self.output

    def postprocess(self, accumulator):
        pass


#####################################################################################################################
def main():
    # client = Client()

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
    save(out, 'output/output.coffea')
    
    elapsed = time.time() - tstart
    print(elapsed)

if __name__ == "__main__":
    main()