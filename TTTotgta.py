import time
import copy

from coffea import processor
from coffea.util import save
from coffea.nanoevents import DelphesSchema

from hist_manager import HistManager
from object_selector import ObjectSelector
from event_selector import EventSelector
from fileset import *

fileset = fileset

class TTPairTotatg(processor.ProcessorABC):
    
    def __init__(self):
        self.hist_manager = HistManager()
        self.hist_manager.define_axes()
        self.hist_manager.define_histograms()
        self.histograms = self.hist_manager.get_histograms()
        self.categories = ["total"] # , "emu", "ee", "mumu"
    
    def define_output_layout(self):
        output = {}
        output["metadata"] = {}
        output["cutflow"] = {}
        output["hists"] = {}
        output["MTT_array"] = {}
        
        for cat in self.categories:
            output["cutflow"][cat] = {}
            output["hists"][cat] = {}
            for hist in self.histograms:
                output["hists"][cat][hist] = {}
                
        return output

    def process(self, events):
        dataset = events.metadata["dataset"]
        self.output = self.define_output_layout()
        self.events = events

        object_selector = ObjectSelector(self.events)
        object_selector.select_good_objects()
        object_selector.count_good_objects()
        
        event_selector = EventSelector(self.events)

        for cat in self.categories:
            selected_events, cutflow, MTT_array = event_selector.select_good_events(cat)
            self.output["cutflow"][cat][dataset] = cutflow
            self.output["MTT_array"][dataset] = MTT_array
            if  len(selected_events) == 0:
                continue
            for name, hist in self.histograms.items():
                # hist_copy = copy.deepcopy(hist)
                hist.fill(selected_events)
                self.output["hists"][cat][name][dataset] = copy.deepcopy(hist.get_histogram())
                hist.reset_histogram()
        
        return self.output

    def postprocess(self, accumulator):
        for dataset, value in fileset.items():
            accumulator["metadata"][dataset] = value["metadata"]


#####################################################################################################################
def main():
    # client = Client()

    tstart = time.time()
    
#     futures_run = processor.Runner(
#         executor = processor.FuturesExecutor(compression=None, workers=4),
#         schema=DelphesSchema,
#         maxchunks=10,
#         chunksize=10000
#     )

#     out = futures_run(
#         fileset,
#         "Delphes",
#         processor_instance=TTPairTotatg() 
#     )
    
    iterative_run = processor.Runner(
       executor = processor.IterativeExecutor(compression=None),
       schema=DelphesSchema,
       chunksize=10000,
        maxchunks=None,
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
