import numpy as np
import hist
import awkward as ak
from dataclasses import dataclass
from typing import Callable, List, Optional

from axis_selection import *

from weight_manager import WeightManager
from weight_variations import VARIATIONS

@dataclass
class Axis:
    name: str 
    label: str 
    bins: int = None
    start: float = None
    stop: float = None
    type: str = "regular"
    function: Optional[Callable] = None
    growth = True
    
    def __post_init__(self):
        if self.function is not None:
            self.get_variable = self.function
    
    def get_variable(self, events):
        raise NotImplementedError(f"Provide a function parameter when creating Axis {self.name}")

class Histogram:
    # Category label used for the un-varied histogram. Kept distinct from
    # the "MUF=1.0_MUR=1.0_PDF=247000" entry in VARIATIONS so the plain
    # nominal weight (no ratio applied) is always available as a cross-check
    # against that variation, which should come out identical.
    NOMINAL_LABEL = "nominal"

    def __init__(self, name, axes:List[Axis], weights=None, apply_variations=False):
        self.name = name
        self.axes = axes
        # variations only make sense for weighted histograms
        self.apply_variations = apply_variations and weights is not None
        hist_axis = []
        for axis in self.axes:
            hist_axis.append(self.get_hist_axis(axis))
        if self.apply_variations:
            hist_axis.append(
                hist.axis.StrCategory(
                    [self.NOMINAL_LABEL] + VARIATIONS,
                    name="variation",
                    label="Scale/PDF systematic variation",
                    growth=False,
                )
            )
        self.weights = weights
        self.histogram = hist.Hist(*hist_axis, name=name, storage="weight")

    def get_hist_axis(self, ax: Axis):
        if ax.name == None:
            ax.name = f"{ax.coll}.{ax.field}"
        if ax.type == "regular" and isinstance(ax.bins, list):
            ax.type = "variable"
        if ax.type == "regular":
            return hist.axis.Regular(
                name=ax.name,
                bins=ax.bins,
                start=ax.start,
                stop=ax.stop,
                label=ax.label,
            )
        elif ax.type == "variable":
            if not isinstance(ax.bins, list):
                raise ValueError(
                    "A list of bins edges is needed as 'bins' parameters for a type='variable' axis"
                )
            return hist.axis.Variable(
                ax.bins,
                name=ax.name,
                label=ax.label
            )
        elif ax.type == "int":
            return hist.axis.Integer(
                name=ax.name,
                start=ax.start,
                stop=ax.stop,
                label=ax.label
            )
        elif ax.type == "intcat":
            return hist.axis.IntCategory(
                ax.bins,
                name=ax.name,
                label=ax.label
            )
        elif ax.type == "strcat":
            return hist.axis.StrCategory(
                ax.bins, name=ax.name, label=ax.label, growth=ax.growth
            )

    def fill(self, events):
        ax = {}
        for axis in self.axes:
            ax[axis.name] = axis.get_variable(events)

        if self.weights is None:
            self.histogram.fill(**ax)
            return

        weight_manager = WeightManager()

        if not self.apply_variations:
            weight = weight_manager.get_weights(events, *self.weights)
            self.histogram.fill(**ax, weight=weight)
            return

        # nominal: usual xsec*luminosity/sum_genweight normalization, no
        # scale/PDF reweighting applied
        nominal_weight = weight_manager.get_weights(events, *self.weights)
        self.histogram.fill(**ax, variation=self.NOMINAL_LABEL, weight=nominal_weight)

        # each scale/PDF variation: same normalization, reweighted by the
        # ratio of that variation's weight to the nominal weight
        for variation_name in VARIATIONS:
            varied_weight = weight_manager.get_weights(
                events, *self.weights, "variation", variation=variation_name
            )
            self.histogram.fill(**ax, variation=variation_name, weight=varied_weight)

    def get_histogram(self):
        return self.histogram
    
    def reset_histogram(self):
        self.histogram.reset()

class HistManager:
    def __init__(self):
        self.axes = {}
        self.histograms = {}
        
    def define_axes(self):
        # self.add_axis("photon_pt",
        #               "$p_T(\gamma) (GeV)$",
        #               [ 20.,  35.,  50.,  70., 100., 130., 165., 200., 250., 300.],
        #               type = "variable",
        #               function = lambda events: ak.flatten(events.GoodPhotons.pt))
        self.add_axis("xsec_photon_pt",
                      "$p_T(\gamma) (GeV)$",
                      [ 20.,  35.,  50.,  70., 130., 200., 300.],
                      type = "variable",
                      function = lambda events: ak.flatten(events.GoodPhotons.pt))
        # self.add_axis("deltaeta_ll",
        #               "$|\Delta\eta(\ell\ell)|$",
        #               [0, 0.5, 1, 1.5, 2, 2.5, 3, 4.5],
        #               type = "variable",
        #               function = lambda events: abs(events.GoodLeptons[:, 0].eta-events.GoodLeptons[:, 1].eta))
        self.add_axis("deltaphi_ll",
                      "$\Delta \phi(\ell\ell)$",
                      [0, 0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8, 3.2],
                      type = "variable",
                      function = lambda events: abs(events.GoodLeptons[:, 0].phi-events.GoodLeptons[:, 1].phi))
        # self.add_axis("ptl1plusptl2",
        #               "$p_{T}(\ell_{1})+p_{T}(\ell_{2})$",
        #               [40, 70, 100, 140, 190, 250, 330, 500],
        #               type = "variable",
        #               function = lambda events: events.GoodLeptons[:, 0].pt+events.GoodLeptons[:, 1].pt)
        # self.add_axis("jet_multiplicity",
        #               "",
        #               ["1,1b", "2,1b", ">=3,1b", "2,2b", ">=3,2b", ">=3,>=3b"],
        #               type = "strcat",
        #               function = get_jet_multiplicity)
        
    def define_histograms(self):
        # self.add_histogram("photon_pt",
        #                    [self.axes["photon_pt"]],
        #                    ["xsec", "luminosity", "sum_genweight"],
        #                    apply_variations=True
        #                   )
        self.add_histogram("diff_xsec_photon_pt",
                           [self.axes["xsec_photon_pt"]],
                           ["xsec", "luminosity", "sum_genweight"],
                           apply_variations=True
                          )
        # self.add_histogram("deltaeta_ll",
        #                    [self.axes["deltaeta_ll"]],
        #                    ["xsec", "luminosity", "sum_genweight"],
        #                    apply_variations=True
        #                   )
        self.add_histogram("deltaphi_ll",
                           [self.axes["deltaphi_ll"]],
                           ["xsec", "luminosity", "sum_genweight"],
                           apply_variations=True
                          )
        # self.add_histogram("ptl1plusptl2",
        #                    [self.axes["ptl1plusptl2"]],
        #                    ["xsec", "luminosity", "sum_genweight"],
        #                    apply_variations=True
        #                   )
        # self.add_histogram("jets_multiplicity",
        #                    [self.axes["jet_multiplicity"]],
        #                    # ["xsec", "luminosity", "sum_genweight"]
        #                   )
        
    def add_axis(self,
                 name,
                 label: str,
                 bins: int = None,
                 start: float = None,
                 stop: float = None,
                 type: str = "regular",
                 function: Optional[Callable] = None
                ):
        self.axes[name] = Axis(name, label, bins, start, stop, type, function)
        
    def add_histogram(self,
                      name,
                      axes:List[str],
                      weights= None,
                      apply_variations = False
                     ):
        self.histograms[name] = Histogram(name, axes, weights, apply_variations)
        
    def get_histogram(self, name):
        return self.histograms[name]
    
    def get_histograms(self):
        return self.histograms
    
