import time
import numpy as np
import json
import hist
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
# import hist.dask as hda
# import dask_awkward as dak

from coffea import processor
from coffea.nanoevents import DelphesSchema

import coffea

# from distributed import Client


#################################################################################################################
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

#################################################################################################################
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
    
    
###################################################################################################################
cms_color = {
    "blue": "#3f90da",
    "orange": "#ffa90e",
    "red": "#bd1f01",
    "gray": "#94a4a2",
    "purple": "#832db6",
    "brown": "#a96b59",
    "dark_prange": "#e76300",
    "beige": "#b9ac70",
    "dark_gray": "#717581",
    "light_blue": "#92dadd"
}

class HistManager:
    def __init__(self):
        pass
        
    def define_histograms(self):
        pt_bins = [ 20.,  35.,  50.,  70., 100., 130., 165., 200., 250., 300.]
        self.hist_pt = hist.Hist(
            hist.axis.Variable(pt_bins, name="pt", label="Photon pT [GeV]"),
            storage="weight",
        )
        self.hist_eta = hist.Hist(
            hist.axis.Regular(15, -1.44, 1.44, name="eta", label="Photon η"),
            storage="weight",
        )
        
    def fill_histogram(self, events, var, weight):
        hist = getattr(self, f"hist_{var}")
        fill_kwargs = {var: ak.flatten(getattr(events.GoodPhotons, var))}
        hist.fill(**fill_kwargs, weight=weight)
       
        return hist
    
    def create_CMS_histograms(self, filename):
        
        with open(filename, 'r') as f:
            file = json.load(f)
            
        group_map = {}
        for i, dictionary in enumerate(file["headers"]):
            if i == 0:
                variable = dictionary["name"]
                continue
            group_map[i-1] = dictionary["name"]

        bin_list = []
        bin_centers = []
        histograms = {}
        values = {}
        errors = {"stat": [], "syst": []}
        for group in group_map.values():
            values[group] = []

        for i, dic in enumerate(file["values"]):
            if i == 0:
                bin_list.append(float(dic['x'][0]["low"]))
            bin_list.append(float(dic["x"][0]["high"]))

            for j,y in enumerate(dic["y"]):
                values[group_map[y["group"]]].append(float(y["value"]))
                if "stat" in y["errors"][0].values():
                    errors["stat"].append(float(y["errors"][0]["symerror"]))
                if "syst" in y["errors"][0].values():
                    errors["syst"].append(float(y["errors"][0]["symerror"]))
        for group in group_map.values():
            histograms[group] = hist.Hist(hist.axis.Variable(list(map(float, bin_list)), name=variable))
            histograms[group][:] = values[group]

        histograms["errors"] = errors

        return histograms

    def extract_hist_data(self, channel, normalize):
        
        self.histograms = self.hist_info[channel]["pt"]
        #Inforamtion from arXiv: 2201.07301v2 and https://www.hepdata.net/record/ins2013377
        cms_info = self.create_CMS_histograms(f"{channel}.json")
        self.histograms.update(cms_info)
        self.define_histograms()
        self.bins = self.hist_pt.axes[0].edges
        self.centers = []
        for i in range(len(self.bins)-1):
            center = (self.bins[i] + self.bins[i+1]) / 2
            self.centers.append(center)
        self.bin_widths = np.diff(self.bins)
        self.errors = self.histograms["errors"]
        self.data_values = np.array(self.histograms["Observed"].values())
        self.mc_values = np.array(self.histograms["Total prediction"].values())
        self.mc_components = {}
        for sample in self.histograms.keys():
            if sample not in ["Observed", "Total prediction", "errors"] and "Signal" not in sample:
                self.mc_components[sample] = np.array(self.histograms[sample].values())
        self.signal_components = {}
        for signal in self.histograms.keys():
            if "Signal" in signal:
                self.signal_components[signal] = np.array(self.histograms[signal].values())
        if normalize:
            self.errors["stat"] = self.errors["stat"]/(np.sum(self.data_values)*self.bin_widths)
            self.errors["syst"] = self.errors["syst"]/(np.sum(self.mc_values)*self.bin_widths)
            self.data_values = self.data_values/(np.sum(self.data_values)*self.bin_widths)
            for sample in self.mc_components.keys():
                self.mc_components[sample] = self.mc_components[sample]/(np.sum(self.mc_values)*self.bin_widths)
            self.mc_values = self.mc_values/(np.sum(self.mc_values)*self.bin_widths)
            for signal in self.histograms.keys():
                if "Signal" in signal:
                    self.signal_components[signal] = self.signal_components[signal]/(np.sum(self.signal_components[signal])*self.bin_widths)
                
        self.colors = [cms_color["orange"], cms_color["purple"], cms_color["red"], cms_color["beige"], cms_color["blue"], cms_color["dark_gray"],]
                
    def define_figure(self):
        self.fig, (self.ax, self.rax) = plt.subplots(
            2, 1, figsize=(10, 8), 
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.0}, 
            sharex=True
        )
    
    def plot_datamc(self, signals, channel):
        bottom = np.zeros(len(self.centers))
        for i, (name, values) in enumerate(self.mc_components.items()):
            self.ax.bar(
                self.centers, values, width=self.bin_widths,
                bottom=bottom, alpha=0.8, label=name, color=self.colors[i],
                edgecolor='black', linewidth=0.5
            )
            bottom += values
            
        self.ax.errorbar(
            self.centers, self.data_values, yerr=self.errors["stat"],
            fmt='o', color='black', markersize=6, capsize=4,
            linewidth=2, label='Data'
        )
        
        for i, signal in enumerate(signals):
            values = self.signal_components[signal]
            self.ax.bar(
                self.centers, values, width=self.bin_widths,
                alpha=1, label=signal, edgecolor=self.colors[-i-1],
                linewidth=1, fill=False
            )
        
        self.ax.fill_between(
            self.centers,
            self.mc_values - self.errors["syst"],
            self.mc_values + self.errors["syst"],
            step='mid', facecolor="None", alpha=0.9, hatch='////',
            label='Syst. Unc.', edgecolor='black', linewidth=0
        )
        
        cats = {"emu": "$e\mu$", "ee": "$ee$", "mumu": "$\mu\mu$"}
        self.ax.text(0.75, 0.45, cats[channel], transform=self.ax.transAxes, 
               fontsize=20, fontweight='bold', va='top')
        
        self.ax.set_ylabel('Events', fontsize=15)
        self.ax.legend(fontsize=10)
        self.ax.grid(True, alpha=0.3)
        self.ax.set_title('Data/MC Comparison', fontsize=16)
        
    def plot_ratio(self):
        ratio = self.data_values / self.mc_values
        ratio_err = self.errors["stat"] / self.mc_values

        # Plot ratio
        self.rax.errorbar(self.centers, ratio, yerr=ratio_err, 
                          fmt='o', color='black', markersize=5, linewidth=1.5,
                          capsize=3, capthick=1.5, label='Data/MC')

        self.rax.axhline(y=1.0, color='red', linestyle='--', linewidth=1.5)

        syst_unc_ratio = self.errors["syst"] / self.mc_values
        self.rax.fill_between(self.centers, 1 - syst_unc_ratio, 1 + syst_unc_ratio,
                 step='mid', facecolor="None", alpha=0.9, hatch='////',
                 label='Syst. Unc.', edgecolor='black', linewidth=0)

        # Set labels and limits
        self.rax.set_xlabel('p$_T$ [GeV]')
        self.rax.set_ylabel('Data/MC')
        self.rax.set_ylim(0.5, 1.5)
        self.rax.set_xlim(self.bins[0], self.bins[-1])
        self.rax.grid(True, alpha=0.3)
        self.rax.set_xlim(self.ax.get_xlim())

        # Add bin edges as x-ticks
        # ax_bottom.set_xticks(bins)

            
    def plot_histograms(self, hist_info, channel="total", signal=[], normalize=False):
        self.hist_info = hist_info
        name = f"photon_pt_{channel}"
        if len(signal):
            name = name + "_Signal"
        if normalize:
            name = name + "_normalized"
            
        self.extract_hist_data(channel, normalize)
        self.define_figure()
        self.plot_datamc(signal, channel)
        self.plot_ratio()
        plt.savefig(f"plots/{name}.png", dpi=300, bbox_inches="tight")
        # plt.savefig(f"plots/{name}.pdf", bbox_inches="tight")
        plt.close()

    def plot_signal_histograms(self, hist_dict, output_dir="plots", normalize = ""):
        
        os.makedirs(output_dir, exist_ok=True)

        for name, h in hist_dict.items():
            if normalize != "":
                bin_widths = np.diff(h.axes[0].edges)
                h = h/(h.sum().value*bin_widths)
            plt.figure()
            h.plot(
                histtype="fill",  
                # edgecolor="black",     
                linewidth=1.5,        
                alpha=0.7,
                color="blue",
            ) 

            if name == "pt":
                plt.xlabel("$p_T(\gamma)$ [GeV]")
                plt.title(f"Photon $p_T$, {normalize}, $T_{{mass}}=1000$")
            elif name == "eta":
                plt.xlabel("$\eta(\gamma)$")
                plt.title(f"Photon $\eta$, {normalize}, $T_{{mass}}=1000$")
            if normalize != "":
                plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)
                plt.gca().yaxis.get_offset_text().set_visible(True)  # Show multiplier
            plt.ylabel("Events")
            plt.grid(True, linestyle="--", alpha=0.5)
            

            if normalize == "":
                x = h.axes[0].centers
                y = h.values()
                plt.fill_between(x, y-np.sqrt(y), y+np.sqrt(y),
                                 step="mid",
                                 facecolor='none',
                                 hatch="////",
                                 edgecolor='black',
                                 linewidth=0,
                                 alpha=0.5
                                )            
                error_patch = Patch(
                    facecolor='none',
                    edgecolor='black',
                    linewidth=0,
                    hatch='////',
                    label='Stat Error'
                )
                plt.legend(handles=[error_patch], loc='best')

            plt.savefig(f"{output_dir}/{normalize}photon_{name}.png", dpi=300, bbox_inches="tight")
            plt.savefig(f"{output_dir}/{normalize}photon_{name}.pdf", bbox_inches="tight")
            plt.close() 
    
#####################################################################################################################
class TTPairTotatg(processor.ProcessorABC):
    
    def __init__(self):
        self.categories = ["total", "emu", "ee", "mumu"]
        self.variables = ["pt", "eta"]
        self.luminosity = 138e3
    
    def define_output_layout(self, dataset):
        self.output = {}
        self.output["nEvents"] = {}
        self.output["hists"] = {}
        for cat in self.categories:
            self.output["hists"][cat] = {}
            for var in self.variables:
                self.output["hists"][cat][var] = {}
        self.output["nEvents"]["primary"] = {}
        self.output["nEvents"]["selected"] = {}

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
            self.output["nEvents"]["selected"][dataset] = len(selected_events)
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