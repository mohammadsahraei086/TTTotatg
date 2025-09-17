import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json

from histogram_plotter import create_CMS_histograms, cms_color

class HistogramXSecPlotter:
    def __init__(self):
        pass
        
    def extract_hist_data(self, hist_name, normalize):
        
        self.histograms = self.hist_info[hist_name]
        #Inforamtion from arXiv: 2201.07301v2 and https://www.hepdata.net/record/ins2013377
        cms_info = create_CMS_histograms(f"json_files/{hist_name}.json")
        self.histograms.update(cms_info)
        self.bins = self.histograms["Observed"].axes[0].edges
        self.centers = []
        for i in range(len(self.bins)-1):
            center = (self.bins[i] + self.bins[i+1]) / 2
            self.centers.append(center)
        self.bin_widths = np.diff(self.bins)
        self.errors = self.histograms["errors"]
        self.errors["stat"] = self.errors["stat"]/self.bin_widths
        self.errors["theory unc."] = self.errors["theory unc."]/self.bin_widths
        self.data_values = np.array(self.histograms["Observed"].values())/self.bin_widths
        self.mc_values = {}
        for sample in self.histograms.keys():
            if sample in ["MG5+PYTHIA8", "MG5+HERWIG7"]:
                self.mc_values[sample] = np.array(self.histograms[sample].values())/self.bin_widths
        self.signal_components = {}
        for signal in self.histograms.keys():
            if "Signal" in signal:
                # epsilon = self.histograms[signal].values()/10000
                self.signal_components[signal] = np.array(self.histograms[signal].values())/(138*self.bin_widths)
        self.x_axis_name = self.histograms["Observed"].axes[0].label
        if normalize:
            self.errors["stat"] = self.errors["stat"]/(np.sum(self.data_values)*self.bin_widths)
            self.errors["theory unc."] = self.errors["theory unc."]/(np.sum(self.mc_values["MG5+PYTHIA8"])*self.bin_widths)
            self.data_values = self.data_values/(np.sum(self.data_values)*self.bin_widths)
            for sample in self.mc_values.keys():
                self.mc_values[sample] = self.mc_values[sample]/(np.sum(self.mc_values[sample])*self.bin_widths)
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
    
    def plot_datamc(self, signals, hist_name, normalize):
        self.ax.step(
            self.bins, np.append(self.mc_values["MG5+PYTHIA8"], self.mc_values["MG5+PYTHIA8"][-1]), where='post',
            alpha=0.8, label="MG5+PYTHIA8", color="black", linewidth=2
        )
        self.ax.step(
            self.bins, np.append(self.mc_values["MG5+HERWIG7"],self.mc_values["MG5+HERWIG7"][-1]), where='post',
            alpha=0.8, label="MG5+HERWIG7", color="blue", linewidth=2
        )
            
        self.ax.errorbar(
            self.centers, self.data_values, yerr=self.errors["stat"],
            fmt='o', color='black', markersize=5, capsize=0,
            linewidth=2, label='Data'
        )
        
        for i, signal in enumerate(signals):
            values = self.signal_components[signal]
            self.ax.step(
                self.bins, np.append(values, values[-1]), where='post',
                alpha=1, linestyle="dashed", label=signal, color=self.colors[i], linewidth=2
            )
        
        lower = self.mc_values["MG5+PYTHIA8"] - self.errors["theory unc."]
        upper = self.mc_values["MG5+PYTHIA8"] + self.errors["theory unc."]
        self.ax.fill_between(
            self.bins,
            np.append(lower, lower[-1]),
            np.append(upper, upper[-1]),
            step='post', facecolor="None", alpha=0.9, hatch='////',
            label="theory unc.", edgecolor='black', linewidth=0
        )
        
        # cats = {"emu": "$e\mu$", "ee": "$ee$", "mumu": "$\mu\mu$"}
        # self.ax.text(0.75, 0.45, cats[channel], transform=self.ax.transAxes, 
        #        fontsize=20, fontweight='bold', va='top')
        
        self.ax.set_ylabel('$d\sigma/dp_T(\gamma)[fb/GeV]$', fontsize=15)
        self.ax.minorticks_on()
        self.ax.legend(fontsize=10)
        self.ax.grid(True, alpha=0.3)
        if normalize:
            self.ax.set_title(f'Normalized differential cross section/{self.x_axis_name}', fontsize=16)
        else:
            self.ax.set_title(f'Differential cross section/{self.x_axis_name}', fontsize=16)
        
    def plot_ratio(self):
        mc_ratio_pythia = self.mc_values["MG5+PYTHIA8"] / self.data_values
        mc_ratio_herwig = self.mc_values["MG5+HERWIG7"] / self.data_values
        data_ratio = self.data_values /self.data_values
        ratio_err = self.errors["stat"] / self.data_values

        self.rax.step(
            self.bins, np.append(mc_ratio_pythia, mc_ratio_pythia[-1]), where='post',
            alpha=0.8, color="black", linewidth=2
        )
        self.rax.step(
            self.bins, np.append(mc_ratio_herwig, mc_ratio_herwig[-1]), where='post',
            alpha=0.8, color="blue", linewidth=2
        )

        # Plot ratio
        self.rax.errorbar(self.centers, data_ratio, yerr=ratio_err, 
                          fmt='o', color='black', markersize=5, linewidth=1.5,
                          capsize=0, capthick=1.5, label='Data/MC')

        self.rax.axhline(y=1.0, color='black', linestyle='--', linewidth=1.5)

        theory_unc_ratio = self.errors["theory unc."] / self.data_values
        lower = mc_ratio_pythia - theory_unc_ratio
        upper = mc_ratio_pythia + theory_unc_ratio
        self.rax.fill_between(self.bins, np.append(lower, lower[-1]), np.append(upper, upper[-1]),
                 step='post', facecolor="None", alpha=0.9, hatch='////',
                 label='Syst. Unc.', edgecolor='k', linewidth=0)

        # Set labels and limits
        self.rax.set_xlabel(self.x_axis_name, fontsize=15)
        self.rax.set_ylabel('Pred./Obs.')
        self.rax.set_ylim(0.5, 1.5)
        self.rax.set_xlim(self.bins[0], self.bins[-1])
        self.rax.minorticks_on()
        self.rax.grid(True, alpha=0.3)
        self.rax.set_xlim(self.ax.get_xlim())

        # Add bin edges as x-ticks
        # ax_bottom.set_xticks(bins)

            
    def plot_histograms(self, hist_info, hist_name, signal=[], normalize=False):
        self.hist_info = hist_info
        name = hist_name
        if len(signal):
            name = name + "_Signal"
        if normalize:
            name = name + "_normalized"
            
        self.extract_hist_data(hist_name, normalize)
        self.define_figure()
        self.plot_datamc(signal, hist_name, normalize)
        self.plot_ratio()
        plt.savefig(f"plots/{name}.png", dpi=300, bbox_inches="tight")
        # plt.savefig(f"plots/{name}.pdf", bbox_inches="tight")
        plt.close()