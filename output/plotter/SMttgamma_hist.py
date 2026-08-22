import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json
import copy
from scipy import stats

from coffea.util import load

from histogram_plotter import create_CMS_histograms, cms_color

class HistogramXSecPlotter:
    def __init__(self, output):
        self.output = output
        
    def compute_chi2(self, values1, values2, sigma, n_params=0):
        """
        values1, values2: arrays of per-bin values (already normalized if desired)
        sigma: combined per-bin uncertainty array (same shape)
        n_params: number of fitted/constrained parameters to subtract from DoF
                  (use 1 if both histograms were normalized to unit area)
        """
        chi2 = np.sum((values1 - values2)**2 / sigma**2)
        dof = len(values1) - n_params
        p_value = 1 - stats.chi2.cdf(chi2, dof)
        return chi2, dof, chi2/dof, p_value
        
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
        self.errors["stat_pythia"] = np.sqrt(self.histograms["MG5+PYTHIA8"].values())
        self.errors["stat_delphes"] = np.sqrt(self.histograms["SMttgamma"].values())
        self.errors["stat_pythia"] = self.errors["stat_pythia"]/(np.sqrt(138)*self.bin_widths)
        self.errors["stat_delphes"] = self.errors["stat_delphes"]/(138*(self.bin_widths))
        self.errors["stat"] = self.errors["stat"]/self.bin_widths
        self.errors["theory unc."] = self.errors["theory unc."]/self.bin_widths
        self.data_values = np.array(self.histograms["Observed"].values())/self.bin_widths
        
        self.mc_values = {}
        for sample in self.histograms.keys():
            if sample in ["MG5+PYTHIA8", "MG5+HERWIG7"]:
                self.mc_values[sample] = np.array(self.histograms[sample].values())/self.bin_widths
                
        self.signal_components = {}
        for signal in self.histograms.keys():
            if (not "MG5" in signal) and (signal != "errors") and (signal != "Observed"):
                # epsilon = np.array(self.histograms[signal].values())/self.output["cutflow"]["total"][signal]["primary"]
                # print(epsilon)
                self.signal_components[signal] = np.array(self.histograms[signal].values())/(138*(self.bin_widths))
                # print(signal, np.sum(np.array(self.histograms[signal].values())/(138)))
                print(signal,np.array(self.histograms[signal].values())/(138*self.bin_widths))
        self.x_axis_name = self.histograms["Observed"].axes[0].label
        
        if normalize:
            self.errors["stat"] = self.errors["stat"]/(np.sum(self.data_values*self.bin_widths))
            self.errors["theory unc."] = self.errors["theory unc."]/(np.sum(self.mc_values["MG5+PYTHIA8"]*self.bin_widths))
            self.data_values = self.data_values/(np.sum(self.data_values*self.bin_widths))
            for sample in self.mc_values.keys():
                self.mc_values[sample] = self.mc_values[sample]/(np.sum(self.mc_values[sample]*self.bin_widths))
            for signal in self.histograms.keys():
                if not "MG5" in signal and signal != "errors" and (signal != "Observed"):
                    self.signal_components[signal] = self.signal_components[signal]/(np.sum(self.signal_components[signal]*self.bin_widths))
                
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
        # self.ax.step(
        #     self.bins, np.append(self.mc_values["MG5+HERWIG7"],self.mc_values["MG5+HERWIG7"][-1]), where='post',
        #     alpha=0.8, label="MG5+HERWIG7", color="blue", linewidth=2
        # )
            
        # self.ax.errorbar(
        #     self.centers, self.data_values, yerr=self.errors["stat"],
        #     fmt='o', color='black', markersize=5, capsize=0,
        #     linewidth=2, label='Data'
        # )
        
        for i, signal in enumerate(signals):
            values = self.signal_components[signal]
            self.ax.step(
                self.bins, np.append(values, values[-1]), where='post',
                alpha=1, linestyle="-", label="MG5+PYTHIA+Delphes", color=self.colors[i], linewidth=2
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
        
        if hist_name == "diff_xsec_photon_pt":
            if normalize == True:
                self.ax.set_ylabel('1/$\sigma$ d$\sigma$/d$p_T$($\gamma$) [1/GeV]', fontsize=20)
            else:
                self.ax.set_ylabel('d$\sigma$/d$p_T$($\gamma$) [fb/GeV]', fontsize=20)
        else:
            self.ax.set_ylabel(f'1/$\sigma$ d$\sigma$/d{self.x_axis_name}', fontsize=20)
        self.ax.minorticks_on()
        
        self.ax.grid(True, alpha=0.3)
        
        self.ax.tick_params(axis='both', which='major', labelsize=14, width=2, length=8)
        self.ax.tick_params(axis='both', which='minor', width=1.5, length=4)
        
        for spine in self.ax.spines.values():
            spine.set_linewidth(2) 
            
        
        sigma = np.sqrt(self.errors["theory unc."]**2 + self.errors["stat_delphes"]**2)  #+ self.errors["stat_pythia"]**2)  
        chi2, dof, chi2_per_dof, pval = self.compute_chi2(
            self.signal_components["SMttgamma"], 
            self.mc_values["MG5+PYTHIA8"], 
            sigma,
            n_params=1 if normalize else 0
        )
        print(f"chi2/dof = {chi2_per_dof:.2f} (chi2={chi2:.2f}, dof={dof}), p={pval:.3f}")
        # legend_title = legend_title = rf'tt$\gamma$' + '\n' + rf'$\chi^2$/DoF = {chi2_per_dof:.2f}'
        # self.ax.legend(title=legend_title, title_fontsize=18, fontsize=14)
        dummy_handle = plt.Line2D([0], [0], color='none', marker='', linestyle='none')
        chi2_label = rf'$\chi^2$/DoF = {chi2_per_dof:.2f}'

        # Get existing legend handles and labels
        handles, labels = self.ax.get_legend_handles_labels()

        # Add the chi2 entry at the end
        handles.append(dummy_handle)
        labels.append(chi2_label)

        # Create legend with all entries including chi2
        self.ax.legend(handles=handles, labels=labels, title="tt$\gamma$", 
                       title_fontsize=18, fontsize=14)
        
        # self.ax.text(0.5, 0.8, f'$\\frac{{Xi^2}}{{DoF}} = {chi2_per_dof:.2f}$', 
        #         transform=self.ax.transAxes,  # Use axes coordinates (0 to 1)
        #         fontsize=18, 
        #         # fontweight='bold',
        #         ha='right',  # horizontal alignment right
        #         va='bottom',  # vertical alignment bottom
        #         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none'))  # optional background
        
        
        self.ax.text(0.98, 1.02, '138 $fb^{-1}$ [13 TeV]', 
                transform=self.ax.transAxes,  # Use axes coordinates (0 to 1)
                fontsize=18, 
                # fontweight='bold',
                ha='right',  # horizontal alignment right
                va='bottom',  # vertical alignment bottom
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none'))  # optional background
        
        
    def plot_ratio(self):
        signal_ratio_pythia = self.signal_components["SMttgamma"] / self.mc_values["MG5+PYTHIA8"]
        # mc_ratio_herwig = self.mc_values["MG5+HERWIG7"] / self.data_values
        # mc_ratio = self.data_values /self.data_values
        # ratio_err = self.errors["stat"] / self.data_values

        self.rax.step(
            self.bins, np.append(signal_ratio_pythia, signal_ratio_pythia[-1]), where='post',
            alpha=0.8, color=self.colors[0], linewidth=2
        )
        # self.rax.step(
        #     self.bins, np.append(mc_ratio_herwig, mc_ratio_herwig[-1]), where='post',
        #     alpha=0.8, color="blue", linewidth=2
        # )

        # Plot ratio
        # self.rax.errorbar(self.centers, data_ratio, yerr=ratio_err, 
        #                   fmt='o', color='black', markersize=5, linewidth=1.5,
        #                   capsize=0, capthick=1.5, label='Data/MC')

        self.rax.axhline(y=1.0, color='black', linestyle='--', linewidth=1.5)

        theory_unc_ratio = self.errors["theory unc."] / self.mc_values["MG5+PYTHIA8"]
        lower = 1 - theory_unc_ratio
        upper = 1 + theory_unc_ratio
        self.rax.fill_between(self.bins, np.append(lower, lower[-1]), np.append(upper, upper[-1]),
                 step='post', facecolor="None", alpha=0.9, hatch='////',
                 label='Syst. Unc.', edgecolor='k', linewidth=0)

        # Set labels and limits
        self.rax.set_xlabel(self.x_axis_name, fontsize=20)
        self.rax.set_ylabel('ratio', fontsize=20)
        self.rax.set_ylim(0.6, 1.4)
        self.rax.set_xlim(self.bins[0], self.bins[-1])
        self.rax.minorticks_on()
        self.rax.grid(True, alpha=0.3)
        self.rax.set_xlim(self.ax.get_xlim())
        
        self.rax.tick_params(axis='both', which='major', labelsize=14, width=2, length=8)
        self.rax.tick_params(axis='both', which='minor', width=1.5, length=4)
        
        for spine in self.rax.spines.values():
            spine.set_linewidth(2) 

        # Add bin edges as x-ticks
        # ax_bottom.set_xticks(bins)

            
    def plot_histograms(self, hist_info, hist_name, signal=[], normalize=False):
        self.hist_info = hist_info
        name = "SMttgama" + hist_name
        if normalize:
            name = name + "_normalized"
            
        self.extract_hist_data(hist_name, normalize)
        self.define_figure()
        self.plot_datamc(signal, hist_name, normalize)
        self.plot_ratio()
        plt.savefig(f"plots/{name}.png", dpi=300, bbox_inches="tight")
        # plt.savefig(f"plots/{name}.pdf", bbox_inches="tight")
        plt.close()


#############################
###### ttg plotter
#############################
if __name__ == "__main__":
    # hist_plotter = HistogramPlotter()
    output = load("../SMttgamma.coffea")
    xsec_hist_plotter = HistogramXSecPlotter(output)
    histograms = {}
    histograms['diff_xsec_photon_pt'] = {}
    histograms['diff_xsec_photon_pt']["SMttgamma"] = None
    for dts in output["hists"]["total"]['diff_xsec_photon_pt']:
        # if "SL" in dts:
        #     continue
        if histograms['diff_xsec_photon_pt']["SMttgamma"] is None:
            histograms['diff_xsec_photon_pt']["SMttgamma"] = copy.deepcopy(output["hists"]["total"]['diff_xsec_photon_pt'][dts])
        else:
            histograms['diff_xsec_photon_pt']["SMttgamma"] += copy.deepcopy(output["hists"]["total"]['diff_xsec_photon_pt'][dts])
        print()
    xsec_hist_plotter.plot_histograms(histograms, 'diff_xsec_photon_pt', signal=["SMttgamma"]) # "Signal_500", "Signal_1000", "Signal_1500", "Signal_2000"
    # xsec_hist_plotter.plot_histograms(histograms, 'diff_xsec_photon_pt', signal=["SMttgamma"], normalize=True) # 
