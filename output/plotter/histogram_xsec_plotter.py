import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json
import copy

from compute_uncertainties import compute_systematic_uncertainties as com_sys

from coffea.util import load

from histogram_plotter import create_CMS_histograms, cms_color

class HistogramXSecPlotter:
    def __init__(self, output):
        self.output = output

    @staticmethod
    def _propagate_normalized_error(values, errors, bin_widths, norm=None):
        """
        Correct error propagation for a "normalize to unit area" transform.

        If f_i = v_i / N with N = sum_j v_j * w_j (w_j = bin width), then N
        itself depends on every bin, so the naive err_i/N is WRONG: it treats
        N as a fixed constant and ignores the anti-correlation the
        normalization introduces between bins (a fluctuation up in one bin
        pulls every other normalized bin down slightly, and pulls the bin
        itself down too, since it also inflates N).

        The correct (first-order / linear error propagation) result, assuming
        the input per-bin errors are uncorrelated with each other, comes from
        the Jacobian:

            d f_i / d v_k = delta_ik / N  -  v_i * w_k / N^2

            var(f_i) = sum_k (df_i/dv_k)^2 * err_k^2
                     = err_i^2 / N^2
                       - 2 * v_i * w_i * err_i^2 / N^3
                       + (v_i^2 / N^4) * sum_k (w_k * err_k)^2

        Sanity check: for a single-bin histogram this correctly gives
        var(f_i) = 0 (normalizing a single bin to unit area leaves no
        freedom, so it can't carry an uncertainty).

        Note: this assumes the errors passed in are uncorrelated bin-to-bin
        (true for independent stat/Poisson errors; an approximation for an
        already-combined symmetric systematic envelope, since the original
        bin-to-bin correlations of the underlying variations aren't
        preserved once everything's been collapsed into a single up/down
        band). If you need to be exact for the signal systematics, normalize
        each MUR/MUF/PDF variation histogram to unit area *before* taking
        differences from nominal in compute_uncertainties.py, rather than
        normalizing the already-combined total_up/total_down band here.
        """
        values = np.asarray(values, dtype=float)
        errors = np.asarray(errors, dtype=float)
        bin_widths = np.asarray(bin_widths, dtype=float)

        if norm is None:
            norm = np.sum(values * bin_widths)

        term1 = errors**2 / norm**2
        term2 = 2 * values * bin_widths * errors**2 / norm**3
        term3 = (values**2 / norm**4) * np.sum((bin_widths * errors) ** 2)

        var = term1 - term2 + term3
        return np.sqrt(np.clip(var, 0, None))

    def extract_hist_data(self, signals, hist_name, normalize):
        
        # self.histograms = self.hist_info[hist_name]
        self.histograms = com_sys(self.hist_info[hist_name])
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
        for signal in signals:
            self.signal_components[signal] = {}
            for val in ['nominal', 'total_up', 'total_down']:
                self.signal_components[signal][val] = np.array(self.histograms[signal][val])/(138*(self.bin_widths))
                
        self.x_axis_name = self.histograms["Observed"].axes[0].label
        
        if normalize:
            # NOTE: order matters here. Each error must be propagated using
            # the *un-normalized* central values (and the norm factor they
            # imply) before those central values get divided down, since the
            # Jacobian above needs v_i and N from the same (pre-normalization)
            # stage.

            data_norm = np.sum(self.data_values * self.bin_widths)
            self.errors["stat"] = self._propagate_normalized_error(
                self.data_values, self.errors["stat"], self.bin_widths, norm=data_norm
            )
            self.data_values = self.data_values / data_norm

            # theory unc. is drawn as a band around MG5+PYTHIA8, so it must be
            # propagated using MG5+PYTHIA8's own normalization, evaluated
            # before MG5+PYTHIA8 itself gets normalized below.
            mc_norm = np.sum(self.mc_values["MG5+PYTHIA8"] * self.bin_widths)
            self.errors["theory unc."] = self._propagate_normalized_error(
                self.mc_values["MG5+PYTHIA8"], self.errors["theory unc."], self.bin_widths, norm=mc_norm
            )

            for sample in self.mc_values.keys():
                self.mc_values[sample] = self.mc_values[sample]/(np.sum(self.mc_values[sample]*self.bin_widths))

            for signal in signals:
                nominal = self.signal_components[signal]["nominal"]
                sig_norm = np.sum(nominal * self.bin_widths)
                up = self._propagate_normalized_error(
                    nominal, self.signal_components[signal]["total_up"], self.bin_widths, norm=sig_norm
                )
                down = self._propagate_normalized_error(
                    nominal, self.signal_components[signal]["total_down"], self.bin_widths, norm=sig_norm
                )
                self.signal_components[signal]["nominal"] = nominal / sig_norm
                self.signal_components[signal]["total_up"] = up
                self.signal_components[signal]["total_down"] = down
                
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
            
        self.ax.errorbar(
            self.centers, self.data_values, yerr=self.errors["stat"],
            fmt='o', color='black', markersize=5, capsize=0,
            linewidth=2, label='Data'
        )
        
        for i, signal in enumerate(signals):
            mass = signal.split("_")[1]
            values = self.signal_components[signal]["nominal"]
            self.ax.step(
                self.bins, np.append(values, values[-1]), where='post',
                alpha=1, linestyle="-", label=f"$M_{{T}}={mass}$", color=self.colors[i], linewidth=2
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

        for i, signal in enumerate(signals):
            lower_sig = self.signal_components[signal]["nominal"] - self.signal_components[signal]["total_down"]
            upper_sig = self.signal_components[signal]["nominal"] + self.signal_components[signal]["total_up"]
            self.ax.fill_between(
                self.bins,
                np.append(lower_sig, lower_sig[-1]),
                np.append(upper_sig, upper_sig[-1]),
                step='post', facecolor="None", alpha=0.9, hatch='////',
                edgecolor=self.colors[i], linewidth=0
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
        self.ax.legend(fontsize=14)
        self.ax.grid(True, alpha=0.3)
        
        self.ax.tick_params(axis='both', which='major', labelsize=14, width=2, length=8)
        self.ax.tick_params(axis='both', which='minor', width=1.5, length=4)
        
        for spine in self.ax.spines.values():
            spine.set_linewidth(2) 
        # if normalize:
        #     self.ax.set_title(f'Normalized differential cross section/{self.x_axis_name}', fontsize=16)
        # else:
        #     self.ax.set_title(f'Differential cross section/{self.x_axis_name}', fontsize=16)
        
#         current_ymax = self.ax.get_ylim()[1]  # Get current upper limit
#         self.ax.set_ylim(bottom=0, top=current_ymax) 
        
#         yticks = self.ax.get_yticks()
#         yticklabels = [str(label) for label in yticks]
#         print(yticks, "#", yticklabels)
#         yticklabels[0] = ""
#         self.ax.set_yticklabels(yticklabels)
        
        
        self.ax.text(0.98, 1.02, '138 $fb^{-1}$ [13 TeV]', 
                transform=self.ax.transAxes,  # Use axes coordinates (0 to 1)
                fontsize=18, 
                # fontweight='bold',
                ha='right',  # horizontal alignment right
                va='bottom',  # vertical alignment bottom
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none'))  # optional background
                
    def plot_ratio(self):
        mc_ratio_pythia = self.mc_values["MG5+PYTHIA8"] / self.data_values
        mc_ratio_herwig = self.mc_values["MG5+HERWIG7"] / self.data_values
        data_ratio = self.data_values /self.data_values
        ratio_err = self.errors["stat"] / self.data_values

        self.rax.step(
            self.bins, np.append(mc_ratio_pythia, mc_ratio_pythia[-1]), where='post',
            alpha=0.8, color="black", linewidth=2
        )
        # self.rax.step(
        #     self.bins, np.append(mc_ratio_herwig, mc_ratio_herwig[-1]), where='post',
        #     alpha=0.8, color="blue", linewidth=2
        # )

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
        self.rax.set_xlabel(self.x_axis_name, fontsize=20)
        self.rax.set_ylabel('Pred./Obs.', fontsize=20)
        self.rax.set_ylim(0.5, 1.4)
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

            
    def plot_histograms(self, hist_info, hist_name, signals=[], normalize=False):
        self.hist_info = hist_info
        name = hist_name
        if len(signals) == 1:
            mass = signals[0].split("_")[1]
            name = name + "_" + mass
        if normalize:
            name = name + "_normalized"
            
        self.extract_hist_data(signals, hist_name, normalize)
        self.define_figure()
        self.plot_datamc(signals, hist_name, normalize)
        self.plot_ratio()
        plt.savefig(f"plots/{name}.png", dpi=300, bbox_inches="tight")
        # plt.savefig(f"plots/{name}.pdf", bbox_inches="tight")
        plt.close()
        
        
if __name__ == "__main__":
    # hist_plotter = HistogramPlotter()
    output = load("../output.coffea")
    xsec_hist_plotter = HistogramXSecPlotter(copy.deepcopy(output))
    histograms = {}
    # for mass in [500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]:
    #     for hist in ['diff_xsec_photon_pt', "deltaphi_ll"]: # , "deltaphi_ll"
    #         xsec_hist_plotter.plot_histograms(copy.deepcopy(output["hists"]["total"]), hist, signals=[f"Signal_{mass}"]) # "Signal_500", "Signal_1000"
    #         xsec_hist_plotter.plot_histograms(copy.deepcopy(output["hists"]["total"]), hist, signals=[f"Signal_{mass}"], normalize=True) #
    for hist in ['diff_xsec_photon_pt', "deltaphi_ll"]: # , "deltaphi_ll"
        xsec_hist_plotter.plot_histograms(copy.deepcopy(output["hists"]["total"]), hist, signals=["Signal_500", "Signal_1000"]) # "Signal_500", "Signal_1000"
        xsec_hist_plotter.plot_histograms(copy.deepcopy(output["hists"]["total"]), hist, signals=["Signal_500", "Signal_1000"], normalize=True) #