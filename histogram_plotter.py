import matplotlib.pyplot as plt
from matplotlib.patches import Patch


class HistogramPlotter:
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