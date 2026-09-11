import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json
import copy

from compute_uncertainties import compute_systematic_uncertainties as com_sys

from coffea.util import load

from histogram_plotter import create_CMS_histograms, cms_color


class HistogramXSecPlotter:
    # Integrated luminosity used to turn event yields into cross sections
    # (fb^-1). Kept as one constant instead of a magic number sprinkled
    # across extract_hist_data() and the plot annotation.
    LUMI = 138.0
    HATCH = '///'

    def __init__(self, output):
        self.output = output

    # ------------------------------------------------------------------
    # small drawing helpers (shared by the main panel and the ratio panel)
    # ------------------------------------------------------------------
    @staticmethod
    def _step_edges(ax, bins, values, **kwargs):
        """Draw a 'post' step line across bin edges (repeats the last
        value so the line extends to the final bin edge)."""
        ax.step(bins, np.append(values, values[-1]), where='post', **kwargs)

    @classmethod
    def _band_edges(cls, ax, bins, lower, upper, edgecolor, label=None):
        """Draw a hatched uncertainty band across bin edges. Centralizing
        this means the 'append last value for the step edge' fix only has
        to be correct in one place."""
        ax.fill_between(
            bins, np.append(lower, lower[-1]), np.append(upper, upper[-1]),
            step='post', facecolor="none", alpha=0.9, hatch=cls.HATCH,
            edgecolor=edgecolor, linewidth=0, label=label
        )

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

    # Default target for 'auto' scaling: how tall (as a fraction of the
    # tallest MC curve) a scaled signal's peak should end up.
    AUTO_SCALE_TARGET_FRACTION = 0.5

    def _resolve_signal_scales(self, signal_scale, signals):
        """
        Turn the user-facing `signal_scale` argument into a per-signal
        {signal: factor} dict, since a single shared factor can leave a
        much smaller signal (e.g. a heavier, steeply-falling resonance)
        invisible even after scaling - see AUTO_SCALE_TARGET_FRACTION.

        Accepts:
          - a number: the same factor applied to every signal (old
            behavior, e.g. signal_scale=10)
          - a dict {signal_name: factor}: explicit per-signal factors,
            e.g. {"Signal_500": 10, "Signal_1000": 100}
          - "auto": for each signal, pick the power-of-ten factor that
            brings its peak bin up to roughly
            AUTO_SCALE_TARGET_FRACTION of the tallest MC curve's peak.
            Rounding to a power of ten keeps the "(x100)" legend tag
            honest and easy to read rather than an arbitrary decimal.
        """
        if signal_scale == "auto":
            mc_peak = max(np.max(v) for v in self.mc_values.values())
            scales = {}
            for signal in signals:
                raw_nominal = (
                    np.array(self.histograms[signal]["nominal"], dtype=float)
                    / (self.LUMI * self.bin_widths)
                )
                sig_peak = np.max(raw_nominal)
                if sig_peak <= 0:
                    scales[signal] = 1
                    continue
                needed = (self.AUTO_SCALE_TARGET_FRACTION * mc_peak) / sig_peak
                # Round down to a whole power of ten; never scale below 1x.
                scales[signal] = 10 ** max(0, int(np.floor(np.log10(needed)))) if needed > 1 else 1
            return scales

        if isinstance(signal_scale, dict):
            return {signal: signal_scale.get(signal, 1) for signal in signals}

        return {signal: signal_scale for signal in signals}

    def _extract_signal_component(self, signal):
        """Pull nominal/up/down for one signal, convert to a cross section
        (divide by luminosity and bin width), and apply that signal's
        entry in self.signal_scale (a per-signal {signal: factor} dict -
        see _resolve_signal_scales).

        The factor k multiplies nominal, total_up AND total_down. This is
        correct linear error propagation for a constant rescale
        (var(k*x) = k^2 * var(x)  =>  err(k*x) = k*err(x)), and it keeps
        total_up/total_down consistent with nominal so the fill_between
        band drawn from them (nominal - total_down, nominal + total_up)
        scales correctly instead of leaving the pre-scale band sitting
        under a post-scale line.
        """
        scale = self.signal_scale[signal]
        component = {}
        for val in ["nominal", "total_up", "total_down"]:
            raw = np.array(self.histograms[signal][val], dtype=float)
            component[val] = (raw / (self.LUMI * self.bin_widths)) * scale
        return component

    def extract_hist_data(self, signals, hist_name, normalize, signal_scale=1):

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

        # Needs mc_values (for 'auto' mode's mc_peak) but must run before
        # _extract_signal_component, which looks up self.signal_scale[signal].
        self.signal_scale = self._resolve_signal_scales(signal_scale, signals)

        self.signal_components = {
            signal: self._extract_signal_component(signal) for signal in signals
        }

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

            
            mc_norm = np.sum(self.mc_values["MG5+PYTHIA8"] * self.bin_widths)
            self.errors["theory unc."] = self._propagate_normalized_error(
                self.mc_values["MG5+PYTHIA8"], self.errors["theory unc."], self.bin_widths, norm=mc_norm
            )

            for sample in self.mc_values.keys():
                self.mc_values[sample] = self.mc_values[sample]/(np.sum(self.mc_values[sample]*self.bin_widths))

            for signal in signals:
                # NOTE: signal_scale (applied in _extract_signal_component)
                # is an overall constant factor on nominal/up/down alike, so
                # it cancels exactly here: sig_norm scales by the same
                # factor as nominal, and every term in
                # _propagate_normalized_error is scale-invariant under a
                # simultaneous rescale of values/errors/norm. So a
                # normalized signal shape looks identical whether
                # signal_scale is 1 or 10 - scaling only matters for the
                # absolute (non-normalized) cross section plot.
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
                
        self.colors = ['#D55E00', '#0072B2', '#009E73', cms_color["red"], '#009E73', cms_color["blue"], cms_color["dark_gray"],]
             
    def define_figure(self, normalize):
        plt.rcParams.update({
            'font.family': 'Times New Roman',      # controls all normal text
            'mathtext.fontset': 'custom',          # 'stix'/'stixsans' etc. ignore
            'mathtext.rm': 'Times New Roman',      # roman (upright) math text
            'mathtext.it': 'Times New Roman:italic',  # italic math text
            'mathtext.bf': 'Times New Roman:bold',    # bold math text
            'font.size': 12,
            'axes.labelsize': 16,
            'axes.titlesize': 16,
        })
        if normalize:
            self.fig, (self.ax) = plt.subplots(
                1, 1, figsize=(10, 8), 
                # gridspec_kw={"height_ratios": [3, 1], "hspace": 0.0}, 
                # sharex=True
            )
        else:           
            self.fig, (self.ax, self.rax) = plt.subplots(
                2, 1, figsize=(10, 8), 
                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.0}, 
                sharex=True
            )
    
    def plot_datamc(self, signals, hist_name, normalize):
        self._step_edges(
            self.ax, self.bins, self.mc_values["MG5+PYTHIA8"],
            alpha=0.8, label=r"$t\bar{t}\gamma$ (SM)", color="black", linewidth=2
        )
        # self._step_edges(
        #     self.ax, self.bins, self.mc_values["MG5+HERWIG7"],
        #     alpha=0.8, label="MG5+HERWIG7", color="blue", linewidth=2
        # )

        if not normalize:
            self.ax.errorbar(
                self.centers, self.data_values, yerr=self.errors["stat"],
                fmt='o', color='black', markersize=5, capsize=0,
                linewidth=2, label='Data'
            )

        # Signals are scaled per-signal by self.signal_scale[signal] (see
        # _resolve_signal_scales / _extract_signal_component); tag each
        # legend entry with its own factor so a "x10" line isn't mistaken
        # for the true cross section, and so two signals scaled by
        # different amounts (e.g. via signal_scale="auto") stay honest.
        # Scaling has no visible effect in the normalized case (it cancels
        # out), so the tag is only added for the absolute plot.
        for i, signal in enumerate(signals):
            mass = signal.split("_")[1]
            values = self.signal_components[signal]["nominal"]
            scale = self.signal_scale[signal]
            label = f"$m_{{T}}={mass}$ GeV"
            if (not normalize) and scale != 1:
                label += rf" ($\times{scale:g}$)"
            self._step_edges(
                self.ax, self.bins, values,
                alpha=1, label=label, color=self.colors[i], linewidth=2
            )

        theory_lower = self.mc_values["MG5+PYTHIA8"] - self.errors["theory unc."]
        theory_upper = self.mc_values["MG5+PYTHIA8"] + self.errors["theory unc."]
        self._band_edges(
            self.ax, self.bins, theory_lower, theory_upper,
            edgecolor='black', label="theory unc."
        )

        if not normalize:
            for i, signal in enumerate(signals):
                component = self.signal_components[signal]
                lower_sig = component["nominal"] - component["total_down"]
                upper_sig = component["nominal"] + component["total_up"]
                self._band_edges(self.ax, self.bins, lower_sig, upper_sig, edgecolor=self.colors[i])

        
        if hist_name == "diff_xsec_photon_pt":
            if normalize:
                self.ax.set_ylabel(r'1/$\sigma$ d$\sigma$/d$p_T$($\gamma$) [1/GeV]', fontsize=20)
            else:
                self.ax.set_ylabel(r'd$\sigma$/d$p_T$($\gamma$) [fb/GeV]', fontsize=20)
        else:
            self.ax.set_ylabel(rf'1/$\sigma$ d$\sigma$/d{self.x_axis_name}', fontsize=20)
        self.ax.minorticks_on()
        self.ax.legend(ncol=1, loc='upper right', fontsize=14,
                                 bbox_to_anchor=(0.97, 0.96),
                                 frameon=True, fancybox=True,
                                 framealpha=0.7, edgecolor="gray", borderpad=0.6)
        self.ax.grid(True, which='major', linestyle='-', linewidth=0.7, alpha=0.3)
        self.ax.set_xlabel(self.x_axis_name, fontsize=20)
        
        self.ax.tick_params(axis='both', which='minor', length=3)
        self.ax.tick_params(axis='both', which='major', length=7)
        self.ax.tick_params(axis='both', which='both' , top=True, right=True, direction='in', labelsize=15)
        self.ax.set_xlim(self.bins[0], self.bins[-1])
        
        for spine in self.ax.spines.values():
            spine.set_linewidth(2) 
        
        
        self.ax.text(0.98, 1.02, f'{self.LUMI:g} $fb^{{-1}}$ [13 TeV]',
                transform=self.ax.transAxes,  # Use axes coordinates (0 to 1)
                fontsize=18, 
                # fontweight='bold',
                ha='right',  # horizontal alignment right
                va='bottom',  # vertical alignment bottom
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none'))  # optional background
                
    def plot_ratio(self):
        mc_ratio_pythia = self.mc_values["MG5+PYTHIA8"] / self.data_values
        data_ratio = self.data_values / self.data_values
        ratio_err = self.errors["stat"] / self.data_values

        self._step_edges(self.rax, self.bins, mc_ratio_pythia, alpha=0.8, color="black", linewidth=2)
        # self._step_edges(self.rax, self.bins, self.mc_values["MG5+HERWIG7"] / self.data_values,
        #                   alpha=0.8, color="blue", linewidth=2)

        # Plot ratio
        self.rax.errorbar(self.centers, data_ratio, yerr=ratio_err, 
                          fmt='o', color='black', markersize=5, linewidth=1.5,
                          capsize=0, capthick=1.5, label='Data/MC')

        self.rax.axhline(y=1.0, color='black', linestyle='--', linewidth=1.5)

        theory_unc_ratio = self.errors["theory unc."] / self.data_values
        ratio_lower = mc_ratio_pythia - theory_unc_ratio
        ratio_upper = mc_ratio_pythia + theory_unc_ratio
        self._band_edges(self.rax, self.bins, ratio_lower, ratio_upper, edgecolor='black', label='Syst. Unc.')

        # Set labels and limits
        self.rax.set_xlabel(self.x_axis_name, fontsize=20)
        self.rax.set_ylabel('Pred./Obs.', fontsize=20)
        self.rax.set_ylim(0.5, 1.4)
        
        self.rax.minorticks_on()
        self.rax.grid(True, which='major', linestyle='-', linewidth=0.7, alpha=0.3)

        self.rax.tick_params(axis='both', which='minor', length=3)
        self.rax.tick_params(axis='both', which='major', length=7)
        self.rax.tick_params(axis='both', which='both' , top=True, right=True, direction='in', labelsize=15)
        
        for spine in self.rax.spines.values():
            spine.set_linewidth(2) 

        # Add bin edges as x-ticks
        # ax_bottom.set_xticks(bins)

            
    def plot_histograms(self, hist_info, hist_name, signals=[], normalize=False, signal_scale=1):
        """
        signal_scale: how much to blow up each signal's nominal value (and
        its total_up/total_down errors, by the same factor - correct
        linear error propagation for a constant rescale) so a tiny signal
        cross section becomes visible next to MG5+PYTHIA8. Has no effect
        on the normalized ('unit area') plots since an overall constant
        factor cancels out of that transform - see the note in
        extract_hist_data(). Accepts:
          - a number, applied to every signal alike (e.g. 10)
          - a dict {signal_name: factor} for per-signal control, e.g.
            {"Signal_500": 10, "Signal_1000": 100} - useful when signals
            span masses with very different cross sections and a single
            shared factor leaves the smaller one invisible
          - "auto" to pick a power-of-ten factor per signal automatically
            (see _resolve_signal_scales)
        """
        self.hist_info = hist_info
        name = hist_name
        if len(signals) == 1:
            mass = signals[0].split("_")[1]
            name = name + "_" + mass
        if normalize:
            name = name + "_normalized"
        elif signal_scale != 1:
            name = name + "_signalscaled"

        self.extract_hist_data(signals, hist_name, normalize, signal_scale)
        self.define_figure(normalize)
        self.plot_datamc(signals, hist_name, normalize)
        if not normalize:
            self.plot_ratio()
        plt.tight_layout()
        # plt.savefig(f"plots/{name}.png", dpi=300, bbox_inches="tight")
        plt.savefig(f"plots/{name}.pdf", bbox_inches="tight")
        plt.close()


if __name__ == "__main__":
    SIGNALS = ["Signal_500"]
    SIGNALS_norm = ["Signal_500", "Signal_1000", "Signal_1500"]
    # "auto" picks a power-of-ten factor per signal so each is comparably
    # visible against MG5+PYTHIA8, even though e.g. Signal_1000's cross
    # section is much smaller than Signal_500's. For explicit control
    # instead, pass a dict, e.g. {"Signal_500": 10, "Signal_1000": 100}.
    SIGNAL_SCALE = "auto"
    HIST_NAMES = ['diff_xsec_photon_pt']

    output = load("../output.coffea")
    xsec_hist_plotter = HistogramXSecPlotter(copy.deepcopy(output))

    for hist in HIST_NAMES:
        xsec_hist_plotter.plot_histograms(
            copy.deepcopy(output["hists"]["total"]), hist, signals=SIGNALS, signal_scale=SIGNAL_SCALE
        )
        xsec_hist_plotter.plot_histograms(
            copy.deepcopy(output["hists"]["total"]), hist, signals=SIGNALS_norm, normalize=True
        )