from scipy.optimize import curve_fit
import scipy
from scipy.stats import chi2

from coffea.util import load
from helper import *


class ComputeLimit:
    def __init__(self, mass, var, in_br = False):
        self.in_br = in_br
        self.mass = mass
        self.var = var
        
        self.width_factor = width_prefactor(mass)
        self.width_wb = gen_val[f"Signal_{mass}"]["width_wb"]
        
        self.a0, self.a1, self.a2 = xsec_factors(mass)
        
        self.acceptance = compute_acceptance(mass, var)
        
        self.V_inv = HepDataParser().get_inverse_covariance_matrix(self.var)
    
    def compute_branching_ratios(self, g3g, g3gamma):
        
        width_tg = self.width_factor * g3g ** 2
        width_tgamma = (3/4) * self.width_factor * g3gamma ** 2
        
        Gamma_total = width_tg + width_tgamma + self.width_wb

        B_g = width_tg / Gamma_total
        B_gamma = width_tgamma / Gamma_total

        return B_g, B_gamma
    
    def compute_xsec(self, g3g):
        
        return self.a0 + self.a1 * g3g ** 2 + self.a2 * g3g ** 4
        
    def get_signal_vector(self, g3g, g3gamma):
        
        xsec_TT = self.compute_xsec(g3g)
        b_g, b_gamma = self.compute_branching_ratios(g3g, g3gamma)
        f1gamma = 2 * b_g * b_gamma        
        
        
        s = xsec_TT * f1gamma * self.acceptance
        
        return s
        
    def chi_square(self, g3g, g3gamma):
        """
        Calculates the chi2 test statistic for given couplings.
        (Eq. 8 in the paper: s^T V^{-1} s)

        Args:
            g3g (float): Dipole coupling to gluons
            g3gamma (float): Dipole coupling to photons

        Returns:
            float: The chi2 value.
        """
        s_vec = self.get_signal_vector(g3g, g3gamma)
        
        chi2_val = s_vec @ self.V_inv @ s_vec
        return chi2_val
    
    def find_contour(self, g3g_range=(0, 10), g3gamma_range=(0, 10), n_points=200):
        """
        Scans the (g3g, g3gamma) parameter space and finds the 95% CL contour.

        Args:
            g3g_range (tuple): (min, max) range for g3g scan.
            g3gamma_range (tuple): (min, max) range for g3gamma scan.
            n_points (int): Number of points along each axis.

        Returns:
            tuple: (X, Y, Z) meshgrid arrays and the chi2 values on the grid.
        """
        g3g_vals = np.linspace(g3g_range[0], g3g_range[1], n_points)
        g3gamma_vals = np.linspace(g3gamma_range[0], g3gamma_range[1], n_points)
        G3g, G3gamma = np.meshgrid(g3g_vals, g3gamma_vals)
        Chi2Grid = np.zeros_like(G3g)

        # Calculate chi2 for each point in the grid
        for i, g3g_val in enumerate(g3g_vals):
            for j, g3gamma_val in enumerate(g3gamma_vals):
                # Avoid (0, 0) to prevent division by zero in BR calculation
                if g3g_val == 0 and g3gamma_val == 0:
                    Chi2Grid[j, i] = 0
                else:
                    Chi2Grid[j, i] = self.chi_square(g3g_val, g3gamma_val)

        return G3g, G3gamma, Chi2Grid
    
    
if __name__ == "__main__":
    # compute_limit = ComputeLimit()
    # compute_limit.chi_square(500, "diff_xsec_photon_pt", 5, 5) 
    
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator

    import scipy
    from scipy.stats import chi2

    import numpy as np

    mass_points = [500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]
    colors = plt.cm.viridis(np.linspace(0, 1, len(mass_points)))
    plt.figure(figsize=(10, 8))
    chi2_68 = chi2.ppf(0.68, df=3)
    chi2_95 = chi2.ppf(0.95, df=3)
    legend_handles = []

    var = "diff_xsec_photon_pt" # "deltaphi_ll"
    for i, mass in enumerate(mass_points): # 
        # compute_limit = ComputeLimit()
        # compute_limit.find_contour(**dic[mass])
        compute_limit = ComputeLimit(mass, var)
        print(f"Start creating grid for mass {mass}...")
        X, Y, Z = compute_limit.find_contour(g3g_range=(0, 1000), g3gamma_range=(0, 1000), n_points=50)
        # X, Y, Z = 0.02 * X, 0.02 * Y, Z
        print("Grid created")
        plt.contourf(X, Y, Z, levels=[chi2_95, Z.max()], colors=[colors[i]], alpha=0.3)
        contour_95 = plt.contour(X, Y, Z, levels=[chi2_95], colors=[colors[i]], linewidths=1.5, linestyles='solid')
        # contour_68 = plt.contour(X, Y, Z, levels=[chi2_68], colors=[colors[i]], linewidths=2, linestyles='dashed')

        patch = Patch(color=colors[i], alpha=0.3, label=f'95% CL region, $m_T$ = {mass} [GeV]')
        # line = Line2D([0], [0], color=colors[i], lw=1.5, linestyle='solid', label=f'95% CL contor, $m_T$ = {mass} [GeV]')
        line = Line2D([0], [0], color=colors[i], lw=2, linestyle='dashed', label=f'68% CL contor, $m_T$ = {mass} [GeV]')

        legend_handles.extend([patch, line])

    # Horizontal line at y = 0
    plt.axhline(0, color='black', linestyle='--', linewidth=1)
    # Vertical line at x = 0
    plt.axvline(0, color='black', linestyle='--', linewidth=1)

    ax = plt.gca() 
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))

    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    ax.tick_params(axis='both', which='minor', length=4)
    ax.tick_params(axis='both', which='major', length=7)

    plt.xlabel('$g_{3g} (\mathrm{TeV^{-1}})$', fontsize = 14)
    plt.ylabel('$g_{3\\gamma} (\mathrm{TeV^{-1}})$', fontsize = 14)
    plt.title(f'95% and 68% CL Exclusion Contours in Coupling Space', fontsize = 18)
    # plt.grid()
    # plt.xlim(0, 1)
    # plt.ylim(0, 1)
    # plt.legend(handles=legend_elements)

    # plt.savefig(f'plots/limit_br_{self.mass}.svg', format='svg', dpi=1200)
    plt.legend(handles=legend_handles, ncol=1, loc='upper right', fontsize=8, bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(f"plots/limit_g.png")
    plt.show()