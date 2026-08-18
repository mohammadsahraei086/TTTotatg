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
