from scipy.optimize import curve_fit, minimize
import scipy
from scipy.stats import chi2
import numpy as np
from coffea.util import load
from ds_helper import *

class ComputeLimit:
    def __init__(self, mass, var, in_br=False):
        self.in_br = in_br
        self.mass = mass
        self.var = var
        
        self.width_factor = width_prefactor(mass)
        self.width_wb = gen_val[f"Signal_{mass}"]["width_wb"]
        
        self.a0, self.a1, self.a2 = xsec_factors(mass)
        
        self.acceptance = compute_acceptance(mass, var)
        
        self.V_inv = HepDataParser().get_inverse_covariance_matrix(self.var)
        
        # Load generation uncertainties (PDF, muR, muF)
        # These should be in a similar format as your gen_val dictionary
        self.load_generation_uncertainties()
        
    def load_generation_uncertainties(self):
        """
        Load generation uncertainties (PDF, scale variations) from file.
        You need to provide these values for your specific mass point.
        """
        # Example structure - you need to fill these with actual values
        # These should be fractional uncertainties (relative) for each photon pt bin
        # For example, if you have 6 bins:
        n_bins = 6
        
        # PDF uncertainty - you might have it from the PDF set variation
        self.pdf_uncertainty = np.array([0.02, 0.025, 0.03, 0.035, 0.04, 0.045])  # Example values
        
        # muR uncertainty - scale variation
        self.muR_uncertainty = np.array([0.015, 0.02, 0.025, 0.03, 0.035, 0.04])  # Example values
        
        # muF uncertainty - scale variation  
        self.muF_uncertainty = np.array([0.015, 0.02, 0.025, 0.03, 0.035, 0.04])  # Example values
        
        # Combine into a diagonal covariance matrix for generation uncertainties
        # Assuming they're independent, but you could add correlations if known
        self.gen_cov_diag = self.pdf_uncertainty**2 + self.muR_uncertainty**2 + self.muF_uncertainty**2
        
    def compute_branching_ratios(self, g3g, g3gamma):
        width_tg = self.width_factor * g3g ** 2
        width_tgamma = (3/4) * self.width_factor * g3gamma ** 2
        
        Gamma_total = width_tg + width_tgamma + self.width_wb

        B_g = width_tg / Gamma_total
        B_gamma = width_tgamma / Gamma_total

        return B_g, B_gamma
    
    def compute_xsec(self, g3g):
        return self.a0 + self.a1 * g3g ** 2 + self.a2 * g3g ** 4
        
    def get_signal_vector(self, g3g, g3gamma, nuisance_params=None):
        """
        Get signal vector with optional nuisance parameters.
        
        Args:
            g3g (float): Dipole coupling to gluons
            g3gamma (float): Dipole coupling to photons
            nuisance_params (np.array): Nuisance parameters for generation uncertainties
        
        Returns:
            np.array: Signal vector
        """
        xsec_TT = self.compute_xsec(g3g)
        b_g, b_gamma = self.compute_branching_ratios(g3g, g3gamma)
        f1gamma = 2 * b_g * b_gamma        
        
        s = xsec_TT * f1gamma * self.acceptance
        
        # Apply nuisance parameters if provided
        if nuisance_params is not None:
            # Nuisance parameters multiply the signal (relative variation)
            # You might want to add them as additive or multiplicative factors
            s = s * (1 + nuisance_params)
        
        return s
    
    def chi_square_with_nuisance(self, g3g, g3gamma, nuisance_params):
        """
        Calculates chi2 with nuisance parameters.
        
        Chi2 = (s(θ) - s_obs)^T V^{-1} (s(θ) - s_obs) + θ^T V_gen^{-1} θ
        
        Where θ are nuisance parameters, V_gen is the generation covariance matrix
        """
        s_vec = self.get_signal_vector(g3g, g3gamma, nuisance_params)
        
        # Main chi2 term
        chi2_main = s_vec @ self.V_inv @ s_vec
        
        # Penalty term for nuisance parameters
        # Assuming nuisance parameters are independent and have diagonal covariance
        chi2_penalty = np.sum(nuisance_params**2 / self.gen_cov_diag)
        
        return chi2_main + chi2_penalty
    
    def profile_chi_square(self, g3g, g3gamma, initial_nuisance=None):
        """
        Profile over nuisance parameters: minimize chi2 with respect to nuisance parameters.
        """
        if initial_nuisance is None:
            initial_nuisance = np.zeros(len(self.gen_cov_diag))
        
        # Define the negative log-likelihood (or chi2) to minimize
        def objective(nuisance):
            return self.chi_square_with_nuisance(g3g, g3gamma, nuisance)
        
        # Minimize with respect to nuisance parameters
        result = minimize(objective, initial_nuisance, method='Nelder-Mead')
        
        return result.fun, result.x
    
    def chi_square_profiled(self, g3g, g3gamma):
        """
        Calculate the profiled chi2 (minimized over nuisance parameters)
        """
        chi2_profiled, _ = self.profile_chi_square(g3g, g3gamma)
        return chi2_profiled
    
    def find_contour(self, g3g_range=(0, 10), g3gamma_range=(0, 10), n_points=200):
        """
        Scans the (g3g, g3gamma) parameter space and finds the 95% CL contour
        with nuisance parameters profiled.
        """
        g3g_vals = np.linspace(g3g_range[0], g3g_range[1], n_points)
        g3gamma_vals = np.linspace(g3gamma_range[0], g3gamma_range[1], n_points)
        G3g, G3gamma = np.meshgrid(g3g_vals, g3gamma_vals)
        Chi2Grid = np.zeros_like(G3g)
        
        # Dictionary to store nuisance parameters at minimum if needed
        self.nuisance_min = {}
        
        # Calculate profiled chi2 for each point in the grid
        for i, g3g_val in enumerate(g3g_vals):
            for j, g3gamma_val in enumerate(g3gamma_vals):
                if g3g_val == 0 and g3gamma_val == 0:
                    Chi2Grid[j, i] = 0
                else:
                    chi2_profiled, nuisance_opt = self.profile_chi_square(g3g_val, g3gamma_val)
                    Chi2Grid[j, i] = chi2_profiled
                    # Store nuisance parameters for diagnostic purposes
                    self.nuisance_min[f"({g3g_val},{g3gamma_val})"] = nuisance_opt

        return G3g, G3gamma, Chi2Grid