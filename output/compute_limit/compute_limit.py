from scipy.optimize import curve_fit, minimize
import scipy
from scipy.stats import chi2
from joblib import Parallel, delayed

from coffea.util import load
from helper import *

def _profile_worker(compute_limit_obj, g3g_val, g3gamma_val):
    """
    Module-level helper (not a method / not a closure) so it can be reliably
    pickled and shipped to worker processes by joblib. Each worker just
    re-runs the same nuisance-parameter minimization your original
    find_contour loop did serially, for one (g3g, g3gamma) grid point.
    """
    if g3g_val == 0 and g3gamma_val == 0:
        return 0.0
    return compute_limit_obj.profile_chi_square(g3g_val, g3gamma_val)


class ComputeLimit:
    def __init__(self, mass, var, in_br=False, kfactor=1.0, from_bin=3, hl_lhc = False):
        self.in_br = in_br
        self.mass = mass
        self.var = var
        self.kfactor = kfactor
        self.hl_lhc = hl_lhc
        self.from_bin = from_bin
        
        self.width_factor = width_prefactor(mass)
        self.width_wb = gen_val[f"Signal_{mass}"]["width_wb"]

        self.a0, self.a1, self.a2 = xsec_factors(mass)

        self.V_inv, self.data_minus_sm = HepDataParser().get_inverse_covariance_matrix_and_data_minus_sm(self.var, from_bin=from_bin, hl_lhc=hl_lhc)

        self.acceptance_gamma, self.acceptance_gammagamma, self.nuisance_deltas = generation_info(mass, var, from_bin=from_bin+1)
        self.nuisance_names = list(self.nuisance_deltas.keys())

    def compute_branching_ratios(self, g3g, g3gamma):

        width_tg = self.width_factor * g3g ** 2
        width_tgamma = (3 / 4) * self.width_factor * g3gamma ** 2

        Gamma_total = width_tg + width_tgamma + self.width_wb

        B_g = width_tg / Gamma_total
        B_gamma = width_tgamma / Gamma_total

        return B_g, B_gamma

    def compute_xsec(self, g3g, g3gamma):
        
        eft_eff1, eft_eff2 = compute_eft_eff(self.mass, g3g, g3gamma, self.from_bin+1)
        xsec1 = self.a0 * np.ones_like(eft_eff1) * self.kfactor + self.a1 * eft_eff1 * g3g ** 2 + self.a2 * eft_eff1 * g3g ** 4
        xsec2 = self.a0 * np.ones_like(eft_eff2) * self.kfactor + self.a1 * eft_eff2 * g3g ** 2 + self.a2 * eft_eff2 * g3g ** 4

        return xsec1 * 1000, xsec2* 1000
    
    def _compute_nominal_signal(self, g3g, g3gamma):
        xsec1, xsec2 = self.compute_xsec(g3g, g3gamma)      # expensive: compute_eft_eff lives here
        b_g, b_gamma = self.compute_branching_ratios(g3g, g3gamma)
        f1gamma, f2gamma = 2 * b_g * b_gamma, b_gamma ** 2
        return xsec1 * f1gamma * self.acceptance_gamma + xsec2 * f2gamma * self.acceptance_gammagamma

    def get_signal_vector(self, s0, theta=None):
        if theta is None:
            s = s0
        else:
            scale = np.ones_like(s0, dtype=float)
            for name, t in zip(self.nuisance_names, theta):
                scale *= (1.0 + t * self.nuisance_deltas[name])
            s = s0 * scale
        return self.data_minus_sm - s if not self.hl_lhc else s

    def chi_square(self, s0, theta=None):
        s_vec = self.get_signal_vector(s0, theta)
        val = s_vec @ self.V_inv @ s_vec
        return val + np.sum(np.square(theta)) if theta is not None else val

    def profile_chi_square(self, g3g, g3gamma):
        s0 = self._compute_nominal_signal(g3g, g3gamma)   # computed ONCE per grid point
        theta0 = np.zeros(len(self.nuisance_names))
        result = minimize(lambda th: self.chi_square(s0, th), theta0, method="BFGS")
        return result.fun

    def create_nonuniform_points(self, start, end, breakpoints, spacings):
        """
        Creates non-uniform grid points with different spacings in different regions.
        
        Args:
            start: Starting value
            end: Ending value
            breakpoints: List of breakpoints where spacing changes
            spacings: List of spacings for each region (must be one more than breakpoints)
        
        Returns:
            numpy array of grid points
        """
        points = []
        current = start
        
        # Handle each region
        for i in range(len(breakpoints) + 1):
            if i == 0:
                region_end = breakpoints[0]
            elif i == len(breakpoints):
                region_end = end
            else:
                region_end = breakpoints[i]
            
            spacing = spacings[i]
            
            # Generate points in this region
            while current <= region_end:
                points.append(current)
                current += spacing
                if current > region_end:
                    break
            
            # Set current to region_end for next region
            current = region_end
        
        # Ensure we have the end point
        if points[-1] != end:
            points.append(end)
        
        return np.array(points)
    
    def find_contour(self, g3g_range=(0, 10), g3gamma_range=(0, 10), n_points=200,
                      n_jobs=-1, verbose=5):
        """
        Scans the (g3g, g3gamma) parameter space and finds the 95% CL contour,
        using the profiled chi2 (nuisance parameters minimized out at every
        grid point). Grid points are independent of each other, so they are
        computed in parallel across CPU cores with joblib.

        Args:
            g3g_range (tuple): (min, max) range for g3g scan.
            g3gamma_range (tuple): (min, max) range for g3gamma scan.
            n_points (int): Number of points along each axis.
            n_jobs (int): Number of worker processes. -1 (default) uses all
                available CPU cores. Set to 1 to fall back to the original
                serial behaviour (e.g. for debugging).
            verbose (int): joblib verbosity (0 = silent, higher = more
                progress output).  gives periodic progress updates, useful
                given how long a full grid can take.

        Returns:5
            tuple: (X, Y, Z) meshgrid arrays and the chi2 values on the grid.
        """
        
        # g3g_vals = np.linspace(g3g_range[0], g3g_range[1], n_points)
        # g3gamma_vals = np.linspace(g3gamma_range[0], g3gamma_range[1], n_points)
        # G3g, G3gamma = np.meshgrid(g3g_vals, g3gamma_vals)
        
        
        breakpoints = [1e-59, 1e-58, 1e-57, 1e-56, 1e-55, 1e-54, 10]  # Where spacing changes , 0.001, 0.01, 0.1, 1, 20, 100
        spacings = [1e-60, 1e-59, 1e-58, 1e-57, 1e-56, 1e-55, 1, 1000] # , 1e-4, 0.0001, 0.001, 0.01, 0.1, 1, 3, 1000
        # breakpoints = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2]  # Where spacing changes , 0.001, 0.01, 0.1, 1, 20, 100
        # spacings = [5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 5e-2, 0.5, 3, 5] # , 1e-4, 0.0001, 0.001, 0.01, 0.1, 1, 3, 1000
        # breakpoints = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2]  # Where spacing changes , 0.001, 0.01, 0.1, 1, 20, 100
        # spacings = [5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 5e-2, 0.5, 3, 100]
        
        g3g_vals = self.create_nonuniform_points(
            g3g_range[0], g3g_range[1], breakpoints, spacings
        )
        g3gamma_vals = self.create_nonuniform_points(
            g3gamma_range[0], g3gamma_range[1], breakpoints, spacings
        )
        
        G3g, G3gamma = np.meshgrid(g3g_vals, g3gamma_vals)

        # Build the flat list of (i, j, g3g_val, g3gamma_val) tasks up front so
        # the parallel results (which come back in the same order they were
        # submitted) can be dropped straight back into the grid.
        tasks = [
            (i, j, g3g_val, g3gamma_val)
            for i, g3g_val in enumerate(g3g_vals)
            for j, g3gamma_val in enumerate(g3gamma_vals)
        ]

        results = Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(_profile_worker)(self, g3g_val, g3gamma_val)
            for (_, _, g3g_val, g3gamma_val) in tasks
        )

        Chi2Grid = np.zeros_like(G3g)
        for (i, j, _, _), chi2_val in zip(tasks, results):
            Chi2Grid[j, i] = chi2_val

        return G3g, G3gamma, Chi2Grid