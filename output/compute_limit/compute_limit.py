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
    def __init__(self, mass, var, in_br=False):
        self.in_br = in_br
        self.mass = mass
        self.var = var

        self.width_factor = width_prefactor(mass)
        self.width_wb = gen_val[f"Signal_{mass}"]["width_wb"]

        self.a0, self.a1, self.a2 = xsec_factors(mass)
        self.kfactor = 1.4

        self.V_inv = HepDataParser().get_inverse_covariance_matrix(self.var)

        self.acceptance, self.nuisance_deltas = generation_info(mass, var)
        self.nuisance_names = list(self.nuisance_deltas.keys())

    def compute_branching_ratios(self, g3g, g3gamma):

        width_tg = self.width_factor * g3g ** 2
        width_tgamma = (3 / 4) * self.width_factor * g3gamma ** 2

        Gamma_total = width_tg + width_tgamma + self.width_wb

        B_g = width_tg / Gamma_total
        B_gamma = width_tgamma / Gamma_total

        return B_g, B_gamma

    def compute_xsec(self, g3g):

        return self.a0 * self.kfactor + self.a1 * g3g ** 2 + self.a2 * g3g ** 4

    def get_signal_vector(self, g3g, g3gamma, theta=None):
        """
        theta: array-like, one entry per nuisance parameter in
        self.nuisance_names (currently ["pdf", "mur", "muf"]), given in
        units of standard deviations. theta=None (or all zeros) reproduces
        the original, nominal signal prediction.

        Each nuisance parameter shifts every bin coherently (fully
        correlated across bins, which is the standard assumption for a
        single PDF/scale variation), with a bin-dependent size given by
        self.nuisance_deltas[name]:

            s_i(theta) = s_i^nominal * prod_k (1 + theta_k * delta_{k,i})
        """
        xsec_TT = self.compute_xsec(g3g)
        b_g, b_gamma = self.compute_branching_ratios(g3g, g3gamma)
        f1gamma = 2 * b_g * b_gamma

        s = xsec_TT * f1gamma * self.acceptance

        if theta is not None:
            scale = np.ones_like(s, dtype=float)
            for name, t in zip(self.nuisance_names, theta):
                scale = scale * (1.0 + t * self.nuisance_deltas[name])
            s = s * scale

        return s

    def chi_square(self, g3g, g3gamma, theta=None):
        """
        Calculates the chi2 test statistic for given couplings and (optionally)
        given nuisance-parameter values:

            chi2(g3g, g3gamma, theta) = s(theta)^T V_inv s(theta) + sum_k theta_k^2

        The sum_k theta_k^2 term is the Gaussian constraint ("penalty") that
        keeps each nuisance parameter close to its nominal value of 0 unless
        the fit to data actually prefers otherwise. With theta=None this is
        identical to the original chi2 (Eq. 8 in the paper).

        Args:
            g3g (float): Dipole coupling to gluons
            g3gamma (float): Dipole coupling to photons
            theta (array-like, optional): nuisance parameter values

        Returns:
            float: The chi2 value.
        """
        s_vec = self.get_signal_vector(g3g, g3gamma, theta)
        chi2_val = s_vec @ self.V_inv @ s_vec

        if theta is not None:
            chi2_val = chi2_val + np.sum(np.square(theta))

        return chi2_val

    def profile_chi_square(self, g3g, g3gamma):
        """
        Profiles out the nuisance parameters at fixed (g3g, g3gamma):

            chi2_prof(g3g, g3gamma) = min_theta chi2(g3g, g3gamma, theta)

        This is the quantity that should now be compared to the chi2_95
        (or chi2_68) threshold -- the number of degrees of freedom used for
        that threshold does NOT change: profiling nuisance parameters out
        does not add degrees of freedom to the test of (g3g, g3gamma), it
        only lets each nuisance pull the prediction within its constraint
        while doing so.

        Returns:
            float: minimized chi2 value.
        """
        n_nuisance = len(self.nuisance_names)
        theta0 = np.zeros(n_nuisance)

        result = minimize(
            lambda theta: self.chi_square(g3g, g3gamma, theta),
            theta0,
            method="BFGS",
        )

        return result.fun

    # def find_contour(self, g3g_range=(0, 10), g3gamma_range=(0, 10), n_points=200):
    #     """
    #     Scans the (g3g, g3gamma) parameter space and finds the 95% CL contour,
    #     using the profiled chi2 (nuisance parameters minimized out at every
    #     grid point).

    #     Args:
    #         g3g_range (tuple): (min, max) range for g3g scan.
    #         g3gamma_range (tuple): (min, max) range for g3gamma scan.
    #         n_points (int): Number of points along each axis.

    #     Returns:
    #         tuple: (X, Y, Z) meshgrid arrays and the chi2 values on the grid.
    #     """
    #     g3g_vals = np.linspace(g3g_range[0], g3g_range[1], n_points)
    #     g3gamma_vals = np.linspace(g3gamma_range[0], g3gamma_range[1], n_points)
    #     G3g, G3gamma = np.meshgrid(g3g_vals, g3gamma_vals)
    #     Chi2Grid = np.zeros_like(G3g)

    #     # Calculate profiled chi2 for each point in the grid
    #     for i, g3g_val in enumerate(g3g_vals):
    #         for j, g3gamma_val in enumerate(g3gamma_vals):
    #             # Avoid (0, 0) to prevent division by zero in BR calculation
    #             if g3g_val == 0 and g3gamma_val == 0:
    #                 Chi2Grid[j, i] = 0
    #             else:
    #                 Chi2Grid[j, i] = self.profile_chi_square(g3g_val, g3gamma_val)

    #     return G3g, G3gamma, Chi2Grid


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
                progress output). 5 gives periodic progress updates, useful
                given how long a full grid can take.

        Returns:
            tuple: (X, Y, Z) meshgrid arrays and the chi2 values on the grid.
        """
        g3g_vals = np.linspace(g3g_range[0], g3g_range[1], n_points)
        g3gamma_vals = np.linspace(g3gamma_range[0], g3gamma_range[1], n_points)
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