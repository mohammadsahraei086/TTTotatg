import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scipy.stats import chi2
from compute_limit import ComputeLimit

# ---------------------------------------------------------------------------
# Caching: each (mass, var, ranges, n_points) grid is computed once and saved
# to disk. As long as those inputs don't change, re-running this script just
# loads the .npz files instead of re-running find_contour. Delete the cache
# file (or pass force_recompute=True) whenever the underlying physics/chi2
# code actually changes.
# ---------------------------------------------------------------------------
CACHE_DIR = "contour_cache"
os.makedirs(CACHE_DIR, exist_ok=True)


def get_contour(mass, var, g3g_range, g3gamma_range, n_points, force_recompute=False):
    cache_file = os.path.join(
        CACHE_DIR,
        f"contour_{var}_m{mass}_n{n_points}_"
        f"g3g{g3g_range[0]}-{g3g_range[1]}_g3gamma{g3gamma_range[0]}-{g3gamma_range[1]}.npz",
    )

    if os.path.exists(cache_file) and not force_recompute:
        data = np.load(cache_file)
        print(f"Loaded cached grid for mass {mass} from {cache_file}")
        return data["X"], data["Y"], data["Z"]

    print(f"No cache found for mass {mass} -- computing (this is the slow part)...")
    compute_limit = ComputeLimit(mass, var)
    X, Y, Z = compute_limit.find_contour(
        g3g_range=g3g_range, g3gamma_range=g3gamma_range, n_points=n_points
    )
    # np.savez(cache_file, X=X, Y=Y, Z=Z)
    print(f"Saved grid for mass {mass} to {cache_file}")
    return X, Y, Z


# ---------------------------------------------------------------------------
# Everything below here is your original plotting code, unchanged except
# that X, Y, Z now come from get_contour() instead of a fresh find_contour()
# call, and it's wrapped in main() behind an __name__ == "__main__" guard.
# That guard is required because find_contour() now uses joblib
# multiprocessing internally: without it, each worker process could end up
# re-importing and re-executing this whole script. Tweak labels/colors/
# titles freely inside main() -- reruns from cache will be fast.
# ---------------------------------------------------------------------------
def main():
    mass_points = [500, 2000]  # , 2250, 2500, 2750, 3000 ,  750, 1000, 1250, 1500, 1750
    for var in ["diff_xsec_photon_pt"]: # , "deltaphi_ll"
        g3g_range = (0, 500)
        g3gamma_range = (0, 250)
        if var == "deltaphi_ll":
            n_points = 200
            chi2_68 = chi2.ppf(0.68, df=6)
            chi2_95 = chi2.ppf(0.95, df=6)
        else:
            n_points = 50
            chi2_68 = chi2.ppf(0.68, df=4)
            chi2_95 = chi2.ppf(0.95, df=4)
    
        colors = plt.cm.viridis(np.linspace(0, 1, len(mass_points)))
        plt.figure(figsize=(10, 8))
        legend_handles = []
    
        for i, mass in enumerate(mass_points):
            X, Y, Z = get_contour(mass, var, g3g_range, g3gamma_range, n_points)
            # keep the unit conversion here (in the plotting stage), not in the cache,
            # so the cached grid stays reusable even if you change this factor later
            X, Y = 0.02 * X, 0.02 * Y
    
            # plt.contourf(Y, X, Z, levels=[chi2_95, Z.max()], colors=[colors[i]], alpha=0.3)
            contour_95 = plt.contour(Y, X, Z, levels=[chi2_95], colors=[colors[i]], linewidths=1.5, linestyles='solid')
            # contour_68 = plt.contour(X, Y, Z, levels=[chi2_68], colors=[colors[i]], linewidths=2, linestyles='dashed')
            patch = Patch(color=colors[i], alpha=0.3, label=f'$m_T$ = {mass} [GeV]')
            # line = Line2D([0], [0], color=colors[i], lw=2, linestyle='dashed', label=f'68% CL contor, $m_T$ = {mass} [GeV]')
            legend_handles.extend([patch])  # , line
    
        plt.axhline(0, color='black', linestyle='--', linewidth=1)
        plt.axvline(0, color='black', linestyle='--', linewidth=1)
        ax = plt.gca()
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.tick_params(axis='both', which='minor', length=4)
        ax.tick_params(axis='both', which='major', length=7)
        plt.xlabel('$c_{t\gamma} [\\mathrm{TeV^{-1}}]$', fontsize=14)
        plt.ylabel('$c_{tg} [\\mathrm{TeV^{-1}}]$', fontsize=14)
        plt.title(f'95% CL Exclusion Contours in Coupling Space ("{var}")', fontsize=18)
        plt.legend(handles=legend_handles, ncol=1, loc='upper right', fontsize=8, bbox_to_anchor=(1, 1))
        plt.grid()
        plt.tight_layout()
        plt.savefig(f"plots/limit_g_{var}.png")


if __name__ == "__main__":
    main()
