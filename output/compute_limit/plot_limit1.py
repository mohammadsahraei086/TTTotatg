import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scipy.stats import chi2
from compute_limit import ComputeLimit
from coffea.util import load

# ---------------------------------------------------------------------------
# Caching: each (mass, var, ranges, n_points, kfactor) grid is computed once
# and saved to disk. As long as those inputs don't change, re-running this
# script just loads the .npz files instead of re-running find_contour.
# Delete the cache file (or pass force_recompute=True) whenever the
# underlying physics/chi2 code actually changes.
# NOTE: np.savez below is commented out, same as you had it -- caching is
# effectively disabled right now (it will still try to load a cache file if
# one happens to exist, but nothing new gets saved). Uncomment it if that
# was just left over from debugging.
# ---------------------------------------------------------------------------
CACHE_DIR = "contour_cache"
os.makedirs(CACHE_DIR, exist_ok=True)


def get_contour(mass, var, g3g_range, g3gamma_range, n_points, kfactor=1.0, force_recompute=False):
    cache_file = os.path.join(
        CACHE_DIR,
        f"contour_{var}_m{mass}_n{n_points}_k{kfactor:.3f}_"
        f"g3g{g3g_range[0]}-{g3g_range[1]}_g3gamma{g3gamma_range[0]}-{g3gamma_range[1]}.npz",
    )

    if os.path.exists(cache_file) and not force_recompute:
        data = np.load(cache_file)
        print(f"Loaded cached grid for mass {mass}, k={kfactor} from {cache_file}")
        return data["X"], data["Y"], data["Z"]

    print(f"No cache found for mass {mass}, k={kfactor} -- computing (this is the slow part)...")
    compute_limit = ComputeLimit(mass, var, kfactor=kfactor)
    X, Y, Z = compute_limit.find_contour(
        g3g_range=g3g_range, g3gamma_range=g3gamma_range, n_points=n_points
    )
    # np.savez(cache_file, X=X, Y=Y, Z=Z)
    print(f"Saved grid for mass {mass}, k={kfactor} to {cache_file}")
    return X, Y, Z

def compute_width(Mu4, gluonflag, gammaflag):
    MHDO = 5000
    fvec3 = 0.1
    MT = 172
    MB = 4.7
    MW = 80.5
    width_T = (((-MT**2 + Mu4**2)*((48*fvec3**2*gammaflag**2*MT**4)/MHDO**2 - (96*fvec3**2*gammaflag**2*MT**2*Mu4**2)/MHDO**2 +
                                   (48*fvec3**2*gammaflag**2*Mu4**4)/MHDO**2))/(96.*cmath.pi*abs(Mu4)**3) + 
               ((-MT**2 + Mu4**2)*((64*fvec3**2*gluonflag**2*MT**4)/MHDO**2 - (128*fvec3**2*gluonflag**2*MT**2*Mu4**2)/MHDO**2 + 
                                   (64*fvec3**2*gluonflag**2*Mu4**4)/MHDO**2))/(96.*cmath.pi*abs(Mu4)**3) + 
               (((48*fvec3**2*MB**4)/MHDO**2 - (96*fvec3**2*MB**2*Mu4**2)/MHDO**2 + (48*fvec3**2*Mu4**4)/MHDO**2 - (24*fvec3**2*MB**2*MW**2)/MHDO**2 - 
                 (144*fvec3**2*MB*Mu4*MW**2)/MHDO**2 - (24*fvec3**2*Mu4**2*MW**2)/MHDO**2 - (24*fvec3**2*MW**4)/MHDO**2)*cmath.sqrt(MB**4 - 2*MB**2*Mu4**2 + Mu4**4 - 2*MB**2*MW**2 -
                                                                                                                                    2*Mu4**2*MW**2 + MW**4))/(96.*cmath.pi*abs(Mu4)**3))
    return width_T/Mu4
    
def get_lambda_eff(mass, output_dir="../output.coffea"):
    output = load(output_dir)
    MTT = np.max(output["MTT_array"][f"Signal_{mass}"].value)
    MTT_TeV = MTT/1000
    
    return MTT_TeV
    
    
    
    
# ---------------------------------------------------------------------------
# Everything below here is your restructured plotting code, unchanged except
# for the k-factor overlay: for every mass, both k-factor variants are drawn
# in the same color, distinguished by line style. Wrapped in main() behind
# __name__ == "__main__" because find_contour() uses joblib multiprocessing
# internally.
# ---------------------------------------------------------------------------
def main():
    mass_points = [500, 2000]  # , 2250, 2500, 2750, 3000 ,  750, 1000, 1250, 1500, 1750

    # First entry is the reference case; further entries are outline-only
    # contours with the next linestyle, same color as their mass.
    kfactor_styles = [
        (1.0, 'solid'),
        (1.4, 'dashed'),
    ]

    for var in ["diff_xsec_photon_pt"]:  # , "deltaphi_ll"
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
        mass_legend_handles = []

        for i, mass in enumerate(mass_points):
            for k_index, (kfactor, linestyle) in enumerate(kfactor_styles):
                X, Y, Z = get_contour(mass, var, g3g_range, g3gamma_range, n_points, kfactor=kfactor)
                # keep the unit conversion here (in the plotting stage), not in the cache,
                # so the cached grid stays reusable even if you change this factor later
                Xp, Yp = 0.02 * X, 0.02 * Y

                # plt.contourf(Yp, Xp, Z, levels=[chi2_95, Z.max()], colors=[colors[i]], alpha=0.3)
                plt.contour(Yp, Xp, Z, levels=[chi2_95], colors=[colors[i]],
                            linewidths=1.5, linestyles=linestyle)
                # contour_68 = plt.contour(Yp, Xp, Z, levels=[chi2_68], colors=[colors[i]], linewidths=2, linestyles='dashed')

            patch = Patch(color=colors[i], alpha=0.3, label=f'$m_T$ = {mass} [GeV]')
            mass_legend_handles.append(patch)

        kfactor_legend_handles = [
            Line2D([0], [0], color='black', lw=1.5, linestyle=linestyle,
                   label=f'k-factor = {kfactor}')
            for kfactor, linestyle in kfactor_styles
        ]

        plt.axhline(0, color='black', linestyle='--', linewidth=1)
        plt.axvline(0, color='black', linestyle='--', linewidth=1)
        ax = plt.gca()
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.tick_params(axis='both', which='minor', length=4)
        ax.tick_params(axis='both', which='major', length=7)
        plt.xlabel('$c_{t\\gamma} [\\mathrm{TeV^{-1}}]$', fontsize=14)
        plt.ylabel('$c_{tg} [\\mathrm{TeV^{-1}}]$', fontsize=14)
        plt.title(f'95% CL Exclusion Contours in Coupling Space ("{var}")', fontsize=18)

        mass_legend = plt.legend(handles=mass_legend_handles, ncol=1, loc='upper right',
                                  fontsize=8, bbox_to_anchor=(1, 1), title='Mass')
        ax.add_artist(mass_legend)
        plt.legend(handles=kfactor_legend_handles, ncol=1, loc='upper right', fontsize=8,
                   bbox_to_anchor=(1, 1 - 0.03 * (len(mass_legend_handles) + 2)), title='K-factor')

        plt.grid()
        plt.tight_layout()
        plt.savefig(f"plots/limit_g_{var}.png")


if __name__ == "__main__":
    main()