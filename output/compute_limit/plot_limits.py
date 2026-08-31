import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, MultipleLocator
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


def get_contour(mass, var, g3g_range, g3gamma_range, n_points, kfactor=1.4, force_recompute=False):
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
    np.savez(cache_file, X=X, Y=Y, Z=Z)
    print(f"Saved grid for mass {mass}, k={kfactor} to {cache_file}")
    return X, Y, Z


# ---------------------------------------------------------------------------
# Physics-validity overlays.
#
# compute_width_ratio(Mu4, gluonflag, gammaflag) is Gamma_T/M_T. Mu4 (the
# mass) is a single scalar per mass point, but gluonflag=g3g and
# gammaflag=g3gamma can be scalars OR numpy arrays (e.g. the raw
# (g3g, g3gamma) meshgrid) -- every operation below is plain +,-,*,/,**,
# which numpy broadcasts elementwise automatically, so passing the full grid
# in one call computes the whole surface at once. No np.vectorize needed.
#
# This is written in the same RAW coupling convention (g3g, g3gamma ~
# O(1-500)) as the rest of the codebase (gen_val's "g" array, the raw X, Y
# grid from find_contour) -- fvec3 and MHDO=Lambda[GeV] are already baked
# into the formula itself, so it should be called with the raw grid, not
# the 0.02-rescaled physical couplings.
# ---------------------------------------------------------------------------
def compute_width_ratio(Mu4, gluonflag, gammaflag):
    MHDO = 5000.0   # Lambda, in GeV (5 TeV)
    fvec3 = 0.1
    MT = 172.0
    MB = 4.7
    MW = 80.5

    gluonflag = np.asarray(gluonflag, dtype=float)
    gammaflag = np.asarray(gammaflag, dtype=float)

    # T -> t gamma partial width (depends on gammaflag = g3gamma)
    term_gamma = (((-MT**2 + Mu4**2) * ((48*fvec3**2*gammaflag**2*MT**4)/MHDO**2
                                         - (96*fvec3**2*gammaflag**2*MT**2*Mu4**2)/MHDO**2
                                         + (48*fvec3**2*gammaflag**2*Mu4**4)/MHDO**2))
                  / (96.*np.pi*abs(Mu4)**3))

    # T -> t g partial width (depends on gluonflag = g3g)
    term_gluon = (((-MT**2 + Mu4**2) * ((64*fvec3**2*gluonflag**2*MT**4)/MHDO**2
                                         - (128*fvec3**2*gluonflag**2*MT**2*Mu4**2)/MHDO**2
                                         + (64*fvec3**2*gluonflag**2*Mu4**4)/MHDO**2))
                  / (96.*np.pi*abs(Mu4)**3))

    # T -> b W partial width: does NOT depend on gluonflag/gammaflag at all
    # (this is the same physical piece as the fixed "width_wb" you already
    # store per mass in gen_val), so it's a single number for this mass.
    sqrt_arg = (MB**4 - 2*MB**2*Mu4**2 + Mu4**4 - 2*MB**2*MW**2 - 2*Mu4**2*MW**2 + MW**4)
    sqrt_arg = max(sqrt_arg, 0.0)  # guard against a below-threshold mass; should not trigger for Mu4 >> MB+MW
    term_bW = (((48*fvec3**2*MB**4)/MHDO**2 - (96*fvec3**2*MB**2*Mu4**2)/MHDO**2
                + (48*fvec3**2*Mu4**4)/MHDO**2 - (24*fvec3**2*MB**2*MW**2)/MHDO**2
                - (144*fvec3**2*MB*Mu4*MW**2)/MHDO**2 - (24*fvec3**2*Mu4**2*MW**2)/MHDO**2
                - (24*fvec3**2*MW**4)/MHDO**2) * np.sqrt(sqrt_arg)
               / (96.*np.pi*abs(Mu4)**3))

    width_T = term_gamma + term_gluon + term_bW
    return width_T / Mu4


def get_lambda_eff_grid(Xp, Yp):
    """
    Lambda_eff(c_tg, c_tgamma) = 1 / sqrt(c_tg^2 + c_tgamma^2), in TeV.
    Takes the PHYSICAL couplings Xp, Yp = 0.02*X, 0.02*Y (c = fvec3/Lambda *
    g = (0.1/5)*g), i.e. the same rescaling already used for the plot axes
    -- NOT the raw grid. (0,0) is replaced with a large finite sentinel
    instead of inf so matplotlib doesn't choke on it.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        lam = 1.0 / np.sqrt(Xp**2 + Yp**2)
    return np.ma.masked_invalid(lam)   # mask, don't fake a value

def get_lambda_eff_threshold(mass, output_dir="../output.coffea"):
    """Max M_TT (in TeV) found in the generated sample for this mass point."""
    try:
        output = load(output_dir)
        MTT = np.max(output["MTT_array"][f"Signal_{mass}"].value)
        return MTT
    except (FileNotFoundError, KeyError, AttributeError) as e:
        print(f"Warning: Could not load MTT for mass {mass}: {e}")
        return 1.0  # fallback value


# ---------------------------------------------------------------------------
# Main plotting code. Overlays, independently toggleable so you don't have
# to comment/uncomment blocks by hand:
#   - SHOW_KFACTOR_CONTOURS: exclusion contour(s), one per kfactor style
#   - SHOW_WIDTH_VALIDITY_BAND: hatched band where 0.1 < Gamma_T/M_T < 0.3,
#     plus the two boundary lines
#   - SHOW_EFT_VALIDITY_BOUNDARY: line where Lambda_eff(c_tg,c_tgamma)
#     equals the sample's max M_TT (EFT no longer valid beyond it)
# Wrapped in main() behind __name__ == "__main__" because find_contour()
# uses joblib multiprocessing internally.
# ---------------------------------------------------------------------------
def main():
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

    
    
    SHOW_KFACTOR_CONTOURS = False
    SHOW_WIDTH_VALIDITY_BAND_0p1 = False
    SHOW_WIDTH_VALIDITY_BAND_0p3 = True
    SHOW_EFT_VALIDITY_BOUNDARY = False

    # 'log' or 'linear' -- switches both x and y axes together.
    AXIS_SCALE = 'log'
    LOG_AXIS_MIN = 1e-5  # only used when AXIS_SCALE == 'log' (log axes can't show 0)

    mass_points = [500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]  # 

    # First entry is the reference case; further entries are outline-only
    # contours with the next linestyle, same color as their mass.
    
    fvec_over_lambda = 0.1/5000

    width_low, width_high = 0.1, 0.3

    for var in ["diff_xsec_photon_pt"]:  # , "deltaphi_ll"
        g3g_range = (0, 1100)
        g3gamma_range = (0, 1100)
        if var == "deltaphi_ll":
            n_points = 0
            chi2_68 = chi2.ppf(0.68, df=6)
            chi2_95 = chi2.ppf(0.95, df=6)
        else:
            n_points = 0
            chi2_68 = chi2.ppf(0.68, df=1)
            chi2_95 = chi2.ppf(0.95, df=1)

        colors = colors = ['#D55E00', '#E69F00', '#0072B2',  '#56B4E9', '#009E73', '#90EE90'
                            , '#BDB76B', '#F0E442', '#9400D3' , '#CC79A7' , '#708090', '#DA70D6',
                            '#800000']
        plt.figure(figsize=(10, 8))
        mass_legend_handles = []

        for i, mass in enumerate(mass_points):
            X = Y = None
            
            if SHOW_KFACTOR_CONTOURS:
                
                X, Y, Z = get_contour(mass, var, g3g_range, g3gamma_range, n_points, kfactor=1)
                # keep the unit conversion here (in the plotting stage), not in the cache,
                # so the cached grid stays reusable even if you change this factor later
                Xp, Yp = fvec_over_lambda * X, fvec_over_lambda * Y

                # plt.contourf(Yp, Xp, Z, levels=[chi2_95, Z.max()], colors=[colors[i]], alpha=0.3)
                plt.contour(Yp, Xp, Z, levels=[chi2_95], colors=[colors[i]],
                            linewidths=2, linestyles='dashed')
                # contour_68 = plt.contour(Yp, Xp, Z, levels=[chi2_68], colors=[colors[i]], linewidths=2, linestyles='dashed')
                

            X, Y, Z = get_contour(mass, var, g3g_range, g3gamma_range, n_points, kfactor=1.4)
            # non_zero_values = Z[Z != 0]
            # print(mass, np.min(non_zero_values))
            # keep the unit conversion here (in the plotting stage), not in the cache,
            # so the cached grid stays reusable even if you change this factor later
            Xp, Yp = fvec_over_lambda * X, fvec_over_lambda * Y

            # plt.contourf(Yp, Xp, Z, levels=[chi2_95, Z.max()], colors=[colors[i]], alpha=0.3)
            plt.contour(Yp, Xp, Z, levels=[chi2_95], colors=[colors[i]],
                        linewidths=2, linestyles='solid')
            # contour_68 = plt.contour(Yp, Xp, Z, levels=[chi2_68], colors=[colors[i]], linewidths=2, linestyles='dashed') 
            

            if SHOW_WIDTH_VALIDITY_BAND_0p1:
                # raw X, Y: compute_width_ratio's fvec3/MHDO already encode the
                # raw -> physical conversion internally
                width_ratio = compute_width_ratio(mass, X, Y)
                plt.contour(Yp, Xp, width_ratio, levels=[width_low], colors=[colors[i]],
                            linewidths=2, linestyles="dotted")
                
            if SHOW_WIDTH_VALIDITY_BAND_0p3:
                # raw X, Y: compute_width_ratio's fvec3/MHDO already encode the
                # raw -> physical conversion internally
                width_ratio = compute_width_ratio(mass, X, Y)
                plt.contour(Yp, Xp, width_ratio, levels=[width_high], colors=[colors[i]],
                            linewidths=2, linestyles="dotted")
                
            if SHOW_EFT_VALIDITY_BOUNDARY:
                mtt_max = get_lambda_eff_threshold(mass)
                lambda_eff_grid = get_lambda_eff_grid(Xp, Yp)
                plt.contour(Yp, Xp, lambda_eff_grid, levels=[mtt_max], colors=[colors[i]],
                            linewidths=2, linestyles='dashdot')

            patch = Line2D([0], [0], color=colors[i], lw=2, label=fr'$m_T = {mass:.0f}\ \mathrm{{GeV}}$')
            mass_legend_handles.append(patch)

        if SHOW_KFACTOR_CONTOURS:
            mass_legend_handles.append(Line2D([0], [0],
                                              color="black",
                                              lw=2,
                                              linestyle='dashed',
                                              label=fr'$K-factor = 1$')
                                      )

        if SHOW_WIDTH_VALIDITY_BAND_0p1:
            mass_legend_handles.append(Line2D([0], [0],
                                              color="black",
                                              lw=2,
                                              linestyle='dotted',
                                              label=fr'$\dfrac{{\Gamma_T}}{{m_T}} = 0.1$')
                                      )
            
        if SHOW_WIDTH_VALIDITY_BAND_0p3:
            mass_legend_handles.append(Line2D([0], [0],
                                              color="black",
                                              lw=2,
                                              linestyle='dotted',
                                              label=fr'$\dfrac{{\Gamma_T}}{{m_T}} = 0.3$')
                                      )

        if SHOW_EFT_VALIDITY_BOUNDARY:
            mass_legend_handles.append(Line2D([0], [0],
                                          color='black',
                                          lw=2,
                                          linestyle='dashdot',
                                          label=fr'$\Lambda_{{eff.}} = \mathrm{{max}}(m_{{T\bar{{T}}}})\ \mathrm{{GeV}}$'
                                         )
                                  )
            
        mass_legend = plt.legend(handles=mass_legend_handles, ncol=1, loc='lower right', fontsize=10,
                                 bbox_to_anchor=(0.97, 0.04),
                                 frameon=True, fancybox=True,
                                 framealpha=0.7, edgecolor="gray", borderpad=0.6)
        
        ax = plt.gca()
        ax.minorticks_on()
        if AXIS_SCALE == 'log':
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(left=LOG_AXIS_MIN, right=0.022)
            ax.set_ylim(bottom=LOG_AXIS_MIN, top=0.022)
            ax.tick_params(axis='both', which='minor', length=3)
            ax.tick_params(axis='both', which='major', length=7)
            ax.tick_params(axis='both', which='both' , top=True, right=True, direction='in', labelsize=15)
            mass_legend = plt.legend(handles=mass_legend_handles, ncol=1, loc='lower right', fontsize=10,
                                 bbox_to_anchor=(0.97, 0.04),
                                 frameon=True, fancybox=True,
                                 framealpha=0.7, edgecolor="gray", borderpad=0.6)
        elif AXIS_SCALE == 'linear':
            ax.set_xscale('linear')
            ax.set_yscale('linear')
            plt.axhline(0, color='black', linestyle='--', linewidth=1)
            plt.axvline(0, color='black', linestyle='--', linewidth=1)
            ax.tick_params(axis='both', which='minor', length=3)
            ax.tick_params(axis='both', which='major', length=7)
            ax.tick_params(axis='both', which='both' , top=True, right=True, direction='in', labelsize=15)
            mass_legend = plt.legend(handles=mass_legend_handles, ncol=1, loc='upper right', fontsize=10,
                                 bbox_to_anchor=(0.97, 0.96),
                                 frameon=True, fancybox=True,
                                 framealpha=0.7, edgecolor="gray", borderpad=0.6)
        else:
            raise ValueError(f"AXIS_SCALE must be 'log' or 'linear', got {AXIS_SCALE!r}")

        ax.set_ylabel(r"$c_{tg}\ [\mathrm{GeV}^{-1}]$"     , fontsize=18)
        ax.set_xlabel(r"$c_{t\gamma}\ [\mathrm{GeV}^{-1}]$", fontsize=18)
        ax.grid(True, which='major', linestyle='-', linewidth=0.7, alpha=0.3)
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
 
        
        ax.text(
            0.052, 1.03,
            #r'$\mathbf{95\%\ CL\ exclusion}$' +
            #r'$, 140\ \mathrm{fb}^{-1} \ $' +
            r'$\mathbf{t\bar{t}\gamma}$' +
            r' , $\mathbf{Untruncated \ Limits}$',
            transform=ax.transAxes,
            fontsize=16,
            ha='left',
            va='bottom',
            bbox=dict(
                boxstyle='round,pad=0.4',
                facecolor='white',
                edgecolor='gray',
                alpha=0.7
            )
        )
        
        if not (SHOW_WIDTH_VALIDITY_BAND_0p1 or SHOW_WIDTH_VALIDITY_BAND_0p3 or SHOW_EFT_VALIDITY_BOUNDARY):
            plt.savefig(f"plots/tta_untruncated_limits_{var}.png")
            plt.savefig(f"plots/tta_untruncated_limits_{var}.pdf")
        if SHOW_WIDTH_VALIDITY_BAND_0p1:
            plt.savefig(f"plots/tta_untruncated_limits_with_Gamma_over_m.png")
            plt.savefig(f"plots/tta_untruncated_limits_with_Gamma_over_m.pdf")
        if SHOW_WIDTH_VALIDITY_BAND_0p3:
            plt.savefig(f"plots/tta_untruncated_limits_with_Gamma_over_m_0p3.png")
            plt.savefig(f"plots/tta_untruncated_limits_with_Gamma_over_m_0p3.pdf")
        if SHOW_EFT_VALIDITY_BOUNDARY:
            plt.savefig(f"plots/tta_untruncated_limits_with_EFT_validity.png")
            plt.savefig(f"plots/tta_untruncated_limits_with_EFT_validity.pdf")

if __name__ == "__main__":
    main()