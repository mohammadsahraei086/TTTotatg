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
    return np.where(np.isfinite(lam), lam, 1e6)


def get_lambda_eff_threshold(mass, output_dir="../output.coffea"):
    """Max M_TT (in TeV) found in the generated sample for this mass point."""
    output = load(output_dir)
    MTT = np.max(output["MTT_array"][f"Signal_{mass}"].value)
    return MTT / 1000


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
    SHOW_KFACTOR_CONTOURS = False
    SHOW_WIDTH_VALIDITY_BAND = True
    SHOW_WIDTH_VALIDITY_BAND_0p1 = False
    SHOW_WIDTH_VALIDITY_BAND_0p3 = False
    SHOW_EFT_VALIDITY_BOUNDARY = False

    # 'log' or 'linear' -- switches both x and y axes together.
    AXIS_SCALE = 'log'
    LOG_AXIS_MIN = 1e-2  # only used when AXIS_SCALE == 'log' (log axes can't show 0)

    mass_points = [500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]  #  ,  

    # First entry is the reference case; further entries are outline-only
    # contours with the next linestyle, same color as their mass.
    kfactor_styles = [
        (1.0, 'dashed'),
        (1.4, 'solid'),
    ]

    width_low, width_high = 0.1, 0.3
    width_low_style = ':'
    width_high_style = '-.'
    width_hatch = '///'

    eft_style = (0, (3, 1, 1, 1))  # densely dash-dotted

    for var in ["diff_xsec_photon_pt"]:  # , "deltaphi_ll"
        g3g_range = (0, 1000)
        g3gamma_range = (0, 1000)
        if var == "deltaphi_ll":
            n_points = 200
            chi2_68 = chi2.ppf(0.68, df=6)
            chi2_95 = chi2.ppf(0.95, df=6)
        else:
            n_points = 300
            chi2_68 = chi2.ppf(0.68, df=4)
            chi2_95 = chi2.ppf(0.95, df=4)

        colors = plt.cm.viridis(np.linspace(0, 1, len(mass_points)))
        plt.figure(figsize=(10, 8))
        mass_legend_handles = []

        for i, mass in enumerate(mass_points):
            X = Y = None

            X, Y, Z = get_contour(mass, var, g3g_range, g3gamma_range, n_points, kfactor=1.4)
            # non_zero_values = Z[Z != 0]
            # print(mass, np.min(non_zero_values))
            # keep the unit conversion here (in the plotting stage), not in the cache,
            # so the cached grid stays reusable even if you change this factor later
            Xp, Yp = 0.02 * X, 0.02 * Y

            # plt.contourf(Yp, Xp, Z, levels=[chi2_95, Z.max()], colors=[colors[i]], alpha=0.3)
            plt.contour(Yp, Xp, Z, levels=[chi2_95], colors=[colors[i]],
                        linewidths=1.5)
            # contour_68 = plt.contour(Yp, Xp, Z, levels=[chi2_68], colors=[colors[i]], linewidths=2, linestyles='dashed')
                    
            if SHOW_KFACTOR_CONTOURS:
                
                X, Y, Z = get_contour(mass, var, g3g_range, g3gamma_range, n_points, kfactor=1)
                # keep the unit conversion here (in the plotting stage), not in the cache,
                # so the cached grid stays reusable even if you change this factor later
                Xp, Yp = 0.02 * X, 0.02 * Y

                # plt.contourf(Yp, Xp, Z, levels=[chi2_95, Z.max()], colors=[colors[i]], alpha=0.3)
                plt.contour(Yp, Xp, Z, levels=[chi2_95], colors=[colors[i]],
                            linewidths=1.5, linestyles='dashed')
                # contour_68 = plt.contour(Yp, Xp, Z, levels=[chi2_68], colors=[colors[i]], linewidths=2, linestyles='dashed')

            if SHOW_WIDTH_VALIDITY_BAND or SHOW_EFT_VALIDITY_BOUNDARY:
                if X is None:
                    # need a raw grid even if kfactor contours are switched off
                    X, Y, _ = get_contour(mass, var, g3g_range, g3gamma_range, n_points, kfactor=1.4)
                Xp, Yp = 0.02 * X, 0.02 * Y

            if SHOW_WIDTH_VALIDITY_BAND_0p1:
                # raw X, Y: compute_width_ratio's fvec3/MHDO already encode the
                # raw -> physical conversion internally
                width_ratio = compute_width_ratio(mass, X, Y)
                plt.contour(Yp, Xp, width_ratio, levels=[width_low], colors=[colors[i]],
                            linewidths=1.5, linestyles=width_low_style)
                
            if SHOW_WIDTH_VALIDITY_BAND_0p3:
                # raw X, Y: compute_width_ratio's fvec3/MHDO already encode the
                # raw -> physical conversion internally
                width_ratio = compute_width_ratio(mass, X, Y)
                plt.contour(Yp, Xp, width_ratio, levels=[width_high], colors=[colors[i]],
                            linewidths=1.5, linestyles=width_high_style)
                
            if SHOW_WIDTH_VALIDITY_BAND:
                # raw X, Y: compute_width_ratio's fvec3/MHDO already encode the
                # raw -> physical conversion internally
                width_ratio = compute_width_ratio(mass, X, Y)
                cs = plt.contourf(Yp, Xp, width_ratio, levels=[width_low, width_high],
                                   colors='none', hatches=[width_hatch])
                for coll in cs.collections:
                    coll.set_edgecolor(colors[i])
                    coll.set_linewidth(0.0)
                plt.contour(Yp, Xp, width_ratio, levels=[width_low], colors=[colors[i]],
                            linewidths=1.0, linestyles=width_low_style)
                plt.contour(Yp, Xp, width_ratio, levels=[width_high], colors=[colors[i]],
                            linewidths=1.0, linestyles=width_high_style)

            if SHOW_EFT_VALIDITY_BOUNDARY:
                # physical Xp, Yp: Lambda_eff is defined in terms of the actual
                # dimensionful Wilson coefficients, squared
                mtt_max = get_lambda_eff_threshold(mass)
                lambda_eff_grid = get_lambda_eff_grid(Xp, Yp)
                plt.contour(Yp, Xp, lambda_eff_grid, levels=[mtt_max], colors=[colors[i]],
                            linewidths=1.2, linestyles=[eft_style])

            patch = Patch(color=colors[i], alpha=0.3, label=f'$m_T$ = {mass} [GeV]')
            mass_legend_handles.append(patch)

        legend_blocks = []

        if SHOW_KFACTOR_CONTOURS:
            kfactor_legend_handles = [
                Line2D([0], [0], color='black', lw=1.5, linestyle=linestyle,
                       label=f'k-factor = {kfactor}')
                for kfactor, linestyle in kfactor_styles
            ]
            legend_blocks.append(('K-factor', kfactor_legend_handles))

        if SHOW_WIDTH_VALIDITY_BAND_0p1:
            width_legend_handles = [
                Line2D([0], [0], color='black', lw=1.3, linestyle=width_low_style,
                       label=f'$\\Gamma_T/M_T$ = {width_low}'),
            ]
            legend_blocks.append(('Width validity', width_legend_handles))
            
        if SHOW_WIDTH_VALIDITY_BAND_0p3:
            width_legend_handles = [
                Line2D([0], [0], color='black', lw=1.3, linestyle=width_high_style,
                       label=f'$\\Gamma_T/M_T$ = {width_high}'),
            ]
            legend_blocks.append(('Width validity', width_legend_handles))
        
        if SHOW_WIDTH_VALIDITY_BAND:
            width_legend_handles = [
                Line2D([0], [0], color='black', lw=1.0, linestyle=width_low_style,
                       label=f'$\\Gamma_T/M_T$ = {width_low}'),
                Line2D([0], [0], color='black', lw=1.0, linestyle=width_high_style,
                       label=f'$\\Gamma_T/M_T$ = {width_high}'),
                Patch(facecolor='none', edgecolor='black', hatch=width_hatch,
                      label=f'{width_low} < $\\Gamma_T/M_T$ < {width_high}'),
            ]
            legend_blocks.append(('Width validity', width_legend_handles))

        if SHOW_EFT_VALIDITY_BOUNDARY:
            eft_legend_handles = [
                Line2D([0], [0], color='black', lw=1.2, linestyle=eft_style,
                       label='$\\Lambda_{eff} = M_{TT}^{max}$'),
            ]
            legend_blocks.append(('EFT validity', eft_legend_handles))

        ax = plt.gca()
        if AXIS_SCALE == 'log':
            # Log axes can't include 0 (log(0) is undefined), so there's no
            # zero reference line here, and ticks use LogLocator instead of
            # a fixed linear spacing. LOG_AXIS_MIN sets an explicit positive
            # lower bound instead of letting the (0, ...) ranges autoscale.
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(left=LOG_AXIS_MIN, right=20)
            ax.set_ylim(bottom=LOG_AXIS_MIN, top=20)
            ax.xaxis.set_major_locator(LogLocator(base=10.0))
            ax.yaxis.set_major_locator(LogLocator(base=10.0))
            ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
            ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
        elif AXIS_SCALE == 'linear':
            ax.set_xscale('linear')
            ax.set_yscale('linear')
            plt.axhline(0, color='black', linestyle='--', linewidth=1)
            plt.axvline(0, color='black', linestyle='--', linewidth=1)
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.yaxis.set_major_locator(MultipleLocator(5))
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
        else:
            raise ValueError(f"AXIS_SCALE must be 'log' or 'linear', got {AXIS_SCALE!r}")

        ax.tick_params(axis='both', which='minor', length=4)
        ax.tick_params(axis='both', which='major', length=7)
        plt.xlabel('$c_{t\\gamma} [\\mathrm{TeV^{-1}}]$', fontsize=14)
        plt.ylabel('$c_{tg} [\\mathrm{TeV^{-1}}]$', fontsize=14)
        plt.title(f'95% CL Exclusion Contours in Coupling Space', fontsize=18)

        mass_legend = plt.legend(handles=mass_legend_handles, ncol=1, loc='upper right',
                                  fontsize=8, bbox_to_anchor=(1, 1), title='Mass')
        ax.add_artist(mass_legend)

        y_offset = 1 - 0.03 * (len(mass_legend_handles) + 2)
        for title, handles in legend_blocks:
            block_legend = plt.legend(handles=handles, ncol=1, loc='upper right', fontsize=8,
                                       bbox_to_anchor=(1, y_offset), title=title)
            ax.add_artist(block_legend)
            y_offset -= 0.03 * (len(handles) + 2)

        plt.grid()
        plt.tight_layout()
        plt.savefig(f"plots/limit_g_{var}.png")


if __name__ == "__main__":
    main()