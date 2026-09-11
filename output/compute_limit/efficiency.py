from coffea.util import load
import matplotlib.pyplot as plt
import contourpy
import numpy as np
from scipy.stats import chi2
from plot_limits import get_contour

data = load("../output.coffea")

plt.figure(figsize=(10, 8))
ax = plt.gca()
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
# plt.rcParams['text.usetex'] = True
colors = ['#D55E00', '#E69F00', '#009E73', '#B8DC70' ,  '#0072B2',  '#56B4E9',
        '#BDB76B', '#F0E442', '#9400D3' , '#CC79A7' , '#708090', '#DA70D6', '#800000']

for i, m in enumerate([500, 750, 1000, 1250]):
    X, Y, Z = get_contour(
                mass=m,
                var="diff_xsec_photon_pt",
                g3g_range=(0, 1100),
                g3gamma_range=(0, 1100),
                n_points=200,
                hl_lhc=False,
                untruncated=False
            )
    Z_min = np.min(Z)
    chi2_95 = chi2.ppf(0.95, df=1)
    
    gen = contourpy.contour_generator(x=X, y=Y, z=Z - Z_min)
    lines = gen.lines(chi2_95)
    
    factor = 0.1/5000
    lam = 1/np.sqrt((factor * lines[0][:,0])**2 + (factor*lines[0][:,1])**2)
    
    mtt = data["arrays"][f"Signal_{m}"]["MTT_array"].value
    eff = []
    for l in lam:
       eff.append(len(mtt[mtt<l])/len(mtt))
    
    plt.plot(lam, eff, label=fr'$m_T = {m:.0f}\ \mathrm{{GeV}}$', color=colors[i],
                            linewidth=2, linestyle='solid')

ax.text(
    0.052, 1.03,
    #r'$\mathbf{95\%\ CL\ exclusion}$' +
    #r'$, 140\ \mathrm{fb}^{-1} \ $' +
    r'$\mathbf{t\bar{t}\gamma}$' + r' , $\mathbf{EFT \ Survival \ Efficiency}$',
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
ax.text(
    0.82, 1.03,
    r'$\mathbf{\mathcal{L} = 138\ fb^{-1}}$',
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
plt.legend(ncol=1, loc='lower right', fontsize=15,
         bbox_to_anchor=(0.97, 0.04),
         frameon=True, fancybox=True,
         framealpha=0.7, edgecolor="gray", borderpad=0.6)

ax.minorticks_on()
ax.tick_params(axis='both', which='minor', length=3)
ax.tick_params(axis='both', which='major', length=7)
ax.tick_params(axis='both', which='both' , top=True, right=True, direction='in', labelsize=15)
plt.xlim(0, 6e3)
ax.set_ylabel(r"Efficiency"    , fontsize=15)
ax.set_xlabel(r"$\Lambda_{eff.} [\mathrm{GeV}]$", fontsize=15)
ax.grid(True, which='major', linestyle='-', linewidth=0.7, alpha=0.3)
plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.savefig(f"plots/efficiency.pdf")