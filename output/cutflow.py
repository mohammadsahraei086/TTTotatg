from coffea.util import load
from collections import defaultdict
import numpy as np

# Load the data from the coffea file
output = load("output.coffea")

# Extract the cutflow information
data = output["cutflow"]["total"]

# Cut names mapping for LaTeX
cut_mapping = {
    'n_lep=2': r'$N_l=2$',
    'leadingLepPT': r'$l_{p_T}^{\text{leading}}>25$',
    'OCLep': r'$\text{OC lepton}$',
    'lepInvariantMass': r'$M_{ll}>20$',
    'onePhoton': r'$N_\gamma=1$',
    'atLeastOneBJet': r'$N_{\text{bjet}}>=1$'
}

# List of cuts in order
cuts = ['n_lep=2', 'leadingLepPT', 'OCLep', 'lepInvariantMass', 'onePhoton', 'atLeastOneBJet']


def get_mass_from_sample(sample_name):
    """Extract mass value from sample name"""
    return int(sample_name.split('_')[-1])

def get_sample_group(sample_name):
    """Get sample group (ttaa or Signal)"""
    if sample_name.startswith('ttaa'):
        return 'ttaa'
    elif sample_name.startswith('Signal'):
        return 'Signal'
    return 'unknown'

def get_sample_label(group):
    """Return the LaTeX label for a sample group"""
    if group == 'ttaa':
        return r'tt$\gamma\gamma$'
    elif group == 'Signal':
        return r'tt$\gamma$'
    return group

def generate_table1(data, sample_group, scale_factor=0.70):
    """
    Table 1: Selection efficiency (%) for each cut compared to primary
    Using yield values for proper normalization
    Transposed: Masses as rows, cuts as columns
    """
    # Get samples for this group
    samples = {name: info for name, info in data.items() if get_sample_group(name) == sample_group}
    
    # Sort by mass
    sorted_samples = sorted(samples.items(), key=lambda x: get_mass_from_sample(x[0]))
    
    # Get the label for this sample group
    sample_label = get_sample_label(sample_group)
    
    # Build table header with cuts as columns - using |c| for lines between columns
    header = "\\begin{table}[h]\n\\centering\n"
    header += f"\\scalebox{{{scale_factor}}}{{\n"
    header += "\\begin{tabular}{|c|" + "|".join(["c"] * len(cuts)) + "|}\n\\hline\n"
    header += "Mass (GeV) & " + " & ".join([cut_mapping[cut] for cut in cuts]) + " \\\\ \\hline\n"
    
    # Build rows for each mass
    rows = []
    for sample_name, info in sorted_samples:
        mass = get_mass_from_sample(sample_name)
        row = f"{mass}"
        primary_yield = info['yield']['primary']
        for cut in cuts:
            cut_yield = info['yield'][cut]
            efficiency = (cut_yield / primary_yield) * 100 if primary_yield > 0 else 0
            row += f" & {efficiency:.2f}\\%"
        row += " \\\\ \\hline\n"
        rows.append(row)
    
    footer = "\\end{tabular}\n}\n\\caption{Selection efficiency for " + sample_label + " samples}\n\\end{table}"
    
    return header + "".join(rows) + footer

def generate_table3(data, scale_factor=0.70):
    """
    Table 3: Two columns showing:
    1. Ratio of ttaa to (ttaa + Signal) after all selections for each mass
    2. Ratio of ttaa to Signal (ttgamma/ttgamma) for same mass
    Using yield values for proper normalization
    Masses as rows, two columns
    """
    # Group samples by mass
    mass_groups = defaultdict(dict)
    
    for sample_name, info in data.items():
        mass = get_mass_from_sample(sample_name)
        group = get_sample_group(sample_name)
        if group in ['ttaa', 'Signal']:
            mass_groups[mass][group] = info['yield']['atLeastOneBJet']
    
    # Sort by mass
    sorted_masses = sorted(mass_groups.keys())
    
    # Build table header - using |c|c|c| for lines between columns
    header = "\\begin{table}[h]\n\\centering\n"
    header += f"\\scalebox{{{scale_factor}}}{{\n"
    header += "\\begin{tabular}{|c|c|c|}\n\\hline\n"
    header += "Mass (GeV) & $\\frac{\\text{tt}\\gamma\\gamma}{\\text{tt}\\gamma\\gamma + \\text{tt}\\gamma}$ (\\%) & $\\frac{\\text{tt}\\gamma\\gamma}{\\text{tt}\\gamma}$ (\\%) \\\\ \\hline\n"
    
    # Build rows
    rows = []
    for mass in sorted_masses:
        if 'ttaa' in mass_groups[mass] and 'Signal' in mass_groups[mass]:
            ttaa = mass_groups[mass]['ttaa']
            signal = mass_groups[mass]['Signal']  # This is ttgamma
            
            # Calculate ratio 1: ttaa / (ttaa + signal)
            ratio1 = (ttaa / (ttaa + signal)) * 100 if (ttaa + signal) > 0 else 0
            
            # Calculate ratio 2: ttaa / signal (now as percentage)
            ratio2 = (ttaa / signal) * 100 if signal > 0 else 0
            
            rows.append(f"{mass} & {ratio1:.2f}\\% & {ratio2:.2f}\\% \\\\ \\hline\n")
    
    footer = "\\end{tabular}\n}\n\\caption{Ratios of tt$\\gamma\\gamma$ events to tt$\\gamma$ events after all selections}\n\\end{table}"
    
    return header + "".join(rows) + footer

def generate_complete_tex_file(data, filename="complete_tables.tex"):
    """
    Generate a complete LaTeX file with all tables in frames
    """
    tex_content = """\\documentclass{beamer}
\\usetheme{default}
\\usepackage{amsmath}
\\usepackage{bm}
\\usepackage{booktabs}
\\usepackage{float}
\\usepackage[caption = false]{subfig}
\\usepackage[]{graphicx}
\\usepackage{adjustbox}

\\title{\\textbf{Analysis Tables}}
\\subtitle{Single Lepton Channel}

\\begin{document}
\\begin{frame}[plain]
    \\maketitle
\\end{frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\\begin{frame}{tt$\\gamma\\gamma$ Samples Selection Efficiency}
"""
    
    # Add Table 1a
    tex_content += generate_table1(data, 'ttaa', scale_factor=0.70)
    
    tex_content += """
\\end{frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\\begin{frame}{tt$\\gamma$ Samples Selection Efficiency}
"""
    
    # Add Table 1b (Signal samples, now labeled as ttgamma)
    tex_content += generate_table1(data, 'Signal', scale_factor=0.7)
    
    tex_content += """
\\end{frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\\begin{frame}{tt$\\gamma\\gamma$ / tt$\\gamma$ Ratios}
"""
    
    # Add Table 3
    tex_content += generate_table3(data, scale_factor=0.75)
    
    tex_content += """
\\end{frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\\end{document}
"""
    
    return tex_content

# Generate and save the complete LaTeX file
tex_content = generate_complete_tex_file(data)
with open("complete_tables.tex", "w") as f:
    f.write(tex_content)

print("Complete LaTeX file has been written to complete_tables.tex")

# Also print individual tables for reference
print("\n" + "=" * 80)
print("TABLE 1a: tt$\\gamma\\gamma$ samples - Selection efficiency")
print("=" * 80)
print(generate_table1(data, 'ttaa', scale_factor=0.75))

print("\n" + "=" * 80)
print("TABLE 1b: tt$\\gamma$ samples - Selection efficiency")
print("=" * 80)
print(generate_table1(data, 'Signal', scale_factor=0.75))

print("\n" + "=" * 80)
print("TABLE 3: tt$\\gamma\\gamma$ / tt$\\gamma$ ratios")
print("=" * 80)
print(generate_table3(data, scale_factor=0.75))