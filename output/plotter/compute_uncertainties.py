import numpy as np
import hist
from typing import Dict
import json
from coffea.util import load


def _values_with_overflow_folded(hist_obj, variation, fold_underflow=False):
    """
    Return the 1D bin contents for a given `variation` slice of hist_obj,
    with the overflow bin folded into the last defined bin.

    hist's/boost-histogram's `.values(flow=True)` returns an array shaped
    like [underflow, bin_0, bin_1, ..., bin_N-1, overflow]. There's no
    single built-in call that folds *only* the overflow into the last bin
    while leaving every other bin untouched: the built-in UHI `sum` slicing
    action (e.g. h[a:hist.overflow:sum]) merges the whole selected range
    into one bin, so it's meant for collapsing several bins together, not
    for topping up one existing edge bin with flow content. So this does it
    by hand on the flow-inclusive array — which is exactly the same
    value/value arithmetic the built-in sum action performs internally.

    By default the underflow bin is simply dropped (not folded anywhere),
    matching "I don't want underflow, but I do want overflow added to the
    last bin". Pass fold_underflow=True if you ever want underflow folded
    into the first bin instead.
    """
    if hist_obj.name == 'diff_xsec_photon_pt':
        values = hist_obj[{'variation': variation}].values(flow=True)
        core = values[1:-1].copy()
        core[-1] += values[-1]
        if fold_underflow:
            core[0] += values[0]
        return core
    else:
        return hist_obj[{'variation': variation}].values()


def compute_systematic_uncertainties(hist_dict: Dict) -> Dict:
    """
    Compute MUR, MUF, and PDF uncertainties from histograms.
    
    Parameters:
    -----------
    hist_dict : dict
        Dictionary containing histograms with variations (from your coffea output)
    
    Returns:
    --------
    dict: Dictionary containing uncertainties for each sample
    """
    
    uncertainties = {}
    
    # Define the variations we need
    nominal_variation = 'nominal'
    mur_up = 'MUF=1.0_MUR=2.0_PDF=247000'
    mur_down = 'MUF=1.0_MUR=0.5_PDF=247000'
    muf_up = 'MUF=2.0_MUR=1.0_PDF=247000'
    muf_down = 'MUF=0.5_MUR=1.0_PDF=247000'
    
    # PDF variations (100 eigenvariations from 247001 to 247100)
    pdf_variations = [f'MUF=1.0_MUR=1.0_PDF=247{i:03d}' for i in range(1, 101)]
    
    for sample_name, hist_obj in hist_dict.items():
        if not isinstance(hist_obj, hist.Hist):
            continue
            
        print(f"Processing {sample_name}...")
        
        # Get nominal values (overflow folded into the last bin, underflow dropped)
        nominal = _values_with_overflow_folded(hist_obj, nominal_variation)
        
        # 1. MUR uncertainty
        mur_up_vals = _values_with_overflow_folded(hist_obj, mur_up)
        mur_down_vals = _values_with_overflow_folded(hist_obj, mur_down)
                
        mur_unc_up = abs(mur_up_vals - nominal)
        mur_unc_down = abs(nominal - mur_down_vals)
        
        # 2. MUF uncertainty
        muf_up_vals = _values_with_overflow_folded(hist_obj, muf_up)
        muf_down_vals = _values_with_overflow_folded(hist_obj, muf_down)
        
        muf_unc_up = abs(muf_up_vals - nominal)
        muf_unc_down = abs(nominal - muf_down_vals)
        
        # 3. PDF uncertainty (RMS of 100 variations)
        pdf_values = []
        for var in pdf_variations:
            try:
                values = _values_with_overflow_folded(hist_obj, var)
                pdf_values.append(values)
            except KeyError:
                continue
        
        if len(pdf_values) > 0:
            pdf_values = np.array(pdf_values)
            # RMS for each bin across all variations
            pdf_rms = np.sqrt(np.mean((pdf_values - nominal)**2, axis=0))
            
            # Symmetric PDF uncertainty (using RMS)
            pdf_unc = pdf_rms  # same for up and down
        else:
            pdf_unc = np.zeros_like(nominal)
        
        # 4. Total uncertainty (combine in quadrature)
        # For up: add the positive variations in quadrature
        # For down: add the negative variations in quadrature
        total_up = np.sqrt(mur_unc_up**2 + muf_unc_up**2 + pdf_unc**2)
        total_down = np.sqrt(mur_unc_down**2 + muf_unc_down**2 + pdf_unc**2)
        
        # Store the results for this sample
        uncertainties[sample_name] = {
            'nominal': nominal,
            # 'mur_up': mur_unc_up.tolist(),
            # 'mur_down': mur_unc_down.tolist(),
            # 'muf_up': muf_unc_up.tolist(),
            # 'muf_down': muf_unc_down.tolist(),
            # 'pdf_unc': pdf_unc.tolist(),  # symmetric
            'total_up': total_up,
            'total_down': total_down,
            # Also store bin information if needed
            'bin_edges': hist_obj.axes[0].edges.tolist(),
        }

    new_uncertainties = {}
    for mass in [500 + 250*i for i in range(11)]:
        new_uncertainties["Signal_" + str(mass)] = {}
        new_uncertainties["Signal_" + str(mass)]['nominal'] = (uncertainties["Signal_" + str(mass)]['nominal'] + uncertainties["ttaa_" + str(mass)]['nominal']).tolist()
        for val in ['total_up', 'total_down']:
            new_uncertainties["Signal_" + str(mass)][val] = np.sqrt(((uncertainties["Signal_" + str(mass)][val])**2)) #+ (uncertainties["ttaa_" + str(mass)][val])**2)).tolist()
        new_uncertainties["Signal_" + str(mass)]['bin_edges'] = uncertainties["Signal_" + str(mass)]['bin_edges']
    
    return new_uncertainties
    
def save_uncertainties(uncertainties: Dict, output_file: str = 'uncertainties.json'):
    """
    Save uncertainties to a JSON file.
    """
    with open(output_file, 'w') as f:
        json.dump(uncertainties, f, indent=2)
    print(f"\nUncertainties saved to {output_file}")

# Example usage
if __name__ == "__main__":
    # Option 1: If you already have the data loaded
    # hist_dict = data['hists']['total']
    # uncertainties = compute_systematic_uncertainties(hist_dict)
    
    # Option 2: Load from file
    # uncertainties = load_and_compute('output.coffea')
    
    # Option 3: If your data is in the format shown in your file
    # Assuming 'hists' is already loaded as shown in your print
    data = load("../output.coffea")
    for var in ["diff_xsec_photon_pt", "deltaphi_ll"]:
        hist_dict = data['hists']['total'][var]  # from your loaded data
    
        uncertainties = compute_systematic_uncertainties(hist_dict)
        
        # Save the results
        save_uncertainties(uncertainties, "json_files/" + var + "_error.json")
