import numpy as np
import hist
from typing import Dict
import json
from coffea.util import load

def compute_systematic_uncertainties(hist_dict: Dict, mass:int) -> Dict:
    uncertainties = {}
    
    # Define the variations we need
    nominal_variation = 'nominal'
    mur_up = 'MUF=1.0_MUR=2.0_PDF=247000'
    mur_down = 'MUF=1.0_MUR=0.5_PDF=247000'
    muf_up = 'MUF=2.0_MUR=1.0_PDF=247000'
    muf_down = 'MUF=0.5_MUR=1.0_PDF=247000'
    
    # PDF variations (100 eigenvariations from 247001 to 247100)
    pdf_variations = [f'MUF=1.0_MUR=1.0_PDF=247{i:03d}' for i in range(1, 101)]
    
    #######################################
    hist_ttag = hist_dict[f"Signal_{mass}"]
    #######################################
    nominal_ttag = hist_ttag[{'variation': nominal_variation}].values()
    # 1. MUR uncertainty
    mur_up_ttag = hist_ttag[{'variation': mur_up}].values()
    mur_down_ttag = hist_ttag[{'variation': mur_down}].values()
    
    mur_unc_ttag = (abs(mur_up_ttag - nominal_ttag) + abs(mur_down_ttag - nominal_ttag))/2
    
    # 2. MUF uncertainty
    muf_up_ttag = hist_ttag[{'variation': muf_up}].values()
    muf_down_ttag = hist_ttag[{'variation': muf_down}].values()
    
    muf_unc_ttag = (abs(muf_up_ttag - nominal_ttag) + abs(muf_down_ttag - nominal_ttag))/2
    
    # 3. PDF uncertainty (RMS of 100 variations)
    pdf_values = []
    for var in pdf_variations:
        values = hist_ttag[{'variation': var}].values()
        pdf_values.append(values)
        
    pdf_values = np.array(pdf_values)
    # RMS for each bin across all variations
    pdf_unc_ttag = np.sqrt(np.mean((pdf_values - nominal_ttag)**2, axis=0))

    #######################################
    hist_ttaa = hist_dict[f"ttaa_{mass}"]
    #######################################
    nominal_ttaa = hist_ttaa[{'variation': nominal_variation}].values()
    # 1. MUR uncertainty
    mur_up_ttaa = hist_ttaa[{'variation': mur_up}].values()
    mur_down_ttaa = hist_ttaa[{'variation': mur_down}].values()
    
    mur_unc_ttaa = (abs(mur_up_ttaa - nominal_ttaa) + abs(mur_down_ttaa - nominal_ttaa))/2
    
    # 2. MUF uncertainty
    muf_up_ttaa = hist_ttaa[{'variation': muf_up}].values()
    muf_down_ttaa = hist_ttaa[{'variation': muf_down}].values()
    
    muf_unc_ttaa = (abs(muf_up_ttaa - nominal_ttaa) + abs(muf_down_ttaa - nominal_ttaa))/2
    
    # 3. PDF uncertainty (RMS of 100 variations)
    pdf_values = []
    for var in pdf_variations:
        values = hist_ttaa[{'variation': var}].values()
        pdf_values.append(values)
        
    pdf_values = np.array(pdf_values)
    # RMS for each bin across all variations
    pdf_unc_ttaa = np.sqrt(np.mean((pdf_values - nominal_ttaa)**2, axis=0))
    
    nominal = (nominal_ttaa + nominal_ttag)
    
    mur = np.sqrt(mur_unc_ttag**2 + mur_unc_ttaa**2)
    muf = np.sqrt(muf_unc_ttag**2 + muf_unc_ttaa**2)
    pdf = np.sqrt(pdf_unc_ttag**2 + pdf_unc_ttaa**2)
    
    total_error = np.sqrt(mur_unc_ttag**2 + muf_unc_ttag**2 + pdf_unc_ttag**2 + mur_unc_ttaa**2 + muf_unc_ttaa**2 + pdf_unc_ttaa**2)
    total_rel_error = total_error/nominal
    
    uncertainties = {
        'mur_rel': (mur/nominal),
        'muf_rel': (muf/nominal),
        'pdf_rel': (pdf/nominal),  # symmetric
        # 'total_error': total_error.tolist(),
        # 'total_rel_error': total_rel_error.tolist(),
    }
    
    return nominal, uncertainties
    
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
        uncertainties_dict = {}
        for mass in [500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000]:
            hist_dict = data['hists']['total'][var]  # from your loaded data
    
            uncertainties = compute_systematic_uncertainties(hist_dict, mass)

            uncertainties_dict[f"Signal_{mass}"] = uncertainties
            
        
        # Save the results
        # save_uncertainties(uncertainties_dict, "json_files/" + var + "_error.json")