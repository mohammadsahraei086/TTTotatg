def load_generation_uncertainties(mass):
    """
    Load generation uncertainties for a given mass point.
    These should be derived from your MC generation.
    
    Args:
        mass (int): Mass of the T particle
    
    Returns:
        dict: Dictionary with PDF, muR, muF uncertainties per bin
    """
    # These should be replaced with actual values from your MC studies
    # For example, from PDF variations (e.g., using NNPDF variations)
    # and scale variations (muR and muF variations by factor 2)
    
    # Example structure - you need to fill these with actual values
    uncertainties = {
        'pdf': {
            'values': [],  # Fractional uncertainty per bin
            'correlation': np.eye(6)  # Correlation matrix for PDF uncertainties
        },
        'muR': {
            'values': [],  # Fractional uncertainty per bin
            'correlation': np.eye(6)
        },
        'muF': {
            'values': [],  # Fractional uncertainty per bin
            'correlation': np.eye(6)
        }
    }
    
    # For now, using placeholder values based on typical top physics
    # You should replace these with actual numbers from your MC studies
    if mass in [500, 750, 1000]:
        uncertainties['pdf']['values'] = np.array([0.02, 0.022, 0.025, 0.028, 0.03, 0.032])
        uncertainties['muR']['values'] = np.array([0.015, 0.018, 0.02, 0.022, 0.025, 0.028])
        uncertainties['muF']['values'] = np.array([0.015, 0.018, 0.02, 0.022, 0.025, 0.028])
    elif mass in [1250, 1500, 1750]:
        uncertainties['pdf']['values'] = np.array([0.025, 0.028, 0.03, 0.033, 0.035, 0.038])
        uncertainties['muR']['values'] = np.array([0.02, 0.022, 0.025, 0.028, 0.03, 0.032])
        uncertainties['muF']['values'] = np.array([0.02, 0.022, 0.025, 0.028, 0.03, 0.032])
    else:  # Higher masses
        uncertainties['pdf']['values'] = np.array([0.03, 0.033, 0.035, 0.038, 0.04, 0.045])
        uncertainties['muR']['values'] = np.array([0.025, 0.028, 0.03, 0.033, 0.035, 0.04])
        uncertainties['muF']['values'] = np.array([0.025, 0.028, 0.03, 0.033, 0.035, 0.04])
    
    # If you have correlations between bins, you can fill the correlation matrices
    # For simplicity, we assume independent bins here (diagonal)
    for key in uncertainties:
        uncertainties[key]['correlation'] = np.eye(len(uncertainties[key]['values']))
    
    return uncertainties