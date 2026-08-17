from scipy.optimize import curve_fit
import numpy as np
from compute_uncertainties import compute_systematic_uncertainties as com_sys

gen_val = {
    "Signal_500": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [0.1817, 0.3563, 0.5890, 0.8798, 1.2289],
        "xsec_TT": [2.325581e+00, 2.355341e+00, 2.431030e+00, 2.516361e+00, 2.638057e+00],
    },
    "Signal_750": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [0.7612, 1.4921, 2.4665, 3.6846, 5.1463],
        "xsec_TT": [2.184231e-01, 2.259208e-01, 2.364962e-01, 2.537361e-01, 2.771376e-01],
    },
    "Signal_1000": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [1.939, 3.8009, 6.2831, 9.3859, 13.1092],
        "xsec_TT": [3.389237e-02, 3.563196e-02, 3.849664e-02, 4.252330e-02, 4.853439e-02],
    },
    "Signal_1250": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [3.913, 7.6707, 12.6802, 18.942, 26.456],
        "xsec_TT": [6.981026e-03, 7.439504e-03, 8.162695e-03, 9.325584e-03, 1.099415e-02],
    },
    "Signal_1500": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [6.883, 13.49100, 22.3014, 33.314, 46.530],
        "xsec_TT": [1.657904e-03, 1.778181e-03, 2.006442e-03, 2.359005e-03, 2.877716e-03],
    },
    "Signal_1750": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [11.0465, 21.6511, 35.7907, 53.4652, 74.6745],
        "xsec_TT": [4.340251e-04, 4.745353e-04, 5.391757e-04, 6.492909e-04, 8.134738e-04],
    },
    "Signal_2000": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [16.602, 32.5411, 53.792, 80.356, 112.233],
        "xsec_TT": [1.207999e-04, 1.310072e-04, 1.518117e-04, 1.855338e-04, 2.372002e-04],
    },
    "Signal_2250": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [23.7503, 46.5507, 76.951, 114.951, 160.552],
        "xsec_TT": [3.454268e-05, 3.757359e-05, 4.363370e-05, 5.329513e-05, 6.961889e-05],
    },
    "Signal_2500": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [32.6886, 64.0697, 105.911, 158.213, 220.975],
        "xsec_TT": [1.017082e-05, 1.108402e-05, 1.277048e-05, 1.576809e-05, 2.057444e-05],
    },
    "Signal_2750": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [43.6164, 85.4882, 141.317, 211.103, 294.847],
        "xsec_TT": [3.027667e-06, 3.277126e-06, 3.737304e-06, 4.596437e-06, 5.970259e-06],
    },
    "Signal_3000": {
        "g": np.array([5, 7, 9, 11, 13]),
        "width_tg": [56.732, 111.1959, 183.8136, 274.585, 383.512],
        "xsec_TT": [9.198671e-07, 9.839057e-07, 1.107198e-06, 1.332289e-06, 1.711181e-06],
    },
}


def br_prefactor(mass):
    g = gen_val[f"Signal_{mass}"]["g"]
    width = gen_val[f"Signal_{mass}"]["width_tg"]
    result = curve_fit(lambda x, a:  a * x ** 2, g, width)
    print("br factor:", result)

    return result[0][0]

def xsec_factors(mass):
    g = gen_val[f"Signal_{mass}"]["g"]
    xsec = gen_val[f"Signal_{mass}"]["xsec_TT"]
    result = curve_fit(lambda x, a, b, c:  a + b * x ** 2 + c * x ** 4, g, xsec)
    print("xsec factors:", result(m)[0])
    print("error xsec factors:","\n", result(m)[1],"\n")
    
    return result[0][0], result[0][1], result[0][2]

def compute_acceptance(mass, var):
    data = load("../output.coffea")
    hist_dict = data['hists']['total'][var]
    info = com_sys(hist_dict, mass)
    
    
    
    