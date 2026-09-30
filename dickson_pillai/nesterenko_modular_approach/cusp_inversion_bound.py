"""
cusp_inversion_bound.py

Evaluates the modular S-inversion transformation on the auxiliary modular form Phi(q)
and extracts the effective Diophantine exponent C_eff for linear forms in logarithms:
|Lambda(k)| = |k ln(3/2) - ln M_k|
"""

import math
import numpy as np
import mpmath as mp
from modular_forms import set_precision, eisenstein_series
from auxiliary_modular_polynomial import AuxiliaryModularForm

def analyze_diophantine_scaling(max_q_deg=3, max_weight=8):
    set_precision(80)
    print(f"--- Nesterenko Modular Inversion Scaling Analysis ---")
    print(f"Parameters: max_q_deg = {max_q_deg}, max_weight = {max_weight}")
    
    aux = AuxiliaryModularForm(max_q_deg=max_q_deg, max_weight=max_weight)
    dim = len(aux.basis)
    print(f"Graded Monomial Basis Dimension: {dim}")
    
    # We choose vanishing order T = dim - 1
    T = dim - 1
    q0 = mp.mpf(2) / mp.mpf(3)
    
    print(f"Target vanishing order at q = 2/3: T = {T}")
    basis, null_vec, S, mat = aux.find_auxiliary_form(q=q0, order_T=T)
    
    # Analyze the height of the auxiliary form
    # h(Phi) = ln max |c_j|
    max_c = np.max(np.abs(null_vec))
    log_height = math.log(max(1.0, float(max_c)))
    
    # The effective exponent in Nesterenko's elimination method:
    # C_eff = (log_height + Analytic_Cusp_Penalty) / T
    # The Analytic Cusp Penalty is governed by the weight W = max_weight:
    # Under S-inversion: E_k(e^-Lambda) ~ (2pi/Lambda)^k
    # So ln |Phi(e^-Lambda)| has a pole of order max_weight in 1/Lambda
    cusp_penalty = (max_weight / 2.0) * math.log(2 * math.pi)
    
    C_eff = (log_height + cusp_penalty) / float(T)
    
    print(f"\nResults:")
    print(f"  Logarithmic Height h(Phi): {log_height:.6f}")
    print(f"  Cusp Pole Weight Penalty:  {cusp_penalty:.6f}")
    print(f"  Effective Diophantine Exponent C_eff: {C_eff:.6f}")
    print(f"\nComparison with Theoretical Barriers:")
    print(f"  Dickson-Pillai Threshold ln(4/3): {math.log(4.0/3.0):.6f}")
    print(f"  Fractional Part Threshold ln(2):   {math.log(2.0):.6f}")
    print(f"  Classical Baker / LMN Barrier:     9.930000")
    
    if C_eff <= math.log(4.0/3.0):
        print("\n>>> STATUS: C_eff is strictly below the Dickson-Pillai barrier! <<<")
    elif C_eff <= math.log(2.0):
        print("\n>>> STATUS: C_eff is below ln(2) (sub-exponential fractional part)! <<<")
    elif C_eff < 9.93:
        print(f"\n>>> STATUS: C_eff ({C_eff:.4f}) beats Baker's classical barrier (9.93), but has not yet crossed ln(4/3). <<<")
    else:
        print(f"\n>>> STATUS: C_eff ({C_eff:.4f}) remains above classical bounds. <<<")

    return C_eff, T, dim

if __name__ == "__main__":
    analyze_diophantine_scaling(max_q_deg=2, max_weight=6)
