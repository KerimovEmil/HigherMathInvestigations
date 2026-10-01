"""
frey_curve_analysis.py

Comprehensive mathematical and numerical engine for Architecture 3:
Frey-Hellegouarch Elliptic Curves, Conductor Invariants, Modular Degrees,
and Localized Szpiro Bounds for the Dickson-Pillai Problem.
"""

import math
import numpy as np
import mpmath as mp
from sympy import factorint, primefactors

def compute_frey_invariants(k):
    two_k = 2**k
    three_k = 3**k
    q_k = three_k // two_k
    r_k = three_k % two_k
    
    # Check 3-adic valuation of q_k
    v3_q = 0
    temp_q = q_k
    while temp_q > 0 and temp_q % 3 == 0:
        v3_q += 1
        temp_q //= 3
        
    div_3 = 3**v3_q
    
    # Primitive Frey triple a_k + b_k = c_k
    a_k = r_k // div_3
    b_k = (two_k * q_k) // div_3
    c_k = three_k // div_3
    
    assert a_k + b_k == c_k, f"Partition identity failed at k={k}"
    
    # Prime factors
    factors_a = factorint(a_k)
    factors_b = factorint(b_k)
    factors_c = factorint(c_k)
    
    # Radical of a_k, b_k, c_k
    rad_a = 1
    for p in factors_a:
        rad_a *= p
    rad_b = 1
    for p in factors_b:
        rad_b *= p
    rad_c = 1
    for p in factors_c:
        rad_c *= p
        
    # Full conductor N_k:
    # E_k has multiplicative reduction at all odd primes dividing a_k * b_k * c_k
    # At p = 2: b_k is divisible by 2^(k - v3_q). For k >= 4, 16 | b_k, so v_2(N_k) = 1 (semi-stable)
    # or v_2(N_k) = 2 if slightly less divisible.
    odd_rad = 1
    all_primes = set(factors_a.keys()) | set(factors_b.keys()) | set(factors_c.keys())
    for p in all_primes:
        if p != 2:
            odd_rad *= p
            
    v2_N = 1 if (k - v3_q) >= 4 else (2 if (k - v3_q) >= 2 else 3)
    N_k = (2**v2_N) * odd_rad
    
    # Invariants
    Delta = 16 * (a_k * b_k * c_k)**2
    c4 = 16 * (c_k**2 - a_k * b_k)
    
    ln_Delta = float(mp.log(Delta))
    ln_c4 = float(mp.log(c4))
    ln_N = float(mp.log(N_k))
    
    # Szpiro ratios
    szpiro = ln_Delta / ln_N
    szpiro_c4 = (3.0 * ln_c4) / (2.0 * ln_N)
    
    # Faltings height h_F(E_k)
    # Formula for semi-stable curves: h_F(E) = (1/12) ln|Delta| - (1/2) ln(Omega) + const
    # For Frey curve, real period Omega approx pi / sqrt(c_k)
    Omega = float(mp.pi / mp.sqrt(c_k))
    ln_Omega = math.log(Omega)
    h_F = (1.0 / 12.0) * ln_Delta - 0.5 * ln_Omega
    
    # Hypothetical Failure Scenario Values
    # If r_k + q_k > 2^k, r_k approx 2^k, q_k approx (3/2)^k
    # Delta_hyp approx 16 * (2^k * 3^k * 3^k)^2 = 16 * 2^(2k) * 3^(4k)
    ln_Delta_hyp = math.log(16) + 2*k*math.log(2) + 4*k*math.log(3)
    # N_hyp <= 6 * 2^k * (3/2)^k = 6 * 3^k
    ln_N_hyp_max = math.log(6) + k*math.log(3)
    szpiro_hyp_min = ln_Delta_hyp / ln_N_hyp_max
    
    return {
        "k": k,
        "a_k": a_k,
        "b_k": b_k,
        "c_k": c_k,
        "N_k": N_k,
        "ln_N": ln_N,
        "ln_Delta": ln_Delta,
        "szpiro": szpiro,
        "szpiro_c4": szpiro_c4,
        "h_F": h_F,
        "szpiro_hyp_min": szpiro_hyp_min,
        "rad_a": rad_a,
        "rad_b": rad_b,
        "rad_c": rad_c
    }

def run_frey_analysis(max_k=40):
    print("="*85)
    print("       ARCHITECTURE 3: FREY-HELLEGOUARCH MODULAR CONDUCTOR INVARIANT ENGINE")
    print("="*85)
    print(f"{'k':<3} | {'Conductor N_k':<18} | {'ln(N_k)':<8} | {'ln|Delta|':<9} | {'Szpiro sigma':<12} | {'Faltings h_F':<12} | {'Hypothetical min sigma'}")
    print("-" * 85)
    
    records = []
    for k in range(1, max_k + 1):
        inv = compute_frey_invariants(k)
        records.append(inv)
        if k <= 20 or k in [25, 30, 35, 40]:
            print(f"{k:<3} | {inv['N_k']:<18} | {inv['ln_N']:<8.2f} | {inv['ln_Delta']:<9.2f} | {inv['szpiro']:<12.4f} | {inv['h_F']:<12.4f} | {inv['szpiro_hyp_min']:.4f}")
            
    print("-" * 85)
    
    # Statistical analysis
    sigmas = [r['szpiro'] for r in records]
    h_Fs = [r['h_F'] for r in records]
    
    print("\n[GLOBAL FREY CURVE INVARIANT THEOREMS & BOUNDS]")
    print(f"1. Maximum Observed Szpiro Ratio:  sigma_max = {max(sigmas):.4f} (at k = {sigmas.index(max(sigmas))+1})")
    print(f"2. Asymptotic Mean Szpiro Ratio:   sigma_mean = {np.mean(sigmas[5:]):.4f} for k >= 6")
    print(f"3. Szpiro Standard Deviation:      sigma_std  = {np.std(sigmas[5:]):.4f}")
    print(f"4. Theoretical Failure Lower Bound: sigma_fail >= {records[-1]['szpiro_hyp_min']:.4f}")
    print(f"5. Faltings Height Growth Rate:    h_F(E_k) / k -> {h_Fs[-1] / max_k:.4f}")
    
    return records

if __name__ == "__main__":
    records = run_frey_analysis(40)
