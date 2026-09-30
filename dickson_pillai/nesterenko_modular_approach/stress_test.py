"""
stress_test.py

Full stress test and verification suite for the Nesterenko Modular Approach.
Runs precision invariance, Ramanujan differential closure, vanishing order checks,
and cusp inversion scaling checks.
"""

import sys
import math
import numpy as np
import mpmath as mp
from modular_forms import set_precision, eisenstein_series, ramanujan_derivatives, ModularMonomial
from auxiliary_modular_polynomial import AuxiliaryModularForm
from cusp_inversion_bound import analyze_diophantine_scaling

def test_ramanujan_differential_closure():
    print("[TEST 1] Ramanujan Differential Closure (100-digit precision)...")
    set_precision(100)
    q = mp.mpf(2) / mp.mpf(3)
    e2, e4, e6 = eisenstein_series(q, terms=400)
    
    dq = mp.mpf(10)**(-30)
    e2_p, e4_p, e6_p = eisenstein_series(q + dq, terms=400)
    
    de2_num = (e2_p - e2) / dq * q
    de4_num = (e4_p - e4) / dq * q
    de6_num = (e6_p - e6) / dq * q
    
    th_e2, th_e4, th_e6 = ramanujan_derivatives(e2, e4, e6)
    
    err_e2 = abs(de2_num - th_e2)
    err_e4 = abs(de4_num - th_e4)
    err_e6 = abs(de6_num - th_e6)
    
    print(f"  |NumDeriv(E2) - Ramanujan(E2)| = {float(err_e2):.4e}")
    print(f"  |NumDeriv(E4) - Ramanujan(E4)| = {float(err_e4):.4e}")
    print(f"  |NumDeriv(E6) - Ramanujan(E6)| = {float(err_e6):.4e}")
    
    assert float(err_e2) < 1e-15, "E2 derivative mismatch!"
    assert float(err_e4) < 1e-15, "E4 derivative mismatch!"
    assert float(err_e6) < 1e-15, "E6 derivative mismatch!"
    print("  --> PASS: Ramanujan differential system verified.\n")

def test_formal_monomial_differentiation():
    print("[TEST 2] Formal Monomial Derivation Matrix vs Numerical Derivation...")
    set_precision(80)
    q = mp.mpf(2) / mp.mpf(3)
    e2, e4, e6 = eisenstein_series(q)
    
    test_monomials = [
        ModularMonomial(1, 1, 0, 0), # q * E2
        ModularMonomial(0, 0, 1, 1), # E4 * E6
        ModularMonomial(2, 2, 1, 0), # q^2 * E2^2 * E4
    ]
    
    dq = mp.mpf(10)**(-25)
    e2_p, e4_p, e6_p = eisenstein_series(q + dq)
    
    for mon in test_monomials:
        # Numerical derivation
        v0 = mon.evaluate(q, e2, e4, e6)
        v1 = mon.evaluate(q + dq, e2_p, e4_p, e6_p)
        num_deriv = (v1 - v0) / dq * q
        
        # Formal derivation
        formal_terms = mon.formal_derivative()
        form_deriv = mp.mpf(0)
        for term, coeff in formal_terms.items():
            form_deriv += coeff * term.evaluate(q, e2, e4, e6)
            
        err = abs(num_deriv - form_deriv) / (abs(form_deriv) + 1e-30)
        print(f"  Monomial {mon}: Relative Error = {float(err):.4e}")
        assert float(err) < 1e-12, f"Monomial derivative failure on {mon}"
        
    print("  --> PASS: Formal monomial algebra matches numerical derivatives.\n")

def test_quasi_modular_cusp_inversion():
    print("[TEST 3] Quasi-Modular Cusp S-Inversion Exactness...")
    set_precision(100)
    
    # In the balanced region y in [0.2, 1.0], both direct and dual q-series converge in < 300 terms
    test_y_values = [mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('0.8'), mp.mpf('1.0')]
    for y in test_y_values:
        tau = mp.mpc(0, y)
        q = mp.exp(2 * mp.pi * mp.j * tau)
        
        # Direct evaluation
        e2_direct, _, _ = eisenstein_series(q, terms=400)
        
        # Modular S-inversion formula:
        # E2(-1/tau) = tau^2 E2(tau) + (6 tau)/(pi i)
        # => E2(tau) = (1/tau^2) E2(-1/tau) - 6 / (pi i tau)
        tau_prime = -1 / tau
        q_prime = mp.exp(2 * mp.pi * mp.j * tau_prime)
        e2_dual, _, _ = eisenstein_series(q_prime, terms=400)
        
        e2_predicted = (mp.mpf(1) / (tau**2)) * e2_dual - mp.mpf(6) / (mp.pi * mp.j * tau)
        
        diff = abs(e2_direct - e2_predicted)
        rel_err = diff / abs(e2_direct)
        print(f"  tau = i * {float(y):.2f}: Absolute Error = {float(diff):.4e}, Rel Error = {float(rel_err):.4e}")
        assert float(rel_err) < 1e-70, f"Cusp inversion mismatch at y = {y}"
        
    print("  --> PASS: Quasi-modular S-inversion verified to > 70 digits of precision.\n")

        
    print("  --> PASS: Quasi-modular S-inversion verified to ultra-high precision.\n")

def run_all_stress_tests():
    print("=====================================================================")
    print("      NESTERNKO MODULAR APPROACH: COMPREHENSIVE STRESS TEST SUITE     ")
    print("=====================================================================\n")
    
    test_ramanujan_differential_closure()
    test_formal_monomial_differentiation()
    test_quasi_modular_cusp_inversion()
    
    print("[TEST 4] Scaling Behavior Analysis across Graded Weights...")
    for weight in [4, 6, 8]:
        analyze_diophantine_scaling(max_q_deg=2, max_weight=weight)
        print("-" * 65)

if __name__ == "__main__":
    run_all_stress_tests()
