"""
auxiliary_modular_polynomial.py

Constructs the auxiliary modular polynomial Phi(q) in Q[q, E2, E4, E6]
and determines integer coefficients c_j to maximize the vanishing order T at q = 2/3.
"""

import math
import numpy as np
import mpmath as mp
from modular_forms import set_precision, eisenstein_series, ramanujan_derivatives, ModularMonomial

class AuxiliaryModularForm:
    def __init__(self, max_q_deg: int = 4, max_weight: int = 12):
        self.max_q_deg = max_q_deg
        self.max_weight = max_weight
        self.basis = []
        self._build_basis()

    def _build_basis(self):
        """Construct the monomial basis {q^j0 E2^j2 E4^j4 E6^j6} with 2j2 + 4j4 + 6j6 <= max_weight."""
        self.basis = []
        for j0 in range(self.max_q_deg + 1):
            for j2 in range(self.max_weight // 2 + 1):
                for j4 in range(self.max_weight // 4 + 1):
                    for j6 in range(self.max_weight // 6 + 1):
                        if 2 * j2 + 4 * j4 + 6 * j6 <= self.max_weight:
                            self.basis.append(ModularMonomial(j0, j2, j4, j6))

    def evaluate_derivative_matrix(self, q, order_T: int):
        """
        Build the evaluation matrix M of size (order_T, len(basis)),
        where M[m, i] = (theta^m basis[i])(q).
        """
        e2, e4, e6 = eisenstein_series(q)
        
        current_polys = [{mon: mp.mpf(1)} for mon in self.basis]
        matrix = []
        
        for m in range(order_T):
            row = []
            for poly in current_polys:
                val = mp.mpf(0)
                for mon, coeff in poly.items():
                    val += coeff * mon.evaluate(q, e2, e4, e6)
                row.append(val)
            matrix.append(row)
            
            # Update to next derivative theta^{m+1}
            next_polys = []
            for poly in current_polys:
                next_poly = {}
                for mon, coeff in poly.items():
                    deriv_terms = mon.formal_derivative()
                    for dmon, dcoeff in deriv_terms.items():
                        if dmon in next_poly:
                            next_poly[dmon] += coeff * dcoeff
                        else:
                            next_poly[dmon] = coeff * dcoeff
                next_polys.append(next_poly)
            current_polys = next_polys

        return matrix

    def find_auxiliary_form(self, q=mp.mpf(2)/mp.mpf(3), order_T: int = None):
        """
        Find integer coefficients c in Z^{dim} such that theta^m Phi(q) = 0 for m < T.
        Uses SVD / nullspace analysis.
        """
        dim = len(self.basis)
        if order_T is None:
            order_T = dim - 1
            
        mat = self.evaluate_derivative_matrix(q, order_T)
        
        # Convert matrix to float numpy array for SVD nullspace
        np_mat = np.zeros((order_T, dim), dtype=np.float64)
        for r in range(order_T):
            for c in range(dim):
                np_mat[r, c] = float(mat[r][c])
                
        # SVD
        U, S, Vt = np.linalg.svd(np_mat)
        
        # The right-singular vector corresponding to smallest singular value
        null_vec = Vt[-1, :]
        
        # Scale to integer approximation
        max_val = np.max(np.abs(null_vec))
        if max_val > 0:
            null_vec = null_vec / max_val
            
        return self.basis, null_vec, S, mat

if __name__ == "__main__":
    set_precision(50)
    print("Testing Auxiliary Modular Form construction...")
    aux = AuxiliaryModularForm(max_q_deg=2, max_weight=6)
    print(f"Basis size: {len(aux.basis)} monomials")
    for i, m in enumerate(aux.basis):
        print(f"  [{i:2d}] {m} (weight={m.weight})")
        
    basis, null_vec, S, mat = aux.find_auxiliary_form(order_T=len(aux.basis)-1)
    print(f"\nSingular values S (smallest = {S[-1]:.4e})")
    print(f"Approximate null vector norm: {np.linalg.norm(null_vec):.4f}")
