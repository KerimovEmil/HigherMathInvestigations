"""
modular_forms.py

High-precision computation of Eisenstein series E2, E4, E6, the Ramanujan derivation
theta = q d/dq, and graded monomial evaluations using mpmath.
"""

import mpmath as mp

def set_precision(dps: int = 100):
    """Set global precision in decimal places."""
    mp.mp.dps = dps

def eisenstein_series(q, terms: int = 300):
    """
    Compute Eisenstein series E2(q), E4(q), E6(q) to high precision.
    Supports both real and complex q inside the unit disk (|q| < 1).
    """
    q_val = q
    is_complex = isinstance(q, mp.mpc) or (hasattr(q, 'imag') and q.imag != 0)
    
    zero = mp.mpc(0) if is_complex else mp.mpf(0)
    one = mp.mpc(1) if is_complex else mp.mpf(1)
    
    e2_sum = zero
    e4_sum = zero
    e6_sum = zero
    
    for n in range(1, terms + 1):
        qn = q_val ** n
        denom = one - qn
        term = qn / denom
        
        e2_sum += n * term
        e4_sum += (n ** 3) * term
        e6_sum += (n ** 5) * term
        
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 5):
            break
            
    e2 = one - mp.mpf(24) * e2_sum
    e4 = one + mp.mpf(240) * e4_sum
    e6 = one - mp.mpf(504) * e6_sum
    
    return e2, e4, e6

def ramanujan_derivatives(e2, e4, e6):
    """
    Compute the Ramanujan derivatives theta(E2), theta(E4), theta(E6)
    where theta = q d/dq.
    theta E2 = (E2^2 - E4) / 12
    theta E4 = (E2 E4 - E6) / 3
    theta E6 = (E2 E6 - E4^2) / 2
    """
    theta_e2 = (e2**2 - e4) / mp.mpf(12)
    theta_e4 = (e2 * e4 - e6) / mp.mpf(3)
    theta_e6 = (e2 * e6 - e4**2) / mp.mpf(2)
    return theta_e2, theta_e4, theta_e6

class ModularMonomial:
    """
    Represents a monomial in Q[q, E2, E4, E6]:
    m = q^{j0} * E2^{j2} * E4^{j4} * E6^{j6}
    """
    def __init__(self, j0: int, j2: int, j4: int, j6: int):
        self.j0 = j0
        self.j2 = j2
        self.j4 = j4
        self.j6 = j6
        self.key = (j0, j2, j4, j6)
        self.weight = 2 * j2 + 4 * j4 + 6 * j6

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return isinstance(other, ModularMonomial) and self.key == other.key

    def evaluate(self, q, e2, e4, e6):
        """Evaluate the monomial at given high-precision values."""
        val = (mp.mpf(q) ** self.j0) * (e2 ** self.j2) * (e4 ** self.j4) * (e6 ** self.j6)
        return val

    def formal_derivative(self):
        """
        Compute theta(m) formally as a linear combination of monomials in Q[q, E2, E4, E6].
        Returns a dict: {ModularMonomial: rational_coefficient}
        """
        terms = {}
        
        def add_term(mon, coeff):
            if mon in terms:
                terms[mon] += coeff
            else:
                terms[mon] = coeff

        # 1. Derivative wrt q: theta(q^{j0}) = j0 * q^{j0}
        if self.j0 != 0:
            add_term(ModularMonomial(self.j0, self.j2, self.j4, self.j6), mp.mpf(self.j0))

        # 2. Derivative wrt E2: j2 * E2^{j2-1} * (E2^2 - E4)/12
        if self.j2 != 0:
            add_term(ModularMonomial(self.j0, self.j2 + 1, self.j4, self.j6), mp.mpf(self.j2) / mp.mpf(12))
            add_term(ModularMonomial(self.j0, self.j2 - 1, self.j4 + 1, self.j6), -mp.mpf(self.j2) / mp.mpf(12))

        # 3. Derivative wrt E4: j4 * E4^{j4-1} * (E2 E4 - E6)/3
        if self.j4 != 0:
            add_term(ModularMonomial(self.j0, self.j2 + 1, self.j4, self.j6), mp.mpf(self.j4) / mp.mpf(3))
            add_term(ModularMonomial(self.j0, self.j2, self.j4 - 1, self.j6 + 1), -mp.mpf(self.j4) / mp.mpf(3))

        # 4. Derivative wrt E6: j6 * E6^{j6-1} * (E2 E6 - E4^2)/2
        if self.j6 != 0:
            add_term(ModularMonomial(self.j0, self.j2 + 1, self.j4, self.j6), mp.mpf(self.j6) / mp.mpf(2))
            add_term(ModularMonomial(self.j0, self.j2, self.j4 + 2, self.j6 - 1), -mp.mpf(self.j6) / mp.mpf(2))

        return terms

    def __repr__(self):
        return f"q^{self.j0} E2^{self.j2} E4^{self.j4} E6^{self.j6}"
