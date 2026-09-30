"""
The Dickson–Pillai Condition (Waring's Problem Sub-Problem):
Condition: 2^k * {(3/2)^k} + floor((3/2)^k) <= 2^k for all k >= 1.

Historical Context & Link to Waring's Problem:
    In 1770, Edward Waring posed the problem of whether every integer is expressible
    as a sum of at most g(k) k-th powers.
    In 1936, Leonard Eugene Dickson and Subbayya Sivasankaranarayana Pillai independently
    proved that if the condition:
        r_k + q_k <= 2^k   (where 3^k = q_k * 2^k + r_k, 0 <= r_k < 2^k)
    holds for all k, then the exact formula for all natural numbers is:
        g(k) = 2^k + floor((3/2)^k) - 2

This module provides:
1. Complete Algebraic Formulations and Equivalence Proofs
2. Theoretical Analysis & Proof Obstructions (Mahler 1957, Baker's linear forms in logarithms)
3. Deep Connections to Related Problems:
   - Mahler's 3/2 Problem (Z-Numbers)
   - Equidistribution of {(3/2)^k} mod 1
   - The Collatz (3x+1) Conjecture (2-adic dynamics)
   - Pisot–Vijayaraghavan Numbers (Pisot vs Rational bases)
4. Formal Implications and Relationships Between These Problems
5. Probabilistic & Heuristic Justification (Borel–Cantelli)
6. Fast Exact-Integer Numerical Verification and Gap Analysis
7. Visualization Tools (Matplotlib)
"""

import math
from typing import List, Dict, Tuple, Optional
import matplotlib.pyplot as plt


# ==============================================================================
# 1. ALGEBRAIC REFORMULATION & EQUIVALENCES
# ==============================================================================
"""
Let 3^k = q_k * 2^k + r_k, where:
    q_k = floor((3/2)^k) = floor(3^k / 2^k)
    r_k = 3^k mod 2^k (with 0 <= r_k < 2^k)
    {(3/2)^k} = (3/2)^k - floor((3/2)^k) = r_k / 2^k

Dickson–Pillai Condition:
    2^k * {(3/2)^k} + floor((3/2)^k) <= 2^k
    <==> 2^k * (r_k / 2^k) + q_k <= 2^k
    <==> r_k + q_k <= 2^k

Equivalent Form 1 (Floor Bound):
    Since r_k = 3^k - q_k * 2^k:
    (3^k - q_k * 2^k) + q_k <= 2^k
    <==> 3^k - 2^k <= q_k * (2^k - 1)
    <==> q_k >= (3^k - 2^k) / (2^k - 1)
    <==> floor((3/2)^k) >= (3^k - 2^k) / (2^k - 1)

Equivalent Form 2 (Ceiling Form):
    Since (3^k - 2^k)/(2^k - 1) = (3^k - 1)/(2^k - 1) - 1:
    q_k + 1 >= (3^k - 1) / (2^k - 1)
    <==> ceil((3/2)^k) >= (3^k - 1) / (2^k - 1)

Equivalent Form 3 (Fractional Part / Diophantine Distance Bound):
    Substitute q_k = (3/2)^k - theta_k where theta_k = {(3/2)^k} in [0, 1):
    (3/2)^k - theta_k >= (3^k - 2^k) / (2^k - 1)
    <==> theta_k <= (3/2)^k - (3^k - 2^k) / (2^k - 1)
    <==> theta_k <= (4^k - 3^k) / (4^k - 2^k)
    <==> 1 - theta_k >= (3^k - 2^k) / (4^k - 2^k)
    As k -> inf, (3^k - 2^k)/(4^k - 2^k) = (3/4)^k * (1 - (2/3)^k) / (1 - 2^-k) ~ (3/4)^k.

    Conclusion: The Dickson–Pillai condition states that {(3/2)^k} cannot approach 1 faster than (3/4)^k:
        1 - {(3/2)^k} >= (3/4)^k * (1 - (2/3)^k) / (1 - 2^-k)
"""


# ==============================================================================
# 2. CONNECTIONS TO RELATED PROBLEMS & IMPLICATION ANALYSIS
# ==============================================================================
"""
A. MAHLER'S 3/2 PROBLEM (Z-NUMBERS)
   - Proposed by Kurt Mahler in 1968.
   - Question: Does there exist any non-zero real xi > 0 such that:
         0 <= { xi * (3/2)^n } < 1/2   for all integers n >= 0?
   - Such a number xi is called a Z-number.
   - Conjecture: No Z-numbers exist. (Status: OPEN).
   - Relationship: Both Mahler's problem and Dickson–Pillai probe whether the fractional
     parts of powers of (3/2) can remain trapped in restricted intervals of [0, 1).

B. EQUIDISTRIBUTION OF {(3/2)^k} MODULO 1
   - Question: Is the sequence {(3/2)^k}_{k=1}^inf uniformly distributed in [0, 1)?
   - Status: OPEN. Even the weaker question of whether {(3/2)^k} is DENSE in [0, 1] remains unsolved.
   - Relationship: If {(3/2)^k} is uniformly distributed, the probability of falling within
     distance (3/4)^k of 1 is exponentially small, making infinite exceptions impossible by Borel–Cantelli.

C. THE COLLATZ (3x + 1) PROBLEM & 2-ADIC DYNAMICS
   - The sequence r_k = 3^k mod 2^k represents the exact binary residue generated when the
     3x+1 map operates over k steps.
   - In the ring of 2-adic integers Z_2, the map x -> 3x is continuous and ergodic with respect
     to Haar measure.
   - The chaotic pseudo-randomness of r_k mod 2^k is the shared core reason why both
     the Collatz conjecture and the Dickson–Pillai condition resist elementary algebraic proofs.

D. PISOT–VIJAYARAGHAVAN NUMBERS (PISOT NUMBERS)
   - A Pisot number is a real algebraic integer alpha > 1 whose all other Galois conjugates
     lie strictly inside the unit disk |z| < 1.
   - For a Pisot number alpha, the distance to the nearest integer ||alpha^n|| -> 0 exponentially fast!
   - Charles Pisot (1938) proved that ONLY algebraic integers can have this property.
   - Because 3/2 is a rational number whose minimal polynomial is 2x - 3 = 0 (not monic over Z),
     3/2 is NOT an algebraic integer.
   - Consequently, (3/2)^k CANNOT approach integers exponentially fast in the manner of Pisot numbers.

--------------------------------------------------------------------------------
LOGICAL IMPLICATIONS: DO THESE PROBLEMS IMPLY EACH OTHER?
--------------------------------------------------------------------------------
1. Strong Equidistribution (with Quantitative Discrepancy Bounds) ==> Dickson–Pillai (for large k):
   - Qualitative equidistribution alone does not rule out finite isolated exceptions.
   - However, an effective discrepancy bound D_N = O(N^(-delta)) on {(3/2)^k} would immediately
     imply Mahler's bound 1 - {(3/2)^k} > (3/4)^k for all large k.

2. Universal Density / Equidistribution for all xi ==> No Z-Numbers Exist:
   - If {xi * (3/2)^n} is dense in [0, 1] for every xi > 0, then {xi * (3/2)^n} must enter [1/2, 1),
     which immediately proves that NO Z-numbers can exist.

3. Pisot's Theorem ==> (3/2)^k has no Pisot-like accumulation:
   - Pisot's structural theorem unconditionally guarantees that powers of 3/2 cannot converge
     exponentially to 0 mod 1.

4. The abc Conjecture ==> Effective Pillai / Diophantine Bounds:
   - The abc conjecture implies that |3^k - q * 2^k| cannot be abnormally small, providing
     strong polynomial-exponential lower bounds on Diophantine approximations.
"""


# ==============================================================================
# 3. THEORETICAL ANALYSIS & PROOF ATTEMPTS
# ==============================================================================
"""
WHY IS THE DICKSON-PILLAI CONDITION HARD TO PROVE?

1. Modular & Algebraic Obstruction:
   - In base 2, 3^k mod 2^k is the low k bits of 3^k.
   - The sequence r_k = 3^k mod 2^k behaves like a pseudo-random bit sequence.
   - Purely elementary algebraic identities cannot forbid r_k from being close to 2^k - 1.

2. Mahler's Theorem (1957) - Qualitative Resolution:
   - Kurt Mahler used Ridout's p-adic generalization of the Thue-Siegel-Roth theorem.
   - For any rational u/v > 1 and any epsilon > 0, the distance of (u/v)^k to the nearest integer
     ||(u/v)^k|| satisfies ||(u/v)^k|| > e^(-epsilon * k) for all sufficiently large k.
   - Taking u/v = 3/2 and epsilon = ln(4/3) - delta ~ 0.28768 - delta > 0:
     Mahler showed that 1 - {(3/2)^k} > (3/4)^k holds for all k >= K_0.
   - Thus, there are at most FINITELY MANY EXCEPTIONAL k!
   - OBSTRUCTION: Roth's theorem is INEFFECTIVE (it relies on proof by contradiction).
     It does not provide an explicit value for K_0.

3. The Conditional Proof via the abc Conjecture (Sinnou David & Michel Waldschmidt):
   - Consider the additive triple:
         A = r_k = 3^k - q_k * 2^k,   B = q_k * 2^k,   C = 3^k
     where gcd(A, B, C) = 1 and A + B = C.
   - The radical of ABC is:
         rad(ABC) = rad(r_k * q_k * 2^k * 3^k) <= 2 * 3 * rad(q_k) * rad(r_k) <= 6 * q_k * r_k
   - If a violation occurs, i.e., r_k + q_k > 2^k:
     Then r_k > 2^k - q_k, and since q_k ~ (3/2)^k, we have:
         rad(ABC) <= 6 * (3/2)^k * 2^k = 6 * 3^k = 6 * C
   - Under the abc Conjecture (C <= K(eps) * rad(ABC)^(1+eps)):
     By analyzing the prime factors and exponents of B = q_k * 2^k, Sinnou David (late 1990s)
     proved that the abc conjecture forces the Dickson–Pillai condition r_k + q_k <= 2^k
     to hold for all sufficiently large k.

4. Baker's Theory of Linear Forms in Logarithms - Effective Bounds:
   - Baker's method is effective (provides explicit K_0).
   - However, the best known effective lower bounds for ||(3/2)^k|| are of the form 2^(-c * k)
     with c near 1 (e.g. Beukers 1981, Dubickas, Bugeaud).
   - To prove our conjecture unconditionally, we need c <= log2(4/3) ~ 0.415037.
   - The gap between 0.415 and the current effective limit ~0.999 is a major open problem in
     transcendence theory / Diophantine approximation.

5. Critical Analysis of Recent Elementary Proof Claims (e.g. arXiv:2508.17950):
   - Recent preprints occasionally claim short elementary proofs by analyzing the continuous
     extensions F_2(x) = 2^x {(3/2)^x} + floor((3/2)^x).
   - For instance, arXiv:2508.17950v1 (2025) considers R_j = (3^j - 1)/(2^j - 1) and correctly
     bounds (3/2)^j < R_j < (3/2)^j + 1.
   - However, it commits a fatal algebraic fallacy by claiming:
         n + 1 <= floor(R_j) ==> n + 1 <= floor((3/2)^j)   [ERROR: dropped the +1!]
   - Retaining the correct bound floor(R_j) <= floor((3/2)^j) + 1 yields floor((3/2)^j) + 1 <= n + 1 <= floor((3/2)^j) + 1,
     which simply means n = floor((3/2)^j)—no contradiction exists.
   - Conclusion: Elementary algebraic bounding alone cannot resolve the pseudo-random bit behavior of 3^k mod 2^k.

6. Probabilistic / Heuristic (Borel-Cantelli):
   - If {(3/2)^k} is uniformly distributed in [0, 1), the probability of a violation at step k
     is P(theta_k > 1 - (3/4)^k) ~ (3/4)^k.
   - The expected number of violations for k >= K is sum_{k=K}^inf (3/4)^k = 4 * (3/4)^K.
   - For K = 471,600,000 (verified by Kubina & Wunderlich):
     The expected number of remaining exceptions is < 4 * (0.75)^(4.716 * 10^8) = 10^(-58,800,000).
   - It is mathematically almost certain that 0 exceptions exist.
"""


# ==============================================================================
# 4. EXACT INTEGER VERIFICATION & METRICS
# ==============================================================================

def check_dickson_pillai_condition(max_k: int = 5000, verbose: bool = True) -> Dict[str, object]:
    """
    Perform exact integer verification of the Dickson–Pillai condition:
        r_k + q_k <= 2^k for all 1 <= k <= max_k.
    Uses pure arbitrary-precision integers with log-ratio tracking to prevent overflow.

    Returns summary dictionary with statistics and worst-case margins.
    """
    min_safety_ratio = float('inf')
    worst_k = 1
    violations = []

    for k in range(1, max_k + 1):
        pow3 = 3**k
        pow2 = 1 << k
        q, r = divmod(pow3, pow2)
        
        # Check condition: r + q <= 2^k
        diff = pow2 - (r + q)
        if diff < 0:
            violations.append((k, q, r, diff))
            continue
        
        # For k >= 2, evaluate safety ratio: (1 - frac) / danger_bound
        # 1 - frac = (2^k - r) / 2^k
        # danger_bound = (3^k - 2^k) / (4^k - 2^k) = (3^k - 2^k) / (2^k * (2^k - 1))
        # safety_ratio = (2^k - r) * (2^k - 1) / (3^k - 2^k)
        if k >= 2:
            num = (pow2 - r) * (pow2 - 1)
            den = pow3 - pow2
            
            log2_num = num.bit_length() - 1 + math.log2(num >> max(0, num.bit_length() - 53)) - (min(53, num.bit_length()) - 1) if num > 0 else float('-inf')
            log2_den = den.bit_length() - 1 + math.log2(den >> max(0, den.bit_length() - 53)) - (min(53, den.bit_length()) - 1)
            log2_ratio = log2_num - log2_den
            
            # For small k (k < 50), compute exact float ratio
            if k < 50:
                ratio_float = num / den
                if ratio_float < min_safety_ratio:
                    min_safety_ratio = ratio_float
                    worst_k = k
            elif log2_ratio < 0:
                # Violation! ratio < 1.0
                min_safety_ratio = math.pow(2, log2_ratio)
                worst_k = k

    is_valid = (len(violations) == 0)
    
    if verbose:
        print(f"=== Dickson-Pillai Condition Verification (k = 1 to {max_k}) ===")
        print(f"All k valid: {is_valid}")
        print(f"Violations found: {len(violations)}")
        print(f"Tightest safety ratio: {min_safety_ratio:.6f} (at k = {worst_k})")
        print("Note: safety_ratio >= 1.0 means condition is strictly satisfied.")
        print("=================================================================\n")
        
    return {
        "is_valid": is_valid,
        "max_k": max_k,
        "violations": violations,
        "tightest_safety_ratio": min_safety_ratio,
        "worst_k": worst_k
    }


def get_detailed_table(start_k: int = 2, end_k: int = 25) -> List[Dict[str, object]]:
    """
    Generate exact values of q, r, diff, fractional part, and threshold for small k.
    """
    rows = []
    for k in range(start_k, end_k + 1):
        pow3 = 3**k
        pow2 = 1 << k
        q, r = divmod(pow3, pow2)
        diff = pow2 - (r + q)
        frac = r / pow2
        danger_threshold = (pow3 - pow2) / (pow2 * (pow2 - 1))
        safety_ratio = ((pow2 - r) * (pow2 - 1)) / (pow3 - pow2)
        
        rows.append({
            "k": k,
            "q": q,
            "r": r,
            "diff": diff,
            "frac": frac,
            "danger_threshold": 1.0 - danger_threshold,
            "safety_ratio": safety_ratio
        })
    return rows


# ==============================================================================
# 5. PLOTTING & VISUALIZATION
# ==============================================================================

def plot_dickson_pillai_analysis(max_k: int = 500, save_path: Optional[str] = None):
    """
    Plot 3 analytical views:
    1. Log2 of Safety Ratio log2((1 - { (3/2)^k }) / (danger threshold)) over k.
    2. The fractional part { (3/2)^k } vs the danger boundary 1 - (3/4)^k.
    3. Distribution (histogram) of { (3/2)^k } modulo 1 demonstrating uniform spread.
    """
    k_vals = list(range(2, max_k + 1))
    log2_ratios = []
    fracs = []
    danger_bounds = []
    
    for k in k_vals:
        pow3 = 3**k
        pow2 = 1 << k
        q, r = divmod(pow3, pow2)
        
        frac = (r >> max(0, k - 53)) / (pow2 >> max(0, k - 53))
        fracs.append(frac)
        
        # Danger bound for fractional part: 1 - (3/4)^k
        danger = 1.0 - math.pow(0.75, k)
        danger_bounds.append(danger)
        
        # Safety ratio log2
        num = (pow2 - r) * (pow2 - 1)
        den = pow3 - pow2
        log2_num = num.bit_length() - 1 + math.log2(num >> max(0, num.bit_length() - 53)) - (min(53, num.bit_length()) - 1) if num > 0 else 0
        log2_den = den.bit_length() - 1 + math.log2(den >> max(0, den.bit_length() - 53)) - (min(53, den.bit_length()) - 1)
        log2_ratios.append(log2_num - log2_den)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Plot 1: Safety Ratio in log2 scale
    axes[0].plot(k_vals, log2_ratios, color='darkblue', lw=1.2)
    axes[0].axhline(y=0.0, color='red', linestyle='--', label='Violation Threshold (log2 Ratio = 0)')
    axes[0].set_title(r"$\log_2(\text{Safety Ratio}) \sim k \log_2(4/3)$")
    axes[0].set_xlabel("k")
    axes[0].set_ylabel(r"$\log_2(\text{Safety Ratio})$")
    axes[0].grid(True, which="both", ls="--", alpha=0.5)
    axes[0].legend()
    
    # Plot 2: Fractional Part vs Danger Curve
    axes[1].scatter(k_vals, fracs, s=6, color='purple', alpha=0.6, label=r"$\{(3/2)^k\}$")
    axes[1].plot(k_vals, danger_bounds, color='red', lw=1.5, label=r"Danger Zone: $1 - (3/4)^k$")
    axes[1].set_title(r"$\{(3/2)^k\}$ vs Danger Curve")
    axes[1].set_xlabel("k")
    axes[1].set_ylabel("Value in [0, 1)")
    axes[1].grid(True, ls="--", alpha=0.5)
    axes[1].legend()

    # Plot 3: Histogram of Fractional Parts
    axes[2].hist(fracs, bins=30, color='teal', edgecolor='black', alpha=0.7, density=True)
    axes[2].axhline(y=1.0, color='crimson', linestyle='--', label='Uniform Dist PDF')
    axes[2].set_title(r"Empirical Distribution of $\{(3/2)^k\}$ mod 1")
    axes[2].set_xlabel(r"$\{(3/2)^k\}$")
    axes[2].set_ylabel("Density")
    axes[2].grid(True, ls="--", alpha=0.5)
    axes[2].legend()
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to {save_path}")
    else:
        plt.show()


if __name__ == "__main__":
    # 1. Run exact check up to k = 5000
    check_dickson_pillai_condition(max_k=5000, verbose=True)
    
    # 2. Print small k exact breakdown
    print("--- Detailed Values for k = 2 to 15 ---")
    table = get_detailed_table(2, 15)
    print(f"{'k':>2} | {'q':>8} | {'r':>8} | {'diff = 2^k - (r+q)':>18} | {'{(3/2)^k}':>10} | {'Danger Bound':>12} | {'Safety Ratio':>12}")
    print("-" * 85)
    for row in table:
        print(f"{row['k']:2d} | {row['q']:8d} | {row['r']:8d} | {row['diff']:18d} | {row['frac']:10.4f} | {row['danger_threshold']:12.4f} | {row['safety_ratio']:12.4f}")
    
    # 3. Plot overview (saving to plots/ directory)
    try:
        import os
        plots_dir = os.path.join(os.path.dirname(__file__), "plots")
        os.makedirs(plots_dir, exist_ok=True)
        plot_path = os.path.join(plots_dir, "dickson_pillai_analysis.png")
        plot_dickson_pillai_analysis(max_k=500, save_path=plot_path)
    except Exception as e:
        print(f"Plotting note: {e}")
